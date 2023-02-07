! <compile=optimized>
#include "copyright.h"
#include "../include/dprec.fh"
#include "../include/assert.fh"
#include "ncsu-config.h"

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ main driver routine to compute energies and forces
subroutine force(xx,ix,ih,ipairs,x,f,ener,vir, &
      fs, rborn, reff, onereff, qsetup,do_list_update,nstep)

#ifndef DISABLE_NCSU
   use ncsu_sander_hooks, only : ncsu_on_force => on_force
#  if !defined(LES) && defined(MPI)
   use remd, only : rem, mdloop
#  endif
#  ifdef MPI
   use ncsu_sander_proxy, only : ncsu_remember_initremd => remember_initremd
#  endif /* MPI */
#endif /* DISABLE_NCSU */
   use file_io_dat
#ifdef LES
   use genbornles
   use pimd_vars, only : nrg_all
#  ifdef MPI
   use les_data, only : elesa, elesb, elesd, elesp
#  endif /* MPI */
#else
   use genborn
   use qmmm_module, only : qmmm_struct, qm2_struct
#endif /* LES */
#ifdef MPI
   use neb_vars, only: ineb, neb_force
   use full_pimd_vars, only: totener
   use softcore, only: sc_ener
#endif /* MPI */
#if defined(LES) && defined(MPI)
   use evb_data, only: nrg_frc
   use pimd_vars, only: equal_part
   use miller, only : dlnQ_dl
   use remd, only : rem ! wasn't used for LES above
#endif /* LES && MPI */
   use poisson_boltzmann, only : pb_force
   use dispersion_cavity, only : npopt, np_force
   use pbtimer_module, only : pbtimer_init, pbtimer_summary
#ifdef RISMSANDER
   use sander_rism_interface, only: rismprm, rism_force
#endif /* RISMSANDER */
#ifdef APBS
   use apbs
#endif /* APBS */
   use trace
   use stack
   use pimd_vars, only: ipimd,nbead,bnd_vir,Epot_spring,Epot_deriv,real_mass,itimass
   use qmmm_module, only : qmmm_nml
   use constants, only : zero, one
   use relax_mat
   use ew_recip
   use parms, only:cn1,cn2,cn6,asol,bsol,pk,rk,tk,numbnd,numang,nptra,nphb,nimprp,&
                   cn3,cn4,cn5!mjhsieh: for another vdwmodel
#ifdef PUPIL_SUPPORT
   use nblist, only:nonbond_list,a,b,c,alpha,beta,gamma,ucell
#else
   use nblist, only:nonbond_list,a,b,c,alpha,beta,gamma
#endif /*PUPIL_SUPPORT*/
#ifdef DSSP
   use dssp, only: fdssp, edssp, idssp
#endif /* DSSP */


   use amoeba_interface,only: AM_VAL_eval,AM_NonBond_eval
   use amoeba_mdin, only : iamoeba,am_nbead

   use amd_mod
   use scaledMD_mod
   use nbips,only: ips,eexips
   use emap,only:temap,emapforce

#ifdef PUPIL_SUPPORT
   use pupildata
#endif /*PUPIL_SUPPORT*/

   use linear_response, only: ilrt, ee_linear_response, energy_m0, energy_w0, &
        energy_vdw0, cn1_lrt, cn2_lrt, crg_m0, crg_w0, do_lrt, f_scratch, &
        lrt_solute_sasa

   use cns_xref
   use xray_interface_module, only: xray_get_derivative, xray_active

!CHARMM Force Field Support
   use charmm_mod, only : charmm_active, charmm_calc_impropers, &
                          charmm_calc_cmap, charmm_calc_urey_bradley,&
                          charmm_dump_gold, &
                          do_charmm_dump_gold
   use ff11_mod, only : cmap_active, calc_cmap
   use decomp, only: &
#ifdef MPI
      collect_dec2,  &
#endif
      init_dec
   use state
   use crg_reloc, only: ifcr, cr_reassign_charge, cr_calc_force
   use sebomd_module, only : sebomd_obj, sebomd_save_forces
   use abfqmmm_module
   use les_data, only : temp0les

   implicit none
  
   integer, intent(in) :: nstep

#ifdef PUPIL_SUPPORT
   character(kind=1,len=5) :: routine="force"
#endif
#if defined(PUPIL_SUPPORT) || defined(MPI)
   integer ierr
#endif
   integer   ipairs(*)
   _REAL_ xx(*)
   integer   ix(*)
   character(len=4) ih(*)
   _REAL_ fs(*),rborn(*),reff(*),dvdl
   _REAL_, intent(out) :: onereff(*)
#include "def_time.h"
#include "ew_frc.h"
#include "ew_cntrl.h"
#include "extra_pts.h"
#include "parallel.h"

#ifdef MPI
#  include "ew_parallel.h"
#  ifdef MPI_DOUBLE_PRECISION
#     undef MPI_DOUBLE_PRECISION
#  endif
   include 'mpif.h'
   integer gb_rad_mpistart, j3, j, i3
   _REAL_ :: temp_amd_totdih
#  ifdef LES
   _REAL_ :: vel0_nrg_sum
#  endif /* LES */
#  ifdef CRAY_PVP
#     define MPI_DOUBLE_PRECISION MPI_REAL8
#  endif
#endif

   logical belly
#include "../include/md.h"
#include "../pbsa/pb_md.h"
#include "box.h"
#include "nmr.h"
#include "../include/memory.h"
#include "extra.h"
#include "tgtmd.h"
#include "multitmd.h"
#include "flocntrl.h"
   integer istart,iend
   _REAL_ evdwex, eelex
   _REAL_ enemap

   logical, intent(inout) :: qsetup
   logical, intent(out) :: do_list_update

   _REAL_  enmr(6),devdis(4),devang(4),devtor(4),devpln(4),devplpt(4),devgendis(4),entr,ecap
   _REAL_  x(*),f(*),vir(4)
   type(state_rec)  ener

   !Local
   _REAL_                     :: ene(30)    !Used locally ONLY
   type(potential_energy_rec) :: pot        !Used locally ONLY

#if defined(LES) && defined(MPI)
   _REAL_  :: nrg_bead(nbead)
#endif /* LES && MPI */

#ifndef LES
   _REAL_ escf
#endif /* LES */

   integer i
   _REAL_  virvsene,eelt,epol,esurf,edisp
#ifdef APBS
   _REAL_ enpol
#endif /* APBS */
#ifdef RISMSANDER
   _REAL_ erism
#endif /*RISMSANDER*/

   _REAL_ ect ! charge transfer

   _REAL_  epolar,aveper,aveind,avetot,emtot,dipiter,dipole_temp
   integer, save :: newbalance
   
!AMD variables
   _REAL_ amd_totdih

!  SEBOMD Gradient test variables
   _REAL_ :: em2, em1, e0, ep1, ep2
   _REAL_ :: fx1, fx2, fxx, xdx, dx

   ect = 0.0

   call trace_enter( 'force' )

   call timer_start(TIME_FORCE)
   if ( idecomp /= 0 .and. icfe == 0) call init_dec
   ene(:) = ZERO 
   call zero_pot_energy(pot)

   belly = ibelly > 0

#ifdef MPI
   if (mpi_orig) then

      !     Check to see if we are done yet in mpi_orig case (tec3).
      !     This is done by monitoring the status of an integer notdone.
      !     If notdone .eq. 1 then we keep going.  notdone is set to zero
      !     when we no longer want to call force().  This perhaps is not the
      !     most efficient means to implement this check...

      call mpi_bcast(notdone,1,mpi_integer,0,commsander,ierr)
      if (notdone /= 1) return

      !       Send copies of xyz coords, setbox common block, vir array
      !       and NTNB value to all nodes from master with a broadcast.

      if (numtasks > 1) then

         call mpi_bcast(box,BC_BOXR,MPI_DOUBLE_PRECISION,0, &
               commsander,ierr)
         call mpi_bcast(ntb,BC_BOXI,mpi_integer,0,commsander,ierr)
         call mpi_bcast(vir,3,MPI_DOUBLE_PRECISION,0, &
               commsander,ierr)
         call mpi_bcast(xx(lcrd),3*natom,MPI_DOUBLE_PRECISION, &
               0,commsander,ierr)
         call mpi_bcast(ntnb,1,mpi_integer,0,commsander,ierr)
         if (iabs(ntb) >= 2) then
            call mpi_bcast(xx(l45),3*natom,MPI_DOUBLE_PRECISION, &
                  0,commsander,ierr)
         end if
      end if
   end if

   istart = iparpt(mytaskid) + 1
   iend = iparpt(mytaskid+1)
#else
   istart = 1
   iend = natom
#endif /* MPI */

   !-----------------------------------
   !   QMMM Variable QM Solvent Scheme
   !-----------------------------------
   if (qmmm_nml%ifqnt .and. qmmm_nml%vsolv > 0) then
     call timer_start(TIME_QMMM)
     call timer_start(TIME_QMMMVARIABLESOLVCALC)
     ! For the moment this is NOT parallel - all threads need to call this.
     call qmmm_vsolv_update(nstep, natom, qsetup, xx(l15), x,  &
                            nbonh, nbona, ntheth, ntheta, nphih, nphia, &
                            ix(iibh),ix(ijbh),ix(iicbh), &
                            ix(iiba),ix(ijba),ix(iicba), &
                            ix(i24),ix(i26),ix(i28),ix(i30), &
                            ix(i32),ix(i34),ix(i36),ix(i38), &
                            ix(i40),ix(i42),ix(i44),ix(i46), &
                            ix(i48),ix(i50),ix(i52),ix(i54), &
                            ix(i56),ix(i58), ix(ibellygp))
     call timer_stop(TIME_QMMMVARIABLESOLVCALC)
     call timer_stop(TIME_QMMM)
   end if 
   !-----------------------------------
   ! END QMMM Variable QM Water Scheme
   !-----------------------------------

   if(iamoeba.eq.1) then
      REQUIRE(am_nbead.eq.ncopy)
   end if

   !     ----- ZERO OUT THE ENERGIES AND FORCES -----
   enoe = 0.d0
   aveper=0.d0
   aveind=0.d0
   avetot=0.d0
   dipiter=0.d0
   dvdl=0.d0
   dipole_temp=0.d0
   enmr(1:6) = 0.d0
   vir(1:4) = 0.d0
   virvsene = 0.d0
   f(1:3*natom+iscale) = 0.d0
#ifdef LES
   if( ipimd>0) nrg_all(1:nbead)=0.d0
#endif

   if (sebomd_obj%do_sebomd) then
     ! step 0: initialize temporary array for saving restraint
     !         (required since SEBOMD forces are stored in a different array)
     call sebomd_save_forces(0,natom,f,xx(gradsebomd),xx(grad1tmp),xx(grad2tmp))
   end if

#ifndef PUPIL_SUPPORT
   if( igb == 0 .and. ipb == 0 .and. iyammp == 0 ) then

      ! (for GB: do all nonbondeds together below)

      call timer_start(TIME_NONBON)
      call timer_start(TIME_LIST)

      if(abfqmmm_param%abfqmmm == 1) qsetup = .true.

      call nonbond_list(x,ix(i04),ix(i06),ix(i08),ix(i10), &
               ntypes,natom/am_nbead,xx,ix,ipairs,ntnb, &
               ix(ibellygp),belly,newbalance, &
               qsetup, &
               do_list_update)
      call timer_stop(TIME_LIST)
      call timer_stop(TIME_NONBON)
   end if
   ! charge reassign here !
   if ( ifcr /= 0 ) then
      call cr_reassign_charge( x, f, pot%ct, xx(l15), natom )
   end if
#endif
   
   if (sebomd_obj%do_sebomd) then
     ! step 1/2 to save forces from ncsu_on_force
     call sebomd_save_forces(1,natom,f,xx(gradsebomd),xx(grad1tmp),xx(grad2tmp))
   end if

#ifndef DISABLE_NCSU
#  ifdef MPI
   call ncsu_remember_initremd(rem.gt.0.and.mdloop.eq.0)
#  endif
   call ncsu_on_force(x, f)
#endif

   if (sebomd_obj%do_sebomd) then
     ! step 2/2 to save forces from ncsu_on_force
     call sebomd_save_forces(2,natom,f,xx(gradsebomd),xx(grad1tmp),xx(grad2tmp))
   end if

   if (sebomd_obj%do_sebomd) then
     ! variables initialization
    
     sebomd_obj%diverror = 0
     if (ntp > 0 .or. ntb == 2) then
        do i=1,3
           sebomd_obj%dbox(i) = box(i)
        enddo
     endif

     ! first call to SEBOMD to initialize the parameters
400  if (sebomd_obj%idcflag == 0) then

       if (ntb == 0) then
          sebomd_obj%dbox(1) = 0.d0
          sebomd_obj%dbox(2) = 0.d0
          sebomd_obj%dbox(3) = 0.d0
       else
          sebomd_obj%dbox(1) = box(1)
          sebomd_obj%dbox(2) = box(2)
          sebomd_obj%dbox(3) = box(3)
!         ! force wrapping (needed by sebomd)
!         call wrap_molecules(nspm,ix(i70),x)
       endif
      
       if (sebomd_obj%chtype > 0) then
!         call sebomd_readatchg(xx(divchg))
       else
          do i = 1, natom
             xx(divchg + i - 1) = xx(l15 + i - 1)/18.2223d0
          enddo
       endif
      
       call sebomd_energy(sebomd_obj%idcflag,sebomd_obj%iflagch, &
                          ntb,natom,vir,sebomd_obj%esebomd, &
                          sebomd_obj%dbox, &
                          x,xx(gradsebomd),nres,ix(i02),xx(divchg),ih(m06), &
                          nspm,ix(i70), &
                          sebomd_obj%pdmx)
     endif ! idcflag = 0

     ! Checking that everything is ok, exiting if not
     if (sebomd_obj%idcflag == -1) then
        write (6,*) 'Error during SEBOMD calculation'
        call mexit(6, 1)
     endif
    
     sebomd_obj%idcflag = 1
    
!    if (ntb.ne.0) then
!       ! force wrapping (needed by sebomd)
!       call wrap_molecules(nspm,ix(i70),x)
!    endif

     ! second call to SEBOMD to calculate the potential energy
     call sebomd_energy(sebomd_obj%idcflag,sebomd_obj%iflagch, &
                        ntb,natom,vir,sebomd_obj%esebomd, &
                        sebomd_obj%dbox, &
                        x,xx(gradsebomd),nres,ix(i02),xx(divchg),ih(m06), &
                        nspm,ix(i70), &
                        sebomd_obj%pdmx)
!    pot%scf = escf
    
     ! Checking that everything is ok, initializing again if not
     if (sebomd_obj%idcflag == -1) then
        if (sebomd_obj%diverror == 0) then
           sebomd_obj%idcflag = 0
           sebomd_obj%diverror = 1
           write (6,*) 'Error during SEBOMD calculation, restarting with a new ', &
                'initial density matrix'
           goto 400
        else
           write (6,*) 'Error during SEBOMD calculation'
           call mexit(6, 1)
        endif
     endif
     ! SEBOMD is currently computing gradient (not forces)
     do i =1,3*natom
        xx(gradsebomd+i-1) = -xx(gradsebomd+i-1)
     enddo

! --------------------------------------------------------------------
     if (sebomd_obj%debugforces /= 0) then
       dx = 1.0d-5
       do i=1,3*natom
         ! central energy E_SEBOMD(x)
         xdx = x(i)
         call sebomd_save_forces(0,natom,f,xx(gradsebomd),xx(grad1tmp),xx(grad2tmp))
         call sebomd_energy(sebomd_obj%idcflag,sebomd_obj%iflagch, &
                            ntb,natom,vir,sebomd_obj%esebomd, &
                            sebomd_obj%dbox, &
                            x,xx(gradsebomd),nres,ix(i02),xx(divchg),ih(m06), &
                            nspm,ix(i70), &
                            sebomd_obj%pdmx)
         e0 = sebomd_obj%esebomd
         fxx=-xx(gradsebomd+i-1)

         ! E_SEBOMD(x-\delta x)
         x(i) = xdx-dx
         call sebomd_save_forces(0,natom,f,xx(gradsebomd),xx(grad1tmp),xx(grad2tmp))
         call sebomd_energy(sebomd_obj%idcflag,sebomd_obj%iflagch, &
                            ntb,natom,vir,sebomd_obj%esebomd, &
                            sebomd_obj%dbox, &
                            x,xx(gradsebomd),nres,ix(i02),xx(divchg),ih(m06), &
                            nspm,ix(i70), &
                            sebomd_obj%pdmx)
         em1 = sebomd_obj%esebomd

         ! E_SEBOMD(x-2\delta x)
         x(i) = xdx-dx-dx
         call sebomd_save_forces(0,natom,f,xx(gradsebomd),xx(grad1tmp),xx(grad2tmp))
         call sebomd_energy(sebomd_obj%idcflag,sebomd_obj%iflagch, &
                            ntb,natom,vir,sebomd_obj%esebomd, &
                            sebomd_obj%dbox, &
                            x,xx(gradsebomd),nres,ix(i02),xx(divchg),ih(m06), &
                            nspm,ix(i70), &
                            sebomd_obj%pdmx)
         em2 = sebomd_obj%esebomd

         ! E_SEBOMD(x+\delta x)
         x(i) = xdx+dx
         call sebomd_save_forces(0,natom,f,xx(gradsebomd),xx(grad1tmp),xx(grad2tmp))
         call sebomd_energy(sebomd_obj%idcflag,sebomd_obj%iflagch, &
                            ntb,natom,vir,sebomd_obj%esebomd, &
                            sebomd_obj%dbox, &
                            x,xx(gradsebomd),nres,ix(i02),xx(divchg),ih(m06), &
                            nspm,ix(i70), &
                            sebomd_obj%pdmx)
         ep1 = sebomd_obj%esebomd

         ! E_SEBOMD(x+2\delta x)
         x(i) = xdx+dx+dx
         call sebomd_save_forces(0,natom,f,xx(gradsebomd),xx(grad1tmp),xx(grad2tmp))
         call sebomd_energy(sebomd_obj%idcflag,sebomd_obj%iflagch, &
                            ntb,natom,vir,sebomd_obj%esebomd, &
                            sebomd_obj%dbox, &
                            x,xx(gradsebomd),nres,ix(i02),xx(divchg),ih(m06), &
                            nspm,ix(i70), &
                            sebomd_obj%pdmx)
         ep2 = sebomd_obj%esebomd

         fx1=-(ep1-em1)/(2.0d0*dx)
         fx2=(1.0d0/ 12.0d0*(ep2-em2) + 2.0d0/  3.0d0*(em1-ep1))/dx

         write(0,'("gradients ",2i4,3f25.16)') (i-1)/3+1,i,fxx,fx2,fx1
         write(0,'("delta grad",2i4,2f25.16)') (i-1)/3+1,i,fxx-fx2,fxx-fx1
         ! restore the current coordinate
         x(i) = xdx
       enddo
     endif
! -------------------------------------------------------------------
   endif ! do_sebomd

   ! ----------------------------------------------------------------
   ! Do weight changes, if requested
   ! ----------------------------------------------------------------

   if (nmropt > 0) &
         call nmrcal(x,f,ih(m04),ih(m02),ix(i02),xx(lwinv),enmr,devdis, &
                     devang,devtor,devplpt,devpln,devgendis,temp0,tautp,cut,&
                     xx(lnmr01),ix(inmr02),xx(l95),31,6,rk,tk,pk,cn1, &
                     cn2,asol,bsol,xx(l15),numbnd,numang,nptra-nimprp, &
                     nimprp,nphb,natom,natom,ntypes,nres, &
                     rad,wel,radhb,welhb,rwell,tgtrmsd,temp0les,-1,'WEIT')
   ! Updated 9/2007 by Matthew Seetin to enable plane-point and plane-plane restraints

   epolar = 0.d0

   ! -----------------------------------------------------------------
   ! EGB if igb>0 and /=6 then we need to calculate the GB radii for this
   ! structure
   ! -----------------------------------------------------------------

   ! If we are doing qm_gb=2 then we need to calculate the GB radii
   ! before calling qm_mm
   if( igb > 0 .and. igb /= 6 .and. igb /= 10 .and. ipb == 0 .and. &
      ( irespa < 2 .or. mod(irespa,nrespai) == 0) ) then
#ifdef MPI
      gb_rad_mpistart = mytaskid+1
#endif
      call timer_start(TIME_EGB)
      call timer_start(TIME_GBRAD1)
      !If qmmm and then this will calculate us the radii for the
      !link atoms and not the mm link pair atoms. The initial radii used are those
      !of the mm link pair's atom type though.
      !The MMlink pair atoms must be the coordinates of the link atoms here.
      if (qmmm_nml%ifqnt) call adj_mm_link_pair_crd(x)

      call egb_calc_radii(igb,natom,x,fs,reff, &
                     onereff,fsmax,rgbmax, rborn, offset, &
                     rbornstat,xx(l188),xx(l189),         &
                     xx(l186),xx(l187), gbneckscale, ncopy, rdt, &
#ifdef LES
                     gbalpha, gbbeta, gbgamma & !Hai Nguyen: keep for sander.LES
#else
                     xx(l2402),xx(l2403),xx(l2404)  & ! Hai Nguyen: gbvalpha,gbvbeta,gbvgamma 
#endif
#ifdef MPI
                     ,gb_rad_mpistart &
#endif
                       )
      if (qmmm_nml%ifqnt) call rst_mm_link_pair_crd(x)

      call timer_stop(TIME_GBRAD1)
      call timer_stop(TIME_EGB)

   end if
 
   ! QM/MM Contributions are now calculated before the NON-Bond info.
   ! ----------------------------------------------------------------
   ! Calculate the qm/mm contributions
   ! ----------------------------------------------------------------

   if(qmmm_nml%ifqnt) then

      ! If we are doing periodic boundaries with QM/MM PME then we need to
      ! do the PME calculation twice. First here to get the potential at
      ! each QM atom due to the PME and then again after all the QM is done
      ! to get the MM-MM potential and all of the gradients.

      if(qmmm_nml%qm_pme) then
         ! Ewald force will put the potential into the qm_ewald%mmpot array.
         call timer_start(TIME_EWALD)
         call ewald_force(x,natom,ix(i04),ix(i06), &
                          xx(l15),cn1,cn2,cn6,eelt,epolar, &
                          f,xx,ix,ipairs,xx(l45),virvsene, xx(lpol), &
#ifdef HAS_10_12
                          xx(lpol2),.true., cn3, cn4, cn5, asol, bsol)
#else
                          xx(lpol2),.true., cn3, cn4, cn5 )
#endif
         call timer_stop(TIME_EWALD)
      endif

      call timer_start(TIME_QMMM)

#ifndef LES
        !========================================================
        !                      REGULAR QMMM
        !========================================================
        !if(qmmm_nml%idc>0)then
        !   call qm_div(x, ix, f, escf, ih(m06))
        !else

         
      call qm_mm(x, natom,qmmm_struct%scaled_mm_charges, &
                 f,escf,periodic,reff,onereff, &
                intdiel,extdiel,Arad, cut,qm2_struct%scf_mchg,ntypes, &
                 ih(m04), ih(m06),xx(lmass), ix(i04),nstep)
        pot%scf = escf
#endif

        !========================================================
        !                  END REGULAR QMMM
        !========================================================
      call timer_stop(TIME_QMMM)
   end if !if(qmmm_nml%ifqnt)

   !---------------------------------------------------------------
   !END qm/mm contributions
   !---------------------------------------------------------------

#ifdef PUPIL_SUPPORT

   !*****************************************************
   !     Getting the Quantum forces with PUPIL package
   !*****************************************************

!  Reconstruct the simulation cell if there is any change
!  call inipupcell(natms,qcell,cell,xxx,yyy,zzz)
   do iPup=1,3    !vector loop
     do jPup=1,3  !Component loop
       qcell((iPup-1)*3+jPup) = ucell(jPup,iPup)
     enddo
   enddo
!  minimum point of the box ..... we assume (0,0,0) ????
   qcell(10) = 0.0d0
   qcell(11) = 0.0d0
   qcell(12) = 0.0d0

!  temporary vector to wrap the real coordinates to pass through
!  PUPIL interface.
!  This stack will be deallocted after recover all PUPIL forces

   call get_stack(l_puptmp,3*natom,routine)
   if(.not. rstack_ok)then
     deallocate(r_stack)
     allocate(r_stack(1:lastrst),stat=alloc_ier)
     call reassign_rstack(routine)
   endif
   REQUIRE(rstack_ok)
   do iPup=1,3*natom
     r_stack(l_puptmp + iPup - 1) = x(iPup)
   end do

   if(ntb > 0) then
     call wrap_molecules(nspm,ix(i70),r_stack(l_puptmp))
     if(ifbox == 2) call wrap_to(nspm,ix(i70),r_stack(l_puptmp),box)
   end if

!  Preparing the coordinates, velocity and classic forces
!  to get quantum force
   do iPup=1,natom
     bs1 = (iPup-1)*9
     bs2 = (iPup-1)*3
     do jPup=1,3
       qcdata(bs1+jPup) = r_stack(l_puptmp + bs2 + jPup - 1)
!        qcdata(bs1+jPup) = x(bs2+jPup)
!          write(6,*) 'Coordinate.',iPup,'==>',realStack(lcrd+bs2+jPup-1),x(bs2+jPup)
!          write(6,*) 'Velocity...',iPup,'==>',realStack(lvel+bs2+jPup-1)
       qcdata(bs1+3+jPup) = realStack(lvel+bs2+jPup-1)
       qcdata(bs1+6+jPup) = f(bs2+jPup)
     enddo
   enddo

!  We are going to use the qmmm_nml and qmmm_struct variables to skip quantum atoms
!  in the force calculation
   qmmm_nml%ifqnt = .true.
   if(pupStep .EQ. 0) then
!!!  To keep initial values from the MD step 1
     ierr = 0
     allocate ( pupnb14(numnb14*3),stat=ierr)
     REQUIRE(ierr == 0)

     allocate ( pupbonh(nbonh*3),stat=ierr)
     REQUIRE(ierr == 0)
     allocate ( pupbona(nbona*3),stat=ierr)
     REQUIRE(ierr == 0)
     allocate ( puptheth(ntheth*4),stat=ierr)
     REQUIRE(ierr == 0)
     allocate ( puptheta(ntheta*4),stat=ierr)
     REQUIRE(ierr == 0)
     allocate ( pupphih(nphih*5),stat=ierr)
     REQUIRE(ierr == 0)
     allocate ( pupphia(nphia*5),stat=ierr)
     REQUIRE(ierr == 0)
     pupnbonh   = nbonh
     pupnbona   = nbona
     pupntheth  = ntheth
     pupntheta  = ntheta
     pupnphih   = nphih
     pupnphia   = nphia
     pupnumnb14 = numnb14
     call copy_14nb(ix(inb_14),pupnb14,numnb14)
     do iPup = 1,nbonh
       bs1 = iPup-1
       bs2 = bs1*3
       pupbonh(bs2+1)  = ix(iibh +bs1)
       pupbonh(bs2+2)  = ix(ijbh +bs1)
       pupbonh(bs2+3)  = ix(iicbh+bs1)
     enddo
     do iPup = 1,nbona
       bs1 = iPup-1
       bs2 = bs1*3
       pupbona(bs2+1)  = ix(iiba +bs1)
       pupbona(bs2+2)  = ix(ijba +bs1)
       pupbona(bs2+3)  = ix(iicba+bs1)
     enddo
     do iPup = 1,ntheth
       bs1 = iPup-1
       bs2 = bs1*4
       puptheth(bs2+1) = ix(i24  +bs1)
       puptheth(bs2+2) = ix(i26  +bs1)
       puptheth(bs2+3) = ix(i28  +bs1)
       puptheth(bs2+4) = ix(i30  +bs1)
     enddo
     do iPup = 1,ntheta
       bs1 = iPup-1
       bs2 = bs1*4
       puptheta(bs2+1) = ix(i32  +bs1)
       puptheta(bs2+2) = ix(i34  +bs1)
       puptheta(bs2+3) = ix(i36  +bs1)
       puptheta(bs2+4) = ix(i38  +bs1)
     enddo
     do iPup = 1,nphih
       bs1 = iPup-1
       bs2 = bs1*5
       pupphih(bs2+1)  = ix(i40  +bs1)
       pupphih(bs2+2)  = ix(i42  +bs1)
       pupphih(bs2+3)  = ix(i44  +bs1)
       pupphih(bs2+4)  = ix(i46  +bs1)
       pupphih(bs2+5)  = ix(i48  +bs1)
     enddo
     do iPup = 1,nphia
       bs1 = iPup-1
       bs2 = bs1*5
       pupphia(bs2+1)  = ix(i50  +bs1)
       pupphia(bs2+2)  = ix(i52  +bs1)
       pupphia(bs2+3)  = ix(i54  +bs1)
       pupphia(bs2+4)  = ix(i56  +bs1)
       pupphia(bs2+5)  = ix(i58  +bs1)
     enddo
   endif

!  Getting the quantum forces for a specific quantum domain
   pupStep  = pupStep + 1
   puperror = 0
   pupLevelData = 3
   call getquantumforces(natom,pupLevelData,pupStep,puperror,qcdata,qcell)
   if (puperror.ne.0) then
     write (6,*) 'Fatal error: Could not obtain quantum forces!'
     call mexit(6,1)
   endif
 ! Quantum energy treatment....
   pot%scf = qmEnergy

!  Deleting interactions CL-QZ if a new list of quantum atoms is given
   if (pupQZchange .ne. 0) then

!    ********** Rebuilding the nonbonding 14 list ***************
!    deleting all connectivity between the QM atoms

!      reinitializing internal nb 14 list structures from the beginning
       numnb14= pupnumnb14
       nbonh  = pupnbonh
       nbona  = pupnbona
       ntheth = pupntheth
       ntheta = pupntheta
       nphih  = pupnphih
       nphia  = pupnphia
       call copy_14nb(pupnb14,ix(inb_14),pupnumnb14)
       do iPup = 1,nbonh
         bs1 = iPup-1
         bs2 = bs1*3
         ix(iibh +bs1) = pupbonh(bs2+1)
         ix(ijbh +bs1) = pupbonh(bs2+2)
         ix(iicbh+bs1) = pupbonh(bs2+3)
       enddo
       do iPup = 1,nbona
         bs1 = iPup-1
         bs2 = bs1*3
         ix(iiba +bs1) = pupbona(bs2+1)
         ix(ijba +bs1) = pupbona(bs2+2)
         ix(iicba+bs1) = pupbona(bs2+3)
       enddo
       do iPup = 1,ntheth
         bs1 = iPup-1
         bs2 = bs1*4
         ix(i24  +bs1) = puptheth(bs2+1)
         ix(i26  +bs1) = puptheth(bs2+2)
         ix(i28  +bs1) = puptheth(bs2+3)
         ix(i30  +bs1) = puptheth(bs2+4)
       enddo
       do iPup = 1,ntheta
         bs1 = iPup-1
         bs2 = bs1*4
         ix(i32  +bs1) = puptheta(bs2+1)
         ix(i34  +bs1) = puptheta(bs2+2)
         ix(i36  +bs1) = puptheta(bs2+3)
         ix(i38  +bs1) = puptheta(bs2+4)
       enddo
       do iPup = 1,nphih
         bs1 = iPup-1
         bs2 = bs1*5
         ix(i40  +bs1) = pupphih(bs2+1)
         ix(i42  +bs1) = pupphih(bs2+2)
         ix(i44  +bs1) = pupphih(bs2+3)
         ix(i46  +bs1) = pupphih(bs2+4)
         ix(i48  +bs1) = pupphih(bs2+5)
       enddo
       do iPup = 1,nphia
         bs1 = iPup-1
         bs2 = bs1*5
         ix(i50  +bs1) = pupphia(bs2+1)
         ix(i52  +bs1) = pupphia(bs2+2)
         ix(i54  +bs1) = pupphia(bs2+3)
         ix(i56  +bs1) = pupphia(bs2+4)
         ix(i58  +bs1) = pupphia(bs2+5)
       enddo

       call deleting_qm_atoms()
       qsetup=.true.

!      Setting as current quantum zone
       pupQZchange = 0

   endif
  ! For PUPIL, rebuild the neighbour list 
  ! and zero the charges on QM atoms at every step
  if ( igb == 0 .and. ipb == 0 .and. iyammp == 0 ) then
         
    ! (for GB: do all nonbondeds together below)
    call timer_start(TIME_NONBON)
    call timer_start(TIME_LIST)
    !do_list_update=.true.         
    call nonbond_list(x,ix(i04),ix(i06),ix(i08),ix(i10), &
                      ntypes,natom/am_nbead,xx,ix,ipairs,ntnb, &
                      ix(ibellygp),belly,newbalance, &
                      qsetup, &
                      do_list_update)
    !call qm_zero_charges(x(L15))
    call timer_stop(TIME_LIST)
    call timer_stop(TIME_NONBON)
  end if
  ! charge reassign here !
  if ( ifcr /= 0 ) then
     call cr_reassign_charge( x, f, pot%ct, xx(l15), natom )
  end if

#endif /*PUPIL_SUPPORT*/


  if(iamoeba.eq.1) then
      vir(1:4)=0.0
   end if
   ! ----------------------------------------------------------------
   ! Calculate the non-bonded contributions
   ! ----------------------------------------------------------------
   call timer_start(TIME_NONBON)


   call timer_start(TIME_EEXIPS)
   if( ips > 0 ) then
      call eexips(evdwex,eelex,istart,iend, ntb,ntypes, &
           ix(i04),ix(i06),ix(i08),ix(i10),xx(l15),cn1,cn2,f,x)
   endif
   call timer_stop(TIME_EEXIPS)

   if( igb == 0 .and. ipb == 0 .and. iyammp == 0 ) then

      ! (for GB: do all nonbondeds together below)

      call timer_start(TIME_EWALD)

      if ( iamoeba == 1 )then
         call AM_NonBond_eval(natom,x,f,vir,xx,ipairs, &
                               evdw,eelt,epolar,&
                               enb14,ee14,diprms,dipiter)
      else

         if ( induced > 0 ) then
              call handle_induced(x,natom,ix(i04),ix(i06), &
                 xx(l15),cn1,cn2,cn6,eelt,epolar,f,xx,ix, &
                 ipairs,xx(lpol),xx(lpol2), &
                 xx(l45),virvsene,ix(i02),ibgwat,nres, &
                 aveper,aveind,avetot,emtot,diprms,dipiter,dipole_temp, &
#ifdef HAS_10_12
                 cn3,cn4,cn5,asol,bsol)
#else
                 cn3,cn4,cn5)
#endif
         else
            if ( ilrt /= 0 ) then
               ! Modifications for computing interaction energy
               ! according to the Linear Response Theory, LIE module
               if ( do_lrt ) then
                  ! call with molecule charges set to zero
                  call ewald_force(x,natom,ix(i04),ix(i06), &
                       crg_m0,cn1,cn2,cn6,energy_m0,epolar, &
                       f_scratch,xx,ix,ipairs,xx(l45),virvsene, xx(lpol), &
#ifdef HAS_10_12
                       xx(lpol2), .false. , cn3, cn4, cn5, asol, bsol)
#else
                       xx(lpol2), .false. , cn3, cn4, cn5)
#endif
                  !write (6,*) 'Em0', energy_m0
                  ! call with water charges set to zero
                  call ewald_force(x,natom,ix(i04),ix(i06), &
                       crg_w0,cn1,cn2,cn6,energy_w0,epolar, &
                       f_scratch,xx,ix,ipairs,xx(l45),virvsene, xx(lpol), &
#ifdef HAS_10_12
                       xx(lpol2), .false. , cn3, cn4, cn5, asol, bsol)
#else
                       xx(lpol2), .false. , cn3, cn4, cn5)
#endif
                  !write (6,*) 'Ew0', energy_w0
                  ! call with full charges but no vdw interaction between solute and solvent
                  call ewald_force(x,natom,ix(i04),ix(i06), &
                       xx(l15),cn1_lrt,cn2_lrt,cn6,eelt,epolar, &
                       f_scratch,xx,ix,ipairs,xx(l45),virvsene, xx(lpol), &
#ifdef HAS_10_12
                       xx(lpol2), .false. , cn3, cn4, cn5, asol, bsol)
#else
                       xx(lpol2), .false. , cn3, cn4, cn5)
#endif
                  energy_vdw0 = evdw
                  call lrt_solute_sasa(x,natom, xx(l165))
               end if
               ! call normal_ewald force this will overwrite everything 
               ! computed above except energy_m0 and energy_w0
               call ewald_force(x,natom,ix(i04),ix(i06), &
                    xx(l15),cn1,cn2,cn6,eelt,epolar, &
                    f,xx,ix,ipairs,xx(l45),virvsene, xx(lpol), &
#ifdef HAS_10_12
                    xx(lpol2), .false. , cn3, cn4, cn5, asol, bsol)
#else
                    xx(lpol2), .false. , cn3, cn4, cn5)
#endif
               energy_vdw0 = evdw - energy_vdw0
               ! count call to ltr, maybe calculate Eee and print it
               call ee_linear_response(eelt, master)
            else ! just call ewald_force normally
               call ewald_force(x,natom,ix(i04),ix(i06),&
                                xx(l15),cn1,cn2,cn6,eelt,epolar, &
                                f,xx,ix,ipairs,xx(l45),virvsene,xx(lpol), &
#ifdef HAS_10_12
                                xx(lpol2), .false. , cn3, cn4, cn5, asol, bsol)
#else
                                xx(lpol2), .false. , cn3, cn4, cn5)
#endif
!              write(0,*) 'EWALD_FORCE; eelt = ', eelt
  
            end if ! ilrt /= 0
         end if ! induced > 0
      end if ! iamoeba == 1

      call timer_stop(TIME_EWALD)

#ifdef MPI
      if(mytaskid == 0)then
#endif
         pot%vdw      = evdw
         pot%elec     = eelt
         pot%hbond    = ehb  !whereis ehb?
#ifdef MPI
      else
         ! energies have already been reduced to the master
         ! node in ewald_force, so here we zero out elements
         ! on non-master nodes:
         pot%vdw      = 0.d0
         pot%elec     = 0.d0
         pot%hbond    = 0.d0

      end if
#endif

      if( ips > 0 )then
         pot%vdw   = pot%vdw   + evdwex
         pot%elec  = pot%elec  + eelex
      endif

   end if  ! ( igb == 0 .and. ipb == 0 .and. iyammp == 0 )

   call timer_stop(TIME_NONBON)

   ! ----------------------------------------------------------------
   ! Calculate the other contributions
   ! ----------------------------------------------------------------

   !     -- when igb==10, all nonbonds are done in routine pb_force, and
   !                      all nonpolar interactions are done in np_force:
   !
   !     -- HG put this part here such that "outflag" is known from a call
   !        of pb_force; outflag is needed in the "bond" routine in the case
   !        of ifcap == 2,5 (i.e., ivcap == 1,5)

#ifdef MPI
   if(mytaskid == 0)then
#endif
      if( igb == 10 .or. ipb /= 0 ) then
         call timer_start(TIME_PBFORCE)
         call pbtimer_init
         call pb_force(natom,nres,ntypes,npdec,ix(i02),ix(i04),ix(i06),ix(i10), &
                 cn1,cn2,xx(l15),x,f,evdw,eelt,epol)
         if ( pbgrid ) pbgrid = .false.
         if ( pbinit ) pbinit = .false.
         pot%vdw  = evdw
         pot%elec = eelt
         pot%pb   = epol
         call timer_stop(TIME_PBFORCE)

         call timer_start(TIME_NPFORCE)
         esurf = 0.0d0; edisp = 0.0d0
         if ( ifcap == 0  .and. npopt /= 0 ) &
            call np_force(natom,nres,ntypes,ix(i02),ix(i04),ix(i06),&
                 cn1,cn2,x,f,esurf,edisp)
         if ( pbprint ) pbprint = .false.
         pot%surf = esurf
         pot%disp = edisp
         call pbtimer_summary
         call timer_stop(TIME_NPFORCE)

      end if  ! ( igb == 10 .or. ipb /= 0 )

#ifdef MPI
   end if
#endif

!  +---------------------------------------------------------------+
!  |  Bonds with H                                                 |
!  +---------------------------------------------------------------+

   call timer_start(TIME_BOND)

   ! initialize bond virial
   if(ipimd>0) bnd_vir = zero

#ifdef MPI /* SOFT CORE */
   ! zero only once, sc bond energy is sum of H and non-H terms
   sc_ener(1) = 0.0d0
#endif

   if( ntf < 2 ) then

      ebdev = 0.d0
      call bond(nbonh,ix(iibh),ix(ijbh),ix(iicbh),x,f,ene(6))
      pot%bond = pot%bond + ene(6)
#ifdef MPI
#  ifdef LES
      if(rem == 2) then
         pot%les = pot%les + elesb
      endif
#  endif
#endif
   end if

!  +---------------------------------------------------------------+
!  |  Bonds without H                                              |
!  +---------------------------------------------------------------+

   if( ntf < 3 ) then

      call bond(nbona+nbper,ix(iiba),ix(ijba),ix(iicba),x,f,ene(7))
      pot%bond = pot%bond + ene(7)
#ifdef MPI
#  ifdef LES
      if(rem == 2) then
         pot%les = pot%les + elesb
      endif
#  endif
#endif
      if (nbonh+nbona > 0) ebdev = sqrt( ebdev/(nbonh+nbona) )
   end if

!  +---------------------------------------------------------------+
!  |  Angles with H                                                |
!  +---------------------------------------------------------------+

   if( ntf < 4 ) then

#ifdef MPI /* SOFT CORE */
      ! zero only once, sc bond energy is sum of H and non-H terms
      sc_ener(2) = 0.0d0
#endif

      eadev = 0.d0

      call angl(ntheth,ix(i24),ix(i26),ix(i28),ix(i30),x,f,ene(8))
      pot%angle = pot%angle + ene(8)
#ifdef MPI
#  ifdef LES
      if(rem == 2) then
         pot%les = pot%les + elesa
      endif
#  endif
#endif
   end if

!  +---------------------------------------------------------------+
!  |  Angles without H                                             |
!  +---------------------------------------------------------------+

   if( ntf < 5 ) then

      call angl(ntheta+ngper,ix(i32),ix(i34),ix(i36),ix(i38),x,f, &
           ene(9)) 
      pot%angle = pot%angle + ene(9)

#ifdef MPI
#  ifdef LES
      if(rem == 2) then
         pot%les = pot%les + elesa
      endif
#  endif
#endif
      if (ntheth+ntheta > 0) eadev = 57.296*sqrt( eadev/(ntheth+ntheta) )
   end if

!  +--------------------------------------------------------------------+
!  | AMD calculate dihedral energy first to estimate dihedral weight,   |
!  | then use it in the regular ephi function. Added by Romelia Salomon |
!  | Fix me later, AMD dihedral weight does NOT support AMOEBA          |
!  +--------------------------------------------------------------------+

   if(iamd .gt. 1)then
!     Dihedrals with H 
      if( ntf < 6 ) then
         call ephi_ene_amd(nphih,ix(i40),ix(i42),ix(i44),ix(i46),ix(i48), &
                           x,amd_dih_noH)
      endif
!     Dihedrals without H 
      if( ntf < 7 ) then
        call ephi_ene_amd(nphia+ndper,ix(i50),ix(i52),ix(i54),ix(i56),ix(i58), &
             x,amd_dih_H)
      endif
#ifdef MPI 
      temp_amd_totdih = amd_dih_noH + amd_dih_H
#  ifdef USE_MPI_IN_PLACE
      call mpi_allreduce(MPI_IN_PLACE,temp_amd_totdih,1,MPI_DOUBLE_PRECISION,mpi_sum,commsander,ierr)
      amd_totdih = temp_amd_totdih
#  else
      call mpi_allreduce(temp_amd_totdih, amd_totdih, &
                         1, MPI_DOUBLE_PRECISION, &
                         mpi_sum, commsander, ierr)
#  endif /* USE_MPI_IN_PLACE */
#else
     amd_totdih = amd_dih_noH + amd_dih_H
#endif /* MPI */
     call calculate_amd_dih_weights(amd_totdih)
   endif

!  +---------------------------------------------------------------+
!  |  Dihedrals with H                                             |
!  +---------------------------------------------------------------+
   ! initialize 14 nb energy virial
   if(ipimd>0) e14vir = zero

   if( ntf < 6 ) then

#ifdef MPI /* SOFT CORE */
      ! zero only once, sc bond energy is sum of H and non-H terms
      sc_ener(3) = 0.0d0
#endif

      call ephi(nphih,ix(i40),ix(i42),ix(i44),ix(i46),ix(i48), &
           xx(l15),ix(i04),x,f,dvdl,ene(10),ene(11),ene(12),xx(l190))

      pot%dihedral = pot%dihedral + ene(10) ! Combine contributions from dihedrals with H
      pot%vdw_14   = pot%vdw_14   + ene(11) ! Combine 1-4 vdw contributions from dihedrals with H
      pot%elec_14  = pot%elec_14  + ene(12) ! Combine 1-4 elec contributions from dihedrals with H


#ifdef MPI
#  ifdef LES
      if(rem == 2) then
         pot%les = pot%les + elesd
      endif
#  endif
#endif
   end if

!  +---------------------------------------------------------------+
!  |  Dihedrals without H                                          |
!  +---------------------------------------------------------------+

   if( ntf < 7 ) then

      call ephi(nphia+ndper,ix(i50),ix(i52),ix(i54),ix(i56),ix(i58), &
           xx(l15),ix(i04),x,f,dvdl,ene(13),ene(14),ene(15),xx(l190))

      pot%dihedral = pot%dihedral + ene(13) ! Combine contributions from dihedrals without H
      pot%vdw_14   = pot%vdw_14   + ene(14) ! Combine 1-4 vdw contributions from dihedrals without H
      pot%elec_14  = pot%elec_14  + ene(15) ! Combine 1-4 elec contributions from dihedrals without H

#ifdef MPI
#  ifdef LES
      if(rem == 2) then
         pot%les = pot%les + elesd
      endif
#  endif
#endif

!  +---------------------------------------------------------------+
!  |  CHARMM IMPROPERS IF CHARMM FORCEFIELD IS IN USE              |
!  +---------------------------------------------------------------+
!  Note: CHARMM does not distinguish between impropers with and 
!        without hydrogen hence it is not possible to strictly
!        conform to the ntf options of sander. Here CHARMM impropers
!        are calculated as long as ntf < 7 - so CHARMM impropers are
!        essentially considered to be in the same set as dihedrals
!        NOT involving hydrogen.
     if (charmm_active) then
       call charmm_calc_urey_bradley(x,pot%angle_ub,f)
       call charmm_calc_impropers(x,pot%imp,f)
       call charmm_calc_cmap(x,pot%cmap,f)
     end if

! --- END CHARMM IMPROPERS IF REQUIRED ---
     if (cmap_active) then
        call calc_cmap(x,pot%cmap,f)
     end if

   end if !ntf < 7

   if(iamoeba==1) then
      call AM_VAL_eval(x,f,vir,ene(6),ene(8),ene(10))
                              !ebond, eangle,etors
      pot%bond     = pot%bond     + ene(6)
      pot%angle    = pot%angle    + ene(8)
      pot%dihedral = pot%dihedral + ene(10)

   end if

   call timer_stop(TIME_BOND)

   if (sebomd_obj%do_sebomd) then
     ! step 1/2 to save emap + entr + ecap forces
     call sebomd_save_forces(1,natom,f,xx(gradsebomd),xx(grad1tmp),xx(grad2tmp))
   end if

   ! --- calculate the EMAP constraint energy ---

   if(temap) then   ! ntr=1 (positional restraints)
       call emapforce(natom,enemap,xx(lmass),x,f )
       pot%emap = enemap
   end if

   ! --- calculate the position constraint energy ---

#ifdef MPI /* SOFT CORE */
      sc_ener(14:19) = 0.d0 ! zero all restraint/constraint energies
#endif
   if(natc > 0 .and. ntr==1) then   ! ntr=1 (positional restraints)
       call xconst(natc,entr,ix(icnstrgp),x,f,xx(lcrdr),xx(l60))
       pot%constraint = entr
   end if

   if ( itgtmd==1 .and. (nattgtfit > 0 .or. nattgtrms > 0) ) then

      ! Calculate rmsd for targeted md (or minimization) if requested.
      ! All nodes do rms fit, could just be master then broadcast.
      ! All nodes need all coordinates for this.

      call rmsfit(xx(lcrdr),x,xx(lmass),ix(itgtfitgp),  &
                  ix(itgtrmsgp),rmsdvalue,nattgtrms,nattgtfit,rmsok)

      if (.not.rmsok) then
         if (master) write (6,*) 'Fatal error: Error calculating RMSD!'
         call mexit(6, 1)
      end if

      call xtgtmd(entr,ix(itgtrmsgp),x,f,xx(lcrdr),xx(lmass),tgtrmsd,tgtmdfrc,rmsdvalue,nattgtrms)
      pot%constraint = entr
   else if(itgtmd == 2) then
      call mtmdcall(entr,xx(lmtmd01),ix(imtmd02),x,f,ih(m04),ih(m02),ix(i02),&
                    ih(m06),xx(lmass),natom,nres,'CALC')
      pot%constraint = entr
   end if

   if(ifcap == 1 .or. ifcap == 2) then
      call capwat(natom,x,f,ecap)
      pot%constraint = pot%constraint + ecap
   else if(ifcap == 3) then
      write(6,*) 'No energy expression for spherical boundary known yet'
      call mexit(6,1)
   else if(ifcap == 4) then
      write(6,*) 'No energy expression for orthorhombic boundary known yet'
      call mexit(6,1)
      !call orth(natom,ix(ibellygp),x,f,eorth)
      !ene(20) = ene(20) + eorth
   end if

   if (sebomd_obj%do_sebomd) then
     ! step 2/2 to save emap + entr + ecap forces
     call sebomd_save_forces(2,natom,f,xx(gradsebomd),xx(grad1tmp),xx(grad2tmp))
   end if

   ! No energy expression for ifcap == 5 given because only
   !    one step of minimization is allowed with this.

   !  (this seems very weird: we have already done an allreduce on molvir
   !  in ewald_force(); this just collects it on processor 0 (with zeroes
   !  on all slave nodes), then later does an allreduce...)

   if( mytaskid == 0 .and. iamoeba == 0 ) then
      vir(1) = vir(1)+0.5d0*molvir(1,1)
      vir(2) = vir(2)+0.5d0*molvir(2,2)
      vir(3) = vir(3)+0.5d0*molvir(3,3)
   end if

   if( igb == 0 .and. ipb == 0 .and. iyammp == 0 ) then
      ener%virvsene    = virvsene
      ener%diprms      = diprms
      ener%dipiter     = dipiter
      ener%dipole_temp = dipole_temp
   end if

   !     ---- get the noesy volume penalty energy: ------

   pot%noe = 0.d0
   if( iredir(4) /= 0 ) then
      call timer_start(TIME_NOE)
      call noecalc(x,f,xx,ix)
      call timer_stop(TIME_NOE)
   end if
   ! Do we need a pot%noe here?  mjw TODO

   !     -- when igb!=0 and igb!=10, all nonbonds are done in routine egb:

   esurf = 0.d0
   if( igb /= 0 .and. igb /= 10 .and. ipb == 0 ) then
      call timer_start(TIME_EGB)
      call egb( x,f,rborn,fs,reff,onereff,xx(l15),ix(i04),ix(i06), &
            ix(i08),ix(i10),xx(l190), &
            cut,ntypes,natom,natbel,epol,eelt,evdw, &
            esurf,dvdl,xx(l165),ix(i82),xx(l170),xx(l175),xx(l180), &
            xx(l185), ncopy &
#ifndef LES
            , xx(l2402),xx(l2403),xx(l2404) &
#endif
            )

      pot%vdw  = evdw
      pot%elec = eelt
      pot%gb   = epol
      pot%surf = esurf
      pot%dvdl = dvdl
      call timer_stop(TIME_EGB)
#ifdef MPI
#  ifdef LES
      if(rem == 2) then
        pot%les = pot%les + elesp
      endif
#  endif
#endif

   end if  ! ( igb /= 0 .and. igb /= 10 .and. ipb == 0 )

#ifdef RISMSANDER
   
   if(rismprm%irism == 1) then
      call timer_start(TIME_RISM)
      call rism_force(x,f,erism,irespa)
      pot%rism = erism
      call timer_stop(TIME_RISM)
   endif
#endif

#ifdef APBS
! APBS forces
      if( mdin_apbs ) then
         if (igb /= 6) then
            write(6, '(a)') '&apbs keyword requires igb=6.'
            call mexit(6,1)
         end if
         call timer_start(TIME_PBFORCE)
! in: coords, radii, charges
! out: updated forces (via apbs_params) and solvation energy (pol + apolar)
         if (sp_apbs) then
            call apbs_spenergy(natom, x, f, eelt, enpol)
         else
            call apbs_force(natom, x, f, pot%vdw, eelt, enpol)
         end if
         pot%pb   = eelt 
         pot%surf = enpol
         call timer_stop(TIME_PBFORCE)

      end if  ! ( mdin_apbs )
#endif /* APBS */

   if (sebomd_obj%do_sebomd) then
     ! step 1/2 to save eshf + epcshf + ealign + ecsa forces
     call sebomd_save_forces(1,natom,f,xx(gradsebomd),xx(grad1tmp),xx(grad2tmp))
   end if

   if( master ) then
      !  These parts of the NMR energies are not parallelized, so only
      !  are done on the master node:
      eshf = 0.d0
      epcshf = 0.d0
      ealign = 0.d0
      ecsa = 0.d0
      if (iredir(5) /= 0) call cshf(natom,x,f)
      if (iredir(7) /= 0) call pcshift(natom,x,f)
      if (iredir(9) /= 0) call csa1(natom,x,f)
      if (iredir(8) /= 0) call align1(natom,x,f,xx(lmass))
   end if

   if (sebomd_obj%do_sebomd) then
     ! step 2/2 to save eshf + epcshf + ealign + ecsa forces
     call sebomd_save_forces(2,natom,f,xx(gradsebomd),xx(grad1tmp),xx(grad2tmp))
   end if

   ! additional force due to charge relocation
   if ( ifcr /= 0 ) then
      call cr_calc_force( f )
   end if

#ifdef MPI

   call timer_barrier( commsander )
   call timer_start(TIME_COLLFRC)

   !     add force, ene, vir, copies from all nodes
   !            also add up newbalance for nonperiodic.

   ! Remember to work on the local instance of the
   ! potential energy array, i.e. pot and NOT the global one,
   ! i.e. ener%pot

   call fdist(f,xx(lfrctmp),pot,vir,newbalance)


   call timer_stop(TIME_COLLFRC)

#endif

   ! ---- at this point, the parallel part of the force calculation is
   !      finished, and the forces have been distributed to their needed
   !      locations.  All forces below here are computed redundantly on
   !      all processors, and added into the force vector.  Hence, below
   !      is the place to put any component of the force calculation that
   !      has not (yet) been parallelized.

   if (sebomd_obj%do_sebomd) then
     ! step 1/2 to save enmr + edssp forces
     call sebomd_save_forces(1,natom,f,xx(gradsebomd),xx(grad1tmp),xx(grad2tmp))
   end if

   ! Calculate the NMR restraint energy contributions, if requested.
   ! (Even though this is not parallelized, it needs to be run on all
   ! threads, since this code is needed for weight changes as well as
   ! for NMR restraint energy analysis.  The whole thing could stand a
   ! major re-write....)

   if (nmropt > 0) &
      call nmrcal(x,f,ih(m04),ih(m02),ix(i02),xx(lwinv),enmr,devdis, &
                  devang,devtor,devplpt,devpln,devgendis,temp0,tautp,cut, &
                  xx(lnmr01),ix(inmr02),xx(l95),31,6,rk,tk,pk,cn1, &
                  cn2,asol,bsol,xx(l15),numbnd,numang,nptra-nimprp, &
                  nimprp,nphb,natom,natom,ntypes,nres, &
                  rad,wel,radhb,welhb,rwell,tgtrmsd,temp0les,-1,'CALC')
#ifdef MPI
   call mpi_reduce(enoe,pot%noe,1,MPI_DOUBLE_PRECISION,mpi_sum,0,commsander,ierr)
   enoe = pot%noe ! so all processors now have the full enoe value
#else
   pot%noe = enoe
#endif

#ifdef DSSP
   if( idssp > 0 ) then
      call fdssp( natom,x,f,edssp )
!     write(6,*) 'edssp = ', edssp
   else
      edssp = 0.d0
   end if
#endif

   if (sebomd_obj%do_sebomd) then
     ! step 2/2 to save enmr + edssp forces
     call sebomd_save_forces(2,natom,f,xx(gradsebomd),xx(grad1tmp),xx(grad2tmp))
   end if

   !     ----- CALCULATE TOTAL ENERGY AND GROUP THE COMPONENTS -----

#ifndef LES
   if( igb == 0 .and. ipb == 0 ) then
      !ene(11) = enb14
      !ene(14) = 0.d0
      !ene(12) = ee14
      !ene(15) = 0.d0 !mjw TODO Grok this...
      pot%vdw_14   = pot%vdw_14   + enb14
      pot%elec_14  = pot%elec_14  + ee14

   endif

#endif

   pot%constraint = pot%constraint+ eshf + epcshf+ pot%noe + &
                    sum(enmr(1:6)) + ealign + ecsa + pot%emap



#ifdef DSSP
   pot%constraint = pot%constraint + edssp
#endif

   pot%polar = epolar

   pot%tot=      pot%vdw        + &
                 pot%elec       + &
                 pot%gb         + &
                 pot%pb         + &
                 pot%bond       + &
                 pot%angle      + &
                 pot%dihedral   + &
                 pot%vdw_14     + &
                 pot%elec_14    + &
                 pot%hbond      + &
                 pot%constraint + &
                 pot%rism       + &
                 pot%ct

   pot%tot = pot%tot + pot%polar + pot%surf + pot%scf + pot%disp

   !Charmm related
   pot%tot = pot%tot + pot%angle_ub + pot%imp + pot%cmap 

   if (sebomd_obj%do_sebomd) then
     ! apply lambda term to forces and energy between full QM (SEBOMD) calculation
     ! and full MM calculation
     ! if lambda .eq. 1.0 (defaults), this means MM calculations is zeroed

     ! step 3: apply restraint forces to SEBOMD forces
     call sebomd_save_forces(3,natom,f,xx(gradsebomd),xx(grad1tmp),xx(grad2tmp))

     ! new potential energy
     ! (to have full restraint, we add pot%constraint to the SEBOMD energy)
     ! (see sebomd_save_forces subroutine for explanation)
     pot%tot = (one-sebomd_obj%lambda)*pot%tot &
                  + sebomd_obj%lambda*(sebomd_obj%esebomd+pot%constraint)

     ! new forces
     do i = 1,3*natom
       f(i) = (one-sebomd_obj%lambda)*f(i) +sebomd_obj%lambda*xx(gradsebomd+i-1)
     end do
   end if

!  +---------------------------------------------------------------+
!  | AMD calculate total potential energy weight, then apply it to |
!  | all the force elements f=f*fwgt. Added by Romelia Salomon     |
!  +---------------------------------------------------------------+

   if(iamd .gt. 0)then
!     call calculate_amd_total_weights(natom,pot%tot,amd_totdih,pot%amd_boost,f,temp0)
     call calculate_amd_total_weights(natom,pot%tot,pot%dihedral,pot%amd_boost,f,temp0)
     ! Update total energy
     pot%tot = pot%tot + pot%amd_boost
   end if

!  +---------------------------------------------------------------+
!  | scaledMD calculate scale the total potential forces. All the  |
!  | force elements f=f*scaledMD_lambda. Added by Romelia Salomon  |
!  +---------------------------------------------------------------+

   if(scaledMD .gt. 0)then
     call scaledMD_scale_frc(natom,pot%tot,f)
     pot%tot = pot%tot *scaledMD_lambda
   end if

   !The handover
   ener%pot = pot


   ener%aveper = aveper
   ener%aveind = aveind
   ener%avetot = avetot



   
   ! This is now historical; MJW Feb 2010
   !
   !    Here is a summary of how the ene array is used.  For parallel runs,
   !    these values get summed then rebroadcast to all nodes (via
   !    mpi_allreduce).

   !    ene(1):      total energy
   !    ene(2):      van der Waals
   !    ene(3):      electrostatic energy
   !    ene(4):      10-12 (hb) energy, or GB energy when igb.gt.0
   !    ene(5):      bond energy
   !    ene(6):      angle energy
   !    ene(7):      torsion angle energy
   !    ene(8):      1-4 nonbonds
   !    ene(9):      1-4 electrostatics
   !    ene(10):     constraint energy
   !    ene(11-19):  used as scratch, but not needed further below
   !    ene(20):     position constraint energy + cap energy
   !    ene(21):     charging free energy result
   !    ene(22):     noe volume penalty
   !    ene(23):     surface-area dependent energy, or cavity energy
   !    ene(24):     potential energy for a subset of atoms
   !    ene(25):     SCF Energy when doing QMMM
   !    ene(26):     implicit solvation dispersion energy


#ifdef PUPIL_SUPPORT
   !*****************************************************
   !     Closing the qmmm structure consideration
   !*****************************************************
!  Adding the quantum forces from last QM calculation
   do iPup=1,pupparticles
     bs1 = (abs(pupqlist(iPup))-1)*3
     !write(6,"(a10,2x,i4,3(2x,e16.10))") 'Classic F:',abs(pupqlist(iPup)), &
     !                                              (f(bs1+jPup),jPup=1,3)
     do jPup=1,3
       bs2    = bs1    + jPup
       f(bs2) = f(bs2) + qfpup(bs2)
     enddo
     !write(6,"(a10,2x,i4,3(2x,e16.10))") 'Quantum F:',abs(pupqlist(iPup)), &
     !                                          (qfpup(bs1+jPup),jPup=1,3)
     !write(6,"(a10,2x,i4,3(2x,e16.10))") 'TOTAL F:',abs(pupqlist(iPup)), &
     !                                          (f(bs1+jPup),jPup=1,3)
   enddo

   ! If there are more that one QM Domain add vdw interaction
   ! among qm particles from different QM Domains
   if (pupnumdomains .gt. 1) then
      call add_vdwqmqm(r_stack(l_puptmp),f,ener,ntypes,ih(m04),ih(m06),ix(i04))
   endif

!  Deallocating temporary stack
   call free_stack(l_puptmp,routine)

!Final forces
!   do iPup=1,natom
!     bs2 = (iPup-1)*3
!     write(6,"(a10,2x,i4,3(2x,e16.10))") 'Total F:',iPup,(f(bs2+jPup),jPup=1,3)
!   enddo

!  Disconnecting qmmmm interactions
   qmmm_nml%ifqnt = .false.

#endif

   ! ----ADD X-RAY TARGET FUNCTION AND GRADIENT
   call cns_xref_run(natom,ih(m04), x,f,ener)

   ! ---- BUILT-IN X-RAY TARGET FUNCTION AND GRADIENT
   if (xray_active) call xray_get_derivative(x,f)

   !     if freezemol has been set, zero out all of the forces for the
   !     real atoms; (no longer necessary to set ibelly).
   if( ifreeze > 0 ) then
      do i=1,3*natom
         f(i) = 0.d0
      end do
   end if

   !     ----- IF BELLY IS ON THEN SET THE BELLY ATOM FORCES TO ZERO -----
   if (belly) call bellyf(natom,ix(ibellygp),f)

!  +---------------------------------------------------------------+
!  |  Interface to EVB                                             |
!  +---------------------------------------------------------------+

#ifdef MPI
#  ifdef LES
!KFW   if( ipimd>0) then
   if( nbead > 0 ) then
      call mpi_allreduce ( nrg_all, nrg_bead, nbead, MPI_DOUBLE_PRECISION &
                         , MPI_SUM, commsander, ierr )
      nrg_all(:) = nrg_bead(:)
   end if
   if( ievb /= 0 ) call evb_ntrfc ( x, f, ener, ix, ipairs, vel0_nrg_sum )
#  else
   if( ievb /= 0 ) call evb_ntrfc ( x, f, ener, xx(lmass), ix, ipairs )
#  endif /* LES */
#endif /* MPI */

   if( ipimd>0 ) then
#ifdef LES
      call pimd_part_spring_force(x,f,real_mass,Epot_spring,Epot_deriv,dvdl)
#  ifdef MPI
      if( ievb /= 0 ) then
         nrg_frc(3)= vel0_nrg_sum
         nrg_frc(2)= equal_part + Epot_deriv
         nrg_frc(1)= nrg_frc(3) + nrg_frc(2)
         dlnQ_dl = dvdl
      endif
#  endif /* MPI */
#else /* NOT LES below */
      ener = ener/nbead 
      f(1:natom*3) = f(1:natom*3)/nbead
      vir(1:3) = vir(1:3) /nbead
      atvir = atvir/nbead
      e14vir = e14vir/nbead
      bnd_vir = bnd_vir/nbead
      call pimd_full_spring_force(x,f,real_mass,Epot_spring,Epot_deriv,dvdl)
#  ifdef MPI
      if (master) call mpi_reduce(ener,totener,state_rec_len,MPI_DOUBLE_PRECISION, &
                      MPI_SUM,0,commmaster,ierr)
#  endif /* MPI */
#endif /* LES */
      ! Pass dvdl = dV/dl for TI w.r.t. mass.
      if (itimass > 0) ener%pot%dvdl = dvdl
   end if
#ifdef MPI

     !CARLOS: ONLY MPI SUPPORTS NEB (MULTISANDER)
   if(ineb>0) then
      if(sanderrank.eq.0) then  !only masters do NEB
         call full_neb_forces( xx(lmass), x, f, ener%pot%tot, ix(itgtfitgp),ix(itgtrmsgp))
      endif

      ! now master will broadcast the neb forces for the rmsgp atoms
      ! all nodes will add this to the current force total
      ! if we wanted all nodes to calculate the neb forces, they
      ! would all need access to the neighbor bead coordinates, which is
      ! probably more expensive than having the master send out the modified
      ! forces

      ! master broadcasts the force update for the rmsgp atoms only

      call mpi_bcast(neb_force,nattgtrms*3,MPI_DOUBLE_PRECISION,0,commsander,ierr)

      do i=1,nattgtrms
         j3 = 3*(i - 1) !pointer into the packed neb force array
         j=ix(itgtrmsgp+i-1) !actual atom # for atom in rms group
         i3 = 3*(j - 1) !pointer into real force array

         f(i3+1)=f(i3+1)+neb_force(j3+1)
         f(i3+2)=f(i3+2)+neb_force(j3+2)
         f(i3+3)=f(i3+3)+neb_force(j3+3)
      enddo

! CARLOS: what is this doing? ener(27) is neb energy
! looks like it reduces entire energy array
! needs MUCH better documentation on details (such as 28)

      if(sanderrank.eq.0) then
         call mpi_reduce(ener,totener,state_rec_len,MPI_DOUBLE_PRECISION,MPI_SUM,0,commmaster,ierr)
      end if

   end if !ineb>0
#endif /*MPI*/


   if (charmm_active) then
     if ( do_charmm_dump_gold == 1 ) then
       call charmm_dump_gold(f,natom,ener)
     endif 
   end if


#ifdef MPI

   if(icfe == 0) then
      if(idecomp == 1 .or. idecomp == 2) then
         call collect_dec2(nres)
      end if
          
      if(idecomp >= 3) then
         call collect_dec2(npdec*npdec)
      end if
   end if

#endif

   call timer_stop(TIME_FORCE)
   call trace_exit( 'force' )
   return


end subroutine force


