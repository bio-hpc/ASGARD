! <compile=optimized>
#include "copyright.h"
#include "../include/assert.fh"
#include "../include/dprec.fh"

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Main routine for Amber's traditional minimization methods
!-----------------------------------------------------------------------
!     --- RUNMIN ---
!-----------------------------------------------------------------------
! Various combinations of steepest descent and conjugate gradient
! optimizations.

subroutine runmin(xx,ix,ih,ipairs,x,fg,w,ib,jb,conp, &
      winv,igrp,skips,ene,carrms, qsetup)

   use fastwt
   use constants, only : zero, one, TEN_TO_MINUS5, TEN_TO_MINUS6
   use qmmm_module, only : qmmm_nml, qmmm_struct, qmmm_mpi, qm2_struct
   use poisson_boltzmann, only: outwat, oution
   use bintraj, only: end_binary_frame
   use file_io_dat

#if defined( MPI )
   use evb_data, only: evb_frc
#endif /* MPI */
   
#ifdef MPI /* SOFT CORE */
   use softcore, only: extra_atoms, sc_ener, sc_dvdl, sc_tot_dvdl, &
                       sc_tot_dvdl_partner, ifsc, sc_mix_sum, sc_print_energies,&
                       sc_dvdl_ee, sc_tot_dvdl_ee, sc_tot_dvdl_partner_ee, &
                       ti_ene_cnt, nsoftcore, nsoftcore_partner, nmixed, nsc
#endif

   use state
   use sebomd_module, only : sebomd_obj
   implicit none

#ifdef MPI
#  include "parallel.h"
#  ifdef MPI_DOUBLE_PRECISION
#    undef MPI_DOUBLE_PRECISION
#  endif
   include 'mpif.h'
   integer ierr
#ifdef CRAY_PVP
#  define MPI_DOUBLE_PRECISION MPI_REAL8
#endif
   integer ist(MPI_STATUS_SIZE), partner
#endif
#include "../include/md.h"
#include "box.h"
#include "../include/memory.h"
#include "nmr.h"
#include "extra.h"
#include "ew_cntrl.h"
#include "../pbsa/pb_md.h"

   ! ------ passed in variables --------------------
   _REAL_   xx(*)
   integer  ix(*), ipairs(*)
   character(len=4) ih(*)
   _REAL_   x(*),fg(*),w(*)
   integer  ib(*),jb(*)
   _REAL_   conp(*),winv(*)
   integer  igrp(*)
   logical  skips(*)
   type(state_rec) ::  ene
   _REAL_   carrms
   logical :: qsetup

   ! ------ External Functions -----------------
   _REAL_ ddot

   ! ------ local variables --------------------
   logical skip,newstr,steep,belly
   _REAL_ betax, ddspln, dfpr, dxst, dxsth
   _REAL_ f, fch, finit, fmin, fnq, fold, gamma, gamden
   _REAL_ ginit, gmin, gnew, gspln, gsqrd, sbound, step
   _REAL_ stepch, stmin, sum, work
   _REAL_ rms,fndfp,swork
   integer i,n,nr,nrx,nct,ndfp,nstcyc,nitp,ier
   integer mstcyc,linmin,iterrs
   integer irsdx,irsdg,iginit,ixopt,igopt,iterc,n_force_calls
   integer iterfm,nfopt,nfbeg,iretry
   logical itdump,ixdump
   logical loutfm
   logical qspatial

   integer maxlin,mxfcon,kstcyc
   ! comments on parameters are guesses by SRB Sep 2003
   parameter ( maxlin = 10 )  ! maximum number of line searches ?
   parameter ( mxfcon =  4 )  ! maximum force con ?
   parameter ( kstcyc =  4 )  ! number of starting cycles ?
   _REAL_  crits, dxstm, dfpred
   parameter ( dxstm  = TEN_TO_MINUS5 )  ! ?
   parameter ( crits  = TEN_TO_MINUS6 )  ! ?
   parameter ( dfpred = ONE )            ! ? in kcal/mol
   logical :: do_list_update=.false.
   logical :: lout

#ifdef MPI
   integer :: ti_ib, ti_jb, ti_nct(2), ti_nct_partner(2), ti_nct_tot, ti_nrp
   _REAL_  :: ti_pot_ene, ti_pot_ene_partner, ti_sum, ti_sum_partner
   _REAL_, allocatable :: frcti(:)
   type(state_rec) :: ecopy
#endif

   ! Zero the state type as done in runmd()
   ene = null_state_rec

   !     ----- EVALUATE SOME CONSTANTS AND INITIALIZE SOME VARIABLES -----

   fmin = 0.0d0
   nr = nrp
   n = 3*nr
   belly = ibelly > 0
   ier = 0
   nct = 0
   if (ntc == 2) nct = nbonh
   if (ntc == 3) nct = nbonh + nbona
   ndfp = n-nct
   if(belly) ndfp = 3*natbel-nct
   ntnb = 1
   fndfp = ndfp
   fnq = sqrt(fndfp)
   sbound = -1.d0
   stmin = 0.d0
   step = 0.d0
   stepch = 0.d0
   nfopt = 0
   nfbeg = 0
   iterrs = 0
   itdump = .false.
   iretry = 0
   gsqrd = 0.d0
   ginit = 0.d0
   gamden = 0.d0
   finit = 0.d0
   dfpr = 0.d0
   ddspln = 0.d0

   if (imin /= 5 .and. master) call amopen(7,mdinfo,'U','F','W')

#ifdef MPI
   ! correct fnq so that it is the same on both nodes
   if (icfe .ne. 0 .and. ifsc .ne. 0) then
     if( master ) then
       partner = ieor(masterrank,1)
       ti_nrp = nmixed + nsoftcore + nsoftcore_partner
       ti_nct(:) = 0
       if (ntc >= 2) then
         do i = 1,nbonh
           ti_ib = ix(iibh + i)/3 + 1
           ti_jb = ix(ijbh + i)/3 + 1
           if (nsc(ti_ib) .eq. 0 .and. nsc(ti_jb) .eq. 0) then
             ti_nct(1) = ti_nct(1) + 1
           else
             ti_nct(2) = ti_nct(2) + 1
           end if
         end do
       end if
       call mpi_sendrecv( ti_nct, 2, MPI_INTEGER, partner, 5, &
                          ti_nct_partner, 2, MPI_INTEGER, partner, 5, &
                          commmaster, ist, ierr )

       ti_nct_tot = ti_nct(1) + ti_nct(2) + ti_nct_partner(2)
       ndfp = 3*ti_nrp-ti_nct_tot
       fndfp = ndfp
       fnq = sqrt(fndfp)

     end if
!      if( numtasks>0 ) then
!         call mpi_bcast(fnq,1,MPI_DOUBLE_PRECISION,0,commsander,ierr)
!      end if  
   end if  
#endif

   rms = 0.0d0
   skip = .false.
   newstr = .false.
   ! determine the number of steepest descent steps
   steep = .false.
   nstcyc = 0
   mstcyc = kstcyc
   if(ntmin == 2) mstcyc = maxcyc
   if(ntmin == 1) mstcyc = ncyc
   if(ntmin > 0) steep = .true.

   ! Ben Roberts: Enable writing to a trajectory file if
   ! requested.
   ! Number of atoms to write to the trajectory
   ! If NTWPRT.NE.0, only print the atoms up to this value
   nrx = n
   if (ntwprt > 0) nrx = ntwprt*3
   ! Trajectory format
   loutfm = (ioutfm <= 0)


#ifdef MPI
   if( icfe /= 0 ) then
      allocate( frcti( n+3*extra_atoms ), stat = ier )
      REQUIRE( ier == 0 )
   end if
#endif

   fold = 0.0d0
   dxst = dx0
   linmin = 0
   gmin = 1.d0
   if (iscale > 0) n = n + iscale

   !     ----- PARTITION THE WORKING ARRAY -----

   irsdx = n
   irsdg = irsdx+n
   iginit = irsdg+n
   ixopt = iginit+n
   igopt = ixopt+n

   !     ----- SET SOME PARAMETERS TO BEGIN THE CALCULATION -----

   iterc = 0
   n_force_calls = 0
   iterfm = iterc
   
   !     ----- LET THE INITIAL SEARCH DIRECTION BE MINUS THE GRADIENT
   !           VECTOR. ITERRS GIVES THE ITERATION NUMBER OF THE MOST
   !           RECENT RESTART , BUT IS SET TO ZERO WHEN STEEPEST DESCENT
   !           DIRECTION IS USED -----

   !====================================================================
   !                    (Here is the beginning of a big loop:)
   20 continue
   !====================================================================

   !     ----- GATHER THE SUBMOLECULES INTO THE BOX -----
   n_force_calls = n_force_calls + 1
   if (mod(n_force_calls,nsnb) == 0) ntnb = 1
   if(ntnb == 1 .and. n_force_calls > 1) steep = .true.

   !====================================================================

   !     ----- CALCULATE THE FORCE AND ENERGY -----

   !     ----- APPLY SHAKE TO CONSTRAIN BONDS IF NECESSARY -----

   !====================================================================

   if(ntc /= 1) then
      fg(1:n) = x(1:n)
      nitp = 0
      qspatial=.false.
      call shake(nr,nbonh,nbona,0,ib,jb,igrp,winv,conp,skips, &
            fg,x,nitp,belly,ix(iifstwt),ix(noshake), &
            qspatial)
      call quick3(fg,x,ix(iifstwr),natom,nres,ix(i02))
      if(nitp <= 0) then
         ! shake failed
         ier = 135
         goto 290
      end if
   end if

   ! reset pb-related flags
   if(master)then
      if ( igb == 10 .or. ipb /= 0 ) then
         !if ( mod(n_force_calls,npbgrid) == 0 .and. n_force_calls /= maxcyc ) pbgrid = .true.
         !if ( mod(n_force_calls,ntpr) == 0 .or. n_force_calls == maxcyc ) pbprint = .true.
         !if ( mod(n_force_calls,nsnbr) == 0 .and. n_force_calls /= maxcyc ) ntnbr = 1
         !if ( mod(n_force_calls,nsnba) == 0 .and. n_force_calls /= maxcyc ) ntnba = 1
         if ( mod(n_force_calls,npbgrid) == 0 ) pbgrid = .true.
         if ( ntpr > 0 .and. mod(n_force_calls,ntpr) == 0 ) pbprint = .true.
         if ( mod(n_force_calls,nsnbr) == 0 ) ntnbr = 1
         if ( mod(n_force_calls,nsnba) == 0 ) ntnba = 1
         npbstep = n_force_calls
        !write(6,*) 'inside runmin', npbgrid, ntpr, nsnbr, nsnba
        !write(6,*) 'inside runmin', n_force_calls, pbgrid, pbprint, ntnbr, ntnba
      end if
   endif

   iprint = 0
   if (n_force_calls == maxcyc .or. n_force_calls == 1) iprint=1
   
   lout = .false.
   if (ntpr > 0) then
      lout = mod(n_force_calls, ntpr) == 0 .or. n_force_calls == 1
   else
      lout = n_force_calls == 1
   end if
   
   ! Ben Roberts: Switches to enable writing of restart files
   ! and trajectory coordinates for this particular step. Will
   ! set the flags to "false" unless the number of steps is right.
   ! Also requires that ntwr and ntwx are non-zero.
   
   ! Restart file
   ixdump = .false.
   if (ntwr /= 0 .and. mod(n_force_calls,ntwr) == 0) ixdump = .true.
   
   ! Trajectory
   ! DRR - Dont write traj during post-processing.
   itdump = .false.
   if (ntwx /= 0 .and. mod(n_force_calls,ntwx) == 0 .and. imin /= 5) itdump = .true.
   
   irespa = n_force_calls

   if (sebomd_obj%do_sebomd) then
     ! write down atomic charges and density matrix if needed
     sebomd_obj%iflagch = 0
     if (sebomd_obj%ntwc /= 0) then
        if (mod(n_force_calls,sebomd_obj%ntwc) == 0) sebomd_obj%iflagch = 1
     endif
!    sebomd_obj%pdmx = 0
!    if (sebomd_obj%pdump /= 0) then
!       if (mod(n_force_calls,ntwr) == 0) sebomd_obj%pdmx = 1
!       if (n_force_calls == maxcyc) sebomd_obj%pdmx = 1
!    endif
   endif
   
   call force(xx,ix,ih,ipairs,x,fg,ene,ene%vir, &
         xx(l96), xx(l97), xx(l98),xx(l99),qsetup, do_list_update,n_force_calls)

#ifdef MPI

   ! If softcore potentials are used, collect their dvdl contribution from the nodes
   ! this is disabled as long as multi-CPU TI minimizations dont work
   if ( ifsc /= 0 ) then
      !call mpi_reduce(sc_dvdl, sc_tot_dvdl, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, commsander, ierr)
      !call mpi_reduce(sc_dvdl_ee, sc_tot_dvdl_ee, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, commsander, ierr)
      sc_tot_dvdl = sc_dvdl
      sc_dvdl=0.0d0 ! zero for next step
      sc_tot_dvdl_ee = sc_dvdl_ee
      sc_dvdl_ee=0.0d0 ! zero for next step
      !call mpi_reduce(sc_ener, sc_ener_tmp, ti_ene_cnt, MPI_DOUBLE_PRECISION, MPI_SUM, 0, commsander, ierr)
      !sc_ener(1:ti_ene_cnt) = sc_ener_tmp(1:ti_ene_cnt)
   end if


   if ( icfe /= 0 )then

      ! ---free energies using thermodynamic integration (icfe /= 0)

      !  --- first, send the forces and energy to your partner:

      if( master ) then
         partner = ieor(masterrank,1)
         call mpi_sendrecv( fg, n, MPI_DOUBLE_PRECISION, partner, 5, &
                            frcti, n+3*extra_atoms, MPI_DOUBLE_PRECISION, partner, 5, &
                            commmaster, ist, ierr )
         call mpi_sendrecv( ene, state_rec_len, MPI_DOUBLE_PRECISION, partner, 5, &
                            ecopy, state_rec_len, MPI_DOUBLE_PRECISION, partner, 5, &
                            commmaster, ist, ierr)
         ! exchange sc-dvdl contributions between masters
         call mpi_sendrecv( sc_tot_dvdl, 1, MPI_DOUBLE_PRECISION, partner, 5, &
                            sc_tot_dvdl_partner, 1, MPI_DOUBLE_PRECISION, partner, 5, &
                            commmaster, ist, ierr )
         call mpi_sendrecv( sc_tot_dvdl_ee, 1, MPI_DOUBLE_PRECISION, partner, 5, &
                            sc_tot_dvdl_partner_ee, 1, MPI_DOUBLE_PRECISION, partner, 5, &
                            commmaster, ist, ierr )
         if( masterrank==0 ) then
            call mix_frcti(frcti,ecopy,fg,ene,n,clambda,klambda)
         else
            call mix_frcti(fg,ene,frcti,ecopy,n,clambda,klambda)
         end if
         ti_pot_ene = sc_ener(12)
         call mpi_sendrecv( ti_pot_ene, 1, MPI_DOUBLE_PRECISION, partner, 5, &
                            ti_pot_ene_partner, 1, MPI_DOUBLE_PRECISION, partner, 5, &
                            commmaster, ist, ierr )
      end if

!      if( numtasks>0 ) then
!         call mpi_bcast(fg,n,MPI_DOUBLE_PRECISION,0,commsander,ierr)
!         call mpi_bcast(ene,51,MPI_DOUBLE_PRECISION,0,commsander,ierr)
!      end if

   end if
#endif

   !f = ene(23)
   f = ene%pot%tot

#if defined( MPI ) 
   if( ievb /= 0 ) f = evb_frc%evb_nrg
#endif /* MPI */

   ntnb = 0
   sum = ddot(n,fg,1,fg,1)
#ifdef MPI /* SOFT CORE */
   if (ifsc == 1) then
      ! This stabilizes softcore minimizations by including the energy of
      ! the decoupled system. This is needed since the forces for the decoupled
      ! system are non-zero, so the corresponding energy should be included here.
      f = f + ti_pot_ene + ti_pot_ene_partner

      ! The forces have already been scaled by lambda, so rather than
      ! repeating the scaling with the sum we add the forces for the
      ! system including the decoupled partner
      if (master) then
        ti_sum = 0.d0     
        do i = 1, nr
          if (nsc(i) .ne. 0) then
             ti_sum = ti_sum + fg(3*(i-1)+1) * fg(3*(i-1)+1)
             ti_sum = ti_sum + fg(3*(i-1)+2) * fg(3*(i-1)+2)
             ti_sum = ti_sum + fg(3*(i-1)+3) * fg(3*(i-1)+3)
           end if
        end do
        call mpi_sendrecv( ti_sum, 1, MPI_DOUBLE_PRECISION, partner, 5, &
                           ti_sum_partner, 1, MPI_DOUBLE_PRECISION, partner, 5, &
                           commmaster, ist, ierr )       
        sum = sum + ti_sum_partner
      end if

      ! Now only called to handle round off error
      call sc_mix_sum(sum)

!      if( numtasks>0 ) then
!         call mpi_bcast(sum,1,MPI_DOUBLE_PRECISION,0,commsander,ierr)
!      end if

   end if
#endif

   if ((ifcap == 2 .or. ifcap == 5) .and. n_force_calls == 1) then  
      ! HG added this to account for waters not being considered if ifcap == 2,5
      fnq = fnq * fnq
      fnq = fnq - 3 * (outwat + oution)
      fnq = sqrt(fnq)
   end if          

   rms = sqrt(sum)/fnq

   !     ----- PRINT THE INTERMEDIATE RESULTS -----

!  x+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++x
!  |  Output EVB data                                              |
!  x+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++x
#ifdef MPI
   if( ievb /= 0 ) call out_evb ( n_force_calls)
#endif

   ! DAN ROE: modified so that during traj post-proc. only final results
   ! are printed.
   if (lout .and. imin /= 5) then
      call report_min_progress( n_force_calls, rms, fg, &
            ene, ih(m04), xx(l15) )  ! ih(m04) = atom names, xx(l15) = charges

#ifdef MPI /* SOFT CORE */
      if (ifsc /= 0) call sc_print_energies(6, sc_ener)
      if (ifsc /= 0) call sc_print_energies(7, sc_ener)
#endif

      !--- Print QM/MM Mulliken Charges if needed ---
      if (qmmm_nml%ifqnt) then
        if (qmmm_nml%printcharges .and. qmmm_mpi%commqmmm_master) then
          call qm2_print_charges(n_force_calls,qmmm_nml%dftb_chg,qmmm_struct%nquant_nlink, &
                                 qm2_struct%scf_mchg,qmmm_struct%iqm_atomic_numbers)
        end if
        if (qmmm_nml%printdipole /= 0) then
          call qmmm_dipole(x,xx(Lmass),ix(i02),ih(m02),nres)
        end if
      end if
!------
   end if

   !====================================================================

   !     ----- DO SOME STEEPEST STEPS BEFORE ENTERING THE CONJUGATE
   !           GRADIENT METHOD -----

   !====================================================================

   if (steep) then
      nstcyc = nstcyc+1
      if (nstcyc <= mstcyc) then

         if(dxst <= crits) dxst = dxstm
         dxst = dxst/2.0d0
         if(f < fold) dxst = dxst*2.4d0
         dxsth = dxst/sqrt(sum)
         if(nstcyc <= 1 .or. f <= fmin) then
            fmin = f
            nfopt = n_force_calls
            w(ixopt+1:ixopt+n) = x(1:n)
            w(igopt+1:igopt+n) = -fg(1:n)
         end if

         !     ----- CHECK FOR CONVERGENCE -----

         if (rms <= drms) then
            goto 300
         end if
         if (n_force_calls >= maxcyc) then
            ier = 131
            goto 290
         end if
         fold = f
         x(1:n) = x(1:n)+dxsth*fg(1:n)
 
         ! Ben Roberts: Write the restart file and/or the trajectory
         ! frame.
         if (master) then
            ! Write a restart file if appropriate
            if (ixdump) call minrit(n_force_calls,nrp,ntxo,x)
            ! Write to the trajectory if appropriate
            if (itdump) then
               call corpac(x,1,nrx,MDCRD_UNIT,loutfm)
               if (ioutfm > 0) call end_binary_frame(MDCRD_UNIT)
            end if
         end if
 
         goto 20

      else
         !                             (arrive here when finished with this
         !                              set of steepest descent cycles)
         steep = .false.
         newstr = .true.
         nstcyc = 0
         mstcyc = kstcyc
      end if
   end if

   !====================================================================

   !     ----- START OF CONJUGATE GRADIENT STEPS -----

   !====================================================================

   fg(1:n) = -fg(1:n)

   if (.not. newstr .and. n_force_calls > 1) goto 82
   70 continue
   w(1:n) = -fg(1:n)
   iterrs = 0
   if(newstr) iterc = 0
   if(iterc > 0) goto 140
   82 continue

   gnew = ddot(n,w,1,fg,1)
   !     ----- STORE THE VALUES OF X, F AND G, IF THEY ARE THE BEST THAT
   !           HAVE BEEN CALCULATED SO FAR. TEST FOR CONVERGENCE ----

   if (newstr .or. n_force_calls == 1) then
      ! artificially set fch less than zero to simplify the code.
      ! fmin will be properly initialized in the nested if statement below
      fch = -ONE
   else
      fch = f - fmin
   end if
   if (fch <= ZERO) then
      if (fch < ZERO .or. gnew/gmin >= -ONE) then
         fmin = f
         gsqrd = sum
         nfopt = n_force_calls
         w(ixopt+1:ixopt+n) = x(1:n)
         w(igopt+1:igopt+n) = fg(1:n)
      end if
      if (rms <= drms) then
         goto 300
      end if
   end if

   !     ----- TEST IF THE VALUE OF MAXCYC ALLOWS ANOTHER CALL OF FUNCT ---

   if (n_force_calls >= maxcyc) then
      ier = 131
      goto 290
   end if
   if (.not.newstr .and. n_force_calls > 1) goto 180

   !    ------ This section is executed at the beginning of a conjugate
   !           gradient set of minimization steps.

   !     ----- SET DFPR TO DX0*GSQRD. DFPR IS THE REDUCTION IN THE FUNCTION
   !           VALUE. STMIN IS USUALLY THE STEP-LENGTH OF THE MOST RECENT
   !           LINE SEARCH THAT GIVES THE LEAST VALUE OF F -----

   !  --- dac change, 10/91:  return to original idea of trying to
   !       go downhill by the absolute amount, DFPRED (which defaults
   !       to 1 kcal/mol, see data statement above).  This can eliminate
   !       very bad initial conjugate gradient steps.

   dfpr = dfpred
   stmin = dfpred/gsqrd

   newstr = .false.

   !====================================================================

   !     ----- Begin the main conjugate gradient iteration -----

   !====================================================================

   140 iterc = iterc+1

   finit = f
   ginit = 0.0d0
   w(iginit+1:iginit+n) = fg(1:n)
   ginit = ddot(n,w,1,fg,1)
   if(ginit >= 0.0d0) goto 260
   gmin = ginit
   sbound = -1.0d0
   nfbeg = n_force_calls
   iretry = -1

   stepch = min(stmin,abs(dfpr/ginit))
   stmin = dxstm

   160 step = stmin+stepch
   dxst = step
   swork = 0.0d0
   
   do i=1,n
      x(i) = w(ixopt+i)+stepch*w(i)
      swork = max(swork,abs(x(i)-w(ixopt+i)))
   end do
   
   ! Ben Roberts: Write the restart file and/or the trajectory frame
   ! if appropriate.
   if(swork > 0.0d0) then
      if (master) then
         ! Write a restart file if appropriate
         if (ixdump) call minrit(n_force_calls,nrp,ntxo,x)
         ! Write to the trajectory if appropriate
         if (itdump) then
         call corpac(x,1,nrx,MDCRD_UNIT,loutfm)
            if (ioutfm > 0) call end_binary_frame(MDCRD_UNIT)
         end if
      end if
      goto 20
   end if
   
   !     "work = swork" may not be needed - wont hurt.  -gls
   work = swork

   !     ----- TERMINATE THE LINE SEARCH IF STEPCH IS EFFECTIVELY ZERO ---

   if(n_force_calls > nfbeg+1 .or. abs(gmin/ginit) > 0.2d0) then
      if (master) write(6,370)
      steep = .true.
      linmin = linmin+1
   end if
   goto 270

   180 work = (fch+fch)/stepch-gnew-gmin
   ddspln = (gnew-gmin)/stepch
   if (n_force_calls > nfopt) then
      sbound = step
   else
      if(gmin*gnew <= 0.0d0) sbound = stmin
      stmin = step
      gmin = gnew
      stepch = -stepch
   end if
   if(fch /= 0.0d0) ddspln = ddspln+(work+work)/stepch

   !     ----- TEST FOR CONVERGENCE OF THE LINE SEARCH, BUT FORCE ATLEAST
   !           TWO STEPS TO BE TAKEN IN ORDER NOT TO LOSE QUADRATIC
   !           TERMINATION -----

   if(gmin == 0.0d0) goto 270
   if(n_force_calls <= nfbeg+1) goto 200
   if(abs(gmin/ginit) <= 0.2d0) goto 270

   !     ----- APPLY THE TEST THAT DEPENDS ON THE PARAMETER MAXLIN -----

   190 if(n_force_calls < nfopt+maxlin) goto 200

   !     ----- POSSIBLE NON BONDED UPDATE. MAKE A RESTART -----

   if (master) write(6,370)
   steep = .true.
   linmin = linmin+1
   goto 270

   200 stepch = 0.5d0*(sbound-stmin)
   if(sbound < -0.5d0) stepch = 9.0d0*stmin
   gspln = gmin+stepch*ddspln
   if(gmin*gspln < 0.0d0) stepch = stepch*gmin/(gmin-gspln)
   goto 160

   !     ----- CALCULATE THE VALUE OF BETAX IN THE NEW DIRECTION -----

   210 sum = ddot(n,fg,1,w(iginit+1),1)
   betax = (gsqrd-sum)/(gmin-ginit)

   !     ----- TEST THAT THE NEW SEARCH DIRECTION CAN BE MADE DOWNHILL.
   !           IF NOT THEN TRY TO IMPROVE THE ACCURACY OF THE LINE
   !           SEARCH -----

   if(abs(betax*gmin) <= 0.2d0*gsqrd) goto 220
   iretry = iretry+1
   if(iretry <= 0) goto 190

   220 if (f < finit) iterfm = iterc
   if (iterc >= iterfm+mxfcon) then
      if (master) write(6,370)
      steep = .true.
      linmin = linmin+1
      goto 270
   end if
   dfpr = stmin*ginit

   !     ----- BRANCH IF A RESTART PROCEDURE IS REQUIRED DUE TO THE
   !           ITERATION NUMBER OR DUE TO THE SCALAR PRODUCT OF
   !           CONSECUTIVE GRADIENTS -----

   if(iretry > 0) goto 70
   if(iterrs == 0) goto 240
   if(iterc-iterrs >= n) goto 240
   if(abs(sum) >= 0.2d0*gsqrd) goto 240

   !     ----- CALCULATE GAMMA IN THE NEW SEARCH DIRECTION. GAMDEN IS
   !           SET BY THE RESTART PROCEDURE -----

   gamma = ddot(n,fg,1,w(irsdg+1),1)
   sum  = ddot(n,fg,1,w(irsdx+1),1)
   gamma = gamma/gamden

   !     ----- RESTART IF THE NEW SEARCH DIRECTION IS NOT SUFFICIENTLY
   !           DOWNHILL ----

   if(abs(betax*gmin+gamma*sum) >= 0.2d0*gsqrd) goto 240

   !     ----- CALCULATE THE NEW SEARCH DIRECTION -----

   w(1:n) = -fg(1:n)+betax*w(1:n)+gamma*w(irsdx+1:irsdx+n)

   !    --- cycle back for more conjugate gradient steps:

   goto 140

   !     ----- APPLY THE RESTART PROCEDURE -----

   240 gamden = gmin-ginit
   do i=1,n
      w(irsdx+i) = w(i)
      w(irsdg+i) = fg(i)-w(iginit+i)
      w(i) = -fg(i)+betax*w(i)
   end do
   iterrs = iterc
   goto 140

   !     ----- SET IER TO INDICATE THAT THE SEARCH DIRECTION IS UPHILL ---

   260 continue
   steep = .true.
   if (master) write(6,370)
   linmin = linmin+1

   !     ----- ENSURE THAT F, X AND G ARE OPTIMAL -----

   270 continue
   if (n_force_calls /= nfopt) then
      f = fmin
      x(1:n)  = w(ixopt+1:ixopt+n)
      fg(1:n) = w(igopt+1:igopt+n)
   end if
   
   if (linmin > 4) then
      ier = 133
      goto 290
   end if
   
   ! Ben Roberts: Write the restart file and/or the trajectory frame
   ! if appropriate.
   if (steep) then
      if (master) then
         ! Write a restart file if appropriate
         if (ixdump) call minrit(n_force_calls,nrp,ntxo,x)
         ! Write to the trajectory if appropriate
         if (itdump) then
            call corpac(x,1,nrx,MDCRD_UNIT,loutfm)
            if (ioutfm > 0) call end_binary_frame(MDCRD_UNIT)
         end if
      end if
      goto 20
   end if
   
   if (ier == 0) goto 210

   290 continue  ! The unconverged minimization terminates 
   if (master) then
      select case ( ier )
      case ( 131 )
         write(6,'(//,a)') '  Maximum number of minimization cycles reached.'
      case ( 133 )
         write(6,'(/5x,a)') '***** REPEATED LINMIN FAILURE *****'
         write(6,'(/5x,a)') '***** SEE http://ambermd.org/Questions/linmin.html FOR MORE INFO *****'
      case ( 135 )
         write(6,'(/5x,a)') '***** ERROR: SHAKE FAILURE *****'
      case default
         ! invalid ier
         ASSERT( .false. )
      end select
   end if

   300 continue  ! The converged minimization terminates 

   if (sebomd_obj%do_sebomd) then
     ! write down atomic charges and density matrix if needed
     sebomd_obj%iflagch = 0
     if (sebomd_obj%ntwc /= 0) then
        if (mod(n_force_calls,sebomd_obj%ntwc) == 0) sebomd_obj%iflagch = 1
     endif
     sebomd_obj%pdmx = 0
     if (sebomd_obj%pdump /= 0) then
        if (mod(n_force_calls,ntwr) == 0) sebomd_obj%pdmx = 1
        if (n_force_calls == maxcyc) sebomd_obj%pdmx = 1
     endif
   endif

   !    do a final force call with iprint=1 to get proper nmr restaint printout:

   iprint=1
   call force(xx,ix,ih,ipairs,x,fg,ene,ene%vir, &
         xx(l96), xx(l97), xx(l98),xx(l99),qsetup, do_list_update,n_force_calls)
#ifdef MPI
   if ( icfe /= 0 )then
      sc_tot_dvdl = sc_dvdl
      sc_dvdl=0.0d0 ! zero for next step
      sc_tot_dvdl_ee = sc_dvdl_ee
      sc_dvdl_ee=0.0d0 ! zero for next step
      ! ---free energies using thermodynamic integration (icfe /= 0)
      !  --- first, send the forces and energy to your partner:
      if( master ) then
         partner = ieor(masterrank,1)
         call mpi_sendrecv( fg, n, MPI_DOUBLE_PRECISION, partner, 5, &
                            frcti, n+3*extra_atoms, MPI_DOUBLE_PRECISION, partner, 5, &
                            commmaster, ist, ierr )
         call mpi_sendrecv( ene, state_rec_len, MPI_DOUBLE_PRECISION, partner, 5, &
                            ecopy, state_rec_len, MPI_DOUBLE_PRECISION, partner, 5, &
                            commmaster, ist, ierr)
         ! exchange sc-dvdl contributions between masters
         call mpi_sendrecv( sc_tot_dvdl, 1, MPI_DOUBLE_PRECISION, partner, 5, &
                            sc_tot_dvdl_partner, 1, MPI_DOUBLE_PRECISION, partner, 5, &
                            commmaster, ist, ierr )
         call mpi_sendrecv( sc_tot_dvdl_ee, 1, MPI_DOUBLE_PRECISION, partner, 5, &
                            sc_tot_dvdl_partner_ee, 1, MPI_DOUBLE_PRECISION, partner, 5, &
                            commmaster, ist, ierr )
         if( masterrank==0 ) then
            call mix_frcti(frcti,ecopy,fg,ene,n,clambda,klambda)
         else
            call mix_frcti(fg,ene,frcti,ecopy,n,clambda,klambda)
         end if
         ti_pot_ene = sc_ener(12)
         call mpi_sendrecv( ti_pot_ene, 1, MPI_DOUBLE_PRECISION, partner, 5, &
                            ti_pot_ene_partner, 1, MPI_DOUBLE_PRECISION, partner, 5, &
                            commmaster, ist, ierr )
      end if
   end if
#endif



   !     ----- WRITE THE FINAL RESULTS -----
   
   ! Ben Roberts:Write restart file here, instead of (as before)
   ! in sander.f, so n_force_calls may be provided.
   if (master) then
      ! Write a restart file if appropriate
      call minrit(n_force_calls,nrp,ntxo,x)
      ! Write to the trajectory if appropriate
      if (itdump) then
         call corpac(x,1,nrx,MDCRD_UNIT,loutfm)
         if (ioutfm > 0) call end_binary_frame(MDCRD_UNIT)
      end if
   end if

   call report_min_results( n_force_calls, rms, x, &
         fg, ene, ih(m04), xx, ix, ih )  ! ih(m04) = atom names
   carrms = rms

#ifdef MPI /* SOFT CORE */
   if (ifsc /= 0) call sc_print_energies(6, sc_ener)
   if (ifsc /= 0) call sc_print_energies(7, sc_ener)
#endif

#ifdef MPI
   if( ievb /= 0 ) then
      call evb_dealloc 
   endif

   if( icfe /= 0 ) then
      deallocate( frcti, stat=ier )
      REQUIRE( ier == 0 )
   end if
#endif

   return

   370 format(/4x,' ... RESTARTED DUE TO LINMIN FAILURE ...')

end subroutine runmin

