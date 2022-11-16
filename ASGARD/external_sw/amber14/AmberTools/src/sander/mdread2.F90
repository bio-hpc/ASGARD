! This file is _not_ meant to be compiled.  It is meant to be included.
#ifdef API
#  define FATAL_ERROR ierr = 1; return
#  define DELAYED_ERROR ierr = 1; return
#else
#  define FATAL_ERROR call mexit(6, 1)
#  define DELAYED_ERROR inerr = 1
#endif /* API */

   use molecule, only: n_iwrap_mask_atoms, iwrap_mask_atoms
   use lmod_driver, only : LMOD_NTMIN_LMOD, LMOD_NTMIN_XMIN, write_lmod_namelist
   use decomp, only : jgroup, indx, irespw
   use findmask
#ifdef LES
   use genbornles, only: isnucat ! gbneck2nu: check if atom belongs to nuc or protein
#else
   use genborn, only: isnucat ! gbneck2nu: check if atom belongs to nuc or protein
#endif
   use constants, only : ZERO, ONE, TWO
   use parms, only: req
   use nbips, only: ips
   use amd_mod, only: iamd,EthreshD,alphaD,EthreshP,alphaP, &
        w_amd,EthreshD_w,alphaD_w,EthreshP_w,alphaP_w
   use nblist, only: a,b,c,alpha,beta,gamma,nbflag,skinnb,sphere,nbtell,cutoffnb
   use amoeba_mdin, only : iamoeba,beeman_integrator
   use amoeba_runmd, only : AM_RUNMD_get_coords
   use nose_hoover_module, only: nchain  ! APJ
   use neb_vars, only: last_neb_atom, ineb
   use cmd_vars, only: adiab_param
   use constantph, only: cnstphread, cnstph_zero, cph_igb, mccycles
   use file_io_dat
   use sander_lib, only: upper
#ifdef LES
   use les_data, only : lestyp, lesfac, nlesty, temp0les
#endif
#ifndef API
   use emap, only : temap, emap_options
   use emil_mod, only : emil_do_calc
   use sebomd_module, only: sebomd_obj, sebomd_write_info, sebomd_check_options
   use qmmm_module, only : qmmm_nml, qmmm_vsolv, qmmm_struct
   use qmmm_vsolv_module, only : print
   use linear_response, only : lrt_interval
   use pimd_vars, only: itimass
   use crg_reloc, only: ifcr, cropt, crcut, crskin, crprintcharges
#else
#  ifdef LES
   use qmmm_module, only : qmmm_nml
#  endif
#endif /* API */
   use pimd_vars, only: ipimd
   use qmmm_module, only: get_atomic_number, qm_gb
#ifdef APBS
   use apbs
#endif /* APBS */
#ifdef EMIL
   use mdin_emil_dat_mod, only : error_hdr, init_emil_dat
#endif /* EMIL */
#ifndef API
   use xray_interface_module, only: xray_active, xray_write_options
#endif /* API */
#if defined(LES) && defined(MPI)
   use evb_pimd, only: evb_pimd_init
#endif /* LES && MPI */
#ifdef MPI
   use softcore, only : ifsc, scalpha, scbeta, dvdl_norest, &
                        sceeorder, logdvdl, dynlmb
   use mbar, only : ifmbar, bar_intervall, bar_l_min, bar_l_max, bar_l_incr
! REMD
   use remd, only : rem, rremd
   use sgld, only : isgld ! for RXSGLD
#endif /* MPI */
   use linear_response, only: ilrt, lrtmask

   implicit none
   _REAL_ x(*)
#  include "../include/memory.h"
   integer ix(lasti)
   character(len=4) ih(*)
   integer ierr, nbond
   integer atom1,atom2
   logical belly,konst
   character(len=2) atype
   integer atomicnumber, hybridization
   integer ngrp,inerr,nr,ir,i,mxresat,j
   integer noshakegp( natom ), natnos
   integer iwrap_maskgp( natom ) , ier
   logical errFlag
   _REAL_ emtmd
#ifndef LES
   logical newstyle
#endif /* LES */
#ifndef API
   integer ntmp
#endif /* API */
   
#ifdef MPI
   !     =========================== AMBER/MPI ===========================
#  ifdef MPI_DOUBLE_PRECISION
#     undef MPI_DOUBLE_PRECISION
#  endif
   include 'mpif.h'
#  include "parallel.h"
   integer ist(MPI_STATUS_SIZE), partner, nbonh_c, num_noshake_c
   integer nquant_c, noshake_overlap_c
   integer crggp( natom )
#  ifdef CRAY_PVP
#     define MPI_DOUBLE_PRECISION MPI_REAL8
#  endif
   !     ========================= END AMBER/MPI =========================
#endif /* MPI */
#  include "../include/md.h"
#  include "box.h"
#  include "nmr.h"
#  include "extra_pts.h"
#  include "ew_cntrl.h"
#  include "ew_pme_recip.h"
#  include "ew_erfc_spline.h"
#  include "ew_mpole.h"
#  include "ew_legal.h"
#  include "def_time.h"
#  include "tgtmd.h"
#  include "multitmd.h"
   
   ! -------------------------------------------------------------------
   !      --- set up resat array, containing string identifying
   !          residue for each atom
   ! -------------------------------------------------------------------
   ierr = 0
   mxresat = min( natom, matom )
   ir = 0
   do i=1,mxresat
      if (i >= ix(ir+i02)) ir=ir+1
      write(resat(i),'(a4,1x,a4,i4)') ih(m04+i-1), &
            ih(ir+m02-1),ir
      !                                    ---null terminator:
      resat(i)(14:14) = char(0)
   end do
   close(unit=8)
   
   ! -------------------------------------------------------------------
   !     ----- SET THE DEFAULT VALUES FOR SOME VARIABLES -----
   ! -------------------------------------------------------------------
   
   nrp = natom
   
#ifndef API
   if (ifbox == 1) write(6, '(/5x,''BOX TYPE: RECTILINEAR'')')
   if (ifbox == 2) write(6, '(/5x,''BOX TYPE: TRUNCATED OCTAHEDRON'')')
   if (ifbox == 3) write(6, '(/5x,''BOX TYPE: GENERAL'')')
#endif
   
   ! For 0< ipimd <= 3, no removal of COM motion
   if (ipimd>0.and.ipimd<=3) then
      nscm = 0
      ndfmin = 0
   endif

   if (ntr.eq.1) &
      nscm = 0

   nsolut =  nrp
   if ( nscm > 0 .and. ntb == 0 ) then
      ndfmin = 6   ! both translation and rotation com motion removed
      if (nsolut == 1) ndfmin = 3
      if (nsolut == 2) ndfmin = 5
   else if ( nscm > 0 ) then
      ndfmin = 3    ! just translation com will be removed
   else
      ndfmin = 0
   end if
   if (ibelly > 0) then   ! No COM Motion Removal, ever.
      nscm = 0
      ndfmin = 0

      !Do not allow ntt=3 with ibelly=1 - this does not make sense and will cause
      !issues.
      if (ntt==3) then
        call sander_bomb("mdread2","ibelly=1 with ntt=3 is not a valid option.", &
                         "Either use a different thermostat or avoid using ibelly.")
      end if
   end if
   if(nscm <= 0) nscm = 0
   if(gamma_ln > 0.0d0)ndfmin=0  ! No COM motion removal for LD simulation
   if(ntt == 4)ndfmin=0  ! No COM motion removal for Nose'-Hoover simulation
   init = 3
   if (irest > 0) init = 4
   if (dielc <= ZERO ) dielc = ONE
   if (tautp <= ZERO ) tautp = 0.2d0
   if (taup <= ZERO ) taup = 0.2d0
   
   !     ----- RESET THE CAP IF NEEDED -----
   
   ! ivcap == 0: Cap will be in effect if it is in the prmtop file (ifcap = 1)

   if(ivcap == 1) then
      ! Cap will be in effect even if not in prmtop file
      !   requires additional information in sander.in file as in the case of ivcap == 3, 4, or 5
      ifcap = 2
   else if(ivcap == 2) then
      ! Inactivate cap
      ifcap = 0
   else if(ivcap == 3) then
      ! Sphere -> not yet implemented
      ifcap = 3
   else if(ivcap == 4) then
      ! Orthorhombus
      ifcap = 4
   else if(ivcap == 5) then
      ! Shell of waters around solute
      ifcap = 5
   end if

   !Support for random seed using time of day in microsecs
   if( ig==-1 ) then
     !Turn on NO_NTT3_SYNC when ig=-1. This means the code no
     !longer synchronized the random numbers between streams when
     !running in parallel giving better scaling.
     no_ntt3_sync = 1
     call microsec(ig)
#ifdef MPI
     write (6, '(a,i8,a)') "Note: ig = -1. Setting random seed to ", ig ," based on wallclock &
                               &time in microseconds"
     write (6, '(a)') "      and disabling the synchronization of random &
                               &numbers between tasks"
     write (6, '(a)') "      to improve performance."
#else 
#  ifndef API
     write (6, '(a,i8,a)') "Note: ig = -1. Setting random seed to ", ig ," based on wallclock &
                                &time in microseconds."
#  endif /* API */
#endif
   end if

#ifdef MPI
   !For some runs using multisander where the coordinates of the two 'replicas' need to be identical,
   !for example TI, it is critical that the random number stream is synchronized between all replicas.
   !Only use the IG value from worldrank=0. Ok to broadcast between the various sander masters since
   !they all call mdread2.
   !Also needs to be synchronized for adaptive QM/MM (qmmm_nml%vsolv > 1)
   if ( (icfe > 0) .or. (qmmm_nml%vsolv > 1) ) then
      ! no_ntt3_sync currently does not work with softcore TI simulations 
      ! see sc_lngdyn in softcore.F90
      ! AWG: I think I also need to set no_ntt3_sync=0 ???
      if ( (ifsc > 0) .or. (qmmm_nml%vsolv > 1) ) no_ntt3_sync = 0 
      call mpi_bcast(ig, 1, MPI_INTEGER, 0, commmaster, ierr)
   end if
#endif
 
   ! -------------------------------------------------------------------
   !     ----- PRINT DATA CHARACTERIZING THE RUN -----
   ! -------------------------------------------------------------------
   
   nr = nrp
#ifndef API
   write(6,9328)
   write(6,9008) title
   write(6,'(/a)') 'General flags:'
   write(6,'(5x,2(a,i8))') 'imin    =',imin,', nmropt  =',nmropt
#endif

   ! Error Checking for REMD
#ifdef MPI
   if (rem/=0) then
      ! Make sure that the number of replicas is even so
      !  that they all have partners in the exchange step
      if (mod(numgroup,2).ne.0) then
         write (6,'(a)') "==================================="
         write (6,'(a)') "REMD requires an even # of replicas"
         write (6,'(a)') "==================================="
         FATAL_ERROR
      endif

      write(6,'(/a)') 'Replica exchange'
      write(6,'(5x,4(a,i8))') 'numexchg=',numexchg,', rem=',rem

      ! REPCRD option temporarily disabled
      if (repcrd == 0) write(6,'(a)') &
         "REMD WARNING: repcrd disabled. Only replica &
         &trajectories/output can be written."

      ! Check for correct number of exchanges
      if (numexchg <= 0) then
         write(6,'(a)') "REMD ERROR: numexchg must be > 0, "
         FATAL_ERROR
      endif

      ! RXSGLD
      if (isgld > 0 .and. rem < 0) then
         write(6, '(a)') 'Multi-D REMD and replica-exchange SGLD are not &
                         &supported yet!'
         FATAL_ERROR
      end if

      ! Hybrid GB
      if (numwatkeep >= 0) then
         write(6,'(5x,4(a,i8))') 'numwatkeep=',numwatkeep,', hybridgb=',hybridgb
         ! Check that user specified GB model for hybrid REMD
         if (hybridgb /= 2 .and. hybridgb /= 5 .and. hybridgb /= 1) then
            write(6,'(a)') "HYBRID REMD ERROR: hybridgb must be 1, 2, or 5."
            FATAL_ERROR
         endif
      else
         !Check that user did not specify GB model if no hybrid run.
         if (hybridgb /= 0) then
            write(6,'(a)') &
            "HYBRID REMD ERROR: numwatkeep must be >= 0 when hybridgb is set."
            FATAL_ERROR
         endif
      endif
      ! RREMD
      if (rremd>0) then
         write(6,'(5x,4(a,i8))') "rremd=",rremd
      endif

      ! ntp > 0 not allowed for remd 
      if (ntp > 0) then
         write(6,'(a,i1)') "ERROR: ntp > 0 not allowed for rem > 0, ntp=", ntp
         FATAL_ERROR
      endif

      ! M-REMD (rem < 0) requires netcdf output.
      if (rem < 0 .and. ioutfm .ne. 1) then
         write(6,'(a)') "ERROR: Multi-D REMD (rem < 0) requires NetCDF &
                        &trajectories (ioutfm=1)"
         FATAL_ERROR
      endif

#  ifdef LES
      ! DAN ROE: Temporarily disable LES REMD until it is verified with new
      !           REMD code
      if (rem==2) then
         write (6,*) "******* LES REM (rem==2) temporarily disabled. Stop. *******"
         FATAL_ERROR
      endif

      if (rem==2 .and. igb/=1) then
         write (6,*) ' partial REM (rem=2) only works with igb=1'
         FATAL_ERROR
      endif
#  else
      if (rem==2) then
         write(6,*) '*******  For rem ==2, partial REM'
         write(6,*) 'use sander.LES with topology created by addles'
         FATAL_ERROR
      endif
#  endif

   endif ! rem>0

#else /* _NO_ MPI */
   ! Check if user set numexchg with no MPI
   if (numexchg>0) write(6,'(a)') &
      "WARNING: numexchg > 0 - for REMD run please recompile sander for &
      &parallel runs."

   ! Check if user set numwatkeep or hybridgb with no MPI - not sensible.
   if (numwatkeep>=0) write(6,'(a)') &
      "WARNING: numwatkeep >= 0 - for hybrid REMD run please recompile &
      &sander for parallel runs."

   if (hybridgb>0) write(6,'(a)') &
      "WARNING: hybridgb > 0 - for hybrid REMD run please recompile &
      &sander for parallel runs."
#endif /* MPI */
   ! End error checking for REMD

#ifndef API
   write(6,'(/a)') 'Nature and format of input:'
   write(6,'(5x,4(a,i8))') 'ntx     =',ntx,', irest   =',irest, &
         ', ntrx    =',ntrx

   write(6,'(/a)') 'Nature and format of output:'
   write(6,'(5x,4(a,i8))') 'ntxo    =',ntxo,', ntpr    =',ntpr, &
         ', ntrx    =',ntrx,', ntwr    =',ntwr
   write(6,'(5x,5(a,i8))') 'iwrap   =',iwrap,', ntwx    =',ntwx, &
         ', ntwv    =',ntwv,', ntwe    =',ntwe
   write(6,'(5x,3(a,i8),a,i7)') 'ioutfm  =',ioutfm, &
         ', ntwprt  =',ntwprt, &
         ', idecomp =',idecomp,', rbornstat=',rbornstat
   if (ntwf > 0) &
      write(6,'(5x, a,i8)') 'ntwf    =',ntwf
   write(6,'(/a)') 'Potential function:'
   write(6,'(5x,5(a,i8))') 'ntf     =',ntf,', ntb     =',ntb, &
         ', igb     =',igb,', nsnb    =',nsnb
   write(6,'(5x,3(a,i8))') 'ipol    =',ipol,', gbsa    =',gbsa, &
         ', iesp    =',iesp
   write(6,'(5x,3(a,f10.5))') 'dielc   =',dielc, &
         ', cut     =',cut,', intdiel =',intdiel

   ! charge relocation
   if ( ifcr /= 0 ) then
      write(6,'(/a)') 'Charge relocation:'
      write(6,'(5x,2(a,i8))') 'cropt   =', cropt, &
                                ', crprintcharges=', crprintcharges
      write(6,'(5x,2(a,f10.5))') 'crcut   =', crcut, ', crskin  =', crskin
   end if

   if (( igb /= 0 .and. igb /= 10 .and. ipb == 0 .and. igb /= 8) &
                                  .or.hybridgb>0.or.icnstph>1) then
      write(6,'(5x,3(a,f10.5))') 'saltcon =',saltcon, &
            ', offset  =',offset,', gbalpha= ',gbalpha
      write(6,'(5x,3(a,f10.5))') 'gbbeta  =',gbbeta, &
            ', gbgamma =',gbgamma,', surften =',surften
      write(6,'(5x,3(a,f10.5))') 'rdt     =',rdt, ', rgbmax  =',rgbmax, &
            '  extdiel =',extdiel
      write(6,'(5x,3(a,i8))') 'alpb  = ',alpb
   end if
   
   ! print output for igb=8
   if ( igb == 8 ) then
      write(6,'(5x,3(a,f10.5))') 'saltcon =',saltcon, &
            ', offset  =',offset,', surften =',surften
      write(6,'(5x,3(a,f10.5))') 'rdt     =',rdt, ', rgbmax  =',rgbmax, &
            '  extdiel =',extdiel
      write(6,'(5x,3(a,i8))') 'alpb  = ',alpb
      write(6,'(5x,3(a,f10.5))') 'gbalphaH  =',gbalphaH, &
            ', gbbetaH   =',gbbetaH,',  gbgammaH  = ',gbgammaH
      write(6,'(5x,3(a,f10.5))') 'gbalphaC  =',gbalphaC, &
            ', gbbetaC   =',gbbetaC,',  gbgammaC  = ',gbgammaC
      write(6,'(5x,3(a,f10.5))') 'gbalphaN  =',gbalphaN, &
            ', gbbetaN   =',gbbetaN,',  gbgammaN  = ',gbgammaN
      write(6,'(5x,3(a,f10.5))') 'gbalphaOS =',gbalphaOS, &
            ', gbbetaOS  =',gbbetaOS,',  gbgammaOS = ',gbgammaOS
      write(6,'(5x,3(a,f10.5))') 'gbalphaP  =',gbalphaP, &
            ', gbbetaP   =',gbbetaP,',  gbgammaP  = ',gbgammaP
      ! gbneck2nu pars
      write(6,'(5x,3(a,f10.5))') 'gb_alpha_hnu  =', gb_alpha_hnu, &
        ', gb_beta_hnu   =', gb_beta_hnu, ',  gb_gamma_hnu  = ', gb_gamma_hnu
      write(6,'(5x,3(a,f10.5))') 'gb_alpha_cnu  =', gb_alpha_cnu, &
        ', gb_beta_cnu   =', gb_beta_cnu, ',  gb_gamma_cnu  = ', gb_gamma_cnu
      write(6,'(5x,3(a,f10.5))') 'gb_alpha_nnu  =', gb_alpha_nnu, &
        ', gb_beta_nnu   =', gb_beta_nnu, ',  gb_gamma_nnu  = ', gb_gamma_nnu
      write(6,'(5x,3(a,f10.5))') 'gb_alpha_onu  =', gb_alpha_osnu, &
        ', gb_beta_onu   =', gb_beta_osnu, ',  gb_gamma_onu  = ', gb_gamma_osnu
      write(6,'(5x,3(a,f10.5))') 'gb_alpha_pnu  =', gb_alpha_pnu, &
        ', gb_beta_pnu   =', gb_beta_pnu, ',  gb_gamma_pnu  = ', gb_gamma_pnu
   end if


   if( alpb /= 0 ) then
      write(6,'(5x,3(a,f10.5))') 'Arad =', Arad
   end if    

   write(6,'(/a)') 'Frozen or restrained atoms:'
   write(6,'(5x,4(a,i8))') 'ibelly  =',ibelly,', ntr     =',ntr
   if( ntr == 1 ) write(6,'(5x,a,f10.5)') 'restraint_wt =', restraint_wt

   if( imin /= 0 ) then
      if( ipimd > 0 ) then
         write(6,'(/a)') 'pimd cannot be used in energy minimization'
         stop
      end if

      write(6,'(/a)') 'Energy minimization:'
      ! print inputable variables applicable to all minimization methods.
      write(6,'(5x,4(a,i8))') 'maxcyc  =',maxcyc,', ncyc    =',ncyc, &
            ', ntmin   =',ntmin
      write(6,'(5x,2(a,f10.5))') 'dx0     =',dx0, ', drms    =',drms

      ! Input flag ntmin determines the method of minimization
      select case ( ntmin )
      case ( 0, 1, 2 )
         ! no specific output
      case ( LMOD_NTMIN_XMIN, LMOD_NTMIN_LMOD )
         call write_lmod_namelist( )
      case default
         ! invalid ntmin
         write(6,'(/2x,a,i3,a)') 'Error: Invalid NTMIN (',ntmin,').'
         stop
      end select
   else
      write(6,'(/a)') 'Molecular dynamics:'
      write(6,'(5x,4(a,i10))') 'nstlim  =',nstlim,', nscm    =',nscm, &
            ', nrespa  =',nrespa
      write(6,'(5x,3(a,f10.5))') 't       =',t, &
            ', dt      =',dt,', vlimit  =',vlimit
      
      if ( ntt == 0 .and. tempi > 0.0d0 .and. irest == 0 ) then
         write(6,'(/a)') 'Initial temperature generation:'
         write(6,'(5x,a,i8)') 'ig      =',ig
         write(6,'(5x,a,f10.5)') 'tempi   =',tempi
      else if( ntt == 1 ) then
         write(6,'(/a)') 'Berendsen (weak-coupling) temperature regulation:'
         write(6,'(5x,3(a,f10.5))') 'temp0   =',temp0, &
               ', tempi   =',tempi,', tautp   =', tautp
#  ifdef LES
         write(6,'(5x,3(a,f10.5))') 'temp0LES   =',temp0les
#  endif
      else if( ntt == 2 ) then
         write(6,'(/a)') 'Anderson (strong collision) temperature regulation:'
         write(6,'(5x,4(a,i8))') 'ig      =',ig, ', vrand   =',vrand
         write(6,'(5x,3(a,f10.5))') 'temp0   =',temp0, ', tempi   =',tempi
      else if( ntt == 3 ) then
         write(6,'(/a)') 'Langevin dynamics temperature regulation:'
         write(6,'(5x,4(a,i8))') 'ig      =',ig
         write(6,'(5x,3(a,f10.5))') 'temp0   =',temp0, &
               ', tempi   =',tempi,', gamma_ln=', gamma_ln
      else if( ntt == 4 ) then
         write(6,'(/a)') 'Nose-Hoover chains'
         write(6,'(5x,(a,f10.5))') 'gamma_ln=', gamma_ln
         write(6,'(5x,(a,i8))') 'number of oscillators=', nchain
      else if( ntt == 5 ) then                                       ! APJ
         write(6,'(/a)') 'Nose-Hoover chains Langevin'               ! APJ
         write(6,'(5x,4(a,i8))') 'ig      =',ig                      ! APJ
         write(6,'(5x,(a,f10.5))') 'gamma_ln=', gamma_ln             ! APJ
         write(6,'(5x,(a,i8))') 'number of oscillators=', nchain     ! APJ
      else if( ntt == 6 ) then                                       ! APJ
         write(6,'(/a)') 'Adaptive Langevin temperature regulation:' ! APJ
         write(6,'(5x,4(a,i8))') 'ig      =',ig                      ! APJ
         write(6,'(5x,3(a,f10.5))') 'temp0   =',temp0, &             ! APJ
               ', tempi   =',tempi,', gamma_ln=', gamma_ln           ! APJ
      else if( ntt == 7 ) then                                       ! APJ
         write(6,'(/a)') 'Adaptive Nose-Hoover chains'               ! APJ
         write(6,'(5x,(a,f10.5))') 'gamma_ln=', gamma_ln             ! APJ
         write(6,'(5x,(a,i8))') 'number of oscillators=', nchain     ! APJ
      else if( ntt == 8 ) then                                       ! APJ
         write(6,'(/a)') 'Adaptive Nose-Hoover chains Langevin'      ! APJ
         write(6,'(5x,4(a,i8))') 'ig      =',ig                      ! APJ
         write(6,'(5x,(a,f10.5))') 'gamma_ln=', gamma_ln             ! APJ
         write(6,'(5x,(a,i8))') 'number of oscillators=', nchain     ! APJ
      else if( ntt == 9) then
         write(6,'(/a)') 'Canonical-isokinetic ensemble regulation:'
         write(6,'(5x,3(a,f10.5))') 'temp0   =',temp0, &
               ', tempi   =',tempi,', gamma_ln=', gamma_ln
         write(6,'(5x,2(a,i10),/)') 'nkija   =',nkija,', idistr  =',idistr
      end if

      if( ntp /= 0 ) then
         write(6,'(/a)') 'Pressure regulation:'
         write(6,'(5x,4(a,i8))') 'ntp     =',ntp
         write(6,'(5x,3(a,f10.5))') 'pres0   =',pres0, &
               ', comp    =',comp,', taup    =',taup
         if (barostat == 2) then
            write(6, '(5x,a)') 'Monte-Carlo Barostat:'
            write(6, '(5x,a,i8)') 'mcbarint  =', mcbarint
         end if
      end if

      if (csurften /= 0) then
         write(6,'(/a)') 'Constant surface tension:'
         write(6,'(5x,a,i8)') 'csurften  =', csurften
         write(6,'(5x,a,f10.5,a,i8)') 'gamma_ten =', gamma_ten, ' ninterface =', ninterface
      end if

   end if

   if( ntc /= 1 ) then
      write(6,'(/a)') 'SHAKE:'
      write(6,'(5x,4(a,i8))') 'ntc     =',ntc,', jfastw  =',jfastw
      write(6,'(5x,3(a,f10.5))') 'tol     =',tol
   end if

   if( ifcap == 1 .or. ifcap == 2 .or. ifcap == 3 ) then
      write(6,'(/a)') 'Water cap:'
      write(6,'(5x,2(a,i8))') 'ivcap   =',ivcap,', natcap  =',natcap
      write(6,'(5x,2(a,f10.5))') 'fcap    =',fcap, ', cutcap  =',cutcap
      write(6,'(5x,3(a,f10.5))') 'xcap    =',xcap, ', ycap    =',ycap,    &
                                 ', zcap    =',zcap
   else if( ifcap == 4 ) then
      write(6,'(/a)') 'Orthorhombus:'
      write(6,'(5x,1(a,i8))')    'ivcap   =',ivcap
      write(6,'(5x,1(a,f10.5))') 'forth   =',forth
      write(6,'(5x,3(a,f10.5))') 'xlorth  =',xlorth,', ylorth  =',ylorth, &
                                 ', zlorth  =',zlorth
      write(6,'(5x,3(a,f10.5))') 'xorth   =',xorth, ', yorth   =',yorth,  &
                                 ', zorth   =',zorth
   else if( ifcap == 5 ) then
      write(6,'(/a)') 'Water shell:'
      write(6,'(5x,(a,i8,a,f10.5))') 'ivcap   =',ivcap,', cutcap  =',cutcap
   endif

   if( nmropt > 0 ) then
      write(6,'(/a)') 'NMR refinement options:'
      write(6,'(5x,4(a,i8))')'iscale  =',iscale,', noeskp  =',noeskp, &
            ', ipnlty  =',ipnlty,', mxsub   =',mxsub
      write(6,'(5x,3(a,f10.5))') 'scalm   =',scalm, &
            ', pencut  =',pencut,', tausw   =',tausw
   end if

   if( numextra > 0 ) then
      write(6,'(/a)') 'Extra-points options:'
      write(6,'(5x,4(a,i8))') 'frameon =',frameon, &
            ', chngmask=',chngmask
   end if

   if( ipol > 0 ) then
      write(6,'(/a)') 'Polarizable options:'
      write(6,'(5x,4(a,i8))') 'indmeth =',indmeth, &
            ', maxiter =',maxiter,', irstdip =',irstdip, &
            ', scaldip =',scaldip
      write(6,'(5x,3(a,f10.5))') &
            'diptau  =',diptau,', dipmass =',dipmass
     if ( ipol > 1 ) then
      write(6,'(5x,3(a,f10.5))') &
            'Default Thole coefficient = ',dipdamp
     end if
   end if

#  ifdef MPI /* SOFT CORE */
   if( icfe /= 0 .or. ifsc/=0) then
      write(6,'(/a)') 'Free energy options:'
      write(6,'(5x,3(a,i8))')       'icfe    =', icfe   , ', ifsc    =', ifsc,    ', klambda =', klambda
      write(6,'(5x,3(a,f8.4))')     'clambda =', clambda, ', scalpha =', scalpha, ', scbeta  =', scbeta
      write(6,'(5x,2(a,i8))')       'sceeorder =', sceeorder, ' dvdl_norest =', dvdl_norest
      write(6,'(5x,a,f8.4,a,i8)')   'dynlmb =', dynlmb, ' logdvdl =', logdvdl
   end if
   if ( ifmbar /= 0 ) then
      write (6,'(/a)') 'FEP MBAR options:'
      write(6,'(5x,a,i8,a,i8)') 'ifmbar  =', ifmbar, ',  bar_intervall = ', bar_intervall
      write(6,'(5x,3(a,f6.4))') 'bar_l_min =', bar_l_min, ',  bar_l_max =', bar_l_max, &
                                ',  bar_l_incr =', bar_l_incr
   end if
#  endif /* MPI */

   if (ilrt /= 0) then
      write (6,*)
      write (6,'(a,i4,a,i4)') ' Linear Response Theory: ilrt =', ilrt, ' lrt_interval =', lrt_interval
      write (6,*)
   end if

   ! Options for TI w.r.t. mass.
   select case (itimass)
   case (0)     ! Default: no TI wrt. mass.
   case (1,2)   ! 1 = use virial est., 2 = use thermodynamic est.
      write(6,'(/a)') 'Isotope effects (thermodynamic integration w.r.t. mass):'
      write(6,'(5x,4(a,i8))') 'itimass =',itimass      
      write(6,'(5x,3(a,f10.5))') 'clambda =',clambda
      if (icfe /= 0) then 
         write(6,'(/2x,a,i2,a,i2,a)') 'Error: Cannot do TI w.r.t. both potential (icfe =', &
            icfe, ') and mass (itimass =', itimass, ').'
         stop
      endif
      if (ipimd == 0) then 
         write(6,'(/2x,a)') 'Error (IPIMD=0): TI w.r.t. mass requires a PIMD run.'
         stop       
      endif
   case default ! Invalid itimass
      write(6,'(/2x,a,i2,a)') 'Error: Invalid ITIMASS (', itimass, ' ).'
      stop
   end select

   if( itgtmd == 1 ) then
      write(6,'(/a)') 'Targeted molecular dynamics:'
      write(6,'(5x,3(a,f10.5))') 'tgtrmsd =',tgtrmsd, &
            ', tgtmdfrc=',tgtmdfrc
   end if
 
   if( icnstph /= 0) then
      write(6, '(/a)') 'Constant pH options:'
      if ( icnstph .ne. 1 ) &
         write(6, '(5x,a,i8)') 'icnstph =', icnstph
      write(6, '(5x,a,i8)') 'ntcnstph =', ntcnstph
      write(6, '(5x,a,f10.5)') 'solvph =', solvph
      if ( icnstph .ne. 1 ) &
         write(6,'(5x,2(a,i8))') 'ntrelax =', ntrelax, ' mccycles =', mccycles
   end if
   if( icnstph /= 2) then
      ntrelax = 0 ! needed for proper behavior of timing
   end if

   if( ntb > 0 ) then
      write(6,'(/a)') 'Ewald parameters:'
      write(6,'(5x,4(a,i8))') 'verbose =',verbose, &
            ', ew_type =',ew_type,', nbflag  =',nbflag, &
            ', use_pme =',use_pme
      write(6,'(5x,4(a,i8))') 'vdwmeth =',vdwmeth, &
            ', eedmeth =',eedmeth,', netfrc  =',netfrc
      write(6, 9002) a, b, c
      write(6, 9003) alpha, beta, gamma
      write(6, 9004) nfft1, nfft2, nfft3
      write(6, 9006) cutoffnb, dsum_tol
      write(6, 9007) ew_coeff
      write(6, 9005) order
      9002 format (5x,'Box X =',f9.3,3x,'Box Y =',f9.3,3x,'Box Z =',f9.3)
      9003 format (5x,'Alpha =',f9.3,3x,'Beta  =',f9.3,3x,'Gamma =',f9.3)
      9004 format (5x,'NFFT1 =',i5  ,7x,'NFFT2 =',i5  ,7x,'NFFT3 =',i5)
      9005 format (5x,'Interpolation order =',i5)
      9006 format (5x,'Cutoff=',f9.3,3x,'Tol   =',e9.3)
      9007 format (5x,'Ewald Coefficient =',f9.5)
   end if

!---- QMMM Options ----

   if( qmmm_nml%ifqnt ) then
      write(6, '(/a)') 'QMMM options:'
      write(6, '(5x,"        ifqnt = True       nquant = ",i8)') &
                 qmmm_struct%nquant
      write(6, '(5x,"         qmgb = ",i8,"  qmcharge = ",i8,"   adjust_q = ",i8)') &
                 qmmm_nml%qmgb, qmmm_nml%qmcharge, qmmm_nml%adjust_q
      write(6, '(5x,"         spin = ",i8,"     qmcut = ",f8.4, "    qmshake = ",i8)') qmmm_nml%spin, &
                 qmmm_nml%qmcut, qmmm_nml%qmshake
      write(6, '(5x,"     qmmm_int = ",i8)') qmmm_nml%qmmm_int
      write(6, '(5x,"lnk_atomic_no = ",i8,"   lnk_dis = ",f8.4," lnk_method = ",i8)') &
                 qmmm_nml%lnk_atomic_no,qmmm_nml%lnk_dis, qmmm_nml%lnk_method
      if ( qmmm_nml%qmtheory%PM3 ) then 
         write(6, '(5x,"     qm_theory =     PM3")',ADVANCE='NO') 
      else if ( qmmm_nml%qmtheory%AM1 ) then
         write(6, '(5x,"     qm_theory =     AM1")',ADVANCE='NO') 
      else if ( qmmm_nml%qmtheory%AM1D ) then
         write(6, '(5x,"     qm_theory =   AM1/d")',ADVANCE='NO') 
      else if ( qmmm_nml%qmtheory%MNDO ) then
         write(6, '(5x,"     qm_theory =    MNDO")',ADVANCE='NO') 
      else if ( qmmm_nml%qmtheory%MNDOD ) then
         write(6, '(5x,"     qm_theory =  MNDO/d")',ADVANCE='NO') 
      else if ( qmmm_nml%qmtheory%PDDGPM3 ) then
         write(6, '(5x,"     qm_theory = PDDGPM3")',ADVANCE='NO') 
      else if ( qmmm_nml%qmtheory%PDDGMNDO ) then
         write(6, '(5x,"     qm_theory =PDDGMNDO")',ADVANCE='NO') 
      else if ( qmmm_nml%qmtheory%PM3CARB1 ) then
         write(6, '(5x,"     qm_theory =PM3CARB1")',ADVANCE='NO') 
      else if ( qmmm_nml%qmtheory%DFTB ) then
         write(6, '(5x,"     qm_theory =    DFTB")',ADVANCE='NO') 
      else if ( qmmm_nml%qmtheory%RM1 ) then
         write(6, '(5x,"     qm_theory =     RM1")',ADVANCE='NO') 
      else if ( qmmm_nml%qmtheory%PDDGPM3_08 ) then
         write(6, '(5x,"     qm_theory = PDDGPM3_08")',ADVANCE='NO') 
      else if ( qmmm_nml%qmtheory%PM6 ) then
         write(6, '(5x,"     qm_theory =     PM6")',ADVANCE='NO') 
      else if ( qmmm_nml%qmtheory%PM3ZNB ) then
         write(6, '(5x,"     qm_theory = PM3/ZnB")',ADVANCE='NO') 
      else if ( qmmm_nml%qmtheory%PM3MAIS ) then
         write(6, '(5x,"     qm_theory = PM3-MAIS")',ADVANCE='NO') 
      else if ( qmmm_nml%qmtheory%EXTERN ) then
         write(6, '(5x,"     qm_theory =     EXTERN")',ADVANCE='NO') 
      else
         write(6, '(5x,"     qm_theory = UNKNOWN!")',ADVANCE='NO') 
      end if
      write (6, '(" verbosity = ",i8)') qmmm_nml%verbosity
      if (qmmm_nml%qmqm_analyt) then
        write(6, '(5x,"       qmqmdx = Analytical")')
      else
        write(6, '(5x,"       qmqmdx = Numerical")')
      end if
      !AWG: if EXTERN in use, skip printing of options that do not apply
      EXTERN: if ( .not. qmmm_nml%qmtheory%EXTERN ) then
         if (qmmm_nml%tight_p_conv) then
            write(6, '(5x," tight_p_conv = True (converge density to SCFCRT)")')
         else
            write(6, '(5x," tight_p_conv = False (converge density to 0.05xSqrt[SCFCRT])")')
         end if
         write(6, '(5x,"      scfconv = ",e9.3,"  itrmax = ",i8)') qmmm_nml%scfconv, qmmm_nml%itrmax
         if (qmmm_nml%printcharges) then
            write(6, '(5x," printcharges = True ")',ADVANCE='NO')
         else
            write(6, '(5x," printcharges = False")',ADVANCE='NO')
         end if
         select case (qmmm_nml%printdipole)
            case (1)
               write(6, '(5x," printdipole = QM   ")',ADVANCE='NO')
            case (2)
               write(6, '(5x," printdipole = QM+MM")',ADVANCE='NO')
            case default
               write(6, '(5x," printdipole = False")',ADVANCE='NO')
         end select 
         if (qmmm_nml%peptide_corr) then
            write(6, '(5x," peptide_corr = True")')
         else
            write(6, '(5x," peptide_corr = False")')
         end if
         if (qmmm_nml%qmmmrij_incore) then
            write(6, '(4x,"qmmmrij_incore = True ")')
         else
            write(6, '(4x,"qmmmrij_incore = False")')
         end if
         if (qmmm_nml%qmqm_erep_incore) then
            write(6, '(2x,"qmqm_erep_incore = True ")')
         else
            write(6, '(2x,"qmqm_erep_incore = False")')
         end if
         if (qmmm_nml%allow_pseudo_diag) then
            write(6, '(7x,"pseudo_diag = True ")',ADVANCE='NO')
            write(6, '("pseudo_diag_criteria = ",f8.4)') qmmm_nml%pseudo_diag_criteria
         else
            write(6, '(7x,"pseudo_diag = False")')
         end if
         write(6,   '(6x,"diag_routine = ",i8)') qmmm_nml%diag_routine
      end if EXTERN
      !If ntb=0 or use_pme =0 then we can't do qm_ewald so overide what the user may
      !have put in the namelist and set the value to false.
      if (qmmm_nml%qm_ewald>0) then
        if (qmmm_nml%qm_pme) then
          write(6, '(10x,"qm_ewald = ",i8, " qm_pme = True ")') qmmm_nml%qm_ewald
        else
          write(6, '(10x,"qm_ewald = ",i8, " qm_pme = False ")') qmmm_nml%qm_ewald
        end if
          write(6, '(10x,"  kmaxqx = ",i4," kmaxqy = ",i4," kmaxqz = ",i4," ksqmaxq = ",i4)') &
                  qmmm_nml%kmaxqx, qmmm_nml%kmaxqy, qmmm_nml%kmaxqz, qmmm_nml%ksqmaxq
      else
        write(6, '(10x,"qm_ewald = ",i8, " qm_pme = False ")') qmmm_nml%qm_ewald
      end if
      !Print the fock matrix prediction params if it is in use.
      if (qmmm_nml%fock_predict>0) then
          write(6, '(6x,"fock_predict = ",i4)') qmmm_nml%fock_predict
          write(6, '(6x,"    fockp_d1 = ",f8.4," fockp_d2 = ",f8.4)') qmmm_nml%fockp_d1, qmmm_nml%fockp_d2
          write(6, '(6x,"    fockp_d2 = ",f8.4," fockp_d4 = ",f8.4)') qmmm_nml%fockp_d3, qmmm_nml%fockp_d4
      end if

      if (qmmm_nml%qmmm_switch) then
         write(6, '(7x,"qmmm_switch = True",3x,"r_switch_lo =",f8.4,3x,"r_switch_hi =",f8.4)') &
               & qmmm_nml%r_switch_lo, qmmm_nml%r_switch_hi
      !else
      !   write(6, '(7x,"qmmm_switch = False")') 
      end if

      if (qmmm_nml%printdipole==2) then
         write(6, '("|",2x,"INFO: To compute MM dipole WAT residues will be stripped")') 
      end if 
   end if

   if (qmmm_nml%vsolv > 0) then
      call print(qmmm_vsolv)
   end if

!---- SEBOMD options ----

   if (sebomd_obj%do_sebomd) then
     call sebomd_write_info()
     call sebomd_check_options()
   endif

!---- XRAY Options ----
   if( xray_active ) call xray_write_options()

! ----EMAP Options-----
   if(temap)call emap_options(5)
! ---------------------

#  ifdef EMIL
! ----EMIL Options-----
   if(emil_do_calc.gt.0)call init_emil_dat(5, 6)
! ---------------------
#  endif


#  ifdef MPI
! --- MPI TIMING OPTIONS ---
      write(6, '(/a)') '| MPI Timing options:'
      write(6, '("|",5x," profile_mpi = ",i8)') profile_mpi
! Sanity check for profile_mpi
      call int_legal_range('profile_mpi',profile_mpi,0,1)
! --------------------------
#  endif /* MPI */
#endif /* API */
 
   cut = cut*cut
   cut_inner = cut_inner*cut_inner
   
   !------------------------------------------------------------------------
   ! If user has requested generalized born electrostatics, set up variables
   !------------------------------------------------------------------------
   
   if( igb == 0 .and. gbsa > 0 ) then
      write(0,*) 'GB/SA calculation is performed only when igb>0'
      FATAL_ERROR
   end if
   if( gbsa == 2 .and. &
       ((imin == 0 .and. nstlim > 1) .or. &
        (imin == 1 .and. maxcyc > 1)) ) then
      write(0,*) 'GBSA=2 only works for single point energy calc'
      FATAL_ERROR
   end if
#ifdef APBS
   if( igb /= 0 .and. igb /= 10 .and. ipb == 0 .and. .not. mdin_apbs) then
#else
   if (( igb /= 0 .and. igb /= 10 .and. ipb == 0 ).or.hybridgb>0.or.icnstph>1) then
#endif /* APBS */
#if defined (LES) && !defined (API)
      write(6,*) 'igb=1,5,7 are working with LES, no SA term included'
#endif
      ! igb7 uses special S_x screening params.
      ! overwrite the tinker values read from the prmtop
      if (igb == 7) then
         do i=1,natom
            if(ix(i100) .eq. 1) then
               atomicnumber = ix(i100+i)
            else
               call get_atomic_number(ih(m04+i-1), x(lmass+i-1), &
                  atomicnumber, errFlag)
#ifdef LES
               if(errFlag) then
                  call get_atomic_number(ih(m04+i-1), &
                     x(lmass+i-1)*lesfac((lestyp(i)-1)*nlesty+lestyp(i)), &
                     atomicnumber, errFlag)
               end if
#endif
            end if

            if (atomicnumber .eq. 6) then
               x(l96+i-1) = 4.84353823306d-1
            else if (atomicnumber .eq. 1) then 
               x(l96+i-1) = 1.09085413633d0
            else if (atomicnumber .eq. 7) then
               x(l96+i-1) = 7.00147318409d-1
            else if (atomicnumber .eq. 8) then
               x(l96+i-1) = 1.06557401132d0
            else if (atomicnumber .eq. 16) then
               x(l96+i-1) = 6.02256336067d-1
            else if (atomicnumber .eq. 15) then
               x(l96+i-1) = 5d-1
            else
               x(l96+i-1) = 5d-1
            end if
         end do
      end if
      
      ! Hai Nguyen: changing S_x screening params for igb = 8
      ! overwrite the tinker values read from the prmtop
      
      if (igb == 8) then
         do i=1,natom
            if (ix(i100) .eq. 1) then
               atomicnumber = ix(i100+i)
            else
               call get_atomic_number(ih(m04+i-1), x(lmass+i-1), &
                  atomicnumber, errFlag)
#ifdef LES
               if(errFlag) then
                  call get_atomic_number(ih(m04+i-1), &
                     x(lmass+i-1)*lesfac((lestyp(i)-1)*nlesty+lestyp(i)), &
                     atomicnumber, errFlag)
               end if
#endif
            end if

            call isnucat(nucat,i,nres,37,ix(i02:i02+nres-1),ih(m02:m02+nres-1)) 
            !Hai Nguyen: check if atom belong to nucleic or protein residue 
            !37 = the total number of different nucleic acid residue type in AMBER 
            !Need to update this number if adding new residue type 
            !update isnucat function in gb_ene.F90 too 
            if (nucat == 1) then 
                ! If atom belongs to nuc, use gbneck2nu pars 
                if (atomicnumber .eq. 6) then 
                  x(l96+i-1) = screen_cnu 
                else if (atomicnumber .eq. 1) then 
                  x(l96+i-1) = screen_hnu 
                else if (atomicnumber .eq. 7) then 
                  x(l96+i-1) = screen_nnu 
                else if (atomicnumber .eq. 8) then 
                  x(l96+i-1) = screen_onu
                else if (atomicnumber .eq. 15) then
                  x(l96+i-1) = screen_pnu
                else
                  x(l96+i-1) = 5.d-1 ! not optimized
                end if
            else
               ! protein
               if (atomicnumber .eq. 6) then
                  x(l96+i-1) = Sc
               else if (atomicnumber .eq. 1) then
                  x(l96+i-1) = Sh
               else if (atomicnumber .eq. 7) then
                  x(l96+i-1) = Sn
               else if (atomicnumber .eq. 8) then
                  x(l96+i-1) = So
               else if (atomicnumber .eq. 16) then
                  x(l96+i-1) = Ss 
               else if (atomicnumber .eq. 15) then
                  x(l96+i-1) = Sp ! We still don't have an optimized Sp parameter
                                  ! if nuc, use gbneck2nu pars 
               else
                  !for atom type Cl,Br,...
                  !These parameters are also not optimized.
                  x(l96+i-1) = 5d-1
               end if 
            end if !check isnucat
         end do !loop i=1, natom
      end if ! igb = 8

      ! Set up for igb == 2, 5, 7, 8
      ! Put gb parameters in arrays 
      if ( igb == 2 .or. igb == 5 .or. igb == 7 .or. &
           hybridgb == 2 .or. hybridgb == 5) then
        do i=1,natom
            x(l2402+i-1) = gbalpha
            x(l2403+i-1) = gbbeta
            x(l2404+i-1) = gbgamma
        end do
      end if
      
      ! IGB = 8, update gb_alpha, gb_beta, gb_gamma
       if ( igb == 8 ) then
           do i=1,natom
              if (ix(i100) .eq. 1) then
                 atomicnumber = ix(i100+i)
              else
                 call get_atomic_number(ih(m04+i-1), x(lmass+i-1), &
                    atomicnumber, errFlag)
#ifdef LES
                 if(errFlag) then
                    call get_atomic_number(ih(m04+i-1), &
                       x(lmass+i-1)*lesfac((lestyp(i)-1)*nlesty+lestyp(i)), &
                       atomicnumber, errFlag)
                 end if
#endif
              end if

          ! check if atom belongs to protein or nucleic acid
          call isnucat(nucat,i,nres,37,ix(i02:i02+nres-1),ih(m02:m02+nres-1))
          if (nucat == 1) then
             !if atom belong to nucleic part, use nuc pars
             if (atomicnumber .eq. 1) then
               x(l2402+i-1) = gb_alpha_hnu
               x(l2403+i-1) = gb_beta_hnu
               x(l2404+i-1) = gb_gamma_hnu
             else if (atomicnumber .eq. 6) then
               x(l2402+i-1) = gb_alpha_cnu
               x(l2403+i-1) = gb_beta_cnu
               x(l2404+i-1) = gb_gamma_cnu
             else if (atomicnumber .eq. 7) then
               x(l2402+i-1) = gb_alpha_nnu
               x(l2403+i-1) = gb_beta_nnu
               x(l2404+i-1) = gb_gamma_nnu
             else if (atomicnumber .eq. 8) then
               x(l2402+i-1) = gb_alpha_osnu
               x(l2403+i-1) = gb_beta_osnu
               x(l2404+i-1) = gb_gamma_osnu
             else if (atomicnumber .eq. 16) then
               x(l2402+i-1) = gb_alpha_osnu
               x(l2403+i-1) = gb_beta_osnu
               x(l2404+i-1) = gb_gamma_osnu
             else if (atomicnumber .eq. 15) then
               x(l2402+i-1) = gb_alpha_pnu
               x(l2403+i-1) = gb_beta_pnu
               x(l2404+i-1) = gb_gamma_pnu
             else
               ! Use GB^OBC (II) (igb = 5) set for other atoms
               x(l2402+i-1) = 1.0d0
               x(l2403+i-1) = 0.8d0
               x(l2404+i-1) = 4.851d0
             end if ! assign nucleic acid pars

          else
             !not nucleic part, use protein pars 
             if (atomicnumber .eq. 1) then
                     x(l2402+i-1) = gbalphaH
                     x(l2403+i-1) = gbbetaH
                     x(l2404+i-1) = gbgammaH
             else if (atomicnumber .eq. 6) then
                     x(l2402+i-1) = gbalphaC
                     x(l2403+i-1) = gbbetaC
                     x(l2404+i-1) = gbgammaC
             else if (atomicnumber .eq. 7) then
                     x(l2402+i-1) = gbalphaN
                     x(l2403+i-1) = gbbetaN
                     x(l2404+i-1) = gbgammaN
             else if (atomicnumber .eq. 8) then
                     x(l2402+i-1) = gbalphaOS
                     x(l2403+i-1) = gbbetaOS
                     x(l2404+i-1) = gbgammaOS
             else if (atomicnumber .eq. 16) then
                     x(l2402+i-1) = gbalphaOS
                     x(l2403+i-1) = gbbetaOS
                     x(l2404+i-1) = gbgammaOS
             else if (atomicnumber .eq. 15) then
                     x(l2402+i-1) = gbalphaP
                     x(l2403+i-1) = gbbetaP
                     x(l2404+i-1) = gbgammaP
             else
                     !use GBOBC set for other atom types
                     x(l2402+i-1) = 1.0d0
                     x(l2403+i-1) = 0.8d0
                     x(l2404+i-1) = 4.851d0
             end if ! assign protein pars
          end if !check nucat
        end do !loop i=1,natom
       end if ! end IGB = 8 section
        
     !       put fs(i)*(rborn(i) - offset) into the "fs" array
      fsmax = 0.d0
      do i=1,natom
         x(l96-1+i) = x(l96-1+i)*( x(l97-1+i) - offset )
         fsmax = max( fsmax, x(l96-1+i) )
         if (rbornstat == 1) then
            x(l186-1+i) = 0.d0
            x(l187-1+i) = 999.d0
            x(l188-1+i) = 0.d0
            x(l189-1+i) = 0.d0
         end if
      end do
      
      !     ---------------------------------------------------------------------
      !       ---get Debye-Huckel kappa (A**-1) from salt concentration (M), assuming:
      !         T = 298.15, epsext=78.5,
      
      kappa = sqrt( 0.10806d0 * saltcon )
      
      !       ---scale kappa by 0.73 to account(?) for lack of ion exclusions:
      
      kappa = 0.73d0* kappa

      !Set kappa for qmmm if needed
      qm_gb%kappa = kappa
      !     ---------------------------------------------------------------------
      
      if ( gbsa == 1 ) then
         
         !     --- assign parameters for calculating SASA according to the
         !         LCPO method ---
         
         do i=1,natom
            ix(i80+i-1)=0
         end do
         
         !         --- get the number of bonded neighbors for each atom:
         
         do i=1,nbona
            atom1=ix(iiba+i-1)/3+1
            atom2=ix(ijba+i-1)/3+1
            ix(i80+atom1-1)=ix(i80+atom1-1)+1
            ix(i80+atom2-1)=ix(i80+atom2-1)+1
         end do

         do i=1,natom
            ix(i80-i)=ix(i80+i-1)
         end do

         do i=1,nbonh
            atom1=ix(iibh+i-1)/3+1
            atom2=ix(ijbh+i-1)/3+1
            ix(i80-atom1)=ix(i80-atom1)+1
            ix(i80-atom2)=ix(i80-atom2)+1
         end do

         
         !         --- construct parameters for SA calculation; note that the
         !             radii stored in L165 are augmented by 1.4 Ang.
         
         do i=1,natom
            write(atype,'(a2)') ih(m06+i-1)
            call upper(atype)
            if (ix(i100) .eq. 1) then
               atomicnumber = ix(i100+i)
            else
               call get_atomic_number(ih(m04+i-1), x(lmass+i-1), &
                  atomicnumber, errFlag)
#ifdef LES
               if(errFlag) then
                  call get_atomic_number(ih(m04+i-1), &
                     x(lmass+i-1)*lesfac((lestyp(i)-1)*nlesty+lestyp(i)), &
                     atomicnumber, errFlag)
               end if
#endif
            end if
            hybridization = ix(i80-i)
            nbond=ix(i80+i-1)
            if (atomicnumber .eq. 6) then
               if (hybridization .eq. 4) then
                  if (nbond == 1) then
                     x(l165-1+i) = 1.70d0 + 1.4d0
                     x(l170-1+i) = 0.77887d0
                     x(l175-1+i) = -0.28063d0
                     x(l180-1+i) = -0.0012968d0
                     x(l185-1+i) = 0.00039328d0
                  else if (nbond == 2) then
                     x(l165-1+i) = 1.70d0 + 1.4d0
                     x(l170-1+i) = 0.56482d0
                     x(l175-1+i) = -0.19608d0
                     x(l180-1+i) = -0.0010219d0
                     x(l185-1+i) = 0.0002658d0
                  else if (nbond == 3) then
                     x(l165-1+i) = 1.70d0 + 1.4d0
                     x(l170-1+i) = 0.23348d0
                     x(l175-1+i) = -0.072627d0
                     x(l180-1+i) = -0.00020079d0
                     x(l185-1+i) = 0.00007967d0
                  else if (nbond == 4) then
                     x(l165-1+i) = 1.70d0 + 1.4d0
                     x(l170-1+i) = 0.00000d0
                     x(l175-1+i) = 0.00000d0
                     x(l180-1+i) = 0.00000d0
                     x(l185-1+i) = 0.00000d0
                  else
                     write(6,*) 'Unusual nbond for CT:', i, nbond, &
                        ' Using default carbon LCPO parameters'
                     x(l165-1+i) = 1.70d0 + 1.4d0
                     x(l170-1+i) = 0.77887d0
                     x(l175-1+i) = -0.28063d0
                     x(l180-1+i) = -0.0012968d0
                     x(l185-1+i) = 0.00039328d0
                  end if
               else
                  if (nbond == 2) then
                     x(l165-1+i) = 1.70d0 + 1.4d0
                     x(l170-1+i) = 0.51245d0
                     x(l175-1+i) = -0.15966d0
                     x(l180-1+i) = -0.00019781d0
                     x(l185-1+i) = 0.00016392d0
                  else if (nbond == 3) then
                     x(l165-1+i) = 1.70d0 + 1.4d0
                     x(l170-1+i) = 0.070344d0
                     x(l175-1+i) = -0.019015d0
                     x(l180-1+i) = -0.000022009d0
                     x(l185-1+i) = 0.000016875d0
                  else
                     write(6,*) 'Unusual nbond for C :', i, nbond, &
                       ' Using default carbon LCPO parameters'
                     x(l165-1+i) = 1.70d0 + 1.4d0
                     x(l170-1+i) = 0.77887d0
                     x(l175-1+i) = -0.28063d0
                     x(l180-1+i) = -0.0012968d0
                     x(l185-1+i) = 0.00039328d0
                  end if
               end if
            else if (atomicnumber .eq. 8) then
               if (atype == 'O ') then
                  x(l165-1+i) = 1.60d0 + 1.4d0
                  x(l170-1+i) = 0.68563d0
                  x(l175-1+i) = -0.1868d0
                  x(l180-1+i) = -0.00135573d0
                  x(l185-1+i) = 0.00023743d0
               else if (atype == 'O2') then
                  x(l165-1+i) = 1.60d0 + 1.4d0
                  x(l170-1+i) = 0.88857d0
                  x(l175-1+i) = -0.33421d0
                  x(l180-1+i) = -0.0018683d0
                  x(l185-1+i) = 0.00049372d0
               else
                  if (nbond == 1) then
                     x(l165-1+i) = 1.60d0 + 1.4d0
                     x(l170-1+i) = 0.77914d0
                     x(l175-1+i) = -0.25262d0
                     x(l180-1+i) = -0.0016056d0
                     x(l185-1+i) = 0.00035071d0
                  else if (nbond == 2) then
                     x(l165-1+i) = 1.60d0 + 1.4d0
                     x(l170-1+i) = 0.49392d0
                     x(l175-1+i) = -0.16038d0
                     x(l180-1+i) = -0.00015512d0
                     x(l185-1+i) = 0.00016453d0
                  else
                     write(6,*) 'Unusual nbond for O:', i, nbond, &
                        ' Using default oxygen LCPO parameters'
                     x(l165-1+i) = 1.60d0 + 1.4d0
                     x(l170-1+i) = 0.77914d0
                     x(l175-1+i) = -0.25262d0
                     x(l180-1+i) = -0.0016056d0
                     x(l185-1+i) = 0.00035071d0
                  end if
               end if
            else if(atomicnumber .eq. 7) then
               if (atype == 'N3') then
                  if (nbond == 1) then
                     x(l165-1+i) = 1.65d0 + 1.4d0
                     x(l170-1+i) = 0.078602d0
                     x(l175-1+i) = -0.29198d0
                     x(l180-1+i) = -0.0006537d0
                     x(l185-1+i) = 0.00036247d0
                  else if (nbond == 2) then
                     x(l165-1+i) = 1.65d0 + 1.4d0
                     x(l170-1+i) = 0.22599d0
                     x(l175-1+i) = -0.036648d0
                     x(l180-1+i) = -0.0012297d0
                     x(l185-1+i) = 0.000080038d0
                  else if (nbond == 3) then
                     x(l165-1+i) = 1.65d0 + 1.4d0
                     x(l170-1+i) = 0.051481d0
                     x(l175-1+i) = -0.012603d0
                     x(l180-1+i) = -0.00032006d0
                     x(l185-1+i) = 0.000024774d0
                  else
                     write(6,*) 'Unusual nbond for N3:', i, nbond, &
                        ' Using default nitrogen LCPO parameters'
                     x(l165-1+i) = 1.65d0 + 1.4d0
                     x(l170-1+i) = 0.078602d0
                     x(l175-1+i) = -0.29198d0
                     x(l180-1+i) = -0.0006537d0
                     x(l185-1+i) = 0.00036247d0
                  end if
               else
                  if (nbond == 1) then
                     x(l165-1+i) = 1.65d0 + 1.4d0
                     x(l170-1+i) = 0.73511d0
                     x(l175-1+i) = -0.22116d0
                     x(l180-1+i) = -0.00089148d0
                     x(l185-1+i) = 0.0002523d0
                  else if (nbond == 2) then
                     x(l165-1+i) = 1.65d0 + 1.4d0
                     x(l170-1+i) = 0.41102d0
                     x(l175-1+i) = -0.12254d0
                     x(l180-1+i) = -0.000075448d0
                     x(l185-1+i) = 0.00011804d0
                  else if (nbond == 3) then
                     x(l165-1+i) = 1.65d0 + 1.4d0
                     x(l170-1+i) = 0.062577d0
                     x(l175-1+i) = -0.017874d0
                     x(l180-1+i) = -0.00008312d0
                     x(l185-1+i) = 0.000019849d0
                  else
                     write(6,*) 'Unusual nbond for N:', i, nbond, &
                        ' Using default nitrogen LCPO parameters'
                     x(l165-1+i) = 1.65d0 + 1.4d0
                     x(l170-1+i) = 0.078602d0
                     x(l175-1+i) = -0.29198d0
                     x(l180-1+i) = -0.0006537d0
                     x(l185-1+i) = 0.00036247d0
                  end if
               end if
            else if(atomicnumber .eq. 16) then
               if (atype == 'SH') then
                  x(l165-1+i) = 1.90d0 + 1.4d0
                  x(l170-1+i) = 0.7722d0
                  x(l175-1+i) = -0.26393d0
                  x(l180-1+i) = 0.0010629d0
                  x(l185-1+i) = 0.0002179d0
               else
                  x(l165-1+i) = 1.90d0 + 1.4d0
                  x(l170-1+i) = 0.54581d0
                  x(l175-1+i) = -0.19477d0
                  x(l180-1+i) = -0.0012873d0
                  x(l185-1+i) = 0.00029247d0
               end if
            else if (atomicnumber .eq. 15) then
               if (nbond == 3) then
                  x(l165-1+i) = 1.90d0 + 1.4d0
                  x(l170-1+i) = 0.3865d0
                  x(l175-1+i) = -0.18249d0
                  x(l180-1+i) = -0.0036598d0
                  x(l185-1+i) = 0.0004264d0
               else if (nbond == 4) then
                  x(l165-1+i) = 1.90d0 + 1.4d0
                  x(l170-1+i) = 0.03873d0
                  x(l175-1+i) = -0.0089339d0
                  x(l180-1+i) = 0.0000083582d0
                  x(l185-1+i) = 0.0000030381d0
               else
                  write(6,*) 'Unusual nbond for P:', i, nbond, &
                     ' Using default phosphorus LCPO parameters'
                  x(l165-1+i) = 1.90d0 + 1.4d0
                  x(l170-1+i) = 0.3865d0
                  x(l175-1+i) = -0.18249d0
                  x(l180-1+i) = -0.0036598d0
                  x(l185-1+i) = 0.0004264d0
               end if
            else if (atype(1:1) == 'Z') then
               x(l165-1+i) = 0.00000d0 + 1.4d0
               x(l170-1+i) = 0.00000d0
               x(l175-1+i) = 0.00000d0
               x(l180-1+i) = 0.00000d0
               x(l185-1+i) = 0.00000d0
            else if (atomicnumber .eq. 1) then
               x(l165-1+i) = 0.00000d0 + 1.4d0
               x(l170-1+i) = 0.00000d0
               x(l175-1+i) = 0.00000d0
               x(l180-1+i) = 0.00000d0
               x(l185-1+i) = 0.00000d0
            else if (atype == 'MG') then
               !  Mg radius = 0.99A: ref. 21 in J. Chem. Phys. 1997, 107, 5422
               !  Mg radius = 1.18A: ref. 30 in J. Chem. Phys. 1997, 107, 5422
               !  Mg radius = 1.45A: Aqvist 1992
               x(l165-1+i) = 1.18d0 + 1.4d0
               !  The following values were taken from O.sp3 with two bonded 
               !  neighbors -> O has the smallest van der Waals radius 
               ! compared to all other elements which had been parametrized
               x(l170-1+i) = 0.49392d0
               x(l175-1+i) = -0.16038d0
               x(l180-1+i) = -0.00015512d0
               x(l185-1+i) = 0.00016453d0
            else if (atype == 'F') then
               x(l165-1+i) = 1.47d0 + 1.4d0
               x(l170-1+i) = 0.68563d0
               x(l175-1+i) = -0.1868d0
               x(l180-1+i) = -0.00135573d0
               x(l185-1+i) = 0.00023743d0 
            else
               ! write( 0,* ) 'bad atom type: ',atype
               ! call mexit( 6,1 )
               x(l165-1+i) = 1.70 + 1.4;
               x(l170-1+i) = 0.51245;
               x(l175-1+i) = -0.15966;
               x(l180-1+i) = -0.00019781;
               x(l185-1+i) = 0.00016392;
               write(6,'(a,a)') 'Using carbon SA parms for atom type', atype 
            end if
         end do  !  i=1,natom
         !
      else if ( gbsa == 2 ) then

         !     --- assign parameters for calculating SASA according to the
         !         ICOSA method; the radii are augmented by 1.4 A ---

         do i=1,natom
            write(atype,'(a2)') ih(m06+i-1)
            if(ix(i100) .eq. 1) then
               atomicnumber = ix(i100+i)
            else
               call get_atomic_number(ih(m04+i-1), x(lmass+i-1), &
                  atomicnumber, errFlag)
#ifdef LES
               if(errFlag) then
                  call get_atomic_number(ih(m04+i-1), &
                     x(lmass+i-1)*lesfac((lestyp(i)-1)*nlesty+lestyp(i)), &
                     atomicnumber, errFlag)
               end if
#endif
            end if
            if (atomicnumber .eq. 7) then
               x(L165-1+i) = 1.55d0 + 1.4d0
            else if (atomicnumber .eq. 6) then
               x(L165-1+i) = 1.70d0 + 1.4d0
            else if (atomicnumber .eq. 1) then
               x(L165-1+i) = 1.20d0 + 1.4d0
            else if (atomicnumber .eq. 8) then
               x(L165-1+i) = 1.50d0 + 1.4d0
            else if (atomicnumber .eq. 15) then
               x(L165-1+i) = 1.80d0 + 1.4d0
            else if (atomicnumber .eq. 16) then
               x(L165-1+i) = 1.80d0 + 1.4d0
            else if (atomicnumber .eq. 17) then
               ! Cl radius
               x(L165-1+i) = 1.70d0 + 1.4d0
            else if (atomicnumber .eq. 9) then
               ! F radius
               x(L165-1+i) = 1.50d0 + 1.4d0
            else if (atomicnumber .eq. 35) then
               ! Br radius
               ! Bondi, J. Phys. Chem. 1964, 68, 441.
               x(L165-1+i) = 1.85d0 + 1.4d0
            else if (atomicnumber .eq. 20) then
               ! Ca radius
               ! Calculated from Aqvist, J. Phys. Chem. 1990, 94, 8021.
               x(L165-1+i) = 1.33d0 + 1.4d0
            else if (atomicnumber .eq. 11) then
               ! Na radius
               ! Calculated from Aqvist, J. Phys. Chem. 1990, 94, 8021.
               x(L165-1+i) = 1.87d0 + 1.4d0
            else if (atomicnumber .eq. 30) then
               ! Zn radius
               ! Hoops, Anderson, Merz, JACS 1991, 113, 8262.
               x(L165-1+i) = 1.10d0 + 1.4d0
            else if (atomicnumber .eq. 12) then
               !  Mg radius = 0.99A: ref. 21 in J. Chem. Phys. 1997, 107, 5422
               !  Mg radius = 1.18A: ref. 30 in J. Chem. Phys. 1997, 107, 5422
               !  Mg radius = 1.45A: Aqvist 1992
               x(L165-1+i) = 1.18d0 + 1.4d0
            else
               write( 0,* ) 'bad atom type: ',atype
               FATAL_ERROR
            end if

            ! dummy LCPO values:
            x(L170-1+i) = 0.0d0
            x(L175-1+i) = 0.0d0
            x(L180-1+i) = 0.0d0
            x(L185-1+i) = 0.0d0
            !  write(6,*) i,' ',atype,x(L165-1+i)
         end do  !  i=1,natom

      end if ! ( gbsa == 1 )
      
   end if  ! ( igb /= 0 .and. igb /= 10 .and. ipb == 0 )

   !-----------------------------------
   ! If a LRT calculation is requested, 
   ! setup the icosa-SASA parameters 
   ! this code is copied from above
   !-----------------------------------

   if ( ilrt /= 0 ) then
      !     --- assign parameters for calculating SASA according to the
      !         ICOSA method; the radii are augmented by 1.4 A ---
      do i=1,natom
         write(atype,'(a2)') ih(m06+i-1)
         if(ix(i100) .eq. 1) then
            atomicnumber = ix(i100+i)
         else
            call get_atomic_number(ih(m04+i-1), x(lmass+i-1), &
               atomicnumber, errFlag)
#ifdef LES
            if(errFlag) then
               call get_atomic_number(ih(m04+i-1), &
                  x(lmass+i-1)*lesfac((lestyp(i)-1)*nlesty+lestyp(i)), &
                  atomicnumber, errFlag)
            end if
#endif
         end if
         if (atomicnumber .eq. 7) then
            x(L165-1+i) = 1.55d0 + 1.4d0
         else if (atomicnumber .eq. 6) then
            x(L165-1+i) = 1.70d0 + 1.4d0
         else if (atomicnumber .eq. 1 .or. &
              ! added for lone pairs
              atype == 'EP') then
            x(L165-1+i) = 1.20d0 + 1.4d0
         else if (atomicnumber .eq. 8) then
            x(L165-1+i) = 1.50d0 + 1.4d0
         else if (atomicnumber .eq. 15) then
            x(L165-1+i) = 1.80d0 + 1.4d0
         else if (atomicnumber .eq. 16) then
            x(L165-1+i) = 1.80d0 + 1.4d0
         else if (atomicnumber .eq. 12) then
            !             Mg radius = 0.99A: ref. 21 in J. Chem. Phys. 1997, 107, 5422
            !             Mg radius = 1.18A: ref. 30 in J. Chem. Phys. 1997, 107, 5422
            !             Mg radius = 1.45A: Aqvist 1992
            x(L165-1+i) = 1.18d0 + 1.4d0
         else
            write( 0,* ) 'bad atom type: ',atype,' cannot perform SASA calculation'
            FATAL_ERROR
         end if  !  atype(1:1) == 'N'
         x(L170-1+i) = 0.0d0
         x(L175-1+i) = 0.0d0
         x(L180-1+i) = 0.0d0
         x(L185-1+i) = 0.0d0
         !write(6,*) i,' ',atype,x(L165-1+i)
      end do  !  i=1,natom
      
   end if ! ( ilrt /= 0 )
   
   !------------------------------------------------------------------------
   ! If user has requested Poisson-Boltzmann electrostatics, set up variables
   !------------------------------------------------------------------------
  
   if ( igb == 10 .or. ipb /= 0 ) then
      call pb_init(ifcap,natom,nres,ntypes,nbonh,nbona,ix(i02),ix(i04),ix(i06),ix(i08),ix(i10),&
                   ix(iibh),ix(ijbh),ix(iiba),ix(ijba),ix(ibellygp),ih(m02),ih(m04),ih(m06),x(l15),x(l97))
   end if  ! ( igb == 10 .or. ipb /= 0 ) 

   if (icnstph /= 0) then
      !  Initialize all constant pH data to 0 and read it in
      call cnstph_zero()
      call cnstphread(x(l15))

      !     Fill proposed charges array from current charges
      do i=1,natom
         x(l190-1+i) = x(l15-1+i)
      end do

      !  If we're doing explicit CpH, fill gbv* arrays
      if ( icnstph .gt. 1 .and. cph_igb == 2 .or. cph_igb == 5) then
        do i=1,natom
            x(l2402+i-1) = gbalpha
            x(l2403+i-1) = gbbeta
            x(l2404+i-1) = gbgamma
        end do
      end if
   end if

!  +---------------------------------------------------------------+
!  |  Read EVB input file                                          |
!  +---------------------------------------------------------------+

   if( ievb /= 0 ) then
#ifdef MPI
      call evb_input
      call evb_init
#  if defined(LES)
!KFW  call evb_pimd_init
#  endif
#else
      write(6,'(/2x,a)') 'Setting ievb>0 requires compilation with MPI'
      FATAL_ERROR
#endif
   endif

   if( iyammp /= 0 ) write( 6, '(a)' ) '  Using yammp non-bonded potential'
   
   ! -------------------------------------------------------------------
   !
   ! -- add check to see if the space in nmr.h is likely to be
   !     too small for this run:
   !     [Note: this check does *not* indicate if MXTAU and MXP are
   !      too small.  There is no easy way to ensure this, since
   !      the experimental intensities are read in a namelist
   !      command: if too many intensities are input, the read
   !      statment may cause a coredump before returning control
   !      to the main program.  Be careful.  sigh....]
   
   if (natom > matom .and. nmropt > 1) then
      write(6,*) 'WARNING: MATOM in nmr.h is smaller than the ', &
            natom,' atoms in this molecule.'
      write(6,*) 'Printout of NMR violations may be compromised.'
   end if
   
   ! -------------------------------------------------------------------
   !     --- checks on bogus data ---
   ! -------------------------------------------------------------------
   
   inerr = 0
      
   if( icfe < 0 .or. icfe > 1 ) then
      write(6,*) 'icfe must be 0 or 1 (icfe=2 is no longer supported)'
      DELAYED_ERROR
   end if
   if( icfe /= 0 .and. numgroup /= 2 ) then
      write(6,*) 'numgroup must be 2 if icfe is set'
      DELAYED_ERROR
   end if
   if (ievb>0) then
#ifdef MPI
!KFW  if( numgroup /= 2 ) then
!KFW     write(6,*) 'numgroup must be 2 if ievb is set'
!KFW     DELAYED_ERROR
!KFW  end if
#else
      write(6,'(/2x,a)') 'Setting ievb>0 requires compilation with MPI'
      DELAYED_ERROR
#endif
   end if
#ifdef PUPIL_SUPPORT
   ! BPR: PUPIL does not work with GB (or, I suppose, PB) for the
   ! time being. It is known to either crash or produce bogus
   ! results.
   if (igb > 0 .or. ipb /= 0) then
      write(6,'(a)') 'Cannot use implicit solvation (GB or PB) with PUPIL'
      DELAYED_ERROR
   end if
#endif /*PUPIL_SUPPORT*/
   if( (igb > 0 .or. ipb /= 0) .and. numextra > 0) then
      write(6,'(a)') 'Cannot use igb>0 with extra-point force fields'
      DELAYED_ERROR
   end if

!AMD validation
   if(iamd.gt.0)then
     if (EthreshD .eq. 0.d0 .and. alphaD .eq. 0.d0 .and. EthreshP .eq. 0.d0 .and. alphaP .eq. 0.d0) then
     write(6,'(a,i3)')'| AMD error all main parameters are 0.0 for Accelerated MD (AMD) or Windowed Accelerated MD (wAMD) '
     DELAYED_ERROR
     endif
    if(w_amd.gt.0)then
     if (EthreshD_w .eq. 0.d0 .and. alphaD_w .eq. 0.d0 .and. EthreshP_w .eq. 0.d0 .and. alphaP_w .eq. 0.d0) then
     write(6,'(a,i3)')'| AMD error all extra parameters are 0.0 for Windowed Accelerated MD (wAMD) LOWERING BARRIERS'
     DELAYED_ERROR
     endif
    endif
   endif
   if (iamd .gt. 0 .and. iamoeba ==1) then
      write(6,*)'amoeba is incompatible with AMD for now'
      inerr=1
   end if

   if (ips < 0 .or. ips > 6) then
      write(6,'(/2x,a,i3,a)') 'IPS (',ips,') must be between 0 and 6'
      DELAYED_ERROR
   end if
   if (ips /= 0 .and. ipol > 0 ) then
      write(6,'(/2x,a)') 'IPS and IPOL are inconsistent options'
      DELAYED_ERROR
   endif
   if (ips /= 0 .and. lj1264 /= 0) then
      write(6, '(/2x,a)') 'IPS and the LJ 12-6-4 potential are incompatible'
      DELAYED_ERROR
   endif
   if ( (igb > 0 .or. ipb /= 0) .and. ips > 0 ) then
      write(6,'(/2x,a,i3,a,i3,a)') 'IGB (',igb,') and ips (',ips, &
          ') cannot both be turned on'
      DELAYED_ERROR
   end if
   if (igb /= 0 .and. igb /= 1 .and. igb /= 2 .and. igb /= 5 &
         .and. igb /= 6 .and. igb /= 7 .and. igb /= 8 .and. igb /= 10) then
      write(6,'(/2x,a,i3,a)') 'IGB (',igb,') must be 0,1,2,5,6,7,8 or 10.'
      DELAYED_ERROR
   end if
   if (alpb /= 0 .and. alpb /= 1 )  then
      write(6,'(/2x,a,i3,a)') 'ALPB (',alpb,') must be 0 or 1.'
      DELAYED_ERROR
   end if
   if (alpb /= 0 .and. igb /= 1 .and. igb /= 2 .and. igb /= 5 .and. igb /=7 )  then
      write(6,'(/2x,a,i3,a)') 'IGB (',igb,') must be 1,2,5, or 7 if ALPB > 0.'
      DELAYED_ERROR
   end if

   if (jar < 0 .or. jar > 1) then
      write(6,'(/2x,a,i3,a)') 'JAR (',jar,') must be 0 or 1'
      DELAYED_ERROR
   end if

#ifdef LES
   if( igb /= 0 .and. igb /= 1 .and. igb /= 5 .and. igb /=7 ) then
      write(6,'(/,a)') 'Error: LES is only compatible with IGB > 0,1,5,7'
      DELAYED_ERROR
   end if
   if( alpb /= 0) then
      write(6,'(/,a)') 'Error: LES is not compatible with ALPB'
      DELAYED_ERROR
   end if
   if( gbsa > 0 ) then
      write(6,'(/,a)') 'Error: LES is not compatible with GBSA > 0'
      DELAYED_ERROR
   end if
   if( qmmm_nml%ifqnt ) then
      write(6,'(/,a)') 'Error: LES is not compatible with QM/MM'
      DELAYED_ERROR
   end if
   if( ipol > 0 ) then
      write(6,'(/,a)') 'Error: LES is not compatible with IPOL > 0'
      DELAYED_ERROR
   end if
   if (temp0les >= 0.d0 .and. iscale > 0 ) then
      write (6,'(/,a)') 'Error: iscale cannot be used with temp0les'
      DELAYED_ERROR
   end if
#else
   rdt = 0
#endif
   if (irest /= 0 .and. irest /= 1) then
      write(6,'(/2x,a,i3,a)') 'IREST (',irest,') must be 0 or 1.'
      DELAYED_ERROR
   end if
   if (ibelly /= 0 .and. ibelly /= 1) then
      write(6,'(/2x,a,i3,a)') 'IBELLY (',ibelly,') must be 0 or 1.'
      DELAYED_ERROR
   end if
   if (imin < 0) then
      write(6,'(/2x,a,i3,a)') 'IMIN (',imin,') must be >= 0.'
      DELAYED_ERROR
   end if
   if (imin == 5) then
      if (ifbox /= 0 .and. ntb == 2) then
         write(6,'(/2x,a)') 'WARNING: IMIN=5 with changing periodic boundaries (NTB=2) can result in'
         write(6,'(/2x,a)') '         odd energies being calculated. Use with caution.'
      endif
#ifdef MPI
      if (sandersize > 1 .and. ntb == 2) then
         write(6,'(/2x,a)') 'ERROR: IMIN=5 and NTB=2 cannot be run with multiple processors.'
         DELAYED_ERROR
      endif
#endif
   end if


   if (iscale > mxvar) then
      write(6,9501) iscale,mxvar
      9501 format('ERROR: ISCALE (',i5,') exceeds MXVAR (',i5, &
            '). See nmr.h')
      DELAYED_ERROR
   end if
   if (ntx < 1 .or. ntx > 7) then
      write(6,'(/2x,a,i3,a)') 'NTX (',ntx,') must be in 1..7'
      DELAYED_ERROR
   end if
   
   if (ntb /= 0 .and. ntb /= 1 .and. ntb /= 2) then
      write(6,'(/2x,a,i3,a)') 'NTB (',ntb,') must be 0, 1 or 2.'
      DELAYED_ERROR
   end if
   if (ntb == 0 .and. iwrap > 0) then
      write(6,'(/2x,a)') 'Error: IWRAP > 0 cannot be used without a periodic box.'
      DELAYED_ERROR
   end if
   
   if (ntt < 0 .or. ntt > 9) then                                      ! APJ
      write(6,'(/2x,a,i3,a)') 'NTT (',ntt,') must be between 0 and 9.' ! APJ
      DELAYED_ERROR
   end if
   if (ntt == 1 .and. tautp < dt) then
      write(6, '(/2x,a,f6.2,a)') 'TAUTP (',tautp,') < DT (step size)'
      DELAYED_ERROR
   end if
   if( ntt < 3 .or. ntt > 9 ) then  ! APJ
      if( gamma_ln > 0.d0 ) then
         write(6,'(a)') 'ntt must be 3 to 9 if gamma_ln > 0' ! APJ
         DELAYED_ERROR
      end if
   end if

   if ( ntt==3 .or. ntt==6 ) nchain = 0 !APJ: Langevin, Adaptive-Langevin must have chain set to zero. ! APJ

   if (ntt == 3 .or. ntt == 4) then
      if ( ntb == 2) then
        !Require gamma_ln > 0.0d0 for ntt=3 and ntb=2 - strange things happen
        !if you run NPT with NTT=3 and gamma_ln = 0.
        if ( gamma_ln <= 0.0d0 ) then
          write(6,'(a)') 'gamma_ln must be > 0 for ntt=3 .or. 4 with ntb=2.'
          DELAYED_ERROR
        end if
      end if
   end if

   ! Isokinetic ensemble checks

   if (ntt == 9) then
      if(gamma_ln <= 0.d0) then
         write(6,'(a)') 'gamma_ln must be > 0 when ntt = 9'
         DELAYED_ERROR
      end if
      if(nkija < 1) then
         write(6,'(a)') 'nkija must be >= 1 when ntt = 9'
         DELAYED_ERROR
      end if
      if(idistr < 0) then
         write(6,'(a)') 'idistr must be >= 0 when ntt = 9'
         DELAYED_ERROR
      end if
      if(ntc /= 1 .or. ntf /= 1) then
         write(6,'(a)') 'ntc and ntf must be = 1 when ntt = 9'
         DELAYED_ERROR
      end if
      if(tempi > temp0) then
         write(6,'(a)') 'tempi must be <= temp0 when ntt = 9'
         DELAYED_ERROR
      end if
   end if

   if (ntp /= 0 .and. ntp /= 1 .and. ntp /= 2 .and. ntp /= 3) then
      write(6,'(/2x,a,i3,a)') 'NTP (',ntp,') must be 0, 1, 2, or 3.'
      DELAYED_ERROR
   end if
   if (ntp == 3 .and. csurften < 1) then
      write(6,'(/2x,a)') 'csurften must be greater than 0 for ntp=3.'
      DELAYED_ERROR
   end if
   if (ntp > 0 .and. taup < dt .and. barostat == 1) then
      write(6, '(/2x,a,f6.2,a)') 'TAUP (',taup,') < DT (step size)'
      DELAYED_ERROR
   end if
   if (npscal < 0 .or. npscal > 1) then
      write(6,'(/2x,a,i3,a)') 'NPSCAL (',npscal,') must be 0 or 1.'
      DELAYED_ERROR
   end if
   if (ntp > 0) then
      if (barostat /= 1 .and. barostat /= 2) then
         write(6, '(/2x,a,i3,a)') 'BAROSTAT (', barostat, ') must be 1 or 2'
         DELAYED_ERROR
      end if
      if (barostat == 2) then
         if (mcbarint <= 0) then
            write(6, '(/2x,a,i3,a)') 'MCBARINT (',mcbarint,') must be positive'
            DELAYED_ERROR
         end if
         if (mcbarint >= nstlim) then
            write(6, '(a)') 'WARNING: mcbarint is greater than the number of &
                     &steps. This is effectively constant volume.'
         end if
      end if
   end if
   
   if (ntc < 1 .or. ntc > 4) then
      write(6,'(/2x,a,i3,a)') 'NTC (',ntc,') must be 1,2,3 or 4.'
      DELAYED_ERROR
   end if
   if (jfastw < 0 .or. jfastw > 4) then
      write(6,'(/2x,a,i3,a)') 'JFASTW (',jfastw,') must be 0->4.'
      DELAYED_ERROR
   end if
   
   if (ntf < 1 .or. ntf > 8) then
      write(6,'(/2x,a,i3,a)') 'NTF (',ntf,') must be in 1..8.'
      DELAYED_ERROR
   end if
   
   if (ioutfm /= 0 .and. ioutfm /= 1) then
      write(6,'(/2x,a,i3,a)') 'IOUTFM (',ioutfm,') must be 0 or 1.'
      DELAYED_ERROR
   end if
   
   if (ntpr < 0) then
      write(6,'(/2x,a,i3,a)') 'NTPR (',ntpr,') must be >= 0.'
      DELAYED_ERROR
   end if
   if (ntwx < 0) then
      write(6,'(/2x,a,i3,a)') 'NTWX (',ntwx,') must be >= 0.'
      DELAYED_ERROR
   end if
   if (ntwv < -1) then
      write(6,'(/2x,a,i3,a)') 'NTWV (',ntwv,') must be >= -1.'
      DELAYED_ERROR
   end if
   if (ntwf < -1) then
      write(6, '(/2x,a,i3,a)') 'NTWF (',ntwf,') must be >= -1.'
      DELAYED_ERROR
   end if
   if (ntwv == -1 .and. ioutfm /= 1) then
      write (6, '(/2x,a)') 'IOUTFM must be 1 for NTWV == -1.'
      DELAYED_ERROR
   end if
   if (ntwf == -1 .and. ioutfm /= 1) then
      write(6, '(/2x,a)') 'IOUTFM must be 1 for NTWF == -1.'
      DELAYED_ERROR
   end if
   if (ntwv == -1 .and. ntwx == 0) then
      write (6, '(/2x,a)') 'NTWX must be > 0 for NTWV == -1.'
      DELAYED_ERROR
   end if
   if (ntwe < 0) then
      write(6,'(/2x,a,i3,a)') 'NTWE (',ntwe,') must be >= 0.'
      DELAYED_ERROR
   end if
   if (ntave < 0) then
      write(6,'(/2x,a,i3,a)') 'NTAVE (',ntave,') must be >= 0.'
      DELAYED_ERROR
   end if
   if (ntr /= 0 .and. ntr /= 1) then
      write(6,'(/2x,a,i3,a)') 'NTR (',ntr,') must be 0 or 1.'
      DELAYED_ERROR
   end if
   if (ntrx /= 0 .and. ntrx /= 1) then
      write(6,'(/2x,a,i3,a)') 'NTRX (',ntrx,') must be 0 or 1.'
      DELAYED_ERROR
   end if
   if (nmropt < 0 .or. nmropt > 2) then
      write(6,'(/2x,a,i3,a)') 'NMROPT (',nmropt,') must be in 0..2.'
      DELAYED_ERROR
   end if
   
   if (idecomp < 0 .or. idecomp > 4) then
      write(6,'(/2x,a)') 'IDECOMP must be 0..4'
      DELAYED_ERROR
   end if
      
   ! check settings related to ivcap

   if(ivcap == 3 .or. ivcap == 4) then
      write(6,'(/2x,a)') 'IVCAP == 3 and IVCAP == 4 currently not implemented'
      DELAYED_ERROR
   endif
   if (ivcap < 0 .and. ivcap > 5) then
      write(6,'(/2x,a)') 'IVCAP must be 0 ... 5'
      DELAYED_ERROR
   end if
   if ((ivcap == 1 .or. ivcap == 5) .and. igb /= 10 .and. ipb == 0) then
      write(6,'(/2x,a)') 'IVCAP == 1,5 only works with Poisson Boltzmann (igb=10 or ipb/=0)'
      DELAYED_ERROR
   end if
   if((ivcap == 1 .or. ivcap == 3 .or. ivcap == 5 ) .and. cutcap <= 0.0d0) then
      write(6,'(/2x,a)') 'For IVCAP == 1,3, or 5, cutcap must be > 0.0'
      DELAYED_ERROR
   endif
   if (ivcap == 4 .and. &
         (xlorth < ZERO .or. ylorth < ZERO .or. zlorth < ZERO .or. &
!  give magic numbers a name  srb aug 2007 !
          xorth > 47114710.0d0 .or. &
!  give magic numbers a name !
          yorth > 47114710.0d0 .or. &
!  give magic numbers a name !
          zorth > 47114710.0d0)) then
      write(6,'(/2x,a)') &
      'For IVCAP == 4, xlorth, ylorth, zlorth, xorth, yorth, zorth must be set'
      DELAYED_ERROR
   end if
   if ((ivcap == 3 .or. ivcap == 4) .and. ibelly == 0) then
      write(6,'(/2x,a,a)') &
         'For IVCAP == 3 or 4, ibelly must be 1 and all atoms', &
         '  not in the spherical or orthorhombic region must be set NOT moving'
      DELAYED_ERROR
   end if
   if (ivcap == 5 .and. (imin /= 1 .or. maxcyc > 1)) then
      write(6,'(/2x,a,a)') &
         'IVCAP == 5 only works for single-point energy calculation'
      DELAYED_ERROR
   end if

   ! check if ifbox variable from prmtop file matches actual angles:

   if ( igb == 0 .and. ipb == 0 .and. ntb /= 0 ) then
      if ( ifbox == 1 ) then
         if ( abs(alpha - 90.0d0) > 1.d-5 .or. &
           abs(beta  - 90.0d0) > 1.d-5 .or. &
           abs(gamma - 90.0d0) > 1.d-5 ) then
           ifbox =3
#ifndef API
           write(6,'(a)') '     Setting ifbox to 3 for non-orthogonal unit cell'
#endif
         end if
      end if

      if ( ifbox == 2 ) then
         if ( abs(alpha - 109.4712190d0) > 1.d-5 .or. &
              abs(beta  - 109.4712190d0) > 1.d-5 .or. &
              abs(gamma - 109.4712190d0) > 1.d-5 ) then
              write(6,'(/2x,a)') &
              'Error: ifbox=2 in prmtop but angles are not correct'
              DELAYED_ERROR
         end if
      end if
   end if

   ! checks for targeted MD
   if (itgtmd /= 0 .and. itgtmd > 2) then
      write(6,'(/2x,a,i3,a)') 'ITGTMD (',itgtmd,') must be 0, 1 or 2.'
      DELAYED_ERROR
   end if
   if (itgtmd == 1 .and. ntr == 1) then
      if (len_trim(tgtfitmask) > 0 .or. len_trim(tgtrmsmask) <= 0) then
         write(6,'(/2x,a)') 'ITGTMD: tgtrmsmask (and not tgtfitmask) ' //  &
                            'should be specified if NTR=1'
         DELAYED_ERROR
      end if
   end if
   ! skip this test until fallback to rgroup() is supported
   !if (itgtmd == 1 .and. ntr == 0) then
   !   if (len_trim(tgtfitmask) == 0 .and. len_trim(tgtrmsmask) == 0) then
   !      write(6,'(/2x,a)')  &
   !        'ITGTMD: both tgtfitmask and tgtrmsmask should be specified if NTR=0'
   !      DELAYED_ERROR
   !   end if
   !end if
   
   !     -- consistency checking
   
   if (imin > 0.and.nrespa > 1)  then
      write(6,'(/2x,a)') 'For minimization, set nrespa,nrespai=1'
      DELAYED_ERROR
   end if
   if (ntp > 0 .and. nrespa > 1) then
      write(6,'(/2x,a)') 'nrespa must be 1 if ntp>0'
      DELAYED_ERROR
   end if
   if  (ntx < 4.and.init /= 3)  then
      write(6,'(/2x,a)') 'NTX / IREST inconsistency'
      DELAYED_ERROR
   end if
   if (ntb == 2 .and. ntp == 0) then
      write(6,'(/2x,a)') 'NTB set but no NTP option (must be 1 or 2)'
      DELAYED_ERROR
   end if
   if (ntp /= 0 .and. ntb /= 2) then
      write(6,'(/,a,a)')' NTP > 0 but not constant pressure P.B.C.', &
            ' (NTB = 2) must be used'
      DELAYED_ERROR
   end if
   if (ntb /= 0 .and. ifbox == 0 .and. ntp /= 0) then
      write(6,'(/,a)') ' (NTB /= 0 and NTP /= 0) but IFBOX == 0'
      write(6,'(/,a)') ' This combination is not supported'
      DELAYED_ERROR
   end if
   if (ntb /= 0 .and. &
         ( box(1) < 1.d0  .or. &
         box(2) < 1.d0  .or. &
         box(3) < 1.d0 ) ) then
      write(6,'(/,a,3f10.3)') ' BOX is too small: ',box(1),box(2),box(3)
      DELAYED_ERROR
   else if (ntb /= 0 .and. &
         (sqrt(cut) >= box(1)*0.5d0 .or. &
         sqrt(cut) >= box(2)*0.5d0 .or. &
         sqrt(cut) >= box(3)*0.5d0) ) then
      write(6,'(/,a)') ' CUT must be < half smallest box dimension'
      DELAYED_ERROR
   end if
   if (ntb /= 0 .and. (igb > 0 .or. ipb /= 0) ) then
      write(6,'(/,a)') ' igb>0 is only compatible with ntb=0'
      DELAYED_ERROR
   end if
#ifdef APBS
   if ( ntb == 0 .and. sqrt(cut) < 8.05 .and. igb /= 10 .and. ipb == 0 .and. &
      .not. mdin_apbs) then
#else
   if ( ntb == 0 .and. sqrt(cut) < 8.05 .and. igb /= 10 .and. ipb == 0 ) then
#endif /* APBS */
      write(6,'(/,a,f8.2)') ' unreasonably small cut for non-periodic run: ', &
         sqrt(cut)
      DELAYED_ERROR
   end if
   if ( rgbmax < 5.d0*fsmax ) then
      write(6,'(/,a,f8.2)') ' rgbmax must be at least ', 5.d0*fsmax
      DELAYED_ERROR
   end if
   if (icfe /= 0 .and. indmeth == 3 ) then
      write(6,'(/,a)') ' indmeth=3 cannot be used with icfe>0'
      DELAYED_ERROR
   end if
   if (icfe /= 0 .and. ibelly /= 0 ) then
      write(6,'(/,a)') ' ibelly cannot be used with icfe'
      DELAYED_ERROR
   end if
#ifdef MPI /* SOFT CORE */
   if (icfe /= 0 .and. dvdl_norest /= 0 ) then
      write(6,'(/,a)') 'dvdl_norest must == 0!'
      write(6,'(/,a)') 'The dvdl_norest option is deprecated.'
      write(6,'(/,a)') 'Restraint energies are seperated, &
        &and do not contribute to dvdl.'
      DELAYED_ERROR
   end if
#endif
   ! Modification done by Ilyas Yildirim
   if (icfe == 1 .and. (klambda < 1 .or. klambda > 6)) then
     write(6,'(/,a)') ' klambda must be between 1 and 6'
     DELAYED_ERROR
   end if
   ! End of modification done by Ilyas Yildirim                                       

   if (clambda < 0.d0 .or. clambda > 1.d0 ) then
      write(6,'(/,a)') ' clambda must be between 0 and 1'
      DELAYED_ERROR
   end if

   if (icfe /= 0 .and. (idecomp == 3 .or. idecomp == 4)) then
      write(6,'(/,a)') ' Pairwise decomposition for thermodynamic integration not implemented'
      DELAYED_ERROR
   end if
   if (icfe /= 0 .and. idecomp /= 0 .and. ipol > 0) then
      write(6,'(/,a)') ' IPOL is incompatible with IDECOMP and ICFE'
      DELAYED_ERROR
   end if
 
#ifdef MPI /* SOFT CORE */
   if (ifsc /= 0) then
      if (ifsc == 2) then
         write (6,'(/,a)') 'The ifsc=2 option is no longer supported. &
               &Internal energies of the soft core region.'
         write (6,'(a)') 'are now handled implicitly and setting ifsc=2 &
               &is no longer needed'
         DELAYED_ERROR
      end if
      if (icfe /= 1 .and. ifsc==1) then
         write (6,'(/,a)') ' Softcore potential requires a standard TI run, set icfe to 1'
         DELAYED_ERROR
      end if
      if ( igb > 0 .or. ipb /= 0 ) then
         write (6,'(/,a)') ' Softcore potential is incompatible with GB (for now)'
         DELAYED_ERROR
      end if
      if ( ntf > 1 ) then
         write (6,'(/,a)') ' Softcore potentials require ntf=1 because SHAKE &
               &constraints on some bonds might be removed'
         DELAYED_ERROR
      end if
      if (clambda > 0.995 .or. clambda < 0.005) then
         write (6,'(/,a)') ' Softcore potentials cannot be used with clambda < 0.005 or > 0.995'
         DELAYED_ERROR
      end if
      if (klambda /= 1) then
         write (6,'(/,a)') ' Softcore potential requires linear mixing, set klambda to 1'
         DELAYED_ERROR
      end if
      if (imin == 1 .and. ntmin /= 2) then
         write (6,'(/,a)') ' Minimizations with ifsc=1 require the steepest descent algorithm.'
         write (6,'(/,a)') ' Set ntmin to 2 and restart'
         DELAYED_ERROR
      end if
   end if
   ! Following is temporary (dec 2014), until/unless problem gets fixed
   if ( logdvdl /= 0 .and. ifsc == 0 ) then
      write (6,'(/,a)') ' Cannot request logdvdl unless ifsc>0'
      DELAYED_ERROR
   end if
   if (ifmbar /= 0) then
      if (icfe /= 1) then
         write (6,'(/,a)') ' MBAR requires a standard TI run, set icfe to 1'
         DELAYED_ERROR
      end if
      if (ifsc /=0 .and. (bar_l_max > 0.995 .or. bar_l_min < 0.005) ) then
         write (6,'(/,a)') ' Softcore potentials cannot be used with &
               &bar_l_min < 0.005 or bar_l_max > 0.995'
         DELAYED_ERROR
      end if
      if (klambda /= 1) then
         write (6,'(/,a)') ' MBAR requires linear mixing, set klambda to 1'
         DELAYED_ERROR
      end if
   end if
#endif

   if (ilrt /= 0 .and. lrtmask == '') then
      write (6,'(a)') 'Linear Response Theory activated, but lrtmask is not set'
      DELAYED_ERROR
   end if
  
   if (idecomp > 0 .and. (ntr > 0 .or. ibelly > 0)) then
      write(6,'(/,a)') 'IDECOMP is not compatible with NTR or IBELLY'
      DELAYED_ERROR
   end if
   if (icnstph /= 0) then

      if ( icnstph < 0 ) then
         write(6, '(/,a)') 'icnstph must be greater than 0'
         DELAYED_ERROR
      end if
      if ( igb == 0 .and. ipb == 0 .and. icnstph == 1 ) then
         write(6, '(/,a)') 'Constant pH using icnstph = 1 requires &
                           &GB implicit solvent'
         DELAYED_ERROR
      end if
      if ( ntb .eq. 0 .and. icnstph .gt. 1 ) then
         write(6, '(/,a)') 'Constant pH using icnstph = 2 requires &
                           &periodic boundary conditions'
         DELAYED_ERROR
      end if
      if (icfe /= 0) then
         write(6, '(/,a)') &
         'Constant pH and thermodynamic integration are incompatable'
         DELAYED_ERROR
      end if

      if (ntcnstph <= 0) then
         write(6, '(/,a)') 'ntcnstph must be a positive integer.'
         DELAYED_ERROR
      end if

      if (icnstph > 1 .and. mccycles <= 0) then
         write(6, '(/,a)') 'mccycles must be a positive integer.'
         DELAYED_ERROR
      end if

   end if ! icnstph

#ifdef noVIRIAL
   if( ntp > 0 .and. barostat == 1 ) then
      write(6,'(/,a)') 'Error: Berendsen barostat is incompatible with noVIRIAL'
      DELAYED_ERROR
   end if
#endif
   
   !-----------------------------------------------------
   !     ----sanity checks for Ewald
   !-----------------------------------------------------
   
   if( igb == 0 .and. ipb == 0 ) then
      call float_legal_range('skinnb: (nonbond list skin) ', &
            skinnb,skinlo,skinhi)
      
      !  --- Will check on sanity of settings after coords are read in
      !      and the extent of the system is determined.
      
      if(periodic == 1)then
         call float_legal_range('skinnb+cutoffnb: (nonbond list cut) ', &
               skinnb+cutoffnb,zero,sphere)
      end if
      if (ntb==0 .and. use_pme/=0) then
         write(6,'(/,a)') &
         'Using PME with a non-periodic simulation does not make sense. Set either ntb>0 of use_pme=0.'
         DELAYED_ERROR
      end if
      call float_legal_range('a: (unit cell size) ',a,boxlo,boxhi)
      call float_legal_range('b: (unit cell size) ',b,boxlo,boxhi)
      call float_legal_range('c: (unit cell size) ',c,boxlo,boxhi)
      call float_legal_range('alpha: (unit cell angle) ', &
            alpha,anglo,anghi)
      call float_legal_range('beta: (unit cell angle) ', &
            beta,anglo,anghi)
      call float_legal_range('gamma: (unit cell angle) ', &
            gamma,anglo,anghi)
      call int_legal_range('order: (interpolation order) ', &
            order,orderlo,orderhi)
      call opt_legal_range('verbose: ',verbose,0,4)
      call opt_legal_range('netfrc: ',netfrc,0,1)
      call opt_legal_range('nbflag: ',nbflag,0,1)
      call opt_legal_range('nbtell: ',nbtell,0,2)
      call opt_legal_range('ew_type: ',ew_type,0,1)
      call opt_legal_range('vdwmeth: ',vdwmeth,0,2)
      call opt_legal_range('eedmeth: ',eedmeth,1,6)
      call opt_legal_range('ee_type: ',ee_type,1,2)
      call opt_legal_range('maxiter: ',maxiter,1,50)
      call opt_legal_range('indmeth: ',indmeth,0,3)
      call opt_legal_range('fix_quad: ',fix_quad,0,1)
      call float_legal_range('eedtbdns: (erfc table density) ', &
            eedtbdns,denslo,denshi)
   end if  ! ( igb == 0 .and. ipb == 0 )

   if( ntb==2 .and. ipimd==1) then
      write(6,*) 'primitive PIMD is incompatible with NTP ensemble'
      inerr=1
   endif

   if( ntb==2 .and. ipimd==3 ) then
      write(6,*) 'CMD is incompatible with NTP ensemble'
      inerr=1
   endif

   if( ipimd==3 .and. adiab_param>=1.0 ) then
      write(6,*) 'For CMD adiab_param must be <=1'
      inerr=1
   endif

   if( ntb==2 .and. ipimd==4 ) then
      write(6,*) 'RPMD is incompatible with NTP ensemble'
      inerr=1
   endif

   if( ntt/=0 .and. ipimd==4 ) then
      write(6,*) 'RPMD is incompatible with NVT ensemble'
      inerr=1
   endif

   if( ntt/=4 .and. ipimd==2 ) then
      write(6,*) 'NMPIMD requires Nose-Hoover chains (ntt=4)'
      inerr=1
   endif

   if( ntt/=4 .and. ipimd==3 ) then
      write(6,*) 'CMD requires Nose-Hoover chains (ntt=4)'
      inerr=1
   endif

#if !defined(MPI)
   if( ineb > 0 ) then
      write(6,*) 'NEB requires MPI'
      inerr=1
   endif
#endif

#ifdef LES
   if( ineb > 0 ) then
      write(6,*) 'NEB no longer works with LES: use multiple groups and a groupfile instead'
      inerr=1
   endif
#endif

   if( ineb > 0 .and. ipimd > 0 ) then
      write(6,*) 'ineb>0 and ipimd>0 are incompatible options'
      inerr=1
   endif

   if( iamoeba == 1 )then
#ifdef LES
      write(6,*)'amoeba is incompatible with LES'
      inerr=1
#endif

      if( ntc > 1 ) then
         write(6,*) 'SHAKE (ntc>1) and amoeba are incompatible options'
         inerr=1
      end if
      if( ntp > 1 .and. beeman_integrator > 0 ) then
         write(6,*) 'ntp>1 is not consistent with the beeman integrator'
         inerr=1
      end if
      if ( igb /= 0 ) then
         write (6,'(/,a)') 'amoeba (iamoeba=1) is incompatible with GB (igb>0)'
         inerr=1
      end if
      if ( ipb /= 0 ) then
         write (6,'(/,a)') 'amoeba (iamoeba=1) is incompatible with PB (ipb>0)'
         inerr=1
      end if
#ifdef MPI
      if (numtasks > 1) then
         write(6, '(/,a)') 'amoeba (iamoeba=1) cannot be run in parallel in &
                           &sander'
         inerr=1
      end if
#endif

   end if

   ! ---WARNINGS:
   
   if ( ibelly == 1 .and. igb == 0 .and. ipb == 0 .and. ntb /= 0 ) then
      write(6,'(/,a,/,a,/,a)') &
            'Warning: Although EWALD will work with belly', &
            '(for equilibration), it is not strictly correct!'
   end if
   
   if (inerr == 1) then
      write(6, '(/,a)') ' *** input error(s)'
      FATAL_ERROR
   end if
   
   ! Load the restrained atoms (ntr=1) or the belly atoms (ibelly=1)
   ! or atoms for targeted md (itgtmd=1). Selections are read from
   ! &cntrl variables or, if these are not defined, it falls back to
   ! the old group input format.
   
   if(mtmd /= 'mtmd') then
     itgtmd = 2
     ntr = 0
     emtmd = 0.0d0
   end if
   konst = ntr > 0
   dotgtmd = itgtmd > 0
   belly = .false.
   natc = 0
   ngrp = 0
   natbel = 0
   nattgtfit = 0  ! number of atoms for tgtmd fitting (=overlap)
   nattgtrms = 0  ! number of atoms for tgtmd rmsd calculation
   nrc = 0
   if(konst.or.dotgtmd) then
      
      ! inserted here to fix the bug that coords are not available
      ! yet when distance based selection (<,>) is requested
#ifdef LES
#  ifdef MPI
      call getcor(nr,x(lcrd),x(lvel),x(lforce),ntx,box,irest,t,temp0les,.FALSE.,solvph)
#  else
      call getcor(nr,x(lcrd),x(lvel),x(lforce),ntx,box,t,.FALSE.)
#  endif
#else
      call AMOEBA_check_newstyle_inpcrd(inpcrd,newstyle)
      if ( newstyle )then
         call AM_RUNMD_get_coords(natom,t,irest,ntb,x(lcrd),x(lvel))
      else
         if( irest == 1 .and. beeman_integrator > 0 ) then
            write(6,*) 'Cannot do a beeman_integrator restart with old-style coordinates'
            call mexit(6,1)
         end if
#  ifdef MPI
         call getcor(nr,x(lcrd),x(lvel),x(lforce),ntx,box,irest,t,temp0,.FALSE.,solvph)
#  else
         call getcor(nr,x(lcrd),x(lvel),x(lforce),ntx,box,t,.FALSE.)
#  endif
      endif
#endif /* LES */
      
      if(itgtmd == 2) then
        call mtmdcall(emtmd,x(lmtmd01),ix(imtmd02),x(lcrd),x(lforce),ih(m04),ih(m02),ix(i02),&
                    ih(m06),x(lmass),natom,nres,'READ')  
      else
#ifndef API
        if (konst) write(6,9408)
        if (dotgtmd) write(6,9409)
#endif
        call rdrest(natom,ntrx,refc,x(lcrdr))
        !close(10)

        ! VH - tgtmd change: preferably call atommask() instead of rgroup()
        if (konst) then
          if( len_trim(restraintmask) <= 0 ) then
              call rgroup(natom,natc,nres,ngrp,ix(i02),ih(m02),ih(m04), &
                          ih(m06),ih(m08),ix(icnstrgp),jgroup,indx,irespw, &
                          npdec,x(l60),konst,dotgtmd,belly,idecomp,5,.true.)
          else
              call atommask( natom, nres, 0, ih(m04), ih(m06), &
                ix(i02), ih(m02), x(lcrd), restraintmask, ix(icnstrgp) )
  
              ! for now, emulate the "GATHER ALL THE CONSTRAINED ATOMS TOGETHER"
              ! section of rgroup(); later, the various masks should be done
              ! differently, i.e. without the "gather", as in the following:
              !     x(l60:l60+natom-1) = restraint_wt
              !     natc = sum(ix(icnstrgp:icnstrgp+natom-1))
  
              natc = 0
              do i=1,natom
                if( ix(icnstrgp-1+i) <= 0 ) cycle
                natc = natc + 1
                ix(icnstrgp-1+natc) = i
                x(l60-1+natc) = restraint_wt
              end do
              write(6,'(a,a,a,i5,a)') '     Mask ', &
              restraintmask(1:len_trim(restraintmask)), ' matches ',natc,' atoms'
          end if
        end if
        nrc = natc
        
        if (itgtmd == 1) then
          if (len_trim(tgtfitmask) <= 0 .and. len_trim(tgtrmsmask) <= 0) then
              ! the following if-endif can be deleted when we stop
              ! supporting rgroup()
              if (konst) then
                ! cannot do both ntr and tgtmd together using old group format
                write(6,'(/2x,a)') 'NTR must be 0 for targeted MD (TGTMD=1)'
                call mexit(6,1)
              else  ! the following only for backward compatibility
                call rgroup(natom,natc,nres,ngrp,ix(i02),ih(m02),ih(m04), &
                            ih(m06),ih(m08),ix(icnstrgp),jgroup,indx,irespw, &
                            npdec,x(l60),konst,dotgtmd,belly,idecomp,5,.true.)
                ! tgtmd atoms are now stored in nattgt, igroup -> icnstrgp
                nattgtfit = natc
                nattgtrms = natc
                do i=1,nattgtfit
                    ix(itgtfitgp-1+i) = ix(icnstrgp-1+i)
                    ix(itgtrmsgp-1+i) = ix(icnstrgp-1+i)
                end do
              end if
          else
              if (ntr == 0) then  ! read tgtfitmask only if ntr=1 
                ! read in atom group for tgtmd fitting (=overlap region)
                call atommask( natom, nres, 0, ih(m04), ih(m06), &
                    ix(i02), ih(m02), x(lcrd), tgtfitmask, ix(itgtfitgp) )
                ! see comments above (for ntr) for the following reduction cycle
                nattgtfit = 0
                do i=1,natom
                  if( ix(itgtfitgp-1+i) <= 0 ) cycle
                  nattgtfit = nattgtfit + 1
                  ix(itgtfitgp-1+nattgtfit) = i
                end do
                write(6,'(a,a,a,i5,a)')  &
                '     Mask "', tgtfitmask(1:len_trim(tgtfitmask)-1),  &
                '" matches ',nattgtfit,' atoms'
              end if
              ! read in atom group for tgtmd rmsd calculation
              call atommask( natom, nres, 0, ih(m04), ih(m06), &
                ix(i02), ih(m02), x(lcrd), tgtrmsmask, ix(itgtrmsgp) )
              nattgtrms = 0
              do i=1,natom
                if( ix(itgtrmsgp-1+i) <= 0 ) cycle
                nattgtrms = nattgtrms + 1
                ix(itgtrmsgp-1+nattgtrms) = i
              end do
              write(6,'(a,a,a,i5,a)')  &
              '     Mask "', tgtrmsmask(1:len_trim(tgtrmsmask)-1),  &
              '" matches ',nattgtrms,' atoms'
          end if
        end if
      
      end if ! (itgtmd == 2)

   end if  ! (konst.or.dotgtmd)

   if (ineb>0) then
      ! carlos: read in fitmask and rmsmask info for NEB, just as done for tgtmd
      ! init last_neb_atom, which is used to determine the limits for the
      ! communication of neighbor coordinates (to reduce size for explicit
      ! water)
      last_neb_atom = 0

      if (ntr /= 0) then
        write(6,'(/2x,a)') 'cannot use NEB with ntr restraints'
! CARLOS: WHY NOT? SHOULD BE OK.
! potential for user error is restrained region overlaps with NEB region. would
! blow up.
        call mexit(6,1)
      else
         ! read in atom group for fitting (=overlap region)
         call atommask( natom, nres, 0, ih(m04), ih(m06), &
            ix(i02), ih(m02), x(lcrd), tgtfitmask, ix(itgtfitgp) )
         ! see comments above (for ntr) for the following reduction cycle
         nattgtfit = 0
         do i=1,natom
           if( ix(itgtfitgp-1+i) <= 0 ) cycle
           nattgtfit = nattgtfit + 1
           ix(itgtfitgp-1+nattgtfit) = i
           if (i.gt.last_neb_atom) last_neb_atom = i
         end do
         write(6,'(a)') "The following selection will be used for NEB structure fitting"
         write(6,'(a,a,a,i5,a)')  &
         '     Mask "', tgtfitmask(1:len_trim(tgtfitmask)-1),  &
         '" matches ',nattgtfit,' atoms'
      end if
      ! read in atom group for tgtmd rmsd calculation
      call atommask( natom, nres, 0, ih(m04), ih(m06), &
         ix(i02), ih(m02), x(lcrd), tgtrmsmask, ix(itgtrmsgp) )
      nattgtrms = 0
      do i=1,natom
        if( ix(itgtrmsgp-1+i) <= 0 ) cycle
        nattgtrms = nattgtrms + 1
        ix(itgtrmsgp-1+nattgtrms) = i
        if (i.gt.last_neb_atom) last_neb_atom = i
      end do
      write(6,'(a)') "The following selection will be used for NEB force application"
      write(6,'(a,a,a,i5,a)')  &
      '     Mask "', tgtrmsmask(1:len_trim(tgtrmsmask)-1),  &
      '" matches ',nattgtrms,' atoms'
      write(6,'(/2x,a,i6)') "Last atom in NEB fitmask or rmsmask is ",last_neb_atom

      if (nattgtrms<=0 .or. nattgtfit <= 0) then
        write(6,'(/2x,a)') 'NEB requires use of tgtfitmask and tgtrmsmask'
        call mexit(6,1)
      endif
   endif
   
   ! dotgtmd may be false here even if doing tgtmd
   ! this is so belly info is read properly? following existing KONST code
   
   dotgtmd=.false.
   konst = .false.
   belly = ibelly > 0
   ngrp = 0
   if(belly) then
      ! inserted here to fix the bug that coords are not available
      ! yet when distance based selection (<,>) is requested
#ifdef LES
#  ifdef MPI
      call getcor(nr,x(lcrd),x(lvel),x(lforce),ntx,box,irest,t,temp0les,.FALSE.,solvph)
#  else
      call getcor(nr,x(lcrd),x(lvel),x(lforce),ntx,box,t,.FALSE.)
#  endif
#else
#  ifdef MPI
      call getcor(nr,x(lcrd),x(lvel),x(lforce),ntx,box,irest,t,temp0,.FALSE.,solvph)
#  else
      call getcor(nr,x(lcrd),x(lvel),x(lforce),ntx,box,t,.FALSE.)
#  endif
#endif /* LES */
#ifndef API
      write(6,9418)
#endif /* API */
      if( len_trim(bellymask) <= 0 ) then
         call rgroup(natom,natbel,nres,ngrp,ix(i02),ih(m02),ih(m04),ih(m06), &
                     ih(m08),ix(ibellygp),jgroup,indx,irespw,npdec, &
                     x(l60),konst,dotgtmd,belly,idecomp,5,.true.)
      else
         call atommask( natom, nres, 0, ih(m04), ih(m06), &
            ix(i02), ih(m02), x(lcrd), bellymask, ix(ibellygp) )
         natbel = sum(ix(ibellygp:ibellygp+natom-1))
         write(6,'(a,a,a,i5,a)') '     Mask ', &
            bellymask(1:len_trim(bellymask)), ' matches ',natbel,' atoms'
      end if
   end if
   call setvar(ix,belly)

   !  see if the user has input a noshakemask string, and process it:
   natnos = 0
   if( len_trim(noshakemask) > 0 ) then
      call atommask( natom, nres, 0, ih(m04), ih(m06), &
         ix(i02), ih(m02), x(lcrd), noshakemask, noshakegp )
      natnos = sum(noshakegp(1:natom))
      write(6,*)
      write(6,'(a,a,a,i5,a)') 'Noshake mask ', &
         noshakemask(1:len_trim(noshakemask)), ' matches ',natnos,' atoms'
      call setnoshake(ix,noshakegp,ntc,num_noshake)
      if( ntf > 1 ) then
         write(6,'(a)') '   Setting ntf to 1'
         ntf = 1
      end if
   end if
   
   ! GMS ------------------------------------
   !  Check for 'iwrap_mask', and process it.
   ! ----------------------------------------
   n_iwrap_mask_atoms = 0
   if( len_trim(iwrap_mask) > 0 .and. iwrap == 2) then
      call atommask( natom, nres, 0, ih(m04), ih(m06), &
         ix(i02), ih(m02), x(lcrd), iwrap_mask, iwrap_maskgp )
      ! iwrap_maskgp is a natom long integer array, with elements:
      ! 0 --> atom is not in iwrap_mask
      ! 1 --> atom is in iwrap_mask
      n_iwrap_mask_atoms = sum(iwrap_maskgp(1:natom))
      write(6,*)
      write(6,'(a,a,a,i5,a)') 'Wrap mask ', &
         iwrap_mask(1:len_trim(iwrap_mask)), ' matches ',n_iwrap_mask_atoms,' atoms:'
      ! Set an array to store the atom numbers of the atoms
      ! in the iwrap_mask
      allocate(iwrap_mask_atoms(n_iwrap_mask_atoms), stat=ier)
      REQUIRE(ier == 0)

      j = 0
      do i=1,natom
        if( iwrap_maskgp(i)>0 ) then
          j = j+1
          iwrap_mask_atoms(j) = i
        end if
      end do
      
      write(6,'(10i5)') (iwrap_mask_atoms(i),i=1,n_iwrap_mask_atoms)
   end if

#ifdef MPI /* SOFT CORE */
   ! lower charges if a crgmask is set
   if ( len_trim(crgmask) > 0 ) then
      call atommask( natom, nres, 0, ih(m04), ih(m06), &
         ix(i02), ih(m02), x(lcrd), crgmask, crggp )
      write(6,'(a,a,a,i5,a)') 'Zero-Charge Mask ',crgmask(1:len_trim(crgmask)), ' matches ',sum(crggp(1:natom)),' atoms'
      call remove_charges(crggp, natom, x(l15))
   end if
#endif

   konst = .false.
   belly = .false.
#ifndef API /* IDECOMP cannot be nonzero in the API */
   if(idecomp > 0) then
      write(6,9428)
      call rgroup(natom,ntmp,nres,ngrp,ix(i02),ih(m02),ih(m04),ih(m06), &
                  ih(m08),ix(ibellygp),jgroup,indx,irespw,npdec, &
                  x(l60),konst,dotgtmd,belly,idecomp,5,.true.)
   end if
#endif
   
   if( ibelly > 0 .and. (igb > 0 .or. ipb /= 0) ) then
      
      !          ---here, the only allowable belly has just the first
      !             NATBEL atoms in the moving part.  Check to see that this
      !             requirement is satisfied:
      
      do i=natbel+1,natom
         if( ix(ibellygp+i-1) /= 0 ) then
            write(6,*) 'When igb>0, the moving part must be at the'
            write(6,*) '   start of the molecule.  This does not seem'
            write(6,*) '   to be the case here.'
            write(6,*) 'natbel,i,igroup(i) = ' &
                  ,natbel,i,ix(ibellygp+i-1)
            call mexit(6,1)
         end if
      end do
   end if
   

   !     ----- CALCULATE THE SQUARE OF THE BOND PARAMETERS FOR SHAKE
   !           THE PARAMETERS ARE PUT SEQUENTIALLY IN THE ARRAY CONP -----
   
   do i=1,nbonh + nbona + nbper
      j = ix(iicbh+i-1)
      x(l50+i-1) = req(j)**2
   end do
   
#ifdef MPI
      if( icfe /= 0 ) then

      !  use the masses of the prmtop file for the first group for both groups:
      !  [only the master nodes communicate here, since non-master nodes
      !   have not yet allocated space]
      ! This leads to problems for dual topology runs, and is therefore skipped
      ! if ifsc is set to one, the masses from both prmtop files are used
         if (ifsc == 0) then
            call mpi_bcast(x(lmass),natom,MPI_DOUBLE_PRECISION,0,commmaster,ierr)
            call mpi_bcast(x(lwinv),natom,MPI_DOUBLE_PRECISION,0,commmaster,ierr)            
            call mpi_bcast(x(l75),natom,MPI_DOUBLE_PRECISION,0,commmaster,ierr)
         end if
         tmass = sum(x(lmass:lmass+natom-1))
         tmassinv = 1.d0/tmass

      !  next, do a minimal sanity check that the SHAKE parameters are
      !  consistent on the two processors:

         ! For Softcore this might be allowed
         ! Put a better check here later
         if( ntc == 2 .and. ifsc == 0) then
            partner = ieor(masterrank,1)
            call mpi_sendrecv( nbonh, 1, MPI_INTEGER, partner, 5, &
                               nbonh_c, 1, MPI_INTEGER, partner, 5, &
                               commmaster, ist, ierr )
            call mpi_sendrecv( num_noshake, 1, MPI_INTEGER, partner, 5, &
                               num_noshake_c, 1, MPI_INTEGER, partner, 5, &
                               commmaster, ist, ierr )
            if (qmmm_nml%ifqnt .and. qmmm_nml%qmshake == 0) then
               ! qtw - if qmshake=0, we need to check the QM atoms
               call mpi_sendrecv( &
                       qmmm_struct%nquant, 1, MPI_INTEGER, partner, 5, &
                       nquant_c, 1, MPI_INTEGER, partner, 5, &
                       commmaster, ist, ierr)
               call mpi_sendrecv( &
                       qmmm_struct%noshake_overlap, 1, MPI_INTEGER,partner, 5, &
                       noshake_overlap_c, 1, MPI_INTEGER, partner, 5, &
                       commmaster, ist, ierr)
               if ( (qmmm_struct%nquant-qmmm_struct%noshake_overlap) /= &
                    (nquant_c-noshake_overlap_c) ) then
                  call sander_bomb('mdread2', &
                       'QMMM: NOSHAKE lists are not match in two groups!', &
                       'try noshakemask in cntrl to match the noshake list')
               end if
            else if( nbonh - num_noshake /= nbonh_c - num_noshake_c ) then
               write(6,*) 'SHAKE lists are not compatible in the two groups!'
               call mexit(6,1)
            end if
         else if( ntc == 3 ) then
            write(6,*) 'ntc = 3 is not compatible with icfe>0'
            call mexit(6,1)
         end if

      end if
#endif

   if ( iamoeba /= 1 )then
      if( igb == 0 .and. ipb == 0 ) &
         call init_extra_pts( &
         ix(iibh),ix(ijbh),ix(iicbh), &
         ix(iiba),ix(ijba),ix(iicba), &
         ix(i24),ix(i26),ix(i28),ix(i30), &
         ix(i32),ix(i34),ix(i36),ix(i38), &
         ix(i40),ix(i42),ix(i44),ix(i46),ix(i48), &
         ix(i50),ix(i52),ix(i54),ix(i56),ix(i58), &
         ih(m06),ix,x,ix(i08),ix(i10), &
         nspm,ix(i70),x(l75),tmass,tmassinv,x(lmass),x(lwinv),req)
   endif
  
   !  DEBUG input; force checking
#ifdef API
   call load_debug
#else
   call load_debug(5)
#endif

   return
   ! -------------------------------------------------------------------------
   ! Standard format statements:
   
#ifndef API
   9328 format(/80('-')/,'   2.  CONTROL  DATA  FOR  THE  RUN',/80('-')/)
   9408 format(/4x,'LOADING THE CONSTRAINED ATOMS AS GROUPS',/)
   9409 format(/4x,'LOADING THE TARGETED MD ATOMS AS GROUPS',/)
   9418 format(/4x,'LOADING THE BELLY ATOMS AS GROUPS',/)
   9428 format(/4x,'LOADING THE DECOMP ATOMS AS GROUPS',/)
   9008 format(a80)
#endif

#undef FATAL_ERROR
#undef DELAYED_ERROR
