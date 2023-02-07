#include "copyright.h"
#include "../include/dprec.fh"

#ifndef PBSA
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Emit the final minimization report 
!-----------------------------------------------------------------------
!     --- REPORT_MIN_RESULTS ---
!-----------------------------------------------------------------------
! Find the maximum component of the gradient (grdmax),
! emit this and other gradient details (printe),
! optionally do pseudocontact shift constraints,
! optionally decompose energies for mm_pbsa,
! and emit nmr related information (nmrptx).

subroutine report_min_results( nstep, gradient_rms, coordinates, &
      forces, energies, igraph, xx, ix, ih )

   use state
   use decomp, only : checkdec, printdec, printpdec
   use file_io_dat
#ifdef RISMSANDER
   use sander_rism_interface, only: rismprm, rism_solvdist_thermo_calc
#endif
   implicit none

   integer, intent(in)             :: nstep
   _REAL_,  intent(in)             :: gradient_rms
   _REAL_,  intent(in)             :: coordinates(*)
   _REAL_,  intent(in)             :: forces(*)
   type(state_rec), intent(in)     :: energies
   character(len=4), intent(in)    :: igraph(*)    ! atom name map
   _REAL_,  intent(inout)          :: xx(*)        ! real dynamic memory
   integer, intent(inout)          :: ix(*)        ! integer dynamic memory
   character(len=4), intent(inout) :: ih(*)        ! hollerith dynamic memory

#  include "box.h"
#  include "extra.h"
#  include "../include/md.h"
#  include "../include/memory.h"
#  include "nmr.h"
#  include "tgtmd.h"
#  include "multitmd.h"

   character(len=4) :: atom_name_of_gmax
   integer          :: atom_number_of_gmax
   _REAL_           :: gradient_max
   _REAL_ emtmd

#ifdef RISMSANDER
   !ensure that RISM thermodynamics are print IFF the user wants them.
   !This part must be run in parallel
   if ( rismprm%irism == 1 .and. rismprm%write_thermo==1)&
        call rism_solvdist_thermo_calc(.false.,0)
#endif /*RISMSANDER*/
   if (master) then
      write(6, '(/ /20x,a,/)' ) 'FINAL RESULTS'
      call grdmax( forces, gradient_max, atom_number_of_gmax )
      atom_name_of_gmax = igraph(atom_number_of_gmax)
      if (imin /= 5) rewind(MDINFO_UNIT)
      call printe( nstep, gradient_rms, gradient_max, energies, &
             atom_number_of_gmax, atom_name_of_gmax )
      if (idecomp > 0) call checkdec(idecomp)
      if (idecomp == 1 .or. idecomp == 2) call printdec(ix)
      if (idecomp == 3 .or. idecomp == 4) call printpdec(ix)
      if (nmropt > 0) then
         if (iredir(7) /= 0) call pcshift(-1,coordinates,forces)
         call nmrptx(6)
         call nmrptx(MDINFO_UNIT)
         call ndvptx(coordinates,forces,ih(m04),ih(m02),ix(i02),nres, &
               xx(l95),natom,xx(lwinv),xx(lnmr01),ix(inmr02),6)
      end if
      if (itgtmd == 2) then
        emtmd = 0.0d0
        call mtmdcall(emtmd,xx(lmtmd01),ix(imtmd02),coordinates,forces,ih(m04),ih(m02),ix(i02),&
                                      ih(m06),xx(lmass),natom,nres,'PRNT')
      end if
      if (imin /= 5) call amflsh(MDINFO_UNIT)
   end if

   return
end subroutine report_min_results


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Emit a minimization progress report 
!-----------------------------------------------------------------------
!     --- REPORT_MIN_PROGRESS ---
!-----------------------------------------------------------------------
! Find the maximum component of the gradient (grdmax),
! emit this and other gradient details (printe),
! and emit nmr related information (nmrptx).

subroutine report_min_progress( nstep, gradient_rms, forces, energies, &
                                igraph, charge )

   use state
   use file_io_dat
   use crg_reloc, only : ifcr, crprintcharges, cr_print_charge
#ifdef RISMSANDER
   use sander_rism_interface, only: rismprm, rism_solvdist_thermo_calc
#endif
   implicit none

   integer, intent(in)           :: nstep
   _REAL_,  intent(in)           :: gradient_rms
   _REAL_,  intent(in)           :: forces(*)
   type(state_rec), intent(in)   :: energies
   character(len=4), intent(in)  :: igraph(*)    ! atom name map
   _REAL_, intent(in)            :: charge(*)

#  include "extra.h"
#  include "../include/md.h"
#  include "nmr.h"

   character(len=4) :: atom_name_of_gmax
   integer          :: atom_number_of_gmax
   _REAL_           :: gradient_max

#ifdef RISMSANDER
   !ensure that RISM thermodynamics are print IFF the user wants them.
   !This part must be run in parallel
   if ( rismprm%irism == 1 .and. rismprm%write_thermo==1)&
        call rism_solvdist_thermo_calc(.false.,0)
#endif /*RISMSANDER*/
   if (master) then
      call grdmax( forces, gradient_max, atom_number_of_gmax )
      atom_name_of_gmax = igraph(atom_number_of_gmax)
      if (imin /= 5) rewind(MDINFO_UNIT)
      call printe( nstep, gradient_rms, gradient_max, energies, &
             atom_number_of_gmax, atom_name_of_gmax )
      if (nmropt > 0) then
         call nmrptx(6)
         call nmrptx(MDINFO_UNIT)
      end if
      if ( ifcr > 0 .and. crprintcharges > 0 ) then
         call cr_print_charge( charge, nstep ) 
      end if
      if (imin /= 5) call amflsh(MDINFO_UNIT)
   end if

   return
end subroutine report_min_progress
#endif /*ifndef PBSA*/

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Compute the maximum gradient component and the corresponding atom
!-----------------------------------------------------------------------
!     --- GRDMAX ---
!-----------------------------------------------------------------------

subroutine grdmax( gradient, gradient_max, atom_number_of_gmax )

   use constants, only : zero
   implicit none
   _REAL_,  intent(in)  :: gradient(*)
   !     This is actually two-dimensional (3,natoms), but to enable
   !     vectorization on IA32 SSE platforms they are treated as
   !     one-dimensional; this may also improve software pipelining !

   _REAL_,  intent(out) :: gradient_max
   integer, intent(out) :: atom_number_of_gmax

#  include "../include/memory.h"

   integer :: i
   _REAL_  :: gi

   gradient_max = ZERO
   atom_number_of_gmax = 1
   do i = 1,3*natom
      gi = abs(gradient(i))
      if (gi > gradient_max) then
         gradient_max = gi
         atom_number_of_gmax = i
      end if
   end do
   atom_number_of_gmax = (atom_number_of_gmax - 1)/3 + 1

   return
end subroutine grdmax


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Print status of a minimization calculation.
!-----------------------------------------------------------------------
!     --- PRINTE ---
!-----------------------------------------------------------------------
! Output the step number, the gradient rms, the max gradient component,
! and the atom label and atom number of the max gradient component.

subroutine printe( nstep, gradient_rms, gradient_max, ene, &
      atom_number_of_gmax, atom_name_of_gmax )
   
#ifdef APBS
   use file_io_dat
#endif
#ifdef RISMSANDER
   use sander_rism_interface, only : rismprm, RISM_NONE, RISM_FULL, RISM_INTERP,&
        rism_calc_type, rism_thermo_print
#endif
   use qmmm_module, only : qmmm_nml
   use cns_xref
   use state ! Access to energy_rec
   use charmm_mod, only : charmm_active
   use crg_reloc, only : ifcr
   use emap,only : temap,scemap
   use ff11_mod, only : cmap_active
   use sebomd_module, only : sebomd_obj

   implicit none
   integer, intent(in)          :: nstep
   _REAL_,  intent(in)          :: gradient_rms
   _REAL_,  intent(in)          :: gradient_max
   type(state_rec), intent(in)  :: ene
   integer, intent(in)          :: atom_number_of_gmax
   character(len=4), intent(in) :: atom_name_of_gmax

#  include "../include/md.h"
#  include "ew_mpole.h"
#  include "ew_cntrl.h"
#  include "tgtmd.h"

   _REAL_  epot,enonb,enele,ehbond,ebond,eangle,edihed,enb14,eel14,egb,epb
   _REAL_  econst,epolar,aveper,aveind,avetot,esurf,edisp,diprms,dipiter
   _REAL_  dipole_temp,escf,dvdl,enemap
#ifdef RISMSANDER
   _REAL_ :: erism
   _REAL_ :: pot_array(potential_energy_rec_len)
#endif /*RISMSANDER*/
   _REAL_ :: ect

   !     ----- Extract Energies. -----

   epot        = ene%pot%tot
   enonb       = ene%pot%vdw
   enele       = ene%pot%elec
   ehbond      = ene%pot%hbond
   egb         = ene%pot%gb 
   epb         = ene%pot%pb
   ebond       = ene%pot%bond
   eangle      = ene%pot%angle
   edihed      = ene%pot%dihedral
   enb14       = ene%pot%vdw_14
   eel14       = ene%pot%elec_14
   econst      = ene%pot%constraint
   epolar      = ene%pot%polar
   dvdl        = ene%pot%dvdl
   aveper      = ene%aveper
   aveind      = ene%aveind
   avetot      = ene%avetot
   esurf       = ene%pot%surf
   diprms      = ene%diprms
   dipiter     = ene%dipiter
   dipole_temp = ene%dipole_temp
   
   escf    = ene%pot%scf
   edisp   = ene%pot%disp
   enemap   = ene%pot%emap
#ifdef RISMSANDER
   erism   = ene%pot%rism
#endif /*RISMSANDER*/
   ect = ene%pot%ct


   

   write(6,9018)
   write(6,9028) nstep, epot, gradient_rms, gradient_max, &
         atom_name_of_gmax, atom_number_of_gmax
   write(6,9038) ebond,eangle,edihed
! CHARMM SPECIFIC ENERGY TERMS
   if (charmm_active) write(6,9039) ene%pot%angle_ub, &
                                    ene%pot%imp,      &
                                    ene%pot%cmap

#ifdef RISMSANDER
   if( igb == 0 .and. ipb == 0 .and. rismprm%irism == 0) then
#else
   if( igb == 0 .and. ipb == 0 ) then
#endif
      write(6,9048) enonb,enele,ehbond
   else if ( igb == 10 .or. ipb /= 0 ) then
      write(6,9050) enonb,enele,epb
#ifdef APBS
   else if ( igb == 6 .and. mdin_apbs ) then
      write(6,9050) enonb,enele,epb
#endif /* APBS */
#ifdef RISMSANDER
   else if(rismprm%irism == 1 )then
      write(6,9051) enonb,enele,erism
#endif
   else
      write(6,9049) enonb,enele,egb
   end if
   write(6,9058) enb14,eel14,econst

   if ( ifcr /= 0 ) write(6,9099) ect

!  wxw: EMAP energy
   if (temap) then
      write (6,9062) enemap,scemap
   endif

   if (qmmm_nml%ifqnt) then
     !write the SCF energy
     if (qmmm_nml%qmtheory%PM3) then
        if (qmmm_nml%qmmm_int == 3) then
           write(6,9090) escf ! PM3-MM*
        else if (qmmm_nml%qmmm_int == 4) then
           write(6,9096) escf ! PM3/MMX2
        else
           write(6,9080) escf
        end if
     else if (qmmm_nml%qmtheory%AM1) then
        write(6,9081) escf
     else if (qmmm_nml%qmtheory%AM1D) then
        write(6,9981) escf
     else if (qmmm_nml%qmtheory%MNDO) then
        write(6,9082) escf
     else if (qmmm_nml%qmtheory%MNDOD) then
        write(6,9982) escf
     else if (qmmm_nml%qmtheory%PDDGPM3) then
        write(6,9083) escf
     else if (qmmm_nml%qmtheory%PDDGMNDO) then
        write(6,9084) escf
     else if (qmmm_nml%qmtheory%PM3CARB1) then
        write(6,9085) escf
     else if (qmmm_nml%qmtheory%DFTB) then
        write(6,9086) escf
     else if (qmmm_nml%qmtheory%RM1) then
        write(6,9087) escf
     else if (qmmm_nml%qmtheory%PDDGPM3_08) then
        write(6,9088) escf
     else if (qmmm_nml%qmtheory%PM6) then
        write(6,9089) escf
     else if (qmmm_nml%qmtheory%PM3ZNB) then
        write(6,9091) escf
     else if (qmmm_nml%qmtheory%EXTERN) then
        write(6,9092) escf
     else if (qmmm_nml%qmtheory%PM3MAIS) then
        write(6,9093) escf
     else
        write(6,'(" ERROR - UNKNOWN QM THEORY")')
     end if
   end if

   if (sebomd_obj%do_sebomd) then
     write(6,9200) sebomd_obj%esebomd
   end if
#ifdef PUPIL_SUPPORT
   write(6,9900) escf
#endif
   if( gbsa > 0 ) write(6,9077) esurf
   if (igb == 10 .or. ipb /= 0) write(6,9074) esurf,edisp
   call write_cns_xref_min_energies( ene )
#ifdef APBS
   if (igb == 6 .and. mdin_apbs ) write(6,9069) esurf
#endif /* APBS */
      if (cmap_active .and. ipol > 0 ) then
          write(6,9066) epolar, ene%pot%cmap
      else
          if (cmap_active) write(6,9067) ene%pot%cmap
          if (epolar /= 0.0) write(6,9068) epolar
      end if
   if (econst /= 0.0) write(6,9078) epot-econst
   if ( dvdl /= 0.d0) write(6,9100) dvdl
   if (induced > 0.and.indmeth < 3) write(6,9190)diprms,dipiter
   if (induced > 0.and.indmeth == 3) write(6,9191)diprms, &
         dipole_temp
   if (itgtmd == 1) then
      write (6,'(a,f8.3)') "Current RMSD from reference: ",rmsdvalue
      write (6,'(a,f8.3)') "Current target RMSD:         ",tgtrmsd
   end if

   call amflsh(6)

   !     ----- SEND IT TO THE INFO FILE -----

   if (imin /= 5) then
      write(7,9018)
      write(7,9028) nstep, epot, gradient_rms, gradient_max, &
            atom_name_of_gmax, atom_number_of_gmax
      write(7,9038) ebond,eangle,edihed
! CHARMM SPECIFIC ENERGY TERMS
   if (charmm_active) write(7,9039) ene%pot%angle_ub, &
                                    ene%pot%imp,      &
                                    ene%pot%cmap


#ifdef RISMSANDER
      if( igb == 0 .and. ipb == 0 .and. rismprm%irism == 0) then
#else
      if( igb == 0 .and. ipb == 0 ) then
#endif
         write(7,9048) enonb,enele,ehbond
      else if ( igb == 10 .or. ipb /= 0 ) then
         write(7,9050) enonb,enele,epb
#ifdef APBS
      else if ( igb == 6 .and. mdin_apbs ) then
         write(7,9050) enonb,enele,epb
#endif /* APBS */
#ifdef RISMSANDER
      else if ( rismprm%irism == 1 ) then
         write(7,9051) enonb,enele,erism
#endif
      else
         write(7,9049) enonb,enele,egb
      end if
      write(7,9058) enb14,eel14,econst

!  wxw: EMAP energy
   if (temap) then
      write (7,9062) enemap,scemap
   endif

      if (qmmm_nml%ifqnt) then
        !write the SCF energy
        if (qmmm_nml%qmtheory%PM3) then
          if (qmmm_nml%qmmm_int == 3) then
             write(7,9090) escf ! PM3-MM*
          else if (qmmm_nml%qmmm_int == 4) then
             write(7,9096) escf ! PM3/MMX2
          else
             write(7,9080) escf
          end if
        else if (qmmm_nml%qmtheory%AM1) then
           write(7,9081) escf
        else if (qmmm_nml%qmtheory%AM1D) then
           write(7,9981) escf
        else if (qmmm_nml%qmtheory%MNDO) then
           write(7,9082) escf
        else if (qmmm_nml%qmtheory%MNDOD) then
           write(7,9982) escf           
        else if (qmmm_nml%qmtheory%PDDGPM3) then
           write(7,9083) escf
        else if (qmmm_nml%qmtheory%PDDGMNDO) then
           write(7,9084) escf
        else if (qmmm_nml%qmtheory%PM3CARB1) then
           write(7,9085) escf
        else if (qmmm_nml%qmtheory%DFTB) then
           write(7,9086) escf
        else if (qmmm_nml%qmtheory%RM1) then
           write(7,9087) escf
        else if (qmmm_nml%qmtheory%PDDGPM3_08) then
           write(7,9088) escf
        else if (qmmm_nml%qmtheory%PM6) then
           write(7,9089) escf
        else if (qmmm_nml%qmtheory%PM3ZNB) then
           write(7,9091) escf
        else if (qmmm_nml%qmtheory%EXTERN) then
           write(7,9092) escf
        else if (qmmm_nml%qmtheory%PM3MAIS) then
           write(7,9093) escf
        else
           write(7,'(" ERROR - UNKNOWN QM THEORY")')
        end if
      end if

      if (sebomd_obj%do_sebomd) then
        write(7,9200) sebomd_obj%esebomd
      end if
#ifdef PUPIL_SUPPORT
      write(7,9900) escf
#endif
      if( gbsa > 0 ) write(7,9077) esurf
      if ( igb == 10 .or. ipb /= 0 ) write(7,9074) esurf,edisp
#ifdef APBS
      if (igb == 6 .and. mdin_apbs ) write(7,9069) esurf
#endif /* APBS */
! FF11 CMAP SPECIFIC ENERGY TERMS
      if (cmap_active .and. epolar /= 0.0 ) then
          write(7,9066) epolar, ene%pot%cmap
      else
          if (cmap_active) write(7,9067) ene%pot%cmap
          if (epolar /= 0.0) write(7,9068) epolar
      end if
      if (econst /= 0.0) write(7,9078) epot-econst
      if ( dvdl /= 0.d0) write(7,9100) dvdl
      if (induced > 0.and.indmeth < 3) write(7,9190)diprms,dipiter
      if (induced > 0.and.indmeth == 3) write(7,9191)diprms, &
            dipole_temp
   end if

#ifdef RISMSANDER
   if(rismprm%irism==1 .and. rismprm%write_thermo==1)then
      if(rism_calc_type(nstep) == RISM_FULL)&
           call rism_thermo_print(.false.,transfer(ene%pot,pot_array))
   end if
#endif /*RISMSANDER*/

   9018 format (/ /,3x,'NSTEP',7x,'ENERGY',10x,'RMS',12x,'GMAX',9x, &
         'NAME',4x,'NUMBER')
   9028 format(1x,i6,2x,3(2x,1pe13.4),5x,a4,2x,i7,/)
   9038 format (1x,'BOND    = ',f13.4,2x,'ANGLE   = ',f13.4,2x, &
         'DIHED      = ',f13.4)
   9039 format (1x, 'UB      = ', f13.4, 2x, 'IMP     = ', f13.4, 2x, &
            'CMAP       = ', f13.4)
   9048 format (1x,'VDWAALS = ',f13.4,2x,'EEL     = ',f13.4,2x, &
         'HBOND      = ',f13.4)
   9049 format (1x,'VDWAALS = ',f13.4,2x,'EEL     = ',f13.4,2x, &
         'EGB        = ',f13.4)
   9050 format (1x,'VDWAALS = ',f13.4,2x,'EEL     = ',f13.4,2x, &
         'EPB        = ',f13.4)

#ifdef RISMSANDER
   9051 format (1x,'VDWAALS = ',f13.4,2x,'EEL     = ',f13.4,2x, &
         'ERISM      = ',f13.4)
#endif

   9058 format (1x,'1-4 VDW = ',f13.4,2x,'1-4 EEL = ',f13.4,2x, &
         'RESTRAINT  = ',f13.4)
   9062 format (1x,'EMAP   = ',f14.4,   '  EMSCORE = ',f8.4)
   9066 format (1x,'EPOLAR  = ',f13.4,2x,'CMAP    = ',f13.4)
   9067 format (1x,'CMAP    = ',f13.4)
   9068 format (1x,'EPOLAR  = ',f13.4)
#ifdef APBS
   9069 format (1x,'ENPOLAR = ',f13.4)
#endif /* APBS */
   9074 format (1x,'ECAVITY = ',f13.4,2x,'EDISPER = ',f13.4)
   9077 format (1x,'ESURF   = ',f13.4)
   9078 format (1x,'EAMBER  = ',f13.4)
   9080 format (1x,'PM3ESCF =',f14.4)
   9081 format (1x,'AM1ESCF =',f14.4)
   9981 format (1x,'AM1DESCF =',f14.4)
   9082 format (1x,'MNDOESCF=',f14.4)
   9982 format (1x,'MNDODESCF=',f14.4)
   9083 format (1x,'PDDGPM3-ESCF=',f14.4)
   9084 format (1x,'PDDGMNDO-ESCF=',f14.4)
   9085 format (1x,'PM3CARB1-ESCF=',f14.4)
   9086 format (1x,'DFTBESCF=',f14.4)
   9087 format (1x,'RM1ESCF =',f14.4)
   9088 format (1x,'PDDGPM3_08-ESCF=',f14.4)
   9089 format (1x,'PM6ESCF =',f14.4)
   9090 format (1x,'PM3MMXESCF =',f14.4)
   9091 format (1x,'PM3ZNBESCF =',f14.4)
   9092 format (1x,'EXTERNESCF =',f14.4)
   9093 format (1x,'PM3MAISESCF =',f14.4)
   9099 format (1x,'ECRG    =',f14.4)
   9096 format (1x,'PM3MMX2ESCF =',f14.4)
   9190 format(1x,'Dipole convergence: rms = ',e10.3,' iters = ',f6.2)
   9191 format(1x,'Dipole convergence: rms = ',e10.3, &
         ' temperature = ',f6.2)
   9100 format (1x,'DV/DL  = ',f14.4)
   9200 format (1x,'ESEBOMD =',f14.4)
#ifdef PUPIL_SUPPORT
   9900 format (1x,'PUPESCF =',f14.4)
#endif
   return
end subroutine printe
