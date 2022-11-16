! <compile=optimized>
#include "copyright.h"
#include "../include/assert.fh"
#include "../include/dprec.fh"

!++++++++++++++++++++++++++++++++++++
!This module contains all variables 
!needed for BAR type FEP calculations
!++++++++++++++++++++++++++++++++++++

  module mbar

    implicit none

#include "parallel.h"
#ifdef MPI
    include 'mpif.h'
    integer,save :: ist(MPI_STATUS_SIZE), partner
#endif

    ! collect FEP energy data for BAR analysis, every bar_intervall steps
    integer,save :: ifmbar, bar_intervall, bar_states, bar_values
    integer :: ier_alloc, ierr, bar_i
    ! minimum, maximum lambda values, increment between those
    _REAL_,save :: bar_l_min, bar_l_max, bar_l_incr         
    ! collects the FEP energy contributions each turn
    _REAL_, dimension(:), allocatable,save :: bar_cont, bar_cont_partner 
    ! collects all final FEP energies     
    ! _REAL_, dimension(:,:), allocatable,save :: bar_ene
    ! contains the lambda values to loop over
    _REAL_, dimension(:), allocatable,save :: bar_lambda
    ! collect mbar contributions on this turn?
    logical, save :: do_mbar
  
#ifdef MPI
  contains
    ! this module contains no subroutines in serial compiles
    ! MBAR functionality is only invoked, setup and used in parallel
    ! so there is no real use for any non-MPI code here
  
!===========================================
    subroutine setup_mbar(nstlim)

      integer, intent(in) :: nstlim
      integer :: i

      if (sanderrank==0) then
         call mpi_barrier( commmaster, ierr)
         partner = ieor(masterrank,1)
      end if

      do_mbar = .false.

      bar_states = 1 + nint( (bar_l_max-bar_l_min) / bar_l_incr )
      bar_values = nstlim / bar_intervall

      allocate (bar_lambda(bar_states),stat=ier_alloc)
      if (ier_alloc /= 0) call sander_bomb('allocate bar_lambda[mbar.f]','cant allocate','')
!      allocate (bar_ene(bar_states,bar_values),stat=ier_alloc)
!      if (ier_alloc /= 0) call sander_bomb('allocate bar_ene[mbar.f]','cant allocate','')
      allocate (bar_cont(bar_states),stat=ier_alloc)
      if (ier_alloc /= 0) call sander_bomb('allocate bar_cont[mbar.f]','cant allocate','')
      allocate (bar_cont_partner(bar_states),stat=ier_alloc)
      if (ier_alloc /= 0) call sander_bomb('allocate bar_cont_partner[mbar.f]','cant allocate','')

      do i=1, bar_states
         bar_lambda(i)=bar_l_min
         bar_l_min = bar_l_min + bar_l_incr
      end do

      if (sanderrank==0) then
         write (6,*)
         write (6,'(a)')           '    MBAR - lambda values considered:'
         write (6,'(a,i4,a,20(1x,f6.4))') '    ', bar_states, ' total: ', bar_lambda(1:bar_states)
         write (6,'(a,i6,a)')      '    Extra energies will be computed ', bar_values, ' times.'
      end if

      bar_cont(1:bar_states) = 0.0d0
      
    end subroutine setup_mbar

!===========================================

    subroutine bar_collect_cont()

      _REAL_, dimension(bar_states) :: bar_cont_tot

      bar_cont_tot(1:bar_states) = 0.0d0

      ! collect contributions from the nodes
      call mpi_reduce(bar_cont, bar_cont_tot, bar_states, MPI_DOUBLE_PRECISION, &
           MPI_SUM, 0, commsander, ierr)

      bar_cont(1:bar_states)=bar_cont_tot(1:bar_states)

    end subroutine bar_collect_cont

!===========================================

    subroutine calc_mbar_energies(etot, etot_partner)

      integer :: i
      _REAL_, intent(in) :: etot, etot_partner
      _REAL_, dimension(bar_states) :: bar_cont_partner
      _REAL_ :: energy

      if (sanderrank==0) then
         ! exchange contributions
         call mpi_sendrecv( bar_cont, bar_states, MPI_DOUBLE_PRECISION, partner, 5, &
              bar_cont_partner, bar_states, MPI_DOUBLE_PRECISION, partner, 5, &
              commmaster, ist, ierr )

         write (6,'(/a)') 'MBAR Energy analysis:'
         do i=1, bar_states
            if (masterrank==0) then
               energy = (1.0d0 - bar_lambda(i)) * (etot + bar_cont(i) ) + &
                    bar_lambda(i) * (etot_partner + bar_cont_partner(i) )
            else
               energy = (1.0d0 - bar_lambda(i)) * (etot_partner + bar_cont_partner(i) ) + &
                    bar_lambda(i) * (etot + bar_cont(i) )
            end if
            write (6,'(a,f6.4,a,f12.4)') 'Energy at ', bar_lambda(i), ' = ', energy
         end do
         bar_cont(1:bar_states) = 0.0d0
      end if

    end subroutine calc_mbar_energies

!===========================================

    subroutine cleanup_mbar()

      if(allocated(bar_lambda)) deallocate (bar_lambda,stat=ier_alloc)
      if (ier_alloc /= 0) call sander_bomb('cleanup_mbar[mbar.f]','cant deallocate bar_lambda','')
!      if(allocated(bar_ene)) deallocate (bar_ene,stat=ier_alloc)
!      if (ier_alloc /= 0) call sander_bomb('cleanup_mbar[mbar.f]','cant deallocate bar_ene','')
      if(allocated(bar_cont)) deallocate (bar_cont,stat=ier_alloc)
      if (ier_alloc /= 0) call sander_bomb('cleanup_mbar[mbar.f]','cant deallocate bar_cont','')
      if(allocated(bar_cont_partner)) deallocate (bar_cont_partner,stat=ier_alloc)
      if (ier_alloc /= 0) call sander_bomb('cleanup_mbar[mbar.f]','cant deallocate bar_cont_partner','')
      
    end subroutine cleanup_mbar

!===========================================
#endif 
  end module mbar
