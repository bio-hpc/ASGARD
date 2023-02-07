! <compile=optimized>
#include "assert.h"

module decomp

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
! --- Module for vibrational entropy decomposition
!
!     Hannes Kopitz
!     18/09/2009
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

integer, public, dimension(:), allocatable  :: ires
integer, public                             :: natm, nrs

integer, private, parameter                 :: ndectype = 6

integer, private, parameter                 :: &

            ! DECOMP contributions for sidechains
            istranl =  1, &  ! DECOMP S_TRN
            isrotnl =  2, &  ! DECOMP S_ROT
            isvibrl =  3, &  ! DECOMP S_VIB
            ! DECOMP contributions for backbones
            ibtranl =  4, &  ! DECOMP S_TRN
            ibrotnl =  5, &  ! DECOMP S_ROT
            ibvibrl =  6     ! DECOMP S_VIB

integer, private, dimension(ndectype)       :: ndecind
!  index into the dec array

double precision,  private, dimension(:), allocatable   :: dec, sumpermode
!  stores decomposition values
double precision,  public, dimension(:,:), allocatable :: sweight
logical, public, dimension(:), allocatable  :: isside

contains

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ allocates scratch space for integer arrays of module decomp
subroutine allocate_intlog_decomp(natom)

   implicit none
   integer, intent(in) :: natom
   integer ier

   allocate(ires(natom), stat = ier)
   REQUIRE(ier == 0)
   ires(:) = 0
   allocate(isside(natom), stat = ier)
   REQUIRE(ier == 0)
   isside(:) = .true.

   return

end subroutine allocate_intlog_decomp


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ allocates scratch space for real array of module decomp
subroutine allocate_real_decomp(natom,nres)

   implicit none
   integer, intent(in) :: natom, nres
   integer i, j, ncnt, ier

   allocate(dec(ndectype * nres), stat = ier)
   REQUIRE(ier == 0)
   allocate(sweight(3*natom-5,natom), stat = ier)
   REQUIRE(ier == 0)
   allocate(sumpermode(3*natom-5), stat = ier)
   REQUIRE(ier == 0)

   ! Init ndecind and dec arrays
   ncnt = 0
   do i=1,ndectype
      ndecind(i) = ncnt
      do j=1,nres
         ncnt = ncnt + 1
         dec(ncnt) = 0.0d0
      end do
   end do

   return

end subroutine allocate_real_decomp


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ deallocates scratch space for integer arrays of module decomp
subroutine deallocate_intlog_decomp( )

   implicit none
   integer ier

   if(allocated(ires)) then
      deallocate(ires, stat = ier)
      REQUIRE(ier == 0)
   else
      ASSERT(.false.)  ! cannot deallocate un-allocated array
   end if
   if(allocated(isside)) then
      deallocate(isside, stat = ier)
      REQUIRE(ier == 0)
   else
      ASSERT(.false.)  ! cannot deallocate un-allocated array
   end if

   return

end subroutine deallocate_intlog_decomp


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ deallocates scratch space for real array of module decomp
subroutine deallocate_real_decomp( )

   implicit none
   integer ier

   if(allocated(dec)) then
      deallocate(dec, stat = ier)
      REQUIRE(ier == 0)
   else
      ASSERT(.false.)  ! cannot deallocate un-allocated array
   end if
   if(allocated(sweight)) then
      deallocate(sweight, stat = ier)
      REQUIRE(ier == 0)
   else
      ASSERT(.false.)  ! cannot deallocate un-allocated array
   end if
   if(allocated(sumpermode)) then
      deallocate(sumpermode, stat = ier)
      REQUIRE(ier == 0)
   else
      ASSERT(.false.)  ! cannot deallocate un-allocated array
   end if

   return

end subroutine deallocate_real_decomp


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ calculate weighting factor for each atom to each normal mode
subroutine build_sweight(natom,nvecs,vec,amass)
   
   implicit double precision(a-h,o-z)

   dimension vec(natom*3,nvecs), amass(*)

   do i=1,nvecs
      sumpermode(i) = 0.0d0
      do j=1,natom
         k = j*3
         sweight(i,j) = amass(j)*(vec(k-2,i)**2 + vec(k-1,i)**2 + vec(k,i)**2)
         sumpermode(i) = sumpermode(i) + sweight(i,j)
      end do
   end do
   i = nvecs+1
   do j=1,natom
      sweight(i,j) = amass(j)
      sumpermode(i) = sumpermode(i) + sweight(i,j)
   end do

end subroutine build_sweight


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ build array which contains the residue no. of each atom
subroutine build_ires(natom,nres,mres)
   
   implicit double precision(a-h,o-z)

   dimension mres(*)

   do i=1,natom
      do j=1,nres
         if(mres(j) .gt. i) then
            ires(i) = j-1
            goto 90
         end if
      end do
      ires(i) = nres
      90 continue
   end do

end subroutine build_ires


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ build array which contains side-chain information for each atom
subroutine build_isside(natom,ianame,iatype)
   
   implicit double precision(a-h,o-z)

   character(len=4) ianame, iatype
   dimension ianame(*), iatype(*)
   character(len=4) igrap(15),isymb(15)
   data igrap /"C   ","O   ","N   ","H   ","CA  ","HA  ","HA2 ","HA3 ", & ! proteins
               "N   ","H1  ","H2  ","H3  ","HA  ","O   ","OXT "/          ! proteins N/C-terminal
   data isymb /"C   ","O   ","N   ","H   ","CT  ","H1  ","H1  ","H1  ", & ! proteins
               "N3  ","H   ","H   ","H   ","HP  ","O2  ","O2  "/          ! proteins N/C-terminal


   do i=1,natom
      do j=1,15
         if(ianame(i) .eq. igrap(j) .and. iatype(i) .eq. isymb(j)) then
            isside(i) = .false.
         end if
      end do
   end do

end subroutine build_isside


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ accumulates contributions fval into S_TOT/VIB/TRN
subroutine decmode(itype,nivec,fval,natom)
   
   implicit none
   
   double precision fval, val
   integer itype
   integer nivec, natm, natom
   integer nind

   do natm=1,natom
      if(itype == 1) then
         !       --- translational entropy contribution
         if(isside(natm)) then
            nind = ndecind(istranl)
         else
            nind = ndecind(ibtranl)
         end if
      else if(itype == 2) then
         !       --- vibrational entropy contribution
         if(isside(natm)) then
            nind = ndecind(isvibrl)
         else
            nind = ndecind(ibvibrl)
         end if
      else
         write(6,*) 'Wrong input for itype: ',itype
         call mexit(6,1)
      end if  ! (itype == 1)

      ! --- Decompose
   
      ! --- Accumulate results
   
      if(nind < 0) then
         write(6,*) 'NIND wrong in DECMODE.'
         write(6,*) 'ITYPE = ',itype
         write(6,*) nivec
         call mexit(6,1)
      end if

      val = fval*sweight(nivec,natm)/sumpermode(nivec)
      dec(nind + ires(natm)) = dec(nind + ires(natm)) + val
   end do
   
   return
end subroutine decmode


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ sums contents in S_TOT/VIB
subroutine checkdec(ndec)
   
   implicit none

   double precision strn, srot, svib
   integer i, ndec
   integer nsstrn, nsbtrn, nssrot, nsbrot, nssvib, nsbvib

   nsstrn = ndecind(istranl)
   nsbtrn = ndecind(ibtranl)
   nssrot = ndecind(isrotnl)
   nsbrot = ndecind(ibrotnl)
   nssvib = ndecind(isvibrl)
   nsbvib = ndecind(ibvibrl)
                                 
   strn = 0.0d0
   srot = 0.0d0
   svib = 0.0d0

   do i=1,ndec
      strn = strn + dec(nsstrn + i) + dec(nsbtrn + i)
      srot = srot + dec(nssrot + i) + dec(nsbrot + i)
      svib = svib + dec(nssvib + i) + dec(nsbvib + i)
   end do
   
   ! --- Output entropies
   
   write(6,300)
   write(6,350) strn
   write(6,360) srot
   write(6,370) svib
   
   write(6,'(/)')
   
   300 format(/ /20x,'CHECK DECOMP - ENTROPIES',/)
   350 format('S_trn = ',f13.4)
   360 format('S_rot = ',f13.4)
   370 format('S_vib = ',f13.4)
   
   return
end subroutine checkdec 


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ outputs contents in S_TRN/ROT/VIB with respect to residues
subroutine printdec(ndec)
   
   implicit none
   
   double precision  entropy
   dimension entropy(3)

   integer i, j
   integer ndec
   integer nsstrn, nsbtrn, nssrot, nsbrot, nssvib, nsbvib
   
   nsstrn = ndecind(istranl)
   nsbtrn = ndecind(ibtranl)
   nssrot = ndecind(isrotnl)
   nsbrot = ndecind(ibrotnl)
   nssvib = ndecind(isvibrl)
   nsbvib = ndecind(ibvibrl)
                                      
   ! --- Output total entropies
   
   write(6,300)
   write(6,350)
   write(6,355)
   do i=1,ndec
      entropy(:) = 0.0d0
      entropy(1) = dec(nsstrn + i) + dec(nsbtrn + i)
      entropy(2) = dec(nssrot + i) + dec(nsbrot + i)
      entropy(3) = dec(nssvib + i) + dec(nsbvib + i)
      write(6,360) i,(entropy(j),j=1,3)
   end do
   
   write(6,310)
   write(6,350)
   write(6,355)
   do i=1,ndec
      entropy(:) = 0.0d0
      entropy(1) = dec(nsbtrn + i)
      entropy(2) = dec(nsbrot + i)
      entropy(3) = dec(nsbvib + i)
      write(6,370) i,(entropy(j),j=1,3)
   end do
   
   write(6,320)
   write(6,350)
   write(6,355)
   do i=1,ndec
      entropy(:) = 0.0d0
      entropy(1) = dec(nsstrn + i)
      entropy(2) = dec(nssrot + i)
      entropy(3) = dec(nssvib + i)
      write(6,380) i,(entropy(j),j=1,3)
   end do
   
   write(6,'(/)')
   
   300 format(/ /20x,'PRINT DECOMP - TOTAL ENTROPIES',/)
   310 format(/ /20x,'PRINT DECOMP - BACKBONE ENTROPIES',/)
   320 format(/ /20x,'PRINT DECOMP - SIDECHAIN ENTROPIES',/)
   350 format(3x,1x,'resid |translational     |rotational        |vibrational       ')
   355 format('==========',3('==================='))
   360 format('TDC',1x,i6,3(1x,f18.9))
   370 format('BDC',1x,i6,3(1x,f18.9))
   380 format('SDC',1x,i6,3(1x,f18.9))
   
   return
end subroutine printdec 


end module decomp
