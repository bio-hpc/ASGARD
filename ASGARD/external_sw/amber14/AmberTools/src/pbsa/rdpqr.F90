#include "copyright.h"
#include "../include/dprec.fh"

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ A mimic of the lite version of Amber rdparm1 for use with pqr reading
subroutine rdpqr1(nf)
   
   use memory_module
   use parms, only : numbnd, numang, nptra, nphb, nimprp, nttyp
   use pbsa_lib, only : get_num_tokens, get_token

   implicit none
   
#  include "../include/md.h"
#  include "box.h"

   integer nf
   integer nhparm
   integer numextra
   integer mbona,mtheta,mphia
   
   ! local variables

   character(len=1024) :: line
   character(len=1024) :: word
   integer :: n = 0
   integer :: i = 0
   integer :: outcome = 0

   natom = 0
   ntypes = 0
   nbonh = 0
   mbona = 0
   ntheth = 0
   mtheta = 0
   nphih = 0
   mphia = 0
   nhparm = 0
   nparm = 0
   nnb = 0
   nres = 0
   nbona = 0
   ntheta = 0
   nphia = 0
   numbnd = 0
   numang = 0
   nptra = 0
   natyp = 0
   nphb = 0
   ifpert = 0
   nbper = 0
   ngper = 0
   ndper = 0
   ifbox = 0
   nmxrs = 0
   ifcap = 0
   numextra = 0
   ncopy = 0

   i = 0 ! atom sequence number counter
   do

      ! read one line into buffer

      read(nf, '(a)', iostat=outcome) line

      ! flawless reading 

      if ( outcome == 0 ) then

         ! Find out how many words are on this field

         call get_num_tokens(line, n)

         ! parsing the line if it is an ATOM and HETATM line

         call get_token(line, 1, word)

         if ( trim(word) .eq. 'ATOM' .or. trim(word) .eq. 'HETATM' ) then
            i = i + 1
         end if

      ! something is wrong

      else if ( outcome > 0 ) then

         write(6,'(a)') 'PBSA Bomb In RDPQR1: Cannot read line'
         write(6,'(a)') trim(line)
         call mexit(6, 1)

      ! end of file reached

      else if ( outcome < 0 ) then

         exit

      end if

   end do ! go to the next line

   ! successfully finishing reading the file

   natom = i

   write(6,8118) natom,ntypes,nbonh,mbona,ntheth,mtheta,nphih,mphia, &
         nhparm,nparm,nnb,nres,nbona,ntheta,nphia,numbnd, &
         numang,nptra,natyp,nphb,ifbox,nmxrs,ifcap,numextra,ncopy
   8118 format(' PQR INPUT', &
         /' NATOM  = ',i7,' NTYPES = ',i7,' NBONH = ',i7,' MBONA  = ',i7, &
         /' NTHETH = ',i7,' MTHETA = ',i7,' NPHIH = ',i7,' MPHIA  = ',i7, &
         /' NHPARM = ',i7,' NPARM  = ',i7,' NNB   = ',i7,' NRES   = ',i7, &
         /' NBONA  = ',i7,' NTHETA = ',i7,' NPHIA = ',i7,' NUMBND = ',i7, &
         /' NUMANG = ',i7,' NPTRA  = ',i7,' NATYP = ',i7,' NPHB   = ',i7, &
         /' IFBOX  = ',i7,' NMXRS  = ',i7,' IFCAP = ',i7,' NEXTRA = ',i7 &
         ,/' NCOPY  = ',i7/)
   
   nttyp = ntypes*(ntypes+1)/2
   
   return 

end subroutine rdpqr1 

!=======================================================================

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ A simple pqr file reader
subroutine rdpqr2(nf, natom, igraph)
!
! PBSA reads a "free-formatted" PQR file with all fields delimited by space.
! It does not require any strict format as in a standard PDB file. This more
! liberal style is to accomodate pqr files of different origins.
! PBSA reads data on a per-line basis using the following format:
! TAG AtomNumber AtomName ResidueName ChainID ResidueNumber X Y Z Charge Radius
!
! TAG: A string specifying either ATOM or HETATM. Lines with other strings are ignored.
! AtomNumber: The sequence no of the atom, which is reset to start from 1.
! AtomName: The atom name.
! ResidueName: The residue name.
! ChainID: The chain ID of the atom, optional, which is ignored.
! ResidueNumber: The sequence no. of the residue, which is ignored.
! X Y Z: The floating numbers representing the atomic coordiantes (in Angstrom).
! Charge: A float number providing the atomic charge (in electron).
! Radius: A float number providing the atomic radius (in Angstrom).
!
! WARNING: It is apparently very very hard to know whether the charge and
! radius fields are swapped as in the Delphi generated pqr file.
! Here we assume the plain P.Q.R. order in the format.
!
   use pbsa_lib, only : get_num_tokens, get_token
   use poisson_boltzmann, only : acrd, acrg
   use solvent_accessibility, only : radi

   implicit none

   ! passed variables

   integer nf
   integer natom
   character(len=4) igraph(natom)
   
   ! local variables

   logical :: chainopt = .false.
   character(len=1024) :: line
   character(len=1024) :: word
   integer :: n = 0
   integer :: i = 0
   integer :: outcome = 0
   character(len=6) :: tag
   character(len=4) resid
   character(len=4) chainid
   integer :: atmnum
   integer :: resnum

   write(6,100)
   100 format(/80(1h-)/,'   3.  ATOMIC COORDINATES, CHARGES, AND RADII FROM PQR INPUT',/80(1h-)/)

   i = 0 ! atom sequence number counter
   do 

      ! read one line into buffer

      read(nf, '(a)', iostat=outcome) line

      ! flawless reading 

      if ( outcome == 0 ) then

         call get_num_tokens(line, n)

         ! Get the first word

         call get_token(line, 1, word)

         if ( trim(word) .eq. 'ATOM' .or. trim(word) .eq. 'HETATM' ) then

            i = i + 1

            ! first check whether we have the correct number of fields

            if ( n .eq. 10 ) then
               chainopt = .false.
            else if ( n .eq. 11 ) then
               chainopt = .true.
            else
               write(6,'(2a,i5)') 'The no. of fields on this line is inconsistent with ', &
                               'the pqr format with or without chain ID.', n
               write(6,'(a)') trim(line)
               call mexit(6, 1)
            end if

            ! next read the line based on whether there is chain id or not
            !
            ! Here is the order of fields on each atom/hetatm entry:
            ! TAG: A string specifying either ATOM or HETATM. Other strings will be ignored.
            ! AtomNumber: The sequence no of the atom, which will be ignored and reset to start from 1.
            ! AtomName: The atom name.
            ! ResidueName: The residue name.
            ! ChainID: The chain ID of the atom, optional, which will be ignored.
            ! ResidueNumber: The sequence no. of the residue, which will be ignored and reset to setart from 1.
            ! X Y Z: The floating numbers representing the atomic coordiantes (in Angstrom).
            ! Charge: A float number providing the atomic charge (in electrons).
            ! Radius: A float number providing the atomic radius (in Angstrom).
            !
            ! WARNING: It is apparently very very hard to know the charge and
            ! radius fields are swapped as in the Delphi generated pqr file.
            ! Here we assume the plain P.Q.R. order in the format.

            if ( chainopt ) read(line, *) tag, atmnum, igraph(i), resid, chainid, resnum, &
                                          acrd(1,i), acrd(2,i), acrd(3,i), acrg(i), radi(i) 
            if ( .not. chainopt ) read(line, *) tag, atmnum, igraph(i), resid, resnum, &
                                          acrd(1,i), acrd(2,i), acrd(3,i), acrg(i), radi(i) 

            write(6,'(a6,i6,a8,a8,i6,f10.3,f10.3,f10.3,f15.6,f15.6)')&
                tag(1:6), i, igraph(i)(1:4), resid(1:4), resnum,&
                acrd(1,i), acrd(2,i), acrd(3,i),acrg(i), radi(i)
         end if

      ! something is wrong

      else if ( outcome > 0 ) then

         write(6,'(a)') 'PBSA Bomb In RDPQR2: Cannot read line'
         write(6,'(a)') trim(line)
         call mexit(6, 1)

      ! end of file reached

      else if ( outcome < 0 ) then

         exit

      end if

   end do ! go to the next line

   write(6,'(a)') 
   write(6,'(a,2i6)') ' The total number of atom entries read is ', i, natom
   write(6,'(a)') ' There is no internal checking imposed.' 
   write(6,'(a)') ' It is assumed that the last two columns are in the order of charges and radii '
   write(6,'(a)') 

   return

end subroutine rdpqr2 
