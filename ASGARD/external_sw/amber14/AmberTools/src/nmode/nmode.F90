#include "assert.h"

!************************************************************************
!                              AMBER                                   **
!                                                                      **
!                  Copyright (c) 1986, 1991, 1995,1997                 **
!             Regents of the University of California                  **
!                       All Rights Reserved.                           **
!                                                                      **
!  This software provided pursuant to a license agreement containing   **
!  restrictions on its disclosure, duplication, and use. This software **
!  contains confidential and proprietary information, and may not be   **
!  extracted or distributed, in whole or in part, for any purpose      **
!  whatsoever, without the express written permission of the authors.  **
!  This notice, and the associated author list, must be attached to    **
!  all copies, or extracts, of this software. Any additional           **
!  restrictions set forth in the license agreement also apply to this  **
!  software.                                                           **
!************************************************************************


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of program nmode here]
program nmode

   use decomp

   implicit double precision (a-h, o-z)
   
   double precision,  dimension(:), allocatable :: x
   integer, dimension(:), allocatable :: ix
   character(len=4), dimension(:), allocatable :: ih

   character(len=1) jobz,uplo
   common /timer/ timtot, tlast, elapst
#  include "files.h"
   
   !     ------------------------------------------------------------------
   
#  include "pointer.h"
#  include "inpdat.h"
#  include "epot.h"
   
   !     ----- Set parameters to be exported via commons
   
   character(len=80) fmt
   character(len=80) fmtin,ifmt,afmt,rfmt,type
   ifmt = '(12I6)'
   afmt = '(20A4)'
   rfmt = '(5E16.8)'
   nf = 51
   
   maxdup = 400
   timtot = 0.0d0
   tlast = 0.0d0
   elapst = 0.0d0
   
   !     ----- Initialize for timing
   
   call timit (timtot, tlast, elapst)
   
   !     ----- Read control data for the job
   
   call nmdfil
   call amopen(5, nmdin, 'O', 'F', 'R')
   if(nmdout /= 'screen') then
      call amopen(6, nmdout, owrite, 'F', 'W')
   end if
   call rdinp
   cut = cut**2
   
   !     ----- Read preliminary info from PARM file
   
   call amopen(51, parm, 'O', 'F', 'R')
   fmtin = afmt
   type = 'TITLE'
   call nxtsec(nf,  6,  0,fmtin,  type,  fmt,  iok)
   read (51,fmt) (ititl(i), i = 1, 20)
   
   fmtin = ifmt
   type = 'POINTERS'
   call nxtsec(nf,  6,  0,fmtin,  type,  fmt,  iok)
   read (51,fmt)   natom, ntypes, nbonh,  mbona,  ntheth, mtheta,      &! 5
         nphih, mphia,  nhparm, nparm,  nnb,    nres,        &! 10
         nbona, ntheta, nphia,  numbnd, numang, nptra,       &! 15
         natyp, nphb, ifpert,idum,idum,idum,idum,idum, &
         idum,idum,idum,idum,idum
   if (ifpert /= 0) then
      write(6,*) ' Error: *** PARM TOPOLOGY CONFIGURED FOR GIBBS'
      call mexit(6,1)
   end if
   nr3 = 3 * natom

   !     kludge for now: set ns3=nr3 (atoms in belly = total), since we
   !     won't read the belly until later:
   ns3 = nr3
   
   write (6,3) (ititl(i),i=1,20)
   3 format (/,3x,' PARM file has the title:',/,5x,20a4)
   
   !     ----- Allocate memory for permanent arrays
   
   call alloc(maxdup, memusd_x, memusd_i, memusd_h)
   
   !     --- dynamic memory allocation:
   
   allocate( x(1:memusd_x), stat = ier )
   REQUIRE( ier == 0 )
   
   allocate( ix(1:memusd_i), stat = ier )
   REQUIRE( ier == 0 )
   
   allocate( ih(1:memusd_h), stat = ier )
   REQUIRE( ier == 0 )
   
   !     ----- Read rest of the parm file
   
   call rdparm (x, ix, ih, x(mchrg), maxdup)
   close (51)
   
   !     ----- Divide into groups if necessary and set group parameters
   
   call setgr (natom,      natc,        nres,        ngroup, &
         ibelly,     natsys,      ix(migres),  ix(mgroup), &
         ix(mipres), ih(mlbres),  ih(mgraph),  ih(msymbl), &
         ih(mitree), icons,       ix(migrp2),  x(mwref), &
         x(mxref) )
   ns3 = natsys*3
   
   !     ----- Read coordinates and minimized reference coordinates
   
   if (ntx == 0) then
      call amopen(53, inpcrd, 'O', 'U', 'R')
   else
      call amopen(53, inpcrd, 'O', 'F', 'R')
   end if
   call getcor(natom,ntx,x(mx))
   if (icons /= 0) then
      if (ntx == 0) then
         call amopen(54, refc, 'O', 'U', 'R')
      else
         call amopen(54, refc, 'O', 'F', 'R')
      end if
      call getref (natom,ntx,x(mxref))
   end if
   
   !     ----- Delete bonds, etc from nonactive part AND set xbel array
   
   call setvar (x, ix, ix(mgroup), ix(mnbel))
   
   call timit (timtot, tlast, elapst)
   tload = elapst
   
   !     ----- Compute second derivative matrix
   
   if (iprr == 1) then
      ipair = 1
   else if (iprw == 1) then
      ipair = 3
   else
      ipair = 2
   end if
   ndrv = 2
   if(ntrun == 3) ndrv=1
   
   call forces (x, ix, ih, ipair, ndrv)
   
   ipair = 4
   
   !     ----- Print energies for reference
   
   call eprnt(0,natom,x(mf))
   
   !     ----- check for job type
   
   if (ntrun == 1) then
      !                             ---get vibrational frequencies
      
      !       ----- first, check to see if the gradient is lower
      !             than drms:
      
      gmax = grdmax(x(mf),rms,3*natom)
      if (rms > drms) then
         write(6,*) 'Root-mean-square gradient of input coords is ', &
               rms
         write(6,*) 'This is greater than the requested maximum:  ', &
               drms
         call mexit(6,1)
      end if
      
      !       ----- Make mass-weighted h by multiplying double
      !       ----- derivatives by appropriate mass factors
      
      call mweight(x(mh),x(mf),dummy,0,ns3,x(mamass),dummy,dummy)
      
      !       ----- if the "level" flag is set, call routine to shift
      !             translations and rotations to higher frequencies
      
      if (ilevel /= 0) call level(x(mh),x(mx),x(mamass),ns3,natsys)
      
      !       ----- Diagonalize h
      
      call timit (timtot, tlast, elapst)
      if( ismem == 0 ) then    ! use the dsyev LAPACK routine
         
         uplo = 'U'
         if (nvect > 0) then
            jobz = 'V'
         else
            jobz = 'N'
         end if
         lwork = 9*ns3
         call fillm (ns3,x(mh),x(mcvec))
         call dsyev(jobz, uplo, ns3, x(mcvec), ns3, x(mcval), &
               x(mcscr),lwork,info)
         if (info /= 0) write(6,*) 'dsyev returns info: ',info
         write(6,'(a10,i8,a22,f8.0)') &
               '| lwork = ',lwork,'; best value would be ',x(mcscr)
         
         !       else   ! save memory with EISPACK
         
         !         call giveis (ns3, nvect, ns3, x(mh),
         !    $       x(mcscr), x(mcval), x(mcvec), ierr)
         !         if (ierr.ne.0) then
         !           write(6,*) 'error: giveis returns ',ierr
         !           call mexit(6,1)
         !         endif
         
      else   ! save memory with dspevd if no eigenvectors are needed
         
         uplo = 'U'
         jobz = 'N'
         lwork = 5*ns3
         lworki = 4*ns3
         call dspevd(jobz, uplo, ns3, x(mh), x(mcval), dummy, &
               1,x(mcscr),lwork, x(mcscr+4*ns3),4*ns3, info)
         if (info /= 0) write(6,*) 'dspevd returns info: ',info
         write(6,'(a10,i8,a22,f8.0)') &
               '| lwork = ',lwork,'; best value would be ',x(mcscr)
         call writlw(lworki,x(mcscr+4*ns3))
         
      end if  ! ( ismem == 0 )
      
      call timit (timtot, tlast, elapst)
      tdiag = elapst
      write(6,'(a,f10.2)') '| Time for diagonalization: ',tdiag
      
      !       ----- Write out the eigenvalues, eigenvectors
      
      if (nvect > 0 .and. idecomp == 0) then
         call mweight(dummy,dummy,dummy,2,ns3,x(mamass),x(mcvec),nvect)
         if (ivform == 0) then
            call amopen(61, vecs, owrite, 'U', 'W')
         else
            call amopen(61, vecs, owrite, 'F', 'W')
         end if
      end if
      call eigout(ns3,nvect,x(mcval),x(mcvec),x(mx),x(mamass),ivform)
      if (ns3 /= natom*3) then
         write(6,*) 'Thermo analysis not supported for belly calc.'
      else
         nvecs =  3*natom - 6
         if(idecomp == 1) then
            call allocate_intlog_decomp(natom)
            call allocate_real_decomp(natom,nres)
            call build_sweight(natom,nvecs,x(mcvec),x(mamass))
            call build_ires(natom,nres,ix(mipres))
            call build_isside(natom,ih(mgraph),ih(msymbl))
         end if
         call thermo(natom,nvecs,ilevel,idecomp,x(mx),x(mamass),x(mcval), &
               x(mh),x(mh+ns3),x(mh+2*ns3), &
               x(mh+3*ns3),t,patm)
         if(idecomp == 1) then
            call checkdec(nres)
            call printdec(nres)
            call deallocate_intlog_decomp()
            call deallocate_real_decomp()
         end if
      end if
      
   else if (ntrun == 2) then
      
      !                             ---find transition state
      
      call tstatf(x,ix,ih,x(mf),x(mh),x(mdd),x(mvect))
      
   else if (ntrun == 3) then
      
      !                             ---conjugate gradient minimization
      
      acc = drms*drms*dble(ns3)
      ndrv = 1
      call zxcgr()
      call putvar(natom,ix(mgroup),x(mxbel),x(mx))
      call savec(natom,x(mx),1)
      write(6,6)
      call eprnt(ifc,natom,x(mf))
      
   else if (ntrun == 4) then
      
      !                             ---modified Newton-Raphson minimization
      
      call nrch(x,ix,ih,x(mh),x(mf),x(mxdir))
      call putvar(natom,ix(mgroup),x(mxbel),x(mx))
      call savec(natom,x(mx),2)
      write(6,6)
      call eprnt(ifc,natom,x(mf))
      
   else if (ntrun == 5) then
      
      !                            ---Langevin normal modes
      
      !       ----- first, check to see if the gradient is lower
      !             than drms:
      
      gmax = grdmax(x(mf),rms,3*natom)
      if (rms > drms) then
         write(6,*) 'Root-mean-square gradient of input coords is ', &
               rms
         write(6,*) 'This is greater than the requested maximum:  ', &
               drms
         call mexit(6,1)
      end if
      call mweight(x(mh),x(mf),x(1),0,nr3,x(mamass),dummy,dummy)
      
      call lmode (natom, x(ma), x(mwr), x(mwi), x(mz), &
            x(mfv1), ix(miv1), x(mh), x(mamass), x(mx), &
            x(mhrad), eta, x(mgam), x(mwinv), nvect, &
            ioseen, hrmax,ih(mgraph))
      
   end if  ! (ntrun == 1)
   
   !     ----- Finish
   
   call timit (timtot, tlast, elapst)
   write (6,5) elapst

   ! deallocate memory:
   deallocate( x, ix, ih, stat = ier )
   REQUIRE( ier == 0 )

   call mexit(6, 0)
   5 format (/ /,'|',5x, 'cpu time = ', f20.2, ' seconds')
   6 format (/ /,30x,'F i n a l   R e s u l t s :'/ /)
end program nmode 

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine fillm here]
subroutine fillm(n,h,a)
   
   !   ---- fill upper triangle of matrix a from aray h:
   
   double precision h(*),a(n,n)
   k = 0
   do 20 i=1,n
      do 10 j=1,i
         k = k + 1
         a(j,i) = h(k)
      10 continue
   20 continue
   return
end subroutine fillm 

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine writlw here]
subroutine writlw(lworki,iw)
   integer lworki,iw
   write(6,'(a10,i8,a22,i8)') &
         '| lworki= ',lworki,'; best value would be ',iw
   return
end subroutine writlw 
