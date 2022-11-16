program nucgen
   
   !************************************************************************
   !                              AMBER                                   **
   !                                                                      **
   !                  Copyright (c) 1986, 1991, 1995                      **
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
   
   ! NUCGEN
   ! copyright @1985,1989 Regents of the University of California
   ! Rev A revision: George Seibel, 1989.
   ! AUTHORS: U.C.Singh, N.Pattabiraman, and S.N.Rao, 1985.
   ! Director: P.A. Kollman
   !
   ! ----- PROGRAM TO GENERATE DOUBLE HELICAL DNA OF DIFFERENT
   !       HELICAL CONFORMATION -----
   
   use nucgen_files_module
   use nucgen_iofile_module
   
   implicit none
   
   integer :: k, nmol, nres, nresm
   character(len=4) ititl(20), kend, iblank, ilbmol
   character(len=4) lbres(1000), lbdum(20)
   data kend,iblank/'END ','    '/
   
   ! --- get file names ---
   call ngfil
   
   ! --- open all files ---
   in = 5
   iout = 6
   ioutb = 10
   !subr amopen(lun,fname,fstat,fform,facc)
   ! -- output
   call amopen(iout, ngout, owrite,'F','W')
   ! -- input
   call amopen(in,   ngin,  'O','F','R')
   ! -- dna data
   call amopen(7,    ngdat, 'O','F','R')
   ! -- pdb output
   call amopen(ioutb,pdbout,owrite,'F','W')
   
   ! ----- READ THE RESIDUE INFORMATION UNTIL AN END CARD
   !       IS FOUND -----
   write(iout,9108)
   write(iout,9118)
   
   nmol = 0
   nres = 0
   
   100 continue
   nmol = nmol+1
   read(in,9008) ititl
   if (ititl(1) .eq. kend) goto 200
   
   write(iout,9008) ititl
   read(in,9008) ilbmol
   
   if (nmol .gt. 2) then
      write(iout,9208)
      call mexit(iout, 1)
   end if
   
   nresm = 0
   
   120 continue
   read(in,9018) (lbdum(k),k=1,16)
   if (lbdum(1) .eq. iblank) goto 180
   
   do k=1,16
      if (lbdum(k) .eq. iblank) goto 120
      nres = nres+1
      nresm = nresm+1
      lbres(nres) = lbdum(k)
   end do
   
   goto 120
   
   180 continue
   
   ! ----- OUTPUT THE MOLECULE INFORMATION -----
   write(iout,9128) nmol,nresm
   write(iout,9158) (lbres(k),k=nres-nresm+1,nres)
   goto 100
   
   200 continue
   
   ! ----- NOW CALL THE GENERATION ROUTINE -----
   call gennuc(nres,lbres,ilbmol)
   
   call mexit(0,0)
   
   9008 format(20a4)
   9018 format(16(a4,1x))
   9108 format(/10x,54(1h-),/10x, &
               'Amber 5   NUCGEN                           UCSF 1997', &
               /10x,54(1h-))
   9118 format(/15x,'INPUT MOLECULES INFORMATION',/)
   9128 format(/5x,'MOLECULE ',i2,' CONTAINS ',i4,' RESIDUES:',/)
   9158 format(5x,12(a4,1x))
   9208 format(/5x,'INPUT CONTAINS MORE THAN TWO MOLECULES AND ASSUMED', &
               /5x,'THAT IT IS NOT A STANDARD DNA OR RNA',/)
   
end program nucgen

!======================================================================

subroutine gennuc(nres,lbres,ilbmol)
   
   !************************************************************************
   !                              AMBER                                   **
   !                                                                      **
   !                  Copyright (c) 1986, 1991, 1995                      **
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
   
   ! ----- ROUTINE TO GENERATE STANDARD DNA AND RNA DOUBLE
   !       HELICES. SIX TYPES ARE AVAILABLE:
   !
   !       ARNA    ...   RIGHT HANDED A-RNA (ARNOTT)
   !       APRNA   ...   RIGHT HANDED A PRIME RNA (ARNOTT)
   !       BDNAL   ...   RIGHT HANDED BDNA (LANGRIDGE)
   !       BDNAA   ...   RIGHT HANDED BDNA (ARNOTT)
   !       SDNA    ...   LEFT HANDED BDNA (SASISEKHARAN)
   !       ADNA    ...   RIGHT HANDED A-DNA (ARNOTT)
   
   use nucgen_iofile_module
   
   implicit none
   
   integer, intent(in) :: nres
   character(len=4), intent(in) :: lbres(1000), ilbmol
   
   real, parameter :: pi = 3.14159265358979323846
   real, parameter :: conv = pi/180.0e0
   
   integer :: i, iatom, iacid, inew, irds(9,2), irde(9,2), ir1(9), ir2(9), &
              ires, ipmol(500), indexd(5,15), indexu(5,15), j, k, &
              nbasp, nf, nresh, ichain, nnseq, idoi, idof, ip, &
              ig, idown, iup
   real :: hxrep(6), hxht(6), urot, uhigh, r(60), phi(60), zz(60), &
           rot, high, hxmul, x, y, z, xrad, yyr
   character(len=8) name(6), nama, jtypm, ispec
   character(len=4) tmpnam, atnam, idna, irna, jseq(9), katn(3), &
                    natn(3), kseq(15), namat(60)
   
   data name/'$ARNA'   ,'$APRNA'  ,'$LBDNA'  ,'$ABDNA'  , &
             '$SBDNA'  ,'$ADNA'   /
   data ispec/'$SPECIAL'/
   data idna,irna/'D   ','R   '/
   data jseq/'HB  ','HE  ','POM ','RIBO','ADE ','GUA ','THY ', &
             'CYT ','URA '/
   data katn/'OA  ','OB  ','O1  '/
   data natn/'O1P ','O2P ','O4  '/
   data irds/1,2,3,6,14,26,40,50,0, 1,2,3,6,15,27,0,41,51/
   data irde/1,2,5,13,25,39,49,59,0, 1,2,5,14,26,40,0,50,59/
   data hxrep/ 32.7, 30.0, 36.0, 36.0, 36.0, 32.7/
   data hxht/ 2.81, 3.00, 3.38, 3.38,-3.38, 2.56/
   data kseq/'A5  ','A   ','A3  ','G5  ','G   ','G3  ','T5  ', &
            'T   ','T3  ','C5  ','C   ','C3  ','U5  ','U   ', &
            'U3  '/
   
   ! ----- EVALUATE SOME CONSTANTS -----
   iatom = 0
   nbasp = 0
   iacid = 0
   nf = 7
   
   ! ----- SET THE CORRESPONDENCE BETWEEN LBRES AND
   !       IPMOL -----
   call setseq(nres,lbres,ipmol,jseq,kseq,inew)
   if (inew .eq. 1) &
      write(iout, '(4x, a)') 'New (1994) residue naming convention'
   
   ! ----- READ THE TYPE OF NUCLEIC ACID NEEDED -----
   read(in,9008) jtypm
   
   ! ----- IACID IS 1 FOR DNA OR 2 FOR RNA -----
   if (ilbmol .eq. idna) iacid = 1
   if (ilbmol .eq. irna) iacid = 2
   if (iacid .eq. 0) goto 900
   
   ! ----- TRANSFERRING THE NUMBER OF ATOMS FOR DEOXY OR RIBOSE SUGAR
   !       AND BASES -----
   
   do i=1,9
      ir1(i) = irds(i,iacid)
      ir2(i) = irde(i,iacid)
   end do
   
   if (inew .eq. 0) goto 139
   
   ! --- SET UP NEW STYLE OF NUCLEIC ACIDS RESIDUES NOMENCLATURE 
   do i=1,5
      do j=1,12
         indexd(i,j) = 0
         indexu(i,j) = 0
      end do
   end do
   
   ! ----- FIND NEW ATOM INDICES FOR a) HB, b) HE, c) SUGARS, -----
   !       d) POM and e) BASES 
   
   ! HB
   do j=1,15,3
      indexd(1,j) = ir1(1)
      indexu(1,j) = ir2(1)
   end do
   
   ! HE
   do j=3,15,3
      indexd(5,j) = ir1(2)
      indexu(5,j) = ir2(2)
   end do
   
   ! POM
   do j=2,15,3
      indexd(2,j) = ir1(3)
      indexu(2,j) = ir2(3)
   end do
   
   do j=3,15,3
      indexd(2,j) = ir1(3)
      indexu(2,j) = ir2(3)
   end do
   
   ! SUGAR
   do j=1,15
      indexd(3,j) = ir1(4)
      indexu(3,j) = ir2(4)
   end do
   
   ! BASES
   j = 1
   do k=5,9
      do i=1,3
         indexd(4,j) = ir1(k)
         indexu(4,j) = ir2(k)
         j = j + 1
      end do
   end do
   
   !do i=1,5
   !   write(6,1010) (indexd(i,j),j=1,15) 
   !   write(6,1010) (indexu(i,j),j=1,15) 
   !end do
   
   !1010 format(1x,15i3)
   
   ! ----- RETRIEVE THE UNIT TWIST AND HEIGHT FROM THE DATA -----
   139 continue
   do i=1,6
      if (jtypm .eq. name(i)) then
         if (jtypm .eq. name(1)) write(6,9218)
         if (jtypm .eq. name(2)) write(6,9228)
         if (jtypm .eq. name(3)) write(6,9238)
         if (jtypm .eq. name(4)) write(6,9248)
         if (jtypm .eq. name(5)) write(6,9258)
         if (jtypm .eq. name(6)) write(6,9268)
         urot = hxrep(i)
         uhigh = hxht(i)
         goto 240
      end if
   end do
   
   if (jtypm .eq. ispec) then
      read(in,9018) urot,uhigh
      write(iout,9108) urot,uhigh
   else
      write(iout,*) 'jtypm does not equal ispec: ', jtypm, ispec
      goto 900
   end if
   
   ! ----- NOW THE HELICAL PARAMETERS ARE AVAILABLE -----
   240 continue
   
   ! ----- READING THE MONOMER UNITS' R ,PHI, ZZ COORDINATES -----
   write(iout,'(2x,a)') 'reading monomer parameters'
   if (inew .eq. 0) then
      
      260 continue
      read(nf,9008,end=900) nama
      if (nama .eq. jtypm) goto 270
      goto 260
      
      270 continue
   else
      ! ----- NEW-AMBER ATOM NAMES TO BE USED -------
      272 continue
      !write(iout,*)' searching for ', jtypm
      read(nf,9008,end=900) nama
      if (nama .eq. jtypm) goto 275
      goto 272
      
      275 continue
   end if
   
   do i=1,60
      read(nf,9028)  namat(i),r(i),phi(i),zz(i)
   end do
   
   rot = 0.0
   high = 0.0
   ires = 0
   nresh = nres/2
   
   ! ----- GENERATE THE COORDINATES OF THE TWO CHAINS -----
   hxmul = -1.0E0
   if (inew .eq. 0) then
      
      ! ---- OLD DNA/RNA NOMENCLATURE ------------
      do ichain=1,2
         if (ichain .eq. 2) hxmul = 1.0E0
         do i=1,nresh
            if (ichain .eq. 1) nnseq = ipmol(i)
            if (ichain .eq. 2) nnseq = ipmol(i+nresh)
            ires = ires+1
            
            ! ----- IF THE RESIDUE IS A POM,HB OR HE SKIP THE SUGAR -----
            if (nnseq .ge. 4) then
               
               ! ----- GENERATE THE SUGAR BACKBONE FIRST -----
               idoi = ir1(4)
               idof = ir2(4)
               
               do ip=idoi,idof
                  xrad = r(ip)
                  yyr = (hxmul*phi(ip)+rot)*conv
                  x = xrad*cos(yyr)
                  y = xrad*sin(yyr)
                  z = hxmul*zz(ip)+high
                  iatom = iatom+1
                  
                  ! --- 3 char names start col. 14, 
                  !     4 char names start col 13. ---
                  atnam = ' '
                  tmpnam = namat(ip)
                  if (tmpnam(4:4) .eq. ' ') then
                     atnam(2:4) = tmpnam
                  else
                     atnam = tmpnam
                  end if
                  write(ioutb,9118) iatom,atnam,jseq(nnseq),ires,x,y,z
               end do
            end if
            
            ! ----- GENERATE THE BASE OR POM ETC COORDINATES -----
            idoi = ir1(nnseq)
            idof = ir2(nnseq)
            
            do ig=idoi,idof
               yyr = (hxmul*phi(ig)+rot)*conv
               xrad = r(ig)
               x = xrad*cos(yyr)
               y = xrad*sin(yyr)
               z = hxmul*zz(ig)+high
               iatom = iatom+1
               
               ! --- 3 char names start col. 14, 
               !     4 char names start col 13. ---
               atnam = ' '
               tmpnam = namat(ig)
               if (tmpnam(4:4) .eq. ' ') then
                  atnam(2:4) = tmpnam
               else
                  atnam = tmpnam
               end if
               
               write(ioutb,9118) iatom,atnam,jseq(nnseq),ires,x,y,z
            end do
            
            ! ----- DO NOT INCREMENT THE HELIX IF HB ETC -----
            if (nnseq .ge. 4) then
               rot = rot+urot
               high = high+uhigh
            end if
            
         ! ----- END OF MOLECULE GENERATION -----
         end do
         
         ! ----- REVERSE THE UNIT HEIGHT AND TWIST IN ORDER TO
         !       DESCEND THE CHAIN -----
         urot = -urot
         rot = rot+urot
         uhigh = -uhigh
         high = high+uhigh
         
      ! ----- END OF CHAIN GENERATION -----
      end do
      
   else
      ! ------- PART USED FOR NEW DNA/RNA NOMENCLATURE ---------
      do ichain=1,2
         if (ichain .eq. 2) hxmul = 1.0E0
         
         do i=1,nresh
            if (ichain .eq. 1) nnseq = ipmol(i)
            if (ichain .eq. 2) nnseq = ipmol(i+nresh)
            ires = ires+1
            
            ! ----- GENERATE NUCLEOSIDE COORDINATES -----
            do k=1,5
               idown = indexd(k,nnseq)
               iup   = indexu(k,nnseq)
               if (idown .ne. 0) then
                  
                  ! --- correct for the last HE atom:
                  if (k .eq. 5.and.i .eq. nresh) then
                     !write(6,1011) k,i
                     !1011 format(1x,'k,i',2i5)
                     high = high + uhigh
                     rot = rot + urot
                  end if
                  
                  do ig=idown,iup
                     yyr = (hxmul*phi(ig)+rot)*conv
                     xrad = r(ig)
                     x = xrad*cos(yyr)
                     y = xrad*sin(yyr)
                     z = hxmul*zz(ig) + high
                     iatom = iatom + 1
                     
                     ! -- terminal hydrogens (HB HE residues in old ff)
                     !    have special goofy names
                     if (i .eq. 1 .and. k .eq. 1 .and. ig .eq. idown) then
                        atnam = ' H5T'
                     else if (i .eq. nresh .and. k .eq. 5 .and. ig .eq. iup) then
                        atnam = ' H3T'
                     else
                        ! --- adjust phosphate & sugar names
                        if (namat(ig) .eq. katn(1)) then
                           namat(ig) = natn(1)
                        else if (namat(ig) .eq. katn(2)) then
                           namat(ig) = natn(2)
                        else if (namat(ig) .eq. katn(3)) then
                           namat(ig) = natn(3)
                        end if
                        
                        ! --- 3 char names start col. 14, 
                        !     4 char names start col 13. ---
                        atnam = ' '
                        tmpnam = namat(ig)
                        if (tmpnam(4:4) .eq. ' ') then
                           atnam(2:4) = tmpnam
                        else
                           atnam = tmpnam
                        end if
                     end if
                     
                     write(ioutb,9118) iatom,atnam,kseq(nnseq),ires,x,y,z
                  end do
               end if
            end do
            
            rot = rot+urot
            high = high+uhigh
            
         ! ----- END OF MOLECULE GENERATION -----
         end do
         
         rot = rot - urot
         high = high - uhigh
         
         ! ----- REVERSE THE UNIT HEIGHT AND TWIST IN ORDER TO
         !       DESCEND THE CHAIN -----
         urot = -urot
         rot = rot+urot
         uhigh = -uhigh
         high = high+uhigh
         
      ! ----- END OF CHAIN GENERATION -----
      end do
   end if
   
   close(unit=nf)
   return
   
   900 write(iout,9408)
   call mexit(iout, 1)
   
   9008 format(5A8)
   9018 format(3F10.4)
   9028 format(A4,3F10.4)
   9108 format(/5X,'SPECIAL HELICAL VALUES', 'urot =',F8.3, &
              ' uhigh =',F8.3,/)
   9118 format('ATOM',2X,I5,1X,A4,1X,A3,2X,I4,4X,3F8.3)
   9218 format(/10X,'GENERATING RIGHT HANDED A-RNA(ARNOTT)',/)
   9228 format(/10X,'GENERATING RIGHT HANDED AP-RNA',/)
   9238 format(/10X,'GENERATING RIGHT HANDED BDNA (LANGRIDGE)',/)
   9248 format(/10X,'GENERATING RIGHT HANDED BDNA (ARNOTT)',/)
   9258 format(/10X,'GENERATING LEFT HANDED BDNA (SASISEKHARAN)',/)
   9268 format(/10X,'GENERATING RIGHT HANDED A-DNA',/)
   9408 format(/10X,'ERROR IN INPUT ... GENNUC')
   
end subroutine gennuc

!======================================================================

subroutine setseq(nres,lbres,ipmol,jseq,kseq,inew)
   
   !************************************************************************
   !                              AMBER                                   **
   !                                                                      **
   !                  Copyright (c) 1986, 1991, 1995                      **
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
   
   use nucgen_iofile_module
   
   implicit none
   
   integer, intent(in) :: nres
   integer, intent(inout) :: ipmol(500)
   integer, intent(out) :: inew
   character(len=4), intent(in) :: lbres(1000), jseq(9), kseq(15)
   
   integer :: i, j, ierr
   
   ! ---- TRY THE OLD NOMENCLATURE STYLE ----
   ierr = 0
   inew = 0 ! BPR 13 Jan 2010: Added to make sure inew has a value
   do i=1,nres
      do j=1,9
         if (lbres(i) .eq. jseq(j)) goto 140
      end do
      
      ierr = 2 
      
      140 continue
      ipmol(i) = j
   end do
   
   if (ierr .eq. 2) then
      
      ! ---- TRY THE NEW NOMENCLATURE STYLE ----
      ierr = 1
      inew = 1
      do i=1,nres
         do j=1,15
            if (lbres(i) .eq. kseq(j)) goto 1400
         end do
         
         ierr = 2 
         write(6,'(a,a)') ' Unknown residue: ', lbres(i)
         
         1400 continue
         ipmol(i) = j
      end do
   end if
   
   if (ierr .lt. 2) return
   
   write(iout,*) ' RESIDUE NAMES ARE NOT UNIFORMLY OLD OR 1994 NUCLEIC ACIDS'
   call mexit(iout, 1)

end subroutine setseq
