!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine rgroup2 here]
subroutine rgroup2(natom,natc,nres,ngrp,ipres,lbres,igraph,isymbl, &
      itree,igroup,weit,xc,konst,belly,nfu)
   
#ifdef DPREC
   implicit double precision (a-h,o-z)
#endif
   
   !     Mods for Rev A by gls:
   !     - cpp selectable precision
   !     - wildcard functionality fixed. (see findp)
   !     - format change to print less garbage
   !     - unreferenced statement labels removed
   !     - changed title write from 12 to 19 A4. (not 20 to avoid wraps)
   !     - changed to read from unit nf -- dac 12/90
   
   logical konst,misc,belly
   
   !     ----- READ IN GROUPS EACH ATOM IS IN -----
   
   !           WILL KEEP READING IN GROUPS OF CARDS UNTIL A BLANK CARD IS
   !           READ
   !           ALL ATOMS IN THIS CARD GROUP WILL BE PUT IN THE SAME GROUP
   !           THE ONLY EXCEPTION IS A RES CARD WHICH STARTS WITH A -
   !           RES NO. THE RESIDUES IN THIS GROUP WILL ALL BE PUT IN INDIV.
   !           GROUPS ANY OTHER GROUPS IN THIS SECTION MUST ALSO START
   !           WITH A  I.E. THEY MUST ALSO BE INDIV. GROUPS
   !           ANY TIME A "FIND" CARD IS READ, ONLY THOSE ATOMS THAT MEET
   !           THE SPECIFICATIONS ON AT LEAST 1 OF THE FIND CARDS
   !           WILL BE INCLUDED IN THE GROUP
   !           THIS OPTION MAY BE ENDED BY READING IN "NFIND",TERMINATING
   !           THE INPUT FOR THE GROUPS, OR BY READING IN ANOTHER "FIND"
   
   common/propf/jgraph(11),jresnm(11),jsymbl(11),jtree(11),isrch
   dimension ititl1(20)
   
   dimension ipres(*),igroup(*),weit(*),igrp(8),jgrp(8)
   dimension lbres(*),igraph(*),isymbl(*),itree(*),xc(*)
   dimension ifld(20),jfld(20),ihol(20),ivar(20),fvar(20)
   
   data ifld/ 2*0, 13*2 , 5*0 /
   data jfld/ 4*1, 16*0 /
   data katn /4hatom /
   data ifind /4hfind /
   data nfind /4hnfin /
   data iiend /4hend /
   data ksear /4hsear/
   
   !     ----- RES CARD LISTS RESIDUE GROUPS IGRP(I) TO JGRP(I) ----
   
   !           IF 1ST VALUE IS NEGATIVE, THEN EVERY RESIDUE FROM IGRP(I)
   !           TO JGRP(I) IS CONSIDERED A SEPARATE GROUP
   !           IF 2ND VALUE   = 0, THEN SUBGROUP CONSISTS OF 1 RESIDUE
   !           ATOM CARD READS IN GROUPS OF ATOMS
   !           IF 2ND ATOM IN EACH PAIR  = 0 , THEN SUBGROUP CONSISTS OF
   !           1 ATOM RES AND ATOM CARDS MAY BE READ IN ANY ORDER
   !           END INPUT WITH AN "END " CARD
   
   !           ZERO NGRP BEFORE CALLING THIS ROUTINE
   !           ROUTINE WILL RETURN WITH THE NUMBER OF THE LAST GROUP READ
   
   ngrp = 0
   nf = nfu
   if (nf <= 0) nf=5
   
   !     ----- INITIALISE THE GROUP ARRAY -----
   
   do 100 i = 1,natom
      igroup(i) = 0
      if (konst) weit(i) = 0.0d0
   100 continue
   22 continue
   itime = 0
   lsign = 0
   izero = 0
   isrch = 0
   natmg = 0
   
   !       ----- READ DIFFERENT GROUPS -----
   
   read(nf,9208) (ititl1(k),k=1,19)
   if(ititl1(1) == iiend) then
      write(6, '(4x,a)') '----- END OF GROUP READ -----'
      goto 900
   end if
   ngrp = ngrp+1
   write(6,'(4x,a,i5,a)') '----- READING GROUP ', &
         ngrp, '; TITLE:'
   write(6,9218) (ititl1(k),k=1,19)
   
   !       ----- IF CONSTRAINED GROUPS READ THE WEIGHT FOR EACH GROUP -----
   
   if (konst) then
      ifld(1) = 3
      ifld(2) = 0
      call rfree(ifld,ihol,ivar,fvar,nf,6)
      wt = fvar(1)
      write(6,9018) ngrp,wt
   end if
   10 continue
   
   !       ----- READ THE GROUP CARDS -----
   
   ifld(1) = 1
   ifld(2) = 2
   call rfree(ifld,ihol,ivar,fvar,nf,6)
   ktypg = ihol(1)
   k = 1
   do 120 i = 1,7
      igrp(i) = ivar(k)
      jgrp(i) = ivar(k+1)
      k = k+2
   120 continue
   if (ktypg == iiend) goto 16
   if (ktypg == nfind) then
      write(6,199)
      isrch = 0
      goto 10
   end if
   if (ktypg == ifind) then
      
      !         ----- FIND OPTION ... READ THE ATOM SPECIFICS -----
      
      write(6,200)
      do 64 iii = 1,12
         call rfree(jfld,ihol,ivar,fvar,nf,6)
         jgraph(iii) = ihol(1)
         jsymbl(iii) = ihol(2)
         jtree(iii) = ihol(3)
         jresnm(iii) = ihol(4)
         if(jgraph(iii) == ksear) goto 65
         write(6,202) jgraph(iii),jsymbl(iii),jtree(iii),jresnm(iii)
      64 continue
      65 continue
      
      isrch = iii-1
      if(isrch > 10) write(6,66) isrch
      
      !         ----- NOW READ IN RES AND ATOMS TO BE SEARCHED -----
      
      goto 10
   end if
   itime = itime+1
   if (ktypg /= katn) then
      
      !         ----- 1ST GROUP OF 1ST CARD MUST BE - IF ANY - NUMBERS ARE
      !               FOUND -----
      
      if(itime == 1.and.igrp(1) < 0) lsign = 1
      do 12 i = 1,7
         i1 = igrp(i)
         if (i1 == 0) goto 10
         i2 = jgrp(i)
         if(i2 > nres) i2 = nres
         if(i1 > 0.and.lsign == 1) goto 36
         if(i1 < 0.and.lsign == 0) goto 36
         if(i1 < 0) i1 = -i1
         if(i2 <= 0) i2 = i1
         if(lsign == 0) write(6,14) ngrp,i1,i2
         do 13 j = i1,i2
            istart = ipres(j)
            iend = ipres(j+1)-1
            do 45 k = istart,iend
               if (isrch > 0) &
                     call findp(k,j,lfind,nres,ipres,lbres,isymbl, &
                     itree,igraph)
               if (isrch > 0.and.lfind == 0) goto 45
               igroup(k) = ngrp
               if(konst) weit(k) = wt
               natmg = natmg+1
            45 continue
            if(lsign == 1) write(6,46) ngrp,j
            if(lsign == 1) ngrp = ngrp+1
         13 continue
      12 continue
      goto 10
   end if
   
   !       ----- ATOM TYPE CONSTRAINTS -----
   
   if(lsign == 1) goto 36
   write(6,51) ngrp
   do 17 i = 1,7
      i1 = igrp(i)
      if(i1 < 0) goto 36
      if(i1 == 0) goto 10
      i2 = jgrp(i)
      if(i2 > natom) i2 = natom
      if(i2 <= 0) i2 = i1
      write(6,52) i1,i2
      do 18 j = i1,i2
         if(isrch > 0) &
               call findp(j,izero,lfind,nres,ipres,lbres,isymbl, &
               itree,igraph)
         if(isrch > 0.and.lfind == 0) goto 18
         natmg = natmg+1
         igroup(j) = ngrp
         if(konst) weit(j) = wt
      18 continue
   17 continue
   goto 10
   
   16 if(lsign == 1) ngrp = ngrp-1
   if(itime == 0) ngrp = ngrp-1
   !       IF(ISRCH.GT.0) WRITE(6,199)
   write(6,222) natmg
   goto 22
   
   36 continue
   write(6,127) ktypg,(igrp(i),jgrp(i),i = 1,7)
   goto 10
   
   !     ----- ALL GROUPS ARE READ RETURN -----
   
   900 continue
   if (konst) then
      
      !       ----- GATHER ALL THE CONSTRAINED ATOMS TOGETHER -----
      
      natc = 0
      do 920 i = 1,natom
         if(igroup(i) <= 0) goto 920
         natc = natc+1
         igroup(natc) = i
         weit(natc) = weit(i)
      920 continue
      
      !       ----- do not PRINT THE HISTORY OF CONSTRAINED ATOMS -----
      
#ifdef debug
      write(6,9108)
      do 940 i = 1,natc
         j = igroup(i)
         j3 = 3*j-3
         write(6,9118) i,j,igraph(j),isymbl(j),itree(j), &
               (xc(j3+k),k=1,3), weit(i)
      940 continue
#endif
   else if (.not.belly) then
      
      !       ----- PUT THE ATOMS WHICH ARE NOT IN THE DEFINED GROUPS
      !             AS THE LAST GROUP -----
      
      misc = .false.
      do 820 i = 1,natom
         if(igroup(i) /= 0) goto 820
         misc = .true.
         igroup(i) = ngrp+1
      820 continue
      if(misc) ngrp = ngrp+1
      !       IF(MISC) WRITE(6,9308) NGRP
   end if
   !  11 FORMAT(A4,1X,7(2I5))
   199 format(1h ,5x,'END OF ATOM SPECIFICATION',/)
   200 format(1h ,5x,'ALL ATOMS THAT MEET 1 OF THE FOLLOWING', &
         ' SPECIFICATIONS WILL BE INCLUDED IN GROUP BELOW',/)
   !  61 FORMAT(A4,A2,A1,A4)
   202 format(1h ,5x,'GRAPH NAME  = ',a4,2x,'SYMBOL  = ',a2,4x, &
         'TREE SYMBOL  = ',a1,5x,'RESIDUE TYPE  = ',a4,/)
   66 format(1h ,5x,'**** NUMBER OF FIND CARDS  = ',i5,2x, &
         'IS TOO BIG ******',/)
   14 format(' GRP',i5,' RES',i5,' TO ',i5)
   46 format(1h ,5x,'GROUP',i5,3x,'CONSISTS OF RESIDUE',i5,/)
   51 format(1h ,5x,'GROUP',i5,3x,'CONSISTS OF ATOMS -',/)
   52 format(1h ,34x,i5,2x,'TO',i5)
   222 format(1h ,5x,'Number of atoms in this group  = ',i5)
   127 format(1h ,5x,'***PROBLEMS WITH GROUP',a4,14i5,'*******',/)
   9018 format(/5x,'GROUP',i5,' HAS HARMONIC CONSTRAINTS',f12.5)
   9118 format(i5,i6,1x,a4,1x,2a4,3f10.4,f12.5)
   9108 format(/ /10x,'HISTORY OF CONSTRAINED ATOMS',/ /)
   9208 format(20a4)
   9218 format(1x,20a4)
   !9308 FORMAT(/5X,'THE GROUP ',I4, ' CONTAINS ALL ATOMS NOT DEFINED',
   !    +       ' AS GROUPS BY THE INPUT',/)
   return
end subroutine rgroup2 
!-----------------------------------------------------------------------
