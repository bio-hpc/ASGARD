C************************************************************************
C           usage: program_name<trajectory_file                         *
C                                                                       *
C this program creates *.crd.N files from the trajectory stream, where  *
C N is a snapshot number                                                *
C************************************************************************
C
C Modified by Holger Gohlke
C --- 4.9.2001
C       Multiple NUMBER_LIG_GROUPS / NUMBER_REC_GROUPS allowed
C
        implicit none   ! real*8(a-h,o-z), integer*4(i-n)
        integer na,ln,lp,nmr,nmo
C max number of atoms
        parameter (NA=500000)
C max length for line
        parameter (LN=80)
C max length for PREFIX and SUFFIX
        parameter (LP=80)
C NMR - max number of mutated residues
        parameter (NMR=40)
C NMO - max number of molecule groups
        parameter (NMO=10)

        character*9 ff
        data ff /'        '/
        character*(LN) line,title,line1,uppercase
c        character*1 line(LN),title(LN),line1(LN),uppercase(LN)
        character*(LP) PREFIX
        character*(LP) SUFFIX_R,SUFFIX_C,SUFFIX_L
        character*(LN) file_r,file_l,file_c,READBOX
        character*3 BOx
        character*(LN) noblank

        integer*4 ind,nn
        integer*4 i,nstart,nstop,nmg,nm,nm1,ma1,ma2,mk,mr
        integer*4 ntotal,nfreq,nprefix,nsuf_r,nsuf_c,nsuf_l
        integer*4 ngr1,ngr2,ngl1,ngl2
        integer*4 NGREC, RSTART(NMO), RSTOP(NMO)
        integer*4 NGLIG, LSTART(NMO), LSTOP(NMO)
        integer*4 NGCOM
        integer*4 RREAL, CREAL, LREAL
        integer*4 MUTANT_ATOM1(NMR),MUTANT_ATOM2(NMR)
        integer*4 MUTANT_KEEP(NMR),MUTANT_REFERENCE(NMR)
        integer*4 l,lfile_r,lfile_l,lfile_c,ifile,mfreq
        integer*4 j,k,ngr,ngl,ngc,mm

        real*8 x(NA), y(NA), z(NA)
        real*8 xR(NA),yR(NA),zR(NA)
        real*8 xC(NA),yC(NA),zC(NA)
        real*8 xL(NA),yL(NA),zL(NA)
        real*8 boxX, boxY, boxZ
        real*8 big, dmach
        real*8 xx,yy,zz,xxx,yyy,zzz,r
C
        print *,' '
        print *,'AmberTools14 Module: make_crg_hg'
        print *,' '
        print *,'usage: make_crd_hg <  trajectory_file '
        print *,' '
      
        call amopen( 10, 'make_crd.in', 'O', 'F', 'R' )
C
C defaults
C
C  On HPUX big=1.d30 was used instead of dmach;
C  If there is still a problem on HPUX then dmach probably can be 
C  rewritten using Fortran 95 intrinsics.  SRB Sep 2009
        big=dmach(3)-10
        NTOTAL=0
        NFREQ=1
        PREFIX=" "
        NPREFIX=0
        NSUF_R=0
        NSUF_C=0
        NSUF_L=0
        SUFFIX_R=" "
        SUFFIX_C=" "
        SUFFIX_L=" "
        file_c=" "
        file_l=" "
        file_r=" "

        NGREC=0
        NGR1=0
        NGR2=0
        NGLIG=0
        NGL1=0
        NGL2=0
        NGCOM=0
        do i=1,NMO
          RSTART(i)=0
          RSTOP(i)=0
          LSTART(i)=0
          LSTOP(i)=0
        end do

        NSTART=1
        NSTOP=999999999
        NMG=0
        NM=0
        do i=1,NMR
          MUTANT_ATOM1(i)=0
          MUTANT_ATOM2(i)=0
          MUTANT_KEEP(i)=0
          MUTANT_REFERENCE(i)=0
        end do
        BOX="NO "
        READBOX=" "
cMFC looks like it is not used        LIGAND="none  "
        NM1=0
        MA1=0
        MA2=0
        MK=0
        MR=0
        RREAL=0
        CREAL=0
        LREAL=0
C
C read settings and files names
C
10      read (10,'(a80)',end=20) line
        line1=uppercase(noblank(line,LN),LN)
        if(line1(1:1).eq.'#') goto 10
        if(line1(1:1).eq.'!') goto 10
        if(line1(1:6).eq.'NTOTAL') then
          ind=index(line,'=')
          read(line((ind+1):LN),*)NTOTAL
        endif
        if(line1(1:5).eq.'NFREQ')then 
          ind=index(line,'=')
          read(line((ind+1):LN),*)NFREQ
        endif
        if(line1(1:6).eq.'PREFIX')then
          ind=index(line,'=')
          PREFIX=noblank(line((ind+1):LN),LN-ind)
          endif
        if(line1(1:8).eq.'SUFFIX_R')then
          ind=index(line,'=')
          SUFFIX_R=noblank(line((ind+1):LN),LN-ind)
          endif
        if(line1(1:8).eq.'SUFFIX_C')then
          ind=index(line,'=')
          SUFFIX_C=noblank(line((ind+1):LN),LN-ind)
          endif
        if(line1(1:8).eq.'SUFFIX_L')then
          ind=index(line,'=')
          SUFFIX_L=noblank(line((ind+1):LN),LN-ind)
          endif
        if(line1(1:20).eq.'NUMBER_MUTANT_GROUPS')then
          ind=index(line,'=')
          read(line((ind+1):LN),*)NMG
          if(NMG.gt.NMR) then
            print *,'only ',NMR,' groups could be mutated'
            print *,'you mutated ',MK,' groups'
            print *,'recompile c.f with the larget NMR number'
            stop
          endif
        endif
        if(line1(1:12).eq.'MUTANT_ATOM1') then
          MA1=MA1+1
          if(MA1.gt.NMR) then
            print *,'only ',NMR,' groups could be mutated'
            print *,'you mutated ',MA1,' groups'
            print *,'recompile c.f with the larget NMR number'
            stop
          endif
          ind=index(line,'=')
          read(line((ind+1):LN),*)MUTANT_ATOM1(MA1)
        endif
        if(line1(1:12).eq.'MUTANT_ATOM2')  then
          MA2=MA2+1
          if(MA2.gt.NMR) then
            print *,'only ',NMR,' groups could be mutated'
            print *,'you mutated ',MA2,' groups'
            print *,'recompile c.f with the larget NMR number'
            stop
          endif
          ind=index(line,'=')
          read(line((ind+1):LN),*)MUTANT_ATOM2(MA2)
        endif
        if(line1(1:11).eq.'MUTANT_KEEP') then
          MK=MK+1
          if(MK.gt.NMR) then
            print *,'only ',NMR,' groups could be mutated'
            print *,'you mutated ',MK,' groups'
            print *,'recompile c.f with the larget NMR number'
            stop
          endif
          ind=index(line,'=')
          read(line((ind+1):LN),*)MUTANT_KEEP(MK)
        endif
        if(line1(1:16).eq.'MUTANT_REFERENCE') then
          MR=MR+1
          if(MR.gt.NMR) then
            print *,'only ',NMR,' groups could be mutated'
            print *,'you mutated ',MR,' groups'
            print *,'recompile c.f with the larget NMR number'
            stop
          endif
          ind=index(line,'=')
          read(line((ind+1):LN),*)MUTANT_REFERENCE(MR)
        endif

        if(index(line,'NSTART').ne.0)then
          ind=index(line,'=')
          read(line((ind+1):LN),*)NSTART
        endif
        if(index(line,'NSTOP').ne.0)then
          ind=index(line,'=')
          read(line((ind+1):LN),*)NSTOP
        endif

        if(line1(1:17).eq.'NUMBER_REC_GROUPS')then
          ind=index(line,'=')
          read(line((ind+1):LN),*)NGREC
          if(NGREC.gt.NMO) then
            print *,'only ',NMO,' groups allowed for NUMBER_REC_GROUPS'
            stop
          endif
        endif
        if(line1(1:17).eq.'NUMBER_LIG_GROUPS')then
          ind=index(line,'=')
          read(line((ind+1):LN),*)NGLIG
          if(NGLIG.gt.NMO) then
            print *,'only ',NMO,' groups allowed for NUMBER_LIG_GROUPS'
            stop
          endif
        endif

        if(line1(1:6).eq.'RSTART') then
          NGR1=NGR1+1
          if(NGR1.gt.NMO) then
            print *,'only ',NMO,' groups allowed for NGR1'
            stop
          endif
          ind=index(line,'=')
          read(line((ind+1):LN),*)RSTART(NGR1)
        endif
        if(line1(1:5).eq.'RSTOP') then
          NGR2=NGR2+1
          if(NGR2.gt.NMO) then
            print *,'only ',NMO,' groups allowed for NGR2'
            stop
          endif
          ind=index(line,'=')
          read(line((ind+1):LN),*)RSTOP(NGR2)
        endif

        if(line1(1:6).eq.'LSTART') then
          NGL1=NGL1+1
          if(NGL1.gt.NMO) then
            print *,'only ',NMO,' groups allowed for NGL1'
            stop
          endif
          ind=index(line,'=')
          read(line((ind+1):LN),*)LSTART(NGL1)
        endif
        if(line1(1:5).eq.'LSTOP') then
          NGL2=NGL2+1
          if(NGL2.gt.NMO) then
            print *,'only ',NMO,' groups allowed for NGL2'
            stop
          endif
          ind=index(line,'=')
          read(line((ind+1):LN),*)LSTOP(NGL2)
        endif

        if(index(line,'BOX').ne.0)then
          ind=index(line,'=')
          READBOX=uppercase(noblank(line((ind+1):LN),LN-ind),LN-ind)
        endif
        goto 10
20      close(10)

        NPREFIX=index(PREFIX,' ')-1
        NSUF_R=index(SUFFIX_R,' ')-1
        NSUF_C=index(SUFFIX_C,' ')-1
        NSUF_L=index(SUFFIX_L,' ')-1
        if(NSUF_L.eq.0)then
          NGLIG=0
        endif
        if(NSUF_C.ne.0)then
          NGCOM=1
        endif
        if(index(READBOX,'Y').ne.0)BOX="YES"
        if(index(READBOX,'N').ne.0)BOX="NO "
C      	print *, NATOM,NTOTAL,NLIGAND,NFREQ,NPREFIX,NSUF_R,NSUF_C,NSUF_L
C      	print 1000, PREFIX, SUFFIX_R,SUFFIX_L,SUFFIX_C,BOX,LIGAND	
C1000	format(a40/)
C
C check sanity
C
        if(NTOTAL.eq.0) then
          print *,"NTOTAL is not determined. See make_crd.in"
          stop
        elseif(NPREFIX.lt.1) then
          print *,
     1  "PREFIX is not determined or have spaces. See make_crd.in"
          stop
        elseif((NGR1.ne.NGR2).or.(NGR1.ne.NGREC)) then
          print *,"Something wrong with receptor groups"
          stop
        elseif((NGL1.ne.NGL2).or.(NGL1.ne.NGLIG)) then
          print *,"Something wrong with ligand groups"
          stop
        endif
C
C check mutant atoms assignment (Thr group has two mutant atoms:
C
        if((MUTANT_ATOM1(1).eq.0).and.(MUTANT_ATOM2(1).eq.0)) then
          MA1=0
          MA2=0
          MK=0
          MR=0
          NMG=0
        endif
        if(MA1.ne.NMG)then
          print *,'MUTANT_ATOM1 was not assigned for some group(s)'
          stop
        endif
        if(MA2.ne.NMG)then
          print *,'MUTANT_ATOM2 was not assigned for some group(s)'
          stop
        endif
        if(MR.ne.NMG)then
          print *,'MUTANT_REFERENCE was not assigned for some group(s)'
          stop
        endif
        if(MK.ne.NMG)then
          print *,'MUTANT_KEEP was not assigned for some group(s)'
          stop
        endif
        do i=1,NMG
        if((MUTANT_ATOM1(i).eq.0).and.(MUTANT_ATOM2(i).gt.0)) then
          MUTANT_ATOM1(i)=MUTANT_ATOM2(i)
          MUTANT_ATOM2(i)=0
        elseif((MUTANT_ATOM2(i).ne.0).and.
     1   (MUTANT_ATOM1(i).ne.0)) then
          if(MUTANT_ATOM1(i).gt.MUTANT_ATOM2(i)) then
            MM=MUTANT_ATOM1(i)
            MUTANT_ATOM1(i)=MUTANT_ATOM2(i)
            MUTANT_ATOM2(i)=MM
          endif
        endif
        enddo
C
C check receptor sizes
C
        if(NGREC.eq.0) then
          NGREC=1
          RSTART(1)=1
          RSTOP(1)=NTOTAL
        endif          
        do i=1,NGREC
          if(RSTART(i).gt.RSTOP(i)) then
            print *,"ERROR in size definition for the receptor"
            print *,"RSTART must be smaller or equal to RSTOP"
            print *,"see make_crd.in"
            stop
          endif
        end do
        do i=1,NGLIG
          if(LSTART(i).gt.LSTOP(i)) then
            print *,"ERROR in size definition for the ligand"
            print *,"LSTART must be smaller or equal to LSTOP"
            print *,"see make_crd.in"
            stop
          endif
        end do
C
C check ligand and mutant atoms definitions
C
        do i=1,NGREC
          RREAL=RREAL+RSTOP(i)-RSTART(i)+1
          CREAL=CREAL+RSTOP(i)-RSTART(i)+1
        end do
        do i=1,NGLIG
          LREAL=LREAL+LSTOP(i)-LSTART(i)+1
          CREAL=CREAL+LSTOP(i)-LSTART(i)+1
        end do

        do i=1,NMG
          if((MUTANT_ATOM1(i).ne.0).and.
     1      (MUTANT_ATOM2(i).ne.0).and.
     1      (MUTANT_KEEP(i).ne.0).and.
     1      (MUTANT_REFERENCE(i).ne.0)) then
            CREAL=CREAL
     1      -(MUTANT_KEEP(i)-MUTANT_ATOM1(i)-2)
            do j=1,NGLIG
              if((MUTANT_ATOM1(i).lt.LSTOP(j)).and.
     1           (MUTANT_ATOM1(i).gt.LSTART(j)))then
                LREAL=LREAL
     1            -(MUTANT_KEEP(i)-MUTANT_ATOM1(i)-2)
              endif
            end do
            do j=1,NGREC
              if((MUTANT_ATOM1(i).lt.RSTOP(j)).and.
     1           (MUTANT_ATOM1(i).gt.RSTART(j)))then
                RREAL=RREAL
     1            -(MUTANT_KEEP(i)-MUTANT_ATOM1(i)-2)
              endif
            end do
          elseif((MUTANT_ATOM1(i).ne.0).and.(MUTANT_ATOM2(i).eq.0))then
            CREAL=CREAL
     1      -(MUTANT_KEEP(i)-MUTANT_ATOM1(i)-1)
            do j=1,NGLIG
              if((MUTANT_ATOM1(i).lt.LSTOP(j)).and.
     1           (MUTANT_ATOM1(i).gt.LSTART(j)))then
                LREAL=LREAL
     1            -(MUTANT_KEEP(i)-MUTANT_ATOM1(i)-1)
              endif
            end do
            do j=1,NGREC
              if((MUTANT_ATOM1(i).lt.RSTOP(j)).and.
     1           (MUTANT_ATOM1(i).gt.RSTART(j)))then
                RREAL=RREAL
     1            -(MUTANT_KEEP(i)-MUTANT_ATOM1(i)-1)
              endif
            end do
          endif
        enddo
C	print *, "CREAL=",CREAL," LREAL=",LREAL," RREAL=",RREAL
C
C start build crd files
C
30      nn=0
        Lfile_r=10
        Lfile_c=11
        Lfile_l=12
        ifile=0
        read (5,'(a80)',end=100) title
40      nn=nn+1
        read (5,'(10f8.3)',end=100)(x(mm),y(mm),z(mm),mm=1,NTOTAL)
        if(BOX.eq.'YES') then
          read (5,'(3f8.3)') boxX, boxY, boxZ
          write(6,*) 'Box info found: ',boxX, boxY, boxZ
        end if
        if ((nn.ge.NSTART).and.(nn.le.NSTOP).and. 
     1      ((mod((nn-NSTART), NFREQ).eq.0).or.(NFREQ.eq.1))) then
          ifile= ifile+1
          if(ifile.lt.10) then
            write(ff,'(i1,8x)')ifile
          elseif(ifile.lt.100) then
            write(ff,'(i2,7x)')ifile
          elseif(ifile.lt.1000) then
            write(ff,'(i3,6x)')ifile
          elseif(ifile.lt.10000) then
            write(ff,'(i4,5x)')ifile
          elseif(ifile.lt.100000) then
            write(ff,'(i5,4x)')ifile
          elseif(ifile.lt.1000000) then
            write(ff,'(i6,3x)')ifile
          elseif(ifile.lt.10000000) then
            write (6,*) 
     1        'Do not know how creat name for file number',ifile
            write (6,*) 
     1        'It is too big, change settings and recompile make_crd.f'
            stop
          endif
C
C claculate the coordinates for the mutated atoms
C
          do i=1,NMG
            xx=0d0
            yy=0d0
            zz=0d0
            xxx=0d0
            yyy=0d0
            zzz=0d0
            r=0
            if((MUTANT_ATOM1(i).ne.0).and.(MUTANT_ATOM2(i).ne.0))then
              xx=x(MUTANT_REFERENCE(i))
              yy=y(MUTANT_REFERENCE(i))
              zz=z(MUTANT_REFERENCE(i))
              xxx=x(MUTANT_ATOM1(i))
              yyy=y(MUTANT_ATOM1(i))
              zzz=z(MUTANT_ATOM1(i))
              r=dsqrt((xxx-xx)*(xxx-xx)+(yyy-yy)*(yyy-yy)+
     1          (zzz-zz)*(zzz-zz))
              x(MUTANT_ATOM1(i))=xx+(xxx-xx)*1.09/r
              y(MUTANT_ATOM1(i))=yy+(yyy-yy)*1.09/r
              z(MUTANT_ATOM1(i))=zz+(zzz-zz)*1.09/r
              xxx=x(MUTANT_ATOM2(i))
              yyy=y(MUTANT_ATOM2(i))
              zzz=z(MUTANT_ATOM2(i))
              r=dsqrt((xxx-xx)*(xxx-xx)+(yyy-yy)*(yyy-yy)+
     1          (zzz-zz)*(zzz-zz))
              x(MUTANT_ATOM2(i))=xx+(xxx-xx)*1.09/r
              y(MUTANT_ATOM2(i))=yy+(yyy-yy)*1.09/r
              z(MUTANT_ATOM2(i))=zz+(zzz-zz)*1.09/r
              if(MUTANT_ATOM1(i)+1.gt.MUTANT_ATOM2(i)-1) goto 500
              do j=MUTANT_ATOM1(i)+1,MUTANT_ATOM2(i)-1
                x(j)=big+10
                y(j)=big+10
                z(j)=big+10
              enddo
500       continue
              if(MUTANT_ATOM2(i)+1.gt.MUTANT_KEEP(i)-1) goto 510
              do j=MUTANT_ATOM2(i)+1,MUTANT_KEEP(i)-1
                x(j)=big+10
                y(j)=big+10
                z(j)=big+10
              enddo
510       continue
            elseif((MUTANT_ATOM1(i).ne.0).and.(MUTANT_ATOM2(i).eq.0))
     1        then
              xx=x(MUTANT_REFERENCE(i))
              yy=y(MUTANT_REFERENCE(i))
              zz=z(MUTANT_REFERENCE(i))
              xxx=x(MUTANT_ATOM1(i))
              yyy=y(MUTANT_ATOM1(i))
              zzz=z(MUTANT_ATOM1(i))
              r=dsqrt((xxx-xx)*(xxx-xx)+(yyy-yy)*(yyy-yy)+
     1          (zzz-zz)*(zzz-zz))
              x(MUTANT_ATOM1(i))=xx+(xxx-xx)*1.09/r
              y(MUTANT_ATOM1(i))=yy+(yyy-yy)*1.09/r
              z(MUTANT_ATOM1(i))=zz+(zzz-zz)*1.09/r
              if(MUTANT_ATOM1(i)+1.gt.MUTANT_KEEP(i)-1) goto 520
              do j=MUTANT_ATOM1(i)+1,MUTANT_KEEP(i)-1
                x(j)=big+10
                y(j)=big+10
                z(j)=big+10
              enddo
520       continue
            endif
          enddo

          j=1
          k=1
          NGR=0
          NGL=0
          NGC=0
          do i=1,NGREC+NGLIG
            if(((NGREC.gt.0).and.(NGREC.ge.j).and.
     1          (RSTART(j).lt.LSTART(k))).or.(k.gt.NGLIG)) then
              do l=RSTART(j),RSTOP(j)
                if(x(l).lt.big) then
                  NGR=NGR+1
                  xR(NGR)=x(l)
                  yR(NGR)=y(l)
                  zR(NGR)=z(l)
                  if(NGCOM.gt.0) then
                    NGC=NGC+1
                    xC(NGC)=x(l)
                    yC(NGC)=y(l)
                    zC(NGC)=z(l)
                  endif
                endif
              end do
              j=j+1
            elseif(((NGLIG.gt.0).and.(NGLIG.ge.k).and.
     1              (LSTART(k).lt.RSTART(j))).or.(j.gt.NGREC)) then
              do l=LSTART(k),LSTOP(k)
                if(x(l).lt.big) then
                  NGL=NGL+1
                  xL(NGL)=x(l)
                  yL(NGL)=y(l)
                  zL(NGL)=z(l)
                  if(NGCOM.gt.0) then
                    NGC=NGC+1
                    xC(NGC)=x(l)
                    yC(NGC)=y(l)
                    zC(NGC)=z(l)
                  endif
                endif
              end do
              k=k+1
            endif
          end do

C
C Check again
C
          if(NSUF_R.gt.0.and.NGR.ne.RREAL) then
            print *,"Something wrong with NGR",NGR," ",RREAL
            stop
          elseif(NSUF_L.gt.0.and.NGL.ne.LREAL) then
            print *,"Something wrong with NGL",NGL," ",LREAL
            stop
          elseif(NSUF_C.gt.0.and.NGC.ne.CREAL) then
            print *,"Something wrong with NGC",NGC," ",CREAL
            stop
          endif
C
C Print crd files
C
          if(NSUF_R.gt.0) then
CCC            if(BOX.eq.'YES')
CCC             --- Image coordinates
CCC     +        call imagecrd(boxX,boxY,boxZ,xR,yR,zR,RREAL)
C The following did not work on OSF1 with Machine compilation keys
C           file_r=PREFIX(1:NPREFIX) // SUFFIX_R(1:NSUF_R) // ".crd." // ff
            file_r(1:NPREFIX)=PREFIX(1:NPREFIX)
            file_r((NPREFIX+1):(NPREFIX+NSUF_R))=SUFFIX_R(1:NSUF_R)
            file_r((NPREFIX+NSUF_R+1):(NPREFIX+NSUF_R+5))=".crd."
            file_r((NPREFIX+NSUF_R+6):(NPREFIX+NSUF_R+15))=ff
            call amopen( Lfile_r, file_r, 'U', 'F', 'W' )
            write(Lfile_r,'(a)') title
            write(Lfile_r,'(i5)') RREAL
            write(Lfile_r,'(6f12.7)')(xR(mm),yR(mm),zR(mm),mm=1,RREAL)
            close(Lfile_r)
          endif

          if(NSUF_L.gt.0) then
CCC            if(BOX.eq.'YES')
CCC             --- Image coordinates
CCC     +        call imagecrd(boxX,boxY,boxZ,xL,yL,zL,LREAL)
C           file_l=PREFIX(1:NPREFIX) // SUFFIX_L(1:NSUF_L) // ".crd." // ff
            file_l(1:NPREFIX)=PREFIX(1:NPREFIX)
            file_l((NPREFIX+1):(NPREFIX+NSUF_L))=SUFFIX_L(1:NSUF_L)
            file_l((NPREFIX+NSUF_L+1):(NPREFIX+NSUF_L+5))=".crd."
            file_l((NPREFIX+NSUF_L+6):(NPREFIX+NSUF_L+15))=ff
            call amopen( Lfile_l, file_l, 'U', 'F', 'W' )
            write(Lfile_l,'(a)') title
            write(Lfile_l,'(i5)') LREAL
            write(Lfile_l,'(6f12.7)')(xL(mm),yL(mm),zL(mm),mm=1,LREAL)
            close(Lfile_l)
          endif

          if(NSUF_C.gt.0) then
CCC            if(BOX.eq.'YES')
CCC             --- Image coordinates
CCC     +        call imagecrd(boxX,boxY,boxZ,xC,yC,zC,CREAL)
C           file_c=PREFIX(1:NPREFIX) // SUFFIX_C(1:NSUF_C) // ".crd." // ff
            file_c(1:NPREFIX)=PREFIX(1:NPREFIX)
            file_c((NPREFIX+1):(NPREFIX+NSUF_C))=SUFFIX_C(1:NSUF_C)
            file_c((NPREFIX+NSUF_C+1):(NPREFIX+NSUF_C+5))=".crd."
            file_c((NPREFIX+NSUF_C+6):(NPREFIX+NSUF_C+15))=ff
            call amopen( Lfile_c, file_c, 'U', 'F', 'W' )
            write(Lfile_c,'(a)') title
            write(Lfile_c,'(i5)') CREAL
            write(Lfile_c,'(6f12.7)')(xC(mm),yC(mm),zC(mm),mm=1,CREAL)
            close(Lfile_c)
          endif
        endif

        goto 40
100     stop
        end
C
C************************************************************************
C
        subroutine imagecrd(boxX,boxY,boxZ,x,y,z,nat)
C
        implicit none    !double precision (a-h,o-z)
C
        integer nat,i
        real*8 boxx,boxy,boxz
        real*8 x(*), y(*), z(*)
C
        do i=1,nat
          x(i) = x(i) - boxX * int(x(i) / boxX)
          y(i) = y(i) - boxY * int(y(i) / boxY)
          z(i) = z(i) - boxZ * int(z(i) / boxZ)
        end do
C
        return
        end
C
C************************************************************************
C
      character*(*) function noblank(line,n)
      implicit none
      integer*4 n,i
      character*(*) line
      do i=1,n
         if((line(i:i).ne.' ').and.(line(i:i).ne.'	'))go to 10
      enddo
 10   if(i.le.n) then
         noblank=line(i:n)
      else
         noblank=' '
      endif 
      return
      end
C
C************************************************************************
C
      real*8 function dmach(job)
      implicit none
      integer job
c     
c     dmach computes machine parameters of floating point
c     arithmetic for use in testing only.  
c     
c     if trouble with automatic computation of these quantities,
c     they can be set by direct assignment statements.
c     assume the computer has
c
c        b = base of arithmetic
c        t = number of base  b  digits
c        l = smallest possible exponent
c        u = largest possible exponent
c
c     then
c
c        eps = b**(1-t)
c        tiny = 100.0*b**(-l+t)
c        huge = 0.01*b**(u-t)
c
c     job is 1, 2 or 3 for epsilon, tiny and huge, respectively.
c
      real*8 eps,tiny,huge,s
c
      eps = 1.0d0
 10   eps = eps/2.0d0
      s = 1.0d0 + eps
      if (s .gt. 1.0d0) go to 10
      eps = 2.0d0*eps
c     
      s = 1.0d0
 20   tiny = s
      s = s/16.0d0
      if (s*1.0 .ne. 0.0d0) go to 20
      tiny = (tiny/eps)*100.0
      huge = 1.0d0/tiny
c     
      if (job .eq. 1) dmach = eps
      if (job .eq. 2) dmach = tiny
      if (job .eq. 3) dmach = huge
      return
      end
C************************************************************************
      character*(*) function uppercase(line,n)
      implicit none
      integer*4 n,i,k
      character*(*) line
      character*26 low,up
      data low /'abcdefgijhklmnopqrstuvwxyz'/
      data up /'ABCDEFGIJHKLMNOPQRSTUVWXYZ'/
      uppercase=line
      do i=1,n
         k=index(low,uppercase(i:i))
         if(k.ne.0) uppercase(i:i)=up(k:k)
      enddo
      return
      end
C************************************************************************

