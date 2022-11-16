      PROGRAM FANTASIAN
c     A program to evaluate magnetic anisotropy tensor parameters
c
c     by
c
c     CERM and  Department of Chemistry, University of Florence
c     Florence, Italy
c    
c     Principal Investigator: Prof. Ivano Bertini
c     Refs:
c  	L. Banci, I. Bertini, K.L. Bren, M.A. Cremonini, H.B. Gray, C.  Luchinat,
c      P. Turano: The use of pseudocontact shifts to refine solution structures
c      of paramagnetic metalloproteins: Met80Ala cyano-cytochrome c as an 
c     example.  JBIC. 1, 117-126, 1996
c	L. Banci, I. Bertini, G. Gori Savellini, A. Romagnoli, P. Turano, M.A.
c      Cremonini, C. Luchinat, H.B. Gray: The pseudocontact shifts as
c      constraints for energy minimization and molecular dynamic calculations
c      on solution structures of paramagnetic metalloproteins.  Proteins
c      Struct.Funct.Genet. 29, 68-76, 1997
c 	L. Banci, I. Bertini, M.A. Cremonini, G. Gori Savellini, C. Luchinat, K.
c      Wüthrich, P. Guentert: PSEUDODYANA for NMR structure calculation of
c      paramagnetic metalloproteins using torsion angle molecular dynamics.
c      J.Biomol.NMR. 12, 553-557, 1998


c     ------------------------------------------------------------
c     Revision made by Gabriele Baldi, October 1999.
c     Subroutine PANG has been rewritten.
c     ------------------------------------------------------------

#include "sup.h"
      character *(*) versn
      parameter (versn='VERSION 2.0 - Nov 1999') ! *** VERSION ***
      
      dimension npt(MAXFE,6), ! Label numbers of ref. systems atoms
     *   ccptx(MAXSTR,MAXFE,6),  ! Coordinates of atoms
     *   ccpty(MAXSTR,MAXFE,6),  !  defining the reference
     *   ccptz(MAXSTR,MAXFE,6),  !  systems.
     *   vx(3), vy(3),vz(3)
      dimension ion(MAXFE)    ! Atom numbers of param. ions
      character char*1        ! Flag for manual input
      CHARACTER varerr*80     ! Error management string
      CHARACTER errmsg*80     ! Error management string
      
      logical InECHO  
      parameter (InECHO=.true.) ! INPUT ECHO FLAG (echo/silent)
      
      gltoll=0.3              ! Default measure tolerance
      
      print *, ' '
      print *, 'Amber11 Module: fantasian'
      print *, ' '
      print *, 'This is Fantasian ',versn ! Version
      
      if(DEBUGP) then ! Output will show debug printing
        print *, ' '
        print *,' =========================='
        print *,' == Debug printing is ON =='
        print *,' =========================='
      endif
      
c	********************
c	*** MANUAL INPUT ***
c	********************
      print *, ' '

      print *, 'Name of observed shifts input file:'
      read (*,'(a)') filename1
      if(InECHO) print '(3x,a)', filename1
      print *, 'Name of Amber pdb input file:'
      read (*,'(a)') filename2
      if(InECHO) print '(3x,a)', filename2
      print *, 'Name for observed out file:'
      read (*,'(a)') filename3 
      if(InECHO) print '(3x,a)', filename3 
      print *, 'Name for output file:'
      read (*,'(a)') fileout
      if(InECHO) print '(3x,a)', fileout
      print *,'Number of atoms in each structure (not more than'
     *,MAXDEF,'):'
      read (*,'(i4)') nat
      if(InECHO) print '(3x,i4)', nat
      if(nat.gt.MAXDEF) then
        errmsg=' Too many atoms in each structure'
	varerr=' code:MAXDEF'
	goto 22
      else if (nat.lt.1) then
        errmsg=' Invalid input parameter'
	varerr=' code:neg'
	goto 22
      endif
      
      print *, 'Number of structures (not more than',MAXSTR,'):'
      read (*,'(i4)') nstr
      if(InECHO) print '(3x,i4)', nstr
      if(nstr.gt.MAXSTR) then
        errmsg=' Too many structures required'
	varerr=' code:MAXSTR'
	goto 22
      else if (nstr.lt.1) then
        errmsg=' Invalid input parameter'
	varerr=' code:neg'
	goto 22
      endif
      
      print *, 'Number of paramagnetic ions (not more than',MAXFE,'):'
      read (*, '(i4)') nfe
      if(InECHO) print '(3x,i4)', nfe
      if(nfe.gt.MAXFE) then
	errmsg=' Too many paramagnetic ions'
	varerr=' code:MAXFE'
	goto 22
      else if (nfe.lt.1) then
        errmsg=' Invalid input parameter'
        varerr=' code:neg'
	goto 22
      endif
      
      do j=1,nfe
        print *,
     *      '  Insert atom number corresponding to paramagnetic ion #',j
        read(*,'(i5)') ion(j)
        if(InECHO) print '(4x,i5)',ion(j)
        if(ion(j).lt.1.or.ion(j).gt.nat) then
          errmsg=' Invalid atom number'
	  varerr=' code:atom'
	  goto 22
        endif
      enddo
      
      print *, 'Value of tolerance (0-1):'
      read (*, '(f8.3)') PERC
      if(InECHO) print '(3x,f8.3)', PERC
      if(PERC.lt.0.or.PERC.gt.1) then
        errmsg=' Value of tolerance out of range'
	varerr=' code:range'
	goto 22
      endif
      
      print *,' '
      print *,'Would you define reference system(s)? (y/n) '
      print *,' '
      read (*,'(a1)') char
      if(InECHO) print '(3x,a)', char
      
      if (char.eq.'y'.or.char.eq.'Y') then
        intsys=1 ! Flag: ref. syst. defined by user
	print *,' '
	print *,'  Unique or multiple reference system(s)? (u/m)'
        print *,' '
	read (*,'(a1)') char
	if(InECHO) print '(4x,a)', char
	if (char.eq.'u'.or.char.eq.'U') then
	  nsystem=1
	else if (char.eq.'m'.or.char.eq.'M') then
	  nsystem=nfe
	else  
	  errmsg=' Invalid input parameter'
          varerr=' code:neg'
	  goto 22
	endif
	 	
        print *,' '
      else
        intsys=0    ! Flag: default reference systems
        nsystem=nfe ! Default nr. of reference systems
      endif
      
c 	Define the required reference systems
c
c	nresno(j) Residual label-number
c	npt(j,i)  Label-numbers of 3 atoms defining the reference system
c		  (this three atoms will be the same for each structure)
      if (intsys.eq.1) then
        do j=1,nsystem
          print *,'Defining reference system nr.',j,'/',nsystem,':'
	  print *,' Insert three atom numbers:'
          print *,'  - first fixes the x direction'
          read(*,'(i5)') npt(j,1)
          if(InECHO) print '(5x,i5)', npt(j,1)
          if(npt(j,1).lt.1.or.npt(j,1).gt.nat) then
            errmsg=' Invalid atom number'
	    varerr=' code:atom'
	    goto 22
          endif
          print *,'  - second is the system origin'
          read(*,'(i5)') npt(j,2)
          if(InECHO) print '(5x,i5)', npt(j,2)
          if(npt(j,2).lt.1.or.npt(j,2).gt.nat) then
            errmsg=' Invalid atom number'
	    varerr=' code:atom'
	    goto 22
          endif
          print *,'  - third is an other point of the plane x-y'
          read(*,'(i5)') npt(j,3)
          if(InECHO) print '(5x,i5)', npt(j,3)
          if(npt(j,3).lt.1.or.npt(j,3).gt.nat) then
            errmsg=' Invalid atom number'
	    varerr=' code:atom'
	    goto 22
          endif
          print *,' '
        enddo
      endif
      
c                                griglia sulle condizioni iniziali
      print *,'Number of computed solutions:'
      read (*,'(i4)') ngrid
      if(InECHO) print '(3x,i4)', ngrid

*     **** END OF MANUAL INPUT ****


* ======================================================================
* ======================================================================
* ======================================================================

c    *****  READ OBSERVED SHIFTS FROM PSEUDO-CONTACTS INPUT FILE  ******
c
c       numres(i)  Residual label number
c       namres(i)  Residual name
c       namat(i)   Name of atom whose shift is observed
c       obs(i)     Observed shift
c       mlprot(i)  Molteplicity (1 or 2 or 3)
c       tolprot(i) Measure tolerance
c       wprot(i)   Weight
c
c     open pseudocontacts input file (assigned: #1)
      print *,' '
      print *,' ------------ S T A R T -------------'
      print *,' '
      open (1,file=filename1,status='old')
      i=0     ! Record counter
      
      print *,'     *****     READING OBSERVED SHIFT:'
      
c     Read-file loop:
      do while (1.gt.0)
        read(1,'(a)',end=99) line ! Read a record from file
        i=i+1
	if(i.gt.MAXOS) then
	  errmsg=' Too many observations in input file'
	  varerr=' code:MAXOS'
	  goto 22
	endif
	read(line,'(i4,1x,a3,1x,a4,1x,f8.3,4x,i3,f8.3,f8.3)') ! Single fields
     *           numres(i),namres(i),namat(i),obs(i),
     *           mlprot(i),tolprot(i),wprot(i)

        if (tolprot(i).lt.EPS) tolprot(i)=gltoll ! Assign default value
        if ((PERC.gt.0).and.(PERC.le.1)) then 
          ttperc=abs(obs(i)*PERC)
          if (ttperc.gt.tolprot(i)) tolprot(i)=ttperc
        end if

        if (mlprot(i).eq.0) mlprot(i)=1 ! Default molteplicity is 1

c	Assigning weight dafaults
        if (wprot(i).eq.0.) wprot(i)=1.0
        if (wprot(i).eq.-1.) wprot(i)=0.0
      enddo
      
 99   ihp=i  ! Nr. of records (lines) read from pseudocontact input file

      print *, '     *****     TOTAL DIPOLAR PROTONS ',ihp
      
* ======================================================================

c    *****  READ STRUCTURE DATA FROM STRUCTURE INPUT FILE  ******
c
c     Open molecule structure input file (assigned: #2)
      open (2,file=filename2,status='old')
c     Open observed otuput file (assigned: #3)
      open (3,file=filename3)
      
c		*****     READ THE  FILE    *****
c
c       num_at    progressive atom label number
c       nam_at    atom name
c       nam_res   residual name
c       num_res   progressive residual label number
c       xp,yp,zp  atom coordinates
c
      print *,'     *****     READING DATA FROM COORDINATE FILE:'
      print *,' '

      do 10 k=1,nstr         ! Loop on each structure
        if(DEBUGP) print *,'Structure',k
      	knat=(k-1)*nat	     ! Index offset for structure k
	
	do 20 i=1,nat        ! Loop on each atom of structure k
	  read(2,'(a)',end=23) line ! Read a record (input file line)
	  
c         Get single fields from record:
          read(line,'(5x,i6,1x,a4,1x,a3,1x,i5,4x,3f8.3)')
     *         num_at(k,i),nam_at(k,i),nam_res(k,i),num_res(k,i),
     *         xp(k,i),yp(k,i),zp(k,i)
     
c    Find out the specified paramagnetic ions; their coordinates are
c    stored in: fx(k,-), fy(k,-), fz(k,-)
c
	  do 13 j=1,nfe
	    if(num_at(k,i).eq.ion(j)+knat) then
	      fx(k,j)=xp(k,i)  ! Store paramagnetic ion coordinates
	      fy(k,j)=yp(k,i)
	      fz(k,j)=zp(k,i)
	      if(DEBUGP.or.k.eq.1)
     *	      print *,' - Paramagn. ion nr.',j,'is ',nam_at(k,i),
     *	              ' (in residue nr.',num_res(k,i),': ',
     *	              nam_res(k,i),')' 
	    end if
 13	  continue
	  
c
c    Find out atoms for the required reference systems; their coordinates are
c    stored in:
c	ccptx(k,j,l)   ccpty(k,j,l)   ccptz(k,j,l)     for l=1,2,3.
c
          if (intsys.eq.1) then   ! Reference systems were required
            do 15 j=1,nsystem     ! Loop on each reference system
            do 15 l=1,3           ! Loop on 3 atoms defining ref. sys. j
              ! Does this atom define a reference system ?
              if (num_at(k,i).eq.npt(j,l)+knat) then
                ccptx(k,j,l)=xp(k,i)
                ccpty(k,j,l)=yp(k,i)  ! Store atom coordinates
                ccptz(k,j,l)=zp(k,i)
              endif
 15         continue
          endif

 20	continue
         if(DEBUGP) print *,' '
 10   continue
      print *,' '
      goto 25

c	**** END OF INPUT FILES READING ****
* ======================================================================

c     Underdimensioned variable error management
 22   print *,' *** ERROR: ',errmsg
      print *,varerr
      print *,' *** PROGRAM STOP.'
      print *,' '
      STOP
      
c     Unexpected end of file while reading coordinates input file
 23   print *,' *** ERROR: Unexpected end of file: ',filename2
      print *,'     (the file was supposed to contain',nstr*nat,
     *'records)'
      print *,' '
      print *,' *** PROGRAM STOP.' 
      print *,' ' 
      STOP

* ======================================================================
 25   continue
c
c	*****     CALCULATE THE REFERENCE SYSTEM	****
c
      if (intsys.eq.1) then ! OK: Reference systems were required
c        print *,'     *****     CALCULATE THE REFERENCE SYSTEM'
c       Medium points (on every structure)
	do 30 j=1,nsystem ! j = index of ref. system
          do 30 m=1,3    ! m = the three points defining the ref. system
            rccptx(j,m) = 0.0 
            rccpty(j,m) = 0.0
            rccptz(j,m) = 0.0
 30     continue
	do 40 j=1,nsystem ! For every reference system:
          do 40 k=1,nstr ! Compute medium points between structures
            do 40 m=1,3
              rccptx(j,m) = rccptx(j,m) + ccptx(k,j,m)/dble(nstr)
              rccpty(j,m) = rccpty(j,m) + ccpty(k,j,m)/dble(nstr)
              rccptz(j,m) = rccptz(j,m) + ccptz(k,j,m)/dble(nstr)
 40     continue

c	For every ref. syst. j, versors Z, X and Y are computed:
c	point 2 is the origin,
c	point 1 lays on x-axis,
c	point 3 defines the x-y plane, so that:
c	(P1-P2)x(P3-P2) gives the Z direction,
c	(P1-P2) gives the X direction,
c	(Z)x(X) gives the Y direction.
c      NOTE: These constraints allow two different solutions:
c	     -> Z=(P1-P2)x(P3-P2), and also
c	     -> Z=(P3-P2)x(P1-P2) (opposite in sign)
c	     This doesn't matter: in both cases final result 
c	     will be the same.

	do 50 j=1,nsystem ! Loop on reference systems
	
c	  Z versor: (P1-P2)x(P3-P2)/module
          vz(1)=(rccpty(j,1)-rccpty(j,2))*
     *              (rccptz(j,3)-rccptz(j,2))-
     *      (rccpty(j,3)-rccpty(j,2))*(rccptz(j,1)-rccptz(j,2))
          vz(2)=(rccptx(j,3)-rccptx(j,2))*
     *              (rccptz(j,1)-rccptz(j,2))-
     *      (rccptx(j,1)-rccptx(j,2))*(rccptz(j,3)-rccptz(j,2))
          vz(3)=(rccptx(j,1)-rccptx(j,2))*
     *               (rccpty(j,3)-rccpty(j,2))-
     *      (rccptx(j,3)-rccptx(j,2))*(rccpty(j,1)-rccpty(j,2))
          rvz=sqrt(vz(1)**2+vz(2)**2+vz(3)**2) ! Module
          vz(1)=vz(1)/rvz 
          vz(2)=vz(2)/rvz
          vz(3)=vz(3)/rvz
	  
c	  X versor: (P1-P2)
          vx(1)=rccptx(j,1)-rccptx(j,2)
          vx(2)=rccpty(j,1)-rccpty(j,2)
          vx(3)=rccptz(j,1)-rccptz(j,2)
          rvx=sqrt(vx(1)**2+vx(2)**2+vx(3)**2) ! Module
          vx(1)=vx(1)/rvx
          vx(2)=vx(2)/rvx
          vx(3)=vx(3)/rvx
	  
c	  Y versor: (Z)x(X)
          vy(1)=vz(2)*vx(3)-vx(2)*vz(3)
          vy(2)=vx(1)*vz(3)-vz(1)*vx(3)
          vy(3)=vz(1)*vx(2)-vx(1)*vz(2)
c
c         Point positions in the global reference system:
c
          center(j,1)=rccptx(j,2) ! center of ref. syst. j (P2)
          center(j,2)=rccpty(j,2)
          center(j,3)=rccptz(j,2) 
	  
          ssx(j,1)=rccptx(j,2)+vx(1) ! X versor of ref. syst. j (P2+X)
          ssx(j,2)=rccpty(j,2)+vx(2)
          ssx(j,3)=rccptz(j,2)+vx(3)
	  
          ssy(j,1)=rccptx(j,2)+vy(1) ! Y versor of ref. syst. j (P2+Y)
          ssy(j,2)=rccpty(j,2)+vy(2)
          ssy(j,3)=rccptz(j,2)+vy(3)
	  
          ssz(j,1)=rccptx(j,2)+vz(1) ! Z versor of ref. syst. j (P2+Z)
          ssz(j,2)=rccpty(j,2)+vz(2)
          ssz(j,3)=rccptz(j,2)+vz(3)

c	  Store director cosines of reference system j:
          axsys(j,1)=vx(1) ! X versor
          axsys(j,2)=vx(2)
          axsys(j,3)=vx(3)
          aysys(j,1)=vy(1) ! Y versor
          aysys(j,2)=vy(2)
          aysys(j,3)=vy(3)
          azsys(j,1)=vz(1) ! Z versor
          azsys(j,2)=vz(2)
          azsys(j,3)=vz(3)

 50     continue
 
      else 
c       Ah! Ref. system was not defined; default ref. systems 
c	are assumed:
c	Number of systems = number of paramagnetic ions
c	Center of system = each paramagnetic ion (medium position 
c						  between structures)
c	Versors = the same of given data (file 'filename2')

C	Calculate the medium centers position (on each structure)
        do 60 j=1,nfe
          do 60 k=1,nstr
            rccptx(j,2) = 0.0
            rccpty(j,2) = 0.0
            rccptz(j,2) = 0.0
 60     continue
        do 65 j=1,nfe
          do 65 k=1,nstr
            rccptx(j,2) = rccptx(j,2) + fx(k,j)/dble(nstr)
            rccpty(j,2) = rccpty(j,2) + fy(k,j)/dble(nstr)
            rccptz(j,2) = rccptz(j,2) + fz(k,j)/dble(nstr)
 65     continue

        do 70 j=1,nfe
c         Points to display:
c
          center(j,1)=rccptx(j,2)
          center(j,2)=rccpty(j,2)
          center(j,3)=rccptz(j,2)
          ssx(j,1)=rccptx(j,2)+1.
          ssx(j,2)=rccpty(j,2)
          ssx(j,3)=rccptz(j,2)
          ssy(j,1)=rccptx(j,2)
          ssy(j,2)=rccpty(j,2)+1.
          ssy(j,3)=rccptz(j,2)
          ssz(j,1)=rccptx(j,2)
          ssz(j,2)=rccpty(j,2)
          ssz(j,3)=rccptz(j,2)+1.

c	  Store director cosines of reference system j:
          axsys(j,1)=1.
          axsys(j,2)=0.
          axsys(j,3)=0.
          aysys(j,1)=0.
          aysys(j,2)=1.
          aysys(j,3)=0.
          azsys(j,1)=0.
          azsys(j,2)=0.
          azsys(j,3)=1.

 70     continue
      endif
      
* ======================================================================
c
c		*****     WRITE OBSERVED OUT FILE	****


      print *,'     *****     WRITING OBSERVED OUT FILE:'
      do 140 k=1,nstr ! For each structure:
* 	    if(DEBUGP) print *,'Structure',k ! Control print
      	ierrflag=0    ! Data mismatch error flag
        icontprot=0   ! Hit counter
        ! For each measured pseodocontact (dipolar proton) j:
        do 150 j=1,ihp
	  do 160 i=1,nat ! Compare data with each atom:
            if (num_res(k,i).eq.numres(j)) then
              if (nam_at(k,i).eq.namat(j)
     *            .or.nam_at(k,i)(2:4).eq.namat(j)(1:3)
     *            .or.nam_at(k,i)(1:3).eq.namat(j)(2:4)) then ! HIT!
                icontprot=icontprot+1
	        write(3,
     *     '(i6,1x,a4,1x,i3,1x,a3,1x,3F8.3,1x,F7.2,1x,i1,2(1x,f7.2))') 
     *          num_at(k,i),nam_at(k,i),num_res(k,i),nam_res(k,i),
     *		xp(k,i),yp(k,i),zp(k,i),obs(j),mlprot(j),
     *          tolprot(j),wprot(j)
		goto 150 ! Breaks i-loop
	      end if
            end if
 160	  continue
          if(ierrflag.eq.0) then ! Control print (data mismatch)
	    ierrflag=1
            print *,' '
	    print *,
     *      ' *** ERROR: following atoms read from file ',filename1
	    print *,'            did not match data in file ',filename2
	  end if 
          print *,namat(j),' in residue',numres(j),'(',namres(j),')'
 150	continue
 140  continue
      print *,'               done.'
 
c========================================================================
c  Data mismatch between input files: Error management
 98   continue
      if (ierrflag.eq.1) then
        print *,' *** ERROR: Data mismatch between input files'
        print *,' '
        print *,' *** PROGRAM STOP'
        print *,' '
        STOP
      endif
c========================================================================
      close(1)
      close(2)
      close(3)
      
      
      nhp=icontprot
      open (1,file=filename3,status='old')
c
c 		*****     RE-READ OBSERVED OUT FILE       *****	
c	coord. of protons		cx(k,i),cy(k,i),cz(k,i)
c       shift parameters values		obs mlprot tolprot wprot
c
      i=0
*      print *,'     *****     READING OBSERVED OUT FILE'
      do 170  k=1,nstr
      	knhp=(k-1)*nhp ! Index offset for structure k
        do 170  i=1,nhp
          read(1,'(A)',end=199) line
          read(line,'(21x,3F8.3,1x,F7.2,i2,2(1x,f7.2))')
     *              cx(k,i),cy(k,i),cz(k,i),
     *              obs(knhp+i),mlprot(knhp+i),
     *              tolprot(knhp+i),wprot(knhp+i)
 170  continue
 199  continue 
      close(1)
      open (3,file=fileout,status='new')
      write(3,'(a)') ' ' 
      write(3,'(1x,a,a)') 'Input file = ',filename3
      write(3,'(1x,a,i6)') 'Total protons ',nhp*nstr
      write(3,'(1x,a,i6)') 'Paramagnetic center',nfe
      write(3,'(a)') ' '
      write(3,'(a)') '   *****   SOLUTIONS   *****' 
c
c		*****     LAUNCH SIMPLEX	******
c
      print *,' '
      print *,' '
      print *,'               RUNNING SIMPLEX ALGORITHM:'
      print *,' '
      print *,' '
      call simplexrun()
      end
 
c==================================================================
c==================================================================

C---------------------------------------------------------------
C
C                        SIMPLEXRUN
C
C   Calculates the susceptibilty tensor of paramagnetic system
C   which minimize the experimental-to-calculated squared shift
C   difference.
C---------------------------------------------------------------
      SUBROUTINE SIMPLEXRUN()
C
      include 'sup.h'
      DIMENSION SIMP(MP,NP),EXTR(NP),VAL(MP)

c     Starting parameters for simplex algorithm
      OLDRESID=1.e+9
      RESID=2.e+9
      iviolation=2000000
      ioldvio=1000000

      delta=3.14159/dble(ngrid+1)
      
c     Simplex Algorithm:
      do 1000 l=1,ngrid
        print *,'Grid #',l,'/',ngrid,':'
            do 999 nk=1,nfe
              SIMP(1,(nk-1)*MPAR+1)=l*delta*(1.-2*rand())
              SIMP(1,(nk-1)*MPAR+2)=l*delta*(1.-2*rand())
              SIMP(1,(nk-1)*MPAR+3)=l*delta*(1.-2*rand())
              simp(1,(nk-1)*MPAR+4)=1000.*(1.-2*rand())
              simp(1,(nk-1)*MPAR+5)=1000.*(1.-2*rand())
 999        continue
          DO 10 I=1,MP
            DO 10 J=1,NP
              SIMP(I,J)=SIMP(1,J)
              if (i.eq.1) then
                SIMP(1,J)=SIMP(1,J)*0.8
              endif
              if (i.eq.j+1) then
                SIMP(I,J)=SIMP(I,J)/0.8*1.2
              endif
  10      CONTINUE
          DO 11 I=1,MP
            DO 12 J=1,NP
  12          EXTR(J)=SIMP(I,J)
            VAL(I)=CASH(EXTR)
  11      CONTINUE
          CALL SIMPCALC(SIMP,VAL)
	  
c         Minimizing optimal solution      
          if (oldresid.ge.resid) then
            OLDRESID=RESID
            ioldvio=iviolation
            do 13 i=1,nfe
              optphi(i)=phi(i)
              optteta(i)=teta(i)
              optomega(i)=omega(i)
              opta1(i)=a1dip(i)
              opta2(i)=a2dip(i)
 13         continue
          endif
        print *,' '
 1000 continue
 
c ------------> End of simplex algorithm <-----------------
c ---------------------------------------------------------

      do 998 nk=1,nfe
        EXTR((nk-1)*MPAR+1)=optphi(nk)
        EXTR((nk-1)*MPAR+2)=optteta(nk)
        EXTR((nk-1)*MPAR+3)=optomega(nk)
        EXTR((nk-1)*MPAR+4)=opta1(nk)
        EXTR((nk-1)*MPAR+5)=opta2(nk)
 998  continue
 
c     Print out the best solution
      print *,' '
      WRITE(*,FMT='(1X,A)') ' *****    BEST SOLUTION    *****'
      print *,' '
      WRITE (*,20) ioldvio, OLDRESID
      print *,' '
      WRITE(3,FMT='(1X,A)') ' '
      WRITE(3,FMT='(1X,A)') ' *****    BEST SOLUTION    *****'
      WRITE(3,FMT='(1X,A)') ' '
      do 15 i=1,nfe
      
        WRITE (*,21) optphi(i),optteta(i),optomega(i)
        WRITE (*,22) opta1(i),opta2(i)
        print *,' '
      
        WRITE (3,21) optphi(i),optteta(i),optomega(i)
        WRITE (3,22) opta1(i),opta2(i)
 15   continue

c     Rotation matrix:
c      if (intsys.eq.1) then
        do 777 i=1,nfe
          P=optphi(i)
          T=optteta(i)
          O=optomega(i)
          asxx(i,1)=cos(P)*cos(O)
          asxx(i,2)=sin(P)*cos(O)
          asxx(i,3)=sin(O)
          asyy(i,1)=-cos(T)*sin(P)-sin(O)*cos(P)*sin(T)
          asyy(i,2)=cos(T)*cos(P)-sin(O)*sin(P)*sin(T)
          asyy(i,3)=sin(T)*cos(O)
          aszz(i,1)=sin(T)*sin(P)-sin(O)*cos(P)*cos(T)
          aszz(i,2)=-sin(T)*cos(P)-sin(O)*sin(P)*cos(T)
          aszz(i,3)=cos(T)*cos(O)
 777    continue
c      endif

c     COMPUTING THE THREE PERMUTATIONS OF TENSOR AXIS:
      CALL PRRES

C     FINAL DATA PRINTING ON OUTPUT FILE:
      ERR=CASH(EXTR)
      call DISPSHIFT
      
 20   FORMAT (1X,' VIOLATIONS ',I6,1X,' RESIDUALS ',F10.3)
 21   FORMAT(1X,'PHI= ',F9.3,1X,'THETA= ',F9.3,1X,'OMEGA= ',F9.3)
 22   FORMAT (1X,'A1= ',F9.3,1X,'A2= ',F9.3)

      RETURN
      END
 
c==================================================================
c==================================================================


C-------------------------------------------------------------
C 
C                        SIMPCALC
C
C   Simplex algorithm kernel
C--------------------------------------------------------------
      SUBROUTINE SIMPCALC(P,Y)
C
      include 'sup.h'
      PARAMETER (AL=1.,BE=.5,GA=2.)
      DIMENSION P(MP,NP),Y(MP),PR(MP),PT(MP),PE(MP),
     *          PC(MP),PK(MP),PBAR(MP),OLDP(NP),EXTR(MP)
c                               inizializzazione dell'algoritmo
      ITER=0
      RNDOM=0
 1    perm=1
      do 1000 while (perm.gt.0)
        perm=0
        do 1100 i=1,np
          if (Y(i).gt.Y(i+1)) then
            app=Y(i+1)
            Y(i+1)=Y(i)
            Y(i)=app
            do 1200 j=1,np
              app=P(i+1,j)
              P(i+1,j)=P(i,j)
              P(i,j)=app
 1200       continue
            perm=perm+1
          endif
 1100   continue  
 1000 continue
      AVERAGE=0
      ASTDEV=0
      DO 100 I=1,MP
 100    AVERAGE=AVERAGE+Y(I)
      AVERAGE=AVERAGE/MP
      DO 101 I=1, MP
 101    ASTDEV=ASTDEV+(Y(I)-AVERAGE)**2.
      ERR=DSQRT(ASTDEV/MP)
      RESID=Y(1)  
      IF (ITER.gt.39999.OR.ERR.lt.1.D-9) THEN
        WRITE (*,200) ITER,IVIOLATION,RESID
        WRITE (3,200) ITER,IVIOLATION,RESID
        WRITE (*,201) 
     *     (P(1,(i-1)*MPAR+1),
     *      P(1,(i-1)*MPAR+2),
     *      P(1,(i-1)*MPAR+3),
     *      P(1,(i-1)*MPAR+4),P(1,(i-1)*MPAR+5),i=1,nfe)
        WRITE (3,201) 
     *     (P(1,(i-1)*MPAR+1),
     *      P(1,(i-1)*MPAR+2),
     *      P(1,(i-1)*MPAR+3),
     *      P(1,(i-1)*MPAR+4),P(1,(i-1)*MPAR+5),i=1,nfe)

 200    FORMAT (1X,'ITERS ',I6,1X,' VIOLATIONS ',I6,
     *          1X,'RESIDUALS= ',F10.3)
 201    FORMAT(1X,'PHI= ',F9.3,1X,'THETA= ',F9.3,1X,'OMEGA= ',F9.3,
     *         1X,'A1= ',F9.3,1X,'A2= ',F9.3)

      do 997 nk=1,nfe
        phi(nk)=P(1,(nk-1)*MPAR+1)
        teta(nk)=P(1,(nk-1)*MPAR+2)
        omega(nk)=P(1,(nk-1)*MPAR+3)
        a1dip(nk)=P(1,(nk-1)*MPAR+4)
        a2dip(nk)=P(1,(nk-1)*MPAR+5)
 997        continue
        RETURN
      ENDIF
C--------------------------------------------------------
      ITER=ITER+1
      DO 12 J=1,NP
 12     PBAR(J)=0
      DO 13 I=1,NP
        DO 13 J=1,NP
 13         PBAR(J)=PBAR(J)+P(I,J)
      DO 14 J=1,NP
        PBAR(J)=PBAR(J)/NP
        PR(J)=(1.+AL)*PBAR(J)-AL*P(MP,J)
 14   CONTINUE
      YR=CASH(PR)
      DO 15 J=1,NP
        PK(J)=PR(J)
 15   CONTINUE
      YK=YR
      IF(YR.lt.Y(NP)) THEN
        IF(YR.lt.Y(1)) THEN
          DO 16 J=1,NP
            PE(J)=GA*PR(J)+(1.-GA)*PBAR(J)
 16       CONTINUE
          YE=CASH(PE)
          iter=iter+1
          IF (YE.LT.Y(1)) THEN
            DO 17 J=1,NP
              PK(J)=PE(J)
 17         CONTINUE
            YK=YE
c            write(*,FMT='(BN,A)') 'expand'
          ENDIF
        ENDIF
      ELSE 
        DO 18 j=1,NP
          PT(j)=P(MP,j)
 18     CONTINUE
        YT=Y(MP)
        IF (YR.LT.YT) THEN
          DO 19 J=1,NP
            PT(J)=PR(J)
 19       CONTINUE
          YT=YR
        ENDIF
        DO 20 J=1,NP
          PC(J)=BE*PT(J)+(1.-BE)*PBAR(J)
 20     CONTINUE
        YC=CASH(PC)
        iter=iter+1
        IF (YC.lt.Y(NP)) THEN
          DO 21 J=1,NP
            PK(J)=PC(J)
 21       CONTINUE
          YK=YC
c          write(*,FMT='(BN,A)') 'contract'
        ELSE
          do 22 i=1,NP
            OLDP(i)=P(2,i)
 22      continue
          DO 24 I=2,NP 
            DO 23 J=1,NP
              P(I,J)=.5*(P(1,J)+P(I,J))
 23         CONTINUE
            PR(I)=P(I,I)  
            Y(I)=CASH(PR)
 24       CONTINUE
          iter=iter+NP-1
          do 25 i=1,NP
            OLDP(i)=abs(OLDP(i)-P(2,i))
 25       continue
          oldpoppa=0 
          do 251 i=1,NP-1
             oldpoppa=oldpoppa+oldp(i)
 251      continue
          if (abs(oldpoppa).lt.1.D-8.and.RNDOM.lt.10) then
c            iter=0
            RNDOM=RNDOM+1
            write(*,FMT='(BN,A,F16.6)') 'random  ',Y(1)
            do 2001 il=1,MP
            do 2002 im=1,NP
              rmr=1+RAND()
              P(il,im)=P(il,im)*rmr
 2002       continue
 2001       continue   
          DO 2011 Il=1,MP
            DO 2012 Jm=1,NP
 2012          EXTR(Jm)=P(Il,Jm)
            Y(Il)=CASH(EXTR)
 2011      CONTINUE
            goto 1 
          endif 
          if (abs(oldpoppa).lt.1.D-9) then
            iter=50000
          endif
          DO 26 J=1,NP
            PK(j)=0.5*(P(1,j)+P(mp,j))
 26       CONTINUE
          YK=CASH(PK)
          iter=iter+1
c          write(*,FMT='(BN,A)') 'shrink' 
        ENDIF
      ENDIF
      DO 27 J=1,NP
        P(MP,J)=PK(J)
 27   CONTINUE
      Y(MP)=YK 
      GO TO 1
      END

 
c==================================================================
c==================================================================

C-------------------------------------------------------------
C
C                       CASH 
C
C--------------------------------------------------------------
      FUNCTION CASH(vett)
C
      include 'sup.h'
      dimension vett(NP)
c     Compute error and violations
      IVIOLATION=0
      TMP1=0
      i=1
      do while (i.le.nhp*nstr)
        shift(i)=0.0
        i=i+1
      enddo
      do 1 n=1,nstr
      do 2 m=1,nfe
        P=vett((m-1)*MPAR+1)
        T=vett((m-1)*MPAR+2)
        O=vett((m-1)*MPAR+3)
        A1D=vett((m-1)*MPAR+4)
        A2D=vett((m-1)*MPAR+5)
        axx=cos(P)*cos(O)
        axy=sin(P)*cos(O)
        axz=sin(O)
        ayx=-cos(T)*sin(P)-sin(O)*cos(P)*sin(T)
        ayy=cos(T)*cos(P)-sin(O)*sin(P)*sin(T)
        ayz=sin(T)*cos(O)
        azx=sin(T)*sin(P)-sin(O)*cos(P)*cos(T)
        azy=-sin(T)*cos(P)-sin(O)*sin(P)*cos(T)
        azz=cos(T)*cos(O)
        tmp2=0.0
        i=1
        ihp=(n-1)*nhp
        do 10 while (i.le.nhp)
          nprotml=mlprot(ihp+i)
          if (nprotml.gt.1) then
            costmult=1/dble(nprotml)
            app=0.0
            do 11 L=1,nprotml
              xapp=CX(n,I)-fx(n,m)
              yapp=CY(n,I)-fy(n,m)
              zapp=CZ(n,I)-fz(n,m)
              r=sqrt(xapp**2+yapp**2+zapp**2)
              scalz=(xapp*azx+yapp*azy+zapp*azz)
              scalx=(xapp*axx+yapp*axy+zapp*axz)
              scaly=(xapp*ayx+yapp*ayy+zapp*ayz)
              g1=(sqrt(3.)*scalz-r)*(sqrt(3.)*scalz+r)
              g2=(scalx-scaly)*(scalx+scaly)
              app=app+(A1D*G1+1.5*A2D*G2)*costmult/r**5
              I=I+1
 11         continue
            do 12 L=1,nprotml
              shift(ihp+i-l)=shift(ihp+i-l)+app
 12         continue
          else
            xapp=CX(n,I)-fx(n,m)
            yapp=CY(n,I)-fy(n,m)
            zapp=CZ(n,I)-fz(n,m)
            r=sqrt(xapp**2+yapp**2+zapp**2)
            scalz=(xapp*azx+yapp*azy+zapp*azz)
            scalx=(xapp*axx+yapp*axy+zapp*axz)
            scaly=(xapp*ayx+yapp*ayy+zapp*ayz)
            g1=(sqrt(3.)*scalz-r)*(sqrt(3.)*scalz+r)
            g2=(scalx-scaly)*(scalx+scaly)
            shift(ihp+i)=shift(ihp+i)+
     *                      (A1D*G1+1.5*A2D*G2)/r**5
            I=I+1
          endif
 10     continue
 2    continue
      do 3 i=1,nhp
          tmp2=abs(shift(ihp+i)-obs(ihp+i))
     *                    -tolprot(ihp+i)
          if (tmp2.gt.0.0) then
            IVIOLATION=IVIOLATION+1
            TMP1=tmp1+tmp2**2*wprot(ihp+i)
     *               /dble(mlprot(ihp+i))
          endif
 3    continue 
 1    continue

      CASH=tmp1
      RETURN
      END

 
c==================================================================
c==================================================================

C--------------------------------------------------------------
C
C                       DISPSHIFT
C
C This routine only prints all final data on the output file
C--------------------------------------------------------------
      SUBROUTINE DISPSHIFT
C
      include 'sup.h'
      character filestr*21, ss*40
      dimension VV(MAXSTR,MAXOS),VMAX(maxos),VMIN(maxos)
      dimension nnmax(maxos),nnmin(maxos)
c                            visualizzazione delle violazioni
      open (1,file=filename3,status='old')
      write(3,'(a)') ' '
      write(3,'(a)')
     * '**********************************************************'
      write(3,'(1x,a)')
     * '    PROTONS            MIN ERR. MAX ERR. OBS.  '
      
      rm=0
      tmp1=0
      smed=0
      stdev=0
      do 2 i=1, nstr
       do 2 j=1, nhp
        tmp2=abs(obs(j+(i-1)*nhp)-SHIFT(j+(i-1)*nhp))
     *       -tolprot(j+(i-1)*nhp)
        if (tmp2.gt.0.0) then
           if (SHIFT(j+(i-1)*nhp).gt.obs(j+(i-1)*nhp)) then
             VV(i,j)= tmp2
           else
             VV(i,j)= -tmp2
           endif
        else
         VV(i,j)= 0.
        endif
 2    continue
      do 44 i=1, nhp
       VMAX(i)=0.
       VMIN(i)=1000.
 44   continue
      do 4 j=1, nhp
       do 4 i=1, nstr
        if (abs(VV(i,j)).gt.abs(VMAX(j))) then 
         VMAX(j)=VV(i,j)
c        nnmax=i
        endif
        if ((abs(VV(i,j)).lt.abs(VMIN(j))).and.(VV(i,j).ne.0.)) then
         VMIN(j)=VV(i,j) 
c        nnmin=i
        endif
 4    continue   
      do 7 j=1, nhp
       do 8 i=1, nstr
        if (abs(VV(i,j)).gt.0.) then
         ss(i:i+1)="*"
c        if (i.eq.nnmax) then
c         ss(i:i+1)="+"
c        endif
c        if (i.eq.nnmin) then
c         ss(i:i+1)="-"
c        endif
        else
         ss(i:i+1)=" "
        endif
 8     continue
       if (VMIN(j).eq.1000.) then
         VMIN(j)=0.
       endif
       read(1,'(A)') filestr
       write(3,'(1X,A,1X,F7.3,1X,F7.3,1X,F7.3,1x,A)') filestr,
     *      VMIN(j), VMAX(j), obs(J), ss 
 7    continue    

      write(3,'(1x)') 
      write(3,'(1x)') 
      write(3,'(1x)') 
      write(3,'(1x,a)')
     * '    PROTONS              MIS.    CALC.    ERR.^2'


      rewind(1)
      do 3 j=1,nhp*nstr
        tmp2=abs(obs(j)-SHIFT(J))-tolprot(J)
        read(1,'(A)') filestr
        if (tmp2.gt.0.0) then
          IVIOLATION=IVIOLATION+1
          rm=rm+1
          smed=smed+tmp2**2*wprot(j)/mlprot(j)
          tmp1=tmp2**2*wprot(j)/mlprot(j)
          write(3,'(1X,A,1X,F7.3,1X,F7.3,1X,F10.5)') filestr,
     *               obs(J),shift(j),tmp1
 15       FORMAT(1X,A,1X,F7.3,1X,F7.3,1X,F7.3)
        else
        write(3,'(1X,A,1X,F7.3,1X,F7.3)') filestr,
     *               obs(J),shift(j)
        endif
 3    continue
      smed=smed/rm
      do 41 j=1,nhp*nstr
        tmp2=abs(obs(j)-SHIFT(J))-tolprot(J)
        if (tmp2.gt.0.0) then
          stdev=stdev+(tmp2**2-smed)**2*wprot(j)/mlprot(j)
        endif
 41   continue
      stdev=stdev/rm
      print *,' '
      WRITE(3,'(A)') ' ' 
      write(*,11) 'mean value of error =',smed
 11   format(1X,A,1X,F10.6)
       write(3,12) 'mean value of error =',smed
 12   format(1X,A,1X,F10.6)
      write(*,13) 'Standard deviation =',stdev
 13   format(1X,A,1X,F10.6)    
       write(3,14) 'Standard deviation =',stdev
 14   format(1X,A,1X,F10.6)         
      close(1)
      close(3)
      RETURN
      END
C     ----------------------------------

 
c==================================================================
c==================================================================

C--------------------------------------------------------------
C
C                      PRRES
C
C  This routine performs the two final permutations of tensor 
c  axis and computes their parameters
C--------------------------------------------------------------
      SUBROUTINE PRRES
C
      include 'sup.h'
c
c	saxx sayy sazz : versors (rotation matrix columns) of tensor 
c
c	saxsys saysys sazsys : versors (r.m.c.) of molecular ref. syst.
c
      dimension saxx(3),sayy(3),sazz(3), ! Rotation matrix of tensor
     *   appx(3),appy(3),appz(3), ! Auxiliary variables for permutations
     *   saxsys(3),saysys(3),sazsys(3) ! Rot. matrix of reference system 
      character namatom*4,namress*3
      PI=3.1415927
c
c	open pdb file
c
      open (2,file=filename2,status='old')
      i=0
      do while (i.le.nstr*nat)
        i=i+1
        read(2,'(A)',END=99) line
        read(line,'(5x,i6,1x,a4,1x,a3,2x,i4,4x,3f8.3)')
     *          numberatom,namatom,namress,nores,xa,ya,za
c               ??? only 'numberatom' & 'nores' are really used?
      enddo
 99   i=0
 
c     For each paramagnetic ion n:
      do 10 n=1,nfe 
      A1=opta1(n)  ! Axial anisotropy
      A2=opta2(n)  ! Rhombic anisotropy
      if (nsystem.gt.1) then
        do 1 i=1,3 ! Rotation matrix of reference system n:
          saxsys(i)=axsys(n,i) ! Components of x-versor (Ex)
          saysys(i)=aysys(n,i) ! Components of y-versor (Ey)
          sazsys(i)=azsys(n,i) ! Components of z-versor (Ez)
 1      continue
      else
        do 2 i=1,3 ! Rotation matrix of the unique reference system:
          saxsys(i)=axsys(1,i) ! Components of x-versor (Ex)
          saysys(i)=aysys(1,i) ! Components of y-versor (Ey)
          sazsys(i)=azsys(1,i) ! Components of z-versor (Ez)
 2      continue
      endif
      do 3 i=1,3 ! Rotation matrix of tensor n:
        saxx(i)=asxx(n,i) ! Components of x-versor (rx)
        sayy(i)=asyy(n,i) ! Components of y-versor (ry)
        sazz(i)=aszz(n,i) ! Components of z-versor (rz)
 3    continue

      print *,
     *'***********************************************************'
      WRITE(3,'(A)') ' '
      WRITE(3,'(A)')
     * '**********************************************************'
      WRITE(3,'(A)') ' '
      WRITE(*,'(A,i3)') ' PARAMAGNETIC CENTER NUMBER ',n

      if(DEBUGP) then ! DEBUG PRINT
        print *,' '
        print *,'Tensor Rotation Matrix: '
        print '(3(2x,f8.4))',saxx(1),sayy(1),sazz(1)
        print '(3(2x,f8.4))',saxx(2),sayy(2),sazz(2)
        print '(3(2x,f8.4))',saxx(3),sayy(3),sazz(3)
        print *,' '
      endif 
      
C     Control: is the tensor really right-hand oriented?
c	[??? But could it be otherwise!?]
c     Right-hand versor rz = (rx)x(ry)
      rzx=saxx(2)*sayy(3)-saxx(3)*sayy(2)
      rzy=saxx(3)*sayy(1)-saxx(1)*sayy(3)
      rzz=saxx(1)*sayy(2)-saxx(2)*sayy(1)
      scalar=sazz(1)*rzx+sazz(2)*rzy+sazz(3)*rzz
      if (scalar.lt.0.) then ! Versor z is left-hand!
c	Switch rx & ry
        do 4 i=1,3
          appx(i)=saxx(i)
          appy(i)=sayy(i)
  4     continue
        do 5 i=1,3
          saxx(i)=appy(i)
          sayy(i)=appx(i)
  5     continue
        A2=-A2 ! Rhombic anisotropy
      endif

      numberatom=numberatom+1
      namatom=' CEN'
      namress='PMC'
      nores=nores+n
      WRITE(3,'(A,i3)') ' PDB COORD. OF PARAMAGNETIC CENTER NUMBER ',n
      WRITE(3,'(A4,1x,i6,1x,a4,1x,a3,2x,i4,4x,3f8.3)')
     * 'ATOM',numberatom,namatom,namress,nores,
     * rccptx(n,2),rccpty(n,2),rccptz(n,2)

c ================ NO PERMUTATION OF AXIS ==================

      WRITE(*,*)
      WRITE(*,'(A)')
     * ' ANGLES BETWEEN REFERENCE AXIS AND TENSOR AXIS (degrees) '
      WRITE(3,'(A)') ' '
      WRITE(3,'(A)') 
     * ' ANGLES BETWEEN REFERENCE AXIS AND TENSOR AXIS (degrees) '
c	(Ex)(rx) Scalar product
      cosa=saxx(1)*saxsys(1)+saxx(2)*saxsys(2)+
     *         saxx(3)*saxsys(3)
c	(Ey)(ry) Scalar product
      cosb=sayy(1)*saysys(1)+sayy(2)*saysys(2)+
     *         sayy(3)*saysys(3)
c	(Ez)(rz) Scalar product
      cosg=sazz(1)*sazsys(1)+sazz(2)*sazsys(2)+
     *         sazz(3)*sazsys(3)
      write (*,23) 
     *acos(cosa)*180./PI,acos(cosb)*180./PI,acos(cosg)*180./PI
      write (3,23) 
     *acos(cosa)*180./PI,acos(cosb)*180./PI,acos(cosg)*180./PI

      WRITE(*,'(A)') ' FIVE TENSOR PARAMETERS '
      WRITE(3,'(A)') ' '
      WRITE(3,'(A)') ' FIVE TENSOR PARAMETERS '
      
      CALL PANG(0,n) ! Compute angles of Euler
      
      WRITE (*,22) 
     * A1,'(',A1*12*PI/1.e+4,' m^3 e+32)  ',
     * A2,'(',A2*12*PI/1.e+4,' m^3 e+32)  '
      WRITE (3,22) 
     * A1,'(',A1*12*PI/1.e+4,' m^3 e+32)  ',
     * A2,'(',A2*12*PI/1.e+4,' m^3 e+32)  '
c      WRITE(3,'(A)') ' '
c      WRITE(3,'(A,i4)') ' ATOMS FOR TENSOR  ',n
c
c      vectors to display:
c
c     Vertex of tensor x-versor in global reference system
      tsx(n,1)=rccptx(n,2)+saxx(1) 
      tsx(n,2)=rccpty(n,2)+saxx(2)
      tsx(n,3)=rccptz(n,2)+saxx(3)
c     Vertex of tensor y-versor in global reference system
      tsy(n,1)=rccptx(n,2)+sayy(1)
      tsy(n,2)=rccpty(n,2)+sayy(2)
      tsy(n,3)=rccptz(n,2)+sayy(3)
c     Vertex of tensor z-versor in global reference system
      tsz(n,1)=rccptx(n,2)+sazz(1)
      tsz(n,2)=rccpty(n,2)+sazz(2)
      tsz(n,3)=rccptz(n,2)+sazz(3)


c    Appending atoms as pdb-records to output file
c
      numberatom=numberatom+1
      namatom=' TEX'
      WRITE(3,'(A)') ' '
      WRITE(3,'(A)') ' PDB COORDINATES OF TENSOR SYSTEM'
      WRITE(3,'(A4,1x,i6,1x,a4,1x,a3,2x,i4,4x,3f8.3)')
     * 'ATOM',numberatom,namatom,namress,nores,
     * tsx(n,1),tsx(n,2),tsx(n,3)
      numberatom=numberatom+1
      namatom=' TEY'
      WRITE(3,'(A4,1x,i6,1x,a4,1x,a3,2x,i4,4x,3f8.3)')
     * 'ATOM',numberatom,namatom,namress,nores,
     * tsy(n,1),tsy(n,2),tsy(n,3)
      numberatom=numberatom+1
      namatom=' TEZ'
      WRITE(3,'(A4,1x,i6,1x,a4,1x,a3,2x,i4,4x,3f8.3)')
     * 'ATOM',numberatom,namatom,namress,nores,
     * tsz(n,1),tsz(n,2),tsz(n,3)
      numberatom=numberatom-3

c ====== FIRST PERMUTATION OF AXIS (x become z; y <- x and z <- y) ========

      WRITE(*,*)
      WRITE(*,'(A,A)')
     * ' ANGLES BETWEEN REFERENCE AXIS AND TENSOR AXIS ' ,
     * '(degrees) 1 PERMUTATION'
      WRITE(3,'(A)') ' '
      WRITE(3,'(A,A)')
     * ' ANGLES BETWEEN REFERENCE AXIS AND TENSOR AXIS ' ,
     * '(degrees) 1 PERMUTATION'
     
      do 6 i=1,3 ! Temporary copy of tensor rotation matrix
        appx(i)=saxx(i) 
        appy(i)=sayy(i)
        appz(i)=sazz(i)
 6     continue 
      do 7 i=1,3 ! Permutation of tensor rotation matrix
        saxx(i)=appz(i) ! New x was z
        sayy(i)=appx(i) ! New y was x
        sazz(i)=appy(i) ! New z was y
 7     continue
 
c	(Ex)(rx) Scalar product
      cosa=saxx(1)*saxsys(1)+saxx(2)*saxsys(2)+
     *         saxx(3)*saxsys(3)
c	(Ey)(ry) Scalar product
      cosb=sayy(1)*saysys(1)+sayy(2)*saysys(2)+
     *         sayy(3)*saysys(3)
c	(Ez)(rz) Scalar product
      cosg=sazz(1)*sazsys(1)+sazz(2)*sazsys(2)+
     *         sazz(3)*sazsys(3)
     
      write (*,23)
     *acos(cosa)*180./PI,acos(cosb)*180./PI,acos(cosg)*180./PI
      write (3,23)
     *acos(cosa)*180./PI,acos(cosb)*180./PI,acos(cosg)*180./PI

      WRITE(*,'(A)') ' FIVE TENSOR PARAMETERS '
      WRITE(3,'(A)') ' '
      WRITE(3,'(A)') ' FIVE TENSOR PARAMETERS '
      
      CALL PANG(1,n)
      
C     Computing new anisotropies for permuted tensor, 
c     starting from the formula:
c	A1=CHIzz-1/2(CHIxx+CHIyy)   axial anisotropy
c	A2=CHIxx-CHIyy              rhombic anisotropy
c
c     Switching x with z, y with x, z with y, we obtain the
c     two new anisotropies:
      b1=-(A1+1.5*A2)/2.0 ! New axial anisotropy
      b2=A1-A2/2.0        ! New rhombic anisotropy
      
      WRITE (*,22) 
     * b1,'(',b1*12*PI/1.e+4,' m^3 e+32)  ',
     * b2,'(',b2*12*PI/1.e+4,' m^3 e+32)  '
      WRITE (3,22) 
     * b1,'(',b1*12*PI/1.e+4,' m^3 e+32)  ',
     * b2,'(',b2*12*PI/1.e+4,' m^3 e+32)  '
c      WRITE(3,'(A,i4)') ' ATOMS FOR TENSOR  ',n

c
c      Vectors to display:
c
c     Vertex of tensor x-versor in global reference system
      tsx(n,1)=rccptx(n,2)+saxx(1)
      tsx(n,2)=rccpty(n,2)+saxx(2)
      tsx(n,3)=rccptz(n,2)+saxx(3)
c     Vertex of tensor y-versor in global reference system
      tsy(n,1)=rccptx(n,2)+sayy(1)
      tsy(n,2)=rccpty(n,2)+sayy(2)
      tsy(n,3)=rccptz(n,2)+sayy(3)
c     Vertex of tensor z-versor in global reference system
      tsz(n,1)=rccptx(n,2)+sazz(1)
      tsz(n,2)=rccpty(n,2)+sazz(2)
      tsz(n,3)=rccptz(n,2)+sazz(3)

c    Appending atoms as pdb-records to output file
c
      numberatom=numberatom+1
      namatom=' TEX'
      WRITE(3,'(A)') ' '
      WRITE(3,'(A)') ' PDB COORDINATES OF TENSOR SYSTEM'
      WRITE(3,'(A4,1x,i6,1x,a4,1x,a3,2x,i4,4x,3f8.3)')
     * 'ATOM',numberatom,namatom,namress,nores,
     * tsx(n,1),tsx(n,2),tsx(n,3)
      numberatom=numberatom+1
      namatom=' TEY'
      WRITE(3,'(A4,1x,i6,1x,a4,1x,a3,2x,i4,4x,3f8.3)')
     * 'ATOM',numberatom,namatom,namress,nores,
     * tsy(n,1),tsy(n,2),tsy(n,3)
      numberatom=numberatom+1
      namatom=' TEZ'
      WRITE(3,'(A4,1x,i6,1x,a4,1x,a3,2x,i4,4x,3f8.3)')
     * 'ATOM',numberatom,namatom,namress,nores,
     * tsz(n,1),tsz(n,2),tsz(n,3)
      numberatom=numberatom-3

c ====== SECOND PERMUTATION OF AXIS (x become z; y <- x and z <- y) ========

      WRITE(*,*)
      WRITE(*,'(A,A)')
     * ' ANGLES BETWEEN REFERENCE AXIS AND TENSOR AXIS ' ,
     * '(degrees) 2 PERMUTATION'
      WRITE(3,'(A)') ' '
      WRITE(3,'(A,A)')
     * ' ANGLES BETWEEN REFERENCE AXIS AND TENSOR AXIS ',
     * '(degrees) 2 PERMUTATION'
     
      do 8 i=1,3 ! Temporary copy of tensor rotation matrix
        appx(i)=saxx(i)
        appy(i)=sayy(i)
        appz(i)=sazz(i)
 8     continue
      do 9 i=1,3 ! Permutation of tensor rotation matrix
        saxx(i)=appz(i) ! New x was z (and initially was y)
        sayy(i)=appx(i) ! New y was x (and initially was z)
        sazz(i)=appy(i) ! New z was y (and initially was x)
 9     continue
 
c	(Ex)(rx) Scalar product
      cosa=saxx(1)*saxsys(1)+saxx(2)*saxsys(2)+
     *         saxx(3)*saxsys(3)
c	(Ey)(ry) Scalar product
      cosb=sayy(1)*saysys(1)+sayy(2)*saysys(2)+
     *         sayy(3)*saysys(3)
c	(Ez)(rz) Scalar product
      cosg=sazz(1)*sazsys(1)+sazz(2)*sazsys(2)+
     *         sazz(3)*sazsys(3)
     
      write (*,23)
     *acos(cosa)*180./PI,acos(cosb)*180./PI,acos(cosg)*180./PI
      write (3,23)
     *acos(cosa)*180./PI,acos(cosb)*180./PI,acos(cosg)*180./PI

      WRITE(*,'(A)') ' FIVE TENSOR PARAMETERS '
      WRITE(3,'(A)') ' '
      WRITE(3,'(A)') ' FIVE TENSOR PARAMETERS '
      
      CALL PANG(2,n)
      
C     Computing new anisotropies for permuted tensor
c     (see note above)
      b1=-(A1-1.5*A2)/2.0 ! New axial anisotropy
      b2=-(A1+A2/2.0)     ! New rhombic anisotropy
      
      WRITE (*,22) 
     * b1,'(',b1*12*PI/1.e+4,' m^3 e+32)  ',
     * b2,'(',b2*12*PI/1.e+4,' m^3 e+32)  '
      WRITE (3,22)
     * b1,'(',b1*12*PI/1.e+4,' m^3 e+32)  ',
     * b2,'(',b2*12*PI/1.e+4,' m^3 e+32)  '
c
c
c      Vectors to display:
c
c     Vertex of tensor x-versor in global reference system
      tsx(n,1)=rccptx(n,2)+saxx(1)
      tsx(n,2)=rccpty(n,2)+saxx(2)
      tsx(n,3)=rccptz(n,2)+saxx(3)
c     Vertex of tensor y-versor in global reference system
      tsy(n,1)=rccptx(n,2)+sayy(1)
      tsy(n,2)=rccpty(n,2)+sayy(2)
      tsy(n,3)=rccptz(n,2)+sayy(3)
c     Vertex of tensor z-versor in global reference system
      tsz(n,1)=rccptx(n,2)+sazz(1)
      tsz(n,2)=rccpty(n,2)+sazz(2)
      tsz(n,3)=rccptz(n,2)+sazz(3)

c    Appending atoms as pdb-records to output file
c
      numberatom=numberatom+1
      namatom=' TEX'
      WRITE(3,'(A)') ' '
      WRITE(3,'(A)') ' PDB COORDINATES OF TENSOR SYSTEM'
      WRITE(3,'(A4,1x,i6,1x,a4,1x,a3,2x,i4,4x,3f8.3)')
     * 'ATOM',numberatom,namatom,namress,nores,
     * tsx(n,1),tsx(n,2),tsx(n,3)
      numberatom=numberatom+1
      namatom=' TEY'
      WRITE(3,'(A4,1x,i6,1x,a4,1x,a3,2x,i4,4x,3f8.3)')
     * 'ATOM',numberatom,namatom,namress,nores,
     * tsy(n,1),tsy(n,2),tsy(n,3)
      numberatom=numberatom+1
      namatom=' TEZ'
      WRITE(3,'(A4,1x,i6,1x,a4,1x,a3,2x,i4,4x,3f8.3)')
     * 'ATOM',numberatom,namatom,namress,nores,
     * tsz(n,1),tsz(n,2),tsz(n,3)

c  = PRINTING VERSORS OF MOLECULAR REFERENCE SYSTEM(S) ON OUTPUT FILE =

      WRITE(3,'(A)') ' '

c     Vertex of molecular system x-versor in global reference system
      saxsys(1)=rccptx(n,2)+saxsys(1)
      saxsys(2)=rccpty(n,2)+saxsys(2)
      saxsys(3)=rccptz(n,2)+saxsys(3)

c     Vertex of molecular system y-versor in global reference system
      saysys(1)=rccptx(n,2)+saysys(1)
      saysys(2)=rccpty(n,2)+saysys(2)
      saysys(3)=rccptz(n,2)+saysys(3)

c     Vertex of molecular system z-versor in global reference system
      sazsys(1)=rccptx(n,2)+sazsys(1)
      sazsys(2)=rccpty(n,2)+sazsys(2)
      sazsys(3)=rccptz(n,2)+sazsys(3)

c    Appending atoms as pdb-records to output file
c
      numberatom=numberatom+1
      namatom=' ISX'
      WRITE(3,'(A)') ' '
      WRITE(3,'(A)') ' PDB COORDINATES OF REFERENCE SYSTEM'
      WRITE(3,'(A4,1x,i6,1x,a4,1x,a3,2x,i4,4x,3f8.3)')
     * 'ATOM',numberatom,namatom,namress,nores,
     * saxsys(1),saxsys(2),saxsys(3)
      numberatom=numberatom+1
      namatom=' ISY'
      WRITE(3,'(A4,1x,i6,1x,a4,1x,a3,2x,i4,4x,3f8.3)')
     * 'ATOM',numberatom,namatom,namress,nores,
     * saysys(1),saysys(2),saysys(3)
      numberatom=numberatom+1
      namatom=' ISZ'
      WRITE(3,'(A4,1x,i6,1x,a4,1x,a3,2x,i4,4x,3f8.3)')
     * 'ATOM',numberatom,namatom,namress,nores,
     * sazsys(1),sazsys(2),sazsys(3)

 10   continue

      close(2)
  
 22   FORMAT
     * (1X,'A1= ',F9.3,2x,a1,F7.3,a,2X,'A2= ',F9.3,2x,a1,F7.3,a)
 23   FORMAT (1X,'arccos x^xt = ',F8.3,4X,'arccos y^yt = ',F8.3,
     *        4X,'arccos z^zt = ',F8.3)

      END


C--------------------------------------------------------------
C
C                      PANG
C
C   This subroutine computes the three new angles of Euler for
C   each permutation of tensor axis
c   [REWRITTEN - October 1999]
C--------------------------------------------------------------
      SUBROUTINE PANG(n,nr)
C
c	n = Number of axis permutation:
c	    0 = none; 1 = first permutation; else = second perm.
c	nr = Number (label) of reference system
c

      include 'sup.h'
      dimension temp(3) ! Auxiliary variable for permutations
      
      PI=3.1415927 
      
c     Composing permuted rotation matrix:
c 	Note: asxx(i,j) = 1st column of rot. matrix of tensor i
c 	      asyy(i,j) = 2nd column of rot. matrix of tensor i
c 	      aszz(i,j) = 3rd column of rot. matrix of tensor i
c 	                  (j=1,2,3; row index)
c-----------------------------------------------
      if (n.eq.0) then ! case 0: No permutation
c       x is x
        sxx=asxx(nr,1)
        sxy=asxx(nr,2)
        sxz=asxx(nr,3)
c       y is y
        syx=asyy(nr,1)
        syy=asyy(nr,2)
        syz=asyy(nr,3)
c       z is z
        szx=aszz(nr,1)
        szy=aszz(nr,2)
        szz=aszz(nr,3)
     
      elseif (n.eq.1) then ! case 1: First permutation of axis
c       x is z
        sxx=aszz(nr,1)
        sxy=aszz(nr,2)
        sxz=aszz(nr,3)
c       y is x
        syx=asxx(nr,1)
        syy=asxx(nr,2)
        syz=asxx(nr,3)
c       z is y
        szx=asyy(nr,1)
        szy=asyy(nr,2)
        szz=asyy(nr,3)
     
      else  ! case default: Second permutation of axis
c       x is y
        sxx=asyy(nr,1)
        sxy=asyy(nr,2)
        sxz=asyy(nr,3)
c       y is z
        syx=aszz(nr,1)
        syy=aszz(nr,2)
        syz=aszz(nr,3)
c       z is x
        szx=asxx(nr,1)
        szy=asxx(nr,2)
        szz=asxx(nr,3)

      endif

c  Computing angles of Euler:
c  let i,j,k = versors of global reference system;
c      sx,sy,sz = versors of tensor;
c  and N = (k) x (sx) the Node line; then:
c
c	phi   is the angle from j to N (right-hand rotation on k)
c	theta is the angle from N to sy (right-hand rotation on sx)
c	omega is the angle from ij plane to sx (left-hand rot. on N).
c
c  In case of sx=k then N is null and:
c
c	phi+theta  is the angle from j to sy (right-hand rotation on k)
c	omega      is +90 (deg)
c       (the single values of phi and theta doesn't matter, so we can
c	 assume theta=0 and phi=phi+theta).
c
c  In case of sx=-k then N is null and:
c	phi-theta  is the angle from j to sy (right-hand rotation on k)
c	omega      is -90 (deg)
c       (the single values of phi and theta doesn't matter, so we can
c	 assume theta=0 and phi=phi-theta).
c
c
C  NOTE WELL: The same rotation matrix is defined by TWO equivalent 
c	      triplet of Euler angles:  (phi,theta,omega)  
c	      but also:           (phi+180,theta+180,180-omega)

c-----------------------------------------------
c     Node line: N=(k)x(sx)= -sxy i + sxx j
      absN=sqrt(sxy**2+sxx**2) ! Node line module
      
      if(absN.lt.EPS) then ! N is null and sx is parallel to k
        t=0
        if(sxz.gt.0) then ! sx=k
	  o=PI/2.
	else              ! sx=-k
	  o=-PI/2.
	endif
	p=acos(syy) ! cos(p)=(sy)(j) [scalar product] = syy
c	Phi sign correction:
c	sign must switch if (j x sy)(k)<0	
	if(syx.gt.0) p=-p ! rotation on k was left-hand
	
      else
c       Phi:
        p=acos(sxx/absN)
c	sign must switch if (j x N)(k)<0	
	if(sxy.lt.0) p=-p ! rotation on k was left-hand

c       Theta:
      	t=acos((syy*sxx-syx*sxy)/absN)
c	sign must switch if (N x sy)(sz)<0	
c	computing (N x sy)(sz):
	scalpr=(sxx*sxx+sxy*sxy)*syz-sxz*(sxy*syy+sxx*syx)
	if(scalpr.lt.0) t=-t ! rotation on sx was left-hand

c       Omega:
	o=acos(absN)
c	sign must switch if (sx x sx')(N)<0 [where sx' is the
c	projection of sx on plane xy], that is: 
c	if (sxx**2+sxy**2)*sxz > 0  --->  sxz > 0
	if(sxz.lt.0) o=-o ! rotation on N was right-hand
      endif
      
      WRITE (*,21) p,t,o 
      WRITE (3,21) p,t,o 

      if(DEBUGP) then ! DEBUG CONTROLS & PRINTS
        WRITE (*,*)
        print *,' Rotation matrix for permutation',n,':'
        print '(3(2x,f8.4))',sxx,syx,szx
        print '(3(2x,f8.4))',sxy,syy,szy
        print '(3(2x,f8.4))',sxz,syz,szz
        print *,' ================================================='
        WRITE (*,*)
c       Control: compute again the rotation matrix from Euler Angles
        sxx=cos(P)*cos(O)
        sxy=sin(P)*cos(O)
        sxz=sin(O)
        syx=-cos(T)*sin(P)-sin(O)*cos(P)*sin(T)
        syy=cos(T)*cos(P)-sin(O)*sin(P)*sin(T)
        syz=sin(T)*cos(O)
        szx=sin(T)*sin(P)-sin(O)*cos(P)*cos(T)
        szy=-sin(T)*cos(P)-sin(O)*sin(P)*cos(T)
        szz=cos(T)*cos(O)
        print *,' Rotation matrix re-computed from Euler Angles:'
        print '(3(2x,f8.4))',sxx,syx,szx
        print '(3(2x,f8.4))',sxy,syy,szy
        print '(3(2x,f8.4))',sxz,syz,szz
        print *,' ================================================='
        WRITE (*,*)
      endif
     
 21   FORMAT(1X,'PHI= ',F7.3,4X,'THETA= ',F7.3,4X,'OMEGA= ',F7.3)
c
      END
