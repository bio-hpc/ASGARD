      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      
      logical DEBUGP
      parameter (DEBUGP=.false.) ! DEBUG PRINT FLAG (print/don't print)

c     Precision threshold:
      parameter (EPS=1.e-6)
      
c*****************************************************************
c      RUNTIME ERROR CODES:
c 
c      code:atom     You referred to an invalid atom number.
c                    The atom numbers must lay between 1 and the
c                    total nr. of atoms in the structure.
c 
c      code:MAXDEF   You required too many atoms in your structure.
c                    If you really need this, please increase the
c                    MAXDEF parameter below and re-compile.
c 
c      code:MAXFE    You required too many paramagnetic ions in 
c                    your structure.
c                    If you really need this, please increase the
c                    MAXFE parameter below and re-compile.
c 
c      code:MAXOS    The number of observations is too high.  
c                    If you really need this, please increase the
c                    MAXOS parameter below and re-compile.
c 
c      code:MAXSTR   You required too many structures.
c                    If you really need this, please increase the
c                    MAXSTR parameter below and re-compile.
c 
c      code:neg      You inserted an invalid integer or character
c                    in input.
c
c      code:range    You inserted a value out of required range.
c
c*****************************************************************

c--------------------------------------------------------------------
      
      parameter (MAXOS=1000, MAXFE=5, MAXSTR=40, MAXDEF=4000,
     * MPAR=5, MP=MAXFE*MPAR+1, NP=MAXFE*MPAR, MH=MAXOS*MAXSTR)
c	MAXOS = Max nr. of observations
c	MAXFE = Max nr. of paramagnetic ions in every structure
c	MAXSTR = Max nr. of structures
c	MAXDEF = Max nr. of atoms in every structure
c	MPAR = Nr. of unknown parameters
c	MP = (see above)
c	NP = (see above)
c	MH = (see above)
c--------------------------------------------------------------------

      character fileout*20 ! Name of output file
      character nameatom*4,namat*4,namres*3,nam_at*4,nam_res*3,
     *          line*80,filename1*60,
     *		filename2*60,filename3*60
c	nameatom = Atom name in the aminoacid (read from input file
c		   'filename2') -- Temporary
c	namat = Name of atom whose pseudoshift has been observed
c		   (read from input file 'filename1')
c	namres = Aminoacid name (read from input file 'filename1')
c	nam_at = Atom name in the aminoacid (read from input file
c		   'filename2') -- Permanent
c	nam_res = Aminoacid name (read from input file 'filename2')
c	line = Buffer employed in file reading 
c	filename1 = Input file containing pseudocontact data
c	filename2 = Input file containg proteine structure data
c	filename3 = Output file name (so called "observed output file")

      dimension xp(MAXSTR,MAXDEF),yp(MAXSTR,MAXDEF),zp(MAXSTR,MAXDEF)
c	xp,yp,zp = Paramagnetic atoms coordinates (read from input 
c		   file 'filename2')

      dimension num_at(MAXSTR,MAXDEF),nam_at(MAXSTR,MAXDEF),
     *          nam_res(MAXSTR,MAXDEF),num_res(MAXSTR,MAXDEF),
     *          namat(MAXDEF),namres(MAXDEF)
c	num_at = Progressive nr. of paramagnetic atom 
c		 (read from input file 'filename2')
c	num_res =  Progressive nr. of aminoacid 
c		   (read from input file 'filename2')

      dimension numres(MH)
c	numres = Progressive nr. of aminoacid 
c		 (read from input file 'filename2')
c--------------------------------------------------------------------

      common /simplex/ 
     *        cx(MAXSTR,MAXOS),cy(MAXSTR,MAXOS),cz(MAXSTR,MAXOS),
     *        fx(MAXSTR,MAXFE),fy(MAXSTR,MAXFE),fz(MAXSTR,MAXFE),
     *        rccptx(MAXFE,3), rccpty(MAXFE,3),rccptz(MAXFE,3),
     *        center(MAXFE,3),ssx(MAXFE,3),ssy(MAXFE,3),ssz(MAXFE,3),
     *        tsx(MAXFE,3),tsy(MAXFE,3),tsz(MAXFE,3),
     *        TOLPROT(MH),WPROT(MH),shift(MH),
     *        obs(MH),mlprot(MH),
     *        phi(MAXFE),teta(MAXFE),omega(MAXFE),
     *        a1dip(MAXFE),a2dip(MAXFE),
     *        optphi(MAXFE),optteta(MAXFE),optomega(MAXFE),
     *        opta1(MAXFE),opta2(MAXFE),
     *        asxx(MAXFE,3), asyy(MAXFE,3), aszz(MAXFE,3),
     *        axsys(MAXFE,3), aysys(MAXFE,3), azsys(MAXFE,3),
     *        toldip,RESID,OLDRESID,ioldvio,IVIOLATION,intsys,
     *        nhp,nfe,nstr,numS,numD,numT,nat,nsystem,ngrid,
     *        line,filename1,filename2,filename3
c    cx,cy,cz = Spatial coordinates of paramagnatic atoms
c    fx,fy,fz = Spatial coordinates of paramagnatic ions
c    rccptx,rccpty,rccptz = coordinates (medium between structures) of
c			    points defining the molecular ref. syst.
c    center = coordinates of center (molecular ref. system)
c    ssx,ssy,ssz = coordinates of three versors x,y,z (molecular ref. system)
c    tsx,tsy,tsz = coordinates of three versors x,y,z (tensor)
c    axsys,aysys,azsys = director cosines of three versors (molecular ref. system)
c    mlprot = proton molteplicity  (read form input file 'filename1')
c    obs = observed shifts (read form input file 'filename1')
c    phi,teta,omega = Euler angles of tensor
c    a1dip,a2dip = Axial & rhombic anisotropy
c    optphi,optteta,optomega = OPTIMAL SOLUTION (Euler angles)
c    opta1,opta2 = OPTIMAL SOLUTION (axial & rhombic anisotropy)
c    asxx,asyy,aszz = Rotation matrix related to optphi,optteta,optomega
