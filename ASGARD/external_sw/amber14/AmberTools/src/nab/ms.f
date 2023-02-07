c molecular surface program
c ms
c
c December 16, 1983
c
c Copyright c 1983
c by Michael Connolly
c
c Written by Michael Connolly
c
c References:
c
c M.L. Connolly, "Solvent-accessible surfaces of proteins
c and nucleic acids", Science, 221, 709-713 (1983).
c
c M.L. Connolly, "Analytical molecular surface calculation",
c Journal of Applied Crystallography, 16, 548-558 (1983).
c
c This program may be freely distributed to anyone.
c
c It is written in fortran 77.
c
c ms calculates the molecular surface of a molecule
c given the coordinates of its atoms.  van der waals radii for
c the atoms and the probe radius must also be specified.
c
c The term molecular surface was introduced by F.M. Richards
c (Annual Reviews of Biophysics and Bioengineering, 1977,
c pages 151-176) with a specific meaning.  The surface
c Richards defined consists of two parts:  (1) the contact
c surface and (2) the reentrant surface.  He defines the
c contact surface to be that part of the van der waals
c surface of each atom which is accessible to a probe sphere
c of a given radius.  He defines the reentrant surface to be
c the inward-facing part of the probe sphere when it is
c simultaneously in contact with more than one atom.
c
c In implementing this definition I have found that there are
c two kinds of reentrant surface:  (1) concave reentrant
c surface, which is generated when the probe sphere
c simultaneously touches three atoms and (2) saddle-shaped
c reentrant surface, which is generated as the probe sphere
c rolls along the crevice between two atoms.
c I have also found that reentrant surface belonging to one
c probe may be contained in the interior volume of an
c overlapping probe and so must be removed.
c
c The input to this program consists of three files.
c
c The first file contains one record specifying the requested
c density of surface points (the number of surface points
c per unit area), the probe radius, the buried surface
c flag and the ascii/binary long/short output flag.
c Normally the buried surface flag will be blank or zero.
c If this flag is equal to 1, then the
c only surface calculated will be the surface of each
c molecule that is buried by the other molecules.
c If this flag is equal to 2, both buried and unburied
c surface are calculated, but the buried surface points
c are flagged by a 1 in the buried output field to the
c right of the surface normal.
c The ascii/binary long/short output flag
c may have one of these four values:
c
c 0      ascii       long
c 1      binary      long
c 2      ascii       short
c 3      binary      short
c
c The various output formats are discussed below.
c
c The format of the parameter record is: (2f10.5,2i5)
c
c The second file contains atomic radius records.
c each record has an integer that is the atom type,
c a van der waals radius for that type,
c and the point density for that type [added by dac 12/95!].
c The atom types need not be contiguous or in
c increasing order. the format is: (i5,2f10.5).
c
c The third file contains the atomic coordinate records.
c Each atomic coordinate record has the x, y and z coordinates
c of the atoms, followed by the atom type, a surface
c request number and a molecule number.
c The format is: (3f10.5,3i5).
c The surface request number may be 0, 1 or 2.
c 0 means that the atom is to be ignored,
c except for occupying a place in the sequence of atom numbers.
c 1 means that no surface is to be generated for this atom,
c but probe positions that collide with this atom will still
c be disallowed.  2 means that we are requested to generate
c surface for this atom.  In most cases 2 should be specified
c for all atoms. The molecule number is an integer that is
c used to divide the atoms into groups so that each group
c may be given its own surface. These groups need not
c correspond to actual molecules, as the program knows
c nothing about bonding. The characters to the right of
c these six fields are not read and may contain anything.
c
c The output from this program consists of two files,
c called 'contact' and 'reentrant', containing the contact
c and reentrant surface points.  All lines in both files have
c the same format.  The first three fields are the atom
c numbers of the atoms the probe was touching when it
c generated the given surface point.  The first number is
c the atom whose van der waals surface the point is closest to.
c The fourth field is the number of atoms the probe was
c touching. If this number is less then three, one or both
c of the second and third fields will be zero. The fifth, sixth
c and seventh fields are the coordinates of the surface point.
c The eighth field is the molecular area associated with
c the surface point.  The ninth, tenth and eleventh fields
c are a unit vector perpendicular to the surface pointing in
c the direction of the probe center.
c If the buried surface flag equals 2, there is a '1' written
c at the end of the record for buried surface points, and
c a '0' for accessible points. The output format is:
c (3i5,i2,3f9.3,4f7.3,i2)
c There is also a short output format: (3i5,i2,3f9.3,i2).
c The corresponding binary records are:
c (4i*2,7r*4,i*2) and (4i*2,3r*4,i*2)
c
c The contact file is generated in atom number order.
c The reentrant file should be sorted
c into atom number order and then merged with the contact file
c using the sort/merge utility programs available at your
c installation.  Then all the surface points belonging to a
c given atom will be in a contiguous series of records.
c
c
c The program has the ability to calculate the
c van der waals surface
c of a molecule simply by specifying a probe radius of zero.
c The reentrant code is bypassed completely and there are no
c before and reentrant files. For a probe radius of zero
c the van der waals surface and the contact surface
c are equivalent.
c
c The flow of the program may be described in general terms
c as follows. First all the input is read. Then the contact
c and reentrant surface is generated. The contact surface
c is in its final form, but the reentrant surface is
c written to a temporary file, called 'before'.
c Each reentrant probe position is written to this file,
c followed by all its surface points.
c After all the contact and reentrant surface has been generated,
c the 'before' file is read and a final reentrant surface file is
c written which contains all the reentrant surface points
c not lying within any reentrant probe.
c
c
c Sometimes it is desirable to calculate the surface of only
c part of a molecule. One cannot simply remove the remaining
c atoms from the input file as this will generate false
c surfaces. One could calculate a surface for the entire
c molecule and then edit out that part of the surface belonging
c to the atoms of interest, but this is a needless
c waste of computer time.  The proper way to accomplish this
c task is to place the probe next to the atoms whose surface
c is requested and to use the remaining atoms only for
c collision checks so that false surfaces will not be
c generated. When the probe is placed next to two or three
c atoms, only one of these need be an atom whose surface
c is requested. only that part of the arc or spherical
c triangle that belongs to the atoms whose surface is
c requested will be written as output.
c
c This program allows the calculation of individual surfaces
c for several interacting molecules. this may be done with
c or without the buried surface option.
c
c This program should be compiled with integer*4 as the default.
c It uses no include files or libraries.
c
      program ms
c
c maxatm     maximum number of atoms
c maxtyp     maximum number of atom types
c maxnbr     maximum number of neighbors an atom may have
c maxsph     maximum number of surface points on a sphere
c maxcir     maximum number of surface points on a circle
c maxarc     maximum number of surface points on an arc
c maxppp     maximum number of surface points per probe
c maxyon     maximum number of yon probes
c maxvic     maximum number of victim probes
c maxeat     maximum number of eaters of a probe's surface
c maxcub     maximum number of cubes in one direction
c
c maxatm must be greater than or equal to maxyon
c because they share the same cubing arrays
c
c
      integer maxatm, maxtyp, maxnbr, maxsph, maxcir, maxarc,
     $    maxppp, maxyon, maxvic, maxeat, maxcub

      parameter (maxatm=15000)
      parameter (maxtyp=6000)
      parameter (maxnbr=2000)
      parameter (maxsph=1000)
      parameter (maxcir=1000)
      parameter (maxarc=1000)
      parameter (maxppp=1000)
      parameter (maxyon=12000)
      parameter (maxvic=6000)
      parameter (maxeat=1000)
      parameter (maxcub=40)



c
c run-time options and parameters
c rp	probe radius, set by -rp VALUE
c ar	additional radius for solvent acc surface with rp=0
C  , set by -ar VALUE
c d	dot density, set by -d VALUE
c ibury	switch for turning on buried surface, set by -b or -bonly
c binary     binary output records instead of ascii, set by -bin
c short      short output records instead of ordinary long, set by -short
c
      real rp,ar
      integer ibury
      integer iabls
      logical binary, short
c
c atom type arrays
c
c itype     atom type number
c rtype     van der waals radius
c nua       number of unit vectors on sphere
c
c dimensioned one more in case input file is too long
      integer itype(maxtyp+1)
      real rtype(maxtyp+1),d(maxtyp+1)
      integer nua(maxtyp)
      integer ntype
      real radmax
c
c arrays for all atoms
c
c co        atomic coordinates
c ias       surface request number
c iat       atom itype
c molnum    molecule number
c rad       radius
c srs       some reentrant surface
c
c dimension arrays 1 more in case input file is too long
      real co(3,maxatm+1)
      integer ias(maxatm+1),molnum(maxatm+1)
      integer iat(maxatm+1)
      real rad(maxatm)
      logical srs(maxatm)
c
c cube arrays
c
c ico       integer cube coordinates
c icuptr    pointer to next atom in cube
c comin     minimum atomic coordinates (cube corner)
c icube     pointer to first atom in list for cube
c scube     some atom with srn=2 in cube
c sscube    some atom with srn=2 in cube or adjoining cubes
c
      integer ico(3,maxatm),icuptr(maxatm)
      real comin(3)
      integer icube(maxcub,maxcub,maxcub)
      logical scube(maxcub,maxcub,maxcub)
      logical sscube(maxcub,maxcub,maxcub)
c
c neighbor arrays
c
c inbr      atom number
c cnbr      coordinates
c rnbr      radius
c snbr      true if srn = 2 for neighbor
c mnbr      mutual neighbor of iatom and jatom
c molnbr    molecule number
c ernbr     expanded radius (rnbr + rp)
c disnbr    distance from neighbor to iatom
c lknbr     link to next farthest out neighbor
c itnl      temporary neighbor list (before sort)
c
      integer inbr(maxnbr)
      real cnbr(3,maxnbr),rnbr(maxnbr)
      logical snbr(maxnbr),mnbr(maxnbr)
      integer molnbr(maxnbr)
      real ernbr(maxnbr)
      real disnbr(maxnbr)
      integer lknbr(maxnbr)
      integer itnl(maxnbr)
c
c circle and sphere unit vector arrays
c
c up        unit vectors for probe
c ua        unit vectors for atom
c eva       extended vectors for atom
c circle    points on a circle
c
      real up(3,maxsph,maxtyp)
      real ua(3,maxsph,maxtyp),eva(3,maxsph,maxtyp)
      real circle(3,maxcir,maxtyp)
      integer nup(maxtyp),ncirc(maxtyp),ptyp
c
c ci        coordinates of atom i
c cj        coordinates of atom j
c ck        coordinates of atom k
c si        srn = 2 for atom i
c sj        srn = 2 for atom j
c sk        srn = 2 for atom k
c
      real ci(3), cj(3), ck(3)
      logical si,sj,sk,sns
c
c geometric construction vectors
c
c vij       vector from atom i to atom j
c uij       unit vector of vij
c q,t       two perpendicular vectors in the saddle plane
c cijk      center of circle of intersection of expanded sphere
c           of atom k with the saddle plane of atoms i and j
c vijk      vector from torus center to cijk
c uijk      unit vector of vijk
c bij       torus center
c aij       starting altitude vector from torus center
c bijk      base point for atoms i, j and k
c aijk      altitude vector for atoms i, j and k
c aijp      altitude vector to probe (rotated)
c a         altitude vector, general, used in orsr
c p         probe coordinates (general, used in orsr)
c pijp      center of probe placed tangent to atoms i and j
c pij       starting center of probe tangent to atoms i and j
c pijk      probe placed tangent to atoms i, j and k
c pipt      probe placed tangent to atom i
c vpi       vector from probe center to contact point with atom i
c vpj       vector from probe center to contact point with atom j
c vpk       vector from probe center to contact point with atom k
c vps0      starting arc points relative to probe center
c vbs0      starting arc points relative to torus center
c vbs       rotated arc point relative to torus center
c arca      area of arc point
c ayon      arc point on yon side of symmetry element
c vector    temporary vector storage
c
c
      real vij(3),uij(3),q(3),t(3),cijk(3),vijk(3),uijk(3)
      real bij(3),aij(3),bijk(3),aijk(3,2),aijp(3,2),a(3)
      real p(3),pij(3),pijp(3,2),pijk(3,2),pipt(3)
      real vpi(3),vpj(3),vpk(3)
      real vps0(3,maxarc),vbs0(3,maxarc,2),vbs(3),arca(maxarc)
      real vector(3)
c
c g         uij, q, t frame for torus
c h         rotation about x-axis
c ghgt      rotation about uij axis
c pow       powers of ghgt
c
c rotation matrices
      real g(3,3),h(3,3),ghgt(3,3),pow(3,3)
c
c s         surface points for probe
c torus     surface points for torus
c n1        atom surface point is on or closest to
c n2,n3     other atoms probe is touching
c yon       whether point lies on yon side of symmetry element
c
c reentrant probe record
      real s(3,maxppp),torus(maxppp)
      integer n1(maxppp),n2(maxppp),n3(maxppp)
      logical yon(maxppp)
      integer imol,np
c
c both      both probe positions free of collisions
c pair      this member of pair free from collisions
c yonprb    probe crosses symmetry element
c found     search flag
c
c logical variables
      logical both,pair(2),ayon(maxarc)
      logical found
      logical yonprb
c
c orsr for non-symmetry-related probes
c the factor of three for the victim arrays
c is based upon experience
c
c py        center of yon probe
c ay        altitude vector of yon probe
c pv        center of victim probe
c av        altitude vector of victim probe
c ivic      list of probe numbers of victims
c ivicp     pointer to next victim with same hash
c molyon    molecule number of yon probe
c molvic    molecule number of victim probe
c eat       coordinates of eaters of probe
c
      real py(3,maxyon),ay(3,maxyon)
      real pv(3,maxvic),av(3,maxvic)
      integer ivic(maxvic),ivicp(maxvic)
      integer molyon(maxyon),molvic(maxvic)
      real eat(3,maxeat)
c
c multiple-molecule variables
c
c bury      probe position is buried
c
      logical bury
c
c counters
c
c nias      number of atoms with given srn
c nshape    number of surface points with given shape
c nlost     number of surface points lost in non-symmetry orsr
c
      integer nias(3)
      integer*4 nshape(3),nlost(3)
c
c output variables
c
c iatnum   atom numbers
c ishape   surface shape number
c ib       buried surface flag 
c 
        integer iatnum(3),ishape,ib
      real outco(3),outvec(3)
c
c
cMAS	common block for optimization
      common/mpck/colpre(3),radpre,r2pre
c	real	colpre(3)
c	real 	radpre,r2pre
c
c
c
c logical functions
c
      logical collid
      logical buried
c
c . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
c
      real pi
      parameter (pi=3.141592654)
c
        character*60 arg
c
c Set default values for run options:
c
 	dens = 8.
      ar = 0.0
      rp = 1.4
      ibury = 0
      iabls = 0
      short = .false.
      binary = .false.
c
c    --- process command-line arguments:
c
      iarg = 0
      indx = iargc()
      if (indx .eq. 0) goto 2
    1 continue
           iarg = iarg + 1
           call getarg(iarg,arg)
           if (arg .eq. '-d') then
                iarg = iarg + 1
                call getarg(iarg,arg)
                read(arg,*) dens
           elseif (arg .eq. '-rp') then
                iarg = iarg + 1
                call getarg(iarg,arg)
                read(arg,*) rp

           elseif (arg .eq. '-ar') then
                iarg = iarg + 1
                call getarg(iarg,arg)
                read(arg,*) ar
           
           

           elseif (arg .eq. '-dfar') then
                iarg = iarg + 1
                call getarg(iarg,arg)
                read(arg,*) dfar
           elseif (arg .eq. '-short') then
                short = .true.
           elseif (arg .eq. '-bin') then
                binary = .true.
           elseif (arg .eq. '-ibury') then
                iarg = iarg + 1
                call getarg(iarg,arg)
                read(arg,*) ibury
           elseif (arg .eq. '-pqr') then
                iarg = iarg + 1
                call getarg(iarg,arg)
                  open(unit=2,file=arg,status='old',form='formatted')
                  rewind(2)
           elseif (arg .eq. '-help') then
                write(6,3)
                stop
           else
                if (arg .eq. ' ') go to 2
                write(6,'(/,5x,a,a)') 'unknown flag: ',arg
                write(6,3)
                stop
           endif
      if (iarg .lt. indx) go to 1
c
    2 continue
    3 format('Usage: ms -d <#> -rp <#> -pqr <file> ')
c
c input value checking
c check for negative probe radius
      if (rp .lt. 0.0) call error(120,0,rp)
c check buried surface flag
      if (ibury .lt. 0 .or. ibury .gt. 2) call error(127,ibury,0.0)
c


c
c initialize srn counters and coordinate minima
      do 350 k = 1,3
         nias(k) = 0
350     continue
      do 400 k = 1,3
         comin(k) = 1000000.0
400     continue



c initialization
      ntype = 1
      radmax = 0.0
      n = 1 
c read atom types with their radii
100	continue
      read (2,150,end=300) itype(ntype),co(1,n),co(2,n),co(3,n),
     $d(ntype),rtype(ntype)

      iat(n) = itype(n)
      ias(n) = 2
      molnum(n) =1
      rtype(ntype) = rtype(ntype) + ar
       


c     
c check for atom overflow
      if (n .gt. maxatm) call error(150,n,0.0)
c check surface request number
      if (ias(n) .lt. 0 .or. ias(n) .gt. 2) call error(160,n,0.0)
c check for new coordinate minima
      do 550 k = 1,3 
         if (co(k,n) .lt. comin(k)) comin(k) = co(k,n)
550     continue
c increment counters for each srn type
      do 600 k = 1,3
         if (ias(n) .eq. k-1) nias(k) = nias(k) + 1
600     continue
c we don't care whether ignored atoms have a radius
      if (ias(n) .eq. 0) go to 700
c     
c look for this atom type number so we may assign a radius
      found = .false.
      do 650 i = 1, ntype
         if (iat(n) .ne. itype(i)) go to 650
         found = .true.
c transfer radius from atom type to atom radius array
         rad(n) = rtype(i)
c end of atom type search loop
650     continue
c check whether atom type found
      if (.not. found) call error(170,n,0.0)
700     continue



150	format(4x,i7,20x,3f8.3,2f8.4)
      if (d(ntype).le.0.0) d(ntype) = dens
c atom radii must be zero or positive
      if (rtype(ntype) .lt. 0.000001) call error(140,ntype,rtype(ntype))
c check for atom type overflow
      if (ntype .gt. maxtyp) call error(130,ntype,0.0)
c check for new maximum radius
      if (rtype(ntype) .gt. radmax) radmax = rtype(ntype)
c number of unit vectors depends on sphere area and input density
      nua(ntype) = (4*pi*rtype(ntype)**2) * d(ntype)
c decrease to array size if too large
      if (nua(ntype) .gt. maxsph) nua(ntype) = maxsph
      if (nua(ntype) .lt. 1) nua(ntype) = 1
c create unit vector arrays
      call genun(ua(1,1,ntype),nua(ntype))
c compute extended vectors for later probe placement
      do 250 isph = 1,nua(ntype)
         do 200 k = 1,3
            eva(k,isph,ntype) = (rtype(ntype) + rp) * ua(k,isph,ntype)
200	   continue
250	continue
c one more atom type
      ntype = ntype + 1
      n = n +1
      natom = n

      go to 100
300	continue
      close(2)
c decrement on end of file
      ntype = ntype - 1
      natom = natom -1

c calculate width of cube from maximum atom radius and probe radius
      width = 2 * (radmax + rp)






c==============================================================================
c
c     set up cube arrays
c     first the integer coordinate arrays
      do 850 i = 1,natom
         do 800 k = 1,3
            ico(k,i) = (co(k,i)-comin(k))/width + 1
            if (ico(k,i) .lt. 1) stop 'cube coordinate too small'
            if (ico(k,i) .gt. maxcub)
     $    stop 'cube coordinate too large'
800	   continue
850	continue
c
c initialize head pointer and srn=2 arrays
      do 1000 k = 1,maxcub
         do 950 j = 1,maxcub
            do 900 i = 1,maxcub
               icube(i,j,k) = 0
               scube(i,j,k) = .false.
               sscube(i,j,k) = .false.
900	      continue
950	   continue
1000	continue
c
c initialize linked list pointers
      do 1050 i = 1,natom
         icuptr(i) = 0
1050	continue
c
c set up head and later pointers for each atom
      do 1250 iatom = 1,natom
c skip atoms with surface request numbers of zero
         if (ias(iatom) .eq. 0) go to 1250
         i = ico(1,iatom)
         j = ico(2,iatom)
         k = ico(3,iatom)
         if (icube(i,j,k) .le. 0) then
c     first atom in this cube
            icube(i,j,k) = iatom
         else
c     add to end of linked list
            iptr = icube(i,j,k)
1100	      continue
c check for duplicate coordinates
            if( molnum(iatom) .eq. molnum(iptr) .and.
     $          dist2(co(1,iatom),co(1,iptr)) .le. 0.0 ) then
               ias(iatom) = 0
               write (6,1150) iatom,iptr
1150	         format(1x,'atom',i5,' dropped (same co as ',i5,')')
               go to 1250
            end if
            if (icuptr(iptr) .le. 0) go to 1200
c move on down the list
            iptr = icuptr(iptr)
            go to 1100
1200	      continue
c store atom number
            icuptr(iptr) = iatom
         end if
c check for surfaced atom
         if (ias(iatom) .eq. 2) scube(i,j,k) = .true.
1250	continue
c
c check for 3 x 3 x 3 with some srn = 2
c
      do 1550 k = 1,maxcub
         do 1500 j = 1,maxcub
            do 1450 i = 1,maxcub
               if (icube(i,j,k) .eq. 0) go to 1450
c check whether this cube or any adjacent cube has srn = 2
               do 1400 k1 = k-1,k+1
                  if (k1 .lt. 1 .or. k1 .gt. maxcub) go to 1400
                  do 1350 j1 = j-1,j+1
                     if (j1 .lt. 1 .or. j1 .gt. maxcub) go to 1350
                     do 1300 i1 = i-1,i+1
                        if (i1 .lt. 1 .or. i1 .gt. maxcub) go to 1300
                        if (scube(i1,j1,k1)) sscube(i,j,k) = .true.
1300	               continue
1350	            continue
1400	         continue
1450	      continue
1500	   continue
1550	continue
c
c initialization
c maximum number of neighbors any atom has
      maxnb = 0
c numbers of surface points
      do 1600 k = 1,3
         nshape(k) = 0
         nlost(k) = 0
1600	continue
c number of yon probes
      ny = 0
c contact and reentrant areas
      areac = 0.0
      arear = 0.0
c
c write out messages
      write (6,1650) natom
1650	format(1x,i5,1x,'atoms')
      write (6,1700) nias(1)
1700	format(1x,i5,1x,'omitted')
      write (6,1750) nias(2)
1750	format(1x,i5,1x,'collision only')
      write (6,1800) nias(3)
1800	format(1x,i5,1x,'surface')
      write (6,1850) dens,rp
1850	format(1x,'surface point density = ',f10.5,5x,
     $   'probe radius = ',f10.5)
      if (ibury .eq. 1) write (6,1900)
1900	format(1x,'buried surface only')
      if (ibury .eq. 2) write (6,1950)
1950	format(1x,'buried surface flagged')
      if (binary) write (6,2000)
2000	format(1x,'binary contact and reentrant files')
      if (short) write (6,2050)
2050	format(1x,'short output records')
c stop if density is not positive
      if (dens .le. 0.0) stop 'non-positive density'
c skip probe and circle setup if van der waals surface
      if (rp .eq. 0.0) go to 2150
      do 2110 ptyp=1,ntype
c
c set up probe sphere and circle
c
        nup(ptyp) = (4 * pi * rp ** 2) * d(ptyp)
        if (nup(ptyp) .lt. 1) nup(ptyp) = 1
        if (nup(ptyp) .gt. maxsph) nup(ptyp) = maxsph
        call genun(up(1,1,ptyp),nup(ptyp))
        ncirc(ptyp) = (2 * pi * rp) * sqrt(d(ptyp))
        if (ncirc(ptyp) .lt. 1) ncirc(ptyp) = 1
        if (ncirc(ptyp) .gt. maxcir) ncirc(ptyp) = maxcir
        do 2100 i = 1, ncirc(ptyp)
           fi = (2 * pi * (i-1))/ncirc(ptyp)
           circle(1,i,ptyp) = rp * cos(fi)
           circle(2,i,ptyp) = rp * sin(fi)
           circle(3,i,ptyp) = 0.0
2100	continue
2110  continue
c
c open before file for writing
      open(4,file='before',form='unformatted',status='unknown')
      rewind(4)
c skip to here if no reentrant surface will be calculated
2150	continue
c
c open contact file for writing

      if (binary) then
         open(7,file='contact',form='unformatted',status='unknown')
      else
         open(7,file='contact',status='unknown')
      end if

      rewind(7)
                 write (7,'(i5)') natom
c
c initialize some reentrant surface to false for each atom
      do 2200 iatom = 1,natom
         srs(iatom) = .false.
2200	continue
c
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c big loop for each atom
      do 5650 iatom = 1, natom
c skip ignored atoms
         if (ias(iatom) .eq. 0) go to 5650
c find the index into the atom type arrays for iatom
         do 2230 idx = 1,ntype
            if (itype(idx) .eq. iat(iatom)) go to 2240
2230	   continue
         stop 'logic error in ms regarding atom types'
2240   ptyp = idx
c
         ici = ico(1,iatom)
         icj = ico(2,iatom)
         ick = ico(3,iatom)
c skip iatom if its cube and adjoining cubes contain only blockers
         if (.not. sscube(ici,icj,ick)) go to 5650
c transfer values from large arrays to iatom variables
         ri = rad(iatom)
         si = ias(iatom) .eq. 2
         do 2250 k = 1,3
            ci(k) = co(k,iatom)
2250	   continue
         imol = molnum(iatom)
c
c gather the neighboring atoms of iatom
c initialize number of neighbors, and number of neighbors in the
c same molecule as atom i
         nnbr = 0
         nimol = 0
c initialize srn = 2 for some neighbor to false
         sns = .false.
c save a little time for distance check
         sumi = 2 * rp + ri
c check iatom cube and adjacent cubes for neighboring atoms
         do 2550 jck = ick-1,ick+1
            if (jck .lt. 1 .or. jck .gt. maxcub) go to 2550
            do 2500 jcj = icj-1,icj+1
               if (jcj .lt. 1 .or. jcj .gt. maxcub) go to 2500
               do 2450 jci = ici-1,ici+1
                  if (jci .lt. 1 .or. jci .gt. maxcub) go to 2450
                  jatom = icube(jci,jcj,jck)
2300	            continue
c check for end of linked list for this cube
                  if (jatom .le. 0) go to 2400
c distance check
                  sum = sumi + rad(jatom)
                  vect1 = abs(co(1,jatom) - ci(1))
                  if (vect1 .ge. sum) go to 2350
                  vect2 = abs(co(2,jatom) - ci(2))
                  if (vect2 .ge. sum) go to 2350
                  vect3 = abs(co(3,jatom) - ci(3))
                  if (vect3 .ge. sum) go to 2350
                  d2 = vect1 ** 2 + vect2 ** 2 + vect3 ** 2
                  if (d2 .ge. sum ** 2) go to 2350
c iatom is not its own neighbor
                  if (iatom .eq. jatom) go to 2350
c we have a new neighbor
                  nnbr = nnbr + 1
c check for neighbor overflow
                  if (nnbr .gt. maxnbr) call error(210,nnbr,0.0)
c save atom number in temporary array
                  itnl(nnbr) = jatom
c check whether surfaced neighbor in same molecule
                  if (ias(jatom) .eq. 2 .and. molnum(jatom) .eq. imol)
     $    sns = .true.
c count the number of atoms in the same molecule as iatom
                  if (imol .eq. molnum(jatom)) nimol = nimol + 1
2350	            continue
c get number of next atom in cube
                  jatom = icuptr(jatom)
                  go to 2300
2400	            continue
2450	         continue
2500	      continue
2550	   continue
c keep track of maximum number of neighbors
c for array-dimensioning purposes
         if (nnbr .gt. maxnb) maxnb = nnbr
c
c no surface for atom i if buried only flag set and
c there are no neighbors from the other molecules
         if (ibury .eq. 1 .and. nimol .eq. nnbr) go to 5650
c
c no surface if iatom and all neighbor atoms
c in the same molecule have surface request numbers < 2
         if (.not. si .and. .not. sns) go to 5650
c
c set up neighbors arrays with jatom in increasing order
c
c initialize minimum neighbor atom number
         jmold = 0
         do 2700 iuse = 1,nnbr
            jmin = natom + 1
            do 2600 jnbr = 1,nnbr
c don't use ones already sorted
               if (itnl(jnbr) .le. jmold) go to 2600
               if (itnl(jnbr) .lt. jmin) then
                  jmin = itnl(jnbr)
                  jminbr = jnbr
               end if
2600	      continue
            jmold = jmin
            jnbr = jminbr
            jatom = itnl(jnbr)
c transfer atom number, coordinates, radius, surface request number,
c molecule number, expanded radius, distance from iatom
            inbr(iuse) = jatom
            do 2650 k = 1,3
               cnbr(k,iuse) = co(k,jatom)
2650	      continue
            rnbr(iuse) = rad(jatom)
            snbr(iuse) = ias(jatom) .eq. 2
            molnbr(iuse) = molnum(jatom)
            ernbr(iuse) = rnbr(iuse) + rp
            disnbr(iuse) = dist2(ci,cnbr(1,iuse))
c initialize link to next farthest out neighbor
            lknbr(iuse) = 0
2700	   continue
c set up a linked list of neighbors in order of
c increasing distance from iatom
c initialize pointer to first neighbor to 0
         lkf = 0
c look for neighbor in same molecule
c we want only atoms in same molecule for collision check
         do 2750 l = 1,nnbr
            if (imol .ne. molnbr(l)) go to 2750
            lkf = l
            go to 2800
2750	   continue
         if (lkf .eq. 0) go to 3000
2800	   continue
c put remaining neighbors in linked list at proper position
         do 2950 l = lkf+1,nnbr
            if (imol .ne. molnbr(l)) go to 2950
            l1 = 0
            l2 = lkf
2850	      continue
            if (disnbr(l) .lt. disnbr(l2)) go to 2900
            l1 = l2
            l2 = lknbr(l2)
            if (l2 .ne. 0) go to 2850
2900	      continue
c add to list
            if (l1 .eq. 0) then
               lkf = l
               lknbr(l) = l2
            else
               lknbr(l1) = l
               lknbr(l) = l2
            end if
2950	   continue
3000	   continue
c
c no reentrant surface will be calculated if we are
c calculating the van der waals surface
c instead of the molecular surface
         if (rp .eq. 0.0) go to 5200
c no reentrant surface if iatom has no neighbors
         if (nimol .le. 0) go to 5200
c
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
c medium loop for each neighbor of iatom
c
         do 5150 jnbr = 1,nnbr
            jatom = inbr(jnbr)
c
c each pair of atoms is considered only once
            if (jatom .le. iatom) go to 5150
c each molecule gets a separate surface
            if (imol .ne. molnbr(jnbr)) go to 5150
c
c tranfer from neighbor arrays to jatom variables
            rj = rnbr(jnbr)
            sj = snbr(jnbr)
            do 3050 k = 1,3
               cj(k) = cnbr(k,jnbr)
3050	      continue
c
c here follow geometric calculations of points, vectors and
c distances used for probe placement in both saddle and
c concave reentrant surface generation
c
c calculate the intersection
c of the expanded spheres of iatom and jatom
c this circle is called the saddle circle
c the plane it lies in is called the saddle plane
c
            do 3100 k = 1,3
               vij(k) = cj(k) - ci(k)
3100	      continue
c create an orthonormal frame
c with uij pointing along the inter-atomic axis
c and q and t defining the saddle plane
            if (anorm(vij) .le. 0.0) then
               write (6,3150) iatom,jatom
3150	         format(1x,'atoms',2i5,' have the same center')
               go to 5150
            end if
            call vnorm(vij,uij)
            call vperp(uij,q)
            call cross(uij,q,t)
c
c calculate the saddle circle center and radius
            dij = anorm(vij)
            f = 0.5*(1.0+((ri+rp)**2-(rj+rp)**2)/dij**2)
c base point
            do 3200 k = 1,3
               bij(k) = ci(k) + f * vij(k)
3200	      continue
            f1 = (ri+rj+2*rp)**2 - dij**2
c skip to bottom of middle loop if atoms are too far apart
            if (f1 .le. 0.0) go to 5150
            f2 = dij**2 - (ri-rj)**2
c skip to bottom of middle loop if one atom inside the other
            if (f2 .le. 0.0) go to 5150
c height (radius of saddle circle)
            hij = sqrt(f1*f2) / (2*dij)
c a starting altitude
            do 3250 k = 1,3
               aij(k) = hij * q(k)
3250	      continue
c
c
c concave reentrant surface
c
c gather mutual neighbors of iatom and jatom
            mutual = 0
            do 3300 knbr = 1, nnbr
               d2 = dist2(cj,cnbr(1,knbr))
               mnbr(knbr) = d2 .lt. (2*rp+rj+rnbr(knbr))**2
     $   .and. knbr .ne. jnbr
               if (mnbr(knbr)) mutual = mutual + 1
3300	      continue
c
c . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
c
c inner loop for each mutual neighbor of iatom and jatom
            ishape = 3
            do 4200 knbr = 1,nnbr
               if (.not. mnbr(knbr)) go to 4200
               katom = inbr(knbr)
c iatom < jatom < katom
               if (katom .le. jatom) go to 4200
               sk = snbr(knbr)
c skip neighbor if all three atom not marked to be surfaced
               if (.not. (si .or. sj .or. sk)) go to 4200
c each molecule gets a separate surface
               if (imol .ne. molnbr(knbr)) go to 4200
c
c tranfer from neighbor array to katom variables
               rk = rnbr(knbr)
               do 3350 k = 1,3
                  ck(k) = cnbr(k,knbr)
3350	         continue
c
c calculate intersection of expanded sphere of katom
c with saddle plane. we will call this the katom circle.
c
c projection of vector,
c from katom to a point on the saddle plane,
c onto iatom-jatom axis,
c in order to get distance katom is from saddle plane
               dk = uij(1) * (bij(1)-ck(1)) + uij(2) * (bij(2)-ck(2)) +
     $   uij(3) * (bij(3)-ck(3))
c
c calculate radius of katom circle
               rijk = (rk+rp) ** 2 - dk ** 2
c skip concave calculation if no intersection
               if (rijk .le. 0.0) go to 4200
               rijk = sqrt(rijk)
c calculate center of katom circle
               do 3400 k = 1,3
                  cijk(k) = ck(k) + dk * uij(k)
3400	         continue
c
c calculate intersection of the katom circle with the saddle circle
               do 3450 k = 1,3
                  vijk(k) = cijk(k) - bij(k)
3450	         continue
               dijk = anorm(vijk)
               if (dijk .le. 0.0) then
                  write (6,3500) iatom,jatom,katom
3500	            format(1x,'atoms',3i5,' have concentric circles')
                  go to 4200
               end if
               f = 0.5 * (1.0+(hij**2-rijk**2)/dijk**2)
c base point bijk is on symmetry plane and saddle plane
               do 3550 k = 1,3
                  bijk(k) = bij(k) + f * vijk(k)
3550	         continue
               f1 = (hij+rijk)**2 - dijk**2
c skip to bottom of inner loop if katom too far away
               if (f1 .le. 0.0) go to 4200
               f2 = dijk**2 - (hij-rijk)**2
c skip to bottom of inner loop if katom circle inside saddle circle
c or vice-versa
               if (f2 .le. 0.0) go to 4200
               hijk = sqrt(f1*f2) / (2*dijk)
               call vnorm(vijk,uijk)
c uij and uijk lie in the symmetry plane passing through the atoms
c so their cross product is perpendicular to this plane
               call cross(uij,uijk,aijk)
c two altitudes
               do 3600 k = 1,3
                  aijk(k,1) = hijk * aijk(k,1)
                  aijk(k,2) = - aijk(k,1)
3600	         continue
c
c probe placement at ends of altitude vectors
               do 3700 ip = 1,2
                  do 3650 k = 1,3
                     pijk(k,ip) = bijk(k) + aijk(k,ip)
3650	            continue
c collision check with mutual neighbors
                  pair(ip) = .not. collid(pijk(1,ip),rp,cnbr,ernbr,mnbr,
     $    nnbr,maxnbr,ishape,jnbr,knbr,molnbr,imol,lkf,lknbr)
3700	         continue
c if neither probe position is allowed, skip to bottom of inner loop
               if (.not. pair(1) .and. .not. pair(2)) go to 4200
               both = pair(1) .and. pair(2)
c some reentrant surface for all three atoms
               srs(iatom) = .true.
               srs(jatom) = .true.
               srs(katom) = .true.
c
c generate surface points
               area = (4 * 3.14159 * rp ** 2)/nup(ptyp)
               do 4150 ip = 1,2
                  if (.not. pair(ip)) go to 4150
c give it some kind of value, in case we don't call buried
                  bury = .false.
c only call buried if we care what the answer is
                  if (ibury .gt. 0) bury = buried(pijk(1,ip),rp,cnbr,rnb
     $   r,mnbr,nnbr,maxnbr,ishape,jnbr,knbr,molnbr,imol)
c skip if not buried and buried surface only flag set
                  if (ibury .eq. 1 .and. .not. bury) go to 4150
c determine whether probe has surface on far side of plane
                  yonprb = hijk .lt. rp .and. .not. both
c calculate vectors defining spherical triangle
c the vectors are given the probe radius as a length
c only for the purpose of making the geometry more clear
                  do 3750 k = 1,3
                     vpi(k) = (ci(k) - pijk(k,ip)) * rp / (ri + rp)
                     vpj(k) = (cj(k) - pijk(k,ip)) * rp / (rj + rp)
                     vpk(k) = (ck(k) - pijk(k,ip)) * rp / (rk + rp)
3750	            continue
                  sign = det(vpi,vpj,vpk)
c initialize number of surface points written
                  np = 1
c gather points on probe sphere lying within triangle
                  do 4000 i = 1,nup(ptyp)
c if the unit vector is pointing away from the symmetry plane
c the surface point cannot lie within the inward-facing triangle
                     if (dot(up(1,i,ptyp),aijk(1,ip)) .gt. 0.0) 
     $                        go to 4000
                     if (sign * det(up(1,i,ptyp),vpj,vpk) .lt. 0) 
     $                        go to 4000
                     if (sign * det(vpi,up(1,i,ptyp),vpk) .lt. 0) 
     $                        go to 4000
                     if (sign * det(vpi,vpj,up(1,i,ptyp)) .lt. 0) 
     $                        go to 4000
                     if (np .gt. maxppp) call error(320,np,0.0)
c calculated whether point is on yon side of plane
               yon(np) = aijk(1,ip) * (aijk(1,ip) + up(1,i,ptyp)) +
     $   aijk(2,ip) * (aijk(2,ip) + up(2,i,ptyp)) +
     $   aijk(3,ip) * (aijk(3,ip) + up(3,i,ptyp)) .lt. 0.0
c overlapping reentrant surface removal
c for symmetry-related probe positions
                     if (yon(np) .and. both) go to 4000
c calculate coordinates of surface point
                     do 3800 k = 1,3
                        s(k,np) = pijk(k,ip) + up(k,i,ptyp) * rp
3800	               continue
c find the closest atom and put the three atom numbers
c in the proper order
c n1 is closest, n2 < n3
                     dsi = dist(s(1,np),ci) - ri
                     dsj = dist(s(1,np),cj) - rj
                     dsk = dist(s(1,np),ck) - rk
                     if (dsi .le. dsj .and. dsi .le. dsk) go to 3850
                     if (dsj .le. dsi .and. dsj .le. dsk) go to 3900
                     if (.not. sk) go to 4000
                     n1(np) = katom
                     n2(np) = iatom
                     n3(np) = jatom
                     go to 3950
3850	               continue
                     if (.not. si) go to 4000
                     n1(np) = iatom
                     n2(np) = jatom
                     n3(np) = katom
                     go to 3950
3900	               continue
                     if (.not. sj) go to 4000
                     n1(np) = jatom
                     n2(np) = iatom
                     n3(np) = katom
3950	               continue
                     np = np + 1
c end of nup loop
4000	            continue
                  np = np - 1
c skip the write if no points
                  if (np .le. 0) go to 4150
c
c
c write the molecule number, shape, number of points,
c probe position and
c the vector from the base to the probe center
                  write (4) imol,ishape,np,(pijk(k,ip),k=1,3),
     $   (aijk(k,ip),k=1,3),yonprb,bury
c save probe in yon probe arrays
                  if (yonprb) then
c check for overflow
                     if (ny .ge. maxyon) call error(720,ny,0.0)
                     ny = ny + 1
                     molyon(ny) = imol
                     do 4050 k = 1,3
                        py(k,ny) = pijk(k,ip)
                        ay(k,ny) = aijk(k,ip)
4050	               continue
                  end if
c
c write surface points for this probe position
                  do 4100 i = 1,np
                     write (4) n1(i),n2(i),n3(i),(s(k,i),k=1,3),area,yon
     $   (i)
4100	            continue
c end of ip loop
4150	         continue
c end of concave reentrant loop
4200	      continue
c
c saddle-shaped reentrant
            ishape = 2
c
c check for neither atom to be surfaces
            if (.not. (si .or. sj)) go to 5150
c
c special check for buried tori
c
c if both atoms are marked to be surface,
c but neither atom has any reentrant surface so far
c (after triangles with all katoms have been checked)
c and if there is some mutual neighbor in the same molecule
c close enough so that the torus cannot be free,
c then we know that this must be a buried torus
c
            if (si .and. sj .and. .not. srs(iatom) .and. .not. srs(jatom
     $   ) .and. mutual .gt. 0) then
               do 4250 knbr = 1,nnbr
                  if (.not. mnbr(knbr)) go to 4250
                  if (imol .ne. molnbr(knbr)) go to 4250
                  d2 = dist2(bij,cnbr(1,knbr))
                  rk2 = ernbr(knbr) ** 2 - hij ** 2
                  if (d2 .lt. rk2) go to 5150
4250	         continue
            end if
c calculate number of rotations of probe pair,
c rotation angle and rotation matrix
            rij = ri/(ri + rp) + rj/(rj + rp)
            avh = (abs(hij-rp) + hij * rij)/3
            nrot = sqrt(d(ptyp)) * pi * avh
            if (nrot .lt. 1) nrot = 1
            angle = 3.14159/nrot
c set up rotation matrix around x-axis
            call imatx(h)
            h(2,2) = cos(angle)
            h(3,3) = h(2,2)
            h(3,2) = sin(angle)
            h(2,3) = - h(3,2)
c calculate matrix to rotate x-axis onto iatom-jatom axis
            do 4300 k = 1,3
               g(k,1) = uij(k)
               g(k,2) = q(k)
               g(k,3) = t(k)
4300	      continue
c make the probe pair rotation matrix be about the iatom-jatom axis
            call conj(h,g,ghgt)
c
c arc generation
            do 4350 k = 1,3
               pij(k) = bij(k) + aij(k)
               vpi(k) = (ci(k) - pij(k)) * rp / (ri + rp)
               vpj(k) = (cj(k) - pij(k)) * rp / (rj + rp)
4350	      continue
c
c rotate circle onto iatom-jatom-probe plane
c and select points between probe-iatom and
c probe-jatom vector to form the arc
            narc = 1
            do 4500 i = 1,ncirc(ptyp)
               if (narc .gt. maxarc) call error(440,narc,0.0)
c rotation
               call multv(circle(1,i,ptyp),g,vps0(1,narc))
c if the vector is pointing away from the symmetry line
c the surface point cannot lie on the inward-facing arc
               if (dot(vps0(1,narc),aij) .gt. 0.0) go to 4500
               call cross(vpi,vps0(1,narc),vector)
               if (dot(g(1,3),vector) .lt. 0.0) go to 4500
               call cross(vps0(1,narc),vpj,vector)
               if (dot(g(1,3),vector) .lt. 0.0) go to 4500
c
c make arc point vectors originate with saddle circle center bij
c rather than probe center because they will be
c rotated around the iatom-jatom axis
               do 4400 k = 1,3
                  vbs0(k,narc,1) = vps0(k,narc) + aij(k)
4400	         continue
c invert arc through line of symmetry
               duij = dot(uij,vbs0(1,narc,1))
               do 4450 k = 1,3
                  vbs0(k,narc,2) = - vbs0(k,narc,1) + 2 * duij * uij(k)
4450	         continue
c
c check whether the arc point crosses the iatom-jatom axis
c and calculate the area associated with the point
               ht = dot(aij,vbs0(1,narc,1))/hij
               ayon(narc) = ht .lt. 0.0
               arca(narc) = (2*3.14159**2*rp*abs(ht))/(ncirc(ptyp)*nrot)
               narc = narc + 1
4500	      continue
            narc = narc - 1
c
c initialize power matrix to identity
            call imatx(pow)
c
c set knbr to zero for collision and buried checks
            knbr = 0
c rotate the probe pair around the pair of atoms
            do 5100 irot = 1,nrot
c multiply altitude vector by power matrix
               call multv(aij,pow,aijp(1,1))
c set up opposing altitude
               do 4550 k = 1,3
                  aijp(k,2) = - aijp(k,1)
4550	         continue
c set up probe sphere positions
               do 4650 ip = 1,2
                  do 4600 k = 1,3
                     pijp(k,ip) = bij(k) + aijp(k,ip)
4600	            continue
c check for collisions with neighboring atoms
                  pair(ip) = .not. collid(pijp(1,ip),rp,cnbr,ernbr,mnbr,
     $   nnbr,maxnbr,ishape,jnbr,knbr,molnbr,imol,lkf,lknbr)
4650	         continue
c no surface generation if neither probe position is allowed
               if (.not. pair(1) .and. .not. pair(2)) go to 5050
               both = pair(1) .and. pair(2)
c some reentrant surface for both atoms
               srs(iatom) = .true.
               srs(jatom) = .true.
c skip to bottom of middle loop if iatom and jatom
c are close enough and the surface point density is
c low enough so that the arc has no points
               if (narc .le. 0) go to 5050
c
c surface generation
               do 5000 ip = 1,2
                  if (.not. pair(ip)) go to 5000
c set default value for bury
                  bury = .false.
c don't check for probe collisions against other molecules
c unless we need to
                  if (ibury .gt. 0) bury = buried(pijp(1,ip),rp,cnbr,rnb
     $   r,mnbr,nnbr,maxnbr,ishape,jnbr,knbr,molnbr,imol)
c skip if not buried and buried surface only flag set
                  if (ibury .eq. 1 .and. .not. bury) go to 5000
c determine whether probe has surface on far side of line
                  yonprb = hij .lt. rp .and. .not. both
                  np = 1
c the saddle-shaped reentrant surface points come from the arc
                  do 4850 i = 1,narc
c overlapping reentrant surface removal
c for symmetry-related probe positions
                     if (both .and. ayon(i)) go to 4850
                     if (np .gt. maxppp) call error(480,np,0.0)
c rotate the arc from the xy plane onto the iatom-jatom-probe plane
                     call multv(vbs0(1,i,ip),pow,vbs)
c make coordinates relative to origin
                     do 4700 k = 1,3
                        s(k,np) = bij(k) + vbs(k)
4700	               continue
c find the closest atom and set up the atom numbers for the point
                     dsi = dist(s(1,np),ci) - ri
                     dsj = dist(s(1,np),cj) - rj
                     if (dsi .le. dsj) go to 4750
                     if (.not. sj) go to 4850
                     n1(np) = jatom
                     n2(np) = iatom
                     n3(np) = 0
                     go to 4800
4750	               continue
                     if (.not. si) go to 4850
                     n1(np) = iatom
                     n2(np) = jatom
                     n3(np) = 0
4800	               continue
c
c we've got a surface point
                     yon(np) = ayon(i)
                     torus(np) = arca(i)
                     np = np + 1
c end of arc point loop
4850	            continue
                  np = np - 1
                  if (np .le. 0) go to 5000
c
c write the molecule number, shape,number of points,
c probe position and the vector from the base to the probe center
                  write (4) imol,ishape,np,(pijp(k,ip),k=1,3),
     $   (aijp(k,ip),k=1,3),yonprb,bury
                  if (yonprb) then
c save probe in yon probe arrays
c check for overflow
                     if (ny .ge. maxyon) call error(720,ny,0.0)
                     ny = ny + 1
                     molyon(ny) = imol
                     do 4900 k = 1,3
                        py(k,ny) = pijp(k,ip)
                        ay(k,ny) = aijp(k,ip)
4900	               continue
                  end if
c
c write surface points for this probe position
                  do 4950 i = 1,np
                     write (4) n1(i),n2(i),n3(i),(s(k,i),k=1,3),
     $   torus(i),yon(i)
c end of arc point loop
4950	            continue
c end of probe pair loop
5000	         continue
c skip to here if both probe positions disallowed or no arc points
5050	         continue
c calculate new power matrix
               call cat(pow,ghgt)
c end of rotation loop
5100	      continue
c end of neighbor loop
5150	   continue
c
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
c skip to here if van der waals surface calculation
5200	   continue
c
c contact surface
         ishape = 1
c skip atom i if marked no surface requested
         if (.not. si) go to 5650
c
c if we are not calculating buried surface
c and the probe radius is greater than zero
c and iatom has at least one neighbor, but no reentrant surface,
c then iatom must be completely inaccessible to the probe
         if (.not. bury .and. rp .gt. 0.0 .and.
     $   nimol .gt. 0 .and. .not. srs(iatom)) go to 5650
c find the index into the atom type arrays for iatom
         do 5250 idx = 1,ntype
            if (itype(idx) .eq. iat(iatom)) go to 5300
5250	   continue
         stop 'logic error in ms regarding atom types'
5300	   continue
         area = (4 * pi * ri ** 2) / nua(idx)
c set jnbr, knbr to zero for collision, buried checks
         jnbr = 0
         knbr = 0
c
c contact probe placement loop
         do 5600 i = 1,nua(idx)
c set up probe coordinates
            do 5350 k = 1,3
               pipt(k) = ci(k) + eva(k,i,idx)
5350	      continue
c check for collision with neighboring atoms
            if (collid(pipt,rp,cnbr,ernbr,mnbr,nnbr,maxnbr,ishape,
     $   jnbr,knbr,molnbr,imol,lkf,lknbr)) go to 5600
c go write it out if we don't care about buried surface
            if (ibury .eq. 0) then
               ib = 0
               go to 5400
            end if
            bury = buried(pipt,rp,cnbr,rnbr,mnbr,nnbr,maxnbr,ishape,
     $   jnbr,knbr,molnbr,imol)
            if (ibury .eq. 1 .and. .not. bury) go to 5600
            if (bury) then
               ib = 1
            else
               ib = 0
            end if
c
5400	      continue
c increment surface point counter for convex surface
            nshape(1) = nshape(1) + 1
c add surface point area to contact area
            areac = areac + area
            iatnum(1) = iatom
            iatnum(2) = 0
            iatnum(3) = 0
            do 5450 k = 1,3
               outco(k) = ci(k) + ri * ua(k,i,idx)
               outvec(k) = ua(k,i,idx)
5450	      continue
c four different output formats
            if (.not. binary .and. .not. short) then
               write (7,5500) iatnum,ishape,outco,area,outvec,ib
5500	         format(3i5,i2,3f9.3,4f7.3,i2)
            else if (binary .and. .not. short) then
               write (7) iatnum,ishape,outco,area,outvec,ib
            else if (.not. binary .and. short) then
c	         write (7,5550) iatnum,ishape,outco,ib
c   ---dac change to write records for DeMon dummy atoms:
             write (7,'(a5,4f8.3)') 'H    ',outco,0.0
5550	         format(3i5,i2,3f9.3,i2)
            else if (binary .and. short) then
               write (7) iatnum,ishape,outco,ib
            end if
c
c end of nua loop
5600	   continue
c end of iatom loop
5650	continue
c
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c write out messages
c
c for array dimensioning
      write(6,5700) maxnb
5700	format(1x,i5,' neighbors maximum')
c
c close contact and before files
      close(7)
c if van der waals surface we are finished
      if (rp .eq. 0.0) go to 8000
      close(4)
c
c orsr    orsr    orsr    orsr    orsr    orsr    orsr    orsr    orsr
c
c
c overlapping reentrant surface removal
c for non-symmetry-related probes
c probe diameter
      dp = 2 * rp
c diameter squared
      dp2 = dp ** 2
c radius squared
      rp2 = rp ** 2
c width for cubing algorithm
c	width = dp
c
c     set up cube arrays
c     first the integer coordinate arrays
      do 5850 iy = 1,ny
         do 5800 k = 1,3
            ico(k,iy) = (py(k,iy)-comin(k)-radmax-rp)/width + 1
            if (ico(k,iy) .lt. 1) ico(k,iy) = 1
            if (ico(k,iy) .gt. maxcub) stop 'cube coordinate too large'
5800	   continue
5850	continue
c
c initialize head pointer array
      do 6000 k = 1,maxcub
         do 5950 j = 1,maxcub
            do 5900 i = 1,maxcub
               icube(i,j,k) = 0
5900	      continue
5950	   continue
6000	continue
c
c initialize linked list pointers
      do 6050 iy = 1,maxyon
         icuptr(iy) = 0
6050	continue
c
c set up head and later pointers for each yon probe
      do 6200 iy = 1,ny
c skip atoms with surface request numbers of zero
         i = ico(1,iy)
         j = ico(2,iy)
         k = ico(3,iy)
         if (icube(i,j,k) .le. 0) then
c     first atom in this cube
            icube(i,j,k) = iy
         else
c     add to end of linked list
            iptr = icube(i,j,k)
6100	      continue
            if (icuptr(iptr) .le. 0) go to 6150
            iptr = icuptr(iptr)
            go to 6100
6150	      continue
            icuptr(iptr) = iy
         end if
6200	continue
c
c reopen before file for reading
      open(4,file='before',form='unformatted',status='old')
      rewind(4)
c
c first pass
c gather victim probes
      nv = 1
c no victim probes if no yon probes
      if (ny .le. 0) go to 6950
      rewind(4)
c initialize victim hashing array
      do 6250 j = 1,maxvic
         ivic(j) = 0
6250	continue
c initialize index of free slot last used
      ifrlst = 0
c initialize probe record number
      i = 1
6300	continue
c check for victim overflow
      if (nv .gt. maxvic) call error(760,nv,0.0)
c read reentrant probe and points
      read (4,end=6950) molvic(nv),ishape,np,
     $   (pv(k,nv),k=1,3),(av(k,nv),k=1,3),yonprb,bury
      do 6350 j = 1,np
         read (4) n1(1),n2(1),n3(1),(s(k,1),k=1,3),area,yon(1)
6350	continue
      if (yonprb) go to 6900
c check if probe too far from symmetry element for possible overlap
      if (anorm(av(1,nv)) .gt. dp) go to 6900
c
c look for overlap with any yon probe in the same molecule
c use cubing algorithm to save time
c
c calculate which cube this probe lies in
      ici = (pv(1,nv)-comin(1)-radmax-rp)/width + 1
      if (ici .lt. 1) ici = 1
      if (ici .gt. maxcub) stop 'cube coordinate too large'
      icj = (pv(2,nv)-comin(2)-radmax-rp)/width + 1
      if (icj .lt. 1) icj = 1
      if (icj .gt. maxcub) stop 'cube coordinate too large'
      ick = (pv(3,nv)-comin(3)-radmax-rp)/width + 1
      if (ick .lt. 1) ick = 1
      if (ick .gt. maxcub) stop 'cube coordinate too large'
c check for overlap with probes in adjoining cubes
      do 6850 jck = ick-1,ick+1
         if (jck .lt. 1 .or. jck .gt. maxcub) go to 6850
         do 6800 jcj = icj-1,icj+1
            if (jcj .lt. 1 .or. jcj .gt. maxcub) go to 6800
            do 6750 jci = ici-1,ici+1
               if (jci .lt. 1 .or. jci .gt. maxcub) go to 6750
               jp = icube(jci,jcj,jck)
c
c
6400	         continue
               if (jp .le. 0) go to 6700
               if (molyon(jp) .ne. molvic(nv)) go to 6650
               x = abs(py(1,jp) - pv(1,nv))
               if (x .ge. dp) go to 6650
               y = abs(py(2,jp) - pv(2,nv))
               if (y .ge. dp) go to 6650
               z = abs(py(3,jp) - pv(3,nv))
               if (z .ge. dp) go to 6650
               d2 = x ** 2 + y ** 2 + z ** 2
               if (d2 .ge. dp2) go to 6650
c check that probes face each other
               if (dot(ay(1,jp),av(1,nv)) .ge. 0.0) go to 6650
c new victim probe
c put into hashing table
               ihash = mod(i,maxvic) + 1
               if (ivic(ihash) .eq. 0) then
c empty slot
                  ivic(ihash) = i
                  ivicp(ihash) = 0
               else
                  iprev = ihash
                  iptr = ivicp(ihash)
6450	            continue
c check for end of linked list
                  if (iptr .eq. 0) go to 6500
                  iprev = iptr
                  iptr = ivicp(iptr)
                  go to 6450
6500	            continue
c look for a free slot
                  do 6550 ifree = ifrlst+1,maxvic
                     if (ivic(ifree) .eq. 0) go to 6600
6550	            continue
                  stop 'victim oveflow'
6600	            continue
c store record number in free slot
                  ivic(ifree) = i
                  ivicp(iprev) = ifree
                  ivicp(ifree) = 0
c new index to last free slot used
                  ifrlst = ifree
               end if
               nv = nv + 1
c one overlap makes this probe a victim
c we don't need to check any more
               go to 6900
6650	         continue
               jp = icuptr(jp)
               go to 6400
6700	         continue
6750	      continue
6800	   continue
6850	continue
c end of yon probe loop
c skip to here if finished with hunt for overlapping probes
6900	continue
      i = i + 1
      go to 6300
c skip to here if there are no yon probes and hence no victims
6950	continue
      nv = nv - 1
c
c open reentrant file for writing

      if (binary) then
         open(8,file='reentrant',form='unformatted',status='unknown')
      else
         open(8,file='reentrant',status='unknown')
      end if

      rewind(8)
c
c second pass
c read, check and write surface points
      rewind(4)
      i = 1
7000	continue
      read (4,end=7850) imol,ishape,np,
     $   (p(k),k=1,3),(a(k),k=1,3),yonprb,bury
c no points can be eaten if this probe is neither yon nor a victim
      neat = 0
      nyeat = 0
      if (ny .le. 0) go to 7450
      ipt = 0
c determine if probe is a yon or victim probe
      if (.not. yonprb) go to 7050
c we've got a yon probe here
      ipt = 2
      go to 7200
7050	continue
      if (nv .le. 0) go to 7450
c hash into table of victim probes
      iptr = mod(i,maxvic) + 1
7100	continue
      if (iptr .eq. 0) go to 7200
      if (ivic(iptr) .eq. 0) go to 7200
      if (ivic(iptr) .eq. i) go to 7150
      iptr = ivicp(iptr)
      go to 7100
7150	continue
c we've got a victim
      ipt = 1
7200	continue
c
      if (ipt .le. 0) go to 7450
c check this victim or yon probe against all yon probes
      do 7300 j = 1,ny
         if (imol .ne. molyon(j)) go to 7300
         if (dist2(p,py(1,j)) .ge. dp2) go to 7300
         if (dot(a,ay(1,j)) .ge. 0.0) go to 7300
c this yon probe could eat some of the probe's points
         neat = neat + 1
         nyeat = nyeat + 1
         if (neat .gt. maxeat) call error(830,neat,0.0)
         do 7250 k = 1,3
            eat(k,neat) = py(k,j)
7250	   continue
c end of yon probe loop
7300	continue
c
c only yon probes can have their points eaten by victims
      if (ipt .le. 1) go to 7450
c check this yon probe against all victim probes
      do 7400 j = 1,nv
         if (imol .ne. molvic(j)) go to 7400
         if (dist2(p,pv(1,j)) .ge. dp2) go to 7400
         if (dot(a,av(1,j)) .ge. 0.0) go to 7400
c this victim probe could eat some of the probe's points
         neat = neat + 1
         if (neat .gt. maxeat) call error(850,neat,0.0)
         do 7350 k = 1,3
            eat(k,neat) = pv(k,j)
7350	   continue
c end of victim probe loop
7400	continue
c
c skip to here if victim or both probe overlap checks omitted
7450	continue
c
c read the surface points belonging to the probe
      do 7750 j = 1,np
         read (4) n1(1),n2(1),n3(1),(s(k,1),k=1,3),area,yon(1)
         if (neat .le. 0) go to 7550
c check surface point against all eaters of this probe
         do 7500 k = 1,neat
c victim probes cannot eat non-yon points of yon probes
            if (yonprb .and. .not. yon(1) .and. k .gt. nyeat) go to 7500
            if (dist2(eat(1,k),s(1,1)) .lt. rp2) go to 7700
7500	   continue
c skip to here if no overlapping probes could eat this point
7550	   continue
         do 7600 k = 1,3
            outvec(k) = (p(k) - s(k,1))/rp
7600	   continue
c reentrant surface point
         nshape(ishape) = nshape(ishape) + 1
         arear = arear + area
c mark whether buried
         ib = 0
         if (ibury .gt. 0 .and. bury) ib = 1
c four possible output formats
         iatnum(1) = n1(1)
         iatnum(2) = n2(1)
         iatnum(3) = n3(1)
         do 7650 k = 1,3
            outco(k) = s(k,1)
7650	   continue
         if (.not. binary .and. .not. short) then
            write (8,5500) iatnum,ishape,outco,area,outvec,ib
         else if (binary .and. .not. short) then
            write (8) iatnum,ishape,outco,area,outvec,ib
         else if (.not. binary .and. short) then
            write (8,5550) iatnum,ishape,outco,ib
         else if (binary .and. short) then
            write (8) iatnum,ishape,outco,ib
         end if
         go to 7750
7700	   continue
         nlost(ishape) = nlost(ishape) + 1
c end of np loop
7750	continue
c end of i loop
7800	continue
      i = i + 1
      go to 7000
7850	continue
c close files
      close(4)
      close(8)
c messages
      write (6,7900) ny,nv
7900	format(1x,i5,' yon and ',i5,' victim probes')
      write(6,7950) nlost(2),nlost(3)
7950	format(1x,i5,' saddle and ',i5,
     $   ' concave surface points removed during non-symmetry orsr')
8000	continue
c write out how many points
      write (6,8050) nshape(1),nshape(2),nshape(3)
8050	format(1x,i5,' contact and ',i5,' saddle and ',
     $   i5,' concave surface points')
      write (6,8100) nshape(1) + nshape(2) + nshape(3)
8100	format(1x,i8,' total surface points')
      write (6,8150) areac,arear,areac+arear
8150	format(1x,'contact area:',f10.3,2x,'reentrant area:',f10.3,
     $   2x,'total area:',f10.3)
      stop
      end
c
c
c subroutines and functions
c
c
c general vector and matrix routines
c
      function dist(a,b)
c distance between a and b
      real a(3)
      real b(3)
      dist = sqrt((a(1)-b(1))**2 + (a(2)-b(2))**2 + (a(3)-b(3))**2)
      return
      end
c
      function dist2(a,b)
c distance between a and b squared
      real a(3)
      real b(3)
      dist2 = (a(1)-b(1))**2 + (a(2)-b(2))**2 + (a(3)-b(3))**2
      return
      end
c
      function anorm(a)
c norm of a
      real a(3)
      anorm = sqrt(a(1)**2 + a(2)**2 + a(3)**2)
      return
      end
c
      function dot(a,b)
c dot product
      real a(3)
      real b(3)
      dot = a(1)*b(1) + a(2)*b(2) + a(3)*b(3)
         return
         end
c
      subroutine cross(a,b,c)
c cross product
      real a(3)
      real b(3)
      real c(3)
      c(1) = a(2) * b(3) - a(3) * b(2)
      c(2) = a(3) * b(1) - a(1) * b(3)
      c(3) = a(1) * b(2) - a(2) * b(1)
      return
      end
c
      subroutine multv(v,a,w)
c multiply v by a giving w
      real a(3,3)
      real v(3)
      real w(3)
      do 50 i = 1, 3
         w(i) = a(i,1)*v(1) + a(i,2)*v(2) + a(i,3)*v(3)
50	continue
      return
      end
c
      subroutine vnorm(a,b)
c normalize a giving b
      real a(3),b(3)
      v = anorm(a)
      do 50 k = 1,3
         b(k) = a(k) / v
50	continue
      return
      end
c
      subroutine vperp(a,b)
c return b perpendicular to a
      real a(3)
      real b(3)
      real p(3)
c find smallest component
      small = 10000.0
      m = 0
      do 50 k = 1,3
         if (abs(a(k)) .ge. small) go to 50
         small = abs(a(k))
         m = k
50	continue
      do 100 k = 1,3
         b(k) = 0.0
         if (k .eq. m) b(k) = 1.0
100	continue
c take projection along a
      dt = a(m) / (a(1)**2 + a(2)**2 + a(3)**2)
      do 150 k = 1, 3
         p(k) = dt * a(k)
c subtract projection from b
         b(k) = b(k) - p(k)
150	continue
c renormalize b
      call vnorm(b,b)
      return
      end
c
      subroutine cat(a,b)
c concatenate matrix b into matrix a
      real a(3,3)
      real b(3,3)
      real temp(3,3)
      do 100 i = 1,3
         do 50 j = 1,3
            temp(i,j) = a(i,1)*b(1,j) + a(i,2)*b(2,j) + a(i,3)*b(3,j)
50	   continue
100	continue
      do 200 i = 1,3
         do 150 j = 1,3
            a(i,j) = temp(i,j)
150	   continue
200	continue
      return
      end
c
      subroutine conj(h,g,ghgt)
c conjugate matrix g with matrix h giving ghgt
      real g(3,3)
      real h(3,3)
      real ghgt(3,3)
      real gt(3,3)
c initialize ghgt matrix to identity
c concatenate g h gt
      call imatx(ghgt)
      call cat(ghgt,g)
      call cat(ghgt,h)
c calculate gt
      do 100 k = 1,3
         do 50 l = 1,3
            gt(k,l) = g(l,k)
50	   continue
100	continue
      call cat(ghgt,gt)
      return
      end
c
      subroutine imatx(a)
c load identity matrix
      real a(3,3)
      do 100 i = 1,3
         do 50 j = 1,3
            a(i,j) = 0.0
50	   continue
         a(i,i) = 1.0
100	continue
      return
      end
c
      function det(a,b,c)
c return triple product of the three vectors
      real a(3)
      real b(3)
      real c(3)
      real ab(3)
      call cross(a,b,ab)
      det = dot(ab,c)
      return
      end
c
c geometric routines
c
c
      logical function collid(p,rp,cnbr,ernbr,mnbr,nnbr,maxnbr,ishape,
     $   jnbr,knbr,molnbr,imol,lkf,lknbr)
c collision check of probe with neighboring atoms
c belonging to the same molecule
cMAS
      common/mpck/colpre(3),radpre,r2pre
c	real	colpre(3)
c	real 	radpre,r2pre
cMAS
      real p(3)
      real cnbr(3,maxnbr)
      real ernbr(maxnbr)
      logical mnbr(maxnbr)
      integer molnbr(maxnbr)
      integer lknbr(maxnbr)
      integer imol,ishape
c
cMAS check whether probe is too close to any neighbor
c
c	check neighbor from previous collision first
      vect1 = abs(p(1) - colpre(1))
      vect2 = abs(p(2) - colpre(2))
      vect3 = abs(p(3) - colpre(3))
      dd2 = vect1 ** 2 + vect2 ** 2 + vect3 ** 2
      if(dd2 .lt. r2pre) then
      	collid = .true.
      	return
      endif
c
cMAS
c
      i = lkf
      go to 100
50	continue
      i = lknbr(i)
100	continue
      if (i .eq. 0) go to 150
      vect1 = abs(p(1) - cnbr(1,i))
      if (vect1 .ge. ernbr(i)) go to 50
      vect2 = abs(p(2) - cnbr(2,i))
      if (vect2 .ge. ernbr(i)) go to 50
      vect3 = abs(p(3) - cnbr(3,i))
      if (vect3 .ge. ernbr(i)) go to 50
      if (i .eq. jnbr .or. i .eq. knbr) go to 50
      sr2 = ernbr(i) ** 2
      dd2 = vect1 ** 2 + vect2 ** 2 + vect3 ** 2
      if (dd2 .ge. sr2) go to 50
      collid = .true.
c
cMAS	found a collision
c
      colpre(1) = cnbr(1,i)
      colpre(2) = cnbr(2,i)
      colpre(3) = cnbr(3,i)
c	radpre = ernbr(i)
      r2pre = sr2
cMAS
      return
150	continue
      collid = .false.
      return
      end
c
      logical function buried(p,rp,cnbr,rnbr,mnbr,nnbr,maxnbr,ishape,
     $   jnbr,knbr,molnbr,imol)
c collision check of probe with neighboring atoms
c belonging to a different molecule
      real p(3)
      real cnbr(3,maxnbr)
      real rnbr(maxnbr)
      logical mnbr(maxnbr)
      integer molnbr(maxnbr)
      integer imol,ishape
c
      if (nnbr .le. 0) go to 100
c check whether probe is too close to any neighbor
      do 50 i = 1, nnbr
         if (imol .eq. molnbr(i)) go to 50
         if (ishape .gt. 1 .and. i .eq. jnbr) go to 50
         if (ishape .eq. 3 .and. (i .eq. knbr .or. .not. mnbr(i)))
     $   go to 50
         sumrad = rp + rnbr(i)
         vect1 = abs(p(1) - cnbr(1,i))
         if (vect1 .ge. sumrad) go to 50
         vect2 = abs(p(2) - cnbr(2,i))
         if (vect2 .ge. sumrad) go to 50
         vect3 = abs(p(3) - cnbr(3,i))
         if (vect3 .ge. sumrad) go to 50
         sr2 = sumrad ** 2
         dd2 = vect1 ** 2 + vect2 ** 2 + vect3 ** 2
         if (dd2 .lt. sr2) go to 150
50	continue
100	continue
      buried = .false.
      go to 200
150	continue
      buried = .true.
200	continue
      return
      end
c
      subroutine genun(u,n)
c generate unit vectors over sphere
      real u(3,n)
      nequat = sqrt(n * 3.14159)
      nvert = 0.5 * nequat
      if (nvert .lt. 1) nvert = 1
      nu = 0
      do 100 i = 0,nvert
         fi = (3.14159 * i) / nvert
         z = cos(fi)
         xy = sin(fi)
         nhor = nequat * xy
         if (nhor .lt. 1) nhor = 1
         do 50 j = 0,nhor-1
            fj = (2 * 3.14159 * j) / nhor
            x = cos(fj) * xy
            y = sin(fj) * xy
            if (nu .ge. n) go to 150
            nu = nu + 1
            u(1,nu) = x
            u(2,nu) = y
            u(3,nu) = z
50	   continue
100	continue
150	continue
      n = nu
      return
      end
c
c error message subroutine
c
      subroutine error(number,int,float)
c
      integer list(15)
      data list/120,127,130,140,150,160,170,
     $   210,320,440,480,720,760,830,850/
c
      do 50 i = 1,15
         if (list(i) .eq. number) go to 150
50	continue
      write (6,100)
100	format(1x,'error of unidentifiable type')
      stop
150	continue
c
      go to (200,300,400,500,600,700,800,900,1000,1100,1200,
     $   1300,1400,1500,1600) i
c
200	write (6,250) number,float
250	format(1x,'error',i5,2x,'negative probe radius:',f10.5)
      stop
300	write (6,350) number,int
350	format(1x,'error',i5,2x,'bad buried surface flag:',i5)
      stop
400	write (6,450) number,int
450	format(1x,'error',i5,2x,'too few or too many atom types:',i5)
      stop
500	write (6,550) number,float,int
550	format(1x,'error',i5,2x,'negative atom radius:',
     $   f10.5,' atom',i5)
      stop
600	write (6,650) number,int
650	format(1x,'error',i5,2x,'too many atoms:',i5)
      stop
700	write (6,750) number,int
750	format(1x,'error',i5,2x,
     $   'invalid surface request number for atom:',i5)
      stop
800	write (6,850) number,int
850	format(1x,'error',i5,2x,'invalid atom type for atom:',i5)
      stop
900	write (6,950) number,int
950	format(1x,'error',i5,2x,'too many neighbors:',i5)
      stop
1000	write (6,1050) number,int
1050	format(1x,'error',i5,2x,'too many points for reentrant probe:',
     $   i5)
      stop
1100	write (6,1150) number,int
1150	format(1x,'error',i5,2x,'too many points for arc:',i5)
      stop
1200	write (6,1250) number,int
1250	format(1x,'error',i5,2x,'too many points for reentrant probe:',
     $   i5)
      stop
1300	write (6,1350) number,int
1350	format(1x,'error',i5,2x,'too many yon probes:',i5)
      stop
1400	write (6,1450) number,int
1450	format(1x,'error',i5,2x,'too many victim probes:',i5)
      stop
1500	write (6,1550) number,int
1550	format(1x,'error',i5,2x,'too many eaters:',i5)
      stop
1600	write (6,1650) number,int
1650	format(1x,'error',i5,2x,'too many eaters:',i5)
      stop
      end
