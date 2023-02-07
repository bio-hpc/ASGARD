#include "../include/dprec.fh"

program TestAmberMask
! Test the syntax and correctness of amber FIND command (from group input).
! VH (Dec 30, 2003)
!
! Command line syntax:
!    ambmask -p prmtop -c inpcrd -prnlev [0|1|2|3] -out [format] -find 'maskstr'
!       prmtop  : (new) amber topology
!       inpcrd  : coordinates
!       format  : the format of output for selected atoms
!                 short .. condensed output with ':' and '@' specifiers
!                 pdb   .. amber-style pdb format
!                 amber .. amber "group input" format (ATOM/RES cards) 
!       maskstr : mask string expression (put it in quotes '..' if spaces
!                 or special shell characters are present)
!       prnlev  : the amount of (debugging) info printed:
!                 0 .. only prints selected atoms in condensed list format
!                 1 .. prints original, tokenized, and RPN form of maskstring,
!                      # of atoms/residues, # of selected atoms, and selected
!                      atoms in amber-style pdb format (this is the default)
!                 2 .. in addition prints mask array after eval() routine,
!                 3 .. and all mask arrays pushed/popped from/into the stack
!
! Notes:
! - Code for reading topology and coordinates was derived from ambpdb.f
! - If no command line arguments are given, program usage is printed out
! - This program is just a driver for mask string parsing. The real work
!   of converting 'maskstr' to mask array is done by the routine findmask()
!   in ../sander/findmask.f
! TODO:
!   complete the documentation accordingly
   use findmask
   
   implicit none
   
   character(80) arg,prmtop,inpcrd,maskstr,title
   character(20) prnlevstr, outformat
   integer iarg, nargc, iargc, i

   integer natom, nres, prnlev
   integer astat, ios
   character(len=4) ititl(20)
   character(len=4), dimension(:), allocatable :: igraph, isymbl, lbres
   integer, dimension(:), allocatable :: ipres
   _REAL_, dimension(:), allocatable :: crd, chg
   integer, dimension(:), allocatable :: mask
   
   ! check argument options
   prmtop  = 'prmtop'  ! default topology file name
   inpcrd  = 'inpcrd'  ! default coordinate file name
   prnlev  = 1         ! default amount of (debugging) information printed
   outformat = 'pdb'   ! default output format is pdb
   maskstr = ':*'      ! if no 'mask' specified, default is everything

   iarg = 0
   nargc = iargc()     ! get number of command line arguments

   if (nargc.eq.0) then
      ! if no cmd line arguments given, print command usage
      write(*,'(A)') "ambmask -p prmtop -c inpcrd -prnlev [0|1|2|3] " //  & 
                "-out [format] -find maskstr"
      write(*,'("   prmtop  : amber topology")')
      write(*,'("   inpcrd  : amber (restrt) coordinates")')
      write(*,'("   prnlev  : amount of (debugging) info printed")')
      write(*,'("   format  : output format: short|pdb|amber")')
      write(*,'("   maskstr : mask string expression (put it in quotes")')
      write(*,'("             if spaces or special shell chars present)"/)')
      stop
   end if

10 continue
      iarg = iarg + 1
      call getarg(iarg,arg)
      if (arg .eq. '-p') then
         iarg = iarg + 1
         call getarg(iarg,prmtop)
      else if (arg .eq. '-c') then
         iarg = iarg + 1
         call getarg(iarg,inpcrd)
      else if (arg .eq. '-prnlev') then
         iarg = iarg + 1
         call getarg(iarg,prnlevstr)
         read(prnlevstr,*,iostat=ios) prnlev
         if (ios > 0) stop "error reading -prnlev"
      else if (arg .eq. '-out') then
         iarg = iarg + 1
         call getarg(iarg,outformat)
      else if (arg .eq. '-find') then
         iarg = iarg + 1
         call getarg(iarg,maskstr)
      else
         write(*,'("Unknown flag: ",A)') arg
         stop
      endif
   if (iarg .lt. nargc) goto 10

   ! open the parm file first to find out array dimensions
   call amopen(10,prmtop,'O','F','R')
   call getdim(natom,nres,ititl,10)
   close(unit=10)

   if (prnlev >= 1) then
      if (natom < 1000000) then
        write(*,'("natom = ",I6)') natom
        write(*,'("nres  = ",I6)') nres
      else
        write(*,'("natom = ",I9)') natom
        write(*,'("nres  = ",I9)') nres
      end if
   end if
   
   ! allocate memory for the topology arrays that we need:
   allocate(igraph(natom), stat=astat)
   if (astat/=0) stop "***igraph allocation error***"

   allocate(isymbl(natom), stat=astat)
   if (astat/=0) stop "***isymbl allocation error***"

   allocate(lbres(nres), stat=astat)
   if (astat/=0) stop "***lbres allocation error***"

   allocate(ipres(nres+1), stat=astat)
   if (astat/=0) stop "***ipres allocation error***"

   allocate(crd(3*natom), stat=astat)
   if (astat/=0) stop "***crd allocation error***"

   allocate(chg(natom), stat=astat)
   if (astat/=0) stop "***chg allocation error***"

   ! open the parm file and load the necessary stuff
   call amopen(10,prmtop,'O','F','R')
   call getnam(natom,nres,igraph,isymbl,ipres,lbres,ititl,10,chg)
   close(unit=10)

   ! read the coordinate file
   call amopen(10,inpcrd,'O','F','R')
   call getcor(natom,crd,10)
   close(unit=10)
   
   allocate(mask(natom), stat=astat)
   if (astat/=0) stop "***mask allocation error***"
   ! initialize 'mask' to no atoms selected
   do i=1,natom
      mask(i) = 0
   end do

   ! character array 'mask' is returned according to mask string 'maskstr'
   call atommask(natom,nres,prnlev,igraph,isymbl,ipres,  &
                 lbres,crd,maskstr,mask)

   ! output the selected atoms in (condensed list) or (pdb) format
   if (index(outformat,"short") > 0) then
      call genatlist(nres,igraph,ipres,6,mask)
   else if (index(outformat,"pdb") > 0) then
      title = maskstr  ! print 'maskstr' into a remark field in pdb
      call genpdb(nres,igraph,isymbl,ipres,lbres,6,crd,title,mask)
   else if (index(outformat,"amber") > 0) then
      call gengroupinput(nres,ipres,mask)
   else
      stop "unknown -format parameter"
   end if

   deallocate(mask)
   ! deallocate topology arrays:
   deallocate(ipres)
   deallocate(lbres)
   deallocate(igraph)
   deallocate(isymbl)
   deallocate(crd)
   deallocate(chg)
   
end program TestAmberMask


!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ get number of atoms/residues from topology (for dynamic allocation)
subroutine getdim(natom,nres,ititl,nf)
   implicit none

   integer natom, nres, nf
   character(len=4) ititl(*)

   integer i, iok
   integer ntypes,nbonh,mbona,ntheth,mtheta,nphih,mphia
   integer nhparm,nparm,nnb,nbona,ntheta,nphia
   integer numbnd,numang,nptra,natyp,nphb,ifpert,nbper,ngper
   integer ndper,mbper,mgper,mdper,ifbox,nmxrs,ifcap
   character(80) fmt,fmtin,ifmt,afmt,rfmt,type

   ifmt = '(12I6)'
   afmt = '(20A4)'
   rfmt = '(5E16.8)'

   ! read title and control variables from the beginning of amber topology
   fmtin = afmt
   type = 'TITLE'
   call nxtsec(nf,  0,  0,fmtin,  type,  fmt,  iok)
   read(nf,fmt) (ititl(i),i=1,20)

   fmtin = ifmt
   type = 'POINTERS'
   call nxtsec(nf,  6,  0,fmtin,  type,  fmt,  iok)
   read(nf,fmt) natom,ntypes,nbonh,mbona,ntheth,mtheta,nphih,mphia,  &
                nhparm,nparm,nnb,nres,nbona,ntheta,nphia,            &
                numbnd,numang,nptra,natyp,nphb,ifpert,nbper,ngper,   &
                ndper,mbper,mgper,mdper,ifbox,nmxrs,ifcap

end subroutine getdim


!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ get atom names, residues names, amber types, etc.
subroutine getnam(natom,nres,igraph,isymbl,ipres,lbres,ititl,nf,chg)
   use constants, only : INV_AMBER_ELECTROSTATIC
   use parms, only : nttyp
   implicit none
   
   integer natom, nres, nf, ipres(*)
   character(len=4) igraph(*), isymbl(*), lbres(*)
   _REAL_ chg(*)
   character(len=4) ititl(*)

   integer i, iok
   integer ntypes,nbonh,mbona,ntheth,mtheta,nphih,mphia
   integer nhparm,nparm,nnb,nbona,ntheta,nphia
   integer numbnd,numang,nptra,natyp,nphb,ifpert,nbper,ngper
   integer ndper,mbper,mgper,mdper,ifbox,nmxrs,ifcap
   character(80) fmt,fmtin,ifmt,afmt,rfmt,type
   ! the following are only needed for reading old-style topology
   integer i6dmax, edmax, idmax, ntype, astat
   integer, dimension(:), allocatable :: i6d, id,jd,kd,ld,ic
   _REAL_, dimension(:), allocatable :: ed

   ifmt = '(12I6)'
   afmt = '(20A4)'
   rfmt = '(5E16.8)'

   ! read needed sections from amber topology

   fmtin = afmt
   type = 'TITLE'
   call nxtsec(nf,  0,  0,fmtin,  type,  fmt,  iok)
   read(nf,fmt) (ititl(i), i=1,20)

   fmtin = ifmt
   type = 'POINTERS'
   call nxtsec(nf,  6,  0,fmtin,  type,  fmt,  iok)
   read(nf,fmt) natom,ntypes,nbonh,mbona,ntheth,mtheta,nphih,mphia,  &
                nhparm,nparm,nnb,nres,nbona,ntheta,nphia,            &
                numbnd,numang,nptra,natyp,nphb,ifpert,nbper,ngper,   &
                ndper,mbper,mgper,mdper,ifbox,nmxrs,ifcap

   ntype = ntypes * ntypes

   ! read atom names, and the charges
   fmtin = afmt
   type = 'ATOM_NAME'
   call nxtsec(nf,  6,  0,fmtin,  type,  fmt,  iok)
   read(nf,fmt) (igraph(i), i=1,natom)

   fmtin = rfmt
   type = 'CHARGE'
   call nxtsec(nf,  6,  0,fmtin,  type,  fmt,  iok)
   read(nf,fmt) (chg(i), i=1,natom)
   chg(1:natom) = chg(1:natom)*INV_AMBER_ELECTROSTATIC

   if( iok.eq.-1 ) then   ! this is an old-style prmtop file
      ! old style prmtop has to be read sequentially, even though
      ! we only need 'lbres', 'ipres', and 'isymbl'
      ! allocate dummy arrays needed for reading of old prmtop
      i6dmax = max(natom, ntype, nres, nnb)
      edmax = max(natom, numbnd, numang, nptra, natyp, nttyp, nphb)
      idmax = max(nbonh, nbona, ntheth, ntheta, nphih, nphia)
      ! debug info
      allocate(i6d(i6dmax), stat=astat)
      if (astat /= 0) stop "***memory allocation error***"

      allocate(ed(edmax), stat=astat)
      if (astat /= 0) stop "***memory allocation error***"

      allocate(id(idmax),jd(idmax),kd(idmax),ld(idmax),ic(idmax), stat=astat)
      if (astat /= 0) stop "***memory allocation error***"
      
      read(nf,rfmt) (ed(i), i=1,natom)      ! amass
      read(nf,ifmt) (i6d(i), i=1,natom)     ! iac
      read(nf,ifmt) (i6d(i), i=1,natom)     ! numex
      read(nf,ifmt) (i6d(i), i=1,ntype)     ! ico
      read(nf,afmt) (lbres(i), i=1,nres)    ! lbres
      read(nf,ifmt) (ipres(i), i=1,nres)    ! ipres
      ipres(nres+1) = natom + 1
      ! read a lot of dummy stuff before you get isymbl
      read(nf,rfmt) (ed(i), i=1,numbnd)     ! rq
      read(nf,rfmt) (ed(i), i=1,numbnd)     ! req
      read(nf,rfmt) (ed(i), i=1,numang)     ! tk
      read(nf,rfmt) (ed(i), i=1,numang)     ! teq
      read(nf,rfmt) (ed(i), i=1,nptra)      ! pk
      read(nf,rfmt) (ed(i), i=1,nptra)      ! pn
      read(nf,rfmt) (ed(i), i=1,nptra)      ! phase
      read(nf,rfmt) (ed(i), i=1,natyp)      ! solty
      read(nf,rfmt) (ed(i), i=1,nttyp)      ! cn1
      read(nf,rfmt) (ed(i), i=1,nttyp)      ! cn2
      read(nf,ifmt) (id(i),jd(i),ic(i), i=1,nbonh)    ! ibh,jbh,icbh
      read(nf,ifmt) (id(i),jd(i),ic(i), i=1,nbona)    ! ib,jb,icb
      ! ith,jth,kth,icth; and  it,jt,kt,ict
      read(nf,ifmt) (id(i),jd(i),kd(i),ic(i), i=1,ntheth)
      read(nf,ifmt) (id(i),jd(i),kd(i),ic(i), i=1,ntheta)
      ! iph,jph,kph,lph,icph; and  ip,jp,kp,lp,icp
      read(nf,ifmt) (id(i),jd(i),kd(i),ld(i),ic(i), i=1,nphih)
      read(nf,ifmt) (id(i),jd(i),kd(i),ld(i),ic(i), i=1,nphia)
      read(nf,ifmt) (i6d(i), i=1,nnb)       ! natex ('next' is 'nnb')
      read(nf,rfmt) (ed(i), i=1,nphb)       ! asol
      read(nf,rfmt) (ed(i), i=1,nphb)       ! bsol
      read(nf,rfmt) (ed(i), i=1,nphb)       ! hbcut
      read(nf,afmt) (isymbl(i), i=1,natom)  ! isymbl
      ! now you have all you need; deallocate dummy arrays
      deallocate(i6d)
      deallocate(ed)
      deallocate(id,jd,kd,ld,ic)

   else                   ! this is new-style prmtop file
      fmtin = afmt
      type = 'RESIDUE_LABEL'
      call nxtsec(nf,  6,  0,fmtin,  type,  fmt,  iok)
      read(nf,fmt) (lbres(i),i=1,nres)
   
      fmtin = ifmt
      type = 'RESIDUE_POINTER'
      call nxtsec(nf,  6,  0,fmtin,  type,  fmt,  iok)
      read(nf,fmt) (ipres(i),i=1,nres)
      ipres(nres+1) = natom+1
   
      fmtin = afmt
      type = 'AMBER_ATOM_TYPE'
      call nxtsec(nf,  6,  0,fmtin,  type,  fmt,  iok)
      read(nf,fmt) (isymbl(i),i=1,natom)
   end if
   
end subroutine getnam


!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ get coordinates from amber coordinate file
subroutine getcor(natom,crd,nf)
   implicit none
   
   _REAL_ crd(*)
   integer natom, nf, i, nat3, ios, matom 
   character(len=4) ititl(20)
   character(len=80) line

   nat3 = 3*natom
   read(nf,'(20A4)') ititl
   read(nf,'(A)') line
   read(line(1:6),*,iostat=ios) matom
   if (ios > 0) stop "error reading natom from inpcrd"

   if(matom.ne.natom) then
      write(6,'("atom no. in prmtop and inpcrd do not match",2I8)') &
            natom,matom
      stop
   end if

   read(nf,'(6F12.7)') (crd(i),i=1,nat3)

end subroutine getcor


!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ generate selected residue/atom list in condensed form
subroutine genatlist(nres,igraph,ipres,nf,mask)
! Generate terse output (it is condensed as much as possible by trying
! to merge contiguous selections). The ':' designates residues and 
! '@' designates atoms.
! Note: string operations which are performed in this routine are
! quite a _nightmare_ in fortran. The logic of the routine is simple
! but it might be quite difficult to follow because of the use of
! multiple logical flags that keep track of the contiguity of residues.

   implicit none
   
   integer nres, nf, ipres(*)
   character(len=4) igraph(*)
   integer mask(*)

   character(80) buffer, flushline
   integer j, j1, j2, k, m
   logical ressel    ! .true. if at least one atom in j-th residue selected
   logical wressel   ! .true. if all atoms in j-th residues selected
   logical startset  ! .true. if residue span starting residue is set

   
   write(nf,'("Condensed format of the atom selection follows:")')
   
   flushline = ''
   startset = .false.
   do j = 1,nres
      ! first residue scan: find out if residue is selected:
      ! a) completely (wressel->true), b) partially (ressel->true)
      ! or c) not at all (ressel->false)
      wressel = .true.
      ressel = .false.
      do k = ipres(j), ipres(j+1)-1
         if (mask(k) == 1) then
            ressel = .true.   ! at least one atom in j-th residue is selected
         else
            wressel = .false. ! all atoms in the residue are not selected
         end if
      end do
      
      ! if whole residue is selected it starts a residue span (which
      ! we define to be at least one residue long)
      if (wressel .and. .not.startset) then
         j1 = j
         startset = .true.  ! this indicates that a residue span is started
      end if
      
      if (ressel) then  ! j-th residue is at least partially selected
         if (wressel) then
            j2 = j   ! if j2 is the last residue, write it here 
                     ! (because it doesnt get written elsewhere)
            if (j2 >= nres) then
               buffer = ''
               if ( j1 == j2 ) then
                  write(buffer,'(":",I8)') j1
               else
                  write(buffer,'(":",I8,"-",I8)') j1,j2
               end if
               call skipSpaceAndFlush(buffer,flushline,nf)
            end if
         else  ! j-th residue only partially selected
            if (startset) then  ! print previous residue span if it existed
               buffer = ''
               if ( j1 == j2 ) then
                  write(buffer,'(":",I8)') j1
               else
                  write(buffer,'(":",I8,"-",I8)') j1,j2
               end if
               call skipSpaceAndFlush(buffer,flushline,nf)
               startset = .false.
            end if
            ! print atom names selected in j-th residue
            buffer = ''
            write(buffer,'(":",I8,"@")') j
            m = index(buffer,'@') + 1
            do k = ipres(j), ipres(j+1)-1
               if (mask(k) == 1) then
                  write(buffer(m:),'(A4,",")') igraph(k)
                  m = m + 5  ! 4 characters + 1 comma
               end if
            end do
            call skipSpaceAndFlush(buffer,flushline,nf)
         end if
      else  ! j-th residue is not selected but it may have ended residue span
         if (startset) then  ! if it did, write the residue span out
            buffer = ''
            if ( j1 == j2 ) then
               write(buffer,'(":",I8)') j1
            else
               write(buffer,'(":",I8,"-",I8)') j1,j2
            end if
            call skipSpaceAndFlush(buffer,flushline,nf)
            startset = .false.
         end if
      end if

   end do ! end the loop over residues (do j=1,nres)

   ! flush whatever is left in the flushline
   write(nf,'(A)') flushline
   
end subroutine genatlist

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ create a buffer for condesed output and flush it when line fills
subroutine skipSpaceAndFlush(buffer,flushline,nf)
   implicit none
   character(*) buffer, flushline
   character(80) condensed
   integer k,i,nf, currlen
   
   i = 1
   condensed = ''
   do k = 1, len(trim(buffer))
      if (buffer(k:k).eq.' ') cycle
      condensed(i:i) = buffer(k:k)
      i = i + 1
   end do

   ! if this was list of atoms, get rid of last comma
   if (condensed(i-1:i-1).eq.',') condensed(i-1:i-1) = ' '

   currlen = len_trim(flushline)
   if (currlen > 50) then   
      write(nf,'(A)') flushline
      flushline(2:) = condensed
   else
      flushline(currlen+2:) = condensed
   end if
   
end subroutine skipSpaceAndFlush 


!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ generate pdb file (amber-style) from selected atoms
subroutine genpdb(nres,igraph,isymbl,ipres,lbres,nf,crd,title,mask)

   implicit none
  
   integer nres, nf, ipres(*)
   character(len=4) igraph(*), isymbl(*), lbres(*)   
   character(*) title
   integer mask(*)
   _REAL_ crd(*)

   integer j, k, jc, m

   ! generate pdb output: put amber atom type into pdb 'element symbol' (77-78)
   write(nf,'("REMARK  ",A)') title(1:72)
   ! write(nf,'(A)') "REMARK   10        20        30        40" //  &
   !                 "        50        60        70        8"
   do j = 1,nres
      do k = ipres(j), ipres(j+1)-1
         ! if atom is selected, write it out
         if (mask(k) == 1) then
            jc = 3*k-3
            write(nf,10) k,igraph(k),lbres(j),j,(crd(jc+m),m=1,3),isymbl(k)
10          format('ATOM',2X,I5,1X,A4,1X,A4,1X,I4,4X,3F8.3,22X,A2)
         end if
      end do
   end do
   write(nf,'("END   ")')
end subroutine genpdb


!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ generate selected residue/atom list in amber 'group input' form
subroutine gengroupinput(nres,ipres,mask)
! This routine is very similar to routine genatlist() above. However,
! it generates selected atoms in a form of ATOM/RES cards that can
! be directly inserted into amber input.
! TODO: detect ATOM number spans in individual residues

   implicit none
   
   integer nres, ipres(*)
   integer mask(*)

   integer j, j1, j2, k
   logical ressel    ! .true. if at least one atom in j-th residue selected
   logical wressel   ! .true. if all atoms in j-th residues selected
   logical startset  ! .true. if residue span starting residue is set

   
   ! write(nf,'("amber group input format of the atom selection:")')
   
   startset = .false.
   do j = 1,nres
      ! first scan through residues: find out if residue is selected:
      ! a) completely (wressel->true), b) partially (ressel->true)
      ! or c) not at all (ressel->false)
      wressel = .true.
      ressel = .false.
      do k = ipres(j), ipres(j+1)-1
         if (mask(k) == 1) then
            ressel = .true.   ! at least one atom in j-th residue is selected
         else
            wressel = .false. ! all atoms in the residue are not selected
         end if
      end do
      
      ! if whole residue is selected it starts a residue span (which
      ! we define to be at least one residue long)
      if (wressel .and. .not.startset) then
         j1 = j
         startset = .true.  ! this indicates that a residue span is started
      end if
      
      if (ressel) then  ! j-th residue is at least partially selected
         if (wressel) then
            j2 = j   ! if j2 is the last residue, write it here 
                     ! (because it doesnt get written elsewhere)
            if (j2 >= nres) write(*,'(" RES ",2I8)') j1,j2
         else  ! j-th residue only partially selected
            if (startset) then  ! print previous residue span if it existed
               write(*,'(" RES ",2I8)') j1,j2
               startset = .false.
            end if
            ! print atom names selected in j-th residue
            ! here is the only difference compared to routine genatlist():
            ! you also need to detect atom span - for now we can just ignore
            ! it (it may produce more ATOM cards but the end result is ok)
            do k = ipres(j), ipres(j+1)-1
               if (mask(k) == 1) then
                  write(*,'(" ATOM ",2I8)') k, k
               end if
            end do
         end if
      else  ! j-th residue is not selected but it may have ended residue span
         if (startset) then  ! if it did, write the residue span out
            write(*,'(" RES ",2I8)') j1,j2
            startset = .false.
         end if
      end if

   end do ! end the loop over residues (do j=1,nres)

end subroutine gengroupinput

