#include "../include/dprec.fh"

module pupildata
   use memory_module, only: ixStack=>ix, realStack=>x, ihStack=>ih
   logical pupactive
   integer puperror,iPup,jPup,bs1,bs2,iresPup,l_puptmp
   integer pupLevelData                 ! Data to pass throught PUPIL Interface
   integer pupStep
   integer pupQZchange                  ! is != 0 if quantum zone has changed
   integer pupparticles                 ! total number of PUPIL particles
   integer pupqmatoms                   ! total number of PUPIL quantum atoms
   integer pupnumdomains                ! total number of PUPIL Quantum Domains
   integer pupnumnb14                   ! total number of initial nonboded 14 pairs
   integer pupnbonh                     ! total number of initial boded H  pairs
   integer pupnbona                     ! total number of initial boded    pairs
   integer pupntheth                    ! total number of initial angled H  tern.
   integer pupntheta                    ! total number of initial angled    tern.
   integer pupnphih                     ! total number of initial dihed. H  quat.
   integer pupnphia                     ! total number of initial dihed.    quat.
   integer, allocatable :: pupnb14  (:) ! initial nonbonded 14 pair list  3,x since it contains ii,jj,ic0
   integer, allocatable :: pupbonh  (:) ! initial boded H  pairs
   integer, allocatable :: pupbona  (:) ! initial boded    pairs
   integer, allocatable :: puptheth (:) ! initial angled H  tern.
   integer, allocatable :: puptheta (:) ! initial angled    tern.
   integer, allocatable :: pupphih  (:) ! initial dihed. H  quat.
   integer, allocatable :: pupphia  (:) ! initial dihed.    quat.
   integer, allocatable :: atmqmdomains (:) ! quantum domain to which belongs each  QM Particle
                                            !  AtmNum1,QMDom1,AtmNum2,QMDom2,... dim=2*pupnumdomains
                                            ! Domain >=1 -> QM Domain ID
   integer, allocatable :: pupqlist (:) ! atom number for quantum list
   integer, allocatable :: pupatm   (:)
   integer, allocatable :: pupres   (:) ! residue pointers
   !integer, allocatable :: ixStack (:) ! Captioning the ix pointer
   _REAL_,  allocatable :: pupchg   (:) ! Initial charges over each atoms in the system
   _REAL_,  allocatable :: qcdata   (:)
   _REAL_,  allocatable :: qfpup    (:)
   _REAL_,  allocatable :: qcell    (:)
   !_REAL_,  allocatable :: realStack(:) ! Captioning the x pointer
   !character(len=4), dimension(:), allocatable :: ihStack  ! Captioning the ih pointer
   character(len=10), dimension(:),allocatable :: keyres   ! to keep the MM key  RESIDUE_NAME
   character(len=20),dimension(:), allocatable :: keyMM    ! to keep the MM key particles RESIDUE_NAME.ATOM_NAME
   character(len=20)  strAux
   _REAL_   qmEnergy
end module pupildata

! ***********************************************************************
! ***********************************************************************

subroutine get_atomic_number_pupil(name,rmass,atomic_number)
!
!     Joan Torras (February 2014)

!     This subroutine assigns atomic number based upon the first and second
!     letter of the atom symbol and checking with the atomic mass 
!
!     name        :   characters containing the atomic name
!     rmass       :   atomic mass which has to fit with the given atom name
!     atomic_number:   atomic number which has to fit with the given atom name 
!
   implicit none

   character(len=4), intent(in)  :: name
   double precision, intent(in)  :: rmass
   integer,          intent(out) :: atomic_number
  
   logical :: check_wmass
   integer :: j, k, coincidence
   character(len=4) :: to_upper
   character(len=4) :: auxname
   character(len=2), dimension(109):: elemnt=(/                &
    'H ','HE',                                            &
    'LI','BE','B ','C ','N ','O ','F ','NE',              &
    'NA','MG','AL','SI','P ','S ','CL','AR',              &
    'K ','CA',                                            &
    'SC','TI','V ','CR','MN','FE','CO','NI','CU','ZN',    &
              'GA','GE','AS','SE','BR','KR',              &
    'RB','SR',                                            & 
    'Y ','ZR','NB','MO','TC','RU','RH','PD','AG','CD',    &
              'IN','SN','SB','TE','I ','XE',              &
    'CS','BA',                                            &
    'LA','CE','PR','ND','PM','SM','EU',                   &
    'GD','TB','DY','HO','ER','TM','YB',                   &
    'LU','HF','TA','W ','RE','OS','IR','PT','AU','HG',    &
              'TL','PB','BI','PO','AT','RN',              &
    'FR','RA',                                            &
    'AC','TH','PA','U ','NP','PU','AM',                   &
    'CM','BK','CF','ES','FM','MD','NO',                   &          
    'LR','RF','DB','SG','BH','HS','MT'/)
!
!      compare values in name to those in elemnt and assign atomic
!      numbers accordingly
!

   auxname = to_upper(name)
   atomic_number = 0
   j = 1
   coincidence = 0
   do while ( (j .le. 109) .and. (coincidence .eq. 0) )
      if(auxname(1:1) .eq. elemnt(j)(1:1)) then
         if(elemnt(j)(2:2) .eq. " ") then
            if(auxname(2:2) .eq. " ") then
               atomic_number = j
               coincidence = coincidence + 1

            else
               do k=j+1,109
                  if((auxname(1:1) .eq. elemnt(k)(1:1)) .and.   &
                     (auxname(2:2) .eq. elemnt(k)(2:2)))then
                     atomic_number = k
                     coincidence = coincidence + 1

                  endif
               enddo
               if (coincidence .eq. 0) then
                     atomic_number = j
                     coincidence = coincidence + 1

               endif
            endif
         else
            if(auxname(2:2) .eq. elemnt(j)(2:2)) then
               atomic_number = j
               coincidence = coincidence + 1

            endif
         endif
      end if
      j = j + 1
   enddo

   ! Checking coincidence between atomic_number and atomic mass
   if (atomic_number .gt. 0 ) then
      if (.not. check_wmass(atomic_number,rmass)) then
         coincidence = 0
      endif
   endif

   ! If not found so far lets shearch by atomic mass
   if(coincidence .eq. 0) then
      j = 1
      do while ( (j .le. 109) .and. (coincidence .eq. 0) )
         if (check_wmass(j,rmass)) then
            atomic_number = j
            coincidence = 1
         endif
         j = j + 1
      enddo
   endif

   if(coincidence .eq. 0) then
      write(6,*) 'PUPIL: Unable to correctly identify element ', name
      call mexit(6,1)
   else
      !write(6,*) 'FOUND: ',elemnt(atomic_number),' related to ',name
   endif

   return

end subroutine get_atomic_number_pupil

! ***********************************************************************
! ***********************************************************************

function to_upper(strIn) result(strOut)
! Adapted from http://www.star.le.ac.uk/~cgp/fortran.html (25 May 2012)
! Original author: Clive Page

     implicit none

     character(len=4), intent(in) :: strIn
     character(len=4) :: strOut
     integer :: i,j

     do i = 1, len(strIn)
          j = iachar(strIn(i:i))
          if (j>= iachar("a") .and. j<=iachar("z") ) then
               strOut(i:i) = achar(iachar(strIn(i:i))-32)
          else
               strOut(i:i) = strIn(i:i)
          end if
     end do

end function to_upper

! ***********************************************************************
! ***********************************************************************

logical function check_wmass (atomic_number, rmass)
!
!       Joan Torras (February 2014)
!
!       This subroutine is checking if the pair made of atomic_number and rmass have
!       coincidence with the values from periodic table
!
!       return a logical value ==> Are the pair of given values:  mass, and atomic 
!       number, in coincidence with the periodic table values with a maximum difference 
!       of 0.5 a.u. of mass ?

   implicit none

   integer,          intent(in) :: atomic_number
   double precision, intent(in) :: rmass

   integer :: j, k, coincidence

   double precision, dimension(109):: wmass=(/  &
  1.0079,  4.0026,6.941   ,9.0122  ,10.811  , 12.0107,14.0067 , 15.9994, 18.9984, 20.1797, &   
 22.9897, 24.305 ,26.9815 ,28.0855 ,30.9738 , 32.065 ,35.453  , 39.948 , 39.0983, 40.078 , &  
 44.9559, 47.867 ,50.9415 ,51.9961 ,54.938  , 55.845 ,58.9332 , 58.6934, 63.546 , 65.39  , & 
 69.723 , 72.64  ,74.9216 ,78.96   ,79.904  , 83.8   ,85.4678 , 87.62  , 88.9059, 91.224 , &
 92.9064, 95.94  ,98.0    ,101.07  ,102.9055,106.42  ,107.8682,112.411 ,114.818 ,118.71  , &
121.76  ,127.6   ,126.9045,131.293 ,132.9055,137.327 ,138.9055,140.116 ,140.9077,144.24  , &
145.0   ,150.36  ,151.964 ,157.25  ,158.9253,162.5   ,164.9303,167.259 ,168.9342,173.04  , &
174.967 ,178.49  ,180.9479,183.84  ,186.207 ,190.23  ,192.217 ,195.078 ,196.9665,200.59  , &
204.3833,207.2   ,208.9804,209.0   ,210.0   ,222.0   ,223.0   ,226.0   ,227.0   ,232.0381, &
231.0359,238.0289,237.0   ,244.0   ,243.0   ,247.0   ,247.0   ,251.0   ,252.0   ,257.0   , &
258.0   ,259.0   ,262.0   ,261.0   ,262.0   ,266.0   ,264.0   ,277.0   ,268.0/)

   check_wmass = abs(rmass-wmass(atomic_number)) .le. 0.5
   !write (6,*) check_wmass, atomic_number, rmass, wmass(atomic_number) 

end function check_wmass

! ***********************************************************************
! ***********************************************************************

subroutine deleting_qm_atoms()
  use pupildata, x=>realStack , ix=>ixStack, ih=>ihStack
  use parms, only: req
  implicit none
   
#include "../include/memory.h"
#include "extra_pts.h"
#include "../include/md.h"

  ! integer i

!  Initializing the list of structures for a new QM zone
  call init_extra_pts( &
         ix(iibh),ix(ijbh),ix(iicbh), &
         ix(iiba),ix(ijba),ix(iicba), &
         ix(i24),ix(i26),ix(i28),ix(i30), &
         ix(i32),ix(i34),ix(i36),ix(i38), &
         ix(i40),ix(i42),ix(i44),ix(i46),ix(i48), &
         ix(i50),ix(i52),ix(i54),ix(i56),ix(i58), &
         ih(m06),ix,x,ix(i08),ix(i10), &
         nspm,ix(i70),x(l75),tmass,tmassinv,x(lmass),x(lwinv),req)

  if(nbonh.gt.0) call setbon(nbonh,ix(iibh),ix(ijbh),ix(iicbh), &
    ix(ibellygp)) ! remove bonds between QM atoms from list

  if(nbona.gt.0) call setbon(nbona,ix(iiba),ix(ijba),ix(iicba), &
    ix(ibellygp)) ! remove bonds between QM atoms from list

  if(ntheth.gt.0) call setang(ntheth,ix(i24),ix(i26),ix(i28),ix(i30), &
    ix(ibellygp)) ! remove angles between QM atoms from list

  if(ntheta.gt.0) call setang(ntheta,ix(i32),ix(i34),ix(i36),ix(i38),&
    ix(ibellygp)) ! remove angles between QM atoms from list

  if(nphih.gt.0) call setdih(nphih,ix(i40),ix(i42),ix(i44),ix(i46), &
    ix(i48), ix(ibellygp)) ! remove dihedrals between QM atoms from list

  if(nphia.gt.0) call setdih(nphia,ix(i50),ix(i52),ix(i54),ix(i56), &
    ix(i58), ix(ibellygp)) ! remove dihedrals between QM atoms from list

  !if(pupQZChange .ne. 0) then
  !  write(6,*) 'NUMBER OF H BONDS            = ',nbonh
  !  write(6,*) 'NUMBER OF   BONDS            = ',nbona
  !  write(6,*) 'NUMBER OF H ANGLES           = ',ntheth
  !  write(6,*) 'NUMBER OF   ANGLES           = ',ntheta
  !  write(6,*) 'NUMBER OF H DIHEDRALS        = ',nphih
  !  write(6,*) 'NUMBER OF   DIHEDRALS        = ',nphia
  !  write(6,*) 'NUMBER OF NONBONDED 14 PAIRS = ',numnb14
  !  do i=1,numnb14
  !     write(6,*) ix(inb_14+(i-1)*3),ix(inb_14+(i-1)*3+1)
  !  enddo
  !endif
  
!  write(6,*) 'nbonh',nbonh,'iibh',iibh,'nbona',nbona,'ntheth',ntheth
!  write(6,*) 'ntheta',ntheta,'nphih',nphih,'nphia',nphia
!  write(6,*) 'i26',i26,'i56',i56
!  write(6,*) 'IGB',igb

  return
end subroutine deleting_qm_atoms

! ***********************************************************************
! ***********************************************************************

subroutine write_energies(str,e)

  implicit none

  _REAL_ e(25)
  character(len=4)  str
  integer i,k
  do i=1,5
    write(6,"(a4,2x,'ENERG',5(2x,i2,2x,d15.8))")   &
            str,( (i-1)*5+k,e((i-1)*5+k), k=1,5 )
  enddo
  write(6,*) " ----------------------------------"
end subroutine write_energies

! ***********************************************************************
! ***********************************************************************

subroutine add_vdwqmqm ( q, f, ener, ntypes, atom_name,     &
                               atom_type, atom_type_index  ) 

!
!       Joan Torras (June 2014)
!
!       This subroutine add vdw forces and energy between those qm particles
!       that belongs to different QM Domains
!

   use parms,         only: cn1, cn2, nttyp
   use pupildata,     only: pupqmatoms,atmqmdomains
   use memory_module, only: natom
   use nblist,        only: cutoffnb, ucell
   use state
 
   implicit none


   _REAL_ , intent(in   ) :: q(natom)
   _REAL_ , intent(inout) :: f(natom)
   type(state_rec),  intent(inout) :: ener
   integer, intent(in   ) :: atom_type_index(natom)
   integer, intent(in   ) :: ntypes
   character(len=4), intent(in) :: atom_name(natom)
   character(len=4), intent(in) :: atom_type(natom)
   
   integer :: m, n, i, j, idx, jdx, ic, idom, jdom, itype, jtype
   integer :: a, b, k, l, p 
   _REAL_  :: dx, dy, dz, fx, fy, fz, rij2, rij2_inv, r6, f6, f12, ff, rcut2
   _REAL_  :: dxi, dyi, dzi, okx, oky, okz, olx, oly, olz, opx, opy, opz
   _REAL_  :: evdw,vdwpot,dvdwpot,sw,dsw

!  +---------------------------------------------------------------------------+
!  |  V = A_ij / R_ij^12 - B_ij / R_ij^6                                       |
!  |  F_Ri = (12 * A_ij / R_ij^14 - 6 * B_ij / R_ij^8 )                        |
!  |       * ( R_i - R_j )                                                     |
!  +---------------------------------------------------------------------------+

   evdw    = 0d0
   rcut2   = cutoffnb*cutoffnb

   !write(6,*) 'CUTOFFNB', cutoffnb

   do m = 1,   pupqmatoms-1
   
      i    = atmqmdomains(m*2 - 1)
      idom = atmqmdomains(m*2    )
      itype= atom_type_index(i)
      
      !write(6,*) 'ATOM NUM',i,atom_name(i),atom_type(i),itype, idom
      
      do n = m+1, pupqmatoms

         j    = atmqmdomains(n*2 - 1)
         jdom = atmqmdomains(n*2    )
         jtype= atom_type_index(j)
 
         !write(6,*) '      ATOM NUM',j,atom_name(j),atom_type(j),jtype, jdom
 
         ! Checking for different QM domain particles
         if (idom .ne. jdom) then
 
            a = min(itype,jtype)
            b = max(itype,jtype)
            ! In case of semi-matrix filling by rows 
            !ic = int((a - 1) * (2*ntypes + 2 - a)/2 + (b - a + 1))
            
            ! In case of semi-matrix filling by columns
            ic = int(b*(b - 1) / 2 + a)

            !write (6,*) '  types',a,b,'pair order ',ic,'total types',ntypes,'A',cn1(ic),'B',cn2(ic)
    
            idx = ( i - 1 ) * 3 
            jdx = ( j - 1 ) * 3

            dxi = q(idx+1) - q(jdx+1) 
            dyi = q(idx+2) - q(jdx+2) 
            dzi = q(idx+3) - q(jdx+3) 

            ! Applying minimum image criterion 
            ! for an Orthorhombic unit cell.
            ! Pupil is dealing with this kind of unit cell so far
            dx = dxi - ucell(1,1) * anint( dxi / ucell(1,1) )
            dy = dyi - ucell(2,2) * anint( dyi / ucell(2,2) )
            dz = dzi - ucell(3,3) * anint( dzi / ucell(3,3) )

            rij2 = dx*dx + dy*dy + dz*dz

            if (rij2 .le. rcut2) then
               call get_switchF(sw,dsw,rij2,rcut2)

               rij2_inv = 1.0d0 / rij2
               r6 = rij2_inv * rij2_inv * rij2_inv

               f6  = cn2(ic) * r6
               f12 = cn1(ic) * r6 * r6 

               vdwpot = (f12 - f6)
               evdw = evdw + sw * vdwpot

               dvdwpot = ( 12.0d0 * f12 - 6.0d0 * f6 ) * rij2_inv 

               ff = dsw * vdwpot + dvdwpot * sw

               fx = ff * dx
               fy = ff * dy
               fz = ff * dz

               f(idx+1) = f(idx+1) + fx
               f(idx+2) = f(idx+2) + fy
               f(idx+3) = f(idx+3) + fz

               f(jdx+1) = f(jdx+1) - fx
               f(jdx+2) = f(jdx+2) - fy
               f(jdx+3) = f(jdx+3) - fz

               !write(6,*) 'CONSIDERING ATOM PAIR',i,atom_name(i),idom,j,atom_name(j), jdom, 'DISTANCE',sqrt(rij2)
               !write(6,*) 'DISTANCE',sqrt(rij2),'FORCE',fx,fy,fz,'PARMS',cn1(ic),cn2(ic)

            endif
         endif
      enddo
   enddo

!  Include the above VDW interactions ... between QM atoms from different
!  quantum Domains (those from the same QM Domain are neglected)         

   !write(6,*) 'VDW ENERGY(previous):',ener%pot%vdw, 'VDW qmqm:',evdw

   ener%pot%vdw = ener%pot%vdw + evdw

end subroutine add_vdwqmqm


!***********************************************************************
subroutine get_switchF ( sw, dsw, rij2, rcut2)

!
!       Joan Torras (September 2014)
!
!       This subroutine calculate the switch function of vdw potential
 
   implicit none


   _REAL_ , intent(inout) :: sw
   _REAL_ , intent(inout) :: dsw
   _REAL_ , intent(in   ) :: rij2
   _REAL_ , intent(in   ) :: rcut2

   _REAL_ :: ron,ron2,den


   ron  = sqrt(rcut2) - 3.0
   ron2 = ron*ron
   den  = (rcut2 - ron2)
   den  = 1.0 / (den*den*den)
   
   if (rij2 .le. ron2) then
       sw  = 1.0
       dsw = 0.0
   else if ( (rij2 .gt. ron2) .and. (rij2 .le. rcut2) ) then
       sw  = ( (rcut2 - rij2)*(rcut2 - rij2)*(rcut2 + 2 * rij2 - 3 * ron2) ) * den
       dsw = 12 * (rcut2 - rij2) * ((ron2/3) - rij2) * den
   else
       sw  = 0.0
       dsw = 0.0
   endif

   !write(6,*) 'SWITCH FUNCT.', sw, 'DERIVATIVE', dsw, 'RIJ', sqrt(rij2)
   
end subroutine get_switchF   
