#include "../include/dprec.fh"

!The 3D-RISM-KH software found here is copyright (c) 2012 by
!Tyler Luchko and David A. Case.
!
!This program is free software: you can redistribute it and/or modify it
!under the terms of the GNU General Public License as published by the Free
!Software Foundation, either version 3 of the License, or (at your option)
!any later version.
!
!This program is distributed in the hope that it will be useful, but
!WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
!or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
!for more details.
!
!You should have received a copy of the GNU General Public License in the
!../../LICENSE file.  If not, see <http://www.gnu.org/licenses/>.
!
!Users of the 3D-RISM capability found here are requested to acknowledge
!use of the software in reports and publications.  Such acknowledgement
!should include the following citations:
!
!1) A. Kovalenko and F. Hirata. J. Chem. Phys., 110:10095-10112  (1999); 
!ibid. 112:10391-10417 (2000).   
!
!2) A. Kovalenko,  in:  Molecular  Theory  of  Solvation,  edited  by  
!F. Hirata  (Kluwer Academic Publishers, Dordrecht, 2003), pp.169-275.  
!
!3) T. Luchko, S. Gusarov, D.R. Roe, C. Simmerling, D.A. Case, J. Tuszynski,
!and  A. Kovalenko, J. Chem. Theory Comput., 6:607-624 (2010). 

!This program takes 3D solvent distribution functions from
!3D-RISM calculations and calculates orientational averages.  

module orave_m
  use safemem
  use getopts_c
  use rism_report_c
  use rism3d_grid_c
  implicit none

  integer, parameter :: CLEN=1024

  !grid  :: grid object
  type(rism3d_grid),save :: grid

  !guvfile : filename for the input guv files
  character(len=CLEN),pointer :: guvfile(:)=>NULL()
  
  !method  : integration method for averaging selected by user
  !outfile : name of the output file, if it exists
  !pdbfile : pdbfile that contains all sample points
  character(len=CLEN) :: method, outfile='', pdbfile='test.pdb'

  !npoint : number of points for integration method
  integer :: npoint

  !guv        :: (product(nr),nsite) solvent distribution function
  _REAL_,pointer :: guv(:,:)=>NULL()

  !dxOrigin :: origin from the DX file.  Tells us how to translate the solute and
  !            how to interpret volume specifications
  !origin :: origin for integration
  !vec    :: orientation vector for polar averaging
  !radius :: radial extent for averaging
  !zlen   :: z-axis extent for averaging.  Centered on origin.
  !dr     :: radial spacing
  !dz     :: z-axis spacing for polar averaging
  _REAL_ :: dxOrigin(3), origin(3),vec(3),radius,zlen,dr,dz

  !abscissa :: list of abscissa (npoints, coordinates and weight)
  !grid2d :: 2d (radial, z, site) grid
  _REAL_, pointer :: abscissa(:,:)=>NULL(), grid2d(:,:,:)=>NULL()

  !z_ave :: perform averaing along the z-axis for polar averaging
  !polar :: polar averaging instead of spherical
  logical :: z_ave, polar=.false., writePDB=.false.

  integer :: nr, nz, zrange(2)
  _REAL_ :: zstart
  contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!get command line options
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine getOptions()
    implicit none
    integer :: err,i, ntoken
    character(len=CLEN),pointer :: extra(:)=>NULL()
    character(len=10*CLEN), pointer :: tempC(:)=>NULL()
    character(len=10*CLEN) :: tempOpt
    logical :: help
    call getopts_add("vec",'',"v",1)
    call getopts_add("origin",'default',"o",1)
    call getopts_add("radius",0d0,"r",1,1)
    call getopts_add("z",-1d0,"z",max=1)
    call getopts_add("dr",0.1d0,max=1)
    call getopts_add("dz",0.1d0,max=1)
    call getopts_add("np",0,max=1,min=1)
    call getopts_add("2d",.false.)
    call getopts_add("guv",'',min=1)
    call getopts_add("method",'',max=1,min=1)
    call getopts_add("pdbout",'',max=1)
    call getopts_add("help",.false.,"h")
    ntoken = getopts_process()
!    call getopts_summary(0)
    call getopts_get("help",help)
    if(help .or. ntoken == 0)then
       call usage()
       stop
    elseif(ntoken <0)then
       call usage()
       call rism_report_error("command line read failed.")
    end if
    
    !for comma separated arguments, read in arguments as strings, then parse
    !arguments into the arrays
    call getopts_get("origin",1,tempopt)
    if(trim(tempopt) .eq. 'default')then
       origin = huge(1d0)
    else
       read(tempopt,*,iostat=err) origin
       if(err/=0) call rism_report_error("Invalid cooridinates for --origin: "&
            //trim(tempopt))
    end if

    if(getopts_narg("vec") == 1)then
       call getopts_get("vec",1,tempopt)
       read(tempopt,*,iostat=err) vec
       if(err/=0) call rism_report_error("Invalid cooridinates for --vec: "&
         //trim(tempopt))
       polar = .true.
       !make it a unit vector
       vec = vec/sqrt(sum(vec**2))
    else
       if(getopts_narg("dz") /=0)&
            call rism_report_error("'--dz' only valid when used with '--vec'.")
       if(getopts_narg("2d") /=0)&
            call rism_report_error("'--2d' only valid when used with '--vec'.")
    end if

    if(getopts_narg("pdbout") == 1)then
       call getopts_get("pdbout",1,pdbfile)
       writePDB=.true.
    end if
    
    call getopts_get("radius",1,radius)
    call getopts_get("z",1,zlen)
    call getopts_get("dr",1,dr)
    call getopts_get("dz",1,dz)
    call getopts_get("np",1,npoint)
    call getopts_get("2d",z_ave)
    z_ave = .not. z_ave
    call getopts_get("method",1,method)

    guvfile=> getopts_getAll("guv",guvfile,len(guvfile))
    extra=> getopts_unread()
    if(size(extra) == 1) outfile = extra(1)
    if( getopts_narg("guv") == 0) then
       call usage()
       call rism_report_error("'--guv' must be specified")
    else if(size(extra) >1)then
       call usage()
       call flush(0)
       !concatenate all extra options into temp
       tempOpt=""
       do i=2, size(extra)
          tempOpt(len_trim(tempOpt)+1:) =  " "//trim(extra(i))
       end do
       call rism_report_error("unknown options:"//trim(tempOpt))
    end if

    if(safemem_dealloc(extra)/=0) &
         call rism_report_error("GETOPTONS: failed to dealloc 'extra'")
    if(safemem_dealloc(tempC)/=0) &
         call rism_report_error("GETOPTONS: failed to dealloc 'tempC'")

    call getopts_cleanup()

  end subroutine getOptions

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!read in input 3D distribution files
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine readDistributionFiles
    use rism3d_opendx
    implicit none
    integer :: i
    _REAL_, pointer :: tmp1(:)
    _REAL_ :: delta(3)
    integer, parameter :: zeroR3(3)=0d0
    integer :: npos(3), nkpos(3), tmppos(3)

    !
    !User supplied distribution functions
    !

    !initialize objects
    call rism3d_grid_new(grid)

    !get grid size
    if(ubound(guvfile,1) >0)then
       call readDXHeader(guvfile(1),dxOrigin,delta,npos)
    end if

    !set grid size
    nkpos = npos
    nkpos(1) = nkpos(1)+2
    call rism3d_grid_resize(grid,delta,npos,nkpos,npos,nkpos,zeroR3,zeroR3)
    !Guv
    if(ubound(guvfile,1) > 0)then
       guv => safemem_realloc(guv,grid%nrTotal,ubound(guvfile,1),.false.)
       do i = 1, ubound(guv,2)
          !check file size
          call readDXHeader(guvfile(i),dxOrigin,delta,tmppos)
          if(sum(abs(tmppos- npos))/=0) call rism_report_error(trim(guvfile(i))//" is the wrong size")
          call readDX(guvfile(i),guv(:,i),grid%nr,dxOrigin,delta)
       end do
    end if
    if(origin(1) == huge(1d0)) origin = dxOrigin + grid%boxlen/2d0
  end subroutine readDistributionFiles


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!allocates memory and sets values for polar abscissa and weights
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine setAbscissa()
    use constants, only : pi
    implicit none
    _REAL_ :: dtheta
    integer :: itheta

    if(polar)then
       if(trim(method) .eq. "uniform")then
          abscissa => safemem_realloc(abscissa,npoint,2)
          dtheta = 2d0*pi/dble(npoint)
          do itheta = 1,npoint
             abscissa(itheta,1) = (itheta-1)*dtheta
             abscissa(itheta,2) = dtheta
          end do
       else
          call usage()
          call rism_report_error("Not a valid method for polar averaging: "//trim(method))
       end if
    elseif(.not. polar)then
       call rism_report_error("Spherical averaging not yet supported.")
    end if
  end subroutine setAbscissa
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!calculates orientational averages
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine average
    implicit none
    if(polar)then
       call averagePolar()
    else
    end if
  end subroutine average

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!calculates polar coordinate orientational averages
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine averagePolar()
    use rism_util, only : freeUnit
    use constants, only : pi
    use quaternion
    implicit none
    integer :: ir, iz, itheta, iguv
    integer :: unit,iostat
    _REAL_ :: z, point(3), quat(4), cross(3), zaxis(3)=(/0,0,1/),&
         angle

    if(.not.inGrid(origin)) call rism_report_error("Origin must be located in solvent grid.")

    if(writePDB)then
       unit = freeUnit()
       open(unit,file=pdbfile,status='replace',iostat=iostat)
       if(iostat/=0)&
            call rism_report_error("failed to open : test.pdb")
    end if

    cross(1) = zaxis(2)*vec(3) - zaxis(3)*vec(2)
    cross(2) = zaxis(3)*vec(1) - zaxis(1)*vec(3)
    cross(3) = zaxis(1)*vec(2) - zaxis(2)*vec(1)
    angle = acos(dot_product(zaxis, vec))
    if(angle == 0d0) cross = (/0,0,1/)
    quat = new_quat(angle,cross)

    nr = dble(radius)/dble(dr)+1
    if(zlen<0d0)then
       !do this quick and dirty for now.  Ideally, we want the longest
       !possible cylinder that will fit in the box with our given
       !radius and orientation.  Rather, we will take the longest line
       !that can fit in the box and test at each point if it fits in
       !the box.  NOTE: the negative value indicates that this has not
       !been set by the user
       zlen = -sqrt(sum(grid%boxlen**2))
       zstart = -sqrt(sum((dxorigin-origin)**2))
       nz = ceiling(abs(zlen)/dble(dz))+1
       zstart = zstart - mod(zstart,dz)
    else
       zstart = -zlen/2d0
       zstart = zstart - mod(zstart,dz)
       nz = ceiling(-2d0*zstart/dble(dz))+1
    end if
    grid2d=>safemem_realloc(grid2d,nr,nz,ubound(guv,2))

    zrange = 0
    if(writePDB)then
       call writepdbline(unit,0,"C","GLY","A",0,origin)
       call writepdbline(unit,0,"H","GLY","A",1,dxorigin)
       call writepdbline(unit,0,"H","GLY","A",1,dxorigin+(/grid%boxlen(1),0d0,0d0/))
       call writepdbline(unit,0,"H","GLY","A",1,dxorigin+(/0d0,grid%boxlen(2),0d0/))
       call writepdbline(unit,0,"H","GLY","A",1,dxorigin+(/0d0,0d0,grid%boxlen(3)/))
       call writepdbline(unit,0,"H","GLY","A",1,dxorigin+(/grid%boxlen(1),grid%boxlen(2),0d0/))
       call writepdbline(unit,0,"H","GLY","A",1,dxorigin+(/0d0,grid%boxlen(3),grid%boxlen(3)/))
       call writepdbline(unit,0,"H","GLY","A",1,dxorigin+(/grid%boxlen(1),0d0,grid%boxlen(3)/))
       call writepdbline(unit,0,"H","GLY","A",1,dxorigin+(/grid%boxlen(1),grid%boxlen(2),grid%boxlen(3)/))
    end if
    do iguv=1,ubound(guv,2)
       do iz=1,nz
          z = zstart+(iz-1)*dz
          grid2d(nr,iz,iguv) = 0d0
          do itheta=1,ubound(abscissa,1)
             point(1) = radius*cos(abscissa(itheta,1))
             point(2) = radius*sin(abscissa(itheta,1))
             point(3) = z
             call rotate_quat(point,quat)
             point = point+origin
             if(inGrid(point))then
                if(zrange(1) == 0)then
                   zrange = iz
                else
                   zrange(2) = iz
                end if
                grid2d(nr,iz,iguv) = grid2d(nr,iz,iguv) + interpVal(point,guv(:,iguv))*abscissa(itheta,2)
                if(writePDB)&
                     call writepdbline(unit,iz,"N","GLY","A",iz+1,point,interpVal(point,guv(:,iguv))*abscissa(itheta,2))
             else
                if(writePDB)then
                   call writepdbline(unit,iz,"O","GLY","A",iz+1,point)
                elseif(zlen>0d0)then
                   write(0,*) "z=", z, point
                   call rism_report_error("z-length too large for origin and radius.")
                end if
                grid2d(:,iz,iguv) = huge(1d0)
                if(.not.writePDB) exit                  
             end if
          end do
          if(grid2d(nr,iz,iguv) == huge(1d0) .and..not.writePDB) cycle
          if(grid2d(nr,iz,iguv) /= huge(1d0)) grid2d(:nr-1,iz,iguv) = 0d0
          do ir=2,nr-1
             do itheta=1,ubound(abscissa,1)
                point(1) = (ir-1)*dr*cos(abscissa(itheta,1))
                point(2) = (ir-1)*dr*sin(abscissa(itheta,1))
                point(3) = z
                call rotate_quat(point,quat)
                point = point+origin
                if(.not.writePDB)then
                   grid2d(ir,iz,iguv) = grid2d(ir,iz,iguv) + interpVal(point,guv(:,iguv))*abscissa(itheta,2)
                else
                   if(grid2d(1,iz,iguv) /= huge(1d0))then
                      grid2d(ir,iz,iguv) = grid2d(ir,iz,iguv) + interpVal(point,guv(:,iguv))*abscissa(itheta,2)
                      call writepdbline(unit,iz,"N","GLY","A",iz+1,point,&
                           interpVal(point,guv(:,iguv))*abscissa(itheta,2))
                   else
                      call writepdbline(unit,iz,"O","GLY","A",iz+1,point)
                   end if
                endif
             end do
          end do

          point(1) = 0
          point(2) = 0
          point(3) = z
          call rotate_quat(point,quat)
          point = point+origin
          if(.not.writePDB)then
             grid2d(1,iz,iguv) = interpVal(point,guv(:,iguv))
          else
             if(grid2d(1,iz,iguv) /= huge(1d0))then
                grid2d(1,iz,iguv) = interpVal(point,guv(:,iguv))
                call writepdbline(unit,iz,"N","GLY","A",iz+1,point,interpVal(point,guv(:,iguv)))
             else
                call writepdbline(unit,iz,"O","GLY","A",iz+1,point)
             end if
          endif
       end do
       !Normalize
       do iz = 1, nz
          do ir = 2, nr
             if(grid2d(ir,iz,iguv) /= huge(1d0))&
                  grid2d(ir,iz,iguv) = grid2d(ir,iz,iguv)/2d0/pi
          end do
       end do
    end do
    if(writePDB)&
         close(unit)

  end subroutine averagePolar

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Determines if a given point is inside the grid
!!!IN:
!!!   point :: xyz point
!!!OUT:
!!!    .true. if it is in the grid
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function inGrid(point)
    implicit none
    _REAL_, intent(in) :: point(3)
    logical :: inGrid
    _REAL_ :: maxExtent(3)
    integer :: id
    inGrid = .true.

    maxExtent = dxOrigin + grid%boxlen - grid%grdspc

    do id=1,3
       if(point(id) < dxOrigin(id) .or. point(id) > maxExtent(id))then
          inGrid = .false.
          return
       end if
    end do
    
  end function inGrid

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Interpolates the value of point from the grid
!!!IN:
!!!   point :: xyz point
!!!   dist  :: distribution to interpolate
!!!OUT:
!!!    interpolated value
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function interpVal(point,dist) result(val)
    implicit none
    _REAL_, intent(in) :: point(3), dist(:)
    _REAL_ :: val
    _REAL_ :: findex(3)
    !nb :: neighbour list
    integer :: nb(8,3), i,j, k
    val = 1d0
    !1) get fraction grid index
    findex = (point - dxOrigin)/grid%grdspc+1d0
    !2) get eight nearest neighbours
    nb(1,:) = floor(findex) !x000
!    if(nb(1,3) == 1) write(0,*) "NB", NB(1,3), point(3), dxorigin(3), findex(3)
    nb(8,:) = nb(1,:)+1 !x111
!    if(nb(8,3) == grid%nr(3)) write(0,*) "NB", NB(8,3), point(3), dxorigin(3), findex(3)
    !catches the case that the point is exactly on the upper boundary
    do i = 1,3
       if(nb(1,i) == grid%nr(i) .and. floor(findex(i)) == ceiling(findex(i)))then
          nb(1,i) = nb(1,i)-1
          nb(8,i) = nb(8,i)-1
       end if
    end do
    !check that both are inbounds
    if(any(nb(1,:) < 1))then
       call rism_report_error("(a,3(f8.3))", "Interpolate lower out-of-bounds error for: ",point)
    end if
    if(nb(8,1) > grid%nr(1) .or. nb(8,2) > grid%nr(2) .or. nb(8,3) > grid%nr(3))then
!       write(0,*) nb(8,:), grid%nr, point(3), (floor(findex(3)) == ceiling(findex(3)))
       call rism_report_error("(a,3(f8.3))", "Interpolate upper out-of-bounds error for: ",point)
    end if

    nb(2,:) = (/nb(1,1),nb(1,2),nb(8,3)/) !x001
    nb(3,:) = (/nb(1,1),nb(8,2),nb(1,3)/) !x010
    nb(4,:) = (/nb(1,1),nb(8,2),nb(8,3)/) !x011
    nb(5,:) = (/nb(8,1),nb(1,2),nb(1,3)/) !x100
    nb(6,:) = (/nb(8,1),nb(1,2),nb(8,3)/) !x101
    nb(7,:) = (/nb(8,1),nb(8,2),nb(1,3)/) !x110
    do i = 1,8
       do j = i+1,8
          if(nb(i,1) == nb(j,1) .and.&
               nb(i,2) == nb(j,2) .and.&
               nb(i,3) == nb(j,3))then
             write(0,*) "neighbours ", i, " and ", j, " are the same:"
             write(0,*) nb(i,:)
             write(0,*) nb(j,:)
             write(0,*) point
             write(0,*) findex
          end if
       end do
    end do
    !turn fractional index to fractional position in voxel
    findex = findex-nb(1,:)
    if(any(findex > 1d0) .or. any(findex < 0d0))then
       write(0,*) "bad findex "
       write(0,*) findex
       write(0,*) point
       
    end if
    !interpolate
    call blend_103(findex(1),findex(2),findex(3),&
         dist(nb(1,1) + (nb(1,2)-1)*grid%nr(1) + (nb(1,3)-1)*grid%nr(1)*grid%nr(2)),&
         dist(nb(2,1) + (nb(2,2)-1)*grid%nr(1) + (nb(2,3)-1)*grid%nr(1)*grid%nr(2)),&
         dist(nb(3,1) + (nb(3,2)-1)*grid%nr(1) + (nb(3,3)-1)*grid%nr(1)*grid%nr(2)),&
         dist(nb(4,1) + (nb(4,2)-1)*grid%nr(1) + (nb(4,3)-1)*grid%nr(1)*grid%nr(2)),&
         dist(nb(5,1) + (nb(5,2)-1)*grid%nr(1) + (nb(5,3)-1)*grid%nr(1)*grid%nr(2)),&
         dist(nb(6,1) + (nb(6,2)-1)*grid%nr(1) + (nb(6,3)-1)*grid%nr(1)*grid%nr(2)),&
         dist(nb(7,1) + (nb(7,2)-1)*grid%nr(1) + (nb(7,3)-1)*grid%nr(1)*grid%nr(2)),&
         dist(nb(8,1) + (nb(8,2)-1)*grid%nr(1) + (nb(8,3)-1)*grid%nr(1)*grid%nr(2)),&
         val)

  end function interpVal

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!write a line to a pdbfile
!!!IN:
!!!   unit : unit to write to
!!!   serial : atom number
!!!   atmname : atom name
!!!   resname : residue name
!!!   chain   : chain name
!!!   resid   : residue number
!!!   xyz     : coordinates
!!!   beta_o  : (optional) beta value
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine writePDBLine(unit, serial,atmname,resname,chain,resid,xyz,beta_o)
    implicit none
    integer, intent(in) :: unit,serial,resid
    character(len=*), intent(in) :: atmname, resname,chain
    _REAL_, intent(in) :: xyz(3)
    _REAL_, optional, intent(in) :: beta_o
    _REAL_ :: beta
    beta = 0
    if(present(beta_o)) beta = beta_o
    write(unit,'(A6,I5,1X,A4,A1,A3,1X,A1,I4,A1,3X,3F8.3,2F6.2,10X,2A2)')&
         "ATOM  ",serial,atmname,"",resname,chain,resid,'',xyz,beta
  end subroutine writePDBLine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!print results
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine writeOutput
    use rism_util, only : freeUnit
    implicit none
    integer :: unit, iostat
    if(len_trim(outfile) > 0)then
       unit = freeUnit()
       open(unit,file=outfile,status='replace',iostat=iostat)
    else
       unit = 6
    end if
    if(iostat/=0)&
         call rism_report_error("failed to open :"//trim(outfile))
    if(z_ave)then
       call write1DOutput(unit)
    else
       call write2DOutput(unit)
    end if
    if(len_trim(outfile) > 0)then
       close(unit)
    end if
  end subroutine writeOutput

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!print 1D radial results
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine write1DOutput(unit)
    use rism_util, only : rmExPrec
    implicit none
    integer, intent(in) :: unit
    character(len=80) :: fmt
    integer :: ir, iv
    write(fmt,'(a,i4,a)') '(1p,',(ubound(grid2d,3)),'(e16.7E3))'
    do ir=1,nr
       write(unit,fmt,advance='no') (ir-1)*dr
       do iv = 1, ubound(grid2d,3)
          write(unit,fmt,advance='no') &
               rmExPrec(sum(grid2d(ir,zrange(1):zrange(2),iv))&
               /(zrange(2)-zrange(1)+1))
       end do
       write(unit,'(a)')
    end do
  end subroutine write1DOutput

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!print 2D radial results
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine write2DOutput(unit)
    use rism_util, only : rmExPrec
    implicit none
    integer, intent(in) :: unit
    character(len=80) :: fmt
    integer :: ir, iz
    write(fmt,'(a,i4,a)') '(1p,',(2+ubound(grid2d,3)),'(e16.7E3))'
    do ir=1,nr
       do iz=1,nz
          if(grid2d(1,iz,1) /= huge(1d0))&
               write(unit,fmt) (ir-1)*dr, (iz-1)*dz + zstart, rmExPrec(grid2d(ir,iz,:))
       end do
       write(unit,'(a)')
    end do
  end subroutine write2DOutput

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Clean up memory
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine cleanup
    implicit none
    integer*8 :: memstats(10)
    integer :: unit , i
    unit = rism_report_getMUnit()
    memstats = memStatus()
    call rism3d_grid_destroy(grid)
    if(safemem_dealloc(guv)/=0) call rism_report_error("Failed to deallocate guv")
    if(safemem_dealloc(guvfile)/=0) call rism_report_error("Failed to deallocate guvfile")
    if(safemem_dealloc(abscissa)/=0) call rism_report_error("Failed to deallocate abscissa")
    if(safemem_dealloc(grid2d)/=0) call rism_report_error("Failed to deallocate grid2d")
    memstats = memStatus()
    write(unit,'(a)') "rism3d.orave memory allocation summary"
    write(unit,'(a)') "Type         Current         Maximum"
    write(unit,'(a,i12,a,f12.5,a)') "Integer  ",memstats(1)," B ",&
         dble(memstats(6))/BYTES_PER_GB," GB"
    write(unit,'(a,i12,a,f12.5,a)') "Real     ",memstats(2)," B ",&
         dble(memstats(7))/BYTES_PER_GB," GB"
    write(unit,'(a,i12,a,f12.5,a)') "Logical  ",memstats(3)," B ",&
         dble(memstats(8))/BYTES_PER_GB," GB"
    write(unit,'(a,i12,a,f12.5,a)') "Character",memstats(4)," B ",&
         dble(memstats(9))/BYTES_PER_GB," GB"
    write(unit,'(a)') "---------------------------------------"
    write(unit,'(a,i12,a,f12.5,a)') "Total    ",memstats(5)," B ",&
         dble(memstats(10))/BYTES_PER_GB," GB"
  end subroutine cleanup

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Usage description
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine usage
    implicit none
    integer :: munit
    munit = rism_report_getMUnit()
    call rism_report_setMUnit(rism_report_getEUnit())
    call rism_report_message("USAGE:   orave (-g|--guv) guvfiles --method name --np npoints")
    call rism_report_message("               (-r|--radius) radius")
    call rism_report_message("               [-o|--origin x,y,z] [--dr dr]  ")
    call rism_report_message("               [-v|--vec x,y,z [-z|--z z-length] [--dz dz] [--2d]]")
    call rism_report_message("               [--pdbout file] [output file]")
    call rism_report_message("")
    call rism_report_message("Orientationally averages the distributions given on the commandline.")
    call rism_report_message("The default averaging is spherical.  If a vector is supplied, then ")
    call rism_report_message("averaging is done with polar coordinates. If '--2d' is specified, then")
    call rism_report_message("averageing is not done over the z-axis and a 2D data set is produced.")
    call rism_report_message("Multiple methods for choosing orientational points to sample at are")
    call rism_report_message("available, though they may be exclusively for spherical or polar")
    call rism_report_message("averaging only. Output is to standout unless a filename is supplied.")
    call rism_report_message("")
    call rism_report_message("--guv     : list of guv files to average.  For spherical and polar")
    call rism_report_message("            averaging, there is one column of output data for each guv file.")
    call rism_report_message("--method  : Use abscissa and weights from method")
    call rism_report_message("              LEB (Lebedev) (spherical)")
    call rism_report_message("              REP (REPULSION) (spherical)")
    call rism_report_message("              ZCW (spherical)")
    call rism_report_message("              SHREWD_REP (SHREWD + REPULSION) (spherical)")
    call rism_report_message("              SHREWD_ZCW (SHREWD + ZCW) (spherical)")
    call rism_report_message("              uniform (spherical or polar)")
    call rism_report_message("            There are only a few sets of valid NP for each method")
    call rism_report_message("            except 'uniform'.")
    call rism_report_message("--np      : Number of sampling points for the averaging method.")
    call rism_report_message("--radius  : Maximum radius to average at.")
    call rism_report_message("--origin  : (center of grid) origin for averaging.")
    call rism_report_message("--dr      : (0.1) radial spacing")
    call rism_report_message("--vec     : orientation vector for polar averaging")
    call rism_report_message("--z       : Maximum length along z-axis.  Otherwise,")
    call rism_report_message("            largest cylinder with radius --radius is used.")
    call rism_report_message("--dz      : (0.1) z-axis spacing for polar averaging")
    call rism_report_message("--2d      : Do not average over z-axis for polar averaging")
    call rism_report_message("--pdbout  : Debugging. Considerably slows calculation. Outputs a")
    call rism_report_message("            PDB format file of all the points sampled. Points in")
    call rism_report_message("            the solvent box are atom type 'N'; those outside are")
    call rism_report_message("            of type 'O'.  Sampled values are placed in the occupancy")
    call rism_report_message("            column.")
    call rism_report_setMUnit(munit)
  end subroutine usage
end module orave_m

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Main program
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program orave
  use orave_m
  implicit none

  call getOptions()
  call readDistributionFiles()
  call setAbscissa()
  call average()
  call writeOutput()
  call cleanup()

end program orave

