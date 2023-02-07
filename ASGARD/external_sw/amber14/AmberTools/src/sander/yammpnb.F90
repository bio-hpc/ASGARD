! <compile=optimized>
#include "copyright.h"
#include "../include/dprec.fh"

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ yammp-style nonbonded potential
subroutine yammpnb(x,f,iac,ico,numex,natex, &
      cut,cn1,cn2,ntypes,natom,evdw)
   
   
   !     Use the yammp-style nonbonded potential: like a lower-bound
   !     distance parabola.
   
   implicit none
   _REAL_  x(*)
   _REAL_  f(*)
   integer iac
   integer ico
   integer numex
   integer natex
   _REAL_  cut
   _REAL_  cn1
   _REAL_  cn2
   integer ntypes
   integer natom
   _REAL_  evdw

#ifdef MPI
#  include "parallel.h"
#ifdef MPI_DOUBLE_PRECISION
#undef MPI_DOUBLE_PRECISION
#endif
   include 'mpif.h'
#ifdef CRAY_PVP
#define MPI_DOUBLE_PRECISION MPI_REAL8
#endif
#endif
#  include "../include/md.h"
   
   _REAL_  de
   _REAL_  diff
   _REAL_  dis
   _REAL_  dedx, dedy, dedz
   _REAL_  dumx, dumy, dumz
   integer i, j
   integer ic
   integer iacI
   integer iexcl, jexcl, nexcl
   integer jexcl_last
   _REAL_  r2
   logical skip
   _REAL_  xi, yi, zi
   _REAL_  xij, yij, zij
   
   evdw = 0.0d0
   
   !    -- big double loop over all pairs of atoms:
   
#ifdef MPI
   iexcl = 1
   do i=1,mytaskid
      iexcl = iexcl + numex(i)
   end do
   do i=mytaskid+1,natom,numtasks
#else
   iexcl = 1
   do i=1,natom
#endif
      xi = x(3*i-2)
      yi = x(3*i-1)
      zi = x(3*i  )
      iaci = ntypes * (iac(i) - 1)
      jexcl = iexcl
      jexcl_last = iexcl + numex(i) -1
      nexcl = natex( jexcl )
      if( jexcl > jexcl_last ) nexcl = 0
      dumx = 0.0d0
      dumy = 0.0d0
      dumz = 0.0d0
      
      do j=i+1,natom
         
         !         -- check the exclusion list:
         
         skip = .false.
         if( j == nexcl ) then
            skip = .true.
            jexcl = jexcl + 1
            nexcl = natex( jexcl )
            if( jexcl > jexcl_last ) nexcl = 0
         end if
         if( .not. skip ) then
            
            de = 0.0d0
            xij = xi - x(3*j-2)
            yij = yi - x(3*j-1)
            zij = zi - x(3*j  )
            r2 = xij*xij + yij*yij + zij*zij
            dis = sqrt(r2)
            
            ic = ico( iaci + iac(j) )
            if( dis < cn2(ic) ) then
               diff = dis - cn2(ic)
               evdw = evdw + cn1(ic)*diff*diff
               de = -2.d0*cn1(ic)*diff/dis
               !             write(6,*) i,j,dis, cn2(ic),evdw
               dedx = de * xij
               dedy = de * yij
               dedz = de * zij
               dumx = dumx + dedx
               dumy = dumy + dedy
               dumz = dumz + dedz
               f(3*j-2) = f(3*j-2) - dedx
               f(3*j-1) = f(3*j-1) - dedy
               f(3*j  ) = f(3*j  ) - dedz
            end if
         end if
      end do  !  j=i+1,natom
      
      f(3*i-2) = f(3*i-2) + dumx
      f(3*i-1) = f(3*i-1) + dumy
      f(3*i  ) = f(3*i  ) + dumz
#ifdef MPI
      do k=1,numtasks
         iexcl = iexcl + numex(i+k-1)
      end do
#else
      iexcl = iexcl + numex(i)
#endif
   end do  !  i=mytaskid+1,natom,numtasks
   
   return
end subroutine yammpnb 
