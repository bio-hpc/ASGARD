      Subroutine se_dipole
      implicit none
#include "sebomd.dim"
#include "se_logic.h"
#include "se_elemnt.h"
#include "se_fulsys.h"
#include "se_xyzcrd.h"
#include "se_glbmat1.h"
#include "se_parproc.h"
#include "se_param.h"
#include "se_param1.h"
#include "se_param2.h"
#include "se_param3.h"
#include "se_param4.h"
#include "se_param5.h"
#include "se_param6.h"
#include "se_param7.h"

!  Dipole.F calculates the dipole of a molecule.  

!  Variable declaration.

      integer counter
      integer i,j
      double precision ptchg(3)
      double precision hybfac(maxatm), hybdip(3)
      double precision x, const, frac1, frac2
      double precision dipol1

!  Initialization block

      do i=1,3
         ptchg(i)=0.d0
         hybdip(i)=0.d0
      enddo
      do j=1,natoms
         hybfac(j)=0.d0
         dipole_atm(1,j) = 0.0d0
         dipole_atm(2,j) = 0.0d0
         dipole_atm(3,j) = 0.0d0
      enddo

      do j=1,natoms
         do i=1,3
            ptchg(i)=ptchg(i)+4.803d0*atchg(j)*xyz(i,j)
            dipole_atm(i,j) = 4.803d0*atchg(j)*xyz(i,j)
         enddo
      enddo

!  Calculate hybridization factors and hybridization contribution to
!  dipole for the mUlliken charge system.

      do j=1,natoms
         x=dble(nquant(iatnum(j)))
         const=((5.0832d0)*(2.d0*x+1.d0)*(2.d0**(2.d0*x+1.d0)))
     .    /sqrt(3.d0)

         if(mndo) then
            frac1=(expnt1(0,iatnum(j))*expnt1(1,iatnum(j)))**
     .       (x+.5d0)
            frac2=(expnt1(0,iatnum(j))+expnt1(1,iatnum(j)))**
     .       (2.d0*x+2.d0)
         elseif (am1) then
            frac1=(expnt2(0,iatnum(j))*expnt2(1,iatnum(j)))**
     .       (x+.5d0)
            frac2=(expnt2(0,iatnum(j))+expnt2(1,iatnum(j)))**
     .       (2.d0*x+2.d0)
         elseif (rm1) then
            frac1=(expnt5(0,iatnum(j))*expnt5(1,iatnum(j)))**
     .       (x+.5d0)
            frac2=(expnt5(0,iatnum(j))+expnt5(1,iatnum(j)))**
     .       (2.d0*x+2.d0)
         elseif (am1d) then
            frac1=(expnt6(0,iatnum(j))*expnt6(1,iatnum(j)))**
     .       (x+.5d0)
            frac2=(expnt6(0,iatnum(j))+expnt6(1,iatnum(j)))**
     .       (2.d0*x+2.d0)
         elseif (pm3pddg) then
            frac1=(expnt7(0,iatnum(j))*expnt7(1,iatnum(j)))**
     .       (x+.5d0)
            frac2=(expnt7(0,iatnum(j))+expnt7(1,iatnum(j)))**
     .       (2.d0*x+2.d0)
         elseif (pm3) then
            frac1=(expnt3(0,iatnum(j))*expnt3(1,iatnum(j)))**
     .       (x+.5d0)
            frac2=(expnt3(0,iatnum(j))+expnt3(1,iatnum(j)))**
     .       (2.d0*x+2.d0)
         else
            write(6,'("Error in dipole.F: no semiempirical method!")')
            stop
         endif
         hybfac(j)=const*frac1/frac2
      enddo

      counter=0
      do j=1,natoms
         if (iatnum(j).eq.1) then  ! s orbital
            counter=counter+1
         else                      ! s,p orbital
            hybdip(1)=hybdip(1)-hybfac(j)*pdiag(counter+2)
            hybdip(2)=hybdip(2)-hybfac(j)*pdiag(counter+4)
            hybdip(3)=hybdip(3)-hybfac(j)*pdiag(counter+7)
            dipole_atm(1,j) = dipole_atm(1,j)-hybfac(j)*pdiag(counter+2)
            dipole_atm(2,j) = dipole_atm(2,j)-hybfac(j)*pdiag(counter+4)
            dipole_atm(3,j) = dipole_atm(3,j)-hybfac(j)*pdiag(counter+7)
            counter=counter+10
         endif
      enddo

      dipol1=sqrt( (ptchg(1)+hybdip(1))**2
     .            +(ptchg(2)+hybdip(2))**2
     .            +(ptchg(3)+hybdip(3))**2)
      if (screen) then
        write (iout,*) ' '
        write (iout,66) dipol1,(ptchg(i)+hybdip(i),i=1,3)
        write (iout,*) ' '
      endif
      return
   66 format(' DIPOLE Total=',f12.6,' X=',f12.6,' Y=',f12.6,' Z=',f12.6)
      end
