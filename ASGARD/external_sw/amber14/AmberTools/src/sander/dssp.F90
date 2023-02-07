#ifdef DSSP
#include "../include/dprec.fh"

module dssp

      integer :: npepc, idssp
      integer, allocatable :: ipepc(:)
      _REAL_  :: edssp

      _REAL_ :: v0 =1.d0      !  Fermi function step locations (as function of V(r) and ron) and widths
      _REAL_ :: wv = 1.d0
      _REAL_ :: ron0 = 1.d0
      _REAL_ :: won = 1.d0
      _REAL_ :: q = 80.d0           !  Q = 332 * dqC * dqN
      _REAL_ :: dssp_wt = 4.5d0     !  weight of the dssp terms

contains


      subroutine fdssp(natom,x,f,edssp)
   !  computes DSSP forces and energy

      implicit none

      integer, intent(in) :: natom       !  number of atoms
      _REAL_, intent(in) :: x(*)           !  coordinates
      _REAL_, intent(inout) :: f(*)        !  force vector, updated here
      _REAL_, intent(out) :: edssp         !  dssp energy term

      integer :: i,j
      _REAL_ :: xc(3),xo(3),xn(3),xh(3),xcn(3),xch(3),xon(3),xoh(3)
      _REAL_ :: rcn2,rch2,ron2,roh2,rcn,rch,ron,roh,rcn3,rch3,ron3,roh3
       
      _REAL_ :: u,v,fv,fr

      _REAL_ :: dvdxc(3),dvdxo(3),dvdxn(3),dvdxh(3),drdxo(3),drdxn(3)
      _REAL_ :: dfvdv,const,dfrdr


      edssp = 0.0

      !  loop over all pairs of peptide groups
         do i=1,npepc-1

         !  get the coordinates needed for the function: carbonyl of peptide "i"
            xc(1:3) = x( 3*(ipepc(i)  )-2 : 3*(ipepc(i)  ) )
            xo(1:3) = x( 3*(ipepc(i)+1)-2 : 3*(ipepc(i)+1) )

            do j=i+1,npepc

            !  NH group of peptide "j"
               xn(1:3) = x( 3*(ipepc(j)+2)-2 : 3*(ipepc(j)+2) )
               xh(1:3) = x( 3*(ipepc(j)+3)-2 : 3*(ipepc(j)+3) )

               xcn(1:3) = xc(1:3) - xn(1:3)
               xch(1:3) = xc(1:3) - xh(1:3)
               xon(1:3) = xo(1:3) - xn(1:3)
               xoh(1:3) = xo(1:3) - xh(1:3)

               rcn2 =  xcn(1)**2 + xcn(2)**2 + xcn(3)**2 
               rch2 =  xch(1)**2 + xch(2)**2 + xch(3)**2 
               ron2 =  xon(1)**2 + xon(2)**2 + xon(3)**2 
               roh2 =  xoh(1)**2 + xoh(2)**2 + xoh(3)**2

               rcn = sqrt(rcn2)
               rch = sqrt(rch2)
               ron = sqrt(ron2)
               roh = sqrt(roh2)

               rcn3 = rcn**3
               rch3 = rch**3
               ron3 = ron**3
               roh3 = roh**3

               call udssp(ron,rch,roh,rcn,u,v,fv,fr)
               edssp = edssp + dssp_wt*u
               if( idssp > 1 .and. abs(u) > 0.01 ) then
                  write(6,'(a,2i5,f10.4)') '      i,j, u = ', i,j,u
               end if

            !  update the forces:  compute DSSP force on all atoms in the peptide pair
            !  (i) C-O...H-N (j)

            !  since U = V * F(V) * F(ron), dU/dx = ( dV/dx ) * F(V) * F(ron) +
            !                                       V * { [ dF(V)/dV ] * [ dV/dx ] } * F(ron) +
            !                                       V * F(V) * [ dF(ron)/dx ] 


            !  add V * [ dF(V)/dx ] * F(r) term
            !  evaluate dV/dx [ = ( dV/dr ) * ( dr/dx ) ]
               dvdxc(1:3) =  q*( -xch( 1 : 3 )/rch3 + xcn( 1 : 3 )/rcn3 )
               dvdxo(1:3) =  q*( -xon( 1 : 3 )/ron3 + xoh( 1 : 3 )/roh3 )
               dvdxn(1:3) =  q*( xon( 1 : 3 )/ron3 - xcn( 1 : 3 )/rcn3 )
               dvdxh(1:3) =  q*( xch( 1 : 3 )/rch3 - xoh( 1 : 3 )/roh3 )

            !  evaluate dF(V)/dV
               dfvdv = ( exp( (v - v0)/wv ) )/wv
               dfvdv = -dfvdv*fv*fv

               const = v*dfvdv*fr*dssp_wt

               f( 3*(ipepc(i)  )-2 : 3*(ipepc(i)  ) ) = f( 3*(ipepc(i)  )-2 : 3*(ipepc(i)  ) ) &
                  - const*dvdxc( 1 : 3 )
               f( 3*(ipepc(i)+1)-2 : 3*(ipepc(i)+1) ) = f( 3*(ipepc(i)+1)-2 : 3*(ipepc(i)+1) ) &
                  - const*dvdxo( 1 : 3 )
               f( 3*(ipepc(j)+2)-2 : 3*(ipepc(j)+2) ) = f( 3*(ipepc(j)+2)-2 : 3*(ipepc(j)+2) ) &
                  - const*dvdxn( 1 : 3 )
               f( 3*(ipepc(j)+3)-2 : 3*(ipepc(j)+3) ) = f( 3*(ipepc(j)+3)-2 : 3*(ipepc(j)+3) ) &
                  - const*dvdxh( 1 : 3 )

            !  add V * F(V) * [ dF(r)/dx ] term:  since F(r) depends only on ron, only xo and xn derivatives
            !  of ron are needed
            !  evaluate dF(r)/dr and dr/dx
               dfrdr = ( exp( (ron - ron0)/won ) )/won
               dfrdr = -dfrdr*fr*fr
               drdxo( 1 : 3 ) = xon( 1 : 3 )/ron
               drdxn( 1 : 3 ) = -xon( 1 : 3 )/ron

               const = v*fv*dfrdr*dssp_wt

               f( 3*(ipepc(i)+1)-2 : 3*(ipepc(i)+1) ) = f( 3*(ipepc(i)+1)-2 : 3*(ipepc(i)+1) ) &
                  - const*drdxo( 1 : 3 )
               f( 3*(ipepc(j)+2)-2 : 3*(ipepc(j)+2) ) = f( 3*(ipepc(j)+2)-2 : 3*(ipepc(j)+2) ) &
                  - const*drdxn( 1 : 3 )

            !  add (dV/dx) * F(V) * F(r) term

               const = fv*fr*dssp_wt

               f( 3*(ipepc(i)  )-2 : 3*(ipepc(i)  ) ) = f( 3*(ipepc(i)  )-2 : 3*(ipepc(i)  ) ) &
                  - const*dvdxc( 1 : 3 )
               f( 3*(ipepc(i)+1)-2 : 3*(ipepc(i)+1) ) = f( 3*(ipepc(i)+1)-2 : 3*(ipepc(i)+1) ) &
                  - const*dvdxo( 1 : 3 )
               f( 3*(ipepc(j)+2)-2 : 3*(ipepc(j)+2) ) = f( 3*(ipepc(j)+2)-2 : 3*(ipepc(j)+2) ) &
                  - const*dvdxn( 1 : 3 )
               f( 3*(ipepc(j)+3)-2 : 3*(ipepc(j)+3) ) = f( 3*(ipepc(j)+3)-2 : 3*(ipepc(j)+3) ) &
                  - const*dvdxh( 1 : 3 )

           end do
        end do

      end subroutine fdssp


      subroutine udssp(ron,rch,roh,rcn,u,v,fv,fr)
      implicit none

      _REAL_, intent(out) :: u,v,fv,fr
      _REAL_, intent(in) :: ron,rch,roh,rcn

!     DSSP potential; separately evaluates cutoff (Fermi) functions in V and ron
!     F(V) and F(ron)

      call vdssp(ron,rch,roh,rcn,v)
      call fermi(v,v0,wv,fv)
      call fermi(ron,ron0,won,fr)

      u = v*fv*fr

      end subroutine udssp


      subroutine vdssp(ron,rch,roh,rcn,v)
      implicit none

      _REAL_,intent(in) :: ron,rch,roh,rcn

      _REAL_,intent(out) :: v

!     DSSP potential (aside from cutoff functions in V and ron)
!     V = 332 * Qc * Qn * ( 1/ron + 1/rch - 1/roh - 1/rcn )

      v = 1.d0/ron + 1.d0/rch - 1.d0/roh - 1.d0/rcn
      v = q*v

      end subroutine vdssp


      subroutine fermi(x,x0,w,f)
      implicit none

      _REAL_, intent(in) :: x,x0,w

      _REAL_, intent(out) :: f

!     Fermi function: { exp[ (x-x0)/w ] + 1 }^-1

      if (w > 0.d0) then
         f = ( exp( (x-x0)/w ) + 1.d0 )
         f = 1.d0/f
      else if (w <= 0.d0) then
         if (x <= x0) then
            f = 1.d0
         else if (x > x0) then
            f = 0.d0
         end if
      end if

      end subroutine fermi

end module dssp

#else

      subroutine dssp
      return
      end subroutine dssp

#endif
