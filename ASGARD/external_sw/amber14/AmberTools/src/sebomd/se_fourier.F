      subroutine se_fourn(data,nn,ndim,isign)
!     Fourier transformations

      implicit none
#include "sebomd.dim"

      integer isign,ndim,nn(ndim)
      double precision data(*)
      integer i1,i2,i2rev,i3,i3rev,ibit,idim,ifp1,ifp2,iip1,iip2,
     *      iip3,k1,k2,n,nprev,nrem,ntot
      double precision tempi,tempr
      double precision theta,wi,wpi,wpr,wr,wtemp

      ntot=1
      do idim=1,ndim
        ntot=ntot*nn(idim)
      end do
      nprev=1
      do idim=1,ndim
        n=nn(idim)
        nrem=ntot/(n*nprev)
        iip1=2*nprev
        iip2=iip1*n
        iip3=iip2*nrem
        i2rev=1
        do i2=1,iip2,iip1
          if(i2.lt.i2rev)then
            do i1=i2,i2+iip1-2,2
              do i3=i1,iip3,iip2
                i3rev=i2rev+i3-i2
                tempr=data(i3)
                tempi=data(i3+1)
                data(i3)=data(i3rev)
                data(i3+1)=data(i3rev+1)
                data(i3rev)=tempr
                data(i3rev+1)=tempi
              end do
            end do
          endif
          ibit=iip2/2
1         if ((ibit.ge.iip1).and.(i2rev.gt.ibit)) then
            i2rev=i2rev-ibit
            ibit=ibit/2
            go to 1
          endif
          i2rev=i2rev+ibit
        end do
        ifp1=iip1
2       if(ifp1.lt.iip2)then
          ifp2=2*ifp1
          theta=isign*pi2/(ifp2/iip1)
          wpr=-2.d0*dsin(0.5d0*theta)**2
          wpi=dsin(theta)
          wr=1.d0
          wi=0.d0
          do i3=1,ifp1,iip1
            do i1=i3,i3+iip1-2,2
              do i2=i1,iip3,ifp2
                k1=i2
                k2=k1+ifp1
                tempr=sngl(wr)*data(k2)-sngl(wi)*data(k2+1)
                tempi=sngl(wr)*data(k2+1)+sngl(wi)*data(k2)
                data(k2)=data(k1)-tempr
                data(k2+1)=data(k1+1)-tempi
                data(k1)=data(k1)+tempr
                data(k1+1)=data(k1+1)+tempi
              end do
            end do
            wtemp=wr
            wr=wr*wpr-wi*wpi+wr
            wi=wi*wpr+wtemp*wpi+wi
          end do 
          ifp1=ifp2
          go to 2
        endif
        nprev=n*nprev
      end do
      return
      end
