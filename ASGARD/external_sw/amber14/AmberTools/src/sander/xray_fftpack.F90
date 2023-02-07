! <compile=optimized>
#include "../include/assert.fh"
module xray_fftpack_module
   !--------------------------------------------------------------------
   ! This is a Fortran90 conversion of FFTPACK 4.1, by Joe Krahn, 2006
   ! I tested FFTPACK 5, which has hand-optimizations, but found it
   ! slower. It probably interferes with modern optimizers, even though
   ! it may have helped when originally written.
   !
   ! Public routines are all wrappers using a derived type FFT descriptor.
   ! Actual FFTPACK routines are all converted to internal subroutines.
   ! Hopefully, some of these internal versions get in-lined.
   !
   ! ERROR return arguments are all optional. If an error occurs without
   ! an error argument, program execution is aborted.
   !
   ! SUBROUTINES:
   !---------------------------------------------------
   ! Allocate a 1D complex FFT descriptor object:
   ! fft1c_new (desc,size,error)
   !   type(fft1c_descriptor), intent(out) :: desc
   !   integer, intent(in) :: size(3)
   !   integer, intent(out), optional :: error
   !
   ! Delete a 1D complex FFT descriptor object:
   ! fft1c_delete (config)
   !
   ! Perform the 1D complex FFT:
   ! fft1c (config,c) ! This is a front-end to FFTPACK CFFTB()
   !
   !---------------------------------------------------
   ! fft1cr_new (config,n,error)
   ! fft1cr_delete (config)
   ! fft1cr (config,c) ! This is a front-end to FFTPACK RFFTB()
   !                   ! Similar to CNS's SFFT1CR
   !
   !---------------------------------------------------
   ! fft3c_new (config,size,error)
   ! fft3c_delete(config)
   ! fft3c(config,mx,my,nx,ny,nz,a,error)
   !
   !--------------------------------------------------------------------
   use xray_globals_module, only: real_kind, rk_, stdout
   implicit none
   public

   integer, parameter :: fft_kind = real_kind

   ! The first 2 typedefs differe only by name, to avoid accidental
   ! mixing of the two FFT types.
   type fft1c_descriptor
      integer :: size=0, nfac
      real(fft_kind), pointer :: wsave(:)=>NULL()
   end type fft1c_descriptor

   type fft1cr_descriptor
      integer :: size=0, nfac
      real(fft_kind), pointer :: wsave(:)=>NULL()
   end type fft1cr_descriptor

   type fft3c_descriptor
      type(fft1c_descriptor) :: fft(3)
   end type fft3c_descriptor

   real(fft_kind), parameter :: pi = &
     3.1415926535897932384626433832795028841971693993_fft_kind
   real(fft_kind), parameter :: twopi = pi * 2.0_fft_kind
   real(fft_kind), parameter :: taur = -0.5_fft_kind
   real(fft_kind), parameter :: taui = &
     0.8660254037844386467637231707529361834714026269_fft_kind ! sqrt(3)/2
   real(fft_kind), parameter :: ti11 = &
     0.9510565162951535721164393333793821434056986341_fft_kind ! sin(2*pi*(1/5))
   real(fft_kind), parameter :: ti12 = &
     0.5877852522924731291687059546390727685976524376_fft_kind ! sin(2*pi*(2/5))
   real(fft_kind), parameter :: tr11 = &
     0.3090169943749474241022934171828190588601545899_fft_kind ! cos(2*pi*(1/5))
   real(fft_kind), parameter :: tr12 = &
     -0.809016994374947424102293417182819058860154590_fft_kind ! cos(2*pi*(2/5))

   type(fft3c_descriptor), save :: fft3c_state

   integer, parameter :: fft_max_prime=5

contains
   !======================================================================
   ! FFT routines derived from FFTPACK 4.1
   ! =====================================================================
   subroutine fft1c (desc,c) ! This is FFTPACK DCFFTB()
      type(fft1c_descriptor) :: desc
      complex(fft_kind) :: c(*)
      if (desc%size == 1) return
      call fftpack_cfftb(desc%size,c,desc%wsave, &
             desc%nfac,desc%wsave(desc%size*2+1))
      return
   end subroutine fft1c

   subroutine fft1c_new (desc,size,error)
      type(fft1c_descriptor), intent(out) :: desc
      integer, intent(in) :: size
      logical, intent(inout), optional :: error
      integer :: alloc_status
      if (present(error)) error=.false.
      if (desc%size>1) call fft1c_delete(desc)
      desc%size=size
      if (size == 1) return
      allocate(desc%wsave(size*2+15),stat=alloc_status)
      REQUIRE(alloc_status==0)
      call fftpack_cffti(size,desc%wsave,desc%nfac,desc%wsave(desc%size*2+1))
      return
   end subroutine fft1c_new

   subroutine fft1c_delete (desc)
      type(fft1c_descriptor) :: desc
      if (desc%size>1) then
         if (.not.associated(desc%wsave)) stop 'BUG in fft1c_delete'
         deallocate(desc%wsave)
      end if
      desc%size=0
   end subroutine fft1c_delete

   !=======================================================================
   ! 1D complex-to-real transform (in-place)
   ! NOTE: complex-to-real requires an even first dimension

   subroutine fft1cr (desc,c)
      type(fft1cr_descriptor) :: desc
      real(fft_kind) :: c(0:*) ! complex in, real out
      external :: fftpack_rfftb
      !local
      integer :: n
      n=desc%size
      if (n <= 1) return
      c(1:n) = c(2:n+1)
      call fftpack_rfftb(desc%size,c,desc%wsave, &
             desc%nfac,desc%wsave(desc%size*2+1))
      return
   end subroutine fft1cr

   subroutine fft1cr_new (desc,n,error)
      type(fft1cr_descriptor), intent(out) :: desc
      integer, intent(in) :: n
      logical, intent(inout), optional :: error
      integer :: alloc_status
      if (present(error)) error=.false.
      if (desc%size>1) call fft1cr_delete(desc)
      desc%size=n
      if (n == 1) return
      allocate(desc%wsave(n*2+15),stat=alloc_status)
      REQUIRE(alloc_status==0)
      call fftpack_crffti(n,desc%wsave,desc%nfac,desc%wsave(desc%size*2+1))
      return
   end subroutine fft1cr_new

   subroutine fft1cr_delete (desc)
      type(fft1cr_descriptor) :: desc
      if (desc%size>1) then
         if (.not.associated(desc%wsave)) stop 'BUG in fft1cr_delete'
         deallocate(desc%wsave)
      end if
      desc%size=0
   end subroutine fft1cr_delete
   !=======================================================================
   ! END OF FFTPACK
   !=======================================================================

   subroutine fft3c_new(desc,size,error)
      type(fft3c_descriptor), intent(out) :: desc
      integer, intent(in) :: size(3)
      logical, intent(inout), optional :: error
      logical :: error_
      integer :: i,j
      error_=.false.
      do i=1,3
         call fft1c_new(desc%fft(i),size(i),error_)
         if (error_) then
            do j=1,i-1
               call fft1c_delete(desc%fft(j))
            end do
            exit
         end if
      end do
      if (present(error)) then
         error=error_
      else if (error_) then
         write(stdout,'(A)') 'ERROR in fft3c_new()'
         call mexit(stdout,1)
      end if
   end subroutine fft3c_new

   subroutine fft3c(desc,cmap,error)
      implicit none
      type(fft3c_descriptor), intent(inout) :: desc
      complex(fft_kind), intent(inout) :: cmap(:,:,:)
      logical, intent(out) :: error
      ! local
      integer :: u,v,w
      ! begin
      error = .false.

      !$OMP parallel default(none)  &
            !$OMP     shared(cmap) private(u,v,w) shared(desc)
      !$OMP do
      do u = lbound(cmap,1),ubound(cmap,1)
         do v = lbound(cmap,2),ubound(cmap,2)
            ! Transform along z (slow direction)
            call fft1c(desc%fft(3), cmap(u,v,:))
         end do
         do w = lbound(cmap,3),ubound(cmap,3)
            ! Transform along y (medium direction)
            call fft1c(desc%fft(2), cmap(u,:,w))
         end do
      end do
      !$OMP end do
      !$OMP end parallel

      !$OMP parallel do default(none) &
            !$OMP     shared(cmap), private(v,w), shared(desc)
      do w = lbound(cmap,3),ubound(cmap,3)
         do v = lbound(cmap,2),ubound(cmap,2)
            ! Transform along x (fast direction)
            call fft1c(desc%fft(1), cmap(:,v,w))
         end do
      end do
      !$OMP end parallel do
   end subroutine fft3c

   subroutine fft3c_delete(desc)
      type(fft3c_descriptor) :: desc
      integer :: i
      do i=1,3
         call fft1c_delete(desc%fft(i))
      end do
      return
   end subroutine fft3c_delete

end module xray_fftpack_module

!=====================================================================

subroutine fftpack_crffti (n,wa,nf,ifac)
   use xray_fftpack_module
   implicit none
   integer :: n
   real(fft_kind) :: wa(n*2)
   integer :: nf,ifac(15)
   ! local variables
   real(fft_kind) :: arg,argh,argld,fi
   integer :: i,ib,ido,ii,ip,is,j,k1,l1,l2,ld,nl,nq,nr,ntry
   integer, parameter :: ntryh(4)=(/4,2,3,5/)
   nl = n
   nf = 0
   j = 0

   LOOP1: do
      j = j+1
      if (j<=4) then
         ntry = ntryh(j)
      else
         ntry = ntry+2
      end if

      LOOP2: do
         nq = nl/ntry
         nr = nl-ntry*nq
         if (nr/=0) cycle LOOP1
         nf = nf+1
         ifac(nf) = ntry
         nl = nq
         if (ntry==2 .and. nf/=1) then
            do i=2,nf
               ib = nf-i+2
               ifac(ib) = ifac(ib-1)
            end do
            ifac(1) = 2
         end if
         if (nl == 1) exit LOOP1
      end do LOOP2
   end do LOOP1

   argh = twopi/float(n)
   is = 0
   l1 = 1
   do k1=1,nf-1
      ip = ifac(k1)
      ld = 0
      l2 = l1*ip
      ido = n/l2
      do j=1,ip-1
         ld = ld+l1
         i = is
         argld = float(ld)*argh
         fi = 0.
         do ii=3,ido,2
            i = i+2
            fi = fi+1.
            arg = fi*argld
            wa(i-1) = cos(arg)
            wa(i) = sin(arg)
         end do
         is = is+ido
      end do
      l1 = l2
   end do  ! k1=1,nf-1
   return
end subroutine fftpack_crffti

subroutine fftpack_cffti(n,wa,nf,ifac)
   use xray_fftpack_module
   implicit none
   integer :: n
   real(fft_kind) :: wa(n*2)
   integer :: nf,ifac(15)
   ! local variables
   real(fft_kind) :: arg,argh,argld,fi
   integer :: i,i1,ib,ido,idot,ii,ip, &
         j,k1,l1,l2,ld,nl,nq,nr,ntry
   integer, parameter :: ntryh(4)=(/3,4,2,5/)
   nl = n
   nf = 0
   j = 0

   LOOP1: do
      j = j+1
      if (j<=4) then
         ntry = ntryh(j)
      else
         ntry = ntry+2
      end if

      LOOP2: do
         nq = nl/ntry
         nr = nl-ntry*nq
         if (nr/=0) cycle LOOP1
         nf = nf+1
         ifac(nf) = ntry
         nl = nq
         if (ntry==2 .and. nf/=1) then
            do i=2,nf
               ib = nf-i+2
               ifac(ib) = ifac(ib-1)
            end do
            ifac(1) = 2
         end if
         if (nl == 1) exit LOOP1
      end do LOOP2
   end do LOOP1

   argh = twopi/real(n,fft_kind)
   i = 2
   l1 = 1
   do k1=1,nf
      ip = ifac(k1)
      ld = 0
      l2 = l1*ip
      ido = n/l2
      idot = ido+ido+2
      do j=1,ip-1
         i1 = i
         wa(i-1) = 1.
         wa(i) = 0.
         ld = ld+l1
         fi = 0.
         argld = real(ld,fft_kind)*argh
         do ii=4,idot,2
            i = i+2
            fi = fi+1.
            arg = fi*argld
            wa(i-1) = cos(arg)
            wa(i) = sin(arg)
         end do
         if (ip > 5) then
            wa(i1-1) = wa(i-1)
            wa(i1) = wa(i)
         end if
      end do
      l1 = l2
   end do ! k1=1,nf
   return
end subroutine fftpack_cffti

subroutine fftpack_cfftb (n,c,wa,nf,ifac)
   use xray_fftpack_module
   implicit none
   integer :: n
   real(fft_kind) :: c(n*2),wa(n*2)
   integer :: nf,ifac(15)
   ! local variables
   real(fft_kind) :: ch(n*2) ! automatic scratch array
   integer :: idl1,ido,idot,ip,iw,ix2,ix3,ix4,k1,l1,l2,nac
   logical :: forward
   forward=.true.
   l1 = 1
   iw = 1
   do k1=1,nf
      ip = ifac(k1)
      l2 = ip*l1
      ido = n/l2
      idot = ido+ido
      idl1 = idot*l1
      select case(ip)
      case(4)
      ix2 = iw+idot
      ix3 = ix2+idot
      if (forward) then
         call passb4 (idot,l1,c,ch,wa(iw),wa(ix2),wa(ix3))
      else
         call passb4 (idot,l1,ch,c,wa(iw),wa(ix2),wa(ix3))
      end if
      forward=.not.forward

      case(2)
      if (forward) then
         call passb2 (idot,l1,c,ch,wa(iw))
      else
         call passb2 (idot,l1,ch,c,wa(iw))
      end if
      forward=.not.forward

      case(3)
      ix2 = iw+idot
      if (forward) then
         call passb3 (idot,l1,c,ch,wa(iw),wa(ix2))
      else
         call passb3 (idot,l1,ch,c,wa(iw),wa(ix2))
      end if
      forward=.not.forward

      case(5)
      ix2 = iw+idot
      ix3 = ix2+idot
      ix4 = ix3+idot
      if (forward) then
         call passb5 (idot,l1,c,ch,wa(iw),wa(ix2),wa(ix3),wa(ix4))
      else
         call passb5 (idot,l1,ch,c,wa(iw),wa(ix2),wa(ix3),wa(ix4))
      end if
      forward=.not.forward

      case default
      if (forward) then
         call passb (nac,idot,ip,l1,idl1,c,c,c,ch,ch,wa(iw))
      else
         call passb (nac,idot,ip,l1,idl1,ch,ch,ch,c,c,wa(iw))
      end if
      if (nac /= 0) forward=.not.forward

      end select

      l1 = l2
      iw = iw+(ip-1)*idot
   end do  ! k1=1,nf
   if (.not.forward) c(1:n*2) = ch(1:n*2)
   return
contains
   subroutine passb (nac,ido,ip,l1,idl1,cc,c1,c2,ch,ch2,wa)
      implicit none
      integer :: nac,ido,ip,l1,idl1
      real(fft_kind) :: cc(ido, ip, l1),c1(ido, l1, ip)
      real(fft_kind) :: c2(idl1, ip),ch(ido, l1, ip)
      real(fft_kind) :: ch2(idl1, ip),wa(*)
      ! local variables
      integer :: i,idij,idj,idl,idlj,idot,idp,ik,inc,ipp2,ipph,j,jc,k,l,lc,nt
      real(fft_kind) :: wai,war
      idot = ido/2
      nt = ip*idl1
      ipp2 = ip+2
      ipph = (ip+1)/2
      idp = ip*ido

      if (ido >= l1) then
         do j=2,ipph
            jc = ipp2-j
            do k=1,l1
               do i=1,ido
                  ch(i  ,k,j ) = cc(i  ,j,k) + cc(i  ,jc,k)
                  ch(i  ,k,jc) = cc(i  ,j,k) - cc(i  ,jc,k)
               end do
            end do
         end do
         do k=1,l1
            do i=1,ido
               ch(i  ,k,1) = cc(i  ,1,k)
            end do
         end do
      else
         do j=2,ipph
            jc = ipp2-j
            do i=1,ido
               do k=1,l1
                  ch(i  ,k,j ) = cc(i  ,j,k) + cc(i  ,jc,k)
                  ch(i  ,k,jc) = cc(i  ,j,k) - cc(i  ,jc,k)
               end do
            end do
         end do
         do i=1,ido
            do k=1,l1
               ch(i  ,k,1) = cc(i  ,1,k)
            end do
         end do
      end if
      idl = 2-ido
      inc = 0
      do l=2,ipph
         lc = ipp2-l
         idl = idl+ido
         do ik=1,idl1
            c2(ik,l ) = ch2(ik,1) + wa(idl-1)*ch2(ik,2)
            c2(ik,lc) = wa(idl)*ch2(ik,ip)
         end do
         idlj = idl
         inc = inc+ido
         do j=3,ipph
            jc = ipp2-j
            idlj = idlj+inc
            if (idlj > idp) idlj = idlj-idp
            war = wa(idlj-1)
            wai = wa(idlj)
            do ik=1,idl1
               c2(ik,l ) = c2(ik,l ) + war*ch2(ik,j )
               c2(ik,lc) = c2(ik,lc) + wai*ch2(ik,jc)
            end do
         end do
      end do  ! l=2,ipph
      do j=2,ipph
         do ik=1,idl1
            ch2(ik,1) = ch2(ik,1) + ch2(ik,j )
         end do
      end do
      do j=2,ipph
         jc = ipp2-j
         do ik=2,idl1,2
            ch2(ik-1,j ) = c2(ik-1,j ) - c2(ik,jc)
            ch2(ik-1,jc) = c2(ik-1,j ) + c2(ik,jc)
            ch2(ik,j ) = c2(ik,j ) + c2(ik-1,jc)
            ch2(ik,jc) = c2(ik,j ) - c2(ik-1,jc)
         end do
      end do
      nac = 1
      if (ido == 2) return
      nac = 0
      c2(1:idl1,1) = ch2(1:idl1,1)
      do j=2,ip
         do k=1,l1
            c1(1,k,j ) = ch(1,k,j )
            c1(2,k,j ) = ch(2,k,j )
         end do
      end do

      if (idot <= l1) then
         idij = 0
         do j=2,ip
            idij = idij+2
            do i=4,ido,2
               idij = idij+2
               do k=1,l1
                  c1(i-1,k,j ) = wa(idij-1)*ch(i-1,k,j ) - wa(idij)*ch(i  ,k,j )
                  c1(i  ,k,j ) = wa(idij-1)*ch(i  ,k,j ) + wa(idij)*ch(i-1,k,j )
               end do
            end do
         end do
      else
         idj = 2-ido
         do j=2,ip
            idj = idj+ido
            do k=1,l1
               idij = idj
               do i=4,ido,2
                  idij = idij+2
                  c1(i-1,k,j ) = wa(idij-1)*ch(i-1,k,j ) - wa(idij)*ch(i  ,k,j )
                  c1(i  ,k,j ) = wa(idij-1)*ch(i  ,k,j ) + wa(idij)*ch(i-1,k,j )
               end do
            end do
         end do
      end if
      return
   end subroutine passb
   !=======================================================================
   subroutine passb2 (ido,l1,cc,ch,wa1)
      implicit none
      integer :: ido,l1
      real(fft_kind) :: cc(ido, 2, l1),ch(ido, l1, 2),wa1(*)
      ! local variables
      integer :: i,k
      real(fft_kind) :: ti2,tr2
      if (ido <= 2) then
         do k=1,l1
            ch(1,k,1) = cc(1,1,k) + cc(1,2,k)
            ch(1,k,2) = cc(1,1,k) - cc(1,2,k)
            ch(2,k,1) = cc(2,1,k) + cc(2,2,k)
            ch(2,k,2) = cc(2,1,k) - cc(2,2,k)
         end do
      else
         do k=1,l1
            do i=2,ido,2
               ch(i-1,k,1) = cc(i-1,1,k) + cc(i-1,2,k)
               tr2 = cc(i-1,1,k) - cc(i-1,2,k)
               ch(i  ,k,1) = cc(i  ,1,k) + cc(i  ,2,k)
               ti2 = cc(i  ,1,k) - cc(i  ,2,k)
               ch(i  ,k,2) = wa1(i-1)*ti2+wa1(i)*tr2
               ch(i-1,k,2) = wa1(i-1)*tr2-wa1(i)*ti2
            end do
         end do
      end if
      return
   end subroutine passb2
   !=======================================================================
   subroutine passb3 (ido,l1,cc,ch,wa1,wa2)
      implicit none
      integer :: ido,l1
      real(fft_kind) :: cc(ido, 3, l1),ch(ido, l1, 3),wa1(*),wa2(*)
      ! local variables
      real(fft_kind) :: ci2,ci3,cr2,cr3,di2,di3,dr2,dr3
      integer :: i,k
      real(fft_kind) :: ti2,tr2
      if (ido == 2) then
         do k=1,l1
            tr2 = cc(1,2,k) + cc(1,3,k)
            cr2 = cc(1,1,k) + taur*tr2
            ch(1,k,1) = cc(1,1,k) + tr2
            ti2 = cc(2,2,k) + cc(2,3,k)
            ci2 = cc(2,1,k) + taur*ti2
            ch(2,k,1) = cc(2,1,k) + ti2
            cr3 = taui*(cc(1,2,k) - cc(1,3,k))
            ci3 = taui*(cc(2,2,k) - cc(2,3,k))
            ch(1,k,2) = cr2-ci3
            ch(1,k,3) = cr2+ci3
            ch(2,k,2) = ci2+cr3
            ch(2,k,3) = ci2-cr3
         end do
      else
         do k=1,l1
            do i=2,ido,2
               tr2 = cc(i-1,2,k) + cc(i-1,3,k)
               cr2 = cc(i-1,1,k) + taur*tr2
               ch(i-1,k,1) = cc(i-1,1,k) + tr2
               ti2 = cc(i  ,2,k) + cc(i  ,3,k)
               ci2 = cc(i  ,1,k) + taur*ti2
               ch(i  ,k,1) = cc(i  ,1,k) + ti2
               cr3 = taui*(cc(i-1,2,k) - cc(i-1,3,k))
               ci3 = taui*(cc(i  ,2,k) - cc(i  ,3,k))
               dr2 = cr2-ci3
               dr3 = cr2+ci3
               di2 = ci2+cr3
               di3 = ci2-cr3
               ch(i  ,k,2) = wa1(i-1)*di2+wa1(i)*dr2
               ch(i-1,k,2) = wa1(i-1)*dr2-wa1(i)*di2
               ch(i  ,k,3) = wa2(i-1)*di3+wa2(i)*dr3
               ch(i-1,k,3) = wa2(i-1)*dr3-wa2(i)*di3
            end do
         end do
      end if
      return
   end subroutine passb3
   !=======================================================================
   subroutine passb4 (ido,l1,cc,ch,wa1,wa2,wa3)
      implicit none
      integer :: ido,l1
      real(fft_kind) :: cc(ido, 4, l1),ch(ido, l1, 4),wa1(*),wa2(*),wa3(*)
      ! local variables
      integer :: i,k
      real(fft_kind) :: ci2,ci3,ci4, cr2,cr3,cr4
      real(fft_kind) :: ti1,ti2,ti3,ti4, tr1,tr2,tr3,tr4
      if (ido == 2) then
         do k=1,l1
            ti1 = cc(2,1,k) - cc(2,3,k)
            ti2 = cc(2,1,k) + cc(2,3,k)
            tr4 = cc(2,4,k) - cc(2,2,k)
            ti3 = cc(2,2,k) + cc(2,4,k)
            tr1 = cc(1,1,k) - cc(1,3,k)
            tr2 = cc(1,1,k) + cc(1,3,k)
            ti4 = cc(1,2,k) - cc(1,4,k)
            tr3 = cc(1,2,k) + cc(1,4,k)
            ch(1,k,1) = tr2+tr3
            ch(1,k,3) = tr2-tr3
            ch(2,k,1) = ti2+ti3
            ch(2,k,3) = ti2-ti3
            ch(1,k,2) = tr1+tr4
            ch(1,k,4) = tr1-tr4
            ch(2,k,2) = ti1+ti4
            ch(2,k,4) = ti1-ti4
         end do
      else
         do k=1,l1
            do i=2,ido,2
               ti1 = cc(i  ,1,k) - cc(i  ,3,k)
               ti2 = cc(i  ,1,k) + cc(i  ,3,k)
               ti3 = cc(i  ,2,k) + cc(i  ,4,k)
               tr4 = cc(i  ,4,k) - cc(i  ,2,k)
               tr1 = cc(i-1,1,k) - cc(i-1,3,k)
               tr2 = cc(i-1,1,k) + cc(i-1,3,k)
               ti4 = cc(i-1,2,k) - cc(i-1,4,k)
               tr3 = cc(i-1,2,k) + cc(i-1,4,k)
               ch(i-1,k,1) = tr2+tr3
               cr3 = tr2-tr3
               ch(i  ,k,1) = ti2+ti3
               ci3 = ti2-ti3
               cr2 = tr1+tr4
               cr4 = tr1-tr4
               ci2 = ti1+ti4
               ci4 = ti1-ti4
               ch(i-1,k,2) = wa1(i-1)*cr2-wa1(i)*ci2
               ch(i  ,k,2) = wa1(i-1)*ci2+wa1(i)*cr2
               ch(i-1,k,3) = wa2(i-1)*cr3-wa2(i)*ci3
               ch(i  ,k,3) = wa2(i-1)*ci3+wa2(i)*cr3
               ch(i-1,k,4) = wa3(i-1)*cr4-wa3(i)*ci4
               ch(i  ,k,4) = wa3(i-1)*ci4+wa3(i)*cr4
            end do  ! i=2,ido,2
         end do  ! k=1,l1
      end if
      return
   end subroutine passb4
   !=======================================================================
   subroutine passb5 (ido,l1,cc,ch,wa1,wa2,wa3,wa4)
      implicit none
      integer :: ido,l1
      real(fft_kind) :: cc(ido,5,l1),ch(ido,l1,5),wa1(*),wa2(*),wa3(*),wa4(*)
      ! local variables
      integer :: i,k
      real(fft_kind) :: ci2,ci3,ci4,ci5, cr2,cr3,cr4,cr5
      real(fft_kind) :: di2,di3,di4,di5, dr2,dr3,dr4,dr5
      real(fft_kind) :: ti2,ti3,ti4,ti5, tr2,tr3,tr4,tr5
      if (ido == 2) then
         do k=1,l1
            ti5 = cc(2,2,k) - cc(2,5,k)
            ti2 = cc(2,2,k) + cc(2,5,k)
            ti4 = cc(2,3,k) - cc(2,4,k)
            ti3 = cc(2,3,k) + cc(2,4,k)
            tr5 = cc(1,2,k) - cc(1,5,k)
            tr2 = cc(1,2,k) + cc(1,5,k)
            tr4 = cc(1,3,k) - cc(1,4,k)
            tr3 = cc(1,3,k) + cc(1,4,k)
            ch(1,k,1) = cc(1,1,k) + tr2+tr3
            ch(2,k,1) = cc(2,1,k) + ti2+ti3
            cr2 = cc(1,1,k) + tr11*tr2+tr12*tr3
            ci2 = cc(2,1,k) + tr11*ti2+tr12*ti3
            cr3 = cc(1,1,k) + tr12*tr2+tr11*tr3
            ci3 = cc(2,1,k) + tr12*ti2+tr11*ti3
            cr5 = ti11*tr5+ti12*tr4
            ci5 = ti11*ti5+ti12*ti4
            cr4 = ti12*tr5-ti11*tr4
            ci4 = ti12*ti5-ti11*ti4
            ch(1,k,2) = cr2-ci5
            ch(1,k,5) = cr2+ci5
            ch(2,k,2) = ci2+cr5
            ch(2,k,3) = ci3+cr4
            ch(1,k,3) = cr3-ci4
            ch(1,k,4) = cr3+ci4
            ch(2,k,4) = ci3-cr4
            ch(2,k,5) = ci2-cr5
         end do  ! k=1,l1
      else
         do k=1,l1
            do i=2,ido,2
               ti5 = cc(i  ,2,k) - cc(i  ,5,k)
               ti2 = cc(i  ,2,k) + cc(i  ,5,k)
               ti4 = cc(i  ,3,k) - cc(i  ,4,k)
               ti3 = cc(i  ,3,k) + cc(i  ,4,k)
               tr5 = cc(i-1,2,k) - cc(i-1,5,k)
               tr2 = cc(i-1,2,k) + cc(i-1,5,k)
               tr4 = cc(i-1,3,k) - cc(i-1,4,k)
               tr3 = cc(i-1,3,k) + cc(i-1,4,k)
               ch(i-1,k,1) = cc(i-1,1,k) + tr2+tr3
               ch(i  ,k,1) = cc(i  ,1,k) + ti2+ti3
               cr2 = cc(i-1,1,k) + tr11*tr2+tr12*tr3
               ci2 = cc(i  ,1,k) + tr11*ti2+tr12*ti3
               cr3 = cc(i-1,1,k) + tr12*tr2+tr11*tr3
               ci3 = cc(i  ,1,k) + tr12*ti2+tr11*ti3
               cr5 = ti11*tr5+ti12*tr4
               ci5 = ti11*ti5+ti12*ti4
               cr4 = ti12*tr5-ti11*tr4
               ci4 = ti12*ti5-ti11*ti4
               dr3 = cr3-ci4
               dr4 = cr3+ci4
               di3 = ci3+cr4
               di4 = ci3-cr4
               dr5 = cr2+ci5
               dr2 = cr2-ci5
               di5 = ci2-cr5
               di2 = ci2+cr5
               ch(i-1,k,2) = wa1(i-1)*dr2-wa1(i)*di2
               ch(i  ,k,2) = wa1(i-1)*di2+wa1(i)*dr2
               ch(i-1,k,3) = wa2(i-1)*dr3-wa2(i)*di3
               ch(i  ,k,3) = wa2(i-1)*di3+wa2(i)*dr3
               ch(i-1,k,4) = wa3(i-1)*dr4-wa3(i)*di4
               ch(i  ,k,4) = wa3(i-1)*di4+wa3(i)*dr4
               ch(i-1,k,5) = wa4(i-1)*dr5-wa4(i)*di5
               ch(i  ,k,5) = wa4(i-1)*di5+wa4(i)*dr5
            end do  ! i=2,ido,2
         end do  ! k=1,l1
      end if
      return
   end subroutine passb5
end subroutine fftpack_cfftb

subroutine fftpack_rfftb(n,c,wa,nf,ifac)
   use xray_fftpack_module
   implicit none
   integer :: n
   real(fft_kind) :: c(n*2),wa(n*2)
   integer :: nf,ifac(15)
   ! local variables
   real(fft_kind) :: ch(n*2) ! automatic scratch array
   integer :: idl1,ido,ip,iw,ix2,ix3,ix4,k1,l1,l2
   logical :: forward ! was na
   forward = .true.
   l1 = 1
   iw = 1
   do k1=1,nf
      ip = ifac(k1)
      l2 = ip*l1
      ido = n/l2
      idl1 = ido*l1
      select case(ip)
      case(4)
      ix2 = iw+ido
      ix3 = ix2+ido
      if (forward) then
         call radb4 (ido,l1,c,ch,wa(iw),wa(ix2),wa(ix3))
      else
         call radb4 (ido,l1,ch,c,wa(iw),wa(ix2),wa(ix3))
      end if
      forward=.not.forward

      case(2)
      if (forward) then
         call radb2 (ido,l1,c,ch,wa(iw))
      else
         call radb2 (ido,l1,ch,c,wa(iw))
      end if
      forward=.not.forward

      case(3)
      ix2 = iw+ido
      if (forward) then
         call radb3 (ido,l1,c,ch,wa(iw),wa(ix2))
      else
         call radb3 (ido,l1,ch,c,wa(iw),wa(ix2))
      end if
      forward=.not.forward

      case(5)
      ix2 = iw+ido
      ix3 = ix2+ido
      ix4 = ix3+ido
      if (forward) then
         call radb5 (ido,l1,c,ch,wa(iw),wa(ix2),wa(ix3),wa(ix4))
      else
         call radb5 (ido,l1,ch,c,wa(iw),wa(ix2),wa(ix3),wa(ix4))
      end if
      forward=.not.forward

      case default
      if (forward) then
         call radbg (ido,ip,l1,idl1,c,c,c,ch,ch,wa(iw))
      else
         call radbg (ido,ip,l1,idl1,ch,ch,ch,c,c,wa(iw))
      end if
      if (ido==1) forward=.not.forward

      end select

      l1 = l2
      iw = iw+(ip-1)*ido
   end do ! k1=1,nf

   if (.not.forward) c(1:n) = ch(1:n)
   return

contains
   subroutine radb2 (ido,l1,cc,ch,wa1)
      integer :: ido,l1
      real(fft_kind) :: cc(ido, 2, l1),ch(ido, l1, 2),wa1(*)
      ! local variables
      integer :: i,ic,idp2,k
      real(fft_kind) :: ti2,tr2
      do k=1,l1
         ch(1,k,1) = cc(1,1,k) + cc(ido,2,k)
         ch(1,k,2) = cc(1,1,k) - cc(ido,2,k)
      end do
      if (ido<2) return
      if (ido/=2) then
         idp2 = ido+2
         do k=1,l1
            do i=3,ido,2
               ic = idp2-i
               ch(i-1,k,1) = cc(i-1,1,k) + cc(ic-1,2,k)
               tr2 = cc(i-1,1,k) - cc(ic-1,2,k)
               ch(i  ,k,1) = cc(i  ,1,k) - cc(ic  ,2,k)
               ti2 = cc(i  ,1,k) + cc(ic  ,2,k)
               ch(i-1,k,2) = wa1(i-2)*tr2-wa1(i-1)*ti2
               ch(i  ,k,2) = wa1(i-2)*ti2+wa1(i-1)*tr2
            end do
         end do
         if (mod(ido,2) == 1) return
      end if
      do k=1,l1
         ch(ido,k,1) = cc(ido,1,k) + cc(ido,1,k)
         ch(ido,k,2) = -(cc(1,2,k) + cc(1,2,k))
      end do
   end subroutine radb2
   !=======================================================================
   subroutine radb3 (ido,l1,cc,ch,wa1,wa2)
      integer :: ido,l1
      real(fft_kind) :: cc(ido, 3, l1),ch(ido, l1, 3),wa1(*),wa2(*)
      ! local variables
      real(fft_kind) :: ci2,ci3,cr2,cr3,di2,di3,dr2,dr3
      integer :: i,ic,idp2,k
      real(fft_kind) :: ti2,tr2
      do k=1,l1
         tr2 = cc(ido,2,k) + cc(ido,2,k)
         cr2 = cc(1,1,k) + taur*tr2
         ch(1,k,1) = cc(1,1,k) + tr2
         ci3 = taui*(cc(1,3,k) + cc(1,3,k))
         ch(1,k,2) = cr2-ci3
         ch(1,k,3) = cr2+ci3
      end do
      if (ido == 1) return
      idp2 = ido+2
      do k=1,l1
         do i=3,ido,2
            ic = idp2-i
            tr2 = cc(i-1,3,k) + cc(ic-1,2,k)
            cr2 = cc(i-1,1,k) + taur*tr2
            ch(i-1,k,1) = cc(i-1,1,k) + tr2
            ti2 = cc(i  ,3,k) - cc(ic  ,2,k)
            ci2 = cc(i  ,1,k) + taur*ti2
            ch(i  ,k,1) = cc(i  ,1,k) + ti2
            cr3 = taui*(cc(i-1,3,k) - cc(ic-1,2,k))
            ci3 = taui*(cc(i  ,3,k) + cc(ic  ,2,k))
            dr2 = cr2-ci3
            dr3 = cr2+ci3
            di2 = ci2+cr3
            di3 = ci2-cr3
            ch(i-1,k,2) = wa1(i-2)*dr2-wa1(i-1)*di2
            ch(i  ,k,2) = wa1(i-2)*di2+wa1(i-1)*dr2
            ch(i-1,k,3) = wa2(i-2)*dr3-wa2(i-1)*di3
            ch(i  ,k,3) = wa2(i-2)*di3+wa2(i-1)*dr3
         end do
      end do  ! k=1,l1
      return
   end subroutine radb3
   !=======================================================================
   subroutine radb4 (ido,l1,cc,ch,wa1,wa2,wa3)
      integer, intent(in) :: ido,l1
      real(fft_kind), intent(in) :: cc(ido, 4, l1),wa1(*),wa2(*),wa3(*)
      real(fft_kind), intent(out) :: ch(ido, l1, 4)
      ! local variables
      integer :: i,ic,idp2,k
      real(fft_kind) :: ci2,ci3,ci4, cr2,cr3,cr4
      real(fft_kind) :: ti1,ti2,ti3,ti4, tr1,tr2,tr3,tr4
      real(fft_kind), parameter :: sqrt2 = &
            1.414213562373095048801688724209698_fft_kind
      do k=1,l1
         tr1 = cc(1,1,k) - cc(ido,4,k)
         tr2 = cc(1,1,k) + cc(ido,4,k)
         tr3 = cc(ido,2,k) + cc(ido,2,k)
         tr4 = cc(1,3,k) + cc(1,3,k)
         ch(1,k,1) = tr2+tr3
         ch(1,k,2) = tr1-tr4
         ch(1,k,3) = tr2-tr3
         ch(1,k,4) = tr1+tr4
      end do
      if (ido<2) return
      if (ido/=2) then
         idp2 = ido+2
         do k=1,l1
            do i=3,ido,2
               ic = idp2-i
               ti1 = cc(i  ,1,k) + cc(ic  ,4,k)
               ti2 = cc(i  ,1,k) - cc(ic  ,4,k)
               ti3 = cc(i  ,3,k) - cc(ic  ,2,k)
               tr4 = cc(i  ,3,k) + cc(ic  ,2,k)
               tr1 = cc(i-1,1,k) - cc(ic-1,4,k)
               tr2 = cc(i-1,1,k) + cc(ic-1,4,k)
               ti4 = cc(i-1,3,k) - cc(ic-1,2,k)
               tr3 = cc(i-1,3,k) + cc(ic-1,2,k)
               ch(i-1,k,1) = tr2+tr3
               cr3 = tr2-tr3
               ch(i  ,k,1) = ti2+ti3
               ci3 = ti2-ti3
               cr2 = tr1-tr4
               cr4 = tr1+tr4
               ci2 = ti1+ti4
               ci4 = ti1-ti4
               ch(i-1,k,2) = wa1(i-2)*cr2-wa1(i-1)*ci2
               ch(i  ,k,2) = wa1(i-2)*ci2+wa1(i-1)*cr2
               ch(i-1,k,3) = wa2(i-2)*cr3-wa2(i-1)*ci3
               ch(i  ,k,3) = wa2(i-2)*ci3+wa2(i-1)*cr3
               ch(i-1,k,4) = wa3(i-2)*cr4-wa3(i-1)*ci4
               ch(i  ,k,4) = wa3(i-2)*ci4+wa3(i-1)*cr4
            end do  ! i=3,ido,2
         end do  ! k=1,l1
         if (mod(ido,2) == 1) return
      end if

      do k=1,l1
         ti1 = cc(1,2,k) + cc(1,4,k)
         ti2 = cc(1,4,k) - cc(1,2,k)
         tr1 = cc(ido,1,k) - cc(ido,3,k)
         tr2 = cc(ido,1,k) + cc(ido,3,k)
         ch(ido,k,1) = tr2+tr2
         ch(ido,k,2) = sqrt2*(tr1-ti1)
         ch(ido,k,3) = ti2+ti2
         ch(ido,k,4) = -sqrt2*(tr1+ti1)
      end do
   end subroutine radb4
   !=======================================================================
   subroutine radb5 (ido,l1,cc,ch,wa1,wa2,wa3,wa4)
      integer :: ido,l1
      real(fft_kind) :: cc(ido,5,l1),ch(ido,l1,5),wa1(*),wa2(*),wa3(*),wa4(*)
      ! local variables
      integer :: i,ic,idp2,k
      real(fft_kind) :: ci2,ci3,ci4,ci5, cr2,cr3,cr4,cr5
      real(fft_kind) :: di2,di3,di4,di5, dr2,dr3,dr4,dr5
      real(fft_kind) :: ti2,ti3,ti4,ti5, tr2,tr3,tr4,tr5
      do k=1,l1
         ti5 = cc(1,3,k) + cc(1,3,k)
         ti4 = cc(1,5,k) + cc(1,5,k)
         tr2 = cc(ido,2,k) + cc(ido,2,k)
         tr3 = cc(ido,4,k) + cc(ido,4,k)
         ch(1,k,1) = cc(1,1,k) + tr2+tr3
         cr2 = cc(1,1,k) + tr11*tr2+tr12*tr3
         cr3 = cc(1,1,k) + tr12*tr2+tr11*tr3
         ci5 = ti11*ti5+ti12*ti4
         ci4 = ti12*ti5-ti11*ti4
         ch(1,k,2) = cr2-ci5
         ch(1,k,3) = cr3-ci4
         ch(1,k,4) = cr3+ci4
         ch(1,k,5) = cr2+ci5
      end do
      if (ido == 1) return
      idp2 = ido+2
      do k=1,l1
         do i=3,ido,2
            ic = idp2-i
            ti5 = cc(i  ,3,k) + cc(ic  ,2,k)
            ti2 = cc(i  ,3,k) - cc(ic  ,2,k)
            ti4 = cc(i  ,5,k) + cc(ic  ,4,k)
            ti3 = cc(i  ,5,k) - cc(ic  ,4,k)
            tr5 = cc(i-1,3,k) - cc(ic-1,2,k)
            tr2 = cc(i-1,3,k) + cc(ic-1,2,k)
            tr4 = cc(i-1,5,k) - cc(ic-1,4,k)
            tr3 = cc(i-1,5,k) + cc(ic-1,4,k)
            ch(i-1,k,1) = cc(i-1,1,k) + tr2+tr3
            ch(i  ,k,1) = cc(i  ,1,k) + ti2+ti3
            cr2 = cc(i-1,1,k) + tr11*tr2+tr12*tr3
            ci2 = cc(i  ,1,k) + tr11*ti2+tr12*ti3
            cr3 = cc(i-1,1,k) + tr12*tr2+tr11*tr3
            ci3 = cc(i  ,1,k) + tr12*ti2+tr11*ti3
            cr5 = ti11*tr5+ti12*tr4
            ci5 = ti11*ti5+ti12*ti4
            cr4 = ti12*tr5-ti11*tr4
            ci4 = ti12*ti5-ti11*ti4
            dr3 = cr3-ci4
            dr4 = cr3+ci4
            di3 = ci3+cr4
            di4 = ci3-cr4
            dr5 = cr2+ci5
            dr2 = cr2-ci5
            di5 = ci2-cr5
            di2 = ci2+cr5
            ch(i-1,k,2) = wa1(i-2)*dr2-wa1(i-1)*di2
            ch(i  ,k,2) = wa1(i-2)*di2+wa1(i-1)*dr2
            ch(i-1,k,3) = wa2(i-2)*dr3-wa2(i-1)*di3
            ch(i  ,k,3) = wa2(i-2)*di3+wa2(i-1)*dr3
            ch(i-1,k,4) = wa3(i-2)*dr4-wa3(i-1)*di4
            ch(i  ,k,4) = wa3(i-2)*di4+wa3(i-1)*dr4
            ch(i-1,k,5) = wa4(i-2)*dr5-wa4(i-1)*di5
            ch(i  ,k,5) = wa4(i-2)*di5+wa4(i-1)*dr5
         end do  ! i=3,ido,2
      end do  ! k=1,l1
      return
   end subroutine radb5
   !=======================================================================
   subroutine radbg (ido,ip,l1,idl1,cc,c1,c2,ch,ch2,wa)
      integer :: ido,ip,l1,idl1
      real(fft_kind) :: cc(ido, ip, l1),c1(ido, l1, ip),c2(idl1, ip)
      real(fft_kind) :: ch(ido, l1, ip),ch2(idl1, ip),wa(*)
      ! local variables
      real(fft_kind) :: ai1,ai2,ar1,ar1h,ar2,ar2h,arg,dc2,dcp,ds2,dsp
      integer :: i,ic,idij,idp2,ik,ipp2,ipph,is,j,j2,jc,k,l,lc,nbd
      arg = twopi/float(ip)
      dcp = cos(arg)
      dsp = sin(arg)
      idp2 = ido+2
      nbd = (ido-1)/2
      ipp2 = ip+2
      ipph = (ip+1)/2
      if (ido >= l1) then
         do k=1,l1
            do i=1,ido
               ch(i,k,1) = cc(i,1,k)
            end do
         end do
      else
         do i=1,ido
            do k=1,l1
               ch(i,k,1) = cc(i,1,k)
            end do
         end do
      end if
      do j=2,ipph
         jc = ipp2-j
         j2 = j+j
         do k=1,l1
            ch(1,k,j ) = cc(ido,j2-2,k) + cc(ido,j2-2,k)
            ch(1,k,jc) = cc(1,j2-1,k) + cc(1,j2-1,k)
         end do
      end do
      if (ido /= 1) then
         if (nbd >= l1) then
            do j=2,ipph
               jc = ipp2-j
               do k=1,l1
                  do i=3,ido,2
                     ic = idp2-i
                     ch(i-1,k,j ) = cc(i-1,2*j-1,k) + cc(ic-1,2*j-2,k)
                     ch(i-1,k,jc) = cc(i-1,2*j-1,k) - cc(ic-1,2*j-2,k)
                     ch(i  ,k,j ) = cc(i  ,2*j-1,k) - cc(ic  ,2*j-2,k)
                     ch(i  ,k,jc) = cc(i  ,2*j-1,k) + cc(ic  ,2*j-2,k)
                  end do
               end do
            end do
         else
            do j=2,ipph
               jc = ipp2-j
               do i=3,ido,2
                  ic = idp2-i
                  do k=1,l1
                     ch(i-1,k,j ) = cc(i-1,2*j-1,k) + cc(ic-1,2*j-2,k)
                     ch(i-1,k,jc) = cc(i-1,2*j-1,k) - cc(ic-1,2*j-2,k)
                     ch(i  ,k,j ) = cc(i  ,2*j-1,k) - cc(ic  ,2*j-2,k)
                     ch(i  ,k,jc) = cc(i  ,2*j-1,k) + cc(ic  ,2*j-2,k)
                  end do
               end do
            end do
         end if
      end if
      ar1 = 1.
      ai1 = 0.
      do l=2,ipph
         lc = ipp2-l
         ar1h = dcp*ar1-dsp*ai1
         ai1 = dcp*ai1+dsp*ar1
         ar1 = ar1h
         do ik=1,idl1
            c2(ik,l ) = ch2(ik,1) + ar1*ch2(ik,2)
            c2(ik,lc) = ai1*ch2(ik,ip)
         end do
         dc2 = ar1
         ds2 = ai1
         ar2 = ar1
         ai2 = ai1
         do j=3,ipph
            jc = ipp2-j
            ar2h = dc2*ar2-ds2*ai2
            ai2 = dc2*ai2+ds2*ar2
            ar2 = ar2h
            do ik=1,idl1
               c2(ik,l ) = c2(ik,l ) + ar2*ch2(ik,j )
               c2(ik,lc) = c2(ik,lc) + ai2*ch2(ik,jc)
            end do
         end do
      end do  ! l=2,ipph
      do j=2,ipph
         do ik=1,idl1
            ch2(ik,1) = ch2(ik,1) + ch2(ik,j )
         end do
      end do
      do j=2,ipph
         jc = ipp2-j
         do k=1,l1
            ch(1,k,j ) = c1(1,k,j ) - c1(1,k,jc)
            ch(1,k,jc) = c1(1,k,j ) + c1(1,k,jc)
         end do
      end do
      if (ido /= 1) then
         if (nbd >= l1) then
            do j=2,ipph
               jc = ipp2-j
               do k=1,l1
                  do i=3,ido,2
                     ch(i-1,k,j ) = c1(i-1,k,j ) - c1(i  ,k,jc)
                     ch(i-1,k,jc) = c1(i-1,k,j ) + c1(i  ,k,jc)
                     ch(i  ,k,j ) = c1(i  ,k,j ) + c1(i-1,k,jc)
                     ch(i  ,k,jc) = c1(i  ,k,j ) - c1(i-1,k,jc)
                  end do
               end do
            end do
         else
            do j=2,ipph
               jc = ipp2-j
               do i=3,ido,2
                  do k=1,l1
                     ch(i-1,k,j ) = c1(i-1,k,j ) - c1(i  ,k,jc)
                     ch(i-1,k,jc) = c1(i-1,k,j ) + c1(i  ,k,jc)
                     ch(i  ,k,j ) = c1(i  ,k,j ) + c1(i-1,k,jc)
                     ch(i  ,k,jc) = c1(i  ,k,j ) - c1(i-1,k,jc)
                  end do
               end do
            end do
         end if
      end if
      if (ido == 1) return
      do ik=1,idl1
         c2(ik,1) = ch2(ik,1)
      end do
      do j=2,ip
         do k=1,l1
            c1(1,k,j ) = ch(1,k,j )
         end do
      end do
      if (nbd <= l1) then
         is = -ido
         do j=2,ip
            is = is+ido
            idij = is
            do i=3,ido,2
               idij = idij+2
               do k=1,l1
                  c1(i-1,k,j) = wa(idij-1)*ch(i-1,k,j) - wa(idij)*ch(i  ,k,j)
                  c1(i  ,k,j) = wa(idij-1)*ch(i  ,k,j) + wa(idij)*ch(i-1,k,j)
               end do
            end do
         end do
      else
         is = -ido
         do j=2,ip
            is = is+ido
            do k=1,l1
               idij = is
               do i=3,ido,2
                  idij = idij+2
                  c1(i-1,k,j) = wa(idij-1)*ch(i-1,k,j) - wa(idij)*ch(i  ,k,j)
                  c1(i  ,k,j) = wa(idij-1)*ch(i  ,k,j) + wa(idij)*ch(i-1,k,j)
               end do
            end do
         end do
      end if
   end subroutine radbg
end subroutine fftpack_rfftb

