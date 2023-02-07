! #+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#
! #|  Compute the inverse of the B matrix                            |#
! #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

#include "../include/dprec.fh"

   subroutine bmat_inv_debug ( bmat, binv, gmat, ginv, natm, ncoord, zero )

   implicit none

   integer, intent(in) :: natm, ncoord 
   _REAL_ , intent(in) :: bmat(ncoord,natm*3), binv(natm*3,ncoord)
   _REAL_ , intent(in) :: gmat(ncoord,ncoord), ginv(ncoord,ncoord)
   _REAL_ , intent(in) :: zero
 
   !  ..........................................................................

   intrinsic :: abs, matmul

   _REAL_  :: BinvBB(natm*3,ncoord), BBinvB(ncoord,natm*3) &
            , BBBinv(natm*3,ncoord), binv_dev(ncoord,natm*3)

   _REAL_  :: GinvGG(ncoord,ncoord), GGinvG(ncoord,ncoord) &
            , GGGinv(ncoord,ncoord), ginv_dev(ncoord,ncoord)

   character(512) :: err_string

!  +---------------------------------------------------------------------------+
!  |  Check inverse of the G matrix                                            |
!  +---------------------------------------------------------------------------+

!  ,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
!  |  [G] = ( [G^-1] [G] ) [G]                                                 |
!  '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

   GinvGG(:,:) = matmul( matmul( ginv, gmat ), gmat)
   ginv_dev(:,:) = gmat - GinvGG

   err_string = "%%%%%%  CHECKING ( [G^-1] [G] ) [G] == [G]  %%%%%%%%%%%%%%%%"
   call inv_error ( gmat, GinvGG, ginv_dev, ncoord, ncoord, zero, err_string )

!  ,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
!  |  [G] = ( [G] [G^-1] ) [G]                                                 |
!  '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

   GGinvG(:,:) = matmul( matmul( gmat, ginv ), gmat)
   ginv_dev(:,:) = gmat - GGinvG

   err_string = "%%%%%%  CHECKING ( [G] [G^-1] ) [G] == [G]  %%%%%%%%%%%%%%%%"
   call inv_error ( gmat, GGinvG, ginv_dev, ncoord, ncoord, zero, err_string )

!  ,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
!  |  [G] = [G] [G^-1] [G]                                                     |
!  '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

   GGinvG(:,:) = matmul( gmat, matmul( ginv, gmat ) )
   ginv_dev(:,:) = gmat - GGinvG

   err_string = "%%%%%%  CHECKING [G] [G^-1] [G] == [G]  %%%%%%%%%%%%%%%%"
   call inv_error ( gmat, GGinvG, ginv_dev, ncoord, ncoord, zero, err_string )

!  ,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
!  |  [G] = [G] [G] [G^-1]                                                     |
!  '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

   GGGinv(:,:) = matmul( gmat, matmul( gmat, ginv ) )
   ginv_dev(:,:) = gmat - GGGinv

   err_string = "%%%%%%  CHECKING [G] [G] [G^-1] == [G]  %%%%%%%%%%%%%%%%"
   call inv_error ( gmat, GGGinv, ginv_dev, ncoord, ncoord, zero, err_string )


   write(6,'(A)') '%%%%%%  DONE CHECKING [G^-1]   %%%%%%%%%%%%%%%%'

!  +---------------------------------------------------------------------------+
!  |  Check inverse of the B matrix                                            |
!  +---------------------------------------------------------------------------+

!  ,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
!  |  [B] = ( [B^-1] [B] ) [B]^t                                               |
!  '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

   BinvBB(:,:) = matmul( matmul( binv, bmat ), transpose(bmat) )
   binv_dev(:,:) = bmat - transpose(BinvBB)

   err_string = "%%%%%%  CHECKING ( [B^-1] [B] ) [B]^t == [B]  %%%%%%%%%%%%%%%%"
   call inv_error ( bmat, transpose(BinvBB), binv_dev, ncoord, natm*3 &
                  , zero, err_string )

!  ,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
!  |  [B] = ( [B] [B^-1] ) [B]                                                 |
!  '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

   BBinvB(:,:) = matmul( matmul( bmat, binv ), bmat)
   binv_dev(:,:) = bmat - BBinvB

   err_string = "%%%%%%  CHECKING ( [B] [B^-1] ) [B] == [B]  %%%%%%%%%%%%%%%%"
   call inv_error ( bmat, BBinvB, binv_dev, ncoord, natm*3, zero, err_string )

!  ,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
!  |  [B] = [B] [B^-1] [B]                                                     |
!  '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

   BBinvB(:,:) = matmul( bmat, matmul( binv, bmat ) )
   binv_dev(:,:) = bmat - BBinvB

   err_string = "%%%%%%  CHECKING [B] [B^-1] [B] == [B]  %%%%%%%%%%%%%%%%"
   call inv_error ( bmat, BBinvB, binv_dev, ncoord, natm*3, zero, err_string )

!  ,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
!  |  [B] = [B]^t [B] [B^-1]                                                   |
!  '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

   BBBinv(:,:) = matmul( transpose(bmat), matmul( bmat, binv ) )
   binv_dev(:,:) = bmat - transpose(BBBinv)

   err_string = "%%%%%%  CHECKING [B] [B] [B^-1] == [B]  %%%%%%%%%%%%%%%%"
   call inv_error ( bmat, transpose(BBBinv), binv_dev, ncoord, natm*3 &
                  , zero, err_string )




   write(6,'(A)') '%%%%%%  DONE CHECKING [B^-1]   %%%%%%%%%%%%%%%%'


!  x+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++x
!  x+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++x

   end subroutine bmat_inv_debug


   subroutine inv_error ( ref, test, mat_dev, mdim, ndim, zero, err_string )

   implicit none

   integer, intent(in) :: mdim, ndim
   _REAL_ , intent(in) :: ref(mdim,ndim), test(mdim,ndim) &
                        , mat_dev(mdim,ndim), zero

   character(512), intent(in) :: err_string
   !  ..........................................................................

   integer :: m, n
   logical :: lerror


   write(6,'(A)')
   write(6,'(A)') trim( adjustl( err_string ) )

   lerror = .false.
   do n = 1, ndim
      do m = 1, mdim
         if( abs( mat_dev(m,n) ) > zero ) then
            write(6,'(2I8,3(1PE16.8))') m, n, ref(m,n), test(m,n) &
                                      , mat_dev(m,n)
            lerror = .true.
         endif
      enddo
   enddo

   if( .not. lerror ) write(6,'(A)') '------  NO ERRORS DETECTED          ----------------'

!  x+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++x
!  x+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++x

   end subroutine inv_error

