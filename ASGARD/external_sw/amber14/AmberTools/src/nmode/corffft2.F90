
!---------------------------------------------------------------------


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine corffft2 here]
subroutine corffft2(ndata,data,table)
   
   !     ----- Calculates autocorrelation function
   !           using the Wiener-Khinchin-Theorem
   !           (s. Comp. Sim. of Liquids, p. 188)
   
   !     ----- ndata is the length of the array data
   !     ----- data is a real array of complex numbers:
   !           data(1)=real(1), data(2)=img(1), ...
   !     ----- it is recommended that the discrete data is appended
   !           with as much zeros to avoid spurious correlations
   !     ----- in addition, data MUST have the dimension of power of 2
   !           (pad the "real" data (plus zeros) with additional
   !            zeros up to the next power of 2)
   !     ----- the result is not yet normalized by (no_of_discrete_data - t)**-1 (!)
   
   !     ----- HG 11/04/2002
   
   implicit double precision (a-h,o-z)
   dimension data(*), table(*)
   
   ntmp = ndata
   do i=1,9999999
      ntmp = ntmp / 2
      if(ntmp <= 1) goto 10
   end do
   10 if(ntmp /= 1) then
      write(0,*) 'ndata is not power of 2: ',ndata
      call mexit(6,1)
   end if
   
   !     --- FFT data
   
   call cfftf(ndata/2, data, table)
   
   !     --- Calc square modulus
   
   do i=1,ndata,2
      data(i) = data(i)*data(i) + data(i+1)*data(i+1)
      data(i+1) = 0.0d0
   end do
   
   !     --- Inverse FFT
   
   call cfftb(ndata/2, data, table)
   
   !     --- Normalize with ndata/2 (since not done in inverse FFT routine)
   
   ddata = 1.0d0 / dble(ndata/2)
   do i=1,ndata
      data(i) = data(i) * ddata
   end do
   
   return
end subroutine corffft2 

!---------------------------------------------------------------------

