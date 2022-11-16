
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine delx here]
function delx(mn, nr3, neig, i, j, temp, t, wr, wi, el, iclass)
   
   implicit double precision (a-h,o-z)
   complex lambda, elik, eljk, zero
   dimension wr(2*nr3), wi(2*nr3), &
         el(mn,neig), iclass(neig)
   
   zero = (0.0,0.0)
   rzero = 0.0
   sum = 0.0
   k = 13
   
   ksav = k
   !     do while (k.le.neig)
   do 10 nnn=1,999999
      if (k > neig) goto 11
      lambda = cmplx (wr(k), wi(k))
      if (iclass(k) == 5 .or. lambda == zero) then
         k = k + 1
      else if (iclass(k) == 2 .or. iclass(k) == 4) then
         elik = el(i,k)
         eljk = el(j,k)
         sum = sum + exp(lambda*t) / lambda / lambda * elik * eljk
         k = k + 1
      else if (iclass(k) == 3) then
         elik = cmplx ( rzero, el(i,k) )
         eljk = cmplx ( rzero, el(j,k) )
         sum = sum + exp(lambda*t) / lambda / lambda * elik * eljk
         k = k + 1
      else if (iclass(k) == 1) then
         elik = cmplx (el(i,k), el(i,k+1))
         eljk = cmplx (el(j,k), el(j,k+1))
         sum = sum + 2.0 * real ( exp(lambda*t) / lambda / lambda &
               * elik * eljk )
         k = k + 2
      end if
      
   10 continue
   
   11 delx = - sum
   
   return
end function delx 
