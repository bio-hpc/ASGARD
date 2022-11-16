
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine count here]
subroutine count (wi, n, kount, nkount, k, del)
   implicit double precision (a-h,o-z)
   dimension wi(n), kount(nkount)
   
   k = 1
   i = 13
   wref = del
   kount(1) = 0
   !     do while (i.le.n .and. k.le.nkount)
   do 20 mmm=1,999999
      if (i > n .or. k > nkount) goto 21
      wii = abs(wi(i))
      if (wii == 0.0) then
         continue
      else if (wii <= wref) then
         kount(k) = kount(k) + 1
      else
         !         do while (wii-wref .gt. del)
         do 10 nnn=1,999999
            if (wii-wref <= del) goto 11
            k = k + 1
            kount(k) = 0
            wref = wref + del
         10 continue
         11 k = k + 1
         kount(k) = 1
         wref = wref + del
      end if
      i = i + 1
   20 continue
   
   21 return
end subroutine count 
