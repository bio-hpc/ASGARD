
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine dsort here]
subroutine dsort (n, dbl)
   implicit double precision (a-h,o-z)
   
   dimension dbl(n)
   
   !     ----- sort keeping track of which element goes where
   
   do 30 i = 1, n-1
      
      !       ----- get pointer to minimum value of remaining entries
      
      ip = i
      do 20 j = i+1, n
         if (dbl(j) > dbl(ip)) ip = j
      20 continue
      
      !       ----- interchange ip-th and i-th dbl
      
      dt = dbl(i)
      dbl(i) = dbl(ip)
      dbl(ip) = dt
      
   30 continue
   
   return
end subroutine dsort 
