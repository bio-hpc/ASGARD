
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine corf here]
subroutine corf(nat3, nvect, freq, vect, x, kup, lup, &
      bose, ibeg, iend, vtemp)
   
   !     ----- to calculate dipole-dipole correlation functions
   
   implicit double precision (a-h,o-z)
   double precision kt
   double precision lkin,llin,lkjn,lljn
   logical bose
   dimension freq(nat3), vect(nat3,nvect), x(nat3), e(3)
   dimension amass(19)
   dimension del(3,3), del0(3,3)
   
   !     functions
   k(i) = 3*(kup-1) + i
   l(i) = 3*(lup-1) + i
   !     end functions
   
   kt=vtemp*0.002
   consq = 0.71942/vtemp
   !                  ----CONSQ = hc/2kT in cm, with T=300K
   !                       (use for quantum, Bose statistics)
   p2 = 1.d0
   
#if 0
                                                                                
         !       --- check orthogonality: (masses here for NMA)
                                                                                
         amass(1) = 1.008
         amass(2) = 12.01
         amass(3) = 1.008
         amass(4) = 1.008
         amass(5) = 12.01
         amass(6) = 16.00
         amass(7) = 14.01
         amass(8) = 1.008
         amass(9) = 12.01
         amass(10) = 1.008
         amass(11) = 1.008
         amass(12) = 1.008
         do ic=1,nat3
           do jc=ic,nat3
             dot = 0.0
             n = 0
             do i=1,nat3
                if( mod(i,3) == 1 ) n = n + 1
                dot = dot + vect(i,ic)*vect(i,jc)*amass(n)
             end do
             write(6,*) 'orthog check: ', ic, jc, dot
           end do
         end do
#endif

   !    ---e(i) is the unit vector along the l-k bond:
   
   req = 0.0
   do 20 i = 1, 3
      e(i) = x(k(i)) - x(l(i))
      req = req + e(i)*e(i)
   20 continue
   req = sqrt(req)
   do i = 1, 3
      e(i) = e(i) / req
   end do
   
   !     ----- loop over the desired modes:
   
   write (6,'(2i5)') kup,lup
   do 35 n=ibeg,iend
      if (freq(n) < 0.5) goto 35
      fre = freq(n)*freq(n)/11791.79
      if (bose) then
         argq = freq(n)*consq
         qcorr = argq/tanh(argq)
      else
         qcorr = 1.0
      end if
      
      !     ----- calculate the correlation matrix for delta
      !        as in Eq. 7.16 of lamm and szabo j. chem. phys. 85, 7334 (1986)
      !        Note that the rhs of this eq. should be multiplied by kT
      
      do i = 1, 3
         ki = k(i)
         li = l(i)
         do j = 1, 3
            kj = k(j)
            lj = l(j)
            del(i,j) = 0.d0
            lkin = vect(ki,n)
            llin = vect(li,n)
            lkjn = vect(kj,n)
            lljn = vect(lj,n)
            delijc = (qcorr/fre) * (lkin - llin) * (lkjn - lljn)
            del(i,j) = kt*delijc
         end do
      end do
      
#if 0
      
      !         ----- correlation in length, Eq. (10.2) of lamm and szabo:
      
      rtr0 = 0.0
      do i = 1, 3
         do j = 1,3
            rtr0 = rtr0 + e(i)*e(j)*del(i,j)
         end do
      end do
      
#endif
      
      !         ----- librational correlation function, using eq. 7.12
      !               of lamm and szabo (without beta on the lhs):
      
      sum = 0.0
      do i = 1, 3
         sum = sum - del(i,i)
         do j = 1,3
            sum = sum + del(i,j)*e(i)*e(j)
         end do
      end do
      sum = (3.0/req**2) * sum
      p2 = p2 + sum
      write (6,'(10x,i5,f10.2,2f10.5)') n,freq(n),sum,p2
   35 continue
   write(6,*) '-----------------------------------------------------'
   
   return
end subroutine corf 
