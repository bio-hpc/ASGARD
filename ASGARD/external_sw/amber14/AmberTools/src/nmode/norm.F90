
!************************************************************************
!                              AMBER                                   **
!                                                                      **
!                  Copyright (c) 1986, 1991, 1995                      **
!             Regents of the University of California                  **
!                       All Rights Reserved.                           **
!                                                                      **
!  This software provided pursuant to a license agreement containing   **
!  restrictions on its disclosure, duplication, and use. This software **
!  contains confidential and proprietary information, and may not be   **
!  extracted or distributed, in whole or in part, for any purpose      **
!  whatsoever, without the express written permission of the authors.  **
!  This notice, and the associated author list, must be attached to    **
!  all copies, or extracts, of this software. Any additional           **
!  restrictions set forth in the license agreement also apply to this  **
!  software.                                                           **
!************************************************************************


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine norm here]
subroutine norm (neig, nm, n, z, wr, wi, gamma)
   
   !     ----- normalize eigenvectors according to eq.(4.2b),
   !     ----- g. lamm and a. szabo, j. chem. phys. 85, 7334 (1986)
   
   implicit double precision (a-h,o-z)
   complex sum, elki, elli, lambda, minusi, zero
   dimension z(nm,n), wr(n), wi(n), gamma(*)
   
   n2 = n/2
   minusi = (0.0, -1.0)
   zero   = (0.0,  0.0)
   rzero  = 0.0
   
   i = 1
#ifdef ENDDO
   do while (i <= neig)
#else
   5 if (i > neig) goto 100
#endif
      lambda = cmplx ( wr(i), wi(i) )
      if (lambda == zero) then
         z(n2+1,i) = 0.0
         i = i + 1
      else if (wi(i) == 0.0) then
         if (wr(i) > 0.0) then
            sum = 0.0
            do 20 k = 1, n2
               elki = z(k,i)
               do 10 l = 1, n2
                  elli = z(l,i)
                  if (k > l) then
                     kl = k*(k-1)/2 + l
                  else
                     kl = l*(l-1)/2 + k
                  end if
                  sum = sum + elki * gamma(kl) * elli
               10 continue
               sum = sum + 2.0 * lambda * elki * elki
            20 continue
            z(n2+1,i) =   1.0
            sbyl = real (sum/lambda)
            const = 1.0 / sqrt(abs(sbyl))
            do 30 k = 1, n2
               z(k,i) = z(k,i) * const
            30 continue
            i = i + 1
         else
            sum = 0.0
            do 50 k = 1, n2
               elki = cmplx (rzero, z(k,i))
               do 40 l = 1, n2
                  elli = cmplx (rzero, z(l,i))
                  if (k > l) then
                     kl = k*(k-1)/2 + l
                  else
                     kl = l*(l-1)/2 + k
                  end if
                  sum = sum + elki * gamma(kl) * elli
               40 continue
               sum = sum + 2.0 * lambda * elki * elki
            50 continue
            sbyl = real (sum/lambda)
            if (sbyl == 0.0) then
               z(n2+1,i) = 0.0
               const = 0.0
            else if (sbyl < 0.0) then
               z(n2+1,i) =   1.0
               const = - 1.0 / sqrt(-sbyl)
            else
               z(n2+1,i) = - 1.0
               const = 1.0 / sqrt(sbyl)
            end if
            do 60 k = 1, n2
               z(k,i) = real ( z(k,i) * const )
            60 continue
            i = i + 1
         end if
      else
         sum = 0.0
         do 80 k = 1, n2
            elki = cmplx ( z(k,i), z(k,i+1) )
            do 70 l = 1, n2
               elli = cmplx ( z(l,i), z(l,i+1) )
               if (k > l) then
                  kl = k*(k-1)/2 + l
               else
                  kl = l*(l-1)/2 + k
               end if
               sum = sum + elki * gamma(kl) * elli
            70 continue
            sum = sum + 2.0 * (lambda * elki * elki)
         80 continue
         sum = sqrt ( sum / lambda )
         z(n2+1,i  ) =   2.0
         z(n2+1,i+1) = - 2.0
         do 90 k = 1, n2
            elki = cmplx ( z(k,i), z(k,i+1) )
            elki = elki / sum
            z(k,i) = real ( elki )
            z(k,i+1) = real ( minusi * elki )
         90 continue
         i = i + 2
      end if
#ifdef ENDDO
   end do
#else
      goto 5
      100 continue
#endif
   
   return
end subroutine norm 
