
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
!+ [Enter a one-line description of subroutine dppfa here]
subroutine dppfa(ap,n,info)
   integer n,info
   double precision ap(1)
   
   !     DPPFA FACTORS A DOUBLE PRECISION SYMMETRIC POSITIVE DEFINITE
   !     MATRIX STORED IN PACKED FORM.
   
   !     DPPFA IS USUALLY CALLED BY DPPCO, BUT IT CAN BE CALLED
   !     DIRECTLY WITH A SAVING IN TIME IF  RCOND  IS NOT NEEDED.
   !     (TIME FOR DPPCO) = (1 + 18/N)*(TIME FOR DPPFA) .
   
   !     ON ENTRY
   
   !        AP      DOUBLE PRECISION (N*(N+1)/2)
   !                THE PACKED FORM OF A SYMMETRIC MATRIX  A .  THE
   !                COLUMNS OF THE UPPER TRIANGLE ARE STORED SEQUENTIALLY
   !                IN A ONE-DIMENSIONAL ARRAY OF LENGTH  N*(N+1)/2 .
   !                SEE COMMENTS BELOW FOR DETAILS.
   
   !        N       INTEGER
   !                THE ORDER OF THE MATRIX  A .
   
   !     ON RETURN
   
   !        AP      AN UPPER TRIANGULAR MATRIX  R , STORED IN PACKED
   !                FORM, SO THAT  A = TRANS(R)*R .
   
   !        INFO    INTEGER
   !                = 0  FOR NORMAL RETURN.
   !                = K  IF THE LEADING MINOR OF ORDER  K  IS NOT
   !                     POSITIVE DEFINITE.
   
   
   !     PACKED STORAGE
   
   !          THE FOLLOWING PROGRAM SEGMENT WILL PACK THE UPPER
   !          TRIANGLE OF A SYMMETRIC MATRIX.
   
   !                K = 0
   !                DO 20 J = 1, N
   !                   DO 10 I = 1, J
   !                      K = K + 1
   !                      AP(K) = A(I,J)
   !             10    CONTINUE
   !             20 CONTINUE
   
   !     LINPACK.  THIS VERSION DATED 08/14/78 .
   !     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB.
   
   !     SUBROUTINES AND FUNCTIONS
   
   !     BLAS DDOT
   !     FORTRAN DSQRT
   
   !     INTERNAL VARIABLES
   
   double precision ddot,t
   double precision s
   integer j,jj,jm1,k,kj,kk
   !     BEGIN BLOCK WITH ...EXITS TO 40
   
   
   jj = 0
   do 30 j = 1, n
      info = j
      s = 0.0d0
      jm1 = j - 1
      kj = jj
      kk = 0
      if (jm1 < 1) goto 20
      do 10 k = 1, jm1
         kj = kj + 1
         t = ap(kj) - ddot(k-1,ap(kk+1),1,ap(jj+1),1)
         kk = kk + k
         t = t/ap(kk)
         ap(kj) = t
         s = s + t*t
      10 continue
      20 continue
      jj = jj + j
      s = ap(jj) - s
      !     ......EXIT
      if (s <= 0.0d0) goto 40
      ap(jj) = dsqrt(s)
   30 continue
   info = 0
   40 continue
   return
end subroutine dppfa 
