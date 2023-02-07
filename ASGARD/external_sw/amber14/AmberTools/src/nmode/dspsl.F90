
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
!+ [Enter a one-line description of subroutine dspsl here]
subroutine dspsl(ap,n,kpvt,b)
   integer n,kpvt(1)
   double precision ap(1),b(1)
   
   !     DSISL SOLVES THE DOUBLE PRECISION SYMMETRIC SYSTEM
   !     A * X = B
   !     USING THE FACTORS COMPUTED BY DSPFA.
   
   !     ON ENTRY
   
   !        AP      DOUBLE PRECISION(N*(N+1)/2)
   !                THE OUTPUT FROM DSPFA.
   
   !        N       INTEGER
   !                THE ORDER OF THE MATRIX  A .
   
   !        KPVT    INTEGER(N)
   !                THE PIVOT VECTOR FROM DSPFA.
   
   !        B       DOUBLE PRECISION(N)
   !                THE RIGHT HAND SIDE VECTOR.
   
   !     ON RETURN
   
   !        B       THE SOLUTION VECTOR  X .
   
   !     ERROR CONDITION
   
   !        A DIVISION BY ZERO MAY OCCUR IF  DSPCO  HAS SET RCOND .EQ. 0.0
   !        OR  DSPFA  HAS SET INFO .NE. 0  .
   
   !     TO COMPUTE  INVERSE(A) * C  WHERE  C  IS A MATRIX
   !     WITH  P  COLUMNS
   !           CALL DSPFA(AP,N,KPVT,INFO)
   !           IF (INFO .NE. 0) GO TO ...
   !           DO 10 J = 1, P
   !              CALL DSPSL(AP,N,KPVT,C(1,J))
   !        10 CONTINUE
   
   !     LINPACK. THIS VERSION DATED 08/14/78 .
   !     JAMES BUNCH, UNIV. CALIF. SAN DIEGO, ARGONNE NAT. LAB.
   
   !     SUBROUTINES AND FUNCTIONS
   
   !     BLAS DAXPY,DDOT
   !     FORTRAN IABS
   
   !     INTERNAL VARIABLES.
   
   double precision ak,akm1,bk,bkm1,ddot,denom,temp
   integer ik,ikm1,ikp1,k,kk,km1k,km1km1,kp
   
   !     LOOP BACKWARD APPLYING THE TRANSFORMATIONS AND
   !     D INVERSE TO B.
   
   k = n
   ik = (n*(n - 1))/2
   10 if (k == 0) goto 80
   kk = ik + k
   if (kpvt(k) < 0) goto 40
   
   !           1 X 1 PIVOT BLOCK.
   
   if (k == 1) goto 30
   kp = kpvt(k)
   if (kp == k) goto 20
   
   !                 INTERCHANGE.
   
   temp = b(k)
   b(k) = b(kp)
   b(kp) = temp
   20 continue
   
   !              APPLY THE TRANSFORMATION.
   
   call daxpy(k-1,b(k),ap(ik+1),1,b(1),1)
   30 continue
   
   !           APPLY D INVERSE.
   
   b(k) = b(k)/ap(kk)
   k = k - 1
   ik = ik - k
   goto 70
   40 continue
   
   !           2 X 2 PIVOT BLOCK.
   
   ikm1 = ik - (k - 1)
   if (k == 2) goto 60
   kp = iabs(kpvt(k))
   if (kp == k - 1) goto 50
   
   !                 INTERCHANGE.
   
   temp = b(k-1)
   b(k-1) = b(kp)
   b(kp) = temp
   50 continue
   
   !              APPLY THE TRANSFORMATION.
   
   call daxpy(k-2,b(k),ap(ik+1),1,b(1),1)
   call daxpy(k-2,b(k-1),ap(ikm1+1),1,b(1),1)
   60 continue
   
   !           APPLY D INVERSE.
   
   km1k = ik + k - 1
   kk = ik + k
   ak = ap(kk)/ap(km1k)
   km1km1 = ikm1 + k - 1
   akm1 = ap(km1km1)/ap(km1k)
   bk = b(k)/ap(km1k)
   bkm1 = b(k-1)/ap(km1k)
   denom = ak*akm1 - 1.0d0
   b(k) = (akm1*bk - bkm1)/denom
   b(k-1) = (ak*bkm1 - bk)/denom
   k = k - 2
   ik = ik - (k + 1) - k
   70 continue
   goto 10
   80 continue
   
   !     LOOP FORWARD APPLYING THE TRANSFORMATIONS.
   
   k = 1
   ik = 0
   90 if (k > n) goto 160
   if (kpvt(k) < 0) goto 120
   
   !           1 X 1 PIVOT BLOCK.
   
   if (k == 1) goto 110
   
   !              APPLY THE TRANSFORMATION.
   
   b(k) = b(k) + ddot(k-1,ap(ik+1),1,b(1),1)
   kp = kpvt(k)
   if (kp == k) goto 100
   
   !                 INTERCHANGE.
   
   temp = b(k)
   b(k) = b(kp)
   b(kp) = temp
   100 continue
   110 continue
   ik = ik + k
   k = k + 1
   goto 150
   120 continue
   
   !           2 X 2 PIVOT BLOCK.
   
   if (k == 1) goto 140
   
   !              APPLY THE TRANSFORMATION.
   
   b(k) = b(k) + ddot(k-1,ap(ik+1),1,b(1),1)
   ikp1 = ik + k
   b(k+1) = b(k+1) + ddot(k-1,ap(ikp1+1),1,b(1),1)
   kp = iabs(kpvt(k))
   if (kp == k) goto 130
   
   !                 INTERCHANGE.
   
   temp = b(k)
   b(k) = b(kp)
   b(kp) = temp
   130 continue
   140 continue
   ik = ik + k + k + 1
   k = k + 2
   150 continue
   goto 90
   160 continue
   return
end subroutine dspsl 
