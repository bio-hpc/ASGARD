
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
!+ [Enter a one-line description of subroutine dppsl here]
subroutine dppsl(ap,n,b)
   integer n
   double precision ap(1),b(1)
   
   !     DPPSL SOLVES THE DOUBLE PRECISION SYMMETRIC POSITIVE DEFINITE
   !     SYSTEM A * X = B
   !     USING THE FACTORS COMPUTED BY DPPCO OR DPPFA.
   
   !     ON ENTRY
   
   !        AP      DOUBLE PRECISION (N*(N+1)/2)
   !                THE OUTPUT FROM DPPCO OR DPPFA.
   
   !        N       INTEGER
   !                THE ORDER OF THE MATRIX  A .
   
   !        B       DOUBLE PRECISION(N)
   !                THE RIGHT HAND SIDE VECTOR.
   
   !     ON RETURN
   
   !        B       THE SOLUTION VECTOR  X .
   
   !     ERROR CONDITION
   
   !        A DIVISION BY ZERO WILL OCCUR IF THE INPUT FACTOR CONTAINS
   !        A ZERO ON THE DIAGONAL.  TECHNICALLY THIS INDICATES
   !        SINGULARITY BUT IT IS USUALLY CAUSED BY IMPROPER SUBROUTINE
   !        ARGUMENTS.  IT WILL NOT OCCUR IF THE SUBROUTINES ARE CALLED
   !        CORRECTLY AND  INFO .EQ. 0 .
   
   !     TO COMPUTE  INVERSE(A) * C  WHERE  C  IS A MATRIX
   !     WITH  P  COLUMNS
   !           CALL DPPCO(AP,N,RCOND,Z,INFO)
   !           IF (RCOND IS TOO SMALL .OR. INFO .NE. 0) GO TO ...
   !           DO 10 J = 1, P
   !              CALL DPPSL(AP,N,C(1,J))
   !        10 CONTINUE
   
   !     LINPACK.  THIS VERSION DATED 08/14/78 .
   !     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB.
   
   !     SUBROUTINES AND FUNCTIONS
   
   !     BLAS DAXPY,DDOT
   
   !     INTERNAL VARIABLES
   
   double precision ddot,t
   integer k,kb,kk
   
   kk = 0
   do 10 k = 1, n
      t = ddot(k-1,ap(kk+1),1,b(1),1)
      kk = kk + k
      b(k) = (b(k) - t)/ap(kk)
   10 continue
   do 20 kb = 1, n
      k = n + 1 - kb
      b(k) = b(k)/ap(kk)
      kk = kk - k
      t = -b(k)
      call daxpy(k-1,t,ap(kk+1),1,b(1),1)
   20 continue
   return
end subroutine dppsl 
