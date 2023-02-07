
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
!+ [Enter a one-line description of subroutine dspfa here]
subroutine dspfa(ap,n, kpvt,     info)
   !       call dspfa(     h,n3,ix(mkpvt),info)

   integer n,kpvt(1),info
   double precision ap(1)

   !     DSPFA FACTORS A DOUBLE PRECISION SYMMETRIC MATRIX STORED IN
   !     PACKED FORM BY ELIMINATION WITH SYMMETRIC PIVOTING.

   !     TO SOLVE  A*X = B , FOLLOW DSPFA BY DSPSL.
   !     TO COMPUTE  INVERSE(A)*C , FOLLOW DSPFA BY DSPSL.
   !     TO COMPUTE  DETERMINANT(A) , FOLLOW DSPFA BY DSPDI.
   !     TO COMPUTE  INERTIA(A) , FOLLOW DSPFA BY DSPDI.
   !     TO COMPUTE  INVERSE(A) , FOLLOW DSPFA BY DSPDI.

   !     ON ENTRY

   !        AP      DOUBLE PRECISION (N*(N+1)/2)
   !                THE PACKED FORM OF A SYMMETRIC MATRIX  A .  THE
   !                COLUMNS OF THE UPPER TRIANGLE ARE STORED SEQUENTIALLY
   !                IN A ONE-DIMENSIONAL ARRAY OF LENGTH  N*(N+1)/2 .
   !                SEE COMMENTS BELOW FOR DETAILS.

   !        N       INTEGER
   !                THE ORDER OF THE MATRIX  A .

   !     OUTPUT

   !        AP      A BLOCK DIAGONAL MATRIX AND THE MULTIPLIERS WHICH
   !                WERE USED TO OBTAIN IT STORED IN PACKED FORM.
   !                THE FACTORIZATION CAN BE WRITTEN  A = U*D*TRANS(U)
   !                WHERE  U  IS A PRODUCT OF PERMUTATION AND UNIT
   !                UPPER TRIANGULAR MATRICES , TRANS(U) IS THE
   !                TRANSPOSE OF  U , AND  D  IS BLOCK DIAGONAL
   !                WITH 1 BY 1 AND 2 BY 2 BLOCKS.

   !        KPVT    INTEGER(N)
   !                AN INTEGER VECTOR OF PIVOT INDICES.

   !        INFO    INTEGER
   !                = 0  NORMAL VALUE.
   !                = K  IF THE K-TH PIVOT BLOCK IS SINGULAR. THIS IS
   !                     NOT AN ERROR CONDITION FOR THIS SUBROUTINE,
   !                     BUT IT DOES INDICATE THAT DSPSL OR DSPDI MAY
   !                     DIVIDE BY ZERO IF CALLED.

   !     PACKED STORAGE

   !          THE FOLLOWING PROGRAM SEGMENT WILL PACK THE UPPER
   !          TRIANGLE OF A SYMMETRIC MATRIX.

   !                K = 0
   !                DO 20 J = 1, N
   !                   DO 10 I = 1, J
   !                      K = K + 1
   !                      AP(K)  = A(I,J)
   !             10    CONTINUE
   !             20 CONTINUE

   !     LINPACK. THIS VERSION DATED 08/14/78 .
   !     JAMES BUNCH, UNIV. CALIF. SAN DIEGO, ARGONNE NAT. LAB.

   !     SUBROUTINES AND FUNCTIONS

   !     BLAS DAXPY,DSWAP,IDAMAX
   !     FORTRAN DABS,DMAX1,DSQRT

   !     INTERNAL VARIABLES

   double precision ak,akm1,bk,bkm1,denom,mulk,mulkm1,t
   double precision absakk,alpha,colmax,rowmax
   integer idamax,ij,ijj,ik,ikm1,im,imax,imaxp1,imim,imj,imk
   integer j,jj,jk,jkm1,jmax,jmim,k,kk,km1,km1k,km1km1,km2,kstep
   logical swap


   !     INITIALIZE

   !     ALPHA IS USED IN CHOOSING PIVOT BLOCK SIZE.

   alpha = (1.0d0 + dsqrt(17.0d0))/8.0d0

   info = 0

   !     MAIN LOOP ON K, WHICH GOES FROM N TO 1.

   k = n
   ik = (n*(n - 1))/2
   10 continue

   !        LEAVE THE LOOP IF K=0 OR K=1.

   !     ...EXIT
   if (k == 0) goto 200
   if (k > 1) goto 20
   kpvt(1) = 1
   if (ap(1) == 0.0d0) info = 1
   !     ......EXIT
   goto 200
   20 continue

   !        THIS SECTION OF CODE DETERMINES THE KIND OF
   !        ELIMINATION TO BE PERFORMED.  WHEN IT IS COMPLETED,
   !        KSTEP WILL BE SET TO THE SIZE OF THE PIVOT BLOCK, AND
   !        SWAP WILL BE SET TO .TRUE. IF AN INTERCHANGE IS
   !        REQUIRED.

   km1 = k - 1
   kk = ik + k
   absakk = dabs(ap(kk))

   !        DETERMINE THE LARGEST OFF-DIAGONAL ELEMENT IN
   !        COLUMN K.

   imax = idamax(k-1,ap(ik+1),1)
   imk = ik + imax
   colmax = dabs(ap(imk))
   if (absakk < alpha*colmax) goto 30
   kstep = 1
   swap = .false.
   goto 90
   30 continue

   !           DETERMINE THE LARGEST OFF-DIAGONAL ELEMENT IN
   !           ROW IMAX.

   rowmax = 0.0d0
   imaxp1 = imax + 1
   im = imax*(imax - 1)/2
   imj = im + 2*imax
   do 40 j = imaxp1, k
      rowmax = dmax1(rowmax,dabs(ap(imj)))
      imj = imj + j
   40 continue
   if (imax == 1) goto 50
   jmax = idamax(imax-1,ap(im+1),1)
   jmim = jmax + im
   rowmax = dmax1(rowmax,dabs(ap(jmim)))
   50 continue
   imim = imax + im
   if (dabs(ap(imim)) < alpha*rowmax) goto 60
   kstep = 1
   swap = .true.
   goto 80
   60 continue
   if (absakk < alpha*colmax*(colmax/rowmax)) goto 70
   kstep = 1
   swap = .false.
   goto 80
   70 continue
   kstep = 2
   swap = imax /= km1
   80 continue
   90 continue
   if (dmax1(absakk,colmax) /= 0.0d0) goto 100

   !           COLUMN K IS ZERO.  SET INFO AND ITERATE THE LOOP.

   kpvt(k) = k
   info = k
   goto 190
   100 continue
   if (kstep == 2) goto 140

   !           1 X 1 PIVOT BLOCK.

   if (.not.swap) goto 120

   !              PERFORM AN INTERCHANGE.

   call dswap(imax,ap(im+1),1,ap(ik+1),1)
   imj = ik + imax
   do 110 jj = imax, k
      j = k + imax - jj
      jk = ik + j
      t = ap(jk)
      ap(jk) = ap(imj)
      ap(imj) = t
      imj = imj - (j - 1)
   110 continue
   120 continue

   !           PERFORM THE ELIMINATION.

   ij = ik - (k - 1)
   do 130 jj = 1, km1
      j = k - jj
      jk = ik + j
      mulk = -ap(jk)/ap(kk)
      t = mulk
      call daxpy(j,t,ap(ik+1),1,ap(ij+1),1)
      ijj = ij + j
      ap(jk) = mulk
      ij = ij - (j - 1)
   130 continue

   !           SET THE PIVOT ARRAY.

   kpvt(k) = k
   if (swap) kpvt(k) = imax
   goto 190
   140 continue

   !           2 X 2 PIVOT BLOCK.

   km1k = ik + k - 1
   ikm1 = ik - (k - 1)
   if (.not.swap) goto 160

   !              PERFORM AN INTERCHANGE.

   call dswap(imax,ap(im+1),1,ap(ikm1+1),1)
   imj = ikm1 + imax
   do 150 jj = imax, km1
      j = km1 + imax - jj
      jkm1 = ikm1 + j
      t = ap(jkm1)
      ap(jkm1) = ap(imj)
      ap(imj) = t
      imj = imj - (j - 1)
   150 continue
   t = ap(km1k)
   ap(km1k) = ap(imk)
   ap(imk) = t
   160 continue

   !           PERFORM THE ELIMINATION.

   km2 = k - 2
   if (km2 == 0) goto 180
   ak = ap(kk)/ap(km1k)
   km1km1 = ikm1 + k - 1
   akm1 = ap(km1km1)/ap(km1k)
   denom = 1.0d0 - ak*akm1
   ij = ik - (k - 1) - (k - 2)
   do 170 jj = 1, km2
      j = km1 - jj
      jk = ik + j
      bk = ap(jk)/ap(km1k)
      jkm1 = ikm1 + j
      bkm1 = ap(jkm1)/ap(km1k)
      mulk = (akm1*bk - bkm1)/denom
      mulkm1 = (ak*bkm1 - bk)/denom
      t = mulk
      call daxpy(j,t,ap(ik+1),1,ap(ij+1),1)
      t = mulkm1
      call daxpy(j,t,ap(ikm1+1),1,ap(ij+1),1)
      ap(jk) = mulk
      ap(jkm1) = mulkm1
      ijj = ij + j
      ij = ij - (j - 1)
   170 continue
   180 continue

   !           SET THE PIVOT ARRAY.

   kpvt(k) = 1 - k
   if (swap) kpvt(k) = -imax
   kpvt(k-1) = kpvt(k)
   190 continue
   ik = ik - (k - 1)
   if (kstep == 2) ik = ik - (k - 2)
   k = k - kstep
   goto 10
   200 continue
   return
end subroutine dspfa 
