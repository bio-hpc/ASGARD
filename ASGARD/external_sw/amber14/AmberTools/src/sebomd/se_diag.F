!------------------------------------------------------------------------
      subroutine se_diag(ndim,a,nevec1,tolera,v,eval1,idegen1,evec1,
     &     ierror)
!
!     driver routine for diagonalization of the real, symmetric,
!     matrix a.
!
!     variables required:
!     ------------------
!
!     ndim = order of the matrix a (i.e., a is ndim by ndim);
!
!     a = real symmetric matrix to be diagonalized.  only the lower
!         half of a need be filled when calling.  a is destroyed by
!         this routine.
!
!     nevec1 = the number of eigenvectors required;
!
!     tolera = tolerance factor used in the qr iteration to determine
!              when off-diagonal entries are essentially zero.  (default
!              is 1.0d-8).
!
!     v = 3 by ndim workspace.
!
!
!     variables returned:
!     ------------------
!
!     eval1 = eigenvalues of a (sorted in increasing algebraic value);
!
!     idegen1 = degeneracies (number of times repeated) for eigenvalues;
!
!     evec1 = eigenvectors of a (in columns of evec1);
!
!     error codes:  ierror=0 - successful calculation.
!                   ierror=1 - no convergence in qr iteration.
!
!
!     programmed by s. l. dixon, oct., 1991.
!
!
      implicit none
#include "sebomd.dim"
#include "se_pass1.h"
      integer ndim, nevec1, ierror
      integer idegen1(*)
      double precision tolera
      double precision a(ndim,*),v(3,*),eval1(*),evec1(ndim,*)

      integer i, j ,k, m
      double precision betah, ancol, anorm, aicol, toltmp, sum, bsum
!
!     flag for whether or not to compute eigenvalues:
!
!
      if(ndim.eq.1)then
        eval1(1) = a(1,1)
        idegen1(1) = 1
        evec1(1,1) = 1.0d0
        return
      endif
!
!-rdc c     tridiagonalize the matrix a.  this will overwrite the diagonal
!     and subdiagonal of a with the tridiagonalized version.  the
!     householder vectors are returned in the rows above the diagonal,
!     and the betahs are returned below the subdiagonal.
!
      call se_tridi(ndim,v,a)
!
!     compute norm of tridiagonal matrix from the "largest" column.
!
      anorm = abs(a(1,1)) + abs(a(2,1))
      if(ndim.gt.2)then
        do i=2,ndim-1
          aicol = abs(a(i-1,i)) + abs(a(i,i)) + abs(a(i+1,i))
          anorm = max(anorm,aicol)
        enddo
      endif
      ancol = abs(a(ndim-1,ndim)) + abs(a(ndim,ndim))
      anorm = max(anorm,ancol)
!
!     get eigenvalues and degeneracies of the tridiagonal matrix a.
!     if the calling routine has not supplied a tolerance factor for
!     off-diagonal entries in the qr iteration, a default of 1.0d-8
!     will be used.
!
      toltmp = tolera
      if(toltmp.le.0.0d0) toltmp = 1.0d-8
      if(doeval) call se_eigval(ndim,a,v,toltmp,anorm,eval1,ierror)
      if(ierror.ne.0) return
!
!     determine degeneracies of eigenvalues.
!
      call se_degen(ndim,eval1,toltmp,anorm,idegen1)
!
!     get eigenvectors of tridiagonalized version of a.
!
      if(nevec1.le.0) return
      call se_eigvec(ndim,nevec1,a,v,toltmp,anorm,eval1,idegen1,evec1)
!
!     premultiply evec1 by the householder matrix used to tridiagonalize
!     a.  this transforms eigenvectors of the tridiagonal matrix to
!     those of the original matrix a.  see subroutine tridi for
!     storage of householder transformation.
!
      if(ndim.gt.2)then
        do k=1,ndim-2
          v(1,k) = a(k+2,k)
        enddo
!
!       swap storage so that the expensive triple loop below doesn't
!       have to jump across columns to get entries of a.
!
        do i=2,ndim
          do j=1,i-1
            a(i,j) = a(j,i)
          enddo
       enddo
        do j=1,nevec1
          do m=2,ndim-1
            k = ndim - m
            betah = v(1,k)
            if(abs(betah).lt.1.0d-50) cycle
            sum = 0.0d0
            do i=k+1,ndim
              sum = sum + a(i,k)*evec1(i,j)
            enddo
            bsum = betah*sum
            do i=k+1,ndim
              evec1(i,j) = evec1(i,j) - a(i,k)*bsum
            enddo
          enddo
        enddo
      endif
      return
      end
!--------------------------------------------------------------------------------
      subroutine se_tridi(ndim,v,a)
!
!     tridiagonalizes a real, symmetric matrix a by the method of
!     householder (j. h. wilkinson, the computer journal, vol. 3,
!     p. 23 (1960)).  ndim is the order of a.  the diagonal and
!     subdiagonal of a are overwritten with the tridiagonalized
!     version of a.  the vectors used in each householder
!     transformation are stored above the diagonal in the first
!     ndim-2 rows of a.  the betahs are returned below the subdiagonal
!     of a.  v is a workspace array.
!
!     programmed by s. l. dixon, oct., 1991.
!
!
      implicit none
      integer ndim
      double precision a(ndim,*),v(3,*)

      integer i, j, k
      double precision thresh, alpha2, aptemp, alpha, betah, sgn
      double precision s, dot12, bhalf, sum, bh12, akv
!
!     thresh will be used as a threshold to determine if a value should
!     be considered to be zero.  this can be changed by the user.
!
      thresh = 1.0d-50
!
!     if a is 2 by 2 or smaller, then it is already tridiagonal -- no
!     need to continue.
!
      if(ndim.le.2) return
      do k=1,ndim-2
!
!       determine the vector v used in the householder transformation p.
!       for each value of k the householder matrix p is defined as:
!
!                   p = i - betah*v*v'
!
!
!       construct a householder transformation only if there is a nonzero
!       off-diagonal element below a(k,k).
!
        alpha2 = 0.0d0
        do i=k+1,ndim
          v(1,i) = a(i,k)
          alpha2 = alpha2 + v(1,i)**2
        enddo
        aptemp = alpha2 - v(1,k+1)**2
        alpha = dsqrt(alpha2)
        if(alpha.ge.thresh)then
          betah = 1.0d0/(alpha*(alpha + abs(v(1,k+1))))
          sgn = sign(1.0d0,v(1,k+1))
          v(1,k+1) = v(1,k+1) + sgn*alpha
!
!-rdc c         now overwrite a with p'*a*p.  the entries below the subdiagonal
!         in the kth column are zeroed by the premultiplication by p'.
!         these entries will be left alone to save time.
!
          akv = aptemp + a(k+1,k)*v(1,k+1)
          s = betah*akv
          a(k+1,k) = a(k+1,k) - s*v(1,k+1)
!
!         now the submatrix consisting of rows k+1,ndim and columns k+1,ndim
!         must be overwritten with the transformation.
!
          dot12 = 0.0d0
          bhalf = betah*0.5d0
          do i=k+1,ndim
            sum = 0.0d0
            do j=k+1,i
              sum = sum + a(i,j)*v(1,j)
            enddo
            if(i.lt.ndim)then
              do j=i+1,ndim
!
!               an upper triangular entry of a will be required.  must use
!               the symmetric entry in the lower triangular part of a.
!
                sum = sum + a(j,i)*v(1,j)
              enddo
            endif
            v(2,i) = betah*sum
            dot12 = dot12 + v(1,i)*v(2,i)
          enddo
          bh12 = bhalf*dot12
          do i=k+1,ndim
            v(2,i) = v(2,i) - bh12*v(1,i)
          enddo
          do j=k+1,ndim
            do i=j,ndim
              a(i,j) = a(i,j) - v(1,i)*v(2,j) - v(2,i)*v(1,j)
            enddo
!
!           store v(1,j) above the diagonal in row k of a
!
            a(k,j) = v(1,j)
          enddo
!
!         store betah below the subdiagonal of a.
!
          a(k+2,k) = betah
        else
!
!         no householder transformation is necessary because the off-
!         diagonals are all essentially zero.
!
          a(k+2,k) = 0.0d0
          do j=k+1,ndim
            a(k,j) = 0.0d0
          enddo
        endif
      enddo
      end
!--------------------------------------------------------------------------------
      subroutine se_eigval(ndim,a,betah,tolera,anorm,eval1,ierror)
!
!     qr routine for the determination of all the eigenvalues
!     of the ndim by ndim symmetric, tridiagonal matrix a.
!
!     input:
!
!     ndim   = size of matrix a.
!     a      = ndim by ndim symmetric tridiagonal matrix.
!     betah   = 3 by ndim workspace.
!     tolera = small number used to determine when off-diagonal
!              elements are essentially zero.
!     anorm  = absolute column norm of tridiagonal matrix a.
!
!
!     returned:
!
!     eval1   = eigenvalues of a in ascending order.
!     ierror = 1 if qr iteration did not converge; 0 otherwise.
!
!     programmed by s. l. dixon.
!
!
      implicit none
      integer ndim, ierror
      double precision a(ndim,*),betah(3,*),eval1(*)
      double precision tolera, anorm

      integer itmax, iter, kupper
      integer i, j, k, l, n, jmax, isort
      double precision epslon, bterm, root1, root2, discr
      double precision signd, p, r, t, w, shift, d, c, s, cc, ss, cs
      double precision csw, ak1, etemp

      ierror = 0
      itmax = 20
!
!     tolerance for off-diagonal elements:
!
      epslon = tolera*anorm
!
!     copy diagonal elements of a to eval1, and subdiagonal elements
!     to betah.
!
      eval1(1) = a(1,1)
      betah(1,1) = a(2,1)
      if(ndim.gt.2)then
        do i=2,ndim-1
          eval1(i) = a(i,i)
          betah(1,i) = a(i+1,i)
        enddo
      endif
      eval1(ndim) = a(ndim,ndim)
!
!     each qr iteration will operate on the unreduced tridiagonal
!     submatrix with upper left element (l,l) and lower right element
!     (n,n).
!
      l = 1
      n = ndim
      iter = 0
!
!     find the smallest unreduced submatrix with lower right corner at
!     (n,n).  i.e., search upward for a betah that is zero.
!
 80   kupper = n-l
      do k=1,kupper
        i = n-k
        if(abs(betah(1,i)).le.epslon)then
          l = i+1
          go to 150
        endif
      enddo
!
!     if we get to the next statement, then there are no zero off-diagonals
!     for the submatrix with upper left a(l,l) and lower right a(n,n).
!     we can still get eigenvalues if the matrix is 2 by 2 or 1 by 1.
!     otherwise, do another qr iteration provided itmax cycles have
!     not occurred.
!
      if(l.eq.n.or.l.eq.n-1)then
        go to 150
      else
        if(iter.eq.itmax)then
          ierror = 1
          return
        else
          go to 200
        endif
      endif
!
!     if we get to 150 then a(l,l-1) is zero and the unreduced submatrix
!     has upper left at a(l,l) and lower right at a(n,n).  we can
!     extract one eigenvalue if this matrix is 1 by 1 and 2 eigenvalues
!     if it is 2 by 2.
!
  150 if(l.eq.n)then
!
!       it's a 1 by 1 and eval1(n) is an eigenvalue.  if l=2 or 1 we are
!       done.  otherwise, update n, reset l and iter, and repeat the
!       search.
!
        if(l.le.2)then
          go to 500
        else
          n = l-1
          l = 1
          iter = 0
          go to 80
        endif
      elseif(l.eq.n-1)then
!
!-rdc c       the unreduced submatrix is a 2 by 2.  overwrite eval1(n-1)
!       and eval1(n) with the eigenvalues of the lower right 2 by 2.
!
        bterm = eval1(n-1) + eval1(n)
        root1 = bterm*0.5d0
        root2 = root1
        discr = bterm**2 - 4.0d0*(eval1(n-1)*eval1(n)-betah(1,n-1)**2)
        if(discr.gt.0.0d0)then
          d = dsqrt(discr)*0.5d0
          root1 = root1 - d
          root2 = root2 + d
        endif
        eval1(n-1) = root1
        eval1(n) = root2
!
!       see if we are done.  if not, reset n, l, and iter and look
!       for next unreduced submatrix.
!
        if(l.le.2)then
          go to 500
        else
          n = l-1
          l = 1
          iter = 0
          go to 80
        endif
      else
!
!       an eigenvalue was found and the new unreduced matrix limits
!       n and l are set.  do a qr iteration on new matrix.
!
        iter = 0
        go to 200
      endif
!
!     qr iteration begins here.
!
 200  iter = iter + 1
!
!     use eigenvalues of the lower right 2 by 2 to compute shift.  shift
!     by the eigenvalue closest to eval1(n).
!
      d = (eval1(n-1) - eval1(n))*0.5d0
      signd = 1.0d0
      if(d.lt.0.0d0) signd = -1.0d0
      shift = eval1(n) + d - signd*dsqrt(d*d + betah(1,n-1)**2)
      p = eval1(l) - shift
      r = betah(1,l)
      t = eval1(l)
      w = betah(1,l)
!
!-rdc c     overwrite a with q'*a*q.
!
      do k=l,n-1
        d = dsqrt(p*p + r*r)
        c = p/d
        s = r/d
        if(k.ne.l) betah(1,k-1) = d
        cc = c*c
        ss = 1.0d0 - cc
        cs = c*s
        csw = 2.0d0*cs*w
        ak1 = eval1(k+1)
        eval1(k) = cc*t + csw + ss*ak1
        p = (cc - ss)*w + cs*(ak1 - t)
        t = ss*t - csw + cc*ak1
        r = s*betah(1,k+1)
        w = c*betah(1,k+1)
      enddo
      betah(1,n-1) = p
      eval1(n) = t
!
!     go back and see if l and n need to be updated.
!
      go to 80
!
!     sort eigenvalues in ascending algebraic order.
!
 500  continue
      do i=2,ndim
        jmax = ndim-i+1
        isort = 0
        do j=1,jmax
          if(eval1(j).gt.eval1(j+1))then
            etemp = eval1(j)
            eval1(j) = eval1(j+1)
            eval1(j+1) = etemp
            isort = 1
          endif
        enddo
        if(isort.eq.0) return
      enddo

      return
      end
!--------------------------------------------------------------------------------
      subroutine se_degen(ndim,eval1,tolera,anorm,idegen1)
!
!     determines degeneracies of the eigenvalues.
!
!     input:
!
!     ndim   = size of matrix being diagonalized.
!     eval1   = sorted eigenvalues (increasing value).
!     tolera = same tolerance used to determine eigenvalues.
!     anorm  = absolute column norm of tridiagonal matrix.
!
!
!     returned:
!
!     idegen1 = degeneracies of eigenvalues.
!
      implicit none
      integer ndim
      integer idegen1(*)
      double precision tolera, anorm
      double precision eval1(*)

      integer nsame, i, j
      double precision dtoler, diff
!
!     determine degeneracies of eigenvalues.  adjacent eigenvalues
!     will be considered to be degenerate when they differ by less
!     than dtoler.
!
      dtoler = max(anorm*dsqrt(tolera),1.0d-8)
      nsame = 1
      do i=2,ndim
        diff = abs(eval1(i-1) - eval1(i))
        if(diff.le.dtoler)then
!
!         eigenvalues i-1 and i are degenerate.
!
          nsame = nsame + 1
          if(i.eq.ndim)then
!
!           we've come to the last requested eigenvalue, and it's time
!           to assign degeneracies for the block ending with the ith
!           eigenvalue.
!
            do j=i-nsame+1,i
              idegen1(j) = nsame
            enddo
          endif
!
!         go to the next eigenvalue (if there are any left) and see if
!         it's degenerate with the nsame eigenvalues we've already
!         found.
!
          cycle
        else
!
!         either eigenvalue i-1 is nondegenerate or it's the last
!         eigenvalue in a degenerate block.  correspondingly, assign the
!         proper degeneracy to i-1 or to each eigenvalue in the block.
!
          do j=i-nsame,i-1
            idegen1(j) = nsame
          enddo
          nsame = 1
!
!         if i=ndim then it must be the case that this last eigenvalue
!         is nondegenerate.
!
          if(i.eq.ndim) idegen1(i) = 1
        endif
      enddo
      return
      end
!--------------------------------------------------------------------------------
      subroutine se_eigvec(ndim,nevec1,a,awork,tolera,anorm,eval1,
     .                  idegen1,evec1)
!
!     inverse iteration routine for eigenvector determination.
!     calculates the eigenvectors of an ndim by ndim symmetric,
!     tridiagonal matrix a.
!
!     input:
!
!     ndim   = size of matrix a.
!     nevec1  = number of eigenvectors required.
!     a      = ndim by ndim tridiagonal matrix.
!     tolera = same tolerance used to determine eigenvalues.
!     anorm  = absolute column norm of tridiagonal matrix a.
!     eval1   = sorted eigenvalues of a.
!     idegen1 = degeneracies of eigenvalues.
!
!
!     returned:
!
!     evec1   = eigenvectors of tridiagonal matrix (in columns).
!
!     programmed by s. l. dixon, oct., 1991.
!
!
      implicit none
      integer ndim, nevec1
      integer idegen1(*)
      double precision tolera, anorm
      double precision a(ndim,*),awork(3,*),eval1(*),evec1(ndim,*)

      integer irand, ivect, nfail, nprtrb, ndegen, mstart
      integer m, i, k, l, jstart
      double precision epslon, z, rndom, adiag, asub, t, etemp, enorm
      double precision atemp, einv
      logical orth

      irand = 13876532
!
!     compute threshold epslon which will be used if the inverse iteration
!     matrix is singular.
!
      epslon = anorm*tolera
!
!     when degeneracies occur, there are rare instances when the
!     degenerate block of eigenvectors are not linearly independent.
!     in these cases, an additional pass through the inverse iteration
!     (with a new set of random numbers) is carried out, i.e., control
!     passes to statement 40.  ivect will keep track of the current
!     starting eigenvector when additional passes are necessary.
!
      ivect = 1
      nfail = 0
      nprtrb = 0
!
!     do one iteration for each eigenvector.
!
 40   orth = .true.
      ndegen = 0
      mstart = ivect
      do m=mstart,nevec1
        if(idegen1(m).gt.1) ndegen = ndegen + 1
        z = eval1(m)
!
!       if the inverse iteration has failed twice due to non-orthogonality
!       of degenerate eigenvectors, perturb the eigenvalue by a small
!       amount.
!
        if(nfail.ge.2)then
          nprtrb = nprtrb + 1
          z = z + 0.001d0*dfloat(nprtrb)*tolera
        endif
!
!       store the tridiagonal entries of the inverse iteration matrix in
!       the 3 by ndim workspace awork.
!
        awork(1,1) = 0.0d0
        awork(2,1) = a(1,1) - z
        awork(3,1) = a(2,1)
        if(ndim.gt.2)then
          do i=2,ndim-1
            awork(1,i) = a(i,i-1)
            awork(2,i) = a(i,i) - z
            awork(3,i) = a(i+1,i)
          enddo
        endif
        awork(1,ndim) = a(ndim,ndim-1)
        awork(2,ndim) = a(ndim,ndim) - z
        awork(3,ndim) = 0.0d0
!
!       assign inverse iteration vector from random numbers.
!
        do i=1,ndim
          call se_random1(irand,rndom)
          rndom = 2.0d0*(rndom - 0.5d0)
          evec1(i,m) = rndom*tolera
        enddo
!
!       carry out forward gaussian elimination with row pivoting
!       on the inverse iteration matrix.
!
        do k=1,ndim-1
          adiag = abs(awork(2,k))
          asub = abs(awork(1,k+1))
          if(adiag.ge.asub)then
!
!           use pivotal element from row k.
!
            if(awork(2,k).eq.0.0d0) awork(2,k) = epslon
            t = awork(1,k+1)/awork(2,k)
            awork(1,k+1) = 0.0d0
            awork(2,k+1) = awork(2,k+1) - t*awork(3,k)
!
!           left-justify equation k so that diagonal entry is stored
!           in awork(1,k).
!
            awork(1,k) = awork(2,k)
            awork(2,k) = awork(3,k)
            awork(3,k) = 0.0d0
!
!           operate on vector as well.
!
            evec1(k+1,m) = evec1(k+1,m) - t*evec1(k,m)
          else
!
!           use pivotal element from row k+1 and swap rows k and k+1.
!
            if(awork(1,k+1).eq.0.0d0) awork(1,k+1) = epslon
            t = awork(2,k)/awork(1,k+1)
            atemp = awork(3,k) - t*awork(2,k+1)
            awork(1,k) = awork(1,k+1)
            awork(2,k) = awork(2,k+1)
            awork(3,k) = awork(3,k+1)
            awork(1,k+1) = 0.0d0
            awork(2,k+1) = atemp
            awork(3,k+1) = -t*awork(3,k+1)
!
!           operate on vector and swap entries.
!
            etemp = evec1(k+1,m)
            evec1(k+1,m) = evec1(k,m) - etemp*t
            evec1(k,m) = etemp
          endif
        enddo
!
!       forward elimination complete.  back substitute to get solution.
!-rdc c       overwrite column m of evec1 with solution.
!
        if(awork(2,ndim).eq.0.0d0) awork(2,ndim) = epslon
        evec1(ndim,m) = evec1(ndim,m)/awork(2,ndim)
        etemp = evec1(ndim-1,m) - awork(2,ndim-1)*evec1(ndim,m)
        evec1(ndim-1,m) = etemp/awork(1,ndim-1)
        enorm = evec1(ndim,m)**2 + evec1(ndim-1,m)**2
        if(ndim.gt.2)then
!
!         caution: problem loop for some ibm rs/6000 compilers.  value
!         of k can get lost when optimize flag is used.
!
          do l=1,ndim-2
            k = ndim-l-1
            etemp = evec1(k,m) - awork(2,k)*evec1(k+1,m)
     .              - awork(3,k)*evec1(k+2,m)
            evec1(k,m) = etemp/awork(1,k)
            enorm = enorm + evec1(k,m)**2
          enddo
        endif
        einv = 1.0d0/dsqrt(enorm)
!
!       normalize eigenvector.
!
        do i=1,ndim
          evec1(i,m) = evec1(i,m)*einv
        enddo
!
!       if we have come to the end of a degenerate block of eigenvectors,
!       orthogonalize the block.
!
        if(ndegen.gt.1)then
          if(ndegen.eq.idegen1(m).or.m.eq.nevec1)then
            jstart = m-ndegen+1
            call se_orthog1(ndim,ndegen,jstart,evec1,orth)
            if(orth)then
              nfail = 0
              nprtrb = 0
!
!             the degenerate vectors were linearly independent and were
!             successfully orthogonalized.
!
              ivect = ivect + ndegen
              ndegen = 0
            else
!
!             the block is apparently not linearly independent.  go back
!             and repeat the inverse iteration for these vectors.  after
!             an independent set has been found, any additional eigenvectors
!             will be determined.
!
              nfail = nfail + 1
              go to 40
            endif
          endif
        endif
!
!       the current eigenvector should be okay if it is nondegenerate.
!
        if(idegen1(m).eq.1) ivect = ivect + 1
      enddo
      return
      end
!--------------------------------------------------------------------------------
      subroutine se_orthog1(ndim,nvect,jstart,vect,orth)
!
!     constructs a set of orthonormal vectors from the nvect linearly
!     independent, normalized vectors in the array vect.  the vectors
!     should be stored columnwise, starting in column jstart.  vect is
!     overwritten with the orthonormal set.  all vectors are ndim by 1.
!     orth is returned with a value of .true. if the set was linearly
!     independent and .false. otherwise.
!
!     programmed by s. l. dixon.
!
!
      implicit none
      integer ndim, nvect, jstart
      double precision vect(ndim,*)
      logical orth

      integer j, jstop, i, jprev
      double precision ortest, dot, vjnorm

!
      orth = .true.
      ortest = 1.0d-8
!
!     begin orthogonalization.
!
      jstop = jstart + nvect - 1
      do j=jstart,jstop
        if(j.gt.jstart)then
!
!         subtract off components of previously determined orthogonal
!         vectors from the vector in column j.
!
          do jprev=jstart,j-1
            dot = 0.0d0
            do i=1,ndim
              dot = dot + vect(i,jprev)*vect(i,j)
            enddo
            do i=1,ndim
              vect(i,j) = vect(i,j) - dot*vect(i,jprev)
            enddo
          enddo
        endif
!
!       normalize column j.
!
        vjnorm = 0.0d0
        do i=1,ndim
          vjnorm = vjnorm + vect(i,j)**2
        enddo
        vjnorm = dsqrt(vjnorm)
!
!       if the norm of this vector is too small then the vectors are
!       not linearly independent.
!
        if(vjnorm.lt.ortest)then
          orth = .false.
          return
        endif
        do i=1,ndim
          vect(i,j) = vect(i,j)/vjnorm
        enddo
      enddo
      return
      end
!--------------------------------------------------------------------------------
      subroutine se_random1(ix,x)
!
!     generates integer*4 and double precision random numbers on
!     the intervals:
!
!                       0 < ix < (2**31)-1
!              and      0.0 < x < 1.0
!
!     on the first call ix should satisfy the top inequality.
!
!     numbers are generated using the relation,
!
!               ix = ix*ic (modulo (2**31)-1),  where ic = 7**5
!
!
      implicit none
      integer ix
      double precision x

      integer i15, i16, i31_1, ic
      integer ix16, ix16ic, ixic31, i16rmd
!
!
!     initialize:     i15 = 2**15
!                     i16 = 2**16
!                   i31_1 = 2**31 - 1
!                      ic = 7**5
!
      data i15 /32768/
      data i16 /65536/
      data i31_1 /2147483647/
      data ic /16807/
!
      save i15,i16,i31_1,ic
!
      ix16 = ix/i16
      i16rmd = ix - ix16*i16
!
!     note that ix = ix16*i16 + i16rmd    (i16rmd = 16-bit remainder)
!
      ix16ic = ix16*ic
      ixic31 = ix16ic/i15
!
!     note:   ix*ic = (ix16*i16 + i16rmd)*ic
!                   = (ix16*ic*i16) + (i16rmd*ic)
!                   = (ix16ic*i16 ) + (i16rmd*ic)
!                   = (  term 1   ) + ( term 2  )
!     and,
!
!     ix16ic = ((ix16ic/i15)*i15)  +  (ix16ic - (ix16ic/i15)*i15))
!             = (  ixic31  )*i15)  +  (ix16ic - (  ixic31  )*i15 )
!                                     (     15-bit remainder     )
!
!     therefore,  term 1 = ((ixic31*i15) + (ix16ic - ixic31*i15))*i16
!
!     then,
!
!              (   term a     )   (        term b          )   ( term c  )
!     ix*ic = ((ixic31*i16*i15) + (ix16ic - ixic31*i15)*i16) + (i16rmd*ic)
!           = (                  term 1                    ) + ( term 2  )
!
!
!     note that term b and term c are both less than 2**31 - 1.  only
!     term a has the possibility of exceeding 2**31 - 1.  but since
!     i16*i15 = 2**31, the factor ixic31 indicates exactly how many times
!     term a "wraps" around the interval (0,2**31 - 1).  thus, for the
!     modulo operation, term a may be replaced by ixic31.  the sum of
!     term a and term b might exceed 2**31 - 1, but we can substract
!     2**31 - 1 from one of them to prevent this from happening.
!
      ix = ixic31 + ((ix16ic-ixic31*i15)*i16 - i31_1) + i16rmd*ic
!
!     add i31_1 back in if the subtraction made ix negative.
!
      if(ix.lt.0) ix = ix + i31_1
!
!     make x random on (0.0,1.0) by multiplying ix by 1.0/i31_1
!
      x = dfloat(ix)*4.6566128752458d-10
      return
      end
