      subroutine se_diagp(norbs,nocc,f,e,fov,v,c)
!
!     pseudo-diagonalization routine for the matrix f.  this routine
!     is used when eigenvalues in scf calculation have sufficiently
!     converged.  their values are used to approximately eliminate
!     significant off-diagonal elements in the occupied-virtual
!     intersection of the mo's with f.  see j.j.p. stewart, c. csaszar,
!     and p. pulay, j. comp. chem, vol. 3, 227-228 (1982).
!
!     norbs = the total number of molecular orbitals.
!
!      nocc = the number of occupied molecular orbitals.
!
!         f = matrix to be pseudodiagonalized.  the storage is linear
!             and requires at least norbs*norbs slots.  the 2-d to 1-d
!             correspondence is (i,j) --> (i+(j-1)*norbs).  only the
!             lower triangle (i.ge.j) is needed and only the upper
!             triangle is altered.
!
!         e = eigenvalues of f.  requires norbs slots.
!
!       fov = temporary storage for partially diagonalized version of f.
!             requires as much storage as f.
!
!         v = temporary storage.  requires norbs slots.
!
!         c = eigenvectors for f.  storage is the same as f.  this routine
!             rotates the eigenvectors to help eliminate the occupied-
!             virtual intersection.
!
!
      implicit none
      double precision f(*),e(*),fov(*),v(*),c(*)
      integer norbs, nocc

      integer i, iterm
      integer j, jterm
      integer k, ki, kj
      double precision epslon, tolera, thresh
      double precision vi, fij, fijmax
      double precision eocc, evirt, d, sigd, w, sine, cosine
      double precision cki, ckj
!

      epslon = 1.0d-8
      tolera = 0.05d0
!
!     copy lower triangle of f to upper triangle for quicker access.
!
      iterm = 0
      do i=1,norbs
        jterm = 0
        do j=1,i
          f(j+iterm) = f(i+jterm)
          jterm = jterm + norbs
        enddo
        iterm = iterm + norbs
      enddo
!
!     form the matrix fov from the occupied-virtual intersection of
!     the molecular orbitals:
!
!                  fov = c(occ)'*f*c(virt)
!
      fijmax = 0.0d0
      jterm = nocc*norbs
      do j=nocc+1,norbs
        iterm = 0
        do i=1,norbs
          vi = 0.0d0
          do k=1,norbs
            vi = vi + f(k+iterm)*c(k+jterm)
          enddo
          v(i) = vi
          iterm = iterm + norbs
        enddo
        iterm = 0
        do i=1,nocc
          fij = 0.0d0
          do k=1,norbs
            fij = fij + c(k+iterm)*v(k)
          enddo
          if(abs(fij).gt.fijmax) fijmax = abs(fij)
          fov(i+jterm) = fij
          iterm = iterm + norbs
        enddo
        jterm = jterm + norbs
      enddo
      thresh = tolera*fijmax
!
!     do a pseudo-jacobi rotation to approximately eliminate significant
!     elements of fov.
!
      iterm = 0
      do i=1,nocc
        jterm = nocc*norbs
        do j=nocc+1,norbs
          if(abs(fov(i+jterm)).lt.thresh) go to 210
          eocc = e(i)
          evirt = e(j)
          fij = fov(i+jterm)
          d = eocc - evirt
          if(d.ne.0.0d0)then
            if(abs(fij/d).lt.epslon) go to 210
          endif
          d = d*0.5d0
          sigd = sign(1.0d0,d)
          w = -sigd*fij/dsqrt(fij**2 + d**2)
          sine = w/dsqrt(2.0d0*(1.0d0 + dsqrt(1.0d0 - w*w)))
          cosine = dsqrt(1.0d0 - sine**2)
!
!         rotate pseudo-eigenvectors.
!
          do k=1,norbs
            ki = k + iterm
            kj = k + jterm
            cki = c(ki)
            ckj = c(kj)
            c(ki) = cki*cosine - ckj*sine
            c(kj) = cki*sine + ckj*cosine
          enddo
 210      jterm = jterm + norbs
        enddo
        iterm = iterm + norbs
      enddo
      return
      end
