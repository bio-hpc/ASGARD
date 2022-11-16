
      double precision qpmec,qpme
      double precision thetapme,bfacxpme
      double precision bfacypme,bfaczpme
      double precision recip1,recip2,recip3,betapme,betapme2,dnspline
      double precision dk1pme,dk2pme,dk3pme
      integer mmax
      integer nspline,k1pme,k2pme,k3pme,k1pmem1,k2pmem1,k3pmem1
      integer k1pmek2pme,k123pme,k123pme2,k1pmenspl,k2pmenspl,k3pmenspl
      common /se_pmeee/
     &     qpmec(maxkpme32),qpme(maxkpme3),
     &     thetapme(maxkpme3),bfacxpme(0:maxkpme-1),
     &     bfacypme(0:maxkpme-1),bfaczpme(0:maxkpme-1),
     &     recip1(3),recip2(3),recip3(3),betapme,betapme2,dnspline,
     &     dk1pme,dk2pme,dk3pme,
     &     mmax,
     &     nspline,k1pme,k2pme,k3pme,k1pmem1,k2pmem1,k3pmem1,
     &     k1pmek2pme,k123pme,k123pme2,k1pmenspl,k2pmenspl,k3pmenspl
