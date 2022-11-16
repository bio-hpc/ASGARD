
      double precision pdiag3, pdiat3 
      double precision pdiag2, pdiat2 
      double precision pdiag1, pdiat1 
      integer ijmat3, ipair3, ip13
      integer ijmat2, ipair2, ip12
      integer ijmat1, ipair1, ip11
      common /se_glbmat4/ 
     $     pdiag3(mxdiag), pdiat3(mxdiat), ijmat3(mxdiat),
     $     ipair3(mbpair), ip13(maxatm),
     $     pdiag2(mxdiag), pdiat2(mxdiat), ijmat2(mxdiat),
     $     ipair2(mbpair), ip12(maxatm),
     $     pdiag1(mxdiag), pdiat1(mxdiat), ijmat1(mxdiat),
     $     ipair1(mbpair), ip11(maxatm) 
