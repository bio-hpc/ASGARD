
      double precision dcbuff1, dcbuff2
      integer nncore
      integer ncore,icorel
      integer icorel1, neighbor
      integer neighborn
      integer ncores
      common /se_subpar/ dcbuff1, dcbuff2,
     &     nncore,
     &     ncore(maxsub),icorel(maxres),
     &     icorel1(maxsub+1), neighbor(maxres+1),
     $     neighborn(maxres+1),
     &     ncores
