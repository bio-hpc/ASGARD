
common /howmny/ natom,  ntypes, nbonh,  ntheth, nphih, &
      nnb,    nres,   nbona,  ntheta, nphia, &
      numbnd, numang, nptra,  natyp,  nphb, &
      natsys, nr3,    ns3,    nb3,    nbmax, &
      nbmax2, mbona,  mtheta, mphia

!     ----- The following are pointers to the first element of
!     ----- various arrays. The pointer name is derived by prefixing
!     ----- array name by m. For example amass begins at x(mamass)
!     ----- Exceptions are the following arrays with six character
!     ----- names:
!     -----          igraph, isymbl, igroup
!     ----- In these cases m replaces the first character

common /pointr/   mchrg,  mamass,  mrk,    mreq,   mtk, &
      mteq,   mfk,     mfpk,   mqeq,   mpk, &
      mpn,    mphase,  msolty, mcn1,   mcn2, &
      masol,  mbsol,   mhbcut, mxref,  msf, &
      momega, mgraph,  miac,   miblo,  mico, &
      mlbres, mipres,  mibh,   mjbh,   micbh, &
      miba,   mjba,    micba,  mith,   mjth, &
      mkth,   micth,   mita,   mjta,   mkta, &
      micta,  miph,    mjph,   mkph,   mlph, &
      micph,  mipa,    mjpa,   mkpa,   mlpa, &
      micpa,  minb,    msymbl, mitree, mgroup, &
      migres, miar1,   mx,     mf,     mh, &
      mnbel,  mxbel,   miar2,  mcscr,  mcval, &
      mcvec,  mdd,     mb,     mroots, mvect, &
      mxdir,  mwref,   migrp2, ma,     mwr, &
      mwi,    mz,      mfv1,   miv1,   mhrad, &
      mgam,   mwinv,   mkpvt,  mxinit, mcn114, &
      mcn214, mjoin,   mrotat, mpol
