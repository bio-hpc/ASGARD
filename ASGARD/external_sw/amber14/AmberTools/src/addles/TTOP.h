c
c alternate topology definitions. 
c
      character*80 tititl
      character*4 tigraph(maxnatom),tlabres(maxnres),tisymbl(maxnatom)
      character*4 titree(maxnatom)
c
      integer tnatom,tntypes,tnbonh,tmbona,tntheth,tmtheta
      integer tnphih,tmphia,tnhparm,tnparm,tnnb,tnres,tnbona,tntheta
      integer tnphia,tnumbnd,tnumang,tnptra,tnatyp,tnphb,tifpret,tnbper
      integer tngper,tndper,tmbper,tmgper,tmdper,tifbox,tnmxrs,tifcap
c
      integer tiac(maxnatom),tnumex(maxnatom),tico(maxntypes*maxntypes)
      integer tipres(maxnres),tibh(maxbnd),tjbh(maxbnd),ticbh(maxbnd) 
      integer tib(maxbnd),tjb(maxbnd),ticb(maxbnd) 
      integer tith(maxang),tjth(maxang),tkth(maxang),ticth(maxang) 
      integer tit(maxang),tjt(maxang),tkt(maxang),tict(maxang) 
      integer tiph(maxdih),tjph(maxdih),tkph(maxdih),tlph(maxdih)
      integer ticph(maxdih),tip(maxdih) 
      integer tjp(maxdih),tkp(maxdih),tlp(maxdih),ticp(maxdih) 
      integer tnatex(maxnext),tjoin(maxnatom),tirotat(maxnatom)
      integer tiptres,tnspm,tnpsol,tnsp(maxnspm),tnatcap
c
c pert stuff
c
      integer tibper(maxbnd),tjbper(maxbnd),ticbper(maxbnd)
      integer t2icbper(maxbnd)
      integer titper(maxang),tjtper(maxang),tktper(maxang)
      integer tictper(maxang),t2ictper(maxang)
      integer tipper(maxdih),tjpper(maxdih),tkpper(maxdih)
      integer ticpper(maxdih),tlpper(maxdih),t2icpper(maxdih)
      integer tiaper(maxnatom),tiacper(maxnatom)
      integer tinegp(maxdih),tjnegp(maxdih),tknegp(maxdih)
      integer tlnegp(maxdih)
      real*8 tcgper(maxnatom),talmper(maxnatom)
      character*4 tlabre2(maxnres),tigrper(maxnatom),tismper(maxnatom)
c
c here is a temporary storage for the SIGN of the torsion atoms.
c make them all positive when read, but need to replace the sign
c upon writing.
c       
      integer tinegh(maxdih),tjnegh(maxdih),tknegh(maxdih)
      integer tlnegh(maxdih)
      integer tineg(maxdih),tjneg(maxdih),tkneg(maxdih),tlneg(maxdih)

c
c ttimass(:) = the perturbed mass for TI in the LES regions
      real*8 tchrg(maxnatom),tamass(maxnatom),trk(maxbnd),
     &       treq(maxbnd),trborn(maxnatom),tfs(maxnatom),
     &       ttimass(maxnatom)
      real*8 ttk(maxang),tteq(maxang),tpk(maxdih),tpn(maxdih),
     &       tphase(maxdih)
      real*8 tsolty(maxnatyp),tcn1(maxntypes*(maxntypes+1)/2)
      real*8 tcn2(maxntypes*(maxntypes+1)/2),tasol(maxnphb),
     &       tbsol(maxnphb)
      real*8 thbcut(maxnphb),tbeta,tbox(3),tcutcap,txcap,tycap,tzcap
c
c all in common
c
      common /ct1/ tititl,tigraph,tlabres,tisymbl,titree
      common /ctp1/ tlabre2,tigrper,tismper
c
      common /it1/ tnatom,tntypes,tnbonh,tmbona,tntheth,tmtheta
      common /it2/ tnphih,tmphia,tnhparm,tnparm,tnnb,tnres,tnbona,
     &             tntheta
      common /it3/ tnphia,tnumbnd,tnumang,tnptra,tnatyp,tnphb,tifpret,
     &             tnbper
      common /it4/ tngper,tndper,tmbper,tmgper,tmdper,tifbox,tnmxrs,
     &             tifcap
      common /it5/ tiac,tnumex,tico,tipres,tibh,tjbh,ticbh
      common /it6/ tib,tjb,ticb,tith,tjth,tkth,ticth
      common /it7/ tit,tjt,tkt,tict,tiph,tjph,tkph,tlph,ticph
      common /it8/ tip,tjp,tkp,tlp,ticp,tnatex,tjoin,tirotat
      common /it9/ tiptres,tnspm,tnpsol,tnsp,tnatcap
      common /it10/ tineg,tjneg,tkneg,tlneg,tinegh,tjnegh,tknegh,tlnegh
      common /itp1/ tibper,tjbper,ticbper,titper,tjtper,tktper,tictper
      common /itp2/ tipper,tjpper,tkpper,tlpper,ticpper,tiaper,tiacper
      common /itp3/ tinegp,tjnegp,tknegp,tlnegp,t2icbper,t2ictper,
     &             t2icpper
c
      common /rt1/ tchrg,tamass,trk,treq,ttk,tteq,tpk,tpn,tphase,
     &             trborn,tfs,ttimass
      common /rt2/ tsolty,tcn1,tcn2,tasol,tbsol,thbcut
      common /rt3/ tbeta,tbox,tcutcap,txcap,tycap,tzcap
      common /rtp1/ tcgper,talmper
      
