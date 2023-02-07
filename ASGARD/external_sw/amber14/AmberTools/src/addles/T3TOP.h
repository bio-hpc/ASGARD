c
c alternate topology definitions. third set.
c
      character*80 t3ititl
      character*4 t3igraph(maxnatom),t3labres(maxnres),
     .        t3isymbl(maxnatom)
      character*4 t3itree(maxnatom)
c
      integer t3natom,t3ntypes,t3nbonh,t3mbona,t3ntheth,t3mtheta
    	integer t3nphih,t3mphia,t3nhparm,t3nparm,t3nnb,t3nres,t3nbona
      integer t3ntheta,t3ifpret,t3nbper,t3nmxrs,t3ifcap
      integer t3nphia,t3numbnd,t3numang,t3nptra,t3natyp,t3nphb
      integer t3ngper,t3ndper,t3mbper,t3mgper,t3mdper,t3ifbox
c
      integer t3iac(maxnatom),t3numex(maxnatom),
     .        t3ico(maxntypes*maxntypes)
      integer t3ipres(maxnres),t3ibh(maxbnd),t3jbh(maxbnd),
     .        t3icbh(maxbnd) 
      integer t3ib(maxbnd),t3jb(maxbnd),t3icb(maxbnd) 
      integer t3ith(maxang),t3jth(maxang),t3kth(maxang),t3icth(maxang) 
      integer t3it(maxang),t3jt(maxang),t3kt(maxang),t3ict(maxang) 
      integer t3iph(maxdih),t3jph(maxdih),t3kph(maxdih),t3lph(maxdih)
      integer t3icph(maxdih),t3ip(maxdih) 
      integer t3jp(maxdih),t3kp(maxdih),t3lp(maxdih),t3icp(maxdih) 
      integer t3natex(maxnext),t3join(maxnatom),t3irotat(maxnatom)
      integer t3iptres,t3nspm,t3npsol,t3nsp(maxnspm),t3natcap
c
c pert stuff
c
        integer t3ibper(maxbnd),t3jbper(maxbnd),t3icbper(maxbnd)
      integer t32icbper(maxbnd)
        integer t3itper(maxang),t3jtper(maxang),t3ktper(maxang)
      integer t3ictper(maxang),t32ictper(maxang)
        integer t3ipper(maxdih),t3jpper(maxdih),t3kpper(maxdih)
        integer t3icpper(maxdih),t3lpper(maxdih),t32icpper(maxdih)
        integer t3iaper(maxnatom),t3iacper(maxnatom)
        integer t3inegp(maxdih),t3jnegp(maxdih),t3knegp(maxdih)
        integer t3lnegp(maxdih)
        real*8 t3cgper(maxnatom),t3almper(maxnatom)
        character*4 t3labre2(maxnres),t3igrper(maxnatom)
      character*4 t3ismper(maxnatom)
c
c here is a temporary storage for the SIGN of the torsion atoms.
c make them all positive when read, but need to replace the sign
c upon writing.
c       
        integer t3inegh(maxdih),t3jnegh(maxdih),t3knegh(maxdih)
      integer t3lnegh(maxdih),t3lneg(maxdih)
        integer t3ineg(maxdih),t3jneg(maxdih),t3kneg(maxdih)

c
      real*8 t3chrg(maxnatom),t3amass(maxnatom),t3rk(maxbnd)
      real*8 t3req(maxbnd),t3phase(maxdih),t3bsol(maxnphb),t3zcap 
      real*8 t3tk(maxang),t3teq(maxang),t3pk(maxdih),t3pn(maxdih)
      real*8 t3solty(maxnatyp),t3cn1(maxntypes*(maxntypes+1)/2)
      real*8 t3cn2(maxntypes*(maxntypes+1)/2),t3asol(maxnphb)
      real*8 t3hbcut(maxnphb),t3beta,t3box(3),t3cutcap,t3xcap,t3ycap
c
c all in common
c
      common /c3t1/ t3ititl,t3igraph,t3labres,t3isymbl,t3itree
      common /c3tp1/ t3labre2,t3igrper,t3ismper
c
      common /i3t1/ t3natom,t3ntypes,t3nbonh,t3mbona,t3ntheth,t3mtheta
    	common /i3t2/ t3nphih,t3mphia,t3nhparm,t3nparm,t3nnb,t3nres
      common /i3t2a/t3nbona,t3ntheta,t3ifpret,t3nbper,t3nmxrs,t3ifcap
      common /i3t3/ t3nphia,t3numbnd,t3numang,t3nptra,t3natyp,t3nphb
      common /i3t4/ t3ngper,t3ndper,t3mbper,t3mgper,t3mdper,t3ifbox
      common /i3t5/ t3iac,t3numex,t3ico,t3ipres,t3ibh,t3jbh,t3icbh
      common /i3t6/ t3ib,t3jb,t3icb,t3ith,t3jth,t3kth,t3icth
      common /i3t7/ t3it,t3jt,t3kt,t3ict,t3iph,t3jph,t3kph,t3lph,t3icph
      common /i3t8/ t3ip,t3jp,t3kp,t3lp,t3icp,t3natex,t3join,t3irotat
      common /i3t9/ t3iptres,t3nspm,t3npsol,t3nsp,t3natcap
      common /i3t10/ t3ineg,t3jneg,t3kneg,t3lneg,t3inegh,t3jnegh,t3knegh
      common /i3t10a/ t3lnegh,t3ictper,t3iacper,t32icpper
      common /i3tp1/ t3ibper,t3jbper,t3icbper,t3itper,t3jtper,t3ktper
        common /i3tp2/ t3ipper,t3jpper,t3kpper,t3lpper,t3icpper,t3iaper
      common /i3tp3/ t3inegp,t3jnegp,t3knegp,t3lnegp,t32icbper,t32ictper
c
      common /r3t1/ t3chrg,t3amass,t3rk,t3req,t3tk,t3teq,t3pk,t3pn
      common /r3t2/ t3solty,t3cn1,t3cn2,t3asol,t3bsol,t3hbcut,t3phase
      common /r3t3/ t3beta,t3box,t3cutcap,t3xcap,t3ycap,t3zcap
      common /r3tp1/ t3cgper,t3almper
      
      

