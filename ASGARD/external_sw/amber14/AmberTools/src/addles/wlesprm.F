c
c write a modifed topology file corresponding to a LES -> noLES pert
c
      subroutine wlesprm(ounit)
#include "SIZE.h"
#include "T3TOP.h"
#include "MISC.h"
#include "TOP.h"
      integer ounit,cnum(maxnatom),count
c

10      format (20a4)
15      format (a)
20      format (12i6)
30      format (5e16.8)

c
c set nparm=1 to make rdparm read the LES info
c
       nparm=1

      write (ounit,15,err=1000) ititl
      write (ounit,20,err=1000) natom,ntypes,t3nbonh,t3mbona,t3ntheth,
     &	t3mtheta,t3nphih,t3mphia,nhparm,nparm,next,nres,t3nbona,t3ntheta,
     &	t3nphia,numbnd,numang,nptra,natyp,nphb,ifpert,t3nbper,t3ngper,
     &	t3ndper,t3mbper,t3mgper,t3mdper,ifbox,nmxrs,ifcap
c
c write atom names
      write (ounit,10,err=1000) (t3igraph(i),i=1,natom)
c charges
      write (ounit,30,err=1000) (t3chrg(i),i=1,natom)
c masses
        if (bigmas) then
        write (6,*) 'modifying masses to ',newm
          do 200 i=1,natom
      	write (6,201) i,origpt(i)
201		format ('Atom ',i5,' modified, orig. atom # ',i5)
                if (imass(origpt(i)).ne.0) t3amass(i)=newm
200       continue
        endif   
      write (ounit,30,err=1000) (t3amass(i),i=1,natom)
c L-J index
      write (ounit,20,err=1000) (t3iac(i),i=1,natom)
c exclusions
      write (ounit,20,err=1000) (numex(i),i=1,natom)
c another L-J index?
      write (ounit,20,err=1000) (ico(i),i=1,ntypes*ntypes)
c residue labels
      write (ounit,10,err=1000) (labres(i),i=1,nres)
c first atom in residue i
      write (ounit,20,err=1000) (ipres(i),i=1,nres)
c bond force constants
      write (ounit,30,err=1000) (rk(i),i=1,numbnd)
c bond equ length
      write (ounit,30,err=1000) (req(i),i=1,numbnd)
c angle force constants
      write (ounit,30,err=1000) (tk(i),i=1,numang)
c angle equ value
      write (ounit,30,err=1000) (teq(i),i=1,numang)
c dihedral force constant
      write (ounit,30,err=1000) (pk(i),i=1,nptra)
c dihedral periodicity
      write (ounit,30,err=1000) (pn(i),i=1,nptra)
c dihedral phase
      write (ounit,30,err=1000) (phase(i),i=1,nptra)
c solty unused
      write (ounit,30,err=1000) (solty(i),i=1,natyp)
c L-J r**12 terms
      write (ounit,30,err=1000) (cn1(i),i=1,ntypes*(ntypes+1)/2)
c L-J r**6 terms
      write (ounit,30,err=1000) (cn2(i),i=1,ntypes*(ntypes+1)/2)
c atoms numbers and index for bonds to hydrogen
      do 100 i=1,t3nbonh
                t3ibh(i)=(t3ibh(i)-1)*3
                t3jbh(i)=(t3jbh(i)-1)*3
100	continue
      write (ounit,20,err=1000) (t3ibh(i),t3jbh(i),t3icbh(i),
     .            i=1,t3nbonh)
c atom numbers and index for bonds not to hydrogen
      do 110 i=1,t3nbona
                t3ib(i)=(t3ib(i)-1)*3
                t3jb(i)=(t3jb(i)-1)*3
110	continue
      write (ounit,20,err=1000) (t3ib(i),t3jb(i),t3icb(i),i=1,t3nbona)
c angle info for hydrogens
      do 120 i=1,t3ntheth
                t3ith(i)=(t3ith(i)-1)*3
                t3jth(i)=(t3jth(i)-1)*3
                t3kth(i)=(t3kth(i)-1)*3
120	continue
      write (ounit,20,err=1000) (t3ith(i),t3jth(i),t3kth(i),t3icth(i),
     &		i=1,t3ntheth)
c angles
      do 130 i=1,t3ntheta
                t3it(i)=(t3it(i)-1)*3
                t3jt(i)=(t3jt(i)-1)*3
                t3kt(i)=(t3kt(i)-1)*3
130	continue
      write (ounit,20,err=1000) (t3it(i),t3jt(i),t3kt(i),t3ict(i),
     &	i=1,t3ntheta)
c dihedrals with hydrogens
      do 140 i=1,t3nphih
                t3iph(i)=(t3iph(i)-1)*3
                t3jph(i)=(t3jph(i)-1)*3
                t3kph(i)=(t3kph(i)-1)*3
                t3lph(i)=(t3lph(i)-1)*3
                if (t3inegh(i).lt.0) then
                        t3iph(i)=-(t3iph(i))
                endif
                if (t3jnegh(i).lt.0) then
                        t3jph(i)=-(t3jph(i))
                endif
                if (t3knegh(i).lt.0) then
                        t3kph(i)=-(t3kph(i))
                endif
                if (t3lnegh(i).lt.0) then
                        t3lph(i)=-(t3lph(i))
                endif
140     continue
      write (ounit,20,err=1000) (t3iph(i),t3jph(i),t3kph(i),t3lph(i),
     &	t3icph(i),i=1,t3nphih)
c dihedrals
      do 150 i=1,t3nphia
                t3ip(i)=(t3ip(i)-1)*3
                t3jp(i)=(t3jp(i)-1)*3
                t3kp(i)=(t3kp(i)-1)*3
                t3lp(i)=(t3lp(i)-1)*3
                if (t3ineg(i).lt.0) then
                        t3ip(i)=-(t3ip(i))
                endif
                if (t3jneg(i).lt.0) then
                        t3jp(i)=-(t3jp(i))
                endif
                if (t3kneg(i).lt.0) then
                        t3kp(i)=-(t3kp(i))
                endif
                if (t3lneg(i).lt.0) then
                        t3lp(i)=-(t3lp(i))
                endif
150     continue

      write (ounit,20,err=1000)(t3ip(i),t3jp(i),t3kp(i),t3lp(i),
     &     t3icp(i),i=1,t3nphia)
c excluded atoms
      write (ounit,20,err=1000) (natex(i),i=1,next)
c r**12 hbond term
      write (ounit,30,err=1000) (asol(i),i=1,nphb)
c r**10 term
      write (ounit,30,err=1000) (bsol(i),i=1,nphb)
c hbcut, not used
      write (ounit,30,err=1000) (hbcut(i),i=1,nphb)
c amber atom types
      write (ounit,10,err=1000) (t3isymbl(i),i=1,natom)
c tree info types
      write (ounit,10,err=1000) (t3itree(i),i=1,natom)
c tree joining info
      write (ounit,20,err=1000) (t3join(i),i=1,natom)
c last atom to move (not used)
      write (ounit,20,err=1000) (t3irotat(i),i=1,natom)
c only if ibox.gt.0)
      if (ifbox.gt.0) then
          write (ounit,20,err=1000) iptres, nspm,nspsol
          write (ounit,20,err=1000) (nsp(i),i=1,nspm)
       	  write (ounit,30,err=1000) beta,box(1),box(2),box(3)
      endif
c cap info
      if (ifcap.gt.0) then
          write (ounit,20,err=1000) natcap
          write (ounit,30,err=1000) cutcap, xcap, ycap, zcap
      endif
c pert info
      if (ifpert.gt.0) then
          do 2100 i=1,t3nbper
                t3ibper(i)=(t3ibper(i)-1)*3
                t3jbper(i)=(t3jbper(i)-1)*3
2100	    continue
          write (ounit,20,err=1000) (t3ibper(i),t3jbper(i),i=1,t3nbper)
          write (ounit,20,err=1000) (t3icbper(i),i=1,2*t3nbper)
          do 2110 i=1,t3ngper
                t3itper(i)=(t3itper(i)-1)*3
                t3jtper(i)=(t3jtper(i)-1)*3
                t3ktper(i)=(t3ktper(i)-1)*3
2110	    continue
          write (ounit,20,err=1000) (t3itper(i),t3jtper(i),t3ktper(i),
     &	    i=1,t3ngper)
          write (ounit,20,err=1000) (t3ictper(i),i=1,2*t3ngper)
          do 2120 i=1,t3ndper
                t3ipper(i)=(t3ipper(i)-1)*3
                t3jpper(i)=(t3jpper(i)-1)*3
                t3kpper(i)=(t3kpper(i)-1)*3
                t3lpper(i)=(t3lpper(i)-1)*3
                if (t3inegp(i).lt.0) then
                        t3ipper(i)=-(t3ipper(i))
                endif
                if (t3jnegp(i).lt.0) then
                        t3jpper(i)=-(t3jpper(i))
                endif
                if (t3knegp(i).lt.0) then
                        t3kpper(i)=-(t3kpper(i))
                endif
                if (t3lnegp(i).lt.0) then
                        t3lpper(i)=-(t3lpper(i))
                endif
2120        continue
          write (ounit,20,err=1000) (t3ipper(i),t3jpper(i),t3kpper(i),
     &	    t3lpper(i),i=1,t3ndper)
          write (ounit,20,err=1000) (t3icpper(i),i=1,2*t3ndper)
          write (ounit,10,err=1000) (labres(i),i=1,nres)
          write (ounit,10,err=1000) (t3igrper(i),i=1,natom)
          write (ounit,10,err=1000) (t3ismper(i),i=1,natom)
          write (ounit,30,err=1000) (t3almper(i),i=1,natom)
          write (ounit,20,err=1000) (t3iaper(i),i=1,natom)
          write (ounit,20,err=1000) (t3iacper(i),i=1,natom)
          write (ounit,30,err=1000) (t3cgper(i),i=1,natom)
      endif
      if (ipol.eq.1) then
      	stop
c	    write (ounit,30,err=1000) (atpol(i),i=1,natom)
      endif
      if (ipol.eq.1.and.ifpert.eq.1) then
      	stop
c	    write (ounit,30,err=1000) (atpol1(i),i=1,natom)
      endif
c
c LES info below will be read by rdparm only if nparm=1
c
      write (ounit,20,err=1000) nlestyp
      write (ounit,20,err=1000) (lestyp(i),i=1,natom)
      write (ounit,30,err=1000) (lesfac(i),i=1,nlestyp*nlestyp)
c
c now write a copy number for each atom
c
      count=1
      do 2000 i=1,natom
      	if (totcop(origpt(i)).eq.1) then
      		cnum(i)=0
      	else
      		cnum(i)=count
      		count=count+1
      	endif
      	if (count.gt.totcop(origpt(i))) count=1
2000	continue
c
c write it
c
      write (ounit,20,err=1000) (cnum(i),i=1,natom)
      
      return
c
1000	write (6,*) 'problem writing parameter file'
      stop
      end
