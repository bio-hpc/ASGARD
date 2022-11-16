

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine rdparm2 here]
subroutine rdparm2()
   
   !      this version for nmanal
   
   implicit double precision (a-h,o-z)
#  include "sizes2.h"
#  include "bad.h"
#  include "files.h"
#  include "infoa.h"
#  include "optf.h"
   
   common/consnb/npair,ntypes,idiel,iyyy,dielc,scnb,scee
   
   common/runhed/ihead(20),ihead1(20)
   common/belly/ibelly,natbel,igroup(maxatom),isymbl(maxatom), &
         itree(maxatom),igres(600)

   !     INCLUDE 'cm2.inc'
   !     INCLUDE (Forcd coordinates common block)
#  include "carts.h"
   
   !     ----- READ THE MOLECULAR TOPOLOGY -----
   
   !       PROCEDURE (Read formatted parm file)
   call amopen(20,parm,'O','F','R')
   read(20,9108) (ihead1(i),i=1,20)
   read(20,9118) natom,ntypes,nbonh,mbona,ntheth,mtheta,nphih, &
         mphia,nhparm,nparm,nnb,nres,nbona,ntheta,nphia, &
         numbnd,numang,nptra,natyp,nphb,ifpert,idum,idum, &
         idum,idum,idum,idum,idum,idum,idum,idum
   
   !     ----- CHECK FOR POSSIBLE OVERFLOWS OF ARRAY BOUNDS -----
   
   if (nbonh > maxbnh .or. nbona > maxbon .or. ntheth > maxanh &
         .or. ntheta > maxang .or. nphih > maxdih &
         .or. nphia > maxdia .or. natom > maxatom) then
      write(6,1000)
      call mexit(6, 1)
   end if
   if (numbnd > 300 .or. numang > 450 .or. nptra > 250) then
      write(6,1001)
      call mexit(6, 1)
   end if
   nbond = nbonh + nbona
   nat3 = 3*natom
   ntype = ntypes*ntypes
   
   !     ----- READ THE SYMBOLS AND THE CHARGES AND THE MASSES -----
   
   read(20,9108) (igrph(i),i = 1,natom)
   read(20,9128) (cg(i),i = 1,natom)
   read(20,9128) (amass(i),i = 1,natom)
   read(20,9118) (iac(i),i = 1,natom)
   read(20,9118) (iblo(i),i = 1,natom)
   read(20,9118) (no(i),i = 1,ntype)
   read(20,9108) (lbres(i),i=1,nres)
   read(20,9118) (ipres(i),i=1,nres)
   ipres(nres+1) = natom+1
   
   !     ----- READ THE PARAMETERS -----
   
   read(20,9128) (rk(i),    i = 1,numbnd)
   read(20,9128) (req(i),   i = 1,numbnd)
   read(20,9128) (tk(i),    i = 1,numang)
   read(20,9128) (teq(i),   i = 1,numang)
   read(20,9128) (pk(i),    i = 1,nptra)
   read(20,9128) (pn(i),    i = 1,nptra)
   read(20,9128) (phase(i), i = 1,nptra)
   read(20,9128) (solty(i), i = 1,natyp)
   
   nttyp = ntypes*(ntypes+1)/2
   
   read(20,9128) (cn1(i),   i = 1,nttyp)
   read(20,9128) (cn2(i),   i = 1,nttyp)
   
   !     ----- READ THE BONDING INFORMATIONS -----
   
   read(20,9118) (ibh(i),jbh(i),icbh(i),i = 1,nbonh)
   read(20,9118) (iba(i),jba(i),icba(i),i = 1,nbona)
   read(20,9118) (ith(i),jth(i),kth(i),icth(i), &
         i = 1,ntheth)
   read(20,9118) (ita(i),jta(i),kta(i),icta(i), &
         i = 1,ntheta)
   read(20,9118) (iph(i),jph(i),kph(i),lph(i),icph(i), &
         i = 1,nphih)
   read(20,9118) (ipa(i),jpa(i),kpa(i),lpa(i),icpa(i), &
         i = 1,nphia)
   read(20,9118) (inb(i),i=1,nnb)
   
   !     ----- READ THE H-BOND PARAMETERS -----
   
   read(20,9128) (asol(i),i=1,nphb)
   read(20,9128) (bsol(i),i=1,nphb)
   read(20,9128) (hbcut(i),i=1,nphb)
   
   !     ----- READ ISYMBL,ITREE ARRAYS -----
   
   read(20,9108) (isymbl(i),i=1,natom)
   read(20,9108) (itree(i),i=1,natom)
   
   9108 format(20a4)
   9118 format(12i6)
   9128 format(5e16.8)
   
   return
   1000 format(/5x,'****** ARRAY BOUNDS OVERFLOW IN RDPARM *****')
   1001 format(/5x,'****** need more parameter space in rdparm**')
   1002 format('In rdparm:',10i5,/10x,10i5)
end subroutine rdparm2 
