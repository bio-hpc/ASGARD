c
c this routine creates the NO-LES -> LES topology file for
c the original lambda=0 state (the state defined in PARM)
c
      subroutine les2prm(uw2prm)

#include "SIZE.h"
#include "TOP.h"
#include "TTOP.h"
#include "T3TOP.h"
#include "TEMP.h"
#include "MISC.h"

        logical per1,per2,per3,per4

10	format (20a4)
15	format (a)
20	format (12i6)
30	format (5e16.8)

c
c first set up atom-specific variables, for lambda=0 case
c where we are going from 1 copy to multiple copies.
c **************************************************************
c WE'LL MAKE LAMBDA=0 THE SINGLE COPY, AND LAMBDA=1 MULT COPIES.
c **************************************************************
c keep in mind this is all for the ORIGINAL lambda=0 topology
c right now the normal topo variables (TOP.h) have the LES info,
c and the TTOP.h has the non-les info
c
      do 100 i=1,natom
c
c name at orig. lam=1 even for lam=0 here
c
      	t3igraph(i)=igrper(i)
      	t3numex(i)=numex(i)
      	t3isymbl(i)=ismper(i)
      	t3itree(i)=itree(i)
      	t3join(i)=join(i)
      	t3irotat(i)=irotat(i)
      	t3igrper(i)=igrper(i)
      	t3ismper(i)=ismper(i)
      	t3almper(i)=almper(i)
c
c set iaper, 1 if atom perturbed (mult copies)
c
      	if (totcop(origpt(i)).gt.1.or.ipert(origpt(i)).eq.1) then
      		t3iaper(i)=1
      	else
      		t3iaper(i)=0
      	endif
c
c charges at lambda=0 are those for mult copies, same for all copies
c use the scaled 'original lambda=0' charges, even though lambda=1 here
c
      	t3chrg(i)=cgper(i)
      	t3iac(i)=iacper(i)
c
c make the params set to the single copy value if it is the 'lreal'
c particle, ie only 1 lreal particle per set of copies.
c
      	if (lreal(i)) then
c
c orig charge (at original lambda=0)
c
      		t3cgper(i)=ocgper(origpt(i))
c
c allow masses to change? can't?
c
      		t3amass(i)=omass(origpt(i))
      		t3iacper(i)=tiacper(i)
      	else
c
c dummy particle, no charge
c
      		t3cgper(i)=0.0d0
      		t3amass(i)=omass(origpt(i))
c
c set dummy particle type!
c
      		t3iacper(i)=dumiac
      	endif

100	continue
c
c now we need to sort the bond list into four types: non-pert, non-pert H,
c pert and boundary. nonpert have all iaper=0. pert have all iaper=1
c boundary have at least 1 of each (for angles & dihedrals too)
c
c keep in mind that 'pert' at this point refers to all LES particles,
c NOT to those actually being perturbed in the original topology
c
c how do we decide if what used to be pert bonds and are now NOT- did they 
c involve hydrogen or not? have to look at names....
c
c also need to keep the constraint bonds separate... what if they involve les?
c and therefore are now pert bonds? put in as pert constraints?
c
      nnop=0
      nnopc=0
      nnoph=0
      nboun=0
      npert=0
c
c now add to these lists every bond in the nbona, nbonh, nbper
c depending on iaper of each of the two atoms.
c
      do 110 i=1,nbonh
          les1=(totcop(origpt(ibh(i))).gt.1)
          les2=(totcop(origpt(jbh(i))).gt.1)
            per1=(ipert(origpt(ibh(i))).eq.1.or.les1)
            per2=(ipert(origpt(jbh(i))).eq.1.or.les2)
c
c bonds with hydrogen. put in nnoph, nboun or npert.
c
          if (.not.per1.and..not.per2) then
c
c doesn't involve les, leave it
c 
      	nnoph=nnoph+1
      	th1(nnoph)=ibh(i)
      	th2(nnoph)=jbh(i)
      	th3(nnoph)=icbh(i)
c 
c decide if it is all les or part les- ie all pert or boundary
c
          else if (les1.and.les2) then
c 
c all les- put in pert. but decide whether it is the 'real' copy or
c a dummy copy at lambda=1, for the lambda=1 params (lambda=0 params
c are the same for all copies)
c
      	    npert=npert+1
      	    tp1(npert)=ibh(i)
      	    tp2(npert)=jbh(i)
c
c tp3 corresponds to lambda=0, tp4 to lambda=1 !
c all copies same at lambda=1. NOT same at lambda=0!
c
c ticbh is the 'full' bond type, icbh is the 'scaled' bond type
c
c assign the value at the LES stage based on user choice of bscale
c
      	    if (bscale.eq.3) then
      	    	tp4(npert)=ticbh(i)
         	    else
      	    	tp4(npert)=icbh(i)
        	    endif
      		
c
      	    if (lreal(ibh(i)).and.lreal(jbh(i))) then
c
c all real, it is the 'only' copy of this bond at single copy end
c 1/N if bscale is 4, otherwise 1
c
      		if (bscale.eq.4) then
      	        	tp3(npert)=icbh(i)
      		else
      	        	tp3(npert)=ticbh(i)
      		endif
c
c THE PROBLEM IS: WHAT FORCE CONSTANT
C TO GIVE THE 'DUMMY' LES PARTICLES (IE WHERE LREAL = FALSE)
c in theory they should have zero force constant since they don't 'exist'
c but that won't maintain topology very well. 
c besides, all this will cancel, right? ;-)
c
c it is also possible to have one real and one not - therefore 
c not a 'real' copy even though both are LES particles.
c for example, bonds across two subspaces might have the 'real' copy for one
c space bonded to a dummy copy from the other space. only ONE out of the
c NaNb bonds should be 'real' - but see comment above. 
c
      	    else
      		if (bscale.eq.1.or.bscale.eq.4) then
      	        	tp3(npert)=icbh(i)
      		else
      	        	tp3(npert)=ticbh(i)
      		endif
      	    endif
          elseif (les1.or.les2) then
c
c here we have mixed LES/non-les, but it may be all pert.
c check to see if it will be put in pert or boundary.
c it will use scaling factors in all cases
c
                if (per1.and.per2) then
c
c pert bond - whether it gets the full force constants depends
c on whether the les particle is the 'real' one or a dummy copy, as above
c
                    npert=npert+1
                    tp1(npert)=ibh(i)
                    tp2(npert)=jbh(i)
                    if (bscale.eq.3) then
                        tp4(npert)=ticbh(i)
                    else
                        tp4(npert)=icbh(i)
                    endif
c
                    if (lreal(ibh(i)).and.lreal(jbh(i))) then
c
c all real, it is the 'only' copy of this bond at lambda =1 (except dummies)
c
                        if (bscale.eq.4) then
                                tp3(npert)=icbh(i)
                        else
                                tp3(npert)=ticbh(i)
                        endif
c
                    else
c
c the non-les MUST be 'real', so the les is a dummy.
c
                        if (bscale.eq.1.or.bscale.eq.4) then
                                tp3(npert)=icbh(i)
                        else
                                tp3(npert)=ticbh(i)
                        endif
c
                    endif
      	else
c
c boundary bond - whether it gets the full force constants depends
c on whether the les particle is the 'real' one or a dummy copy, as above
c
      	    nboun=nboun+1
      	    tb1(nboun)=ibh(i)
      	    tb2(nboun)=jbh(i)
      	    if (bscale.eq.3) then
      	    	tb4(nboun)=ticbh(i)
         	    else
      	    	tb4(nboun)=icbh(i)
        	    endif
c
      	    if (lreal(ibh(i)).and.lreal(jbh(i))) then
c
c all real, it is the 'only' copy of this bond at lambda =1 (except dummies)
c
      		if (bscale.eq.4) then
      	        	tb3(nboun)=icbh(i)
      		else
      	        	tb3(nboun)=ticbh(i)
      		endif
c
      	    else
c 
c the non-les MUST be 'real', so the les is a dummy.
c
      		if (bscale.eq.1.or.bscale.eq.4) then
      	        	tb3(nboun)=icbh(i)
      		else
      	        	tb3(nboun)=ticbh(i)
      		endif
c
      	    endif
      	endif
c
c last case is the possibility of no LES but yet some user-chosen pert atoms.
c already ruled out no pert, some les, all les so all that's left
c is no les with some pert or no les with all pert
c
            elseif (per1.and.per2) then
c
c pert but no scaling
c
                    npert=npert+1
                    tp1(npert)=ibh(i)
                    tp2(npert)=jbh(i)
                    tp4(npert)=icbh(i)
                    tp3(npert)=icbh(i)
            else
c
c boundary and no scaling
c
                    nboun=nboun+1
                    tb1(nboun)=ibh(i)
                    tb2(nboun)=jbh(i)
                    tb4(nboun)=icbh(i)
                    tb3(nboun)=icbh(i)
          endif
110	continue
c	
c now the normal bonds, but leave the constraint bonds for later
c
      do 120 i=1,mbona
          les1=(totcop(origpt(ib(i))).gt.1)
          les2=(totcop(origpt(jb(i))).gt.1)
            per1=(ipert(origpt(ib(i))).eq.1.or.les1)
            per2=(ipert(origpt(jb(i))).eq.1.or.les2)
          if (.not.per1.and..not.per2) then
      	nnop=nnop+1
      	tn1(nnop)=ib(i)
      	tn2(nnop)=jb(i)
      	tn3(nnop)=icb(i)
          else if (les1.and.les2) then
      	    npert=npert+1
      	    tp1(npert)=ib(i)
      	    tp2(npert)=jb(i)
      	    if (bscale.eq.3) then
      	    	tp4(npert)=ticb(i)
         	    else
      	    	tp4(npert)=icb(i)
        	    endif
      	    if (lreal(ib(i)).and.lreal(jb(i))) then
      		if (bscale.eq.4) then
      	        	tp3(npert)=icb(i)
      		else
      	        	tp3(npert)=ticb(i)
      		endif
      	    else
      		if (bscale.eq.1.or.bscale.eq.4) then
      	        	tp3(npert)=icb(i)
      		else
      	        	tp3(npert)=ticb(i)
      		endif
      	    endif
          elseif (les1.or.les2) then
                if (per1.and.per2) then
                    npert=npert+1
                    tp1(npert)=ib(i)
                    tp2(npert)=jb(i)
                    if (bscale.eq.3) then
                        tp4(npert)=ticb(i)
                    else
                        tp4(npert)=icb(i)
                    endif
                    if (lreal(ib(i)).and.lreal(jb(i))) then
                        if (bscale.eq.4) then
                                tp3(npert)=icb(i)
                        else
                                tp3(npert)=ticb(i)
                        endif
                    else
                        if (bscale.eq.1.or.bscale.eq.4) then
                                tp3(npert)=icb(i)
                        else
                                tp3(npert)=ticb(i)
                        endif
                    endif
      	else
      	    nboun=nboun+1
      	    tb1(nboun)=ib(i)
      	    tb2(nboun)=jb(i)
      	    if (bscale.eq.3) then
      	    	tb4(nboun)=ticb(i)
         	    else
      	    	tb4(nboun)=icb(i)
        	    endif
      	    if (lreal(ib(i)).and.lreal(jb(i))) then
      		if (bscale.eq.4) then
      	        	tb3(nboun)=icb(i)
      		else
      	        	tb3(nboun)=ticb(i)
      		endif
      	    else
      		if (bscale.eq.1.or.bscale.eq.4) then
      	        	tb3(nboun)=icb(i)
      		else
      	        	tb3(nboun)=ticb(i)
      		endif
      	    endif
      	endif
c
            elseif (per1.and.per2) then
                    npert=npert+1
                    tp1(npert)=ib(i)
                    tp2(npert)=jb(i)
                    tp4(npert)=icb(i)
                    tp3(npert)=icb(i)
            else
                    nboun=nboun+1
                    tb1(nboun)=ib(i)
                    tb2(nboun)=jb(i)
                    tb4(nboun)=icb(i)
                    tb3(nboun)=icb(i)
          endif
120	continue
c	
c now the pert bonds (what were pert in the original topology)
c this part differes from the OTHER les perturbation in that
c the original lambda=0 params are used for all lambda here (some scaled)
c
      do 130 i=1,nbper
          les1=(totcop(origpt(ibper(i))).gt.1)
          les2=(totcop(origpt(jbper(i))).gt.1)
            per1=(ipert(origpt(ibper(i))).eq.1.or.les1)
            per2=(ipert(origpt(jbper(i))).eq.1.or.les2)
          if (.not.per1.and..not.per2) then
c
c NOT PERT ANYMORE! decide where to put it based on whether either name 
c  begins with H or DH - and use the atom names at lambda = 1!
c ******* THIS WILL AFFECT WHAT TO SHAKE, NOTHING ELSE**************
c
      	if (isymbl(ibper(i))(1:1).eq.'H'.or.
     &		isymbl(jbper(i))(1:1).eq.'H'.or.
     &		isymbl(ibper(i))(1:2).eq.'DH'.or.
     &		isymbl(jbper(i))(1:2).eq.'DH') then
c
      		nnoph=nnoph+1
      		th1(nnoph)=ibper(i)
      		th2(nnoph)=jbper(i)
c use lambda=0 type
      		th3(nnoph)=ticbper(i+nbper)
              else
c no H
      		nnop=nnop+1
      		tn1(nnop)=ibper(i)
      		tn2(nnop)=jbper(i)
      		tn3(nnop)=ticbper(i+nbper)
c
      	endif
          elseif (les1.and.les2) then
c
c ok, still a pert
c
      	    npert=npert+1
      	    tp1(npert)=ibper(i)
      	    tp2(npert)=jbper(i)
      	    if (bscale.eq.3) then
      	    	tp4(npert)=ticbper(i+nbper)
      	    else
      	    	tp4(npert)=icbper(i+nbper)
      	    endif
      	    if (lreal(ibper(i)).and.lreal(jbper(i))) then
c
c use the original lambda=0 type for both lambda=0 and lambda=1
c use the lambda=1 params for the OTHER LES pert topo file.
c
      		if (bscale.eq.4) then
      	        	tp3(npert)=icbper(i+nbper)
      		else
      	        	tp3(npert)=ticbper(i+nbper)
      		endif
      	    else
      		if (bscale.eq.1.or.bscale.eq.4) then
      	        	tp3(npert)=icbper(i+nbper)
      		else
      	        	tp3(npert)=ticbper(i+nbper)
      		endif
      	    endif
          elseif (les1.or.les2) then
                if (per1.and.per2) then
                    npert=npert+1
                    tp1(npert)=ibper(i)
                    tp2(npert)=jbper(i)
                    if (bscale.eq.3) then
                        tp4(npert)=ticbper(i+nbper)
                    else
                        tp4(npert)=icbper(i+nbper)
                    endif
                    if (lreal(ibper(i)).and.lreal(jbper(i))) then
                        if (bscale.eq.4) then
                                tp3(npert)=icbper(i+nbper)
                        else
                                tp3(npert)=ticbper(i+nbper)
                        endif
                    else
                        if (bscale.eq.1.or.bscale.eq.4) then
                                tp3(npert)=icbper(i+nbper)
                        else
                                tp3(npert)=ticbper(i+nbper)
                        endif
                    endif
                else
c boundary
      	    nboun=nboun+1
      	    tb1(nboun)=ibper(i)
      	    tb2(nboun)=jbper(i)
      	    if (bscale.eq.3) then
      	    	tb4(nboun)=ticbper(i+nbper)
      	    else
      	    	tb4(nboun)=icbper(i+nbper)
      	    endif
      	    if (lreal(ibper(i)).and.lreal(jbper(i))) then
      		if (bscale.eq.4) then
      	        	tb3(nboun)=icbper(i+nbper)
      		else
      	        	tb3(nboun)=ticbper(i+nbper)
      		endif
      	    else
      		if (bscale.eq.1.or.bscale.eq.4) then
      	        	tb3(nboun)=icbper(i+nbper)
      		else
      	        	tb3(nboun)=ticbper(i+nbper)
      		endif
      	    endif
      	endif
            elseif (per1.and.per2) then
                    npert=npert+1
                    tp1(npert)=ibper(i)
                    tp2(npert)=jbper(i)
                    tp4(npert)=icbper(i+nbper)
                    tp3(npert)=icbper(i+nbper)
            else
                    nboun=nboun+1
                    tb1(nboun)=ibper(i)
                    tb2(nboun)=jbper(i)
                    tb4(nboun)=icbper(i+nbper)
                    tb3(nboun)=icbper(i+nbper)
          endif
130	continue
c
c now do constraint bonds- but what if LES is involved?
c how do you do a perturbed constraint bond? what are they?
c
c ok to add at the end of the bond list now. non-constraint bonds done.
c
c set the number of current bonds in the new lists that aren't constraints
c
      nnopc=nnop
c
      do 140 i=mbona+1,nbona
          les1=(totcop(origpt(ib(i))).gt.1)
          les2=(totcop(origpt(jb(i))).gt.1)
            per1=(ipert(origpt(ib(i))).eq.1.or.les1)
            per2=(ipert(origpt(jb(i))).eq.1.or.les2)
          if (.not.per1.and..not.per2) then
      	nnop=nnop+1
      	tn1(nnop)=ib(i)
      	tn2(nnop)=jb(i)
      	tn3(nnop)=icb(i)
          else if (les1.and.les2) then
c
c pert a constraint with LES?
c
      	write (6,*) 'Warning- using LES for a constraint- untested'
      	    npert=npert+1
      	    tp1(npert)=ib(i)
      	    tp2(npert)=jb(i)
      	    if (bscale.eq.3) then
      	    	tp4(npert)=ticb(i)
      	    else
      	    	tp4(npert)=icb(i)
      	    endif
      	    if (lreal(ib(i)).and.lreal(jb(i))) then
      		if (bscale.eq.4) then
      	        	tp3(npert)=icb(i)
      		else
      	        	tp3(npert)=ticb(i)
      		endif
      	    else
      		if (bscale.eq.4.or.bscale.eq.1) then
      	        	tp3(npert)=icb(i)
      		else
      	        	tp3(npert)=ticb(i)
      		endif
      	    endif
          elseif (les1.or.les2) then
                write (6,*) 
     &		 'Warning- using LES for a constraint- untested'
                if (per1.and.per2) then
                    npert=npert+1
                    tp1(npert)=ib(i)
                    tp2(npert)=jb(i)
                    if (bscale.eq.3) then
                        tp4(npert)=ticb(i)
                    else
                        tp4(npert)=icb(i)
                    endif
                    if (lreal(ib(i)).and.lreal(jb(i))) then
                        if (bscale.eq.4) then
                                tp3(npert)=icb(i)
                        else
                                tp3(npert)=ticb(i)
                        endif
                    else
                        if (bscale.eq.1.or.bscale.eq.4) then
                                tp3(npert)=icb(i)
                        else
                                tp3(npert)=ticb(i)
                        endif
                    endif
                else
      	write (6,*) 'Warning- using LES for a constraint- untested'
      	    nboun=nboun+1
      	    tb1(nboun)=ib(i)
      	    tb2(nboun)=jb(i)
      	    if (bscale.eq.3) then
      	    	tb4(nboun)=ticb(i)
      	    else
      	    	tb4(nboun)=icb(i)
      	    endif
      	    if (lreal(ib(i)).and.lreal(jb(i))) then
      		if (bscale.eq.4) then
      	        	tb3(nboun)=icb(i)
      		else
      	        	tb3(nboun)=ticb(i)
      		endif
      	    else
      		if (bscale.eq.1.or.bscale.eq.4) then
      	        	tb3(nboun)=icb(i)
      		else
      	        	tb3(nboun)=ticb(i)
      		endif
      	    endif
      	endif
            elseif (per1.and.per2) then
                    npert=npert+1
                    tp1(npert)=ib(i)
                    tp2(npert)=jb(i)
                    tp4(npert)=icb(i)
                    tp3(npert)=icb(i)
            else
                    nboun=nboun+1
                    tb1(nboun)=ib(i)
                    tb2(nboun)=jb(i)
                    tb4(nboun)=icb(i)
                    tb3(nboun)=icb(i)
          endif
140	continue
c
c now move all this to the t3 arrays to free up the temp ones for the angles
c
      do 150 i=1,nnop
      	t3ib(i)=tn1(i)
      	t3jb(i)=tn2(i)
      	t3icb(i)=tn3(i)
150	continue
      do 155 i=1,nnoph
      	t3ibh(i)=th1(i)
      	t3jbh(i)=th2(i)
      	t3icbh(i)=th3(i)
155	continue
      do 160 i=1,npert
      	t3ibper(i)=tp1(i)
      	t3jbper(i)=tp2(i)
      	t3icbper(i)=tp4(i)
      	t3icbper(i+npert+nboun)=tp3(i)
160	continue
      do 165 i=1,nboun
      	j=i+npert
      	t3ibper(j)=tb1(i)
      	t3jbper(j)=tb2(i)
      	t3icbper(j)=tb4(i)
      	t3icbper(j+npert+nboun)=tb3(i)
165	continue
      t3nbonh=nnoph
      t3nbona=nnop
      t3mbona=nnopc
      t3nbper=npert+nboun
      t3mbper=npert
c
c now do the angles essentially the same
c
      nnop=0
      nnopc=0
      nnoph=0
      nboun=0
      npert=0
c
      do 200 i=1,ntheth
          les1=(totcop(origpt(ith(i))).gt.1)
          les2=(totcop(origpt(jth(i))).gt.1)
          les3=(totcop(origpt(kth(i))).gt.1)
            per1=(ipert(origpt(ith(i))).eq.1.or.les1)
            per2=(ipert(origpt(jth(i))).eq.1.or.les2)
            per3=(ipert(origpt(kth(i))).eq.1.or.les3)
          if (.not.per1.and..not.per2.and..not.per3) then
      	nnoph=nnoph+1
      	th1(nnoph)=ith(i)
      	th2(nnoph)=jth(i)
      	th3(nnoph)=kth(i)
      	th4(nnoph)=icth(i)
          else if (les1.and.les2.and.les3) then
      	    npert=npert+1
      	    tp1(npert)=ith(i)
      	    tp2(npert)=jth(i)
      	    tp3(npert)=kth(i)
      	    if (ascale.eq.3) then
      	    	tp5(npert)=ticth(i)
      	    else
      	    	tp5(npert)=icth(i)
      	    endif
      	    if ((lreal(ith(i)).and.lreal(jth(i)).and.
     &		    lreal(kth(i)))) then
      		if (ascale.eq.4) then
      	        	tp4(npert)=icth(i)
      		else
      	        	tp4(npert)=ticth(i)
      		endif
      	    else
      		if (ascale.eq.1.or.ascale.eq.4) then
      	        	tp4(npert)=icth(i)
      		else
      	        	tp4(npert)=ticth(i)
      		endif
      	    endif
          elseif (les1.or.les2.or.les3) then
                if (per1.and.per2.and.per3) then
c pert
                    npert=npert+1
                    tp1(npert)=ith(i)
                    tp2(npert)=jth(i)
                    tp3(npert)=kth(i)
                    if (ascale.eq.3) then
                        tp5(npert)=ticth(i)
                    else
                        tp5(npert)=icth(i)
                    endif
                    if ((lreal(ith(i)).and.lreal(jth(i)).and.
     &              lreal(kth(i)))) then
                        if (ascale.eq.4) then
                                tp4(npert)=icth(i)
                        else
                                tp4(npert)=ticth(i)
                        endif
                    else
                        if (ascale.eq.1.or.ascale.eq.4) then
                                tp4(npert)=icth(i)
                        else
                                tp4(npert)=ticth(i)
                        endif
                    endif
                else
      	    nboun=nboun+1
      	    tb1(nboun)=ith(i)
      	    tb2(nboun)=jth(i)
      	    tb3(nboun)=kth(i)
      	    if (ascale.eq.3) then
      	    	tb5(nboun)=ticth(i)
      	    else
      	    	tb5(nboun)=icth(i)
      	    endif
      	    if ((lreal(ith(i)).and.lreal(jth(i)).and.
     &		    lreal(kth(i)))) then
      		if (ascale.eq.4) then
      	        	tb4(nboun)=icth(i)
      		else
      	        	tb4(nboun)=ticth(i)
      		endif
      	    else
      		if (ascale.eq.1.or.ascale.eq.4) then
      	        	tb4(nboun)=icth(i)
      		else
      	        	tb4(nboun)=ticth(i)
      		endif
      	    endif
      	endif
           elseif (per1.and.per2.and.per3) then
c all pert no les
                    npert=npert+1
                    tp1(npert)=ith(i)
                    tp2(npert)=jth(i)
                    tp3(npert)=kth(i)
                    tp4(npert)=icth(i)
                    tp5(npert)=icth(i)
           else
c some pert no les
                    nboun=nboun+1
                    tb1(nboun)=ith(i)
                    tb2(nboun)=jth(i)
                    tb3(nboun)=kth(i)
                    tb4(nboun)=icth(i)
                    tb5(nboun)=icth(i)
          endif
200	continue
      do 210 i=1,mtheta
          les1=(totcop(origpt(it(i))).gt.1)
          les2=(totcop(origpt(jt(i))).gt.1)
          les3=(totcop(origpt(kt(i))).gt.1)
            per1=(ipert(origpt(it(i))).eq.1.or.les1)
            per2=(ipert(origpt(jt(i))).eq.1.or.les2)
            per3=(ipert(origpt(kt(i))).eq.1.or.les3)
          if (.not.per1.and..not.per2.and..not.per3) then
      	nnop=nnop+1
      	tn1(nnop)=it(i)
      	tn2(nnop)=jt(i)
      	tn3(nnop)=kt(i)
      	tn4(nnop)=ict(i)
          else if (les1.and.les2.and.les3) then
      	    npert=npert+1
      	    tp1(npert)=it(i)
      	    tp2(npert)=jt(i)
      	    tp3(npert)=kt(i)
      	    if (ascale.eq.3) then
      	    	tp5(npert)=tict(i)
      	    else
      	    	tp5(npert)=ict(i)
      	    endif
      	    if ((lreal(it(i)).and.lreal(jt(i)).and.
     &		    lreal(kt(i)))) then
      		if (ascale.eq.4) then
      	        	tp4(npert)=ict(i)
      		else
      	        	tp4(npert)=tict(i)
      		endif
      	    else
      		if (ascale.eq.1.or.ascale.eq.4) then
      	        	tp4(npert)=ict(i)
      		else
      	        	tp4(npert)=tict(i)
      		endif
      	    endif
          elseif (les1.or.les2.or.les3) then
                if (per1.and.per2.and.per3) then
                    npert=npert+1
                    tp1(npert)=it(i)
                    tp2(npert)=jt(i)
                    tp3(npert)=kt(i)
                    if (ascale.eq.3) then
                        tp5(npert)=tict(i)
                    else
                        tp5(npert)=ict(i)
                    endif
                    if ((lreal(it(i)).and.lreal(jt(i)).and.
     &              lreal(kt(i)))) then
                        if (ascale.eq.4) then
                                tp4(npert)=ict(i)
                        else
                                tp4(npert)=tict(i)
                        endif
                    else
                        if (ascale.eq.1.or.ascale.eq.4) then
                                tp4(npert)=ict(i)
                        else
                                tp4(npert)=tict(i)
                        endif
                    endif
                else
      	    nboun=nboun+1
      	    tb1(nboun)=it(i)
      	    tb2(nboun)=jt(i)
      	    tb3(nboun)=kt(i)
      	    if (ascale.eq.3) then
      		tb5(nboun)=tict(i)
      	    else
      		tb5(nboun)=ict(i)
      	    endif
      	    if ((lreal(it(i)).and.lreal(jt(i)).and.
     &		    lreal(kt(i)))) then
      		if (ascale.eq.4) then
      	        	tb4(nboun)=ict(i)
      		else
      	        	tb4(nboun)=tict(i)
      		endif
      	    else
      		if (ascale.eq.1.or.ascale.eq.4) then
      	        	tb4(nboun)=ict(i)
      		else
      	        	tb4(nboun)=tict(i)
      		endif
      	    endif
      	endif
           elseif (per1.and.per2.and.per3) then
c all pert no les
                    npert=npert+1
                    tp1(npert)=it(i)
                    tp2(npert)=jt(i)
                    tp3(npert)=kt(i)
                    tp4(npert)=ict(i)
                    tp5(npert)=ict(i)
           else
c some pert no les
                    nboun=nboun+1
                    tb1(nboun)=it(i)
                    tb2(nboun)=jt(i)
                    tb3(nboun)=kt(i)
                    tb4(nboun)=ict(i)
                    tb5(nboun)=ict(i)
          endif
210	continue
      do 220 i=1,ngper
          les1=(totcop(origpt(itper(i))).gt.1)
          les2=(totcop(origpt(jtper(i))).gt.1)
          les3=(totcop(origpt(ktper(i))).gt.1)
            per1=(ipert(origpt(itper(i))).eq.1.or.les1)
            per2=(ipert(origpt(jtper(i))).eq.1.or.les2)
            per3=(ipert(origpt(ktper(i))).eq.1.or.les3)
          if (.not.per1.and..not.per2.and..not.per3) then
      	if (isymbl(itper(i))(1:1).eq.'H'.or.
     &		isymbl(jtper(i))(1:1).eq.'H'.or.
     &		isymbl(ktper(i))(1:1).eq.'H'.or.
     &		isymbl(itper(i))(1:2).eq.'DH'.or.
     &		isymbl(jtper(i))(1:2).eq.'DH'.or.
     &		isymbl(ktper(i))(1:2).eq.'DH') then
      		nnoph=nnoph+1
      		th1(nnoph)=itper(i)
      		th2(nnoph)=jtper(i)
      		th3(nnoph)=ktper(i)
      		th4(nnoph)=ictper(i+ngper)
              else
      		nnop=nnop+1
      		tn1(nnop)=itper(i)
      		tn2(nnop)=jtper(i)
      		tn3(nnop)=ktper(i)
      		tn4(nnop)=ictper(i+ngper)
      	endif
          elseif (les1.and.les2.and.les3) then
      	    npert=npert+1
      	    tp1(npert)=itper(i)
      	    tp2(npert)=jtper(i)
      	    tp3(npert)=ktper(i)
      	    if (ascale.eq.3) then
      		tp5(npert)=tictper(i+ngper)
      	    else
      		tp5(npert)=ictper(i+ngper)
      	    endif
      	    if ((lreal(itper(i)).and.lreal(jtper(i)).and.
     &		    lreal(ktper(i)))) then
      		if (ascale.eq.4) then
      	        	tp4(npert)=ictper(i+ngper)
      		else
      	        	tp4(npert)=tictper(i+ngper)
      		endif
      	    else
      		if (ascale.eq.1.or.ascale.eq.4) then
      	        	tp4(npert)=ictper(i+ngper)
      		else
      	        	tp4(npert)=tictper(i+ngper)
      		endif
      	    endif
          elseif (les1.or.les2.or.les3) then
                if (per1.and.per2.and.per3) then
                    npert=npert+1
                    tp1(npert)=itper(i)
                    tp2(npert)=jtper(i)
                    tp3(npert)=ktper(i)
                    if (ascale.eq.3) then
                        tp5(npert)=tictper(i+ngper)
                    else
                        tp5(npert)=ictper(i+ngper)
                    endif
                    if ((lreal(itper(i)).and.lreal(jtper(i)).and.
     &              lreal(ktper(i)))) then
                        if (ascale.eq.4) then
                                tp4(npert)=ictper(i+ngper)
                        else
                                tp4(npert)=tictper(i+ngper)
                        endif
                    else
                        if (ascale.eq.1.or.ascale.eq.4) then
                                tp4(npert)=ictper(i+ngper)
                        else
                                tp4(npert)=tictper(i+ngper)
                        endif
                    endif
                else
      	    nboun=nboun+1
      	    tb1(nboun)=itper(i)
      	    tb2(nboun)=jtper(i)
      	    tb3(nboun)=ktper(i)
      	    if (ascale.eq.3) then
      		tb5(nboun)=tictper(i+ngper)
      	    else
      		tb5(nboun)=ictper(i+ngper)
      	    endif
      	    if ((lreal(itper(i)).and.lreal(jtper(i)).and.
     &		    lreal(ktper(i)))) then
      		if (ascale.eq.4) then
      	        	tb4(nboun)=ictper(i+ngper)
      		else
      	        	tb4(nboun)=tictper(i+ngper)
      		endif
      	    else
      		if (ascale.eq.1.or.ascale.eq.4) then
      	        	tb4(nboun)=ictper(i+ngper)
      		else
      	        	tb4(nboun)=tictper(i+ngper)
      		endif
      	    endif
      	endif
           elseif (per1.and.per2.and.per3) then
c all pert no les
                    npert=npert+1
                    tp1(npert)=itper(i)
                    tp2(npert)=jtper(i)
                    tp3(npert)=ktper(i)
                    tp4(npert)=ictper(i+ngper)
                    tp5(npert)=ictper(i+ngper)
           else
c some pert no les
                    nboun=nboun+1
                    tb1(nboun)=itper(i)
                    tb2(nboun)=jtper(i)
                    tb3(nboun)=ktper(i)
                    tb4(nboun)=ictper(i+ngper)
                    tb5(nboun)=ictper(i+ngper)
          endif
220	continue
      nnopc=nnop
      do 230 i=mtheta+1,ntheta
          les1=(totcop(origpt(it(i))).gt.1)
          les2=(totcop(origpt(jt(i))).gt.1)
          les3=(totcop(origpt(kt(i))).gt.1)
            per1=(ipert(origpt(it(i))).eq.1.or.les1)
            per2=(ipert(origpt(jt(i))).eq.1.or.les2)
            per3=(ipert(origpt(kt(i))).eq.1.or.les3)
          if (.not.per1.and..not.per2.and..not.per3) then
      	nnop=nnop+1
      	tn1(nnop)=it(i)
      	tn2(nnop)=jt(i)
      	tn3(nnop)=kt(i)
      	tn4(nnop)=ict(i)
          else if (les1.and.les2.and.les3) then
      	write (6,*) 'Warning- using LES for a constraint- untested'
      	    npert=npert+1
      	    tp1(npert)=it(i)
      	    tp2(npert)=jt(i)
      	    tp3(npert)=kt(i)
      	    if (ascale.eq.3) then
      		tp5(npert)=tict(i)
      	    else
      		tp5(npert)=ict(i)
      	    endif
      	    if ((lreal(it(i)).and.lreal(jt(i)).and.
     &		    lreal(kt(i)))) then
      		if (ascale.eq.4) then
      	        	tp4(npert)=ict(i)
      		else
      	        	tp4(npert)=tict(i)
      		endif
      	    else
      		if (ascale.eq.1.or.ascale.eq.4) then
      	        	tp4(npert)=ict(i)
      		else
      	        	tp4(npert)=tict(i)
      		endif
      	    endif
          elseif (les1.or.les2.or.les3) then
                if (per1.and.per2.and.per3) then
                    write (6,*) 
     &		    'Warning- using LES for a constraint- untested'
                    npert=npert+1
                    tp1(npert)=it(i)
                    tp2(npert)=jt(i)
                    tp3(npert)=kt(i)
                    if (ascale.eq.3) then
                        tp5(npert)=tict(i)
                    else
                        tp5(npert)=ict(i)
                    endif
                    if ((lreal(it(i)).and.lreal(jt(i)).and.
     &              lreal(kt(i)))) then
                        if (ascale.eq.4) then
                                tp4(npert)=ict(i)
                        else
                                tp4(npert)=tict(i)
                        endif
                    else
                        if (ascale.eq.1.or.ascale.eq.4) then
                                tp4(npert)=ict(i)
                        else
                                tp4(npert)=tict(i)
                        endif
                    endif
                else
                    write (6,*) 
     &		    'Warning- using LES for a constraint- untested'
                    nboun=nboun+1
                    tb1(nboun)=it(i)
                    tb2(nboun)=jt(i)
                    tb3(nboun)=kt(i)
                    if (ascale.eq.3) then
                        tb5(nboun)=tict(i)
                    else
                        tb5(nboun)=ict(i)
                    endif
                    if ((lreal(it(i)).and.lreal(jt(i)).and.
     &              lreal(kt(i)))) then
                        if (ascale.eq.4) then
                                tb4(nboun)=ict(i)
                        else
                                tb4(nboun)=tict(i)
                        endif
                    else
                        if (ascale.eq.1.or.ascale.eq.4) then
                                tb4(nboun)=ict(i)
                        else
                                tb4(nboun)=tict(i)
                        endif
                    endif
      	endif
           elseif (per1.and.per2.and.per3) then
c all pert no les
                    npert=npert+1
                    tp1(npert)=it(i)
                    tp2(npert)=jt(i)
                    tp3(npert)=kt(i)
                    tp4(npert)=ict(i)
                    tp5(npert)=ict(i)
           else
c some pert no les
                    nboun=nboun+1
                    tb1(nboun)=it(i)
                    tb2(nboun)=jt(i)
                    tb3(nboun)=kt(i)
                    tb4(nboun)=ict(i)
                    tb5(nboun)=ict(i)
          endif
230	continue
c
c now move all this to the t3 arrays to free up the temp ones for the dihedrals
c
      do 240 i=1,nnop
      	t3it(i)=tn1(i)
      	t3jt(i)=tn2(i)
      	t3kt(i)=tn3(i)
      	t3ict(i)=tn4(i)
240	continue
      do 245 i=1,nnoph
      	t3ith(i)=th1(i)
      	t3jth(i)=th2(i)
      	t3kth(i)=th3(i)
      	t3icth(i)=th4(i)
245	continue
      do 250 i=1,npert
      	t3itper(i)=tp1(i)
      	t3jtper(i)=tp2(i)
      	t3ktper(i)=tp3(i)
      	t3ictper(i)=tp5(i)
      	t3ictper(i+npert+nboun)=tp4(i)
250	continue
      do 255 i=1,nboun
      	j=i+npert
      	t3itper(j)=tb1(i)
      	t3jtper(j)=tb2(i)
      	t3ktper(j)=tb3(i)
      	t3ictper(j)=tb5(i)
      	t3ictper(j+npert+nboun)=tb4(i)
255	continue
      t3ntheth=nnoph
      t3ntheta=nnop
      t3mtheta=nnopc
      t3ngper=npert+nboun
      t3mgper=npert
c
c now dihedrals
c
      nnop=0
      nnopc=0
      nnoph=0
      nboun=0
      npert=0
c
      do 260 i=1,nphih
          les1=(totcop(origpt(iph(i))).gt.1)
          les2=(totcop(origpt(jph(i))).gt.1)
          les3=(totcop(origpt(kph(i))).gt.1)
          les4=(totcop(origpt(lph(i))).gt.1)
            per1=(ipert(origpt(iph(i))).eq.1.or.les1)
            per2=(ipert(origpt(jph(i))).eq.1.or.les2)
            per3=(ipert(origpt(kph(i))).eq.1.or.les3)
            per4=(ipert(origpt(lph(i))).eq.1.or.les4)
          if (.not.per1.and..not.per2.and..not.per3.and.
     &	    .not.per4) then
      	nnoph=nnoph+1
      	th1(nnoph)=iph(i)
      	th2(nnoph)=jph(i)
      	th3(nnoph)=kph(i)
      	th4(nnoph)=lph(i)
      	th5(nnoph)=icph(i)
c need the ineg info too...
      	th7(nnoph)=inegh(i)
      	th8(nnoph)=jnegh(i)
      	th9(nnoph)=knegh(i)
      	th10(nnoph)=lnegh(i)
          else if (les1.and.les2.and.les3.and.les4) then
      	    npert=npert+1
      	    tp1(npert)=iph(i)
      	    tp2(npert)=jph(i)
      	    tp3(npert)=kph(i)
      	    tp4(npert)=lph(i)
      	    if (dscale.eq.3) then
      		tp6(npert)=ticph(i)
      	    else
      		tp6(npert)=icph(i)
      	    endif
      	    if ((lreal(iph(i)).and.lreal(jph(i)).and.
     &		    lreal(kph(i)).and.lreal(lph(i)))) then
      		if (dscale.eq.4) then
      	        	tp5(npert)=icph(i)
      		else
      	        	tp5(npert)=ticph(i)
      		endif
      	    else
c
c set it to no barrier - free rotation  
c no- that converges poorly. leave it and hope it cancels.
c
c		        tp5(npert)=ticph(i)+realtors
      		if (dscale.eq.1.or.dscale.eq.4) then
      	        	tp5(npert)=icph(i)
      		else
      	        	tp5(npert)=ticph(i)
      		endif
      	    endif
      	    tp7(npert)=inegh(i)
      	    tp8(npert)=jnegh(i)
      	    tp9(npert)=knegh(i)
      	    tp10(npert)=lnegh(i)
          elseif (les1.or.les2.or.les3.or.les4) then
                if (per1.and.per2.and.per3.and.per4) then
                    npert=npert+1
                    tp1(npert)=iph(i)
                    tp2(npert)=jph(i)
                    tp3(npert)=kph(i)
                    tp4(npert)=lph(i)
                    if (dscale.eq.3) then
                        tp6(npert)=ticph(i)
                    else
                        tp6(npert)=icph(i)
                    endif
                    if ((lreal(iph(i)).and.lreal(jph(i)).and.
     &              lreal(kph(i)).and.lreal(lph(i)))) then
                        if (dscale.eq.4) then
                                tp5(npert)=icph(i)
                        else
                                tp5(npert)=ticph(i)
                        endif
                    else
c                       tp5(npert)=ticph(i)+realtors
                        if (dscale.eq.1.or.dscale.eq.4) then
                                tp5(npert)=icph(i)
                        else
                                tp5(npert)=ticph(i)
                        endif
                    endif
                    tp7(npert)=inegh(i)
                    tp8(npert)=jnegh(i)
                    tp9(npert)=knegh(i)
                    tp10(npert)=lnegh(i)
                else
      	    nboun=nboun+1
      	    tb1(nboun)=iph(i)
      	    tb2(nboun)=jph(i)
      	    tb3(nboun)=kph(i)
      	    tb4(nboun)=lph(i)
      	    if (dscale.eq.3) then
      		tb6(nboun)=ticph(i)
      	    else
      		tb6(nboun)=icph(i)
      	    endif
      	    if ((lreal(iph(i)).and.lreal(jph(i)).and.
     &		    lreal(kph(i)).and.lreal(lph(i)))) then
      		if (dscale.eq.4) then
      	        	tb5(nboun)=icph(i)
      		else
      	        	tb5(nboun)=ticph(i)
      		endif
      	    else
c		        tb5(nboun)=ticph(i)+realtors
      		if (dscale.eq.1.or.dscale.eq.4) then
      	        	tb5(nboun)=icph(i)
      		else
      	        	tb5(nboun)=ticph(i)
      		endif
      	    endif
      	    tb7(nboun)=inegh(i)
      	    tb8(nboun)=jnegh(i)
      	    tb9(nboun)=knegh(i)
      	    tb10(nboun)=lnegh(i)
      	endif
            elseif (per1.and.per2.and.per3.and.per4) then
                    npert=npert+1
                    tp1(npert)=iph(i)
                    tp2(npert)=jph(i)
                    tp3(npert)=kph(i)
                    tp4(npert)=lph(i)
                    tp5(npert)=icph(i)
                    tp6(npert)=icph(i)
                    tp7(npert)=inegh(i)
                    tp8(npert)=jnegh(i)
                    tp9(npert)=knegh(i)
                    tp10(npert)=lnegh(i)
            else
                    nboun=nboun+1
                    tb1(nboun)=iph(i)
                    tb2(nboun)=jph(i)
                    tb3(nboun)=kph(i)
                    tb4(nboun)=lph(i)
                    tb5(nboun)=icph(i)
                    tb6(nboun)=icph(i)
                    tb7(nboun)=inegh(i)
                    tb8(nboun)=jnegh(i)
                    tb9(nboun)=knegh(i)
                    tb10(nboun)=lnegh(i)
          endif
260	continue
      do 270 i=1,mphia
          les1=(totcop(origpt(ip(i))).gt.1)
          les2=(totcop(origpt(jp(i))).gt.1)
          les3=(totcop(origpt(kp(i))).gt.1)
          les4=(totcop(origpt(lp(i))).gt.1)
            per1=(ipert(origpt(ip(i))).eq.1.or.les1)
            per2=(ipert(origpt(jp(i))).eq.1.or.les2)
            per3=(ipert(origpt(kp(i))).eq.1.or.les3)
            per4=(ipert(origpt(lp(i))).eq.1.or.les4)
          if (.not.per1.and..not.per2.and..not.per3.and.
     &	    .not.per4) then
      	nnop=nnop+1
      	tn1(nnop)=ip(i)
      	tn2(nnop)=jp(i)
      	tn3(nnop)=kp(i)
      	tn4(nnop)=lp(i)
      	tn5(nnop)=icp(i)
      	tn7(nnop)=ineg(i)
      	tn8(nnop)=jneg(i)
      	tn9(nnop)=kneg(i)
      	tn10(nnop)=lneg(i)
          else if (les1.and.les2.and.les3.and.les4) then
      	    npert=npert+1
      	    tp1(npert)=ip(i)
      	    tp2(npert)=jp(i)
      	    tp3(npert)=kp(i)
      	    tp4(npert)=lp(i)
      	    if (dscale.eq.3) then
      		tp6(npert)=ticp(i)
      	    else
      		tp6(npert)=icp(i)
      	    endif
      	    tp7(npert)=ineg(i)
      	    tp8(npert)=jneg(i)
      	    tp9(npert)=kneg(i)
      	    tp10(npert)=lneg(i)
      	    if ((lreal(ip(i)).and.lreal(jp(i)).and.
     &		    lreal(kp(i)).and.lreal(lp(i)))) then
      		if (dscale.eq.4) then
      	        	tp5(npert)=icp(i)
      		else
      	        	tp5(npert)=ticp(i)
      		endif
      	    else
c		        tp5(npert)=ticp(i)+realtors
      		if (dscale.eq.1.or.dscale.eq.4) then
      	        	tp5(npert)=icp(i)
      		else
      	        	tp5(npert)=ticp(i)
      		endif
      	    endif
          elseif (les1.or.les2.or.les3.or.les4) then
                if (per1.and.per2.and.per3.and.per4) then
                    npert=npert+1
                    tp1(npert)=ip(i)
                    tp2(npert)=jp(i)
                    tp3(npert)=kp(i)
                    tp4(npert)=lp(i)
                    if (dscale.eq.3) then
                        tp6(npert)=ticp(i)
                    else
                        tp6(npert)=icp(i)
                    endif
                    tp7(npert)=ineg(i)
                    tp8(npert)=jneg(i)
                    tp9(npert)=kneg(i)
                    tp10(npert)=lneg(i)
                    if ((lreal(ip(i)).and.lreal(jp(i)).and.
     &              lreal(kp(i)).and.lreal(lp(i)))) then
                        if (dscale.eq.4) then
                                tp5(npert)=icp(i)
                        else
                                tp5(npert)=ticp(i)
                        endif
                    else
c                       tp5(npert)=ticp(i)+realtors
                        if (dscale.eq.1.or.dscale.eq.4) then
                                tp5(npert)=icp(i)
                        else
                                tp5(npert)=ticp(i)
                        endif
                    endif
                else
      	    nboun=nboun+1
      	    tb1(nboun)=ip(i)
      	    tb2(nboun)=jp(i)
      	    tb3(nboun)=kp(i)
      	    tb4(nboun)=lp(i)
      	    if (dscale.eq.3) then
      		tb6(nboun)=ticp(i)
      	    else
      		tb6(nboun)=icp(i)
      	    endif
      	    tb7(nboun)=ineg(i)
      	    tb8(nboun)=jneg(i)
      	    tb9(nboun)=kneg(i)
      	    tb10(nboun)=lneg(i)
      	    if ((lreal(ip(i)).and.lreal(jp(i)).and.
     &		    lreal(kp(i)).and.lreal(lp(i)))) then
      	        if (dscale.eq.4) then
      			tb5(nboun)=icp(i)
      		else
      			tb5(nboun)=ticp(i)
      		endif
      	    else
c		        tb5(nboun)=ticp(i)+realtors
      	        if (dscale.eq.1.or.dscale.eq.4) then
      			tb5(nboun)=icp(i)
      		else
      			tb5(nboun)=ticp(i)
      		endif
      	    endif
      	endif
            elseif (per1.and.per2.and.per3.and.per4) then
                    npert=npert+1
                    tp1(npert)=ip(i)
                    tp2(npert)=jp(i)
                    tp3(npert)=kp(i)
                    tp4(npert)=lp(i)
                    tp5(npert)=icp(i)
                    tp6(npert)=icp(i)
                    tp7(npert)=ineg(i)
                    tp8(npert)=jneg(i)
                    tp9(npert)=kneg(i)
                    tp10(npert)=lneg(i)
            else
                    nboun=nboun+1
                    tb1(nboun)=ip(i)
                    tb2(nboun)=jp(i)
                    tb3(nboun)=kp(i)
                    tb4(nboun)=lp(i)
                    tb5(nboun)=icp(i)
                    tb6(nboun)=icp(i)
                    tb7(nboun)=ineg(i)
                    tb8(nboun)=jneg(i)
                    tb9(nboun)=kneg(i)
                    tb10(nboun)=lneg(i)
          endif
270	continue
      do 280 i=1,ndper
          les1=(totcop(origpt(ipper(i))).gt.1)
          les2=(totcop(origpt(jpper(i))).gt.1)
          les3=(totcop(origpt(kpper(i))).gt.1)
          les4=(totcop(origpt(lpper(i))).gt.1)
            per1=(ipert(origpt(ipper(i))).eq.1.or.les1)
            per2=(ipert(origpt(jpper(i))).eq.1.or.les2)
            per3=(ipert(origpt(kpper(i))).eq.1.or.les3)
            per4=(ipert(origpt(lpper(i))).eq.1.or.les4)
          if (.not.per1.and..not.per2.and..not.per3.and.
     &	    .not.per4) then
      	if (isymbl(ipper(i))(1:1).eq.'H'.or.
     &		isymbl(jpper(i))(1:1).eq.'H'.or.
     &		isymbl(kpper(i))(1:1).eq.'H'.or.
     &		isymbl(lpper(i))(1:1).eq.'H'.or.
     &		isymbl(ipper(i))(1:2).eq.'DH'.or.
     &		isymbl(jpper(i))(1:2).eq.'DH'.or.
     &		isymbl(kpper(i))(1:2).eq.'DH'.or.
     &		isymbl(lpper(i))(1:2).eq.'DH') then
      		nnoph=nnoph+1
      		th1(nnoph)=ipper(i)
      		th2(nnoph)=jpper(i)
      		th3(nnoph)=kpper(i)
      		th4(nnoph)=lpper(i)
      		th5(nnoph)=icpper(i+ndper)
      	        th7(nnoph)=inegp(i)
      	        th8(nnoph)=jnegp(i)
      	        th9(nnoph)=knegp(i)
      	        th10(nnoph)=lnegp(i)
              else
      		nnop=nnop+1
      		tn1(nnop)=ipper(i)
      		tn2(nnop)=jpper(i)
      		tn3(nnop)=kpper(i)
      		tn4(nnop)=lpper(i)
      		tn5(nnop)=icpper(i+ndper)
      	        tn7(nnop)=inegp(i)
      	        tn8(nnop)=jnegp(i)
      	        tn9(nnop)=knegp(i)
      	        tn10(nnop)=lnegp(i)
      	endif
          elseif (les1.and.les2.and.les3.and.les4) then
      	    npert=npert+1
      	    tp1(npert)=ipper(i)
      	    tp2(npert)=jpper(i)
      	    tp3(npert)=kpper(i)
      	    tp4(npert)=lpper(i)
      	    if (dscale.eq.3) then
      		tp6(npert)=ticpper(i+ndper)
      	    else
      		tp6(npert)=icpper(i+ndper)
      	    endif
      	    tp7(npert)=inegp(i)
      	    tp8(npert)=jnegp(i)
      	    tp9(npert)=knegp(i)
      	    tp10(npert)=lnegp(i)
      	    if ((lreal(ipper(i)).and.lreal(jpper(i)).and.
     &		    lreal(kpper(i)).and.lreal(lpper(i)))) then
      	        if (dscale.eq.4) then
      			tp5(npert)=icpper(i+ndper)
      		else
      			tp5(npert)=ticpper(i+ndper)
      		endif
      	    else
c		        tp5(npert)=ticpper(i+ndper)+realtors
      	        if (dscale.eq.1.or.dscale.eq.4) then
      			tp5(npert)=icpper(i+ndper)
      		else
      			tp5(npert)=ticpper(i+ndper)
      		endif
      	    endif
          elseif (les1.or.les2.or.les3.or.les4) then
                if (per1.and.per2.and.per3.and.per4) then
                    npert=npert+1
                    tp1(npert)=ipper(i)
                    tp2(npert)=jpper(i)
                    tp3(npert)=kpper(i)
                    tp4(npert)=lpper(i)
                    if (dscale.eq.3) then
                        tp6(npert)=ticpper(i+ndper)
                    else
                        tp6(npert)=icpper(i+ndper)
                    endif
                    tp7(npert)=inegp(i)
                    tp8(npert)=jnegp(i)
                    tp9(npert)=knegp(i)
                    tp10(npert)=lnegp(i)
                    if ((lreal(ipper(i)).and.lreal(jpper(i)).and.
     &              lreal(kpper(i)).and.lreal(lpper(i)))) then
                        if (dscale.eq.4) then
                                tp5(npert)=icpper(i+ndper)
                        else
                                tp5(npert)=ticpper(i+ndper)
                        endif
                    else
c                       tp5(npert)=ticpper(i)+realtors
                        if (dscale.eq.1.or.dscale.eq.4) then
                                tp5(npert)=icpper(i+ndper)
                        else
                                tp5(npert)=ticpper(i+ndper)
                        endif
                    endif
                else
      	    nboun=nboun+1
      	    tb1(nboun)=ipper(i)
      	    tb2(nboun)=jpper(i)
      	    tb3(nboun)=kpper(i)
      	    tb4(nboun)=lpper(i)
      	    if (dscale.eq.3) then
      		tb6(nboun)=ticpper(i+ndper)
      	    else
      		tb6(nboun)=icpper(i+ndper)
      	    endif
      	    tb7(nboun)=inegp(i)
      	    tb8(nboun)=jnegp(i)
      	    tb9(nboun)=knegp(i)
      	    tb10(nboun)=lnegp(i)
      	    if ((lreal(ipper(i)).and.lreal(jpper(i)).and.
     &		    lreal(kpper(i)).and.lreal(lpper(i)))) then
      		if (dscale.eq.4) then
      	        	tb5(nboun)=icpper(i+ndper)
      		else
      	        	tb5(nboun)=ticpper(i+ndper)
      		endif
      	    else
c		        tb5(nboun)=ticpper(i+ndper)+realtors
      		if (dscale.eq.1.or.dscale.eq.4) then
      	        	tb5(nboun)=icpper(i+ndper)
      		else
      	        	tb5(nboun)=ticpper(i+ndper)
      		endif
      	    endif
      	endif
            elseif (per1.and.per2.and.per3.and.per4) then
                    npert=npert+1
                    tp1(npert)=ipper(i)
                    tp2(npert)=jpper(i)
                    tp3(npert)=kpper(i)
                    tp4(npert)=lpper(i)
                    tp5(npert)=icpper(i+ndper)
                    tp6(npert)=icpper(i+ndper)
                    tp7(npert)=inegp(i)
                    tp8(npert)=jnegp(i)
                    tp9(npert)=knegp(i)
                    tp10(npert)=lnegp(i)
            else
                    nboun=nboun+1
                    tb1(nboun)=ipper(i)
                    tb2(nboun)=jpper(i)
                    tb3(nboun)=kpper(i)
                    tb4(nboun)=lpper(i)
                    tb5(nboun)=icpper(i+ndper)
                    tb6(nboun)=icpper(i+ndper)
                    tb7(nboun)=inegp(i)
                    tb8(nboun)=jnegp(i)
                    tb9(nboun)=knegp(i)
                    tb10(nboun)=lnegp(i)
          endif
280	continue
      nnopc=nnop
      do 290 i=mphia+1,nphia
c
c constraints
c
      	write (6,*) 'Constraint torsions... untested!'
          les1=(totcop(origpt(ip(i))).gt.1)
          les2=(totcop(origpt(jp(i))).gt.1)
          les3=(totcop(origpt(kp(i))).gt.1)
          les4=(totcop(origpt(lp(i))).gt.1)
            per1=(ipert(origpt(ip(i))).eq.1.or.les1)
            per2=(ipert(origpt(jp(i))).eq.1.or.les2)
            per3=(ipert(origpt(kp(i))).eq.1.or.les3)
            per4=(ipert(origpt(lp(i))).eq.1.or.les4)
          if (.not.per1.and..not.per2.and..not.per3.and.
     &	    .not.per4) then
      	nnop=nnop+1
      	tn1(nnop)=ip(i)
      	tn2(nnop)=jp(i)
      	tn3(nnop)=kp(i)
      	tn4(nnop)=lp(i)
      	tn5(nnop)=icp(i)
      	tn7(nnop)=ineg(i)
      	tn8(nnop)=jneg(i)
      	tn9(nnop)=kneg(i)
      	tn10(nnop)=lneg(i)
          else if (les1.and.les2.and.les3.and.les4) then
      	    npert=npert+1
      	    tp1(npert)=ip(i)
      	    tp2(npert)=jp(i)
      	    tp3(npert)=kp(i)
      	    tp4(npert)=lp(i)
      	    if (dscale.eq.3) then
      		tp6(npert)=ticp(i)
      	    else
      		tp6(npert)=icp(i)
      	    endif
      	    tp7(npert)=ineg(i)
      	    tp8(npert)=jneg(i)
      	    tp9(npert)=kneg(i)
      	    tp10(npert)=lneg(i)
      	    if ((lreal(ip(i)).and.lreal(jp(i)).and.
     &		    lreal(kp(i)).and.lreal(lp(i)))) then
      	        if (dscale.eq.4) then
      			tp5(npert)=icp(i)
      		else
      			tp5(npert)=ticp(i)
      		endif
      	    else
c		        tp5(npert)=ticp(i)+realtors
      	        if (dscale.eq.1.or.dscale.eq.4) then
      			tp5(npert)=icp(i)
      		else
      			tp5(npert)=ticp(i)
      		endif
      	    endif
          elseif (les1.or.les2.or.les3.or.les4) then
                if (per1.and.per2.and.per3.and.per4) then
                    npert=npert+1
                    tp1(npert)=ip(i)
                    tp2(npert)=jp(i)
                    tp3(npert)=kp(i)
                    tp4(npert)=lp(i)
                    if (dscale.eq.3) then
                        tp6(npert)=ticp(i)
                    else
                        tp6(npert)=icp(i)
                    endif
                    tp7(npert)=ineg(i)
                    tp8(npert)=jneg(i)
                    tp9(npert)=kneg(i)
                    tp10(npert)=lneg(i)
                    if ((lreal(ip(i)).and.lreal(jp(i)).and.
     &              lreal(kp(i)).and.lreal(lp(i)))) then
                        if (dscale.eq.4) then
                                tp5(npert)=icp(i)
                        else
                                tp5(npert)=ticp(i)
                        endif
                    else
c                       tp5(npert)=ticp(i)+realtors
                        if (dscale.eq.1.or.dscale.eq.4) then
                                tp5(npert)=icp(i)
                        else
                                tp5(npert)=ticp(i)
                        endif
                    endif
                else
      	    nboun=nboun+1
      	    tb1(nboun)=ip(i)
      	    tb2(nboun)=jp(i)
      	    tb3(nboun)=kp(i)
      	    tb4(nboun)=lp(i)
      	    if (dscale.eq.3) then
      		tb6(nboun)=ticp(i)
      	    else
      		tb6(nboun)=icp(i)
      	    endif
      	    tb7(nboun)=ineg(i)
      	    tb8(nboun)=jneg(i)
      	    tb9(nboun)=kneg(i)
      	    tb10(nboun)=lneg(i)
      	    if ((lreal(ip(i)).and.lreal(jp(i)).and.
     &		    lreal(kp(i)).and.lreal(lp(i)))) then
      	        if (dscale.eq.4) then
      			tb5(nboun)=icp(i)
      		else
      			tb5(nboun)=ticp(i)
      		endif
      	    else
c		        tb5(nboun)=ticp(i)+realtors
      	        if (dscale.eq.1.or.dscale.eq.4) then
      			tb5(nboun)=icp(i)
      		else
      			tb5(nboun)=ticp(i)
      		endif
      	    endif
      	endif
            elseif (per1.and.per2.and.per3.and.per4) then
                    npert=npert+1
                    tp1(npert)=ip(i)
                    tp2(npert)=jp(i)
                    tp3(npert)=kp(i)
                    tp4(npert)=lp(i)
                    tp5(npert)=icp(i)
                    tp6(npert)=icp(i)
                    tp7(npert)=ineg(i)
                    tp8(npert)=jneg(i)
                    tp9(npert)=kneg(i)
                    tp10(npert)=lneg(i)
            else
                    nboun=nboun+1
                    tb1(nboun)=ip(i)
                    tb2(nboun)=jp(i)
                    tb3(nboun)=kp(i)
                    tb4(nboun)=lp(i)
                    tb5(nboun)=icp(i)
                    tb6(nboun)=icp(i)
                    tb7(nboun)=ineg(i)
                    tb8(nboun)=jneg(i)
                    tb9(nboun)=kneg(i)
                    tb10(nboun)=lneg(i)
          endif
290	continue
c
c now move all this to the t3 arrays to free up the temp ones for the dihedrals
c
      do 300 i=1,nnop
      	t3ip(i)=tn1(i)
      	t3jp(i)=tn2(i)
      	t3kp(i)=tn3(i)
      	t3lp(i)=tn4(i)
      	t3icp(i)=tn5(i)
      	t3ineg(i)=tn7(i)
      	t3jneg(i)=tn8(i)
      	t3kneg(i)=tn9(i)
      	t3lneg(i)=tn10(i)
300	continue
      do 305 i=1,nnoph
      	t3iph(i)=th1(i)
      	t3jph(i)=th2(i)
      	t3kph(i)=th3(i)
      	t3lph(i)=th4(i)
      	t3icph(i)=th5(i)
      	t3inegh(i)=th7(i)
      	t3jnegh(i)=th8(i)
      	t3knegh(i)=th9(i)
      	t3lnegh(i)=th10(i)
305	continue
      do 310 i=1,npert
      	t3ipper(i)=tp1(i)
      	t3jpper(i)=tp2(i)
      	t3kpper(i)=tp3(i)
      	t3lpper(i)=tp4(i)
      	t3icpper(i)=tp6(i)
      	t3icpper(i+npert+nboun)=tp5(i)
      	t3inegp(i)=tp7(i)
      	t3jnegp(i)=tp8(i)
      	t3knegp(i)=tp9(i)
      	t3lnegp(i)=tp10(i)
310	continue
      do 315 i=1,nboun
      	j=i+npert
      	t3ipper(j)=tb1(i)
      	t3jpper(j)=tb2(i)
      	t3kpper(j)=tb3(i)
      	t3lpper(j)=tb4(i)
      	t3icpper(j)=tb6(i)
      	t3icpper(j+npert+nboun)=tb5(i)
      	t3inegp(j)=tb7(i)
      	t3jnegp(j)=tb8(i)
      	t3knegp(j)=tb9(i)
      	t3lnegp(j)=tb10(i)
315	continue
      t3nphih=nnoph
      t3nphia=nnop
      t3mphia=nnopc
      t3ndper=npert+nboun
      t3mdper=npert
c
c now write it....
c
      call wlesprm(uw2prm)
      
      return

      end
