#
#************************************************************************
#                             AMBER                                    **
#                                                                      **
#                Copyright (c) 1986, 1991, 1995                        **
#             Regents of the University of California                  **
#                       All Rights Reserved.                           ** 
#                                                                      **
#  This software provided pursuant to a license agreement containing   **
#  restrictions on its disclosure, duplication, and use. This software **
#  contains confidential and proprietary information, and may not be   **
#  extracted or distributed, in whole or in part, for any purpose      **
#  whatsoever, without the express written permission of the authors.  **
#  This notice, and the associated author list, must be attached to    **
#  all copies, or extracts, of this software. Any additional           **
#  restrictions set forth in the license agreement also apply to this  **
#  software.                                                           **
#************************************************************************
#
SHELL=/bin/sh
VPATH=
include ../config.h
.SUFFIXES:
.SUFFIXES: .F90 .c .o

# Nmode Fortran source files are free format
LOCALFLAGS = $(FREEFORMAT_FLAG)

# NMODE is always double precision

SRC=	nmode.F90 alloc.F90 angl.F90 assert.F90 \
	bond.F90 bresna.F90 cdiv.F90 chck.F90 chdir.F90 chlamba.F90 \
	debug.F90 decomp.F90 difang.F90 difbon.F90 dihdup.F90 dihed.F90 \
	dirch.F90 dsort.F90 \
	eigout.F90 ephi.F90 eprnt.F90 errstp.F90 \
	findl.F90 forces.F90 \
	getcor.F90 getf.F90 getm.F90 getref.F90 grdmax.F90 \
	level.F90 lmdout.F90 lmode.F90 \
	makel.F90 movecm.F90 mweight.F90 mweit.F90 \
	new_stuff.F90 \
	nmdfil.F90 nonbon.F90 nonbon1.F90 norm.F90 norml.F90 nrch.F90 \
	oseen.F90 pairr.F90 pairw.F90 percon.F90 pload.F90 putvar.F90 \
	rdinp.F90 rdparm.F90 remrot.F90 resdis.F90 \
	savec.F90 savei.F90  savet.F90 setang.F90 setatm.F90 setbon.F90 \
	setcor.F90 setdih.F90 setgam.F90 setgr.F90 setvar.F90 spsort.F90 \
	thermo.F90 timit.F90 traco.F90  tsearch.F90 tstate.F90 \
	xconst.F90 zerout.F90 zxcgr.F90 mexit.F90

OBJ =	nmode.o alloc.o angl.o assert.o \
	bond.o bresna.o cdiv.o chck.o chdir.o chlamba.o \
	debug.o decomp.o difang.o difbon.o dihdup.o dihed.o \
	dirch.o dsort.o \
	eigout.o ephi.o eprnt.o errstp.o \
	findl.o forces.o \
	getcor.o getf.o getm.o getref.o grdmax.o \
	level.o lmdout.o lmode.o \
	makel.o movecm.o mweight.o mweit.o \
	new_stuff.o \
	nmdfil.o nonbon.o nonbon1.o norm.o norml.o nrch.o \
	oseen.o pairr.o pairw.o percon.o pload.o putvar.o \
	rdinp.o rdparm.o remrot.o resdis.o \
	savec.o savei.o  savet.o setang.o setatm.o setbon.o \
	setcor.o setdih.o setgam.o setgr.o setvar.o spsort.o \
	thermo.o timit.o traco.o  tsearch.o tstate.o \
	xconst.o zerout.o zxcgr.o mexit.o

SRC1=   nmanal.F90 mdread.F90 rdvect.F90 proj.F90 stigrp.F90 setgrp.F90\
        group.F90 grpprj.F90 rdparm2.F90 rmsf.F90 excise.F90 dotref.F90\
        prtma.F90 reori.F90 movecm.F90 bmat.F90 belly.F90 quasip.F90\
        fluct.F90 nmdfil.F90 corf.F90 mexit.F90 lmdin.F90 rgroup2.F90 \
        corfl.F90 corffft.F90 corffft2.F90 ../sander/pubfft.F90 fileio.F90

OBJ1=   mdread.o nmanal.o rdvect.o proj.o stigrp.o setgrp.o\
        group.o grpprj.o rdparm2.o rmsf.o excise.o dotref.o\
        prtma.o reori.o movecm.o bmat.o belly.o quasip.o\
        fluct.o nmdfil.o corf.o mexit.o lmdin.o \
        corfl.o corffft.o corffft2.o pubfft.o fileio.o rgroup2.o

OBJ2=   lmanal.o distr.o lmfil.o mexit.o count.o \
        fluct.o lmdin.o dsort.o getref.o fileio.o \
        corfl.o corffft.o corffft2.o pubfft.o 

LIBSRC= ../lib/rfree.F90   ../lib/matinv.F90  ../lib/nxtsec.F90 \
	../lib/amopen.F90 
LIBOBJ= ../lib/rfree.o   ../lib/matinv.o ../lib/nxtsec.o \
	../lib/amopen.o ../lib/rgroup.o
LIBOBJNOPATH= rfree.o   matinv.o nxtsec.o \
	amopen.o rgroup.o

$(BINDIR)/nmode$(SFX):	$(OBJ) libobj syslib netlib
	$(FC) $(FPPFLAGS) $(FFLAGS) $(AMBERFFLAGS) -o $(BINDIR)/nmode$(SFX) \
		$(OBJ) $(LIBOBJ) \
		../lib/sys.a -L$(LIBDIR) $(FLIBSF) $(LDFLAGS) $(AMBERLDFLAGS)

nmanal$(SFX):	$(OBJ1) libobj netlib binpos.o
	$(FC) $(FPPFLAGS) $(FFLAGS) $(AMBERFFLAGS) -o nmanal$(SFX) \
		$(OBJ1) $(LIBOBJ) binpos.o \
		-L$(LIBDIR) $(FLIBSF) $(LDFLAGS) $(AMBERLDFLAGS)

lmanal$(SFX):	$(OBJ2) libobj binpos.o
	$(FC) $(FPPFLAGS) $(FFLAGS) $(AMBERFFLAGS) -o lmanal$(SFX) \
		$(OBJ2) $(LIBOBJ) binpos.o \
        -L$(LIBDIR) $(FLIBSF) $(LDFLAGS) $(AMBERLDFLAGS)

quasih$(SFX):     quasih.o fileio.o thermo.o mexit.o \
				syslib netlib libobj binpos.o
	$(FC) $(FPPFLAGS) $(FFLAGS) $(AMBERFFLAGS) -o quasih$(SFX) \
		quasih.o fileio.o thermo.o mexit.o \
		$(LIBOBJ) binpos.o ../lib/sys.a -L$(LIBDIR) $(FLIBSF) \
		$(LDFLAGS) $(AMBERLDFLAGS)

calink$(SFX):     calink.o fileio.o netlib libobj syslib
	$(FC) $(FPPFLAGS) $(FFLAGS) $(AMBERFFLAGS) -o calink$(SFX) \
		calink.o fileio.o $(LIBOBJ) -L$(LIBDIR) $(FLIBSF) \
		../lib/sys.a $(LDFLAGS) $(AMBERLDFLAGS)

vecs_to_binpos$(SFX):     vecs_to_binpos.o fileio.o libobj syslib binpos.o
	$(FC) $(FPPFLAGS) $(FFLAGS) $(AMBERFFLAGS) -o vecs_to_binpos$(SFX) \
		vecs_to_binpos.o fileio.o \
		$(LIBOBJ) binpos.o ../lib/sys.a \
		$(LDFLAGS) $(AMBERLDFLAGS)

libobj::
	cd ../lib && $(MAKE) $(LIBOBJNOPATH)

netlib:  
	cd ../blas && $(MAKE) $(BLAS)
	cd ../lapack && $(MAKE) $(LAPACK)
	cd ../arpack && $(MAKE) install

syslib::
	cd ../lib && $(MAKE) sys.a

binpos.o:  ../lib/binpos.c
	$(CC) -c $(CNOOPTFLAGS) $(CFLAGS) $(AMBERCFLAGS) \
		-o binpos.o ../lib/binpos.c

#-----------

calink.o:	calink.F90
	$(FC) $(FPPFLAGS) -c $(FOPTFLAGS) $(FFLAGS) $(AMBERFFLAGS) -o $@ $<

rdparm.o:	rdparm.F90
	$(FC) $(FPPFLAGS) -c $(FOPTFLAGS) $(FFLAGS) $(AMBERFFLAGS) -o $@ $<

rdparm2.o:	rdparm2.F90
	$(FC) $(FPPFLAGS) -c $(FOPTFLAGS) $(FFLAGS) $(AMBERFFLAGS) -o $@ $<

nonbon.o:	nonbon.F90
	$(FC) $(FPPFLAGS) -c $(FOPTFLAGS) $(FFLAGS) $(AMBERFFLAGS) -o $@ $<

nonbon1.o:   	nonbon1.F90
	$(FC) $(FPPFLAGS) -c $(FOPTFLAGS) $(FFLAGS) $(AMBERFFLAGS) -o $@ $<

bond.o:   	bond.F90
	$(FC) $(FPPFLAGS) -c $(FOPTFLAGS) $(FFLAGS) $(AMBERFFLAGS) -o $@ $<

angl.o:		angl.F90
	$(FC) $(FPPFLAGS) -c $(FOPTFLAGS) $(FFLAGS) $(AMBERFFLAGS) -o $@ $<

bresna.o:	bresna.F90
	$(FC) $(FPPFLAGS) -c $(FOPTFLAGS) $(FFLAGS) $(AMBERFFLAGS) -o $@ $<

resdis.o:	resdis.F90
	$(FC) $(FPPFLAGS) -c $(FOPTFLAGS) $(FFLAGS) $(AMBERFFLAGS) -o $@ $<

decomp.o:	decomp.F90
	$(FC) $(FPPFLAGS) -c $(FOPTFLAGS) $(FFLAGS) $(AMBERFFLAGS) -o $@ $<

pubfft.o:	../sander/pubfft.F90
	$(FC) $(FPPFLAGS) -c $(FREEFORMAT_FLAG) $(FOPTFLAGS) $(FFLAGS) \
		$(AMBERFFLAGS) -o $@ ../sander/pubfft.F90
	
#----- header dependencies:
alloc.o: pointer.h
alloc.o: inpdat.h
belly.o: sizes2.h
belly.o: bad.h
belly.o: carts.h
belly.o: infoa.h
belly.o: optf.h
belly.o: opta.h
bmat.o: sizes2.h
bmat.o: bad.h
chck.o: pointer.h
chck.o: epot.h
chck.o: inpdat.h
chdir.o: inpdat.h
chlamba.o: inpdat.h
corfl.o: files.h
corffft.o: sizes2.h
corffft.o: anal.h
corffft.o: files.h
corffft.o: infoa.h
distr.o: files.h
dotref.o: sizes2.h
eprnt.o: epot.h
excise.o: sizes2.h
excise.o: infoa.h
fileio.o: fileio.h
fluct.o: sizes2.h
fluct.o: anal.h
forces.o: pointer.h
forces.o: inpdat.h
forces.o: epot.h
forces.o: pol.h
forces.o: sizes2.h
group.o: sizes2.h
group.o: infoa.h
grpprj.o: anal.h
lmdin.o: files.h
lmfil.o: files.h
lmode.o: files.h
mdread.o: sizes2.h
mdread.o: anal.h
mdread.o: infoa.h
mdread.o: optf.h
mdread.o: opta.h
mdread.o: files.h
mweight.o: sizes2.h
nmanal.o: sizes2.h
nmanal.o: anal.h
nmanal.o: bad.h
nmanal.o: files.h
nmanal.o: infoa.h
nmanal.o: carts.h
nmanal.o: optf.h
nmanal.o: opta.h
nmdfil.o: files.h
nmode.o: files.h
nmode.o: pointer.h
nmode.o: inpdat.h
nmode.o: epot.h
nmode.o: decomp.o
nrch.o: pointer.h
nrch.o: inpdat.h
pairw.o: files.h
proj.o: sizes2.h
proj.o: anal.h
proj.o: infoa.h
quasip.o: sizes2.h
quasip.o: anal.h
quasip.o: files.h
rdinp.o: inpdat.h
rdparm.o: pointer.h
rdparm.o: inpdat.h
rdparm.o: pol.h
rdparm2.o: sizes2.h
rdparm2.o: bad.h
rdparm2.o: files.h
rdparm2.o: infoa.h
rdparm2.o: optf.h
rdparm2.o: carts.h
rdvect.o: sizes2.h
reori.o: sizes2.h
rmsf.o: sizes2.h
rmsf.o: anal.h
rmsf.o: infoa.h
savec.o: files.h
savec.o: inpdat.h
savet.o: files.h
savet.o: inpdat.h
setgam.o: files.h
setvar.o: pointer.h
setvar.o: inpdat.h
tsearch.o: pointer.h
tsearch.o: inpdat.h
tstate.o: pointer.h
tstate.o: inpdat.h
calink.o: calink.h

#-----

install:	$(BINDIR)/nmode$(SFX)

clean:
		/bin/rm -f *.o  nmanal$(SFX) *.mod

uninstall:
		-rm -f $(BINDIR)/nmode$(SFX)

#  default rules for Fortran and C compilation:

.F90.o:   $<
	$(FC) $(FPPFLAGS) -c $(FNOOPTFLAGS) $(FFLAGS) $(AMBERFFLAGS) -o $@ $<

.c.o:
	$(CC) -c $(CPPFLAGS) $(CNOOPTFLAGS) $(CFLAGS) $(AMBERCFLAGS) -o $@ $<
