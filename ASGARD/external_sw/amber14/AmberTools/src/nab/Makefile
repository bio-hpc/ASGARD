include ../config.h

.SUFFIXES:
.SUFFIXES:  .nab .c .o .f

%.o: %.nab $(BINDIR)/$(MPI)nab$(SFX) $(BINDIR)/$(MPI)nab2c$(SFX)
	$(BINDIR)/$(MPI)nab$(SFX) -c $(NABFLAGS) $<
	@rm -f $*.c

NAB2COBJS=	\
	nab2c.o	\
	cgen.o	\
	checkexpr.o	\
	dumpnode.o\
	errormsg.o	\
	fixexpr.o	\
	genavscode.o	\
	node.o	\
	symbol.o	\
	y.tab.o

LIBNABOBJS=	\
	allatom_to_dna3.o \
	axis2frame.o	\
	bonded_atoms.o \
	boundsutil.o	\
	chirvol.o	\
	database.o  \
	deform_ener.o \
	dna3.o		\
	dna3_to_allatom.o \
	dump.o		\
	embed.o		\
	fd_helix.o  \
	frameutil.o	\
	getdist.o	\
	gsub.o	\
	hashutil.o	\
	jacobi.o	\
	match.o		\
	matop.o		\
	molio.o		\
	molsurf.o   \
	molutil.o	\
	nab2avs.o	\
	nab2geom.o	\
	pointutil.o	\
	prm.o \
	regexp.o	\
	reslib.o	\
	select_atoms.o	\
	traceback.o \
	mme_init.o \
	stringutil.o	\
	temputil.o	\
	tm_malloc.o	\
	bdna.o	\
	dg_helix.o \
	getchivol.o \
	getres.o	\
	getseq_from_pdb.o \
	getxyz_from_pdb.o \
	getpdb_prm.o \
	getpdb_rlb.o \
	linkprot.o	\
	link_na.o	\
	na_anal.o 	\
	setboundsfromdb.o \
	setchiplane.o	\
	wc_basepair.o	\
	wc_complement.o	\
	wc_helix.o \
	yammp.o

EXAMPLES=	\
	../examples/program_*.nab

#===========================================================================

install: $(BINDIR)/$(MPI)nab$(SFX) $(BINDIR)/$(MPI)nab2c$(SFX) libnab.a
	mv -f  libnab.a $(LIBDIR)
	rm -f $(INCDIR)/nabcode.h $(INCDIR)/defreal.h $(INCDIR)/nab.h
	cp nabcode.h $(INCDIR)
	cp nabtypes.h $(INCDIR)
	cp defreal.h $(INCDIR)
	cp nab.h $(INCDIR)
	cp xmin_opt.h lmod_opt.h nab_netcdf.h $(INCDIR)
#	cp lmodprmtop $(BINDIR)

#  Force rebuilding of nab, since there is no easy to to see if FLIBS has
#  changed.
$(BINDIR)/$(MPI)nab$(SFX)::
	$(CC) -DCC='"$(CC)"' -DCPP='"$(CPP)"' -DFLIBS='"$(FLIBS)"' \
		$(CNOOPTFLAGS) $(CFLAGS) $(AMBERCFLAGS) \
		$(LDFLAGS) $(AMBERLDFLAGS) -o $(BINDIR)/$(MPI)nab$(SFX) nab.c

$(BINDIR)/$(MPI)nab2c$(SFX):	$(NAB2COBJS)
	$(CC) $(CFLAGS) $(AMBERCFLAGS) $(LDFLAGS) $(AMBERLDFLAGS) \
		-o $(BINDIR)/$(MPI)nab2c$(SFX) $(NAB2COBJS)

#===========================================================================

libnab.a: $(LIBNABOBJS) $(BINDIR)/$(MPI)nab$(SFX) $(BINDIR)/$(MPI)nab2c$(SFX)
	-rm -f libnab.a
	$(AR) libnab.a $(LIBNABOBJS)
	-ranlib libnab.a

ms$(SFX):  ms.o
	$(CC) $(CFLAGS) $(AMBERCFLAGS) $(LDFLAGS) $(AMBERLDFLAGS) \
		-o ms$(SFX)  ms.o ../lib/f2c.a -lm

y.tab.h: y.tab.c

y.tab.c:  nabgrm.y lex.yy.c
	$(YACC) -d nabgrm.y

y.tab.o: y.tab.c
	$(CC) -c -D$(LEX) $(CFLAGS) $(AMBERCFLAGS) y.tab.c

lex.yy.c: nablex.l
	$(LEX) nablex.l

lex.dg_options.c: dg_options.l
	$(LEX) -t dg_options.l | sed -e 's/{stdout}/{DONTchangeSTDOUThere}/' -e 's/yy/dgo/g' -e 's/stdout/nabout/' -e 's/DONTchangeSTDOUThere/stdout/' > lex.dg_options.c

checkexpr.c: y.tab.o
	(cd ../semantics; $(MAKE))
	mv ../semantics/checkexpr.c .


#    Compiles with optimized code:

molio.o:  molio.c
	$(CC) -c $(COPTFLAGS) $(CFLAGS) $(AMBERCFLAGS) molio.c

mme_init.o:  mme_init.c
	$(CC) -c -D$(LEX) $(COPTFLAGS) $(CFLAGS) $(AMBERCFLAGS) mme_init.c

boundsutil.o:  boundsutil.c
	$(CC) -c -DAMBERHOME='"$(AMBERHOME)"' $(COPTFLAGS) $(CFLAGS) $(AMBERCFLAGS) boundsutil.c

embed.o:  embed.c
	$(CC) -c -D$(LEX) $(COPTFLAGS) $(CFLAGS) $(AMBERCFLAGS) embed.c

bonded_atoms.o:  bonded_atoms.c
	$(CC) -c $(COPTFLAGS) $(CFLAGS) $(AMBERCFLAGS) bonded_atoms.c

ms.o:  ms.c
	$(CC) -c $(COPTFLAGS) $(CFLAGS) $(AMBERCFLAGS) ms.c

sasad.o:  sasad.c
	$(CC) -c $(COPTFLAGS) $(CFLAGS) $(AMBERCFLAGS) sasad.c

#===========================================================================

clean:
	-rm -f *.o nab$(SFX) nab2c$(SFX) libnab.a checkexpr.c
	-rm -f y.tab.c y.tab.h lex.dg_options.c lex.yy.c
	(cd ../semantics && $(MAKE) clean )

uninstall: clean
	-rm -f \
	$(BINDIR)/nab$(SFX) \
        $(BINDIR)/nab2c$(SFX) \
	$(BINDIR)/mpinab$(SFX) \
        $(BINDIR)/mpinab2c$(SFX) \
	$(LIBDIR)/libnab.a \
	$(INCDIR)/nabcode.h \
	$(INCDIR)/nabtypes.h \
	$(INCDIR)/defreal.h \
        $(INCDIR)/nab.h \
        $(INCDIR)/nab_netcdf.h \
	$(INCDIR)/xmin_opt.h \
	$(INCDIR)/lmod_opt.h \
        $(BINDIR)/lmodprmtop

examples: serial 
	@for i in $(EXAMPLES); do $(BINDIR)/nab $(NABFLAGS) -o `ls $$i | sed 's/.nab/$(SFX)/'` $$i; done
	(cd ../examples; $(MAKE) examples )

#===========================================================================

axis2frame.o:	\
	nab.h
binpos.o:  \
	memutil.h
bonded_atoms.o: \
	nab.h
boundsutil.o:	\
	nab.h  \
	memutil.h \
    traceback.h \
	errormsg.h
cgen.o:	\
	nab.h	\
	errormsg.h	\
	y.tab.h
chirvol.o:	\
	nab.h
database.o:  \
	database.h
dump.o:		\
	nab.h
dumpnode.o:	\
	nab.h	\
	y.tab.h
embed.o:	\
	nab.h	\
	errormsg.h  \
	memutil.h  \
    traceback.h \
    chirvol.h \
	lex.dg_options.c
errormsg.o:	\
	errormsg.h
fixexpr.o:	\
	nab.h	\
	errormsg.h	\
	y.tab.h
frameutil.o:	\
	nab.h	\
	memutil.h
genavscode.o:	\
	nab.h	\
	y.tab.h
gsub.o:    \
	nab.h
hashutil.o:	\
	nab.h
jacobi.o:	\
	memutil.h
matop.o: \
	nabcode.h
match.o:	\
	nab.h
memutil.o:	\
	nab.h
molio.o:	\
	nab.h	\
	errormsg.h
molsurf.o:  \
	nab.h   \
	molsurf.h
molutil.o:	\
	nab.h	\
	errormsg.h	\
	memutil.h
nab2avs.o:	\
	nab.h
nab2c.o:	\
	nab.h
nab2geom.o:	\
	nab.h
node.o:	\
	nab.h	\
	errormsg.h	\
	y.tab.h
pointutil.o: \
	nab.h
prm.o: \
	nab.h
reslib.o: \
	nab.h \
	errormsg.h \
	memutil.h  \
    chirvol.h  \
	database.h
select_atoms.o:	\
	nab.h
mme_init.o: \
	nab.h 
stringutil.o: \
	nab.h
symbol.o:	\
	nab.h	\
	errormsg.h
temputil.o: \
	nab.h
traceback.o: \
	errormsg.h
y.tab.o:	\
	nab.h	\
	errormsg.h \
	lex.yy.c

