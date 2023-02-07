# -*-Makefile-*-

# This should point to the location of config.h for your
# amber/nab installation.
include $(AMBERHOME)/src/config.h

PROG=		peptide
OBJS=		lprot.o peptide.o get_the_opts.o
HDRS=		$(PROG).h get_the_opts.h
CFLAGS+=	-g -pipe -I/usr/local/include -I$(AMBERHOME)/include
LDFLAGS+=	-L$(AMBERHOME)/lib -L/usr/local/lib
NABC=		nab
CLEANFILES=	lprot.c peptide.c A12.pdb
PREFIX?=	$(HOME)
INSTALL_USER=	mldodson
INSTALL_GRP=	mldodson
DESTDIR=	$(PREFIX)/
BINDIR=		bin
MANDIR=		man
MAN1=           $(PROG).1

.SUFFIXES: .gz .o .3 .1 .c .nab .h .inc .txt

.nab.o:
	$(NABC) -c $(NABFLAGS) $(CFLAGS) $(AMBERCFLAGS) $<
	@rm -f $*.c

.txt.1:
	txt2man -s 1 -t $(PROG) $*.txt > $@

all: $(PROG) $(MAN1)

$(PROG): $(OBJS)
	$(NABC) -g -o $@ $(OBJS) $(LDFLAGS) $(AMBERLDFLAGS)

$(OBJS): $(HDRS)

$(MAN1).gz: $(MAN1)
	gzip -cn $(MAN1) > $@

test: $(PROG)
	./$(PROG) -o A12.pdb -t ALPHA -s AAAAAAAAAAAA
	diff A12.pdb save/A12.pdb

beforeinstall: all $(MAN1).gz
	-@mkdir -p $(DESTDIR)bin
	-@mkdir -p $(DESTDIR)man/man1

afterinstall:

install:  beforeinstall
	./install.sh -c -s -o $(INSTALL_USER) -g $(INSTALL_GRP) -m 555 $(PROG) \
		$(DESTDIR)bin
	./install.sh -c -o $(INSTALL_USER) -g $(INSTALL_GRP) -m 444 $(PROG).1.gz \
		$(DESTDIR)man/man1
	@$(MAKE) afterinstall

realclean: clean
	rm -f *~ *core

clean:
	rm -f $(OBJS) $(CLEANFILES) $(PROG) *.gz

.PHONY: all test clean realclean beforeinstall afterinstall install
