include ../config.h

# Makefile for building and installing AMBER 10 Amgauss

ESPTDOC = ESSfile.pdf Gfchk.pdf Glib.pdf Aprmtop.pdf Aprmcrd.pdf
FIXEDPWD=$(subst /,\/,$(PWD))
PROG=amgauss
SHELL=/bin/sh
VPATH=lib/ESPT

install: fixlib chmod exelink 

chmod:
	-chmod a+x $(PROG).pl

clean: uninstall
	-rm -f doc/*.pdf

uninstall:
	-rm -f $(BINDIR)/$(PROG)


# create documentation in pdf format
# from embedded perldoc.
doc: $(PROG).pdf
	mv $(PROG).pdf doc/

# link to the development version of ESPT
esptlink:
	ln -s ../../ESPT lib/ESPT

# make a link in AMBER's binaries directory
exelink:
	ln -sf ../src/$(PROG)/$(PROG).pl $(BINDIR)/$(PROG)

# make ESPT documentation
extradoc: $(ESPTDOC)
	mv *.pdf doc/

# update the use lib statement to reflect where amgauss actually resides
fixlib:
	perl -pe 's/^use\slib\s".*$(PROG)\/lib";/use lib "$(FIXEDPWD)\/lib";/g' $(PROG).plx > amgauss.pl  

%.ps : %.pl
	pod2man $< | groff -man -Tps -t >$@

%.ps : %.pm
	pod2man $< | groff -man -Tps -t >$@

%.pdf : %.ps
	ps2pdf $<
