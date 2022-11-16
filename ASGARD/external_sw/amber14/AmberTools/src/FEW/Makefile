SHELL=/bin/sh
include ../config.h

install:
	/bin/sh setup.sh
	(cd $(BINDIR) && /bin/ln -sf ../AmberTools/src/FEW/FEW.pl . )

clean:

uninstall:
	(cd $(BINDIR) && /bin/rm -f FEW.pl)
	/bin/rm -rf $(AMBERHOME)/AmberTools/src/FEW/additional_libs/PerlMol
