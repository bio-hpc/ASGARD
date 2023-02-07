# installation for amberlite:

include ../config.h

install: $(BINDIR)/nab$(SFX) $(BINDIR)/ffgbsa$(SFX) $(BINDIR)/minab$(SFX) \
      $(BINDIR)/mdnab$(SFX)
	./setup.sh $(BINDIR) $(PYTHON)

$(BINDIR)/nab$(SFX): ../nab/nab.c
	cd ../nab && $(MAKE) $@

$(BINDIR)/ffgbsa$(SFX): ffgbsa.nab
	$(BINDIR)/nab -o $(BINDIR)/ffgbsa$(SFX) ffgbsa.nab

$(BINDIR)/minab$(SFX): minab.nab
	$(BINDIR)/nab -o $(BINDIR)/minab$(SFX) minab.nab

$(BINDIR)/mdnab$(SFX): mdnab.nab
	$(BINDIR)/nab -o $(BINDIR)/mdnab$(SFX) mdnab.nab

clean:
	/bin/rm -f *.c

uninstall:
	(cd $(BINDIR) && /bin/rm -f ffgbsa minab mdnab pytleap pymdpbsa pdb4amber)
