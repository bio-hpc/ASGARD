include ../config.h

install::	
	cd prepare_input && $(MAKE) clean && $(MAKE) install
	cd fantasian && $(MAKE) clean && $(MAKE) install
	cp senergy sviol sviol2 $(BINDIR)

clean::	
	cd prepare_input && $(MAKE) clean 
	cd fantasian && $(MAKE) clean

uninstall:
	-(cd prepare_input && $(MAKE) uninstall)
	-(cd fantasian && $(MAKE) uninstall)
	-rm -f $(BINDIR)/senergy
	-rm -f $(BINDIR)/sviol
	-rm -f $(BINDIR)/sviol2
