if DIFFTOOL
TESTDIR = tests
endif

# if DOC_ON
# DOCDIR = doc
# endif

# data scripts $(DOCDIR) examples
SUBDIRS = src tools $(TESTDIR)

include_HEADERS = config.h

EXTRA_DIST = autogen.sh

#if DOC_ON
#docs:
#	$(MAKE) -C doc disthtml
#endif

#if DIFFTOOL
#test:
#	$(MAKE) -C tests check
#endif

MAINTAINERCLEANFILES = configure aclocal.m4 acinclude.m4 autom4te.cache \
        Makefile.in config.h.in
