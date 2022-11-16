include ../../config.h

LIBRARY		= pdb

LIBARCH		= lib$(LIBRARY).a
OBJS		= pdbrun.o pdb_read.o pdb_sprntf.o pdb_sscanf.o pdb_write.o ms.o
SRCS		= pdbrun.c pdb_read.c pdb_sprntf.c pdb_sscanf.c pdb_write.c ms.c

install:	$(LIBARCH)

$(LIBARCH):     $(OBJS)
		$(AR) $(LIBARCH) $(OBJS)
		$(RANLIB) $(LIBARCH)
		mv $(LIBARCH) $(LIBDIR)

clean:
	-rm -f $(OBJS) $(LIBDIR)/$(LIBARCH)

