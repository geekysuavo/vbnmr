
# basenames and filenames of headers in vbnmr
H=vbnmr quad vfgp
HDR=$(addsuffix .h,$(H))

# installation directories.
DESTDIR=
PREFIX=$(DESTDIR)/usr/local
INCDIR=$(PREFIX)/include/vbnmr

# non-file targets and accepted filename extensions.
.PHONY: all install
.SUFFIXES: .h

# global, default (nilpotent) make target.
all:

# installation target.
install:
	@echo " INSTALL vbnmr"
	@install -d $(INCDIR)
	@install -d $(INCDIR)/util
	@install $(HDR) $(INCDIR)

