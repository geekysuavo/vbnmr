
# executable filenames.
CD=cd
WC=wc -l
RM=rm -f
TAR=tar czf
GREP=grep -RHni --color

# directory and date for tarball creation.
DIR=vbnmr
DATE=$(shell date +%Y%m%d)

# non-file targets.
.PHONY: all clean again lines fixme dist

# global, default make target.
all:

# intermediate file cleanup target.
clean:
	@$(MAKE) -sC lib clean
	@$(MAKE) -sC tests/ans1d clean
	@$(MAKE) -sC tests/linear clean
	@$(MAKE) -sC tests/mf-slice clean
	@$(MAKE) -sC tests/mf-toy1d clean
	@$(MAKE) -sC tests/mf-toy2d clean
	@$(MAKE) -sC tests/slice clean
	@$(MAKE) -sC tests/toy1d clean

# full recompilation target.
again: clean all

# line-count reporting target.
lines:
	@echo " WC lib"
	@$(WC) lib/*.[ch]
	@echo " WC vbnmr"
	@$(WC) vbnmr/*.[ch]

# fixme statement reporting target.
fixme:
	@echo " FIXME lib"
	@$(GREP) fixme lib/*.[ch] || echo " None found"
	@echo " FIXME vbnmr"
	@$(GREP) fixme vbnmr/*.[ch] || echo " None found"

# tarball creation target.
dist: clean
	@echo " DIST $(DATE)"
	@$(RM) ../$(DIR)-$(DATE).tgz
	@$(CD) .. && $(TAR) $(DIR)-$(DATE).tgz $(DIR)/ && $(CD) $(DIR)

