
# executable filenames.
CC=gcc
RM=rm -f
WC=wc -l
GREP=grep -RHni --color

# compilation and linkage flags.
CFLAGS=-fPIC -ggdb -O3 -std=gnu99 -I.. -Wall -Wextra
LDFLAGS=-shared
LIBS=-lfftw3 -lvfl -lm

# output library filename.
LIB=libvbnmr.so

# filenames of C sources.
SRC=vbnmr.c quad.c vfgp.c

# filenames of compiled objects.
OBJ=$(SRC:.c=.o)

# non-file targets and accepted filename extensions.
.PHONY: all clean again lines fixme
.SUFFIXES: .c .o

# global, default make target.
all: $(LIB)

# LIB: library linkage make target.
$(LIB): $(OBJ)
	@echo " LD $@"
	@$(CC) $(LDFLAGS) $^ -o $@ $(LIBS)

# source compilation target.
.c.o:
	@echo " CC $^"
	@$(CC) $(CFLAGS) -c $^ -o $@

# intermediate file cleanup target.
clean:
	@echo " CLEAN"
	@$(RM) $(LIB) $(OBJ)

# full recompilation target.
again: clean all

# line-count report target.
lines:
	@echo " WC"
	@$(WC) *.[ch]

# fixme report target.
fixme:
	@echo " FIXME"
	@$(GREP) fixme *.[ch] || echo " None found"

