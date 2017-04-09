
CC=gcc
RM=rm -rf
WC=wc -l
GREP=grep -RHni --color

CFLAGS=-ggdb -O3 -std=gnu99 -I../.. -Wall -Wextra
LDFLAGS=-L../../lib
LIBS=-lvbnmr -lvfl

ENV=env DYLD_LIBRARY_PATH=../../lib $(TESTFLAGS)

OBJ=$(SRC:.c=.o)

all: $(BIN)

$(BIN): $(OBJ)
	@echo " LD $@"
	@$(CC) $(CFLAGS) $(LDFLAGS) $^ -o $@ $(LIBS)

.c.o:
	@echo " CC $^"
	@$(CC) $(CFLAGS) -c $^ -o $@

clean:
	@echo " CLEAN"
	@$(RM) $(BIN) $(BIN).dSYM $(BIN).out $(BIN).err $(OBJ) $(DAT)

again: clean all

test: $(BIN)
	@echo " TEST $^"
	@$(ENV) ./$(BIN) 1>$(BIN).out 2>$(BIN).err

lines:
	@echo " WC"
	@$(WC) *.[ch]

fixme:
	@echo " FIXME"
	@$(GREP) fixme *.[ch] || echo " None found"

