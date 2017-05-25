
# executable filenames.
CC=gcc
RM=rm -rf
WC=wc -l
VFL=vflang
GREP=grep -RHni --color

# compilation and linkage flags.
CFLAGS+= -ggdb -O3 -std=gnu99 -Wall -Wextra
LIBS=-lvbnmr -lvfl

# environment for executing compiled tests.
ENV=env $(TESTFLAGS)

# filenames of compiled objects.
OBJ=$(SRC:.c=.o)

# global, default make target.
all: $(BIN)

# test linkage target.
$(BIN): $(OBJ)
	@echo " LD $@"
	@$(CC) $(CFLAGS) $(LDFLAGS) $^ -o $@ $(LIBS)

# source compilation target.
.c.o:
	@echo " CC $^"
	@$(CC) $(CFLAGS) -c $^ -o $@

# intermediate file cleanup target.
clean:
	@echo " CLEAN"
	@$(RM) $(BIN) $(BIN).dSYM $(BIN).out $(BIN).err $(BIN).log $(OBJ) $(DAT)

# full recompilation target.
again: clean all

# interpreted testing target (default).
test: $(BIN).vfl
	@echo " TEST $^"
	@$(VFL) $^

# compiled testing target (alternative).
test-c: $(BIN)
	@echo " TEST $^"
	@$(ENV) ./$(BIN) 1>$(BIN).out 2>$(BIN).err

