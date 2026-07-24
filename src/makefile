CC      ?= gcc
# -ffp-contract=off: gcc and clang fuse a*b+c into FMA by default where the
# target has it (ARM64, x86-64-v3); that shifts last digits in mat5, mdsm,
# icovops and nlregops against the .ref files.
CFLAGS   = -std=c99 -O2 -ffp-contract=off

ifeq ($(OS),Windows_NT)
  EXE = .exe
endif

SRCS    = $(wildcard *.c)
OBJS    = $(SRCS:.c=.o)
TARGET  = tda$(EXE)

.PHONY: all clean smoke check

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CC) $(CFLAGS) -o $@ $^ -lm

%.o: %.c
	$(CC) $(CFLAGS) -c -o $@ $<

clean:
	rm -f $(OBJS) $(TARGET)

smoke: $(TARGET)
	cd .. && python3 tests/smoke.py src/$(TARGET)
	cd .. && python3 tests/smoke2.py src/$(TARGET)

check: $(TARGET)
	cd .. && python3 tests/check.py src/$(TARGET) examples/exam examples/ehhnew examples/tests examples/coverage
