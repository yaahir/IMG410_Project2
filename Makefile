CC=gcc
CFLAGS=-std=c11 -Wall -Wextra -pedantic -O2
LDFLAGS=-lm

all: v3test

v3test: v3test.o v3math.o
	$(CC) $(CFLAGS) -o v3test v3test.o v3math.o $(LDFLAGS)

v3test.o: v3test.c v3math.h
	$(CC) $(CFLAGS) -c v3test.c

v3math.o: v3math.c v3math.h
	$(CC) $(CFLAGS) -c v3math.c

clean:
	rm -f *.o v3test
