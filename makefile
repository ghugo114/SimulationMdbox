
CC = gcc
CFLAGS = -O2 -g
LDLIBS = -lm

OBJS = box.o nrutil.o rk4.o gnuplot_i.o

all: box

box: $(OBJS)
	$(CC) $(CFLAGS) $(OBJS) -o box $(LDLIBS)

%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -f $(OBJS) box particle_*.dat

.PHONY: all clean

# Explicit dependencies (optional)
gnuplot_i.o: gnuplot_i.c gnuplot_i.h
rk4.o: rk4.c nrutil.h
nrutil.o: nrutil.c nrutil.h
box.o: box.c nrutil.h gnuplot_i.h







