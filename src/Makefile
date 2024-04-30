CC = mpicc
CFLAGS = -O3 -g -Wall 


# HDF5 = /usr/local/hdf5

# BIN = .
# LIB = -L$(HDF5)/lib -lhdf5 -lm -lmpi #-lfftw3
# INC = -I. -I$(HDF5)/include 
# SRC = $(wildcard *.c)
# OBJ = $(SRC:.c=.o)

BIN = ../bin
LIB = -lm -lmpi #-lfftw3
INC = -I../include
SRC = $(wildcard *.c)
OBJ = $(SRC:.c=.o)



all: clean main

main:	$(OBJ)
	$(CC) $(CFLAGS) -o $(BIN)/main $(OBJ) $(LIB)

%.o: %.c 
	$(CC) $(CFLAGS) -c $^ -o $@ $(INC) $(LIB)

clean:
	find . -name "*.o"   -exec rm {} \;
	find . -name "*.c%"  -exec rm {} \;
	find . -name "*.bck" -exec rm {} \;
	find . -name "*~"    -exec rm {} \;
	find . -name "\#*"   -exec rm {} \;
	rm -f $(OBJ) main


