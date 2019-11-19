CC = mpiCC
CFLAGS = -O3 -std=c++11

all: laplace

laplace: src/tp2AminMansourElham.cc
	$(CC) src/tp2AminMansourElham.cc $(CFLAGS) -o $@


clean: rm laplace
