MPICC?=mpicc

all: inertia.so

inertia.so: inertia.c
	gcc -Wall -shared -Wl,-soname,inertia -o inertia.so -fPIC inertia.c -lm -lgsl -lgslcblas -fopenmp

clean:
	rm *.so
	rm *~
