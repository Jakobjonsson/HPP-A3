galsim: galsim_pthreads.o
	gcc -Wall -O3 -ffast-math -march=native -fopt-info-vec -o galsim galsim_pthreads.o -lm -lpthread
galsim_pthreads.o: galsim_pthreads.c
	gcc -O3 -ffast-math -march=native -fopt-info-vec -c galsim_pthreads.c -o galsim_pthreads.o -lpthread
clean:
	rm -f galsim galsim_pthreads.o
