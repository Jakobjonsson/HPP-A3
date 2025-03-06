galsim: galsim_openmp.o
	gcc -Wall -O3 -fopenmp -ffast-math -march=native -fopt-info-omp -fopt-info-vec -o galsim galsim_openmp.o -lm
galsim_openmp.o: galsim_openmp.c
	gcc -O3 -fopenmp -ffast-math -march=native -fopt-info-vec -c galsim_openmp.c -o galsim_openmp.o
clean:
	rm -f galsim galsim_openmp.o
