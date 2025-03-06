galsim: galsim.o
	gcc -Wall -O3 -ffast-math -march=native -fopt-info-vec -o galsim galsim.o -lm
galsim.o: galsim.c
	gcc -O3 -ffast-math -march=native -fopt-info-vec -c galsim.c -o galsim.o
clean:
	rm -f galsim galsim.o