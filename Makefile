galsim: galsim.o
	gcc -O3 -o galsim galsim.o
galsim.o: galsim.c
	gcc -O3 -c galsim.c -o galsim.o
clean:
	rm -f galsim galsim.o
