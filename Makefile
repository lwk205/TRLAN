LIB2= -Wl,--start-group /opt/intel/mkl/lib/intel64/libmkl_intel_lp64.a /opt/intel/mkl/lib/intel64/libmkl_core.a /opt/intel/mkl/lib/intel64/libmkl_gnu_thread.a -Wl,--end-group -lpthread -lm -ldl -lgfortran

all: clean dlasru.o dbdsqru.o CGS2.o DGEBRDG_4_BISIDE.o DGEBRDG_LP1.o ERR.o RESGKL.o doqds1.o doqds3.o dlartg6.o dlartg7.o dfma0.o RESGKL_MAIN.o fileinput_gkl

%.o: src/%.f90
	gfortran -fopenmp -Wall -O3 -mtune=native -march=native -mcmodel=medium -c -o src/$@ $< 

%.o: lib/%.f
	gfortran -fopenmp -Wall -O3 -mtune=native -march=native -mcmodel=medium  -c -o lib/$@ $<

dlartg7.o:
	gfortran -fopenmp -Wall -mcmodel=medium -O2 -c -o lib/dlartg7.o lib/dlartg7.f

dfma0.o:
	gcc -O3 -fopenmp -Wall -mtune=native -march=native -mcmodel=medium -c -o lib/dfma0.o lib/dfma0.c

fileinput_gkl:
	gcc -fopenmp -Wall -O3 -mtune=native -march=native -mcmodel=medium -o fileinput_gkl fileinput_gkl.c lib/dlasru.o lib/dbdsqru.o src/CGS2.o src/DGEBRDG_4_BISIDE.o src/DGEBRDG_LP1.o src/ERR.o src/RESGKL.o lib/doqds1.o lib/doqds3.o lib/dlartg6.o lib/dlartg7.o lib/dfma0.o src/RESGKL_MAIN.o ${LIB2}

clean:
	rm -rf lib/*.o src/*.o fileinput_gkl

