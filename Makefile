libpath=/usr/local/Thermo-Calc/4.1/
libf=libcalculation-engine-linux-x86_64-gfortran44-7.1.0.7670.so

all:
	clear
	rm -rf mechanics/*.dat kinetics.dat YStress.dat composition.dat *.LOG fitness.dat time.dat
	python GA_opt.py

compilec:
	clear
	rm -rf *.o
	gcc -c -I $(libpath) -fPIC TC_OPT.c
	gcc TC_OPT.o $(libpath)$(libf) -shared -o TC_OPT.so -lm

fortran:
	clear
	f2py -c -m irreverisble irreverisble.f90

