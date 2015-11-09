libpath=/usr/local/Thermo-Calc/4.1/
libf=libcalculation-engine-linux-x86_64-gfortran44-7.1.0.7670.so

run_fname=TC_CALC

all:
	clear
	rm -rf mechanics/*.dat kinetics.dat YStress.dat composition.dat *.LOG fitness.dat
	python GA.py

opt:
	clear
	rm -rf mechanics/*.dat kinetics.dat YStress.dat composition.dat *.LOG fitness.dat time.dat
	python GA_opt.py

compilec:
	clear
	rm -rf *.o
	gfortran -ffree-form -fPIC -c anns_Gprime.f95
	gfortran -ffree-form -fPIC -c anns_Gamma.f95
	gcc -c -I $(libpath) -fPIC $(run_fname).c
	gcc $(run_fname).o anns_Gprime.o anns_Gamma.o -lgfortran $(libpath)$(libf) -shared -o $(run_fname).so -lm
	gcc -c -I $(libpath) -fPIC TC_OPT.c
	gcc TC_OPT.o anns_Gprime.o anns_Gamma.o -lgfortran $(libpath)$(libf) -shared -o TC_OPT.so -lm

fortran:
	clear
	f2py -c -m irreverisble irreverisble.f90

