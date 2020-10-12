#The compiler:
FC=gfortran #-O -g -Wall -Wextra -Warray-temporaries -Wconversion -fimplicit-none -fbacktrace -ffree-line-length-0 -fcheck=all 
CPP=/lib/cpp -ffreestanding -C -ansi
CC=gcc -DHAVE_CONFIG_H -I. -I.. -g -O0 -MT aux.o -MD -MP -MF .deps/aux.Tpo -c 
FCFLAGS=   -O3 -I/opt/etsf/include/ -lxcf90 -lxc 
CPPFLAGS= -I.
LDFLAGS= -O3 -lblas -llapack path/libnag.a -L/opt/etsf/lib -lxcf90 -lxc 

#Executable
PROGRAM=hippo

OBJECT= common.o auxil.o main.o vary_occ.o \
       vary_occ_NAG.o construct_f.o  non_loc_eff_pot.o  \
       vary_occ_theta.o OEP.o aux_OEP.o vary_occ_theta_n.o  \
       readgamess.o functionals.o total_energy.o vary_orbs_direct.o \
       grid_module.o TotalSpinSq.o FermiDirac.o Model_OEP.o \
       Effective_orbital.o invtest_analyt.o \

#"make" for all:
all: $(PROGRAM)

hippo: $(OBJECT)
	$(FC) $(FCFLAGS) -o $@ $^ $(LDFLAGS)

%.o: %.f90
	$(FC) $(FCFLAGS) -c $<

clean:
	rm -f *.o *.mod *.MOD
