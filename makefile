FFLAGS = -O -fno-emulate-complex -ffast-math
FC=mpif77
#LIBS=-L/home/dow/bin/Helas -ldhelas
LIBS=
prog=wz
OBJ=$(prog).o kern.o config.o cuts.o smear.o integration.o reno25.o alphas2.o random.o winput.o lib.o intlib.o boxlib.o dfunc.o libweak.o libqed.o mrst2004qed.o dfunc_new.o graphs.o Cteq6Pdf-2008.o
OBJF=$(prog).f kern.f config.f cuts.f smear.f integration.f reno25.f alphas2.f random.f winput.f lib.f intlib.f boxlib.f dfunc.f libweak.f libqed.f mrst2004qed.f dfunc_new.f graphs.f Cteq6Pdf-2008.f

$(prog).x: $(OBJ) 
	$(FC) $(FFLAGS) -o $@ $(OBJ) $(LIBS)
