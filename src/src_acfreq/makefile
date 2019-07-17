include ../make.inc 
.SUFFIXES:      .f90

OBJS= acrgn.o  acpatrd.o calcacfreq.o gaussj.o getsac.o \
      init_c.o  laguer.o mrqcof.o  mrqmin.o  nllsq.o \
      ratfun.o setpatrd.o setsac.o setwpa.o stdesc.o   
LIB=../libacfreq.a

$(LIB): $(OBJS) 
	$(AR) $(ARFLAGS) $(LIB) $(OBJS) 
	$(RANLIB) $(LIB)

clean:
	rm -f *.o *__genmod.* 

cleanall:
	rm -f *.o $(LIB)

%.o : %.f90
	$(FC) $(FFLAGS) -c $< -o $@
#	$(AR) $(ARFLAGS) $(TARGET) $@
.f90.o:
	$(FC) $(FFLAGS) -c $< -o $@ 
#.o.a:
##	$(AR) $(ARFLAGS)  $@ $< 

