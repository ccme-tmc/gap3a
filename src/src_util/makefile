include ../make.inc 
.SUFFIXES:      .f90

OBJS=  bessj.o combin.o covsrt.o derfc.o eixb.o factr.o fint.o gammln.o gaulag.o gauleg.o \
       higam.o polynom.o shelsort.o sphbes.o ylm.o veclen.o vecprojlen.o ten_rvctrv.o
LIB=../libutil.a

$(LIB): $(OBJS) 
	$(AR) $(ARFLAGS) $(LIB) $(OBJS) 
	$(RANLIB) $(LIB)

clean:
	rm -f *.o 

cleanall:
	rm -f *.o $(LIB)

%.o : %.f90
	$(FC) $(FFLAGS) -c $< -o $@
#	$(AR) $(ARFLAGS) $(TARGET) $@
.f90.o:
	$(FC) $(FFLAGS) -c $< -o $@ 
#.o.a:
##	$(AR) $(ARFLAGS)  $@ $< 

