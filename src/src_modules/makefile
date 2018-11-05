include ../make.inc 
.SUFFIXES:      .f90

OBJS= task.o constants.o lebedev_laikov.o reallocate.o recipvec.o \
	  acfd.o angintegrals.o bands.o barcoul.o bzinteg.o anisotropy.o  \
	  core.o crpa.o dielmat.o eigenvec.o fouri.o freq.o \
      hfexch.o kmeshintp.o  kpoints.o lapack.o  lapwlo.o minmmat.o mixbasis.o  modmpi.o  modules.o  mommat.o    \
      param.o prodfun.o radwf.o rspevec.o selfenergy.o struk.o  \
      xcpot.o   
      
# objects that are needed in both sequential and parallel version 
SP_OBJS= modmpi.o dielmat.o


MODFILES = $(OBJS:.o=.f90) $(SP_OBJS:.o=.f90)
DIST_FILE = modules.tgz
SEQMODS = ../modules.a
MPIMODS = ../modules_mpi.a


seq: keep_seq_files 
	$(MAKE) $(SEQMODS) 

para:  keep_para_files 
	$(MAKE) $(MPIMODS) FC=$(MPIFC) FFLAGS='$(MPIFFLAGS)'



$(SEQMODS): $(OBJS) 
	$(AR) $(ARFLAGS) $(SEQMODS) $(OBJS) 
	$(RANLIB) $(SEQMODS)

$(MPIMODS): $(OBJS)  
	$(AR) $(ARFLAGS) $(MPIMODS) $(OBJS)
	$(RANLIB) $(MPIMODS)

keep_seq_files:
	if [ -f .para ]; then \
	  rm -f .para $(SP_OBJS) $(SP_OBJS:.o=.mod); \
	fi
	touch .seq

keep_para_files:
	if [ -f .seq ]; then \
	  rm -f .seq $(SP_OBJS) $(SP_OBJS:.o=.mod); \
	fi
	touch .para

bands.o : kpoints.o 
barcoul.o: constants.o kpoints.o mixbasis.o struk.o
bzinteg.o : bands.o core.o freq.o kpoints.o struk.o constants.o  lebedev_laikov.o
core.o: radwf.o
dielmat.o: constants.o mixbasis.o recipvec.o 
eigenvec.o: task.o bands.o struk.o lapwlo.o recipvec.o
minmmat.o : task.o
mixbasis.o: recipvec.o reallocate.o struk.o  core.o 
modmpi.o : constants.o
moments.o: lapwlo.o
radwf.o :  lapwlo.o bands.o 
recipvec.o : kpoints.o
selfenergy.o : freq.o 
xcpot.o: selfenergy.o bands.o kpoints.o struk.o
struk.o : reallocate.o

clean:
	rm -f $(OBJS) *.mod .para .seq
#
#  generate distribution file
#
dist:
	tar -cPzvf $(DIST_FILE) makefile $(MODFILES)

%.o : %.f90
	$(FC) $(FFLAGS) -c $< -o $@
#	$(AR) $(ARFLAGS) $(TARGET) $@
.f90.o:
	$(FC) $(FFLAGS) -c $< -o $@ 
#.o.a:
##	$(AR) $(ARFLAGS)  $@ $< 

