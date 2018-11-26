include ../make.inc 
.SUFFIXES:      .f90

OBJS= acfd.o lapwlo.o modules.o angintegrals.o lapack.o fouri.o \
	  freq.o task.o constants.o kpoints.o lebedev_laikov.o \
	  reallocate.o hfexch.o kmeshintp.o \
	  modmpi.o  recipvec.o bands.o \
      radwf.o  eigenvec.o \
      barcoul.o mixbasis.o  anisotropy.o  bzinteg.o \
      core.o crpa.o dielmat.o\
      minmmat.o mommat.o \
      param.o prodfun.o rspevec.o selfenergy.o struk.o  \
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
	$(AR) $(ARFLAGS) $(SEQMODS) $(OBJS);
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

clean:
	rm -f $(OBJS) *.mod .para .seq
#
#  generate distribution file
#
dist:
	tar -cPzvf $(DIST_FILE) makefile $(MODFILES)

# independent modules and modules directly dependent on them
#    acfd.o
#    constants.o
#    freq.o
#    kpoints.o
bands.o : kpoints.o 
recipvec.o : kpoints.o
#    lapwlo.o
#    reallocate.o
rotmat.o : reallocate.o
struk.o : reallocate.o
#    task.o
modmpi.o : constants.o task.o

# ======================================
radwf.o :  lapwlo.o \
           bands.o struk.o
eigenvec.o: task.o constants.o kpoints.o lapwlo.o \
            bands.o struk.o recipvec.o

# ======================================
core.o : \
         \
         radwf.o

mommat.o: lapwlo.o constants.o \
          struk.o bands.o \
          \
          core.o
mixbasis.o: reallocate.o lapwlo.o \
            struk.o bands.o recipvec.o \
            \
            core.o
bzinteg.o : lebedev_laikov.o freq.o kpoints.o \
            struk.o bands.o \
            \
            core.o

# ======================================
minmmat.o : task.o bands.o kpoints.o \
            \
            \
            core.o \
            mixbasis.o 
dielmat.o: task.o constants.o kpoints.o \
           bands.o \
           \
           core.o \
           mixbasis.o
anisotropy.o: task.o constants.o \
              modmpi.o \
			  \
              \
		      mixbasis.o bzinteg.o
crpa.o: freq.o kpoints.o \
        bands.o \
        \
        \
        mixbasis.o
barcoul.o: constants.o kpoints.o \
           struk.o \
           \
           \
           mixbasis.o 

# ======================================
selfenergy.o : freq.o task.o kpoints.o \
          bands.o \
          \
          \
          \
          crpa.o
selfenergy-nov15.o : freq.o task.o kpoints.o \
          	         bands.o \
          	         \
          	         \
          	         \
	                 crpa.o
xcpot.o: kpoints.o \
         bands.o struk.o \
         \
         \
         \
         \
         selfenergy.o
# ======================================

%.o : %.f90
	$(FC) $(FFLAGS) -c $< -o $@
#	$(AR) $(ARFLAGS) $(TARGET) $@
.f90.o:
	$(FC) $(FFLAGS) -c $< -o $@ 
#.o.a:
##	$(AR) $(ARFLAGS)  $@ $< 

