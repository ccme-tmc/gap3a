include ../make.inc
# Add file types to the known list of suffixes (file extensions)
.SUFFIXES: .f .exe .obj

# The fortran compiler and its options:
PROG=./dmftproj.x
#FC = ifort
#FFLAGS = -g -traceback
#FFLAGS = -O2
#LFLAGS = $(LAPACK) 

# Define object file names
objs=modules.o dmftproj.o readcomline.o set_ang_trans.o setsym.o \
     set_rotloc.o timeinv.o read_k_list.o set_projections.o orthogonal.o density.o \
     rot_projectmat.o symmetrize_mat.o rot_dens.o orthogonal_wannier.o \
     outputqmc.o outbwin.o outband.o outputcrpa.o set_projections_silicon.o newbasis.o 
     
.f.o:
	- $(F77) $(FFLAGS) -c $*.f

$(PROG): $(objs)
	-$(F77) $(objs) $(LDFLAGS) -o $(PROG)

clean:
	rm -f *.o *.mod *__genmod.* 

cleanall: 
	rm -f $(PROG) *.o *.mod

