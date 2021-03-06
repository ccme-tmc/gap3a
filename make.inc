#------------------------------------------------------------------------------#
#                          make.inc for Linux                                  #
#------------------------------------------------------------------------------#


# use ifort+MKL is used for lapack
#  LAPACK= -L/export/opt/intel/compiler/11.0/074/mkl/lib/em64t -lmkl_lapack -lmkl_em64t -lguide
# LAPACK=  -L/home/jiang/intel/mkl/10.1.3.027/lib/em64t -lmkl_lapack  -lmkl_em64t -lguide -lpthread
LAPACK = -L$(MKLROOT)/lib/intel64 -lmkl_intel_lp64 -lmkl_core -lmkl_sequential -lpthread -lm

# use a local version of lapack + gotoblas
# LAPACK= -L./src_numlibs -llapack -lgotoblas
#LAPACK = -L./src_numlibs -llapack -L/export/opt/intel/mkl/10.1.1.019/lib/em64t -lmkl_em64t -lguide

 CC = icc
 CFLAGS = -O2

 #options for debug version 
 FC = ifort -FR
 F77 = ifort -fixed 
 DOPTS =  -O0 -nogen-interfaces -warn interfaces -debug-parameters all -traceback -g -fpe0 -warn all -check uninit -check bounds -ftrapuv  -check noarg_temp_created -warn nounused -DDEBUG
 OOPTS =  -O3 -mp1 -w -prec_div -pc80 -nopad -ip -fpe0 
 FPP= -fpp -DINTEL
 FFLAGS = -r8  $(FPP) $(OOPTS) # -convert big_endian
 LDFLAGS = $(LAPACK) 
 
 MPIFC=mpif90
 MPITAG = mpi
 MPIFFLAGS = -fc=ifort  $(FFLAGS) -DMPI 

# this is to support that different compilers may use different way to include *.mod 
# in ifort, it is just -I, but in SunOS f95 compiler, it uses -M 
 USEMOD = -I

 YACC = bison -y  
 AR  = ar 
 ARFLAGS= crv 
 RANLIB  = ar s

.SUFFIXES:      .F90
.SUFFIXES:      .f90

