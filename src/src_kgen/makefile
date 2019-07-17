###############################################################################
#                                                                             #
#  make              ... generate executable                                  #
#  make clean        ... delete unnecessary files                             #
#  make check        ... check static semantics of the FORTRAN routines       #
#  make xref         ... generate a cross reference listing                   #
#                                                                             #
###############################################################################

#------------------------------------------------------------------------------
#
#                          IBM RS6000  AIX 3.2
#
# FC ......... compiler name
# FOPT ....... optimization flags
# FGEN ....... code generation flags (flags not related to optimization)
# LDFLAGS .... linker flags
# LIBS ....... libraries needed to build the executable
# DESTDIR .... directory where the executable should go (incl. trailing '/'!)
# EXECNAME ... name of the resulting executable (without prefixed path!)
#
#------------------------------------------------------------------------------
include ../make.inc
#FC = $(FC) 
#FOPT = -r8 $(FPP) $(DOPTS)
#FOPT     = -O3 -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML
#FGEN     = 
#LDFLAGS  = -Bstatic
LIBS     =
DESTDIR  = ./
EXECNAME = kgen

###############################################################################
###############################################################################

#FFLAGS = $(FGEN) $(FOPT)
EXEC   = $(DESTDIR)$(EXECNAME)

#..............................................................................
#
#  Routines which can be compiled unmodified
#
OBJS1 =	main.o sdefl.o matmm.o sdef.o groups.o bravai.o arbmsh.o basdiv.o \
	reduz.o zuord.o tetdiv.o tetcnt.o ord1.o indexx.o gbass.o def0.o \
	drval.o initi.o initdr.o divisi.o addinv.o

#..............................................................................
#
#  Routines which may require preprocessing before compiling on some machines
#  (Sourcecode-level optimization)
#
PREOBJS =

#..............................................................................
#
OBJS = $(OBJS1) $(PREOBJS)
PRJS = $(OBJS1:.o=.prj) $(PREOBJS:.o=.prj)

#..............................................................................
#
.SUFFIXES:	.prj

#..............................................................................
#
#  Build executable
#
$(EXEC):	$(OBJS)
	$(FC) -o $(EXEC) $(OBJS) $(LDFLAGS) $(LIBS)

#..............................................................................
#
#  All routines depend upon an include file (contains common PARAMETERS)
#

#..............................................................................
#
#  perform static semantics checking (variables used before set, type
#  consistency, correct usage of common blocks, etc.)
#
check:	$(PRJS)
	ftnchek $(PRJS)

#..............................................................................
#
#  remove object files, preprocessed source files and files needed for
#  static semantics checking
#
clean:
	rm  -f $(OBJS) $(PREOBJS:.o=.P) $(OBJS1:.o=.P) $(PRJS) \
		$(EXECNAME).xref $(EXEC)  *__genmod.* 

#..............................................................................
#
#  build a cross reference table (file '$(EXECNAME).xref')
#
xref:
	fxref $(OBJS1:.o=.f) $(PREOBJS:.o=.f) >$(EXECNAME).xref

#..............................................................................
#
#  define rules for generating files needed to check the static semantics
#  of the FORTRAN source code and for generating object files from source files
#
.f.prj:
	ftnchek -project -f77 -declare -sixchar -novice=2 -noextern -library $<
.f.o:
	$(FC) $(FFLAGS) -c $<
.f90.o:
	$(FC) $(FFLAGS) -c $<
