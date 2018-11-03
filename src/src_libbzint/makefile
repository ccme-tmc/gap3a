include ../make.inc
.SUFFIXES:	.f90 

DESTDIR  = ./
LIBNAME = ../libbzint.a
DIST_FILE =$(DESTDIR)/src_libbzint.tgz
LIB   = $(DESTDIR)$(LIBNAME)

#
#  Files defining modules
#
MODOBJS = bzint.o order.o kgen_internals.o tetra_internal.o polyhedron.o 

#
#  Files for generation of k-point grids
#
KGENFILES = kgen.o kgenq.o kqgen.o

#
# Subroutines for use of kgenfiles only
#

KGENSUBS =  coorskp.o divisi.o jget.o idkp.o jset.o  sym2int.o     \
	    tetinit.o tgen.o tgenq.o
            
#
#  Files for BZ integrations and convolutions
#

INTFILES = cart2int.o cart2int0.o tetcw.o tetcorecw.o tetiw.o tetiwsurf.o tetraint.o        \
		   tetraqint.o calcdos.o calcidosve.o dostet.o idos.o        \
		   fermi.o

#
# Subroutine for use of intfiles only
#

INTSUBS = bloechlcor.o convw1t.o convw.o convcorew.o dos1t.o         \
	  generictetra.o genericfunf.o genericprism.o intdos1t.o       \
	  intweight1t.o intw.o intwsurf.o ksurf.o relnodes.o         \
	  setnodes.o sortnodes.o sortsurf.o surfnodes.o tlinap.o     \
	  unrepnodes.o sorteq.o stweight_imag.o stweight_real.o        \
          stweight_itaylor.o stweight_rtaylor.o reduk.o redusym.o  \
          stweight_numeric.o
      
#
#  Files that can be called independently (and are also used by the others programs)
#
FREEFILES = cartezian.o gbass.o intern.o rbass.o factorize.o         \
            average_degen_convw.o average_degen_weights.o 

#
# All files
#

OBJS = $(KGENFILES) $(INTFILES) $(FREEFILES) $(KGENSUBS) $(INTSUBS)
#

#

#..............................................................................
#
#  Build executable
#
$(LIB):	$(OBJS) $(MODOBJS)
	$(AR) $(ARFLAGS) $(LIB) $(OBJS) $(MODOBJS)
	$(RANLIB) $(LIB)

#..............................................................................
#
#  All routines depend upon an include file (contains common PARAMETERS)
#

$(OBJS): $(MODOBJS) makefile  

#..............................................................................
#
#  remove object files, preprocessed source files and files needed for
#  static semantics checking
#
clean:  
	rm  -f $(OBJS) $(MODOBJS) *.a *.mod  *__genmod.* 

#..............................................................................
#
#  generate latex documentation
doc:
	cp tetra.tex $(LIB:.a=.tex)
	@echo '\newpage' >> $(LIB:.a=.tex)
	@echo '\section{Routine/Function Prologues}' >> $(LIB:.a=.tex)
	protex -b -n $(MODOBJS) >> $(LIB:.a=.tex)
	@echo '\newpage' >> $(LIB:.a=.tex)
	@echo '\subsection{Subroutines for the generation of k-point grids}' >> $(LIB:.a=.tex)
	protex -b -n $(KGENFILES) >> $(LIB:.a=.tex)
	@echo '\newpage' >> $(LIB:.a=.tex)
	@echo '\subsection{Subroutines used by the k-point generation subroutines}' >> $(LIB:.a=.tex)
	protex -b -n $(KGENSUBS) >> $(LIB:.a=.tex)
	@echo '\newpage' >> $(LIB:.a=.tex)
	@echo '\subsection{Subroutines for the Brillouin Zone integration}' >> $(LIB:.a=.tex)
	protex -b -n $(INTFILES) >> $(LIB:.a=.tex)
	@echo '\newpage' >> $(LIB:.a=.tex)
	@echo '\subsection{Subroutines used by the integration subroutines}' >> $(LIB:.a=.tex)
	protex -b -n $(INTSUBS) >> $(LIB:.a=.tex)
	@echo '\newpage' >> $(LIB:.a=.tex)
	@echo '\subsection{Other useful subroutines that can also be called separatedly}' >> $(LIB:.a=.tex)
	protex -b -n $(FREEFILES) >> $(LIB:.a=.tex)
	@echo '\newpage' >> $(LIB:.a=.tex)
	@echo '\input{tests}' >> $(LIB:.a=.tex)
	@echo '\end{document}' >> $(LIB:.a=.tex)

#..............................................................................
#
#  generate distribution file
#
dist:   
	tar -cvf $(DIST_FILE) makefile $(FILES)
        
#.f90.o: 
#	$(FC) $(FFLAGS) -c $<
       
%.o : %.f90
	$(FC) $(FFLAGS) -c $< -o $@
 
f90.o:
	$(FC) $(FFLAGS) -c $< -o $@        
