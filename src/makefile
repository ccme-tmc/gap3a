############################################################################
#                                                                             #
#                          Generic Makefile for gw                            #
#                                                                             #
#  make           ... generate executable for the REAL sequential version     #
#  make real      ... generate executable for the REAL sequential version     #
#  make complex   ... generate executable for the COMPLEX sequential version  #
#  make rp        ... generate executable for the REAL parallel version       #
#  make cp        ... generate executable for the COMPLEX parallel version    #
#  make all       ... generate real, complex, rp, cp                          #
#  make last      ... repeat latest make (real, complex, rp, or cp)           #
#  make clean     ... delete unnecessary files                                # 
#  make bzint     ... generate the bzint library                              #
#                                                                             #
#  CAUTION: This makefile must not be called via the '-f' parameter of 'make'.#
#           (Don't use 'make -f Makefile.aix32' !!)                           #
#           (Always (symbolic) link the appropriate Makefile to 'Makefile'!!) #
#                                                                             #
###############################################################################
#
# FC ........... compiler name
# MPIFC ......... compiler name for parallel compilation
# FGEN ......... code generation flags (flags not related to optimization)
# LDFLAGS ...... linker flags
# R_LIBS ....... libraries needed to build the REAL executable
# C_LIBS ....... libraries needed to build the COMPLEX executable
# RP_LIBS ...... libraries needed to build the REAL parallel executable
# CP_LIBS ...... libraries needed to build the COMPLEX parallel executable
# DESTDIR ...... dir. where the executable should go (without trailing '/'!)
# R_EXECNAME ... name of the resulting REAL executable (without prefixed path!)
# C_EXECNAME ... name of the resulting COMPLEX executable ( -"- )
# RP_EXECNAME ... name of the resulting REAL parallel executable ( -"- )
# CP_EXECNAME ... name of the resulting COMPLEX parallel executable ( -"- )
#
include ./make.inc

vernum = 3
verdiv = a
version = $(vernum)$(verdiv)
PACKNAME = gap$(version)
S_EXECNAME = gap$(version).x
P_EXECNAME = gap$(version)-mpi.x  
DESTDIR = .
S_EXEC = $(DESTDIR)/$(S_EXECNAME)
P_EXEC = $(DESTDIR)/$(P_EXECNAME)

SCRIPTS = gap${version}_init gap$(vernum)_analy gap$(vernum)_bandanaly gap$(vernum)_gwnvf gap$(vernum)_initwann gap$(vernum)_kgen gap$(vernum)_newin1 gap$(vernum)_save gap_utils.py gap$(vernum)_w2kset gap$(vernum)_w2w gap$(vernum)_scaninp 

BZDIR= ./src_libbzint
LIBBZ = ./libbzint.a

PARSERDIR = ./src_libparser
LIBPARSER = ./libparser.a

ACFREQDIR = ./src_acfreq
LIBACFREQ = ./libacfreq.a

UTILDIR = ./src_util
LIBUTIL = ./libutil.a

KGENDIR = ./src_kgen  
DMFTPRDIR = ./src_dmftproj

LIBS = $(LIBBZ) $(LIBPARSER) $(LIBACFREQ) $(LIBUTIL)
INCLUDES = $(USEMOD)$(MODDIR) $(USEMOD)$(PARSERDIR) $(USEMOD)$(BZDIR)

MODDIR = ./src_modules
SEQMODS = modules.a 
MPIMODS = modules_mpi.a

###############################################################################

DOC_FILE =    ./doc/gw.tex
DOCTOT_FILE = ./doc/gwtot.tex

#..............................................................................
#

# files that use MPI library 
MPIOBJS =  

# files used by both sequential and parallel code 
SP_OBJS = calcemac_nlf.o calceps.o calceps_2d.o calcgw0.o calcscgw0.o calcescgw.o \
          calcselfx.o calcselfc.o  calc_mwm3.o calcselfcmn.o calcselfxmn.o \
		  coul_coef.o \
          io_cleanup.o io_initfiles.o io_sxcmn.o io_utils.o main.o \
          readingw.o calc_sxcmn.o calcselfcmn_mwm.o \
          scgw_check_conv.o scgw_set_restart.o scgw_update_bande.o \
          task_acfd.o task_crpa.o task_emac.o task_eps.o task_gw.o task_gwsc.o task_gw2wann.o task_mwm.o  

OBJS =  $(SP_OBJS) \
        bandanaly.o \
        bz_calcqdepw.o bz_setiksym.o bz_setkiw.o bz_setkmesh.o \
        calcecrpa.o calceffmass.o calceqp.o   \
        calcexhf.o calcgcoef.o calchead.o calcjlam.o  calcloabc.o\
        calcmmatvv.o calcmmatcc.o calcmmatcv.o \
        calcminm.o calcminc.o calcmicm.o calcmicc.o calcmommat.o  calcmwm.o \
        calcplasmon.o calcqpwf_dmat.o calcrwf.o  coul_calcsing.o      \
        calcvhmn.o  calcvorbnn.o calcudot.o calcu.o    \
        calc_spectfun_nk.o calc_spectfun.o \
        coul_barcq0.o coul_barc.o coul_setvm0.o coul_setvm1.o coul_mpwmix.o coul_mpwipw.o coul_setev.o coul_strcnst.o coul_wmix0.o  \
        coul_barcq0_cutoff.o coul_barc_cutoff.o coul_setvm0_cutoff.o coul_setvm1_cutoff.o coul_strcnst_cutoff.o \
        crpa_calcmill.o crpa_calcumat.o crpa_calcvmat.o crpa_check_unitary.o crpa_output.o \
        crpa_readpln.o crpa_readw2w.o crpa_setmask.o \
        diagipw.o \
        expand_evec.o \
        fourintp.o freq_convl.o freq_intpl_ac.o freq_set_fgrid.o \
        genauxf.o gendjmm.o getcgcoef.o getdjmm.o get_minm.o      \
        int1ipw.o intipw.o intstipw.o \
        io_eps.o io_eqp.o io_mwm.o io_mwm3.o io_vxcmn.o \
        k2cart.o kfrac2cart.o kp_k0index.o kp_interp.o kip_qpeintp.o kip_readeqp.o kip_readenk.o \
        latgen.o locdef.o lohns.o  \
        mb_calcrlamint.o mb_calcs3r.o mb_setumix.o mb_setuprod.o momradintc.o   momradintv.o     \
        orthog_corewf.o outwin.o           \
        prep_ang_int.o \
        radmesh.o readvector.o rint13.o rot1tosph.o rotate.o rotdef.o     \
        scgw_herm_sxc.o \
        set_lapwcoef.o \
        set_ipw.o      set_ksref.o    set_lapwlo.o  \
        set_minm.o     set_mixbasis.o set_struct.o  \
        setkk0index.o  setkpg.o       setlocmixind.o \
        setrindex.o   symmom.o \
        task_acont.o task_chkbz.o \
        task_coul.o \
        task_ldau.o task_nvf.o  task_ppgw.o \
        testbarc.o  threeylm.o \
        trans_minm.o   trans_vector.o trans_vmn.o trilinear.o \
        unrotate.o update_vector.o \
        w2k_calcuxcu.o w2k_genvectorf.o  w2k_readenergy.o  w2k_readstruct.o  \
        w2k_readvxc.o  w2k_writeklist.o  w2k_calcvxcmn.o   w2k_calcvxcnn.o  w2k_readeband.o w2k_readcore.o \
        w2k_readin1.o  w2k_readmommat.o  w2k_readvsp.o     w2k_writevector.o w2k_vxccub.o \
        writeqgen.o    writesym.o  write_sxc_wann.o   

#
#..............................................................................
#
#  Object files for sequential and parallel versions
#
S_OBJS = $(OBJS) $(SEQMODS) $(LIBS)  
P_OBJS = $(OBJS) $(MPIMODS) $(LIBS)

#..............................................................................
#
#  Build executable (either REAL or COMPLEX versions)
#

# Documented files


default: seq 

all: seq para  

seq:  keep_seq_files kgen dmftpr bin
	$(MAKE) $(S_EXEC) 
	cp $(S_EXEC) ../bin/$(S_EXECNAME)

para: keep_para_files kgen dmftpr bin
	$(MAKE) $(P_EXEC) FC=$(MPIFC) FFLAGS='$(MPIFFLAGS)'
	cp $(P_EXEC) ../bin/$(P_EXECNAME)

bin: $(SCRIPTS)
	mkdir -p ../bin/
	cp $(SCRIPTS) ../bin/

$(S_EXEC): SEQMOD $(LIBS) $(OBJS) 
	$(FC) $(FFLAGS) -o $(S_EXEC) $(S_OBJS) $(LDFLAGS) 

$(P_EXEC): MPIMOD $(LIBS) $(OBJS) $(MPIOBJS) 
	$(FC) $(FFLAGS) -o $(P_EXEC) $(P_OBJS) $(LDFLAGS) 

kgen: 
	cd $(KGENDIR); make; cd ..

dmftpr:
	cd $(DMFTPRDIR); make; cd ..

SEQMOD:   
	cd $(MODDIR); make seq; cd ..

MPIMOD: 
	cd $(MODDIR); make para; cd ..

$(LIBBZ): $(BZDIR)/*.f90  
	cd $(BZDIR); make; cd .. 

$(LIBACFREQ): $(ACFREQDIR)/*.f90 
	cd $(ACFREQDIR); make; cd ..

$(LIBUTIL): $(UTILDIR)/*.f90
	cd $(UTILDIR); make; cd ..

$(LIBPARSER): $(PARSERDIR)/*.f90 $(PARSERDIR)/*.c  
	cd $(PARSERDIR); make ; cd ..

keep_seq_files:
	if [ -f .para ]; then \
	rm -f .para $(SP_OBJS)  ; \
	fi
	touch .seq

keep_para_files:
	if [ -f .seq ]; then \
	rm -f .seq $(SP_OBJS) ; \
	fi
	touch .para

#..............................................................................
#FILES = $(OBJS:.o=.f90)
#
#  generate distribution file
#
pack: 
	rm -f ../$(PACKNAME).tgz
	mkdir -p ../$(PACKNAME)
	$(MAKE) prepsrc
	$(MAKE) prepdoc
	$(MAKE) preptest
	$(MAKE) prepexamples
	tar -czvf ../$(PACKNAME).tgz ../$(PACKNAME)
	rm -rf ../$(PACKNAME)

prepsrc:
	mkdir -p ../$(PACKNAME)/src
	cp -rf ../makefile ../README.md  ../$(PACKNAME)/
	cp -rf *.f90 makefile make.inc-* $(SCRIPTS) ../$(PACKNAME)/src
	mkdir -p ../$(PACKNAME)/src/$(KGENDIR); cd $(KGENDIR); cp *.f makefile README ../../$(PACKNAME)/src/$(KGENDIR); cd .. 
	mkdir -p ../$(PACKNAME)/src/$(DMFTPRDIR); cd $(DMFTPRDIR); cp *.f makefile ../../$(PACKNAME)/src/$(DMFTPRDIR); cd .. 
	mkdir -p ../$(PACKNAME)/src/$(MODDIR); cd $(MODDIR); cp *.f90 makefile ../../$(PACKNAME)/src/$(MODDIR); cd ..
	mkdir -p ../$(PACKNAME)/src/$(BZDIR); cd $(BZDIR); cp *.f90 *.f makefile ../../$(PACKNAME)/src/$(BZDIR); cd ..
	mkdir -p ../$(PACKNAME)/src/$(ACFREQDIR); cd $(ACFREQDIR); cp *.f90 *.f makefile ../../$(PACKNAME)/src/$(ACFREQDIR); cd ..
	mkdir -p ../$(PACKNAME)/src/$(UTILDIR); cd $(UTILDIR); cp *.f90 *.f makefile ../../$(PACKNAME)/src/$(UTILDIR); cd .. 
	mkdir -p ../$(PACKNAME)/src/$(PARSERDIR); cd $(PARSERDIR); cp *.f90 *.c *.h *.in *.y makefile ../../$(PACKNAME)/src/$(PARSERDIR); cd .. 

prepdoc:
	mkdir -p ../$(PACKNAME)/doc
	cp ../doc/README* ../$(PACKNAME)/doc

preptest:
	mkdir -p ../$(PACKNAME)/test
	cp ../test/* ../$(PACKNAME)/test


prepexamples:
	mkdir -p ../$(PACKNAME)/examples
	cp ../examples/* ../$(PACKNAME)/examples

#..............................................................................
#
#  remove unnecessary files (executable(s) are not removed)
#
clean:
	rm -f *.o *__genmod.* .seq .para  

cleanall:
	make clean
	rm -f ../$(PACKNAME).tgz
	rm -f *.mod
	rm -f $(S_EXEC) $(P_EXEC) $(LIBS)  *.swp  $(SEQMODS) $(MPIMODS)
	cd ../bin; rm -f gap*;cd ../src/
	cd $(MODDIR); make clean; cd .. 
	cd $(ACFREQDIR); make clean; cd ..
	cd $(BZDIR); make clean; cd .. 
	cd $(PARSERDIR); make clean; cd ..
	cd $(UTILDIR); make clean; cd ..
	cd $(KGENDIR); make clean; cd ..
	cd $(DMFTPRDIR); make cleanall; cd .. 

cleansvn:
	rm -rf .svn */.svn 
#..............................................................................
#
#  define inference rules to generate object files from source files
#  be chosen.)
#
.f.o:
	$(FC) $(INCLUDES) $(FFLAGS)  -c $<
.c.o:
	$(CC) -c $<
.f90.o:
	$(FC) $(INCLUDES) $(FFLAGS) -c $<
