############################################################################
#                                                                             #
#                          Generic Makefile for gw                            #
#                                                                             #
#  make           ... generate executable for the COMPLEX sequential version  #
#  make seq       ... generate sequential binary                              #
#  make para      ... generate parallel binary                                #
#  make all       ... generate sequential and parallel build                  #
#  make clean     ... delete unnecessary files                                # 
#  make pack      ... pack all necessary files for distribution               # 
#  make cleanall  ... delete all generated files                              # 
#                                                                             #
###############################################################################
#
#
SRC = "./src/"

default: seq
	cd $(SRC); make; cd ..

all: 
	cd $(SRC); make all; cd ..

seq: 
	cd $(SRC); make seq; cd ..

para:
	cd $(SRC); make para; cd ..

pack:  
	cd $(SRC); make pack; cd ..
         
#..............................................................................
#
#  remove unnecessary files (executable(s) are not removed)
#
clean:
	cd $(SRC); make clean; cd ..

cleanall:
	cd $(SRC); make cleanall; cd ..

