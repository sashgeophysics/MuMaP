SOURCES=  main.f90 

MAINOBJ = main.o

MODULES =  global.mod microgeodynamics.mod regional.mod

MODOBJ = global.o microgeodynamics.o regional.o

EXEC = mumap.exe

FC = ifort
FLAG = -O3  -g  -debug inline_debug_info -Bdynamic -auto -o
FLAG2= -c -g  -debug inline_debug_info -Bdynamic -auto -free
#MKLPATH for cluster
MKLPATH=-L$(MKLROOT)/lib/intel64 -lmkl_lapack95_lp64  -Wl,--start-group -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -Wl,--end-group -openmp -lpthread
DATAPATH=../data/
#MKLPATH for laptop
#MKLPATH=-L$(MKLROOT)/lib/ia32  -lmkl_blas95 -lmkl_lapack95 -lmkl_intel -lmkl_sequential -lmkl_core -lpthread
#MKLINCLUDE for cluster
MKLINCLUDE=-I/opt/intel/composerxe-2011.4.191/mkl/include/intel64
#MKLINCLUDE for laptop
#MKLINCLUDE=-I$(F95ROOT)/include/ia32 
MKLFLAG= -L$(MKLPATH) $(MKLINCLUDE) 
#####
# No need to edit beyond here unless modules
# are being added or subtracted
#####
OBJS=$(SOURCES:.f90=.o)
MOD_FILES=$(MODULES:.f90=.o)
####################
# Compile Executable
####################
$(EXEC):object_files
	$(FC) $(FLAG) $(EXEC) $(OBJS)  $(MODOBJ) $(MKLFLAG)
####################

####################
# Compile Objects
####################
object_files:phase3
	$(FC) $(FLAG2) $(SOURCES)
####################
# Compile rules for modules.  These must be done by hand 
# to account for dependencies
phase3:phase2
	$(FC) $(FLAG2) regional.f90 $(MKLFLAG)
####################
phase2:phase1
	$(FC) $(FLAG2) microgeodynamics.f90 $(MKLFLAG)
####################
phase1:
	$(FC) $(FLAG2) global.f90 $(MKLFLAG)
##############
clean:
	rm -rf *.mod *o *.exe *.dat *~ $(DATAPATH)*.vtk $(DATAPATH)*.csv $(DATAPATH)*.dat




