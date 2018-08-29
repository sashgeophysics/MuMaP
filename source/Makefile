FC = ftn #ifort
FLAG = -O3  -g -o # -debug inline_debug_info -Bdynamic -auto -o
FLAG2= -c -g #-debug inline_debug_info -Bdynamic -auto -free
#MKLROOT=/opt/intel/composer_xe_2013_sp1.1.106/mkl
MKLPATH=-L$(MKLROOT)/lib/intel64 -lmkl_lapack95_lp64  -Wl,--start-group -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -Wl,--end-group -openmp -lpthread
HIMORB=../data/HI/MORB/
HICARB=../data/HI/carbonate/
CSMORB=../data/CS/MORB/
CSCARB=../data/CS/carbonate/
MKLINCLUDE=-I/opt/intel/composerxe-2011.4.191/mkl/include/intel64


#MKLPATH for ARCHER
MKLPATH= ${MKLROOT}/lib/intel64/libmkl_blas95_lp64.a ${MKLROOT}/lib/intel64/libmkl_lapack95_lp64.a -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a ${MKLROOT}/lib/intel64/libmkl_core.a ${MKLROOT}/lib/intel64/libmkl_sequential.a -Wl,--end-group -lpthread -lm -fPIE
#MKLinclude for ARCHER
MKLINCLUDE= -I${MKLROOT}/include/intel64/lp64 -I${MKLROOT}/include

MKLFLAG= -L$(MKLPATH) $(MKLINCLUDE) 
SOURCES=  hmkc13.f90 LAB.f90 hmt.f90 agius_hi.f90
MODULES =  global.mod microgeodynamics.mod regional.mod
MODOBJ = global.o microgeodynamics.o regional.o
MOD_FILES=$(MODULES:.f90=.o)

AGIUS_SOURCE=agius_hi.f90
AGIUS_OBJECT=$(AGIUS_SOURCE:.f90=.o)
AGIUS_EXEC=agius_hi.exe

HMT_SOURCE=hmt.f90
HMT_OBJECT=$(HMT_SOURCE:.f90=.o)
HMT_EXEC=hmt.exe

HMKC_SOURCE=hmkc13.f90
HMKC_OBJECT=$(HMKC_SOURCE:.f90=.o)
HMKC_EXEC=hmkc13.exe

LAB_SOURCE=LAB.f90
LAB_OBJECT=$(LAB_SOURCE:.f90=.o)
LAB_EXEC=lab.exe
#OBJS=$(SOURCES:.f90=.o)

#####
# No need to edit beyond here unless modules
# are being added or subtracted
#####
####################
# Compile Executable
####################
.PHONY:all HMKC LAB HMT AGIUS

all: HMKC LAB HMT AGIUS

AGIUS:$(AGIUS_EXEC)

HMT:$(HMT_EXEC)

HMKC:$(HMKC_EXEC)

LAB:$(LAB_EXEC)

$(AGIUS_EXEC): object_files 
	$(FC) $(FLAG) $(AGIUS_EXEC) $(AGIUS_OBJECT)  $(MODOBJ) $(MKLFLAG)

$(HMT_EXEC): object_files 
	$(FC) $(FLAG) $(HMT_EXEC) $(HMT_OBJECT)  $(MODOBJ) $(MKLFLAG)


$(HMKC_EXEC): object_files 
	$(FC) $(FLAG) $(HMKC_EXEC) $(HMKC_OBJECT)  $(MODOBJ) $(MKLFLAG)

$(LAB_EXEC):  object_files 
	$(FC) $(FLAG) $(LAB_EXEC) $(LAB_OBJECT)  $(MODOBJ) $(MKLFLAG)
####################
#lab.exe:object_files
#	$(FC) $(FLAG) lab.exe $(OBJS)  $(MODOBJ) $(MKLFLAG)
####################
# Compile Objects
####################
object_files:phase3
	$(FC) $(FLAG2) $(SOURCES)
####################
# Compile rules for modules.  These must be done by hand 
# to account for dependencies
phase3:phase2
	$(FC) $(FLAG2) regional.f90 $(MKLINCLUDE)
####################
phase2:phase1
	$(FC) $(FLAG2) microgeodynamics.f90 $(MKLINCLUDE)
####################
phase1:
	$(FC) $(FLAG2) global.f90 $(MKLINCLUDE)
##############
clean:
	rm -rf *.mod *o *.exe *.dat *~ 
dataclean:
	rm -rf $(HIMORB)*.* $(HICARB)*.* $(CSMORB)*.* $(CSCARB)*.* 



