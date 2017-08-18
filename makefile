MA?= mine
path?= obj/
BINDIR?= mod

ifeq ($(MA), mine)
 FLINKER = g++
 FLAGS = -02 -pthread
 LIB =   -I /usr/lib/openmpi/include/ -I /usr/include -I /usr/local/lib/ -larmadillo -llapack -lm
else
 FLINKER = ftn
 FLAGS = -w
 LIB =  -L$(FFTW_LIB) -lfftw3  
endif



# -lpacklib  -lpawlib -lpacklib -lmathlib -lgraflib -lgrafX11 \-lkernlib
 
include $(SLEPC_DIR)/conf/slepc_common 
PETSC_INCLUDE = -I$(PETSC_DIR)/include -I /usr/include/lam/
PETSC_ARCH_INCLUDE = -I$(PETSC_DIR)/$(PETSC_ARCH)/include


Dicky: 	$(FLINKER) $(FLAGS) -o $@ src/dicky.c++ $(LIB)
	
cDicky: $(FLINKER) $(FLAGS) -o $@ src/dicky_coherant.c++ $(LIB)