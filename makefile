MA?= mine
path?= obj/
BINDIR?= mod

ifeq ($(MA), mine)
 cc = g++
 FLAGS = -fopenmp
 LIB =   -I ~/code/armadillo/include -DARMA_DONT_USE_WRAPPER -lblas -llapack -larpack -lm
endif
ifeq ($(MA),213)
 cc = g++
 FLAGS = -fopenmp -DMKL_ILP64 -m64 -I${MKLROOT}/include
 LIB =   -I ~/armadillo/include -DARMA_DONT_USE_WRAPPER  -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_ilp64.a ${MKLROOT}/lib/intel64/libmkl_gnu_thread.a ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group -lgomp -lpthread -lm -ldl
endif


All: Dicky cDicky

Dicky: 
	$(cc) src/dicky.cpp -o $@ $(FLAGS)  $(LIB)

cDicky: 
	$(cc) src/dicky_coherant.cpp -o $@ $(FLAGS)  $(LIB)
	
clean: 
	rm Dicky cDicky