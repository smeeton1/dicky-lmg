MA?= mine
path?= obj/
BINDIR?= mod

ifeq ($(MA), mine)
 cc = g++
 FLAGS = -fopenmp
 LIB =    -larmadillo -llapack -lm
ifeq ($(MA),213)
 cc = g++
 FLAGS = -fopenmp
 LIB =   -I ~/armadillo/include -DARMA_DONT_USE_WRAPPER  -DMKL_ILP64 -m64 -I${MKLROOT}/include -lm
endif


All: Dicky cDicky

Dicky: 
	$(cc) src/dicky.cpp -o $@ $(FLAGS)  $(LIB)

cDicky: 
	$(cc) src/dicky_coherant.cpp -o $@ $(FLAGS)  $(LIB)
	
clean: 
	rm Dicky cDicky