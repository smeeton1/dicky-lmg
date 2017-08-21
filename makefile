MA?= mine
path?= obj/
BINDIR?= mod

ifeq ($(MA), mine)
 cc = g++
 FLAGS = -O2 -pthread
 LIB =    -larmadillo -llapack -lm
endif


All: Dicky cDicky

Dicky: Dicky
	$(cc) src/dicky.cpp -o $@ $(FLAGS)  $(LIB)

cDicky: cDicky
	$(cc) src/dicky_coherant.cpp -o $@ $(FLAGS)  $(LIB)
	
clean: 
	rm Dicky cDicky