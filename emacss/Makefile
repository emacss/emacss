CPP=g++
CFLAGS=-g -fast
LM=-lm -std=c99

EMACSS=emacss.cpp input.cpp info.cpp params.cpp cluster.cpp \
	stellar_evolution.cpp dynamics.cpp
EMACSS_OBJ=$(EMACSS:.cpp=.o)


all: emacss

emacss: $(EMACSS_OBJ) 
	$(CPP) $(CFLAGS) $(EMACSS_OBJ) $(LM) -o $@

clean:
	rm *.o
 
