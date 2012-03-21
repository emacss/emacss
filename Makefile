CPP=g++
CFLAGS=-g -fast
LM=-lm -std=c99

EMACSS=emacss.cpp params.cpp pre-collapse.cpp input.cpp output.cpp
EMACSS_OBJ=$(EMACSS:.cpp=.o)


all: emacss

emacss: $(EMACSS_OBJ) 
	$(CPP) $(CFLAGS) $(EMACSS_OBJ) $(LM) -o $@

clean:
	rm *.o
 
