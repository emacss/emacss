CPP=g++
CFLAGS=-g -fast
LM=-lm -std=c99

EMACSS=emacss.cpp input.cpp info.cpp params.cpp cluster.cpp output.cpp
DYN=dynamics.cpp
SE=stellar_evolution.cpp 
EMACSS_OBJ=$(EMACSS:.cpp=.o)
DYN_OBJ=$(DYN:.cpp=.o)
SE_OBJ=$(SE:.cpp=.o)

all: emacss_dev

emacss_dev: $(EMACSS_OBJ) $(DYN_OBJ) $(SE_OBJ)
	$(CPP) -fast $(CFLAGS) $(EMACSS_OBJ) $(DYN_OBJ) $(SE_OBJ) $(LM) -o $@

clean:
	rm *.o
 
emacss.o: emacss.cpp
	$(CPP) -g -c $^ -o $@

input.o: emacss/input.cpp
	$(CPP) -g -c $^ -o $@

info.o: emacss/info.cpp
	$(CPP) -g -c $^ -o $@

params.o: emacss/params.cpp
	$(CPP) -g -c $^ -o $@

cluster.o: emacss/cluster.cpp
	$(CPP) -g -c $^ -o $@

output.o: emacss/output.cpp
	$(CPP) -g -c $^ -o $@

stellar_evolution.o: stevo/stellar_evolution.cpp
	$(CPP) -g -c $^ -o $@

dynamics.o: dynamics/dynamics.cpp
	$(CPP) -g -c $^ -o $@
