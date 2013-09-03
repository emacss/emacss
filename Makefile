CPP=g++
CFLAGS=-g -fast
LM=-lm -std=c99

EMACSS=emacss.cpp input.cpp info.cpp params.cpp cluster.cpp output.cpp
DYN=dynamics.cpp
EMACSS_OBJ=$(EMACSS:.cpp=.o)
DYN_OBJ=$(DYN:.cpp=.o)

all: emacss

emacss: $(EMACSS_OBJ) $(DYN_OBJ) 
	$(CPP) $(CFLAGS) $(EMACSS_OBJ) $(DYN_OBJ) $(LM) -o $@

clean:
	rm *.o
 
emacss.o: emacss.cpp
	$(CPP) -g -c $^ -o $@

input.o: emacss_dir/input.cpp
	$(CPP) -g -c $^ -o $@

info.o: emacss_dir/info.cpp
	$(CPP) -g -c $^ -o $@

params.o: emacss_dir/params.cpp
	$(CPP) -g -c $^ -o $@

cluster.o: emacss_dir/cluster.cpp
	$(CPP) -g -c $^ -o $@

output.o: emacss_dir/output.cpp
	$(CPP) -g -c $^ -o $@

dynamics.o: dynamics/dynamics.cpp
	$(CPP) -g -c $^ -o $@
