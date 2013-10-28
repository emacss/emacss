CPP=g++
CFLAGS = -fPIC -c -g
LM=-lm -std=c99
LIB = -shared -W1,-soname,lib

EMACSS=emacss.cpp input.cpp info.cpp params.cpp cluster.cpp output.cpp
DYN=dynamics.cpp
SE=stellar_evolution.cpp 
EMACSS_OBJ=$(EMACSS:.cpp=.o)
DYN_OBJ=$(DYN:.cpp=.o)
SE_OBJ=$(SE:.cpp=.o)

all: libemacss

libemacss: $(EMACSS_OBJ) $(DYN_OBJ) $(SE_OBJ)
	$(CPP) $(LIB) -o $@.so $(EMACSS_OBJ) $(DYN_OBJ) $(SE_OBJ)
	
emacss_exe: $(EMACSS_OBJ) $(DYN_OBJ) $(SE_OBJ)
	$(CPP) -fast -g $(EMACSS_OBJ) $(DYN_OBJ) $(SE_OBJ) $(LM) -o $@

clean:
	rm *.o
 
emacss.o: emacss.cpp
	$(CPP) $(CFLAGS) $^ -o $@

input.o: emacss/input.cpp
	$(CPP) $(CFLAGS) $^ -o $@

info.o: emacss/info.cpp
	$(CPP) $(CFLAGS) $^ -o $@

params.o: emacss/params.cpp
	$(CPP) $(CFLAGS) $^ -o $@

cluster.o: emacss/cluster.cpp
	$(CPP) $(CFLAGS) $^ -o $@

output.o: emacss/output.cpp
	$(CPP) $(CFLAGS) $^ -o $@

stellar_evolution.o: stevo/stellar_evolution.cpp
	$(CPP) $(CFLAGS) $^ -o $@

dynamics.o: dynamics/dynamics.cpp
	$(CPP) $(CFLAGS) $^ -o $@
