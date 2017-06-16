XLIBS = -L/usr/X11R6/lib -lX11
XINCLUDE = -I/usr/X11R6/include/

all: bin colloidal

bin:
	mkdir -p bin
	mkdir -p results
colloidal: source/colloidal.c source/random_numrec.c source/initialization.c source/running.c source/timing.c source/global_variables.c
	gcc source/colloidal.c source/random_numrec.c source/initialization.c source/running.c source/timing.c source/global_variables.c -o build/colloidal -lm -O3
iceplot: source/iceplot.c 
	gcc source/iceplot.c -o plotting/iceplot -lm  $(XLIBS) $(XINCLUDE)
