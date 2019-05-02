all:	rms45.cpp
	`root-config --cxx --cflags` -o rms45 rms45.cpp `root-config --glibs`

bpm:	bpm45.cpp
	`root-config --cxx --cflags` -o bpm45 bpm45.cpp `root-config --glibs`

view:	trajectory_display.cpp
	`root-config --cxx --cflags` -o view trajectory_display.cpp `root-config --glibs`

emittance:	emittance45.cpp
	`root-config --cxx --cflags` -o emittance45 emittance45.cpp `root-config --glibs`

clean:	bpm45 view
	rm ./bpm45 ./view
