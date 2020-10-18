## SET FILEPATHS HERE ##
#CSVPATH = ./../ ## FOR LOCAL TESTING
CSVPATH	= ./../model_20200625/
########################

all:	rms45.cpp
	`root-config --cxx --cflags` -o rms45 rms45.cpp `root-config --glibs`

bpm:	bpm45.cpp
	`root-config --cxx --cflags` -o bpm45 bpm45.cpp `root-config --glibs`

view:	trajectory_display.cpp
	`root-config --cxx --cflags` -o view trajectory_display.cpp `root-config --glibs`

emittance:	emittance45.cpp
	`root-config --cxx --cflags` -o emittance45 emittance45.cpp `root-config --glibs`

emitscan:	emitscan45.cpp
	`root-config --cxx --cflags` -o emitscan emitscan45.cpp `root-config --glibs`

duplicatecsv:	$(CSVPATH)testplane-emittance.csv
	cat $(CSVPATH)testplane-emittance.csv | cut -d "," -f 1-4 > $(CSVPATH)testplane-bpm.csv

clean:	$(CSVPATH)testplane-emittance.csv
	rm $(CSVPATH)testplane-emittance.csv
	rm $(CSVPATH)testplane-bpm.csv
