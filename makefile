#compiler
CXX=g++

#make
MAKC=make
MAKD=tabix-0.2.6/
MAKF=tabix-0.2.6/Makefile

#compiler flags
COPTI= -O3
CDEBG= -g
CWARN= -Wall -Wextra -Wno-sign-compare
CVERP= -D_VERBOSE

#linker flags
LOPTI= -O3
LDEBG= -g
LSTDD= -lm -lpthread -lboost_iostreams -lboost_program_options -lgsl -lblas -lz -I/gscmnt/gc2719/halllab/src/boost_1_57_0/include -L/gscmnt/gc2719/halllab/src/boost_1_57_0/lib 
LTABX= tabix-0.2.6/libtabix.a

#executable file
EFILE= bin/fastQTL

#header files
HFILE= $(shell find src -name *.h)

#source files
CFILE= $(shell find src -name *.cpp)

#source path
VPATH= $(shell for file in `find src -name *.cpp`; do echo $$(dirname $$file); done)

#include path
ISTDP= -Isrc -Itabix-0.2.6

#object files
OFILE= $(shell for file in `find src -name *.cpp`; do echo obj/$$(basename $$file .cpp).o; done)

#default
all: dynamic

#dynamic release
dynamic: CFLAG=$(COPTI) $(CWARN) 
dynamic: LFLAG=$(LOPTI) $(LSTDD)
dynamic: IFLAG=$(ISTDP)
dynamic: $(EFILE)

#dynamic release $(CVERP)
verbose: CFLAG=$(COPTI) $(CWARN) $(CVERP) 
verbose: LFLAG=$(LOPTI) $(LSTDD)
verbose: IFLAG=$(ISTDP)
verbose: $(EFILE)

#debug release
debug: CFLAG=$(CDEBG) $(CWARN) $(CVERP)
debug: LFLAG=$(LDEBG) $(LSTDD)
debug: IFLAG=$(ISTDP)
debug: $(EFILE)

#tabix
tabix-0.2.6/libtabix.a:
	cd $(MAKD) && $(MAKC) && cd ../..

$(EFILE): $(OFILE) $(LTABX)
	$(CXX) $^ -o $@ $(LFLAG)

obj/%.o: %.cpp $(HFILE)
	$(CXX) -o $@ -c $< $(CFLAG) $(IFLAG)

clean: 
	rm -f obj/*.o $(EFILE) && cd $(MAKD) && $(MAKC) clean && cd ../..

oxford:
	cp $(EFILE) ~/bin/.

install:
	cp $(EFILE) /usr/local/bin/.
