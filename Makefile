ifndef CC
  CC = gcc	
endif

## Advanced compiling options
#DEBUG = 1
#MAC = 1
#64BITS = 1 #UNCOMMENT TO COMPILE FOR 64-bits

IDIR_BASIC =include/basic
IDIR_ALIGNMENT_TESTS =include/test/alignment

BIN = bin

#main program includes
IDIR_ALIGNMENT =include/alignment
IDIR_BASIC =include/basic

CFLAGS_CUNIT          = -L/usr/local/lib 

ifdef MAC
 MACFLAG    = -fnested-functions
 CC = gcc
#We need to force gcc as nested funcitons are not available on other compilers
endif 

ARCH = 
#ifndef 32_BITS
#   ARCH=-m64
#endif

OPT  		      = $(ARCH) -Wall -O3 $(MACFLAG) -DNO_BOUNDS_CHECK  

ifdef DEBUG
 OPT                  = $(ARCH)  -Wall  -O0 -g3 $(MACFLAG)  -D__DEBUG 
endif



PyroClean: OPT += -DALIGN -g3


ifdef TUNE
 OPT +=  -mtune=$(TUNE)
endif

LINKOPT               = -lm

#we should make this an option, to make or not make a multithreaded compilation. 
ifdef THREADS
 OPT += -DTHREADS=$(THREADS)
 LINKOPT += -lpthread -lrt
endif




CFLAGS_BASIC          = -I$(IDIR_BASIC) 
CFLAGS_ALIGNMENT	  = -I$(IDIR_BASIC) -I$(IDIR_ALIGNMENT)
CFLAGS_ALIGNMENT_TESTS  = -I$(IDIR_ALIGNMENT_TESTS)  -I$(IDIR_ALIGNMENT)  -I$(IDIR_CUNIT)  

LINKOPT += $(OPT) 


ALIGN_OBJ	  =  obj/binary_tree.o obj/peptide.o obj/alignment.o obj/seq.o obj/flags.o obj/file_reader.o obj/logger.o

BASIC_TESTS_OBJ =  obj/element.o obj/flags.o obj/path.o obj/binary_kmer.o obj/seq.o obj/test_binary_kmer.o obj/test_seq.o obj/run_basic_tests.o


ALIGNMENT_TESTS_OBJ = $(ALIGN_OBJ) obj/test_alignment.o obj/run_alignment_tests.o 

all: remove_objects apps

apps: remove_objects PyroClean

tests: remove_objects run_basic_tests 

PyroClean: remove_objects $(ALIGN_OBJ) obj/aligner.o
	mkdir -p $(BIN); $(CC) $(LINKOPT)  -o $(BIN)/PyroClean $(ALIGN_OBJ) obj/aligner.o

run_alignment_tests : clean $(ALIGNMENT_TESTS_OBJ)
	mkdir -p $(BIN); $(CC) $(LINKOPT) $(CFLAGS_CUNIT) -o $(BIN)/run_alignment_tests $(ALIGNMENT_TESTS_OBJ) -lcunit

.PHONY : clean cortex_rain_opts
clean :
	rm -rf $(BIN)
	rm -rf obj

.PHONY : remove_objects	
remove_objects:
	rm -fr obj

#pattern rules

obj/%.o : src/basic/%.c include/basic/%.h
	mkdir -p obj/; $(CC) $(CFLAGS_BASIC) $(OPT) -c $< -o $@


obj/%.o :  src/test/alignment/%.c 
	mkdir -p obj/; $(CC) $(CFLAGS_ALIGNMENT) $(CFLAGS_ALIGNMENT_TESTS) $(OPT) -c $< -o $@
obj/%.o :  src/alignment/%.c 
	mkdir -p obj/; $(CC) $(CFLAGS_ALIGNMENT) $(OPT) -c $< -o $@
obj/%.o : src/test/basic/%.c
	mkdir -p obj/; $(CC) $(CFLAGS_BASIC_TESTS) $(OPT) -c $< -o $@	
obj/%.o : src/test/alignment/%.c 
	mkdir -p obj/; $(CC) $(CFLAGS_ALIGNMENT_TESTS) $(OPT) -c $< -o $@	
