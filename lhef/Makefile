OBJ_names := particle LHEFReader LHEF2Wqq

LIB_names := lhef

DIRS := lib

CFLAGS := -Wall -O3 $(shell root-config --cflags)

BAK = backup`date +%s`.tar.gz

.PHONY: all clean backup

OBJ := $(addprefix lib/, $(addsuffix .o, $(OBJ_names)))
LIB := $(addprefix lib/lib, $(addsuffix .a, $(LIB_names)))

all: $(DIRS) $(LIB)

# directories rule
$(DIRS):
	@mkdir -p $@

# object rule
$(OBJ): lib/%.o: src/%.cxx src/%.h Makefile
	@echo "Compiling "$@" ... "
	@g++ $(CFLAGS) -c $(filter %.cxx,$^) -o $@

# executable rule
$(LIB): $(OBJ) Makefile
	@echo "Creating "$@" ... "
	@ar rvs $@ $(filter %.o,$^)

lib/LHEFReader.o: src/LHEF.h
lib/LHEF2Wqq.o: src/LHEF.h src/LHEFReader.h

clean:
	@rm -f $(OBJ) $(LIB)
	@echo "Cleaned ..."

backup:
	@echo "Creating "$(BAK)
	@mkdir -p BAK
	@tar cvzfh BAK/$(BAK) $(wildcard src/*.c*) $(wildcard src/*.h) Makefile
