SHELL := /bin/bash
CPP := g++
#~/local/gcc/bin/g++

EXE_names := plot_theory plot_data \
  fit_polariz \
  cut_polariz_theory cut_polariz_data \
  weight_dist

OBJ_names := functions inverse fcnp fcnpf TMinuitObj \
  chisq_polariz_base chisq_polariz_op chisq_polariz_tilt chisq_polariz_ediff chisq_polariz_ediff_abs

DIRS := lib bin

LHEF_PATH := lhef

CFLAGS := -std=c++0x -Wall -O3 \
          $(shell root-config --cflags) \
          -I$(LHEF_PATH)/src
GLIBS := -L$(LHEF_PATH)/lib -llhef -lparse -Llib -lwqq \
         $(shell root-config --glibs) -lMinuit \
         -lboost_program_options

BAK = backup`date +%s`.tar.gz

.PHONY: all clean backup

OBJ := $(addprefix lib/, $(addsuffix .o, $(OBJ_names)))
OBJX := $(addprefix lib/, $(addsuffix .o, $(EXE_names)))
EXE := $(addprefix bin/, $(EXE_names))

all: $(DIRS) $(EXE)

# directories rule
$(DIRS):
	@mkdir -p $@

# object rule
$(OBJ): lib/%.o: src/%.cxx src/%.h Makefile
	@echo -e "Compiling \E[0;49;96m"$@"\E[0;0m ... "
	@$(CPP) $(CFLAGS) -c $(filter %.cxx,$^) -o $@

# object rule for those containing main
$(OBJX): lib/%.o: src/%.cxx Makefile
	@echo -e "Compiling \E[0;49;94m"$@"\E[0;0m ... "
	@$(CPP) $(CFLAGS) -c $(filter %.cxx,$^) -o $@

# archive library rule
lib/libwqq.a: $(OBJ) Makefile
	@echo -e "Archiving \E[0;49;93m"$@"\E[0;0m ... "
	@ar rvs $@ $(filter %.o,$^) > /dev/null

# executable rule
$(EXE): bin/%: lib/%.o Makefile
	@echo -e "Linking \E[0;49;92m"$@"\E[0;0m ... "
	@$(CPP) -O3 $(filter %.o,$^) -o $@ $(GLIBS)

# OBJX dependencies
lib/plot_theory.o: src/functions.h src/inverse.h
lib/plot_data.o: $(wildcard $(LHEF_PATH)/src/*.h)
lib/cut_polariz_theory.o: src/functions.h src/inverse.h
lib/cut_polariz_data.o: src/inverse.h src/fcnp.h $(wildcard $(LHEF_PATH)/src/*.h)
lib/weight_dist.o: src/functions.h src/inverse.h $(wildcard $(LHEF_PATH)/src/*.h)

lib/chisq_polariz_base.o: src/TMinuitObj.h
lib/chisq_polariz_op.o: src/functions.h src/TMinuitObj.h src/chisq_polariz_base.h
lib/chisq_polariz_tilt.o: src/functions.h src/TMinuitObj.h src/chisq_polariz_base.h
lib/chisq_polariz_ediff.o: src/functions.h src/TMinuitObj.h src/chisq_polariz_base.h
lib/chisq_polariz_ediff_abs.o: src/functions.h src/TMinuitObj.h src/chisq_polariz_base.h

lib/fit_polariz.o: src/TMinuitObj.h src/chisq_polariz_base.h src/chisq_polariz_op.h src/chisq_polariz_tilt.h src/chisq_polariz_ediff.h src/chisq_polariz_ediff_abs.h $(wildcard $(LHEF_PATH)/src/*.h)

# EXE dependencies
bin/plot_theory: lib/libwqq.a
bin/plot_data: $(LHEF_PATH)/lib/liblhef.a
bin/fit_polariz: $(LHEF_PATH)/lib/liblhef.a lib/libwqq.a
bin/cut_polariz_theory: lib/libwqq.a
bin/cut_polariz_data: $(LHEF_PATH)/lib/liblhef.a lib/libwqq.a
bin/weight_dist: $(LHEF_PATH)/lib/liblhef.a lib/libwqq.a

clean:
	@rm -f $(OBJ) $(OBJX) $(EXE)
	@echo "Cleaned ..."

backup:
	@echo -e "Creating \E[0;49;93m"$(BAK)"\E[0;0m"
	@mkdir -p BAK
	@tar cvzfh BAK/$(BAK) $(wildcard src/*.cxx) $(wildcard src/*.h) $(wildcard scripts/*.sh) Makefile
