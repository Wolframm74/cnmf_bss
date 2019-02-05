MEXSUFFIX = mexa64
IX= mexa64
#MATLABHOME = /usr/local/MATLAB/R2013a
MATLABHOME = /usr/local/MATLAB/MATLAB_Production_Server/R2013a
MEX=g++
MEXCXX=echo
CXX=g++

CFLAGS = -fPIC -pthread -DMATLAB_MEX_FILE -ansi -D_GNU_SOURCE -fno-omit-frame-pointer -pthread -O3 -DNDEBUG

CXX=g++

ARMA_DIR=/home/tung/Documents/School/thesis2/mycode/libs

MEX_DIR=../../mex
#I think this will get the name of the folder the Makefile is currently in. 

CURR_DIR=$(MEX_DIR)/debug_compilation_1

INCLUDE = -I $(ARMA_DIR)/armadillo-4.500.1/include/ -I $(ARMA_DIR)/armadillo-4.500.1/mex_interface -I $(CURR_DIR)/include
#-I../openblas_install/include/

INCLUDE_MEX = -I$(MATLABHOME)/extern/include
#

MEXFLAGS =  -shared -Wl,--no-undefined -Wl,-rpath-link,$(MATLABHOME)/bin/glnxa64 -L$(MATLABHOME)/bin/glnxa64 -lmx -lmex -lmat -lm -lut

BIN_OUTPUT_FLAGS =  -Wl,--no-undefined -Wl,-rpath-link,$(MATLABHOME)/bin/glnxa64 -L$(MATLABHOME)/bin/glnxa64 -lmx -lmex -lmat -lm -lut

#

CXXFLAGS = -c -Wall -g -DARMA_BLAS_LONG -std=c++11 -mcmodel=large
#-v

#$(DEBUG) $(FINAL) $(OPT) $(EXTRA_OPT)

LINKFLAGS = -pthread
#-v
#-u mexFunction -Wl,-e,mexFunction

STATIC_LIBS = -lopenblas

#LIBS = -llapack

#try setting LIBS to what I think are MATLAB's internal LAPACK and BLAS libs
LIBS= -lmwlapack -lmwblas -lpthread

SEARCH_DIRS = -L /home/tung/Documents/openblas_install/lib/

UPDATE_MODULES_DIR=$(MEX_DIR)/separate_cpp_files

all: nguyen_2015_bin nguyen_2015_mex

INCLUDE_FILES= $(CURR_DIR)/include/local_inc.hpp $(CURR_DIR)/include/matrix_dimensions.hpp 
#$(CURR_DIR)/include/matrix_objs.hpp 

NGUYEN_2015_MEX_OBJS = Xhat_obj.o W_obj.o Phi_S_obj.o T_obj_mex.o V_obj_mex.o Z_obj.o Y_obj.o update_engines_mex.o algorithm_entry_mex.o helper_funs.o normalize_funs.o plot_and_update_mex.o 
#pre_calculation_modules_outside.o

NGUYEN_2015_BIN_OBJS = algorithm_entry_bin.o 
#Xhat_obj.o T_obj_bin.o V_obj_bin.o update_engines_bin.o helper_funs.o normalize_funs.o plot_and_update_bin.o 
#pre_calculation_modules_outside.o
#W_obj.o Phi_S_obj.o Z_obj.o Y_obj.o 

nguyen_2015_mex: $(NGUYEN_2015_MEX_OBJS)
	$(CXX) $(LINKFLAGS) -o $@.$(MEXSUFFIX) $(NGUYEN_2015_MEX_OBJS) $(LIBS) $(MEXFLAGS) 

algorithm_entry_mex.o: algorithm_entry_hybrid.cpp $(INCLUDE_FILES)
	$(CXX) -DMEX_OUTPUT_FLAG $(CXXFLAGS) $(CFLAGS) -o $@ $< $(INCLUDE) $(INCLUDE_MEX) 

nguyen_2015_bin: $(NGUYEN_2015_BIN_OBJS)
	$(CXX) $(LINKFLAGS) -o $@ $(NGUYEN_2015_BIN_OBJS) $(LIBS) $(BIN_OUTPUT_FLAGS) 

algorithm_entry_bin.o: algorithm_entry_bin.cpp $(INCLUDE_FILES)
	$(CXX) -DBINARY_OUTPUT_FLAG $(CXXFLAGS) $(CFLAGS) -o $@ $< $(INCLUDE) $(INCLUDE_MEX) 

Xhat_obj.o: $(UPDATE_MODULES_DIR)/Xhat_update.cpp $(INCLUDE_FILES)
	$(CXX) $(CXXFLAGS) $(CFLAGS) -o $@ $< $(INCLUDE) $(INCLUDE_MEX) 

W_obj.o: $(UPDATE_MODULES_DIR)/W_update.cpp $(INCLUDE_FILES)
	$(CXX) $(CXXFLAGS) $(CFLAGS) -o $@ $< $(INCLUDE) $(INCLUDE_MEX) 

Phi_S_obj.o: $(UPDATE_MODULES_DIR)/Phi_S_update.cpp $(INCLUDE_FILES)
	$(CXX) $(CXXFLAGS) $(CFLAGS) -o $@ $< $(INCLUDE) $(INCLUDE_MEX) 	

T_obj_mex.o: $(UPDATE_MODULES_DIR)/T_update.cpp $(INCLUDE_FILES)
	$(CXX) -DMEX_OUTPUT_FLAG $(CXXFLAGS) $(CFLAGS) -o $@ $< $(INCLUDE) $(INCLUDE_MEX) 

T_obj_bin.o: $(UPDATE_MODULES_DIR)/T_update.cpp $(INCLUDE_FILES)
	$(CXX) -DBINARY_OUTPUT_FLAG $(CXXFLAGS) $(CFLAGS) -o $@ $< $(INCLUDE) $(INCLUDE_MEX) 

V_obj_mex.o: $(UPDATE_MODULES_DIR)/V_update.cpp $(INCLUDE_FILES)
	$(CXX) -DMEX_OUTPUT_FLAG $(CXXFLAGS) $(CFLAGS) -o $@ $< $(INCLUDE) $(INCLUDE_MEX) 

V_obj_bin.o: $(UPDATE_MODULES_DIR)/V_update.cpp $(INCLUDE_FILES)
	$(CXX) -DBINARY_OUTPUT_FLAG $(CXXFLAGS) $(CFLAGS) -o $@ $< $(INCLUDE) $(INCLUDE_MEX) 

Z_obj.o: $(UPDATE_MODULES_DIR)/Z_update.cpp $(INCLUDE_FILES)
	$(CXX) $(CXXFLAGS) $(CFLAGS) -o $@ $< $(INCLUDE) $(INCLUDE_MEX) 

Y_obj.o: $(UPDATE_MODULES_DIR)/Y_update.cpp $(INCLUDE_FILES)
	$(CXX) $(CXXFLAGS) $(CFLAGS) -o $@ $< $(INCLUDE) $(INCLUDE_MEX) 

update_engines_mex.o: update_engines_hybrid.cpp $(INCLUDE_FILES)
	$(CXX) -DMEX_OUTPUT_FLAG $(CXXFLAGS) $(CFLAGS) -o $@ $< $(INCLUDE) $(INCLUDE_MEX) 

update_engines_bin.o: update_engines_hybrid.cpp $(INCLUDE_FILES)
	$(CXX) -DBINARY_OUTPUT_FLAG $(CXXFLAGS) $(CFLAGS) -o $@ $< $(INCLUDE) $(INCLUDE_MEX) 	

helper_funs.o: $(MEX_DIR)/helper_funs.cpp $(INCLUDE_FILES)
	$(CXX) $(CXXFLAGS) $(CFLAGS) -o $@ $< $(INCLUDE) $(INCLUDE_MEX) 

normalize_funs.o: $(MEX_DIR)/normalize_funs.cpp $(INCLUDE_FILES)
	$(CXX) $(CXXFLAGS) $(CFLAGS) -o $@ $< $(INCLUDE) $(INCLUDE_MEX) 	

plot_and_update_mex.o: plot_and_update_hybrid.cpp $(INCLUDE_FILES)
	$(CXX) -DMEX_OUTPUT_FLAG $(CXXFLAGS) $(CFLAGS) -o $@ $< $(INCLUDE) $(INCLUDE_MEX) 	

plot_and_update_bin.o: plot_and_update_hybrid.cpp $(INCLUDE_FILES)
	$(CXX) -DBINARY_OUTPUT_FLAG $(CXXFLAGS) $(CFLAGS) -o $@ $< $(INCLUDE) $(INCLUDE_MEX) 	


#pre_calculation_modules_outside.o: $(PRE_CALC_DIR)/pre_calculation_modules.cpp $(INCLUDE_FILES)
#	$(CXX) $(CXXFLAGS) $(CFLAGS) -o $@ $< $(INCLUDE) $(INCLUDE_MEX) 	

.PHONY: clean

clean:
	rm -f top_down_threaded_mex.$(MEXSUFFIX)
	rm -f *.o	