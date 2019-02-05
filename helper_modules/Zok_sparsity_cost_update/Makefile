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

MEX_DIR=../mex
#I think this will get the name of the folder the Makefile is currently in. 

HELPER_MODULES_DIR=$(MEX_DIR)/helper_modules

INCLUDE = -I $(ARMA_DIR)/armadillo-4.500.1/include/ -I $(ARMA_DIR)/armadillo-4.500.1/mex_interface -I $(MEX_DIR)/include -I $(HELPER_MODULES_DIR)/include
#-I../openblas_install/include/

INCLUDE_MEX = -I$(MATLABHOME)/extern/include
#

MEXFLAGS =  -shared -Wl,--no-undefined -Wl,-rpath-link,$(MATLABHOME)/bin/glnxa64 -L$(MATLABHOME)/bin/glnxa64 -lmx -lmex -lmat -lm -lut -leng
#

CXXFLAGS = -c -Wall -g -DARMA_BLAS_LONG -std=c++11 -mcmodel=medium
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

all: nguyen_2015_mex

INCLUDE_FILES= $(MEX_DIR)/include/local_inc.hpp $(MEX_DIR)/include/matrix_dimensions.hpp $(HELPER_MODULES_DIR)/include/Phi_S_update_extra.hpp $(HELPER_MODULES_DIR)/include/Vkn_update_extra.hpp $(HELPER_MODULES_DIR)/include/Zok_update_extra.hpp $(HELPER_MODULES_DIR)/include/Zok_sparsity_cost_update.hpp
#$(HELPER_MODULES_DIR)/include/covariance_cost_update.hpp

#$(HELPER_MODULES_DIR)/include/Phi_S_update_extra.hpp

algorithm_entry_mex.o: algorithm_entry_mex.cpp $(INCLUDE_FILES)
	$(CXX) $(CXXFLAGS) $(CFLAGS) -o $@ $< $(INCLUDE) $(INCLUDE_MEX) 

Xhat_obj.o: $(UPDATE_MODULES_DIR)/Xhat_update.cpp $(INCLUDE_FILES)
	$(CXX) $(CXXFLAGS) $(CFLAGS) -o $@ $< $(INCLUDE) $(INCLUDE_MEX) 

W_obj.o: $(UPDATE_MODULES_DIR)/W_update.cpp $(INCLUDE_FILES)
	$(CXX) $(CXXFLAGS) $(CFLAGS) -o $@ $< $(INCLUDE) $(INCLUDE_MEX) 

Phi_S_obj.o: $(UPDATE_MODULES_DIR)/Phi_S_update.cpp $(INCLUDE_FILES)
	$(CXX) $(CXXFLAGS) $(CFLAGS) -o $@ $< $(INCLUDE) $(INCLUDE_MEX) 	

T_obj.o: $(UPDATE_MODULES_DIR)/T_update.cpp $(INCLUDE_FILES)
	$(CXX) $(CXXFLAGS) $(CFLAGS) -o $@ $< $(INCLUDE) $(INCLUDE_MEX) 

V_obj.o: $(UPDATE_MODULES_DIR)/V_update.cpp $(INCLUDE_FILES)
	$(CXX) $(CXXFLAGS) $(CFLAGS) -o $@ $< $(INCLUDE) $(INCLUDE_MEX) 

Z_obj.o: $(UPDATE_MODULES_DIR)/Z_update.cpp $(INCLUDE_FILES)
	$(CXX) $(CXXFLAGS) $(CFLAGS) -o $@ $< $(INCLUDE) $(INCLUDE_MEX) 

Y_obj.o: $(UPDATE_MODULES_DIR)/Y_update.cpp $(INCLUDE_FILES)
	$(CXX) $(CXXFLAGS) $(CFLAGS) -o $@ $< $(INCLUDE) $(INCLUDE_MEX) 

update_engines.o: update_engines.cpp $(INCLUDE_FILES)
	$(CXX) $(CXXFLAGS) $(CFLAGS) -o $@ $< $(INCLUDE) $(INCLUDE_MEX) 

helper_funs.o: helper_funs.cpp $(INCLUDE_FILES)
	$(CXX) $(CXXFLAGS) $(CFLAGS) -o $@ $< $(INCLUDE) $(INCLUDE_MEX) 

normalize_funs.o: normalize_funs.cpp $(INCLUDE_FILES)
	$(CXX) $(CXXFLAGS) $(CFLAGS) -o $@ $< $(INCLUDE) $(INCLUDE_MEX) 	

plot_and_update.o: plot_and_update.cpp $(INCLUDE_FILES)
	$(CXX) $(CXXFLAGS) $(CFLAGS) -o $@ $< $(INCLUDE) $(INCLUDE_MEX) 	

Z_Y_cluster_obj.o: $(UPDATE_MODULES_DIR)/Z_Y_cluster_update.cpp $(INCLUDE_FILES)
	$(CXX) $(CXXFLAGS) $(CFLAGS) -o $@ $< $(INCLUDE) $(INCLUDE_MEX) 	

VVt_ZtZ_YYt_update_wrappers_obj.o: $(UPDATE_MODULES_DIR)/VVt_ZtZ_YYt_update_wrappers.cpp $(INCLUDE_FILES)
	$(CXX) $(CXXFLAGS) $(CFLAGS) -o $@ $< $(INCLUDE) $(INCLUDE_MEX) 	

V_kn_weighted_scaling_cols.o: $(MEX_DIR)/helper_modules/V_kn_weighted_scaling_cols.cpp $(INCLUDE_FILES)
	$(CXX) $(CXXFLAGS) $(CFLAGS) -o $@ $< $(INCLUDE) $(INCLUDE_MEX) 

##HELPER MODULES##############################################################

#Phi_S_update_extra
Phi_S_update_extra.o: $(HELPER_MODULES_DIR)/Phi_S_update_extra/Phi_S_update_extra.cpp $(INCLUDE_FILES)	
	$(CXX) $(CXXFLAGS) $(CFLAGS) -o $@ $< $(INCLUDE) $(INCLUDE_MEX) 

Vnk_update_Phi_S.o: $(HELPER_MODULES_DIR)/Phi_S_update_extra/Vnk_update_Phi_S.cpp $(INCLUDE_FILES)	
	$(CXX) $(CXXFLAGS) $(CFLAGS) -o $@ $< $(INCLUDE) $(INCLUDE_MEX) 	

#Phi_S_update_extra
Vkn_update_extra.o: $(HELPER_MODULES_DIR)/Vkn_update_extra/Vkn_update_extra.cpp $(INCLUDE_FILES)	
	$(CXX) $(CXXFLAGS) $(CFLAGS) -o $@ $< $(INCLUDE) $(INCLUDE_MEX) 	

#Zok_update_extra
Zok_Ylk_update_extra.o: $(HELPER_MODULES_DIR)/Zok_update_extra/Zok_Ylk_update_extra.cpp $(INCLUDE_FILES)	
	$(CXX) $(CXXFLAGS) $(CFLAGS) -o $@ $< $(INCLUDE) $(INCLUDE_MEX) 		

Zok_Zol_update_extra.o: $(HELPER_MODULES_DIR)/Zok_update_extra/Zok_Zol_update_extra.cpp $(INCLUDE_FILES)	
	$(CXX) $(CXXFLAGS) $(CFLAGS) -o $@ $< $(INCLUDE) $(INCLUDE_MEX) 			

#covariance_cost_update

# COVARIANCE_COST_UPDATE_OBJS = covariance_common.o covariance_cost_update_engine.o covariance_Y_update.o covariance_Z_update.o 

# covariance_common.o: $(HELPER_MODULES_DIR)/covariance_cost_update/covariance_common.cpp $(INCLUDE_FILES)	
# 	$(CXX) $(CXXFLAGS) $(CFLAGS) -o $@ $< $(INCLUDE) $(INCLUDE_MEX) 			

# covariance_cost_update_engine.o: $(HELPER_MODULES_DIR)/covariance_cost_update/covariance_cost_update_engine.cpp $(INCLUDE_FILES)	
# 	$(CXX) $(CXXFLAGS) $(CFLAGS) -o $@ $< $(INCLUDE) $(INCLUDE_MEX) 			

# covariance_Y_update.o: $(HELPER_MODULES_DIR)/covariance_cost_update/covariance_Y_update.cpp $(INCLUDE_FILES)	
# 	$(CXX) $(CXXFLAGS) $(CFLAGS) -o $@ $< $(INCLUDE) $(INCLUDE_MEX) 			

# covariance_Z_update.o: $(HELPER_MODULES_DIR)/covariance_cost_update/covariance_Z_update.cpp $(INCLUDE_FILES)	
# 	$(CXX) $(CXXFLAGS) $(CFLAGS) -o $@ $< $(INCLUDE) $(INCLUDE_MEX) 						

#target output file: nguyen_2015_mex

#Zok_sparsity_cost_update

ZOK_SPARSITY_COST_UPDATE_OBJS = Zok_sparsity_common.o Zok_sparsity_Ylk_update.o Zok_sparsity_Zol_update.o

Zok_sparsity_common.o: $(HELPER_MODULES_DIR)/Zok_sparsity_cost_update/Zok_sparsity_common.cpp $(INCLUDE_FILES)	
	$(CXX) $(CXXFLAGS) $(CFLAGS) -o $@ $< $(INCLUDE) $(INCLUDE_MEX) 		

Zok_sparsity_Ylk_update.o: $(HELPER_MODULES_DIR)/Zok_sparsity_cost_update/Zok_sparsity_Ylk_update.cpp $(INCLUDE_FILES)	
	$(CXX) $(CXXFLAGS) $(CFLAGS) -o $@ $< $(INCLUDE) $(INCLUDE_MEX) 			

Zok_sparsity_Zol_update.o: $(HELPER_MODULES_DIR)/Zok_sparsity_cost_update/Zok_sparsity_Zol_update.cpp $(INCLUDE_FILES)	
	$(CXX) $(CXXFLAGS) $(CFLAGS) -o $@ $< $(INCLUDE) $(INCLUDE_MEX) 				

NGUYEN_2015_MEX_OBJS = Xhat_obj.o W_obj.o Phi_S_obj.o T_obj.o V_obj.o Z_obj.o Y_obj.o update_engines.o algorithm_entry_mex.o helper_funs.o normalize_funs.o plot_and_update.o Z_Y_cluster_obj.o VVt_ZtZ_YYt_update_wrappers_obj.o V_kn_weighted_scaling_cols.o Phi_S_update_extra.o Vnk_update_Phi_S.o Vkn_update_extra.o Zok_Ylk_update_extra.o Zok_Zol_update_extra.o $(ZOK_SPARSITY_COST_UPDATE_OBJS)
#$(COVARIANCE_COST_UPDATE_OBJS)

nguyen_2015_mex: $(NGUYEN_2015_MEX_OBJS)
	$(CXX) $(LINKFLAGS) -o $@.$(MEXSUFFIX) $(NGUYEN_2015_MEX_OBJS) $(LIBS) $(MEXFLAGS) 

.PHONY: clean

clean:
	rm -f top_down_threaded_mex.$(MEXSUFFIX)
	rm -f *.o	