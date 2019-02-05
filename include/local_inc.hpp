#ifndef LOCAL_INC_H_
#define LOCAL_INC_H_

#include "stdint.h" /*http://stackoverflow.com/questions/13362084/difference-between-uint32-and-uint32-t*/
#include "math.h" /*http://www.cplusplus.com/reference/cmath/*/
#include "assert.h"
#include <queue>
#include <vector>
#include <armaMex.hpp>
#include "matrix.h"
#include <pthread.h>
#include <stdlib.h>
#include <semaphore.h>
/*#include "local_helper_templates.hpp"*/

/*#define MEX_UNIT_TESTING*/
/*#define REGULAR_USE*/

/*#define DEBUG*/
/*#define DEBUG1*/
/*#define DEBUG34534*/

#define COMPUTE_L_SWITCH
#define N_ITERATIONS_TV 40

#include "matrix_dimensions.hpp"

/*typedef struct {
	int Xhat_queue_array[N_static*M_static][2];
} Xhat_struct_t; 

typedef struct {
	int H_queue_array[F_static*L_static][2];
} H_struct_t; 

typedef struct {
	int Phi_H_queue_array[L_static*M_static][2];
} Phi_H_struct_t; 

typedef struct {
	int S_queue_array[N_static*L_static][2];
} S_struct_t; */

typedef struct {

	uint16_t M;
	uint16_t F;
	uint16_t N;
	uint16_t K;
	uint16_t L;
	uint16_t O; 

	mxArray *T_fk_mxArray_p;
	mxArray *Y_lk_mxArray_p;
	mxArray *Z_ol_mxArray_p;
	mxArray *V_nk_mxArray_p;


	cx_cube* Xhat_fnm_p;
	cube* Xhat_low_fnm_p;
	cx_cube* E_conj_fnm_p;
	cx_cube* Xtilde_fnm_p;
	cube* W_fom_p;
	cx_cube* W_fom_cx_p;
	mat* Z_ol_p;
	mat* Y_lk_p; 
	mat* T_fk_p; 
	mat* V_nk_p; 
	cube* Phi_W_fom_p; 
	/*cube* Phi_S_fnk_p;*/ 
	cube* Phi_S_nkf_p; 
	cube* Phi_S_fkn_p; 		
	cx_cube* expj_Phi_S_nkf_p;
	cx_cube* expj_Phi_S_fkn_p;
	cx_cube* expj_Phi_W_fom_p;
	cx_mat* U_N_mat_fnk_p; 		
/*	cx_cube* S_fnl_p;
	cube::fixed<F_static, N_static, K_static> Phi_S_fnkl_arr[L_static];*/	/*pass in the array name*/

	/*std::queue<int*> queue_intpz;*/
	std::queue<int> queue_intz; 
	int work_queue_size; 
	int which_update;
	char update_string[80];

	FILE* Xhat_pFile; 

#ifdef BOTTOM_UP_MEX_FLAG

	/*Bottom Up / Pairwise Merge related variables: -----------------------------------------------------------------------------------------------------------------------------------------------------------------*/
	int Lcurrent;
	mat* frob_norm_results_vector_p;
	//uint16_t bu_queue_num_elements;		

	/*Thread related info*/
	/*bool thread_identity;*/ /*false: thread1, true: thread2*/
	int mid_point;
	int frob_final_size;
	int num_frob_worker_threads;

#endif	

} arg_struct_t; 

typedef struct {
	int thread_iter; 
	arg_struct_t* argStruct_p;
} thread_arg_t;

#ifdef BOTTOM_UP_MEX_FLAG
extern void W_callMATLAB_wrapper(void);
#endif

#ifdef ENABLE_MATLAB_ENGINE_PLOTTING
extern Engine* mlEngine_p;
#endif

/*Preprocessing serial functions*/
extern void W_0th_sum_preprocess(mat* Z_ol_p, mat* Y_lk_p);
extern void V_0th_sum_preprocess(mat* V_nk_p);

/*Normalizing functions*/
extern void normalize_V_nk_frob_norm(mat* V_nk_p, mat* T_fk_p);
extern void normalize_Z_ol_frob_norm(mat* Z_ol_p);
extern void normalize_Y_lk_compensate_Z_Lstatic_eq_Ltarget(mat* Y_lk_p, mat* Z_ol_p);
extern void normalize_Z_ok_compensate_TV(mat* Z_ol_p, mat* Y_lk_p, mat* T_fk_p, mat* V_nk_p);
extern void normalize_Z_ol(mat* Z_ol_p, mat* Y_lk_p);
extern void normalize_V_nk(mat* V_nk_p, mat* T_fk_p);
extern void normalize_V_nk_2(mat* V_nk_p, mat* Y_lk_p);
extern void normalize_Y_lk(mat* Y_lk_p, mat* T_fk_p);
extern void normalize_Y_lk_2(mat* Y_lk_p, mat* V_nk_p);
extern void normalize_T_fk(mat* T_fk_p, mat* V_nk_p);
extern void normalize_T_fk_arth_mean(mat* T_fk_p, mat* V_nk_p);
extern void normalize_T_fk_arth_mean_2(mat* T_fk_p, mat* Y_lk_p);

/*Normalizing global variables*/
extern double accum_T_fk;
extern double accum_V_nk;
extern double accum_Z_ol;
extern double accum_W_fom;
extern cube* W_fom_local_p;

/*Also need function prototypes for 'update_engines.cpp' module */
extern void update_engine(arg_struct_t* argStruct_p, int which_update);

/*plot_and_update functions*/
extern void create_lhs_matrices(mxArray *plhs[], arg_struct_t* argStruct_p);
extern void plot_and_update_Xhat(mxArray *plhs[], arg_struct_t* argStruct_p, int n);
extern void plot_and_update_Tfk(mxArray *plhs[], const mxArray *prhs[], arg_struct_t* argStruct_p, int n);
extern void plot_and_update_Vnk(mxArray *plhs[], arg_struct_t* argStruct_p, int n);
extern void plot_and_update_Zol(mxArray *plhs[], const mxArray *prhs[], arg_struct_t* argStruct_p, int n);
extern void plot_and_update_Ylk(mxArray *plhs[], arg_struct_t* argStruct_p, int n);
extern void plot_and_update_Phi_S(mxArray *plhs[], arg_struct_t* argStruct_p, int n);
extern void plot_and_update_Wfom(mxArray *plhs[], arg_struct_t* argStruct_p, int n);

/*Helper functions*/
extern void populate_queue_wrt_one_index(arg_struct_t* argStruct_p, int dim1_len);
extern void populate_queue_wrt_pair_indices(arg_struct_t* argStruct_p, int dim1_len, int dim2_len); 
extern bool calculate_queue_index(arg_struct_t* argStruct_p, int* int_ptr, bool* last_element_p); 
extern bool calculate_pair_queue_indices(arg_struct_t* argStruct_p, int* index1_p, int* index2_p, bool* last_element_p); 
extern void populate_W_fom_cx(cube* W_fom_p, cx_cube* W_fom_cx_p, cube* Phi_W_fom_p); 
extern void clear_queue_intz(std::queue<int> &queue_intz);
extern void threads_exit_and_signal(int num_threads);

extern void projfunc_wrapper_Zol(mat* Z_ol_p, mxArray *Z_ol_mxArray_p);
extern void projfunc_wrapper_Ylk(mat* Y_lk_p, mxArray *Y_lk_mxArray_p);
extern void projfunc_wrapper_Tfk(mat* T_fk_p, mxArray *T_fk_mxArray_p);

extern void orthogonalize_Ylk_square_wrapper(mat* Y_lk_p, mxArray *Y_lk_mxArray_p);
extern void orthogonalize_Vnk_wrapper(mat* V_nk_p, mxArray *V_nk_mxArray_p);

extern void YYt_update_wrapper(mat* Y_lk_p, mxArray *Y_lk_mxArray_p, mxArray *plhs[]);
extern void ZtZ_update_wrapper(mat* Z_ol_p, mxArray *Z_ol_mxArray_p, mxArray *plhs[]);
extern void VVt_update_wrapper(mat* V_nk_p, mxArray *V_nk_mxArray_p, mxArray *plhs[]);
extern void TtT_update_wrapper(mat* T_fk_p, mxArray *T_fk_mxArray_p, mxArray *plhs[]); 

extern void TDOA_update_module1_wrapper(mxArray *plhs[], const mxArray *prhs[], arg_struct_t* argStruct_p);
extern void Plot_significant_L_wrapper(mxArray *plhs[], const mxArray *prhs[], arg_struct_t* argStruct_p);
extern void Plot_significant_L_wrapper_3(mxArray *plhs[], const mxArray *prhs[], arg_struct_t* argStruct_p);

/*Inside W_update*/
extern void magnitude_square_rooted_processing_W_wrapper(arg_struct_t* argStruct_p, mxArray* plhs[]);

/*Inside Phi_S_update*/
extern void Phi_S_L1_norm_V_checker(cx_cube* expj_Phi_S_nkf_p, cx_cube* expj_Phi_S_fkn_p, mat* V_nk_p);

extern void V_nk_MA_filter_rows(mat* V_nk_p);
void populate_window_vec(void);

/*"Cluster" update functions*/
void Z_cluster_update_wrapper(arg_struct_t* argStruct_p);
void Y_cluster_update_wrapper(arg_struct_t* argStruct_p);

/*Function Signatures*/
extern void* Xhat_primary_auxfun_start(void* arg);
extern void* Z_primary_auxfun_start(void* arg);
extern void* Y_primary_auxfun_start(void* arg);
extern void* T_primary_auxfun_start(void* arg);
extern void* V_primary_auxfun_start(void* arg);
extern void* W_primary_auxfun_start(void* arg);
extern void* Phi_S_primary_auxfun_start(void* arg);

extern void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]); /*Keep this here. Can't remember exactly why I need(ed) this, but I think it may have been a result of the gcc compiler not knowing which function was the entry point of the .mex object file.*/

#define MAX_NUM_SUMS 3
#define NUM_WORKER_THREADS 4
/*#define NUM_WORKER_THREADS 8*/

#define mexPrintf(...) /*Define this if you want to define mexPrintf() as being NULL/compiler replaces all instances of it with empty space*/
							/*Undefine it (Comment it out) if you want it to be defined as normal */

#define fprintf(...) /*Define this if you want to define mexPrintf() as being NULL/compiler replaces all instances of it with empty space*/
							/*Undefine it (Comment it out) if you want it to be defined as normal */

/*Semaphores*/
extern sem_t check_queue_global_sem; 	/*Double check you init'd this*/
extern sem_t threads_asleep_sem; 		/*Double check you init'd this*/
extern sem_t queue_sem;
extern sem_t sleep_sem;
extern sem_t threads_while_condition_0th_sum_sem;
extern sem_t threads_while_condition_1st_sum_sem;
extern sem_t threads_while_condition_2nd_sum_sem; 
extern sem_t last_element_sleep_sem;
extern sem_t work_actually_computed_sem;
extern sem_t checkpoint_sem;
extern sem_t wait_till_queue_populated_sem;
/*extern sem_t threads_while_condition_3rd_sum_sem;*/

/*Global variables*/
extern bool threads_while_condition_0th_sum_flag;
extern bool threads_while_condition_1st_sum_flag;
extern bool threads_while_condition_2nd_sum_flag; 
/*extern bool threads_while_condition_3rd_sum_flag;*/
extern bool last_element_found_global_flag; 		/*No corresponding semaphore for this obj*/
extern bool check_queue_global_flag; 
extern bool last_element_sleep_flag;
extern bool work_actually_completed_flag; 
extern bool Xhat_pop_MKF_for_Phi_S_flag;

extern int threads_asleep_ctr; 
extern int work_actually_computed_ctr;

/*Debug variables*/
#ifdef DEBUG
extern sem_t print_sem;
#endif 

#ifdef DEBUG1
extern bool Xhat_after_V_flag;
#endif

/*Global variables*/
/*extern mat::fixed<F_static, O_static> W_fo_m[M_static];		
extern cx_mat::fixed<F_static, O_static> W_fo_m_cx[M_static];	*/

extern cube::fixed<M_static, N_static, F_static> Intensor_Xhat_low_mnf;
extern cx_cube::fixed<M_static, N_static, F_static> Intensor_E_conj_mnf;

/*extern cx_cube::fixed<F_static, K_static, M_static> Xhat_outtensor_cx_FKM; 
extern cube::fixed<F_static, K_static, M_static> Xhat_outtensor_real_FKM;*/

extern cx_cube::fixed<M_static, K_static, F_static> Xhat_outtensor_cx_MKF; 
extern cube::fixed<M_static, K_static, F_static> Xhat_outtensor_real_MKF;

extern cx_cube::fixed<F_static, L_static, M_static> Xhat_outtensor_cx_FLM; 
extern cube::fixed<F_static, L_static, M_static> Xhat_outtensor_real_FLM;

extern colvec::fixed<N_static> ones_col_Nx1; /*This needs to be set to .ones() just once. ever, before either Z update or W update */
extern colvec::fixed<O_static> ones_col_Ox1;	/*Inside normalize_funs.cpp*/
extern colvec::fixed<L_static> ones_col_Lx1;	/*Inside normalize_funs.cpp*/
extern mat::fixed<F_static, O_static> ones_mat_FxO;	/*Inside W_update.cpp*/
extern rowvec::fixed<K_static> ones_row_1xK; /*Z update*/
extern colvec::fixed<F_static> ones_col_Fx1;	/*Inside normalize_funs.cpp*/
/*extern cx_cube::fixed<F_static, K_static, N_static> reference_phase_tensor_cx_FKN;
extern cx_cube::fixed<N_static, K_static, F_static> reference_phase_tensor_cx_NKF;
extern cube::fixed<F_static, K_static, N_static> ones_tensor_real_FKN;
extern cube::fixed<N_static, K_static, F_static> ones_tensor_real_NKF;
extern cx_cube::fixed<F_static, K_static, N_static> phase_shift_tensor_cx_FKN;
extern cx_cube::fixed<N_static, K_static, F_static> phase_shift_tensor_cx_NKF;*/
extern mat::fixed<N_static, K_static> ones_mat_NxK;	/*Inside W_update.cpp*/

extern mat::fixed<K_static, L_static> Z_inmat_real_Ykl;

/*helper fun: threads_exit_and_signal() related stuff*/
extern int threads_exit_ctr; /*ctr*/
extern sem_t threads_exit_lock_sem; /*lock semaphore*/
extern sem_t threads_exit_sleep_sem; /*sleep semaphore*/

extern int which_update_helper;
extern FILE* Xhat_pFile_helper;

extern bool last_thread_to_sleep_flag_array[MAX_NUM_SUMS][NUM_WORKER_THREADS];
extern bool last_thread_to_sleep_flag_local[2][M_static];

extern bool N_ITERATIONS_TV_bool;

extern mxArray* scalar_global_p[3];
extern colvec::fixed<1> scalar_global[3];

extern double L1_global;
extern double L2_global;

/*extern cube::fixed<F_static, N_static, M_static> X_tilde_fnm_abs_cube;
extern cube::fixed<F_static, N_static, M_static> X_tilde_fnm_binary_mask_cube;*/

#include "Phi_S_update_extra.hpp"
#include "Vkn_update_extra.hpp"
#include "Tfk_update_extra.hpp"
#include "Zok_update_extra.hpp"
#include "Ylk_update_extra.hpp"
#include "Zol_update_extra.hpp"
/*#include "covariance_cost_update.hpp"*/
/*#include "Zok_sparsity_cost_update.hpp"*/
/*#include "hfk_Af_LS_estimator_update.hpp"*/
#include "cluster_pullup.hpp"
#include "k_means_clustering.hpp"
#include "tdoa_update_mx.hpp"
#include "tdoa_update_m1.hpp"
#include "tdoa_update_m2.hpp"
#include "tdoa_update_m3.hpp"
#include "Fevotte_cost_function_m1.hpp"
#include "Fevotte_cost_function_m2_helper_funs.hpp"
#include "Fevotte_EM_update_m1.hpp"
#include "beamforming_update_m1.hpp"
#include "complex_argument_costfun_m1.hpp"
#include "complex_argument_costfun_m2.hpp"
#include "complex_argument_costfun_m3.hpp"
#include "complex_argument_costfun_m4.hpp"
#include "multiplicative_outputs_api.hpp"
#include "separate_cpp_files.hpp"
#include "complex_argument_interchannel_clustering_m1.hpp"
#include "Sawada_source_reconstruction_m1.hpp"
#include "Compute_deltas.hpp"

#endif
/*#ifndef LOCAL_INC_H_*/