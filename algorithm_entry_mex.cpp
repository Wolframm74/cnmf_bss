#include "local_inc.hpp"

#ifdef ENABLE_MATLAB_ENGINE_PLOTTING
Engine* mlEngine_p;
#endif

bool N_ITERATIONS_TV_bool;

/*Global cx_mat's or in some cases cubes */
/*Take the convention to set the columns as whatever dimension has the longest length*/

/*3-dim tensors*/
cx_cube::fixed<F_static, N_static, M_static> Xtilde_fnm;

cx_cube::fixed<F_static, N_static, M_static> Xhat_fnm;	/*work queue: f=1:F*/
cube::fixed<F_static, N_static, M_static> Xhat_low_fnm;
cx_cube::fixed<F_static, N_static, M_static> E_conj_fnm;

cube::fixed<M_static, N_static, F_static> Intensor_Xhat_low_mnf;	/*Like in Phi_S_update, a re-mapping should take place of elements from Xhat_low_fnm to Xhat_low_mnf */
cx_cube::fixed<M_static, N_static, F_static> Intensor_E_conj_mnf;

/*cx_cube::fixed<F_static, M_static, N_static> Xhat_fmn;	
cube::fixed<F_static, M_static, N_static> Xhat_low_fmn;
cx_cube::fixed<F_static, M_static, N_static> E_conj_fmn;*/

/*cube::fixed<F_static, L_static, M_static> H_flm; 
cube::fixed<F_static, L_static, M_static> Phi_H_flm; */

cube::fixed<F_static, O_static, M_static> W_fom;		/*work queue: f=1:F*/
cx_cube::fixed<F_static, O_static, M_static> W_fom_cx;	/*work queue: f=1:F*/
mat::fixed<O_static, L_static> Z_ol;					/*work queue: l=1:L*/
mat::fixed<L_static, K_static> Y_lk;					/*work queue: k=1:K*/

mat::fixed<F_static, K_static> T_fk; 				/*work queue: k=1:K*/
mat::fixed<N_static, K_static> V_nk; 				/*work queue: k=1:K*/

cube::fixed<F_static, O_static, M_static> Phi_W_fom; 							/*work queue: f=1:F*/
/*cube::fixed<F_static, N_static, K_static> Phi_S_fnk;*/ 							/*work queue: k=1:K*/
cube::fixed<N_static, K_static, F_static> Phi_S_nkf;
cube::fixed<F_static, K_static, N_static> Phi_S_fkn;

/* Set these to exp(Phi_W_fom) and exp(Phi_S_fnk) respectively, upon entry of parse_arguments, as well as after updating Phi_W_fom and Phi_S_fnk */
cx_cube::fixed<F_static, O_static, M_static> Phase_W_fom; 			/*Not sure if either of these are even useful. */
/*cx_cube::fixed<F_static, N_static, K_static> expj_Phi_S_fnk;*/ 					
cx_cube::fixed<N_static, K_static, F_static> expj_Phi_S_nkf; 					
cx_cube::fixed<F_static, K_static, N_static> expj_Phi_S_fkn; 

cx_cube::fixed<F_static, O_static, M_static> expj_Phi_W_fom;					

cx_cube::fixed<F_static, N_static, L_static> S_fnl; 

/*Semaphores, followed by the stuff that they protect */
sem_t b_exit_cond_sem; 	/*Can probably comment this out?*/

/*Global variables*/
bool b_exit_cond_global;

/*Global variables*/

/*Semaphores*/
sem_t check_queue_global_sem; 	/*Double check you init'd this*/
sem_t threads_asleep_sem; 		/*Double check you init'd this*/
sem_t queue_sem;	/*the queue itself is a member of the arg_struct_t struct */
sem_t sleep_sem;
sem_t threads_while_condition_0th_sum_sem;
sem_t threads_while_condition_1st_sum_sem;
sem_t threads_while_condition_2nd_sum_sem; 
sem_t last_element_sleep_sem;
sem_t work_actually_computed_sem;
sem_t checkpoint_sem;
sem_t wait_till_queue_populated_sem;
/*extern sem_t threads_while_condition_3rd_sum_sem;*/

/*Global variables*/
bool threads_while_condition_0th_sum_flag;
bool threads_while_condition_1st_sum_flag;
bool threads_while_condition_2nd_sum_flag; 
/*extern bool threads_while_condition_3rd_sum_flag;*/
bool last_element_found_global_flag; 		/*No corresponding semaphore for this obj*/
bool check_queue_global_flag; 
bool last_element_sleep_flag;
bool work_actually_completed_flag; 
int threads_asleep_ctr; 
int work_actually_computed_ctr;
/*bool Xhat_pop_MKF_for_Phi_S_flag;*/
bool Xhat_pop_NKFM_for_T_or_Phi_S_flag;
bool Xhat_pop_FKNM_for_V_flag;
bool Xhat_pop_W_cx_outtensor_flag;

bool last_thread_to_sleep_flag_array[MAX_NUM_SUMS][NUM_WORKER_THREADS];

arg_struct_t argStruct; 

#ifdef DEBUG1
bool Xhat_after_V_flag;
#endif 

#ifdef DEBUG
sem_t print_sem;
#endif 



/*static void parse_4th_order_tensor(arg_struct_t* argStruct_p, const mxArray *prhs[]){

	int l_iter;
	double* head_p=mxGetPr(prhs[11]);
	int offset;

	double* front_p;

	//Need to loop through 1:L, to populate Phi_S_fnkl 

	for (l_iter=0; l_iter<L_static; l_iter++){

		offset=l_iter*(F_static*N_static*K_static);

		front_p=(head_p+offset);

		(argStruct_p->Phi_S_fnkl_arr[l_iter])=Cube<double>(front_p, F_static, N_static, K_static, false, true);

		//(argStruct_p->Phi_S_fnkl_arr[l_iter]).print();

	}

}*/


static void parse_arguments(arg_struct_t* argStruct_p, const mxArray *prhs[]){

	mxArray* output_arg_p;	/*My application doesn't allocate any memory related to the output argument. */
							/*MATLAB will allocate the memory and upon returning from mexCallMATLAB, I think this pointer will be pointing to it */

	mxArray* input_arg_p;	/*My application DOES need to allocate memory for the input argument. */

	input_arg_p=armaCreateMxMatrix(2,1, mxDOUBLE_CLASS, mxREAL);
	mat some_mat=zeros<mat>(2,1);

	/*Iterators*/
	int l_iter, k_iter; 

	/*Argument list: 	prhs[0] 	Xtilde_fnm, 
						prhs[1] 	Xhat_fnm, 
						prhs[2] 	W_fom, 
						prhs[3]		Z_ol, 
						prhs[4]		Y_lk, 
						prhs[5] 	T_fk, 
						prhs[6] 	V_nk, 
						prhs[7]		Phi_W_fom,
						prhs[8]		Phi_S_fnk,
						prhs[9] 	M, 
						prhs[10] 	F, 
						prhs[11] 	N, 
						prhs[12] 	K, 
						prhs[13] 	L, 				
						prhs[14]	O, 
						prhs[15]	Phi_S_fkn,
						prhs[16]	expj_Phi_S_nkf,
						prhs[17]	expj_Phi_S_fkn,

						*/

	Xtilde_fnm.set_real(armaGetCubePr(prhs[0], true)); 	
	Xtilde_fnm.set_imag(armaGetCubePi(prhs[0], true)); 	
	argStruct_p->Xtilde_fnm_p=&Xtilde_fnm;

	/*Xhat_fnm = armaGetCx(prhs[1], true); 	*/
	/*may not need to read this from MATLAB*/
	argStruct_p->Xhat_fnm_p=&Xhat_fnm;		/*good.*/

	W_fom.set_real(armaGetCubePr(prhs[2], true)); 	
	argStruct_p->W_fom_p=&W_fom;

	/*Project V_nk onto the nonnegative orthant. ie: set all negative elements to zero. */
	(W_fom).elem(find((W_fom)<0)).zeros();

	Z_ol.set_real(armaGetPr(prhs[3], true)); 	
	argStruct_p->Z_ol_p=&Z_ol; 

	/*Project Z_ol onto the nonnegative orthant. ie: zero all negative elements*/
	(Z_ol).elem(find((Z_ol)<0)).zeros();	

	Y_lk.set_real(armaGetPr(prhs[4], true)); 	
	argStruct_p->Y_lk_p=&Y_lk; 

	T_fk.set_real(armaGetPr(prhs[5], true)); 	
	argStruct_p->T_fk_p=&T_fk; 

	/*Set negative elements to zero*/
	(T_fk).elem(find((T_fk)<0)).zeros();

	V_nk.set_real(armaGetPr(prhs[6], true)); 	
	argStruct_p->V_nk_p=&V_nk; 

	Phi_W_fom.set_real(armaGetCubePr(prhs[7],true));
	argStruct_p->Phi_W_fom_p=&Phi_W_fom; 		

	Phi_S_nkf.set_real(armaGetCubePr(prhs[8],true));
	argStruct_p->Phi_S_nkf_p=&Phi_S_nkf; 		

	argStruct_p->M = armaGetScalar<int>(prhs[9]); 	
	argStruct_p->F = armaGetScalar<int>(prhs[10]); 	
	argStruct_p->N = armaGetScalar<int>(prhs[11]); 	
	argStruct_p->K = armaGetScalar<int>(prhs[12]); 	
	argStruct_p->L = armaGetScalar<int>(prhs[13]); 	
	argStruct_p->O = armaGetScalar<int>(prhs[14]); 	

	Phi_S_fkn.set_real(armaGetCubePr(prhs[15],true));
	argStruct_p->Phi_S_fkn_p=&Phi_S_fkn; 		

	expj_Phi_S_nkf.set_real(armaGetCubePr(prhs[16],true));
	expj_Phi_S_nkf.set_imag(armaGetCubePi(prhs[16],true));	
	argStruct_p->expj_Phi_S_nkf_p=&expj_Phi_S_nkf; 		

	expj_Phi_S_fkn.set_real(armaGetCubePr(prhs[17],true));
	expj_Phi_S_fkn.set_imag(armaGetCubePi(prhs[17],true));	
	argStruct_p->expj_Phi_S_fkn_p=&expj_Phi_S_fkn; 		

	expj_Phi_W_fom.set_real(armaGetCubePr(prhs[18],true));
	expj_Phi_W_fom.set_imag(armaGetCubePi(prhs[18],true));
	argStruct_p->expj_Phi_W_fom_p=&expj_Phi_W_fom;

	/*Test out the parsing of the 4th order tensor*/
	/*parse_4th_order_tensor(argStruct_p, prhs);*/

	/*Lastly, set the following pointer to its corresponding object */
	/*argStruct_p->S_fnl_p=&S_fnl; */

	Xhat_low_fnm.zeros();
	argStruct_p->Xhat_low_fnm_p=&Xhat_low_fnm;	/*good*/

	E_conj_fnm.zeros();
	argStruct_p->E_conj_fnm_p=&E_conj_fnm;		/*good*/

	argStruct_p->W_fom_cx_p=&W_fom_cx;

	populate_W_fom_cx(argStruct_p->W_fom_p, argStruct_p->W_fom_cx_p, argStruct_p->Phi_W_fom_p);

}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){

double X_L1_before;
double X_L1_after;
double T_L1_before;
double T_L1_after;
double V_L1_before;
double V_L1_after;
double W_L1_before;
double W_L1_after;
double PS_L1_before;
double PS_L1_after;

double X_L2_before;
double X_L2_after;
double T_L2_before;
double T_L2_after;
double V_L2_before;
double V_L2_after;
double W_L2_before;
double W_L2_after;
double PS_L2_before;
double PS_L2_after;

N_ITERATIONS_TV_bool=true;

colvec::fixed<1> local_scalar;

/*double sparsity_fill_value=0.7;

populate_sigma_vec(sparsity_fill_value);*/

/*Call this function*/
populate_window_vec();

populate_6x1_possible_permutations();

/*Set this to ones*/
ones_col_Nx1.ones();
ones_col_Ox1.ones();
ones_col_Lx1.ones();
ones_mat_FxO.ones();
/*ones_mat_NxK.ones();*/
ones_row_1xK.ones();
ones_col_Fx1.ones();
/*ones_tensor_real_FKN.ones();
ones_tensor_real_NKF.ones();*/

#ifdef ENABLE_MATLAB_ENGINE_PLOTTING

local_scalar.print("line 258\n");

	if (!(mlEngine_p = engOpen(""))) {
		mexErrMsgTxt("Can't start MATLAB engine\n");
		/*return EXIT_FAILURE;*/
	}

/*maybe addpath() to where the plot functions reside*/
engEvalString(mlEngine_p, "working_directory=pwd();");	

engEvalString(mlEngine_p, "variable1_=zeros(5,1);");	

engEvalString(mlEngine_p, "figure, plot(variable1_);");	
engEvalString(mlEngine_p, "title(working_directory);");	

engEvalString(mlEngine_p, "addpath('../../../../libs/getComputerName');");	

engEvalString(mlEngine_p, "mycomputer=getComputerName();");
engEvalString(mlEngine_p, "if (strcmp(mycomputer, 'tung-gc670aar-aba-a6120n')); codepath_prefix='/mnt/asus_share/'; elseif (strcmp(mycomputer, 'tn-k53e')); codepath_prefix='/home/tung/Documents/'; elseif (strcmp(mycomputer, 'tung-all-series-4770k')); codepath_prefix='/mnt/asus_share/'; end");

/*engEvalString(mlEngine_p, "codepath_prefix='/mnt/asus_share/';");
engEvalString(mlEngine_p, "elseif (strcmp(mycomputer, 'tn-k53e'))");   
engEvalString(mlEngine_p, "codepath_prefix='/home/tung/Documents/';");
engEvalString(mlEngine_p, "end");*/

engEvalString(mlEngine_p, "addpath([codepath_prefix, 'School/thesis2/mycode/nguyen_2016/code_vectorised_60_3_orthogonal_speech/top_down/mex']);");
engEvalString(mlEngine_p, "addpath([codepath_prefix, 'School/thesis2/mycode/nguyen_2016/code_vectorised_60_3_orthogonal_speech/top_down/mex/separate_cpp_files']);");
engEvalString(mlEngine_p, "addpath([codepath_prefix, 'School/thesis2/mycode/codegen_common/bin']);");

/*addpaths:'matlab_common'*/
engEvalString(mlEngine_p, "addpath([codepath_prefix, 'School/thesis2/mycode/matlab_common/frob_norm_processing_W_module']);");
engEvalString(mlEngine_p, "addpath([codepath_prefix, 'School/thesis2/mycode/matlab_common/plot_significant_clusters_module']);");
engEvalString(mlEngine_p, "addpath([codepath_prefix, 'School/thesis2/mycode/nguyen_2016/code_vectorised_60_3_orthogonal_speech/top_down/mex/helper_modules/TDOA_update_module1']);");

engEvalString(mlEngine_p, "addpath([codepath_prefix, 'School/thesis2/mycode/matlab_common/plot_P_ln_and_C_l']);");
engEvalString(mlEngine_p, "addpath([codepath_prefix, 'School/thesis2/mycode/matlab_common/plot_update_rule_deltas']);");

local_scalar.print("line 279\n");

#endif	

/*BUGFIX: this was set to set_imag in code_vectorised_12*/
/*phase_shift_tensor_cx_FKN.set_real((-1)*ones_tensor_real_FKN);
phase_shift_tensor_cx_NKF.set_real((-1)*ones_tensor_real_NKF);

reference_phase_tensor_cx_FKN.set_imag(ones_tensor_real_FKN);
reference_phase_tensor_cx_NKF.set_imag(ones_tensor_real_NKF);*/

/*Set W_fom_local_p local to normalize_funs.cpp*/
W_fom_local_p=&W_fom;

/*Create array of pointer(s) to mxArrays*/
mxArray *array_p[4];

int n;
double L_value;

arg_struct_t* argStruct_p=&argStruct; 

/*open file pointer for Xhat_update*/
argStruct_p->Xhat_pFile=fopen("Xhat_update_log_file.txt", "w");
Xhat_pFile_helper=argStruct_p->Xhat_pFile;

/*Init flag*/
/*Xhat_pop_MKF_for_Phi_S_flag=false;*/
Xhat_pop_NKFM_for_T_or_Phi_S_flag=false;

/*mexCallMATLAB stuff*/

parse_arguments(argStruct_p, prhs);

/*Do abs and mask stuff here*/
/*X_tilde_fnm_abs_cube=abs(Xtilde_fnm);
X_tilde_fnm_binary_mask_cube.ones();*/
/*(X_tilde_fnm_binary_mask_cube).elem(find(X_tilde_fnm_binary_mask_cube<=0)).fill(0.00000001);;*/
/*(X_tilde_fnm_binary_mask_cube).elem(find(X_tilde_fnm_binary_mask_cube<=0.05)).fill(0);;*/

mexPrintf("HELLOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO");


/*Plot T, V, H and Phi_H beforehand as a spot check */

#ifdef DEBUG34534

armaSetPr(plhs[4], T_fk);
mexCallMATLAB(0, NULL, 1, &plhs[4], "plot_T_fk");

armaSetPr(plhs[5], V_nk);
mexCallMATLAB(0, NULL, 1, &plhs[5], "plot_V_nk");

mexPrintf("line 302");

armaSetCubePr(plhs[1], W_fom);
mexCallMATLAB(0, NULL, 1, &plhs[1], "plot_W_fom");

mexPrintf("line 307");

/*armaSetCubePr(plhs[4], Phi_W_fom);
mexCallMATLAB(0, NULL, 1, &plhs[4], "plot_Phi_H_flm");*/

#endif

create_lhs_matrices(plhs, argStruct_p);

argStruct_p->Z_ol_mxArray_p=plhs[2];
argStruct_p->Y_lk_mxArray_p=plhs[3];
argStruct_p->V_nk_mxArray_p=plhs[5];
argStruct_p->T_fk_mxArray_p=plhs[4];

/*Start While*/
/*for (n=1; n<(N_ITERATIONS_TV+1); n++){

plot_and_update_Xhat(plhs, argStruct_p, n);

plot_and_update_Tfk(plhs, prhs, argStruct_p, n);

plot_and_update_Xhat(plhs, argStruct_p, n);

plot_and_update_Vnk(plhs, argStruct_p, n);

}*/

N_ITERATIONS_TV_bool=false;

complex_argument_interchannel_clustering_m1_entry(plhs, argStruct_p, 1, 0);
complex_argument_interchannel_clustering_m1_lookdirs_search_entry(plhs, argStruct_p, 1, 0);

Xhat_pop_FKNM_for_V_flag=false;

compute_deltas_init_rowvecs();


/*Start While*/
/*for (n=1; n<(2+1); n++){*/
for (n=1; n<(N_ITERATIONS+1); n++){

/*Xhat_pop_NKFM_for_T_or_Phi_S_flag=true;*/	//Set flag
plot_and_update_Xhat(plhs, argStruct_p, n);
/*Xhat_pop_NKFM_for_T_or_Phi_S_flag=false;*/	//Unset flag

NMF_auxfun1_rowvec(n-1)=L1_global;
NMF_auxfun2_rowvec(n-1)=L2_global;
NMF_rowvec_combined(n-1)=L1_global+L2_global;

compute_deltas_NMF_auxfun_rowvecs_send_data_to_Matlab_Eng_and_plot(plhs);

X_L1_before=L1_global;
X_L2_before=L2_global;

T_L1_before=L1_global;
T_L2_before=L2_global;

plot_and_update_Tfk(plhs, prhs, argStruct_p, n);

Xhat_pop_FKNM_for_V_flag=true;	//Set flag
plot_and_update_Xhat(plhs, argStruct_p, n);
Xhat_pop_FKNM_for_V_flag=false;	//Unset flag

T_L1_after=L1_global;
T_L2_after=L2_global;

update_T_update_and_plot_rowvecs(plhs, n-1, T_L1_after-T_L1_before, T_L2_after-T_L2_before, (T_L1_after+T_L2_after)-(T_L1_before+T_L2_before));

V_L1_before=L1_global;
V_L2_before=L2_global;

plot_and_update_Vnk(plhs, argStruct_p, n);

/*Xhat_pop_NKFM_for_T_or_Phi_S_flag=true;*/	//Set flag
plot_and_update_Xhat(plhs, argStruct_p, n);
/*Xhat_pop_NKFM_for_T_or_Phi_S_flag=false;*/	//Unset flag

V_L1_after=L1_global;
V_L2_after=L2_global;

update_V_update_and_plot_rowvecs(plhs, n-1, V_L1_after-V_L1_before, V_L2_after-V_L2_before, (V_L1_after+V_L2_after)-(V_L1_before+V_L2_before));

PS_L1_before=L1_global;
PS_L2_before=L2_global;

plot_and_update_Phi_S(plhs, argStruct_p, n);

Xhat_pop_W_cx_outtensor_flag =true;
plot_and_update_Xhat(plhs, argStruct_p, n);
Xhat_pop_W_cx_outtensor_flag=false; 

PS_L1_after=L1_global;
PS_L2_after=L2_global;

update_PS_update_and_plot_rowvecs(plhs, n-1, PS_L1_after-PS_L1_before, PS_L2_after-PS_L2_before, (PS_L1_after+PS_L2_after)-(PS_L1_before+PS_L2_before));

W_L1_before=L1_global;
W_L2_before=L2_global;
 
plot_and_update_Wfom(plhs, argStruct_p, n);

plot_and_update_Xhat(plhs, argStruct_p, n);

W_L1_after=L1_global;
W_L2_after=L2_global;

update_W_update_and_plot_rowvecs(plhs, n-1, W_L1_after-W_L1_before, W_L2_after-W_L2_before, (W_L1_after+W_L2_after)-(W_L1_before+W_L2_before));

// Populate this for Z update
Z_inmat_real_Ykl=trans(*(argStruct_p->Y_lk_p));

plot_and_update_Zol(plhs, prhs, argStruct_p, n);

/*plot_and_update_Xhat(plhs, argStruct_p, n);

plot_and_update_Ylk(plhs, argStruct_p, n);*/

X_L1_after=L1_global;
X_L2_after=L2_global;
update_X_update_and_plot_rowvecs(plhs, n-1, X_L1_after-X_L1_before, X_L2_after-X_L2_before, (X_L1_after+X_L2_after)-(X_L1_before+X_L2_before));

Plot_significant_L_wrapper(plhs, prhs, argStruct_p);
/*Plot_significant_L_wrapper_3(plhs, prhs, argStruct_p);*/
/*Sawada_source_reconstruction_m1_entry(plhs, argStruct_p);*/

}	/*End While*/	

/* Finally, compute S_fnl */
/*update_engine(argStruct_p, 6);

#ifdef 	MEX_UNIT_TESTING
armaSetCubeCx(plhs[5], S_fnl);
mexCallMATLAB(0, NULL, 1, &plhs[5], "plot_S_fnl");
#endif*/


armaSetCubeCx(plhs[0], Xhat_fnm);
armaSetCubePr(plhs[1], W_fom); 
armaSetPr(plhs[2], Z_ol); 
armaSetPr(plhs[3], Y_lk); 
armaSetPr(plhs[4], T_fk); 
armaSetPr(plhs[5], V_nk);
armaSetCubeCx(plhs[6], expj_Phi_S_nkf);
armaSetCubeCx(plhs[7], expj_Phi_S_fkn);

/*close file pointer for Xhat_update*/
fclose(argStruct_p->Xhat_pFile);

#ifdef ENABLE_MATLAB_ENGINE_PLOTTING
engClose(mlEngine_p);
#endif

return;

}

