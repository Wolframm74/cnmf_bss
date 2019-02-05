#include "local_inc.hpp"

colvec::fixed<N_static> ones_col_Nx1; /*This needs to be set to .ones() just once. ever, before either Z update or W update */

/*0th sum corresponds to basic matrix multiplication*/
mat::fixed<O_static, K_static> W_fom_outmat_0th_real_OK; 

/*1st sum stuff*/
/*cube::fixed<N_static, O_static, F_static> outtensor_1st_real_NOF;
cx_cube::fixed<N_static, O_static, F_static> outtensor_1st_cx_NOF;
static cx_mat::fixed<N_static, K_static> dummy_mat_cx_NK[NUM_WORKER_THREADS];
static mat::fixed<K_static, K_static> dummy_mat_outerprod_KK[NUM_WORKER_THREADS];
static mat::fixed<N_static, K_static> dummy_mat_real_NK[NUM_WORKER_THREADS];*/

/*2nd sum stuff*/
cube::fixed<F_static, O_static, M_static> outtensor_2nd_real_FOM;
cx_cube::fixed<F_static, O_static, M_static> outtensor_2nd_cx_FOM;

/*Numerator tensor*/
cube::fixed<F_static, O_static, M_static> numerator_tensor_real_FOM;

/*cube::fixed<M_static, O_static, F_static> W_outtensor_2nd_real_MOF;
cx_cube::fixed<M_static, O_static, F_static> W_outtensor_2nd_cx_MOF;*/

void W_0th_sum_preprocess(mat* Z_ol_p, mat* Y_lk_p){

W_fom_outmat_0th_real_OK=((*Z_ol_p)*(*Y_lk_p));

}

/*
Argument List:
arrayp[0] outtensor_1st_real_NOF.
arrayp[1] outtensor_1st_cx_NOF.
arrayp[2] outtensor_2nd_real_FOM.
arrayp[3] outtensor_2nd_cx_FOM. */

#ifdef BOTTOM_UP_MEX_FLAG

void W_callMATLAB_wrapper(void){

mxArray *array_p[4];

array_p[0]=armaCreateMxMatrix(outtensor_1st_real_NOF.n_rows, outtensor_1st_real_NOF.n_cols, outtensor_1st_real_NOF.n_slices, mxDOUBLE_CLASS, mxREAL);
armaSetCubePr(array_p[0], outtensor_1st_real_NOF);

array_p[1]=armaCreateMxMatrix(outtensor_1st_cx_NOF.n_rows, outtensor_1st_cx_NOF.n_cols, outtensor_1st_cx_NOF.n_slices, mxDOUBLE_CLASS, mxCOMPLEX);
armaSetCubeCx(array_p[1], outtensor_1st_cx_NOF);

/*array_p[2]=armaCreateMxMatrix(outtensor_2nd_real_FOM.n_rows, outtensor_2nd_real_FOM.n_cols, outtensor_2nd_real_FOM.n_slices, mxDOUBLE_CLASS, mxREAL);
armaSetCubePr(array_p[2], outtensor_2nd_real_FOM);

array_p[3]=armaCreateMxMatrix(outtensor_2nd_cx_FOM.n_rows, outtensor_2nd_cx_FOM.n_cols, outtensor_2nd_cx_FOM.n_slices, mxDOUBLE_CLASS, mxCOMPLEX);
armaSetCubeCx(array_p[3], outtensor_1st_cx_NOF);*/

mexCallMATLAB(0, NULL, 2, &array_p[0], "Phi_W_unit_test");

}

#endif

/*Remember that this function should be called after first calling populate_W_fom_cx() at least once in order to have an up to date W_fom_cx */
void magnitude_square_rooted_processing_W_wrapper(arg_struct_t* argStruct_p, mxArray* plhs[]){

/*Z_ol, 			plhs[2]*/	
armaSetPr(plhs[2], (*(argStruct_p->Z_ol_p)));
engPutVariable(mlEngine_p, "Z_ol", plhs[2]);

/*T_fk, 			plhs[4]*/	
armaSetPr(plhs[4], (*(argStruct_p->T_fk_p)));
engPutVariable(mlEngine_p, "T_fk", plhs[4]);

/*W_fom_cx*/
armaSetCubeCx(plhs[8], (*(argStruct_p->W_fom_cx_p)));

/*engPutVariable(mlEngine_p, "Z_ol", plhs[2]);*/
engPutVariable(mlEngine_p, "W_fom_cx", plhs[8]);
/*engEvalString(mlEngine_p, "[W_fom, Z_ol, T_fk]=frob_norm_processing_W(W_fom_cx, Z_ol, T_fk);");*/
engEvalString(mlEngine_p, "[W_fom]=frob_norm_processing_W(W_fom_cx);");

plhs[1]=engGetVariable(mlEngine_p, "W_fom");
/*armaSetCubePr(plhs[1], (*(argStruct_p->W_fom_p)));*/

(*(argStruct_p->W_fom_p)).set_real(armaGetCubePr(plhs[1], true)); 	

populate_W_fom_cx(argStruct_p->W_fom_p, argStruct_p->W_fom_cx_p, argStruct_p->Phi_W_fom_p);

/*Also get Z_ol*/
/*plhs[2]=engGetVariable(mlEngine_p, "Z_ol");
(*(argStruct_p->Z_ol_p)).set_real(armaGetPr(plhs[2], true)); */

/*Also get T_fk*/
/*plhs[4]=engGetVariable(mlEngine_p, "T_fk");
(*(argStruct_p->T_fk_p)).set_real(armaGetPr(plhs[4], true)); */

/*engEvalString(mlEngine_p, "[V_nk]=projfunc_wrapper_matlab_Vnk_mex(V_nk, N_static, K_static);");
plhs[5]=engGetVariable(mlEngine_p, "V_nk");
(*(argStruct_p->V_nk_p)).set_real(armaGetPr(plhs[5], true)); 	

// Clip the negative elements again..
(*(argStruct_p->V_nk_p)).elem(find((*(argStruct_p->V_nk_p))<=0)).fill(0.00000001);*/

}

static mat::fixed<N_static, K_static> accu_kronmat_real_local_NK[NUM_WORKER_THREADS];

static void W_primary_auxfun_compute_output_at_pair_indices_do_work(cx_cube* Xhat_fnm_p, cube* Xhat_low_fnm_p, mat* Z_ol_p, mat* Y_lk_p, mat* T_fk_p, mat* V_nk_p, cx_cube* expj_Phi_W_fom_p, cx_cube* expj_Phi_S_nkf_p, int f_index, int o_index, int thread_iter){

int m_index;

/*populate the denominator tensor outtensor_2nd_real_FOM*/
for (m_index=0; m_index<M_static; m_index++){

accu_kronmat_real_local_NK[thread_iter]=(*V_nk_p)%kron(ones_col_Nx1%trans((*Xhat_low_fnm_p).slice(m_index).row(f_index)), ((*T_fk_p).row(f_index))%(W_fom_outmat_0th_real_OK.row(o_index)));

outtensor_2nd_real_FOM(f_index, o_index, m_index)=as_scalar(accu(accu_kronmat_real_local_NK[thread_iter]));

}

}

/*W_primary_auxfun_compute_outstanding_work_last_thread(argStruct_p->W_fom_p, argStruct_p->expj_Phi_W_fom_p);*/
static void W_primary_auxfun_compute_outstanding_work_last_thread(cube* W_fom_p, cx_cube* expj_Phi_W_fom_p){

mat mean_W=zeros<mat>(1);

/*outtensor_2nd_cx_FOM.print("outtensor_2nd_cx_FOM:");
(*expj_Phi_W_fom_p).print("expj_Phi_W_fom:");*/

/*numerator_tensor_real_FOM=outtensor_2nd_real_FOM+real(outtensor_2nd_cx_FOM%(*expj_Phi_W_fom_p));*/

outtensor_2nd_cx_FOM=Xhat_outtensor_cx_FOM;

numerator_tensor_real_FOM=outtensor_2nd_real_FOM+real(outtensor_2nd_cx_FOM);

/*numerator_tensor_real_FOM.print("numerator_tensor_real_FOM:");
outtensor_2nd_real_FOM.print("outtensor_2nd_real_FOM:");*/

/*Normalize by l1-norms*/
/*numerator_tensor_real_FOM=numerator_tensor_real_FOM/accu(numerator_tensor_real_FOM);

outtensor_2nd_real_FOM=outtensor_2nd_real_FOM/accu(outtensor_2nd_real_FOM);*/

/*POSTPONE THIS~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
/*ok. should be fine. auxiliary function 2 stuff should accumulate numerator stuff to numerator_tensor_real_FOM and denominator stuff to outtensor_2nd_real_FOM, going forward*/
/*(*W_fom_p)=(*W_fom_p)%(numerator_tensor_real_FOM/outtensor_2nd_real_FOM);*/

/*NVM DO IT NOW~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/


/*Project V_nk onto the nonnegative orthant. ie: set all negative elements to zero. */
/*(*W_fom_p).elem(find((*W_fom_p)<0)).zeros();*/
(*W_fom_p).elem(find((*W_fom_p)<0)).fill(0.00000001);

(*W_fom_p).elem(find((*W_fom_p)>100)).fill(100);

/*(*W_fom_p).elem(find((*W_fom_p)<0)).fill(0.0000000000001);*/

/*Make W_fom_p mean (1/(O*mean_Z_ol)) */
accum_W_fom=accu(*W_fom_p);

/*Subtract off the current mean*/
/**W_fom_p=*W_fom_p-(accum_W_fom/(((double)F_static)*((double)O_static)*((double)M_static)))*ones(F_static, O_static, M_static);*/

/*Set a new mean forcibly*/
/**W_fom_p=(*W_fom_p)+sqrt(1/((double)O_static))*ones(F_static, O_static, M_static);*/

mean_W=accu(*W_fom_p)/(F_static*O_static*M_static);

mean_W.print("W update: line 121: mean_W=");

}

/*real double matrix for which each f,o'th element is the frobenius norm for that Mx1 column vector*/
mat::fixed<F_static, O_static> FrobNormMat;
mat::fixed<F_static, O_static> ones_mat_FxO;		/*Set this to ones in the algorithm mex entry point*/

/*Compute the Frobenius norm across columns of 1:M for all f=1:F and o=1:O */
static void normalize_W_local(cube* W_fom_p){

int m;

FrobNormMat.zeros();

/*for 1:M; accumulate the squares*/ 
for (m=0; m<M_static; m++){

	FrobNormMat=FrobNormMat+((*W_fom_p).slice(m)%(*W_fom_p).slice(m));

}

/*sqrt(of an FxO matrix)*/
FrobNormMat=sqrt(FrobNormMat);

/*for 1:M; elementwise multiply each slice of W_fom with ones(F_static, O_static)./FrobNormMat */	
for (m=0; m<M_static; m++){

	(*W_fom_p).slice(m)=((*W_fom_p).slice(m))%(ones_mat_FxO/FrobNormMat);

}

}	


static void W_primary_auxfun_compute_output_at_pair_indices_complete_work_wrapper(arg_struct_t* argStruct_p, int thread_iter){


int f_index;
int o_index;

bool last_element_queue_flag=false;

bool compute_flag; 

while (threads_while_condition_0th_sum_flag&&(!(last_thread_to_sleep_flag_array[0][thread_iter]))){

	if ((!last_element_found_global_flag)&&(check_queue_global_flag)){

	/*Semaphore or mutex lock */
	sem_wait(&queue_sem);

	compute_flag=calculate_pair_queue_indices(argStruct_p, &f_index, &o_index, &last_element_queue_flag); 

	/*Semaphore or mutex unlock */
	sem_post(&queue_sem);

	}

	if (compute_flag){

	W_primary_auxfun_compute_output_at_pair_indices_do_work(argStruct_p->Xhat_fnm_p, argStruct_p->Xhat_low_fnm_p, argStruct_p->Z_ol_p, argStruct_p->Y_lk_p, argStruct_p->T_fk_p, argStruct_p->V_nk_p, argStruct_p->expj_Phi_W_fom_p, argStruct_p->expj_Phi_S_nkf_p, f_index, o_index, thread_iter);

	sem_wait(&work_actually_computed_sem);
	work_actually_computed_ctr++;

	if (work_actually_computed_ctr==argStruct_p->work_queue_size){

		/*Wake up the last element thread*/
		/*sem_post(&last_element_sleep_sem);*/
		work_actually_completed_flag=true; 

		/*FLAG THIS thread_iter'th THREAD AS THE LAST THREAD TO SLEEP. IE: IT DOESN'T SLEEP. IT GETS TO BREAK THE LOOP ON ITS WAY OUT.*/
		last_thread_to_sleep_flag_array[0][thread_iter]=true;

		/*Set this flag on so that threads no longer go into the check queue block. May not be necessary if last_element_found_global is already set*/
		sem_wait(&check_queue_global_sem);
		check_queue_global_flag=true;
		sem_post(&check_queue_global_sem);

	}

	sem_post(&work_actually_computed_sem);

	compute_flag=false;

	}

	/*Code block that follows is a trap. Last thread to sleep should be able to bypass this trap. 
	Other threads than this thread should only enter the trap once the work_actually_completed_flag has been set to true. 
	*/	

	if (work_actually_completed_flag&&(!(last_thread_to_sleep_flag_array[0][thread_iter]))){

	sem_wait(&threads_asleep_sem);
	threads_asleep_ctr++;
	sem_post(&threads_asleep_sem);

	sem_wait(&sleep_sem);

	sem_wait(&threads_asleep_sem);
	threads_asleep_ctr--;	
	sem_post(&threads_asleep_sem);

	}

}	

}

void* W_primary_auxfun_start(void* arg){

int i_iter;
int thread_iter;

thread_arg_t* threadArg_p;
arg_struct_t* argStruct_p;

threadArg_p=(thread_arg_t*)arg;
argStruct_p=threadArg_p->argStruct_p;
thread_iter=threadArg_p->thread_iter; 

W_primary_auxfun_compute_output_at_pair_indices_complete_work_wrapper(argStruct_p, thread_iter);

if (last_thread_to_sleep_flag_array[0][thread_iter]){

W_primary_auxfun_compute_outstanding_work_last_thread(argStruct_p->W_fom_p, argStruct_p->expj_Phi_W_fom_p);

sem_wait(&threads_while_condition_0th_sum_sem);
threads_while_condition_0th_sum_flag=false;	
sem_post(&threads_while_condition_0th_sum_sem);

/*Sem post to wake NUM_WORKER_THREADS-1 of the sleeping threads*/

for (i_iter=0; i_iter<(NUM_WORKER_THREADS-1); i_iter++){
	sem_post(&sleep_sem);
}

/*Sem wait NUM_WORKER_THREADS-1 times so that you only proceed from this point once that many of the sleeping flags have high fived and successfully made it out of the "while traps" */

for (i_iter=0; i_iter<(NUM_WORKER_THREADS-1); i_iter++){
	sem_wait(&checkpoint_sem);
}

}

if (!last_thread_to_sleep_flag_array[0][thread_iter]){

sem_post(&checkpoint_sem);

}

threads_exit_and_signal(NUM_WORKER_THREADS);
     
}

