#include "local_inc.hpp"

static void inspect_Zol_scale_down_Ylk(mat* Z_ol_p, mat* Y_lk_p){

int l_index;

double lth_l1_norm_in_Z;

double threshold=0.2*pow(10, -9);

for (l_index=0; l_index<L_static; l_index++){

lth_l1_norm_in_Z=as_scalar(sum((*Z_ol_p).col(l_index)));

	if (lth_l1_norm_in_Z<(threshold)){

	(*Y_lk_p).row(l_index)=0.00001*(*Y_lk_p).row(l_index);

	}

}

}

/*0th sum stuff*/
/*cube::fixed<F_static, N_static, L_static> outtensor_0th_real_FNL;
cx_cube::fixed<F_static, N_static, L_static> outtensor_0th_cx_FNL;*/
mat::fixed<K_static, L_static> Z_inmat_real_Ykl;

/*cube::fixed<M_static, N_static, F_static> Intensor_Xhat_low_mnf;
cx_cube::fixed<M_static, N_static, F_static> Intensor_E_conj_mnf;*/

/*cube::fixed<N_static, L_static, F_static> outtensor_0th_real_NLF;
cx_cube::fixed<N_static, L_static, F_static> outtensor_0th_cx_NLF;

static cx_mat::fixed<N_static, K_static> dummy_mat_cx_NK[NUM_WORKER_THREADS]; 
static mat::fixed<K_static, K_static> dummy_mat_outerprod_KK[NUM_WORKER_THREADS]; 

static mat::fixed<N_static, K_static> dummy_mat_real_NK[NUM_WORKER_THREADS]; */

/*1st sum stuff*/
/*cube::fixed<L_static, F_static, M_static> outtensor_1st_real_LFM;
cx_cube::fixed<L_static, F_static, M_static> outtensor_1st_cx_LFM;*/

/*2nd sum stuff*/
/*cube::fixed<O_static, L_static, M_static> outtensor_2nd_real_OLM;
cx_cube::fixed<O_static, L_static, M_static> outtensor_2nd_cx_OLM;*/

/*3rd sum stuff*/
mat::fixed<O_static, L_static> outmat_3rd_real_OL_den;
mat::fixed<O_static, L_static> outmat_3rd_real_OL_num;
/*cx_mat::fixed<O_static, L_static> outmat_3rd_cx_OL;*/

static void Z_primary_auxfun_compute_output_at_pair_indices_do_work(cx_cube* Xhat_fnm_p, cube* Xhat_low_fnm_p, mat* Y_lk_p, mat* T_fk_p, mat* V_nk_p, cube* W_fom_p, cx_cube* W_fom_cx_p, cx_cube* expj_Phi_S_nkf_p, int o_index, int l_index, int thread_iter){

}

void projfunc_wrapper_Zol(mat* Z_ol_p, mxArray *Z_ol_mxArray_p){

mxArray* output_arg_p[1];

armaSetPr(Z_ol_mxArray_p, *Z_ol_p);

/*mexCallMATLAB(1, &output_arg_p, 1, &plhs[2], "projfunc_wrapper_matlab");*/
/*mexCallMATLAB(0, NULL, 1, &Z_ol_mxArray_p, "projfunc_wrapper_matlab");*/
mexCallMATLAB(1, &output_arg_p[0], 1, &Z_ol_mxArray_p, "projfunc_wrapper_matlab_Zol");

/*memcpy(mxGetPr(Z_ol_mxArray_p), mxGetPr(output_arg_p), O_static*L_static*sizeof(double));*/

/*arrayops::copy((*Z_ol_p).memptr(), mxGetPr(output_arg_p[0]), O_static*L_static);*/

(*Z_ol_p).set_real(armaGetPr(output_arg_p[0], true)); 	

}

/*With the exception of the W_fom update, the computation of the 3rd sum involves only M components, and therefore no need to set up a work queue for the computations for this*/
/*Z_primary_auxfun_compute_outstanding_work_last_thread(argStruct_p->Z_ol_p, argStruct_p->Y_lk_p, argStruct_p);*/
static void Z_primary_auxfun_compute_outstanding_work_last_thread(mat* Z_ol_p, mat* Y_lk_p, arg_struct_t* argStruct_p){

}

static void Z_primary_auxfun_compute_output_at_pair_indices_complete_work_wrapper(arg_struct_t* argStruct_p, int thread_iter){


int o_index;
int l_index;

bool last_element_queue_flag=false;

bool compute_flag; 

while (threads_while_condition_0th_sum_flag&&(!(last_thread_to_sleep_flag_array[0][thread_iter]))){

	if ((!last_element_found_global_flag)&&(check_queue_global_flag)){

	/*Semaphore or mutex lock */
	sem_wait(&queue_sem);

	compute_flag=calculate_pair_queue_indices(argStruct_p, &o_index, &l_index, &last_element_queue_flag); 

	/*Semaphore or mutex unlock */
	sem_post(&queue_sem);

	}

	if (compute_flag){

	Z_primary_auxfun_compute_output_at_pair_indices_do_work(argStruct_p->Xhat_fnm_p, argStruct_p->Xhat_low_fnm_p, argStruct_p->Y_lk_p, argStruct_p->T_fk_p, argStruct_p->V_nk_p, argStruct_p->W_fom_p, argStruct_p->W_fom_cx_p, argStruct_p->expj_Phi_S_nkf_p, o_index, l_index, thread_iter);

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


void* Z_primary_auxfun_start(void* arg){

int i_iter;
int thread_iter;

thread_arg_t* threadArg_p;
arg_struct_t* argStruct_p;

threadArg_p=(thread_arg_t*)arg;
argStruct_p=threadArg_p->argStruct_p;
thread_iter=threadArg_p->thread_iter; 

Z_primary_auxfun_compute_output_at_pair_indices_complete_work_wrapper(argStruct_p, thread_iter);

if (last_thread_to_sleep_flag_array[0][thread_iter]){

Z_primary_auxfun_compute_outstanding_work_last_thread(argStruct_p->Z_ol_p, argStruct_p->Y_lk_p, argStruct_p);

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
