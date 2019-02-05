#include "local_inc.hpp"

/*
int complex_argument_costfun_m4_m_index_b;
int complex_argument_costfun_m4_m_index_a;
*/

static rowvec::fixed<K_static> W_channel_b_complex_argument_m4_num_spreadrow_1xK[NUM_WORKER_THREADS_COMPLEX_ARGUMENT_COSTFUN_M4];
static rowvec::fixed<K_static> W_channel_b_complex_argument_m4_den_spreadrow_1xK[NUM_WORKER_THREADS_COMPLEX_ARGUMENT_COSTFUN_M4];

mat::fixed<F_static, O_static> W_channel_b_complex_argument_m4_num_outmat_FO;
mat::fixed<F_static, O_static> W_channel_b_complex_argument_m4_den_outmat_FO;

static void W_channel_b_complex_argument_m4_compute_output_at_pair_indices_do_work(cx_cube* Xhat_fnm_p, cube* Xhat_low_fnm_p, mat* Z_ol_p, mat* Y_lk_p, mat* T_fk_p, mat* V_nk_p, cx_cube* expj_Phi_W_fom_p, cx_cube* expj_Phi_S_nkf_p, int f_index, int o_index, int thread_iter){

int l_index, n_index;


/*num stuff: Formula: real(( 2nd term )%( 3rd term )) */	

// 2nd term


// 3rd term
/*conj(expj_Phi_U_flnm[complex_argument_costfun_m4_m_index_b](f_index, l_index, n_index))*( ((complex_argument_costfun_m4_Error_mat_FN(f_index, n_index))*((*Xhat_fnm_p)(f_index, n_index, complex_argument_costfun_m4_m_index_b))) )*conj((*expj_Phi_W_fom_p)(f_index, o_index, complex_argument_costfun_m4_m_index_b))*((*Z_ol_p)(o_index, l_index))*((*V_nk_p).row(n_index)%conj((*expj_Phi_S_nkf_p).slice(f_index).row(n_index)))%( (*Y_lk_p).row(l_index)%(*T_fk_p).row(f_index) )*/


/*den stuff: Formula: ( 2nd term )%( 3rd term ) */	

// 2nd term	

// 3rd term	
/*(((complex_argument_costfun_m4_Magnitude_model_mat_FN(f_index, n_index))*((*Xhat_low_fnm_p)(f_index, n_index, complex_argument_costfun_m4_m_index_b)))*((*Z_ol_p)(o_index, l_index))*((*V_nk_p).row(n_index))%((*Y_lk_p).row(l_index))%((*T_fk_p).row(f_index)))*/


W_channel_b_complex_argument_m4_num_spreadrow_1xK[thread_iter].zeros();
W_channel_b_complex_argument_m4_den_spreadrow_1xK[thread_iter].zeros();

for (l_index=0; l_index<L_static; l_index++){

for (n_index=0; n_index<N_static; n_index++){	

W_channel_b_complex_argument_m4_num_spreadrow_1xK[thread_iter]=W_channel_b_complex_argument_m4_num_spreadrow_1xK[thread_iter]+real(conj(expj_Phi_U_flnm[complex_argument_costfun_m4_m_index_b](f_index, l_index, n_index))*( ((complex_argument_costfun_m4_Error_mat_FN(f_index, n_index))*((*Xhat_fnm_p)(f_index, n_index, complex_argument_costfun_m4_m_index_b))) )*conj((*expj_Phi_W_fom_p)(f_index, o_index, complex_argument_costfun_m4_m_index_b))*((*Z_ol_p)(o_index, l_index))*((*V_nk_p).row(n_index)%conj((*expj_Phi_S_nkf_p).slice(f_index).row(n_index)))%( (*Y_lk_p).row(l_index)%(*T_fk_p).row(f_index) ));
W_channel_b_complex_argument_m4_den_spreadrow_1xK[thread_iter]=W_channel_b_complex_argument_m4_den_spreadrow_1xK[thread_iter]+(((complex_argument_costfun_m4_Magnitude_model_mat_FN(f_index, n_index))*((*Xhat_low_fnm_p)(f_index, n_index, complex_argument_costfun_m4_m_index_b)))*((*Z_ol_p)(o_index, l_index))*((*V_nk_p).row(n_index))%((*Y_lk_p).row(l_index))%((*T_fk_p).row(f_index)));

}

}

W_channel_b_complex_argument_m4_num_outmat_FO(f_index, o_index)=2*accu(W_channel_b_complex_argument_m4_num_spreadrow_1xK[thread_iter]);
W_channel_b_complex_argument_m4_den_outmat_FO(f_index, o_index)=2*accu(W_channel_b_complex_argument_m4_den_spreadrow_1xK[thread_iter]);

W_channel_b_complex_argument_m4_num_outmat_FO(f_index, o_index)=W_channel_b_complex_argument_m4_num_outmat_FO(f_index, o_index)+W_channel_b_complex_argument_m4_den_outmat_FO(f_index, o_index);

}

static void W_channel_b_complex_argument_m4_compute_outstanding_work_last_thread(cube* W_fom_p, arg_struct_t* argStruct_p){

/*(*W_fom_p).slice(complex_argument_costfun_m4_m_index_b)=(*W_fom_p).slice(complex_argument_costfun_m4_m_index_b)%(W_channel_b_complex_argument_m4_num_outmat_FO/W_channel_b_complex_argument_m4_den_outmat_FO);*/

numerator_tensor_real_FOM.slice(complex_argument_costfun_m4_m_index_b)=numerator_tensor_real_FOM.slice(complex_argument_costfun_m4_m_index_b)+W_channel_b_complex_argument_m4_num_outmat_FO;
outtensor_2nd_real_FOM.slice(complex_argument_costfun_m4_m_index_b)=outtensor_2nd_real_FOM.slice(complex_argument_costfun_m4_m_index_b)+W_channel_b_complex_argument_m4_den_outmat_FO;

}

static void W_channel_b_complex_argument_m4_compute_output_at_pair_indices_complete_work_wrapper(arg_struct_t* argStruct_p, int thread_iter){

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

	W_channel_b_complex_argument_m4_compute_output_at_pair_indices_do_work(argStruct_p->Xhat_fnm_p, argStruct_p->Xhat_low_fnm_p, argStruct_p->Z_ol_p, argStruct_p->Y_lk_p, argStruct_p->T_fk_p, argStruct_p->V_nk_p, argStruct_p->expj_Phi_W_fom_p, argStruct_p->expj_Phi_S_nkf_p, f_index, o_index, thread_iter);

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

/*Split this mondule into two modules. One for the a channel, One for the b channel. */
void* W_channel_b_complex_argument_m4_start(void* arg){

int i_iter;
int thread_iter;

thread_arg_t* threadArg_p;
arg_struct_t* argStruct_p;

threadArg_p=(thread_arg_t*)arg;
argStruct_p=threadArg_p->argStruct_p;
thread_iter=threadArg_p->thread_iter; 

W_channel_b_complex_argument_m4_compute_output_at_pair_indices_complete_work_wrapper(argStruct_p, thread_iter);

if (last_thread_to_sleep_flag_array[0][thread_iter]){

W_channel_b_complex_argument_m4_compute_outstanding_work_last_thread(argStruct_p->W_fom_p, argStruct_p);

sem_wait(&threads_while_condition_0th_sum_sem);
threads_while_condition_0th_sum_flag=false;	
sem_post(&threads_while_condition_0th_sum_sem);

/*Sem post to wake NUM_WORKER_THREADS_COMPLEX_ARGUMENT_COSTFUN_M4-1 of the sleeping threads*/

for (i_iter=0; i_iter<(NUM_WORKER_THREADS_COMPLEX_ARGUMENT_COSTFUN_M4-1); i_iter++){
	sem_post(&sleep_sem);
}

/*Sem wait NUM_WORKER_THREADS_COMPLEX_ARGUMENT_COSTFUN_M4-1 times so that you only proceed from this point once that many of the sleeping flags have high fived and successfully made it out of the "while traps" */

for (i_iter=0; i_iter<(NUM_WORKER_THREADS_COMPLEX_ARGUMENT_COSTFUN_M4-1); i_iter++){
	sem_wait(&checkpoint_sem);
}

}

if (!last_thread_to_sleep_flag_array[0][thread_iter]){

sem_post(&checkpoint_sem);

}

threads_exit_and_signal(NUM_WORKER_THREADS_COMPLEX_ARGUMENT_COSTFUN_M4);
    
}