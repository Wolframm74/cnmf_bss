#include "local_inc.hpp"

/*
int complex_argument_costfun_m4_m_index_b;
int complex_argument_costfun_m4_m_index_a;
*/

/*static mat::fixed<N_static, O_static> T_complex_argument_m4_num_spreadmat_NO[NUM_WORKER_THREADS_COMPLEX_ARGUMENT_COSTFUN_M4];*/
static mat::fixed<N_static, O_static> T_complex_argument_m4_den_spreadmat_NO[NUM_WORKER_THREADS_COMPLEX_ARGUMENT_COSTFUN_M4];

mat::fixed<F_static, K_static> T_complex_argument_m4_num_outmat_FK;
mat::fixed<F_static, K_static> T_complex_argument_m4_den_outmat_FK;

static void T_complex_argument_m4_compute_output_at_pair_indices_do_work(cx_cube* Xhat_fnm_p, cube* Xhat_low_fnm_p, mat* Z_ol_p, mat* Y_lk_p, mat* V_nk_p, cube* W_fom_p, cx_cube* W_fom_cx_p, cx_cube* expj_Phi_S_nkf_p, int f_index, int k_index, int thread_iter){

int l_index;

/*num stuff: Formula: real(( 1st term + 2nd term )%( 3rd term )) */	

/*den stuff: Formula: ( 1st term + 2nd term )%( 3rd term ) */	

/*
// 1st term
kron( trans((*Xhat_low_fnm_p).slice(complex_argument_costfun_m4_m_index_a).row(f_index)) , (*W_fom_p).slice(complex_argument_costfun_m4_m_index_b).row(f_index) )

// 2nd term	
kron( trans((*Xhat_low_fnm_p).slice(complex_argument_costfun_m4_m_index_b).row(f_index)) , (*W_fom_p).slice(complex_argument_costfun_m4_m_index_a).row(f_index) )

// 3rd term	
(*Y_lk_p)(l_index, k_index)*kron( trans(complex_argument_costfun_m4_Magnitude_model_mat_FN.row(f_index))%trans((*V_nk_p).row(k_index)) , trans((*Z_ol_p).col(l_index)) )*/

/*T_complex_argument_m4_num_spreadmat_NO[thread_iter].zeros();*/
T_complex_argument_m4_den_spreadmat_NO[thread_iter].zeros();

for (l_index=0; l_index<L_static; l_index++){

/*T_complex_argument_m4_num_spreadmat_NO[thread_iter]=T_complex_argument_m4_num_spreadmat_NO[thread_iter]+real(( ( kron( strans( (*Xhat_fnm_p).slice(complex_argument_costfun_m4_m_index_a).row(f_index) )%conj( (*expj_Phi_S_nkf_p).slice(f_index).col(k_index)) , conj((*W_fom_cx_p).slice(complex_argument_costfun_m4_m_index_b).row(f_index)) ) ) + ( kron( conj(trans( (*Xhat_fnm_p).slice(complex_argument_costfun_m4_m_index_b).row(f_index) ))%(*expj_Phi_S_nkf_p).slice(f_index).col(k_index) , (*W_fom_cx_p).slice(complex_argument_costfun_m4_m_index_a).row(f_index) ) ) )%( (*Y_lk_p)(l_index, k_index)*kron( trans(complex_argument_costfun_m4_Error_mat_FN.row(f_index))%((*V_nk_p).col(k_index)) , trans((*Z_ol_p).col(l_index)) ) ));*/
T_complex_argument_m4_den_spreadmat_NO[thread_iter]=T_complex_argument_m4_den_spreadmat_NO[thread_iter]+( ( kron( trans((*Xhat_low_fnm_p).slice(complex_argument_costfun_m4_m_index_a).row(f_index)) , (*W_fom_p).slice(complex_argument_costfun_m4_m_index_b).row(f_index) ) ) + ( kron( trans((*Xhat_low_fnm_p).slice(complex_argument_costfun_m4_m_index_b).row(f_index)) , (*W_fom_p).slice(complex_argument_costfun_m4_m_index_a).row(f_index) ) ) )%( (*Y_lk_p)(l_index, k_index)*kron( trans(complex_argument_costfun_m4_Magnitude_model_mat_FN.row(f_index))%((*V_nk_p).col(k_index)) , trans((*Z_ol_p).col(l_index)) ) );

}	

/*compute the numerator part*/
T_complex_argument_m4_num_outmat_FK(f_index, k_index)=2*real(accu( (*V_nk_p).col(k_index)%strans(complex_argument_costfun_m4_Error_mat_FN.row(f_index))%(  strans((*Xhat_fnm_p).slice(complex_argument_costfun_m4_m_index_a).row(f_index))%conj(Xhat_out_4way_tensor_cx_nkfm[complex_argument_costfun_m4_m_index_b].slice(f_index).col(k_index))%conj((*expj_Phi_S_nkf_p).slice(f_index).col(k_index)) + trans((*Xhat_fnm_p).slice(complex_argument_costfun_m4_m_index_b).row(f_index))%(Xhat_out_4way_tensor_cx_nkfm[complex_argument_costfun_m4_m_index_a].slice(f_index).col(k_index))%((*expj_Phi_S_nkf_p).slice(f_index).col(k_index)) ) ));

/*T_complex_argument_m4_num_outmat_FK(f_index, k_index)=2*accu(T_complex_argument_m4_num_spreadmat_NO[thread_iter]);*/
T_complex_argument_m4_den_outmat_FK(f_index, k_index)=2*accu(T_complex_argument_m4_den_spreadmat_NO[thread_iter]);

T_complex_argument_m4_num_outmat_FK(f_index, k_index)=T_complex_argument_m4_num_outmat_FK(f_index, k_index)+T_complex_argument_m4_den_outmat_FK(f_index, k_index);

}

static void T_complex_argument_m4_compute_outstanding_work_last_thread(mat* T_fk_p, arg_struct_t* argStruct_p){

(*T_fk_p)=(*T_fk_p)%((T_complex_argument_m4_num_outmat_FK+outmat_3rd_real_FK_num)/(T_complex_argument_m4_den_outmat_FK+outmat_3rd_real_FK_den));

}

static void T_complex_argument_m4_compute_output_at_pair_indices_complete_work_wrapper(arg_struct_t* argStruct_p, int thread_iter){

int f_index;
int k_index;

bool last_element_queue_flag=false;

bool compute_flag; 

while (threads_while_condition_0th_sum_flag&&(!(last_thread_to_sleep_flag_array[0][thread_iter]))){

	if ((!last_element_found_global_flag)&&(check_queue_global_flag)){

	/*Semaphore or mutex lock */
	sem_wait(&queue_sem);

	compute_flag=calculate_pair_queue_indices(argStruct_p, &f_index, &k_index, &last_element_queue_flag); 

	/*Semaphore or mutex unlock */
	sem_post(&queue_sem);

	}

	if (compute_flag){

	T_complex_argument_m4_compute_output_at_pair_indices_do_work(argStruct_p->Xhat_fnm_p, argStruct_p->Xhat_low_fnm_p, argStruct_p->Z_ol_p, argStruct_p->Y_lk_p, argStruct_p->V_nk_p, argStruct_p->W_fom_p, argStruct_p->W_fom_cx_p, argStruct_p->expj_Phi_S_nkf_p, f_index, k_index, thread_iter);

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

void* T_complex_argument_m4_start(void* arg){

int i_iter;
int thread_iter;

thread_arg_t* threadArg_p;
arg_struct_t* argStruct_p;

threadArg_p=(thread_arg_t*)arg;
argStruct_p=threadArg_p->argStruct_p;
thread_iter=threadArg_p->thread_iter; 

T_complex_argument_m4_compute_output_at_pair_indices_complete_work_wrapper(argStruct_p, thread_iter);

if (last_thread_to_sleep_flag_array[0][thread_iter]){

T_complex_argument_m4_compute_outstanding_work_last_thread(argStruct_p->T_fk_p, argStruct_p);

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