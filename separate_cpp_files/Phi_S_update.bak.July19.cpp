#include "local_inc.hpp"

/*cube::fixed<F_static, K_static, N_static> ones_tensor_real_FKN;
cube::fixed<N_static, K_static, F_static> ones_tensor_real_NKF;*/

/*cube::fixed<F_static, K_static, N_static> outtensor_real_FKN;*/
cube::fixed<N_static, K_static, F_static> outtensor_real_NKF;

cx_cube::fixed<F_static, K_static, N_static> outtensor_target_cx_FKN;	/*let this one be related to the real tensor*/
cx_cube::fixed<N_static, K_static, F_static> outtensor_target_cx_NKF; 

void Phi_S_update_auxfun1_zero_internal_quantities(void){

outtensor_real_NKF.zeros();
outtensor_target_cx_FKN.zeros();
outtensor_target_cx_NKF.zeros();

}

static void Phi_S_primary_auxfun_compute_output_at_pair_indices_do_work(mat* T_fk_p, mat* V_nk_p, cube* Xhat_low_fnm_p, cx_cube* E_conj_fnm_p, cx_cube* expj_Phi_S_nkf_p, cx_cube* expj_Phi_S_fkn_p, int f_index, int n_index, int thread_iter, cube* W_fom_p, mat* Z_ol_p, mat* Y_lk_p){

int m_index, k_index; 

for (m_index=0; m_index<M_static; m_index++){

outtensor_real_NKF.slice(f_index).row(n_index)=outtensor_real_NKF.slice(f_index).row(n_index)+((*Xhat_low_fnm_p)(f_index, n_index, m_index))*((*T_fk_p).row(f_index)%(*V_nk_p).row(n_index))%((Xhat_outtensor_real_FLM.slice(m_index).row(f_index))*(*Y_lk_p));

/*outtensor_target_cx_NKF.slice(f_index).row(n_index)=outtensor_target_cx_NKF.slice(f_index).row(n_index)+conj((*E_conj_fnm_p)(f_index, n_index, m_index))*conj(Xhat_out_4way_tensor_cx_nkfm[m_index].slice(f_index).row(n_index));*/
outtensor_target_cx_NKF.slice(f_index).row(n_index)=outtensor_target_cx_NKF.slice(f_index).row(n_index)+((*E_conj_fnm_p)(f_index, n_index, m_index))*(Xhat_out_4way_tensor_cx_nkfm[m_index].slice(f_index).row(n_index));

}

outtensor_target_cx_NKF.slice(f_index).row(n_index)=outtensor_target_cx_NKF.slice(f_index).row(n_index)%((*T_fk_p).row(f_index)%(*V_nk_p).row(n_index));

/*outtensor_target_cx_NKF.slice(f_index).row(n_index)=-((conj((*expj_Phi_S_fkn_p).slice(n_index).row(f_index)))%outtensor_target_cx_NKF.slice(f_index).row(n_index)+2*outtensor_real_NKF.slice(f_index).row(n_index));*/

outtensor_target_cx_NKF.slice(f_index).row(n_index)=(outtensor_target_cx_NKF.slice(f_index).row(n_index)+(conj((*expj_Phi_S_fkn_p).slice(n_index).row(f_index)))%outtensor_real_NKF.slice(f_index).row(n_index));

outtensor_target_cx_FKN.slice(n_index).row(f_index)=outtensor_target_cx_NKF.slice(f_index).row(n_index);

}

static void Phi_S_primary_auxfun_compute_output_at_pair_indices_complete_work_wrapper(arg_struct_t* argStruct_p, int thread_iter){

int f_index;
int n_index;

bool last_element_queue_flag=false;

bool compute_flag; 

while (threads_while_condition_0th_sum_flag&&(!(last_thread_to_sleep_flag_array[0][thread_iter]))){

	if ((!last_element_found_global_flag)&&(check_queue_global_flag)){

	/*Semaphore or mutex lock */
	sem_wait(&queue_sem);

	compute_flag=calculate_pair_queue_indices(argStruct_p, &f_index, &n_index, &last_element_queue_flag); 

	/*Semaphore or mutex unlock */
	sem_post(&queue_sem);

	}

	if (compute_flag){

	/*Xhat_primary_auxfun_compute_output_at_pair_indices_do_work(cx_cube* Xhat_fnm_p, cube* Xhat_low_fnm_p, mat* Z_ol_p, mat* Y_lk_p, mat* T_fk_p, mat* V_nk_p, cube* W_fom_p, cx_cube* W_fom_cx_p, cx_cube* expj_Phi_S_nkf_p, cx_cube* E_conj_fnm_p, cx_cube* Xtilde_fnm_p, cx_cube* expj_Phi_W_fom_p, int f_index, n_index, thread_iter);		*/
	Phi_S_primary_auxfun_compute_output_at_pair_indices_do_work(argStruct_p->T_fk_p, argStruct_p->V_nk_p, argStruct_p->Xhat_low_fnm_p, argStruct_p->E_conj_fnm_p, argStruct_p->expj_Phi_S_nkf_p, argStruct_p->expj_Phi_S_fkn_p, f_index, n_index, thread_iter, argStruct_p->W_fom_p, argStruct_p->Z_ol_p, argStruct_p->Y_lk_p);

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

static cube::fixed<F_static, K_static, N_static> absval_tensor_FKN;
static cube::fixed<N_static, K_static, F_static> absval_tensor_NKF;

void* Phi_S_primary_auxfun_start(void* arg){

double epsilon_value;

double mu_value;

mu_value=0.1;

int i_iter;
int thread_iter;
cx_double local_cx_double;

thread_arg_t* threadArg_p;
arg_struct_t* argStruct_p;

threadArg_p=(thread_arg_t*)arg;
argStruct_p=threadArg_p->argStruct_p;
thread_iter=threadArg_p->thread_iter; 

epsilon_value=0.000000000001;

Phi_S_primary_auxfun_compute_output_at_pair_indices_complete_work_wrapper(argStruct_p, thread_iter);

if (last_thread_to_sleep_flag_array[0][thread_iter]){

/*POSTPONE THIS~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

absval_tensor_NKF=abs(outtensor_target_cx_NKF);
absval_tensor_FKN=abs(outtensor_target_cx_FKN);

outtensor_target_cx_NKF.elem(find(absval_tensor_NKF==0)).fill(0.00000001);
outtensor_target_cx_FKN.elem(find(absval_tensor_FKN==0)).fill(0.00000001);

(absval_tensor_NKF).elem(find(absval_tensor_NKF<=0)).fill(0.00000001);
(absval_tensor_FKN).elem(find(absval_tensor_FKN<=0)).fill(0.00000001);

*(argStruct_p->expj_Phi_S_nkf_p)=conj(outtensor_target_cx_NKF/absval_tensor_NKF); //%phase_shift_tensor_cx_NKF;
*(argStruct_p->expj_Phi_S_fkn_p)=conj(outtensor_target_cx_FKN/absval_tensor_FKN); //%phase_shift_tensor_cx_FKN;

outtensor_target_cx_NKF.zeros();
outtensor_target_cx_FKN.zeros();

/*(*(argStruct_p->expj_Phi_S_nkf_p))=(*(argStruct_p->expj_Phi_S_nkf_p))+2*mu_value*(outtensor_target_cx_NKF);
(*(argStruct_p->expj_Phi_S_fkn_p))=(*(argStruct_p->expj_Phi_S_fkn_p))+2*mu_value*(outtensor_target_cx_FKN);*/

/*outtensor_target_cx_NKF=-conj((*(argStruct_p->expj_Phi_S_nkf_p))%outtensor_target_cx_NKF);
outtensor_target_cx_FKN=-conj((*(argStruct_p->expj_Phi_S_fkn_p))%outtensor_target_cx_FKN);*/

/*output tensors set*/

sem_wait(&threads_while_condition_0th_sum_sem);
threads_while_condition_0th_sum_flag=false;	
sem_post(&threads_while_condition_0th_sum_sem);

/*Sem post to wake NUM_WORKER_THREADS-1 of the sleeping threads*/
/*sem_post(&sleep_sem);
sem_post(&sleep_sem);
sem_post(&sleep_sem);*/

for (i_iter=0; i_iter<(NUM_WORKER_THREADS-1); i_iter++){
	sem_post(&sleep_sem);
}

/*Sem wait NUM_WORKER_THREADS-1 times so that you only proceed from this point once that many of the sleeping flags have high fived and successfully made it out of the "while traps" */
/*sem_wait(&checkpoint_sem);
sem_wait(&checkpoint_sem);
sem_wait(&checkpoint_sem);*/

for (i_iter=0; i_iter<(NUM_WORKER_THREADS-1); i_iter++){
	sem_wait(&checkpoint_sem);
}

}

if (!last_thread_to_sleep_flag_array[0][thread_iter]){

sem_post(&checkpoint_sem);

}

threads_exit_and_signal(NUM_WORKER_THREADS);
     
}


static rowvec::fixed<K_static> L1_norm_rowvec;

void Phi_S_L1_norm_V_checker(cx_cube* expj_Phi_S_nkf_p, cx_cube* expj_Phi_S_fkn_p, mat* V_nk_p){

L1_norm_rowvec=sum(*V_nk_p, 0);

int k_index, f_index, n_index;

for (k_index=0; k_index<K_static; k_index++){

	if (L1_norm_rowvec(k_index)<0.3){

	for (f_index=0; f_index<F_static; f_index++){

		for (n_index=0; n_index<N_static; n_index++){

			(*expj_Phi_S_nkf_p)(n_index, k_index, f_index)=1;

			(*expj_Phi_S_fkn_p)(f_index, k_index, n_index)=1;

		}

	}

	}

}


}