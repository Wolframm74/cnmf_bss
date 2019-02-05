#include "local_inc.hpp"

/*cube::fixed<F_static, N_static, M_static> X_tilde_fnm_abs_cube;
cube::fixed<F_static, N_static, M_static> X_tilde_fnm_binary_mask_cube;*/

/*0th sum stuff*/
cx_cube::fixed<F_static, L_static, M_static> Xhat_outtensor_cx_FLM; 
cube::fixed<F_static, L_static, M_static> Xhat_outtensor_real_FLM;

cx_cube::fixed<F_static, K_static, N_static> Xhat_out_4way_tensor_cx_fknm[M_static];
cx_cube::fixed<N_static, K_static, F_static> Xhat_out_4way_tensor_cx_nkfm[M_static];

cx_cube::fixed<F_static, O_static, M_static> Xhat_outtensor_cx_FOM; 
cx_cube::fixed<F_static, O_static, M_static> Xhat_dummytensor_cx_FOM[NUM_WORKER_THREADS]; 

/*1st sum stuff*/
/*cx_cube::fixed<F_static, M_static, K_static> Xhat_outtensor_cx_FMK;
cube::fixed<F_static, M_static, K_static> Xhat_outtensor_real_FMK;*/
/*cx_cube::fixed<F_static, K_static, M_static> Xhat_outtensor_cx_FKM; 
cube::fixed<F_static, K_static, M_static> Xhat_outtensor_real_FKM;

cx_cube::fixed<M_static, K_static, F_static> Xhat_outtensor_cx_MKF; 	
cube::fixed<M_static, K_static, F_static> Xhat_outtensor_real_MKF;*/

void Xhat_serial_function(mat* Z_ol_p, cx_cube* W_fom_cx_p, cube* W_fom_p){

int m_index;

for (m_index=0; m_index<M_static; m_index++){

Xhat_outtensor_cx_FLM.slice(m_index)=(*W_fom_cx_p).slice(m_index)*(*Z_ol_p);

Xhat_outtensor_real_FLM.slice(m_index)=(*W_fom_p).slice(m_index)*(*Z_ol_p);

}

}

/*spreadmats in LxK*/


/*2nd sum stuff*/

/*queue elements are a pair of indices n, m*/
/*static void compute_Xhat_fnm(cx_cube* Xhat_fnm_p, cube* W_fom_p, mat* Z_ol_p, mat* Y_lk_p, mat* T_fk_p, mat* V_nk_p, cube* Phi_W_fom_p, cube* Phi_S_fnk_p, cube* Xhat_low_fnm_p, int f, int M, int N, int L, int K, int O)*/

/*Need to think if all the other updates, not just Xhat_update can live with indeed shifting M to the far right for Xhat, Xhat_low, and E_conj*/
static void Xhat_primary_auxfun_compute_output_at_pair_indices_do_work(cx_cube* Xhat_fnm_p, cube* Xhat_low_fnm_p, mat* Z_ol_p, mat* Y_lk_p, mat* T_fk_p, mat* V_nk_p, cube* W_fom_p, cx_cube* W_fom_cx_p, cx_cube* expj_Phi_S_nkf_p, cx_cube* E_conj_fnm_p, cx_cube* Xtilde_fnm_p, cx_cube* expj_Phi_W_fom_p, int f_index, int n_index, int thread_iter){

int m_index;

int l_index;

for (m_index=0; m_index<M_static; m_index++){

/*Integrate out L*/

if (Xhat_pop_FKNM_for_V_flag){

Xhat_out_4way_tensor_cx_fknm[m_index].slice(n_index).row(f_index)=((Xhat_outtensor_cx_FLM.slice(m_index).row(f_index))%(expj_Phi_U_flnm[m_index].slice(n_index).row(f_index)))*(*Y_lk_p);	/*assign 1xK=(1xL)*(LxK) into 1xK th output row vector*/

/*Integrate out K, fnth output= (1xK)*(Kx1) */

(*Xhat_fnm_p)(f_index, n_index, m_index)=as_scalar(((Xhat_out_4way_tensor_cx_fknm[m_index].slice(n_index).row(f_index))%((*T_fk_p).row(f_index))%((*expj_Phi_S_nkf_p).slice(f_index).row(n_index)))*trans((*V_nk_p).row(n_index)));

}

/*if (Xhat_pop_NKFM_for_T_or_Phi_S_flag){*/
/*I obsoleted Xhat_pop_NKFM_for_T_or_Phi_S_flag because for example if you called the function with no arguments but the flag Xhat_pop_NKFM_for_T_or_Phi_S_flag was off
then it wouldn't even for example populate Xhat_fnm */

/*This is now the default*/
/*Note that the code will only run as it should if Xhat_pop_FKNM_for_V_flag is set to false!*/
if (!Xhat_pop_FKNM_for_V_flag){

Xhat_out_4way_tensor_cx_nkfm[m_index].slice(f_index).row(n_index)=((Xhat_outtensor_cx_FLM.slice(m_index).row(f_index))%(expj_Phi_U_flnm[m_index].slice(n_index).row(f_index)))*(*Y_lk_p);	/*assign 1xK=(1xL)*(LxK) into 1xK th output row vector*/

/*Integrate out K, fnth output= (1xK)*(Kx1) */

(*Xhat_fnm_p)(f_index, n_index, m_index)=as_scalar(((Xhat_out_4way_tensor_cx_nkfm[m_index].slice(f_index).row(n_index))%((*T_fk_p).row(f_index))%((*expj_Phi_S_nkf_p).slice(f_index).row(n_index)))*trans((*V_nk_p).row(n_index)));

}

/*}*/

/*Compute magnitude model Xhat_low_fnm_p*/

(*Xhat_low_fnm_p)(f_index, n_index, m_index)=accu( (*Y_lk_p)%kron( trans(Xhat_outtensor_real_FLM.slice(m_index).row(f_index)) , ((*T_fk_p).row(f_index))%( (*V_nk_p).row(n_index) ) ) );

(*E_conj_fnm_p)(f_index, n_index, m_index)=conj((*Xtilde_fnm_p)(f_index, n_index, m_index)-(*Xhat_fnm_p)(f_index, n_index, m_index));	


if (Xhat_pop_W_cx_outtensor_flag){

for (l_index=0; l_index<L_static; l_index++){

/*Xhat_dummytensor_cx_FOM[thread_iter].slice(m_index).row(f_index)=Xhat_dummytensor_cx_FOM[thread_iter].slice(m_index).row(f_index)+((*E_conj_fnm_p)(f_index, n_index, m_index))*(expj_Phi_U_flnm[m_index](f_index, l_index, n_index))*sum( kron( trans( (*V_nk_p).row(n_index)%(*T_fk_p).row(f_index)%(*Y_lk_p).row(l_index) ) , ( trans( (*Z_ol_p).col(l_index) )%(*expj_Phi_W_fom_p).slice(m_index).row(f_index) ) ) , 0);*/
Xhat_dummytensor_cx_FOM[thread_iter].slice(m_index).row(f_index)=Xhat_dummytensor_cx_FOM[thread_iter].slice(m_index).row(f_index)+((*E_conj_fnm_p)(f_index, n_index, m_index))*(expj_Phi_U_flnm[m_index](f_index, l_index, n_index))*sum( kron( strans( (*V_nk_p).row(n_index)%(*T_fk_p).row(f_index)%(*Y_lk_p).row(l_index)%((*expj_Phi_S_nkf_p).slice(f_index).row(n_index)) ) , ( trans( (*Z_ol_p).col(l_index) )%(*expj_Phi_W_fom_p).slice(m_index).row(f_index) ) ) , 0);

/*copy the procedure for if( Xhat_pop_NKFM_for_T_or_Phi_S_flag )*/

/*Xhat_out_4way_tensor_cx_nkfm[m_index].slice(f_index).row(n_index)=((Xhat_outtensor_cx_FLM.slice(m_index).row(f_index))%(expj_Phi_U_flnm[m_index].slice(n_index).row(f_index)))*(*Y_lk_p);*/	/*assign 1xK=(1xL)*(LxK) into 1xK th output row vector*/

/*Integrate out K, fnth output= (1xK)*(Kx1) */

/*(*Xhat_fnm_p)(f_index, n_index, m_index)=as_scalar(((Xhat_out_4way_tensor_cx_nkfm[m_index].slice(f_index).row(n_index))%((*T_fk_p).row(f_index))%((*expj_Phi_S_nkf_p).slice(f_index).row(n_index)))*trans((*V_nk_p).row(n_index)));*/

}


}

}



}

static void Xhat_primary_auxfun_compute_outstanding_work_last_thread(void){

int iter_ctr;

Xhat_outtensor_cx_FOM.zeros();

for (iter_ctr=0; iter_ctr<NUM_WORKER_THREADS; iter_ctr++){

Xhat_outtensor_cx_FOM=Xhat_outtensor_cx_FOM+Xhat_dummytensor_cx_FOM[iter_ctr];

}

}

static void Xhat_primary_auxfun_compute_output_at_pair_indices_complete_work_wrapper(arg_struct_t* argStruct_p, int thread_iter){

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
	Xhat_primary_auxfun_compute_output_at_pair_indices_do_work(argStruct_p->Xhat_fnm_p, argStruct_p->Xhat_low_fnm_p, argStruct_p->Z_ol_p, argStruct_p->Y_lk_p, argStruct_p->T_fk_p, argStruct_p->V_nk_p, argStruct_p->W_fom_p, argStruct_p->W_fom_cx_p, argStruct_p->expj_Phi_S_nkf_p, argStruct_p->E_conj_fnm_p, argStruct_p->Xtilde_fnm_p, argStruct_p->expj_Phi_W_fom_p, f_index, n_index, thread_iter);

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


void* Xhat_primary_auxfun_start(void* arg){

int i_iter;
int thread_iter;

thread_arg_t* threadArg_p;
arg_struct_t* argStruct_p;

threadArg_p=(thread_arg_t*)arg;
argStruct_p=threadArg_p->argStruct_p;
thread_iter=threadArg_p->thread_iter; 

Xhat_primary_auxfun_compute_output_at_pair_indices_complete_work_wrapper(argStruct_p, thread_iter);

if (last_thread_to_sleep_flag_array[0][thread_iter]){

/*Xhat_primary_auxfun_compute_outstanding_work_last_thread(argStruct_p->V_nk_p, argStruct_p->Y_lk_p, argStruct_p);*/

if (Xhat_pop_W_cx_outtensor_flag){

Xhat_primary_auxfun_compute_outstanding_work_last_thread();

}

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