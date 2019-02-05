#include "local_inc.hpp"

/*0th sum stuff*/
/*cx_cube::fixed<F_static, M_static, L_static> Xhat_outtensor_cx_FML; 
cube::fixed<F_static, M_static, L_static> Xhat_outtensor_real_FML;*/
cx_cube::fixed<F_static, L_static, M_static> Xhat_outtensor_cx_FLM; 
cube::fixed<F_static, L_static, M_static> Xhat_outtensor_real_FLM;

/*1st sum stuff*/
/*cx_cube::fixed<F_static, M_static, K_static> Xhat_outtensor_cx_FMK;
cube::fixed<F_static, M_static, K_static> Xhat_outtensor_real_FMK;*/
cx_cube::fixed<F_static, K_static, M_static> Xhat_outtensor_cx_FKM; 
cube::fixed<F_static, K_static, M_static> Xhat_outtensor_real_FKM;

cx_cube::fixed<M_static, K_static, F_static> Xhat_outtensor_cx_MKF; 	
cube::fixed<M_static, K_static, F_static> Xhat_outtensor_real_MKF;

/*2nd sum stuff*/

/*queue elements are a pair of indices n, m*/
/*static void compute_Xhat_fnm(cx_cube* Xhat_fnm_p, cube* W_fom_p, mat* Z_ol_p, mat* Y_lk_p, mat* T_fk_p, mat* V_nk_p, cube* Phi_W_fom_p, cube* Phi_S_fnk_p, cube* Xhat_low_fnm_p, int f, int M, int N, int L, int K, int O)*/

/*static void compute_0th_sum(cx_cube* W_fom_cx_p, cube* W_fom_p, cx_mat* Z_ol_p, int m_index){*/
static void compute_0th_sum(cube* W_fom_p, cx_cube* W_fom_cx_p, mat* Z_ol_p, int m_index){

	/*Integrate out o index. FxL=FxO*OxL */
	(Xhat_outtensor_cx_FLM.slice(m_index))=((*W_fom_cx_p).slice(m_index))*(*Z_ol_p);

	(Xhat_outtensor_real_FLM.slice(m_index))=((*W_fom_p).slice(m_index))*(*Z_ol_p);

}

static void compute_1st_sum(mat* Y_lk_p, int m_index){

	/*Integrate out l index. FxK=FxL*LxK */	/*I think obvious that you should shift M right for Xhat_outtensor_cx_FML */
	Xhat_outtensor_cx_FKM.slice(m_index)=(Xhat_outtensor_cx_FLM.slice(m_index))*(*Y_lk_p);

	Xhat_outtensor_real_FKM.slice(m_index)=(Xhat_outtensor_real_FLM.slice(m_index))*(*Y_lk_p);

}

/*Need to think if all the other updates, not just Xhat_update can live with indeed shifting M to the far right for Xhat, Xhat_low, and E_conj*/
static void compute_2nd_sum(cx_cube* Xhat_fnm_p, cx_cube* Xtilde_fnm_p, cube* Xhat_low_fnm_p, cx_cube* E_conj_fnm_p, mat* T_fk_p, mat* V_nk_p, cx_cube* expj_Phi_S_nkf_p, int m_index, int f_index){

	/*Integrate out k index. */	/*This also is benefitted from shifting M right. */
	/*(*Xhat_fnm_p).subcube(f_index,0,m_index,f_index,N_static-1,m_index)=(Xhat_outtensor_cx_FKM.subcube(f_index, 0, m_index, f_index, K_static-1, m_index)%(*T_fk_p).submat(f_index, 0, f_index, K_static-1))*strans( (*expj_Phi_S_fnk_p).subcube(f_index,0,0,f_index, N_static-1, K_static-1) %(*V_nk_p) );*/
	(*Xhat_fnm_p).slice(m_index).row(f_index)=(Xhat_outtensor_cx_FKM.slice(m_index).row(f_index)%(*T_fk_p).row(f_index))*strans( (*expj_Phi_S_nkf_p).slice(f_index)%(*V_nk_p) );

	(*Xhat_low_fnm_p).slice(m_index).row(f_index)=(Xhat_outtensor_real_FKM.slice(m_index).row(f_index)%(*T_fk_p).row(f_index))*trans(*V_nk_p);

	(*E_conj_fnm_p).slice(m_index).row(f_index)=conj((*Xtilde_fnm_p).slice(m_index).row(f_index)-(*Xhat_fnm_p).slice(m_index).row(f_index));

}

static void compute_2nd_sum_pop_MKF(cx_cube* Xhat_fnm_p, cx_cube* Xtilde_fnm_p, cube* Xhat_low_fnm_p, cx_cube* E_conj_fnm_p, mat* T_fk_p, mat* V_nk_p, cx_cube* expj_Phi_S_nkf_p, int m_index, int f_index){
	
	/*Populate MKF tensors*/
	Xhat_outtensor_cx_MKF.slice(f_index).row(m_index)=Xhat_outtensor_cx_FKM.slice(m_index).row(f_index);

	Xhat_outtensor_real_MKF.slice(f_index).row(m_index)=Xhat_outtensor_real_FKM.slice(m_index).row(f_index);

	/*Integrate out k index. */	/*This also is benefitted from shifting M right. */
	/*(*Xhat_fnm_p).subcube(f_index,0,m_index,f_index,N_static-1,m_index)=(Xhat_outtensor_cx_FKM.subcube(f_index, 0, m_index, f_index, K_static-1, m_index)%(*T_fk_p).submat(f_index, 0, f_index, K_static-1))*strans( (*expj_Phi_S_fnk_p).subcube(f_index,0,0,f_index, N_static-1, K_static-1) %(*V_nk_p) );*/
	(*Xhat_fnm_p).slice(m_index).row(f_index)=(Xhat_outtensor_cx_FKM.slice(m_index).row(f_index)%(*T_fk_p).row(f_index))*strans( (*expj_Phi_S_nkf_p).slice(f_index)%(*V_nk_p) );

	(*Xhat_low_fnm_p).slice(m_index).row(f_index)=(Xhat_outtensor_real_FKM.slice(m_index).row(f_index)%(*T_fk_p).row(f_index))*trans(*V_nk_p);

	(*E_conj_fnm_p).slice(m_index).row(f_index)=conj((*Xtilde_fnm_p).slice(m_index).row(f_index)-(*Xhat_fnm_p).slice(m_index).row(f_index));

}


static void Xhat_0th_sum_start(arg_struct_t* argStruct_p, int thread_iter){

int m_index;
bool last_element_queue_flag=false;		/*this may not be asbolutely necessary as before in the legacy code */

bool compute_flag;

while (threads_while_condition_0th_sum_flag&&(!(last_thread_to_sleep_flag_local[0][thread_iter]))){

	if ((!last_element_found_global_flag)&&(check_queue_global_flag)){

	/*Semaphore or mutex lock */
	sem_wait(&queue_sem);

	compute_flag=calculate_queue_index(argStruct_p, &m_index, &last_element_queue_flag);  	/*Check that ourside world is populating queue wrt f=1:F and not n=1:N*/

	/*Semaphore or mutex unlock */
	sem_post(&queue_sem);

	}

	if (compute_flag){

	compute_0th_sum(argStruct_p->W_fom_p, argStruct_p->W_fom_cx_p, argStruct_p->Z_ol_p, m_index);

	sem_wait(&work_actually_computed_sem);
	work_actually_computed_ctr++;

	if (work_actually_computed_ctr==argStruct_p->work_queue_size){

		/*Wake up the last element thread*/
		/*sem_post(&last_element_sleep_sem);*/
		work_actually_completed_flag=true; 

		/*FLAG THIS thread_iter'th THREAD AS THE LAST THREAD TO SLEEP. IE: IT DOESN'T SLEEP. IT GETS TO BREAK THE LOOP ON ITS WAY OUT.*/
		last_thread_to_sleep_flag_local[0][thread_iter]=true;

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

	if (work_actually_completed_flag&&(!(last_thread_to_sleep_flag_local[0][thread_iter]))){

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

static void Xhat_1st_sum_start(arg_struct_t* argStruct_p, int thread_iter){
int m_index;

bool last_element_queue_flag=false;

bool compute_flag;

while (threads_while_condition_1st_sum_flag&&(!(last_thread_to_sleep_flag_local[1][thread_iter]))){

	if ((!last_element_found_global_flag)&&(check_queue_global_flag)){

	/*Semaphore or mutex lock */
	sem_wait(&queue_sem);

	compute_flag=calculate_queue_index(argStruct_p, &m_index, &last_element_queue_flag); 

	/*Semaphore or mutex unlock */
	sem_post(&queue_sem);

	}

	if (compute_flag){

	compute_1st_sum(argStruct_p->Y_lk_p, m_index);

	sem_wait(&work_actually_computed_sem);
	work_actually_computed_ctr++;

	if (work_actually_computed_ctr==argStruct_p->work_queue_size){

		/*Wake up the last element thread*/
		/*sem_post(&last_element_sleep_sem);*/
		work_actually_completed_flag=true; 

		/*FLAG THIS thread_iter'th THREAD AS THE LAST THREAD TO SLEEP. IE: IT DOESN'T SLEEP. IT GETS TO BREAK THE LOOP ON ITS WAY OUT.*/
		last_thread_to_sleep_flag_local[1][thread_iter]=true;

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

	if (work_actually_completed_flag&&(!(last_thread_to_sleep_flag_local[1][thread_iter]))){

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

static void Xhat_2nd_sum_start(arg_struct_t* argStruct_p, int thread_iter){

int f_index;	
int m_index;

bool last_element_queue_flag=false;

bool compute_flag; 

while (threads_while_condition_1st_sum_flag&&(!(last_thread_to_sleep_flag_array[0][thread_iter]))){

	if ((!last_element_found_global_flag)&&(check_queue_global_flag)){

	/*Semaphore or mutex lock */
	sem_wait(&queue_sem);

	compute_flag=calculate_pair_queue_indices(argStruct_p, &m_index, &f_index, &last_element_queue_flag); 

	/*Semaphore or mutex unlock */
	sem_post(&queue_sem);

	}

	if (compute_flag){

	if (Xhat_pop_MKF_for_Phi_S_flag){

	compute_2nd_sum_pop_MKF(argStruct_p->Xhat_fnm_p, argStruct_p->Xtilde_fnm_p, argStruct_p->Xhat_low_fnm_p, argStruct_p->E_conj_fnm_p, argStruct_p->T_fk_p, argStruct_p->V_nk_p, argStruct_p->expj_Phi_S_nkf_p, m_index, f_index);

	} else {

	compute_2nd_sum(argStruct_p->Xhat_fnm_p, argStruct_p->Xtilde_fnm_p, argStruct_p->Xhat_low_fnm_p, argStruct_p->E_conj_fnm_p, argStruct_p->T_fk_p, argStruct_p->V_nk_p, argStruct_p->expj_Phi_S_nkf_p, m_index, f_index);

	}

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

void* Xhat_start_M_threads(void* arg){

thread_arg_t* threadArg_p;
arg_struct_t* argStruct_p;

bool work_not_yet_completed_flag=true;
bool note_to_sleep_flag=false;

threadArg_p=(thread_arg_t*)arg;
argStruct_p=threadArg_p->argStruct_p;

Xhat_0th_sum_start(argStruct_p, thread_iter);

if (last_thread_to_sleep_flag_local[0][thread_iter]){

sem_wait(&threads_while_condition_0th_sum_sem);
threads_while_condition_0th_sum_flag=false;	
sem_post(&threads_while_condition_0th_sum_sem);

/*Sem post to wake M_static-1 of the sleeping threads*/
sem_post(&sleep_sem);

/*Sem wait M_static-1 times so that you only proceed from this point once that many of the sleeping flags have high fived and successfully made it out of the "while traps" */
sem_wait(&checkpoint_sem);

}

if (!last_thread_to_sleep_flag_local[0][thread_iter]){

sem_post(checkpoint_sem);
sem_wait(&wait_till_queue_populated);

}

if (last_thread_to_sleep_flag_local[0][thread_iter]){

/*Repopulate queue*/
populate_queue_wrt_one_index(argStruct_p, M_static);

/*I think also reset: last_element_found_global_flag to false */
/*Consider creating a semaphore for this and wrapping it, but for now do without*/
last_element_found_global_flag=false;

/*Reset*/
work_actually_completed_flag=false;
work_actually_computed_ctr=0;	

sem_wait(&check_queue_global_sem);
check_queue_global_flag=true;
sem_post(&check_queue_global_sem);

/*Let the threads that are waiting until queue is populated proceed: call sem_post() three times*/	
sem_post(&wait_till_queue_populated);

}

Xhat_1st_sum_start(argStruct_p, thread_iter);

if (last_thread_to_sleep_flag_local[1][thread_iter]){

sem_wait(&threads_while_condition_1st_sum_sem);
threads_while_condition_1st_sum_flag=false;	
sem_post(&threads_while_condition_1st_sum_sem);

/*Sem post to wake M_static-1 of the sleeping threads*/
sem_post(&sleep_sem);

/*Sem wait M_static-1 times so that you only proceed from this point once that many of the sleeping flags have high fived and successfully made it out of the "while traps" */
sem_wait(&checkpoint_sem);

}

if (!last_thread_to_sleep_flag_local[1][thread_iter]){

sem_post(checkpoint_sem);
sem_wait(&wait_till_queue_populated);

}

if (last_thread_to_sleep_flag_local[1][thread_iter]){

/*Repopulate queue*/
populate_queue_wrt_pair_indices(argStruct_p, M_static, F_static);

/*I think also reset: last_element_found_global_flag to false */
/*Consider creating a semaphore for this and wrapping it, but for now do without*/
last_element_found_global_flag=false;

/*Reset*/
work_actually_completed_flag=false;
work_actually_computed_ctr=0;	

sem_wait(&check_queue_global_sem);
check_queue_global_flag=true;
sem_post(&check_queue_global_sem);

/*Let the threads that are waiting until queue is populated proceed: call sem_post() three times*/	
sem_post(&wait_till_queue_populated);

}

}

void* Xhat_start(void* arg){

int thread_iter;

thread_arg_t* threadArg_p;
arg_struct_t* argStruct_p;

threadArg_p=(thread_arg_t*)arg;
argStruct_p=threadArg_p->argStruct_p;
thread_iter=threadArg_p->thread_iter; 

Xhat_2nd_sum_start(argStruct_p, thread_iter);

if (last_thread_to_sleep_flag_array[0][thread_iter]){

sem_wait(&threads_while_condition_2nd_sum_sem);
threads_while_condition_2nd_sum_flag=false;	
sem_post(&threads_while_condition_2nd_sum_sem);

/*Sem post to wake NUM_WORKER_THREADS-1 of the sleeping threads*/
sem_post(&sleep_sem);

/*Sem wait NUM_WORKER_THREADS-1 times so that you only proceed from this point once that many of the sleeping flags have high fived and successfully made it out of the "while traps" */
sem_wait(&checkpoint_sem);

}

if (!last_thread_to_sleep_flag_array[0][thread_iter]){

sem_post(checkpoint_sem);

}

threads_exit_and_signal(NUM_WORKER_THREADS);

}