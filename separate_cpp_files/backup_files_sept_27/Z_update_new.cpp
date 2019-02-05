#include "local_inc.hpp"

/*0th sum stuff*/
/*cube::fixed<F_static, N_static, L_static> outtensor_0th_real_FNL;
cx_cube::fixed<F_static, N_static, L_static> outtensor_0th_cx_FNL;*/
mat::fixed<K_static, L_static> Z_inmat_real_Ykl;

/*cube::fixed<M_static, N_static, F_static> Intensor_Xhat_low_mnf;
cx_cube::fixed<M_static, N_static, F_static> Intensor_E_conj_mnf;*/

cube::fixed<N_static, L_static, F_static> outtensor_0th_real_NLF;
cx_cube::fixed<N_static, L_static, F_static> outtensor_0th_cx_NLF;

static cx_mat::fixed<N_static, K_static> dummy_mat_cx_NK[NUM_WORKER_THREADS]; 
static mat::fixed<K_static, K_static> dummy_mat_outerprod_KK[NUM_WORKER_THREADS]; 

static mat::fixed<N_static, K_static> dummy_mat_real_NK[NUM_WORKER_THREADS]; 

/*1st sum stuff*/
cube::fixed<L_static, F_static, M_static> outtensor_1st_real_LFM;
cx_cube::fixed<L_static, F_static, M_static> outtensor_1st_cx_LFM;

/*2nd sum stuff*/
cube::fixed<O_static, L_static, M_static> outtensor_2nd_real_OLM;
cx_cube::fixed<O_static, L_static, M_static> outtensor_2nd_cx_OLM;

/*3rd sum stuff*/
mat::fixed<O_static, L_static> outmat_3rd_real_OL_den;
mat::fixed<O_static, L_static> outmat_3rd_real_OL_num;
cx_mat::fixed<O_static, L_static> outmat_3rd_cx_OL;

static void compute_0th_sum(mat* T_fk_p, mat* V_nk_p, cx_cube* expj_Phi_S_nkf_p, int f_index, int thread_iter){

#ifdef DEBUG	
mexPrintf("Z_update: pthread_self(): %d, compute_0th_sum(), f_index=%d, work_actually_computed_ctr=%d \n", (int)pthread_self(), f_index, work_actually_computed_ctr);	
#endif

/*KxK outer product computation*/	/*fine, I think*/
dummy_mat_outerprod_KK[thread_iter]=trans((*T_fk_p).row(f_index))*((*T_fk_p).row(f_index));

/*outtensor_0th_real_FNL computation*/	/*fine, I think*/
/*outtensor_0th_real_FNL.subcube(0, n_index, 0, F_static-1, n_index, L_static-1)=(*T_fk_p)*dummy_mat_outerprod_KK[thread_iter]*trans(*Y_lk_p);*/
outtensor_0th_real_NLF.slice(f_index)=(*V_nk_p)*dummy_mat_outerprod_KK[thread_iter]*Z_inmat_real_Ykl;

/*Prepping dummy_mat_real_NK*/
dummy_mat_real_NK[thread_iter]=kron(ones_col_Nx1, (*T_fk_p).row(f_index)); 

/*outtensor_0th_cx_FNL computation*/	/*Might be ok, since you're only doing elementwise stuff on expj_Phi_S_fnk_p*/

dummy_mat_cx_NK[thread_iter]=((*expj_Phi_S_nkf_p).slice(f_index))%(*V_nk_p);

/*dummy_mat_cx_NK[thread_iter]=dummy_mat_cx_NK[thread_iter].each_row()%(*T_fk_p).row(f_index);*/
dummy_mat_cx_NK[thread_iter]=(dummy_mat_cx_NK[thread_iter])%(dummy_mat_real_NK[thread_iter]);

/*outtensor_0th_cx_FNL.subcube(0, n_index, 0, F_static-1, n_index, L_static-1)=dummy_mat_cx_NK[thread_iter]*trans(*Y_lk_p);*/
outtensor_0th_cx_NLF.slice(f_index)=dummy_mat_cx_NK[thread_iter]*Z_inmat_real_Ykl;

}

/*static void compute_1st_sum(cx_cube* Xhat_low_fnm_p, cx_cube* E_conj_fnm_p, int m_index, int f_index){*/
static void compute_1st_sum(cube* Xhat_low_fnm_p, cx_cube* E_conj_fnm_p, int m_index, int f_index){	

/*outtensor_1st_real_LFM computation*/

/*Integrate out n. 1xL=1xN*NxL */	/*Xhat_low_fnm_p: represent differently*/	/*If possible, F should be moved to the very far right */
/*outtensor_1st_real_LFM.subcube(f_index, m_index, 0, f_index, m_index, L_static-1)=(*Xhat_low_fnm_p).subcube(f_index, 0, m_index, f_index, N_static-1, m_index)*outtensor_0th_real_FNL.subcube(f_index,0,0,f_index,N_static-1,L_static-1);	*/
/*outtensor_1st_real_LFM.slice(m_index).col(f_index)=trans(Intensor_Xhat_low_mnf.slice(f_index).row(m_index)*outtensor_0th_real_NLF.slice(f_index));*/
#ifdef DEBUG	
mexPrintf("Z_update: pthread_self(): %d, compute_1st_sum(), m_index=%d, f_index=%d, work_actually_computed_ctr=%d \n", (int)pthread_self(), m_index, f_index, work_actually_computed_ctr);	
#endif

outtensor_1st_real_LFM.slice(m_index).col(f_index)=trans(((*Xhat_low_fnm_p).slice(m_index).row(f_index))*outtensor_0th_real_NLF.slice(f_index));

/*outtensor_1st_cx_LFM computation*/	/*E_conj_fnm_p represent differently*/	/*Likewise, shift F */
/*outtensor_1st_cx_LFM.subcube(f_index, m_index, 0, f_index, m_index, L_static-1)=(*E_conj_fnm_p).subcube(f_index, 0, m_index, f_index, N_static-1, m_index)*outtensor_0th_cx_FNL.subcube(f_index,0,0,f_index,N_static-1,L_static-1);	*/
/*outtensor_1st_cx_LFM.slice(m_index).col(f_index)=strans(Intensor_E_conj_mnf.slice(f_index).row(m_index)*outtensor_0th_cx_NLF.slice(f_index));	*/
outtensor_1st_cx_LFM.slice(m_index).col(f_index)=strans(((*E_conj_fnm_p).slice(m_index).row(f_index))*outtensor_0th_cx_NLF.slice(f_index));	

}

static void compute_2nd_sum(cube* W_fom_p, cx_cube* W_fom_cx_p, int m_index, int l_index){

#ifdef DEBUG	
mexPrintf("Z_update: pthread_self(): %d, compute_2nd_sum(), m_index=%d, l_index=%d, work_actually_computed_ctr=%d \n", (int)pthread_self(), m_index, l_index, work_actually_computed_ctr);	
#endif

/*Integrate out f, put the result in the real output tensor. 1xO=1xF*FxO*/	/*outtensor_1st_real_LFM might be ok as is. W_fom needs to be rep'd differently */
/*outtensor_2nd_real_MOL.subcube(m_index,0,l_index,m_index, O_static-1,l_index)=outtensor_1st_real_LFM.subcube(0,m_index,l_index, F_static-1,m_index,l_index)*(*W_fom_p).subcube(0,0,m_index,F_static-1, O_static-1,m_index);*/
outtensor_2nd_real_OLM.slice(m_index).col(l_index)=trans(outtensor_1st_real_LFM.slice(m_index).row(l_index)*(*W_fom_p).slice(m_index));

/*Integrate out f for the cx output tensor. */	
/*outtensor_2nd_cx_MOL.subcube(m_index,0,l_index,m_index, O_static-1,l_index)=outtensor_1st_cx_LFM.subcube(0,m_index,l_index, F_static-1,m_index,l_index)*(*W_fom_cx_p).subcube(0,0,m_index,F_static-1, O_static-1,m_index);*/
outtensor_2nd_cx_OLM.slice(m_index).col(l_index)=strans(outtensor_1st_cx_LFM.slice(m_index).row(l_index)*(*W_fom_cx_p).slice(m_index));

}


/*With the exception of the W_fom update, the computation of the 3rd sum involves only M components, and therefore no need to set up a work queue for the computations for this*/
static void compute_3rd_sum_last_thread(mat* Z_ol_p, mat* Y_lk_p){

#ifdef Z_UNIT_DEBUG

outtensor_1st_cx_LFM.slice(0).print("outtensor_1st_cx_LFM.slice(0):");
outtensor_1st_cx_LFM.slice(1).print("outtensor_1st_cx_LFM.slice(1):");

outtensor_2nd_cx_OLM.slice(0).print("outtensor_2nd_cx_OLM.slice(0):");
outtensor_2nd_cx_OLM.slice(1).print("outtensor_2nd_cx_OLM.slice(1):");

#endif

/*Do the sum for m=0:M_static-1*/
/*outmat_3rd_real_OL_den=outtensor_2nd_real_MOL.subcube(0,0,0,0,O_static-1,L_static-1)+outtensor_2nd_real_MOL.subcube(1,0,0,1,O_static-1,L_static-1);*/
outmat_3rd_real_OL_den=outtensor_2nd_real_OLM.slice(0)+outtensor_2nd_real_OLM.slice(1);	

#ifdef Z_UNIT_DEBUG

outmat_3rd_real_OL_den.print("outmat_3rd_real_OL_den:");

outtensor_2nd_cx_OLM.slice(0).print("outtensor_2nd_cx_OLM.slice(0):");

outtensor_2nd_cx_OLM.slice(1).print("outtensor_2nd_cx_OLM.slice(1):");

#endif

outmat_3rd_cx_OL=outtensor_2nd_cx_OLM.slice(0)+outtensor_2nd_cx_OLM.slice(1);


outmat_3rd_real_OL_num=outmat_3rd_real_OL_den+real(outmat_3rd_cx_OL);

#ifdef Z_UNIT_DEBUG

outmat_3rd_cx_OL.print("outmat_3rd_cx_OL:");

outmat_3rd_real_OL_num.print("outmat_3rd_real_OL_num:");

#endif

(*Z_ol_p)=(*Z_ol_p)%(outmat_3rd_real_OL_num/outmat_3rd_real_OL_den);

/*Project Z_ol onto the nonnegative orthant. ie: zero all negative elements*/
(*Z_ol_p).elem(find((*Z_ol_p)<0)).zeros();

/*Call normalize_Z_ol() */
normalize_Z_ol(Z_ol_p, Y_lk_p);

}

static void Z_0th_sum_start(arg_struct_t* argStruct_p, int thread_iter){

int f_index;
bool last_element_queue_flag=false;		/*this may not be asbolutely necessary as before in the legacy code */

bool compute_flag;

while (threads_while_condition_0th_sum_flag&&(!(last_thread_to_sleep_flag_array[0][thread_iter]))){

	if ((!last_element_found_global_flag)&&(check_queue_global_flag)){

	/*Semaphore or mutex lock */
	sem_wait(&queue_sem);

	compute_flag=calculate_queue_index(argStruct_p, &f_index, &last_element_queue_flag); 	/*Check that ourside world is populating queue wrt f=1:F and not n=1:N*/

	/*Semaphore or mutex unlock */
	sem_post(&queue_sem);

	}

	if (compute_flag){

	compute_0th_sum(argStruct_p->T_fk_p, argStruct_p->V_nk_p, argStruct_p->expj_Phi_S_nkf_p, f_index, thread_iter);

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

static void Z_1st_sum_start(arg_struct_t* argStruct_p, int thread_iter){
int m_index;
int f_index;

bool last_element_queue_flag=false;

bool compute_flag;

while (threads_while_condition_1st_sum_flag&&(!(last_thread_to_sleep_flag_array[1][thread_iter]))){

	if ((!last_element_found_global_flag)&&(check_queue_global_flag)){

	/*Semaphore or mutex lock */
	sem_wait(&queue_sem);

	compute_flag=calculate_pair_queue_indices(argStruct_p, &m_index, &f_index, &last_element_queue_flag); 

	/*Semaphore or mutex unlock */
	sem_post(&queue_sem);

	}

	if (compute_flag){

	compute_1st_sum(argStruct_p->Xhat_low_fnm_p, argStruct_p->E_conj_fnm_p, m_index, f_index);	

	sem_wait(&work_actually_computed_sem);
	work_actually_computed_ctr++;

	if (work_actually_computed_ctr==argStruct_p->work_queue_size){

		/*Wake up the last element thread*/
		/*sem_post(&last_element_sleep_sem);*/
		work_actually_completed_flag=true; 

		/*FLAG THIS thread_iter'th THREAD AS THE LAST THREAD TO SLEEP. IE: IT DOESN'T SLEEP. IT GETS TO BREAK THE LOOP ON ITS WAY OUT.*/
		last_thread_to_sleep_flag_array[1][thread_iter]=true;

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

	if (work_actually_completed_flag&&(!(last_thread_to_sleep_flag_array[1][thread_iter]))){

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

static void Z_2nd_sum_start(arg_struct_t* argStruct_p, int thread_iter){

int l_index;	
int m_index;

bool last_element_queue_flag=false;

bool compute_flag; 

while (threads_while_condition_1st_sum_flag&&(!(last_thread_to_sleep_flag_array[2][thread_iter]))){

	if ((!last_element_found_global_flag)&&(check_queue_global_flag)){

	/*Semaphore or mutex lock */
	sem_wait(&queue_sem);

	compute_flag=calculate_pair_queue_indices(argStruct_p, &m_index, &l_index, &last_element_queue_flag); 

	/*Semaphore or mutex unlock */
	sem_post(&queue_sem);

	}

	if (compute_flag){

	compute_2nd_sum(argStruct_p->W_fom_p, argStruct_p->W_fom_cx_p, m_index, l_index);

	sem_wait(&work_actually_computed_sem);
	work_actually_computed_ctr++;

	if (work_actually_computed_ctr==argStruct_p->work_queue_size){

		/*Wake up the last element thread*/
		/*sem_post(&last_element_sleep_sem);*/
		work_actually_completed_flag=true; 

		/*FLAG THIS thread_iter'th THREAD AS THE LAST THREAD TO SLEEP. IE: IT DOESN'T SLEEP. IT GETS TO BREAK THE LOOP ON ITS WAY OUT.*/
		last_thread_to_sleep_flag_array[2][thread_iter]=true;

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

	if (work_actually_completed_flag&&(!(last_thread_to_sleep_flag_array[2][thread_iter]))){

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

void* Z_start(void* arg){

int thread_iter;

thread_arg_t* threadArg_p;
arg_struct_t* argStruct_p;

threadArg_p=(thread_arg_t*)arg;
argStruct_p=threadArg_p->argStruct_p;
thread_iter=threadArg_p->thread_iter; 

Z_0th_sum_start(argStruct_p, thread_iter);

if (last_thread_to_sleep_flag_array[0][thread_iter]){

sem_wait(&threads_while_condition_0th_sum_sem);
threads_while_condition_0th_sum_flag=false;	
sem_post(&threads_while_condition_0th_sum_sem);

/*Sem post to wake NUM_WORKER_THREADS-1 of the sleeping threads*/
sem_post(&sleep_sem);
sem_post(&sleep_sem);
sem_post(&sleep_sem);

/*Sem wait NUM_WORKER_THREADS-1 times so that you only proceed from this point once that many of the sleeping flags have high fived and successfully made it out of the "while traps" */
sem_wait(&checkpoint_sem);
sem_wait(&checkpoint_sem);
sem_wait(&checkpoint_sem);

}

if (!last_thread_to_sleep_flag_array[0][thread_iter]) {

sem_post(&checkpoint_sem);
sem_wait(&wait_till_queue_populated);

}

if (last_thread_to_sleep_flag_array[0][thread_iter]){

/*Repopulate queue wrt M, F*/
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
sem_post(&wait_till_queue_populated);
sem_post(&wait_till_queue_populated);

}

/*non last thread to sleep threads will wait as a first step once they get into the function below*/
Z_1st_sum_start(argStruct_p, &thread_iter);

if (last_thread_to_sleep_flag_array[1][thread_iter]){

sem_wait(&threads_while_condition_1st_sum_sem);
threads_while_condition_1st_sum_flag=false;	
sem_post(&threads_while_condition_1st_sum_sem);

/*Sem post to wake NUM_WORKER_THREADS-1 of the sleeping threads*/
sem_post(&sleep_sem);
sem_post(&sleep_sem);
sem_post(&sleep_sem);

/*Sem wait NUM_WORKER_THREADS-1 times so that you only proceed from this point once that many of the sleeping flags have high fived and successfully made it out of the "while traps" */
sem_wait(&checkpoint_sem);
sem_wait(&checkpoint_sem);
sem_wait(&checkpoint_sem);

}

if (!last_thread_to_sleep_flag_array[1][thread_iter]){

sem_post(&checkpoint_sem);
sem_wait(&wait_till_queue_populated);

}

if (last_thread_to_sleep_flag_array[1][thread_iter]){

/*Repopulate queue wrt M, L */
populate_queue_wrt_pair_indices(argStruct_p, M_static, L_static);

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
sem_post(&wait_till_queue_populated);
sem_post(&wait_till_queue_populated);

}

Z_2nd_sum_start(argStruct_p, &thread_iter);
	
if (last_thread_to_sleep_flag_array[2][thread_iter]){

sem_wait(&threads_while_condition_2nd_sum_sem);
threads_while_condition_2nd_sum_flag=false;	
sem_post(&threads_while_condition_2nd_sum_sem);

/*Sem post to wake NUM_WORKER_THREADS-1 of the sleeping threads*/
sem_post(&sleep_sem);
sem_post(&sleep_sem);
sem_post(&sleep_sem);

/*Sem wait NUM_WORKER_THREADS-1 times so that you only proceed from this point once that many of the sleeping flags have high fived and successfully made it out of the "while traps" */
sem_wait(&checkpoint_sem);
sem_wait(&checkpoint_sem);
sem_wait(&checkpoint_sem);

}

if (!last_thread_to_sleep_flag_array[2][thread_iter]){

sem_post(&checkpoint_sem);

}

threads_exit_and_signal(NUM_WORKER_THREADS);

}
