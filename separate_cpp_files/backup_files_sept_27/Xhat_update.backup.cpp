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

#ifdef DEBUG	
mexPrintf("Xhat_update: pthread_self(): %d, compute_0th_sum(), m_index=%d \n", (int)pthread_self(), m_index);		
#endif	

	/*Integrate out o index. FxL=FxO*OxL */
	(Xhat_outtensor_cx_FLM.slice(m_index))=((*W_fom_cx_p).slice(m_index))*(*Z_ol_p);

	(Xhat_outtensor_real_FLM.slice(m_index))=((*W_fom_p).slice(m_index))*(*Z_ol_p);

}

static void compute_1st_sum(mat* Y_lk_p, int m_index){

#ifdef DEBUG	
mexPrintf("Xhat_update: pthread_self(): %d, compute_1st_sum(), m_index=%d \n", (int)pthread_self(), m_index);		
#endif	
	
	/*Integrate out l index. FxK=FxL*LxK */	/*I think obvious that you should shift M right for Xhat_outtensor_cx_FML */
	Xhat_outtensor_cx_FKM.slice(m_index)=(Xhat_outtensor_cx_FLM.slice(m_index))*(*Y_lk_p);

	Xhat_outtensor_real_FKM.slice(m_index)=(Xhat_outtensor_real_FLM.slice(m_index))*(*Y_lk_p);

}

/*Need to think if all the other updates, not just Xhat_update can live with indeed shifting M to the far right for Xhat, Xhat_low, and E_conj*/
static void compute_2nd_sum(cx_cube* Xhat_fnm_p, cx_cube* Xtilde_fnm_p, cube* Xhat_low_fnm_p, cx_cube* E_conj_fnm_p, mat* T_fk_p, mat* V_nk_p, cx_cube* expj_Phi_S_nkf_p, int m_index, int f_index){

#ifdef DEBUG	
mexPrintf("Xhat_update: pthread_self(): %d, compute_2nd_sum(), m_index=%d, f_index=%d \n", (int)pthread_self(), m_index, f_index);		
#endif	
	
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


static void Xhat_0th_sum_start(arg_struct_t* argStruct_p, bool* last_element_p){

/*Local variables*/
/*int l_index;*/
int m_index;
bool last_thread_to_sleep=false;

/*bool b_exit_cond_local=false; */

bool compute_flag; 
/*bool last_element=false;*/

while(threads_while_condition_0th_sum_flag&&(!(*last_element_p))) {		/*Change b_exit_cond_global to b_exit_cond_global_0th_sum?? */

	if (!(*last_element_p)) {	/*Assuming you make the outer while a function of last_element_p, then you can remove this. Can re-enable it if and when you're investigating bugs/debugging. */

	if ((!last_element_found_global_flag)&&(check_queue_global_flag)){

	/*Semaphore or mutex lock */
	sem_wait(&queue_sem);

	compute_flag=calculate_queue_index(argStruct_p, &m_index, last_element_p); 

	mexPrintf("Xhat_0th_sum_start: pthread_self(): %d returned from calculate_queue_index() with m_index= %d, compute_flag=%d, threads_asleep_ctr=%d, last_element flag = %d, work_actually_computed_ctr=%d \n", (int) pthread_self(), m_index, (int)compute_flag, threads_asleep_ctr, (int)(*last_element_p), work_actually_computed_ctr);

	/*Semaphore or mutex unlock */
	sem_post(&queue_sem);

	}

	/*if (b_exit_cond_g||last_element*/	
	if (compute_flag){

	mexPrintf("Xhat_0th_sum_start: pthread_self(): %d about to enter from compute_0th_sum() with m_index= %d, compute_flag=%d, threads_asleep_ctr=%d, last_element flag = %d, work_actually_computed_ctr=%d \n", (int) pthread_self(), m_index, (int)compute_flag, threads_asleep_ctr, (int)(*last_element_p), work_actually_computed_ctr);

	compute_0th_sum(argStruct_p->W_fom_p, argStruct_p->W_fom_cx_p, argStruct_p->Z_ol_p, m_index);

	sem_wait(&work_actually_computed_sem);
	work_actually_computed_ctr++;

	mexPrintf("Xhat_0th_sum_start: pthread_self(): %d returned from compute_0th_sum() with m_index= %d, compute_flag=%d, threads_asleep_ctr=%d, last_element flag = %d, work_actually_computed_ctr=%d \n", (int) pthread_self(), m_index, (int)compute_flag, threads_asleep_ctr, (int)(*last_element_p), work_actually_computed_ctr);

	if ((work_actually_computed_ctr==argStruct_p->work_queue_size)&&last_element_sleep_flag){

		/*Wake up the last element thread*/
		/*sem_post(&last_element_sleep_sem);*/
		work_actually_completed_flag=true; 

	}

	sem_post(&work_actually_computed_sem);

	compute_flag=false;

	}
	/*end if*/

	mexPrintf("Xhat_0th_sum_start: pthread_self(): %d with m_index= %d, compute_flag=%d, threads_asleep_ctr=%d TRYING TO ENTER SLEEP BLOCK!!!!!, last_element flag = %d, work_actually_computed_ctr=%d, last_element_found_global_flag=%d, work_actually_completed_flag=%d \n", (int) pthread_self(), m_index, (int)compute_flag, threads_asleep_ctr, (int)(*last_element_p), work_actually_computed_ctr, last_element_found_global_flag, work_actually_completed_flag);

	if ((!compute_flag)&&(last_element_found_global_flag)&&(work_actually_completed_flag)){

	sem_wait(&threads_asleep_sem);

	threads_asleep_ctr++;

	if (threads_asleep_ctr==NUM_WORKER_THREADS-1){

		sem_wait(&check_queue_global_sem);
		check_queue_global_flag=false;
		sem_post(&check_queue_global_sem);

		last_thread_to_sleep=true;

	}

	sem_post(&threads_asleep_sem);

	mexPrintf("Xhat_0th_sum_start: pthread_self(): %d with m_index= %d, compute_flag=%d, threads_asleep_ctr=%d ABOUT TO GO TO SLEEP!!!!!, last_element flag = %d, work_actually_computed_ctr=%d \n", (int) pthread_self(), m_index, (int)compute_flag, threads_asleep_ctr, (int)(*last_element_p), work_actually_computed_ctr);

	if(last_thread_to_sleep&&last_element_sleep_flag){
		sem_post(&last_element_sleep_sem);
	}

	/*sem_wait(&sleep_sem). aka sleep. when sem_post() is called from the outside world three times by the last_element==true thread, the three threads that were asleep/blocked, should awake and continue from this point. */
	sem_wait(&sleep_sem);

	mexPrintf("Xhat_0th_sum_start: pthread_self(): %d with m_index= %d, compute_flag=%d, threads_asleep_ctr=%d JUST WOKE UP FROM SLEEP!!!!!, last_element flag = %d, work_actually_computed_ctr \n", (int) pthread_self(), m_index, (int)compute_flag, threads_asleep_ctr, (int)(*last_element_p), work_actually_computed_ctr);

	sem_wait(&threads_asleep_sem);

	threads_asleep_ctr--;

	if (threads_asleep_ctr==0){

		sem_wait(&check_queue_global_sem);
		check_queue_global_flag=true;
		sem_post(&check_queue_global_sem);
	}

	sem_post(&threads_asleep_sem);

	}

	}

	/*If there's anything that would be better off being done, resetting, here as opposed to outside this function, by the last_element==true thread, can catch the thread with an if(){} and make it do something. */

} /*end while*/

}

static void Xhat_1st_sum_start(arg_struct_t* argStruct_p, bool* last_element_p){
/*Local variables*/
/*int l_index;*/
int m_index;
bool last_thread_to_sleep=false;
/*bool b_exit_cond_local=false; */

bool compute_flag; 
/*bool last_element=false;*/

while(threads_while_condition_1st_sum_flag&&(!(*last_element_p))) {		/*Change b_exit_cond_global to b_exit_cond_global_0th_sum?? */

	if (!(*last_element_p)) {	/*Assuming you make the outer while a function of last_element_p, then you can remove this. Can re-enable it if and when you're investigating bugs/debugging. */

	if ((!last_element_found_global_flag)&&(check_queue_global_flag)){

	/*Semaphore or mutex lock */
	sem_wait(&queue_sem);

	compute_flag=calculate_queue_index(argStruct_p, &m_index, last_element_p); 

	mexPrintf("Xhat_1st_sum_start: pthread_self(): %d returned from calculate_queue_index() with m_index= %d, compute_flag=%d, threads_asleep_ctr=%d, last_element flag = %d, work_actually_computed_ctr=%d \n", (int) pthread_self(), m_index, (int)compute_flag, threads_asleep_ctr, (int)(*last_element_p), work_actually_computed_ctr);

	/*Semaphore or mutex unlock */
	sem_post(&queue_sem);

	}

	/*if (b_exit_cond_g||last_element*/	
	if (compute_flag){

	mexPrintf("Xhat_1st_sum_start: pthread_self(): %d about to enter from compute_0th_sum() with m_index= %d, compute_flag=%d, threads_asleep_ctr=%d, last_element flag = %d, work_actually_computed_ctr=%d \n", (int) pthread_self(), m_index, (int)compute_flag, threads_asleep_ctr, (int)(*last_element_p), work_actually_computed_ctr);

	compute_1st_sum(argStruct_p->Y_lk_p, m_index);

	sem_wait(&work_actually_computed_sem);
	work_actually_computed_ctr++;

	mexPrintf("Xhat_1st_sum_start: pthread_self(): %d returned from compute_0th_sum() with m_index= %d, compute_flag=%d, threads_asleep_ctr=%d, last_element flag = %d, work_actually_computed_ctr=%d \n", (int) pthread_self(), m_index, (int)compute_flag, threads_asleep_ctr, (int)(*last_element_p), work_actually_computed_ctr);

	if ((work_actually_computed_ctr==argStruct_p->work_queue_size)&&last_element_sleep_flag){

		/*Wake up the last element thread*/
		/*sem_post(&last_element_sleep_sem);*/
		work_actually_completed_flag=true; 

	}

	sem_post(&work_actually_computed_sem);

	compute_flag=false;

	}
	/*end if*/

	mexPrintf("Xhat_1st_sum_start: pthread_self(): %d with m_index= %d, compute_flag=%d, threads_asleep_ctr=%d TRYING TO ENTER SLEEP BLOCK!!!!!, last_element flag = %d, work_actually_computed_ctr=%d, last_element_found_global_flag=%d, work_actually_completed_flag=%d \n", (int) pthread_self(), m_index, (int)compute_flag, threads_asleep_ctr, (int)(*last_element_p), work_actually_computed_ctr, last_element_found_global_flag, work_actually_completed_flag);

	if ((!compute_flag)&&(last_element_found_global_flag)&&(work_actually_completed_flag)){

	sem_wait(&threads_asleep_sem);

	threads_asleep_ctr++;

	if (threads_asleep_ctr==NUM_WORKER_THREADS-1){

		sem_wait(&check_queue_global_sem);
		check_queue_global_flag=false;
		sem_post(&check_queue_global_sem);

		last_thread_to_sleep=true;

	}

	sem_post(&threads_asleep_sem);

	mexPrintf("Xhat_1st_sum_start: pthread_self(): %d with m_index= %d, compute_flag=%d, threads_asleep_ctr=%d ABOUT TO GO TO SLEEP!!!!!, last_element flag = %d, work_actually_computed_ctr=%d \n", (int) pthread_self(), m_index, (int)compute_flag, threads_asleep_ctr, (int)(*last_element_p), work_actually_computed_ctr);

	if(last_thread_to_sleep&&last_element_sleep_flag){
		sem_post(&last_element_sleep_sem);
	}

	/*sem_wait(&sleep_sem). aka sleep. when sem_post() is called from the outside world three times by the last_element==true thread, the three threads that were asleep/blocked, should awake and continue from this point. */
	sem_wait(&sleep_sem);

	mexPrintf("Xhat_1st_sum_start: pthread_self(): %d with m_index= %d, compute_flag=%d, threads_asleep_ctr=%d JUST WOKE UP FROM SLEEP!!!!!, last_element flag = %d, work_actually_computed_ctr \n", (int) pthread_self(), m_index, (int)compute_flag, threads_asleep_ctr, (int)(*last_element_p), work_actually_computed_ctr);

	sem_wait(&threads_asleep_sem);

	threads_asleep_ctr--;

	if (threads_asleep_ctr==0){

		sem_wait(&check_queue_global_sem);
		check_queue_global_flag=true;
		sem_post(&check_queue_global_sem);
	}

	sem_post(&threads_asleep_sem);

	}

	}

	/*If there's anything that would be better off being done, resetting, here as opposed to outside this function, by the last_element==true thread, can catch the thread with an if(){} and make it do something. */

} /*end while*/
	
}

static void Xhat_2nd_sum_start(arg_struct_t* argStruct_p, bool* last_element_p){
	
/*Local variables*/
/*int l_index;*/
int f_index;	
int m_index;
bool last_thread_to_sleep=false;

/*bool b_exit_cond_local=false; */

bool compute_flag; 
/*bool last_element=false;*/

while(threads_while_condition_2nd_sum_flag&&(!(*last_element_p))) {		/*Change b_exit_cond_global to b_exit_cond_global_0th_sum?? */

	if (!(*last_element_p)) {	/*Assuming you make the outer while a function of last_element_p, then you can remove this. Can re-enable it if and when you're investigating bugs/debugging. */

	if ((!last_element_found_global_flag)&&(check_queue_global_flag)){

	/*Semaphore or mutex lock */
	sem_wait(&queue_sem);

	/*You need to create a different function to calculate a pair of queue indices. And therefore should also pass in f_index here to be updated by the function. */
	compute_flag=calculate_pair_queue_indices(argStruct_p, &m_index, &f_index, last_element_p); 

	mexPrintf("Xhat_2nd_sum_start: pthread_self(): %d returned from calculate_pair_queue_indices() with m_index= %d, f_index=%d, threads_asleep_ctr=%d, compute_flag=%d, last_element_found_global_flag=%d \n", (int) pthread_self(), m_index, f_index, threads_asleep_ctr, (int)compute_flag, (int)last_element_found_global_flag);

	/*Semaphore or mutex unlock. */
	sem_post(&queue_sem);

	}

	/*if (b_exit_cond_g||last_element*/	
	if (compute_flag){

	if (Xhat_pop_MKF_for_Phi_S_flag){

	compute_2nd_sum_pop_MKF(argStruct_p->Xhat_fnm_p, argStruct_p->Xtilde_fnm_p, argStruct_p->Xhat_low_fnm_p, argStruct_p->E_conj_fnm_p, argStruct_p->T_fk_p, argStruct_p->V_nk_p, argStruct_p->expj_Phi_S_nkf_p, m_index, f_index);

	} else {

	compute_2nd_sum(argStruct_p->Xhat_fnm_p, argStruct_p->Xtilde_fnm_p, argStruct_p->Xhat_low_fnm_p, argStruct_p->E_conj_fnm_p, argStruct_p->T_fk_p, argStruct_p->V_nk_p, argStruct_p->expj_Phi_S_nkf_p, m_index, f_index);

	}

	sem_wait(&work_actually_computed_sem);
	work_actually_computed_ctr++;

	if ((work_actually_computed_ctr==argStruct_p->work_queue_size)&&last_element_sleep_flag){

		/*Wake up the last element thread*/
		/*sem_post(&last_element_sleep_sem);*/
		work_actually_completed_flag=true;

	}

	sem_post(&work_actually_computed_sem);

	mexPrintf("Xhat_2nd_sum_start: pthread_self(): %d returned from compute_2nd_sum() with m_index= %d, f_index=%d, threads_asleep_ctr=%d, compute_flag=%d, last_element_found_global_flag=%d \n", (int) pthread_self(), m_index, f_index, threads_asleep_ctr, (int)compute_flag, (int) last_element_found_global_flag);	

	compute_flag=false;

	}
	/*end if*/

	if ((!compute_flag)&&(last_element_found_global_flag)&&(work_actually_completed_flag)) {

	sem_wait(&threads_asleep_sem);

	threads_asleep_ctr++;

	if (threads_asleep_ctr==NUM_WORKER_THREADS-1){

		sem_wait(&check_queue_global_sem);
		check_queue_global_flag=false;
		sem_post(&check_queue_global_sem);

		last_thread_to_sleep=true; 
	}

	sem_post(&threads_asleep_sem);

	if(last_thread_to_sleep&&last_element_sleep_flag){
		sem_post(&last_element_sleep_sem);
	}

	/*sem_wait(&sleep_sem). aka sleep. when sem_post() is called from the outside world three times by the last_element==true thread, the three threads that were asleep/blocked, should awake and continue from this point. */
	sem_wait(&sleep_sem);

	sem_wait(&threads_asleep_sem);

	threads_asleep_ctr--;

	if (threads_asleep_ctr==0){

		sem_wait(&check_queue_global_sem);
		check_queue_global_flag=true;
		sem_post(&check_queue_global_sem);
	}

	sem_post(&threads_asleep_sem);

	}

	}

	/*If there's anything that would be better off being done, resetting, here as opposed to outside this function, by the last_element==true thread, can catch the thread with an if(){} and make it do something. */

} /*end while*/
	
}

void* Xhat_start(void* arg){

thread_arg_t* threadArg_p;
arg_struct_t* argStruct_p;

bool last_element_flag=false; 
/*int thread_iter;*/

/*Cast the void* arg into a arg_struct_t**/
threadArg_p=(thread_arg_t*)arg;
argStruct_p=threadArg_p->argStruct_p; 
/*thread_iter=threadArg_p->thread_iter; */	/*Xhat_update doesn't require thread_iter, however, some other updates do. */

/*threads_while_condition_0th_sum_flag (a global flag) needs to be true upon any and all threads entering this point. */
Xhat_0th_sum_start(argStruct_p, &last_element_flag);
	
/*Allow the last thread to get to this point, and charge the last thread with the task of re-populating the queue. */	

if (last_element_flag){

/*sem_wait(&work_actually_computed_sem);*/

mexPrintf("After Xhat_0th_sum_start(), work_queue_size=%d !!!!!!!!!!!!!!!!!!!!!!!!!!! \n", argStruct_p->work_queue_size);

if (work_actually_computed_ctr<(argStruct_p->work_queue_size)){

	mexPrintf("After Xhat_0th_sum_start(), checking if work_actually_computed_ctr<argStruct_p->work_queue_size=%d!!!!!!!!!!!!\n", argStruct_p->work_queue_size);

	/*Sleep until the last thread to sleep gets to a point as close as possible to when it's about to sleep*/
	/*The intent is for the thread to wake you right before it goes to sleep.*/	
	last_element_sleep_flag=true; 

	sem_wait(&last_element_sleep_sem);

	/*Reset last_element_sleep_flag*/
	last_element_sleep_flag=false;	

}

/*sem_post(&work_actually_computed_sem);*/

mexPrintf("pthread_self():%d returned from Xhat_0th_sum_start() AND ABOUT TO REPOPULATE QUEUE FROM THE OUTSIDE WORLD!!!!!!!!!!!!!!!!!!!1\n", (int)pthread_self());

/*If queue isn't empty at this point*/
/*If you ever see weird behaviour where a certain update is consistently jumping into this code maybe set a breakpoint and investigate*/
/*if (argStruct_p->queue_intz.size()>0){
  clear_queue_intz(argStruct_p->queue_intz);
}*/

/*Repopulate queue*/
populate_queue_wrt_one_index(argStruct_p, M_static);

/*I think also reset: last_element_found_global_flag to false */
/*Consider creating a semaphore for this and wrapping it, but for now do without*/
last_element_found_global_flag=false;

/*Set the threads_while_condition_0th_sum_flag to a state that allows the just awoken threads to break their while loops. */
/*Wrap this write operation with a semaphore lock.*/
sem_wait(&threads_while_condition_0th_sum_sem);
threads_while_condition_0th_sum_flag=false;	
sem_post(&threads_while_condition_0th_sum_sem);

/*Wake the other threads by calling sem_post() three times. */
	sem_post(&sleep_sem);
	sem_post(&sleep_sem);
	sem_post(&sleep_sem);

	mexPrintf("After Xhat_0th_sum_start() SEM_POST() TIMES THREE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!");

/*Reset*/
work_actually_completed_flag=false;
work_actually_computed_ctr=0;	

} /*endif*/

/*reset last_element_flag to false for all threads*/
last_element_flag=false;

/*threads_while_condition_1st_sum_flag (a global flag) needs to be true upon any and all threads entering this point. */
Xhat_1st_sum_start(argStruct_p, &last_element_flag);

if (last_element_flag){

/*sem_wait(&work_actually_computed_sem);*/

if (work_actually_computed_ctr<(argStruct_p->work_queue_size)){

	/*Sleep until the last thread to sleep gets to a point as close as possible to when it's about to sleep*/
	/*The intent is for the thread to wake you right before it goes to sleep.*/	
	last_element_sleep_flag=true; 

	sem_wait(&last_element_sleep_sem);

	/*Reset last_element_sleep_flag*/
	last_element_sleep_flag=false;	

}

/*sem_post(&work_actually_computed_sem);*/

/*Repopulate queue*/
populate_queue_wrt_pair_indices(argStruct_p, M_static, F_static);

/*I think also reset: last_element_found_global_flag to false */
/*Consider creating a semaphore for this and wrapping it, but for now do without*/
last_element_found_global_flag=false;

/*Set the threads_while_condition_1st_sum_flag to a state that allows the just awoken threads to break their while loops. */
/*Wrap this write operation with a semaphore lock.*/
sem_wait(&threads_while_condition_1st_sum_sem);
threads_while_condition_1st_sum_flag=false;	
sem_post(&threads_while_condition_1st_sum_sem);

/*Wake the other threads by calling sem_post() three times. */
	sem_post(&sleep_sem);
	sem_post(&sleep_sem);
	sem_post(&sleep_sem);

/*Reset*/
work_actually_completed_flag=false;
work_actually_computed_ctr=0;		

} /*endif*/

/*reset last_element_flag to false for all threads*/
last_element_flag=false; 	

/*threads_while_condition_2nd_sum_flag (a global flag) needs to be true upon any and all threads entering this point. */
Xhat_2nd_sum_start(argStruct_p, &last_element_flag);

/*if(last_element): wake sleeping threads by only calling sem_post() three times, but no need to repopulate queue. */
if   (last_element_flag){

/*sem_wait(&work_actually_computed_sem);*/

if (work_actually_computed_ctr<(argStruct_p->work_queue_size)){

	/*Sleep until the last thread to sleep gets to a point as close as possible to when it's about to sleep*/
	/*The intent is for the thread to wake you right before it goes to sleep.*/	
	last_element_sleep_flag=true; 

	sem_wait(&last_element_sleep_sem);

	/*Reset last_element_sleep_flag*/
	last_element_sleep_flag=false;	

}

/*sem_post(&work_actually_computed_sem);*/

/*Don't populate the queue at this point! No updates follow. */

/*Since no updates follow, this probably isn't absolutely necessary*/
last_element_found_global_flag=false;

/*Set the threads_while_condition_1st_sum_flag to a state that allows the just awoken threads to break their while loops. */
/*Wrap this write operation with a semaphore lock.*/
sem_wait(&threads_while_condition_2nd_sum_sem);
threads_while_condition_2nd_sum_flag=false;	
sem_post(&threads_while_condition_2nd_sum_sem);

/*Wake the other threads by calling sem_post() three times. */
	sem_post(&sleep_sem);
	sem_post(&sleep_sem);
	sem_post(&sleep_sem);		

/*Reset*/
work_actually_completed_flag=false;
work_actually_computed_ctr=0;		

} /*endif*/

/*This probably isn't necessary either as threads will be destroyed and re-created but do it anyways*/
last_element_flag=false; 	

}