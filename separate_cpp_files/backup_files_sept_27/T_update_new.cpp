#include "local_inc.hpp"

/*2nd sum stuff*/
cube::fixed<F_static, K_static, M_static> outtensor_2nd_real_FKM;
cx_cube::fixed<F_static, K_static, M_static> outtensor_2nd_cx_FKM;

/*3rd sum stuff*/
mat::fixed<F_static, K_static> outmat_3rd_real_FK_den;
mat::fixed<F_static, K_static> outmat_3rd_real_FK_num;
cx_mat::fixed<F_static, K_static> outmat_3rd_cx_FK;

static void compute_2nd_sum(cube* Xhat_low_fnm_p, cx_cube* E_conj_fnm_p, mat* V_nk_p, cx_cube* expj_Phi_S_nkf_p, int m_index, int f_index){

#ifdef DEBUG	
mexPrintf("T_update: pthread_self(): %d, compute_2nd_sum(), m_index=%d, f_index=%d \n", (int)pthread_self(), m_index, f_index);		
#endif	

/*1xK=(1xN)*(NxK) */
outtensor_2nd_real_FKM.slice(m_index).row(f_index)=(*Xhat_low_fnm_p).slice(m_index).row(f_index)*(*V_nk_p);

/*If you notice that the T_update is slow, you can try allocating a fixed NxK dummy cx_mat to store the intermediate result of the elementwise multiplication */
/*expj_Phi_S_fnk needs to be implemented as two tensors:
one that is NxKxF (for this particular update).
+maybe a few others ways to accomodate other updates. 
*/
outtensor_2nd_cx_FKM.slice(m_index).row(f_index)=(*E_conj_fnm_p).slice(m_index).row(f_index)*((*expj_Phi_S_nkf_p).slice(f_index)%(*V_nk_p));

}

static void compute_3rd_sum_last_thread(mat* T_fk_p, mat* V_nk_p){

outmat_3rd_real_FK_den.zeros();
outmat_3rd_cx_FK.zeros();

/*Integrate out m*/
for (int m_iter=0; m_iter<M_static; m_iter++){	
outmat_3rd_real_FK_den=outmat_3rd_real_FK_den+Xhat_outtensor_real_FKM.slice(m_iter)%outtensor_2nd_real_FKM.slice(m_iter);

outmat_3rd_cx_FK=outmat_3rd_cx_FK+Xhat_outtensor_cx_FKM.slice(m_iter)%outtensor_2nd_cx_FKM.slice(m_iter);
}

outmat_3rd_real_FK_num=outmat_3rd_real_FK_den+real(outmat_3rd_cx_FK);

(*T_fk_p)=(*T_fk_p)%(outmat_3rd_real_FK_num/outmat_3rd_real_FK_den);

/*Set negative elements to zero*/
(*T_fk_p).elem(find((*T_fk_p)<0)).zeros();

normalize_T_fk(T_fk_p, V_nk_p);

}

static void T_2nd_sum_start(arg_struct_t* argStruct_p, bool* last_element_p){

int m_index;
int f_index;
bool last_thread_to_sleep=false; 

bool compute_flag; 
/*bool last_element=false;*/

while(threads_while_condition_2nd_sum_flag&&(!(*last_element_p))) {		/*Change b_exit_cond_global to b_exit_cond_global_0th_sum?? */

	if (!(*last_element_p)) {	/*Assuming you make the outer while a function of last_element_p, then you can remove this. Can re-enable it if and when you're investigating bugs/debugging. */

	if ((!last_element_found_global_flag)&&(check_queue_global_flag)){

	/*Semaphore or mutex lock */
	sem_wait(&queue_sem);

	/*You need to create a different function to calculate a pair of queue indices. And therefore should also pass in f_index here to be updated by the function. */
	compute_flag=calculate_pair_queue_indices(argStruct_p, &m_index, &f_index, last_element_p); 

	/*Semaphore or mutex unlock.*/
	sem_post(&queue_sem);

	}

	/*if (b_exit_cond_g||last_element*/	
	if (compute_flag){

	compute_2nd_sum(argStruct_p->Xhat_low_fnm_p, argStruct_p->E_conj_fnm_p, argStruct_p->V_nk_p, argStruct_p->expj_Phi_S_nkf_p, m_index, f_index);

	sem_wait(&work_actually_computed_sem);
	work_actually_computed_ctr++;

	if (work_actually_computed_ctr==argStruct_p->work_queue_size){

		/*Wake up the last element thread*/
		/*sem_post(&last_element_sleep_sem);*/
		work_actually_completed_flag=true; 

	}

	sem_post(&work_actually_computed_sem);

	compute_flag=false;

	}
	/*end if*/

	/*Fix: I think this absolutely needs to be dependent on !(*last_element_p). ie: the thread that finds the last element MUSTN'T EVEN ENTER THIS CODE BLOCK
	AND EXPOSE ITSELF TO THE POSSIBILITY OF GOING TO SLEEP.....*/
	if ((!compute_flag)&&(last_element_sleep_flag)&&(work_actually_completed_flag)&&(!(*last_element_p))){

	sem_wait(&threads_asleep_sem);

	threads_asleep_ctr++;

	if (threads_asleep_ctr==NUM_WORKER_THREADS-1){

		sem_wait(&check_queue_global_sem);
		check_queue_global_flag=false;
		sem_post(&check_queue_global_sem);

		last_thread_to_sleep=true;

	}

	sem_post(&threads_asleep_sem);

	if(last_thread_to_sleep){
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

void* T_start(void* arg){

thread_arg_t* threadArg_p;
arg_struct_t* argStruct_p;

bool last_element_flag=false; 
/*int thread_iter;*/

/*Cast the void* arg into a arg_struct_t**/
threadArg_p=(thread_arg_t*)arg;
argStruct_p=threadArg_p->argStruct_p; 
/*thread_iter=threadArg_p->thread_iter; */

T_2nd_sum_start(argStruct_p, &last_element_flag);

if (last_element_flag){

if (threads_asleep_ctr<(NUM_WORKER_THREADS-1)){

	last_element_sleep_flag=true;

	sem_wait(&last_element_sleep_sem);

	last_element_sleep_flag=false;

}	

compute_3rd_sum_last_thread(argStruct_p->T_fk_p, argStruct_p->V_nk_p);

last_element_found_global_flag=false; 

sem_wait(&threads_while_condition_2nd_sum_sem);
threads_while_condition_2nd_sum_flag=false;
sem_post(&threads_while_condition_2nd_sum_sem);

/*Wake up other threads*/
	sem_post(&sleep_sem);
	sem_post(&sleep_sem);
	sem_post(&sleep_sem);

/*Reset*/	
work_actually_completed_flag=false;
work_actually_computed_ctr=0;	

}

last_element_flag=false;

threads_exit_and_signal(NUM_WORKER_THREADS);
		
}