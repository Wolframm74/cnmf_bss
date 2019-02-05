#include "local_inc.hpp"

/*2nd sum stuff*/
cube::fixed<N_static, K_static, M_static> outtensor_2nd_real_NKM;
cx_cube::fixed<N_static, K_static, M_static> outtensor_2nd_cx_NKM;

mat::fixed<F_static, K_static> dummy_mat_real_FK;
cx_mat::fixed<F_static, K_static> dummy_mat_cx_FK;

/*3rd sum stuff*/
mat::fixed<N_static, K_static> outmat_3rd_real_NK_den;
mat::fixed<N_static, K_static> outmat_3rd_real_NK_num;
cx_mat::fixed<N_static, K_static> outmat_3rd_cx_NK;

static void compute_2nd_sum(cube* Xhat_low_fnm_p, cx_cube* E_conj_fnm_p, mat* T_fk_p, cx_cube* expj_Phi_S_fkn_p, int m_index, int n_index){

#ifdef DEBUG	
mexPrintf("V_update: pthread_self(): %d, compute_2nd_sum(), m_index=%d, n_index=%d \n", (int)pthread_self(), m_index, n_index);		
#endif	

/*populate intermediate matrices*/
dummy_mat_real_FK=Xhat_outtensor_real_FKM.slice(m_index)%(*T_fk_p);

dummy_mat_cx_FK=Xhat_outtensor_cx_FKM.slice(m_index)%(*T_fk_p)%(*expj_Phi_S_fkn_p).slice(n_index);

/*outtensor_2nd_real_NKM: 1xK=1xF*FxK*/
outtensor_2nd_real_NKM.slice(m_index).row(n_index)=trans((*Xhat_low_fnm_p).slice(m_index).col(n_index))*dummy_mat_real_FK;

/*outtensor_2nd_cx_NKM*/	
outtensor_2nd_cx_NKM.slice(m_index).row(n_index)=strans((*E_conj_fnm_p).slice(m_index).col(n_index))*dummy_mat_cx_FK;

}

static void compute_3rd_sum_last_thread(mat* V_nk_p){

/*The function calculates the nonnegative ratios. Should be necessary to zero the dummy/memory matrices at least once as follows: */
/*Only on thread should end up calling this function. Not a bad thing to zero in case you didn't zero at least once already. */	
outmat_3rd_real_NK_den.zeros();
outmat_3rd_cx_NK.zeros();

/*Integrate out m*/
for (int m_iter=0; m_iter<M_static; m_iter++){	
outmat_3rd_real_NK_den=outmat_3rd_real_NK_den+outtensor_2nd_real_NKM.slice(m_iter);

outmat_3rd_cx_NK=outmat_3rd_cx_NK+outtensor_2nd_cx_NKM.slice(m_iter);

}

outmat_3rd_real_NK_num=outmat_3rd_real_NK_den+real(outmat_3rd_cx_NK);

(*V_nk_p)=(*V_nk_p)%(outmat_3rd_real_NK_num/outmat_3rd_real_NK_den);

/*Project V_nk onto the nonnegative orthant. ie: set all negative elements to zero. */
(*V_nk_p).elem(find((*V_nk_p)<0)).zeros();

}

static void V_2nd_sum_start(arg_struct_t* argStruct_p, bool* last_element_p){

int m_index;
int n_index;
bool last_thread_to_sleep=false; 

bool compute_flag; 
/*bool last_element=false;*/

while(threads_while_condition_2nd_sum_flag&&(!(*last_element_p))) {		/*Change b_exit_cond_global to b_exit_cond_global_0th_sum?? */

	if (!(*last_element_p)) {	/*Assuming you make the outer while a function of last_element_p, then you can remove this. Can re-enable it if and when you're investigating bugs/debugging. */

	if ((!last_element_found_global_flag)&&(check_queue_global_flag)){

	/*Semaphore or mutex lock */
	sem_wait(&queue_sem);

	/*You need to create a different function to calculate a pair of queue indices. And therefore should also pass in f_index here to be updated by the function. */
	compute_flag=calculate_pair_queue_indices(argStruct_p, &m_index, &n_index, last_element_p); 

	/*Semaphore or mutex unlock.*/
	sem_post(&queue_sem);

	}

	/*if (b_exit_cond_g||last_element*/	
	if (compute_flag){

	compute_2nd_sum(argStruct_p->Xhat_low_fnm_p, argStruct_p->E_conj_fnm_p, argStruct_p->T_fk_p, argStruct_p->expj_Phi_S_fkn_p, m_index, n_index);

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


void* V_start(void* arg){

thread_arg_t* threadArg_p;
arg_struct_t* argStruct_p;

bool last_element_flag=false; 
/*int thread_iter;*/

/*Cast the void* arg into a arg_struct_t**/
threadArg_p=(thread_arg_t*)arg;
argStruct_p=threadArg_p->argStruct_p; 
/*thread_iter=threadArg_p->thread_iter; */

V_2nd_sum_start(argStruct_p, &last_element_flag);

if (last_element_flag){

if (threads_asleep_ctr<(NUM_WORKER_THREADS-1)){

	last_element_sleep_flag=true;

	sem_wait(&last_element_sleep_sem);

	last_element_sleep_flag=false;

}	

compute_3rd_sum_last_thread(argStruct_p->V_nk_p);

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