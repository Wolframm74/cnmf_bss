#include "local_inc.hpp"

static mat::fixed<O_static, L_static> Partial_wrt_Zol;

static cx_cube::fixed<F_static, L_static, M_static> FLM_tensor_local;
static cx_mat::fixed<F_static*L_static, M_static> FLxM_mat_local;
static cx_mat::fixed<M_static, F_static*L_static> MxFL_mat_local;
static cx_cube::fixed<M_static, F_static, L_static> MFL_tensor_quantity1_local;

static cx_cube::fixed<M_static, L_static, F_static> MLF_tensor_quantity2_local;

void hfk_Af_LS_estimator_Zol_update_API_serial_preprocess(mat* Y_lk_p, cx_cube* W_fom_cx_p){

int m_index, f_index;

/*populate quantity 1*/
for (m_index=0; m_index<M_static; m_index++){

FLM_tensor_local.slice(m_index)=h_fkm_shared_quantity.slice(m_index)*trans(*Y_lk_p);

}

/*rotate it*/

/*copy*/
arrayops::copy(FLxM_mat_local.memptr() , FLM_tensor_local.memptr(), F_static*L_static*M_static);

/*strans()*/
MxFL_mat_local=strans(FLxM_mat_local);

/*copy into target*/
arrayops::copy(MFL_tensor_quantity1_local.memptr(), MxFL_mat_local.memptr(), M_static*F_static*L_static);


/*populate quantity 2*/
for (f_index=0; f_index<F_static; f_index++){

MLF_tensor_quantity2_local.slice(f_index)=hfk_Af_LS_estimator_update_quantity_A_mkf.slice(f_index)*trans(*Y_lk_p);
	
}

rotate_W_fom_cx(W_fom_cx_p);

}

static void compute_Partial_hfk_Af_Cost_wrt_Zol(int o_index, int l_index, int thread_iter){

int f_index;

double output_val=0;

for (f_index=0; f_index<F_static; f_index++){

output_val=output_val+real(as_scalar(trans(MFL_tensor_quantity1_local.slice(l_index).col(f_index))*(W_mfo_cx_shared.slice(o_index).col(f_index))));

output_val=output_val+real(as_scalar(trans(W_mfo_cx_shared.slice(o_index).col(f_index))*(MFL_tensor_quantity1_local.slice(l_index).col(f_index))));

output_val=output_val-real(as_scalar(trans(MLF_tensor_quantity2_local.slice(f_index).col(l_index))*(W_mfo_cx_shared.slice(o_index).col(f_index))));

output_val=output_val-real(as_scalar(trans(W_mfo_cx_shared.slice(o_index).col(f_index))*(MLF_tensor_quantity2_local.slice(f_index).col(l_index))));

}	

Partial_wrt_Zol(o_index, l_index)=output_val;

}

static void compute_stuff_last_thread(mat* Z_ol_p){

double aeta=0.000001;

/*exponentiated gradient multiplicative update*/
(*Z_ol_p)=(*Z_ol_p)%exp(-aeta*Partial_wrt_Zol);

}


static void compute_Partial_hfk_Af_Cost_wrt_Zol_start(arg_struct_t* argStruct_p, int thread_iter){

int o_index, l_index;
bool last_element_queue_flag=false;		/*this may not be asbolutely necessary as before in the legacy code */

bool compute_flag;

while (threads_while_condition_0th_sum_flag&&(!(last_thread_to_sleep_flag_array[0][thread_iter]))){

	if ((!last_element_found_global_flag)&&(check_queue_global_flag)){

	/*Semaphore or mutex lock */
	sem_wait(&queue_sem);

	compute_flag=calculate_pair_queue_indices(argStruct_p, &o_index, &l_index, &last_element_queue_flag); 	/*Check that ourside world is populating queue wrt f=1:F and not n=1:N*/

	/*Semaphore or mutex unlock */
	sem_post(&queue_sem);

	}

	if (compute_flag){

	compute_Partial_hfk_Af_Cost_wrt_Zol(o_index, l_index, thread_iter);

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

/*static void Z_1st_sum_start(arg_struct_t* argStruct_p, int thread_iter){

}*/


void* compute_Partial_hfk_Af_Cost_wrt_Zol_thread_start(void* arg){

int thread_iter;
int i_iter;

thread_arg_t* threadArg_p;
arg_struct_t* argStruct_p;

threadArg_p=(thread_arg_t*)arg;
argStruct_p=threadArg_p->argStruct_p;
thread_iter=threadArg_p->thread_iter; 

compute_Partial_hfk_Af_Cost_wrt_Zol_start(argStruct_p, thread_iter);

if (last_thread_to_sleep_flag_array[0][thread_iter]){

compute_stuff_last_thread(argStruct_p->Z_ol_p);

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

if (!last_thread_to_sleep_flag_array[0][thread_iter]) {

sem_post(&checkpoint_sem);
/*sem_wait(&wait_till_queue_populated_sem);*/

}

threads_exit_and_signal(NUM_WORKER_THREADS);

}
