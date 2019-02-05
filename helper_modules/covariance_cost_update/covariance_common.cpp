#include "local_inc.hpp"

/*have to assign all these matrices and tensors their proper dimensions at some point */

static cx_mat::fixed<F_static*N_static, M_static> Xhat_fn_m;
static cx_mat::fixed<F_static*N_static, M_static> Xtilde_fn_m;	
static cx_mat::fixed<F_static*O_static, M_static> W_fo_m_cx;
static cx_mat::fixed<N_static*K_static, F_static> expj_Phi_S_nk_f;

static cx_mat::fixed<M_static, F_static*N_static> Xhat_m_fn;
static cx_mat::fixed<M_static, F_static*N_static> Xtilde_m_fn;	
static cx_mat::fixed<M_static, F_static*O_static> W_m_fo_cx;
static cx_mat::fixed<F_static, N_static*K_static> expj_Phi_S_f_nk;

/*xtilde_fnm, xhat_fnm, w_fom rotated */
static cx_cube::fixed<M_static, F_static, N_static> Xhat_mfn;
static cx_cube::fixed<M_static, F_static, N_static> Xtilde_mfn;	
static cx_cube::fixed<M_static, F_static, O_static> W_mfo_cx;
cx_cube::fixed<F_static, N_static, K_static> expj_Phi_S_fnk;

/*Just need to call this once ever, but don't forget to call it at least once*/
static void covariance_common_rotate_Xtilde(cx_cube* Xtilde_fnm_p){

	/*Copy into the MxFN destination*/
	arrayops::copy(Xtilde_fn_m.memptr(), (*Xtilde_fnm_p).memptr(), F_static*N_static*M_static );

	Xtilde_m_fn=strans(Xtilde_fn_m);

	arrayops::copy(Xtilde_mfn.memptr(), Xtilde_m_fn.memptr(), F_static*N_static*M_static );	

}

/*need to call this once every time in between updating any of the covariance related things. */
static void covariance_common_rotate_tensors(cx_cube* Xhat_fnm_p, cx_cube* W_fom_cx_p, cx_cube* expj_Phi_S_nkf_p){

	/*Copy into the destination*/
	arrayops::copy(Xhat_m_fn.memptr(), (*Xhat_fnm_p).memptr(), F_static*N_static*M_static );

	Xhat_m_fn=strans(Xhat_m_fn);

	arrayops::copy(Xhat_mfn.memptr(), Xhat_m_fn.memptr(), F_static*N_static*M_static );	

	/*Copy into the destination*/
	arrayops::copy(W_fo_m_cx.memptr(), (*W_fom_cx_p).memptr(), F_static*O_static*M_static );

	W_m_fo_cx=strans(W_fo_m_cx);

	arrayops::copy(W_mfo_cx.memptr(), W_m_fo_cx.memptr(), F_static*O_static*M_static );	

	/*do it again*/
	arrayops::copy(expj_Phi_S_nk_f.memptr(), (*expj_Phi_S_nkf_p).memptr(), F_static*N_static*K_static);

	expj_Phi_S_f_nk=strans(expj_Phi_S_nk_f);

	arrayops::copy(expj_Phi_S_fnk.memptr(), expj_Phi_S_f_nk.memptr(), F_static*N_static*K_static);

}

cube::fixed<F_static, N_static, O_static> tr_Efn_Efn_xw;
cube::fixed<F_static, N_static, O_static> tr_Efn_Efn_wx;

/*local variables but allocated globally instead of locally, all are a function of thread_iter */
static cx_mat::fixed<M_static, M_static> tempmat_local_1[NUM_WORKER_THREADS];
static cx_mat::fixed<M_static, M_static> tempmat_local_2[NUM_WORKER_THREADS];

/* int thread_iter*/
static void compute_tr_E_E_xw(int f_index, int n_index, int thread_iter){

int o_index;

/*Note: trans(cx_mat) corresponds to the conjugate transpose of the cx_mat*/
tempmat_local_1[thread_iter]=Xtilde_mfn.slice(n_index).col(f_index)*trans(Xtilde_mfn.slice(n_index).col(f_index));

/*Note: trans(cx_mat) corresponds to the conjugate transpose of the cx_mat*/
tempmat_local_1[thread_iter]=tempmat_local_1[thread_iter]+trans(tempmat_local_1[thread_iter]);

for (o_index=0; o_index<O_static; o_index++){

tempmat_local_2[thread_iter]=(tempmat_local_1[thread_iter])*(Xhat_mfn.slice(n_index).col(f_index))*trans(W_mfo_cx.slice(o_index).col(f_index));

tr_Efn_Efn_xw(f_index, n_index, o_index)=real(trace(tempmat_local_2[thread_iter]));

tempmat_local_2[thread_iter]=(tempmat_local_1[thread_iter])*W_mfo_cx.slice(o_index).col(f_index)*trans(Xhat_mfn.slice(n_index).col(f_index));

tr_Efn_Efn_wx(f_index, n_index, o_index)=real(trace(tempmat_local_2[thread_iter]));

}

}

/*static void compute_1st_sum(cx_cube* Xhat_low_fnm_p, cx_cube* E_conj_fnm_p, int m_index, int f_index){*/
/*static void compute_1st_sum(cube* Xhat_low_fnm_p, cx_cube* E_conj_fnm_p, int m_index, int f_index){	

}

static void compute_2nd_sum(cube* W_fom_p, cx_cube* W_fom_cx_p, int m_index, int l_index){

}

void projfunc_wrapper_Zol(mat* Z_ol_p, mxArray *Z_ol_mxArray_p){

}*/

/*With the exception of the W_fom update, the computation of the 3rd sum involves only M components, and therefore no need to set up a work queue for the computations for this*/
static void compute_tr_E_E_xw_last_thread(void){

}

static void compute_tr_E_E_xw_start(arg_struct_t* argStruct_p, int thread_iter){

int f_index, n_index;
bool last_element_queue_flag=false;		/*this may not be asbolutely necessary as before in the legacy code */

bool compute_flag;

while (threads_while_condition_0th_sum_flag&&(!(last_thread_to_sleep_flag_array[0][thread_iter]))){

	if ((!last_element_found_global_flag)&&(check_queue_global_flag)){

	/*Semaphore or mutex lock */
	sem_wait(&queue_sem);

	compute_flag=calculate_pair_queue_indices(argStruct_p, &f_index, &n_index, &last_element_queue_flag); 	/*Check that ourside world is populating queue wrt f=1:F and not n=1:N*/

	/*Semaphore or mutex unlock */
	sem_post(&queue_sem);

	}

	if (compute_flag){

	compute_tr_E_E_xw(f_index, n_index, thread_iter);

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

}

static void Z_2nd_sum_start(arg_struct_t* argStruct_p, int thread_iter){

}*/

void* compute_tr_E_E_xw_thread_start(void* arg){

int thread_iter;
int i_iter;

thread_arg_t* threadArg_p;
arg_struct_t* argStruct_p;

threadArg_p=(thread_arg_t*)arg;
argStruct_p=threadArg_p->argStruct_p;
thread_iter=threadArg_p->thread_iter; 

compute_tr_E_E_xw_start(argStruct_p, thread_iter);

if (last_thread_to_sleep_flag_array[0][thread_iter]){

compute_tr_E_E_xw_last_thread();

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
