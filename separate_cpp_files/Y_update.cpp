#include "local_inc.hpp"

/*1st sum stuff*/
cube::fixed<F_static, K_static, M_static> outtensor_1st_real_FKM;
cx_cube::fixed<F_static, K_static, M_static> outtensor_1st_cx_FKM;

/*2nd sum stuff*/
cube::fixed<L_static, K_static, M_static> outtensor_2nd_real_LKM;
cx_cube::fixed<L_static, K_static, M_static> outtensor_2nd_cx_LKM;

/*3rd sum stuff*/
mat::fixed<L_static, K_static> outmat_3rd_real_LK_den;
mat::fixed<L_static, K_static> outmat_3rd_real_LK_num;
cx_mat::fixed<L_static, K_static> outmat_3rd_cx_LK;

static void compute_1st_sum(cube* Xhat_low_fnm_p, cx_cube* E_conj_fnm_p, mat* V_nk_p, cx_cube* expj_Phi_S_nkf_p, int m_index, int f_index){

#ifdef DEBUG	
mexPrintf("Y_update: pthread_self(): %d, compute_1st_sum(), m_index=%d, f_index=%d \n", (int)pthread_self(), m_index, f_index);		
#endif

/*1xK=(1xN)*(NxK) */
outtensor_1st_real_FKM.slice(m_index).row(f_index)=(*Xhat_low_fnm_p).slice(m_index).row(f_index)*(*V_nk_p);

/*If you notice that the T_update is slow, you can try allocating a fixed NxK dummy cx_mat to store the intermediate result of the elementwise multiplication */
outtensor_1st_cx_FKM.slice(m_index).row(f_index)=(*E_conj_fnm_p).slice(m_index).row(f_index)*((*expj_Phi_S_nkf_p).slice(f_index)%(*V_nk_p));

}

/*This computation should involve Xhat_outtensor_cx_FLM and Xhat_outtensor_real_FLM*/
/*Dimensions: 1xL=(1xF)*(FxL) */
static void compute_2nd_sum(mat* T_fk_p, int m_index, int k_index){

#ifdef DEBUG	
mexPrintf("Y_update: pthread_self(): %d, compute_2nd_sum(), m_index=%d, k_index=%d \n", (int)pthread_self(), m_index, k_index);		
#endif	

/*outtensor real*/	/*Quick glance: shift M over to the right for every subcube() */	/*this is a row by matrix multiplication if i'm not mistaken, try compiling and see what it says. */	/*If not you'll have to think of a clever implementation*/
/*outtensor_2nd_real_LKM.subcube(k_index,m_index,0,k_index,m_index,L_static-1)=(outtensor_1st_real_FKM.subcube(0,m_index,k_index,F_static-1,m_index,k_index)%(*T_fk_p).submat(0,k_index,F_static-1,k_index))*Xhat_outtensor_real_FLM.subcube(0,m_index,0,F_static-1,m_index,L_static-1);*/
outtensor_2nd_real_LKM.slice(m_index).col(k_index)=trans(trans(outtensor_1st_real_FKM.slice(m_index).col(k_index)%(*T_fk_p).col(k_index))*Xhat_outtensor_real_FLM.slice(m_index));

/*outtensor cx*/	
/*outtensor_2nd_cx_LKM.subcube(k_index,m_index,0,k_index,m_index,L_static-1)=(outtensor_1st_cx_FKM.subcube(0,m_index,k_index,F_static-1,m_index,k_index)%(*T_fk_p).submat(0,k_index,F_static-1,k_index))*Xhat_outtensor_cx_FLM.subcube(0,m_index,0,F_static-1,m_index,L_static-1);*/
outtensor_2nd_cx_LKM.slice(m_index).col(k_index)=strans(strans(outtensor_1st_cx_FKM.slice(m_index).col(k_index)%(*T_fk_p).col(k_index))*Xhat_outtensor_cx_FLM.slice(m_index));

}

void projfunc_wrapper_Ylk(mat* Y_lk_p, mxArray *Y_lk_mxArray_p){

mxArray* output_arg_p[1];

armaSetPr(Y_lk_mxArray_p, *Y_lk_p);

/*mexCallMATLAB(1, &output_arg_p, 1, &plhs[2], "projfunc_wrapper_matlab");*/
/*mexCallMATLAB(0, NULL, 1, &Z_ol_mxArray_p, "projfunc_wrapper_matlab");*/
mexCallMATLAB(1, &output_arg_p[0], 1, &Y_lk_mxArray_p, "projfunc_wrapper_matlab_Ylk");

/*memcpy(mxGetPr(Z_ol_mxArray_p), mxGetPr(output_arg_p), O_static*L_static*sizeof(double));*/

/*arrayops::copy((*Z_ol_p).memptr(), mxGetPr(output_arg_p[0]), O_static*L_static);*/

(*Y_lk_p).set_real(armaGetPr(output_arg_p[0], true)); 	

}

void orthogonalize_Ylk_square_wrapper(mat* Y_lk_p, mxArray *Y_lk_mxArray_p){

mxArray* output_arg_p[1];

armaSetPr(Y_lk_mxArray_p, *Y_lk_p);

/*mexCallMATLAB(1, &output_arg_p, 1, &plhs[2], "projfunc_wrapper_matlab");*/
/*mexCallMATLAB(0, NULL, 1, &Z_ol_mxArray_p, "projfunc_wrapper_matlab");*/
mexCallMATLAB(1, &output_arg_p[0], 1, &Y_lk_mxArray_p, "orthogonalize_Y_lk_square");

/*memcpy(mxGetPr(Z_ol_mxArray_p), mxGetPr(output_arg_p), O_static*L_static*sizeof(double));*/

/*arrayops::copy((*Z_ol_p).memptr(), mxGetPr(output_arg_p[0]), O_static*L_static);*/

(*Y_lk_p).set_real(armaGetPr(output_arg_p[0], true)); 	

}

static void compute_3rd_sum_last_thread(mat* Y_lk_p, mat* compensation_mat_p, mat* Z_ol_p){

outmat_3rd_real_LK_den.zeros();
outmat_3rd_cx_LK.zeros();

int m_iter;

/*Integrate out m*/
for (m_iter=0; m_iter<M_static; m_iter++){

outmat_3rd_real_LK_den=outmat_3rd_real_LK_den+outtensor_2nd_real_LKM.slice(m_iter);

outmat_3rd_cx_LK=outmat_3rd_cx_LK+outtensor_2nd_cx_LKM.slice(m_iter);

}

outmat_3rd_real_LK_num=outmat_3rd_real_LK_den+real(outmat_3rd_cx_LK);

outmat_3rd_real_LK_num.print("outmat_3rd_real_LK_num:");
outmat_3rd_real_LK_den.print("outmat_3rd_real_LK_den:");

/*Normalize by l1-norms*/
/*outmat_3rd_real_LK_num=outmat_3rd_real_LK_num/accu(outmat_3rd_real_LK_num);

outmat_3rd_real_LK_den=outmat_3rd_real_LK_den/accu(outmat_3rd_real_LK_den);*/

/*Zok_Ylk_update_extra_sparsity(Z_ol_p, Y_lk_p);*/

/*Zok_sparsity_cost_update stuff*/
/*compute_Zok_sparsity_common_Quantities(Z_ol_p, Y_lk_p);

Zok_sparsity_Ylk_update(Z_ol_p);*/

/*(*Y_lk_p)=(*Y_lk_p)%((outmat_3rd_real_LK_num+nve_part_Partial_wrt_Ylk)/(outmat_3rd_real_LK_den+pve_part_Partial_wrt_Ylk)); */

/*POSTPONE THIS~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
/*(*Y_lk_p)=(*Y_lk_p)%((outmat_3rd_real_LK_num)/(outmat_3rd_real_LK_den)); */

//+Zok_sparsity_outmat_LK_den));

/*Set nonnegative elements to zero*/
/*(*Y_lk_p).elem(find((*Y_lk_p)<0)).zeros();*/
(*Y_lk_p).elem(find((*Y_lk_p)<=0)).fill(0.00000001);
/*(*Y_lk_p).elem(find((*Y_lk_p)<0)).fill(0.0000000000001);*/

/*normalize_Y_local(Y_lk_p);*/
/*normalize_Y_lk(Y_lk_p, T_fk_p);*/
/*normalize_Y_lk(Y_lk_p, compensation_mat_p);*/

/*Populate this for Z update: move this to Y_lk_post_update.!!*/
/*Z_inmat_real_Ykl=trans(*Y_lk_p);*/

}

static void Y_1st_sum_start(arg_struct_t* argStruct_p, int thread_iter){

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

	compute_1st_sum(argStruct_p->Xhat_low_fnm_p, argStruct_p->E_conj_fnm_p, argStruct_p->V_nk_p, argStruct_p->expj_Phi_S_nkf_p, m_index, f_index);

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

static void Y_2nd_sum_start(arg_struct_t* argStruct_p, int thread_iter){

int m_index;
int k_index;

bool last_element_queue_flag=false;

bool compute_flag; 

while (threads_while_condition_2nd_sum_flag&&(!(last_thread_to_sleep_flag_array[2][thread_iter]))){

	if ((!last_element_found_global_flag)&&(check_queue_global_flag)){

	/*Semaphore or mutex lock */
	sem_wait(&queue_sem);

	compute_flag=calculate_pair_queue_indices(argStruct_p, &m_index, &k_index, &last_element_queue_flag); 

	/*Semaphore or mutex unlock */
	sem_post(&queue_sem);

	}

	if (compute_flag){

	compute_2nd_sum(argStruct_p->T_fk_p, m_index, k_index);

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

void* Y_primary_auxfun_start(void* arg){

int i_iter;
int thread_iter;

thread_arg_t* threadArg_p;
arg_struct_t* argStruct_p;

threadArg_p=(thread_arg_t*)arg;
argStruct_p=threadArg_p->argStruct_p;
thread_iter=threadArg_p->thread_iter; 

Y_1st_sum_start(argStruct_p, thread_iter);

if (last_thread_to_sleep_flag_array[1][thread_iter]){

sem_wait(&threads_while_condition_1st_sum_sem);
threads_while_condition_1st_sum_flag=false;	
sem_post(&threads_while_condition_1st_sum_sem);

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

if (!last_thread_to_sleep_flag_array[1][thread_iter]){

sem_post(&checkpoint_sem);
sem_wait(&wait_till_queue_populated_sem);

}

if (last_thread_to_sleep_flag_array[1][thread_iter]){

populate_queue_wrt_pair_indices(argStruct_p, M_static, K_static);

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
/*sem_post(&wait_till_queue_populated_sem);
sem_post(&wait_till_queue_populated_sem);
sem_post(&wait_till_queue_populated_sem);*/

for (i_iter=0; i_iter<(NUM_WORKER_THREADS-1); i_iter++){
	sem_post(&wait_till_queue_populated_sem);
}

}

Y_2nd_sum_start(argStruct_p, thread_iter);
	
if (last_thread_to_sleep_flag_array[2][thread_iter]){

compute_3rd_sum_last_thread(argStruct_p->Y_lk_p, argStruct_p->V_nk_p, argStruct_p->Z_ol_p);

sem_wait(&threads_while_condition_2nd_sum_sem);
threads_while_condition_2nd_sum_flag=false;	
sem_post(&threads_while_condition_2nd_sum_sem);

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

if (!last_thread_to_sleep_flag_array[2][thread_iter]){

sem_post(&checkpoint_sem);

}

threads_exit_and_signal(NUM_WORKER_THREADS);

}

