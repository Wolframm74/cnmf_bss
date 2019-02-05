#include "local_inc.hpp"

colvec::fixed<N_static> ones_col_Nx1; /*This needs to be set to .ones() just once. ever, before either Z update or W update */

/*0th sum corresponds to basic matrix multiplication*/
mat::fixed<K_static, O_static> W_fom_outmat_0th_real_KO; 

/*1st sum stuff*/
cube::fixed<N_static, O_static, F_static> outtensor_1st_real_NOF;
cx_cube::fixed<N_static, O_static, F_static> outtensor_1st_cx_NOF;
static cx_mat::fixed<N_static, K_static> dummy_mat_cx_NK[NUM_WORKER_THREADS];
static mat::fixed<K_static, K_static> dummy_mat_outerprod_KK[NUM_WORKER_THREADS];
static mat::fixed<N_static, K_static> dummy_mat_real_NK[NUM_WORKER_THREADS];

/*2nd sum stuff*/
cube::fixed<F_static, O_static, M_static> outtensor_2nd_real_FOM;
cx_cube::fixed<F_static, O_static, M_static> outtensor_2nd_cx_FOM;

/*Numerator tensor*/
cube::fixed<F_static, O_static, M_static> numerator_tensor_real_FOM;

/*cube::fixed<M_static, O_static, F_static> W_outtensor_2nd_real_MOF;
cx_cube::fixed<M_static, O_static, F_static> W_outtensor_2nd_cx_MOF;*/

void W_0th_sum_preprocess(mat* Z_ol_p, mat* Y_lk_p){

W_fom_outmat_0th_real_KO=trans((*Z_ol_p)*(*Y_lk_p));

}

static void compute_1st_sum(mat* T_fk_p, mat* V_nk_p, cx_cube* expj_Phi_S_nkf_p, int f_index, int thread_iter){

#ifdef DEBUG	
mexPrintf("W_update: pthread_self(): %d, compute_1st_sum(), f_index=%d, work_actually_computed_ctr=%d, check_queue_global_flag=%d\n", (int)pthread_self(), f_index, work_actually_computed_ctr, check_queue_global_flag);		
#endif

/*dummy_mat_outerprod_KK[thread_iter]=trans((*V_nk_p).submat(n_index, 0, n_index, K_static-1))*((*V_nk_p).submat(n_index,0,n_index,K_static-1));*/
dummy_mat_outerprod_KK[thread_iter]=trans((*T_fk_p).row(f_index))*((*T_fk_p).row(f_index));

/*outtensor_1st_real_NOF*/
/*outtensor_1st_real_NOF.subcube(0,n_index,0,F_static-1,n_index,O_static-1)=(*T_fk_p)*dummy_mat_outerprod_KK[thread_iter]*W_fom_outmat_0th_real_KO;*/
outtensor_1st_real_NOF.slice(f_index)=(*V_nk_p)*dummy_mat_outerprod_KK[thread_iter]*W_fom_outmat_0th_real_KO;

/*Prepping dummy_mat_real_NK*/
dummy_mat_real_NK[thread_iter]=kron(ones_col_Nx1, (*T_fk_p).row(f_index));

/*Prepping dummy_mat_cx_FK*/
/*dummy_mat_cx_FK[thread_iter]=((*expj_Phi_S_fnk_p).subcube(0,n_index,0,F_static-1,n_index,K_static-1))%(*T_fk_p);*/
dummy_mat_cx_NK[thread_iter]=((*expj_Phi_S_nkf_p).slice(f_index))%(*V_nk_p);

/*dummy_mat_cx_FK[thread_iter]=dummy_mat_cx_FK[thread_iter].each_row()%(*V_nk_p).submat(n_index,0,n_index,K_static-1);*/
/*dummy_mat_cx_NK[thread_iter]=(dummy_mat_cx_NK[thread_iter].each_row())%((*T_fk_p).row(f_index));*/
dummy_mat_cx_NK[thread_iter]=(dummy_mat_cx_NK[thread_iter])%(dummy_mat_real_NK[thread_iter]);

/*outtensor_1st_cx_NOF*/
/*outtensor_1st_cx_NOF.subcube(0,n_index,0,F_static-1,n_index,O_static-1)=dummy_mat_cx_FK[thread_iter]*W_fom_outmat_0th_real_KO;*/
outtensor_1st_cx_NOF.slice(f_index)=dummy_mat_cx_NK[thread_iter]*W_fom_outmat_0th_real_KO;

}

static void compute_2nd_sum(cube* Xhat_low_fnm_p, cx_cube* E_conj_fnm_p, int m_index, int f_index){

#ifdef DEBUG	
mexPrintf("W_update: pthread_self(): %d, compute_2nd_sum(), m_index=%d, f_index=%d, work_actually_computed_ctr=%d, check_queue_global_flag=%d \n", (int)pthread_self(), m_index, f_index, work_actually_computed_ctr, check_queue_global_flag);		
#endif

/*outtensor_1st_real_FML computation*/

/*Integrate out n. 1xL=1xN*NxO */
/*outtensor_2nd_real_FOM.subcube(f_index, m_index, 0, f_index, m_index, O_static-1)=(*Xhat_low_fnm_p).subcube(f_index, 0, m_index, f_index, N_static-1, m_index)*outtensor_1st_real_NOF.submat(f_index,0,0,f_index,N_static-1,O_static-1);	*/
/*outtensor_2nd_real_FOM.slice(m_index).row(f_index)=(Intensor_Xhat_low_mnf.slice(f_index).row(m_index))*(outtensor_1st_real_NOF.slice(f_index));*/
outtensor_2nd_real_FOM.slice(m_index).row(f_index)=((*Xhat_low_fnm_p).slice(m_index).row(f_index))*(outtensor_1st_real_NOF.slice(f_index));

/*outtensor_2nd_real_MOF.slice(f_index).row(m_index)=outtensor_2nd_real_FOM.slice(m_index).row(f_index);*/

/*outtensor_1st_cx_FML computation*/	/*No clue/not clear what i'm doing dimensionswise. could be missing a transpose. */
/*outtensor_2nd_cx_FOM.subcube(f_index, m_index, 0, f_index, m_index, O_static-1)=(*E_conj_fnm_p).subcube(f_index, 0, m_index, f_index, N_static-1, m_index)*outtensor_1st_cx_NOF.submat(f_index,0,0,f_index,N_static-1,O_static-1);	*/
/*outtensor_2nd_cx_FOM.slice(m_index).row(f_index)=(Intensor_E_conj_mnf.slice(f_index).row(m_index))*(outtensor_1st_cx_NOF.slice(f_index));*/
outtensor_2nd_cx_FOM.slice(m_index).row(f_index)=((*E_conj_fnm_p).slice(m_index).row(f_index))*(outtensor_1st_cx_NOF.slice(f_index));

/*outtensor_2nd_cx_MOF.slice(f_index).row(m_index)=outtensor_2nd_cx_FOM.slice(m_index).row(f_index);*/

}

static void compute_stuff_last_thread(cube* W_fom_p, cx_cube* expj_Phi_W_fom_p){

/*outtensor_2nd_cx_FOM.print("outtensor_2nd_cx_FOM:");
(*expj_Phi_W_fom_p).print("expj_Phi_W_fom:");*/

numerator_tensor_real_FOM=outtensor_2nd_real_FOM+real(outtensor_2nd_cx_FOM%(*expj_Phi_W_fom_p));

/*numerator_tensor_real_FOM.print("numerator_tensor_real_FOM:");
outtensor_2nd_real_FOM.print("outtensor_2nd_real_FOM:");*/

/*(*W_fom_p)=(*W_fom_p)%(numerator_tensor_real_FOM/outtensor_2nd_real_FOM);*/

}

/*real double matrix for which each f,o'th element is the frobenius norm for that Mx1 column vector*/
mat::fixed<F_static, O_static> FrobNormMat;
mat::fixed<F_static, O_static> ones_mat_FxO;		/*Set this to ones in the algorithm mex entry point*/

/*Compute the Frobenius norm across columns of 1:M for all f=1:F and o=1:O */
static void normalize_W_local(cube* W_fom_p){

int m;

FrobNormMat.zeros();

/*for 1:M; accumulate the squares*/ 
for (m=0; m<M_static; m++){

	FrobNormMat=FrobNormMat+((*W_fom_p).slice(m)%(*W_fom_p).slice(m));

}

/*sqrt(of an FxO matrix)*/
FrobNormMat=sqrt(FrobNormMat);

/*for 1:M; elementwise multiply each slice of W_fom with ones(F_static, O_static)./FrobNormMat */	
for (m=0; m<M_static; m++){

	(*W_fom_p).slice(m)=((*W_fom_p).slice(m))%(ones_mat_FxO/FrobNormMat);

}

}	


static void W_1st_sum_start(arg_struct_t* argStruct_p, int thread_iter){

int f_index;

bool last_element_queue_flag=false;

bool compute_flag;

while (threads_while_condition_1st_sum_flag&&(!(last_thread_to_sleep_flag_array[1][thread_iter]))){

	if ((!last_element_found_global_flag)&&(check_queue_global_flag)){

	/*Semaphore or mutex lock */
	sem_wait(&queue_sem);

	compute_flag=calculate_queue_index(argStruct_p, &f_index, &last_element_queue_flag); 

	/*Semaphore or mutex unlock */
	sem_post(&queue_sem);

	}

	if (compute_flag){

	compute_1st_sum(argStruct_p->T_fk_p, argStruct_p->V_nk_p, argStruct_p->expj_Phi_S_nkf_p, f_index, thread_iter);

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

static void W_2nd_sum_start(arg_struct_t* argStruct_p, int thread_iter){

int f_index;	
int m_index;

bool last_element_queue_flag=false;

bool compute_flag; 

while (threads_while_condition_2nd_sum_flag&&(!(last_thread_to_sleep_flag_array[2][thread_iter]))){

	if ((!last_element_found_global_flag)&&(check_queue_global_flag)){

	/*Semaphore or mutex lock */
	sem_wait(&queue_sem);

	compute_flag=calculate_pair_queue_indices(argStruct_p, &m_index, &f_index, &last_element_queue_flag); 

	/*Semaphore or mutex unlock */
	sem_post(&queue_sem);

	}

	if (compute_flag){

	compute_2nd_sum(argStruct_p->Xhat_low_fnm_p, argStruct_p->E_conj_fnm_p, m_index, f_index);

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

void* W_start(void* arg){

int thread_iter;

thread_arg_t* threadArg_p;
arg_struct_t* argStruct_p;

threadArg_p=(thread_arg_t*)arg;
argStruct_p=threadArg_p->argStruct_p;
thread_iter=threadArg_p->thread_iter; 

W_1st_sum_start(argStruct_p, thread_iter);

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
sem_wait(&wait_till_queue_populated_sem);

}

if (last_thread_to_sleep_flag_array[1][thread_iter]){

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
sem_post(&wait_till_queue_populated_sem);
sem_post(&wait_till_queue_populated_sem);
sem_post(&wait_till_queue_populated_sem);

}

W_2nd_sum_start(argStruct_p, thread_iter);
	
if (last_thread_to_sleep_flag_array[2][thread_iter]){



/*Call serial functions*/

compute_stuff_last_thread(argStruct_p->W_fom_p, argStruct_p->expj_Phi_W_fom_p);

normalize_W_local(argStruct_p->W_fom_p);

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

