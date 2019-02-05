#include "local_inc.hpp"

static cx_rowvec::fixed<L_static> Phi_S_complex_argument_m4_fnk_spreadcol_1xL[NUM_WORKER_THREADS_COMPLEX_ARGUMENT_COSTFUN_M4];
static cx_rowvec::fixed<L_static> Phi_S_complex_argument_m4_fng_spreadcol_1xL[NUM_WORKER_THREADS_COMPLEX_ARGUMENT_COSTFUN_M4];

cx_cube::fixed<N_static, K_static, F_static> Phi_S_complex_argument_m4_outmat_NKF;
cx_cube::fixed<F_static, K_static, N_static> Phi_S_complex_argument_m4_outmat_FKN;
/*static cx_cube::fixed<F_static, N_static, K_static> Phi_S_complex_argument_m4_outmat_FNK;*/
/*static mat::fixed<F_static, N_static, K_static> Phi_S_complex_argument_m4_fng_outmat_FNK;*/

/*
int complex_argument_costfun_m4_m_index_b;
int complex_argument_costfun_m4_m_index_a;
*/

static void Phi_S_complex_argument_m4_compute_output_at_pair_indices_do_work(cx_cube* Xhat_fnm_p, cube* Xhat_low_fnm_p, mat* Z_ol_p, mat* Y_lk_p, mat* T_fk_p, mat* V_nk_p, cube* W_fom_p, cx_cube* W_fom_cx_p, cx_cube* expj_Phi_S_nkf_p, int f_index, int n_index, int thread_iter){

int k_index;

/*Phi_S(f,n,k) stuff: Formula: real(( 1st term + 2nd term )%(3rd term))*/	

// 1st term	
/*conj(complex_argument_costfun_m4_Error_mat_FN(f_index, n_index))*conj((*Xhat_fnm_p)(f_index, n_index, complex_argument_costfun_m4_m_index_a))*kron( trans((*W_fom_cx_p).slice(complex_argument_costfun_m4_m_index_b).row(f_index)) , trans(ones_col_Lx1) )*/

// 2nd term	
/*complex_argument_costfun_m4_Magnitude_model_mat_FN(f_index, n_index)*((*Xhat_low_fnm_p)(f_index, n_index, complex_argument_costfun_m4_m_index_a))*conj((*expj_Phi_S_nkf_p)(n_index, k_index, f_index))*kron( trans((W_fom_p).slice(complex_argument_costfun_m4_m_index_b).row(f_index)) , trans(ones_col_Lx1) )*/

// 3rd term	
/*((*V_nk_p)(n_index, k_index))*((*T_fk_p)(f_index, k_index))*(*Z_ol_p)%kron(ones_col_Ox1, trans((*Y_lk_p).col(k_index)) )*/

/*Phi_S(f,n,g) den stuff: Formula: (  1st term + 2nd term )%(3rd term)+4th term */	

// 1st term	
/*(complex_argument_costfun_m4_Error_mat_FN(f_index, n_index))*conj((*Xhat_fnm_p)(f_index, n_index, complex_argument_costfun_m4_m_index_b))*kron( trans((*W_fom_cx_p).slice(complex_argument_costfun_m4_m_index_a).row(f_index)) , trans(ones_col_Lx1) )*/

// 2nd term	
/*complex_argument_costfun_m4_Magnitude_model_mat_FN(f_index, n_index)*((*Xhat_low_fnm_p)(f_index, n_index, complex_argument_costfun_m4_m_index_b))*conj((*expj_Phi_S_nkf_p)(n_index, k_index, f_index))*kron( trans((W_fom_p).slice(complex_argument_costfun_m4_m_index_a).row(f_index)) , trans(ones_col_Lx1) )*/

// 3rd term	
/*((*V_nk_p)(n_index, k_index))*((*T_fk_p)(f_index, k_index))*(*Z_ol_p)%kron(ones_col_Ox1, trans((*Y_lk_p).col(k_index)) )*/

Phi_S_complex_argument_m4_fnk_spreadcol_1xL[thread_iter].zeros();
Phi_S_complex_argument_m4_fng_spreadcol_1xL[thread_iter].zeros();

for (k_index=0; k_index<K_static; k_index++){

Phi_S_complex_argument_m4_fnk_spreadcol_1xL[thread_iter]=((*expj_Phi_S_nkf_p)(n_index, k_index, f_index))*((*Xhat_low_fnm_p)(f_index, n_index, complex_argument_costfun_m4_m_index_a))*(complex_argument_costfun_m4_Magnitude_model_mat_FN(f_index, n_index))*((*T_fk_p)(f_index, k_index))*((*V_nk_p)(n_index, k_index))*(trans((*Y_lk_p).col(k_index))%Xhat_outtensor_real_FLM.slice(complex_argument_costfun_m4_m_index_b).row(f_index));
Phi_S_complex_argument_m4_fng_spreadcol_1xL[thread_iter]=((*expj_Phi_S_nkf_p)(n_index, k_index, f_index))*((*Xhat_low_fnm_p)(f_index, n_index, complex_argument_costfun_m4_m_index_b))*(complex_argument_costfun_m4_Magnitude_model_mat_FN(f_index, n_index))*((*T_fk_p)(f_index, k_index))*((*V_nk_p)(n_index, k_index))*(trans((*Y_lk_p).col(k_index))%Xhat_outtensor_real_FLM.slice(complex_argument_costfun_m4_m_index_a).row(f_index));

Phi_S_complex_argument_m4_outmat_NKF(n_index, k_index, f_index)=(((*T_fk_p)(f_index, k_index))*((*V_nk_p)(n_index, k_index))*(complex_argument_costfun_m4_Error_mat_FN(f_index, n_index))*((*Xhat_fnm_p)(f_index, n_index, complex_argument_costfun_m4_m_index_a))*conj(Xhat_out_4way_tensor_cx_nkfm[complex_argument_costfun_m4_m_index_b](n_index, k_index, f_index)))+(((*T_fk_p)(f_index, k_index))*((*V_nk_p)(n_index, k_index))*(conj(complex_argument_costfun_m4_Error_mat_FN(f_index, n_index)))*((*Xhat_fnm_p)(f_index, n_index, complex_argument_costfun_m4_m_index_b))*conj(Xhat_out_4way_tensor_cx_nkfm[complex_argument_costfun_m4_m_index_a](n_index, k_index, f_index)));

Phi_S_complex_argument_m4_outmat_NKF(n_index, k_index, f_index)=Phi_S_complex_argument_m4_outmat_NKF(n_index, k_index, f_index)+accu(2*Phi_S_complex_argument_m4_fnk_spreadcol_1xL[thread_iter])+accu(2*Phi_S_complex_argument_m4_fng_spreadcol_1xL[thread_iter]);

Phi_S_complex_argument_m4_outmat_NKF(n_index, k_index, f_index)=-Phi_S_complex_argument_m4_outmat_NKF(n_index, k_index, f_index)*conj((*expj_Phi_S_nkf_p)(n_index, k_index, f_index));

Phi_S_complex_argument_m4_outmat_FKN(f_index, k_index, n_index)=Phi_S_complex_argument_m4_outmat_NKF(n_index, k_index, f_index);

}



}

static cube::fixed<F_static, K_static, N_static> absval_tensor_FKN;	
static cube::fixed<N_static, K_static, F_static> absval_tensor_NKF; 


static void Phi_S_complex_argument_m4_compute_outstanding_work_last_thread_main(cx_cube* expj_Phi_S_nkf_p, cx_cube* expj_Phi_S_fkn_p, arg_struct_t* argStruct_p){

mat* T_fk_p=argStruct_p->T_fk_p;
mat* V_nk_p=argStruct_p->V_nk_p;
double mu_value;

/*mu_value=10;*/
double mu_value1=0.01;
/*double mu_value2=0.001;*/
double mu_value2=0.0;

/*(*expj_Phi_S_nkf_p)=(*expj_Phi_S_nkf_p)-2*(Phi_S_complex_argument_m4_outmat_NKF+outtensor_target_cx_NKF);*/

/*FKN update below*/

/*(*expj_Phi_S_nkf_p)=(*expj_Phi_S_nkf_p)+2*mu_value*(outtensor_target_cx_NKF);*/

/*(*expj_Phi_S_nkf_p)=(*expj_Phi_S_nkf_p)+2*mu_value*(Phi_S_complex_argument_m4_outmat_NKF+outtensor_target_cx_NKF);*/
(*expj_Phi_S_nkf_p)=(*expj_Phi_S_nkf_p)+2*(mu_value2*Phi_S_complex_argument_m4_outmat_NKF+mu_value1*outtensor_target_cx_NKF);

/*(*expj_Phi_S_nkf_p)=(*expj_Phi_S_nkf_p)-2*mu_value*(Phi_S_complex_argument_m4_outmat_NKF+outtensor_target_cx_NKF);*/

/*(*expj_Phi_S_nkf_p)=(*expj_Phi_S_nkf_p)+2*mu_value*(Phi_S_complex_argument_m4_outmat_NKF);*/

/*FKN update below*/

/*(*expj_Phi_S_fkn_p)=(*expj_Phi_S_fkn_p)+2*mu_value*(outtensor_target_cx_FKN);*/

/*(*expj_Phi_S_fkn_p)=(*expj_Phi_S_fkn_p)+2*mu_value*(Phi_S_complex_argument_m4_outmat_FKN+outtensor_target_cx_FKN);*/

(*expj_Phi_S_fkn_p)=(*expj_Phi_S_fkn_p)+2*(mu_value2*Phi_S_complex_argument_m4_outmat_FKN+mu_value1*outtensor_target_cx_FKN);

/*(*expj_Phi_S_fkn_p)=(*expj_Phi_S_fkn_p)-2*mu_value*(Phi_S_complex_argument_m4_outmat_FKN+outtensor_target_cx_FKN);*/

/*(*expj_Phi_S_fkn_p)=(*expj_Phi_S_fkn_p)+2*mu_value*(Phi_S_complex_argument_m4_outmat_FKN);*/

// Can later modify the code below to instead compute the rank one approximation at each k bin and spreak it into T_fk and V_nk 


// rotate
// Rotating function defined in complex_argument_costfun_m1.hpp
rotate_expj_Phi_S_nkf_p_to_FNK(argStruct_p->expj_Phi_S_nkf_p);
// Can now safely access an up-to-date rotated_expj_Phi_S_FNK, from this point on

//On the way out compute a rank 1 approximation at each k_index and update T_fk and V_nk at each k_index
mat U_svd;
vec s_svd;
mat V_svd;

double max_T_kth;
double max_V_kth;

int k_index; 

for (k_index=0; k_index<(K_static); k_index++){

	svd(U_svd, s_svd, V_svd, abs(rotated_expj_Phi_S_FNK.slice(k_index)) );

	if (max(abs(s_svd(0)*U_svd.col(0)))>1){

	max_T_kth=max(abs(s_svd(0)*U_svd.col(0)));

	} else {

	max_T_kth=1;

	}

	if (max(abs(V_svd.col(0)))>1){

	max_V_kth=max(abs(V_svd.col(0)));

	} else {

	max_V_kth=1;

	}

	(*T_fk_p).col(k_index)=((abs(s_svd(0)*U_svd.col(0)))/max_T_kth)%(*T_fk_p).col(k_index);

	(*V_nk_p).col(k_index)=((abs(V_svd.col(0)))/max_V_kth)%(*V_nk_p).col(k_index);

}

// With the updated rank1 approximation updated into T_fk and V_nk, can now Divide out the absolute value to obtain unit magnitude in each FxNxK bin. 
absval_tensor_NKF=abs((*expj_Phi_S_nkf_p));
absval_tensor_FKN=abs((*expj_Phi_S_fkn_p));

(*expj_Phi_S_nkf_p).elem(find(absval_tensor_NKF==0)).fill(1);
(*expj_Phi_S_fkn_p).elem(find(absval_tensor_FKN==0)).fill(1);



(*expj_Phi_S_nkf_p)=(*expj_Phi_S_nkf_p)/abs((*expj_Phi_S_nkf_p));

//FKN update below
(*expj_Phi_S_fkn_p)=(*expj_Phi_S_fkn_p)/abs((*expj_Phi_S_fkn_p));

/*clean up*/
outtensor_target_cx_NKF.zeros();
outtensor_target_cx_FKN.zeros();

}

static void Phi_S_complex_argument_m4_compute_outstanding_work_last_thread_alternate(cx_cube* expj_Phi_S_nkf_p, arg_struct_t* argStruct_p){

/*FKN update below*/

}

static void Phi_S_complex_argument_m4_compute_output_at_pair_indices_complete_work_wrapper(arg_struct_t* argStruct_p, int thread_iter){


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

	Phi_S_complex_argument_m4_compute_output_at_pair_indices_do_work(argStruct_p->Xhat_fnm_p, argStruct_p->Xhat_low_fnm_p, argStruct_p->Z_ol_p, argStruct_p->Y_lk_p, argStruct_p->T_fk_p, argStruct_p->V_nk_p, argStruct_p->W_fom_p, argStruct_p->W_fom_cx_p, argStruct_p->expj_Phi_S_nkf_p, f_index, n_index, thread_iter);

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

void* Phi_S_complex_argument_m4_start(void* arg){

int i_iter;
int thread_iter;

thread_arg_t* threadArg_p;
arg_struct_t* argStruct_p;

threadArg_p=(thread_arg_t*)arg;
argStruct_p=threadArg_p->argStruct_p;
thread_iter=threadArg_p->thread_iter; 

Phi_S_complex_argument_m4_compute_output_at_pair_indices_complete_work_wrapper(argStruct_p, thread_iter);

if (last_thread_to_sleep_flag_array[0][thread_iter]){

/*Phi_S_complex_argument_m4_compute_outstanding_work_last_thread_main(argStruct_p->expj_Phi_S_nkf_p, argStruct_p);*/
Phi_S_complex_argument_m4_compute_outstanding_work_last_thread_main(argStruct_p->expj_Phi_S_nkf_p, argStruct_p->expj_Phi_S_fkn_p, argStruct_p);

sem_wait(&threads_while_condition_0th_sum_sem);
threads_while_condition_0th_sum_flag=false;	
sem_post(&threads_while_condition_0th_sum_sem);

/*Sem post to wake NUM_WORKER_THREADS_COMPLEX_ARGUMENT_COSTFUN_M4-1 of the sleeping threads*/

for (i_iter=0; i_iter<(NUM_WORKER_THREADS_COMPLEX_ARGUMENT_COSTFUN_M4-1); i_iter++){
	sem_post(&sleep_sem);
}

/*Sem wait NUM_WORKER_THREADS_COMPLEX_ARGUMENT_COSTFUN_M4-1 times so that you only proceed from this point once that many of the sleeping flags have high fived and successfully made it out of the "while traps" */

for (i_iter=0; i_iter<(NUM_WORKER_THREADS_COMPLEX_ARGUMENT_COSTFUN_M4-1); i_iter++){
	sem_wait(&checkpoint_sem);
}

}

if (!last_thread_to_sleep_flag_array[0][thread_iter]){

sem_post(&checkpoint_sem);

}

threads_exit_and_signal(NUM_WORKER_THREADS_COMPLEX_ARGUMENT_COSTFUN_M4);
     

}