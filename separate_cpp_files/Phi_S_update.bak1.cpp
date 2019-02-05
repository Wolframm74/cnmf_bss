#include "local_inc.hpp"

cx_cube::fixed<F_static, K_static, N_static> phase_shift_tensor_cx_FKN;
cx_cube::fixed<N_static, K_static, F_static> phase_shift_tensor_cx_NKF;

cx_cube::fixed<F_static, K_static, N_static> reference_phase_tensor_cx_FKN;
cx_cube::fixed<N_static, K_static, F_static> reference_phase_tensor_cx_NKF;
cube::fixed<F_static, K_static, N_static> ones_tensor_real_FKN;
cube::fixed<N_static, K_static, F_static> ones_tensor_real_NKF;

/*2nd sum stuff*/
/*cube::fixed<F_static, N_static, K_static> outtensor_2nd_real_FNK;
cx_cube::fixed<F_static, N_static, K_static> outtensor_2nd_cx_FNK;*/

cube::fixed<F_static, K_static, N_static> outtensor_2nd_real_FKN;
cx_cube::fixed<F_static, K_static, N_static> outtensor_2nd_cx_FKN;

cube::fixed<N_static, K_static, F_static> outtensor_2nd_real_NKF;
cx_cube::fixed<N_static, K_static, F_static> outtensor_2nd_cx_NKF;

/*output f,n,k related stuff*/
/*mat::fixed<F_static, N_static> dummy_mat_real_FN;
cx_mat::fixed<F_static, N_static> dummy_mat_cx_FN;
cx_cube::fixed<F_static, N_static, K_static> outtensor_target_cx_FNK;*/	/*Do a mexCallMATLAB() on this */

cx_cube::fixed<F_static, K_static, N_static> tensor_dummy_cx_FKN;	/*let this one be related to the cx tensor*/
cx_cube::fixed<N_static, K_static, F_static> tensor_dummy_cx_NKF; 

cx_cube::fixed<F_static, K_static, N_static> outtensor_target_cx_FKN;	/*let this one be related to the real tensor*/
cx_cube::fixed<N_static, K_static, F_static> outtensor_target_cx_NKF; 

static cube::fixed<F_static, K_static, N_static> absval_tensor_FKN;	/*let this one be related to the real tensor*/
static cube::fixed<N_static, K_static, F_static> absval_tensor_NKF; 

/*Phase_FKN and Phase_NKF's relationship to the outtensor_target's is conj(outtensor_target_cx*___) */

cube::fixed<F_static, M_static, N_static> Xhat_low_fmn_local;
cx_cube::fixed<F_static, M_static, N_static> E_conj_fmn_local;

static mat::fixed<O_static, L_static> dummy_mat_OL_local1[NUM_WORKER_THREADS];
static cube::fixed<N_static, K_static, F_static> Partial_wrt_Phi_S_nkf_out;
static cube::fixed<F_static, K_static, N_static> Partial_wrt_Phi_S_fkn_out;

static void compute_2nd_sum(mat* T_fk_p, mat* V_nk_p, cube* Xhat_low_fnm_p, cx_cube* E_conj_fnm_p, cx_cube* expj_Phi_S_nkf_p, cx_cube* expj_Phi_S_fkn_p, int n_index, int thread_iter, cube* W_fom_p, mat* Z_ol_p, mat* Y_lk_p){

double outvalue_other_thing;

#ifdef DEBUG	
mexPrintf("Phi_S_update: pthread_self(): %d, compute_2nd_sum(), n_index=%d\n", (int)pthread_self(), n_index);		
#endif	

int f_index; 
int m_index; 

int k_index;

for (m_index=0; m_index<M_static; m_index++){

Xhat_low_fmn_local.slice(n_index).col(m_index)=(*Xhat_low_fnm_p).slice(m_index).col(n_index);

E_conj_fmn_local.slice(n_index).col(m_index)=(*E_conj_fnm_p).slice(m_index).col(n_index);

}


for (f_index=0; f_index<F_static; f_index++){

/*Dimensions of matrix product: 1xK=(1xM)*(MxK) */

/*outtensor_2nd_real_FNK.subcube(f_index,n_index,0,f_index,n_index,K_static-1)=(*Xhat_low_fnm_p).subcube(f_index,n_index,0,f_index,n_index,M_static-1)*Xhat_outtensor_real_FMK.subcube(f_index,0,0,f_index,M_static-1,K_static-1);*/
/*outtensor_2nd_real_FKN.slice(n_index).row(f_index)=(*Xhat_low_fmn_p).slice(n_index).row(f_index)*(Xhat_outtensor_real_MKF).slice(f_index);*/
outtensor_2nd_real_FKN.slice(n_index).row(f_index)=Xhat_low_fmn_local.slice(n_index).row(f_index)*(Xhat_outtensor_real_MKF).slice(f_index);	
outtensor_2nd_real_NKF.slice(f_index).row(n_index)=outtensor_2nd_real_FKN.slice(n_index).row(f_index);

/*outtensor_2nd_cx_FNK.subcube(f_index,n_index,0,f_index,n_index,K_static-1)=(*E_conj_fmn_p).subcube(f_index,n_index,0,f_index,n_index,M_static-1)*Xhat_outtensor_cx_FMK.subcube(f_index,0,0,f_index,M_static-1,K_static-1);*/
/*outtensor_2nd_cx_FKN.slice(n_index).row(f_index)=(*E_conj_fmn_p).slice(n_index).row(f_index)*(Xhat_outtensor_cx_MKF).slice(f_index);*/
outtensor_2nd_cx_FKN.slice(n_index).row(f_index)=E_conj_fmn_local.slice(n_index).row(f_index)*(Xhat_outtensor_cx_MKF).slice(f_index);
outtensor_2nd_cx_NKF.slice(f_index).row(n_index)=outtensor_2nd_cx_FKN.slice(n_index).row(f_index);

/*Do a few more steps to compute the dummy and output tensors while you're here*/
tensor_dummy_cx_FKN.slice(n_index).row(f_index)=(outtensor_2nd_cx_FKN.slice(n_index).row(f_index))%((*T_fk_p).row(f_index))%((*V_nk_p).row(n_index));

tensor_dummy_cx_NKF.slice(f_index).row(n_index)=tensor_dummy_cx_FKN.slice(n_index).row(f_index);

outtensor_target_cx_FKN.slice(n_index).row(f_index)=conj((*expj_Phi_S_fkn_p).slice(n_index).row(f_index))%((*T_fk_p).row(f_index))%((*V_nk_p).row(n_index));	/*Should be noted that this step should be properly related to the Phi_S=-phase(arg) step */

outtensor_target_cx_FKN.slice(n_index).row(f_index)=(outtensor_target_cx_FKN.slice(n_index).row(f_index))%(outtensor_2nd_real_FKN.slice(n_index).row(f_index));

outtensor_target_cx_NKF.slice(f_index).row(n_index)=outtensor_target_cx_FKN.slice(n_index).row(f_index);

outtensor_target_cx_FKN.slice(n_index).row(f_index)=outtensor_target_cx_FKN.slice(n_index).row(f_index)+tensor_dummy_cx_FKN.slice(n_index).row(f_index);

outtensor_target_cx_NKF.slice(f_index).row(n_index)=outtensor_target_cx_NKF.slice(f_index).row(n_index)+tensor_dummy_cx_NKF.slice(f_index).row(n_index);

/*OTHER THING*/
	for (k_index=0; k_index<K_static; k_index++){

		dummy_mat_OL_local1[thread_iter].zeros();

		for (m_index=0; m_index<M_static; m_index++){

			/*accumulate something*/
			dummy_mat_OL_local1[thread_iter]=dummy_mat_OL_local1[thread_iter]+((*Xhat_low_fnm_p)(f_index, n_index, m_index))*(*Z_ol_p)%kron( trans((*W_fom_p).slice(m_index).row(f_index)), trans((*Y_lk_p).col(k_index)) );

		}

		outvalue_other_thing=accu(dummy_mat_OL_local1[thread_iter]);

		/*populate the fkn'th and nkf'th phases*/
		Partial_wrt_Phi_S_nkf_out(n_index, k_index, f_index)=-2*real(outtensor_target_cx_NKF(n_index, k_index, f_index)*outvalue_other_thing);

		Partial_wrt_Phi_S_fkn_out(f_index, k_index, n_index)=-2*real(outtensor_target_cx_FKN(f_index, k_index, n_index)*outvalue_other_thing);

	}

/*END OTHER THING*/	

}

}

static void Phi_S_2nd_sum_start(arg_struct_t* argStruct_p, int thread_iter){

int n_index;

bool last_element_queue_flag=false;

bool compute_flag; 

while (threads_while_condition_2nd_sum_flag&&(!(last_thread_to_sleep_flag_array[0][thread_iter]))){

	if ((!last_element_found_global_flag)&&(check_queue_global_flag)){

	/*Semaphore or mutex lock */
	sem_wait(&queue_sem);

	compute_flag=calculate_queue_index(argStruct_p, &n_index, &last_element_queue_flag);

	/*Semaphore or mutex unlock */
	sem_post(&queue_sem);

	}

	if (compute_flag){

	compute_2nd_sum(argStruct_p->T_fk_p, argStruct_p->V_nk_p, argStruct_p->Xhat_low_fnm_p, argStruct_p->E_conj_fnm_p, argStruct_p->expj_Phi_S_nkf_p, argStruct_p->expj_Phi_S_fkn_p, n_index, thread_iter, argStruct_p->W_fom_p, argStruct_p->Z_ol_p, argStruct_p->Y_lk_p);

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

} /*end while*/


}

static cx_cube::fixed<N_static, K_static, F_static> expj_Partial_wrt_Phi_S_nkf_out;
static cx_cube::fixed<F_static, K_static, N_static> expj_Partial_wrt_Phi_S_fkn_out;

void* Phi_S_start(void* arg){

double epsilon_value;

int i_iter;
int thread_iter;
cx_double local_cx_double;

thread_arg_t* threadArg_p;
arg_struct_t* argStruct_p;

threadArg_p=(thread_arg_t*)arg;
argStruct_p=threadArg_p->argStruct_p;
thread_iter=threadArg_p->thread_iter; 

epsilon_value=0.000000000001;

Phi_S_2nd_sum_start(argStruct_p, thread_iter);

if (last_thread_to_sleep_flag_array[0][thread_iter]){

/*POSTPONE THIS~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*absval_tensor_NKF=abs(outtensor_target_cx_NKF);
absval_tensor_FKN=abs(outtensor_target_cx_FKN);

outtensor_target_cx_NKF.elem(find(absval_tensor_NKF==0)).fill(0.00000001);
outtensor_target_cx_FKN.elem(find(absval_tensor_FKN==0)).fill(0.00000001);

(absval_tensor_NKF).elem(find(absval_tensor_NKF<=0)).fill(0.00000001);
(absval_tensor_FKN).elem(find(absval_tensor_FKN<=0)).fill(0.00000001);

*(argStruct_p->expj_Phi_S_nkf_p)=conj(outtensor_target_cx_NKF/absval_tensor_NKF); //%phase_shift_tensor_cx_NKF;
*(argStruct_p->expj_Phi_S_fkn_p)=conj(outtensor_target_cx_FKN/absval_tensor_FKN); //%phase_shift_tensor_cx_FKN;

outtensor_target_cx_NKF.zeros();
outtensor_target_cx_FKN.zeros();*/

/*Phase shift the output, using the newly calculated partial derivatives*/
Partial_wrt_Phi_S_nkf_out;
Partial_wrt_Phi_S_fkn_out;

expj_Partial_wrt_Phi_S_nkf_out.set_real(cos(epsilon_value*Partial_wrt_Phi_S_nkf_out));
expj_Partial_wrt_Phi_S_nkf_out.set_imag(sin(epsilon_value*Partial_wrt_Phi_S_nkf_out));

expj_Partial_wrt_Phi_S_fkn_out.set_real(cos(epsilon_value*Partial_wrt_Phi_S_fkn_out));
expj_Partial_wrt_Phi_S_fkn_out.set_imag(sin(epsilon_value*Partial_wrt_Phi_S_fkn_out));

(*(argStruct_p->expj_Phi_S_nkf_p))=(*(argStruct_p->expj_Phi_S_nkf_p))%conj(expj_Partial_wrt_Phi_S_nkf_out);
(*(argStruct_p->expj_Phi_S_fkn_p))=(*(argStruct_p->expj_Phi_S_fkn_p))%conj(expj_Partial_wrt_Phi_S_fkn_out);

outtensor_target_cx_NKF.zeros();
outtensor_target_cx_FKN.zeros();

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

if (!last_thread_to_sleep_flag_array[0][thread_iter]){

sem_post(&checkpoint_sem);

}

threads_exit_and_signal(NUM_WORKER_THREADS);
     
}


static rowvec::fixed<K_static> L1_norm_rowvec;

void Phi_S_L1_norm_V_checker(cx_cube* expj_Phi_S_nkf_p, cx_cube* expj_Phi_S_fkn_p, mat* V_nk_p){

L1_norm_rowvec=sum(*V_nk_p, 0);

int k_index, f_index, n_index;

for (k_index=0; k_index<K_static; k_index++){

	if (L1_norm_rowvec(k_index)<0.3){

	for (f_index=0; f_index<F_static; f_index++){

		for (n_index=0; n_index<N_static; n_index++){

			(*expj_Phi_S_nkf_p)(n_index, k_index, f_index)=1;

			(*expj_Phi_S_fkn_p)(f_index, k_index, n_index)=1;

		}

	}

	}

}


}