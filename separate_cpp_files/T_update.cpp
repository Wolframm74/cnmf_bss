#include "local_inc.hpp"

/*2nd sum stuff*/
/*cube::fixed<F_static, K_static, M_static> outtensor_2nd_real_FKM;
cx_cube::fixed<F_static, K_static, M_static> outtensor_2nd_cx_FKM;*/

/*3rd sum stuff*/
mat::fixed<F_static, K_static> outmat_3rd_real_FK_den;
mat::fixed<F_static, K_static> outmat_3rd_real_FK_num;
/*cx_mat::fixed<F_static, K_static> outmat_3rd_cx_FK;*/

static void T_primary_auxfun_compute_output_at_pair_indices_do_work(cx_cube* Xhat_fnm_p, cube* Xhat_low_fnm_p, mat* Z_ol_p, mat* Y_lk_p, mat* V_nk_p, cube* W_fom_p, cx_cube* W_fom_cx_p, cx_cube* expj_Phi_S_nkf_p, cx_cube* E_conj_fnm_p, int f_index, int k_index, int thread_iter){

int m_index, l_index;

/*double den_value_local;
double num_value_local;*/

double real_accu_value_local=0;
cx_double cx_accu_value_local=0;

for (m_index=0; m_index<M_static; m_index++){

	cx_accu_value_local=cx_accu_value_local+accu(strans((*E_conj_fnm_p).slice(m_index).row(f_index))%(*V_nk_p).col(k_index)%((*expj_Phi_S_nkf_p).slice(f_index).col(k_index))%(Xhat_out_4way_tensor_cx_nkfm[m_index].slice(f_index).col(k_index)));

	for (l_index=0; l_index<L_static; l_index++){

	real_accu_value_local=real_accu_value_local+(Xhat_outtensor_real_FLM(f_index, l_index, m_index))*((*Y_lk_p)(l_index, k_index))*accu((*V_nk_p).col(k_index)%trans((*Xhat_low_fnm_p).slice(m_index).row(f_index)));

	}

}

outmat_3rd_real_FK_den(f_index, k_index)=real_accu_value_local;
outmat_3rd_real_FK_num(f_index, k_index)=outmat_3rd_real_FK_den(f_index, k_index)+real(cx_accu_value_local);

/*(*T_fk_p)(f_index, k_index)=((*T_fk_p)(f_index, k_index))*(num_value_local/den_value_local);*/

}

void projfunc_wrapper_Tfk(mat* T_fk_p, mxArray *T_fk_mxArray_p){

mxArray* output_arg_p[1];

armaSetPr(T_fk_mxArray_p, *T_fk_p);

/*mexCallMATLAB(1, &output_arg_p, 1, &plhs[2], "projfunc_wrapper_matlab");*/
/*mexCallMATLAB(0, NULL, 1, &T_fk_mxArray_p, "projfunc_wrapper_matlab");*/
mexCallMATLAB(1, &output_arg_p[0], 1, &T_fk_mxArray_p, "projfunc_wrapper_matlab_Tfk");

/*memcpy(mxGetPr(T_fk_mxArray_p), mxGetPr(output_arg_p), O_static*L_static*sizeof(double));*/

/*arrayops::copy((*T_fk_p).memptr(), mxGetPr(output_arg_p[0]), O_static*L_static);*/

(*T_fk_p).set_real(armaGetPr(output_arg_p[0], true)); 	

}

/*T_primary_auxfun_compute_outstanding_work_last_thread(argStruct_p->T_fk_p);*/
static void T_primary_auxfun_compute_outstanding_work_last_thread(mat* T_fk_p){

/*outmat_3rd_real_FK_den.zeros();
outmat_3rd_cx_FK.zeros();

// Integrate out m
for (int m_iter=0; m_iter<M_static; m_iter++){	
outmat_3rd_real_FK_den=outmat_3rd_real_FK_den+Xhat_outtensor_real_FKM.slice(m_iter)%outtensor_2nd_real_FKM.slice(m_iter);

outmat_3rd_cx_FK=outmat_3rd_cx_FK+Xhat_outtensor_cx_FKM.slice(m_iter)%outtensor_2nd_cx_FKM.slice(m_iter);
}

outmat_3rd_real_FK_num=outmat_3rd_real_FK_den+real(outmat_3rd_cx_FK);*/

/*outmat_3rd_real_FK_num.print("outmat_3rd_real_FK_num");
outmat_3rd_real_FK_den.print("outmat_3rd_real_FK_den");*/

/*Normalize by l1-norms*/
/*outmat_3rd_real_FK_num=outmat_3rd_real_FK_num/accu(outmat_3rd_real_FK_num);

outmat_3rd_real_FK_den=outmat_3rd_real_FK_den/accu(outmat_3rd_real_FK_den);*/

/*POSTPONE THIS~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
/*(*T_fk_p)=(*T_fk_p)%(outmat_3rd_real_FK_num/outmat_3rd_real_FK_den);*/

/*(*T_fk_p).print("T_fk INSIDE MEX, LINE 49:");*/

/*Set negative elements to zero*/
/*(*T_fk_p).elem(find((*T_fk_p)<0)).zeros();*/
(*T_fk_p).elem(find((*T_fk_p)<=0)).fill(0.00000001);

/*(*T_fk_p).elem(find((*T_fk_p)>10000)).fill(10000);*/

/*Subtract off the current mean*/
/*Store the current mean as a reference so that you can shift it back to this mean after doing and T_fk post processing. ex: sparsity post processing*/
accum_T_fk=accu(*T_fk_p);

/**T_fk_p=(*T_fk_p)-(accum_T_fk/(((double)F_static)*((double)K_static)))*ones(F_static, K_static);*/

/*Set a new mean forcibly*/
/**T_fk_p=(*T_fk_p)+sqrt(1/((double)K_static))*ones(F_static, K_static);*/

/*(*T_fk_p).print("T_fk INSIDE MEX: LINE 61");*/

}

static void T_primary_auxfun_compute_output_at_pair_indices_complete_work_wrapper(arg_struct_t* argStruct_p, int thread_iter){

int f_index;
int k_index;

bool last_element_queue_flag=false;

bool compute_flag; 

while (threads_while_condition_0th_sum_flag&&(!(last_thread_to_sleep_flag_array[0][thread_iter]))){

	if ((!last_element_found_global_flag)&&(check_queue_global_flag)){

	/*Semaphore or mutex lock */
	sem_wait(&queue_sem);

	compute_flag=calculate_pair_queue_indices(argStruct_p, &f_index, &k_index, &last_element_queue_flag); 

	/*Semaphore or mutex unlock */
	sem_post(&queue_sem);

	}

	if (compute_flag){

	/*T_primary_auxfun_compute_output_at_pair_indices_do_work(cx_cube* Xhat_fnm_p, cube* Xhat_low_fnm_p, mat* Z_ol_p, mat* Y_lk_p, mat* V_nk_p, cube* W_fom_p, cx_cube* W_fom_cx_p, cx_cube* expj_Phi_S_nkf_p, cx_cube* E_conj_fnm_p, int f_index, int k_index, int thread_iter);*/
	T_primary_auxfun_compute_output_at_pair_indices_do_work(argStruct_p->Xhat_fnm_p, argStruct_p->Xhat_low_fnm_p, argStruct_p->Z_ol_p, argStruct_p->Y_lk_p, argStruct_p->V_nk_p, argStruct_p->W_fom_p, argStruct_p->W_fom_cx_p, argStruct_p->expj_Phi_S_nkf_p, argStruct_p->E_conj_fnm_p, f_index, k_index, thread_iter);

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


void* T_primary_auxfun_start(void* arg){

int i_iter;
int thread_iter;

thread_arg_t* threadArg_p;
arg_struct_t* argStruct_p;

threadArg_p=(thread_arg_t*)arg;
argStruct_p=threadArg_p->argStruct_p;
thread_iter=threadArg_p->thread_iter; 

T_primary_auxfun_compute_output_at_pair_indices_complete_work_wrapper(argStruct_p, thread_iter);

if (last_thread_to_sleep_flag_array[0][thread_iter]){

T_primary_auxfun_compute_outstanding_work_last_thread(argStruct_p->T_fk_p);

sem_wait(&threads_while_condition_0th_sum_sem);
threads_while_condition_0th_sum_flag=false;	
sem_post(&threads_while_condition_0th_sum_sem);

/*Sem post to wake NUM_WORKER_THREADS-1 of the sleeping threads*/

for (i_iter=0; i_iter<(NUM_WORKER_THREADS-1); i_iter++){
	sem_post(&sleep_sem);
}

/*Sem wait NUM_WORKER_THREADS-1 times so that you only proceed from this point once that many of the sleeping flags have high fived and successfully made it out of the "while traps" */

for (i_iter=0; i_iter<(NUM_WORKER_THREADS-1); i_iter++){
	sem_wait(&checkpoint_sem);
}

}

if (!last_thread_to_sleep_flag_array[0][thread_iter]){

sem_post(&checkpoint_sem);

}

threads_exit_and_signal(NUM_WORKER_THREADS);
     
}

