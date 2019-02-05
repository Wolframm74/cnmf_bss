#include "local_inc.hpp"

/*Local variables (to this module) definitions, but globally visible */

/*populate_W_fom_cx() related variables: */
cube::fixed<F_static, O_static, M_static> dummy_tensor_fom;

/*threads_exit_and_signal()*/
int threads_exit_ctr; /*ctr*/
sem_t threads_exit_lock_sem; /*lock semaphore*/
sem_t threads_exit_sleep_sem; /*sleep semaphore*/

int which_update_helper;
FILE* Xhat_pFile_helper;

void threads_exit_and_signal(int num_threads){

bool last_thread_to_exit=false; 

mexPrintf("threads_exit_and_signal(), JUST ENTERED threads_exit_and_signal()!!!!!!!!!!!: pthread_self():%d, check_queue_global_flag=%d, threads_exit_ctr=%d, \n", (int)pthread_self(), (int)check_queue_global_flag, threads_exit_ctr);

/*Increments a global ctr specifying how many threads have already passed by this point*/
/*Counter needs to have been set (or reset) prior to this point between algorithm updates. Should be okay to set it to zero every time upon entry to update_engines() */

/*sem_wait() / "lock" the ctr (binary) semaphore */	
sem_wait(&threads_exit_lock_sem);

	/*Increment the counter*/
	threads_exit_ctr++;

mexPrintf("threads_exit_and_signal(), INSIDE threads_exit_and_signal()!!!!!!!!!!!: JUST INCREMENTED THREADS_EXIT_CTR pthread_self():%d, check_queue_global_flag=%d, threads_exit_ctr=%d, \n", (int)pthread_self(), (int)check_queue_global_flag, threads_exit_ctr);

	/*If the counter equals some specific value signifying that the current thread is the last thread to exit, set the last_thread_to_exit bool to true. */
	if (threads_exit_ctr==num_threads){
		last_thread_to_exit=true;
	}

/*sem_post() / "unlock" the binary semaphore*/
sem_post(&threads_exit_lock_sem);

if (last_thread_to_exit){

mexPrintf("threads_exit_and_signal(), INSIDE threads_exit_and_signal()!!!!!!!!!!!: FOUND THE LAST_THREAD_TO_EXIT THREAD!!!!!!!!!!!!!! pthread_self():%d, check_queue_global_flag=%d, threads_exit_ctr=%d, \n", (int)pthread_self(), (int)check_queue_global_flag, threads_exit_ctr);

/*Signal ( sem_post() ) w/ the sleep semaphore to the main thread to wake up*/
/*Main thread should then join all threads and safely wait until they have shut down to resume via pthread_join()'s' */	

	sem_post(&threads_exit_sleep_sem);

}

mexPrintf("threads_exit_and_signal(), EXITING threads_exit_and_signal() NOW!!!!!!!!!!!: pthread_self():%d, check_queue_global_flag=%d, threads_exit_ctr=%d, \n", (int)pthread_self(), (int)check_queue_global_flag, threads_exit_ctr);

}

void populate_queue_wrt_one_index(arg_struct_t* argStruct_p, int dim1_len){

int dim_iter; 

  /*Can maybe wrap this with a binary semaphore or mutex*/
	sem_wait(&queue_sem);

	argStruct_p->work_queue_size=dim1_len;

  /*If you suspect that different updates and their threads are overlapping / competing to access the "queue_intz" queue you can first assert that it should have its size reset to zero (ie: be empty) at this point */
	assert(argStruct_p->queue_intz.size()==0);

	  for (dim_iter=0; dim_iter<dim1_len; dim_iter++){

	    argStruct_p->queue_intz.push(dim_iter);

		mexPrintf("populate_queue_wrt_one_index: pthread_self(): %d, queue_intz.size():%d, dim_iter=%d, threads_asleep_ctr= %d \n", (int)pthread_self(), (int)(argStruct_p->queue_intz.size()), dim_iter, threads_asleep_ctr );

	  }

  	sem_post(&queue_sem);

  /*End wrap. */

}

void populate_queue_wrt_pair_indices(arg_struct_t* argStruct_p, int dim1_len, int dim2_len){

int d1, d2; 

  /*Can maybe wrap this with a binary semaphore or mutex*/
	sem_wait(&queue_sem);

	argStruct_p->work_queue_size=dim1_len*dim2_len;

  /*If you suspect that different updates and their threads are overlapping / competing to access the "queue_intz" queue you can first assert that it should have its size reset to zero (ie: be empty) at this point */
	assert(argStruct_p->queue_intz.size()==0);

	  for (d1=0; d1<dim1_len; d1++){

	  	for (d2=0; d2<dim2_len; d2++){

	  	/*Need to take note of these next two lines when popping elements off the front of the queue. */
		argStruct_p->queue_intz.push(d2);	  		
	    argStruct_p->queue_intz.push(d1);

	    mexPrintf("populate_queue_wrt_pair_indices; %s ;", argStruct_p->update_string);
		mexPrintf(" pthread_self(): %d, queue_intz.size():%d, d1=%d, d2=%d, threads_asleep_ctr=%d \n", (int)pthread_self(), (int)(argStruct_p->queue_intz.size()), d1, d2, threads_asleep_ctr );

		}

	  }

  	sem_post(&queue_sem);

  /*End wrap. */

}

bool calculate_queue_index(arg_struct_t* argStruct_p, int* int_ptr, bool* last_element_p){

	mexPrintf("calculate_queue_index; %s; pthread_self(): %d, queue_intz.size():%d, threads_asleep_ctr=%d \n", argStruct_p->update_string, (int)pthread_self(), (int)(argStruct_p->queue_intz.size()), threads_asleep_ctr );

	/*If the size is greater than 1 */
	if((argStruct_p->queue_intz.size())>1){

		/*Read from the front element of the queue */	
		*int_ptr=argStruct_p->queue_intz.front();

		/*Pop off the front element - you don't need it anymore */
		argStruct_p->queue_intz.pop();

		return true; 

	} else if((argStruct_p->queue_intz.size())==1){	/*Endif - else if queue size == 1 : */

		/*Set the b_update_last_queue_element=true; */
		*last_element_p=true;

		/*Set the global last queue element flag*/		
		last_element_found_global_flag=true; 

		mexPrintf("pthread_self(): %d FOUND THE LAST ELEMENT!!!!!!!!!", (int)pthread_self());

		/*Read from the front element of the queue */	
		*int_ptr=argStruct_p->queue_intz.front();

		/*Pop off the front element - you don't need it anymore */
		argStruct_p->queue_intz.pop();

		return true; 

	} else if ((argStruct_p->queue_intz.size())<=0)	{ /*End else if */

		return false;

	}

}

bool calculate_pair_queue_indices(arg_struct_t* argStruct_p, int* index1_p, int* index2_p, bool* last_element_p){

	mexPrintf("inside calculate_pair_queue_indices(): pthread_self(): %d  queue size: %d, threads_asleep_ctr=%d \n", (int) pthread_self(), argStruct_p->queue_intz.size(), threads_asleep_ctr );

	if (argStruct_p->which_update==1){
	fprintf(argStruct_p->Xhat_pFile, "inside calculate_pair_queue_indices(): pthread_self(): %d  queue size: %d, threads_asleep_ctr=%d \n", (int) pthread_self(), argStruct_p->queue_intz.size(), threads_asleep_ctr);		
	}

	/*If the size is greater than 1 */
	if((argStruct_p->queue_intz.size())>2){

		/*Read from the front element of the queue */	
		*index2_p=argStruct_p->queue_intz.front();

		/*Pop off the front element - you don't need it anymore */
		argStruct_p->queue_intz.pop();

		/*Read from the front element of the queue */	
		*index1_p=argStruct_p->queue_intz.front();

		/*Pop off the front element - you don't need it anymore */
		argStruct_p->queue_intz.pop();

		return true; 

	} else if((argStruct_p->queue_intz.size())==2){	/*Endif - else if queue size == 1 : */

		mexPrintf("inside calculate_pair_queue_indices(): line 128; pthread_self(): %d  queue size: %d, threads_asleep_ctr=%d \n", (int) pthread_self(), argStruct_p->queue_intz.size(), threads_asleep_ctr );
		
		if (argStruct_p->which_update==1){
		fprintf(argStruct_p->Xhat_pFile, "inside calculate_pair_queue_indices(): line 128; pthread_self(): %d  queue size: %d, threads_asleep_ctr=%d \n", (int) pthread_self(), argStruct_p->queue_intz.size(), threads_asleep_ctr );		
		}

		/*Set the b_update_last_queue_element=true; */
		*last_element_p=true;

		/*Set the global last queue element flag*/		
		last_element_found_global_flag=true; 

		mexPrintf("pthread_self(): %d FOUND THE LAST ELEMENT!!!!!!!!!", (int)pthread_self());

		if (argStruct_p->which_update==1){
		fprintf(argStruct_p->Xhat_pFile, "pthread_self(): %d FOUND THE LAST ELEMENT!!!!!!!!!", (int)pthread_self());		
		}


		/*Read from the front element of the queue */	
		*index2_p=argStruct_p->queue_intz.front();

		/*Pop off the front element - you don't need it anymore */
		argStruct_p->queue_intz.pop();

		/*Read from the front element of the queue */	
		*index1_p=argStruct_p->queue_intz.front();

		/*Pop off the front element - you don't need it anymore */
		argStruct_p->queue_intz.pop();

		mexPrintf("inside calculate_pair_queue_indices(): line 150; pthread_self(): %d  queue size: %d, threads_asleep_ctr=%d \n", (int) pthread_self(), argStruct_p->queue_intz.size(), threads_asleep_ctr );

		if (argStruct_p->which_update==1){
		fprintf(argStruct_p->Xhat_pFile, "inside calculate_pair_queue_indices(): line 150; pthread_self(): %d  queue size: %d, threads_asleep_ctr=%d \n", (int) pthread_self(), argStruct_p->queue_intz.size(), threads_asleep_ctr );		
		}


		return true; 

	} else if ((argStruct_p->queue_intz.size())<=1)	{ /*End else if */

		mexPrintf("inside calculate_pair_queue_indices(): line 156; pthread_self(): %d  queue size: %d, threads_asleep_ctr=%d \n", (int) pthread_self(), argStruct_p->queue_intz.size(), threads_asleep_ctr );

		if (argStruct_p->which_update==1){
		fprintf(argStruct_p->Xhat_pFile, "inside calculate_pair_queue_indices(): line 156; pthread_self(): %d  queue size: %d, threads_asleep_ctr=%d \n", (int) pthread_self(), argStruct_p->queue_intz.size(), threads_asleep_ctr );		
		}


		return false;

	}

}

void populate_W_fom_cx(cube* W_fom_p, cx_cube* W_fom_cx_p, cube* Phi_W_fom_p){

/*populate the real part*/
dummy_tensor_fom=(*W_fom_p)%cos(*Phi_W_fom_p);

(*W_fom_cx_p).set_real(dummy_tensor_fom);

/*populate the imag part*/	
dummy_tensor_fom=(*W_fom_p)%sin(*Phi_W_fom_p);

(*W_fom_cx_p).set_imag(dummy_tensor_fom);

}

void forward_mapping_Phase_S(){



}

void backward_mapping_Phase_S(){


	
}

void clear_queue_intz(std::queue<int> &queue_intz){

std::queue<int> dummy_queue;
std::swap(queue_intz, dummy_queue);

}

/*Argument list:

Xtilde_fnm, 		prhs[0]
Z_ol, 			plhs[2]	
Y_lk, 			plhs[3]	
T_fk, 			plhs[4]	
V_nk, 			plhs[5]	
W_fom, 			plhs[1]
expj_Phi_W_fom, 	prhs[18]
expj_Phi_S_fkn		plhs[7]

*/

void TDOA_update_module1_wrapper(mxArray *plhs[], const mxArray *prhs[], arg_struct_t* argStruct_p){

/*Xtilde_fnm, 		prhs[0]*/
engPutVariable(mlEngine_p, "Xtilde_fnm", prhs[0]);

/*Z_ol, 			plhs[2]*/	
armaSetPr(plhs[2], (*(argStruct_p->Z_ol_p)));
engPutVariable(mlEngine_p, "Z_ol", plhs[2]);

/*Y_lk, 			plhs[3]*/	
armaSetPr(plhs[3], (*(argStruct_p->Y_lk_p)));
engPutVariable(mlEngine_p, "Y_lk", plhs[3]);

/*T_fk, 			plhs[4]*/	
armaSetPr(plhs[4], (*(argStruct_p->T_fk_p)));
engPutVariable(mlEngine_p, "T_fk", plhs[4]);

/*V_nk, 			plhs[5]*/	
armaSetPr(plhs[5], (*(argStruct_p->V_nk_p)));
engPutVariable(mlEngine_p, "V_nk", plhs[5]);

/*W_fom, 		plhs[1]*/
armaSetCubePr(plhs[1], (*(argStruct_p->W_fom_p)));
engPutVariable(mlEngine_p, "W_fom", plhs[1]);

/*expj_Phi_W_fom, 	prhs[18]*/
engPutVariable(mlEngine_p, "expj_Phi_W_fom", prhs[18]);

/*expj_Phi_S_fkn		plhs[7]*/
armaSetCubeCx(plhs[7], *(argStruct_p->expj_Phi_S_fkn_p));
engPutVariable(mlEngine_p, "expj_Phi_S_fkn", plhs[7]);

/*example code*/
/*armaSetPr(plhs[5], (*(argStruct_p->V_nk_p)));
engPutVariable(mlEngine_p, "V_nk", plhs[5]);*/

engEvalString(mlEngine_p, "[Z_ol, W_fom] = tdoa_update_m1(Xtilde_fnm, Z_ol, Y_lk, T_fk, V_nk, W_fom, expj_Phi_W_fom, expj_Phi_S_fkn);");

plhs[2]=engGetVariable(mlEngine_p, "Z_ol");
(*(argStruct_p->Z_ol_p)).set_real(armaGetPr(plhs[2], true)); 

plhs[1]=engGetVariable(mlEngine_p, "W_fom");
(*(argStruct_p->W_fom_p)).set_real(armaGetCubePr(plhs[1], true)); 	

/*nonnegative clipping*/
(*(argStruct_p->Z_ol_p)).elem(find((*(argStruct_p->Z_ol_p))<=0)).fill(0.00000001);

(*(argStruct_p->Z_ol_p)).elem(find((*(argStruct_p->Z_ol_p))>1)).fill(1);

(*(argStruct_p->W_fom_p)).elem(find((*(argStruct_p->W_fom_p))<=0)).fill(0.00000001);

}

void Plot_significant_L_wrapper(mxArray *plhs[], const mxArray *prhs[], arg_struct_t* argStruct_p){

/*Xtilde_fnm, 		prhs[0]*/
/*engPutVariable(mlEngine_p, "Xtilde_fnm", prhs[0]);*/

/*Z_ol, 			plhs[2]*/	
armaSetPr(plhs[2], (*(argStruct_p->Z_ol_p)));
engPutVariable(mlEngine_p, "Z_ol", plhs[2]);

/*Y_lk, 			plhs[3]*/	
armaSetPr(plhs[3], (*(argStruct_p->Y_lk_p)));
engPutVariable(mlEngine_p, "Y_lk", plhs[3]);

/*T_fk, 			plhs[4]*/	
armaSetPr(plhs[4], (*(argStruct_p->T_fk_p)));
engPutVariable(mlEngine_p, "T_fk", plhs[4]);

/*V_nk, 			plhs[5]*/	
armaSetPr(plhs[5], (*(argStruct_p->V_nk_p)));
engPutVariable(mlEngine_p, "V_nk", plhs[5]);

/*W_fom, 		plhs[1]*/
armaSetCubePr(plhs[1], (*(argStruct_p->W_fom_p)));
engPutVariable(mlEngine_p, "W_fom", plhs[1]);

/*expj_Phi_W_fom, 	prhs[18]*/
engPutVariable(mlEngine_p, "expj_Phi_W_fom", prhs[18]);

/*expj_Phi_S_fkn		plhs[7]*/
armaSetCubeCx(plhs[7], *(argStruct_p->expj_Phi_S_fkn_p));
engPutVariable(mlEngine_p, "expj_Phi_S_fkn", plhs[7]);

/*example code*/
/*armaSetPr(plhs[5], (*(argStruct_p->V_nk_p)));
engPutVariable(mlEngine_p, "V_nk", plhs[5]);*/

engEvalString(mlEngine_p, "plot_significant_l_clusters_2(Z_ol, Y_lk, T_fk, V_nk, W_fom, expj_Phi_W_fom, expj_Phi_S_fkn);");


}

void Plot_significant_L_wrapper_3(mxArray *plhs[], const mxArray *prhs[], arg_struct_t* argStruct_p){

/*Xtilde_fnm, 		prhs[0]*/
/*engPutVariable(mlEngine_p, "Xtilde_fnm", prhs[0]);*/

/*Z_ol, 			plhs[2]*/	
armaSetPr(plhs[2], (*(argStruct_p->Z_ol_p)));
engPutVariable(mlEngine_p, "Z_ol", plhs[2]);

/*Y_lk, 			plhs[3]*/	
armaSetPr(plhs[3], (*(argStruct_p->Y_lk_p)));
engPutVariable(mlEngine_p, "Y_lk", plhs[3]);

/*T_fk, 			plhs[4]*/	
armaSetPr(plhs[4], (*(argStruct_p->T_fk_p)));
engPutVariable(mlEngine_p, "T_fk", plhs[4]);

/*V_nk, 			plhs[5]*/	
armaSetPr(plhs[5], (*(argStruct_p->V_nk_p)));
engPutVariable(mlEngine_p, "V_nk", plhs[5]);

/*W_fom, 		plhs[1]*/
armaSetCubePr(plhs[1], (*(argStruct_p->W_fom_p)));
engPutVariable(mlEngine_p, "W_fom", plhs[1]);

/*expj_Phi_W_fom, 	prhs[18]*/
engPutVariable(mlEngine_p, "expj_Phi_W_fom", prhs[18]);

/*expj_Phi_S_fkn		plhs[7]*/
armaSetCubeCx(plhs[7], *(argStruct_p->expj_Phi_S_fkn_p));
engPutVariable(mlEngine_p, "expj_Phi_S_fkn", plhs[7]);

engPutVariable(mlEngine_p, "Xtilde_fnm", prhs[0]);

/*example code*/
/*armaSetPr(plhs[5], (*(argStruct_p->V_nk_p)));
engPutVariable(mlEngine_p, "V_nk", plhs[5]);*/

engEvalString(mlEngine_p, "plot_significant_l_clusters_3(Z_ol, Y_lk, T_fk, V_nk, W_fom, expj_Phi_W_fom, expj_Phi_S_fkn, Xtilde_fnm);");


}