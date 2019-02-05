#include "local_inc.hpp"

/*Shared variables*/

int complex_argument_costfun_m4_m_index_b;
int complex_argument_costfun_m4_m_index_a;

cx_mat::fixed<F_static, N_static> complex_argument_costfun_m4_Error_mat_FN;
mat::fixed<F_static, N_static> complex_argument_costfun_m4_Magnitude_model_mat_FN;

/*End// shared variables*/

void update_engine_complex_argument_costfun_m4(arg_struct_t* argStruct_p, int which_update, int m_index_b, int m_index_a){

/*Empty queue to copy from if necessary*/
/*std::queue<int> dummy_queue;*/   /*Make it a local variable on a function's stack to be destroyed every time the function returns */

/*Iterators*/
int thread_iter; 

/*Threads*/
pthread_t worker_threads[NUM_WORKER_THREADS_COMPLEX_ARGUMENT_COSTFUN_M4];
/*pthread_t worker_threads_M[M_static];*/

/*Thread args*/
thread_arg_t thread_args[NUM_WORKER_THREADS_COMPLEX_ARGUMENT_COSTFUN_M4]; 
/*thread_arg_t thread_args_M[M_static];*/

thread_arg_t* thread_arg_pz[NUM_WORKER_THREADS_COMPLEX_ARGUMENT_COSTFUN_M4];
/*thread_arg_t* thread_arg_pz_M[M_static];*/

/*Dimensions related to how to populate the queue for the first time*/

int dim1_len=0;
int dim2_len=0;
bool dim_flag=false; 

int ind1;
int ind2;

/*Set up the function pointer*/

void *(*start_rtn)(void*);

/*Set the indices!*/
complex_argument_costfun_m4_m_index_b=m_index_b;
complex_argument_costfun_m4_m_index_a=m_index_a;

/*Compute the Error matrix*/

complex_argument_costfun_m4_Error_mat_FN=((*(argStruct_p->Xtilde_fnm_p)).slice(m_index_b)%conj((*(argStruct_p->Xtilde_fnm_p)).slice(m_index_a)))-((*(argStruct_p->Xhat_fnm_p)).slice(m_index_b)%conj((*(argStruct_p->Xhat_fnm_p)).slice(m_index_a)));

complex_argument_costfun_m4_Magnitude_model_mat_FN=((*(argStruct_p->Xhat_low_fnm_p)).slice(m_index_b))%((*(argStruct_p->Xhat_low_fnm_p)).slice(m_index_a));

/*Init the relevant semaphores*/
sem_init(&check_queue_global_sem, 1, 1);
sem_init(&threads_asleep_sem, 1, 1);
sem_init(&queue_sem, 1, 1);
sem_init(&sleep_sem, 1, 0);
sem_init(&checkpoint_sem, 1, 0);
sem_init(&wait_till_queue_populated_sem, 1, 0);
/*sem_init(&last_element_sleep_sem, 1, 0);*/
sem_init(&threads_while_condition_0th_sum_sem, 1, 1);
sem_init(&threads_while_condition_1st_sum_sem, 1, 1);
sem_init(&threads_while_condition_2nd_sum_sem, 1, 1);
sem_init(&work_actually_computed_sem, 1, 1);
sem_init(&threads_exit_lock_sem, 1, 1);
sem_init(&threads_exit_sleep_sem, 1, 0);

/*Parse which update variable*/

switch (which_update){

    case 1:
      start_rtn=T_complex_argument_m4_start;
      dim1_len=(double)F_static;
      dim2_len=(double)K_static;

      dim_flag=true;

      strcpy(argStruct_p->update_string,"T_complex_argument_m4_start");      

      break;
    case 2:
      start_rtn=V_complex_argument_m4_start;
      dim1_len=(double)N_static;
      dim2_len=(double)K_static;

      dim_flag=true;

      strcpy(argStruct_p->update_string,"V_complex_argument_m4_start");      

      break;                  
    case 3:
	  start_rtn=W_channel_a_complex_argument_m4_start;
      dim1_len=(double)F_static;
      dim2_len=(double)O_static;

      dim_flag=true;

      strcpy(argStruct_p->update_string,"W_channel_a_complex_argument_m4_start");      

      /*accum_Z_ol=accu(*(argStruct_p->Z_ol_p));*/

      break;      
    case 4:
      start_rtn=W_channel_b_complex_argument_m4_start;
      dim1_len=(double)F_static;
      dim2_len=(double)O_static;

      dim_flag=true;

      strcpy(argStruct_p->update_string,"W_channel_b_complex_argument_m4_start");      

      /*accum_Z_ol=accu(*(argStruct_p->Z_ol_p));*/

      break;
    case 5:
      start_rtn=Phi_S_complex_argument_m4_start;
      dim1_len=(double)F_static;
      dim2_len=(double)N_static;

      dim_flag=true;

      strcpy(argStruct_p->update_string,"Phi_S_complex_argument_m4_start");      

      break;
    case 6:

      /*Rotating function defined in complex_argument_costfun_m1.hpp*/
      rotate_expj_Phi_S_nkf_p_to_FNK(argStruct_p->expj_Phi_S_nkf_p);

      /*Can now safely access an up-to-date rotated_expj_Phi_S_FNK, from this point on*/

      start_rtn=Z_complex_argument_m4_start;
      dim1_len=(double)O_static;
      dim2_len=(double)L_static;

      dim_flag=true;
      strcpy(argStruct_p->update_string,"Z_complex_argument_m4_start");      

      break;
    case 7:

      rotate_expj_Phi_S_nkf_p_to_FNK(argStruct_p->expj_Phi_S_nkf_p);

      /*Can now safely access an up-to-date rotated_expj_Phi_S_FNK, from this point on*/

      start_rtn=Y_complex_argument_m4_start;
      dim1_len=(double)L_static;
      dim2_len=(double)K_static;

      dim_flag=true;

      strcpy(argStruct_p->update_string,"Y_complex_argument_m4_start");      

      break;      
    default: 
      assert(0);
      break;	
}

mexPrintf("update_engines, line 110: pthread_self(): %d, which_update: %d\n", (int)pthread_self(), which_update);

/*Set this*/
argStruct_p->which_update=which_update;

if (argStruct_p->queue_intz.size()>0){
  clear_queue_intz(argStruct_p->queue_intz);
}

if (dim_flag==false){

	/*Init queue wrt single index. */
	populate_queue_wrt_one_index(argStruct_p, dim1_len);

}

if (dim_flag==true){

	/*Init queue wrt a pair of indices. */
	populate_queue_wrt_pair_indices(argStruct_p, dim1_len, dim2_len);

}

/*Enable all while flags*/
threads_while_condition_0th_sum_flag=true;
threads_while_condition_1st_sum_flag=true;
threads_while_condition_2nd_sum_flag=true; 

/*Init to 0*/
sem_wait(&threads_exit_lock_sem);
threads_exit_ctr=0;
sem_post(&threads_exit_lock_sem);

/*Init this to false*/
last_element_found_global_flag=false;

/*Init this to 0*/
threads_asleep_ctr=0;

/*Init this to true*/
sem_wait(&check_queue_global_sem);
check_queue_global_flag=true;
sem_post(&check_queue_global_sem);

/*Reset*/
work_actually_completed_flag=false;
work_actually_computed_ctr=0;

/*Reset this to all false for every entry*/
for (ind1=0; ind1<MAX_NUM_SUMS; ind1++){

  for (ind2=0; ind2<NUM_WORKER_THREADS_COMPLEX_ARGUMENT_COSTFUN_M4; ind2++){

    last_thread_to_sleep_flag_array[ind1][ind2]=false;

  }

}

/*Start worker threads */
for (thread_iter=0; thread_iter<NUM_WORKER_THREADS_COMPLEX_ARGUMENT_COSTFUN_M4; thread_iter++){

	thread_args[thread_iter].thread_iter=thread_iter;

	thread_arg_pz[thread_iter]=&thread_args[thread_iter];

  thread_arg_pz[thread_iter]->argStruct_p=argStruct_p;

	/*pthread_create(&worker_threads[thread_iter], NULL, start_rtn, (void*) argStruct_p);*/
  pthread_create(&worker_threads[thread_iter], NULL, start_rtn, (void*) thread_arg_pz[thread_iter]);

}

  mexPrintf("update_engines(): MAIN THREAD DONE CREATING THREAD. (LINE 202) CALLING SEM_WAIT() pthread_self(): %d, \n", (int)pthread_self());

/*MAIN THREAD SLEEP!!!!*/
sem_wait(&threads_exit_sleep_sem);

  mexPrintf("update_engines(): CREATING THREAD. (LINE 207) RETURNED FROM SEM_WAIT() pthread_self(): %d, \n", (int)pthread_self());

/*Join the worker threads */

for (thread_iter=0; thread_iter<NUM_WORKER_THREADS_COMPLEX_ARGUMENT_COSTFUN_M4; thread_iter++){

	pthread_join(worker_threads[thread_iter], NULL);

}

  mexPrintf("update_engines(): (LINE 217) RETURNED FROM ALL PTHREAD_JOIN()'S' pthread_self(): %d, \n", (int)pthread_self());

/*Destroy semaphores*/
sem_destroy(&queue_sem);
sem_destroy(&sleep_sem);
sem_destroy(&checkpoint_sem);
sem_destroy(&wait_till_queue_populated_sem);
/*sem_destroy(&last_element_sleep_sem);*/
sem_destroy(&check_queue_global_sem);
sem_destroy(&threads_asleep_sem);
sem_destroy(&threads_while_condition_0th_sum_sem);
sem_destroy(&threads_while_condition_1st_sum_sem);
sem_destroy(&threads_while_condition_2nd_sum_sem);
sem_destroy(&work_actually_computed_sem);
sem_destroy(&threads_exit_lock_sem);
sem_destroy(&threads_exit_sleep_sem);

}

