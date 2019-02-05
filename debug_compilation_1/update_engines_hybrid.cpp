#include "local_inc.hpp"

bool update_W_fom_cx_flag; 	/*Xhat_update (or the update parsing) should check against this with an if statement of some sorts. If true it should then reset it to an off state*/	/*W_update should switch this from an off state to an on state. */

void update_engine(arg_struct_t* argStruct_p, int which_update){

/*Empty queue to copy from if necessary*/
/*std::queue<int> dummy_queue;*/   /*Make it a local variable on a function's stack to be destroyed every time the function returns */

/*Iterators*/
int thread_iter; 

/*Threads*/
pthread_t worker_threads[NUM_WORKER_THREADS];
/*pthread_t worker_threads_M[M_static];*/

/*Thread args*/
thread_arg_t thread_args[NUM_WORKER_THREADS]; 
/*thread_arg_t thread_args_M[M_static];*/

thread_arg_t* thread_arg_pz[NUM_WORKER_THREADS];
/*thread_arg_t* thread_arg_pz_M[M_static];*/

/*Dimensions related to how to populate the queue for the first time*/

int dim1_len=0;
int dim2_len=0;
bool dim_flag=false; 

int ind1;
int ind2;

/*Set up the function pointer*/

void *(*start_rtn)(void*);

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
      start_rtn=Xhat_start;

      dim1_len=argStruct_p->M;  
      dim2_len=argStruct_p->F;

      dim_flag=true; 

      strcpy(argStruct_p->update_string,"Xhat_start");

      /*Can disable this if you don't need it*/
      which_update_helper=which_update;

      break;
    case 2:
      start_rtn=T_start;
      dim1_len=argStruct_p->M;  
      dim2_len=argStruct_p->F;

      dim_flag=true; 

      strcpy(argStruct_p->update_string,"T_start");      

      break;
    case 3:
      start_rtn=V_start;
      dim1_len=argStruct_p->M; 
      dim2_len=argStruct_p->N;

      dim_flag=true; 

      strcpy(argStruct_p->update_string,"V_start");      

      break;                  
    case 4:
      start_rtn=W_start;
      dim1_len=argStruct_p->F; 

      strcpy(argStruct_p->update_string,"W_start");      

      accum_Z_ol=accu(*(argStruct_p->Z_ol_p));

      break;
    case 5:
      start_rtn=Phi_S_start;
      dim1_len=argStruct_p->N; 

      strcpy(argStruct_p->update_string,"Phi_S_start");      

      break;
    case 6:

      start_rtn=Z_start;
      dim1_len=argStruct_p->F;

      strcpy(argStruct_p->update_string,"Z_start");      

      break;
    case 7:

      start_rtn=Y_start;
      dim1_len=argStruct_p->M;
      dim2_len=argStruct_p->F;

      dim_flag=true;

      strcpy(argStruct_p->update_string,"Y_start");      

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

  for (ind2=0; ind2<NUM_WORKER_THREADS; ind2++){

    last_thread_to_sleep_flag_array[ind1][ind2]=false;

  }

}

/*Start worker threads */
for (thread_iter=0; thread_iter<NUM_WORKER_THREADS; thread_iter++){

	thread_args[thread_iter].thread_iter=thread_iter;

	thread_arg_pz[thread_iter]=&thread_args[thread_iter];

  thread_arg_pz[thread_iter]->argStruct_p=argStruct_p;

	/*pthread_create(&worker_threads[thread_iter], NULL, start_rtn, (void*) argStruct_p);*/
  pthread_create(&worker_threads[thread_iter], NULL, start_rtn, (void*) thread_arg_pz[thread_iter]);

}

/*MAIN THREAD SLEEP!!!!*/
sem_wait(&threads_exit_sleep_sem);

/*Join the worker threads */

for (thread_iter=0; thread_iter<NUM_WORKER_THREADS; thread_iter++){

	pthread_join(worker_threads[thread_iter], NULL);

}

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


void update_engine_M_threads(arg_struct_t* argStruct_p){

/*Empty queue to copy from if necessary*/
/*std::queue<int> dummy_queue;*/   /*Make it a local variable on a function's stack to be destroyed every time the function returns */

int ind1;
int ind2;

/*Iterators*/
int thread_iter; 

/*Threads*/
/*pthread_t worker_threads[NUM_WORKER_THREADS];*/
pthread_t worker_threads_M[M_static];

/*Thread args*/
/*thread_arg_t thread_args[NUM_WORKER_THREADS]; */
thread_arg_t thread_args_M[M_static];

/*thread_arg_t* thread_arg_pz[NUM_WORKER_THREADS];*/
thread_arg_t* thread_arg_pz_M[M_static];

/*Dimensions related to how to populate the queue for the first time*/

int dim1_len=0;
int dim2_len=0;
bool dim_flag=false; 

/*Set up the function pointer*/

void *(*start_rtn)(void*);

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


      start_rtn=Xhat_start_M_threads;
      dim1_len=argStruct_p->M; 
      strcpy(argStruct_p->update_string,"Xhat_start_M_threads");

/*Set this*/
argStruct_p->which_update=1;

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
/*threads_while_condition_3rd_sum_flag=true;*/

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
for (ind1=0; ind1<2; ind1++){

  for (ind2=0; ind2<M_static; ind2++){

    last_thread_to_sleep_flag_local[ind1][ind2]=false;

  }

}

      /*Start worker threads */
      for (thread_iter=0; thread_iter<M_static; thread_iter++){

        thread_args_M[thread_iter].thread_iter=thread_iter;

        thread_arg_pz_M[thread_iter]=&thread_args_M[thread_iter];

        thread_arg_pz_M[thread_iter]->argStruct_p=argStruct_p;

        pthread_create(&worker_threads_M[thread_iter], NULL, start_rtn, (void*) thread_arg_pz_M[thread_iter]);

      }

      /*MAIN THREAD SLEEP!!!!*/
      sem_wait(&threads_exit_sleep_sem);

      /*Join the worker threads */

      for (thread_iter=0; thread_iter<M_static; thread_iter++){

        pthread_join(worker_threads_M[thread_iter], NULL);

      }

      /*Reset stuff here:*/

/*Think/verify/double check what needs to be reset cleanly here for a clean start to Xhat_start()*/

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