void* Xhat_start_M_threads(void* arg){

thread_arg_t* threadArg_p;
arg_struct_t* argStruct_p;

bool last_element_flag=false;
bool work_not_yet_completed_flag=true;
bool note_to_sleep_flag=false;

threadArg_p=(thread_arg_t*)arg;
argStruct_p=threadArg_p->argStruct_p;

Xhat_0th_sum_start(argStruct_p, &last_element_flag);

if (last_element_flag){

while(work_not_yet_completed_flag||all_threads_not_yet_asleep){

	sem_wait(&work_actually_computed_sem);

	if (work_actually_completed_flag){

	/*set the while flag off so that you can freely proceed*/		
	work_not_yet_completed_flag=false;		
	note_to_sleep_flag=false;

	} else {

	/*send a note to yourself to sleep after you've posted*/
	note_to_sleep_flag=true;		

	}

	sem_post(&work_actually_computed_sem);

	sem_wait(&threads_asleep_sem);

	if (threads_asleep_ctr<M_static-1){

	note_to_sleep_flag=true;

	} else {

	/*All threads are somehow asleep*/		
	all_threads_not_yet_asleep=false;
	note_to_sleep_flag=false;

	}

	sem_post(&threads_asleep_sem);

	if(note_to_sleep_flag){

	/*sleep until the right thread for the job calls sem_post() in order to wake you!!!*/		
	sem_wait(&last_element_sleep_sem);

	}

}

}




}