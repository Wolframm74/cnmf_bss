#include "local_inc.hpp"

#define SHIFTMAT_SCALING_FACTOR 1

mat::fixed<N_static, K_static> ones_mat_NxK;

/*0th sum: not really related to the 0th sum, but compute the shift matrix*/
mat::fixed<N_static, K_static> shiftmat_real_NK;

/*2nd sum stuff*/
cube::fixed<N_static, K_static, M_static> outtensor_2nd_real_NKM;
cx_cube::fixed<N_static, K_static, M_static> outtensor_2nd_cx_NKM;

mat::fixed<F_static, K_static> dummy_mat_real_FK[NUM_WORKER_THREADS];
cx_mat::fixed<F_static, K_static> dummy_mat_cx_FK[NUM_WORKER_THREADS];

/*3rd sum stuff*/
mat::fixed<N_static, K_static> outmat_3rd_real_NK_den;
mat::fixed<N_static, K_static> outmat_3rd_real_NK_num;
cx_mat::fixed<N_static, K_static> outmat_3rd_cx_NK;

void V_0th_sum_preprocess(mat* V_nk_p){

int n_index;

for (n_index=N_static-1; n_index>0; n_index--){

shiftmat_real_NK.row(n_index)=(*V_nk_p).row(n_index-1);

}


}

static void compute_2nd_sum(cube* Xhat_low_fnm_p, cx_cube* E_conj_fnm_p, mat* T_fk_p, cx_cube* expj_Phi_S_fkn_p, int m_index, int n_index, int thread_iter){

#ifdef DEBUG	
mexPrintf("V_update: pthread_self(): %d, compute_2nd_sum(), m_index=%d, n_index=%d \n", (int)pthread_self(), m_index, n_index);		
#endif	

/*populate intermediate matrices*/
dummy_mat_real_FK[thread_iter]=Xhat_outtensor_real_FKM.slice(m_index)%(*T_fk_p);	/*needs to be a function of a thread_iter*/

dummy_mat_cx_FK[thread_iter]=Xhat_outtensor_cx_FKM.slice(m_index)%(*T_fk_p)%(*expj_Phi_S_fkn_p).slice(n_index);	/*needs to be a function of a thread_iter too!!!!*/

/*outtensor_2nd_real_NKM: 1xK=1xF*FxK*/
outtensor_2nd_real_NKM.slice(m_index).row(n_index)=trans((*Xhat_low_fnm_p).slice(m_index).col(n_index))*dummy_mat_real_FK[thread_iter];

/*outtensor_2nd_cx_NKM*/	
outtensor_2nd_cx_NKM.slice(m_index).row(n_index)=strans((*E_conj_fnm_p).slice(m_index).col(n_index))*dummy_mat_cx_FK[thread_iter];

}

void orthogonalize_Vnk_wrapper(mat* V_nk_p, mxArray *V_nk_mxArray_p){

mxArray* output_arg_p[1];

armaSetPr(V_nk_mxArray_p, *V_nk_p);

mexCallMATLAB(1, &output_arg_p[0], 1, &V_nk_mxArray_p, "orthogonalize_V_nk");

(*V_nk_p).set_real(armaGetPr(output_arg_p[0], true)); 	

}

static void compute_3rd_sum_last_thread(mat* V_nk_p, mat* compensation_mat_p){

double gamma=0.5;
double lamda1=0.5;
double lamda2=0.25;
double p_value=0.5;

mat mean_V=zeros<mat>(1);

/*The function calculates the nonnegative ratios. Should be necessary to zero the dummy/memory matrices at least once as follows: */
/*Only on thread should end up calling this function. Not a bad thing to zero in case you didn't zero at least once already. */	
outmat_3rd_real_NK_den.zeros();
outmat_3rd_cx_NK.zeros();

/*Integrate out m*/
for (int m_iter=0; m_iter<M_static; m_iter++){	
outmat_3rd_real_NK_den=outmat_3rd_real_NK_den+outtensor_2nd_real_NKM.slice(m_iter);

outmat_3rd_cx_NK=outmat_3rd_cx_NK+outtensor_2nd_cx_NKM.slice(m_iter);

}

outmat_3rd_real_NK_num=outmat_3rd_real_NK_den+real(outmat_3rd_cx_NK);

/*outmat_3rd_real_NK_den=outmat_3rd_real_NK_den+(lamda2)*((*V_nk_p)-shiftmat_real_NK);*/

outmat_3rd_real_NK_den=outmat_3rd_real_NK_den+(lamda1)*p_value*pow((*V_nk_p), p_value-1);

(*V_nk_p)=(*V_nk_p)+(1/gamma)*(*V_nk_p)-((1/gamma)*(*V_nk_p))%(outmat_3rd_real_NK_den/outmat_3rd_real_NK_num);

/*(*V_nk_p)=(*V_nk_p)%(outmat_3rd_real_NK_num/outmat_3rd_real_NK_den);*/

(*V_nk_p).elem(find((*V_nk_p)<=0)).fill(0.00000001);

(*V_nk_p).elem(find((*V_nk_p)>50)).fill(50);

}

static void V_2nd_sum_start(arg_struct_t* argStruct_p, int thread_iter){

int m_index;
int n_index;

bool last_element_queue_flag=false;

bool compute_flag; 

while (threads_while_condition_2nd_sum_flag&&(!(last_thread_to_sleep_flag_array[0][thread_iter]))){

	if ((!last_element_found_global_flag)&&(check_queue_global_flag)){

	/*Semaphore or mutex lock */
	sem_wait(&queue_sem);

	compute_flag=calculate_pair_queue_indices(argStruct_p, &m_index, &n_index, &last_element_queue_flag); 

	/*Semaphore or mutex unlock */
	sem_post(&queue_sem);

	}

	if (compute_flag){

	compute_2nd_sum(argStruct_p->Xhat_low_fnm_p, argStruct_p->E_conj_fnm_p, argStruct_p->T_fk_p, argStruct_p->expj_Phi_S_fkn_p, m_index, n_index, thread_iter);

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

void* V_start(void* arg){

int thread_iter;

thread_arg_t* threadArg_p;
arg_struct_t* argStruct_p;

threadArg_p=(thread_arg_t*)arg;
argStruct_p=threadArg_p->argStruct_p;
thread_iter=threadArg_p->thread_iter; 

V_2nd_sum_start(argStruct_p, thread_iter);

if (last_thread_to_sleep_flag_array[0][thread_iter]){

compute_3rd_sum_last_thread(argStruct_p->V_nk_p, argStruct_p->Y_lk_p);

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

if (!last_thread_to_sleep_flag_array[0][thread_iter]){

sem_post(&checkpoint_sem);

}

threads_exit_and_signal(NUM_WORKER_THREADS);
     
}

/*[0]
[1]
[2]
[3]
[4]
[5] window center
[6]
[7]
[8]
[9]
[10]*/


bool overlap1;
bool overlap2;
bool overlap3;


#define HALF_WINDOW_LEN 5
#define WINDOW_LEN (HALF_WINDOW_LEN*2)+1

rowvec::fixed<WINDOW_LEN> storage_vec;	/*storage for the row vector to be dotted with the window*/
				/*Indexing the correct "center" sample is key here. */


rowvec::fixed<WINDOW_LEN> window_vec;	/*odd numbered filter*/

rowvec::fixed<WINDOW_LEN> window_dummy;


rowvec::fixed<N_static> outvec_1xN;

/*Indices that index the vector to copy from (V_kn for some row index k) */
int storage_index1;
int storage_index2;

/*Indices tha index the vector to be copied to (*/
int window_index1;
int window_index2;	

void populate_window_vec(void){

int i;

double DC_gain;
double x;

for (i=0; i<WINDOW_LEN; i++){

x=(double)(i-HALF_WINDOW_LEN);

window_vec(i)=(1/(((double)HALF_WINDOW_LEN)*sqrt(2*(datum::pi))))*exp(-(pow(x,2))/(2*((double)HALF_WINDOW_LEN)));

}

DC_gain=sum(window_vec);

window_vec=(1/DC_gain)*window_vec;

}

static void V_filter_dot_op(int n_index){

outvec_1xN(n_index)=dot(window_dummy, storage_vec);

}

static void V_filter_row(rowvec& V_kn_kth_row){

colvec::fixed<1> scalar_local;
int n_index;

window_dummy.zeros();
storage_vec.zeros();

/*for index=0:N_static-1*/
for (n_index=0; n_index<N_static; n_index++){

/*	scalar_local=(double)n_index;

	scalar_local.print("V_filter_row: n_index=");*/

	if (n_index<HALF_WINDOW_LEN){

		storage_index1=0;
		storage_index2=(WINDOW_LEN-1)-(HALF_WINDOW_LEN-n_index);
		window_index1=(HALF_WINDOW_LEN-n_index);
		window_index2=WINDOW_LEN-1;

		storage_vec(span(storage_index1,storage_index2))=V_kn_kth_row(span(storage_index1,storage_index2));
		window_dummy(span(storage_index1,storage_index2))=window_vec(span(window_index1,window_index2));

	} else if (n_index<(N_static-HALF_WINDOW_LEN)) {

		storage_index1=n_index-HALF_WINDOW_LEN;
		storage_index2=n_index+HALF_WINDOW_LEN;
		window_index1=0;
		window_index2=WINDOW_LEN-1;

		storage_vec(span(window_index1,window_index2))=V_kn_kth_row(span(storage_index1,storage_index2));
		window_dummy=window_vec;


	} else {	/*assume n_index >= (N_static-WINDOW_LEN) */

		window_dummy.zeros();
		storage_vec.zeros();

		storage_index1=n_index-HALF_WINDOW_LEN;

/*		scalar_local=(double)storage_index1;
		scalar_local.print("storage_index1");*/

		storage_index2=N_static-1;

/*		scalar_local=(double)storage_index2;
		scalar_local.print("storage_index2");*/

		window_index1=0;

/*		scalar_local=(double)window_index1;
		scalar_local.print("window_index1");*/

		window_index2=(N_static-1)-storage_index1;	

/*		scalar_local=(double)window_index2;
		scalar_local.print("window_index2");		*/

		storage_vec(span(window_index1,window_index2))=V_kn_kth_row(span(storage_index1,storage_index2));
		window_dummy(span(window_index1,window_index2))=window_vec(span(window_index1,window_index2));

	}

	V_filter_dot_op(n_index);


}

}

void V_nk_MA_filter_rows(mat* V_nk_p){

int k_index;

rowvec::fixed<N_static> V_kn_kth_row;

for (k_index=0; k_index<K_static; k_index++){

	V_kn_kth_row=trans((*V_nk_p).col(k_index));

	V_filter_row(V_kn_kth_row);

	(*V_nk_p).col(k_index)=trans(outvec_1xN);

}

}


/*: Expected output:

n_index=0

	overlap1=true;

storage_index1=0;
storage_index2=10-5=5;
window_index1=5;
window_index2=10;	

window_dummy(0:5)=window_vec(5:10);

storage_vec(0:5)=V_kn_nth_row(0:5)

n_index=1

	overlap1=true;

storage_index1=0;
storage_index2=10-4=6;
window_index1=4;
window_index2=10;	

window_dummy(0:6)=window_vec(4:10);

storage_vec(0:6)=V_kn_nth_row(0:6)

n_index=4

	overlap1=true;

storage_index1=0;
storage_index2=9;
window_index1=1;
window_index2=10;	

window_dummy(0:9)=window_vec(1:10);

storage_vec(0:9)=V_kn_nth_row(0:9)

n_index=5

	overlap2=true;

storage_index1=0;
storage_index2=10;
window_index1=0;
window_index2=10;	

window_dummy(0:10)=window_vec(0:10);

storage_vec(0:10)=V_kn_nth_row(0:10)

n_index=6

	overlap2=true;

storage_index1=1;
storage_index2=11;
window_index1=0;
window_index2=10;	

window_dummy(0:10)=window_vec(0:10);

storage_vec(0:10)=V_kn_nth_row(1:11)

n_index=778

	overlap2=true;

storage_index1=773;
storage_index2=783;
window_index1=0;
window_index2=10;	

window_dummy(0:10)=window_vec(0:10);
storage_vec(0:10)=V_kn_nth_row(773:783)

n_index=779

	overlap3=true;

storage_index1=774;
storage_index2=783;
window_index1=0;
window_index2=9;	

window_dummy.zeros(), storage_vec.zeros();
window_dummy(0:9)=window_vec(0:9);
storage_vec(0:9)=V_kn_nth_row(774:783)

n_index=783

	overlap3=true;

storage_index1=777;
storage_index2=783;
window_index1=0;
window_index2=6;	

window_dummy.zeros(), storage_vec.zeros();
window_dummy(0:6)=window_vec(0:6);
storage_vec(0:6)=V_kn_nth_row(777:783)

n_index=783

	overlap3=true;

storage_index1=778;
storage_index2=783;
window_index1=0;
window_index2=5;	

window_dummy.zeros(), storage_vec.zeros();
window_dummy(0:5)=window_vec(0:5);
storage_vec(0:5)=V_kn_nth_row(778:783)

*/