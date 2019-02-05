#include "local_inc.hpp"

/*0th sum stuff*/
cube::fixed<F_static, N_static, L_static> outcube_0th_sum_FNL_real; 		/*Call reshape on the kronecker product output matrix to populate this tensor */	
mat::fixed<F_static, N_static> kronmat_0th_sum_dummy_real[NUM_WORKER_THREADS];			/*This is FxN intermediate result, 0th sum */
mat::fixed<F_static, N_static*L_static> kronmat_fat_0th_sum_dummy_real[NUM_WORKER_THREADS];		/*This is an Fx(N*L) intermediate result, 0th sum*/
mat::fixed<F_static, N_static*L_static> kronmat_fat_0th_sum_local_real[NUM_WORKER_THREADS];		/*The thread_iter'th array element is an Fx(N*L) matrix that stores the result of the accumulation for a certain thread. */
mat::fixed<F_static, N_static*L_static> kronmat_fat_0th_sum_output_real; 						/*kronecker product output matrix. All the local results need to be added to this output matrix. */		

cx_cube::fixed<F_static, N_static, L_static> outcube_0th_sum_FNL_cx; /*[NUM_WORKER_THREADS-1]; */
mat::fixed<F_static, N_static> kronmat_0th_sum_dummy_cx[NUM_WORKER_THREADS];			/*This is FxN intermediate result, 0th sum */
mat::fixed<F_static, N_static*L_static> kronmat_fat_0th_sum_dummy_cx[NUM_WORKER_THREADS];		/*This is an Fx(N*L) intermediate result, 0th sum*/
mat::fixed<F_static, N_static*L_static> kronmat_fat_0th_sum_local_cx[NUM_WORKER_THREADS];		/*The thread_iter'th array element is an Fx(N*L) matrix that stores the result of the accumulation for a certain thread. */
mat::fixed<F_static, N_static*L_static> kronmat_fat_0th_sum_output_cx; 						/*kronecker product output matrix. All the local results need to be added to this output matrix. */		

/*1st sum stuff*/


/*2nd sum stuff*/

cube::fixed<F_static, M_static, L_static> dummy_real_FML[NUM_WORKER_THREADS-1];
cube::fixed<M_static, O_static, L_static> dummy_real_MOL[NUM_WORKER_THREADS-1];
mat::fixed<O_static, L_static> dummy_real_OL[NUM_WORKER_THREADS-1];					/*m sum output (real) */

cx_cube::fixed<F_static, M_static, L_static> dummy_cx_FML[NUM_WORKER_THREADS-1];
cx_cube::fixed<M_static, O_static, L_static> dummy_cx_MOL[NUM_WORKER_THREADS-1];
cx_mat::fixed<O_static, L_static> dummy_cx_OL[NUM_WORKER_THREADS-1];				/*m sum output (cx) */

/*local sums*/


/*Each thread accumulates its local sum/tensors. Once there is no work left on the work queue, each thread can add its accumulated local sum to the global sum/tensor */


/*global sum(s)*/

/*queue elements are just the k index*/
/*static void compute_Z_ol(cx_cube* Xhat_fnm_p, cx_cube* Xtilde_fnm_p, cube* Xhat_low_fnm_p, cube* W_fom_p, mat* Z_ol_p, mat* Y_lk_p, mat* T_fk_p, mat* V_nk_p, cube* Phi_W_fom_p, cube* Phi_S_fnk_p, int l, int F, int N, int M, int K, int O){*/

/*Pass in pointers to Phase_W_fom and Phase_S_fnk as opposed to pointers to Phi_W_fom and Phi_S_fnk to whatever functions require them */

static void compute_0th_sum(mat* Y_lk_p, mat* T_fk_p, mat* V_nk_p, cx_cube* Phase_S_fnk_p, k_index, thread_iter){

/*This (and any other) static function(s)'s only responsibilty is to accumulate the dummy kronmat to the local kronmat. Z_start can accumulate the local kronmat to the output kronmat based on some logical states as seen by the calling function(s) */

/*Compute the FxN outer product component and write the result to kronmat_0th_sum_dummy_real[thread_iter] */
kronmat_0th_sum_dummy_real[thread_iter]=kron( trans((*V_nk_p).col(k_index)) , (*T_fk_p).col(k_index) );

/*Take the kth slice of Phase_S_fnk and do elementwise muliplication on the kth FxN component (outer product). Write the result to kronmat_0th_sum_dummy_cx[thread_iter]*/
kronmat_0th_sum_dummy_cx[thread_iter]=kronmat_0th_sum_dummy_real[thread_iter]%((*Phase_S_fnk).subcube(0,0,k_index, F_static-1,N_static-1,k_index));

/*Add kronmat_fat_dummy to kronmat_fat_local for the real and cx cases */	
kronmat_fat_0th_sum_local_real=kronmat_fat_0th_sum_local_real+kronmat_fat_0th_sum_dummy_real;

kronmat_fat_0th_sum_local_cx=kronmat_fat_0th_sum_local_cx+kronmat_fat_0th_sum_dummy_cx;

/*Expect that by some other mechanism (controlled by calling functions or otherwise) that the local kronmat_fat's will be added to the output kronmat_fat, and reshaping and assignment to outcube_0th_sum_FNL_real and outcube_0th_sum_FNL_cx will occur */	
/*ie: everything regarding the 0th sum is prepped prior to compute_1st_sum() being called */	

}

static void compute_1st_sum(){


}

static void compute_2nd_sum(){


}

/*With the exception of the W_fom update, the computation of the 3rd sum involves only M components, and therefore no need to set up a work queue for the computations for this*/
static void compute_3rd_sum(){

}


static void normalize_across_an_index(){


}


void* Z_start(void* arg){

arg_struct_t* argStruct_p;

/*Local variables*/
/*int l_index;*/
int l_index;

bool b_exit_cond_local=false; 
bool b_update_last_queue_element=false;

/*Cast the void* arg into a arg_struct_t**/
argStruct_p=(arg_struct_t*)arg;

/*int L=(int)argStruct_p->L;
int K=(int)argStruct_p->K;*/

while((!b_exit_cond_global)&&(!b_exit_cond_local)){

	/*Semaphore or mutex lock */
	sem_wait(&queue_sem);

	/*If the size is greater than 1 */
	if((argStruct_p->queue_intz.size())>1){

		/*Read from the front element of the queue */	
		l_index=argStruct_p->queue_intz.front();

		/*Pop off the front element - you don't need it anymore */
		argStruct_p->queue_intz.pop();

	} else if((argStruct_p->queue_intz.size())==1){	/*Endif - else if queue size == 1 : */

		/*Set the b_update_last_queue_element=true; */
		b_update_last_queue_element=true;

		/*Read from the front element of the queue */	
		l_index=argStruct_p->queue_intz.front();

		/*Pop off the front element - you don't need it anymore */
		argStruct_p->queue_intz.pop();

	} else if ((argStruct_p->queue_intz.size())<=0)	{ /*End else if */

	/*Endif - else if : */

		/*Exit the while loop? set the local exit condition.  */
		b_exit_cond_local=true; 

	}

	/*End else if */

	/*Semaphore or mutex unlock */
	sem_post(&queue_sem);

	/*if (b_exit_cond_g||b_update_last_queue_element*/	/*also make this a function of b_exit_cond_local */
	if (((!b_exit_cond_global)&&(!b_exit_cond_local))){

	/*Update should be outside the sem-lock/unlock so that other threads don't have to wait until the entireee update is complete to access the queue. */
	compute_Z_ol(argStruct_p->Xhat_fnm_p, argStruct_p->Xtilde_fnm_p, argStruct_p->Xhat_low_fnm_p, argStruct_p->W_fom_p, argStruct_p->Z_ol_p, argStruct_p->Y_lk_p, argStruct_p->T_fk_p, argStruct_p->V_nk_p, argStruct_p->Phi_W_fom_p, argStruct_p->Phi_S_fnk_p, l_index, argStruct_p->F, argStruct_p->N, argStruct_p->M, argStruct_p->K, argStruct_p->O);
	/*compute_V_nkl(argStruct_p->H_flm_p, argStruct_p->T_fkl_p, argStruct_p->V_nkl_p, argStruct_p->Phi_H_flm_p, k_index, argStruct_p->F, argStruct_p->N, argStruct_p->M, argStruct_p->L);*/
	/*argStruct_p->Phi_S_fnkl_arr, argStruct_p->C_flkmn_arr, argStruct_p->B_flkmn_arr,*/

		/*if b_update_last_queue_element*/
		if (b_update_last_queue_element==true){

		/*Set b_exit_cond_global */
		b_exit_cond_global=true;

		}

	}
	/*end if*/

}
/*end while*/
	
}