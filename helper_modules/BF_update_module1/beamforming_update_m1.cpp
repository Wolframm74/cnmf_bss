#include "local_inc.hpp"

#define DESIRED_NO_ITERATIONS_BEAMFORMING_KMEANS 50

static cx_cube::fixed<F_static, N_static, K_static> expj_Phi_S_fnk_local;
static cx_mat::fixed<N_static*K_static, F_static> entry_point1_fun1_mat1_NKxF;
static cx_mat::fixed<F_static, N_static*K_static> entry_point1_fun1_mat2_FxNK;

static void entry_point1_fun1_rotate_Phi_S(cx_cube* expj_Phi_S_nkf_p){

/*copy the nkf input into NKxF */
arrayops::copy(entry_point1_fun1_mat1_NKxF.memptr(), (*expj_Phi_S_nkf_p).memptr(), F_static*N_static*K_static);

/*transpose NKxF into FxNK*/	
entry_point1_fun1_mat2_FxNK=strans(entry_point1_fun1_mat1_NKxF);

/*copy FxNK into output FNK*/	
arrayops::copy(expj_Phi_S_fnk_local.memptr(), entry_point1_fun1_mat2_FxNK.memptr(), F_static*N_static*K_static);

}

static rowvec::fixed<K_static> frob_norm_rowvec_norm_V;

static void normalize_V_nk_local(mat* V_nk_p){

frob_norm_rowvec_norm_V=sqrt(sum( (*V_nk_p)%(*V_nk_p) , 0));

/*Divide out each col by its frob norm:*/
(*V_nk_p)=(*V_nk_p)/kron(ones_col_Nx1, frob_norm_rowvec_norm_V);

}

static mat::fixed<N_static, K_static> V_nk_complement; 

static void entry_point1_fun2_compute_V_complement(mat* V_nk_p, mat* T_fk_p){

mat* V_nk_complement_p=&V_nk_complement;

/*first normalize V*/
normalize_V_nk_frob_norm(V_nk_p, T_fk_p);

/*Set the complement to ones*/

V_nk_complement.ones();

/*Subtract off V from the complement*/
V_nk_complement=V_nk_complement-(*V_nk_p);

/*There should only be ones left in positions where V was originally zero. This is what you want*/

/*Normalize the complement to bring the average energy down a bit*/
normalize_V_nk_local(V_nk_complement_p);

}


static colvec::fixed<O_static> kith_dummy_sum_vec;
static mat::fixed<O_static, K_static> Z_ok_target_BF;
static mat::fixed<O_static, K_static> Z_ok_sum_matrix_BF;
static cx_mat::fixed<F_static, N_static> dummy_mat_FxN_entry_point1_fun3; 

void entry_point1_fun3_compute_Zko_target(mat* T_fk_p, mat* V_nk_p, cx_cube* W_fom_cx_p){

Z_ok_target_BF.zeros();

int k_i_index, k_index;

int m_index, f_index, n_index, o_index; 

double sum_min_o_index_value, avg_less_one;

uword min_o_index;

for (k_i_index=0; k_i_index<K_static; k_i_index++){

	kith_dummy_sum_vec.zeros();

	for (o_index=0; o_index<O_static; o_index++){


		for (m_index=0; m_index<M_static; m_index++){

					for (k_index=0; k_index<K_static; k_index++){

						if (k_index!=k_i_index){

						dummy_mat_FxN_entry_point1_fun3=expj_Phi_S_fnk_local.slice(k_index);

						dummy_mat_FxN_entry_point1_fun3=dummy_mat_FxN_entry_point1_fun3%kron( (*T_fk_p).col(k_index) , trans(ones_col_Nx1) );

						dummy_mat_FxN_entry_point1_fun3=dummy_mat_FxN_entry_point1_fun3%kron( (*W_fom_cx_p).slice(m_index).col(o_index) , trans(ones_col_Nx1) );

						dummy_mat_FxN_entry_point1_fun3=dummy_mat_FxN_entry_point1_fun3%kron( ones_col_Fx1, trans((*V_nk_p).col(k_index)) );

						dummy_mat_FxN_entry_point1_fun3=dummy_mat_FxN_entry_point1_fun3%kron( ones_col_Fx1, trans(V_nk_complement.col(k_index)) );

						kith_dummy_sum_vec(o_index)=kith_dummy_sum_vec(o_index)+real(accu(dummy_mat_FxN_entry_point1_fun3%conj(dummy_mat_FxN_entry_point1_fun3)));

						}

					}

		}

	}

	Z_ok_sum_matrix_BF.col(k_i_index)=kith_dummy_sum_vec;

	kith_dummy_sum_vec.min(min_o_index);

	sum_min_o_index_value=kith_dummy_sum_vec(min_o_index);

	kith_dummy_sum_vec(min_o_index)=0; 

	avg_less_one=as_scalar(sum(kith_dummy_sum_vec))/(((double)O_static)-1);

	/*Z_ok_target_BF((int)min_o_index, k_i_index)=1;*/

	Z_ok_target_BF((int)min_o_index, k_i_index)=(avg_less_one-sum_min_o_index_value)/avg_less_one;

}


}

static mat::fixed<O_static, L_static> chat_oq_estimate;
static colvec::fixed<L_static> kth_rand_value_for_qth_starting_cluster_center;

static void entry_point1_fun4_fun0_random_initialization(void){

int q_index;

kth_rand_value_for_qth_starting_cluster_center.randu();
kth_rand_value_for_qth_starting_cluster_center=(((double)K_static)-1)*kth_rand_value_for_qth_starting_cluster_center;
floor(kth_rand_value_for_qth_starting_cluster_center);

for (q_index=0; q_index<L_static; q_index++){

/*kth_rand_value_for_qth_starting_cluster_center(q_index)=randu(1, distr_param(0, K_static-1));*/

chat_oq_estimate.col(q_index)=Z_ok_target_BF.col((int)(kth_rand_value_for_qth_starting_cluster_center(q_index)))+0.1*randu<colvec>(O_static);

}

}

static mat::fixed<L_static, K_static> indicators_P_qk;

static mat::fixed<L_static, K_static> indicators_P_qk_frob_den_mat;

static rowvec::fixed<3> index_display_rowvec;
static rowvec::fixed<3> frob_display_rowvec_kq_kr_frobden;

#define STABILITY_FACTOR_KMEANS 10000000

static void entry_point1_fun4_fun1_compute_indicators(void){

int k_index, q_index, r_index;

double frob_norm_kq;
double frob_norm_kr;

double frob_den; 

Z_ok_target_BF.print("line 181, Z_ok_target_BF");
chat_oq_estimate.print("line 182, chat_oq_estimate");

for (q_index=0; q_index<L_static; q_index++){

	index_display_rowvec(0)=q_index;

	/*for (k_index=0; k_index<(L_static*F_static-L_static); k_index++){*/
	for (k_index=0; k_index<(K_static); k_index++){

	index_display_rowvec(1)=k_index;

	/*index_display_rowvec.print("display indices: q_index, k_index");*/

		frob_den=0;

		frob_norm_kq=sqrt(as_scalar(sum(square(abs( Z_ok_target_BF.col(k_index)-chat_oq_estimate.col(q_index) )))));
		/*frob_norm_kq=(as_scalar(sum(square(abs( Z_ok_target_BF.col(k_index)-chat_oq_estimate.col(q_index) )))));*/

		frob_display_rowvec_kq_kr_frobden(0)=frob_norm_kq;

		/*frob_display_rowvec_kq_kr_frobden.print("frob_display_rowvec_kq_kr_frobden, frob_norm_kq changed:");*/

		for (r_index=0; r_index<L_static; r_index++){

			index_display_rowvec(2)=r_index;

			/*index_display_rowvec.print("display indices: q_index, k_index, r_index");*/

			frob_norm_kr=sqrt(as_scalar(sum(square(abs( Z_ok_target_BF.col(k_index)-chat_oq_estimate.col(r_index) )))));
			/*frob_norm_kr=(as_scalar(sum(square(abs( Z_ok_target_BF.col(k_index)-chat_oq_estimate.col(r_index) )))));*/

			frob_display_rowvec_kq_kr_frobden(1)=frob_norm_kr;

			frob_den=frob_den+pow((frob_norm_kq/frob_norm_kr), 2);

			frob_display_rowvec_kq_kr_frobden(2)=frob_den;

			/*frob_display_rowvec_kq_kr_frobden.print("frob_display_rowvec_kq_kr_frobden, frob_norm_kr, frobden changed");*/

		}

		/*indicators_P_qk(q_index, k_index)=1/frob_den;*/

		indicators_P_qk(q_index, k_index)=1/frob_den;

		indicators_P_qk_frob_den_mat(q_index, k_index)=frob_den;

	}

}

/*indicators_P_qk.print("indicators_P_qk, line 234:");*/

}

static colvec::fixed<O_static> means_computation_dummy_Ox1_2;
static double means_computation_accu_double_dummy_2;

static void entry_point1_fun4_fun2_compute_means(void){

int q_index, k_index;

for (q_index=0; q_index<L_static; q_index++){

	means_computation_dummy_Ox1_2.zeros();
	means_computation_accu_double_dummy_2=0;

	for (k_index=0; k_index<K_static; k_index++){

	means_computation_dummy_Ox1_2=means_computation_dummy_Ox1_2+pow(indicators_P_qk(q_index, k_index), 2)*Z_ok_target_BF.col(k_index);

	means_computation_accu_double_dummy_2=means_computation_accu_double_dummy_2+pow(indicators_P_qk(q_index, k_index), 2);

	}

	chat_oq_estimate.col(q_index)=means_computation_dummy_Ox1_2/means_computation_accu_double_dummy_2;

}

}

static void entry_point1_fun4_compute_k_means(void){

entry_point1_fun4_fun0_random_initialization();

chat_oq_estimate.print("chat_oq_estimate, line 232:");

int iter_ctr;

/*Z_ok_target_BF=((double)STABILITY_FACTOR_KMEANS)*Z_ok_target_BF;

chat_oq_estimate=((double)STABILITY_FACTOR_KMEANS)*chat_oq_estimate;*/

for (iter_ctr=0; iter_ctr<DESIRED_NO_ITERATIONS_BEAMFORMING_KMEANS; iter_ctr++){

entry_point1_fun4_fun1_compute_indicators();

entry_point1_fun4_fun2_compute_means();

}

chat_oq_estimate.print("chat_oq_estimate, line 242");

}

rowvec::fixed<L_static> entry_point1_fun5_dummy_vec1;
colvec::fixed<NUM_PERMUTATIONS> entry_point1_fun5_dummy_vec2_sum_of_squares_vector;

static mat::fixed<O_static, L_static> chat_oq_estimate_permutation_aligned;
static mat::fixed<L_static, K_static> indicators_P_qk_permutation_aligned;

static void entry_point1_fun5_compute_minimum_sum_squares_permutation(mat* Z_ol_p){

int i_iter;

int l_index, q_index;

uword q_index_uword, min_permutation_output_index_uword;

int min_permutation_output_index;

double i_iterth_permutation_sum_of_squares_dummy;

/*chat_oq_estimate and rhat_ML*/
for (i_iter=0; i_iter<NUM_PERMUTATIONS; i_iter++){

	i_iterth_permutation_sum_of_squares_dummy=0; 

	for (l_index=0; l_index<L_static; l_index++){

	entry_point1_fun5_dummy_vec1=permutation_list.slice(i_iter).row(l_index);

	/*q_index=(int)entry_point1_fun5_dummy_vec1.index_max();*/
	entry_point1_fun5_dummy_vec1.max(q_index_uword);
	q_index=(int)q_index_uword;

	i_iterth_permutation_sum_of_squares_dummy=i_iterth_permutation_sum_of_squares_dummy+as_scalar(sum(square(abs((*(Z_ol_p)).col(l_index)-chat_oq_estimate.col(q_index)))));

	}

	entry_point1_fun5_dummy_vec2_sum_of_squares_vector(i_iter)=i_iterth_permutation_sum_of_squares_dummy;

}

/*min_permutation_output_index=entry_point1_fun5_dummy_vec2_sum_of_squares_vector.index_min();*/
entry_point1_fun5_dummy_vec2_sum_of_squares_vector.min(min_permutation_output_index_uword);
min_permutation_output_index=(int)min_permutation_output_index_uword;

chat_oq_estimate_permutation_aligned=chat_oq_estimate*permutation_list.slice(min_permutation_output_index);

indicators_P_qk_permutation_aligned=permutation_list.slice(min_permutation_output_index)*indicators_P_qk;

}

static colvec::fixed<L_static> frob_norm_colvec_norm_Y_local;

/*Frobenius norm normalize both Y_lk and P_qk */
static void frob_norm_normalize_Y_lk_local(mat* Y_lk_p){

int l_iter;

/*compute each frobenius norm in Y_lk for l=1:L*/

frob_norm_colvec_norm_Y_local=sqrt(sum((*Y_lk_p)%(*Y_lk_p), 1));

/*Divide out each row by its frob norm:*/
(*Y_lk_p)=(*Y_lk_p)/kron(frob_norm_colvec_norm_Y_local, ones_row_1xK);

}

static void entry_point1_fun6_compute_midpoints(mat* Z_ol_p, mat* Y_lk_p){

mat* Z_ol_target_p=&chat_oq_estimate_permutation_aligned; 

mat* Y_lk_target_p=&indicators_P_qk_permutation_aligned;

normalize_Z_ol_frob_norm(Z_ol_p);
normalize_Z_ol_frob_norm(Z_ol_target_p);

(*Z_ol_p)=(0.2*(*Z_ol_p)+0.8*(*Z_ol_target_p));

/*frob_norm_normalize_Y_lk_local(Y_lk_p);
frob_norm_normalize_Y_lk_local(Y_lk_target_p);

(*Y_lk_p)=(0.2*(*Y_lk_p)+0.8*(*Y_lk_target_p));*/

}

void Beamforming_update_module1_entry_point1(mat* Z_ol_p, mat* Y_lk_p, mat* T_fk_p, mat* V_nk_p, cx_cube* W_fom_cx_p, cx_cube* expj_Phi_S_nkf_p){

/*Rotate Phi_S*/
entry_point1_fun1_rotate_Phi_S(expj_Phi_S_nkf_p);

/*Rotate function might be useful if you wanted to create a code that used accu(FxN) instead of so many for loops as done above. */	

/*Compute the complement of V*/	
entry_point1_fun2_compute_V_complement(V_nk_p, T_fk_p);

/*Iterate for k=1:K, o=1:O*/	
entry_point1_fun3_compute_Zko_target(T_fk_p, V_nk_p, W_fom_cx_p);

Z_ok_sum_matrix_BF.print("Z_ok_sum_matrix_BF, line 367:");

/*This should compute a target matrix Z_ko target*/	

/*Integrate Z_ko target into the current model*/	
entry_point1_fun4_compute_k_means();

/*compute the permutation that minimizes the sum of squares*/
entry_point1_fun5_compute_minimum_sum_squares_permutation(Z_ol_p);

entry_point1_fun6_compute_midpoints(Z_ol_p, Y_lk_p);

}