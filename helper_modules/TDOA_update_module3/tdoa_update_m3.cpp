#include "local_inc.hpp"

#define DESIRED_NO_ITERATIONS_TDOAU_M3_KMEANS 50

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

static cx_cube::fixed<F_static, N_static, K_static> Y_est_FxNxK;
static cx_mat::fixed<F_static, N_static> entry_point1_fun2_mat1_FN;

static void entry_point1_fun2_compute_Y_est(mat* T_fk_p, mat* V_nk_p){

int k_index;

for (k_index=0; k_index<K_static; k_index++){

	Y_est_FxNxK.slice(k_index)=expj_Phi_S_fnk_local.slice(k_index)%kron((*T_fk_p).col(k_index), trans((*V_nk_p).col(k_index)));

}

}

static cx_cube::fixed<M_static, F_static, N_static> Xtilde_mfn_local;
static cx_mat::fixed<F_static*N_static, M_static> entry_point1_fun3_fun0_mat1_FNxM;
static cx_mat::fixed<M_static, F_static*N_static> entry_point1_fun3_fun0_mat2_MxFN;

static cx_cube::fixed<K_static, F_static, N_static> Y_est_KxFxN_local;
static cx_mat::fixed<F_static*N_static, K_static> entry_point1_fun3_fun0_mat3_FNxK;
static cx_mat::fixed<K_static, F_static*N_static> entry_point1_fun3_fun0_mat4_KxFN;


static void entry_point1_fun3_fun0_rotate_X_Y(cx_cube* Xtilde_fnm_p){

/*copy the nkf input into NKxF */
arrayops::copy(entry_point1_fun3_fun0_mat1_FNxM.memptr(), (*Xtilde_fnm_p).memptr(), F_static*N_static*M_static);

/*transpose NKxF into FxNK*/	
entry_point1_fun3_fun0_mat2_MxFN=strans(entry_point1_fun3_fun0_mat1_FNxM);

/*copy FxNK into output FNK*/	
arrayops::copy(Xtilde_mfn_local.memptr(), entry_point1_fun3_fun0_mat2_MxFN.memptr(), F_static*N_static*M_static);

/*copy the nkf input into NKxF */
arrayops::copy(entry_point1_fun3_fun0_mat3_FNxK.memptr(), Y_est_FxNxK.memptr(), F_static*N_static*K_static);

/*transpose NKxF into FxNK*/	
entry_point1_fun3_fun0_mat4_KxFN=strans(entry_point1_fun3_fun0_mat3_FNxK);

/*copy FxNK into output FNK*/	
arrayops::copy(Y_est_KxFxN_local.memptr(), entry_point1_fun3_fun0_mat4_KxFN.memptr(), F_static*N_static*K_static);


}

static cx_cube::fixed<M_static, K_static, F_static_admissible> A_pseudoinv_MKF; 
static cx_mat::fixed<M_static, K_static> entry_point1_fun3_fun1_mat1_MK;
static cx_mat::fixed<K_static, K_static> entry_point1_fun3_fun1_mat2_KK;

static void entry_point1_fun3_fun1_compute_A_mkf(void){

int f_index, n_index; 

int f_index_admissible;

for (f_index=0; f_index<F_static_admissible; f_index++){

	f_index_admissible=f_index+1; 

	entry_point1_fun3_fun1_mat1_MK.zeros();
	entry_point1_fun3_fun1_mat2_KK.zeros();

	for (n_index=0; n_index<N_static; n_index++){

	entry_point1_fun3_fun1_mat1_MK=entry_point1_fun3_fun1_mat1_MK+Xtilde_mfn_local.slice(n_index).col(f_index_admissible)*trans(Y_est_KxFxN_local.slice(n_index).col(f_index_admissible));

	entry_point1_fun3_fun1_mat2_KK=entry_point1_fun3_fun1_mat2_KK+Y_est_KxFxN_local.slice(n_index).col(f_index_admissible)*trans(Y_est_KxFxN_local.slice(n_index).col(f_index_admissible));

	}

/*	A_pseudoinv_MKF.slice(f_index)=(1/((double)N_static))*entry_point1_fun3_fun1_mat1_MK*pinv(entry_point1_fun3_fun1_mat2_KK);*/
	A_pseudoinv_MKF.slice(f_index)=(1/((double)N_static))*entry_point1_fun3_fun1_mat1_MK*pinv((1/((double)N_static))*entry_point1_fun3_fun1_mat2_KK);

}

}

static cube::fixed<M_static, K_static, F_static_admissible> rhat_MKF; 
static mat::fixed<M_static, K_static> rhat_MK; 

static void entry_point1_fun3_fun1_compute_r_k_hat(void){

int f_index, k_index, m_index;

double freq_value;

/*populate rhat*/
for (k_index=0; k_index<K_static; k_index++){

	rhat_MK.col(k_index).zeros();

	for (f_index=0; f_index<F_static_admissible; f_index++){

		for (m_index=0; m_index<M_static; m_index++){

			freq_value=(((double)f_index+1)*((double)FS_CONSTANT))/((double)N_STFT_CONSTANT);

			rhat_MKF(m_index, k_index, f_index)=-compute_angle_tdoa_update(A_pseudoinv_MKF(m_index, k_index, f_index)/A_pseudoinv_MKF((int)J_VALUE, k_index, f_index))/(2*(datum::pi)*freq_value);

		}

	rhat_MK.col(k_index)=rhat_MK.col(k_index)+rhat_MKF.slice(f_index).col(k_index);

	}

	rhat_MK.col(k_index)=(1/((double)F_static_admissible))*rhat_MK.col(k_index);

}

}

static void entry_point1_fun3_compute_A_rhat(cx_cube* Xtilde_fnm_p){

entry_point1_fun3_fun0_rotate_X_Y(Xtilde_fnm_p);

entry_point1_fun3_fun1_compute_A_mkf();

/*A_pseudoinv_MKF.slice(0).print("A_pseudoinv_MKF.slice(0) line 137;");*/

entry_point1_fun3_fun1_compute_r_k_hat();

rhat_MK.print("rhat_MK, line 141:");

}

static mat::fixed<M_static, L_static> chat_mq_estimate;
static colvec::fixed<L_static> kth_rand_value_for_qth_starting_cluster_center;

static void entry_point1_fun4_fun0_random_initialization(void){

int q_index;

kth_rand_value_for_qth_starting_cluster_center.randu();
kth_rand_value_for_qth_starting_cluster_center=(((double)K_static)-1)*kth_rand_value_for_qth_starting_cluster_center;
floor(kth_rand_value_for_qth_starting_cluster_center);

for (q_index=0; q_index<L_static; q_index++){

/*kth_rand_value_for_qth_starting_cluster_center(q_index)=randu(1, distr_param(0, K_static-1));*/

chat_mq_estimate.col(q_index)=rhat_MK.col((int)(kth_rand_value_for_qth_starting_cluster_center(q_index)))+0.0001*randn<colvec>(M_static);

}

}

static mat::fixed<L_static, K_static> indicators_P_qk;

static mat::fixed<L_static, K_static> indicators_P_qk_frob_den_mat;

static rowvec::fixed<3> index_display_rowvec;
static rowvec::fixed<3> frob_display_rowvec_kq_kr_frobden;

static void entry_point1_fun4_fun1_compute_indicators(void){

int k_index, q_index, r_index;

double frob_norm_kq;
double frob_norm_kr;

double frob_den; 

/*rhat_MK.print("line 181, rhat_MK");
chat_mq_estimate.print("line 182, chat_mq_estimate");*/

for (q_index=0; q_index<L_static; q_index++){

	index_display_rowvec(0)=q_index;

	/*for (k_index=0; k_index<(L_static*F_static-L_static); k_index++){*/
	for (k_index=0; k_index<(K_static); k_index++){

	index_display_rowvec(1)=k_index;

	/*index_display_rowvec.print("display indices: q_index, k_index");*/

		frob_den=0;

		frob_norm_kq=sqrt(as_scalar(sum(square(abs( rhat_MK.col(k_index)-chat_mq_estimate.col(q_index) )))));
		/*frob_norm_kq=(as_scalar(sum(square(abs( rhat_MK.col(k_index)-chat_mq_estimate.col(q_index) )))));*/

		frob_display_rowvec_kq_kr_frobden(0)=frob_norm_kq;

		/*frob_display_rowvec_kq_kr_frobden.print("frob_display_rowvec_kq_kr_frobden, frob_norm_kq changed:");*/

		for (r_index=0; r_index<L_static; r_index++){

			index_display_rowvec(2)=r_index;

			/*index_display_rowvec.print("display indices: q_index, k_index, r_index");*/

			frob_norm_kr=sqrt(as_scalar(sum(square(abs( rhat_MK.col(k_index)-chat_mq_estimate.col(r_index) )))));
			/*frob_norm_kr=(as_scalar(sum(square(abs( rhat_MK.col(k_index)-chat_mq_estimate.col(r_index) )))));*/

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

static colvec::fixed<M_static> means_computation_dummy_Mx1_2;
static double means_computation_accu_double_dummy_2;

static void entry_point1_fun4_fun2_compute_means(void){

int q_index, k_index;

for (q_index=0; q_index<L_static; q_index++){

	means_computation_dummy_Mx1_2.zeros();
	means_computation_accu_double_dummy_2=0;

	for (k_index=0; k_index<K_static; k_index++){

	means_computation_dummy_Mx1_2=means_computation_dummy_Mx1_2+pow(indicators_P_qk(q_index, k_index), 2)*rhat_MK.col(k_index);

	means_computation_accu_double_dummy_2=means_computation_accu_double_dummy_2+pow(indicators_P_qk(q_index, k_index), 2);

	}

	chat_mq_estimate.col(q_index)=means_computation_dummy_Mx1_2/means_computation_accu_double_dummy_2;

}

}

static void entry_point1_fun4_compute_k_means(void){

entry_point1_fun4_fun0_random_initialization();

chat_mq_estimate.print("chat_mq_estimate, line 232:");

int iter_ctr;





for (iter_ctr=0; iter_ctr<DESIRED_NO_ITERATIONS_TDOAU_M3_KMEANS; iter_ctr++){

entry_point1_fun4_fun1_compute_indicators();

entry_point1_fun4_fun2_compute_means();

}

/*indicators_P_qk.print("indicators_P_qk, line 236");

indicators_P_qk_frob_den_mat.print("indicators_P_qk_frob_den_mat, line 252"); */

chat_mq_estimate.print("chat_mq_estimate, line 242");

}

static cx_cube::fixed<M_static, L_static, F_static_admissible> h_mlf_local;

static void entry_point1_fun5_compute_h_flm(mat* Z_ol_p, cube* W_fom_p, cx_cube* expj_Phi_W_fom_p){

int f_index, l_index, m_index, o_index, f_index_admissible;
cx_double dummy_value_local=0; 

for (f_index=0; f_index<F_static_admissible; f_index++){

	f_index_admissible=f_index+1;

	for (l_index=0; l_index<L_static; l_index++){

		for (m_index=0; m_index<M_static; m_index++){

			dummy_value_local=0;

			for (o_index=0; o_index<O_static; o_index++){

				dummy_value_local=dummy_value_local+((*W_fom_p)(f_index_admissible, o_index, m_index))*((*expj_Phi_W_fom_p)(f_index_admissible, o_index, m_index))*((*Z_ol_p)(o_index, l_index));

			}

			h_mlf_local(m_index, l_index, f_index)=dummy_value_local;

		}

	}

}

}

static cube::fixed<M_static, L_static, F_static_admissible> rhat_MLF; 
static mat::fixed<M_static, L_static> rhat_ML; 									

static void entry_point1_fun6_compute_rhat(void){

int f_index, l_index, m_index;

double freq_value;

/*populate rhat*/
for (l_index=0; l_index<L_static; l_index++){

	rhat_ML.col(l_index).zeros();

	for (f_index=0; f_index<F_static_admissible; f_index++){

		for (m_index=0; m_index<M_static; m_index++){

			freq_value=(((double)f_index+1)*((double)FS_CONSTANT))/((double)N_STFT_CONSTANT);

			rhat_MLF(m_index, l_index, f_index)=-compute_angle_tdoa_update(h_mlf_local(m_index, l_index, f_index)/h_mlf_local((int)J_VALUE, l_index, f_index))/(2*(datum::pi)*freq_value);

		}

	rhat_ML.col(l_index)=rhat_ML.col(l_index)+rhat_MLF.slice(f_index).col(l_index);

	}

	rhat_ML.col(l_index)=(1/((double)F_static_admissible))*rhat_ML.col(l_index);

}

}

/*#define NUM_PERMUTATIONS 6*/

cube::fixed<L_static, L_static, NUM_PERMUTATIONS> permutation_list;
rowvec::fixed<L_static> v1_permutation_basis_vec;
rowvec::fixed<L_static> v2_permutation_basis_vec;
rowvec::fixed<L_static> v3_permutation_basis_vec;

/*This function needs to be called serially at least once before calling any other functions provided by this module*/
void populate_6x1_possible_permutations(void){

int iter_ctr=0;

v1_permutation_basis_vec(0)=1;

v2_permutation_basis_vec(1)=1;

v3_permutation_basis_vec(2)=1;

permutation_list.slice(iter_ctr).row(0)=v1_permutation_basis_vec;
permutation_list.slice(iter_ctr).row(1)=v2_permutation_basis_vec;
permutation_list.slice(iter_ctr).row(2)=v3_permutation_basis_vec;

iter_ctr=1;

permutation_list.slice(iter_ctr).row(0)=v1_permutation_basis_vec;
permutation_list.slice(iter_ctr).row(1)=v3_permutation_basis_vec;
permutation_list.slice(iter_ctr).row(2)=v2_permutation_basis_vec;

iter_ctr=2;

permutation_list.slice(iter_ctr).row(0)=v2_permutation_basis_vec;
permutation_list.slice(iter_ctr).row(1)=v1_permutation_basis_vec;
permutation_list.slice(iter_ctr).row(2)=v3_permutation_basis_vec;

iter_ctr=3;

permutation_list.slice(iter_ctr).row(0)=v2_permutation_basis_vec;
permutation_list.slice(iter_ctr).row(1)=v3_permutation_basis_vec;
permutation_list.slice(iter_ctr).row(2)=v1_permutation_basis_vec;

iter_ctr=4;

permutation_list.slice(iter_ctr).row(0)=v3_permutation_basis_vec;
permutation_list.slice(iter_ctr).row(1)=v1_permutation_basis_vec;
permutation_list.slice(iter_ctr).row(2)=v2_permutation_basis_vec;

iter_ctr=5;

permutation_list.slice(iter_ctr).row(0)=v3_permutation_basis_vec;
permutation_list.slice(iter_ctr).row(1)=v2_permutation_basis_vec;
permutation_list.slice(iter_ctr).row(2)=v1_permutation_basis_vec;

}

rowvec::fixed<L_static> entry_point1_fun7_fun1_dummy_vec1;
colvec::fixed<NUM_PERMUTATIONS> entry_point1_fun7_fun1_dummy_vec2_sum_of_squares_vector;

static mat::fixed<M_static, L_static> chat_mq_estimate_permutation_aligned;
static mat::fixed<L_static, K_static> indicators_P_qk_permutation_aligned;

static void entry_point1_fun7_compute_minimum_sum_squares_permutation(void){

int i_iter;

int l_index, q_index;

uword q_index_uword, min_permutation_output_index_uword;

int min_permutation_output_index;

double i_iterth_permutation_sum_of_squares_dummy;

/*chat_mq_estimate and rhat_ML*/
for (i_iter=0; i_iter<NUM_PERMUTATIONS; i_iter++){

	i_iterth_permutation_sum_of_squares_dummy=0; 

	for (l_index=0; l_index<L_static; l_index++){

	entry_point1_fun7_fun1_dummy_vec1=permutation_list.slice(i_iter).row(l_index);

	/*q_index=(int)entry_point1_fun7_fun1_dummy_vec1.index_max();*/
	entry_point1_fun7_fun1_dummy_vec1.max(q_index_uword);
	q_index=(int)q_index_uword;

	i_iterth_permutation_sum_of_squares_dummy=i_iterth_permutation_sum_of_squares_dummy+as_scalar(sum(square(abs(rhat_ML.col(l_index)-chat_mq_estimate.col(q_index)))));

	}

	entry_point1_fun7_fun1_dummy_vec2_sum_of_squares_vector(i_iter)=i_iterth_permutation_sum_of_squares_dummy;

}

/*min_permutation_output_index=entry_point1_fun7_fun1_dummy_vec2_sum_of_squares_vector.index_min();*/
entry_point1_fun7_fun1_dummy_vec2_sum_of_squares_vector.min(min_permutation_output_index_uword);
min_permutation_output_index=(int)min_permutation_output_index_uword;

chat_mq_estimate_permutation_aligned=chat_mq_estimate*permutation_list.slice(min_permutation_output_index);

indicators_P_qk_permutation_aligned=permutation_list.slice(min_permutation_output_index)*indicators_P_qk;

}

static colvec::fixed<L_static> frob_norm_colvec_norm_Y_comp_Z;

/*Frobenius norm normalize both Y_lk and P_qk */
static void frob_norm_normalize_Y_lk_local(mat* Y_lk_p){

int l_iter;

/*compute each frobenius norm in Y_lk for l=1:L*/

frob_norm_colvec_norm_Y_comp_Z=sqrt(sum((*Y_lk_p)%(*Y_lk_p), 1));

/*Divide out each row by its frob norm:*/
(*Y_lk_p)=(*Y_lk_p)/kron(frob_norm_colvec_norm_Y_comp_Z, ones_row_1xK);

}

static void entry_point1_fun8_fun1_update_Y_lk(mat* Y_lk_p){

mat* P_qk_p=&indicators_P_qk_permutation_aligned;

/*normalize Y*/
frob_norm_normalize_Y_lk_local(Y_lk_p);

/*normalize P*/
frob_norm_normalize_Y_lk_local(P_qk_p);

/*(*Y_lk_p)=((*Y_lk_p)+(*P_qk_p))/2;*/

(*Y_lk_p)=(0.1*(*Y_lk_p)+0.9*(*P_qk_p));

/*normalize Y again*/
frob_norm_normalize_Y_lk_local(Y_lk_p);

}

static mat::fixed<M_static, L_static> rhat_ML_target; 	

static cx_cube::fixed<M_static, L_static, F_static_admissible> bhat_MLF; 
static cube::fixed<M_static, L_static, F_static_admissible> bhat_MLF_real; 
static cube::fixed<M_static, L_static, F_static_admissible> bhat_MLF_imag; 

static void entry_point1_fun8_fun2_compute_bhat(void){

int f_index, l_index, m_index;

double freq_value;	

for (l_index=0; l_index<L_static; l_index++){

	for (f_index=0; f_index<F_static_admissible; f_index++){

		freq_value=(((double)f_index+1)*((double)FS_CONSTANT))/((double)N_STFT_CONSTANT);

		bhat_MLF_real.slice(f_index).col(l_index)=real(cos(-2*(datum::pi)*freq_value*rhat_ML_target.col(l_index)));

		bhat_MLF_imag.slice(f_index).col(l_index)=imag(sin(-2*(datum::pi)*freq_value*rhat_ML_target.col(l_index)));

	}

}

bhat_MLF.set_real(bhat_MLF_real);
bhat_MLF.set_imag(bhat_MLF_imag);

}

static mat::fixed<O_static, L_static> partial_d_wrt_Zol_1;
static mat::fixed<O_static, L_static> partial_d_wrt_Zol_2;
static cx_colvec::fixed<F_static_admissible> entry_point1_fun5_colvec1_Fx1;
static cx_colvec::fixed<F_static_admissible> entry_point1_fun5_colvec2_Fx1;
static cx_colvec::fixed<F_static_admissible> entry_point1_fun5_colvec3_Fx1;

static void entry_point1_fun8_fun3_update_Z_ol(mat* Z_ol_p, cube* W_fom_p, cx_cube* expj_Phi_W_fom_p){

int o_index, l_index, f_index, f_index_admissible;

int m_index=1; 

/*double gamma_local=50;*/

double gamma_local=10;

double epsilon_value; 

for (o_index=0; o_index<O_static; o_index++){

	for (l_index=0; l_index<L_static; l_index++){

		for (f_index=0; f_index<F_static_admissible; f_index++){

			f_index_admissible=f_index+1; 

			entry_point1_fun5_colvec1_Fx1(f_index)=conj(bhat_MLF(m_index, l_index, f_index));

			entry_point1_fun5_colvec3_Fx1(f_index)=h_mlf_local((int)J_VALUE, l_index, f_index)*((*W_fom_p)(f_index_admissible, o_index, m_index))*((*expj_Phi_W_fom_p)(f_index_admissible, o_index, m_index))-h_mlf_local(m_index, l_index, f_index)*((*W_fom_p)(f_index_admissible, o_index, (int)J_VALUE))*((*expj_Phi_W_fom_p)(f_index_admissible, o_index, (int)J_VALUE));

			entry_point1_fun5_colvec3_Fx1(f_index)=entry_point1_fun5_colvec3_Fx1(f_index)/(pow(abs(h_mlf_local(m_index, l_index, f_index)), 2));

			entry_point1_fun5_colvec1_Fx1(f_index)=entry_point1_fun5_colvec1_Fx1(f_index)*entry_point1_fun5_colvec3_Fx1(f_index);

			entry_point1_fun5_colvec2_Fx1(f_index)=h_mlf_local(m_index, l_index, f_index)/h_mlf_local((int)J_VALUE, l_index, f_index);

			entry_point1_fun5_colvec2_Fx1(f_index)=entry_point1_fun5_colvec2_Fx1(f_index)/abs(entry_point1_fun5_colvec2_Fx1(f_index));

			entry_point1_fun5_colvec2_Fx1(f_index)=entry_point1_fun5_colvec2_Fx1(f_index)*entry_point1_fun5_colvec3_Fx1(f_index);

		}

		partial_d_wrt_Zol_1(o_index, l_index)=-2*as_scalar(sum(real(entry_point1_fun5_colvec1_Fx1)));

		partial_d_wrt_Zol_2(o_index, l_index)=as_scalar(sum(real(entry_point1_fun5_colvec2_Fx1)));

		epsilon_value=(1/gamma_local)*(*Z_ol_p)(o_index, l_index)/partial_d_wrt_Zol_2(o_index, l_index);

		(*Z_ol_p)(o_index, l_index)=(*Z_ol_p)(o_index, l_index)-epsilon_value*(partial_d_wrt_Zol_1(o_index, l_index)+partial_d_wrt_Zol_2(o_index, l_index));

	}

}

}

static void entry_point1_fun8_compute_targets_and_midpoints(mat* Y_lk_p, mat* Z_ol_p, cube* W_fom_p, cx_cube* expj_Phi_W_fom_p){

/*rhat_ML_target=(rhat_ML+chat_mq_estimate_permutation_aligned)/2;*/
rhat_ML_target=(0.1*rhat_ML+0.9*chat_mq_estimate_permutation_aligned);	

rhat_ML_target.print("rhat_ML_target, line 592:");

entry_point1_fun8_fun1_update_Y_lk(Y_lk_p);

/*(*Y_lk_p).print("Y_lk: line 536");*/

entry_point1_fun8_fun2_compute_bhat();

entry_point1_fun8_fun3_update_Z_ol(Z_ol_p, W_fom_p, expj_Phi_W_fom_p);

/*(*Z_ol_p).print("Y_lk: line 542");*/

}

void TDOA_update_module3_entry_point1(cx_cube* Xtilde_fnm_p, mat* Z_ol_p, mat* Y_lk_p, mat* T_fk_p, mat* V_nk_p, cube* W_fom_p, cx_cube* expj_Phi_W_fom_p, cx_cube* expj_Phi_S_nkf_p){

entry_point1_fun1_rotate_Phi_S(expj_Phi_S_nkf_p);

entry_point1_fun2_compute_Y_est(T_fk_p, V_nk_p);

entry_point1_fun3_compute_A_rhat(Xtilde_fnm_p);

entry_point1_fun4_compute_k_means();	/*q information should be available by this point*/

entry_point1_fun5_compute_h_flm(Z_ol_p, W_fom_p, expj_Phi_W_fom_p);

/*Compute rhat_MLF and rhat_ML from h_flm instead of from A_pseudoinv_MLF */
/*Might have to rotate h_flm first. */

/*need to resolve permutations between q and l from this point on. What exactly does this mean?*/

/*Recall that our target quantity is rhat_ML
Need to find a mapping (permutation) from q=1:L to l=1:L that minimizes a sum of squares criterion...
between chat_mq_estimate and rhat_ML.
*/

/*This permutation should not only tell us how to re-compute rhat_ML and h_flm, 
but also should tell us how to use indicators_P_qk to update Y_lk

For re-computing rhat_ML and h_flm, use the permutation obtained from chat_mq_estimate, to select the q mean that best matches a certain l mean, then take the midpoint between the two means (of the TDOA vectors),
use the new midpoints indexed by l to generate a new set of b values. and call Z_ol update given the new b values. This should effectively and properly update the L TDOA vectors for all the classes rhat_ML; in the correct direction. 

ex: may simply need to do a rearrangment of the rows of P_qk first. Could then take a midpoint between Y_lk and P_qk (q rearranged) to update Y_lk
 */

/*Compute rhat_ML*/
entry_point1_fun6_compute_rhat();

/*Compute the 6x1 sequence of permutations of 3x3 permutation matrices. Don't need to use Heap's algorithm since 3!=6, however, could use it for 4!=24 or greater*/
/*Function to populate it should be called serially, not here, but at least once before calling this function frequently */

/*compute the permutation that minimizes the sum of squares*/
entry_point1_fun7_compute_minimum_sum_squares_permutation();

/*Compute midpoints according to permutation aligned things*/
entry_point1_fun8_compute_targets_and_midpoints(Y_lk_p, Z_ol_p, W_fom_p, expj_Phi_W_fom_p);

/*Given the midpoints/targets, update Z_ol and Y_lk */

}

