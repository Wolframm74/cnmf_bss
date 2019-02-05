#include "local_inc.hpp"

/*

Plain and simple task is to find the optimal match per 3=L x O=27 look directions. This should correspond to a min in each of the lth class bins. 

/*
1. Index into expj_Phi_Wfom for the b,ath pair. Generate a column vector in R^Fx1

2. Index into chat_fq; for q=1:L

3. populate an LxO matrix

4. populate a 2xL matrix, in the left hand column is l=1:L, in the right hand column is the corresponding minimum value best matching look direction index. 

5. populate Z_ol accordingly. 

6. 

*/

static cx_colvec::fixed<F_static> entry_point1_fun1_lth_dummy_colvec;

static cx_colvec::fixed<F_static> entry_point1_fun1_oth_dummy_colvec;

static rowvec::fixed<O_static> entry_point1_fun1_lth_dummy_rowvec_real;

static int minimum_indices_colvec[L_static];

static mat::fixed<L_static, O_static> entry_point1_fun1_optimal_match_LxO_mat;

static void entry_point1_fun1_lookdirs_search(cx_cube* expj_Phi_W_fom_p, mat* Z_ol_p, int m_index_b, int m_index_a){

int l_index, o_index;

uword o_index_uword;

(*Z_ol_p).zeros();

for (l_index=0; l_index<L_static; l_index++){

	entry_point1_fun1_lth_dummy_colvec.set_real(cos(chat_fq_estimate.col(l_index)));
	entry_point1_fun1_lth_dummy_colvec.set_imag(sin(chat_fq_estimate.col(l_index)));

	for (o_index=0; o_index<O_static; o_index++){

/*	entry_point1_fun1_oth_dummy_colvec_channel_b.set_real();
	entry_point1_fun1_oth_dummy_colvec_channel_b.set_real();

	entry_point1_fun1_oth_dummy_colvec_channel_a.set_real();
	entry_point1_fun1_oth_dummy_colvec_channel_a.set_real();	*/

	entry_point1_fun1_optimal_match_LxO_mat(l_index, o_index)=accu(square(abs(entry_point1_fun1_lth_dummy_colvec-((*expj_Phi_W_fom_p).slice(m_index_b).col(o_index)%conj((*expj_Phi_W_fom_p).slice(m_index_a).col(o_index))))));

	}

	entry_point1_fun1_lth_dummy_rowvec_real=entry_point1_fun1_optimal_match_LxO_mat.row(l_index);

	entry_point1_fun1_lth_dummy_rowvec_real.min(o_index_uword);

	minimum_indices_colvec[l_index]=(int)o_index_uword;

	(*Z_ol_p)( minimum_indices_colvec[l_index] , l_index )=1; 

}

}	

static void entry_point1_fun2_populate_expj_Phi_U_flnm(cx_cube* expj_Phi_W_fom_p, mat* Z_ol_p){

int m_index;
int n_index;
int l_index;

for (m_index=0; m_index<M_static; m_index++){

	for (l_index=0; l_index<L_static; l_index++){

		for (n_index=0; n_index<N_static; n_index++){

			if (indicators_P_qn(l_index, n_index)<((double)0.15)){

			expj_Phi_U_flnm[m_index].slice(n_index).col(l_index)=conj((*expj_Phi_W_fom_p).slice(m_index).col(minimum_indices_colvec[l_index]));

			}

		}

	}

}

}

static mat::fixed<K_static, N_static> indicators_P_kn;	

static cx_mat::fixed<N_static, K_static*F_static> mat_NxKF;
static cx_mat::fixed<K_static, F_static> dummy_cx_mat_KF;
static cx_rowvec::fixed<K_static*F_static> dummy_cx_rowvec_KF;

void complex_argument_interchannel_clustering_m1_lookdirs_search_populate_ones_expj_Phi_S_fkn(cx_cube* expj_Phi_S_fkn_p, cx_cube* expj_Phi_S_nkf_p){

int n_index, k_index; 

/*indicators_P_kn=trans(*Y_lk_p)*indicators_P_qn;*/

for (n_index=0; n_index<N_static; n_index++){

	for (k_index=0; k_index<K_static; k_index++){

	if (indicators_P_kn(k_index, n_index)<((double)0.15)){

		(*expj_Phi_S_fkn_p).slice(n_index).col(k_index).ones();

	}

	}

	dummy_cx_mat_KF=strans((*expj_Phi_S_fkn_p).slice(n_index));

	/*arrayops::copy(entry_point1_fun1_mat1_NKxF.memptr(), (*expj_Phi_S_nkf_p).memptr(), F_static*N_static*K_static);*/
	arrayops::copy(dummy_cx_rowvec_KF.memptr(), dummy_cx_mat_KF.memptr(), K_static*F_static);

	mat_NxKF.row(n_index)=dummy_cx_rowvec_KF;

}

arrayops::copy((*expj_Phi_S_nkf_p).memptr(), mat_NxKF.memptr(), N_static*K_static*F_static);

}

static mat::fixed<N_static, K_static> V_nk_mask;

static void complex_argument_interchannel_clustering_m1_lookdirs_search_boost_and_suppress_V_kn_compute_mask(void){

int n_index, k_index; 

V_nk_mask.ones();

/*V_nk_mask.elem(find( trans(indicators_P_kn) > 0.55 )).fill(1.5);*/

/*V_nk_mask.elem(find( trans(indicators_P_kn) > 0.2 )).fill(1.01);*/

V_nk_mask.elem(find( trans(indicators_P_kn) < 0.05 )).fill(0.8);

V_nk_mask.elem(find( trans(indicators_P_kn) < 0.1 )).fill(0.9);

V_nk_mask.elem(find( trans(indicators_P_kn) > 0.2 )).fill(1.2);

V_nk_mask.elem(find( trans(indicators_P_kn) > 0.4 )).fill(1.5);


}

void complex_argument_interchannel_clustering_m1_lookdirs_search_boost_and_suppress_V_kn(mat* V_nk_p){

(*V_nk_p)=(*V_nk_p)%V_nk_mask;

}

/*mat::fixed<L_static, N_static> indicators_P_qn_balanced_outmat;*/
rowvec::fixed<N_static> row_energy_accumed_ep1_fun3;

static void entry_point1_fun3_balance_rows_indicator_matrix(cx_cube* Xtilde_fnm_p, int m_index_b, int m_index_a){

int n_index;

int l_index;

rowvec::fixed<N_static> frob_rowvec;
colvec::fixed<L_static> frob_colvec;
colvec::fixed<L_static> frob_colvec_2;
double frob_avg;

frob_colvec=sum(indicators_P_qn%indicators_P_qn , 1);

frob_colvec=sqrt(frob_colvec);

frob_avg=accu(frob_colvec)/((double)L_static);

for (l_index=0; l_index<L_static; l_index++){

/*assign all rows to have the target average frobenius norm*/

/*Divide out the old frobenius norm by having it on the den, assign the new one by having it on the num*/
indicators_P_qn.row(l_index)=indicators_P_qn.row(l_index)*(frob_avg/frob_colvec(l_index));

}

/*indicators_P_qn_balanced_outmat can be used as a debug variable. If you comment out the function it is useless
and the algorithm depends only on indicators_P_qn as it would typically, otherwise.
*/

/*indicators_P_qn=indicators_P_qn;*/

/*first calculate row_energy_accumed_ep1_fun3*/
row_energy_accumed_ep1_fun3=sum( abs((*Xtilde_fnm_p).slice(m_index_b))+abs((*Xtilde_fnm_p).slice(m_index_a))  ,0);

frob_rowvec=sum(indicators_P_qn%indicators_P_qn, 0);

frob_rowvec=sqrt(frob_rowvec);

/*now to across the rows for n=1:N*/
/*for each n_index index the row vector and treat that element as the target frobenius norm*/

for (n_index=0; n_index<N_static; n_index++){

indicators_P_qn.col(n_index)=indicators_P_qn.col(n_index)*( row_energy_accumed_ep1_fun3(n_index) / frob_rowvec(n_index));

}

frob_colvec=sum(indicators_P_qn%indicators_P_qn , 1);

frob_colvec=sqrt(frob_colvec);

/*reset the indicator matrix to have each row have frob norm frob_avg*/
for (l_index=0; l_index<L_static; l_index++){

/*assign all rows to have the target average frobenius norm*/

/*Divide out the old frobenius norm by having it on the den, assign the new one by having it on the num*/
indicators_P_qn.row(l_index)=indicators_P_qn.row(l_index)*(frob_avg/frob_colvec(l_index));

}

}

void complex_argument_interchannel_clustering_m1_lookdirs_search_entry(mxArray *plhs[], arg_struct_t* argStruct_p, int m_index_b, int m_index_a){

entry_point1_fun1_lookdirs_search(argStruct_p->expj_Phi_W_fom_p, argStruct_p->Z_ol_p, m_index_b, m_index_a);

entry_point1_fun1_optimal_match_LxO_mat.print("entry_point1_fun1_optimal_match_LxO_mat:");

expj_Phi_U_flnm[0].ones();
expj_Phi_U_flnm[1].ones();

/*Beware!! this function modifies indicators_P_qn. Turn it off if you would like to have just the untouched output of k-means.
In which a particular class may be dominant; have stronger indicators than two or more other classes. 
*/
entry_point1_fun3_balance_rows_indicator_matrix(argStruct_p->Xtilde_fnm_p, m_index_b, m_index_a);

/*Do another plot check*/
send_data_to_Matlab_Eng_and_plot(plhs);

/*populate expj_Phi_U_flnm[0] with the conjugate values in a class-dependent way. */
entry_point1_fun2_populate_expj_Phi_U_flnm(argStruct_p->expj_Phi_W_fom_p, argStruct_p->Z_ol_p);

/*Populate this once to spread the results for l=1:L into k=1:K*/
indicators_P_kn=trans(*(argStruct_p->Y_lk_p))*indicators_P_qn;

/*Call it for the first time here, and once every iteration after updating Phi_S to ensure that ones are populated in the correct locations, even after updating Phi_S. Thereby neglecting/overwriting anything in bins where ones should be. */
complex_argument_interchannel_clustering_m1_lookdirs_search_populate_ones_expj_Phi_S_fkn(argStruct_p->expj_Phi_S_fkn_p, argStruct_p->expj_Phi_S_nkf_p);

V_nk_mask.ones();

complex_argument_interchannel_clustering_m1_lookdirs_search_boost_and_suppress_V_kn_compute_mask();

}