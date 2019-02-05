#include "local_inc.hpp"

int iter_ctr;

/*#define Q_static 6*/
/*#define Q_static L_static*/

#define DESIRED_NO_ITERATIONS_COMPLEX_ARGUMENT_INTERCH_CLUSTER_M1 200

/*Premise of this is under the guess that the estimated phase shift values occuring at k bins of disinterest to a class, 
and at time bins n=1:N (still at a particular k bin) of disinterest.. yes.. the phase shift values occuring at these locations 
are most likely WRONG/IRRELEVANT and introduce error into the convergence of the data model towards the correct target/solution. 
And therefore we take the conjugate of (H_fl, l=whatever the class of interest is) in order to yes, make the output phase at 
these time frequency (FxN) bins have output phase zero. exp(j*smth)=1 or arg(exp(j(smth))=0
With the help of again Y_lk, the membership indicator pre defined partition*/

/*
1. Cluster over the argument of phase difference between channels into l=1:L=3 class bins. 

2. This should reveal 3 significant look directions. Need to figure out and narrow down from 27 possible 
look directions what those look directions are.

3. Start with 24 k bin partition into 3 possible partitions with 8 k bins per partition like you 
BEEN doin as similar as to Fevotte 2009. Y_lk is the partition matrix.

4. Scanning across time bins, when you compare the arg(X2*conj(X1)), where we can say this is
 the b, a-th phase difference FxN matrix where b=2, a=1 (ones based indexing). 
 You need to supress the time activation columns of regions not in your class. 

5. Do k means clustering with vectors in R^Fx1, over data index n=1:N,
to reduce the true rank of the data space to target size L=3. Obtain the L=3 phase difference mixing vectors
 which reveal the interchannel phase differences each of which is in R^Fx1, which reveal the time difference
  between the b, a th channels as described in step 4. 

6. Populate expj_Phi_S(f,n,k) in such a way - using the membership matrix as to ensure that time activation
 bins for time regions of "disinterest" to a particular class have zero (green- in matlab colomap domain,
  for arg(X2*conj(X1)) ) output phase!!!!! Regions of disinterest should be entirely determined defined,
 I think by the Y_lk, membership matrix. 

Ex: Say partition for l=1 includes k bins {1, 2, .., 8}

Then k bins {9, 10, .... , 24} should be modified in such a way... by the conj(H_fl, l=1}. And this can and 
must be done within expj_Phi_S(f,n,k) at the correct n bins!! These must be probably saved from the k means
 clustering step, I'm thinking.
Laymens terms: 

"Turn the clocks to make output phase of time activation bins not in your class have output complex argument zero. 
Aka phase = exp(j*argument) = 1"
Regions of disinterest for a particular class are time regions for that class such that
 the b-a'th interchannel difference pattern discovered by k-means for that class are INACTIVE (GREEN).

*/

static mat::fixed<F_static, N_static> complex_interchannel_arg_G_ba_quantity;  
mat::fixed<F_static, Q_static> chat_fq_estimate;
static colvec::fixed<Q_static> nth_rand_value_for_qth_starting_cluster_center;
static colvec::fixed<F_static> dummy_colvec_1;

static void complex_argument_interchannel_clustering_m1_entry_fun0_initialize_means_svd(void){

int q_index;

mat U_agba;
vec s_agba;
mat V_agba;

double neg_min;
double pos_max;

/*call svd*/
svd(U_agba, s_agba, V_agba, complex_interchannel_arg_G_ba_quantity );

for (q_index=0; q_index<Q_static; q_index++){

/*index the qth svd thing. save the result into column vector. Don't forget to switch the sign*/

chat_fq_estimate.col(q_index)=-s_agba(q_index)*U_agba.col(q_index);

dummy_colvec_1=chat_fq_estimate.col(q_index);

neg_min=dummy_colvec_1.min();

pos_max=dummy_colvec_1.max();

if (abs(neg_min)>abs(pos_max)){

/*find and clip any negative values greater than it*/
(dummy_colvec_1).elem(find(dummy_colvec_1<=-pos_max)).fill(-pos_max);

/*The vector should now be balanced in terms of the abs(max)=abs(min) should be equal*/
dummy_colvec_1=((datum::pi)*dummy_colvec_1)/abs(pos_max);

} else if (abs(pos_max)>abs(neg_min)){

/*find and clip any positve values greater than it*/
(dummy_colvec_1).elem(find(dummy_colvec_1>=abs(neg_min))).fill(abs(neg_min));

dummy_colvec_1=((datum::pi)*dummy_colvec_1)/abs(neg_min);

}

chat_fq_estimate.col(q_index)=dummy_colvec_1;

}

}

static void entry_point1_fun1_fun0_random_initialization(void){

int q_index;

nth_rand_value_for_qth_starting_cluster_center.randu();
nth_rand_value_for_qth_starting_cluster_center=(((double)N_static)-1)*nth_rand_value_for_qth_starting_cluster_center;
floor(nth_rand_value_for_qth_starting_cluster_center);

for (q_index=0; q_index<Q_static; q_index++){

/*nth_rand_value_for_qth_starting_cluster_center(q_index)=randu(1, distr_param(0, N_static-1));*/

chat_fq_estimate.col(q_index)=complex_interchannel_arg_G_ba_quantity.col((int)(nth_rand_value_for_qth_starting_cluster_center(q_index)))+0.0001*randn<colvec>(F_static);

}

}

static colvec::fixed<Q_static> entry_point1_fun1_fun1_dummy_vec1;

mat::fixed<Q_static, N_static> indicators_P_qn;
static mat::fixed<N_static, Q_static> indicators_P_nq;

static mat::fixed<Q_static, N_static> indicators_P_qn_frob_den_mat;

static rowvec::fixed<3> index_display_rowvec;
static rowvec::fixed<3> frob_display_rowvec_nq_nr_frobden;

static void entry_point1_fun1_fun1_compute_indicators(void){

int n_index, q_index, r_index;

double frob_norm_nq;
double frob_norm_nr;

double frob_den; 

double outsum;

uword q_index_uword;

/*complex_interchannel_arg_G_ba_quantity.print("line 181, complex_interchannel_arg_G_ba_quantity");
chat_fq_estimate.print("line 182, chat_fq_estimate");*/

for (n_index=0; n_index<(N_static); n_index++){

	outsum=sum(abs(complex_interchannel_arg_G_ba_quantity.col(n_index)));

	if (outsum<(((double)F_static)*1.1)){

	indicators_P_qn.col(n_index).zeros();

	} else {

	for (q_index=0; q_index<Q_static; q_index++){

		frob_den=0;

		/*frob_norm_nq=sqrt(as_scalar(sum(square(abs( complex_interchannel_arg_G_ba_quantity.col(n_index)-chat_fq_estimate.col(q_index) )))));*/
		frob_norm_nq=(as_scalar(sum(square(abs( complex_interchannel_arg_G_ba_quantity.col(n_index)-chat_fq_estimate.col(q_index) )))));

		frob_display_rowvec_nq_nr_frobden(0)=frob_norm_nq;

		/*frob_display_rowvec_nq_nr_frobden.print("frob_display_rowvec_nq_nr_frobden, frob_norm_nq changed:");*/

		for (r_index=0; r_index<Q_static; r_index++){

			index_display_rowvec(2)=r_index;

			/*index_display_rowvec.print("display indices: q_index, n_index, r_index");*/

			/*frob_norm_nr=sqrt(as_scalar(sum(square(abs( complex_interchannel_arg_G_ba_quantity.col(n_index)-chat_fq_estimate.col(r_index) )))));*/
			frob_norm_nr=(as_scalar(sum(square(abs( complex_interchannel_arg_G_ba_quantity.col(n_index)-chat_fq_estimate.col(r_index) )))));

			frob_display_rowvec_nq_nr_frobden(1)=frob_norm_nr;

			frob_den=frob_den+pow((frob_norm_nq/frob_norm_nr), 2);

			frob_display_rowvec_nq_nr_frobden(2)=frob_den;

			/*frob_display_rowvec_nq_nr_frobden.print("frob_display_rowvec_nq_nr_frobden, frob_norm_nr, frobden changed");*/

		}

		/*indicators_P_qn(q_index, n_index)=1/frob_den;*/

		indicators_P_qn(q_index, n_index)=1/frob_den;

		indicators_P_qn_frob_den_mat(q_index, n_index)=frob_den;

	}

	}

entry_point1_fun1_fun1_dummy_vec1=indicators_P_qn.col(n_index);

entry_point1_fun1_fun1_dummy_vec1.max(q_index_uword);

/*indicators_P_qn.col(n_index).zeros();*/

/*if (indicators_P_qn((int)q_index_uword, n_index)>0.2){

indicators_P_qn((int)q_index_uword, n_index)=1;

}*/

}

/*indicators_P_nq=trans(indicators_P_qn);

V_nk_MA_filter_rows(&indicators_P_nq, (int)Q_static);

indicators_P_qn=trans(indicators_P_nq);*/

/*indicators_P_qn.print("indicators_P_qn, line 234:");*/

}

static rowvec::fixed<Q_static> entry_point1_fun1_fun2_dummy_vec1;

static colvec::fixed<F_static> means_computation_dummy_Fx1_2;
static double means_computation_accu_double_dummy_2;

static void entry_point1_fun1_fun2_compute_means(void){

int q_index, n_index;

double current_ratio;
double difference_boost_ratio=0.75; 

uword q_index_strongest_uword;

for (q_index=0; q_index<Q_static; q_index++){

	means_computation_dummy_Fx1_2.zeros();
	means_computation_accu_double_dummy_2=0;

	for (n_index=0; n_index<N_static; n_index++){

	means_computation_dummy_Fx1_2=means_computation_dummy_Fx1_2+pow(indicators_P_qn(q_index, n_index), 2)*complex_interchannel_arg_G_ba_quantity.col(n_index);

	means_computation_accu_double_dummy_2=means_computation_accu_double_dummy_2+pow(indicators_P_qn(q_index, n_index), 2);

	}

	chat_fq_estimate.col(q_index)=means_computation_dummy_Fx1_2/means_computation_accu_double_dummy_2;

}

}



static void complex_argument_interchannel_clustering_m1_entry_fun1_k_means_clustering(void){

/*entry_point1_fun1_fun0_random_initialization();*/



for (iter_ctr=0; iter_ctr<DESIRED_NO_ITERATIONS_COMPLEX_ARGUMENT_INTERCH_CLUSTER_M1; iter_ctr++){

entry_point1_fun1_fun1_compute_indicators();

entry_point1_fun1_fun2_compute_means();

}

}

static cx_mat::fixed<F_static, N_static> dummy_mat_cx_FN_common_local;

static mat::fixed<F_static, N_static> outmat_FN_local;
static mat::fixed<F_static, N_static> nummat_FN_local;
static mat::fixed<F_static, N_static> denmat_FN_local;

/*Function local data*/
static mat& compute_arg_X_FN_local(void){

int f_index, n_index; 

mat& outmat_FN_ref=outmat_FN_local;

/*dummy_mat_cx_FN_common_local=(*X_fn_p);*/

dummy_mat_cx_FN_common_local.elem(find(abs(dummy_mat_cx_FN_common_local)==0)).fill(0.00000001);

nummat_FN_local=imag(dummy_mat_cx_FN_common_local);

denmat_FN_local=sqrt(square(real(dummy_mat_cx_FN_common_local))+square(imag(dummy_mat_cx_FN_common_local)))+real(dummy_mat_cx_FN_common_local);

outmat_FN_local=nummat_FN_local/denmat_FN_local;

outmat_FN_local=2*atan(outmat_FN_local);

/*find nan and set to zero. Everything else should be non-nan */
outmat_FN_local.elem( find_nonfinite(outmat_FN_local )).zeros();

return outmat_FN_ref;

}

/*static mat::fixed<Q_static, N_static> indicators_P_qn;
static mat::fixed<F_static, Q_static> chat_fq_estimate;*/

static bool send_data_to_Matlab_Eng_and_plot_init_flag=false;

static void send_data_to_Matlab_Eng_and_plot_init(mxArray *plhs[]){

plhs[9]=armaCreateMxMatrix(indicators_P_qn.n_rows, indicators_P_qn.n_cols, mxDOUBLE_CLASS, mxREAL);

plhs[10]=armaCreateMxMatrix(chat_fq_estimate.n_rows, chat_fq_estimate.n_cols, mxDOUBLE_CLASS, mxREAL);

send_data_to_Matlab_Eng_and_plot_init_flag=true;

}

void send_data_to_Matlab_Eng_and_plot(mxArray *plhs[]){

if (!send_data_to_Matlab_Eng_and_plot_init_flag){

send_data_to_Matlab_Eng_and_plot_init(plhs);

}

/*plhs[9]=armaCreateMxMatrix(indicators_P_qn.n_rows, indicators_P_qn.n_cols, mxDOUBLE_CLASS, mxREAL);

plhs[10]=armaCreateMxMatrix(chat_fq_estimate.n_rows, chat_fq_estimate.n_cols, mxDOUBLE_CLASS, mxREAL);*/

armaSetPr(plhs[9], indicators_P_qn);

armaSetPr(plhs[10], chat_fq_estimate);

engPutVariable(mlEngine_p, "indicators_P_qn", plhs[9]);

engPutVariable(mlEngine_p, "chat_fq_estimate", plhs[10]);

engEvalString(mlEngine_p, "plot_P_ln_and_C_l(indicators_P_qn, chat_fq_estimate);");

}

void complex_argument_interchannel_clustering_m1_entry(mxArray *plhs[], arg_struct_t* argStruct_p, int m_index_b, int m_index_a){

int q_index;

double current_ratio;
double difference_boost_ratio=0.9; 

uword q_index_strongest_uword;

/*Set the dummy mat*/
dummy_mat_cx_FN_common_local=(*(argStruct_p->Xtilde_fnm_p)).slice(m_index_b)%conj((*(argStruct_p->Xtilde_fnm_p)).slice(m_index_a));

/*arg() function provided as part of the complex_argument_costfun_module1 module */
complex_interchannel_arg_G_ba_quantity=compute_arg_X_FN_local();

/*compute svd as a better initialization of the mean vectors*/
complex_argument_interchannel_clustering_m1_entry_fun0_initialize_means_svd();

send_data_to_Matlab_Eng_and_plot(plhs);

complex_argument_interchannel_clustering_m1_entry_fun1_k_means_clustering();

/*Filter indicators and boost means a bit on the way out*/

/*indicators_P_nq=trans(indicators_P_qn);

V_nk_MA_filter_rows(&indicators_P_nq, (int)Q_static);

indicators_P_qn=trans(indicators_P_nq);*/

/*indicators_P_qn.print("indicators_P_qn, line 234:");*/

/*Boost other means*/

/*entry_point1_fun1_fun2_dummy_vec1=sum(abs(chat_fq_estimate), 0);*/
entry_point1_fun1_fun2_dummy_vec1=sqrt(sum((chat_fq_estimate)%(chat_fq_estimate), 0));

entry_point1_fun1_fun2_dummy_vec1.max(q_index_strongest_uword);

for (q_index=0; q_index<Q_static; q_index++){

	if (q_index!=((int)q_index_strongest_uword)){

	current_ratio=(entry_point1_fun1_fun2_dummy_vec1(q_index))/(entry_point1_fun1_fun2_dummy_vec1.max(q_index_strongest_uword));

	// boost the mean
	chat_fq_estimate.col(q_index)=((current_ratio+(1-current_ratio)*difference_boost_ratio)/(current_ratio))*chat_fq_estimate.col(q_index);

	}

}

/*(*(argStruct_p->V_nk_p)).elem(find((*(argStruct_p->V_nk_p))>1)).fill(1);*/
chat_fq_estimate.elem(find(chat_fq_estimate>(datum::pi))).fill(datum::pi);
chat_fq_estimate.elem(find(chat_fq_estimate<(-datum::pi))).fill(-datum::pi);

/*	chat_fq_estimate.col(1)=1.2*chat_fq_estimate.col(1);
	chat_fq_estimate.col(2)=1.2*chat_fq_estimate.col(2);*/

send_data_to_Matlab_Eng_and_plot(plhs);

}