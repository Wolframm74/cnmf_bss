#include "local_inc.hpp"

/*STFT CONSTANTS*/

#define DESIRED_NO_ITERATIONS 50

static cx_mat::fixed<F_static*N_static, M_static> rotation_dummymat_FNxM;
static cx_mat::fixed<M_static, F_static*N_static> rotation_dummymat_MxFN;

static cx_cube::fixed<M_static, F_static, N_static> Xtilde_mfn_rotated;

/*memcpy() example usage:
arrayops::copy(dest, source, numel);
*/

static void Xtilde_fnm_rotation_fun1(cx_cube* Xtilde_fnm_p){

/*copy (probably memcpy) from FNM into output FNxM*/
arrayops::copy(rotation_dummymat_FNxM.memptr(), (*Xtilde_fnm_p).memptr(), F_static*N_static*M_static);

/*transpose the result and place the MxFN output somewhere*/	
rotation_dummymat_MxFN=strans(rotation_dummymat_FNxM);

/*memcpy the transposed result into an MFN output*/	
arrayops::copy(Xtilde_mfn_rotated.memptr(), rotation_dummymat_MxFN.memptr(), F_static*N_static*M_static);

}

static cx_cube::fixed<M_static, F_static, N_static> Xtilde_mfn_normalized;
/*static cx_cube::fixed<F_static, N_static, M_static> Xtilde_fnm_normalized;*/

/*static void Xtilde_fnm_normalization_fun1(cx_cube* Xtilde_fnm_p){*/
static void Xtilde_fnm_normalization_fun1(void){

int f_index, n_index, m_index;

cx_double phase_X_m1_fn;
double frob_norm_X_fn;

/*Xtilde_fnm_normalized(f_index, n_index, m_index)=;*/

for (f_index=0; f_index<F_static; f_index++){

	for (n_index=0; n_index<N_static; n_index++){

		/*phase_X_m1_fn=conj(((*Xtilde_fnm_p)(f_index, n_index, 1))/abs((*Xtilde_fnm_p)(f_index, n_index, 1)))*/

		/*to get just the phase of m=1th element divide out its absolute value. */
		phase_X_m1_fn=conj((Xtilde_mfn_rotated(1, f_index, n_index))/abs(Xtilde_mfn_rotated(1, f_index, n_index)));

		/*frob_norm_X_fn=0;*/

		/*Xtilde_fnm_normalized(f_index, n_index, m_index)=((*Xtilde_fnm_p)(f_index, n_index, m_index))*phase_X_m1_fn;*/

		frob_norm_X_fn=sqrt(as_scalar(sum(abs(Xtilde_mfn_rotated.slice(n_index).col(f_index))%abs(Xtilde_mfn_rotated.slice(n_index).col(f_index)))));

		Xtilde_mfn_normalized.slice(n_index).col(f_index)=(phase_X_m1_fn*Xtilde_mfn_rotated.slice(n_index).col(f_index))/frob_norm_X_fn;

	}

}

}

colvec::fixed<L_static> nth_rand_value_for_f_lth_starting_cluster_center;

cx_cube::fixed<M_static, L_static, F_static> h_mlf_kmeans_estimate_output;

static void step1_fun_random_initialization(int f_index){

int l_index;

nth_rand_value_for_f_lth_starting_cluster_center.randu();
nth_rand_value_for_f_lth_starting_cluster_center=(((double)N_static)-1)*nth_rand_value_for_f_lth_starting_cluster_center;
floor(nth_rand_value_for_f_lth_starting_cluster_center);

/*call randi() to return random integer values*/

for (l_index=0; l_index<L_static; l_index++){

/*nth_rand_value_for_f_lth_starting_cluster_center(l_index)=randi(1, distr_param(0, N_static-1));*/

h_mlf_kmeans_estimate_output.slice(f_index).col(l_index)=Xtilde_mfn_normalized.slice((int)(nth_rand_value_for_f_lth_starting_cluster_center(l_index))).col(f_index);

}

}

mat::fixed<L_static, N_static> indicators_B_ln; 

static void step2_fun_indicators_computation(int f_index){

int l_index, n_index, r_index;

double frob_norm_nl;
double frob_norm_nr;


double frob_den; 

for (l_index=0; l_index<L_static; l_index++){

	for (n_index=0; n_index<N_static; n_index++){

		frob_den=0;

		frob_norm_nl=sqrt(as_scalar(sum(square(abs(Xtilde_mfn_normalized.slice(n_index).col(f_index)-h_mlf_kmeans_estimate_output.slice(f_index).col(l_index))))));

		for (r_index=0; r_index<L_static; r_index++){

		frob_norm_nr=sqrt(as_scalar(sum(square(abs(Xtilde_mfn_normalized.slice(n_index).col(f_index)-h_mlf_kmeans_estimate_output.slice(f_index).col(r_index))))));

		frob_den=frob_den+pow((frob_norm_nl/frob_norm_nr), 2);

		}

	indicators_B_ln(l_index, n_index)=1/frob_den;

	}

}

}

cx_colvec::fixed<M_static> means_computation_dummy_Mx1;
static double means_computation_accu_double_dummy;

static void step3_fun_means_computation(int f_index){

int l_index, n_index;

for (l_index=0; l_index<L_static; l_index++){

	means_computation_dummy_Mx1.zeros();
	means_computation_accu_double_dummy=0;

	for (n_index=0; n_index<N_static; n_index++){

	means_computation_dummy_Mx1=means_computation_dummy_Mx1+pow(indicators_B_ln(l_index, n_index), 2)*Xtilde_mfn_normalized.slice(n_index).col(f_index);

	means_computation_accu_double_dummy=means_computation_accu_double_dummy+pow(indicators_B_ln(l_index, n_index), 2);

	}

h_mlf_kmeans_estimate_output.slice(f_index).col(l_index)=means_computation_dummy_Mx1/means_computation_accu_double_dummy;

}

}



void TDOA_update_module2_entry_point1(mxArray *plhs[], const mxArray *prhs[], arg_struct_t* argStruct_p){

int f_index=0; 

int iter_ctr;

/*step 0: rotation so that m fits into columns */
Xtilde_fnm_rotation_fun1(argStruct_p->Xtilde_fnm_p);

/*step 0.5: normalization*/
Xtilde_fnm_normalization_fun1();

for (f_index=0; f_index<F_static; f_index++){

/*step 1: randomly initialize h_l=1(f) ... h_l=L(f)*/
step1_fun_random_initialization(f_index);

	/*loop over 2 and 3 for a fixed number of satisfactory iterations*/

	for (iter_ctr=0; iter_ctr<DESIRED_NO_ITERATIONS; iter_ctr++){

	/*step 2: compute assignment of fuzzy indicators b_ln */
	step2_fun_indicators_computation(f_index);

	/*step 3: compute assignment of means h_l(f) for l=1:L_static */
	step3_fun_means_computation(f_index);

	}

}

}


double compute_angle_tdoa_update(cx_double input){

double input_real;
double input_im;

double angle;

int case_int; 

input_real=real(input);
input_im=imag(input);

/*Case 1: Check if 'input' is >= 0. return angle=0;*/

if (input_im==0){

if (input_real>=0){

angle=0;

}

/*Case 3: Check if 'input' is exactly real and negative. return angle=pi*/	
else if (input_real<0){

angle=datum::pi; 

}

} else {

/*Otherwise simply call atan2 and return the 'angle' of 'input' */

angle=2*atan(input_im/(input_real+sqrt(pow(input_real, 2)+pow(input_im, 2))));

}

return angle; 

}

cx_cube::fixed<M_static, L_static, F_static_admissible> rhat_mlf_kmeans_estimate_output;	/*subtract off 1 from F_SPATIAL_ALIASING_MAX_INDEX in order to disregard the frequency 0 set of values. */
cx_mat::fixed<M_static, E_static> rhat_me_kmeans_estimate_output;

/*This function whose name lacks description does several things the most important of which is to correct compute rhat_mlf_kmeans_estimate_output from the complex domain valued h_mlf_kmeans_estimate_output.

Note that we skip over accessing the frequency 0th frequency of h_mlf_kmeans_estimate_output (ie: this is the case wherever you see (f_index+1) instead of just f_index, 
because this frequency is inadmissible, along with the frequencies above the spatial aliasing frequency, which we leave as part of the frequency set for now. 

*/
static void entry_point2_fun1(void){

int f_index, l_index, m_index;

double freq_value; 

int J_value=0; 

/*Could later modify this code to include only nonambiguous frequencies below the spatial aliasing frequency. */
/*Would need to re-define E, thus consequently resizing the container rhat_me_kmeans_estimate_output, for example */
for (f_index=0; f_index<(F_static_admissible); f_index++){

	for (l_index=0; l_index<L_static; l_index++){

		for (m_index; m_index<M_static; m_index++){

			freq_value=(((double)f_index+1)*((double)FS_CONSTANT))/((double)N_STFT_CONSTANT);

			rhat_mlf_kmeans_estimate_output(m_index, l_index, f_index)=-compute_angle_tdoa_update(h_mlf_kmeans_estimate_output(m_index, l_index, f_index+1)/h_mlf_kmeans_estimate_output(J_value, l_index, f_index+1))/(2*(datum::pi)*freq_value);

		}	

	}

}

arrayops::copy(rhat_me_kmeans_estimate_output.memptr(), rhat_mlf_kmeans_estimate_output.memptr(), M_static*E_static );

}

cx_mat::fixed<M_static, L_static> chat_mq_estimate;
colvec::fixed<L_static> eth_rand_value_for_qth_starting_cluster_center;

static void entry_point2_fun2_random_initialization(void){

int q_index;

eth_rand_value_for_qth_starting_cluster_center.randu();
eth_rand_value_for_qth_starting_cluster_center=(((double)E_static)-1)*eth_rand_value_for_qth_starting_cluster_center;
floor(eth_rand_value_for_qth_starting_cluster_center);

for (q_index=0; q_index<L_static; q_index++){

/*eth_rand_value_for_qth_starting_cluster_center(q_index)=randu(1, distr_param(0, E_static-1));*/

chat_mq_estimate.col(q_index)=rhat_me_kmeans_estimate_output.col((int)(eth_rand_value_for_qth_starting_cluster_center(q_index)));

}

}

mat::fixed<L_static, E_static> indicators_G_qe;

static void entry_point2_fun3_indicators_computation(void){

int e_index, q_index, r_index;

double frob_norm_eq;
double frob_norm_er;

double frob_den; 

for (q_index=0; q_index<L_static; q_index++){

	/*for (e_index=0; e_index<(L_static*F_static-L_static); e_index++){*/
	for (e_index=0; e_index<(E_static); e_index++){

		frob_den=0;

		frob_norm_eq=sqrt(as_scalar(sum(square(abs(rhat_me_kmeans_estimate_output.col(e_index)-chat_mq_estimate.col(q_index))))));

		for (r_index=0; r_index<L_static; r_index++){

			frob_norm_er=sqrt(as_scalar(sum(square(abs(rhat_me_kmeans_estimate_output.col(e_index)-chat_mq_estimate.col(r_index))))));

			frob_den=frob_den+pow((frob_norm_eq/frob_norm_er), 2);

		}

		indicators_G_qe(q_index, e_index)=1/frob_den;

	}

}

}

cx_colvec::fixed<M_static> means_computation_dummy_Mx1_2;
static double means_computation_accu_double_dummy_2;

static void entry_point2_fun4_mean_computation(void){

int q_index, e_index;

for (q_index=0; q_index<L_static; q_index++){

	means_computation_dummy_Mx1_2.zeros();
	means_computation_accu_double_dummy_2=0;

	for (e_index=0; e_index<E_static; e_index++){

	means_computation_dummy_Mx1_2=means_computation_dummy_Mx1_2+pow(indicators_G_qe(q_index, e_index), 2)*rhat_me_kmeans_estimate_output.col(e_index);

	means_computation_accu_double_dummy_2=means_computation_accu_double_dummy_2+pow(indicators_G_qe(q_index, e_index), 2);

	}

	chat_mq_estimate.col(q_index)=means_computation_dummy_Mx1_2/means_computation_accu_double_dummy_2;

}

}


void TDOA_update_module2_entry_point2(mxArray *plhs[], const mxArray *prhs[], arg_struct_t* argStruct_p){

/*convert the h_mlf information into r_mlf (ie: complex TDOA info into time domain TDOA information). populate rhat_mlf_kmeans_estimate_output */
/*This function is going to reshape rhat_mlf_kmeans_estimate_output into MxLF (roughly), more precisely reshapes it into MxE, where E=(F-1)*L */
entry_point2_fun1();

entry_point2_fun2_random_initialization();

/*Dont forget to wrap the following two functions in a loop for some desired No. of iterations */

entry_point2_fun3_indicators_computation();

entry_point2_fun4_mean_computation();

}