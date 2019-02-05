#include "local_inc.hpp"

static colvec::fixed<L_static> Z_colvec_before_01;

void Z_pull_vector_elements_towards_mean_wrapper(mat* Z_ol_p){

colvec_dynamic_01=trans(sum(*Z_ol_p, 0));

/*Save*/
Z_colvec_before_01=colvec_dynamic_01;

pull_vector_elements_towards_mean((int) L_static, 0.01);

(*Z_ol_p)=(*Z_ol_p)%kron(ones_col_Ox1, trans(colvec_dynamic_01/Z_colvec_before_01));

}

static colvec::fixed<L_static> Z_colvec_before_02;

void Z_pull_up_vector_elements_generalized_gaussian_wrapper(mat* Z_ol_p){

colvec_dynamic_02=trans(sum(*Z_ol_p, 0));

/*Save*/
Z_colvec_before_02=colvec_dynamic_02;

pull_up_vector_elements_generalized_gaussian(8, 0.6, (int) L_static , 0.15, 0.9);

(*Z_ol_p)=(*Z_ol_p)%kron(ones_col_Ox1, trans(colvec_dynamic_02/Z_colvec_before_02));

}

static rowvec::fixed<L_static> Z_L1_norm_rowvec;

void check_Z_scale_down_Y_low_threshold(mat* Z_ol_p, mat* Y_lk_p){

int l_index;
double L1_norm_threshold=0.00000032875;

Z_L1_norm_rowvec=sum(*Z_ol_p, 0);

for (l_index=0; l_index<L_static; l_index++){

	if (Z_L1_norm_rowvec(l_index)<L1_norm_threshold){

		(*Y_lk_p).row(l_index).zeros();

	}

}

(*Y_lk_p).elem( find_nonfinite((*Y_lk_p)) ).zeros();

// Clip the negative elements again..
(*Y_lk_p).elem(find((*Y_lk_p)<=0)).fill(0.00000001);

(*Y_lk_p).elem(find((*Y_lk_p)>10)).fill(10);

}

static colvec::fixed<L_static> Y_L1_norm_colvec;
static colvec::fixed<L_static> Y_L1_norm_scaling_colvec;

void check_Z_scale_down_Y_high_threshold(mat* Z_ol_p, mat* Y_lk_p){

int l_index;
double L1_norm_threshold=0.003;
double max_L1_norm_in_Y;

/*Input info*/
Z_L1_norm_rowvec=sum(*Z_ol_p, 0);

Z_L1_norm_rowvec.print("Inside check_Z_scale_down_Y_high_threshold, Z_L1_norm_rowvec=");

/*Needs to be scaled elementwise*/
Y_L1_norm_colvec=sum(*Y_lk_p, 1);

max_L1_norm_in_Y=max(Y_L1_norm_colvec);

/*Default thing is to scale the output by 1*/
Y_L1_norm_scaling_colvec.ones();

for (l_index=0; l_index<L_static; l_index++){

	if ((Z_L1_norm_rowvec(l_index))<=L1_norm_threshold){

		Y_L1_norm_scaling_colvec(l_index)=(0.047/L1_norm_threshold)*Z_L1_norm_rowvec(l_index)*(max_L1_norm_in_Y);

		/*To adjust the L1 norm of a vector to a target l1 norm multiply elementwise by the target, divide elementwise by the current L1 norm*/
		(*Y_lk_p).row(l_index)=(Y_L1_norm_scaling_colvec(l_index)/Y_L1_norm_colvec(l_index))*(*Y_lk_p).row(l_index);

	}

}

(*Y_lk_p).elem( find_nonfinite((*Y_lk_p)) ).zeros();

// Clip the negative elements again..
(*Y_lk_p).elem(find((*Y_lk_p)<=0)).fill(0.00000001);

(*Y_lk_p).elem(find((*Y_lk_p)>10)).fill(10);

}