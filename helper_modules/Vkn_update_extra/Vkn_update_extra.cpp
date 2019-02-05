#include "local_inc.hpp"

static rowvec::fixed<K_static> T_fk_p_l1_norm_rowvec;
colvec::fixed<K_static> V_kn_l1_norm_colvec;


//This should give a small boost to the rows of V_kn that have higher L1 norm. The max_L1_norm_in_V value should get a max boost of 1.2 and every other kth L1 norm in V should get a smaller boost relative to 1.2
static void check_significant_V_k_boost_T_k(mat* T_fk_p, double max_L1_norm_in_V, int k_index){

//linear/affine scale_factor generating function
double scale_factor1, scale_factor2;
double M_value;
double offset=1;

//erf() scale_factor generating function
double mu;
double sigma;

//M_value=(0.2/max_L1_norm_in_V);

//M_value=0.2/exp(max_L1_norm_in_V);
M_value=0.02/exp(max_L1_norm_in_V);

//scale_factor=M_value*(V_kn_l1_norm_colvec(k_index))+offset;

scale_factor1=M_value*exp(V_kn_l1_norm_colvec(k_index))+offset;

mu=0.15*max_L1_norm_in_V;
sigma=0.25*mu;

//pass the scale factor again through a gaussian CDF erf function
scale_factor2=(0.5*(1+erf((V_kn_l1_norm_colvec(k_index)-mu)/(sigma*sqrt(2)))));

(*T_fk_p).col(k_index)=scale_factor1*scale_factor2*(*T_fk_p).col(k_index);

}


static void check_significant_V_k_boost_T_k_generalized_gaussian(mat* T_fk_p, double max_L1_norm_in_V, int k_index){

double mu_value, alpha_value;
double mu, alpha, beta; 

double scale_factor;

beta=8;

mu_value=0.6;

alpha_value=(1-mu_value)/2;

mu=mu_value*max_L1_norm_in_V;

alpha=alpha_value*max_L1_norm_in_V;

scale_factor=0.005*exp(  -pow((abs(V_kn_l1_norm_colvec(k_index)-mu))/alpha, beta)  )+1;

(*T_fk_p).col(k_index)=scale_factor*(*T_fk_p).col(k_index);

}	

void check_V_scale_down_T(mat* T_fk_p, mat* V_nk_p){

int k_index;

double T_kth_l1_norm, V_kth_l1_norm, T_kth_target_l1_norm, T_kth_target_l1_norm_final;

double max_L1_norm_in_V;

T_fk_p_l1_norm_rowvec=sum(*T_fk_p, 0);

V_kn_l1_norm_colvec=sum(trans(*V_nk_p), 1);

max_L1_norm_in_V=max(V_kn_l1_norm_colvec);

for (k_index=0; k_index<K_static; k_index++){

V_kth_l1_norm=V_kn_l1_norm_colvec(k_index);

T_kth_l1_norm=T_fk_p_l1_norm_rowvec(k_index);

if (V_kth_l1_norm<=0.3){

T_kth_target_l1_norm=(2000/0.3)*V_kth_l1_norm;

T_kth_target_l1_norm_final=min(T_kth_target_l1_norm, T_kth_l1_norm);

(*T_fk_p).col(k_index)=(T_kth_target_l1_norm/T_kth_l1_norm)*(*T_fk_p).col(k_index);

}

if ((V_kth_l1_norm>=0.3)&&(V_kth_l1_norm<=0.6)){

/*For the adjacent interval, double the slope to as allow a proportionally double amount of freedom in terms of scaling*/
T_kth_target_l1_norm=2*(2000/0.3)*V_kth_l1_norm;

T_kth_target_l1_norm_final=min(T_kth_target_l1_norm, T_kth_l1_norm);

(*T_fk_p).col(k_index)=(T_kth_target_l1_norm/T_kth_l1_norm)*(*T_fk_p).col(k_index);

}

/*check_significant_V_k_boost_T_k(T_fk_p, max_L1_norm_in_V, k_index);*/

check_significant_V_k_boost_T_k_generalized_gaussian(T_fk_p, max_L1_norm_in_V, k_index);

}



}

/*void check_V_scale_down_Y(mat* Y_lk_p, mat* V_nk_p){


}*/

/*only normalize the columns of Y_lk that are of "interest". The cols are indexed by k=1:K. Only normalize the k such that the L1 norm of the kth row of V_kn is above some L1 norm threshold. These are the k's of interest. */
/*Need to have kept a record of this from during when you most recently computed the L1 norm for all rows of V_kn. */
/*void normalize_Ylk_second(){




}*/


/*void V_nk_erf_threshold(mat* V_nk_p){

double mu=0.02;
double sigma=0.005;
double scale_factor;

int n_index;
int k_index=0;

for (n_index=0; n_index<N_static; n_index++){

	for (n_index=0; n_index<N_static; n_index++){

	scale_factor=0.5*(1+erf(((*V_nk_p)(n_index, k_index)-mu)/(sigma*sqrt(2))));

	(*V_nk_p)(n_index, k_index)=scale_factor*((*V_nk_p)(n_index, k_index));

	}

}

}*/

/*void V_nk_erf_threshold(mat* V_nk_p){

double mu=0.02;
double sigma=0.005;
double scale_factor;

int n_index;
int k_index;

for (n_index=0; n_index<N_static; n_index++){

	for (n_index=0; n_index<N_static; n_index++){

	scale_factor=0.5*(1+erf(((*V_nk_p)(n_index, k_index)-mu)/(sigma*sqrt(2))));

	(*V_nk_p)(n_index, k_index)=scale_factor*((*V_nk_p)(n_index, k_index));

	}

}

}*/

void V_nk_erf_threshold(mat* V_nk_p){

double mu=0.002;
double sigma=0.0005;
double scale_factor;

int n_index;
int k_index;

for (n_index=0; n_index<N_static; n_index++){

	for (k_index=0; k_index<K_static; k_index++){

	scale_factor=0.5*(1+erf(((*V_nk_p)(n_index, k_index)-mu)/(sigma*sqrt(2))));

	(*V_nk_p)(n_index, k_index)=scale_factor*((*V_nk_p)(n_index, k_index));

	}

}

}