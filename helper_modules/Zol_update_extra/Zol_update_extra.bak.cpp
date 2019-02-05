#include "local_inc.hpp"

static rowvec::fixed<K_static> Z_ol_p_l1_norm_rowvec;

static void check_significant_Z_l_pull_up_weak_clusters_generalized_gaussian(mat* Z_ol_p, double max_L1_norm_in_Z, int l_index){

double mu_value, alpha_value;
double mu, alpha, beta; 

double scale_factor;

beta=8;

mu_value=0.6;

alpha_value=(1-mu_value)/2;

mu=mu_value*max_L1_norm_in_Z;

alpha=alpha_value*max_L1_norm_in_Z;

scale_factor=0.15*exp(  -pow((abs(Z_ol_p_l1_norm_rowvec(l_index)-mu))/alpha, beta)  )+0.9;

(*Z_ol_p).col(l_index)=scale_factor*(*Z_ol_p).col(l_index);

}	

void check_Z_pull_up_weak_clusters(mat* Z_ol_p){

int l_index;

double Z_lth_l1_norm, V_kth_l1_norm, T_kth_target_l1_norm, T_kth_target_l1_norm_final;

double max_L1_norm_in_Z;

Z_ol_p_l1_norm_rowvec=sum(*Z_ol_p, 0);

max_L1_norm_in_Z=max(Z_ol_p_l1_norm_rowvec);

for (l_index=0; l_index<L_static; l_index++){

/*check_significant_V_k_boost_T_k(Z_ol_p, max_L1_norm_in_Z, l_index);*/

check_significant_Z_l_pull_up_weak_clusters_generalized_gaussian(Z_ol_p, max_L1_norm_in_Z, l_index);

}

}
