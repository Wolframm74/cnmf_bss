#include "local_inc.hpp"

#define L_target 3

static colvec::fixed<L_static> Y_colvec_before_02;

void Y_pull_up_vector_elements_generalized_gaussian_wrapper(mat* Y_lk_p){

colvec_dynamic_02=(sum(*Y_lk_p, 1));

/*Save*/
Y_colvec_before_02=colvec_dynamic_02;

pull_up_vector_elements_generalized_gaussian(8, 0.5, (int) L_static , 0.05, 1);

(*Y_lk_p)=(*Y_lk_p)%kron(ones_row_1xK, (colvec_dynamic_02/Y_colvec_before_02));

}

/*Represent the L1 norms*/
static rowvec::fixed<L_static> Z_L1_norm_rowvec;
static colvec::fixed<L_static> Y_L1_norm_colvec;

/*Scaled relative to their respective L1 norms*/
static rowvec::fixed<L_static> Z_L1_norm_rowvec_scaledvec;
static colvec::fixed<L_static> Y_L1_norm_colvec_scaledvec;

static colvec::fixed<L_static> dot_outvec;

void relate_significant_Zl_clusters_to_Yl_Clusters(mat* Z_ol_p, mat* Y_lk_p){

(*Y_lk_p).print("relate_significant_Zl_clusters_to_Yl_Clusters, Y_lk, before:");

double sum_L1_norms_Z;

double sum_L1_norms_Y;

Z_L1_norm_rowvec=sum(*Z_ol_p, 0);

Y_L1_norm_colvec=sum(*Y_lk_p, 1);

sum_L1_norms_Z=as_scalar(sum(Z_L1_norm_rowvec));

sum_L1_norms_Y=as_scalar(sum(Y_L1_norm_colvec));

Z_L1_norm_rowvec_scaledvec=Z_L1_norm_rowvec/(sum_L1_norms_Z/((double)L_static));

Y_L1_norm_colvec_scaledvec=Y_L1_norm_colvec/(sum_L1_norms_Y/((double)L_static));

dot_outvec=trans(Z_L1_norm_rowvec_scaledvec)%Y_L1_norm_colvec_scaledvec;

dot_outvec.print("dot_outvec, before:");

/*normalize dot_outvec so that it sums to sum_L1_norms_Y*/

dot_outvec=dot_outvec*(sum_L1_norms_Y/as_scalar(sum(dot_outvec)));

(*Y_lk_p)=(*Y_lk_p)%kron(dot_outvec/Y_L1_norm_colvec, ones_row_1xK);

Z_L1_norm_rowvec_scaledvec.print("Z_L1_norm_rowvec_scaledvec");

Y_L1_norm_colvec_scaledvec.print("Y_L1_norm_colvec_scaledvec");

dot_outvec.print("dot_outvec");

(*Y_lk_p).print("relate_significant_Zl_clusters_to_Yl_Clusters, Y_lk, after:");

}

/*static mat::fixed<L_target, 2> Ltargetx2_info_colmatrix;*/
static colvec::fixed<L_target> Ltarget_amplitudes_col;
static uvec::fixed<L_target> Ltarget_indices_col;
static colvec::fixed<L_static> Y_L1_norm_colvec_secondary;

void balance_energies_of_Ltarget_most_significant_clusters(mat* Y_lk_p){

int l_ctr;
double returned_index;

double total_energy, avg_energy;

Y_L1_norm_colvec=sum(*Y_lk_p, 1);
Y_L1_norm_colvec_secondary=Y_L1_norm_colvec;


for (l_ctr=0; l_ctr<L_target; l_ctr++){

Ltarget_amplitudes_col(l_ctr)=Y_L1_norm_colvec.max(Ltarget_indices_col(l_ctr));

/*Ltarget_amplitudes_col(l_ctr, 2)=returned_index;*/

Y_L1_norm_colvec(Ltarget_indices_col(l_ctr))=0;

}

total_energy=as_scalar(sum(Ltarget_amplitudes_col));

avg_energy=total_energy/((double)L_target);

Y_L1_norm_colvec.elem(Ltarget_indices_col).fill(avg_energy);

/*Kronecker spreading*/
(*Y_lk_p)=(*Y_lk_p)%kron(Y_L1_norm_colvec/Y_L1_norm_colvec_secondary, ones_row_1xK);

}