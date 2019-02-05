#include "local_inc.hpp"

/*mat::fixed<O_static, K_static> Zok_sparsity_common_sharedmat_Zok; 
rowvec::fixed<K_static> Q1_1xK;
rowvec::fixed<K_static> Q2_1xK;
mat::fixed<O_static, K_static> Q3_OxK;
rowvec::fixed<K_static> Q4_1xK;
rowvec::fixed<K_static> Q5_1xK;
rowvec::fixed<K_static> Q6_1xK;
mat::fixed<O_static, K_static> Q7_OxK;

extern colvec::fixed<K_static> sigma_vec;

*/

mat::fixed<L_static, K_static> pve_part_Partial_wrt_Ylk;
mat::fixed<L_static, K_static> nve_part_Partial_wrt_Ylk;

mat::fixed<L_static, K_static> outmat1_Q3_Zol_LxK;
mat::fixed<L_static, K_static> outmat2_Zok_Zol_LxK;

static void compute_outmat1(mat* Z_ol_p){

outmat1_Q3_Zol_LxK=trans(*Z_ol_p)*Q3_OxK;

}

static void compute_outmat2(mat* Z_ol_p){

outmat2_Zok_Zol_LxK=trans(*Z_ol_p)*Zok_sparsity_common_sharedmat_Zok;

}

void Zok_sparsity_Ylk_update(mat* Z_ol_p){

int k_index;

compute_outmat1(Z_ol_p);

compute_outmat2(Z_ol_p);

for (k_index=0; k_index<K_static; k_index++){

pve_part_Partial_wrt_Ylk.col(k_index)=Zok_sparsity_lamda*Q1_1xK(k_index)*(outmat1_Q3_Zol_LxK.col(k_index)*sigma_vec(k_index)*Q2_1xK(k_index)+(outmat2_Zok_Zol_LxK.col(k_index))*((Q4_1xK(k_index)*Q5_1xK(k_index))/(Q2_1xK(k_index))));

nve_part_Partial_wrt_Ylk.col(k_index)=Zok_sparsity_lamda*Q1_1xK(k_index)*((outmat2_Zok_Zol_LxK.col(k_index))*sigma_vec(k_index)*((Q5_1xK(k_index))/(Q2_1xK(k_index)))+outmat1_Q3_Zol_LxK.col(k_index)*(Q4_1xK(k_index)*Q2_1xK(k_index)));

}

}