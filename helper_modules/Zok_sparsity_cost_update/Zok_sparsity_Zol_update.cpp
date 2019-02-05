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

mat::fixed<O_static, L_static> pve_part_Partial_wrt_Zol;
mat::fixed<O_static, L_static> pve_dummy_mat1_ol;
mat::fixed<O_static, L_static> pve_dummy_mat2_ol;

mat::fixed<O_static, L_static> nve_part_Partial_wrt_Zol;
mat::fixed<O_static, L_static> nve_dummy_mat1_ol;
mat::fixed<O_static, L_static> nve_dummy_mat2_ol;


void Zok_sparsity_Zol_update(mat* Y_lk_p){

int k_index;

for (k_index=0; k_index<K_static; k_index++){

pve_dummy_mat1_ol=(Q3_OxK.col(k_index)*trans((*Y_lk_p).col(k_index)))*Q2_1xK(k_index);

pve_dummy_mat2_ol=(Zok_sparsity_common_sharedmat_Zok.col(k_index)*trans((*Y_lk_p).col(k_index)))*((Q5_1xK(k_index))/(Q2_1xK(k_index)));

nve_dummy_mat1_ol=pve_dummy_mat2_ol;

nve_dummy_mat2_ol=pve_dummy_mat1_ol;

pve_dummy_mat1_ol=pve_dummy_mat1_ol*sigma_vec(k_index);

pve_dummy_mat2_ol=pve_dummy_mat2_ol*Q4_1xK(k_index);

nve_dummy_mat1_ol=nve_dummy_mat1_ol*sigma_vec(k_index);

nve_dummy_mat2_ol=nve_dummy_mat2_ol*Q4_1xK(k_index);

pve_part_Partial_wrt_Zol=Zok_sparsity_lamda*Q1_1xK(k_index)*(pve_dummy_mat1_ol+pve_dummy_mat2_ol);

nve_part_Partial_wrt_Zol=Zok_sparsity_lamda*Q1_1xK(k_index)*(nve_dummy_mat1_ol+nve_dummy_mat2_ol);

}

}