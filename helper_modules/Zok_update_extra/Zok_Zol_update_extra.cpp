#include "local_inc.hpp"

mat::fixed<O_static, K_static> Zok_shared_mat;

double accu_Zok_square;

void compute_Zok(mat* Z_ol_p, mat* Y_lk_p){

	Zok_shared_mat=(*Z_ol_p)*(*Y_lk_p);

	accu_Zok_square=accu(Zok_shared_mat%Zok_shared_mat);

}

mat::fixed<O_static, L_static> Zok_sparsity_outmat_OL_den;


void Zok_Zol_update_extra_sparsity(mat* Z_ol_p, mat* Y_lk_p){

double exponent_value;

exponent_value=(1)+(((double)PVALUE_FIXED)-2);

Zok_sparsity_outmat_OL_den=(accu_Zok_square*10000)*((double)PVALUE_FIXED)*pow(Zok_shared_mat, exponent_value)*trans(*Y_lk_p);

}