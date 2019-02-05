#include "local_inc.hpp"

mat::fixed<L_static, K_static> Zok_sparsity_outmat_LK_den;

void Zok_Ylk_update_extra_sparsity(mat* Z_ol_p, mat* Y_lk_p){

compute_Zok(Z_ol_p, Y_lk_p);

double exponent_value;

exponent_value=(1)+(((double)PVALUE_FIXED)-2);

Zok_sparsity_outmat_LK_den=(accu_Zok_square*10000)*((double)PVALUE_FIXED)*trans(*Z_ol_p)*(pow(Zok_shared_mat, exponent_value));

}