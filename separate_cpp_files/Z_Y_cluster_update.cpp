#include "local_inc.hpp"

/*common matrices*/

mat::fixed<L_static, L_static> dummy_mat1_Y_Yt_LL;
mat::fixed<O_static, L_static> ones_mat_OxL;
mat::fixed<L_static, K_static> ones_mat_LxK;

/*Z_cluster_update matrices*/

mat::fixed<O_static, L_static> Z_outmat_num_OL;
mat::fixed<O_static, L_static> Z_outmat_den_OL;

/*Y_cluster_update matrices*/

mat::fixed<L_static, L_static> dummy_mat2_Zt_Z_LL;

mat::fixed<L_static, K_static> Y_outmat_den_LK;
mat::fixed<L_static, K_static> Y_outmat_num_LK;

static void Z_cluster_update(mat* Z_ol_p, mat* Y_lk_p){

/*Compute Y*Y' */
dummy_mat1_Y_Yt_LL=(*Y_lk_p)*trans(*Y_lk_p);

/*Compute Z*Y*Y' */
Z_outmat_den_OL=(*Z_ol_p)*dummy_mat1_Y_Yt_LL;

/*Compute numerator as Z*Y*Y' * Y*Y' */
Z_outmat_num_OL=Z_outmat_den_OL*dummy_mat1_Y_Yt_LL;

ones_mat_OxL.ones();

/*Z_ol=Z_ol%(ones(O, L)+num/den)*/
(*Z_ol_p)=(*Z_ol_p)%(ones_mat_OxL+(Z_outmat_num_OL/Z_outmat_den_OL));

ones_mat_OxL.fill(0.5);

/*Z_ol=0.5*Z_ol;*/	
(*Z_ol_p)=ones_mat_OxL%(*Z_ol_p);

}

static void Y_cluster_update(mat* Y_lk_p, mat* Z_ol_p){

/*Compute Z'Z*/
dummy_mat2_Zt_Z_LL=trans(*Z_ol_p)*(*Z_ol_p);

/*Compute Z'Z*Y */	
Y_outmat_den_LK=dummy_mat2_Zt_Z_LL*(*Y_lk_p);

/*Compute YY' */	/*Maybe don't even need to re-compute this assuming you call Z first (always), before Y*/
dummy_mat1_Y_Yt_LL=(*Y_lk_p)*trans(*Y_lk_p);

/*Compute numerator term 1*/	
Y_outmat_num_LK=dummy_mat1_Y_Yt_LL*Y_outmat_den_LK;

/*Numerator = numerator term 1 + numerator term 2 */	
Y_outmat_num_LK=Y_outmat_num_LK+(dummy_mat2_Zt_Z_LL*dummy_mat1_Y_Yt_LL*(*Y_lk_p));

ones_mat_LxK.ones();

/*Y_lk=Y_lk%(ones(L, K)+num/den)*/
(*Y_lk_p)=(*Y_lk_p)%(ones_mat_LxK+(Y_outmat_num_LK/Y_outmat_den_LK));

ones_mat_LxK.fill(0.5);

/*Y_lk=0.5*Y_lk;*/		
(*Y_lk_p)=ones_mat_LxK%(*Y_lk_p);

}

/*Expose this function to outside world*/
void Z_cluster_update_wrapper(arg_struct_t* argStruct_p){

double accum_Z_ol_local; 

Z_cluster_update(argStruct_p->Z_ol_p, argStruct_p->Y_lk_p);

// Subtract off the current mean
accum_Z_ol_local=accu(*(argStruct_p->Z_ol_p));

*(argStruct_p->Z_ol_p)=(*(argStruct_p->Z_ol_p))-(accum_Z_ol_local/(((double)O_static)*((double)L_static)))*ones(O_static, L_static);

// Set new mean forcibly
*(argStruct_p->Z_ol_p)=(*(argStruct_p->Z_ol_p))+(1)*(accum_Z_ol/(((double)O_static)*((double)L_static)))*ones(O_static, L_static);

/*Optionally clip the negative elements*/
(*(argStruct_p->Z_ol_p)).elem(find((*(argStruct_p->Z_ol_p))<=0)).fill(0.00000001);

}

/*Expose this function to outside world*/
void Y_cluster_update_wrapper(arg_struct_t* argStruct_p){

Y_cluster_update(argStruct_p->Y_lk_p, argStruct_p->Z_ol_p);

/*Optionally clip the negative elements*/
(*(argStruct_p->Y_lk_p)).elem(find((*(argStruct_p->Y_lk_p))<=0)).fill(0.00000001);

}