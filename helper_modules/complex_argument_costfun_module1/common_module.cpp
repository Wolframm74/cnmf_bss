#include "local_inc.hpp"

mat::fixed<F_static, N_static> second_dummy_mat_real_FN;
mat::fixed<F_static, N_static> third_dummy_mat_real_FN;

colvec::fixed<1> scalar_global_real;
cx_colvec::fixed<1> scalar_global_cx;

/*Ones stuff*/
mat::fixed<F_static, N_static> ones_mat_FxN;
cube::fixed<F_static, N_static, K_static> ones_cube_FNK;	/*Inside Phi_S_update2.cpp*/
rowvec::fixed<N_static> ones_row_1xN; /*Phi_S_update2.cpp*/
mat::fixed<F_static, K_static> ones_cube_FK;

/*Function local data: m=1:2 F_staticxK_static matrices*/
cx_mat::fixed<F_static, K_static> dummy_mat_H_fk_1;

cx_mat::fixed<F_static, L_static> dummy_mat_H_fl_1;

void populate_dummy_mat_H_fl_1(int m_index, cx_cube* W_fom_cx_p, mat* Z_ol_p){

/*Populate dummy_mat_H_fl*/
dummy_mat_H_fl_1=((*W_fom_cx_p).slice(m_index))*(*Z_ol_p);

}

void populate_dummy_mat_H_fk_1(int m_index, cx_cube* W_fom_cx_p, mat* Z_ol_p, mat* Y_lk_p){

/*Populate dummy_mat_H_fl*/
dummy_mat_H_fl_1=((*W_fom_cx_p).slice(m_index))*(*Z_ol_p);

/*Populate dummy_mat_H_fk*/
dummy_mat_H_fk_1=dummy_mat_H_fl_1*(*Y_lk_p);

/*dummy_mat_H_fk_1.print("dummy_mat_H_fk_1: Phi_S_update2, common_modue: line 38:");*/

}

cx_mat::fixed<F_static, N_static> dummy_mat_cx_FN_common;

mat::fixed<F_static, N_static> outmat_FN;
mat::fixed<F_static, N_static> nummat_FN;
mat::fixed<F_static, N_static> denmat_FN;

/*Function local data*/
mat& compute_arg_X_FN(cx_cube* X_fnm, int m_index){

int f_index, n_index; 

mat& outmat_FN_ref=outmat_FN;

dummy_mat_cx_FN_common=(*X_fnm).slice(m_index);

dummy_mat_cx_FN_common.elem(find(abs(dummy_mat_cx_FN_common)==0)).fill(0.00000001);

nummat_FN=imag(dummy_mat_cx_FN_common);

denmat_FN=sqrt(square(real(dummy_mat_cx_FN_common))+square(imag(dummy_mat_cx_FN_common)))+real(dummy_mat_cx_FN_common);

outmat_FN=nummat_FN/denmat_FN;

outmat_FN=2*atan(outmat_FN);

/*find nan and set to zero. Everything else should be non-nan */
outmat_FN.elem( find_nonfinite(outmat_FN )).zeros();

return outmat_FN_ref;

}

/*Function local data*/

mat::fixed<F_static, N_static> G1_outmat_fn;
mat::fixed<F_static, N_static> H1_outmat_fn;

mat::fixed<F_static, N_static> sqrt_denmat_1_fn;

void populate_G_H_fn_channel1(int m_index, cx_cube* Xhat_fnm_p){

H1_outmat_fn=real((*Xhat_fnm_p).slice(m_index));

G1_outmat_fn=imag((*Xhat_fnm_p).slice(m_index));

/*Compute  for the denominator function*/
sqrt_denmat_1_fn=sqrt(square(G1_outmat_fn)+square(H1_outmat_fn));
/*May use the above result again later*/

H1_outmat_fn=sqrt_denmat_1_fn+H1_outmat_fn;

}


cx_mat::fixed<N_static*K_static, F_static> temptensor_NKxF;
cx_mat::fixed<F_static, N_static*K_static> temptensor_FxNK;
cx_cube::fixed<F_static, N_static, K_static> rotated_expj_Phi_S_FNK;

void rotate_expj_Phi_S_nkf_p_to_FNK(cx_cube* expj_Phi_S_nkf_p){

/*copy to temptensor_NKxF*/
arrayops::copy(temptensor_NKxF.memptr(), (*expj_Phi_S_nkf_p).memptr(), F_static*N_static*K_static);

/*transpose temptensor_NKxF into output result FxNK*/	
temptensor_FxNK=trans(temptensor_NKxF);

/*copy result to rotated_expj_Phi_S_FNK*/
arrayops::copy(rotated_expj_Phi_S_FNK.memptr(), temptensor_FxNK.memptr(), F_static*N_static*K_static);

}