#include "local_inc.hpp"

/*mat::fixed<F_static, N_static> second_dummy_mat_real_FN;
mat::fixed<F_static, N_static> third_dummy_mat_real_FN;*/

/*Function local data: m=1:2 F_staticxK_static matrices*/
/*cx_mat::fixed<F_static, K_static> dummy_mat_H_fk_1;*/
static cx_mat::fixed<F_static, L_static> dummy_mat_H_fl_1_local;

void cxarg_m2_populate_dummy_mat_H_fl_1(cx_mat& dummy_mat_H_fl_1_ref, int m_index, cx_cube* W_fom_cx_p, mat* Z_ol_p){

/*Populate dummy_mat_H_fl*/
dummy_mat_H_fl_1_ref=((*W_fom_cx_p).slice(m_index))*(*Z_ol_p);

}

void cxarg_m2_populate_dummy_mat_H_fk_1(cx_mat& dummy_mat_H_fk_1_ref, int m_index, cx_cube* W_fom_cx_p, mat* Z_ol_p, mat* Y_lk_p){

/*Populate dummy_mat_H_fl*/
dummy_mat_H_fl_1_local=((*W_fom_cx_p).slice(m_index))*(*Z_ol_p);

/*Populate dummy_mat_H_fk*/
dummy_mat_H_fk_1_ref=dummy_mat_H_fl_1_local*(*Y_lk_p);

}

/*Function local data: m=1:2 F_staticxK_static matrices*/
/*cx_mat::fixed<F_static, K_static> dummy_mat_conj_H_fk_1;*/
static cx_mat::fixed<F_static, L_static> dummy_mat_conj_H_fl_1_local;

void cxarg_m2_populate_dummy_mat_conj_H_fl_1(cx_mat& dummy_mat_conj_H_fl_1_ref, int m_index, cx_cube* W_fom_cx_p, mat* Z_ol_p){

/*Populate dummy_mat_H_fl*/
dummy_mat_conj_H_fl_1_ref=conj((*W_fom_cx_p).slice(m_index))*(*Z_ol_p);

}

void cxarg_m2_populate_dummy_mat_conj_H_fk_1(cx_mat& dummy_mat_conj_H_fk_1_ref, int m_index, cx_cube* W_fom_cx_p, mat* Z_ol_p, mat* Y_lk_p){

/*Populate dummy_mat_H_fl*/
dummy_mat_conj_H_fl_1_local=conj((*W_fom_cx_p).slice(m_index))*(*Z_ol_p);

/*Populate dummy_mat_H_fk*/
dummy_mat_conj_H_fk_1_ref=dummy_mat_conj_H_fl_1_local*(*Y_lk_p);

}

/*Already defined in common_module.cpp under module 1*/
/*static cx_mat::fixed<N_static*K_static, F_static> temptensor_NKxF;
static cx_mat::fixed<F_static, N_static*K_static> temptensor_FxNK;
cx_cube::fixed<F_static, N_static, K_static> rotated_expj_Phi_S_FNK;

void rotate_expj_Phi_S_nkf_p_to_FNK(cx_cube* expj_Phi_S_nkf_p){

arrayops::copy(temptensor_NKxF.memptr(), (*expj_Phi_S_nkf_p).memptr(), F_static*N_static*K_static);

temptensor_FxNK=trans(temptensor_NKxF);

arrayops::copy(rotated_expj_Phi_S_FNK.memptr(), temptensor_FxNK.memptr(), F_static*N_static*K_static);

}*/