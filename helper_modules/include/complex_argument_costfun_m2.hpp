extern void complex_argument_costfun_m2_V_update3_entry(arg_struct_t* argStruct_p);
extern void complex_argument_costfun_m2_Phi_S_update3_entry(arg_struct_t* argStruct_p);
extern void complex_argument_costfun_m2_W_update3_entry(arg_struct_t* argStruct_p);

extern cx_mat::fixed<F_static, N_static> G_ba_FxN_cx_mat; 

/*common module*/
extern void cxarg_m2_populate_dummy_mat_H_fl_1(cx_mat& dummy_mat_H_fl_1_ref, int m_index, cx_cube* W_fom_cx_p, mat* Z_ol_p);
extern void cxarg_m2_populate_dummy_mat_H_fk_1(cx_mat& dummy_mat_H_fk_1_ref, int m_index, cx_cube* W_fom_cx_p, mat* Z_ol_p, mat* Y_lk_p);
extern void cxarg_m2_populate_dummy_mat_conj_H_fl_1(cx_mat& dummy_mat_conj_H_fl_1_ref, int m_index, cx_cube* W_fom_cx_p, mat* Z_ol_p);
extern void cxarg_m2_populate_dummy_mat_conj_H_fk_1(cx_mat& dummy_mat_conj_H_fk_1_ref, int m_index, cx_cube* W_fom_cx_p, mat* Z_ol_p, mat* Y_lk_p);

/*Rotating function already defined in complex_argument_costfun_m1.hpp*/
/*extern cx_cube::fixed<F_static, N_static, K_static> rotated_expj_Phi_S_FNK;
extern void rotate_expj_Phi_S_nkf_p_to_FNK(cx_cube* expj_Phi_S_nkf_p);*/