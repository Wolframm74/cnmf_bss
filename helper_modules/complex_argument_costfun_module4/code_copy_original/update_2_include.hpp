/*----------------------------------------------------------------------------------------------------------------------------------------------------------------*/

/*Exposed Function local data*/

extern cx_mat::fixed<F_static, K_static> dummy_mat_H_fk_1;
extern cx_mat::fixed<F_static, L_static> dummy_mat_H_fl_1;

/*Exposed Function(s)*/

extern void populate_dummy_mat_H_fl_1(int m_index, cx_cube* W_fom_cx_p, mat* Z_ol_p);
extern void populate_dummy_mat_H_fk_1(int m_index, cx_cube* W_fom_cx_p, mat* Z_ol_p, mat* Y_lk_p);

/*----------------------------------------------------------------------------------------------------------------------------------------------------------------*/

/*Exposed Function local data*/

extern mat::fixed<F_static, N_static> outmat_FN;
extern mat::fixed<F_static, N_static> nummat_FN;
extern mat::fixed<F_static, N_static> denmat_FN;

/*Exposed Function(s)*/

extern mat& compute_arg_X_FN(cx_cube* X_fnm, int m_index);

/*----------------------------------------------------------------------------------------------------------------------------------------------------------------*/

/*Exposed Function local data*/

extern mat::fixed<F_static, N_static> G1_outmat_fn;
extern mat::fixed<F_static, N_static> H1_outmat_fn;
extern mat::fixed<F_static, N_static> sqrt_denmat_1_fn;

/*Exposed Function(s)*/

extern void populate_G_H_fn_channel1(int m_index, cx_cube* Xhat_fnm_p);

/*----------------------------------------------------------------------------------------------------------------------------------------------------------------*/

extern cx_cube::fixed<F_static, N_static, K_static> rotated_expj_Phi_S_FNK;

extern void rotate_expj_Phi_S_nkf_p_to_FNK(cx_cube* expj_Phi_S_nkf_p);

/*----------------------------------------------------------------------------------------------------------------------------------------------------------------*/


extern mat::fixed<F_static, N_static> dummy_mat_real_FN;	/*Phi_S_update_2*/
extern mat::fixed<F_static, N_static> second_dummy_mat_real_FN;
extern mat::fixed<F_static, N_static> third_dummy_mat_real_FN;

/*----------------------------------------------------------------------------------------------------------------------------------------------------------------*/

/*Ones objects*/

extern mat::fixed<F_static, N_static> ones_mat_FxN;		/*Z_update_2*/
extern cube::fixed<F_static, N_static, K_static> ones_cube_FNK;	/*Inside Phi_S_update2.cpp*/
extern rowvec::fixed<N_static> ones_row_1xN; /*Phi_S_update2.cpp*/
extern mat::fixed<F_static, K_static> ones_cube_FK;

extern void T_update2_entry(arg_struct_t* argStruct_p);
extern void V_update2_entry(arg_struct_t* argStruct_p);
extern void Y_update2_entry(arg_struct_t* argStruct_p);
extern void Z_update2_entry(arg_struct_t* argStruct_p);

extern colvec::fixed<1> scalar_global_real;
extern cx_colvec::fixed<1> scalar_global_cx;