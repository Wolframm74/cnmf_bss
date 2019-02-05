/*AUXILIARY FUNCTION 1:*/

/*V_nk output*/

extern mat::fixed<N_static, K_static> outmat_3rd_real_NK_den;
extern mat::fixed<N_static, K_static> outmat_3rd_real_NK_num;

/*T_fk output*/
extern mat::fixed<F_static, K_static> outmat_3rd_real_FK_den;
extern mat::fixed<F_static, K_static> outmat_3rd_real_FK_num;

/*Z_ol output*/
extern mat::fixed<O_static, L_static> outmat_3rd_real_OL_den;
extern mat::fixed<O_static, L_static> outmat_3rd_real_OL_num;

/*W_fom output*/
extern cube::fixed<F_static, O_static, M_static> numerator_tensor_real_FOM;
extern cube::fixed<F_static, O_static, M_static> outtensor_2nd_real_FOM;	/*denominator tensor*/

/*Y_lk output*/
extern mat::fixed<L_static, K_static> outmat_3rd_real_LK_den;
extern mat::fixed<L_static, K_static> outmat_3rd_real_LK_num;

/*Phi S*/
extern cx_cube::fixed<F_static, K_static, N_static> outtensor_target_cx_FKN;	/*let this one be related to the real tensor*/
extern cx_cube::fixed<N_static, K_static, F_static> outtensor_target_cx_NKF; 

/*AUXILIARY FUNCTION 2: INTERCHANNEL*/

/*V_nk output*/

extern mat::fixed<N_static, K_static> V_complex_argument_m4_num_outmat_NK;
extern mat::fixed<N_static, K_static> V_complex_argument_m4_den_outmat_NK;

/*T_fk output*/
extern mat::fixed<F_static, K_static> T_complex_argument_m4_num_outmat_FK;
extern mat::fixed<F_static, K_static> T_complex_argument_m4_den_outmat_FK;

/*Z_ol output*/
extern mat::fixed<O_static, L_static> Z_complex_argument_m4_num_outmat_OL;
extern mat::fixed<O_static, L_static> Z_complex_argument_m4_den_outmat_OL;

/*W_fom output*/
extern mat::fixed<F_static, O_static> W_channel_a_complex_argument_m4_num_outmat_FO;
extern mat::fixed<F_static, O_static> W_channel_a_complex_argument_m4_den_outmat_FO;

extern mat::fixed<F_static, O_static> W_channel_b_complex_argument_m4_num_outmat_FO;
extern mat::fixed<F_static, O_static> W_channel_b_complex_argument_m4_den_outmat_FO;

/*Y_lk output*/
extern mat::fixed<L_static, K_static> Y_complex_argument_m4_num_outmat_LK;
extern mat::fixed<L_static, K_static> Y_complex_argument_m4_den_outmat_LK;

/*Phi S*/
extern cx_cube::fixed<N_static, K_static, F_static> Phi_S_complex_argument_m4_outmat_NKF;
extern cx_cube::fixed<F_static, K_static, N_static> Phi_S_complex_argument_m4_outmat_FKN;