extern cx_cube::fixed<F_static, L_static, N_static> expj_Phi_U_flnm[M_static];

extern void Xhat_serial_function(mat* Z_ol_p, cx_cube* W_fom_cx_p, cube* W_fom_p);

/*Xhat_update*/
extern cx_cube::fixed<F_static, L_static, M_static> Xhat_outtensor_cx_FLM; 
extern cube::fixed<F_static, L_static, M_static> Xhat_outtensor_real_FLM;

extern cx_cube::fixed<F_static, K_static, N_static> Xhat_out_4way_tensor_cx_fknm[M_static];
extern cx_cube::fixed<N_static, K_static, F_static> Xhat_out_4way_tensor_cx_nkfm[M_static];

extern cx_cube::fixed<F_static, O_static, M_static> Xhat_outtensor_cx_FOM; 
extern cx_cube::fixed<F_static, O_static, M_static> Xhat_dummytensor_cx_FOM[NUM_WORKER_THREADS]; 

extern bool Xhat_pop_NKFM_for_T_or_Phi_S_flag;
extern bool Xhat_pop_FKNM_for_V_flag;
extern bool Xhat_pop_W_cx_outtensor_flag;

/*Phi_S update*/
extern void Phi_S_update_auxfun1_zero_internal_quantities(void);