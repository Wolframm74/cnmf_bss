/*hfk_Af_LS_estimator_update_common*/
extern mat::fixed<O_static, K_static> hfk_Af_LS_estimator_update_quantity_Z_ok;

extern cx_cube::fixed<F_static, K_static, M_static> h_fkm_shared_quantity;

extern cx_cube::fixed<M_static, F_static, K_static> hfk_Af_LS_estimator_update_quantity_h_mfk;
extern cx_cube::fixed<K_static, F_static, N_static> hfk_Af_LS_estimator_update_quantity_s_hat_kfn;
extern cx_cube::fixed<M_static, K_static, F_static> hfk_Af_LS_estimator_update_quantity_A_mkf;

extern cx_cube::fixed<M_static, F_static, O_static> W_mfo_cx_shared;

extern void hfk_Af_LS_estimator_update_compute_Z_ok(mat* Z_ol_p, mat* Y_lk_p);

extern void hfk_Af_LS_estimator_update_compute_h_mfk(cx_cube* W_fom_cx_p);

extern void hfk_Af_LS_estimator_update_compute_s_hat(mat* T_fk_p, mat* V_nk_p, cx_cube* expj_Phi_S_nkf_p);

extern void hfk_Af_LS_estimator_update_compute_A(cx_cube* Xtilde_fnm_p);

extern void rotate_W_fom_cx(cx_cube* W_fom_cx_p);

extern void compute_h_mfk_s_hat_Amkf_common(arg_struct_t* argStruct_p);

extern void compute_J(void);
extern double J_value_secondary_objective;

/*hfk_Af_LS_estimator_update_engine*/
extern void update_engine_hfk_Af_LS_estimator(arg_struct_t* argStruct_p, int which_update);

/*hfk_Af_LS_estimator_Wfom_update*/
extern void* compute_Partial_hfk_Af_Cost_wrt_Wfom_thread_start(void* arg);

/*hfk_Af_LS_estimator_Ylk_update*/
extern void hfk_Af_LS_estimator_Ylk_update_compute_h_mfl(cx_cube* W_fom_cx_p);
extern void* compute_Partial_hfk_Af_Cost_wrt_Ylk_thread_start(void* arg);

/*hfk_Af_LS_estimator_Zol_update*/
extern void hfk_Af_LS_estimator_Zol_update_API_serial_preprocess(mat* Y_lk_p, cx_cube* W_fom_cx_p);
extern void* compute_Partial_hfk_Af_Cost_wrt_Zol_thread_start(void* arg);