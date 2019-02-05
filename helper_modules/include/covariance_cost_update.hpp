/*covariance cost update engine*/
extern void update_engine_covariance_cost(arg_struct_t* argStruct_p, int which_update);

/*covariance common*/
extern cx_cube::fixed<F_static, N_static, K_static> expj_Phi_S_fnk;
extern void* compute_tr_E_E_xw_thread_start(void* arg);

extern cube::fixed<F_static, N_static, O_static> tr_Efn_Efn_xw;
extern cube::fixed<F_static, N_static, O_static> tr_Efn_Efn_wx;

/*covariance z update*/
extern void* compute_Partial_Covariance_Cost_wrt_Zol_thread_start(void* arg);

/*covariance y update*/
extern void* compute_Partial_Covariance_Cost_wrt_Ylk_thread_start(void* arg);