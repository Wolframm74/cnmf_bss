#define NUM_WORKER_THREADS_COMPLEX_ARGUMENT_COSTFUN_M4 4

extern void update_engine_complex_argument_costfun_m4(arg_struct_t* argStruct_p, int which_update, int m_index_b, int m_index_a);

extern int complex_argument_costfun_m4_m_index_b;
extern int complex_argument_costfun_m4_m_index_a;

extern cx_mat::fixed<F_static, N_static> complex_argument_costfun_m4_Error_mat_FN;
extern mat::fixed<F_static, N_static> complex_argument_costfun_m4_Magnitude_model_mat_FN;

/*Thread start functions*/
extern void* Phi_S_complex_argument_m4_start(void* arg);
extern void* V_complex_argument_m4_start(void* arg);
extern void* T_complex_argument_m4_start(void* arg);
extern void* W_channel_a_complex_argument_m4_start(void* arg);
extern void* W_channel_b_complex_argument_m4_start(void* arg);
extern void* Z_complex_argument_m4_start(void* arg);
extern void* Y_complex_argument_m4_start(void* arg);