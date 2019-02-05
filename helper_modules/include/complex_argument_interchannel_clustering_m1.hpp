extern void complex_argument_interchannel_clustering_m1_entry(mxArray *plhs[], arg_struct_t* argStruct_p, int m_index_b, int m_index_a);
extern void complex_argument_interchannel_clustering_m1_lookdirs_search_entry(mxArray *plhs[], arg_struct_t* argStruct_p, int m_index_b, int m_index_a);

#define Q_static L_static

extern mat::fixed<F_static, Q_static> chat_fq_estimate;
extern mat::fixed<Q_static, N_static> indicators_P_qn;

extern void complex_argument_interchannel_clustering_m1_lookdirs_search_populate_ones_expj_Phi_S_fkn(cx_cube* expj_Phi_S_fkn_p, cx_cube* expj_Phi_S_nkf_p);
extern void complex_argument_interchannel_clustering_m1_lookdirs_search_boost_and_suppress_V_kn(mat* V_nk_p); 

extern void send_data_to_Matlab_Eng_and_plot(mxArray *plhs[]);