extern void compute_deltas_init_rowvecs(void);

extern void update_W_update_and_plot_rowvecs(mxArray *plhs[], int n_iter, double delta1, double delta2, double delta3);
extern void update_T_update_and_plot_rowvecs(mxArray *plhs[], int n_iter, double delta1, double delta2, double delta3);
extern void update_V_update_and_plot_rowvecs(mxArray *plhs[], int n_iter, double delta1, double delta2, double delta3);
extern void update_PS_update_and_plot_rowvecs(mxArray *plhs[], int n_iter, double delta1, double delta2, double delta3);
extern void update_X_update_and_plot_rowvecs(mxArray *plhs[], int n_iter, double delta1, double delta2, double delta3);

extern rowvec::fixed<N_ITERATIONS> NMF_auxfun1_rowvec;
extern rowvec::fixed<N_ITERATIONS> NMF_auxfun2_rowvec;
extern rowvec::fixed<N_ITERATIONS> NMF_rowvec_combined;

extern void compute_deltas_NMF_auxfun_rowvecs_send_data_to_Matlab_Eng_and_plot(mxArray *plhs[]);