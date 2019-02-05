#include "local_inc.hpp"

static rowvec::fixed<N_ITERATIONS> V_update_delta_rowvec_auxfun1;
static rowvec::fixed<N_ITERATIONS> T_update_delta_rowvec_auxfun1;
static rowvec::fixed<N_ITERATIONS> PS_update_delta_rowvec_auxfun1;
static rowvec::fixed<N_ITERATIONS> W_update_delta_rowvec_auxfun1;

static rowvec::fixed<N_ITERATIONS> X_update_delta_rowvec_auxfun1;

static rowvec::fixed<N_ITERATIONS> V_update_delta_rowvec_auxfun2;
static rowvec::fixed<N_ITERATIONS> T_update_delta_rowvec_auxfun2;
static rowvec::fixed<N_ITERATIONS> PS_update_delta_rowvec_auxfun2;
static rowvec::fixed<N_ITERATIONS> W_update_delta_rowvec_auxfun2;

static rowvec::fixed<N_ITERATIONS> X_update_delta_rowvec_auxfun2;

static rowvec::fixed<N_ITERATIONS> V_update_delta_rowvec_combined;
static rowvec::fixed<N_ITERATIONS> T_update_delta_rowvec_combined;
static rowvec::fixed<N_ITERATIONS> PS_update_delta_rowvec_combined;
static rowvec::fixed<N_ITERATIONS> W_update_delta_rowvec_combined;

static rowvec::fixed<N_ITERATIONS> X_update_delta_rowvec_combined;

rowvec::fixed<N_ITERATIONS> NMF_auxfun1_rowvec;
rowvec::fixed<N_ITERATIONS> NMF_auxfun2_rowvec;
rowvec::fixed<N_ITERATIONS> NMF_rowvec_combined;

void compute_deltas_init_rowvecs(void){

V_update_delta_rowvec_auxfun1.zeros();
T_update_delta_rowvec_auxfun1.zeros();
PS_update_delta_rowvec_auxfun1.zeros();
W_update_delta_rowvec_auxfun1.zeros();

V_update_delta_rowvec_auxfun2.zeros();
T_update_delta_rowvec_auxfun2.zeros();
PS_update_delta_rowvec_auxfun2.zeros();
W_update_delta_rowvec_auxfun2.zeros();

V_update_delta_rowvec_combined.zeros();
T_update_delta_rowvec_combined.zeros();
PS_update_delta_rowvec_combined.zeros();
W_update_delta_rowvec_combined.zeros();

NMF_auxfun1_rowvec.zeros();
NMF_auxfun2_rowvec.zeros();
NMF_rowvec_combined.zeros();

}

static rowvec::fixed<N_ITERATIONS> dummy_rowvec_auxfun1;

static bool compute_deltas_send_data_to_Matlab_Eng_and_plot_init_flag=false;

static void compute_deltas_send_data_to_Matlab_Eng_and_plot_init(mxArray *plhs[]){

plhs[12]=armaCreateMxMatrix(1, dummy_rowvec_auxfun1.n_cols, mxDOUBLE_CLASS, mxREAL);

plhs[13]=armaCreateMxMatrix(1, dummy_rowvec_auxfun1.n_cols, mxDOUBLE_CLASS, mxREAL);

plhs[14]=armaCreateMxMatrix(1, dummy_rowvec_auxfun1.n_cols, mxDOUBLE_CLASS, mxREAL);

compute_deltas_send_data_to_Matlab_Eng_and_plot_init_flag=true;

}

static void compute_deltas_W_send_data_to_Matlab_Eng_and_plot(mxArray *plhs[]){

if (!compute_deltas_send_data_to_Matlab_Eng_and_plot_init_flag){

compute_deltas_send_data_to_Matlab_Eng_and_plot_init(plhs);

}

armaSetPr(plhs[12], W_update_delta_rowvec_auxfun1);

armaSetPr(plhs[13], W_update_delta_rowvec_auxfun2);

armaSetPr(plhs[14], W_update_delta_rowvec_combined);

engPutVariable(mlEngine_p, "W_update_delta_rowvec_auxfun1", plhs[12]);

engPutVariable(mlEngine_p, "W_update_delta_rowvec_auxfun2", plhs[13]);

engPutVariable(mlEngine_p, "W_update_delta_rowvec_combined", plhs[14]);

engEvalString(mlEngine_p, "plot_W_update_deltas(W_update_delta_rowvec_auxfun1, W_update_delta_rowvec_auxfun2, W_update_delta_rowvec_combined);");

}

void update_W_update_and_plot_rowvecs(mxArray *plhs[], int n_iter, double delta1, double delta2, double delta3){

W_update_delta_rowvec_auxfun1(n_iter)=delta1;

W_update_delta_rowvec_auxfun2(n_iter)=delta2;

W_update_delta_rowvec_combined(n_iter)=delta3;

compute_deltas_W_send_data_to_Matlab_Eng_and_plot(plhs);

}

static void compute_deltas_T_send_data_to_Matlab_Eng_and_plot(mxArray *plhs[]){

if (!compute_deltas_send_data_to_Matlab_Eng_and_plot_init_flag){

compute_deltas_send_data_to_Matlab_Eng_and_plot_init(plhs);

}

armaSetPr(plhs[12], T_update_delta_rowvec_auxfun1);

armaSetPr(plhs[13], T_update_delta_rowvec_auxfun2);

armaSetPr(plhs[14], T_update_delta_rowvec_combined);

engPutVariable(mlEngine_p, "T_update_delta_rowvec_auxfun1", plhs[12]);

engPutVariable(mlEngine_p, "T_update_delta_rowvec_auxfun2", plhs[13]);

engPutVariable(mlEngine_p, "T_update_delta_rowvec_combined", plhs[14]);

engEvalString(mlEngine_p, "plot_T_update_deltas(T_update_delta_rowvec_auxfun1, T_update_delta_rowvec_auxfun2, T_update_delta_rowvec_combined);");

}

void update_T_update_and_plot_rowvecs(mxArray *plhs[], int n_iter, double delta1, double delta2, double delta3){

T_update_delta_rowvec_auxfun1(n_iter)=delta1;

T_update_delta_rowvec_auxfun2(n_iter)=delta2;

T_update_delta_rowvec_combined(n_iter)=delta3;

compute_deltas_T_send_data_to_Matlab_Eng_and_plot(plhs);
	
}

static void compute_deltas_V_send_data_to_Matlab_Eng_and_plot(mxArray *plhs[]){

if (!compute_deltas_send_data_to_Matlab_Eng_and_plot_init_flag){

compute_deltas_send_data_to_Matlab_Eng_and_plot_init(plhs);

}

armaSetPr(plhs[12], V_update_delta_rowvec_auxfun1);

armaSetPr(plhs[13], V_update_delta_rowvec_auxfun2);

armaSetPr(plhs[14], V_update_delta_rowvec_combined);

engPutVariable(mlEngine_p, "V_update_delta_rowvec_auxfun1", plhs[12]);

engPutVariable(mlEngine_p, "V_update_delta_rowvec_auxfun2", plhs[13]);

engPutVariable(mlEngine_p, "V_update_delta_rowvec_combined", plhs[14]);

engEvalString(mlEngine_p, "plot_V_update_deltas(V_update_delta_rowvec_auxfun1, V_update_delta_rowvec_auxfun2, V_update_delta_rowvec_combined);");

}

void update_V_update_and_plot_rowvecs(mxArray *plhs[], int n_iter, double delta1, double delta2, double delta3){

V_update_delta_rowvec_auxfun1(n_iter)=delta1;

V_update_delta_rowvec_auxfun2(n_iter)=delta2;

V_update_delta_rowvec_combined(n_iter)=delta3;

compute_deltas_V_send_data_to_Matlab_Eng_and_plot(plhs);
	
}

static void compute_deltas_PS_send_data_to_Matlab_Eng_and_plot(mxArray *plhs[]){

if (!compute_deltas_send_data_to_Matlab_Eng_and_plot_init_flag){

compute_deltas_send_data_to_Matlab_Eng_and_plot_init(plhs);

}

armaSetPr(plhs[12], PS_update_delta_rowvec_auxfun1);

armaSetPr(plhs[13], PS_update_delta_rowvec_auxfun2);

armaSetPr(plhs[14], PS_update_delta_rowvec_combined);

engPutVariable(mlEngine_p, "PS_update_delta_rowvec_auxfun1", plhs[12]);

engPutVariable(mlEngine_p, "PS_update_delta_rowvec_auxfun2", plhs[13]);

engPutVariable(mlEngine_p, "PS_update_delta_rowvec_combined", plhs[14]);

engEvalString(mlEngine_p, "plot_PS_update_deltas(PS_update_delta_rowvec_auxfun1, PS_update_delta_rowvec_auxfun2, PS_update_delta_rowvec_combined);");

}

void update_PS_update_and_plot_rowvecs(mxArray *plhs[], int n_iter, double delta1, double delta2, double delta3){

PS_update_delta_rowvec_auxfun1(n_iter)=delta1;

PS_update_delta_rowvec_auxfun2(n_iter)=delta2;

PS_update_delta_rowvec_combined(n_iter)=delta3;
	
compute_deltas_PS_send_data_to_Matlab_Eng_and_plot(plhs);

}

static void compute_deltas_X_send_data_to_Matlab_Eng_and_plot(mxArray *plhs[]){

if (!compute_deltas_send_data_to_Matlab_Eng_and_plot_init_flag){

compute_deltas_send_data_to_Matlab_Eng_and_plot_init(plhs);

}

armaSetPr(plhs[12], X_update_delta_rowvec_auxfun1);

armaSetPr(plhs[13], X_update_delta_rowvec_auxfun2);

armaSetPr(plhs[14], X_update_delta_rowvec_combined);

engPutVariable(mlEngine_p, "X_update_delta_rowvec_auxfun1", plhs[12]);

engPutVariable(mlEngine_p, "X_update_delta_rowvec_auxfun2", plhs[13]);

engPutVariable(mlEngine_p, "X_update_delta_rowvec_combined", plhs[14]);

engEvalString(mlEngine_p, "plot_X_update_deltas(X_update_delta_rowvec_auxfun1, X_update_delta_rowvec_auxfun2, X_update_delta_rowvec_combined);");

}

void update_X_update_and_plot_rowvecs(mxArray *plhs[], int n_iter, double delta1, double delta2, double delta3){

X_update_delta_rowvec_auxfun1(n_iter)=delta1;

X_update_delta_rowvec_auxfun2(n_iter)=delta2;

X_update_delta_rowvec_combined(n_iter)=delta3;

compute_deltas_X_send_data_to_Matlab_Eng_and_plot(plhs);
	
}

void compute_deltas_NMF_auxfun_rowvecs_send_data_to_Matlab_Eng_and_plot(mxArray *plhs[]){

if (!compute_deltas_send_data_to_Matlab_Eng_and_plot_init_flag){

compute_deltas_send_data_to_Matlab_Eng_and_plot_init(plhs);

}

armaSetPr(plhs[12], NMF_auxfun1_rowvec);

armaSetPr(plhs[13], NMF_auxfun2_rowvec);

armaSetPr(plhs[14], NMF_rowvec_combined);

engPutVariable(mlEngine_p, "NMF_auxfun1_rowvec", plhs[12]);

engPutVariable(mlEngine_p, "NMF_auxfun2_rowvec", plhs[13]);

engPutVariable(mlEngine_p, "NMF_rowvec_combined", plhs[14]);

engEvalString(mlEngine_p, "plot_NMF_auxfun_rowvecs(NMF_auxfun1_rowvec, NMF_auxfun2_rowvec, NMF_rowvec_combined);");

}