#include "local_inc.hpp"

static cx_cube::fixed<F_static, N_static, K_static> expj_Phi_S_fnk_local;
static cx_mat::fixed<N_static*K_static, F_static> entry_point1_fun1_mat1_NKxF;
static cx_mat::fixed<F_static, N_static*K_static> entry_point1_fun1_mat2_FxNK;

static void entry_point1_fun1_rotate_Phi_S(cx_cube* expj_Phi_S_nkf_p){

/*copy the nkf input into NKxF */
arrayops::copy(entry_point1_fun1_mat1_NKxF.memptr(), (*expj_Phi_S_nkf_p).memptr(), F_static*N_static*K_static);

/*transpose NKxF into FxNK*/	
entry_point1_fun1_mat2_FxNK=strans(entry_point1_fun1_mat1_NKxF);

/*copy FxNK into output FNK*/	
arrayops::copy(expj_Phi_S_fnk_local.memptr(), entry_point1_fun1_mat2_FxNK.memptr(), F_static*N_static*K_static);

}

static cx_cube::fixed<F_static, N_static, L_static> Y_est_FxNxL;
static cx_mat::fixed<F_static, N_static> entry_point1_fun2_mat1_FN;

static void entry_point1_fun2_compute_Y_est(mat* Y_lk_p, mat* T_fk_p, mat* V_nk_p){

int l_index, k_index;

for (l_index=0; l_index<L_static; l_index++){

	entry_point1_fun2_mat1_FN.zeros();

	for (k_index=0; k_index<K_static; k_index++){

		entry_point1_fun2_mat1_FN=entry_point1_fun2_mat1_FN+((*Y_lk_p)(l_index, k_index))*expj_Phi_S_fnk_local.slice(k_index)%kron((*T_fk_p).col(k_index), trans((*V_nk_p).col(k_index)));

	}


Y_est_FxNxL.slice(l_index)=entry_point1_fun2_mat1_FN;

}

}

static cx_cube::fixed<M_static, F_static, N_static> Xtilde_mfn_local;
static cx_mat::fixed<F_static*N_static, M_static> entry_point1_fun3_fun0_mat1_FNxM;
static cx_mat::fixed<M_static, F_static*N_static> entry_point1_fun3_fun0_mat2_MxFN;

static cx_cube::fixed<L_static, F_static, N_static> Y_est_LxFxN_local;
static cx_mat::fixed<F_static*N_static, L_static> entry_point1_fun3_fun0_mat3_FNxL;
static cx_mat::fixed<L_static, F_static*N_static> entry_point1_fun3_fun0_mat4_LxFN;


static void entry_point1_fun3_fun0_rotate_X_Y(cx_cube* Xtilde_fnm_p){

/*copy the nkf input into NKxF */
arrayops::copy(entry_point1_fun3_fun0_mat1_FNxM.memptr(), (*Xtilde_fnm_p).memptr(), F_static*N_static*M_static);

/*transpose NKxF into FxNK*/	
entry_point1_fun3_fun0_mat2_MxFN=strans(entry_point1_fun3_fun0_mat1_FNxM);

/*copy FxNK into output FNK*/	
arrayops::copy(Xtilde_mfn_local.memptr(), entry_point1_fun3_fun0_mat2_MxFN.memptr(), F_static*N_static*M_static);

/*copy the nkf input into NKxF */
arrayops::copy(entry_point1_fun3_fun0_mat3_FNxL.memptr(), Y_est_FxNxL.memptr(), F_static*N_static*L_static);

/*transpose NKxF into FxNK*/	
entry_point1_fun3_fun0_mat4_LxFN=strans(entry_point1_fun3_fun0_mat3_FNxL);

/*copy FxNK into output FNK*/	
arrayops::copy(Y_est_LxFxN_local.memptr(), entry_point1_fun3_fun0_mat4_LxFN.memptr(), F_static*N_static*L_static);


}

static cx_cube::fixed<M_static, L_static, F_static_admissible> A_pseudoinv_MLF; 
static cx_mat::fixed<M_static, L_static> entry_point1_fun3_fun1_mat1_ML;
static cx_mat::fixed<L_static, L_static> entry_point1_fun3_fun1_mat2_LL;

static void entry_point1_fun3_fun1_compute_A(void){

int f_index, n_index; 

int f_index_admissible;

for (f_index=0; f_index<F_static_admissible; f_index++){

	f_index_admissible=f_index+1; 

	entry_point1_fun3_fun1_mat1_ML.zeros();
	entry_point1_fun3_fun1_mat2_LL.zeros();

	for (n_index=0; n_index<N_static; n_index++){

	entry_point1_fun3_fun1_mat1_ML=entry_point1_fun3_fun1_mat1_ML+Xtilde_mfn_local.slice(n_index).col(f_index_admissible)*trans(Y_est_LxFxN_local.slice(n_index).col(f_index_admissible));

	entry_point1_fun3_fun1_mat2_LL=entry_point1_fun3_fun1_mat2_LL+Y_est_LxFxN_local.slice(n_index).col(f_index_admissible)*trans(Y_est_LxFxN_local.slice(n_index).col(f_index_admissible));

	}

	/*A_pseudoinv_MLF.slice(f_index)=(1/((double)N_static))*entry_point1_fun3_fun1_mat1_ML*pinv(entry_point1_fun3_fun1_mat2_LL);*/
	A_pseudoinv_MLF.slice(f_index)=(1/((double)N_static))*entry_point1_fun3_fun1_mat1_ML*pinv((1/((double)N_static))*entry_point1_fun3_fun1_mat2_LL);

}

}

static cube::fixed<M_static, L_static, F_static_admissible> rhat_MLF; 
mat::fixed<M_static, L_static> rhat_ML; 										/*expose this*/
static cx_cube::fixed<M_static, L_static, F_static_admissible> bhat_MLF; 
static cube::fixed<M_static, L_static, F_static_admissible> bhat_MLF_real; 
static cube::fixed<M_static, L_static, F_static_admissible> bhat_MLF_imag; 

static void entry_point1_fun3_fun1_compute_rhat_bhat(void){

int f_index, l_index, m_index;

double freq_value;

/*populate rhat*/
for (l_index=0; l_index<L_static; l_index++){

	rhat_ML.col(l_index).zeros();

	for (f_index=0; f_index<F_static_admissible; f_index++){

		for (m_index=0; m_index<M_static; m_index++){

			freq_value=(((double)f_index+1)*((double)FS_CONSTANT))/((double)N_STFT_CONSTANT);

			rhat_MLF(m_index, l_index, f_index)=-compute_angle_tdoa_update(A_pseudoinv_MLF(m_index, l_index, f_index)/A_pseudoinv_MLF((int)J_VALUE, l_index, f_index))/(2*(datum::pi)*freq_value);

		}

	rhat_ML.col(l_index)=rhat_ML.col(l_index)+rhat_MLF.slice(f_index).col(l_index);

	}

	rhat_ML.col(l_index)=(1/((double)F_static_admissible))*rhat_ML.col(l_index);

}

/*populate bhat*/
for (l_index=0; l_index<L_static; l_index++){

	for (f_index=0; f_index<F_static_admissible; f_index++){

		freq_value=(((double)f_index+1)*((double)FS_CONSTANT))/((double)N_STFT_CONSTANT);

		/*bhat_MLF.slice(f_index).col(l_index)=exp(-2*(datum::pi)*freq_value*rhat_MLF.slice(f_index).col(l_index));*/

		/*bhat_MLF_real.slice(f_index).col(l_index)=real(cos(-2*(datum::pi)*freq_value*rhat_MLF.slice(f_index).col(l_index)));

		bhat_MLF_imag.slice(f_index).col(l_index)=imag(sin(-2*(datum::pi)*freq_value*rhat_MLF.slice(f_index).col(l_index)));*/

		bhat_MLF_real.slice(f_index).col(l_index)=real(cos(-2*(datum::pi)*freq_value*rhat_ML.col(l_index)));

		bhat_MLF_imag.slice(f_index).col(l_index)=imag(sin(-2*(datum::pi)*freq_value*rhat_ML.col(l_index)));

	}

}

bhat_MLF.set_real(bhat_MLF_real);
bhat_MLF.set_imag(bhat_MLF_imag);

}

static void entry_point1_fun3_compute_A_rhat_bhat(cx_cube* Xtilde_fnm_p){

entry_point1_fun3_fun0_rotate_X_Y(Xtilde_fnm_p);

entry_point1_fun3_fun1_compute_A();

entry_point1_fun3_fun1_compute_rhat_bhat();

}

static cx_cube::fixed<M_static, L_static, F_static_admissible> h_mlf_local;

static void entry_point1_fun4_compute_h_flm(mat* Z_ol_p, cube* W_fom_p, cx_cube* expj_Phi_W_fom_p){

int f_index, l_index, m_index, o_index, f_index_admissible;
cx_double dummy_value_local=0; 

for (f_index=0; f_index<F_static_admissible; f_index++){

	f_index_admissible=f_index+1;

	for (l_index=0; l_index<L_static; l_index++){

		for (m_index=0; m_index<M_static; m_index++){

			dummy_value_local=0;

			for (o_index=0; o_index<O_static; o_index++){

				dummy_value_local=dummy_value_local+((*W_fom_p)(f_index_admissible, o_index, m_index))*((*expj_Phi_W_fom_p)(f_index_admissible, o_index, m_index))*((*Z_ol_p)(o_index, l_index));

			}

			h_mlf_local(m_index, l_index, f_index)=dummy_value_local;

		}

	}

}

}

static mat::fixed<F_static, O_static> partial_d_wrt_Wfom_1;
static mat::fixed<F_static, O_static> partial_d_wrt_Wfom_2;
static cx_colvec::fixed<L_static> entry_point1_fun5_colvec1_Lx1;
static cx_colvec::fixed<L_static> entry_point1_fun5_colvec2_Lx1;

/*static mat::fixed<F_static_admissible, O_static> epsilon_mat_FO;*/

static void entry_point1_fun5_update_W_fom(mat* Z_ol_p, cube* W_fom_p, cx_cube* expj_Phi_W_fom_p){

int o_index, f_index, f_index_admissible; 

int m_index=1;

/*double gamma_local=50;*/

/*double gamma_local=10;*/
double gamma_local=1;

double epsilon_value; 

partial_d_wrt_Wfom_1.zeros();
partial_d_wrt_Wfom_2.zeros();

for (o_index=0; o_index<O_static; o_index++){

	for (f_index=0; f_index<F_static_admissible; f_index++){	

	f_index_admissible=f_index+1;

	entry_point1_fun5_colvec1_Lx1=conj(strans(bhat_MLF.slice(f_index).row(m_index)));

	entry_point1_fun5_colvec1_Lx1=entry_point1_fun5_colvec1_Lx1/(strans(h_mlf_local.slice(f_index).row(m_index)));

	entry_point1_fun5_colvec1_Lx1=((*expj_Phi_W_fom_p)(f_index_admissible, o_index, m_index))*entry_point1_fun5_colvec1_Lx1;

	entry_point1_fun5_colvec1_Lx1=entry_point1_fun5_colvec1_Lx1%trans(((*Z_ol_p).row(o_index)));

	partial_d_wrt_Wfom_1(f_index_admissible, o_index)=-2*as_scalar(sum(real(entry_point1_fun5_colvec1_Lx1)));

	entry_point1_fun5_colvec2_Lx1=strans(h_mlf_local.slice(f_index).row(m_index));

	entry_point1_fun5_colvec2_Lx1=entry_point1_fun5_colvec2_Lx1/strans(h_mlf_local.slice(f_index).row((int)J_VALUE));

	entry_point1_fun5_colvec2_Lx1=2*entry_point1_fun5_colvec2_Lx1/abs(entry_point1_fun5_colvec2_Lx1);

	entry_point1_fun5_colvec2_Lx1=entry_point1_fun5_colvec2_Lx1/strans(h_mlf_local.slice(f_index).row((int)J_VALUE));

	entry_point1_fun5_colvec2_Lx1=((*expj_Phi_W_fom_p)(f_index_admissible, o_index, m_index))*entry_point1_fun5_colvec2_Lx1;

	entry_point1_fun5_colvec2_Lx1=entry_point1_fun5_colvec2_Lx1%trans((*Z_ol_p).row(o_index));

	partial_d_wrt_Wfom_2(f_index_admissible, o_index)=as_scalar(sum(real(entry_point1_fun5_colvec2_Lx1)));

	epsilon_value=(1/gamma_local)*(*W_fom_p)(f_index_admissible, o_index, m_index)/partial_d_wrt_Wfom_2(f_index_admissible, o_index);

	(*W_fom_p)(f_index_admissible, o_index, m_index)=(*W_fom_p)(f_index_admissible, o_index, m_index)-epsilon_value*(partial_d_wrt_Wfom_1(f_index_admissible, o_index)+partial_d_wrt_Wfom_2(f_index_admissible, o_index));

	}

}

}

static mat::fixed<O_static, L_static> partial_d_wrt_Zol_1;
static mat::fixed<O_static, L_static> partial_d_wrt_Zol_2;
static cx_colvec::fixed<F_static_admissible> entry_point1_fun5_colvec1_Fx1;
static cx_colvec::fixed<F_static_admissible> entry_point1_fun5_colvec2_Fx1;
static cx_colvec::fixed<F_static_admissible> entry_point1_fun5_colvec3_Fx1;

static void entry_point1_fun6_update_Z_ol(mat* Z_ol_p, cube* W_fom_p, cx_cube* expj_Phi_W_fom_p){

int o_index, l_index, f_index, f_index_admissible;

int m_index=1; 

/*double gamma_local=50;*/

double gamma_local=10;

double epsilon_value; 

for (o_index=0; o_index<O_static; o_index++){

	for (l_index=0; l_index<L_static; l_index++){

		for (f_index=0; f_index<F_static_admissible; f_index++){

			f_index_admissible=f_index+1; 

			entry_point1_fun5_colvec1_Fx1(f_index)=conj(bhat_MLF(m_index, l_index, f_index));

			entry_point1_fun5_colvec3_Fx1(f_index)=h_mlf_local((int)J_VALUE, l_index, f_index)*((*W_fom_p)(f_index_admissible, o_index, m_index))*((*expj_Phi_W_fom_p)(f_index_admissible, o_index, m_index))-h_mlf_local(m_index, l_index, f_index)*((*W_fom_p)(f_index_admissible, o_index, (int)J_VALUE))*((*expj_Phi_W_fom_p)(f_index_admissible, o_index, (int)J_VALUE));

			entry_point1_fun5_colvec3_Fx1(f_index)=entry_point1_fun5_colvec3_Fx1(f_index)/(pow(abs(h_mlf_local(m_index, l_index, f_index)), 2));

			entry_point1_fun5_colvec1_Fx1(f_index)=entry_point1_fun5_colvec1_Fx1(f_index)*entry_point1_fun5_colvec3_Fx1(f_index);

			entry_point1_fun5_colvec2_Fx1(f_index)=h_mlf_local(m_index, l_index, f_index)/h_mlf_local((int)J_VALUE, l_index, f_index);

			entry_point1_fun5_colvec2_Fx1(f_index)=entry_point1_fun5_colvec2_Fx1(f_index)/abs(entry_point1_fun5_colvec2_Fx1(f_index));

			entry_point1_fun5_colvec2_Fx1(f_index)=entry_point1_fun5_colvec2_Fx1(f_index)*entry_point1_fun5_colvec3_Fx1(f_index);

		}

		partial_d_wrt_Zol_1(o_index, l_index)=-2*as_scalar(sum(real(entry_point1_fun5_colvec1_Fx1)));

		partial_d_wrt_Zol_2(o_index, l_index)=as_scalar(sum(real(entry_point1_fun5_colvec2_Fx1)));

		epsilon_value=(1/gamma_local)*(*Z_ol_p)(o_index, l_index)/partial_d_wrt_Zol_2(o_index, l_index);

		(*Z_ol_p)(o_index, l_index)=(*Z_ol_p)(o_index, l_index)-epsilon_value*(partial_d_wrt_Zol_1(o_index, l_index)+partial_d_wrt_Zol_2(o_index, l_index));

	}

}

}

void TDOA_update_module1_entry_point1(cx_cube* Xtilde_fnm_p, mat* Z_ol_p, mat* Y_lk_p, mat* T_fk_p, mat* V_nk_p, cube* W_fom_p, cx_cube* expj_Phi_W_fom_p, cx_cube* expj_Phi_S_nkf_p){

entry_point1_fun1_rotate_Phi_S(expj_Phi_S_nkf_p);

entry_point1_fun2_compute_Y_est(Y_lk_p, T_fk_p, V_nk_p);

entry_point1_fun3_compute_A_rhat_bhat(Xtilde_fnm_p);

entry_point1_fun4_compute_h_flm(Z_ol_p, W_fom_p, expj_Phi_W_fom_p);

/*entry_point1_fun5_update_W_fom(Z_ol_p, W_fom_p, expj_Phi_W_fom_p);*/

entry_point1_fun6_update_Z_ol(Z_ol_p, W_fom_p, expj_Phi_W_fom_p);

rhat_ML.print("rhat_ML: line 340:");

}

