#include "local_inc.hpp"

colvec::fixed<F_static> Fx1_dummy_vec;
cx_colvec::fixed<F_static> Fx1_dummy_cx_vec;
cube::fixed<F_static, K_static, N_static> Phi_S_fkn_local;

static void compute_Phi_S_fkn(cx_cube* expj_Phi_S_fkn_p){

Phi_S_fkn_local=2*atan((imag(*expj_Phi_S_fkn_p))/(sqrt(square(real(*expj_Phi_S_fkn_p))+square(imag(*expj_Phi_S_fkn_p)))+real(*expj_Phi_S_fkn_p)));

}

#define MU_LOCAL 0.15
#define SIGMA_LOCAL 0.04

void Vnk_update_Phi_S(mat* V_nk_p, cx_cube* expj_Phi_S_fkn_p, arg_struct_t* argStruct_p){

compute_Phi_S_fkn(expj_Phi_S_fkn_p);

int n_index, k_index; 

double V_nk_mask_scalar;

for (n_index=0; n_index<N_static; n_index++){

	for (k_index=0; k_index<K_static; k_index++){

	V_nk_mask_scalar=((1/2)*(1+erf(((*V_nk_p)(n_index, k_index)-((double)MU_LOCAL))/(((double)SIGMA_LOCAL)*sqrt(2)))));

	Fx1_dummy_vec=Phi_S_fkn_local.slice(n_index).col(k_index);

	Fx1_dummy_vec=V_nk_mask_scalar*Fx1_dummy_vec;

	Fx1_dummy_cx_vec.set_imag(sin(Fx1_dummy_vec));

	Fx1_dummy_cx_vec.set_real(cos(Fx1_dummy_vec));

	(*expj_Phi_S_fkn_p).slice(n_index).col(k_index)=Fx1_dummy_cx_vec;

	

	}

}

populate_expj_Phi_S_nkf(argStruct_p);

}