#include "local_inc.hpp"

/*2nd sum stuff*/
/*cube::fixed<F_static, N_static, K_static> outtensor_2nd_real_FNK;
cx_cube::fixed<F_static, N_static, K_static> outtensor_2nd_cx_FNK;*/

cube::fixed<F_static, K_static, N_static> outtensor_2nd_real_FKN;
cx_cube::fixed<F_static, K_static, N_static> outtensor_2nd_cx_FKN;

cube::fixed<N_static, K_static, F_static> outtensor_2nd_real_NKF;
cx_cube::fixed<N_static, K_static, F_static> outtensor_2nd_cx_NKF;

/*output f,n,k related stuff*/
/*mat::fixed<F_static, N_static> dummy_mat_real_FN;
cx_mat::fixed<F_static, N_static> dummy_mat_cx_FN;
cx_cube::fixed<F_static, N_static, K_static> outtensor_target_cx_FNK;*/	/*Do a mexCallMATLAB() on this */

cx_cube::fixed<F_static, K_static, N_static> tensor_dummy_cx_FKN;	/*let this one be related to the cx tensor*/
cx_cube::fixed<N_static, K_static, F_static> tensor_dummy_cx_NKF; 

cx_cube::fixed<F_static, K_static, N_static> outtensor_target_cx_FKN;	/*let this one be related to the real tensor*/
cx_cube::fixed<N_static, K_static, F_static> outtensor_target_cx_NKF; 

/*Phase_FKN and Phase_NKF's relationship to the outtensor_target's is conj(outtensor_target_cx*___) */

cube::fixed<F_static, M_static, N_static> Xhat_low_fmn_local;
cx_cube::fixed<F_static, M_static, N_static> E_conj_fmn_local;

static void compute_2nd_sum(mat* T_fk_p, mat* V_nk_p, cube* Xhat_low_fnm_p, cx_cube* E_conj_fnm_p, int n_index){

#ifdef DEBUG	
mexPrintf("Phi_S_update: pthread_self(): %d, compute_2nd_sum(), n_index=%d\n", (int)pthread_self(), n_index);		
#endif	

int f_index; 
int m_index; 

for (m_index=0; m_index<M_static; m_index++){

Xhat_low_fmn_local.slice(n_index).col(m_index)=(*Xhat_low_fnm_p).slice(m_index).col(n_index);

E_conj_fmn_local.slice(n_index).col(m_index)=(*E_conj_fnm_p).slice(m_index).col(n_index);

}


for (f_index=0; f_index<F_static; f_index++){

/*Dimensions of matrix product: 1xK=(1xM)*(MxK) */

/*outtensor_2nd_real_FNK.subcube(f_index,n_index,0,f_index,n_index,K_static-1)=(*Xhat_low_fnm_p).subcube(f_index,n_index,0,f_index,n_index,M_static-1)*Xhat_outtensor_real_FMK.subcube(f_index,0,0,f_index,M_static-1,K_static-1);*/
/*outtensor_2nd_real_FKN.slice(n_index).row(f_index)=(*Xhat_low_fmn_p).slice(n_index).row(f_index)*(Xhat_outtensor_real_MKF).slice(f_index);*/
outtensor_2nd_real_FKN.slice(n_index).row(f_index)=Xhat_low_fmn_local.slice(n_index).row(f_index)*(Xhat_outtensor_real_MKF).slice(f_index);	
outtensor_2nd_real_NKF.slice(f_index).row(n_index)=trans(outtensor_2nd_real_FKN.slice(n_index).row(f_index));

/*outtensor_2nd_cx_FNK.subcube(f_index,n_index,0,f_index,n_index,K_static-1)=(*E_conj_fmn_p).subcube(f_index,n_index,0,f_index,n_index,M_static-1)*Xhat_outtensor_cx_FMK.subcube(f_index,0,0,f_index,M_static-1,K_static-1);*/
/*outtensor_2nd_cx_FKN.slice(n_index).row(f_index)=(*E_conj_fmn_p).slice(n_index).row(f_index)*(Xhat_outtensor_cx_MKF).slice(f_index);*/
outtensor_2nd_cx_FKN.slice(n_index).row(f_index)=E_conj_fmn_local.slice(n_index).row(f_index)*(Xhat_outtensor_cx_MKF).slice(f_index);
outtensor_2nd_cx_NKF.slice(f_index).row(n_index)=strans(outtensor_2nd_cx_FKN.slice(n_index).row(f_index));

/*Do a few more steps to compute the dummy and output tensors while you're here*/
tensor_dummy_cx_FKN.slice(n_index).row(f_index)=(outtensor_2nd_cx_FKN.slice(n_index).row(f_index))%((*T_fk_p).row(f_index))%((*V_nk_p).row(n_index));

tensor_dummy_cx_NKF.slice(f_index).row(n_index)=tensor_dummy_cx_FKN.slice(n_index).row(f_index);

outtensor_target_cx_FKN.slice(n_index).row(f_index)=conj(outtensor_target_cx_FKN.slice(n_index).row(f_index))%((*T_fk_p).row(f_index))%((*V_nk_p).row(n_index));	/*Should be noted that this step should be properly related to the Phi_S=-phase(arg) step */

outtensor_target_cx_NKF.slice(f_index).row(n_index)=outtensor_target_cx_FKN.slice(n_index).row(f_index);

}

}

/*static void compute_output_fnk(mat* T_fk_p, mat* V_nk_p, cx_cube* expj_Phi_S_fnk_p, int n_index, int f_index){

//tensor_dummy stuff
tensor_dummy_cx_FKN.slice(n_index).row(f_index)=outtensor_2nd_cx_FKN.slice(n_index).row(f_index)

tensor_dummy_cx_NKF=

//outtensor_target stuff

outtensor_target_cx_FKN=

outtensor_target_cx_NKF=

}*/

void* Phi_S_start(void* arg){

arg_struct_t* argStruct_p;

	
}