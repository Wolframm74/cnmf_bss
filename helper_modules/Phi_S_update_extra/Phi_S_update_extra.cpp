#include "local_inc.hpp"

/*static void T_fk_update_Phi_S_wrapper(){



}*/

cx_mat::fixed<K_static, F_static> dummy_mat_cx_KxF;
cx_rowvec::fixed<K_static*F_static> dummy_row_cx_KF;
cx_mat::fixed<N_static, K_static*F_static> dummy_mat_cx_NxKF;

void populate_expj_Phi_S_nkf(arg_struct_t* argStruct_p){

int n_index;

for (n_index=0; n_index<N_static; n_index++){

dummy_mat_cx_KxF=strans((*(argStruct_p->expj_Phi_S_fkn_p)).slice(n_index));

/*Copy the matrix into the row*/
/*memcpy*/
arrayops::copy(dummy_row_cx_KF.memptr(), dummy_mat_cx_KxF.memptr(), K_static*F_static);

dummy_mat_cx_NxKF.row(n_index)=dummy_row_cx_KF;

}

arrayops::copy((*(argStruct_p->expj_Phi_S_nkf_p)).memptr(), dummy_mat_cx_NxKF.memptr(), N_static*K_static*F_static);

}

void V_nk_update_Phi_S_wrapper(mxArray *plhs[], arg_struct_t* argStruct_p){

(scalar_global[0])(0)=(double)F_static;
(scalar_global[1])(0)=(double)K_static; 
(scalar_global[2])(0)=(double)N_static; 

armaSetPr(scalar_global_p[0], scalar_global[0]);
armaSetPr(scalar_global_p[1], scalar_global[1]);
armaSetPr(scalar_global_p[2], scalar_global[2]);

engPutVariable(mlEngine_p, "F_static", scalar_global_p[0]);
engPutVariable(mlEngine_p, "K_static", scalar_global_p[1]);
engPutVariable(mlEngine_p, "N_static", scalar_global_p[2]);

//Set V, put it into the Matlab workspace 
armaSetPr(plhs[5], (*(argStruct_p->V_nk_p)));
engPutVariable(mlEngine_p, "V_nk", plhs[5]);

//Set expj_Phi_S, put it into the Matlab workspace
armaSetCubeCx(plhs[7], (*(argStruct_p->expj_Phi_S_fkn_p)));
engPutVariable(mlEngine_p, "expj_Phi_S_fkn", plhs[7]);

//Call the function
engEvalString(mlEngine_p, "[expj_Phi_S_fkn_out]=Vnk_update_Phi_S_mex(V_nk, F_static, K_static, N_static, expj_Phi_S_fkn);");

plhs[7]=engGetVariable(mlEngine_p, "expj_Phi_S_fkn_out");

//Write the result into 
(*(argStruct_p->expj_Phi_S_fkn_p)).set_real(armaGetCubePr(plhs[7],true));
(*(argStruct_p->expj_Phi_S_fkn_p)).set_imag(armaGetCubePi(plhs[7],true));	

populate_expj_Phi_S_nkf(argStruct_p);
	
}


