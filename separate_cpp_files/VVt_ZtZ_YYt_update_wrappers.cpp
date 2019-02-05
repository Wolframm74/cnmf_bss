#include "local_inc.hpp"

/*#define TURN_ON_LOCAL_FUNCTIONS*/

void YYt_update_wrapper(mat* Y_lk_p, mxArray *Y_lk_mxArray_p, mxArray *plhs[]){

#ifndef TURN_ON_LOCAL_FUNCTIONS

#ifdef ENABLE_MATLAB_ENGINE_PLOTTING
armaSetPr(plhs[3], (*Y_lk_p));
engPutVariable(mlEngine_p, "Y_lk", plhs[3]);
engEvalString(mlEngine_p, "[Y_lk]=YYt_update(Y_lk);");
plhs[3]=engGetVariable(mlEngine_p, "Y_lk");
(*Y_lk_p).set_real(armaGetPr(plhs[3], true)); 	
// Clip the negative elements again..
(*Y_lk_p).elem(find((*Y_lk_p)<=0)).fill(0.00000001);
#endif

#endif

#ifdef TURN_ON_LOCAL_FUNCTIONS

mxArray* output_arg_p[1];

armaSetPr(Y_lk_mxArray_p, *Y_lk_p);

mexCallMATLAB(1, &output_arg_p[0], 1, &Y_lk_mxArray_p, "YYt_update");

(*Y_lk_p).set_real(armaGetPr(output_arg_p[0], true)); 	

#endif

}

void ZtZ_update_wrapper(mat* Z_ol_p, mxArray *Z_ol_mxArray_p, mxArray *plhs[]){

#ifndef TURN_ON_LOCAL_FUNCTIONS

#ifdef ENABLE_MATLAB_ENGINE_PLOTTING
armaSetPr(plhs[2], (*Z_ol_p));
engPutVariable(mlEngine_p, "Z_ol", plhs[2]);
engEvalString(mlEngine_p, "[Z_ol]=ZtZ_update(Z_ol);");
plhs[2]=engGetVariable(mlEngine_p, "Z_ol");
(*Z_ol_p).set_real(armaGetPr(plhs[2], true)); 	
// Clip the negative elements again..
(*Z_ol_p).elem(find((*Z_ol_p)<=0)).fill(0.00000001);
#endif

#endif

#ifdef TURN_ON_LOCAL_FUNCTIONS

mxArray* output_arg_p[1];

armaSetPr(Z_ol_mxArray_p, *Z_ol_p);

mexCallMATLAB(1, &output_arg_p[0], 1, &Z_ol_mxArray_p, "ZtZ_update");

(*Z_ol_p).set_real(armaGetPr(output_arg_p[0], true)); 	

#endif

}

void VVt_update_wrapper(mat* V_nk_p, mxArray *V_nk_mxArray_p, mxArray *plhs[]){

#ifndef TURN_ON_LOCAL_FUNCTIONS

#ifdef ENABLE_MATLAB_ENGINE_PLOTTING
armaSetPr(plhs[5], (*V_nk_p));
engPutVariable(mlEngine_p, "V_nk", plhs[5]);
engEvalString(mlEngine_p, "[V_nk]=VVt_update(V_nk);");
plhs[5]=engGetVariable(mlEngine_p, "V_nk");
(*V_nk_p).set_real(armaGetPr(plhs[5], true)); 	
// Clip the negative elements again..
(*V_nk_p).elem(find((*V_nk_p)<=0)).fill(0.00000001);
#endif

#endif

#ifdef TURN_ON_LOCAL_FUNCTIONS

mxArray* output_arg_p[1];

armaSetPr(V_nk_mxArray_p, *V_nk_p);

mexCallMATLAB(1, &output_arg_p[0], 1, &V_nk_mxArray_p, "VVt_update");

(*V_nk_p).set_real(armaGetPr(output_arg_p[0], true)); 	

#endif

}

void TtT_update_wrapper(mat* T_fk_p, mxArray *T_fk_mxArray_p, mxArray *plhs[]){

#ifndef TURN_ON_LOCAL_FUNCTIONS

#ifdef ENABLE_MATLAB_ENGINE_PLOTTING
armaSetPr(plhs[4], (*T_fk_p));
engPutVariable(mlEngine_p, "T_fk", plhs[4]);
engEvalString(mlEngine_p, "[T_fk]=TtT_update(T_fk);");
plhs[4]=engGetVariable(mlEngine_p, "T_fk");
(*T_fk_p).set_real(armaGetPr(plhs[4], true)); 	
// Clip the negative elements again..
(*T_fk_p).elem(find((*T_fk_p)<=0)).fill(0.00000001);
#endif

#endif

#ifdef TURN_ON_LOCAL_FUNCTIONS

mxArray* output_arg_p[1];

armaSetPr(T_fk_mxArray_p, *T_fk_p);

mexCallMATLAB(1, &output_arg_p[0], 1, &T_fk_mxArray_p, "TtT_update");

(*T_fk_p).set_real(armaGetPr(output_arg_p[0], true)); 	

#endif

}