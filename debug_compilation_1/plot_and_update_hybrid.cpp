#include "local_inc.hpp"

/*globals local to this module*/
mxArray* input_arg_p;
colvec::fixed<3> some_mat;

#ifdef COMPUTE_L_SWITCH
static double compute_L(arg_struct_t* argStruct_p){

int f, n, m;

double L_value=0;
double absvalue;

	for (f=0; f<(argStruct_p->F); f++){

		for (n=0; n<(argStruct_p->N); n++){

			for (m=0; m<(argStruct_p->M); m++){

				absvalue=abs( (*(argStruct_p->Xtilde_fnm_p))(f,n,m)-(*(argStruct_p->Xhat_fnm_p))(f,n,m));

				L_value=L_value+absvalue*absvalue;
				mexPrintf("f:%d, n:%d, m:%d, L_value:%f ........\n", f, n, m, L_value);

			}


		}


	}

return L_value;

}

#endif 

void create_lhs_matrices(mxArray *plhs[], arg_struct_t* argStruct_p){

plhs[0]=armaCreateMxMatrix((*(argStruct_p->Xhat_fnm_p)).n_rows, (*(argStruct_p->Xhat_fnm_p)).n_cols, (*(argStruct_p->Xhat_fnm_p)).n_slices, mxDOUBLE_CLASS, mxCOMPLEX);

plhs[1]=armaCreateMxMatrix((*(argStruct_p->W_fom_p)).n_rows, (*(argStruct_p->W_fom_p)).n_cols, (*(argStruct_p->W_fom_p)).n_slices, mxDOUBLE_CLASS, mxREAL);

plhs[2]=armaCreateMxMatrix((*(argStruct_p->Z_ol_p)).n_rows, (*(argStruct_p->Z_ol_p)).n_cols, mxDOUBLE_CLASS, mxREAL);

plhs[3]=armaCreateMxMatrix((*(argStruct_p->Y_lk_p)).n_rows, (*(argStruct_p->Y_lk_p)).n_cols, mxDOUBLE_CLASS, mxREAL);

plhs[4]=armaCreateMxMatrix((*(argStruct_p->T_fk_p)).n_rows, (*(argStruct_p->T_fk_p)).n_cols, mxDOUBLE_CLASS, mxREAL);

plhs[5]=armaCreateMxMatrix((*(argStruct_p->V_nk_p)).n_rows, (*(argStruct_p->V_nk_p)).n_cols, mxDOUBLE_CLASS, mxREAL);

plhs[6]=armaCreateMxMatrix((*(argStruct_p->expj_Phi_S_nkf_p)).n_rows, (*(argStruct_p->expj_Phi_S_nkf_p)).n_cols, (*(argStruct_p->expj_Phi_S_nkf_p)).n_slices, mxDOUBLE_CLASS, mxCOMPLEX);

plhs[7]=armaCreateMxMatrix((*(argStruct_p->expj_Phi_S_fkn_p)).n_rows, (*(argStruct_p->expj_Phi_S_fkn_p)).n_cols, (*(argStruct_p->expj_Phi_S_fkn_p)).n_slices, mxDOUBLE_CLASS, mxCOMPLEX);

/*Additionally create the "input_arg" matrix*/
input_arg_p=armaCreateMxMatrix(3,1, mxDOUBLE_CLASS, mxREAL);

/*Init this dummy vector to zeros*/
some_mat.zeros();

}

/*update_engine: case 1
plhs[0]
*/
void plot_and_update_Xhat(mxArray *plhs[], arg_struct_t* argStruct_p, int n){

double L_value;

/*Compute Xhat_fnm, B_flkmn, C_flkmn ---------------------------------------------------------------- */
update_engine_M_threads(argStruct_p);
update_engine(argStruct_p, 1);

#ifdef 	MEX_UNIT_TESTING
/*armaSetCubeCx(array_p[0], Xhat_outtensor_cx_FLM);
armaSetCubePr(array_p[1], Xhat_outtensor_real_FLM);
armaSetCubeCx(array_p[2], Xhat_outtensor_cx_FKM);
armaSetCubePr(array_p[3], Xhat_outtensor_real_FKM);
mexCallMATLAB(0, NULL, 4, array_p, "plot_X_matrices");*/

armaSetCubeCx(plhs[0], (*(argStruct_p->Xhat_fnm_p)));
mexCallMATLAB(0, NULL, 1, &plhs[0], "plot_Xhat_fnm");

	#ifdef COMPUTE_L_SWITCH
	L_value=compute_L(argStruct_p);
	some_mat(2)=L_value;
	#endif

#endif

/* Print status to MATLAB */
some_mat(0)=(double)n; 
some_mat(1)=(double)1;
armaSetPr(input_arg_p, some_mat);
mexCallMATLAB(0, NULL, 1, &input_arg_p, "nguyen_2015_print");

}

/*update_engine: case 2
plhs[4]
*/
void plot_and_update_Tfk(mxArray *plhs[], arg_struct_t* argStruct_p, int n){

mexPrintf("About to compute T_fk \n");

/*Compute T_fk  ------------------------------------------------------------------------------------ */	
update_engine(argStruct_p, 2);

#ifdef 	MEX_UNIT_TESTING
armaSetPr(plhs[4], (*(argStruct_p->T_fk_p)));
mexCallMATLAB(0, NULL, 1, &plhs[4], "plot_T_fk");
#endif

/* Print status to MATLAB */
some_mat(0)=(double)n; 
some_mat(1)=(double)2;
armaSetPr(input_arg_p, some_mat);
mexCallMATLAB(0, NULL, 1, &input_arg_p, "nguyen_2015_print");

}

/*update_engine: case 3
plhs[5]
*/
void plot_and_update_Vnk(mxArray *plhs[], arg_struct_t* argStruct_p, int n){	

/*Call the function to preprocess and output a required matrix*/
V_0th_sum_preprocess(argStruct_p->V_nk_p);

/*Compute V_nk  ------------------------------------------------------------------------------------ */
update_engine(argStruct_p, 3);

#ifdef 	MEX_UNIT_TESTING
armaSetPr(plhs[5], (*(argStruct_p->V_nk_p)));
mexCallMATLAB(0, NULL, 1, &plhs[5], "plot_V_nk");
#endif

/* Print status to MATLAB */
some_mat(0)=(double)n; 
some_mat(1)=(double)3;
armaSetPr(input_arg_p, some_mat);
mexCallMATLAB(0, NULL, 1, &input_arg_p, "nguyen_2015_print");

}

/*update_engine: case 6
plhs[2]
*/
void plot_and_update_Zol(mxArray *plhs[], arg_struct_t* argStruct_p, int n){	

double accum_Z_ol_local; 

mexPrintf("About to compute Z_ol \n");

/*Populate this for Z update*/
Z_inmat_real_Ykl=trans(*(argStruct_p->Y_lk_p));

/*Compute Z_ol  ------------------------------------------------------------------------- */	
update_engine(argStruct_p, 6);

/*Call projfunc_wrapper() */
projfunc_wrapper_Zol(argStruct_p->Z_ol_p, argStruct_p->Z_ol_mxArray_p);

/*Subtract off the current mean*/
accum_Z_ol_local=accu(*(argStruct_p->Z_ol_p));

*(argStruct_p->Z_ol_p)=(*(argStruct_p->Z_ol_p))-(accum_Z_ol_local/(((double)O_static)*((double)L_static)))*ones(O_static, L_static);

/*Set new mean forcibly*/
*(argStruct_p->Z_ol_p)=(*(argStruct_p->Z_ol_p))+(0.95)*(accum_Z_ol/(((double)O_static)*((double)L_static)))*ones(O_static, L_static);

/*Clip the negative elements again..*/
(*(argStruct_p->Z_ol_p)).elem(find((*(argStruct_p->Z_ol_p))<0)).fill(0.00000001);

#ifdef 	MEX_UNIT_TESTING
armaSetPr(plhs[2], (*(argStruct_p->Z_ol_p)));
mexCallMATLAB(0, NULL, 1, &plhs[2], "plot_Z_ol");
#endif

/* Print status to MATLAB */
some_mat(0)=(double)n; 
some_mat(1)=(double)6;
armaSetPr(input_arg_p, some_mat);
mexCallMATLAB(0, NULL, 1, &input_arg_p, "nguyen_2015_print");

}

/*update_engine: case 7
plhs[3]
*/
void plot_and_update_Ylk(mxArray *plhs[], arg_struct_t* argStruct_p, int n){	

mexPrintf("About to compute Y_lk \n");

/*Compute Y_lk  ------------------------------------------------------------------------- */	
update_engine(argStruct_p, 7);

projfunc_wrapper_Ylk(argStruct_p->Y_lk_p, argStruct_p->Y_lk_mxArray_p);

normalize_Y_lk(argStruct_p->Y_lk_p, argStruct_p->V_nk_p);

#ifdef 	MEX_UNIT_TESTING
armaSetPr(plhs[3], (*(argStruct_p->Y_lk_p)));
mexCallMATLAB(0, NULL, 1, &plhs[3], "plot_Y_lk");
#endif

mexPrintf("algorithm_entry_mex.cpp: pthread_self(): %d, RETURNED FROM UPDATE ENGINE FOR Y_LK, threads_asleep_ctr: %d\n", (int)pthread_self(), threads_asleep_ctr);

/* Print status to MATLAB */
some_mat(0)=(double)n; 
some_mat(1)=(double)7;
armaSetPr(input_arg_p, some_mat);
mexCallMATLAB(0, NULL, 1, &input_arg_p, "nguyen_2015_print");

}

/*update_engine: case 4
plhs[1]
*/
void plot_and_update_Wfom(mxArray *plhs[], arg_struct_t* argStruct_p, int n){	

mexPrintf("About to compute W_fom \n");

/*Call the function to preprocess and output a required matrix*/
W_0th_sum_preprocess(argStruct_p->Z_ol_p, argStruct_p->Y_lk_p);

/*Compute W_fom  ------------------------------------------------------------------------- */	
update_engine(argStruct_p, 4);

#ifdef BOTTOM_UP_MEX_FLAG

W_callMATLAB_wrapper();

#endif

populate_W_fom_cx(argStruct_p->W_fom_p, argStruct_p->W_fom_cx_p, argStruct_p->Phi_W_fom_p);

#ifdef 	MEX_UNIT_TESTING
armaSetCubePr(plhs[1], (*(argStruct_p->W_fom_p)));
mexCallMATLAB(0, NULL, 1, &plhs[1], "plot_W_fom");
#endif

/* Print status to MATLAB */
some_mat(0)=(double)n; 
some_mat(1)=(double)4;
armaSetPr(input_arg_p, some_mat);
mexCallMATLAB(0, NULL, 1, &input_arg_p, "nguyen_2015_print");

}

/*Argument list:
array_p[0] expj_Phi_S_fkn, 
array_p[1] T_fk, 
array_p[2] V_nk, 
array_p[3] Xhat_low_fnm, 
array_p[4] E_conj_fnm, 
array_p[5] Xhat_outtensor_real_MKF, 
array_p[6] Xhat_outtensor_cx_MKF, 
array_p[7] Y_lk, 
array_p[8] Z_ol, 
array_p[9] W_fom, 
array_p[10] expj_Phi_W_fom
array_p[11] dim_mat
*/
static void mexCallMATLAB_Phi_S(mxArray *plhs[], arg_struct_t* argStruct_p){

colvec::fixed<6> dim_mat;
mxArray *array_p[12];

/*array_p[0] expj_Phi_S_fkn, */
array_p[0]=plhs[7];
armaSetCubeCx(array_p[0], *(argStruct_p->expj_Phi_S_fkn_p));

/*array_p[1] T_fk, */
array_p[1]=plhs[4];
armaSetPr(array_p[1], *(argStruct_p->T_fk_p));

/*array_p[2] V_nk, */
array_p[2]=plhs[5];
armaSetPr(array_p[2], *(argStruct_p->V_nk_p));

/*array_p[3] Xhat_low_fnm, */
array_p[3]=armaCreateMxMatrix((*(argStruct_p->Xhat_low_fnm_p)).n_rows, (*(argStruct_p->Xhat_low_fnm_p)).n_cols, (*(argStruct_p->Xhat_low_fnm_p)).n_slices, mxDOUBLE_CLASS, mxREAL);
armaSetCubePr(array_p[3], *(argStruct_p->Xhat_low_fnm_p));

/*array_p[4] E_conj_fnm, */
array_p[4]=armaCreateMxMatrix((*(argStruct_p->E_conj_fnm_p)).n_rows, (*(argStruct_p->E_conj_fnm_p)).n_cols, (*(argStruct_p->E_conj_fnm_p)).n_slices, mxDOUBLE_CLASS, mxCOMPLEX);
armaSetCubeCx(array_p[4], *(argStruct_p->E_conj_fnm_p));

/*array_p[5] Xhat_outtensor_real_MKF, */
array_p[5]=armaCreateMxMatrix(Xhat_outtensor_real_MKF.n_rows, Xhat_outtensor_real_MKF.n_cols, Xhat_outtensor_real_MKF.n_slices, mxDOUBLE_CLASS, mxREAL);
armaSetCubePr(array_p[5], Xhat_outtensor_real_MKF);

/*array_p[6] Xhat_outtensor_cx_MKF, */
array_p[6]=armaCreateMxMatrix(Xhat_outtensor_cx_MKF.n_rows, Xhat_outtensor_cx_MKF.n_cols, Xhat_outtensor_cx_MKF.n_slices, mxDOUBLE_CLASS, mxCOMPLEX);
armaSetCubeCx(array_p[6], Xhat_outtensor_cx_MKF);

/*array_p[7] Y_lk, */
array_p[7]=plhs[3];
armaSetPr(array_p[7], *(argStruct_p->Y_lk_p));

/*array_p[8] Z_ol, */
array_p[8]=plhs[2];
armaSetPr(array_p[8], *(argStruct_p->Z_ol_p));

/*array_p[9] W_fom, */
array_p[9]=plhs[1];
armaSetCubePr(array_p[9], *(argStruct_p->W_fom_p));

/*array_p[10] expj_Phi_W_fom*/
array_p[10]=armaCreateMxMatrix((*(argStruct_p->expj_Phi_W_fom_p)) .n_rows, (*(argStruct_p->expj_Phi_W_fom_p)) .n_cols, (*(argStruct_p->expj_Phi_W_fom_p)) .n_slices, mxDOUBLE_CLASS, mxCOMPLEX);
armaSetCubeCx(array_p[10], *(argStruct_p->expj_Phi_W_fom_p));

/*array_p[11] dim_mat*/
array_p[11]=armaCreateMxMatrix(6,1, mxDOUBLE_CLASS, mxREAL);

dim_mat(0)=M_static;
dim_mat(1)=F_static;
dim_mat(2)=N_static;
dim_mat(3)=K_static;
dim_mat(4)=L_static;
dim_mat(5)=O_static;

armaSetPr(array_p[11], dim_mat);

mexCallMATLAB(0, NULL, 12, &array_p[0], "Phi_S_unit_test");

}

/*update_engine: case 5
plhs[6], plhs[7]
*/
void plot_and_update_Phi_S(mxArray *plhs[], arg_struct_t* argStruct_p, int n){

/*mexCallMATLAB_Phi_S(plhs, argStruct_p);*/

/*Compute Phi_S_fnk  ------------------------------------------------------------------------------- */	
update_engine(argStruct_p, 5); 

/*Plot Phi_S_fnk??*/
#ifdef 	MEX_UNIT_TESTING
armaSetCubeCx(plhs[7], (*(argStruct_p->expj_Phi_S_fkn_p)));
mexCallMATLAB(0, NULL, 1, &plhs[7], "plot_Phi_S_fnk");
#endif

/* Print status to MATLAB */
some_mat(0)=(double)n; 
some_mat(1)=(double)5;
armaSetPr(input_arg_p, some_mat);
mexCallMATLAB(0, NULL, 1, &input_arg_p, "nguyen_2015_print");

}	