#include "local_inc.hpp"

/*globals local to this module*/
mxArray* input_arg_p;
colvec::fixed<3> some_mat;

mxArray* scalar_global_p[3];
colvec::fixed<1> scalar_global[3];

double L_value_main_objective;
double L_value_threshold_local;
double L_value_threshold_local2;

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

/*L_value_threshold_local2=2.25*pow(10, 4);*/
L_value_threshold_local2=2.5*pow(10, 4);	
L_value_threshold_local=2.5*pow(10, 4);
L_value_main_objective=L_value;

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

plhs[8]=armaCreateMxMatrix((*(argStruct_p->W_fom_cx_p)).n_rows, (*(argStruct_p->W_fom_cx_p)).n_cols, (*(argStruct_p->W_fom_cx_p)).n_slices, mxDOUBLE_CLASS, mxCOMPLEX);

/*Additionally create the "input_arg" matrix*/
input_arg_p=armaCreateMxMatrix(3,1, mxDOUBLE_CLASS, mxREAL);

scalar_global_p[0]=armaCreateMxMatrix(1,1, mxDOUBLE_CLASS, mxREAL);
scalar_global_p[1]=armaCreateMxMatrix(1,1, mxDOUBLE_CLASS, mxREAL);
scalar_global_p[2]=armaCreateMxMatrix(1,1, mxDOUBLE_CLASS, mxREAL);

/*Init this dummy vector to zeros*/
some_mat.zeros();

}

static void plot_Zok_wrapper(mxArray* plhs[]){

engPutVariable(mlEngine_p, "Z_ol", plhs[2]);
engPutVariable(mlEngine_p, "Y_lk", plhs[3]);
engEvalString(mlEngine_p, "plot_Z_ok(Z_ol, Y_lk);");

}

/*To do list: 

*/
static void post_update_Tfk(mxArray *plhs[], arg_struct_t* argStruct_p, int n){

#ifdef ENABLE_MATLAB_ENGINE_PLOTTING
/*armaSetPr(plhs[4], (*(argStruct_p->T_fk_p)));
engPutVariable(mlEngine_p, "T_fk", plhs[4]);
engEvalString(mlEngine_p, "plot_T_fk(T_fk,1);");

engEvalString(mlEngine_p, "[T_fk]=projfunc_wrapper_matlab_Tfk(T_fk);");
*/
#endif

/*Call projfunc_wrapper() */
/*projfunc_wrapper_Tfk(argStruct_p->T_fk_p, argStruct_p->T_fk_mxArray_p);*/

#ifdef ENABLE_MATLAB_ENGINE_PLOTTING
/*armaSetPr(plhs[4], (*(argStruct_p->T_fk_p)));
engPutVariable(mlEngine_p, "T_fk", plhs[4]);
engEvalString(mlEngine_p, "plot_T_fk(T_fk,1);");*/
#endif

// Subtract off the current mean
/*accum_T_fk_local=accu(*(argStruct_p->T_fk_p));

*(argStruct_p->T_fk_p)=(*(argStruct_p->T_fk_p))-(accum_T_fk_local/(((double)F_static)*((double)K_static)))*ones(F_static, K_static);

// Set new mean forcibly
*(argStruct_p->T_fk_p)=(*(argStruct_p->T_fk_p))+(1)*(accum_T_fk/(((double)F_static)*((double)K_static)))*ones(F_static, K_static);*/

/*T_pull_up_vector_elements_generalized_gaussian_wrapper(argStruct_p->T_fk_p);*/

(*(argStruct_p->T_fk_p)).elem( find_nonfinite((*(argStruct_p->T_fk_p))) ).zeros();

(*(argStruct_p->T_fk_p)).elem(find((*(argStruct_p->T_fk_p))<=0)).fill(0.00000001);

(*(argStruct_p->T_fk_p)).elem(find((*(argStruct_p->T_fk_p))>100)).fill(100);

/*if (!N_ITERATIONS_TV_bool){

TtT_update_wrapper(argStruct_p->T_fk_p, argStruct_p->T_fk_mxArray_p, plhs);

}*/

}	

/*To do list: 

*/
static void post_update_Vnk(mxArray *plhs[], arg_struct_t* argStruct_p, int n){


/*orthogonalize_Vnk_wrapper(argStruct_p->V_nk_p, argStruct_p->V_nk_mxArray_p);

(*(argStruct_p->V_nk_p)).elem(find((*(argStruct_p->V_nk_p))<=0)).fill(0.00000001);*/

VVt_update_wrapper(argStruct_p->V_nk_p, argStruct_p->V_nk_mxArray_p, plhs);

if (L_value_main_objective>(L_value_threshold_local)){

V_nk_MA_filter_rows(argStruct_p->V_nk_p);

V_weighted_scaling(argStruct_p->V_nk_p);

// Clip the negative elements again..
(*(argStruct_p->V_nk_p)).elem(find((*(argStruct_p->V_nk_p))<=0)).fill(0.00000001);

(*(argStruct_p->V_nk_p)).elem(find((*(argStruct_p->V_nk_p))>10)).fill(10);

}

#ifdef ENABLE_MATLAB_ENGINE_PLOTTING

(scalar_global[0])(0)=(double)N_static;
(scalar_global[1])(0)=(double)K_static; 

armaSetPr(scalar_global_p[0], scalar_global[0]);
armaSetPr(scalar_global_p[1], scalar_global[1]);

armaSetPr(plhs[5], (*(argStruct_p->V_nk_p)));

/*engPutVariable(mlEngine_p, "N_static", scalar_global_p[0]);
engPutVariable(mlEngine_p, "K_static", scalar_global_p[1]);
engPutVariable(mlEngine_p, "V_nk", plhs[5]);

engEvalString(mlEngine_p, "[V_nk]=projfunc_wrapper_matlab_Vnk_mex(V_nk, N_static, K_static);");
plhs[5]=engGetVariable(mlEngine_p, "V_nk");
(*(argStruct_p->V_nk_p)).set_real(armaGetPr(plhs[5], true)); 	

// Clip the negative elements again..
(*(argStruct_p->V_nk_p)).elem(find((*(argStruct_p->V_nk_p))<=0)).fill(0.00000001);*/
#endif

if (!N_ITERATIONS_TV_bool){

/*V_nk_update_Phi_S_wrapper(plhs, argStruct_p);*/
/*Vnk_update_Phi_S(argStruct_p->V_nk_p, argStruct_p->expj_Phi_S_fkn_p, argStruct_p);*/

/*V_nk_erf_threshold(argStruct_p->V_nk_p);

// Clip the negative elements again..
(*(argStruct_p->V_nk_p)).elem(find((*(argStruct_p->V_nk_p))<=0)).fill(0.00000001);

(*(argStruct_p->V_nk_p)).elem(find((*(argStruct_p->V_nk_p))>10)).fill(10);*/

/*Call the new function to properly scale the cols of T_fk */
/*check_V_scale_down_T(argStruct_p->T_fk_p, argStruct_p->V_nk_p);

// Clip the negative elements again..
(*(argStruct_p->V_nk_p)).elem(find((*(argStruct_p->V_nk_p))<=0)).fill(0.00000001);

(*(argStruct_p->V_nk_p)).elem(find((*(argStruct_p->V_nk_p))>10)).fill(10);*/

/*normalize_V_nk_frob_norm(argStruct_p->V_nk_p, argStruct_p->T_fk_p);

(*(argStruct_p->V_nk_p)).elem(find((*(argStruct_p->V_nk_p))<=0)).fill(0.00000001);

(*(argStruct_p->V_nk_p)).elem(find((*(argStruct_p->V_nk_p))>10)).fill(10);*/

}

/*NaN protection*/
(*(argStruct_p->V_nk_p)).elem( find_nonfinite((*(argStruct_p->V_nk_p))) ).zeros();


}	

/*To do list: 

comment out : orthogonalize_Z_ol_skinny();
enable/comment in: ZtZ_update_wrapper(). Within this function can play with the SVD of the Z_ol. Then compute the ZtZ of it. Then follows the update to make Z_ol converge accordingly; result: Z_ol is more orthogonal than before the update. 

*/
static void post_update_Zol(mxArray *plhs[], arg_struct_t* argStruct_p, int n){



/*Call projfunc_wrapper() */
/*projfunc_wrapper_Zol(argStruct_p->Z_ol_p, argStruct_p->Z_ol_mxArray_p);*/

// Subtract off the current mean
/*accum_Z_ol_local=accu(*(argStruct_p->Z_ol_p));

*(argStruct_p->Z_ol_p)=(*(argStruct_p->Z_ol_p))-(accum_Z_ol_local/(((double)O_static)*((double)L_static)))*ones(O_static, L_static);

// Set new mean forcibly
*(argStruct_p->Z_ol_p)=(*(argStruct_p->Z_ol_p))+0.95*(accum_Z_ol/(((double)O_static)*((double)L_static)))*ones(O_static, L_static);*/

// Clip the negative elements again..
(*(argStruct_p->Z_ol_p)).elem(find((*(argStruct_p->Z_ol_p))<=0)).fill(0.00000001);

(*(argStruct_p->Z_ol_p)).elem(find((*(argStruct_p->Z_ol_p))>1)).fill(1);

/*Perform symmetric NMF style clustering */
/*Z_cluster_update_wrapper(argStruct_p);*/

if (L_value_main_objective>(L_value_threshold_local)){

ZtZ_update_wrapper(argStruct_p->Z_ol_p, argStruct_p->Z_ol_mxArray_p, plhs);

// Clip the negative elements again..
(*(argStruct_p->Z_ol_p)).elem(find((*(argStruct_p->Z_ol_p))<=0)).fill(0.00000001);

(*(argStruct_p->Z_ol_p)).elem(find((*(argStruct_p->Z_ol_p))>1)).fill(1);

}

/*Call the new function to properly scale the cols of T_fk */
/*check_Z_pull_up_weak_clusters(argStruct_p->Z_ol_p);*/
/*Z_pull_up_vector_elements_generalized_gaussian_wrapper(argStruct_p->Z_ol_p);

// Clip the negative elements again..
(*(argStruct_p->Z_ol_p)).elem(find((*(argStruct_p->Z_ol_p))<=0)).fill(0.00000001);

(*(argStruct_p->Z_ol_p)).elem(find((*(argStruct_p->Z_ol_p))>100)).fill(100);*/

Z_weighted_scaling(argStruct_p->Z_ol_p);

(*(argStruct_p->Z_ol_p)).elem(find((*(argStruct_p->Z_ol_p))<=0)).fill(0.00000001);

(*(argStruct_p->Z_ol_p)).elem(find((*(argStruct_p->Z_ol_p))>1)).fill(1);

/*check_Z_scale_down_Y_low_threshold(argStruct_p->Z_ol_p, argStruct_p->Y_lk_p);

(*(argStruct_p->Z_ol_p)).elem(find((*(argStruct_p->Z_ol_p))<=0)).fill(0.00000001);

(*(argStruct_p->Z_ol_p)).elem(find((*(argStruct_p->Z_ol_p))>1)).fill(1);

check_Z_scale_down_Y_high_threshold(argStruct_p->Z_ol_p, argStruct_p->Y_lk_p);

(*(argStruct_p->Z_ol_p)).elem(find((*(argStruct_p->Z_ol_p))<=0)).fill(0.00000001);

(*(argStruct_p->Z_ol_p)).elem(find((*(argStruct_p->Z_ol_p))>1)).fill(1);*/

normalize_Z_ol_frob_norm(argStruct_p->Z_ol_p);

(*(argStruct_p->Z_ol_p)).elem(find((*(argStruct_p->Z_ol_p))<=0)).fill(0.00000001);

(*(argStruct_p->Z_ol_p)).elem(find((*(argStruct_p->Z_ol_p))>1)).fill(1);

#ifdef ENABLE_MATLAB_ENGINE_PLOTTING
/*armaSetPr(plhs[2], (*(argStruct_p->Z_ol_p)));
engPutVariable(mlEngine_p, "Z_ol", plhs[2]);
engEvalString(mlEngine_p, "[Z_ol]=orthogonalize_Z_ol_skinny(Z_ol);");
plhs[2]=engGetVariable(mlEngine_p, "Z_ol");
(*(argStruct_p->Z_ol_p)).set_real(armaGetPr(plhs[2], true)); 	

// Clip the negative elements again..
(*(argStruct_p->Z_ol_p)).elem(find((*(argStruct_p->Z_ol_p))<=0)).fill(0.00000001);*/
#endif


}	

/*To do list: 

*/
static void post_update_Ylk(mxArray *plhs[], arg_struct_t* argStruct_p, int n){




/*projfunc_wrapper_Ylk(argStruct_p->Y_lk_p, argStruct_p->Y_lk_mxArray_p);*/

/*orthogonalize_Ylk_square_wrapper(argStruct_p->Y_lk_p, argStruct_p->Y_lk_mxArray_p);*/

/*relate_significant_Zl_clusters_to_Yl_Clusters(argStruct_p->Z_ol_p, argStruct_p->Y_lk_p);*/

/*Perform symmetric NMF style clustering */
/*Y_cluster_update_wrapper(argStruct_p);*/

/*Can enable the .fill here or after normalize Y_lk if things get numerically unstable. */
(*(argStruct_p->Y_lk_p)).elem(find((*(argStruct_p->Y_lk_p))<=0)).fill(0.00000001);

(*(argStruct_p->Y_lk_p)).elem(find((*(argStruct_p->Y_lk_p))>10)).fill(10);

/*if (L_value_main_objective>(L_value_threshold_local)){*/

/*(*(argStruct_p->Y_lk_p)).elem(find((*(argStruct_p->Y_lk_p))>100)).fill(100);*/

/*normalize_Y_lk(argStruct_p->Y_lk_p, argStruct_p->V_nk_p);*/

if (L_value_main_objective>(L_value_threshold_local)){

YYt_update_wrapper(argStruct_p->Y_lk_p, argStruct_p->Y_lk_mxArray_p, plhs);

(*(argStruct_p->Y_lk_p)).elem( find_nonfinite((*(argStruct_p->Y_lk_p))) ).zeros();

// Clip the negative elements again..
(*(argStruct_p->Y_lk_p)).elem(find((*(argStruct_p->Y_lk_p))<=0)).fill(0.00000001);

(*(argStruct_p->Y_lk_p)).elem(find((*(argStruct_p->Y_lk_p))>10)).fill(10);

}

/*(*(argStruct_p->Y_lk_p)).print("Y_lk after YYt_update_wrapper:");*/

/*Y_weighted_scaling(argStruct_p->Y_lk_p);

// Clip the negative elements again..
(*(argStruct_p->Y_lk_p)).elem(find((*(argStruct_p->Y_lk_p))<=0)).fill(0.00000001);

(*(argStruct_p->Y_lk_p)).elem(find((*(argStruct_p->Y_lk_p))>1)).fill(1);*/

/*Y_pull_up_vector_elements_generalized_gaussian_wrapper(argStruct_p->Y_lk_p);

// Clip the negative elements again..
(*(argStruct_p->Y_lk_p)).elem(find((*(argStruct_p->Y_lk_p))<=0)).fill(0.00000001);

(*(argStruct_p->Y_lk_p)).elem(find((*(argStruct_p->Y_lk_p))>10)).fill(10);

balance_energies_of_Ltarget_most_significant_clusters(argStruct_p->Y_lk_p);

// Clip the negative elements again..
(*(argStruct_p->Y_lk_p)).elem(find((*(argStruct_p->Y_lk_p))<=0)).fill(0.00000001);

(*(argStruct_p->Y_lk_p)).elem(find((*(argStruct_p->Y_lk_p))>10)).fill(10);*/

/*}*/

/*if (L_value_main_objective<(L_value_threshold_local2)){*/

/*K MEANS CLUSTERING*/
k_means_clustering_entry(argStruct_p->Z_ol_p, argStruct_p->Y_lk_p);

// Clip the negative elements again..
(*(argStruct_p->Y_lk_p)).elem(find((*(argStruct_p->Y_lk_p))<=0)).fill(0.00000001);

(*(argStruct_p->Y_lk_p)).elem(find((*(argStruct_p->Y_lk_p))>10)).fill(10);

(*(argStruct_p->Z_ol_p)).elem(find((*(argStruct_p->Z_ol_p))<=0)).fill(0.00000001);

(*(argStruct_p->Z_ol_p)).elem(find((*(argStruct_p->Z_ol_p))>1)).fill(1);

//Frobenius norm normalization
normalize_Z_ok_compensate_TV(argStruct_p->Z_ol_p, argStruct_p->Y_lk_p, argStruct_p->T_fk_p, argStruct_p->V_nk_p);

// Clip the negative elements again..
(*(argStruct_p->Y_lk_p)).elem(find((*(argStruct_p->Y_lk_p))<=0)).fill(0.00000001);

(*(argStruct_p->Y_lk_p)).elem(find((*(argStruct_p->Y_lk_p))>10)).fill(10);

(*(argStruct_p->Y_lk_p)).elem( find_nonfinite((*(argStruct_p->Y_lk_p))) ).zeros();

(*(argStruct_p->Z_ol_p)).elem(find((*(argStruct_p->Z_ol_p))<=0)).fill(0.00000001);

(*(argStruct_p->Z_ol_p)).elem(find((*(argStruct_p->Z_ol_p))>1)).fill(1);

(*(argStruct_p->Z_ol_p)).elem( find_nonfinite((*(argStruct_p->Z_ol_p))) ).zeros();

(*(argStruct_p->T_fk_p)).elem(find((*(argStruct_p->T_fk_p))<=0)).fill(0.00000001);

(*(argStruct_p->T_fk_p)).elem(find((*(argStruct_p->T_fk_p))>100)).fill(100);

(*(argStruct_p->T_fk_p)).elem( find_nonfinite((*(argStruct_p->T_fk_p))) ).zeros();

normalize_Y_lk_compensate_Z_Lstatic_eq_Ltarget(argStruct_p->Y_lk_p, argStruct_p->Z_ol_p);

// Clip the negative elements again..
(*(argStruct_p->Y_lk_p)).elem(find((*(argStruct_p->Y_lk_p))<=0)).fill(0.00000001);

(*(argStruct_p->Y_lk_p)).elem(find((*(argStruct_p->Y_lk_p))>1)).fill(1);

(*(argStruct_p->Z_ol_p)).elem(find((*(argStruct_p->Z_ol_p))<=0)).fill(0.00000001);

(*(argStruct_p->Z_ol_p)).elem(find((*(argStruct_p->Z_ol_p))>1)).fill(1);

/*}*/

/*Populate this for Z update*/
Z_inmat_real_Ykl=trans((*(argStruct_p->Y_lk_p)));

}	


/*To do list: 

*/
static void post_update_Phi_S(mxArray *plhs[], arg_struct_t* argStruct_p, int n){

if (!N_ITERATIONS_TV_bool){

/*V_nk_update_Phi_S_wrapper(plhs, argStruct_p);*/

/*Vnk_update_Phi_S(argStruct_p->V_nk_p, argStruct_p->expj_Phi_S_fkn_p, argStruct_p);*/

}

}	

static void plot_Xhat_fnm_wrapper(mxArray* plhs[]){

engPutVariable(mlEngine_p, "Xhat_fnm", plhs[0]);
engEvalString(mlEngine_p, "plot_Xhat_fnm(Xhat_fnm);");

}

static void plot_Xhat_Xtilde_fnm_wrapper(mxArray* plhs[], const mxArray *prhs[]){

engPutVariable(mlEngine_p, "Xtilde_fnm", prhs[0]);
engPutVariable(mlEngine_p, "Xhat_fnm", plhs[0]);
engEvalString(mlEngine_p, "plot_Xhat_Xtilde_fnm(Xhat_fnm, Xtilde_fnm);");

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

#ifndef ENABLE_MATLAB_ENGINE_PLOTTING
mexCallMATLAB(0, NULL, 1, &plhs[0], "plot_Xhat_fnm");
#endif

#ifdef ENABLE_MATLAB_ENGINE_PLOTTING
engPutVariable(mlEngine_p, "Xhat_fnm", plhs[0]);
#endif

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
void plot_and_update_Tfk(mxArray *plhs[], const mxArray *prhs[], arg_struct_t* argStruct_p, int n){

double accum_T_fk_local; 

/*Temp array of pointers of mxArray*/
mxArray* plot_array_local[2];

/*colvec::fixed<1> scalar_local;*/

(scalar_global[0])(0)=1;

armaSetPr(scalar_global_p[0], scalar_global[0]);

plot_array_local[0]=plhs[4];
plot_array_local[1]=scalar_global_p[0];

mexPrintf("About to compute T_fk \n");

/*Compute T_fk  ------------------------------------------------------------------------------------ */	
update_engine(argStruct_p, 2);

post_update_Tfk(plhs, argStruct_p, n);

#ifdef 	MEX_UNIT_TESTING
armaSetPr(plhs[4], (*(argStruct_p->T_fk_p)));
/*mexCallMATLAB(0, NULL, 1, &plhs[4], "plot_T_fk");*/

#ifndef ENABLE_MATLAB_ENGINE_PLOTTING
mexCallMATLAB(0, NULL, 2, plot_array_local, "plot_T_fk");
#endif

#ifdef ENABLE_MATLAB_ENGINE_PLOTTING

/*plot_Xhat_fnm_wrapper(plhs);*/
plot_Xhat_Xtilde_fnm_wrapper(plhs, prhs);

engPutVariable(mlEngine_p, "T_fk", plhs[4]);
engEvalString(mlEngine_p, "plot_T_fk(T_fk,1);");
#endif

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

post_update_Vnk(plhs, argStruct_p, n);

#ifdef 	MEX_UNIT_TESTING
armaSetPr(plhs[5], (*(argStruct_p->V_nk_p)));

#ifndef ENABLE_MATLAB_ENGINE_PLOTTING
mexCallMATLAB(0, NULL, 1, &plhs[5], "plot_V_nk");
#endif

#ifdef ENABLE_MATLAB_ENGINE_PLOTTING
engPutVariable(mlEngine_p, "V_nk", plhs[5]);
engEvalString(mlEngine_p, "plot_V_nk(V_nk);");

plot_Xhat_fnm_wrapper(plhs);

#endif

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

post_update_Zol(plhs, argStruct_p, n);

/*#ifdef 	MEX_UNIT_TESTING
armaSetPr(plhs[2], (*(argStruct_p->Z_ol_p)));

#ifndef ENABLE_MATLAB_ENGINE_PLOTTING
mexCallMATLAB(0, NULL, 1, &plhs[2], "plot_Z_ol");
#endif

#ifdef ENABLE_MATLAB_ENGINE_PLOTTING
engPutVariable(mlEngine_p, "Z_ol", plhs[2]);
engEvalString(mlEngine_p, "plot_Z_ol(Z_ol);");
#endif

#endif

update_engine_hfk_Af_LS_estimator(argStruct_p, 2);

post_update_Zol(plhs, argStruct_p, n);*/

/*update_engine_covariance_cost(argStruct_p, 1);

update_engine_covariance_cost(argStruct_p, 2);*/

#ifdef 	MEX_UNIT_TESTING
armaSetPr(plhs[2], (*(argStruct_p->Z_ol_p)));

#ifndef ENABLE_MATLAB_ENGINE_PLOTTING
mexCallMATLAB(0, NULL, 1, &plhs[2], "plot_Z_ol");
#endif

#ifdef ENABLE_MATLAB_ENGINE_PLOTTING
engPutVariable(mlEngine_p, "Z_ol", plhs[2]);
engEvalString(mlEngine_p, "plot_Z_ol(Z_ol);");
#endif

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

post_update_Ylk(plhs, argStruct_p, n);

/*#ifdef 	MEX_UNIT_TESTING
armaSetPr(plhs[3], (*(argStruct_p->Y_lk_p)));

#ifndef ENABLE_MATLAB_ENGINE_PLOTTING
mexCallMATLAB(0, NULL, 1, &plhs[3], "plot_Y_lk");
#endif

#ifdef ENABLE_MATLAB_ENGINE_PLOTTING
engPutVariable(mlEngine_p, "Y_lk", plhs[3]);
engEvalString(mlEngine_p, "plot_Y_lk(Y_lk);");

plot_Xhat_fnm_wrapper(plhs);

#endif

#endif

update_engine_hfk_Af_LS_estimator(argStruct_p, 3);

post_update_Ylk(plhs, argStruct_p, n);*/

/*update_engine_covariance_cost(argStruct_p, 1);

update_engine_covariance_cost(argStruct_p, 3);*/

plot_Zok_wrapper(plhs);

#ifdef 	MEX_UNIT_TESTING
armaSetPr(plhs[3], (*(argStruct_p->Y_lk_p)));

#ifndef ENABLE_MATLAB_ENGINE_PLOTTING
mexCallMATLAB(0, NULL, 1, &plhs[3], "plot_Y_lk");
#endif

#ifdef ENABLE_MATLAB_ENGINE_PLOTTING
engPutVariable(mlEngine_p, "Y_lk", plhs[3]);
engEvalString(mlEngine_p, "plot_Y_lk(Y_lk);");

plot_Xhat_fnm_wrapper(plhs);

#endif

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

/*populate_W_fom_cx(argStruct_p->W_fom_p, argStruct_p->W_fom_cx_p, argStruct_p->Phi_W_fom_p);

update_engine_hfk_Af_LS_estimator(argStruct_p, 1);*/

#ifdef BOTTOM_UP_MEX_FLAG

/*W_callMATLAB_wrapper();*/

#endif

populate_W_fom_cx(argStruct_p->W_fom_p, argStruct_p->W_fom_cx_p, argStruct_p->Phi_W_fom_p);

/*magnitude_square_rooted_processing_W_wrapper(argStruct_p, plhs);*/

#ifdef 	MEX_UNIT_TESTING
armaSetCubePr(plhs[1], (*(argStruct_p->W_fom_p)));

#ifndef ENABLE_MATLAB_ENGINE_PLOTTING
mexCallMATLAB(0, NULL, 1, &plhs[1], "plot_W_fom");
#endif

#ifdef ENABLE_MATLAB_ENGINE_PLOTTING
engPutVariable(mlEngine_p, "W_fom", plhs[1]);
engEvalString(mlEngine_p, "plot_W_fom(W_fom);");

#endif

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

/*Temp array of pointers of mxArray*/
mxArray* plot_array_local[3];

/*Make this global, and reuse it?*/
/*colvec::fixed<1> scalar_global[2];*/

(scalar_global[0])(0)=1;
(scalar_global[1])(0)=1; 

armaSetPr(scalar_global_p[0], scalar_global[0]);
armaSetPr(scalar_global_p[1], scalar_global[1]);

plot_array_local[0]=plhs[7];
plot_array_local[1]=scalar_global_p[0];
plot_array_local[2]=scalar_global_p[1];

/*mexCallMATLAB_Phi_S(plhs, argStruct_p);*/

/*Compute Phi_S_fnk  ------------------------------------------------------------------------------- */	
update_engine(argStruct_p, 5); 

post_update_Phi_S(plhs, argStruct_p, n);

/*Phi_S_L1_norm_V_checker(argStruct_p->expj_Phi_S_nkf_p, argStruct_p->expj_Phi_S_fkn_p, argStruct_p->V_nk_p);*/

/*Plot Phi_S_fnk??*/
#ifdef 	MEX_UNIT_TESTING
armaSetCubeCx(plhs[7], (*(argStruct_p->expj_Phi_S_fkn_p)));
/*mexCallMATLAB(0, NULL, 1, &plhs[7], "plot_Phi_S_fnk");*/

#ifndef ENABLE_MATLAB_ENGINE_PLOTTING
mexCallMATLAB(0, NULL, 3, plot_array_local, "plot_Phi_S_fnk");
#endif

#ifdef ENABLE_MATLAB_ENGINE_PLOTTING
engPutVariable(mlEngine_p, "expj_Phi_S_fkn", plhs[7]);
engEvalString(mlEngine_p, "plot_Phi_S_fnk(expj_Phi_S_fkn, 1, 1);");
#endif

#endif

/* Print status to MATLAB */
some_mat(0)=(double)n; 
some_mat(1)=(double)5;
armaSetPr(input_arg_p, some_mat);
mexCallMATLAB(0, NULL, 1, &input_arg_p, "nguyen_2015_print");

}	