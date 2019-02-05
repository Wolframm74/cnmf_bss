#include "local_inc.hpp"

/*cx_cube::fixed<F_static, N_static, O_static> P_Xh1_wrt_Zol_fno_l[L_static];
cx_cube::fixed<F_static, N_static, O_static> P_Xh2_wrt_Zol_fno_l[L_static];*/

/*Function and data moved to common_module.cpp*/
/*cx_mat::fixed<N_static*K_static, F_static> temptensor_NKxF;
cx_cube::fixed<F_static, N_static, K_static> rotated_expj_Phi_S_FNK;

static void rotate_expj_Phi_S_nkf_p_to_FNK(cx_cube* expj_Phi_S_nkf_p){

//copy to temptensor_NKxF

//trans temptensor_NKxF into output result FxNK	

//copy result to rotated_expj_Phi_S_FNK


}*/


cx_mat::fixed<F_static, N_static> P_Xh1_wrt_Zol_fn;

mat::fixed<F_static, N_static> P_G1_wrt_Zol_fn;

mat::fixed<F_static, N_static> P_H1_wrt_Zol_fn;

mat::fixed<F_static, N_static> P_H1_wrt_Zol_num_fn;

mat::fixed<F_static, N_static> P_arg_Xh1_wrt_Zol_fn;

mat::fixed<O_static, L_static> Pwrt_C_Zol_ol;

void Z_update2_Zol(mat* Z_ol_p){

double epsilon_Zol;

epsilon_Zol=0.01;

(*Z_ol_p)=(*Z_ol_p)-(epsilon_Zol)*(Pwrt_C_Zol_ol);

}

static void compute_Pwrt_C_Zol_ol(cx_cube* W_fom_cx_p, mat* Y_lk_p, mat* T_fk_p, mat* V_nk_p, mat* Z_ol_p, cx_cube* Xhat_fnm_p, cx_cube* Xtilde_fnm_p, cx_cube* expj_Phi_S_nkf_p){

int l_index, o_index;

int k_index;

int m_index;

Pwrt_C_Zol_ol.zeros();

rotate_expj_Phi_S_nkf_p_to_FNK(expj_Phi_S_nkf_p);

for (m_index=0; m_index<M_static; m_index++){

/*Calling this will populate G1_outmat_fn, H1_outmat_fn, sqrt_denmat_1_fn*/
populate_G_H_fn_channel1(m_index, Xhat_fnm_p);

/*Elementise wise FxNxK product of the thing you computed with the previous function vs the difference between Pwrt_arg_Xh1_fnk and Pwrt_arg_Xh2_fnk */
dummy_mat_real_FN.zeros();

dummy_mat_real_FN=compute_arg_X_FN(Xhat_fnm_p, m_index);

/*dummy_mat_real_FN=dummy_mat_real_FN+compute_arg_X_FN(Xtilde_fnm_p, m_index_2);*/

dummy_mat_real_FN=dummy_mat_real_FN-compute_arg_X_FN(Xtilde_fnm_p, m_index);

/*dummy_mat_real_FN=dummy_mat_real_FN-compute_arg_X_FN(Xhat_fnm_p, m_index_2);*/

dummy_mat_real_FN=2*dummy_mat_real_FN;

second_dummy_mat_real_FN=real((*Xhat_fnm_p).slice(m_index));

third_dummy_mat_real_FN=2*ones_mat_FxN/(square(H1_outmat_fn)+square(G1_outmat_fn));

for (o_index=0; o_index<O_static; o_index++){

	for (l_index=0; l_index<L_static; l_index++){

		scalar_global_real=(double)o_index;

		scalar_global_real.print("Z_update_2, line 83: o_index=");

		scalar_global_real=(double)l_index;		

		scalar_global_real.print("Z_update_2, line 87: l_index=");

		P_Xh1_wrt_Zol_fn.zeros();

		/*P_Xh2_wrt_Zol_fn.zeros();*/

		for (k_index=0; k_index<K_static; k_index++){

		P_Xh1_wrt_Zol_fn=P_Xh1_wrt_Zol_fn+((*Y_lk_p)(l_index, k_index))*(((*W_fom_cx_p).slice(m_index).col(o_index)%(*T_fk_p).col(k_index))*(trans((*V_nk_p).col(k_index))))%(rotated_expj_Phi_S_FNK.slice(k_index));

		}		

		P_G1_wrt_Zol_fn=imag(P_Xh1_wrt_Zol_fn);

		P_H1_wrt_Zol_num_fn=second_dummy_mat_real_FN%real(P_Xh1_wrt_Zol_fn)+G1_outmat_fn%imag(P_Xh1_wrt_Zol_fn);

		P_H1_wrt_Zol_fn=P_H1_wrt_Zol_num_fn/sqrt_denmat_1_fn;
		
		P_arg_Xh1_wrt_Zol_fn=third_dummy_mat_real_FN%(P_G1_wrt_Zol_fn%H1_outmat_fn-G1_outmat_fn%P_H1_wrt_Zol_fn);

		Pwrt_C_Zol_ol(o_index, l_index)=Pwrt_C_Zol_ol(o_index, l_index)+accu(dummy_mat_real_FN%(P_arg_Xh1_wrt_Zol_fn));

	}

}

}

Z_update2_Zol(Z_ol_p);

}


void Z_update2_entry_print(void){


/*m=1st channel*/

scalar_global_cx=as_scalar(accu(P_Xh1_wrt_Zol_fn));
scalar_global_cx.print("P_Xh1_wrt_Zol_fn");

scalar_global_real=as_scalar(accu(P_arg_Xh1_wrt_Zol_fn));
scalar_global_real.print("P_arg_Xh1_wrt_Zol_fn");


/*Function local data: m=1:2 F_staticxK_static matrices*/
scalar_global_cx=as_scalar(accu(dummy_mat_H_fk_1));
scalar_global_real.print("dummy_mat_H_fk_1");

scalar_global_cx=as_scalar(accu(dummy_mat_H_fl_1));
scalar_global_cx.print("dummy_mat_H_fl_1");

scalar_global_real=as_scalar(accu(P_H1_wrt_Zol_fn));	
scalar_global_real.print("P_H1_wrt_Zol_fn");

scalar_global_real=as_scalar(accu(P_H1_wrt_Zol_num_fn)); 	/*This should come in as real(*Xhat_fnm_p.slice(m_index)) spreaded over k*/
scalar_global_real.print("P_H1_wrt_Zol_num_fn");

/*Function local data*/

scalar_global_real=as_scalar(accu(G1_outmat_fn));
scalar_global_real.print("G1_outmat_fn");

scalar_global_real=as_scalar(accu(H1_outmat_fn));
scalar_global_real.print("H1_outmat_fn");

/*Function local data*/

scalar_global_real=as_scalar(accu(P_G1_wrt_Zol_fn));
scalar_global_real.print("P_G1_wrt_Zol_fn");

scalar_global_real=as_scalar(accu(outmat_FN));
scalar_global_real.print("outmat_FN");

scalar_global_real=as_scalar(accu(nummat_FN));
scalar_global_real.print("nummat_FN");

scalar_global_real=as_scalar(accu(denmat_FN));
scalar_global_real.print("denmat_FN");

/*Function local data*/
scalar_global_real=as_scalar(accu(Pwrt_C_Zol_ol));
scalar_global_real.print("Pwrt_C_Zol_ol");

}

void Z_update2_entry(arg_struct_t* argStruct_p){

compute_Pwrt_C_Zol_ol(argStruct_p->W_fom_cx_p, argStruct_p->Y_lk_p, argStruct_p->T_fk_p, argStruct_p->V_nk_p, argStruct_p->Z_ol_p, argStruct_p->Xhat_fnm_p, argStruct_p->Xtilde_fnm_p, argStruct_p->expj_Phi_S_nkf_p);

Z_update2_entry_print();

}