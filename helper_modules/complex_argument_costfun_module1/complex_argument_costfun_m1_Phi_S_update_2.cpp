#include "local_inc.hpp"

/*Dummy mats*/
static cx_mat::fixed<F_static, N_static> dummy_mat_cx_FN;
mat::fixed<F_static, N_static> dummy_mat_real_FN;

/*m=1st channel*/

static cx_cube::fixed<F_static, N_static, K_static> Pwrt_Xh1_fnk;
static cube::fixed<F_static, N_static, K_static> Pwrt_arg_Xh1_fnk;

/*m=2nd channel*/

/*Function local data: m=1:2 F_staticxK_static matrices*/
/*MOVED TO COMMON_MODULE.CPP*/
/*cx_mat::fixed<F_static, K_static> dummy_mat_H_fk_1;
cx_mat::fixed<F_static, K_static> dummy_mat_H_fk_2;

cx_mat::fixed<F_static, L_static> dummy_mat_H_fl_1;
cx_mat::fixed<F_static, L_static> dummy_mat_H_fl_2;*/


static void populate_Pwrt_Xh1_fnk(int m_index, cx_cube* W_fom_cx_p, mat* Z_ol_p, mat* Y_lk_p, cx_cube* expj_Phi_S_fkn_p, mat* T_fk_p, mat* V_nk_p){
/*To answer the question why would you not just hardcore m_index to 1 for this function and 2 for the other: Because sometimes you might need to use m_indices not equal to just 1 and 2.
In the case that  you have >=2 microphones. 
*/

int k_index, n_index; 
cx_double fill_value;

fill_value.imag()=1; 

populate_dummy_mat_H_fk_1(m_index, W_fom_cx_p, Z_ol_p, Y_lk_p);

/*Populate dummy_mat_H_fl*/
/*dummy_mat_H_fl_1=((*W_fom_cx_p).slice(m_index))*(*Z_ol_p);*/

/*Populate dummy_mat_H_fk*/
/*dummy_mat_H_fk_1=dummy_mat_H_fl_1*(*Y_lk_p);*/

/*Spread it across n=0:N-1, by doing a kron() with a ones_row_1xN */	

/*Do the elementwise multiplication of it with T_fk, V_nk, and expj_Phi_S */


for (k_index=0; k_index<K_static; k_index++){

	dummy_mat_cx_FN.fill(fill_value);

	dummy_mat_cx_FN=dummy_mat_cx_FN%(kron(dummy_mat_H_fk_1.col(k_index), ones_row_1xN));

	dummy_mat_cx_FN=dummy_mat_cx_FN%((*T_fk_p).col(k_index)*trans((*V_nk_p).col(k_index)));

	Pwrt_Xh1_fnk.slice(k_index)=dummy_mat_cx_FN;

	for (n_index=0; n_index<N_static; n_index++){

	Pwrt_Xh1_fnk.slice(k_index).col(n_index)=Pwrt_Xh1_fnk.slice(k_index).col(n_index)%(*expj_Phi_S_fkn_p).slice(n_index).col(k_index);

	}

}


}


/*Function local data for functions populate_Pwrt_H1_fnk and populate_Pwrt_H2_fnk*/
/*This local function data will need to be moved up, otherwise will generate a compiler error*/
static cube::fixed<F_static, N_static, K_static> Pwrt_H1_fnk;	/*Can be used as a dummy variable while you are computing the num and den functions at first*/	/*This should come in as imag(*Xhat_fnm_p.slice(m_index)) spreaded over k*/
static cube::fixed<F_static, N_static, K_static> Pwrt_H1_num_fnk; 	/*This should come in as real(*Xhat_fnm_p.slice(m_index)) spreaded over k*/
/*cube::fixed<F_static, N_static, K_static> Pwrt_H1_den_fnk;*/

/*cube::fixed<F_static, N_static, K_static> Pwrt_H2_den_fnk;*/

/*Function local data*/

static cube::fixed<F_static, N_static, K_static> G1_fn_spreaded_fnk;
static cube::fixed<F_static, N_static, K_static> H1_fn_spreaded_fnk;

static cube::fixed<F_static, N_static, K_static> sqrt_denmat_1_fnk;

static void populate_G_H_fn_channel1_local(int m_index, cx_cube* Xhat_fnm_p){

int k_index;

Pwrt_H1_fnk.slice(0)=imag((*Xhat_fnm_p).slice(m_index));	
Pwrt_H1_num_fnk.slice(0)=real((*Xhat_fnm_p).slice(m_index));

G1_fn_spreaded_fnk.slice(0)=imag((*Xhat_fnm_p).slice(m_index));

/*Compute .slice(0) for the denominator function*/
dummy_mat_real_FN=sqrt(square(Pwrt_H1_num_fnk.slice(0))+square(Pwrt_H1_fnk.slice(0)));
/*May use the above result again later*/

sqrt_denmat_1_fnk.slice(0)=dummy_mat_real_FN;

H1_fn_spreaded_fnk.slice(0)=dummy_mat_real_FN+Pwrt_H1_num_fnk.slice(0);

/*now that you have the .slice(0) functions for both num and den, spread the result across the k=1:K-1 slices*/
for (k_index=0; k_index<K_static-1; k_index++){

G1_fn_spreaded_fnk.slice(k_index)=G1_fn_spreaded_fnk.slice(0);

H1_fn_spreaded_fnk.slice(k_index)=H1_fn_spreaded_fnk.slice(0);

sqrt_denmat_1_fnk.slice(k_index)=sqrt_denmat_1_fnk.slice(0);

/*Precalculate stuff for one of the next serial functions. Also need to spread these across k*/
Pwrt_H1_fnk.slice(k_index)=Pwrt_H1_fnk.slice(0);

Pwrt_H1_num_fnk.slice(k_index)=Pwrt_H1_num_fnk.slice(0);

}

scalar_global_real=as_scalar(accu(Pwrt_H1_fnk));	
scalar_global_real.print("Pwrt_H1_fnk");

}

/*Function local data*/

cube::fixed<F_static, N_static, K_static> Pwrt_G1_fnk;

static void populate_Pwrt_G1_fnk(void){

/*Take the imaginary part of Pwrt_Xh1_fnk*/
Pwrt_G1_fnk=imag(Pwrt_Xh1_fnk);


}

/*cube::fixed<F_static, N_static, K_static> Pwrt_H1_fnk_finite_ceiling;*/

static void populate_Pwrt_H1_fnk(int m_index, cx_cube* Xhat_fnm_p){

/*Take the imaginary and real parts of the various matrices, tensors to compute this*/

/*Can use Pwrt_H1_fnk, Pwrt_H2_fnk as dummy vars to temp store things */

/*Compute the num function*/	
Pwrt_H1_fnk=Pwrt_H1_fnk%imag(Pwrt_Xh1_fnk);
Pwrt_H1_num_fnk=Pwrt_H1_num_fnk%real(Pwrt_Xh1_fnk);

/*Compute the den function => It's already been computed and saved within Hm_fn_spreaded_fnk */	

Pwrt_H1_fnk=Pwrt_H1_num_fnk+Pwrt_H1_fnk;

/*ones_cube_FNK=sign(H1_fn_spreaded_fnk);*/

H1_fn_spreaded_fnk.elem(find(abs(H1_fn_spreaded_fnk)==0)).fill(0.00000001);

sqrt_denmat_1_fnk.elem(find(abs(sqrt_denmat_1_fnk)==0)).fill(0.00000001);

/*H1_fn_spreaded_fnk.elem(find(abs(H1_fn_spreaded_fnk)>=100000000)).fill(100000000);*/

/*H1_fn_spreaded_fnk=abs(H1_fn_spreaded_fnk)%ones_cube_FNK;

ones_cube_FNK.ones();*/

Pwrt_H1_fnk=Pwrt_H1_fnk/sqrt_denmat_1_fnk;

Pwrt_H1_fnk.elem(find_nonfinite(Pwrt_H1_fnk)).fill(100000000);

Pwrt_H1_fnk=Pwrt_H1_fnk+real(Pwrt_Xh1_fnk);

}


/*Function local data*/

static void populate_Pwrt_arg_Xh1_fnk(void){

/*Compute (2/(h^2+g^2)), followed by the rest */

/*H1_fn_spreaded_fnk.elem(find(abs(H1_fn_spreaded_fnk)<=0)).fill(0.00000001);

G1_fn_spreaded_fnk.elem(find(abs(G1_fn_spreaded_fnk)<=0)).fill(0.00000001);*/

Pwrt_arg_Xh1_fnk=2*ones_cube_FNK/(square(H1_fn_spreaded_fnk)+square(G1_fn_spreaded_fnk));

Pwrt_arg_Xh1_fnk=Pwrt_arg_Xh1_fnk%(Pwrt_G1_fnk%H1_fn_spreaded_fnk-G1_fn_spreaded_fnk%Pwrt_H1_fnk);

}

/*mat& compute_arg_X_FN(cx_cube* X_fnm, int m_index): moved to common_module.cpp*/

/*Function local data*/
cube::fixed<F_static, N_static, K_static> Pwrt_C_phi_fnk;
cube::fixed<F_static, N_static, K_static> Pwrt_C_phi_fnk_mth_component;

static cube& compute_Pwrt_C_phi_fnk(int m_index, cx_cube* Xhat_fnm_p, cx_cube* Xtilde_fnm_p){

int k_index;

cube& Pwrt_C_phi_fnk_mth_component_ref=Pwrt_C_phi_fnk_mth_component;

Pwrt_C_phi_fnk_mth_component.zeros();

/*Elementise wise FxNxK product of the thing you computed with the previous function vs the difference between Pwrt_arg_Xh1_fnk and Pwrt_arg_Xh2_fnk */
dummy_mat_real_FN.zeros();

dummy_mat_real_FN=compute_arg_X_FN(Xhat_fnm_p, m_index);

/*dummy_mat_real_FN=dummy_mat_real_FN+compute_arg_X_FN(Xtilde_fnm_p, m_index_2);*/

/*dummy_mat_real_FN=dummy_mat_real_FN-compute_arg_X_FN(Xtilde_fnm_p, m_index);*/

/*dummy_mat_real_FN=dummy_mat_real_FN-compute_arg_X_FN(Xhat_fnm_p, m_index_2);*/

dummy_mat_real_FN=2*dummy_mat_real_FN;

/*Spead the result over k=0:K-1*/
for (k_index=0; k_index<K_static; k_index++){

Pwrt_C_phi_fnk_mth_component.slice(k_index)=dummy_mat_real_FN;

}	

Pwrt_C_phi_fnk_mth_component=Pwrt_C_phi_fnk_mth_component%(Pwrt_arg_Xh1_fnk);

return Pwrt_C_phi_fnk_mth_component_ref;

}

static cx_cube::fixed<F_static, N_static, K_static> phase_diff_cx_FNK;
static double epsilon_Phi_S_value;

static void populate_phase_diff_FNK(cx_cube* expj_Phi_S_fkn_p, cx_cube* expj_Phi_S_nkf_p){

int f_index, n_index, k_index;

epsilon_Phi_S_value=1000;

phase_diff_cx_FNK.set_real(cos(-epsilon_Phi_S_value*Pwrt_C_phi_fnk));

phase_diff_cx_FNK.set_imag(sin(-epsilon_Phi_S_value*Pwrt_C_phi_fnk));

for (f_index=0; f_index<F_static; f_index++){

	for (n_index=0; n_index<N_static; n_index++){

		for (k_index=0; k_index<K_static; k_index++){

		((*expj_Phi_S_fkn_p)(f_index, k_index, n_index))=((*expj_Phi_S_fkn_p)(f_index, k_index, n_index))*phase_diff_cx_FNK(f_index, n_index, k_index);

		((*expj_Phi_S_nkf_p)(n_index, k_index, f_index))=((*expj_Phi_S_nkf_p)(n_index, k_index, f_index))*phase_diff_cx_FNK(f_index, n_index, k_index);

		}

	}		

}



}

void Phi_S_update2_entry(void){


/*m=1st channel*/

scalar_global_cx=as_scalar(accu(Pwrt_Xh1_fnk));
scalar_global_cx.print("Pwrt_Xh1_fnk");

scalar_global_real=as_scalar(accu(Pwrt_arg_Xh1_fnk));
scalar_global_real.print("Pwrt_arg_Xh1_fnk");


/*Function local data: m=1:2 F_staticxK_static matrices*/
scalar_global_cx=as_scalar(accu(dummy_mat_H_fk_1));
scalar_global_real.print("dummy_mat_H_fk_1");

scalar_global_cx=as_scalar(accu(dummy_mat_H_fl_1));
scalar_global_cx.print("dummy_mat_H_fl_1");

scalar_global_real=as_scalar(accu(Pwrt_H1_fnk));	
scalar_global_real.print("Pwrt_H1_fnk");

scalar_global_real=as_scalar(accu(Pwrt_H1_num_fnk)); 	/*This should come in as real(*Xhat_fnm_p.slice(m_index)) spreaded over k*/
scalar_global_real.print("Pwrt_H1_num_fnk");

/*Function local data*/

scalar_global_real=as_scalar(accu(G1_fn_spreaded_fnk));
scalar_global_real.print("G1_fn_spreaded_fnk");

scalar_global_real=as_scalar(accu(H1_fn_spreaded_fnk));
scalar_global_real.print("H1_fn_spreaded_fnk");

/*Function local data*/

scalar_global_real=as_scalar(accu(Pwrt_G1_fnk));
scalar_global_real.print("Pwrt_G1_fnk");

scalar_global_real=as_scalar(accu(outmat_FN));
scalar_global_real.print("outmat_FN");

scalar_global_real=as_scalar(accu(nummat_FN));
scalar_global_real.print("nummat_FN");

scalar_global_real=as_scalar(accu(denmat_FN));
scalar_global_real.print("denmat_FN");

/*Function local data*/
scalar_global_real=as_scalar(accu(Pwrt_C_phi_fnk));
scalar_global_real.print("Pwrt_C_phi_fnk");

}

void complex_argument_costfun_m1_Phi_S_update2_entry(arg_struct_t* argStruct_p){

int m_index;

Pwrt_C_phi_fnk.zeros();

/*for (m_index=0; m_index<M_static; m_index++){*/

m_index=0; 

populate_Pwrt_Xh1_fnk(m_index, argStruct_p->W_fom_cx_p, argStruct_p->Z_ol_p, argStruct_p->Y_lk_p, argStruct_p->expj_Phi_S_fkn_p, argStruct_p->T_fk_p, argStruct_p->V_nk_p);

populate_G_H_fn_channel1_local(m_index, argStruct_p->Xhat_fnm_p);

populate_Pwrt_G1_fnk();

populate_Pwrt_H1_fnk(m_index, argStruct_p->Xhat_fnm_p);

populate_Pwrt_arg_Xh1_fnk();

Pwrt_C_phi_fnk=Pwrt_C_phi_fnk+compute_Pwrt_C_phi_fnk(m_index, argStruct_p->Xhat_fnm_p, argStruct_p->Xtilde_fnm_p);

/*}*/

populate_phase_diff_FNK(argStruct_p->expj_Phi_S_fkn_p, argStruct_p->expj_Phi_S_nkf_p);

/*Phi_S_update2_entry();*/

}


/*Supplementary modules for computing secondary cost function*/

/*Function 1: populate the partial derivative for m=1st channel. Need the caller to specify the m index. The function will use this m_index to index the correct index of W_fom
and to populate the "Hfk,1" dummy variable*/

/*Function 2: do the same for the m=2nd channel*/

/*Call the function that will take the difference of the two partial derivative tensors, as well as the compute the difference of arg() functions;
the F_staticxN_static result I think needs to be spread across a dimension of k=1:K slices.  Then dotted with the difference of partial derivative tensors.
After elementwise multiplication by 2, the result can be used for gradient descent update of Phi_S_fnk.
*/

/*In theory any of the functions in this module can be used within an external for loop so long as their call sequence up to and including 
for example the gradient descent update of Phi_S is preserved.

For example: basic use: call them for the m=1th,2nd channels

Call again for m=2, 3rd channels.

Call again for the m=1th, 3rd channels. 

*/

/*Define the tensors Pwrt_Xhat_fnm_1, Pwrt_Xhat_fnm_2
Pwrt_arg_Xhat_fnm_1, Pwrt_arg_Xhat_fnm_2
*/