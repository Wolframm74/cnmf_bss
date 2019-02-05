#include "local_inc.hpp"

static cx_cube::fixed<F_static, N_static, L_static> Ytilde_out_fnl_1_channel_m_eq1;

static bool send_data_to_Matlab_Eng_and_plot_Sawada_source_reconstruction_init_flag=false;

static void send_data_to_Matlab_Eng_and_plot_Sawada_source_reconstruction_init(mxArray *plhs[]){

plhs[11]=armaCreateMxMatrix(Ytilde_out_fnl_1_channel_m_eq1.n_rows, Ytilde_out_fnl_1_channel_m_eq1.n_cols, Ytilde_out_fnl_1_channel_m_eq1.n_slices , mxDOUBLE_CLASS, mxCOMPLEX);

send_data_to_Matlab_Eng_and_plot_Sawada_source_reconstruction_init_flag=true;

}

void send_data_to_Matlab_Eng_and_plot_Sawada_source_reconstruction(mxArray *plhs[]){

if (!send_data_to_Matlab_Eng_and_plot_Sawada_source_reconstruction_init_flag){

send_data_to_Matlab_Eng_and_plot_Sawada_source_reconstruction_init(plhs);

}

/*plhs[11]=armaCreateMxMatrix(Ytilde_out_fnl_1_channel_m_eq1.n_rows, Ytilde_out_fnl_1_channel_m_eq1.n_cols, Ytilde_out_fnl_1_channel_m_eq1.n_slices , mxDOUBLE_CLASS, mxCOMPLEX);*/

armaSetCubeCx(plhs[11], Ytilde_out_fnl_1_channel_m_eq1);

engPutVariable(mlEngine_p, "Ytilde_out_fnl_1_channel_m_eq1", plhs[11]);

engEvalString(mlEngine_p, "plot_significant_l_clusters_4(Ytilde_out_fnl_1_channel_m_eq1);");

}

static cx_cube::fixed<F_static, L_static, M_static> h_flm_local;

static cx_cube::fixed<M_static, F_static, N_static> Ytilde_out_mfnl[L_static];

static cx_colvec::fixed<M_static> dummy_mixing_vec_Mx1;
static cx_colvec::fixed<M_static> Xhat_colvec_Mx1;
static cx_colvec::fixed<M_static> Xtilde_colvec_Mx1;

void Sawada_source_reconstruction_m1_entry_source_reconstruction(cx_cube* Xtilde_fnm_p, cx_cube* Xhat_fnm_p, mat* V_nk_p, mat* T_fk_p, mat* Z_ol_p, mat* Y_lk_p, cx_cube* W_fom_cx_p){

int f_index, n_index, m_index, l_index;

for (m_index=0; m_index<M_static; m_index++){

/*Populate the three way tensor h_flm*/

	h_flm_local.slice(m_index)=((*W_fom_cx_p).slice(m_index))*(*Z_ol_p);

}

for (n_index=0; n_index<N_static; n_index++){

	for (f_index=0; f_index<F_static; f_index++){

		for (m_index=0; m_index<M_static; m_index++){

			/*Populate the three pertinent Mx1 column vectors*/
			
			Xhat_colvec_Mx1(m_index)=(*Xhat_fnm_p)(f_index, n_index, m_index);
			Xtilde_colvec_Mx1(m_index)=(*Xtilde_fnm_p)(f_index, n_index, m_index);

		}

		for (l_index=0; l_index<L_static; l_index++){

			for (m_index=0; m_index<M_static; m_index++){

			dummy_mixing_vec_Mx1(m_index)=h_flm_local(f_index, l_index, m_index)*(expj_Phi_U_flnm[m_index](f_index, l_index, n_index));

			}

			(Ytilde_out_mfnl[l_index]).slice(n_index).col(f_index)=sqrt(as_scalar(accu(((*Y_lk_p).row(l_index))%((*T_fk_p).row(f_index))%((*V_nk_p).row(n_index)))))*(dummy_mixing_vec_Mx1*trans(dummy_mixing_vec_Mx1))*pinv(Xhat_colvec_Mx1*trans(Xhat_colvec_Mx1))*Xtilde_colvec_Mx1;

			Ytilde_out_fnl_1_channel_m_eq1(f_index, n_index, l_index)=Ytilde_out_mfnl[l_index](0, f_index, n_index);

		}

	}

}

}



void Sawada_source_reconstruction_m1_entry(mxArray *plhs[], arg_struct_t* argStruct_p){

Sawada_source_reconstruction_m1_entry_source_reconstruction(argStruct_p->Xtilde_fnm_p, argStruct_p->Xhat_fnm_p, argStruct_p->V_nk_p, argStruct_p->T_fk_p, argStruct_p->Z_ol_p, argStruct_p->Y_lk_p, argStruct_p->W_fom_cx_p);

send_data_to_Matlab_Eng_and_plot_Sawada_source_reconstruction(plhs);


}