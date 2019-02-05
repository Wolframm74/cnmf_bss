#include "local_inc.hpp"

static colvec::fixed<K_static> T_colvec_before_01;

void T_pull_vector_elements_towards_mean_wrapper(mat* T_fk_p){

colvec_dynamic_01=trans(sum(*T_fk_p, 0));

/*Save*/
T_colvec_before_01=colvec_dynamic_01;

pull_vector_elements_towards_mean((int) K_static, 0.01);

(*T_fk_p)=(*T_fk_p)%kron(ones_col_Fx1, trans(colvec_dynamic_01/T_colvec_before_01));

}

static colvec::fixed<K_static> T_colvec_before_02;

void T_pull_up_vector_elements_generalized_gaussian_wrapper(mat* T_fk_p){

colvec_dynamic_02=trans(sum(*T_fk_p, 0));

/*Save*/
T_colvec_before_02=colvec_dynamic_02;

pull_up_vector_elements_generalized_gaussian(8, 0.5, (int) K_static , 0.05, 1);

(*T_fk_p)=(*T_fk_p)%kron(ones_col_Fx1, trans(colvec_dynamic_02/T_colvec_before_02));

}