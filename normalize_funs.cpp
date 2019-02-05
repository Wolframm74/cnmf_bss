#include "local_inc.hpp"

double accum_T_fk;
double accum_V_nk;
double accum_Z_ol;
double accum_W_fom;
cube* W_fom_local_p;

mat::fixed<O_static, L_static> dummy_mat_cumsum_OL; 
colvec::fixed<O_static> ones_col_Ox1;
mat::fixed<O_static, L_static> kronmat_real_OL;
mat::fixed<L_static, K_static> compensation_mat_LK; 
rowvec::fixed<K_static> ones_row_1xK;

/*L2 norm*/
void normalize_Z_ol(mat* Z_ol_p, mat* Y_lk_p){

/*First take the squares*/
dummy_mat_cumsum_OL=(*Z_ol_p)%(*Z_ol_p);

/*Ok. I think this is now fully correct. */
dummy_mat_cumsum_OL=cumsum(dummy_mat_cumsum_OL, 0);

kronmat_real_OL=kron(ones_col_Ox1, dummy_mat_cumsum_OL.row(O_static-1) );

kronmat_real_OL=sqrt(kronmat_real_OL);

#ifdef Z_UNIT_DEBUG

kronmat_real_OL.print("kronmat_real_OL:");

#endif

/*(*Z_ol_p).each_row()/dummy_mat_cumsum_OL.submat(O_static-1,0,O_static-1,L_static-1);*/

(*Z_ol_p)=(*Z_ol_p)/kronmat_real_OL;

/*Compensate for the Z_ol scaling by scaling Y_lk*/

compensation_mat_LK=kron(ones_row_1xK, trans(dummy_mat_cumsum_OL.row(O_static-1)));

compensation_mat_LK=sqrt(compensation_mat_LK);

/*(*Y_lk_p)=(*Y_lk_p)%compensation_mat_LK;*/

}


mat::fixed<N_static, K_static> dummy_mat_cumsum_NK; 
mat::fixed<N_static, K_static> kronmat_real_NK;
mat::fixed<F_static, K_static> V_compensation_mat_FK;

colvec::fixed<F_static> ones_col_Fx1;

/*L2 norm*/
void normalize_V_nk(mat* V_nk_p, mat* T_fk_p){

/*First take the squares*/
dummy_mat_cumsum_NK=(*V_nk_p)%(*V_nk_p);

/*cumsum along dim=0*/
dummy_mat_cumsum_NK=cumsum(dummy_mat_cumsum_NK, 0);

kronmat_real_NK=kron(ones_col_Nx1, dummy_mat_cumsum_NK.row(N_static-1));	/*ones_col_Nx1 is apparantly defined somewhere already: W_update*/

kronmat_real_NK=sqrt(kronmat_real_NK);

(*V_nk_p)=(*V_nk_p)/kronmat_real_NK;

/*Compensate*/
V_compensation_mat_FK=kron(ones_col_Fx1, dummy_mat_cumsum_NK.row(N_static-1));

V_compensation_mat_FK=sqrt(V_compensation_mat_FK);

/*(*T_fk_p)=(*T_fk_p)%V_compensation_mat_FK;*/

}

/*L2 norm*/
void normalize_V_nk_2(mat* V_nk_p, mat* Y_lk_p){

/*First take the squares*/
dummy_mat_cumsum_NK=(*V_nk_p)%(*V_nk_p);

/*cumsum along dim=0*/
dummy_mat_cumsum_NK=cumsum(dummy_mat_cumsum_NK, 0);

kronmat_real_NK=kron(ones_col_Nx1, dummy_mat_cumsum_NK.row(N_static-1));	/*ones_col_Nx1 is apparantly defined somewhere already: W_update*/

kronmat_real_NK=sqrt(kronmat_real_NK);

(*V_nk_p)=(*V_nk_p)/kronmat_real_NK;

/*Compensate*/
compensation_mat_LK=kron(ones_col_Lx1, dummy_mat_cumsum_NK.row(N_static-1));

compensation_mat_LK=sqrt(compensation_mat_LK);

(*Y_lk_p)=(*Y_lk_p)%compensation_mat_LK;

}

mat::fixed<L_static, K_static> dummy_mat_cumsum_LK;
colvec::fixed<L_static> ones_col_Lx1;
mat::fixed<L_static, K_static> kronmat_real_LK;
mat::fixed<F_static, K_static> Y_compensation_mat_FK;

/*Arithmetic mean*/
void normalize_Y_lk(mat* Y_lk_p, mat* T_fk_p){

/*Cumsum along dim=0*/
dummy_mat_cumsum_LK=cumsum((*Y_lk_p), 0);

kronmat_real_LK=kron(ones_col_Lx1, dummy_mat_cumsum_LK.row(L_static-1));	/*Check if ones_col_Lx1 is defined elsewhere already */

(*Y_lk_p)=(*Y_lk_p)/kronmat_real_LK;

/*Compensate*/
Y_compensation_mat_FK=kron(ones_col_Fx1, dummy_mat_cumsum_LK.row(L_static-1));

/*(*T_fk_p)=(*T_fk_p)%Y_compensation_mat_FK;*/

}

#ifdef BOTTOM_UP_MEX_FLAG

void normalize_Y_lk_given_L_current(mat* Y_lk_p, mat* T_fk_p){

int l;

kronmat_real_LK.zeros();

/*Cumsum along dim=0*/
dummy_mat_cumsum_LK=cumsum((*Y_lk_p), 0);

/*kronmat_real_LK=kron(ones_col_Lx1, dummy_mat_cumsum_LK.row(L_static-1));*/	/*Check if ones_col_Lx1 is defined elsewhere already */

for (l=0; l<L_current; l++){

(*Y_lk_p).row(l)=(*Y_lk_p).row(l)/dummy_mat_cumsum_LK.row(L_static-1);

}

/*(*Y_lk_p)=(*Y_lk_p)/kronmat_real_LK;*/

/*Compensate*/
/*Y_compensation_mat_FK=kron(ones_col_Fx1, dummy_mat_cumsum_LK.row(L_static-1));*/

/*(*T_fk_p)=(*T_fk_p)%Y_compensation_mat_FK;*/

}

#endif

mat::fixed<N_static, K_static> compensation_mat_NK;

/*Arithmetic mean*/
void normalize_Y_lk_2(mat* Y_lk_p, mat* V_nk_p){

/*Cumsum along dim=0*/
dummy_mat_cumsum_LK=cumsum((*Y_lk_p), 0);

kronmat_real_LK=kron(ones_col_Lx1, dummy_mat_cumsum_LK.row(L_static-1));	/*Check if ones_col_Lx1 is defined elsewhere already */

(*Y_lk_p)=(*Y_lk_p)/kronmat_real_LK;

/*Compensate*/
compensation_mat_NK=kron(ones_col_Nx1, dummy_mat_cumsum_LK.row(L_static-1));

(*V_nk_p)=(*V_nk_p)%compensation_mat_NK;

}

mat::fixed<F_static, K_static> dummy_mat_cumsum_FK;
mat::fixed<F_static, K_static> kronmat_real_FK;

/*L2 norm*/
void normalize_T_fk(mat* T_fk_p, mat* V_nk_p){

/*First take the squares*/
dummy_mat_cumsum_FK=(*T_fk_p)%(*T_fk_p);

/*Cumsum along dim=0*/
dummy_mat_cumsum_FK=cumsum(dummy_mat_cumsum_FK, 0);

kronmat_real_FK=kron(ones_col_Fx1, dummy_mat_cumsum_FK.row(F_static-1));

kronmat_real_FK=sqrt(kronmat_real_FK);

(*T_fk_p)=(*T_fk_p)/kronmat_real_FK;

/*Compensate*/
compensation_mat_NK=kron(ones_col_Nx1, dummy_mat_cumsum_FK.row(F_static-1));

compensation_mat_NK=sqrt(compensation_mat_NK);

/*(*V_nk_p)=(*V_nk_p)%compensation_mat_NK;*/

}

/*Arithmetic mean*/
void normalize_T_fk_arth_mean(mat* T_fk_p, mat* V_nk_p){

/*Cumsum along dim=0*/
dummy_mat_cumsum_FK=cumsum((*T_fk_p), 0);

kronmat_real_FK=kron(ones_col_Fx1, dummy_mat_cumsum_FK.row(F_static-1));

(*T_fk_p)=(*T_fk_p)/kronmat_real_FK;

/*Compensate*/
compensation_mat_NK=kron(ones_col_Nx1, dummy_mat_cumsum_FK.row(F_static-1));

/*(*V_nk_p)=(*V_nk_p)%compensation_mat_NK;*/

}

/*Arithmetic mean*/
void normalize_T_fk_arth_mean_2(mat* T_fk_p, mat* Y_lk_p){

/*Cumsum along dim=0*/
dummy_mat_cumsum_FK=cumsum((*T_fk_p), 0);

kronmat_real_FK=kron(ones_col_Fx1, dummy_mat_cumsum_FK.row(F_static-1));

(*T_fk_p)=(*T_fk_p)/kronmat_real_FK;

/*Compensate*/
compensation_mat_LK=kron(ones_col_Lx1, dummy_mat_cumsum_FK.row(F_static-1));

/*(*Y_lk_p)=(*Y_lk_p)%compensation_mat_LK;*/

}

static double V_frob_avg;
static rowvec::fixed<K_static> V_frob_avgvec;
static rowvec::fixed<K_static> V_frobnormvec;	/*compute simply the frob norm*/
static rowvec::fixed<K_static> V_ratiovec;	/*want the avg to map to 1*/

static void normalize_Z_ok_dep_V(mat* V_nk_p){

/*first compute frob norm: sum(V, 0)*/
V_frobnormvec=sum(*V_nk_p, 0);

/*next compute the avg*/	
V_frob_avg=sum(V_frobnormvec)/((double)K_static);

V_frob_avgvec.fill(V_frob_avg);

/*ratio=frob/avg*/	
V_ratiovec=V_frobnormvec/V_frob_avgvec;

}

static mat::fixed<O_static, K_static> Z_ok_tempmat_local;
static rowvec::fixed<K_static> checkvec_local_1xK;
/*mat::fixed<F_static, K_static> kronmat_real_FK;*/

void normalize_Z_ok_compensate_TV(mat* Z_ol_p, mat* Y_lk_p, mat* T_fk_p, mat*  V_nk_p){

/*normalize_Z_ok_dep_V(V_nk_p);*/

Z_ok_tempmat_local=(*Z_ol_p)*(*Y_lk_p);

/*Modify Y_lk in terms of its K bins*/
(*Y_lk_p)=(*Y_lk_p)/kron(ones_col_Lx1, sqrt(sum(Z_ok_tempmat_local%Z_ok_tempmat_local, 0)) );
/*(*Y_lk_p)=(*Y_lk_p)%kron(ones_col_Lx1, V_ratiovec/(sqrt(sum(Z_ok_tempmat_local%Z_ok_tempmat_local, 0))) );*/

/*Can do a simple check to see if Z_ok is now normalized, if you verify that 

sqrt(sum(Z_ok_out%Z_ok_out, 0))

is in fact a 1xK row vector of 1's

*/

Z_ok_tempmat_local=(*Z_ol_p)*(*Y_lk_p);

checkvec_local_1xK=sqrt(sum(Z_ok_tempmat_local%Z_ok_tempmat_local, 0));

checkvec_local_1xK.print("normalize_Z_ok_compensate_T: checkvec_local_1xK:");

/*Compensate onto T_fk*/
(*T_fk_p)=(*T_fk_p)%kron(ones_col_Fx1, sqrt(sum(Z_ok_tempmat_local%Z_ok_tempmat_local, 0)));

/*Also compensate T_fk using V_ratiovec*/
/*(*T_fk_p)=(*T_fk_p)%kron(ones_col_Fx1, V_ratiovec);*/

}

/*Assumption normalize_Y_lk_compensate_Z: L_static=L_target*/
static colvec::fixed<L_static> frob_norm_colvec_norm_Y_comp_Z;

/*If L_static>L_target, you need to first sort and find the L_target most significant clusters*/

void normalize_Y_lk_compensate_Z_Lstatic_eq_Ltarget(mat* Y_lk_p, mat* Z_ol_p){

int l_iter;

/*compute each frobenius norm in Y_lk for l=1:L*/

frob_norm_colvec_norm_Y_comp_Z=sqrt(sum((*Y_lk_p)%(*Y_lk_p), 1));

/*Divide out each row by its frob norm:*/
(*Y_lk_p)=(*Y_lk_p)/kron(frob_norm_colvec_norm_Y_comp_Z, ones_row_1xK);

/*Compensate by scaling each col of Z_ol by the transpose of the column vector*/
/*(*Z_ol_p)=(*Z_ol_p)%kron(ones_col_Ox1, trans(frob_norm_colvec_norm_Y_comp_Z));*/

}

static rowvec::fixed<L_static> frob_norm_rowvec_norm_Z;

void normalize_Z_ol_frob_norm(mat* Z_ol_p){

(*Z_ol_p).print("Z_ol before:");

frob_norm_rowvec_norm_Z=sqrt(sum( (*Z_ol_p)%(*Z_ol_p) , 0));

frob_norm_rowvec_norm_Z.print("inside normalize_Z_ol_frob_norm, line 318:, normalize_Z_ol_frob_norm row vector:");

/*Divide out each col by its frob norm:*/
(*Z_ol_p)=(*Z_ol_p)/kron(ones_col_Ox1, frob_norm_rowvec_norm_Z);

(*Z_ol_p).print("Z_ol after:");

frob_norm_rowvec_norm_Z=sqrt(sum( (*Z_ol_p)%(*Z_ol_p) , 0));

frob_norm_rowvec_norm_Z.print("inside normalize_Z_ol_frob_norm, line 318:, normalize_Z_ol_frob_norm row vector, AFTER!!!!!!!:");

}

static rowvec::fixed<K_static> frob_norm_rowvec_norm_V;

void normalize_V_nk_frob_norm(mat* V_nk_p, mat* T_fk_p){

frob_norm_rowvec_norm_V=sqrt(sum( (*V_nk_p)%(*V_nk_p) , 0));

/*Divide out each col by its frob norm:*/
(*V_nk_p)=(*V_nk_p)/kron(ones_col_Nx1, frob_norm_rowvec_norm_V);

/*Compensate onto T_fk*/
(*T_fk_p)=(*T_fk_p)%kron(ones_col_Fx1, frob_norm_rowvec_norm_V);

}