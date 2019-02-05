#include "local_inc.hpp"

mat::fixed<O_static, K_static> Zok_sparsity_common_sharedmat_Zok; 
rowvec::fixed<K_static> Q1_1xK;
rowvec::fixed<K_static> Q2_1xK;
mat::fixed<O_static, K_static> Q3_OxK;
rowvec::fixed<K_static> Q4_1xK;
rowvec::fixed<K_static> Q5_1xK;
rowvec::fixed<K_static> Q6_1xK;
mat::fixed<O_static, K_static> Q7_OxK;

static void compute_Zok_sparsity_common_sharedmat_Zok(mat* Z_ol_p, mat* Y_lk_p){

	Zok_sparsity_common_sharedmat_Zok=(*Z_ol_p)*(*Y_lk_p);

}

static void compute_Zok_sparsity_common_Quantity1(void){

Q1_1xK=(1/(sqrt((double)O_static)-1))*(ones_row_1xK/Q6_1xK);

}

static void compute_Zok_sparsity_common_Quantity2(void){

Q2_1xK=sqrt(Q6_1xK);

}

static void compute_Zok_sparsity_common_Quantity3(void){

Q3_OxK=Q7_OxK/Zok_sparsity_common_sharedmat_Zok;

}

/*This needs to be augmented with sqrt((double)O_static): Did it */
static void compute_Zok_sparsity_common_Quantity4(void){

Q4_1xK=Q5_1xK/Q2_1xK;

Q4_1xK=sqrt((double)O_static)*ones_row_1xK-Q4_1xK;

Q4_1xK=(1/(sqrt((double)O_static)-1))*Q4_1xK;

}

static void compute_Zok_sparsity_common_Quantity5(void){

Q5_1xK=sum(Q7_OxK, 0);

}

static void compute_Zok_sparsity_common_Quantity6(void){

Q6_1xK=sum(square(Zok_sparsity_common_sharedmat_Zok), 0);

}

static void compute_Zok_sparsity_common_Quantity7(void){

Q7_OxK=abs(Zok_sparsity_common_sharedmat_Zok);

}

void compute_Zok_sparsity_common_Quantities(mat* Z_ol_p, mat* Y_lk_p){

/*compute Z_ok*/
compute_Zok_sparsity_common_sharedmat_Zok(Z_ol_p, Y_lk_p);

/*compute abs(Z_ok): only depends on Z_ok*/
compute_Zok_sparsity_common_Quantity7();

/*compute the sum of squares of Z_ok: only depends on Z_ok */
compute_Zok_sparsity_common_Quantity6();

/*Depends on 7*/
compute_Zok_sparsity_common_Quantity3();

/*Depends on 6*/
compute_Zok_sparsity_common_Quantity2();

/*Depends on 7*/
compute_Zok_sparsity_common_Quantity5();

/*Depends on 2, 5*/
compute_Zok_sparsity_common_Quantity4();

/*Depends on 6*/
compute_Zok_sparsity_common_Quantity1();

}

colvec::fixed<K_static> sigma_vec;
double Zok_sparsity_lamda;

void populate_sigma_vec(double fill_value){

Zok_sparsity_lamda=pow(10, 10);

sigma_vec.fill(fill_value);

}