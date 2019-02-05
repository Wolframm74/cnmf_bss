/*Zok_sparsity_common.cpp*/

extern mat::fixed<O_static, L_static> Zok_sparsity_common_sharedmat_Zok; 
extern rowvec::fixed<K_static> Q1_1xK;
extern rowvec::fixed<K_static> Q2_1xK;
extern mat::fixed<O_static, K_static> Q3_OxK;
extern rowvec::fixed<K_static> Q4_1xK;
extern rowvec::fixed<K_static> Q5_1xK;
extern rowvec::fixed<K_static> Q6_1xK;
extern mat::fixed<O_static, K_static> Q7_OxK;

extern void compute_Zok_sparsity_common_Quantities(mat* Z_ol_p, mat* Y_lk_p);

extern colvec::fixed<K_static> sigma_vec;
extern double Zok_sparsity_lamda;

extern void populate_sigma_vec(double fill_value);

/*Zok_sparsity_Ylk_update.cpp*/
extern mat::fixed<L_static, K_static> pve_part_Partial_wrt_Ylk;
extern mat::fixed<L_static, K_static> nve_part_Partial_wrt_Ylk;
extern void Zok_sparsity_Ylk_update(mat* Z_ol_p);

/*Zok_sparsity_Zol_update.cpp*/
extern mat::fixed<O_static, L_static> pve_part_Partial_wrt_Zol;
extern mat::fixed<O_static, L_static> nve_part_Partial_wrt_Zol;
extern void Zok_sparsity_Zol_update(mat* Y_lk_p);


