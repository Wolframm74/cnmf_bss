#define PVALUE_FIXED 1.2

extern double accu_Zok_square;

extern mat::fixed<O_static, K_static> Zok_shared_mat;

extern mat::fixed<O_static, L_static> Zok_sparsity_outmat_OL_den;
extern mat::fixed<L_static, K_static> Zok_sparsity_outmat_LK_den;

extern void compute_Zok(mat* Z_ol_p, mat* Y_lk_p);

extern void Zok_Zol_update_extra_sparsity(mat* Z_ol_p, mat* Y_lk_p);

extern void Zok_Ylk_update_extra_sparsity(mat* Z_ol_p, mat* Y_lk_p);