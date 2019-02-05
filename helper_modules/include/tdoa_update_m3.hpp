extern void TDOA_update_module3_entry_point1(cx_cube* Xtilde_fnm_p, mat* Z_ol_p, mat* Y_lk_p, mat* T_fk_p, mat* V_nk_p, cube* W_fom_p, cx_cube* expj_Phi_W_fom_p, cx_cube* expj_Phi_S_nkf_p);

extern void populate_6x1_possible_permutations(void);

#define NUM_PERMUTATIONS 6

extern cube::fixed<L_static, L_static, NUM_PERMUTATIONS> permutation_list;

