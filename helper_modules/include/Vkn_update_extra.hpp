extern void check_V_scale_down_T(mat* T_fk_p, mat* V_nk_p);
extern void V_nk_erf_threshold(mat* V_nk_p);

extern void V_weighted_scaling(mat* V_nk_p); 
extern void V_weighted_scaling_over_classes(mat* V_nk_p, mat* Y_lk_p);