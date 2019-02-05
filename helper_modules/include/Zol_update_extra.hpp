/*extern void check_Z_pull_up_weak_clusters(mat* Z_ol_p);*/

extern void Z_pull_vector_elements_towards_mean_wrapper(mat* Z_ol_p);

extern void Z_pull_up_vector_elements_generalized_gaussian_wrapper(mat* Z_ol_p);

extern void check_Z_scale_down_Y_low_threshold(mat* Z_ol_p, mat* Y_lk_p);
extern void check_Z_scale_down_Y_high_threshold(mat* Z_ol_p, mat* Y_lk_p);

/*Z_ol_weighted_scaling_rows.cpp*/
extern void Z_weighted_scaling(mat* Z_ol_p);