extern void populate_expj_Phi_S_nkf(arg_struct_t* argStruct_p);
extern void V_nk_update_Phi_S_wrapper(mxArray *plhs[], arg_struct_t* argStruct_p);

extern void Vnk_update_Phi_S(mat* V_nk_p, cx_cube* expj_Phi_S_fkn_p, arg_struct_t* argStruct_p);