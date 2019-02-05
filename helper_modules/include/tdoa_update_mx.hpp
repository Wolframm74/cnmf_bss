
#define F_SPATIAL_ALIASING_MAX_INDEX 218 /*1701Hz / 8000 Hz =~ f_index=218 */

#define F_static_admissible (F_SPATIAL_ALIASING_MAX_INDEX-1)

#define E_static L_static*F_static_admissible


#define WLEN_CONSTANT 4*512
#define HLEN_CONSTANT WLEN_CONSTANT/4
#define NFFT_CONSTANT WLEN_CONSTANT
#define FS_CONSTANT 16000

#define N_STFT_CONSTANT NFFT_CONSTANT

#define J_VALUE 0 

extern double compute_angle_tdoa_update(cx_double input);