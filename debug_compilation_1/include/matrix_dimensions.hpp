/*#define NUM_WORKER_THREADS 4*/

#define MEX_UNIT_TESTING
/*#define REGULAR_USE*/

#ifdef MEX_UNIT_TESTING

	#define N_ITERATIONS 201

	/*Dimensions for static matrix allocation*/

/*	#define M_static 2
	#define L_static 3
	#define F_static 20
	#define N_static 10
	#define K_static 8
	#define O_static 12*/

	#define M_static 2
	#define L_static 30
	#define F_static 513
	#define N_static 2128		/*Allocate an array of "1052364096" doubles in order to store the Phi_S_Struct data type. As a double is approx 8 bytes, this amounts to approx 8 gigs.  */
	/*#define N_static 393*/
	#define K_static 60
	#define O_static 110

#endif

#ifdef REGULAR_USE

	#define N_ITERATIONS 201

	/*Dimensions for static matrix allocation*/
	#define M_static 2
	#define L_static 30
	#define F_static 513
	/*#define N_static 393*/
	#define N_static 2128
	#define K_static 60
	#define O_static 110	/*This is from resolutions of 36 azimuth directions over 360 degrees and 8 elevation directions over 180 degrees. */

#endif