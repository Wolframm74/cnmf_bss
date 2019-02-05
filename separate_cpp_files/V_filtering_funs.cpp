
/*[0]
[1]
[2]
[3]
[4]
[5] window center
[6]
[7]
[8]
[9]
[10]*/


bool overlap1;
bool overlap2;
bool overlap3;


#define HALF_WINDOW_LEN 5
#define WINDOW_LEN (HALF_WINDOW_LEN*2)+1

rowvec::fixed<WINDOW_LEN> storage_vec;	/*storage for the row vector to be dotted with the window*/
				/*Indexing the correct "center" sample is key here. */


rowvec::fixed<WINDOW_LEN> window_vec;	/*odd numbered filter*/

rowvec::fixed<WINDOW_LEN> window_dummy;


rowvec::fixed<N_static> outvec_1xN;

/*Indices that index the vector to copy from (V_kn for some row index k) */
copy_from_index1;
copy_from_index2;

/*Indices tha index the vector to be copied to (*/
copy_to_index1;
copy_to_index2;	

static void V_filter_dot_op(int n_index){

outvec_1xN(n_index)=dot(window_dummy, storage_vec);

}

static void V_filter_row(rowvec& V_kn_kth_row){

int n_index;

window_dummy.zeros();
storage_vec.zeros();

/*for index=0:N_static-1*/
for (n_index=0; n_index<N_static; n_index++){

	if (n_index<HALF_WINDOW_LEN){

		storage_index1=0;
		storage_index2=(WINDOW_LEN-1)-(HALF_WINDOW_LEN-n_index);
		window_index1=(HALF_WINDOW_LEN-n_index);
		window_index2=WINDOW_LEN-1;

		storage_vec(storage_index1:storage_index2)=V_kn_kth_row(storage_index1:storage_index2);
		window_dummy(storage_index1:storage_index2)=window_vec(window_index1:window_index2);

	} else if (n_index<(N_static-HALF_WINDOW_LEN)) {

		storage_index1=n_index-HALF_WINDOW_LEN;
		storage_index2=n_index+HALF_WINDOW_LEN;
		window_index1=0;
		window_index2=WINDOW_LEN-1;

		storage_vec(window_index1:window_index2)=V_kn_kth_row(storage_index1:storage_index2);
		window_dummy=window_vec;


	} else {	/*assume n_index >= (N_static-WINDOW_LEN) */

		window_dummy.zeros();
		storage_vec.zeros();

		storage_index1=n_index-HALF_WINDOW_LEN;
		storage_index2=N_static-1;
		window_index1=0;
		window_index2=(N_static-1)-storage_index1;	

		storage_vec(window_index1:window_index2)=V_kn_kth_row(window_index1:window_index2);
		window_dummy(window_index1:window_index2)=window_vec(storage_index1:storage_index2);

	}

	V_filter_dot_op(n_index);


}

}

void V_nk_MA_filter_rows(*V_nk_p){

int k_index;

rowvec::fixed<N_static> V_kn_kth_row;

for (k_index=0; k_index<K_static; k_index++){

	V_kn_kth_row=transpose((*V_nk_p).col(k_index));

	V_filter_row(V_kn_kth_row);

}

}


/*: Expected output:

n_index=0

	overlap1=true;

storage_index1=0;
storage_index2=10-5=5;
window_index1=5;
window_index2=10;	

window_dummy(0:5)=window_vec(5:10);

storage_vec(0:5)=V_kn_nth_row(0:5)

n_index=1

	overlap1=true;

storage_index1=0;
storage_index2=10-4=6;
window_index1=4;
window_index2=10;	

window_dummy(0:6)=window_vec(4:10);

storage_vec(0:6)=V_kn_nth_row(0:6)

n_index=4

	overlap1=true;

storage_index1=0;
storage_index2=9;
window_index1=1;
window_index2=10;	

window_dummy(0:9)=window_vec(1:10);

storage_vec(0:9)=V_kn_nth_row(0:9)

n_index=5

	overlap2=true;

storage_index1=0;
storage_index2=10;
window_index1=0;
window_index2=10;	

window_dummy(0:10)=window_vec(0:10);

storage_vec(0:10)=V_kn_nth_row(0:10)

n_index=6

	overlap2=true;

storage_index1=1;
storage_index2=11;
window_index1=0;
window_index2=10;	

window_dummy(0:10)=window_vec(0:10);

storage_vec(0:10)=V_kn_nth_row(1:11)

n_index=778

	overlap2=true;

storage_index1=773;
storage_index2=783;
window_index1=0;
window_index2=10;	

window_dummy(0:10)=window_vec(0:10);
storage_vec(0:10)=V_kn_nth_row(773:783)

n_index=779

	overlap3=true;

storage_index1=774;
storage_index2=783;
window_index1=0;
window_index2=9;	

window_dummy.zeros(), storage_vec.zeros();
window_dummy(0:9)=window_vec(0:9);
storage_vec(0:9)=V_kn_nth_row(774:783)

n_index=783

	overlap3=true;

storage_index1=777;
storage_index2=783;
window_index1=0;
window_index2=6;	

window_dummy.zeros(), storage_vec.zeros();
window_dummy(0:6)=window_vec(0:6);
storage_vec(0:6)=V_kn_nth_row(777:783)

n_index=783

	overlap3=true;

storage_index1=778;
storage_index2=783;
window_index1=0;
window_index2=5;	

window_dummy.zeros(), storage_vec.zeros();
window_dummy(0:5)=window_vec(0:5);
storage_vec(0:5)=V_kn_nth_row(778:783)

*/