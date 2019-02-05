global HALF_WINDOW_LEN
global WINDOW_LEN (HALF_WINDOW_LEN*2)+1

HALF_WINDOW_LEN=5
WINDOW_LEN=(HALF_WINDOW_LEN*2)+1

storage_vec=zeros(1, WINDOW_LEN);	%storage for the row vector to be dotted with the window
				%Indexing the correct "center" sample is key here. 


window_vec=zeros(1, WINDOW_LEN);	%odd numbered filter

outvec_1xN=zeros(1, N_static);

%Indices that index the vector to copy from (V_kn for some row index k) 

%Indices tha index the vector to be copied to (

static void V_filter_dot_op(int n_index){

outvec_1xN(n_index)=dot(window_vec, storage_vec);

}

static void V_filter_row(){

int n_index;

%for index=0:N_static-1
for (n_index=0; n_index<N_static; n_index++){

	storage_vec()=

	%window_vec()=

	V_filter_dot_op(n_index);


}

}

void V_nk_MA_filter_rows(){

}