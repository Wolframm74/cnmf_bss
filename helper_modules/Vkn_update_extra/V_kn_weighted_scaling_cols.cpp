#include "local_inc.hpp"

/*1st function: for col=0:N_static-1; 
find the max and min for all k=0:K_static-1.
compute the difference
subtract off the min from all elements

start going through swtich/case statements

check each "min-adjusted element" and compare it against the max/min-adjusted element.

funnel each min adjusted element into a case where depending on its ratio compared to the max min adjusted element it will be scaled down in value.

Scaled down more severely for small valued elements. 

Large valued elements comparable to the max valued element should be assigned a gain closer to 1. 

Return the column vector of scaled values.

Entry point function can set it to the corresponding col of V_kn. 

*/

colvec::fixed<K_static> nth_colvec_dummy;
colvec::fixed<K_static> nth_colvec_of_V_kn;



void V_scale_nth_column(void){

double min_value;
double max_value;
double max_value_adj;
double dummy_value;
int k_index; 
int i_iter;

min_value=min(nth_colvec_of_V_kn);

max_value=max(nth_colvec_of_V_kn);

//the cu rrent max element of nth_colvec_dummy should be equal to the max_value_adj, calculated below. 
nth_colvec_dummy=nth_colvec_of_V_kn-(min_value)*nth_colvec_dummy.ones();

max_value_adj=max_value-min_value;

for (k_index=0; k_index<K_static; k_index++){

dummy_value=nth_colvec_dummy(k_index);



//Since max_value_adj represents the max() of nth_colvec_dummy, begin with it as a first case, and proceed logically to cover the other cases, parsing the cases in descending order (in terms of dummy value's magnitude) 
if(dummy_value>=max_value_adj){

	nth_colvec_of_V_kn(k_index)=dummy_value*1.01;

	}
else {


	for (i_iter=0; i_iter<10; i_iter++){

	if( (dummy_value<(10-i_iter)*0.1*max_value_adj) && (dummy_value>=((10-(i_iter+1))*0.1*max_value_adj)) ){

	nth_colvec_of_V_kn(k_index)=dummy_value*(1.01-exp(-0.5*(10-i_iter)));

	} /*end-if*/

	} /*end-for*/

}	/*end-else*/

}

}
/*void V_scale_nth_column(void){

double min_value;
double max_value;
double max_value_adj;
double dummy_value;
int k_index; 


min_value=min(nth_colvec_of_V_kn);

max_value=max(nth_colvec_of_V_kn);

//the cu rrent max element of nth_colvec_dummy should be equal to the max_value_adj, calculated below. 
nth_colvec_dummy=nth_colvec_of_V_kn-(min_value)*nth_colvec_dummy.ones();

max_value_adj=max_value-min_value;

for (k_index=0; k_index<K_static; k_index++){

dummy_value=nth_colvec_dummy(k_index);



//Since max_value_adj represents the max() of nth_colvec_dummy, begin with it as a first case, and proceed logically to cover the other cases, parsing the cases in descending order (in terms of dummy value's magnitude) 
if(dummy_value>=max_value_adj){

	nth_colvec_of_V_kn(k_index)=dummy_value;

	}
else if(dummy_value<max_value_adj && dummy_value>=(0.9*max_value_adj)){

	nth_colvec_of_V_kn(k_index)=dummy_value*0.9;

	}
else if(dummy_value<(0.9*max_value_adj) && dummy_value>=(0.8*max_value_adj)){

	nth_colvec_of_V_kn(k_index)=dummy_value*0.8;

	}
else if(dummy_value<(0.8*max_value_adj) && dummy_value>=(0.7*max_value_adj)){

	nth_colvec_of_V_kn(k_index)=dummy_value*0.7;

	}
else if(dummy_value<(0.7*max_value_adj) && dummy_value>=(0.6*max_value_adj)){

	nth_colvec_of_V_kn(k_index)=dummy_value*0.6;

	}
else if(dummy_value<(0.6*max_value_adj) && dummy_value>=(0.5*max_value_adj)){

	nth_colvec_of_V_kn(k_index)=dummy_value*0.5;

	}
else if(dummy_value<(0.5*max_value_adj) && dummy_value>=(0.4*max_value_adj)){

	nth_colvec_of_V_kn(k_index)=dummy_value*0.4;

	}
else if(dummy_value<(0.4*max_value_adj) && dummy_value>=(0.3*max_value_adj)){

	nth_colvec_of_V_kn(k_index)=dummy_value*0.3;

	}		
else if(dummy_value<(0.3*max_value_adj) && dummy_value>=(0.2*max_value_adj)){

	nth_colvec_of_V_kn(k_index)=dummy_value*0.2;

	}
else if(dummy_value<(0.2*max_value_adj) && dummy_value>=(0.1*max_value_adj)){

	nth_colvec_of_V_kn(k_index)=dummy_value*0.1;

	}
else if(dummy_value<(0.1*max_value_adj) && dummy_value>=(0*max_value_adj)){

	nth_colvec_of_V_kn(k_index)=dummy_value*0.05;

	}


}

}*/


/*entry point: required arg: mat* V_nk_p, transpose it*/
void V_weighted_scaling(mat* V_nk_p){

int n_index;

for (n_index=0; n_index<N_static; n_index++){

nth_colvec_of_V_kn=trans((*V_nk_p).row(n_index));

V_scale_nth_column();

(*V_nk_p).row(n_index)=trans(nth_colvec_of_V_kn);

}

}


/*------------------------------------------------------------------------------------------------------------------------*/

static colvec::fixed<L_static> nth_colvec_of_V_kn_2;


static mat::fixed<L_static, K_static> Y_lk_spreadmat; 
static rowvec::fixed<K_static> output_rowvec_local_1xK;

void V_scale_nth_column_over_classes(mat* Y_lk_p){

double max_value;
int l_index; 
int i_iter;

max_value=max(nth_colvec_of_V_kn_2);

for (l_index=0; l_index<L_static; l_index++){


nth_colvec_of_V_kn_2(l_index)=1.05-exp(-3*((nth_colvec_of_V_kn_2(l_index))/max_value));


}

/*nth_colvec_of_V_kn_2.print("nth_colvec_of_V_kn_2: line 206:");*/

/*spread the info over Y_lk using kronecker product*/
Y_lk_spreadmat=kron(nth_colvec_of_V_kn_2, ones_row_1xK)%(*Y_lk_p);

output_rowvec_local_1xK=sum(Y_lk_spreadmat, 0);

}



/*entry point: required arg: mat* V_nk_p, transpose it*/
void V_weighted_scaling_over_classes(mat* V_nk_p, mat* Y_lk_p){

int n_index, l_index;

for (n_index=0; n_index<N_static; n_index++){

/*Partition the K columns into their L classes via Y_lk*/

/*Search for the most dominant class l*/
for (l_index=0; l_index<L_static; l_index++){

nth_colvec_of_V_kn_2(l_index)=as_scalar(sum( ((*V_nk_p).row(n_index))%(*Y_lk_p).row(l_index) ));

}

/*The most dominant class l should persist its value where as the other two classes need to be attenuated relative to their respective amplitudes*/
/*Collapse a Kx1 into an Lx1 using the info provided by Y_lk*/
/*assign it an Lx1 column vector*/

V_scale_nth_column_over_classes(Y_lk_p);

(*V_nk_p).row(n_index)=(*V_nk_p).row(n_index)%output_rowvec_local_1xK;

}

}