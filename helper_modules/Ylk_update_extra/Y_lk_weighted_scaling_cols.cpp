#include "local_inc.hpp"

/*1st function: for col=0:K_static-1; 
find the max and min for all k=0:L_static-1.
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

colvec::fixed<L_static> kth_colvec_dummy;
colvec::fixed<L_static> kth_colvec_of_Y_lk;



void Y_scale_kth_column(void){

double min_value;
double max_value;
double max_value_adj;
double dummy_value;
int l_index; 
int i_iter;

min_value=min(kth_colvec_of_Y_lk);

max_value=max(kth_colvec_of_Y_lk);

//the cu rrent max element of kth_colvec_dummy should be equal to the max_value_adj, calculated below. 
kth_colvec_dummy=kth_colvec_of_Y_lk-(min_value)*kth_colvec_dummy.ones();

max_value_adj=max_value-min_value;

for (l_index=0; l_index<L_static; l_index++){

dummy_value=kth_colvec_dummy(l_index);



//Since max_value_adj represents the max() of kth_colvec_dummy, begin with it as a first case, and proceed logically to cover the other cases, parsing the cases in descending order (in terms of dummy value's magnitude) 
if(dummy_value>=max_value_adj){

	kth_colvec_of_Y_lk(l_index)=dummy_value*1.01;

	}
else {


	for (i_iter=0; i_iter<10; i_iter++){

	if( (dummy_value<(10-i_iter)*0.1*max_value_adj) && (dummy_value>=((10-(i_iter+1))*0.1*max_value_adj)) ){

	kth_colvec_of_Y_lk(l_index)=dummy_value*(1.01-exp(-0.5*(10-i_iter)));

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
int l_index; 


min_value=min(kth_colvec_of_Y_lk);

max_value=max(kth_colvec_of_Y_lk);

//the cu rrent max element of kth_colvec_dummy should be equal to the max_value_adj, calculated below. 
kth_colvec_dummy=kth_colvec_of_Y_lk-(min_value)*kth_colvec_dummy.ones();

max_value_adj=max_value-min_value;

for (l_index=0; l_index<L_static; l_index++){

dummy_value=kth_colvec_dummy(l_index);



//Since max_value_adj represents the max() of kth_colvec_dummy, begin with it as a first case, and proceed logically to cover the other cases, parsing the cases in descending order (in terms of dummy value's magnitude) 
if(dummy_value>=max_value_adj){

	kth_colvec_of_Y_lk(l_index)=dummy_value;

	}
else if(dummy_value<max_value_adj && dummy_value>=(0.9*max_value_adj)){

	kth_colvec_of_Y_lk(l_index)=dummy_value*0.9;

	}
else if(dummy_value<(0.9*max_value_adj) && dummy_value>=(0.8*max_value_adj)){

	kth_colvec_of_Y_lk(l_index)=dummy_value*0.8;

	}
else if(dummy_value<(0.8*max_value_adj) && dummy_value>=(0.7*max_value_adj)){

	kth_colvec_of_Y_lk(l_index)=dummy_value*0.7;

	}
else if(dummy_value<(0.7*max_value_adj) && dummy_value>=(0.6*max_value_adj)){

	kth_colvec_of_Y_lk(l_index)=dummy_value*0.6;

	}
else if(dummy_value<(0.6*max_value_adj) && dummy_value>=(0.5*max_value_adj)){

	kth_colvec_of_Y_lk(l_index)=dummy_value*0.5;

	}
else if(dummy_value<(0.5*max_value_adj) && dummy_value>=(0.4*max_value_adj)){

	kth_colvec_of_Y_lk(l_index)=dummy_value*0.4;

	}
else if(dummy_value<(0.4*max_value_adj) && dummy_value>=(0.3*max_value_adj)){

	kth_colvec_of_Y_lk(l_index)=dummy_value*0.3;

	}		
else if(dummy_value<(0.3*max_value_adj) && dummy_value>=(0.2*max_value_adj)){

	kth_colvec_of_Y_lk(l_index)=dummy_value*0.2;

	}
else if(dummy_value<(0.2*max_value_adj) && dummy_value>=(0.1*max_value_adj)){

	kth_colvec_of_Y_lk(l_index)=dummy_value*0.1;

	}
else if(dummy_value<(0.1*max_value_adj) && dummy_value>=(0*max_value_adj)){

	kth_colvec_of_Y_lk(l_index)=dummy_value*0.05;

	}


}

}*/


/*entry point: required arg: mat* Y_lk_p, transpose it*/
void Y_weighted_scaling(mat* Y_lk_p){

int k_index;

for (k_index=0; k_index<K_static; k_index++){

kth_colvec_of_Y_lk=(*Y_lk_p).col(k_index);

Y_scale_kth_column();

(*Y_lk_p).col(k_index)=kth_colvec_of_Y_lk;

}


}

