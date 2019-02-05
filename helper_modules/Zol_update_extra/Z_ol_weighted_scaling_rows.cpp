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

rowvec::fixed<L_static> oth_rowvec_dummy;
rowvec::fixed<L_static> oth_rowvec_of_Z_ol;



void Z_scale_oth_row(void){

double min_value;
double max_value;
double max_value_adj;
double dummy_value;
int l_index; 
int i_iter;

min_value=min(oth_rowvec_of_Z_ol);

max_value=max(oth_rowvec_of_Z_ol);

//the cu rrent max element of kth_colvec_dummy should be equal to the max_value_adj, calculated below. 
oth_rowvec_dummy=oth_rowvec_of_Z_ol-(min_value)*oth_rowvec_dummy.ones();

max_value_adj=max_value-min_value;

for (l_index=0; l_index<L_static; l_index++){

dummy_value=oth_rowvec_dummy(l_index);



//Since max_value_adj represents the max() of kth_colvec_dummy, begin with it as a first case, and proceed logically to cover the other cases, parsing the cases in descending order (in terms of dummy value's magnitude) 
if(dummy_value>=max_value_adj){

	oth_rowvec_of_Z_ol(l_index)=dummy_value*1.01;

	}
else {


	for (i_iter=0; i_iter<10; i_iter++){

	if( (dummy_value<(10-i_iter)*0.1*max_value_adj) && (dummy_value>=((10-(i_iter+1))*0.1*max_value_adj)) ){

	oth_rowvec_of_Z_ol(l_index)=dummy_value*(1.01-exp(-0.5*(10-i_iter)));

	} /*end-if*/

	} /*end-for*/

}	/*end-else*/

}

}


/*entry point: required arg: mat* Y_lk_p, transpose it*/
void Z_weighted_scaling(mat* Z_ol_p){

int o_index;

for (o_index=0; o_index<O_static; o_index++){

oth_rowvec_of_Z_ol=(*Z_ol_p).row(o_index);

Z_scale_oth_row();

(*Z_ol_p).row(o_index)=oth_rowvec_of_Z_ol;

}


}
