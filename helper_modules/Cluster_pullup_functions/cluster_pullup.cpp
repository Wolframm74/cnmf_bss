#include "local_inc.hpp"

/*Pass in a vector sequence: boost the elements below the max l1 norm in the sequence, using a generalized gaussian */

/*Assumed that the caller passes in a vector of L1 norms*/

colvec colvec_dynamic_02;
colvec colvec_dynamic_02_scaling_vector;
colvec colvec_dynamic_02_ones;

void pull_up_vector_elements_generalized_gaussian(double beta_value, double mu_value, int Vec_Size, double boost_value, double offset_value){

double mu, alpha, alpha_value;

double max_element;

max_element=max(colvec_dynamic_02);

mu=mu_value*max_element;

colvec_dynamic_02_ones=colvec_dynamic_02;

colvec_dynamic_02_ones.ones();

alpha_value=(1-mu_value)/2;

mu=mu_value*max_element;

/*colvec_dynamic_02_ones=mu*colvec_dynamic_02_ones;*/

alpha=alpha_value*max_element;

colvec_dynamic_02_scaling_vector=colvec_dynamic_02-colvec_dynamic_02_ones;

colvec_dynamic_02_scaling_vector=boost_value*exp(  -pow((abs(colvec_dynamic_02-mu*colvec_dynamic_02_ones))/alpha, beta_value)  )+offset_value*colvec_dynamic_02_ones;

colvec_dynamic_02=colvec_dynamic_02_scaling_vector%colvec_dynamic_02;

}

/*Would like to allocate a global object here*/
colvec colvec_dynamic_01;
colvec colvec_dynamic_01_ones;

static colvec partial_derivative_vec;

static double mean_value;

/*Pass in a vector sequence: pull the elements towards the mean*/
void pull_vector_elements_towards_mean(int Vec_Size, double aeta_value){

mean_value=(as_scalar(sum(colvec_dynamic_01)))/((double)Vec_Size);

colvec_dynamic_01_ones=colvec_dynamic_01;

colvec_dynamic_01_ones.ones();

colvec_dynamic_01_ones=mean_value*colvec_dynamic_01_ones;

/*Would like for the global object to dynamically change its size depending on the size of the intended vector to be operated on*/

/*Function should do same simple elementwise operations regardless of vector size. */	
 
partial_derivative_vec=2*(colvec_dynamic_01-colvec_dynamic_01_ones);

colvec_dynamic_01=colvec_dynamic_01%exp((aeta_value)*partial_derivative_vec);

}

/*Thoughts:

Confusion always stems back to the fact that you ingorant to some of the details that occur when an assignment (=) allows one to copy the data from one matrx / auxiliary memory to another. 

*/

