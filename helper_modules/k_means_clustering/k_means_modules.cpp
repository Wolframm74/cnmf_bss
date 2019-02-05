#include "local_inc.hpp"

#define L_target 3

#define K_MEANS_N_ITERATIONS 5

colvec::fixed<L_target> indices_vec;
static mat::fixed<O_static, K_static> Z_ok_local;


rowvec::fixed<L_static> collapse_Z_ol_summed_rowvec;

/*1. Search for the L_target most significant clusters. Save the indices of these clusters. */
void search_L_target_significant_clusters(mat* Z_ol_p, mat* Y_lk_p){

int l_index;

uword max_index;

collapse_Z_ol_summed_rowvec=sum(*Z_ol_p, 0);

for (l_index=0; l_index<L_target; l_index++){

collapse_Z_ol_summed_rowvec.max(max_index);

collapse_Z_ol_summed_rowvec(max_index)=0;

/*Save the index*/
indices_vec(l_index)=max_index;

}

/*indices_vec should be populated with L_target entries by here.*/

/*And compute Z_ok while you're here. */

Z_ok_local=(*Z_ol_p)*(*Y_lk_p);

Z_ok_local.print("search_L_target_significant_clusters: Z_ok_local:");

}

static colvec::fixed<O_static> Ox1_colvec_quantity1_local;
static colvec::fixed<O_static> Ox1_colvec_quantity2_local;

static rowvec::fixed<K_static> frob_norm_Y_rowvec_1xK;

/*2. Compute Y_lk */
void kmeans_compute_Y_lk(mat* Y_lk_p, mat* Z_ol_p){

double den_value;

int l_index, k_index, r_index;

int l_index_value, r_index_value;

/*collapse Y into a rowvec*/
frob_norm_Y_rowvec_1xK=sqrt(sum((*Y_lk_p)%(*Y_lk_p), 0));

for (l_index=0; l_index<L_target; l_index++){

l_index_value=indices_vec(l_index);

	for (k_index=0; k_index<K_static; k_index++){

		den_value=0;

		for (r_index=0; r_index<L_target; r_index++){

			r_index_value=indices_vec(r_index);

			Ox1_colvec_quantity1_local=Z_ok_local.col(k_index)-(*Z_ol_p).col(l_index_value);

			Ox1_colvec_quantity2_local=Z_ok_local.col(k_index)-(*Z_ol_p).col(r_index_value);

			den_value=den_value+ pow(( sqrt( as_scalar(trace(Ox1_colvec_quantity1_local*trans(Ox1_colvec_quantity1_local))) )/ sqrt( as_scalar( trace(Ox1_colvec_quantity2_local*trans(Ox1_colvec_quantity2_local)) ) ) ), 2);

		}

		(*Y_lk_p)(l_index_value,k_index)=1/den_value;
		/*(*Y_lk_p)(l_index_value,k_index)=frob_norm_Y_rowvec_1xK(k_index)/den_value;*/

	}

}


}

/*3. Compute Z_ol */
void kmeans_compute_Z_ol(mat* Z_ol_p, mat* Y_lk_p){

int l_index;

int l_index_value;

for (l_index=0; l_index<L_target; l_index++){

l_index_value=indices_vec(l_index);

(*Z_ol_p).col(l_index_value)=sum( kron( ones_col_Ox1 , square((*Y_lk_p).row(l_index_value)) )%Z_ok_local , 1)/as_scalar(sum( square((*Y_lk_p).row(l_index_value)) ));

}	

}



void k_means_clustering_entry(mat* Z_ol_p, mat* Y_lk_p){

int n_iter_kmeans;

search_L_target_significant_clusters(Z_ol_p, Y_lk_p);

/*Do roughly 5 iterations then exit.*/

for (n_iter_kmeans=0; n_iter_kmeans<K_MEANS_N_ITERATIONS; n_iter_kmeans++){

kmeans_compute_Y_lk(Y_lk_p, Z_ol_p);

kmeans_compute_Z_ol(Z_ol_p, Y_lk_p);

}

}