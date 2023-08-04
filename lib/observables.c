/**
* \class observables
* \brief library of functions for the calculation of several observables
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <my_malloc.h>
#include <observables.h>
#include <io.h>
#include <time.h>
#include <geometry.h>
#include <alignment.h>
#include <time.h>
#include <geometry.h>

void compute_coupling_matrix(double *coupling_mat, traj *Trajectory, int fr_id, float sigma){
    /** 
    * routine that computes the coupling matrix over a frame
    * 
    * `coupling_mat` : coupling matrix
    * 
    * `Trajectory` : traj object
    * 
    * `fr_id` : frame ID
    * 
    * `sigma` : sigma parameter
    */ 
    double gauss_factor = 4*pow(sigma,2);
    int i,j;
    int knt; // check that is equal to pairs at the end!
    double dist_sq;
    knt = 0;
    for (i=0; i < Trajectory->n_at-1; i++){
        for (j=i+1; j < Trajectory->n_at; j++){
            dist_sq = pow(Trajectory->traj_coords[fr_id][3*i]-Trajectory->traj_coords[fr_id][3*j],2) + pow(Trajectory->traj_coords[fr_id][3*i+1]-Trajectory->traj_coords[fr_id][3*j+1],2) + pow(Trajectory->traj_coords[fr_id][3*i+2]-Trajectory->traj_coords[fr_id][3*j+2],2);
            coupling_mat[knt] = exp(-dist_sq/gauss_factor);
            knt++;
        }
    }
}

double compute_atomistic_coord_number(double *coupling_mat,traj *Trajectory,FILE *f_out_l){
    /** 
    * routine that computes the atomistic coordination number for a certain coupling matrix. 
    * Double counting is necessary to ensure proper normalisation to norm and scalar product
    * 
    * `coupling_mat` : coupling matrix
    * 
    * `Trajectory` : traj object
    * 
    * `f_out_l` : output filename
    */ 
    int i;
    double n_coord_at = Trajectory->n_at; // diagonal terms
    for (i=0; i < Trajectory->n_at * (Trajectory->n_at - 1) / 2; i++){
      n_coord_at += coupling_mat[i]*2;
    }
    fprintf(f_out_l,"not normalized n_coord_at: %lf\n", n_coord_at);
    n_coord_at = n_coord_at/(1.0*Trajectory->n_at);
    fprintf(f_out_l,"normalized n_coord_at: %lf\n", n_coord_at);
    return n_coord_at;
}

void compute_norm(cg_mapping *mapping, double *coupling_mat, double n_coord_at, int fr_id, FILE *f_out_l){
  /** 
  * routine that computes the norm of a mapping over a frame of a trajectory
  *
  * Parameters
  * ----------
  *
  * `mapping` : cg_mapping object
  * 
  * `coupling_mat` : coupling matrix
  * 
  * `n_coord_at` : atomistic coordination number
  * 
  * `fr_id` : frame index
  * 
  * `f_out_l` : output filename
  */ 
  int i,j,knt;
  double norm = 0.0;
  knt = 0;
  for (i=0; i < mapping->n_at; i++){
    norm += mapping->mapping[i]; // self-interaction term
    for (j=i+1; j < mapping->n_at; j++){
      norm += 2*coupling_mat[knt]*mapping->mapping[i]*mapping->mapping[j];
      knt +=1;
    }
  }
  fprintf(f_out_l, "norm[%d] = %lf\n",fr_id, norm);
  norm = norm/n_coord_at;
  fprintf(f_out_l,"normalized_norm[%d] = %lf\n",fr_id, norm);
  mapping->norms[fr_id] = norm;
}

double compute_distance(cg_mapping *mapping, cg_mapping *mapping_prime, double *coupling_mat, double n_coord_at, int fr_id, FILE *f_out_l){
    /** 
    * routine that computes the distance and cosine between a pair of cg mappings
    *
    * Parameters
    * ----------
    * 
    * `mapping`, `mapping_prime` : cg_mapping objects
    * 
    * `coupling_mat` : coupling matrix
    * 
    * `n_coord_at` : atomistic coordination number
    * 
    * `fr_id` : frame index
    * 
    * `f_out_l` : output filename
    */
    int i,j,knt;
    double dot_product = 0.0;
    knt = 0;
    for (i=0; i < mapping->n_at; i++){
      dot_product+=mapping->mapping[i]*mapping_prime->mapping[i];// self-interaction term
      for (j=i+1; j < mapping_prime->n_at; j++){
        dot_product += coupling_mat[knt]*mapping->mapping[i]*mapping_prime->mapping[j];
        knt += 1;
      }
    }
    //intermediate dot_product
    knt = 0;
    for (j=0; j < mapping_prime->n_at; j++){
      for (i=j+1; i < mapping->n_at; i++){
        dot_product += coupling_mat[knt]*mapping->mapping[i]*mapping_prime->mapping[j];
        knt += 1;
      }
    }
    fprintf(f_out_l,"dot product[%d] = %lf\n", fr_id,dot_product);
    //printf("normalized dot_product = %lf\n", dot_product/n_at);
    double cos = dot_product/(n_coord_at*sqrt(mapping->norms[fr_id])*sqrt(mapping_prime->norms[fr_id]));
    fprintf(f_out_l,"cosine[%d] %lf\n",fr_id, cos);
    double sq_dist = mapping->norms[fr_id] + mapping_prime->norms[fr_id] - 2*dot_product/n_coord_at;
    double distance = 0.0;
    fprintf(f_out_l,"sq_dist[%d] %lf\n",fr_id, sq_dist);
    if (sq_dist > 0.0){
        distance = sqrt(sq_dist);
    }
    else if (sq_dist < 0.0){
        if (-sq_dist < 0.000001){
            fprintf(f_out_l,"roundoff error, negative squared distance = %8.10lf\n", sq_dist);
            distance = 0.0;
        }
        else{
            fprintf(f_out_l,"negative squared distance %8.10lf..\n",sq_dist);
            distance = -1.0;
        }
    }
    fprintf(f_out_l,"distance[%d] %lf\n",fr_id, distance);
    return distance;
}

void compute_mapping_norms(traj *Trajectory, cg_mapping *mapping,FILE *f_out_l){
    /** 
    * routine that computes the norm of a mapping over a MD trajectory
    *
    * Parameters
    * ----------
    *
    * `Trajectory` : traj object
    * 
    * `mapping` : cg_mapping object
    * 
    * `f_out_l` : output filename
    */
    // declaring variables
    double *coupling_mat; // coupling matrix
    coupling_mat = d1t(Trajectory->n_at * (Trajectory->n_at - 1) / 2);
    mapping->norms = d1t(Trajectory->frames);
    int fr_id;
    double n_coord_at = -1.0;
    // iterating over frames
    for(fr_id=0; fr_id < Trajectory->frames; fr_id++) {
        compute_coupling_matrix(coupling_mat, Trajectory, fr_id, 1.9); // sigma is hard coded to half the distance between alpha carbons
        if (fr_id == 0){
            // compute atomistic coordination number to be used throughout the trajectory for normalisation
            n_coord_at = compute_atomistic_coord_number(coupling_mat,Trajectory,f_out_l);
        }
        compute_norm(mapping, coupling_mat, n_coord_at, fr_id, f_out_l);
    }
    free_d1t(coupling_mat);
}

void compute_mapping_distances(traj *Trajectory, cg_mapping *mapping,cg_mapping *mapping_prime, FILE *f_out_l){
    /** 
    * routine that computes the distances and cosines between two mappings provided by the user over a MD trajectory
    *
    * Parameters
    * ----------
    *
    * `Trajectory` : traj object
    * 
    * `mapping`, `mapping_prime` : cg_mapping objects
    * 
    * `f_out_l` : output filename
    */
    // declaring variables
    double *coupling_mat; // coupling matrix
    coupling_mat = d1t(Trajectory->n_at * (Trajectory->n_at - 1) / 2);
    mapping->norms = d1t(Trajectory->frames);
    mapping_prime->norms = d1t(Trajectory->frames);
    int fr_id;
    double n_coord_at = -1.0;
    // iterating over frames
    for(fr_id=0; fr_id < Trajectory->frames; fr_id++) {
        compute_coupling_matrix(coupling_mat, Trajectory, fr_id, 1.9); // sigma is hard coded to half the distance between alpha carbons
        if (fr_id == 0){
            // compute atomistic coordination number to be used throughout the trajectory for normalisation
            n_coord_at = compute_atomistic_coord_number(coupling_mat,Trajectory,f_out_l);
        }
        compute_norm(mapping, coupling_mat, n_coord_at, fr_id, f_out_l);
        compute_norm(mapping_prime, coupling_mat, n_coord_at, fr_id, f_out_l);
        fprintf(f_out_l, "norm(mapping)[%d] = %lf, norm(mapping_prime)[%d] = %lf\n", fr_id, mapping->norms[fr_id],fr_id,mapping_prime->norms[fr_id]);
        compute_distance(mapping,mapping_prime,coupling_mat,n_coord_at,fr_id,f_out_l);
    }
    free_d1t(coupling_mat);
}
void compute_mapping_distmat(traj *Trajectory, cg_mapping *mapping_matrix[], int nmaps,  FILE *f_out_l, char *distmat_filename){
    /** 
    * routine that computes the distance matrix between a set of mappings over a single structure
    *
    * Parameters
    * ----------
    *
    * `Trajectory` : traj object
    * 
    * `mappings_filename` : filename with the chosen mappings
    * 
    * `namps` : number of mappings
    * 
    * `f_out_l` : output filename
    */
    // coupling matrix
    double *coupling_mat; // coupling matrix
    coupling_mat = d1t(Trajectory->n_at * (Trajectory->n_at - 1) / 2);
    double n_coord_at;
    compute_coupling_matrix(coupling_mat, Trajectory, 0, 1.9); // sigma is hard coded to half the distance between alpha carbons
    printf("coupling_matrix[0] = %lf\n", coupling_mat[0]);
    n_coord_at = compute_atomistic_coord_number(coupling_mat,Trajectory,f_out_l);
    // mapping norms
    int map_idx;
    for (map_idx = 0; map_idx < nmaps; map_idx++){
        compute_norm(mapping_matrix[map_idx], coupling_mat, n_coord_at, 0, f_out_l);
    }
    // distance matrix calculation
    int map_jdx;
    int mapping_pairs = nmaps*(nmaps-1)/2;
    printf("%d mapping pairs\n", mapping_pairs);
    double *distances;
    distances = d1t(mapping_pairs);
    int knt = 0; 
    time_t seconds;
    time_t seconds_ref = time(NULL);
    for (map_idx=0; map_idx < nmaps; map_idx++) {
        for (map_jdx=map_idx+1; map_jdx < nmaps; map_jdx++) {
            distances[knt] = compute_distance(mapping_matrix[map_idx],mapping_matrix[map_jdx],coupling_mat,n_coord_at,0,f_out_l);
            knt += 1;
      }
    }
    // open distance matrix file TODO: routinize it
    FILE *f_distmat;
    f_distmat = open_file_w(distmat_filename);
    int sdt, sum_range;
    knt = 0;
    for (map_idx=0; map_idx < nmaps; map_idx++) {
        for (map_jdx=0; map_jdx < nmaps; map_jdx++) {
            if (map_jdx < map_idx){
                if (map_jdx%2 == 0){sum_range = (1+map_jdx)*map_jdx/2;}
                else{sum_range = (1+map_jdx)*(map_jdx/2) + (map_jdx/2 + 1);}
                sdt = map_jdx*(nmaps) - sum_range + (map_idx - map_jdx - 1);
                fprintf(f_distmat,"%lf ",distances[sdt]);
            }
            else if (map_jdx > map_idx){
                fprintf(f_distmat,"%lf ",distances[knt]);
                knt += 1;
                if (map_jdx == nmaps-1){fprintf(f_distmat,"\n");}
            }
            else{
                fprintf(f_distmat,"0.0 ");
            }
        }
    }
    fclose(f_distmat);
    free_d1t(distances);
    free_d1t(coupling_mat);
}

void compute_variances(int nclust, double *variances, int *cluster_list, int *cluster_list_idx, double *energies) {
    /** 
    * routine that computes the variance of the energies
    *
    * Parameters
    * ----------
    *
    * `nclust` : number of macrostates
    * 
    * `variances` : vector of variances
    * 
    * `cluster_list` : list of cluster IDs
    * 
    * `cluster_list_idx` : list of cluster indices
    * 
    * `energies` : array of energies
    */
    int i, cl, pop;
    double sum, sum1, average;
    for (cl = 0; cl < nclust; cl++) {
        average = 0.0;
        pop = cluster_list_idx[cl + 1] - cluster_list_idx[cl];
        if (pop == 0) {
            printf("POPULATION ZERO: cluster %d\n", cl);
            variances[cl] = 0.0;
            continue;
        }
        sum = 0.0;
        sum1 = 0.0;
        for (i = cluster_list_idx[cl]; i < cluster_list_idx[cl + 1]; i++) {
            sum = sum + energies[cluster_list[i]];
        }
        average = sum / (float) pop;
        /*  Compute  variance  */
        for (i = cluster_list_idx[cl]; i < cluster_list_idx[cl + 1]; i++) {
            sum1 = sum1 + pow((energies[cluster_list[i]] - average), 2);
        }
        variances[cl] = sum1 / (float) pop;
    }
}

void compute_pR(int nclust, double *p_R, int *cluster_list, int *cluster_list_idx, double *energies) {
    /** 
    * routine that computes the variance of the energies
    *
    * Parameters
    * ----------
    *
    * `nclust` : number of macrostates
    * 
    * `variances` : vector of variances
    * 
    * `cluster_list` : list of cluster IDs
    * 
    * `cluster_list_idx` : list of cluster indices
    * 
    * `energies` : array of energies
    */
    int i, cl, pop;
    double sum, average;
    for (cl = 0; cl < nclust; cl++) {
        average = 0.0;
        pop = cluster_list_idx[cl + 1] - cluster_list_idx[cl];
        if (pop == 0) {
            printf("POPULATION ZERO: cluster %d\n", cl);
        }
        sum = 0.0;
        for (i = cluster_list_idx[cl]; i < cluster_list_idx[cl + 1]; i++) {
            sum = sum + energies[cluster_list[i]];
        }
        //printf("pop[%d] = %d\n", cl, pop);
        /*  Compute  p_R  */
        p_R[cl] = sum / (float) pop;
        //printf("current p_R[%d] = %lf\n", cl, p_R[cl]);
    }
}

double get_smap(int frames, int curr_nclust, int *clusters, double *energies) {
    /** 
    * routine that computes the observable given the current `nclust` and the current `clusters`
    *
    * Parameters
    * ----------
    *
    * `frames` : number of frames
    * 
    * `curr_nclust` : current index of CG macrostate
    * 
    * `clusters` : list of cluster IDs
    * 
    * `energies` : array of energies
    */
    int *cluster_list;//,*cluster_list_idx;
    int cluster_list_idx[curr_nclust + 1];
    cluster_list = i1t(frames);
    compute_clusters_list(clusters, cluster_list, cluster_list_idx, frames, curr_nclust);
    // gives cluster_list [7 8 9 1 2 3 4 5 6 10 11 13 14 15 12] and cluster_list_idx 0 3 9 14.
    // you know (variable frames) that you have 15 configurations, so the last cluster ends there and has population
    // frames - cluster_list_idx[nclust-1]
    //////////////////////////////////////////////////////////////////////////////////
    // variances
    double *variances;
    variances = d1t(curr_nclust);
    compute_variances(curr_nclust, variances, cluster_list, cluster_list_idx, energies);
    double *smap_R;
    smap_R = d1t(curr_nclust);
    int cl;
    double smap = 0.0;
    //printf("computed variances\n");
    for (cl = 0; cl < curr_nclust; cl++) {
        smap_R[cl] = 0.5 * variances[cl] * (cluster_list_idx[cl + 1] - cluster_list_idx[cl]) / (float) frames;
        //printf("smap[%d] = %lf\n", cl, smap_R[cl]);
        smap += smap_R[cl];
    }
    //printf("smap(ncl=%d) = %lf\n", curr_nclust, smap);
    free_d1t(variances);
    free_d1t(smap_R);
    free_i1t(cluster_list);
    return smap;
}

double get_kl(int frames, int curr_nclust, int *clusters, double *energies) {
    /** 
    * routine that computes the observable given the current `nclust` and the current `clusters`
    *
    * Parameters
    * ----------
    *
    * `frames` : number of frames
    * 
    * `curr_nclust` : current index of CG macrostate
    * 
    * `clusters` : list of cluster IDs
    * 
    * `energies` : array of energies
    */
    int *cluster_list;//,*cluster_list_idx;
    int cluster_list_idx[curr_nclust + 1];
    cluster_list = i1t(frames);
    compute_clusters_list(clusters, cluster_list, cluster_list_idx, frames, curr_nclust);
    // probabilities
    double *p_R;
    p_R = d1t(curr_nclust);
    compute_pR(curr_nclust, p_R, cluster_list, cluster_list_idx, energies);
    int fr;
    //for (cl = 0; cl < curr_nclust; cl++) {printf("p_R[%d] = %lf\n", cl, p_R[cl]);}
    double smap = 0.0;
    double smap_r ;
    //printf("computed variances\n");
    for (fr = 0; fr < frames; fr++) {
        //printf("clusters[%d] = %d, p_r[%d] = %lf p_R[%d] = %lf\n",fr, clusters[fr], fr, energies[fr], clusters[fr], p_R[clusters[fr]-1]);
        smap_r = energies[fr] * log(energies[fr]/p_R[clusters[fr]-1]);
        //printf("smap_r[%d] = %lf\n", fr, smap_r);
        smap += smap_r;
    }
    //printf("smap(ncl=%d) = %lf\n", curr_nclust, smap);
    //free_d1t(variances);
    free_d1t(p_R);
    free_i1t(cluster_list);
    return smap;
}

void compute_smap_spins(spin_traj *Trajectory, cg_mapping *mapping){

    int ncl = spin_clustering(Trajectory->traj_coords, mapping->mapping, mapping->n_at, Trajectory->frames, mapping->clusters, mapping->n_cg);
    mapping->smap = get_kl(Trajectory->frames, ncl, mapping->clusters, Trajectory->energies);
    printf("curr. smap = %lf\n", mapping->smap);
}

void overall_compute_smap(alignments *align, clust_params *clustering, traj *Trajectory, cg_mapping *mapping, int verbose, int kl_flag){
    /** 
    * routine that calls `get_smap` with the correct parameters
    * 
    * Parameters
    * ----------
    *
    * `rmsd_mat` : condensed matrix of pairwise RMSDs
    * 
    * `clustering` : clust_params object
    * 
    * `Trajectory` : traj object
    * 
    * `mapping` : cg_mapping object
    * 
    * `verbose` : tunes the level of verbosity
    * 
    * `f_out` : output filename 
    */
    double **Z;
    mapping->smap = 0.0;
    int k;
    if (clustering->crit != 1) {
        Z = d2t(Trajectory->frames - 1, 4);
        // linkage matrix
        hierarchical_clustering(align->rmsd_mat, Trajectory->frames, Trajectory->pairs, mapping->size, Z);
        // observable computation
        if (clustering->crit == 0) {
            if (verbose == 1) {
                printf("criterion %d (maxclust): clustering of dist_mat (len %d, %d frames) into %d clusters \n", clustering->crit,
                       Trajectory->pairs, Trajectory->frames, clustering->ncl);
            }
            cluster_maxclust_dist(Z, mapping->clusters, Trajectory->frames,  clustering->ncl);
            if (kl_flag == 0){mapping->smap = get_smap(Trajectory->frames, clustering->ncl, mapping->clusters, Trajectory->energies);}
            else{mapping->smap = get_kl(Trajectory->frames, clustering->ncl, mapping->clusters, Trajectory->energies);}
            if (verbose == 1) {
                int cl_id;
                printf("clusters with nclust = %d\n", clustering->ncl);
                for (cl_id = 0; cl_id < Trajectory->frames; cl_id++) { printf("%d ", mapping->clusters[cl_id]); }
                printf("\n");
            }
        } else if (clustering->crit == 3) {
            // It means it's distance-based clustering
            if (verbose == 1) {
                printf("Warning: criterion %d for computing observable, distance-based clustering\nNOT COMPATIBLE WITH Giulini et. al (2020)\n",
                       clustering->crit);
                printf("criterion %d (distance): clustering of dist_mat (len %d, %d frames) using cophenetic distance %lf \n",
                       clustering->crit, Trajectory->pairs, Trajectory->frames, clustering->c_distance);
            }
            cluster_dist(Z, mapping->clusters, clustering->c_distance, Trajectory->frames);
            int ncl = 1;
            int snt = 0;
            int ncls[Trajectory->frames];
            int present = 1;
            ncls[0] = mapping->clusters[0];
            for (k = 1; k < Trajectory->frames; k++) {
                for (snt = 0; snt < ncl; snt++) {
                    if (ncls[snt] == mapping->clusters[k]) { present = 0; }
                }
                if (present == 1) {
                    ncls[ncl] = mapping->clusters[k];
                    ncl += 1;
                }
                present = 1;
            }
            if (verbose == 1) {
                int cl_id;
                printf("clusters with nclust = %d\n", ncl);
                for (cl_id = 0; cl_id < Trajectory->frames; cl_id++) { printf("%d ", mapping->clusters[cl_id]); }
                printf("\n");
            }
            if (kl_flag == 0){mapping->smap = get_smap(Trajectory->frames, ncl, mapping->clusters, Trajectory->energies);}
            else{mapping->smap = get_kl(Trajectory->frames, ncl, mapping->clusters, Trajectory->energies);}
        } else if (clustering->crit == 2) {
            int cl_inside = 3;
            double stride = (clustering->max_ncl - clustering->min_ncl) / ((float) cl_inside + 1);
            int clusts[cl_inside + 2];
            int dummy, d;
            clusts[0] = clustering->min_ncl;
            for (dummy = 1; dummy < cl_inside + 1; dummy++) { clusts[dummy] = (int) clusts[dummy - 1] + stride; }
            clusts[cl_inside + 1] = clustering->max_ncl;
            for (dummy = 0; dummy < cl_inside + 2; dummy++) {
                cluster_maxclust_dist(Z, mapping->clusters, Trajectory->frames, clusts[dummy]);
                if (kl_flag == 0){mapping->smap += get_smap(Trajectory->frames, clusts[dummy], mapping->clusters, Trajectory->energies);}
                else{mapping->smap += get_kl(Trajectory->frames, clusts[dummy], mapping->clusters, Trajectory->energies);}
            }
            mapping->smap = mapping->smap / ((float) cl_inside + 2);
        }
    }
    else{
        // fast clustering
        Z = d2t(Trajectory->eff_frames - 1 , 4);
        hierarchical_clustering(align->rmsd_mat, Trajectory->eff_frames, Trajectory->pairs, mapping->size, Z);
        cluster_maxclust_dist(Z, mapping->clusters, Trajectory->eff_frames,  clustering->ncl);
        int cl_id;
        // reordering clusters. last frame has to be treated differently
        mapping->clusters[Trajectory->frames-1] = mapping->clusters[Trajectory->eff_frames-1];
        mapping->clusters[Trajectory->eff_frames-1] = 0;
        for (cl_id = Trajectory->eff_frames - 2; cl_id > 0; cl_id--){
            mapping->clusters[cl_id*Trajectory->stride] = mapping->clusters[cl_id];
            mapping->clusters[cl_id] = 0;
        }
        // assigning missing points
        int prev = Trajectory->frames -1;
        int next = Trajectory->frames - 1 - (Trajectory->frames-1)%Trajectory->stride;
        //double rmsd_next , rmsd_prev;
        for (cl_id = Trajectory->frames -2; cl_id > 0; cl_id--){
            if (Trajectory->strides[cl_id] == 0){
                if (mapping->clusters[prev] == mapping->clusters[next]){mapping->clusters[cl_id] = mapping->clusters[prev];}
                else{
                    //printf("cl_id = %d pivot clusters (ids %d %d) are not equal (%d %d), choose the closer\n", cl_id, prev, next, mapping->clusters[prev] , mapping->clusters[next]);
                    //printf("rmsd_prev = %8.6lf rmsd_next %8.6lf\n", align->rmsd_vector[cl_id*2+1] , align->rmsd_vector[cl_id*2]);
                    if (align->rmsd_vector[cl_id*2+1] < align->rmsd_vector[cl_id*2]){mapping->clusters[cl_id] = mapping->clusters[prev];}
                    else{mapping->clusters[cl_id] = mapping->clusters[next];}
                }
            }
            else if(Trajectory->strides[cl_id] == 1){
                prev = next;
                next = prev - Trajectory->stride;
            }
        }
        if (verbose == 1) {
            printf("clusters with nclust = %d\n", clustering->ncl);
            for (cl_id = 0; cl_id < Trajectory->frames; cl_id++) { printf("%d ", mapping->clusters[cl_id]); }
            printf("\n");
        }
        //cluster_maxclust_dist(Z, mapping->clusters, Trajectory->frames,  clustering->ncl);
        if (kl_flag == 0){mapping->smap = get_smap(Trajectory->frames, clustering->ncl, mapping->clusters, Trajectory->energies);}
        else{mapping->smap = get_kl(Trajectory->frames, clustering->ncl, mapping->clusters, Trajectory->energies);}
    }
    free_d2t(Z);
}
