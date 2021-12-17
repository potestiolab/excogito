#ifndef HDR_OBS
#define HDR_OBS

#include <stdio.h>
#include <mapping.h>
#include <hierarchical_clustering.h>
#include <traj.h>
#include <alignment.h>
// observables.c

void compute_coupling_matrix(double *coupling_mat, traj *Trajectory, int fr_id, float sigma);

double compute_atomistic_coord_number(double *coupling_mat,traj *Trajectory,FILE *f_out_l);

void compute_norm(cg_mapping *mapping, double *coupling_mat, double n_coord_at, int fr_id, FILE *f_out_l);

double compute_distance(cg_mapping *mapping, cg_mapping *mapping_prime, double *coupling_mat, double n_coord_at, int fr_id, FILE *f_out_l);

void compute_mapping_norms(traj *Trajectory, cg_mapping *mapping,FILE *f_out_l);

void compute_mapping_distances(traj *Trajectory, cg_mapping *mapping,cg_mapping *mapping_prime, FILE *f_out_l);

void compute_mapping_distmat(traj *Trajectory, cg_mapping *mapping_matrix[], int nmaps, FILE *f_out_l, char *distmat_filename);

void compute_variances(int Nclust, double *variances, int *cluster_list, int *cluster_list_idx, double *energies);

double get_smap(int frames, int curr_nclust, int *clusters, double *energies);

//void overall_compute_smap(double *rmsd_mat, clust_params *clustering, traj *Trajectory, cg_mapping *mapping, int verbose);

void overall_compute_smap(alignments *align, clust_params *clustering, traj *Trajectory, cg_mapping *mapping, int verbose, int kl_flag);


#endif
