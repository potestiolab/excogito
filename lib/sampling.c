/**
* \file sampling
* \brief library of functions containing the available sampling strategies
*/

#include <stdio.h>
#include <stdlib.h>
#include <my_malloc.h>
#include <math.h>
#include <sampling.h>
#include <alignment.h>
#include <mapping.h>
#include <observables.h>

void my_make_a_move(cg_mapping *old_mapping, cg_mapping *new_mapping, int rem_add[2]) {
    /**
    * function that swaps two atoms inside a CG mapping
    * Parameters
    * ----------
    *
    * `old_mapping` : cg_mapping object
    * 
    * `new_mapping` : cg_mapping object
    * 
    * `rem_add` : vector of length 2 containing the removed and added atom index
    */
    int i, rd, new_rd;
    for (i = 0; i < old_mapping->n_at; i++) { new_mapping->mapping[i] = old_mapping->mapping[i]; }
    rd = rand() % old_mapping->n_at;
    if (new_mapping->mapping[rd] == 1) {
        rem_add[0] = rd;
        new_mapping->mapping[rd] = 0;
        while (1) {
            new_rd = rand() % old_mapping->n_at;
            if (new_mapping->mapping[new_rd] == 0 & new_rd != rem_add[0]) {
                rem_add[1] = new_rd;
                new_mapping->mapping[new_rd] = 1;
                break;
            }
        }
    } else {
        rem_add[1] = rd;
        new_mapping->mapping[rd] = 1;
        while (1) {
            new_rd = rand() % old_mapping->n_at;
            if (new_mapping->mapping[new_rd] == 1 & new_rd != rem_add[1]) {
                rem_add[0] = new_rd;
                new_mapping->mapping[new_rd] = 0;
                break;
            }
        }
    }
}

void accept_move(int rem_add[2], cg_mapping *mapping, cg_mapping *new_mapping, alignments *align, alignments *new_align, traj *Trajectory, clust_params *clustering, double **new_coefficients_matrix){
    /**
    * routine that accepts a Simulated Annealing move. It updates all the relevant observables.
    *
    * Parameters
    * ----------
    *
    * `rem_add` : vector of length 2 containing the removed and added atom index
    * 
    * `alignments` : align object
    * 
    * `mapping` : cg_mapping object
    * 
    * `verbose` : tunes the level of verbosity
    */
    int dummy, k;
    mapping->mapping[rem_add[0]] = 0;
    mapping->mapping[rem_add[1]] = 1;
    mapping->smap = new_mapping->smap;
    for (dummy = 0; dummy < Trajectory->pairs; dummy++) { align->rmsd_mat[dummy] = new_align->rmsd_mat[dummy]; }
    if (clustering->crit==3){
        for (dummy = 0; dummy < Trajectory->frames*2; dummy++){align->rmsd_vector[dummy] = new_align->rmsd_vector[dummy];}
    }   
}

void accept_move_spins(int rem_add[2], cg_mapping *mapping, cg_mapping *new_mapping){
    /**
    * routine that accepts a Simulated Annealing Spins move. It updates all the relevant observables.
    *
    * Parameters
    * ----------
    *
    * `rem_add` : vector of length 2 containing the removed and added atom index
    * 
    * `mapping` : cg_mapping object
    * 
    * `new_mapping` : cg_mapping object
    */
    mapping->mapping[rem_add[0]] = 0;
    mapping->mapping[rem_add[1]] = 1;
    mapping->smap = new_mapping->smap;
}

double tzero_estimation_spins(spin_traj *Trajectory, int cgnum, FILE *f_out_l) {
    /**
    * routine that makes *nrun* unbiased moves. It is used to estimate the starting temperature for Simulated Annealing.
    *
    * Parameters
    * ----------
    *
    * `Trajectory` : traj object
    * 
    * `alignments` : align object
    * 
    * `mapping` : cg_mapping object
    * 
    * `verbose` : tunes the level of verbosity
    */
    printf("tzero_estimation\n");
    int n_moves = 20;
    double *deltas; // vector of delta_smap values
    deltas = d1t(n_moves); // fixed number of moves generated
    int dummy, nr;
    int rem_add[2]; // vector containing removed and added atom
    // mapping
    printf("allocating mappings\n");
    cg_mapping *mapping = malloc (sizeof(cg_mapping));
    mapping->n_at = Trajectory->n_at;
    mapping->n_cg = cgnum;
    mapping->mapping = i1t(Trajectory->n_at);
    mapping->smap = -1.0;
    mapping->clusters = i1t(Trajectory->frames);
    mapping->size = i1t(Trajectory->frames);
    cg_mapping *new_mapping = malloc(sizeof(cg_mapping));
    new_mapping->n_at = Trajectory->n_at;
    new_mapping->n_cg = cgnum;
    new_mapping->smap = -1.0;
    new_mapping->mapping = i1t(Trajectory->n_at);
    new_mapping->clusters = i1t(Trajectory->frames);
    new_mapping->size = i1t(Trajectory->frames);
    int d_idx = 0;
    generate_random_mapping(mapping,f_out_l);
    convert_mapping(mapping, f_out_l);
    compute_smap_spins(Trajectory, mapping, verbose);
    printf("computing start smap\n");
    //overall_compute_smap(align, clustering, Trajectory, mapping, verbose, kl_flag);
    for (nr = 0; nr < n_moves; nr++) {
        my_make_a_move(mapping, new_mapping, rem_add);
        convert_mapping(new_mapping, f_out_l);
        // if (clustering->crit==3){correct_rmsd_fastclust(align, Trajectory, align, mapping->n_cg, rem_add[0],rem_add[1]);}
        // else if (kl_flag!=2){correct_rmsd(align, Trajectory, align, mapping->n_cg, rem_add[0],rem_add[1]);}
        // compute new smap
        compute_smap_spins(Trajectory, new_mapping, verbose);
        // add delta
        deltas[d_idx] = fabs(new_mapping->smap - mapping->smap);
        // updating mapping with new smap and new site
        mapping->smap = new_mapping->smap;
        mapping->mapping[rem_add[0]] = 0;
        mapping->mapping[rem_add[1]] = 1;
        // increasing index
        d_idx = d_idx + 1;
    }
    double avg_delta = 0.0;
    fprintf(f_out_l, "collected deltas:\n");
    for (d_idx = 0; d_idx < n_moves; d_idx++){
        fprintf(f_out_l, "%lf ", deltas[d_idx]);
        avg_delta += deltas[d_idx];
    }
    avg_delta = avg_delta/(n_moves*1.0);
    fprintf(f_out_l, "\naverage delta = %8.6lf\n", avg_delta);
    double corr_t_zero = -avg_delta/(log(0.75));
    free_mapping(mapping);
    free_mapping(new_mapping);
    free_d1t(deltas);
    return corr_t_zero;
}

double tzero_estimation(traj *Trajectory, clust_params *clustering, int cgnum, int rsd, int verbose, int kl_flag, FILE *f_out_l) {
    /**
    * routine that makes *nrun* unbiased moves. It is used to estimate the starting temperature for Simulated Annealing.
    *
    * Parameters
    * ----------
    *
    * `Trajectory` : traj object
    * 
    * `alignments` : align object
    * 
    * `mapping` : cg_mapping object
    * 
    * `verbose` : tunes the level of verbosity
    */
    printf("tzero_estimation");
    double *deltas; // vector of delta_smap values
    deltas = d1t(20); // fixed number of moves generated
    int dummy, nr;
    int rem_add[2]; // vector containing removed and added atom
    // alignment object
    alignments *align = malloc (sizeof(alignments));
    align->rmsd_mat = d1t(Trajectory->pairs);
    align->rotation_matrices = d2t(Trajectory->pairs, 9);
    align->coms = d2t(Trajectory->frames, 3);
    align->rsd = rsd;
    align->rmsd_vector = d1t(Trajectory->frames * 2);
    align->rotation_matrices_vector = d2t(Trajectory->frames * 2, 9);
    double **new_coefficients_matrix;
    // mapping
    cg_mapping *mapping = malloc (sizeof(cg_mapping));
    mapping->n_at = Trajectory->n_at;
    mapping->n_cg = cgnum;
    mapping->mapping = i1t(Trajectory->n_at);
    mapping->smap = -1.0;
    mapping->clusters = i1t(Trajectory->frames);
    mapping->size = i1t(Trajectory->frames);
    cg_mapping *new_mapping = malloc(sizeof(cg_mapping));
    new_mapping->n_at = Trajectory->n_at;
    new_mapping->n_cg = cgnum;
    new_mapping->smap = -1.0;
    new_mapping->mapping = i1t(Trajectory->n_at);
    new_mapping->clusters = i1t(Trajectory->frames);
    new_mapping->size = i1t(Trajectory->frames);
    int d_idx = 0;
    generate_random_mapping(mapping,f_out_l);
    if (clustering->crit==3){cycle_alignment_fastclust(Trajectory, align, mapping);}
    else{cycle_alignment(Trajectory, align, mapping);}
    overall_compute_smap(align, clustering, Trajectory, mapping, verbose, kl_flag);
    for (nr = 0; nr < 20; nr++) {
        my_make_a_move(mapping, new_mapping, rem_add);
        convert_mapping(new_mapping, f_out_l);
        if (clustering->crit==3){correct_rmsd_fastclust(align, Trajectory, align, mapping->n_cg, rem_add[0],rem_add[1]);}
        else{correct_rmsd(align, Trajectory, align, mapping->n_cg, rem_add[0],rem_add[1]);}
        // compute new smap
        overall_compute_smap(align, clustering, Trajectory, new_mapping, verbose, kl_flag);
        // add delta
        deltas[d_idx] = fabs(new_mapping->smap - mapping->smap);
        // updating mapping with new smap and new site
        mapping->smap = new_mapping->smap;
        mapping->mapping[rem_add[0]] = 0;
        mapping->mapping[rem_add[1]] = 1;
        // increasing index
        d_idx = d_idx + 1;
    }
    double avg_delta = 0.0;
    fprintf(f_out_l, "collected deltas:\n");
    for (d_idx = 0; d_idx < 20; d_idx++){
        fprintf(f_out_l, "%lf ", deltas[d_idx]);
        avg_delta += deltas[d_idx];
    }
    avg_delta = avg_delta/(20*1.0);
    fprintf(f_out_l, "\naverage delta = %8.6lf\n", avg_delta);
    double corr_t_zero = -avg_delta/(log(0.75));
    free_mapping(mapping);
    free_mapping(new_mapping);
    free_alignment(align);
    free_d1t(deltas);
    return corr_t_zero;
}

void simulated_annealing(traj *Trajectory, clust_params *clustering, MC_params *SA_params, int cgnum, int rsd, int verbose, int kl_flag, FILE *f_out_l){
    /**
    * simulated annealing optimisation
    * Parameters
    * ----------
    *
    * `Trajectory` : traj object
    * 
    * `alignments` : align object
    * 
    * `mapping` : cg_mapping object
    * 
    * `SA_params` : set of Monte Carlo parameters
    * 
    * `verbose` : tunes the level of verbosity
    * 
    * `f_out_l` : output filename
    */
    fprintf(f_out_l, "simulated_annealing optimisation\n");
    fprintf(f_out_l, "kl_flag is %d: ", kl_flag);
    if (kl_flag == 1){
        fprintf(f_out_l,"explicitly calculating Kullback-Leibler divergences\n");
    }
    else if (kl_flag == 0){
        fprintf(f_out_l,"calculating approximated expressions\n");
    }
    // setting verbose to 0 if MC_steps*frames > 10^6
    if (SA_params->MC_steps*Trajectory->frames > 1000000){
        printf("MC_steps*frames (%d) > 10^6: setting verbose to 0 to avoid wasting too much disk space\n", SA_params->MC_steps*Trajectory->frames);
        verbose = 0;
    }
    cg_mapping *mapping = malloc (sizeof(cg_mapping));
    cg_mapping *new_mapping = malloc (sizeof(cg_mapping));
    cg_mapping *lowest_mapping = malloc (sizeof(cg_mapping));
    // variables
    int mc; // indexes of mc step and nrun
    int k; // counter
    double temp;
    int dummy;
    int rem_add[2];
    double **new_coefficients_matrix;
    //new_rmsd_mat = d1t(Trajectory->pairs);
    // MC
    double r, p;
    alignments *align = malloc (sizeof(alignments));
    align->rmsd_mat = d1t(Trajectory->pairs);
    align->rotation_matrices = d2t(Trajectory->pairs, 9);
    align->coms = d2t(Trajectory->frames, 3);
    align->rsd = rsd;
    align->rmsd_vector = d1t(Trajectory->frames * 2);
    align->rotation_matrices_vector = d2t(Trajectory->frames * 2, 9);
    // new alignment
    alignments *new_align =  malloc (sizeof(alignments));
    new_align->rmsd_mat = d1t(Trajectory->pairs);
    new_align->rmsd_vector = d1t(Trajectory->frames * 2);
    // mapping
    mapping->n_at = Trajectory->n_at;
    mapping->n_cg = cgnum;
    mapping->mapping = i1t(Trajectory->n_at);
    mapping->smap = -1.0;
    mapping->clusters = i1t(Trajectory->frames);
    mapping->size = i1t(Trajectory->frames);
    new_mapping->n_at = mapping->n_at;
    new_mapping->n_cg = cgnum;
    new_mapping->mapping = i1t(mapping->n_at);
    new_mapping->clusters = i1t(Trajectory->frames);
    new_mapping->size = i1t(Trajectory->frames);
    new_mapping->smap = -1.0;
    lowest_mapping->n_at = mapping->n_at;
    lowest_mapping->n_cg = cgnum;
    lowest_mapping->mapping = i1t(mapping->n_at);
    lowest_mapping->clusters = i1t(Trajectory->frames);
    lowest_mapping->size = i1t(Trajectory->frames);
    lowest_mapping->smap = -1.0;
    printf("allocated variables\n");
    // random mapping and rotation matrices (and first rmsd_mat)
    printf("generating random mapping\n");
    generate_random_mapping(mapping, f_out_l);
    fprintf(f_out_l, "Initial mapping\n");
    convert_mapping(mapping, f_out_l);
    // if fast clustering
    if (clustering->crit==3){
        cycle_alignment_fastclust(Trajectory, align, mapping);            
    }
    else{cycle_alignment(Trajectory, align, mapping);}
    // calculate start mapping entropy
    overall_compute_smap(align, clustering, Trajectory, mapping, verbose,kl_flag);
    fprintf(f_out_l, "start observable is %lf\n", mapping->smap);
    fprintf(f_out_l, "starting opt with T0 = %lf and decay_time %lf\n", SA_params->t_zero, SA_params->decay_time);
    update_mapping(mapping, lowest_mapping, Trajectory->frames);
    // optimisation
    for (mc = 0; mc < SA_params->MC_steps; mc++) {
        temp = SA_params->t_zero * exp(-((double) mc) / SA_params->decay_time);
        fprintf(f_out_l, "MC_step %d temp = %lf \n", mc, temp);
        // print mapping every 250 SA steps
        if (mc % 50 == 0 & mc != 0 & verbose == 1) {
            fprintf(f_out_l, "mapping\n");
            convert_mapping(mapping, f_out_l);
        }
        // check if we have to update rotation matrices and coms
        if ((mc % SA_params->rotmats_period == 0) & (mc != 0) & clustering->crit != 4) {
            fprintf(f_out_l, "Updating rotation matrices and coms\n");
            if (clustering->crit==3){cycle_alignment_stride(Trajectory, align, mapping);}
            else{cycle_alignment(Trajectory, align, mapping);}
        }
        my_make_a_move(mapping, new_mapping, rem_add);
        fprintf(f_out_l, "move made: removed %d added %d\n", rem_add[0], rem_add[1]);
        // new rmsd matrix
        if (clustering->crit==3){
            correct_rmsd_fastclust(new_align, Trajectory, align, mapping->n_cg, rem_add[0],rem_add[1]);}
        else{correct_rmsd(new_align, Trajectory, align, mapping->n_cg, rem_add[0],rem_add[1]);}
        // observable
        overall_compute_smap(new_align, clustering, Trajectory, new_mapping, verbose, kl_flag);
        fprintf(f_out_l, "new_smap %lf\n", new_mapping->smap);
        // MC rule
        if ((new_mapping->smap - mapping->smap) < 0) {
            fprintf(f_out_l, "move accepted\n");
            accept_move(rem_add,mapping,new_mapping,align,new_align,Trajectory,clustering,new_coefficients_matrix);
            // check lowest
            if (mapping->smap < lowest_mapping->smap) {
                fprintf(f_out_l, "Lowest mapping entropy found\n");
                update_mapping(mapping, lowest_mapping,Trajectory->frames);
                fprintf(f_out_l, "Printing lowest_smap mapping\n");
                convert_mapping(lowest_mapping, f_out_l);
            }
        } else {
            r = rand() / (double) RAND_MAX;
            p = exp((mapping->smap - new_mapping->smap) / temp);
            if (verbose == 1){fprintf(f_out_l, "r= %lf p = %lf\n", r, p);}
            if (p > r) {
                fprintf(f_out_l, "p > r move accepted\n");
                accept_move(rem_add,mapping,new_mapping,align,new_align,Trajectory,clustering,new_coefficients_matrix);
                }
            else {
                fprintf(f_out_l, "move rejected\n");
            }
        }
        fprintf(f_out_l, "SMAP=%lf\n", mapping->smap);
    }
    fprintf(f_out_l, "END OF SIMULATED ANNEALING OPTIMISATION\n");
    fprintf(f_out_l, "lowest_smap = %lf\n", lowest_mapping->smap);
    fprintf(f_out_l, "lowest_mapping = \n");
    for (dummy = 0; dummy < lowest_mapping->n_at; dummy++) { fprintf(f_out_l, "%d ", lowest_mapping->mapping[dummy]); }
    convert_mapping(lowest_mapping, f_out_l);
    fprintf(f_out_l, "last_smap %lf\n", mapping->smap);
    fprintf(f_out_l, "last_mapping\n");
    for (dummy = 0; dummy < mapping->n_at; dummy++) { fprintf(f_out_l, "%d ", mapping->mapping[dummy]); }
    convert_mapping(mapping, f_out_l);
    // freeing variables
    free_mapping(mapping);
    free_mapping(new_mapping);
    free_mapping(lowest_mapping);
    free_alignment(align);
    free_new_alignment(new_align);
}

void simulated_annealing_spins(traj *Trajectory, MC_params *SA_params, int cgnum, int verbose, FILE *f_out_l){
    /**
    * simulated annealing optimisation of a spin system
    * Parameters
    * ----------
    *
    * `Trajectory` : traj object
    * 
    * `mapping` : cg_mapping object
    * 
    * `SA_params` : set of Monte Carlo parameters
    * 
    * `verbose` : tunes the level of verbosity
    * 
    * `f_out_l` : output filename
    */
    fprintf(f_out_l, "simulated_annealing_spins optimisation\n");
    // setting verbose to 0 if MC_steps*frames > 10^6
    if (SA_params->MC_steps*Trajectory->frames > 1000000){
        printf("MC_steps*frames (%d) > 10^6: setting verbose to 0 to avoid wasting too much disk space\n", SA_params->MC_steps*Trajectory->frames);
        verbose = 0;
    }
    cg_mapping *mapping = malloc (sizeof(cg_mapping));
    cg_mapping *new_mapping = malloc (sizeof(cg_mapping));
    cg_mapping *lowest_mapping = malloc (sizeof(cg_mapping));
    // variables
    int mc; // indexes of mc step and nrun
    int k; // counter
    double temp;
    int dummy;
    int rem_add[2];
    // MC
    double r, p;
    // mapping
    mapping->n_at = Trajectory->n_at;
    mapping->n_cg = cgnum;
    mapping->mapping = i1t(Trajectory->n_at);
    mapping->smap = -1.0;
    mapping->clusters = i1t(Trajectory->frames);
    mapping->size = i1t(Trajectory->frames);
    new_mapping->n_at = mapping->n_at;
    new_mapping->n_cg = cgnum;
    new_mapping->mapping = i1t(mapping->n_at);
    new_mapping->clusters = i1t(Trajectory->frames);
    new_mapping->size = i1t(Trajectory->frames);
    new_mapping->smap = -1.0;
    lowest_mapping->n_at = mapping->n_at;
    lowest_mapping->n_cg = cgnum;
    lowest_mapping->mapping = i1t(mapping->n_at);
    lowest_mapping->clusters = i1t(Trajectory->frames);
    lowest_mapping->size = i1t(Trajectory->frames);
    lowest_mapping->smap = -1.0;
    printf("allocated variables\n");
    // random mapping and rotation matrices (and first rmsd_mat)
    printf("generating random mapping\n");
    generate_random_mapping(mapping, f_out_l);
    fprintf(f_out_l, "Initial mapping\n");
    convert_mapping(mapping, f_out_l);
    // calculate start mapping entropy
    compute_smap_spins(Trajectory, mapping, verbose);
    fprintf(f_out_l, "start observable is %lf\n", mapping->smap);
    fprintf(f_out_l, "starting opt with T0 = %lf and decay_time %lf\n", SA_params->t_zero, SA_params->decay_time);
    update_mapping(mapping, lowest_mapping, Trajectory->frames);
    // optimisation
    for (mc = 0; mc < SA_params->MC_steps; mc++) {
        temp = SA_params->t_zero * exp(-((double) mc) / SA_params->decay_time);
        fprintf(f_out_l, "MC_step %d temp = %lf \n", mc, temp);
        // print mapping every 50 SA steps
        if (mc % 50 == 0 & mc != 0 & verbose == 1) {
            fprintf(f_out_l, "mapping\n");
            convert_mapping(mapping, f_out_l);
        }
        my_make_a_move(mapping, new_mapping, rem_add);
        fprintf(f_out_l, "move made: removed %d added %d\n", rem_add[0], rem_add[1]);
        // observable
        compute_smap_spins(Trajectory, new_mapping, verbose);
        fprintf(f_out_l, "new_smap %lf\n", new_mapping->smap);
        // MC rule
        if ((new_mapping->smap - mapping->smap) < 0) {
            fprintf(f_out_l, "move accepted\n");
            accept_move_spins(rem_add,mapping,new_mapping);
            // check lowest
            if (mapping->smap < lowest_mapping->smap) {
                fprintf(f_out_l, "Lowest mapping entropy found\n");
                update_mapping(mapping, lowest_mapping,Trajectory->frames);
                fprintf(f_out_l, "Printing lowest_smap mapping\n");
                convert_mapping(lowest_mapping, f_out_l);
            }
        } else {
            r = rand() / (double) RAND_MAX;
            p = exp((mapping->smap - new_mapping->smap) / temp);
            if (verbose == 1){fprintf(f_out_l, "r= %lf p = %lf\n", r, p);}
            if (p > r) {
                fprintf(f_out_l, "p > r move accepted\n");
                accept_move_spins(rem_add,mapping,new_mapping);
                }
            else {
                fprintf(f_out_l, "move rejected\n");
            }
        }
        fprintf(f_out_l, "SMAP=%lf\n", mapping->smap);
    }
    fprintf(f_out_l, "END OF SIMULATED ANNEALING OPTIMISATION\n");
    fprintf(f_out_l, "lowest_smap = %lf\n", lowest_mapping->smap);
    fprintf(f_out_l, "lowest_mapping = \n");
    for (dummy = 0; dummy < lowest_mapping->n_at; dummy++) { fprintf(f_out_l, "%d ", lowest_mapping->mapping[dummy]); }
    convert_mapping(lowest_mapping, f_out_l);
    fprintf(f_out_l, "last_smap %lf\n", mapping->smap);
    fprintf(f_out_l, "last_mapping\n");
    for (dummy = 0; dummy < mapping->n_at; dummy++) { fprintf(f_out_l, "%d ", mapping->mapping[dummy]); }
    convert_mapping(mapping, f_out_l);
    // freeing variables
    free_mapping(mapping);
    free_mapping(new_mapping);
    free_mapping(lowest_mapping);
}