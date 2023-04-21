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
#include <time.h>

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

void accept_move(int rem_add[2], cg_mapping *mapping, cg_mapping *new_mapping, traj *Trajectory){
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
    mapping->res = new_mapping->res;
    /*for (dummy = 0; dummy < Trajectory->pairs; dummy++) { align->rmsd_mat[dummy] = new_align->rmsd_mat[dummy]; }
    if (clustering->crit==3){
        for (dummy = 0; dummy < Trajectory->frames*2; dummy++){align->rmsd_vector[dummy] = new_align->rmsd_vector[dummy];}
    }   */										//(!)
}

double tzero_estimation(traj *Trajectory, int cgnum, int rsd, int verbose, int kl_flag, FILE *f_out_l) {		//(!)
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
    double *deltas; // vector of delta_smap values
    deltas = d1t(20); // fixed number of moves generated
    int dummy, nr;
    int rem_add[2]; // vector containing removed and added atom
    // mapping
    cg_mapping *mapping = malloc (sizeof(cg_mapping));
    mapping->n_at = Trajectory->n_at;
    mapping->n_cg = cgnum;
    mapping->mapping = i1t(Trajectory->n_at);
    mapping->smap = -1.0;
    mapping->clusters = i1t(Trajectory->frames);
    mapping->size = i1t(Trajectory->frames);
    mapping->idx_cluster = i1t(Trajectory->frames); //ADDED				//(!)
    mapping->omega = i1t(Trajectory->frames);//ADDED					//(!)
    cg_mapping *new_mapping = malloc(sizeof(cg_mapping));
    new_mapping->n_at = Trajectory->n_at;
    new_mapping->n_cg = cgnum;
    new_mapping->smap = -1.0;
    new_mapping->mapping = i1t(Trajectory->n_at);
    new_mapping->clusters = i1t(Trajectory->frames);
    new_mapping->size = i1t(Trajectory->frames);
    new_mapping->idx_cluster = i1t(Trajectory->frames);//ADDED				//(!)
    new_mapping->omega = i1t(Trajectory->frames);//ADDED				//(!)
    int d_idx = 0;
    generate_random_mapping(mapping,f_out_l);
    clustering_and_smap(Trajectory, mapping); //ADDED			//(!)
    for (nr = 0; nr < 20; nr++) {
        my_make_a_move(mapping, new_mapping, rem_add);
        convert_mapping(new_mapping, f_out_l);
        clustering_and_smap(Trajectory, new_mapping); //ADDED			//(!)
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
    free_d1t(deltas);

    //corr_t_zero = 0.7;								//(!!!)
    
    return corr_t_zero;
}

void simulated_annealing(traj *Trajectory, MC_params *SA_params, int cgnum, int rsd, int verbose, int kl_flag, FILE *f_out_l, int core){
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

    char out_filename[1000];   
    /*                                         
    sprintf(out_filename, "./MC_CONVERGENCE/MC_converge_N%d_Ncg%d_MC%d_c%d.csv", Trajectory->n_at, cgnum, SA_params->MC_steps, core);
    FILE* f_smap_mc;
    f_smap_mc = open_file_w(out_filename);
    fprintf(f_smap_mc, "Smap\n");
    */

    sprintf(out_filename, "./MC_RESOLUTION/MC_resolution_N%d_Ncg%d_MC%d_c%d.csv", Trajectory->n_at, cgnum, SA_params->MC_steps, core);	//(!)
    FILE* f_res_mc;
    f_res_mc = open_file_w(out_filename);
    fprintf(f_res_mc, "Resolution\n");

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
    time_t seconds;
    time_t seconds_ref = time(NULL);
    // MC
    double r, p;
    // mapping
    mapping->n_at = Trajectory->n_at;
    mapping->n_cg = cgnum;
    mapping->mapping = i1t(Trajectory->n_at);
    mapping->smap = -1.0;
    mapping->res  = -1.0;
    mapping->clusters = i1t(Trajectory->frames);
    mapping->size = i1t(Trajectory->frames);
    mapping->idx_cluster = i1t(Trajectory->frames); //ADDED				//(!)
    mapping->omega = i1t(Trajectory->frames);//ADDED					//(!)
    new_mapping->n_at = mapping->n_at;
    new_mapping->n_cg = cgnum;
    new_mapping->mapping = i1t(mapping->n_at);
    new_mapping->clusters = i1t(Trajectory->frames);
    new_mapping->size = i1t(Trajectory->frames);
    new_mapping->smap = -1.0;
    new_mapping->res  = -1.0;
    lowest_mapping->res = -1.0;								//(!)
    new_mapping->idx_cluster = i1t(Trajectory->frames);//ADDED				//(!)
    new_mapping->omega = i1t(Trajectory->frames);//ADDED				//(!)
    lowest_mapping->n_at = mapping->n_at;
    lowest_mapping->n_cg = cgnum;
    lowest_mapping->mapping = i1t(mapping->n_at);
    lowest_mapping->clusters = i1t(Trajectory->frames);
    lowest_mapping->size = i1t(Trajectory->frames);
    lowest_mapping->smap = -1.0;
    lowest_mapping->idx_cluster = i1t(Trajectory->frames);//ADDED			//(!)
    lowest_mapping->omega = i1t(Trajectory->frames);//ADDED				//(!)
    printf("allocated variables\n");
    // random mapping and rotation matrices (and first rmsd_mat)
    printf("generating random mapping\n");
    generate_random_mapping(mapping, f_out_l);
    fprintf(f_out_l, "Initial mapping\n");
    convert_mapping(mapping, f_out_l);
    clustering_and_smap(Trajectory, mapping); //ADDED			//(!)
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
        my_make_a_move(mapping, new_mapping, rem_add);
        fprintf(f_out_l, "move made: removed %d added %d\n", rem_add[0], rem_add[1]);
        // observable
        clustering_and_smap(Trajectory, new_mapping); //ADDED			//(!)
        //fprintf(f_out_l, "new_smap %lf\n", new_mapping->smap);
        
        // MC rule
        if ((new_mapping->smap - mapping->smap) < 0) {
            fprintf(f_out_l, "move accepted\n");
            accept_move(rem_add,mapping,new_mapping,Trajectory);
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
                accept_move(rem_add,mapping,new_mapping,Trajectory);
                }
            else {
                fprintf(f_out_l, "move rejected\n");
            }
        }
        fprintf(f_out_l, "SMAP=%lf\n", mapping->smap);
        
        if (mc % 100 == 0 & mc != 0)
        {
            seconds = time(NULL);
            printf("MC step number: %d\n", mc);
            printf("time elapsed: %ld seconds\n", seconds - seconds_ref);
        }

	//fprintf(f_smap_mc, "%lf\n", lowest_mapping->smap);			//(!)
        fprintf(f_res_mc, "%lf\n", mapping->res);                        	//(!)
        //printf("\n RES = %lf\n", mapping->res);                        		//(!)
    }
    fprintf(f_out_l, "END OF SIMULATED ANNEALING OPTIMISATION\n");
    fprintf(f_out_l, "lowest_smap = %lf\n", lowest_mapping->smap);
    fprintf(f_out_l, "lowest_mapping = \n");
    
    // SAVING MINIMUM MAPPING                                            		
    sprintf(out_filename, "./RESULTS/MappingMin_N%d_ncg%d_MC%d_c%d.dat", mapping->n_at, cgnum, SA_params->MC_steps, core);
    FILE* f_smap;
    f_smap = open_file_w(out_filename);
    fprintf(f_smap, "n_at = %d\n", mapping->n_at);
    fprintf(f_smap, "n_cg = %d\n", cgnum);
    fprintf(f_smap, "MC steps = %d\n", SA_params->MC_steps);
    fprintf(f_smap, "decay time = %lf\n", SA_params->decay_time);
    fprintf(f_smap, "lowest_smap = %lf\n", lowest_mapping->smap);
    fprintf(f_smap, "lowest_mapping =\n");                           
    convert_mapping(lowest_mapping, f_smap);
    fprintf(f_smap, "\nexecution_time = %ld\n", seconds - seconds_ref); 
    
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
