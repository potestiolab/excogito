#include <io.h>
#include <stdio.h>
#include <stdlib.h>
#include <my_malloc.h>
#include <traj.h>
#include <hierarchical_clustering.h>
#include <alignment.h>
#include <mapping.h>
#include <sampling.h>
#include <math.h>
#include <omp.h>
#include <stdbool.h>

int spin_clustering(int **traj_coords, int *mapping, int atomnum, int frames, int *mapping_clusters, int cgnum){
    /**
     * spin clustering
     * 
     */

    int i, j,k;
    int confs[frames][cgnum];
    int conf_pops[frames];
    for (j = 0; j<frames; j++){
        conf_pops[j] = 0;
    }
    int cg_idx;
    
    int cg_conf_num = 0;
    int cg_conf_idx = 0;
    int cg_conf[cgnum];
    char char_cg_conf[cgnum];
    for (j = 0; j<frames; j++){
        cg_idx = 0;
        for (k = 0; k<atomnum; k++){
            if (mapping[k] == 1){
                cg_conf[cg_idx] = traj_coords[j][k];
                cg_idx += 1;
            }
        }
        // looping over previous CG confs
        bool equals = 1;
        for (i = 0; i<cg_conf_num; i++){
            equals = 0;
            for (k = 0; k<cgnum; k++){
                if (confs[i][k] != cg_conf[k]){
                    equals = 1;
                    break;
                }
            }
            if (equals == 0){
                // found configuration
                cg_conf_idx = i;
                conf_pops[i] += 1;
                break;
            }
        }
        if (equals == 1){
            // new configuration
            for (k = 0; k<cgnum; k++){
                confs[cg_conf_num][k] = cg_conf[k];
                conf_pops[cg_conf_num] = 1;
                cg_conf_idx = cg_conf_num;
            }
            cg_conf_num += 1;
        }
        // filling mapping->clusters
        mapping_clusters[j] = cg_conf_idx + 1;
    }
    printf("returning cg_conf_num = %d\n", cg_conf_num);
    return cg_conf_num;
}

void optimize_spins(arguments *arguments, parameters *cc){
    /**
    * subprogram to optimize the coarse-grained mapping by minimising the KL divergence version of its mapping entropy
    * 
    * Parameters
    * ----------
    * 
    * `arguments` : `arguments` object (command line arguments)
    * 
    * `parameters` : `parameters` object
    */
    FILE *fe; // declaration of error file 
    printf("subprogram optimize_spins\n");  
    char out_filename[1000];
    // spin trajectory
    spin_traj *Trajectory = malloc (sizeof(spin_traj));
    Trajectory->frames = cc->frames;
    Trajectory->n_at = cc->atomnum;
    Trajectory->traj_coords = i2t(cc->frames, cc->atomnum);
    Trajectory->energies = d1t(cc->frames);
    Trajectory->stride = cc->stride;
    // OMP calculation
    int nthreads, tid;
    // read files
    printf("reading trajectory\n");
    read_spinFile(arguments->trajectory_file, Trajectory);
    printf("reading probabilities\n");
    read_EnergyFile(arguments->probability_file, Trajectory);
    int i;
    for (i = 0; i < Trajectory->frames; i++){printf("energies %d = %lf\n", i, Trajectory->energies[i]);}
    int isprob = check_probabilities(Trajectory->energies, Trajectory->frames);
    if (isprob != 1){
        fe = fopen("error.dat", "w");
        fprintf(fe, "Error. The probability file %s is not normalized to 1\nAborting.\n", arguments->probability_file);
        printf("Error. The probability file %s is not normalized to 1\nAborting.\n", arguments->probability_file);
        fclose(fe);
        exit(EXIT_FAILURE);
    }
    #pragma omp parallel private(tid)
        {
            tid = omp_get_thread_num();
            if (tid == 0) {
                nthreads = omp_get_num_threads();
                printf("Number of available threads: %d\n", nthreads);
            }
        }
        if (cc->Ncores > nthreads) {
            fe = fopen("error.dat", "w");
            fprintf(fe, "Error. Larger number of threads set (%d) than available (%d)\n", cc->Ncores, nthreads);
            printf("Error. Larger number of threads set (%d) than available (%d)\n", cc->Ncores, nthreads);
            fclose(fe);
            exit(EXIT_FAILURE);
        } else {
            nthreads = cc->Ncores;
            printf("Number of used threads: %d\n", nthreads);
        }
        int q;
        if (cc->Flag_t_zero != 1){
            // if temperature is not defined
            double *tzeros; // vector of delta_smap values
            tzeros = d1t(cc->Ncores);
            printf("estimating simulated annealing start temperature\n");
            #pragma omp parallel shared(Trajectory) private(q, out_filename)
            {
            #pragma omp for schedule(static, 1)
                for (q = 0; q < nthreads; q++) {
                    sprintf(out_filename, "./%skls_fast_delta_N%d_%d.dat", arguments->prot_code, cc->cgnum, q); //
                    FILE *f_out_l;
                    f_out_l = open_file_w(out_filename);
                    printf("launching tzero_estimation\n");
                    tzeros[q] = tzero_estimation_spins(Trajectory, cc->cgnum, f_out_l);
                    fprintf(f_out_l, "t_zero[%d] estimation concluded: starting temperature for simulated annealing = %8.6lf", q, tzeros[q]);
                    fclose(f_out_l);
                }
            #pragma omp barrier
            }
            // calculating mean t_zero
            cc->t_zero = 0.0;
            for (q = 0; q < nthreads; q++) {cc->t_zero += tzeros[q];}
            cc->t_zero = cc->t_zero/nthreads;
            free_d1t(tzeros);
        }
        printf("t_zero %lf decay_time %lf\n", cc->t_zero, cc->decay_time);
        #pragma omp parallel shared(Trajectory) private(q, out_filename)
        {
        #pragma omp for schedule(static, 1)
            for (q = 0; q < nthreads; q++) {
                sprintf(out_filename, "./%skls_SA_N%d_%d.dat", arguments->prot_code, cc->cgnum, q); //
                FILE *f_out_l;
                f_out_l = open_file_w(out_filename);
                MC_params *SA_params = malloc (sizeof(MC_params));
                SA_params->t_zero = cc->t_zero;
                SA_params->MC_steps = cc->MC_steps;
                if (cc->Flag_decay_time != 1){
                    double d_decay_time = -SA_params->MC_steps*1.0/log(0.001); // so that to have 1/1000 of the start temperature at the end. TODO: convert to float
                    fprintf(f_out_l,"calculating decay parameter for simulated annealing\nUsing decay parameter = %8.6lf\n", d_decay_time);
                    SA_params->decay_time = d_decay_time;
                }
                else{SA_params->decay_time = cc->decay_time;}
                SA_params->rotmats_period = cc->rotmats_period;
                printf("rsd = %d\n", cc->rsd);
                simulated_annealing_spins(Trajectory, SA_params, cc->cgnum, arguments->verbose, f_out_l);
                fclose(f_out_l);
                // free SA_params
                free(SA_params);
            }
        #pragma omp barrier
                }
        printf("finished pragma part\n");
        printf("finished mapping optimisation");
        // free trajectory
        free_i2t(Trajectory->traj_coords);
        free_d1t(Trajectory->energies);
        free(Trajectory);
}
