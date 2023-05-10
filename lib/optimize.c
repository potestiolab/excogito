#include <io.h>
#include <stdlib.h>
#include <my_malloc.h>
#include <traj.h>
#include <hierarchical_clustering.h>
#include <alignment.h>
#include <mapping.h>
#include <sampling.h>
#include <math.h>
#include <omp.h>

void optimize(arguments *arguments, parameters *cc){
    /**
    * subprogram to optimize the coarse-grained mapping by minimising its mapping entropy
    * 
    * Parameters
    * ----------
    * 
    * `arguments` : `arguments` object (command line arguments)
    * 
    * `parameters` : `parameters` object
    */ 
    FILE *fe; // declaration of error file 
    printf("subprogram optimize\n");  
    char out_filename[1000];
    int spins = 0;
    // clustering
    clust_params *clustering = malloc (sizeof(clust_params));
    clustering->crit = cc->criterion;
    clustering->ncl = cc->nclust;
    clustering->max_ncl = cc->max_nclust;
    clustering->min_ncl = cc->min_nclust;
    clustering->c_distance = cc->distance;
    traj *Trajectory = malloc (sizeof(traj));
    Trajectory->frames = cc->frames;
    Trajectory->n_at = cc->atomnum;
    Trajectory->traj_coords = d2t(cc->frames, 3 * cc->atomnum);
    Trajectory->energies = d1t(cc->frames);
    Trajectory->stride = cc->stride;
    if (clustering->crit == 3){
        printf("criterion = 3\n");
        printf("cc->frames/cc->stride = %d\n",cc->frames/cc->stride);
        if ((Trajectory->frames-1)%Trajectory->stride == 0){
            Trajectory->eff_frames = cc->frames/cc->stride + 1;
        }
        else{
            Trajectory->eff_frames = (cc->frames-1)/cc->stride + 2;
        }
        printf("effective frames = %d\n", Trajectory->eff_frames);
        Trajectory->pairs = Trajectory->eff_frames * (Trajectory->eff_frames - 1)/ 2;
        Trajectory->strides = i1t(cc->frames);
        int idx = 0;
        for (idx = 0; idx < cc->frames; idx++){
            if (idx%cc->stride == 0) {Trajectory->strides[idx] = 1;}
            else{Trajectory->strides[idx] = 0;}
        }
        Trajectory->strides[cc->frames-1] = 1;
        printf("strides\n");
        for (idx = 0; idx < cc->frames; idx++){printf("%d ", Trajectory->strides[idx]);}
    }
    else{Trajectory->pairs = cc->frames * (cc->frames - 1) / 2;}
    printf("frames = %d\n", Trajectory->frames);
    // computing parameters and allocating objects
    printf("overall pairs = %d\n", Trajectory->pairs);
    // OMP calculation
    int nthreads, tid;
    // read files
    printf("reading trajectory\n");
    read_TrajectoryFile(arguments->trajectory_file, Trajectory, spins);                                                             //(!) MODIFIED
    printf("reading energies\n");
    read_EnergyFile(arguments->energy_file, Trajectory);
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
            #pragma omp parallel shared(Trajectory, clustering) private(q, out_filename)
            {
            #pragma omp for schedule(static, 1)
                for (q = 0; q < nthreads; q++) {
                    sprintf(out_filename, "./%s_fast_delta_N%d_%d.dat", arguments->prot_code, cc->cgnum, q); //
                    FILE *f_out_l;
                    f_out_l = open_file_w(out_filename);
                    //tzeros[q] = tzero_estimation(Trajectory, clustering, cc->cgnum, cc->rsd, arguments->verbose, 0, f_out_l);     //(!) ELIMINATE because we change it for optimize_kl
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
        #pragma omp parallel shared(Trajectory, clustering) private(q, out_filename)
        {
#pragma omp for schedule(static, 1)
            for (q = 0; q < nthreads; q++) {
                sprintf(out_filename, "./%s_SA_N%d_%d.dat", arguments->prot_code, cc->cgnum, q); //
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
                //simulated_annealing(Trajectory, clustering, SA_params, cc->cgnum, cc->rsd, arguments->verbose, 0, f_out_l);   //(!) ELIMINATE because we change it for optimize_kl
                fclose(f_out_l);
                // free SA_params
                free(SA_params);
            }
#pragma omp barrier
        }
        printf("finished pragma part\n");
        printf("finished mapping optimisation");
        // free trajectory
        free_d2t(Trajectory->traj_coords);
        free_d1t(Trajectory->energies);
        if (clustering->crit == 3){free_i1t(Trajectory->strides);}
        free(Trajectory);
        // free clustering
        free(clustering);
}
