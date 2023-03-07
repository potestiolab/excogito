#include <io.h>
#include <stdlib.h>
#include <my_malloc.h>
#include <traj.h>
#include <hierarchical_clustering.h>
#include <alignment.h>
#include <mapping.h>
#include <sampling.h>
#include <math.h>
#include <geometry.h>
#include <observables.h>

void random_mapping(arguments *arguments, parameters *cc){
    /**
    * subprogram to randomly generate coarse-grained representations and measure the associated mapping entropies
    * 
    * Parameters
    * ----------
    * 
    * `arguments` : `arguments` object (command line arguments)
    * 
    * `parameters` : `parameters` object
    */
    FILE *fe; // declaration of error file 
    printf("subprogram random\n");   
    if (cc->Ncores != 1) { // todo: parallelize
        fe = fopen("error.dat", "w");
        fprintf(fe, "Error. Task random requires Ncores == 1. here Ncores == %d\nAborting.\n", cc->Ncores);
        printf("Error. Task random requires Ncores == 1. here Ncores == %d\nAborting.\n", cc->Ncores);
        fclose(fe);
        exit(EXIT_FAILURE);
    }
    char out_filename[1000];
    // clustering
    clust_params *clustering = malloc (sizeof(clust_params));
    clustering->crit = cc->criterion;
    clustering->ncl = cc->nclust;
    clustering->max_ncl = cc->max_nclust;
    clustering->min_ncl = cc->min_nclust;
    clustering->c_distance = cc->distance;
    // trajectory
    traj *Trajectory = malloc (sizeof(traj));
    Trajectory->frames = cc->frames;
    Trajectory->n_at = cc->atomnum;
    Trajectory->traj_coords = d2t(cc->frames, 3 * cc->atomnum);
    Trajectory->energies = d1t(cc->frames);
    Trajectory->stride = cc->stride;
    if (clustering->crit == 1){
        printf("criterion = %d\n", clustering->crit);
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
    // trajectory was the last object to be declared "in common"
    printf("reading energies\n");
    // read trajectory
    printf("reading trajectory\n");
    read_TrajectoryFile(arguments->trajectory_file, Trajectory);
    read_EnergyFile(arguments->energy_file, Trajectory);
    // declaring variables
    cg_mapping *mapping = malloc (sizeof(cg_mapping));
    mapping->n_at = cc->atomnum;
    mapping->n_cg = cc->cgnum;
    mapping->mapping = i1t(mapping->n_at);
    zero_vec_i(mapping->mapping,mapping->n_at);
    mapping->smap = -1.0;
    mapping->clusters = i1t(Trajectory->frames);
    mapping->size = i1t(Trajectory->frames);
    alignments *align = malloc (sizeof(alignments));
    align->rmsd_mat = d1t(Trajectory->pairs);
    align->rotation_matrices = d2t(Trajectory->pairs, 9);
    align->coms = d2t(cc->frames, 3);
    align->rsd = cc->rsd;
    if (clustering->crit == 1){
        align->rmsd_vector = d1t(Trajectory->frames * 2);
        align->rotation_matrices_vector = d2t(Trajectory->frames * 2, 9);
    }
    int rd_map;
    sprintf(out_filename, "./%s_random_N%d.dat", arguments->prot_code, cc->cgnum);
    FILE *f_out_l;
    f_out_l = open_file_w(out_filename);
    for (rd_map = 0; rd_map < cc->n_mappings; rd_map++) {
        fprintf(f_out_l, "generating random mapping number %d\n", rd_map);
        generate_random_mapping(mapping, f_out_l);
        fprintf(f_out_l, "computing rotmats\n");
        if (clustering->crit==1){
            cycle_alignment_fastclust(Trajectory, align, mapping);            
        }
        else{cycle_alignment(Trajectory, align, mapping);}
        fprintf(f_out_l, "computing observable\n");
        overall_compute_smap(align, clustering, Trajectory, mapping, arguments->verbose, 0);
        fprintf(f_out_l, "random_smap = %8.6lf\n", mapping->smap);
    }
    fclose(f_out_l);
    printf("finished random mappings generation");
    // free mapping
    free_mapping(mapping);
    // free alignment
    free_d1t(align->rmsd_mat);
    free_d2t(align->rotation_matrices);
    free_d2t(align->coms);
    if (clustering->crit == 1){
        free_d1t(align->rmsd_vector);
        free_d2t(align->rotation_matrices_vector);
    }
    free(align);
    // free trajectory
    free_d2t(Trajectory->traj_coords);
    free_d1t(Trajectory->energies);
    if (clustering->crit == 1){free_i1t(Trajectory->strides);}
    free(Trajectory);
    // free clustering
    free(clustering);
}
