#include <io.h>
#include <stdlib.h>
#include <my_malloc.h>
#include <traj.h>
#include <mapping.h>
#include <math.h>
#include <geometry.h>
#include <observables.h>

void distance(arguments *arguments, parameters *cc){
    /**
    * subprogram to calculate the distance matrix between a data set of mappings (provided by the user) over a single conformation
    * 
    * Parameters
    * ----------
    * 
    * `arguments` : `arguments` object (command line arguments)
    * 
    * `parameters` : `parameters` object
    */ 
    FILE *fe; // declaration of error file 
    printf("subprogram norm\n");
    if (cc->Ncores != 1) { // todo: parallelize
        fe = fopen("error.dat", "w");
        fprintf(fe, "Error. Task distance requires Ncores == 1. here Ncores == %d\nAborting.\n", cc->Ncores);
        printf("Error. Task distance requires Ncores == 1. here Ncores == %d\nAborting.\n", cc->Ncores);
        fclose(fe);
        exit(EXIT_FAILURE);
    }
    char out_filename[200];
    int spins = 0;
    traj *Trajectory = malloc (sizeof(traj));
    Trajectory->frames = cc->frames;
    Trajectory->n_at = cc->atomnum;
    Trajectory->traj_coords = d2t(cc->frames, 3 * cc->atomnum);
    Trajectory->pairs = cc->frames * (cc->frames - 1) / 2;
    printf("frames = %d\n", Trajectory->frames);
    // computing parameters and allocating objects
    printf("overall pairs = %d\n", Trajectory->pairs);
    // read trajectory
    printf("reading trajectory\n");
    read_TrajectoryFile(arguments->trajectory_file, Trajectory, spins);
    printf("distance matrix calculation for %d mappings\n", cc->n_mappings);
    sprintf(out_filename, "./%s_distmat_N%d.dat", arguments->prot_code, cc->cgnum);
    FILE *f_out_l;
    f_out_l = open_file_w(out_filename);
    fprintf(f_out_l,"computing distance matrix for a set of mappings over a configuration\n");
    if (Trajectory->frames != 1){
        printf("number of frames %d != 1: not valid for task distance\nAborting.\n", Trajectory->frames);
        fe = open_file_w("error.dat");
        fprintf(fe, "number of frames %d != 1: not valid for task distance\nAborting.\n", Trajectory->frames);
        fclose(fe);
        exit(EXIT_FAILURE);
    }
    char distmat_filename[200];
    sprintf(distmat_filename, "./%s_distmat_N%d.csv", arguments->prot_code, cc->cgnum);
    // declaring mapping matrix
    int map_idx;
    cg_mapping *mapping_matrix[cc->n_mappings];
    for (map_idx = 0; map_idx < cc->n_mappings; map_idx++){
        printf("allocating struct %d\n", map_idx);
        mapping_matrix[map_idx] = malloc(sizeof(cg_mapping));
        mapping_matrix[map_idx]->n_at = cc->atomnum;
        mapping_matrix[map_idx]->n_cg = cc->cgnum;
        mapping_matrix[map_idx]->mapping = i1t(cc->atomnum);
        mapping_matrix[map_idx]->norms = d1t(Trajectory->frames); // only one value
        zero_vec_i(mapping_matrix[map_idx]->mapping,mapping_matrix[map_idx]->n_at);
    }
    // loading mappings
    printf("loading mapping matrix\n");
    ////load_mapping_matrix(arguments->mapping_matrix, f_out_l, mapping_matrix, cc->n_mappings);  //OLD MARCO
    read_mapping_matrix(arguments->mapping_matrix, f_out_l, mapping_matrix, cc->n_mappings);
    for (map_idx = 0; map_idx < cc->n_mappings; map_idx++){
        convert_mapping(mapping_matrix[map_idx],f_out_l);
    }
    compute_mapping_distmat(Trajectory,mapping_matrix,cc->n_mappings,f_out_l,distmat_filename); 
}
