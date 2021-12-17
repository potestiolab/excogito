#include <io.h>
#include <stdlib.h>
#include <my_malloc.h>
#include <traj.h>
#include <mapping.h>
#include <math.h>
#include <geometry.h>
#include <observables.h>

void norm(arguments *arguments, parameters *cc){
    /**
    * subprogram to To calculate the norm of a mapping (provided by the user) throughout a trajectory
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
        fprintf(fe, "Error. Task norm requires Ncores == 1. here Ncores == %d\nAborting.\n", cc->Ncores);
        printf("Error. Task norm requires Ncores == 1. here Ncores == %d\nAborting.\n", cc->Ncores);
        fclose(fe);
        exit(EXIT_FAILURE);
    }
    char out_filename[1000];
    traj *Trajectory = malloc (sizeof(traj));
    Trajectory->frames = cc->frames;
    Trajectory->n_at = cc->atomnum;
    Trajectory->traj_coords = d2t(cc->frames, 3 * cc->atomnum);
    //Trajectory->energies = d1t(cc->frames);
    Trajectory->pairs = cc->frames * (cc->frames - 1) / 2;
    printf("frames = %d\n", Trajectory->frames);
    // performance measurement
    //time_t seconds;
    //time_t seconds_ref = time(NULL);
    // computing parameters and allocating objects
    printf("overall pairs = %d\n", Trajectory->pairs);
    // read trajectory
    printf("reading trajectory\n");
    read_TrajectoryFile(arguments->trajectory_file, Trajectory);
    cg_mapping *mapping = malloc (sizeof(cg_mapping));
    //mapping = malloc(sizeof(cg_mapping));
    mapping->n_at = cc->atomnum;
    mapping->n_cg = cc->cgnum;
    mapping->mapping = i1t(cc->atomnum);
    zero_vec_i(mapping->mapping,mapping->n_at);
    sprintf(out_filename, "./%s_norms_N%d.dat", arguments->prot_code, cc->cgnum);
    FILE *f_out_l;
    f_out_l = open_file_w(out_filename);
    fprintf(f_out_l,"computing norms for provided mapping\n");
    read_MappingFile(arguments->mapping_file, f_out_l, mapping);
    compute_mapping_norms(Trajectory,mapping,f_out_l);
    fprintf(f_out_l, "finished computation of mapping norms");
    fclose(f_out_l);
    // free mapping
    free_i1t(mapping->mapping);
    free(mapping);
    // free trajectory
    free_d2t(Trajectory->traj_coords);
    //free_d1t(Trajectory->energies);
    free(Trajectory);
}
