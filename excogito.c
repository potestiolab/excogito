/*! \mainpage EXCOGITO
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <my_malloc.h>
#include <time.h>
#include <omp.h>
#include <sys/types.h>
#include <unistd.h>
#include <limits.h>
#include <stdbool.h>
// include called headers
#include <sampling.h>
#include <io.h>
#include <geometry.h>
#include <mapping.h>
#include <observables.h>
#include <alignment.h>
#include <ini.h>
#include <traj.h>
#include <optimize.h>
#include <random_mapping.h>
#include <measure.h>
#include <norm.h>
#include <cosine.h>
#include <distance.h>
#include <argp.h>
#include <measure_kl.h>
#include <optimize_kl.h>
#include <random_mapping_kl.h>
#include <optimize_spins.h>


int main(int argc, char *argv[]) {
    /**
    * main file of the program
    */
    time_t seconds;
    time_t seconds_ref = time(NULL);

    char * tasks [] = {"optimize", "optimize_kl", "random","random_kl","measure","measure_kl", "norm", "cosine", "distance", "optimize_spins"};
    int n_tasks = sizeof(tasks)/sizeof(tasks[0]);

    arguments *args = malloc(sizeof(arguments));

    args->task             = "-1";
    args->silent           = 0   ;
    args->verbose          = 0   ;
    args->parameter_file   = "-1";
    args->energy_file      = "-1";
    args->mapping_file     = "-1";
    args->mapping_file2    = "-1"; 
    args->trajectory_file  = "-1";
    args->prot_code        = "-1"; 
    args->mapping_matrix   = "-1"; 
    args->probability_file = "-1"; 

    FILE *fe; // declaration of error file 

    /* PARSING ARGUMENTS: every option seen by parse_opt will be reflected in args */
    if(argc == 1){
        print_usage_main(argv);
        exit(EXIT_FAILURE); 
    }

    int tk = 0;
    bool found = false;
    for (tk = 0 ; tk < n_tasks; tk++){
        if (strcmp(argv[1], tasks[tk]) == 0){
            found = true;
        }
    }
    if (found == false){
        fe = fopen("error.dat", "w");
        if(strcmp(argv[1], "--usage") == 0){
            print_usage_main(argv);
            exit(EXIT_FAILURE); 
        }
        if(strcmp(argv[1], "--help") == 0 || strcmp(argv[1], "-h") == 0){
            print_help_main(argv);
            exit(EXIT_FAILURE);  
        }
        if(strcmp(argv[1], "--verbose") == 0 || strcmp(argv[1], "-v") == 0){
            printf("Error. %s cannot be the first argument. %s not in the list of accepted tasks\n", argv[1], argv[1]);
            print_usage_main(argv);
            exit(EXIT_FAILURE);
        }  
        fprintf(fe, "Error. %s not in the list of accepted tasks\n", argv[1]);
        printf("Error. %s not in the list of accepted tasks\n", argv[1]);
        fclose(fe); 
        print_usage_main(argv); 
        exit(EXIT_FAILURE);
    }

    check_argv_errors(argv, argc);

    argp_parse (&argp, argc, argv, ARGP_IN_ORDER, 0, args);

    /* Check if the mandatory files are present. Mapping file is optional */
    mandatory_files_present(args, argv);

    srand(time(NULL) + getpid()); // randomizing srand
    printf("parameter file is %s\n", args->parameter_file);   
    printf("trajectory file is %s\n", args->trajectory_file); 
    printf("task is %s\n", args->task); 
    printf("verbose = %d\n", args->verbose);
    printf("Setting Parameters\n");

    parameters *cc = malloc(sizeof(parameters));
    
    read_ParameterFile(args, cc);
    
    // subprograms
    if (strcmp(args->task , "optimize") == 0){
        optimize(args, cc);
    }
    if (strcmp(args->task , "optimize_kl") == 0){
        optimize_kl(args, cc);
    }
    if (strcmp(args->task , "optimize_spins") == 0){
        optimize_spins(args, cc);
    }
    if (strcmp(args->task , "random") == 0){
        random_mapping(args, cc);
    }
    if (strcmp(args->task , "random_kl") == 0){
        random_mapping_kl(args, cc);
    }
    if (strcmp(args->task , "measure") == 0){
        measure(args, cc);
    }
    if (strcmp(args->task , "measure_kl") == 0){
        measure_kl(args, cc);
    }
    if (strcmp(args->task , "norm") == 0){
        norm(args, cc);
    }
    if (strcmp(args->task , "cosine") == 0){
        cosine(args, cc);
    }
    if (strcmp(args->task , "distance") == 0){
        distance(args, cc);
    }
    free(args);
    free(cc);
    seconds = time(NULL);
    printf("\nOverall execution time: %ld seconds", seconds-seconds_ref);
    return 0;
}
