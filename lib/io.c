/**
* \class io
* \brief library of functions for all input-output operations
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <io.h>
#include <argp.h>
#include <ini.h>
#include <traj.h>
#include <mapping.h>

//const char *argp_program_version     = "MeTool Program";

const char *argp_program_bug_address = "<raffaele.fiorentini@unitn.it>"; 

void check_int_string(const char *str, int row, char *fname){
    /**
    * routine that checks if the string token in account reading a generic FILE is an INTEGER number
    *
    * Parameters
    * ----------
    *
    * `str`   : string token in account
    *  
    * `row`   : number of row where the string is found. 
    *
    * `fname` : filename read 
    */
 
    FILE *fe; 
    int j;

    for(j=0; j<strlen(str); j++){    // Reading each character of the string 
   
        if(isalpha(str[j])){
            fe = fopen("error.dat", "w");
            fprintf(fe, "Error. The line must be integer. In %s, the character '%c' is present in '%s' at %dth row\n", fname, str[j], str, row+1);
            printf("Error. The line must be integer. In %s, the character '%c' is present in '%s' at %dth row\n", fname, str[j], str, row+1);
            fclose(fe); 
            exit(EXIT_FAILURE);
        }

        else if(str[j]>= '0' && str[j] <= '9')
            continue;

        else if(str[j] == ' ' || str[j] =='\n' || str[j] =='\t' || str[j] =='\v' || str[j] =='\r')
            continue;

        else if(str[j] =='.'){
            fe = fopen("error.dat", "w");
            fprintf(fe, "Error. The line must be integer. In %s, DOT is present in %s at %dth row\n", fname, str, row+1);
            printf("Error. The line must be integer. In %s, DOT is present in %s at %dth row\n", fname, str, row+1);
            fclose(fe); 
            exit(EXIT_FAILURE);
        }

        else{
            fe = fopen("error.dat", "w");
            fprintf(fe, "Error. The line must be integer. In %s, the special character %c is present in %s at %dth row\n", fname, str[j], str, row+1);
            printf("Error. The line must be integer. In %s, the special character %c is present in %s at %dth row\n", fname, str[j], str, row+1);
            fclose(fe);
            exit(EXIT_FAILURE);
        }
    }
}

void check_int_string_iniFile(const char *str, char *fname, char *name){
    /**
    * routine that checks if the string token in account is an INTEGER number. It works only for ini Files 
    *
    * Parameters
    * ----------
    *
    * `str`   : string token in account
    *  
    * `fname` : parameter filename  
    *
    * `name`  : name of each parameter in the file 
    */ 

    FILE *fe; 
    int j;

    for(j=0; j<strlen(str); j++){    // Reading each character of the string 
   
        if(isalpha(str[j])){
            fe = fopen("error.dat", "w");
            fprintf(fe,"Error. The line must be integer. In %s, the character '%c' is present in '%s'(%s = %s)\n", fname, str[j], str, name, str);
            printf("Error. The line must be integer. In %s, the character '%c' is present in '%s'(%s = %s)\n", fname, str[j], str, name, str);
            fclose(fe); 
            exit(EXIT_FAILURE);
        }

        else if(str[j]>= '0' && str[j] <= '9')
            continue;

        else if(str[j] == ' ' || str[j] =='\n' || str[j] =='\t' || str[j] =='\v' || str[j] =='\r')
            continue;

        else if(str[j] =='.'){
            fe = fopen("error.dat", "w");
            fprintf(fe,"Error. The line must be integer. In %s, DOT is present in %s (%s = %s)\n", fname, str, name, str);
            printf("Error. The line must be integer. In %s, DOT is present in %s (%s = %s)\n", fname, str, name, str);
            fclose(fe); 
            exit(EXIT_FAILURE);
        }

        else{
            fe = fopen("error.dat", "w");
            fprintf(fe,"Error. The line must be integer. In %s, the special character %c is present in %s (%s = %s)\n", fname, str[j], str, name, str);
            printf("Error. The line must be integer. In %s, the special character %c is present in %s (%s = %s)\n", fname, str[j], str, name, str);
            fclose(fe);
            exit(EXIT_FAILURE);
        }
    }
}

/* INI FILE ...... ... . . .*/

int handler(void* config, const char* section, const char* name, const char* value){

    // Config instance for filling in the values.
    parameters* pconfig = (parameters*)config;
    // Define a macro for checking Sections and keys under the sections.
    #define MATCH(s, n) strcmp(section, s) == 0 && strcmp(name, n) == 0

    // Check that each VALUE is an INTEGER and, in that case, associate the value at each pconfig
    if(MATCH("Parameters", "atomnum")){
        check_int_string_iniFile(value, "parameter file", "atomnum"); 
        pconfig->atomnum = atoi(value);
        pconfig->Flag_atomnum = 1; 
    }
    else if(MATCH("Parameters", "frames")){
        check_int_string_iniFile(value, "parameter file", "frames");
        pconfig->frames = atoi(value);
        pconfig->Flag_frames = 1;
    } 
    else if(MATCH("Parameters", "cgnum")){
        check_int_string_iniFile(value, "parameter file", "cgnum");
        pconfig->cgnum = atoi(value);  
        pconfig->Flag_cgnum = 1;     
    }
    else if(MATCH("Parameters", "nclust")){
        check_int_string_iniFile(value, "parameter file", "nclust");
        pconfig->nclust = atoi(value);
        pconfig->Flag_nclust = 1;
    }
    else if(MATCH("Parameters", "n_mappings")){
        check_int_string_iniFile(value, "parameter file", "n_mappings");
        pconfig->n_mappings = atoi(value); 
        pconfig->Flag_n_mappings = 1; 
    }
    else if(MATCH("Parameters", "MC_steps")){
        check_int_string_iniFile(value, "parameter file", "MC_steps");
        pconfig->MC_steps = atoi(value);
        pconfig->Flag_MC_steps = 1; 
    }
    else if(MATCH("Parameters", "rotmats_period")){
        check_int_string_iniFile(value, "parameter file", "rotmats_period");
        pconfig->rotmats_period = atoi(value); 
        pconfig->Flag_rotmats_period = 1; 
    }
    else if(MATCH("Parameters", "t_zero")){
        check_float_string_iniFile(value, "parameter file", "t_zero");
        pconfig->t_zero = atof(value);
        pconfig->Flag_t_zero = 1; 
    }
    else if(MATCH("Parameters", "distance")){
        check_float_string_iniFile(value, "parameter file", "distance");
        pconfig->distance = atof(value);
        pconfig->Flag_distance = 1; 
    }
    else if(MATCH("Parameters", "criterion")){
        check_int_string_iniFile(value, "parameter file", "criterion");
        pconfig->criterion = atoi(value);
        pconfig->Flag_criterion = 1; 
    }
    else if(MATCH("Parameters", "max_nclust")){
        check_int_string_iniFile(value, "parameter file", "max_nclust");
        pconfig->max_nclust = atoi(value);
        pconfig->Flag_max_nclust = 1;
    }
    else if(MATCH("Parameters", "min_nclust")){ 
        check_int_string_iniFile(value, "parameter file", "min_nclust");
        pconfig->min_nclust = atoi(value);
        pconfig->Flag_min_nclust = 1; 
    }
    else if(MATCH("Parameters", "Ncores")){
        check_int_string_iniFile(value, "parameter file", "Ncores");
        pconfig->Ncores = atoi(value);
        pconfig->Flag_Ncores = 1; 
    }
    else if(MATCH("Parameters", "decay_time")){
        check_float_string_iniFile(value, "parameter file", "decay_time");
        pconfig->decay_time = atof(value);
        pconfig->Flag_decay_time = 1;
    }
    else if(MATCH("Parameters", "rsd")){
        check_int_string_iniFile(value, "parameter file", "rsd");
        pconfig->rsd = atoi(value);
        pconfig->Flag_rsd = 1; 
    }
    else if(MATCH("Parameters", "stride")){
        check_int_string_iniFile(value, "parameter file", "stride");
        pconfig->stride = atoi(value);
        pconfig->Flag_stride = 1; 
    }
    else
        return 0;

    return 1; 
}

/* Function parameters pp_config that reads each value of ini file... */
parameters pp_config(parameters config){

int atomnum             = config.atomnum;
int frames              = config.frames;
int cgnum               = config.cgnum;
int nclust              = config.nclust;
int n_mappings          = config.n_mappings;
int MC_steps            = config.MC_steps;
int rotmats_period      = config.rotmats_period;
float t_zero            = config.t_zero;
float distance          = config.distance;
int criterion           = config.criterion;
int max_nclust          = config.max_nclust;
int min_nclust          = config.min_nclust;
int Ncores              = config.Ncores;
int decay_time          = config.decay_time;
int rsd                 = config.rsd;
int stride              = config.stride;

// Flags of each variable
int Flag_atomnum           = config.Flag_atomnum;
int Flag_frames 	       = config.Flag_frames; 	     
int Flag_cgnum 	           = config.Flag_cgnum; 	          
int Flag_nclust 	       = config.Flag_nclust; 	     
int Flag_n_mappings        = config.Flag_n_mappings;
int Flag_MC_steps	       = config.Flag_MC_steps;	  
int Flag_rotmats_period    = config.Flag_rotmats_period;   
int Flag_t_zero 	       = config.Flag_t_zero; 	  
int Flag_distance 	       = config.Flag_distance; 	  
int Flag_criterion         = config.Flag_criterion;        
int Flag_max_nclust        = config.Flag_max_nclust;       
int Flag_min_nclust        = config.Flag_min_nclust;       
int Flag_Ncores 	       = config.Flag_Ncores; 	  
int Flag_decay_time        = config.Flag_decay_time;       
int Flag_rsd               = config.Flag_rsd;  
int Flag_stride            = config.Flag_stride; 
 
return config;
}


void check_empty_file(FILE *f, char *filename){
    /**
    * routine that checks if the file required exists. If it is the case, check if it is empty or not. 
    *
    * Parameters
    * ----------
    *
    * `f`        : FILE structure that represents the file opened.
    *
    * `filename` : filename read  
    */
 
    int size;    

    FILE *fe;

    if(f == NULL){
        fe = fopen("error.dat", "w");
        fprintf(fe, "Error while opening the file. '%s' does not exist\n", filename);
        printf("Error while opening the file. '%s' does not exist\n", filename);
        fclose(fe);
        exit(EXIT_FAILURE);
    }

    fseek(f, 0, SEEK_END);    
    size = ftell(f);
 
    if(size == 0){
        fe = fopen("error.dat", "w");
        fprintf(fe, "Error. Your file is empty. Please, fill out %s with significant data or use another file\n", filename); 
        printf("Error. Your file is empty. Please, fill out %s with significant data or use another file\n", filename); 
        fclose(fe);
        exit(EXIT_FAILURE); 
    }
}

int n_rows (FILE *f){
    /**
    * routine that returns the number of rows in a file. It counts correctly this number even if the last row does not present \n
    *
    * Parameter
    * ----------
    *
    * `f` : FILE structure that represents the file opened.
    */
 
    size_t rows_i = 0, n = 0;
    char *buf = NULL;

    fseek(f, 0, SEEK_SET);

    while(getline(&buf, &n, f) != -1) rows_i++;

    free(buf);

    return rows_i;
}

void check_empty_rows(char *str){
    /**
    * routine that checks if a generic line is empty or not
    *
    * Parameter
    * ----------
    *
    * `str`   : string token in account 
    */

    FILE *fe;

    if(str[strspn(str, " \t\v\r\n")] == '\0'){   
        fe = fopen("error.dat", "w");
        fprintf(fe, "In your file there is, at least, an empty row.\n");
        printf("In your file there is, at least, an empty row.\n");
        exit(EXIT_FAILURE);
    }
}

void check_argv_errors(char *argv[], int argc){
    /** 
    * routine that checks the correctness of command line arguments
    *
    * Parameter
    * ----------
    *
    * `argv[]` : array of command line arguments
    * 
    * `argc`   : number of command line arguments
    */
    FILE *fe;
    int i; 

    for(i=1; i<argc; i++){     // avoiding i = 0 beacuse we know in advance that argv[0] == ./smap_tool

        if(i%2 == 0){          // if even...
 
            if(argv[i][0] == '-' && strlen(argv[i]) == 1){
                fe = fopen("error.dat", "w");
                fprintf(fe, "Error. '-' is not accepted as flag. Use, for example, '-p' instead of '-'. Look below for further help.\n"); 
                printf("Error. '-' flag is not accepted as flag. Use, for example, '-p' instead of '-'. Look below for further help.\n");
                fclose(fe); 
                print_help(argv);
                exit(EXIT_FAILURE); 
            }

            if(argv[i][0] == '-' && argv[i][1] != '-' &&  strlen(argv[i]) == 2)
                continue; 

            if(argv[i][0] == '-' && argv[i][1] != '-' &&  strlen(argv[i]) > 2){
                fe = fopen("error.dat", "w");
                fprintf(fe, "Error. Each flag must contain '-' plus ONLY ONE letter. Example: -p. Look below for further help.\n");
                printf("Error. Each flag must contain '-' plus ONLY ONE letter. Example: -p. Look below for further help.\n");
                fclose(fe); 
                print_help(argv); 
                exit(EXIT_FAILURE);  
            }

            if(argv[i][0] == '-' && argv[i][1] == '-' &&  strlen(argv[i]) == 2){
                fe = fopen("error.dat", "w");
                fprintf(fe, "Error. '--' is not accepted as flag. Use, for example, '--p' instead of '--'. Look below for further help.\n");
                printf("Error. '--' is not accepted as flag. Use, for example, '--p' instead of '--'. Look below for further help.\n");
                fclose(fe); 
                print_help(argv);
                exit(EXIT_FAILURE);
            }

            if(argv[i][0] == '-' && argv[i][1] == '-' &&  strlen(argv[i]) > 2)
                continue;          
                
        }

    }

}

void check_float_string(char *str, int row, char *fname){
    /**
    * routine that checks if the string token in account reading a generic FILE is an Float number
    *
    * Parameters
    * ----------
    *
    * `str`   : string token in account
    *
    * `row`   : number of row where the string is found.
    *
    * `fname` : filename read  
    */
 
    FILE *fe; 
    int j; 
    int count_dots = 0;  
    int count_minus = 0; 

    for(j=0; j<strlen(str); j++){
        if(isalpha(str[j])){
            fe = fopen("error.dat", "w"); 
            fprintf(fe, "Error. Each line must be FLOAT. In %s, the character %c is present in %s at %dth row\n", fname, str[j], str, row+1);
            printf("Error. Each line must be FLOAT. In %s, the character %c is present in %s at %dth row\n", fname, str[j], str, row+1);
            fclose(fe); 
            exit(EXIT_FAILURE);
        }

        else if(str[j]>= '0' && str[j] <= '9')
            continue;

        else if(str[j] == ' ' || str[j] =='\n' || str[j] =='\t' || str[j] =='\v' || str[j] =='\r')
            continue; 

        else if(j ==0 && str[j] == '-' )  // Only the 1st letter of the string can be "-" that indicates that the string is negative number. 
            continue; 

        else if(str[j] == '.'){
            count_dots++; 
            if(count_dots>1){             // FLOAT numbers require that each string must have ONLY one dot.
                fe = fopen("error.dat", "w");
                fprintf(fe, "Error. Each line must be FLOAT. In %s, more than one dot is present in %s at %dth row\n", fname, str, row+1); 
                printf("Error. Each line must be FLOAT. In %s, more than one dot is present in %s at %dth row\n", fname, str, row+1); 
                fclose(fe);
                exit(EXIT_FAILURE);
            }
	}   		
        else{
            fe = fopen("error.dat", "w");
            fprintf(fe, "Error. Each line must be FLOAT. In %s, the special character %c is present in %s at %dth row\n", fname, str[j], str, row+1);
            printf("Error. Each line must be FLOAT. In %s, the special character %c is present in %s at %dth row\n", fname, str[j], str, row+1);
            fclose(fe); 
            exit(EXIT_FAILURE);
        }
    }
    if(count_dots<1){                      // Avoid the possibility of having INTEGER numbers.
        fe = fopen("error.dat", "w");
        fprintf(fe,     "Error. Each line must be FLOAT. In %s, the %dth row contains an integer number--> %s\n", fname, row+1, str);
        printf("Error. Each line must be FLOAT. In %s, the %dth row contains an integer number--> %s\n", fname, row+1, str);
        fclose(fe); 
        exit(EXIT_FAILURE); 
    }   
}

void check_float_string_iniFile(const char *str, char *fname, char *name){
    /**
    * routine that checks if the string token in account is an Float number. It works only for ini Files 
    *
    * Parameters
    * ----------
    *
    * `str`   : string token in account
    *  
    * `fname` : parameter filename  
    *
    * `name`  : name of each parameter in the file 
    */ 

    FILE *fe; 
    int j; 
    int count_dots = 0;  
    int count_minus = 0; 

    for(j=0; j<strlen(str); j++){
     
        if(isalpha(str[j])){
            fe = fopen("error.dat", "w");
            fprintf(fe,"Error. The line must be Float. In %s, the character '%c' is present in '%s'(%s = %s)\n", fname, str[j], str, name, str);
            printf("Error. The line must be Float. In %s, the character '%c' is present in '%s'(%s = %s)\n", fname, str[j], str, name, str); 
            fclose(fe); 
            exit(EXIT_FAILURE);
        }

        else if(str[j]>= '0' && str[j] <= '9')
            continue;

        else if(str[j] == ' ' || str[j] =='\n' || str[j] =='\t' || str[j] =='\v' || str[j] =='\r')
            continue; 

        else if(j ==0 && str[j] == '-' )  // Only the 1st letter of the string can be "-" that indicates that the string is negative number. 
            continue; 

        else if(str[j] == '.'){
            count_dots++; 
            if(count_dots>1){             // FLOAT numbers require that each string must have ONLY one dot.
                fe = fopen("error.dat", "w");
                fprintf(fe,"Error. The line must be Float. In %s, more than one dot is present in %s (%s = %s)\n", fname, str, name, str); 
                printf("Error. The line must be Float. In %s, more than one dot is present in %s (%s = %s)\n", fname, str, name, str); 
                fclose(fe);
                exit(EXIT_FAILURE);
            }
	}   		
        else{
            fe = fopen("error.dat", "w");
            fprintf(fe,"Error. The line must be Float. In %s, the special character %c is present in %s (%s = %s)\n", fname, str[j], str, name, str);
            printf("Error. The line must be Float. In %s, the special character %c is present in %s (%s = %s)\n", fname, str[j], str, name, str);
            fclose(fe); 
            exit(EXIT_FAILURE);
        }
    }
    if(count_dots<1){                      // Avoid the possibility of having INTEGER numbers.
        fe = fopen("error.dat", "w");
        fprintf(fe, "Error. The line must be Float. '%s' contains an integer number (%s = %s) --> %s\n", fname, name, str, str);
        printf("Error. The line must be Float. '%s' contains an integer number (%s = %s) --> %s\n", fname, name, str, str);
        fclose(fe); 
        exit(EXIT_FAILURE); 
    }   
}

int columns(char* string){
    /**
    * routine that returns the number of columns for each row inside the file chosen.
    *
    * Parameter
    * ----------
    *
    * `string`   : string token in account
    */
 
    int counter = 0;

    while(string){
        counter++;
        string = strtok(NULL, " \t\v\r\n");
    }
    return counter; 
}

void print_usage_main(char *argv[]){
    /**
    * routine that prints the usage of the program
    *
    * Parameter
    * ----------
    *
    * `argv[]` : array of command line arguments
    */
    char * tasks [] = {"optimize", "random", "measure", "norm", "cosine", "distance", "optimize_kl", "measure_kl", "random_kl", "optimize_spins"};
    int n_tasks = sizeof(tasks)/sizeof(tasks[0]);

    int tk = 0;
    
    printf("\n"); 
    for(tk = 0 ; tk < n_tasks; tk++){
        printf("Usage: %s %s [OPTIONS]\n", argv[0], tasks[tk]);
    }

    printf("\nTry %s -h or %s --help for more information.\n\n", argv[0], argv[0]);
}

void print_help_main(char *argv[]){
    /**
    * routine that prints detailed information about the program
    *
    * Parameter
    * ----------
    *
    * `argv[]` : array of command line arguments
    */
    char * tasks [] = {"optimize", "random", "measure", "norm", "cosine", "distance","optimize_kl", "measure_kl", "random_kl", "optimize_spins"};
    int n_tasks = sizeof(tasks)/sizeof(tasks[0]);

    int tk = 0;
    
    printf("\n"); 
    for(tk = 0 ; tk < n_tasks; tk++){
        printf("Usage: %s %s [OPTIONS]\n", argv[0], tasks[tk]);
    }

    printf("\n---------------------------------------------------------------------------\n");
    printf("Please, choose one of the following tasks:\n\n\
   *random*                  To randomly generate coarse-grained representations\n\
                             and measure the associated mapping entropies\n\n\
   *optimize*                To optimize the coarse-grained mapping by minimising\n\
                             its mapping entropy\n\n\
   *measure*                 To measure the mapping entropy of a mapping\n\
                             provided by the user\n\n\
   *norm*                    To calculate the norm of a mapping (provided by the user)\n\
                             throughout a trajectory\n\n\
   *cosine*                  To calculate pairwise distance and cosine between a pair\n\
                             of mappings (provided by the user) throughout a trajectory\n\n\
   *distance*                To calculate the distance matrix between a data set\n\
                             of mappings (provided by the user) over a single conformation\n\n\
   *measure_kl*              To measure the KL divergence version of the mapping entropy for a mapping\n\
                             provided by the user\n\n\
   *optimize_kl*             To optimize the coarse-grained mapping by minimising\n\
                             the KL divergence version of its mapping entropy\n\n\
   *random_kl*               To randomly generate coarse-grained representations\n\
                             and measure KL divergence version of the associated mapping entropies\n\
----------------------------------------------------------------------------\n\n");
 
    printf("Hereafter the list of OPTIONS:\n\n");

    printf("  -p,  --p           FILE       Parameter File (.ini format).\n");
    printf("  -t   --t           FILE       Trajectory File (.xyz format).\n");
    printf("  -e   --e           FILE       Energy File (.txt or .dat format).\n");
    printf("  -c   --code        STR        String that identifies your structure.\n");
    printf("  -m   --m1          FILE       1st Mapping File (.txt or .dat format).\n");
    printf("  -n   --m2          FILE       2nd Mapping File (.txt or .dat format).\n");
    printf("  -x   --matrix      FILE       Mapping Matrix File (.txt or .dat format).\n");
    printf("  -r   --probs       FILE       Probability File (.txt or .dat format).\n");
    printf(" [-h] [--help]                  Give this help list\n");
    printf(" [-v] [--verbose]               Produce verbose output\n\n");           

    printf("Try %s <TASK> for more information about the mandatory options of a specific task\n\n", argv[0]); 
    
    printf("Report bugs to <raffaele.fiorentini@unitn.it>\n"); 
}

void print_help(char *argv[]){
    /**
    * routine that prints some help
    *
    * Parameter
    * ----------
    *
    * `argv[]` : array of command line arguments
    */
    printf("\n");
    printf("Usage: %s %s -p <Parameter FILE> -t <Trajectory FILE> ", argv[0], argv[1]);
    if(strcmp(argv[1], "random") == 0 || strcmp(argv[1], "optimize") == 0 || strcmp(argv[1], "measure") == 0)
        printf("-e <Energy FILE> ");

    printf("-c <Protein Code> "); 

    if(strcmp(argv[1], "measure") == 0 || strcmp(argv[1], "norm") == 0 || strcmp(argv[1], "cosine") == 0 || strcmp(argv[1], "measure_kl") == 0)
       printf("-m <Mapping FILE> ");
    if(strcmp(argv[1], "cosine") == 0)
       printf("-n <2nd Mapping FILE> ");
    if(strcmp(argv[1], "distance") == 0)
       printf("-x <Mapping Matrix FILE>");
    if(strcmp(argv[1], "measure_kl") == 0 || strcmp(argv[1], "optimize_kl") == 0 || strcmp(argv[1], "random_kl") == 0 || strcmp(argv[1], "optimize_spins") == 0)
       printf("-r <Probability FILE> ");
    
    printf("\n");

    printf("   or: %s %s --p <Parameter FILE> --t <Trajectory FILE> ", argv[0], argv[1]);
    if(strcmp(argv[1], "random") == 0 || strcmp(argv[1], "optimize") == 0 || strcmp(argv[1], "measure") == 0)
        printf("--e <Energy FILE> ");

    printf("--code <Protein Code> ");

    if(strcmp(argv[1], "measure") == 0 || strcmp(argv[1], "norm") == 0 || strcmp(argv[1], "cosine") == 0  || strcmp(argv[1], "measure_kl") == 0)
        printf("--m1 <Mapping File> ");
    if(strcmp(argv[1], "cosine") == 0) 
        printf("--m2 <2nd Mapping File> ");
    if(strcmp(argv[1], "distance") == 0)
        printf("--matrix <Mapping Matrix FILE> ");
    if(strcmp(argv[1], "measure_kl") == 0 || strcmp(argv[1], "optimize_kl") == 0 || strcmp(argv[1], "random_kl") == 0 || strcmp(argv[1], "optimize_spins") == 0)
        printf("--prob <Probability FILE> ");

    printf("\n");  

    printf("\n---------------------------------------------------------------------------\n");
    printf("The *%s* task requires the following inputs:\n\n", argv[1]);
    printf("   Parameter FILE            File with the list of parameters in .ini format;\n");
    printf("   Trajectory FILE           File with the trajectory coordinates in .xyz format\n");

    if(strcmp(argv[1], "random") == 0 || strcmp(argv[1], "optimize") == 0 || strcmp(argv[1], "measure") == 0)
        printf("   Energy FILE               File with the value of energy for each frame;\n");

    printf("   Protein Code              String that identifies your structure.\n");

    if(strcmp(argv[1], "measure") == 0 || strcmp(argv[1], "norm") == 0 || strcmp(argv[1], "cosine") == 0 || strcmp(argv[1], "measure_kl") == 0)
        printf("   Mapping FILE              File showing the retained CG sites.\n");
    if(strcmp(argv[1], "cosine") == 0)
        printf("   Mapping FILE 2            2nd File with the retained CG sites.\n");
    if(strcmp(argv[1], "distance") == 0)
        printf("   Mapping Matrix FILE       File with a set of CG mappings (one per row).\n");    
    if(strcmp(argv[1], "measure_kl") == 0 || strcmp(argv[1], "optimize_kl") == 0 || strcmp(argv[1], "random_kl") == 0 || strcmp(argv[1], "optimize_spins") == 0)
        printf("   Probability FILE          File with the probability value for each frame.\n");

    printf("----------------------------------------------------------------------------\n\n"); 
    printf("Hereafter the list of flags:\n\n");
    printf("  -p,  --p           FILE       Parameter File (.ini format).\n");
    printf("  -t   --t           FILE       Trajectory File (.xyz format).\n");

    if(strcmp(argv[1], "random") == 0 || strcmp(argv[1], "optimize") == 0 || strcmp(argv[1], "measure") == 0)
        printf("  -e   --e           FILE       Energy File (.txt or .dat format).\n");
 
    printf("  -c   --code        STR        String that identifies your structure.\n");

    if(strcmp(argv[1], "measure") == 0 || strcmp(argv[1], "norm") == 0 || strcmp(argv[1], "cosine") == 0 || strcmp(argv[1], "measure_kl") == 0)
        printf("  -m   --m1          FILE       Mapping File (.txt or .dat format)\n");
    if(strcmp(argv[1], "cosine") == 0)
        printf("  -n   --m2          FILE       2nd Mapping File (.txt or .dat format)\n");
    if(strcmp(argv[1], "distance") == 0)
        printf("  -x   --matrix      FILE       Mapping Matrix File (.txt or .dat format)\n");
    if(strcmp(argv[1], "measure_kl") == 0 || strcmp(argv[1], "optimize_kl") == 0 || strcmp(argv[1], "random_kl") == 0 || strcmp(argv[1], "optimize_spins") == 0)
        printf("  -r   --probs        FILE      Probability File (.txt or .dat format).\n");

    printf(" [-h] [--help]                  Give this help list\n");
    printf(" [-v] [--verbose]               Produce verbose output\n\n");           

    printf("\n\n");
    printf("Report bugs to <raffaele.fiorentini@unitn.it>\n\n");
}


void check_files(char **pars, char **pars_names, int n_pars, char *argv[]){
    /**
    * routine that checks if all command line arguments are correctly provided
    *
    * Parameter
    * ----------
    *
    * `pars`       : parameters
    * 
    * `pars_names` : names of parameters
    * 
    * `n_pars`     : number of parameters
    * 
    * `argv[]`     : array of command line arguments
    */
    FILE *fe;     // error file 
    int par = 0;
    for (par = 0 ; par < n_pars; par++){
        if(strcmp(pars[par], "-1") == 0){
            fe = fopen("error.dat", "w");
            fprintf(fe, "\nError. The %s is missing. It is required for the task chosen. Look below for further help.\n\n", pars_names[par]);
            printf("\nError. The %s is missing. It is required for the task chosen. Look below for further help.\n\n", pars_names[par]);
            fclose(fe);
            print_help(argv);
            exit(EXIT_FAILURE);
        }
    }
}

void check_parameters(int *pars, char **pars_names, int n_pars){
    /**
    * routine that checks if all mandatory parameters are correctly provided
    *
    * Parameter
    * ----------
    *
    * `pars`       : parameters
    * 
    * `pars_names` : names of parameters
    * 
    * `n_pars`     : number of parameters
    */
    FILE *fe;     // error file 
    int par = 0;
    for (par = 0 ; par < n_pars; par++){
        printf("flag for %s is %d\n", pars_names[par], pars[par]);
        if(pars[par] != 1){
            fe = fopen("error.dat", "w");
            fprintf(fe, "Error. Parameter %s is missing.\n", pars_names[par]);
            printf("Error. Parameter %s is missing.\n", pars_names[par]);
            fclose(fe);
            exit(EXIT_FAILURE);
        }
    }
}

void check_optional_parameters(parameters *cc){
    /**
    * routine that checks optional parameters for the tasks that need them
    *
    * Parameters
    * -----------
    *
    * `cc` : parameters
    */
    FILE *fe;     // error file 
    if (cc->Flag_rsd != 1){
        // rsd not defined
        printf("Parameter rsd not specified. Using default RMSD.\n");
        cc->rsd = 0;
    }
    else{
        if (cc->rsd != 0 && cc->rsd != 1){
            fe = fopen("error.dat", "w");
            fprintf(fe, "Error. Parameter rsd (%d) must be 0 or 1.\n",cc->rsd);
            printf("Error. Parameter rsd (%d) must be 0 or 1.\n",cc->rsd);
            fclose(fe);
            exit(EXIT_FAILURE);
        }
    }
    // clustering parameters
    if (cc->Flag_criterion != 1){
        // clustering criterion not specified
        printf("Warning. Parameters for clustering not specified. Using default values.\n");
        cc->criterion = 0;
        cc->nclust = cc->frames/100;
        cc->rotmats_period = 1;
        printf("nclust = %d\n", cc->nclust);
    }
    else{
        if (cc->criterion == 0){
            // fixed number of clusters
            if (cc->Flag_nclust != 1){
                fe = fopen("error.dat", "w");
                fprintf(fe, "Error. Parameter nclust is missing.\n");
                printf("Error. Parameter nclust is missing.\n");
                fclose(fe);
                exit(EXIT_FAILURE);
            }
            else{
                if (cc->nclust < 2 || cc->nclust >= cc->frames){
                    fe = fopen("error.dat", "w");
                    fprintf(fe, "Error. Parameter nclust (%d) must be between 2 and frames-1 (%d).\n", cc->nclust, cc->frames-1);
                    printf("Error. Parameter nclust (%d) must be between 2 and frames-1 (%d).\n", cc->nclust, cc->frames-1);
                    fclose(fe);
                    exit(EXIT_FAILURE);
                }
            }
        }
        else if(cc->criterion == 3){
            // cophenetic distance clustering
            if (cc->Flag_distance != 1){
                fe = fopen("error.dat", "w");
                fprintf(fe, "Error. Parameter distance is missing.\n");
                printf("Error. Parameter distance is missing.\n");
                fclose(fe);
                exit(EXIT_FAILURE);
            }
            else{
                if (cc->distance < 0.0){
                    fe = fopen("error.dat", "w");
                    fprintf(fe, "Error. Parameter distance must be higher than zero.\n");
                    printf("Error. Parameter distance must be higher than zero.\n");
                    fclose(fe);
                    exit(EXIT_FAILURE);
                }
            }
        }
        else if (cc->criterion == 2){
            // multiple numbers of clusters
            if (cc->Flag_min_nclust != 1){
                fe = fopen("error.dat", "w");
                fprintf(fe, "Error. Parameter min_nclust is missing.\n");
                printf("Error. Parameter min_nclust is missing.\n");
                fclose(fe);
                exit(EXIT_FAILURE);
            }
            else if (cc->Flag_max_nclust != 1){
                fe = fopen("error.dat", "w");
                fprintf(fe, "Error. Parameter max_nclust is missing.\n");
                printf("Error. Parameter max_nclust is missing.\n");
                fclose(fe);
                exit(EXIT_FAILURE);
            }
            else{
                if (cc->min_nclust > (cc->max_nclust-4)){
                    fe = fopen("error.dat", "w");
                    fprintf(fe, "Error. Parameter min_clust (%d) cannot be higher than max_nclust-4 (%d).\n", cc->min_nclust, cc->max_nclust-4);
                    printf("Error. Parameter min_clust (%d) cannot be higher than max_nclust-4 (%d).\n", cc->min_nclust, cc->max_nclust-4);
                    fclose(fe);
                    exit(EXIT_FAILURE);
                }
            }
        }
        else if (cc->criterion == 1){
            // fast clustering
            if (cc->Flag_nclust != 1){
                fe = fopen("error.dat", "w");
                fprintf(fe, "Error. Parameter nclust is missing.\n");
                printf("Error. Parameter nclust is missing.\n");
                fclose(fe);
                exit(EXIT_FAILURE);
            }
            if (cc->Flag_stride != 1){
                fe = fopen("error.dat", "w");
                fprintf(fe, "Error. Parameter stride is missing.\n");
                printf("Error. Parameter stride is missing.\n");
                fclose(fe);
                exit(EXIT_FAILURE);
            }
            else{
                if (cc->frames/cc->stride < cc->nclust){
                    fe = fopen("error.dat", "w");
                    fprintf(fe, "Error. Criterion is %d: parameter nclust (%d) must be larger than frames/stride (%d).\n", cc->criterion, cc->nclust, cc->frames/cc->stride);
                    printf("Error. Criterion is %d: parameter nclust (%d) must be larger than frames/stride (%d).\n", cc->criterion, cc->nclust, cc->frames/cc->stride);
                    fclose(fe);
                    exit(EXIT_FAILURE);
                }
            }
        }
        else{
            fe = fopen("error.dat", "w");
            fprintf(fe, "Error. Parameter criterion (%d) for hierarchical clustering not acceptable\n", cc->criterion);
            printf("Error. Parameter criterion (%d) for hierarchical clustering not acceptable.\n", cc->criterion);
            fclose(fe);
            exit(EXIT_FAILURE);
        }
    }
}

void mandatory_files_present(arguments *arguments, char *argv[]){    
    /**
    * routine that checks if the mandatory files are present
    *
    * Parameters
    * -----------
    *
    * `arguments` : command line arguments
    *
    * `argv[]`    : array of command line arguments
    */
    int n_pars = -1;
    if (strcmp(arguments->task , "optimize") == 0){
        char * pars [] = {arguments->parameter_file, arguments->trajectory_file, arguments->energy_file, arguments->prot_code };
        char * pars_names [] = { "parameter file", "trajectory file", "energy file", "protein code" };
        n_pars = sizeof(pars)/sizeof(pars[0]);
        check_files(pars,pars_names,n_pars,argv);
    }
    else if (strcmp(arguments->task , "random") == 0){
        char * pars [] = {arguments->parameter_file, arguments->trajectory_file, arguments->energy_file, arguments->prot_code };
        char * pars_names [] = { "parameter file", "trajectory file", "energy file", "protein code" };
        n_pars = sizeof(pars)/sizeof(pars[0]);
        check_files(pars,pars_names,n_pars,argv);
    }
    else if (strcmp(arguments->task , "measure") == 0){
        char * pars [] = {arguments->parameter_file, arguments->trajectory_file, arguments->energy_file, arguments->prot_code, arguments->mapping_file };
        char * pars_names [] = { "parameter file", "trajectory file", "energy file", "protein code", "mapping file" };
        n_pars = sizeof(pars)/sizeof(pars[0]);
        check_files(pars,pars_names,n_pars,argv);
    }
    else if (strcmp(arguments->task , "norm") == 0){
        char * pars [] = {arguments->parameter_file, arguments->trajectory_file, arguments->prot_code, arguments->mapping_file };
        char * pars_names [] = { "parameter file", "trajectory file", "protein code", "mapping file" };
        n_pars = sizeof(pars)/sizeof(pars[0]);
        check_files(pars,pars_names,n_pars,argv);
    }
    else if (strcmp(arguments->task , "cosine") == 0){
        char * pars [] = {arguments->parameter_file, arguments->trajectory_file, arguments->prot_code, arguments->mapping_file, arguments->mapping_file2 };
        char * pars_names [] = { "parameter file", "trajectory file", "protein code", "mapping file", "second mapping file" };
        n_pars = sizeof(pars)/sizeof(pars[0]);
        check_files(pars,pars_names,n_pars,argv);
    }
    else if (strcmp(arguments->task , "distance") == 0){
        char * pars [] = {arguments->parameter_file, arguments->trajectory_file, arguments->prot_code, arguments->mapping_matrix };
        char * pars_names [] = { "parameter file", "trajectory file", "protein code", "mapping matrix" };
        n_pars = sizeof(pars)/sizeof(pars[0]);
        check_files(pars,pars_names,n_pars,argv);
    }
    else if (strcmp(arguments->task , "measure_kl") == 0){
        char * pars [] = {arguments->parameter_file, arguments->trajectory_file, arguments->prot_code, arguments->mapping_file, arguments->probability_file };
        char * pars_names [] = { "parameter file", "trajectory file", "protein code", "mapping file", "probability file" };
        n_pars = sizeof(pars)/sizeof(pars[0]);
        check_files(pars,pars_names,n_pars,argv);
    }
    else if (strcmp(arguments->task , "optimize_kl") == 0){
        char * pars [] = {arguments->parameter_file, arguments->trajectory_file, arguments->prot_code, arguments->probability_file };
        char * pars_names [] = { "parameter file", "trajectory file", "protein code", "probability file" };
        n_pars = sizeof(pars)/sizeof(pars[0]);
        check_files(pars,pars_names,n_pars,argv);
    }
    else if (strcmp(arguments->task , "random_kl") == 0){
        char * pars [] = {arguments->parameter_file, arguments->trajectory_file, arguments->prot_code, arguments->probability_file };
        char * pars_names [] = { "parameter file", "trajectory file", "protein code", "probability file" };
        n_pars = sizeof(pars)/sizeof(pars[0]);
        check_files(pars,pars_names,n_pars,argv);
    }
    else if (strcmp(arguments->task , "optimize_spins") == 0){
        char * pars [] = {arguments->parameter_file, arguments->trajectory_file, arguments->prot_code, arguments->probability_file };
        char * pars_names [] = { "parameter file", "trajectory file", "protein code", "probability file" };
        n_pars = sizeof(pars)/sizeof(pars[0]);
        check_files(pars,pars_names,n_pars,argv);
    }

}

void init_parameters(parameters *cc){
    /**
    * routine that initialises the parameters
    * 
    * Parameters
    * -----------
    * 
    * `cc` : parameters object
    */
    cc->atomnum = -1;
    cc->frames = -1;
    cc->cgnum = -1;
    cc->nclust = -1;
    cc->n_mappings = -1;
    cc->MC_steps = -1;
    cc->rotmats_period = -1;
    cc->t_zero = -1.0;
    cc->distance = -1.0;
    cc->criterion = 0;
    cc->max_nclust = -1;
    cc->min_nclust = -1;
    cc->Ncores = 1; // predefined value for Ncores
    cc->decay_time = -1.0;
    cc->rsd = 0;
    cc->stride = -1;
    // init flags
    cc->Flag_atomnum = -1;      
    cc->Flag_frames = -1;
    cc->Flag_cgnum = -1;
    cc->Flag_nclust = -1;
    cc->Flag_n_mappings = -1;
    cc->Flag_MC_steps = -1;
    cc->Flag_rotmats_period = -1;
    cc->Flag_t_zero = -1;
    cc->Flag_distance = -1;
    cc->Flag_criterion = -1;
    cc->Flag_max_nclust = -1;
    cc->Flag_min_nclust = -1;
    cc->Flag_Ncores = -1;
    cc->Flag_decay_time = -1;
    cc->Flag_rsd = -1;
    cc->Flag_stride = -1;
}

void read_ParameterFile(arguments *arguments, parameters *cc){
    /** 
    * routine that reads the input parameter file
    *
    * Parameter
    * ----------
    *
    * `ParameterFileName` : parameter filename
    */
 
    FILE *fp; 

    /* Opening and check if it exist and if it empty */
    fp = fopen(arguments->parameter_file, "r");
    check_empty_file(fp, arguments->parameter_file);
    printf("Initialize Parameters\n");
    init_parameters(cc);
    /* Associate each value of parameter file (cc.VALUE --> e.g. cc.atonmum is the number of atoms in ini file). */
    printf("Reading Parameter FILE...\n");

    /* Parsing the .ini file, Check that each value is an INTEGER */
    ini_parse(arguments->parameter_file, handler, cc);
    //printf("cc_Ncores %d\n", cc->Ncores);
    int n_pars = -1;
    if (strcmp(arguments->task , "optimize") == 0){
        int pars [] = { cc->Flag_atomnum, cc->Flag_frames, cc->Flag_cgnum, cc->Flag_MC_steps };
        char * pars_names [] = { "atomnum", "frames", "cgnum", "MC_steps"};
        n_pars = sizeof(pars)/sizeof(pars[0]);
        check_parameters(pars, pars_names, n_pars);
        check_optional_parameters(cc);
    }
    else if (strcmp(arguments->task , "random") == 0){
        int pars [] = { cc->Flag_atomnum, cc->Flag_frames, cc->Flag_cgnum, cc->Flag_n_mappings };
        char * pars_names [] = { "atomnum", "frames", "cgnum", "n_mappings"};
        n_pars = sizeof(pars)/sizeof(pars[0]);
        check_parameters(pars, pars_names, n_pars);
        check_optional_parameters(cc);
    }
    else if (strcmp(arguments->task , "measure") == 0){
        int pars [] = { cc->Flag_atomnum, cc->Flag_frames, cc->Flag_cgnum };
        char * pars_names [] = { "atomnum", "frames", "cgnum"};
        n_pars = sizeof(pars)/sizeof(pars[0]);
        check_parameters(pars, pars_names, n_pars);
        check_optional_parameters(cc);
    }
    else if (strcmp(arguments->task , "norm") == 0){
        int pars [] = { cc->Flag_atomnum, cc->Flag_frames, cc->Flag_cgnum };
        char * pars_names [] = { "atomnum", "frames", "cgnum"};
        n_pars = sizeof(pars)/sizeof(pars[0]);
        check_parameters(pars, pars_names, n_pars);
    }
    else if (strcmp(arguments->task , "cosine") == 0){
        int pars [] = { cc->Flag_atomnum, cc->Flag_frames, cc->Flag_cgnum };
        char * pars_names [] = { "atomnum", "frames", "cgnum"};
        n_pars = sizeof(pars)/sizeof(pars[0]);
        check_parameters(pars, pars_names, n_pars);
    }
    else if (strcmp(arguments->task , "distance") == 0){
        int pars [] = { cc->Flag_atomnum, cc->Flag_frames, cc->Flag_cgnum, cc->Flag_n_mappings };
        char * pars_names [] = { "atomnum", "frames", "cgnum", "n_mappings"};
        n_pars = sizeof(pars)/sizeof(pars[0]);
        check_parameters(pars, pars_names, n_pars);
    }
    else if (strcmp(arguments->task , "measure_kl") == 0){
        int pars [] = { cc->Flag_atomnum, cc->Flag_frames, cc->Flag_cgnum };
        char * pars_names [] = { "atomnum", "frames", "cgnum"};
        n_pars = sizeof(pars)/sizeof(pars[0]);
        check_parameters(pars, pars_names, n_pars);
        check_optional_parameters(cc);
    }
    else if (strcmp(arguments->task , "optimize_kl") == 0){
        int pars [] = { cc->Flag_atomnum, cc->Flag_frames, cc->Flag_cgnum, cc->Flag_MC_steps};
        char * pars_names [] = { "atomnum", "frames", "cgnum", "MC_steps"};
        n_pars = sizeof(pars)/sizeof(pars[0]);
        check_parameters(pars, pars_names, n_pars);
        check_optional_parameters(cc);
    }
    else if (strcmp(arguments->task , "random_kl") == 0){
        int pars [] = { cc->Flag_atomnum, cc->Flag_frames, cc->Flag_cgnum, cc->Flag_n_mappings};
        char * pars_names [] = { "atomnum", "frames", "cgnum", "n_mappings"};
        n_pars = sizeof(pars)/sizeof(pars[0]);
        check_parameters(pars, pars_names, n_pars);
        check_optional_parameters(cc);
    }
    else if (strcmp(arguments->task , "optimize_spins") == 0){
        int pars [] = { cc->Flag_atomnum, cc->Flag_frames, cc->Flag_cgnum, cc->Flag_MC_steps};
        char * pars_names [] = { "atomnum", "frames", "cgnum", "MC_steps"};
        n_pars = sizeof(pars)/sizeof(pars[0]);
        check_parameters(pars, pars_names, n_pars);
        check_optional_parameters(cc);
    }
    fclose(fp);
}

FILE * open_file_w(char *filename) {
    /**
    * routine that opens a file in write mode
    */

    FILE *fp;
    FILE *fe;

    if ((fp = fopen(filename, "w")) == NULL) {
        printf("Error. Could not open file %s.\n.Exiting.\n", filename);
        fe = fopen("error.dat", "w");
        fprintf(fe, "Error. Could not open file %s.\n.Exiting.\n", filename);
        fclose(fe);
        exit(1);
    }
    return (fp);
}

/**********/
FILE * open_file_r(char *filename) {
    /**
    * routine that opens a file in read mode
    */
    FILE *fp;
    FILE *fe;

    if ((fp = fopen(filename, "r")) == NULL) {
        printf("Error. Could not open file %s.\n.Exiting.\n", filename);
        fe = open_file_w("error.dat");
        fprintf(fe, "Error. Could not open file %s.\n.Exiting.\n", filename);
        fclose(fe);
        exit(1);
    }
    return (fp);
}

/**********/
FILE * open_file_a(char *filename) {
    /**
    * routine that opens a file in append mode
    */
    FILE *fp;
    FILE *fe;

    if ((fp = fopen(filename, "a")) == NULL) {
        printf("Error. Could not open file %s.\n.Exiting.\n", filename);
        fe = open_file_w("error.dat");
        fprintf(fe, "Error. Could not open file %s.\n.Exiting.\n", filename);
        fclose(fe);
        exit(1);
    }
    return (fp);
}

/**********/

void close_file(FILE *fp) {
    /**
    * routine that closes a file
    */
    fclose(fp);
}
