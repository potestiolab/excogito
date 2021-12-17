#ifndef HDR_IO
#define HDR_IO

#include <stdio.h>
#include <string.h>
#include <argp.h>
#include <ini.h>
#include <stdlib.h>

FILE *open_file_w(char *filename);

FILE *open_file_r(char *filename);

FILE *open_file_a(char *filename);

void close_file(FILE *fp);

void print_help_main(char * argv[]); 

void print_help(char *argv[]);
 
/* Program documentation.*/ 
static char doc_main[] =
"\n---------------------------------------------------------------------------\n\
Please, choose one of the following tasks:\n\n\
   *random*		     To randomly generate coarse-grained representations\n\
                             and measure the associated mapping entropies;\n\n\
   *optmize*  		     To optimize the coarse-grained mapping by minimising\n\
                             its mapping entropy\n\n\
   *measure*  	 	     To measure the mapping entropy of a mapping\n\
                             provided by the user (in the form of a .txt file)\n\n\
   *norm*		     To calculate the norm of a mapping (provided by the user)\n\
                             throughout a trajectory\n\n\
   *cosine*		     To calculate pairwise distance and cosine between a pair\n\
			     of mappings (provided by the user) throughout a trajectory\n\n\
   *distance*		     To calculate the distance matrix between a data set\n\
			     of mappings (provided by the user) over a single conformation\n\
----------------------------------------------------------------------------\n\n\
Hereafter the list of OPTIONS:";

/* A description of the arguments we accept. */
static char args_doc_main[] = "random\noptimize\nmeasure\nnorm\ncosine\ndistance";

static struct argp_option options_main[] = {
   {"verbose",      'v',          0,  0,  "Produce verbose output" },
   {"quiet",        'q',          0,  OPTION_HIDDEN,  "Don't produce any output" },
   {"help",         'h',          0,  0,  "Give this help list"},
   {"p",            'p',      "FILE", OPTION_HIDDEN,  "Parameter file in ini.format (mandatory)"},
   {"e",            'e',      "FILE", OPTION_HIDDEN,  "Energy file (mandatory for tasks 0, 1, 2, 3)"},
   {"m1",           'm',      "FILE", OPTION_HIDDEN,  "Mapping file1 (mandatory for tasks 2, 4, 5)" },
   {"m2",           'n',      "FILE", OPTION_HIDDEN,  "Mapping file2 (mandatory for task 5)" },
   {"t",  	    't',      "FILE", OPTION_HIDDEN,  "Trajectory file in .xyz format (mandatory)"},
   {"code",         'c',      "STR",  OPTION_HIDDEN,  "String that identifies your structure (mandatory)"},
   {"matrix",       'x',      "STR",  OPTION_HIDDEN,  "mapping_matrix"},
   {"prob",         'r',      "FILE", OPTION_HIDDEN,  "Probability file"},
   { 0 }
};

/* Used by main to communicate with parse_opt. This structure is used to hold all of the arguments */
typedef struct {
    int   silent, verbose;            
    char *parameter_file;             /*!< input parameter file      */ 
    char *energy_file;                /*!< input energy file 	   */
    char *mapping_file;               /*!< input first mapping file  */
    char *mapping_file2;              /*!< input second mapping file */
    char *trajectory_file;            /*!< input trajectory file 	   */
    char *prot_code;                  /*!< protein code 		   */
    char *task;                       /*!< task: (optimize, random, measure, norm, cosine, distance) */
    char *mapping_matrix;             /*!< input mapping matrix */
    char *probability_file;            /*!< input probability file */
} arguments; 


static error_t parse_opt (int key, char *arg, struct argp_state *state){

    // Get the input argument from argp_parse, which we know is a pointer to our arguments structure.
    arguments *arguments = state->input; 

    switch (key)
      {
      case ARGP_KEY_ARG:
        if (state->arg_num == 0)      // we impose that arguments.task is the first argument after ./smap_tool 
            arguments->task = arg;
        else{                         // if command-line arguments are more than one...
            printf("\nError. Only one task is allowed. Check that each flag (e.g. -p) is followed by its specific argument\n");
            printf("Check also that each optional flag (-h --usage or -v) preceedes the other ones\n"); 
            printf("Example:\n\n ./smap_tool <Task> -v -p <parameter FILE> is accepted\n"); 
            printf(" ./smap_tool <Task> -p -v <parameter FILE> is not accepted.\n\n");
            printf("Look below for more information\n");  
            print_help(state->argv);
            exit(EXIT_FAILURE);    
        } 
        break;
      case 'q': case 's':
        arguments->silent = 1;
        break;
      case 'h':
        if (state->arg_num == 0){     // cosi evitiamo che stampi il SUO help se mettiamo -h dopo il task (ex. /smap TASK -h)
            print_help_main(state->argv);
            exit(EXIT_FAILURE);
        }
        else{
            print_help(state->argv); 
            exit(EXIT_FAILURE);
        } 
        break; 
      case 'v':
        arguments->verbose = 1;
        break;
      case 'p':
        arguments->parameter_file = arg;
        break;
      case 'e':
        arguments->energy_file = arg;
        break;
      case 'm':
        arguments->mapping_file = arg;
        break;
      case 'n': 
        arguments->mapping_file2 = arg;
        break; 
      case 'x':
        arguments->mapping_matrix = arg;
        break; 
      case 't':
        arguments->trajectory_file = arg;
        break;
      case 'c': 
        arguments->prot_code = arg; 
        break; 
      case 'r': 
        arguments->probability_file = arg; 
        break; 

      default:
        return ARGP_ERR_UNKNOWN;
      }
    return 0;
}

// Our argp parser.
static struct argp argp = { options_main, parse_opt, args_doc_main, doc_main };


// define a structure for holding the values in "Parameters".
typedef struct{
    int   atomnum;
    int   frames;
    int   cgnum;
    int   nclust;
    int   n_mappings;
    int   MC_steps;
    int   rotmats_period;
    float t_zero;
    float distance;
    int	  criterion;
    int   max_nclust;
    int   min_nclust;
    int   Ncores;
    float   decay_time;
    int   rsd;
    int   stride;
    int   Flag_atomnum;
    int   Flag_frames;
    int   Flag_cgnum;
    int   Flag_nclust;
    int   Flag_n_mappings;
    int   Flag_MC_steps;
    int   Flag_rotmats_period;
    int   Flag_t_zero;
    int   Flag_distance;
    int	  Flag_criterion;
    int   Flag_max_nclust;
    int   Flag_min_nclust;
    int   Flag_Ncores;
    int   Flag_decay_time;
    int   Flag_rsd;
    int   Flag_stride;
} parameters; 

int handler(void* config, const char* section, const char* name, const char* value);  // eliminate static (before int)  

parameters pp_config(parameters config);   // ppconfig is a function... 

void print_usage_main(char *argv[]);



void check_files(char **pars, char **pars_names, int n_pars, char *argv[]); 

void check_empty_file(FILE *f, char *filename);

int n_rows(FILE *f); 

void check_empty_rows(char *str); 

void check_int_string(const char *str, int row, char *fname);

void check_int_string_iniFile(const char *str, char *fname, char *name); 

void check_argv_errors(char *argv[], int argc); 

void check_float_string(char *str, int row, char *fname); 

void check_float_string_iniFile(const char *str, char *fname, char *name); 

int columns(char *string); 

void mandatory_files_present(arguments *arguments, char *argv[]);

void read_ParameterFile(arguments *arguments, parameters *cc);

void check_optional_parameters(parameters *cc);

void check_parameters(int *pars, char **pars_names, int n_pars);

struct cg_mapping; 
//void load_mapping_matrix(char *mappings_filename, FILE *f_out_l, struct cg_mapping *mapping_matrix[], int nmaps); 
void read_mapping_matrix(char *mappings_filename, FILE *f_out_l, struct cg_mapping *mapping_matrix[], int nmaps);
#endif
