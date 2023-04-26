#ifndef HDR_MAPPING
#define HDR_MAPPING

#include <stdio.h>

/**
* \class cg_mapping
* \brief structure that defines a cg mapping
*/
typedef struct cg_mapping{
    int n_at;       /*!< number of atoms in the atomistic structure */
    int n_cg;       /*!< number of CG sites */
    int *mapping;   /*!< binary array defining the CG mapping */
    int *idx_cluster;    /*!< 1D array of indices that tells you in what cluster the config belongs. */           //(!)
    int *omega;     /*!< 1D array of omega_1. */        				//(!)  
    double smap;    /*!< value of mapping entropy */
    double res;     /*!< value of resolution */
    int *clusters;  /*!< array CG macrostates */
    int *size;      /*!< sizes of CG macrostates */
    double *norms;  /*!< moduli of CG mapping over the trajectory */
} cg_mapping;

struct arguments;

void free_mapping(cg_mapping *mapping);

void convert_mapping(cg_mapping *mapping, FILE *f_out);

void generate_random_mapping(cg_mapping *mapping, FILE *f_out);

void update_mapping(cg_mapping *curr_mapping, cg_mapping *old_mapping, int frames);

void read_MappingFile(char *MappingFileName, FILE *f_out_l, cg_mapping *mapping); 

void read_mapping_matrix(char *mappings_filename, FILE *f_out_l, cg_mapping *mapping_matrix[], int nmaps); 

///void load_mapping_matrix(char *mappings_filename, FILE *f_out_l, cg_mapping *mapping_matrix[], int nmaps);


#endif
