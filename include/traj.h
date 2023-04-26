#ifndef HDR_TRAJ
#define HDR_TRAJ

#include <stdio.h>
#include <string.h>
#include <argp.h>
#include <ini.h>
#include <io.h>

/**
* \class traj
* \brief structure that defines a MD trajectory
*/
typedef struct traj{
    int frames;             /*!< number of frames in the trajectory */
    double **traj_coords;   /*!< 2D array of xyz coordinates */
    double *energies;       /*!< 1D array of energies. One value per frame. */
    double *energies_cg;    /*!< 1D array of CG energies. One value per frame. */  	//(!)
    int n_at;               /*!< number of atoms in the atomistic structure */
    int pairs;              /*!< number of possible pairs of structures */
    int *strides;           /*!< vector of configurations to consider (criterion 3)*/
    int stride;             /*!< number of configurations between each pivot for clustering (criterion 3)*/ 
    int eff_frames;         /*!< number of effective frames in the trajectory (criterion 3)*/ 
}traj;

struct arguments; 

int check_probabilities(double *probabilities, int prob_length);

void read_EnergyFile(char *EnergyFileName, traj *Trajectory); 

void read_TrajectoryFile(char *TrajFileName, traj *Trajectory);   
   
#endif
