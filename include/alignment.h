#ifndef HDR_ALIGN
#define HDR_ALIGN

#include <mapping.h>
#include <traj.h>

/**
* \class alignments
* \brief structure that defines the current alignments stored in memory
*/
typedef struct alignments{
    double *rmsd_mat;           /*!< condensed pairwise RMSD matrix */
    double **rotation_matrices; /*!< condensed matrix of pairwise rotation matrices */
    double **coms;              /*!< array of centers of mass */
    int rsd;                    /*!< RSD parameter. {0: use the RMSD, 1: use the RSD} */
    double *rmsd_vector;        /*!< RMSD vector for fast, 1D, clustering */
    double **rotation_matrices_vector; /*!< vector of pairwise rotation matrices */
}alignments;

// alignment.c

void free_new_alignment(alignments *new_align);

void free_alignment(alignments *align);

void align_two_frames(double *frame_ref, double *frame_middle, int ref_id, int middle_id, cg_mapping *mapping, alignments *align);

double optimal_alignment(double **x, double **y, int mapping_length, double u[][3]);

void correct_rmsd(alignments *new_align, traj *Trajectory, alignments *prev_align, int cgnum, int removed, int added);

void cycle_alignment(traj *Trajectory, alignments *align, cg_mapping *mapping);

void cycle_alignment_fastclust(traj *Trajectory, alignments *align, cg_mapping *mapping);

void correct_rmsd_fastclust(alignments *new_align, traj *Trajectory, alignments *prev_align, int cgnum, int removed, int added);

void cycle_alignment_stride(traj *Trajectory, alignments *align, cg_mapping *mapping);

void align_traj_to_reference(traj *Trajectory, int ref_id);

#endif
