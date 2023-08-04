#ifndef HDR_OPTIMIZE_SPINS
#define HDR_OPTIMIZE_SPINS

#include <stdio.h>
#include <io.h>

int spin_clustering(int **traj_coords, int *mapping, int atomnum, int frames, int *mapping_clusters, int cgnum); 

void optimize_spins(arguments *arguments, parameters *cc);

#endif