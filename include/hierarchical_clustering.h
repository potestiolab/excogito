#ifndef HDR_CLUSTER
#define HDR_CLUSTER

/**
* \class clust_params
* \brief structure that defines the parameters for hierarchical clustering
*/
typedef struct clust_params{
    int crit;           /*!< criterion for clustering structures. {0: single nclust, 1: distance-based, 2: multiple nclust, 3: fast clustering} */
    int ncl;            /*!< number of clusters (if crit is 0) */
    int max_ncl;        /*!< maximum number of clusters (if crit is 2) */
    int min_ncl;        /*!< minimum number of clusters (if crit is 2) */
    double c_distance;  /*!< maximum cophenetic distance (if crit is 1) */
}clust_params;


void mergesort_merge(double **arr, int l, int m, int r, int dim, int dims);
void my_mergesort(double **arr, int l, int r, int dim, int dims);
int condensed_index(int frames, int i, int j);
double new_dist(double d_xi, double d_yi, double d_xy, int size_x, int size_y, int size_i);
int is_visited(unsigned char *bitset, int i);
void set_visited(unsigned char *bitset, int i);
void get_max_dist_for_each_cluster(double **Z, double *MD, int frames);
void cluster_monocrit(double **Z, double *MC, int *T, double cutoff, int frames);
void cluster_maxclust_monocrit(double **Z, double *MC, int *T, int n, int max_nc);
void cluster_maxclust_dist(double **Z, int *T, int frames, int Nclust);
void cluster_dist(double **Z, int *T, double cutoff, int frames);
int find(int x, int *self_parent);
int merge(int *self_parent, int *self_size, int next_label, int x, int y);
void label(double **Z, int frames);
void hierarchical_clustering(double *rmsd_mat, int n, int couples, int *size, double **Z);
void compute_clusters_list(int *clusters, int *cluster_list, int *cluster_list_idx, int frames, int Nclust);

#endif