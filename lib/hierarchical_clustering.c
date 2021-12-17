/**
* \class hierarchical_clustering
* \brief library of functions that perform hierarchical clustering
*
* Credits to scipy authors:
*
* Copyright (C) Damian Eads, 2007-2008. New BSD License.
*
* hierarchy.py (derived from cluster.py, http://scipy-cluster.googlecode.com)
*
* Author: Damian Eads
* Date:   September 22, 2007
*
* Copyright (c) 2007, 2008, Damian Eads
*
* All rights reserved.
*
* Redistribution and use in source and binary forms, with or without
* modification, are permitted provided that the following conditions
* are met:
*   - Redistributions of source code must retain the above
*     copyright notice, this list of conditions and the
*     following disclaimer.
*   - Redistributions in binary form must reproduce the above copyright
*     notice, this list of conditions and the following disclaimer
*     in the documentation and/or other materials provided with the
*     distribution.
*   - Neither the name of the author nor the names of its
*     contributors may be used to endorse or promote products derived
*     from this software without specific prior written permission.
*
* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
* "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
* LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
* A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
* OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
* SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
* LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
* DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
* THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
* (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
* OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#include <stdlib.h>
#include <math.h>
#include <my_malloc.h>
#include <string.h> /* memset */
#include <hierarchical_clustering.h>

void mergesort_merge(double **arr, int l, int m, int r, int dim, int dims) {
    
    int d;
    int i, j, k;
    int n1 = m - l + 1;
    int n2 = r - m;

    /* create temp arrays */
    double L[n1][dims], R[n2][dims];

    /* Copy data to temp arrays L[] and R[] */
    for (i = 0; i < n1; i++) {
        for (d = 0; d < dims; d++) {
            L[i][d] = arr[l + i][d];
        }
    }
    for (j = 0; j < n2; j++) {
        for (d = 0; d < dims; d++) {
            R[j][d] = arr[m + 1 + j][d];
        }
    }
    /* Merge the temp arrays back into arr[l..r]*/
    i = 0; // Initial index of first subarray 
    j = 0; // Initial index of second subarray 
    k = l; // Initial index of merged subarray 
    while (i < n1 && j < n2) {
        if (L[i][dim] <= R[j][dim]) {
            for (d = 0; d < dims; d++) { arr[k][d] = L[i][d]; }
            i++;
        } else {
            for (d = 0; d < dims; d++) { arr[k][d] = R[j][d]; }
            j++;
        }
        k++;
    }
    /* Copy the remaining elements of L[], if there 
       are any */
    while (i < n1) {
        for (d = 0; d < dims; d++) { arr[k][d] = L[i][d]; }
        i++;
        k++;
    }
    /* Copy the remaining elements of R[], if there 
       are any */
    while (j < n2) {
        for (d = 0; d < dims; d++) { arr[k][d] = R[j][d]; }
        j++;
        k++;
    }
}

/*l is for left index and r is right index of the 
sub-array of arr to be sorted */
void my_mergesort(double **arr, int l, int r, int dim, int dims) {
    if (l < r) {
        // Same as (l+r)/2, but avoids overflow for 
        // large l and h 
        int m = l + (r - l) / 2;
        // Sort first and second halves 
        my_mergesort(arr, l, m, dim, dims);
        my_mergesort(arr, m + 1, r, dim, dims);
        mergesort_merge(arr, l, m, r, dim, dims);
    }
}

int condensed_index(int frames, int i, int j) {
    /**
    * 
    * `frames` : number of observations
    * 
    * `i` : node
    * 
    * `j` : node
    */
    int condensed_index = -1;
    if (i < j) {
        condensed_index = frames * i - (i * (i + 1) / 2) + (j - i - 1);
    } else if (i > j) {
        condensed_index = frames * j - (j * (j + 1) / 2) + (i - j - 1);
    }
    else{
        FILE *fe;     // error file 
        fe = fopen("error.dat", "w");
        fprintf(fe, "Error. arguments to condensed_index cannot be equal (i==j==%d)\n",i);
        printf("Error. arguments to condensed_index cannot be equal (i==j==%d)\n",i);
        fclose(fe);
        exit(EXIT_FAILURE);
    }
    return condensed_index;
}

double new_dist(double d_xi, double d_yi, double d_xy, int size_x, int size_y, int size_i) {
    return (size_x * d_xi + size_y * d_yi) / (size_x + size_y);
}

int is_visited(unsigned char *bitset, int i) {
    /**
    * routine that checks if node i was visited.
    * 
    * Parameters
    * ----------
    *
    * `bitset` : char defining visits
    * 
    * `i` : node
    */
    return bitset[i >> 3] & (1 << (i & 7));
}

void set_visited(unsigned char *bitset, int i) {
    /**
    * routine that marks node i as visited.
    * 
    * Parameters
    * ----------
    *
    * `bitset` : char defining visits
    * 
    * `i` : node
    */
    bitset[i >> 3] |= 1 << (i & 7);
}

void get_max_dist_for_each_cluster(double **Z, double *MD, int frames) {
    /**
    * Get the maximum inconsistency coefficient for each non-singleton cluster.
    * 
    * Parameters
    * ----------
    * `Z` : linkage matrix.
    *
    * `MD` : array to store the result.
    *
    * `frames` : number of observations.
    */
    int k, i_lc, i_rc, root;
    double max_dist, max_l, max_r;
    int *curr_node;
    curr_node = i1t(frames);
    int visited_size = (((frames * 2) - 1) >> 3) + 1;
    unsigned char *visited = (unsigned char *) malloc(visited_size);
    memset(visited, 0, visited_size);
    k = 0;
    curr_node[0] = 2 * frames - 2;
    while (k >= 0) {
        root = curr_node[k] - frames;
        i_lc = (int) Z[root][0];
        i_rc = (int) Z[root][1];
        if (i_lc >= frames) {
            if (!is_visited(visited, i_lc)) {
                set_visited(visited, i_lc);
                k += 1;
                curr_node[k] = i_lc;
                continue;
            }
        }

        if (i_rc >= frames) {
            if (!is_visited(visited, i_rc)) {
                set_visited(visited, i_rc);
                k += 1;
                curr_node[k] = i_rc;
                continue;
            }
        }
        max_dist = Z[root][2];
        if (i_lc >= frames) {
            max_l = MD[i_lc - frames];
            if (max_l > max_dist) {
                max_dist = max_l;
            }
        }
        if (i_rc >= frames) {
            max_r = MD[i_rc - frames];
            if (max_r > max_dist) {
                max_dist = max_r;
            }
        }
        MD[root] = max_dist;
        k -= 1;
    }
    free(visited);
    free_i1t(curr_node);
}

void cluster_monocrit(double **Z, double *MC, int *T, double cutoff, int frames) {
    /**
    * Form flat clusters by monocrit criterion.
    * Parameters
    * ----------
    * `Z` : linkage matrix.
    * 
    * `MC` : monotonic criterion array.
    *
    * `T` : array to store the cluster numbers. The i'th observation belongs to cluster `T[i]`.
    *
    * `cutoff` : Clusters are formed when the MC values are less than or equal to `cutoff`.
    *
    * `frames` : number of observations
    */
    int k, i_lc, i_rc, root, n_cluster = 0, cluster_leader = -1;
    int *curr_node;
    curr_node = i1t(frames);

    int visited_size = (((frames * 2) - 1) >> 3) + 1;
    unsigned char *visited = (unsigned char *) malloc(visited_size);
    //if not visited:
    //    raise MemoryError
    memset(visited, 0, visited_size);
    k = 0;
    curr_node[0] = 2 * frames - 2;
    while (k >= 0) {
        root = curr_node[k] - frames;
        i_lc = (int) Z[root][0];
        i_rc = (int) Z[root][1];
        if (cluster_leader == -1 && MC[root] <= cutoff) {  // found a cluster
            cluster_leader = root;
            n_cluster += 1;
        }
        if (i_lc >= frames && !is_visited(visited, i_lc)) {
            set_visited(visited, i_lc);
            k += 1;
            curr_node[k] = i_lc;
            continue;
        }
        if (i_rc >= frames && !is_visited(visited, i_rc)) {
            set_visited(visited, i_rc);
            k += 1;
            curr_node[k] = i_rc;
            continue;
        }
        if (i_lc < frames) {
            if (cluster_leader == -1) { n_cluster += 1; } // singleton cluster
            T[i_lc] = n_cluster;
        }
        if (i_rc < frames) {
            if (cluster_leader == -1) { n_cluster += 1; } // singleton cluster
            T[i_rc] = n_cluster;
        }
        if (cluster_leader == root) { cluster_leader = -1; }  //back to the leader
        k -= 1;
    }
    free(visited);
    free_i1t(curr_node);
}

void cluster_maxclust_monocrit(double **Z, double *MC, int *T, int frames, int max_nc) {
    /**
    * Form flat clusters by maxclust_monocrit criterion.
    * Parameters
    * ----------
    * `Z` : linkage matrix
    *
    * `MC` : monotonic criterion array
    *
    * `T` : array to store the cluster numbers. The i'th observation belongs to cluster `T[i]`
    *
    * `frames` : number of observations
    *
    * `max_nc` : The maximum number of clusters
    */
    int i, k, i_lc, i_rc, root, nc, lower_idx, upper_idx;
    double thresh;
    int *curr_node;
    curr_node = i1t(frames);
    int visited_size = (((frames * 2) - 1) >> 3) + 1;
    unsigned char *visited = (unsigned char *) malloc(visited_size);
    lower_idx = 0;
    upper_idx = frames - 1;
    while (upper_idx - lower_idx > 1) {
        i = (lower_idx + upper_idx) >> 1;
        thresh = MC[i];
        memset(visited, 0, visited_size);
        nc = 0;
        k = 0;
        curr_node[0] = 2 * frames - 2;
        while (k >= 0) {
            root = curr_node[k] - frames;
            i_lc = (int) Z[root][0];
            i_rc = (int) Z[root][1];
            if (MC[root] <= thresh) {  // this subtree forms a cluster
                nc += 1;
                if (nc > max_nc) { break; }  // illegal
                k -= 1;
                set_visited(visited, i_lc);
                set_visited(visited, i_rc);
                continue;
            }
            if (!is_visited(visited, i_lc)) {
                set_visited(visited, i_lc);
                if (i_lc >= frames) {
                    k += 1;
                    curr_node[k] = i_lc;
                    continue;
                } else {  // singleton cluster
                    nc += 1;
                    if (nc > max_nc) { break; }
                }
            }
            if (!is_visited(visited, i_rc)) {
                set_visited(visited, i_rc);
                if (i_rc >= frames) {
                    k += 1;
                    curr_node[k] = i_rc;
                    continue;
                } else {   // singleton cluster
                    nc += 1;
                    if (nc > max_nc) { break; }
                }
            }
            k -= 1;
        }
        if (nc > max_nc) { lower_idx = i; }
        else { upper_idx = i; }
    }
    free(visited);
    free_i1t(curr_node);
    cluster_monocrit(Z, MC, T, MC[upper_idx], frames);
}

void cluster_maxclust_dist(double **Z, int *T, int frames, int nclust) {
    /**
    * routine that converts the dendrogram into nclust clusters
    *
    * Parameters
    * ----------
    * `Z` : linkage matrix.
    *
    * `T` : array to store the cluster numbers. The i'th observation belongs to cluster `T[i]`.
    *
    * `frames` : number of observations.
    *
    * `nclust` : number of desired clusters.
    */
    double *max_dists;
    max_dists = d1t(frames);
    get_max_dist_for_each_cluster(Z, max_dists, frames);
    cluster_maxclust_monocrit(Z, max_dists, T, frames, nclust);
    free_d1t(max_dists);
}

void cluster_dist(double **Z, int *T, double cutoff, int frames) {
    /*
    * routine that generates clusters from the dendrogram according to a certain cophenetic distance
    *
    * Parameters
    * ----------
    * `Z`: linkage matrix.
    *
    * `T` : The array to store the cluster numbers
    *
    * `cutoff` : the cophenetic distance
    *
    * `frames` : number of observations
    */
    double *max_dists;
    max_dists = d1t(frames);
    get_max_dist_for_each_cluster(Z, max_dists, frames);
    cluster_monocrit(Z, max_dists, T, cutoff, frames);
    free_d1t(max_dists);
}

int find(int x, int *self_parent) {
    int p = x;
    while (self_parent[x] != x) {
        x = self_parent[x];
    }
    while (self_parent[p] != x) {
        p = self_parent[p];
        self_parent[p] = x;
    }
    return x;
}

int merge(int *self_parent, int *self_size, int next_label, int x, int y) {
    self_parent[x] = next_label;
    self_parent[y] = next_label;
    int size = self_size[x] + self_size[y];
    self_size[next_label] = size;
    next_label += 1;
    return size, next_label;
}

void label(double **Z, int frames) {
    /**
    * routine that correctly labels clusters in the unsorted dendrogram
    * 
    * Parameters
    * ----------
    * 
    * `Z` : linkage matrix
    * 
    * `frames` : number of observations
    */
    int dummy;
    int x, y, x_root, y_root;
    int i, k;
    // declaring vectors
    int *self_parent;
    int *self_size;
    int next_label = frames;
    // initialising vectors
    self_parent = i1t(2 * frames - 1);
    self_size = i1t(2 * frames - 1);
    for (i = 0; i < 2 * frames - 1; i++) {
        self_parent[i] = i;
        self_size[i] = 1;
    }
    // labelling
    for (i = 0; i < frames - 1; i++) {
        x = (int) Z[i][0];
        y = (int) Z[i][1];
        x_root = find(x, self_parent);
        y_root = find(y, self_parent);
        if (x_root < y_root) {
            Z[i][0] = x_root;
            Z[i][1] = y_root;
        } else {
            Z[i][0] = y_root;
            Z[i][1] = x_root;
        }
        Z[i][3], next_label = merge(self_parent, self_size, next_label, x_root, y_root);
    }
    free(self_parent);
    free(self_size);
}

void hierarchical_clustering(double *rmsd_mat, int frames, int pairs, int *size, double **Z) {
    /**
    * overall routine for hierarchical clustering
    * 
    * Parameters
    * ----------
    *
    * `rmsd_mat` : condensed pairwise RMSD matrix
    * 
    * `frames` : number of observations
    *
    * `pairs` : possible pairs of structures
    *
    * `size` : size of the clusters (it is nclust long)
    *
    * `Z` : linkage matrix 
    */
    int i, j, k;
    int dummy;
    double *D;
    double current_min, dist;
    // neighbors
    int cluster_chain[frames];
    int chain_length = 0;
    // clusters
    int x, y;
    int nx, ny, ni; //sizes

    D = d1t(pairs);
    // initialisation
    for (i = 0; i < frames; i++) { cluster_chain[i] = 0; }
    for (i = 0; i < frames; i++) { size[i] = 1; } // initial size is 1
    for (i = 0; i < pairs; i++) { D[i] = rmsd_mat[i]; } // copying dist mat
    // linkage matrix
    for (k = 0; k < (frames - 1); k++) {
        if (chain_length == 0) {
            chain_length = 1;
            for (i = 0; i < frames; i++) {
                if (size[i] > 0) {
                    cluster_chain[0] = i;
                    break;
                }
            }
        }
        // Go through chain of neighbors until two mutual neighbors are found.
        while (1) {
            x = cluster_chain[chain_length - 1];
            // We want to prefer the previous element in the chain as the
            // minimum, to avoid potentially going in cycles.
            if (chain_length > 1) {
                y = cluster_chain[chain_length - 2];
                current_min = D[condensed_index(frames, x, y)];
            } else { current_min = INFINITY; }
            for (i = 0; i < frames; i++) {
                if (size[i] == 0 || x == i) { continue; }
                dist = D[condensed_index(frames, x, i)];
                if (dist < current_min) {
                    current_min = dist;
                    y = i;
                }
            }

            if (chain_length > 1) {
                if (y == cluster_chain[chain_length - 2]) {
                    break;
                }
            }
            cluster_chain[chain_length] = y;
            chain_length += 1;
        }
        // Merge clusters x and y and pop them from stack.
        chain_length -= 2;
        // This is a convention used in fastcluster.
        //if x > y:
        //    x, y = y, x
        // get the original numbers of points in clusters x and y
        nx = size[x];
        ny = size[y];
        // Record the new node.
        Z[k][0] = x;
        Z[k][1] = y;
        Z[k][2] = current_min;
        Z[k][3] = nx + ny;
        size[x] = 0; // Cluster x will be dropped.
        size[y] = nx + ny; // Cluster y will be replaced with the new cluster
        // Update the distance matrix.
        for (i = 0; i < frames; i++) {
            ni = size[i];
            if (ni == 0 || i == y) {
                continue;
            }
            D[condensed_index(frames, i, y)] = new_dist(
                    D[condensed_index(frames, i, x)],
                    D[condensed_index(frames, i, y)],
                    current_min, nx, ny, ni);
        }
    }
    // sorting data
    my_mergesort(Z, 0, frames - 2, 2, 4);
    // labelling data
    label(Z, frames);
    // freeing dist_mat
    free_d1t(D);
}

void compute_clusters_list(int *clusters, int *cluster_list, int *cluster_list_idx, int frames, int nclust) {
    /** 
    * routine that computes the list of cluster IDs
    *
    * Parameters
    * ----------
    * 
    * `clusters` : list of labels (one for each frame)
    * 
    * `cluster_list` : ordered list of labels 
    * 
    * `cluster_list_idx` :  is an index vector that stores the sum of populations up to each index 
    * 
    * `frames` : number of observations
    * 
    * `nclust` : number of clusters
    */
    int k, index;
    int counter_vector[nclust]; /*!< array that stores the population of each clusters */
    int cluster_list_idx_cp[nclust];

    for (k = 0; k < nclust; k++) { counter_vector[k] = 0; }
    for (k = 0; k < frames; k++) { counter_vector[clusters[k] - 1] += 1; } // -1 because clusters id are one-based

    cluster_list_idx[0] = 0;
    for (k = 1; k < nclust + 1; k++) {
        cluster_list_idx[k] = cluster_list_idx[k - 1] + counter_vector[k - 1];
    }
    // copy this vector
    for (k = 0; k < nclust + 1; k++) { cluster_list_idx_cp[k] = cluster_list_idx[k]; }
    // now let's build  
    for (k = 0; k < frames; k++) {
        index = cluster_list_idx_cp[clusters[k] - 1]; // again 1-based clusters
        cluster_list[index] = k;
        cluster_list_idx_cp[clusters[k] - 1] += 1;
    }
}
