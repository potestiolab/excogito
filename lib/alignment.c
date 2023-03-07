/**
* \class alignment
* \brief library of functions that perform alignments of pairs of structures
*/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <my_malloc.h>
#include <alignment.h>
#include <geometry.h>

void free_new_alignment(alignments *new_align){
    /**
    * routine that frees an alignments object used in criterion 1
    * 
    * Parameters
    * ----------
    * 
    * `new_align`: alignments object
    */
    free_d1t(new_align->rmsd_mat);
    free_d1t(new_align->rmsd_vector);
}

void free_alignment(alignments *align){
    /**
    * routine that frees an alignments object
    * 
    * Parameters
    * ----------
    * 
    * `align`: alignments object
    */
    free_d1t(align->rmsd_mat);
    free_d2t(align->rotation_matrices);
    free_d2t(align->coms);
    free_d1t(align->rmsd_vector);
    free_d2t(align->rotation_matrices_vector);
    free(align);
}

double optimal_alignment(double **x, double **y, int cgnum, double u[][3]) {
    /**
    * routine that computes the Kabsch alignment and the rmsd between two configurations
    * 
    *
    * Parameters
    * ----------
    * 
    * `x`, `y` : CG structures
    * 
    * `cgnum` : length of CG mapping
    * 
    * `u` : rotation matrix
    */
    void myjacobi(double a[][3], int n, double *d, double v[][3], int *nrot);
    int i, j, k, sign[3], order[3], nrot;
    double e, e0;
    double r[3][3], rt[3][3], temp; //, **x, **y;
    double a[3][3], eval[3], evec[3][3];
    double eigenvalues[3], eigenvectors[3][3], b[3][3];
    int speak = 1;

    e0 = 0.0;
    for (i = 0; i < cgnum; i++) {
        e0 += 0.5 * norm_d(x[i], 3) * norm_d(x[i], 3);
        e0 += 0.5 * norm_d(y[i], 3) * norm_d(y[i], 3);
    }

    for (i = 0; i < 3; i++) {
        for (j = 0; j < 3; j++) {
            r[i][j] = 0.0;
            for (k = 0; k < cgnum; k++) {
                r[i][j] += y[k][i] * x[k][j];
            }
            rt[j][i] = r[i][j];
        }
    }

    if (isnan(e0) == 1) printf("Found a NaN in Kabsch alignment at checkpoint 1\n");

    for (i = 0; i < 3; i++) {
        for (j = 0; j < 3; j++) {
            a[i][j] = 0;
            for (k = 0; k < 3; k++) {
                a[i][j] += rt[i][k] * r[k][j];
            }
        }
    }

    myjacobi(a, 3, eval, evec, &nrot);
    /* we add small quantities in order to remove potentially dangerous degeneracies */
    if (isnan(eval[0]) == 1) printf("Found a NaN in Kabsch alignment at checkpoint 2\n");
    if (isnan(eval[1]) == 1) printf("Found a NaN in Kabsch alignment at checkpoint 3\n");
    if (isnan(eval[2]) == 1) printf("Found a NaN in Kabsch alignment at checkpoint 4\n");

    eval[0] += 0.0000000000000001;
    eval[1] += 0.00000000000000001;
    eval[2] += 0.00000000000000002;

    if ((eval[0] < eval[1]) && (eval[0] < eval[2])) {
        order[0] = 1;
        order[1] = 2;
        order[2] = 0;
    }

    if ((eval[1] < eval[0]) && (eval[1] < eval[2])) {
        order[0] = 0;
        order[1] = 2;
        order[2] = 1;
    }

    if ((eval[2] < eval[0]) && (eval[2] < eval[1])) {
        order[0] = 0;
        order[1] = 1;
        order[2] = 2;
    }

    for (i = 0; i < 3; i++) {
        eigenvalues[i] = eval[order[i]];
        for (j = 0; j < 3; j++) {
            eigenvectors[i][j] = evec[j][order[i]];
        }
    }

    normalize_d(eigenvectors[0], 3);
    normalize_d(eigenvectors[1], 3);
    vecprod_d(eigenvectors[0], eigenvectors[1], eigenvectors[2]);
    normalize_d(eigenvectors[2], 3);

    for (i = 0; i < 2; i++) {
        for (j = 0; j < 3; j++) {
            b[i][j] = 0;
            for (k = 0; k < 3; k++) {
                b[i][j] += r[j][k] * eigenvectors[i][k];
            }
        }
        normalize_d(b[i], 3);
    }

    vecprod_d(b[0], b[1], b[2]);
    normalize_d(b[2], 3);

    temp = 0.0;
    for (i = 0; i < 3; i++) {
        for (j = 0; j < 3; j++) {
            temp += b[2][i] * r[i][j] * eigenvectors[2][j];
        }
    }
    sign[2] = +1;
    if (temp < 0)
        sign[2] = -1;
    sign[0] = sign[1] = 1;

    if (fabs(eigenvalues[2]) < 0.0000001) {
        e = e0 - sqrt(eigenvalues[0]) - sqrt(eigenvalues[1]);
    } else {
        e = e0 - sqrt(eigenvalues[0]) - sqrt(eigenvalues[1]) - sign[2] * sqrt(eigenvalues[2]);
    }

    if (isnan(e) == 1) {
        printf("Found a NaN in Kabsch alignment at checkpoint 5 | \n");
        printf("e %lf e0 %lf e1 %lf e2 %lf e3 %lf\n", e, e0, eigenvalues[0], eigenvalues[1], eigenvalues[2]);
    }
/********************/
    e = 2.0 * e / cgnum;
    if (e < 0.0) {
        if (fabs(e) < 1.0e-3) {
            if (speak == 1) {
                printf("Warning. In Kabsch alignment found slightly negative value of e (%e). Roundoff error? I will set it equal to zero.\n",
                       e);
                e = 0.0;
            }
        }
            /* occasionally, when dealing with two practically identical configurations
                the value of e may be slightly negative due to the small offsets and roundoff errors.
                In this case we set it equal to zero. */
        else {
            FILE *fe; 
            fe = fopen("error.dat", "w");
            fprintf(fe, "Error. In Kabsch alignment found negative value of e: (%lf)\n", e);
            printf("Error. In Kabsch alignment found negative value of e: (%lf)\n", e);
            fclose(fe);
            exit(EXIT_FAILURE);
        }
    }
    e = sqrt(e);
    if (isnan(e) == 1) printf("Found a NaN in Kabsch alignment at checkpoint 6\n");
/********************/
    // filling rotation_matrix
    for (i = 0; i < 3; i++) {
        for (j = 0; j < 3; j++) {
            u[i][j] = 0.0;
            for (k = 0; k < 3; k++) {
                u[i][j] += b[k][i] * eigenvectors[k][j];
            }
        }
    }
    return (e);
}

void align_two_frames(double *frame_ref, double *frame_middle, int ref_id, int middle_id, cg_mapping *mapping, alignments *align){
     /**
    * routine that aligns a pair of frames in a trajectory, calling `optimal_alignment`
    * 
    * Parameters
    * ----------
    * 
    * `frame_ref` : reference frame
    * 
    * `frame_middle` : frame in between two pivot clusters
    * 
    * `ref_id` : id (index) of frame_ref in the trajectory
    * 
    * `middle_id` : id (index) of frame_middle in the trajectory
    * 
    * `mapping` : cg_mapping object
    * 
    * `align` : alignments object
    */ 
    double **x, **y;
    double u[3][3], cm1[3], cm2[3];
    int vector_idx, rot_one, rot_two;
    double rmsd;
    zero_vec_d(cm1, 3);
    zero_vec_d(cm2, 3);
    x = d2t(mapping->n_cg, 3);
    y = d2t(mapping->n_cg, 3);
    int i, j;
    for (i = 0; i < mapping->n_at; i++) {
        if (mapping->mapping[i] == 1) {
            for (j = 0; j < 3; j++) {
                cm1[j] += frame_ref[i * 3 + j] / mapping->n_cg;
                cm2[j] += frame_middle[i * 3 + j] / mapping->n_cg;
            }
        }
    }
    align->coms[middle_id][0] = cm2[0];
    align->coms[middle_id][1] = cm2[1];
    align->coms[middle_id][2] = cm2[2];
    int idx = 0;
    for (i = 0; i < mapping->n_at; i++) {
        if (mapping->mapping[i] == 1) {
            for (j = 0; j < 3; j++) {
                x[idx][j] = frame_ref[i * 3 + j] - cm1[j];
                y[idx][j] = frame_middle[i * 3 + j] - cm2[j];
            }
            idx += 1;
        }
    }
    
    // extracting rmsd
    rmsd = optimal_alignment(x, y, mapping->n_cg, u);
    if (ref_id < middle_id){vector_idx = middle_id*2;}
    else if(ref_id > middle_id){vector_idx = middle_id*2 + 1;}
    else{printf("error");}
    align->rmsd_vector[vector_idx] = rmsd;
    //printf("%d vs %d : rmsd = %lf\n", ref_id, middle_id, rmsd);
    for (rot_one = 0; rot_one < 3; rot_one++) {
        for (rot_two = 0; rot_two < 3; rot_two++) {
            align->rotation_matrices_vector[vector_idx][rot_one * 3 + rot_two] = u[rot_one][rot_two];
        }
    }
    //printf("%d vs %d : rot_vect = %lf %lf %lf %lf %lf %lf %lf %lf %lf\n", ref_id, middle_id, align->rotation_matrices_vector[vector_idx][0], align->rotation_matrices_vector[vector_idx][1], align->rotation_matrices_vector[vector_idx][2],align->rotation_matrices_vector[vector_idx][3],align->rotation_matrices_vector[vector_idx][4],align->rotation_matrices_vector[vector_idx][5],align->rotation_matrices_vector[vector_idx][6], align->rotation_matrices_vector[vector_idx][7],align->rotation_matrices_vector[vector_idx][8]);
    //printf("result %lf\n", rmsd);
    free_d2t(x);
    free_d2t(y);
}

void cycle_alignment_stride(traj *Trajectory, alignments *align, cg_mapping *mapping) {
    /**
    * routine that cycles over all pairs of frames in a trajectory, calling `optimal_alignment`
    * 
    * Parameters
    * ----------
    * 
    * `Trajectory` : traj object
    * 
    * `align` : alignments object
    * 
    * `mapping` : cg_mapping object
    */ 
    int idx; // index for the CG structures 
    int n, s;
    int knt = 0;
    int i, j;
    double **x, **y;
    double u[3][3], cm1[3], cm2[3];
    int rot_one, rot_two, rot_idx;
    printf("calculating rotation matrices..\n");
    // allocating vectors
    x = d2t(mapping->n_cg, 3);
    y = d2t(mapping->n_cg, 3);
    int st_id = 0; // structure ID
    //printf("rmsd_mat = \n");
    for (n = 0; n < Trajectory->frames; n++) {
        // check if this structure has to be considered
        if (Trajectory->strides[n] == 1){
            // insert structure and RMSD computation
            //printf("structure n %d\n",n);
            zero_vec_d(cm1, 3);
            for (i = 0; i < mapping->n_at; i++) {
                if (mapping->mapping[i] == 1) {
                    for (j = 0; j < 3; j++) {
                        cm1[j] += Trajectory->traj_coords[n][i * 3 + j] / mapping->n_cg;
                    }
                }
            }
            // updating coms vector
            align->coms[n][0] = cm1[0];
            align->coms[n][1] = cm1[1];
            align->coms[n][2] = cm1[2];
            // filling cg structure
            idx = 0;
            for (i = 0; i < mapping->n_at; i++) {
                if (mapping->mapping[i] == 1) {
                    for (j = 0; j < 3; j++) {
                        x[idx][j] = Trajectory->traj_coords[n][i * 3 + j] - cm1[j];
                    }
                    idx += 1;
                }
            }
            // cycling over all the other structures
            for (s = n + 1; s < Trajectory->frames; s++) {
                if (Trajectory->strides[s] == 1){
                    //printf("structure s %d\n",s);
                    zero_vec_d(cm2, 3);
                    for (i = 0; i < mapping->n_at; i++) {
                        if (mapping->mapping[i] == 1) {
                            for (j = 0; j < 3; j++) {
                                cm2[j] += Trajectory->traj_coords[s][i * 3 + j] / mapping->n_cg;
                            }
                        }
                    }
                    idx = 0; // once again idx is set to zero prior to the creation of the structures
                    for (i = 0; i < mapping->n_at; i++) {
                        if (mapping->mapping[i] == 1) {
                            for (j = 0; j < 3; j++) {
                                y[idx][j] = Trajectory->traj_coords[s][i * 3 + j] - cm2[j];
                            }
                            idx += 1;
                        }
                    }
                    // Alignment AND RMSD computation
                    align->rmsd_mat[knt] = optimal_alignment(x, y, mapping->n_cg, u);
                    if (align->rsd != 0) {
                        if (align->rsd == 1) {
                            align->rmsd_mat[knt] = align->rmsd_mat[knt] * sqrt(mapping->n_cg);
                        }
                    }
                    for (rot_one = 0; rot_one < 3; rot_one++) {
                        for (rot_two = 0; rot_two < 3; rot_two++) {
                            align->rotation_matrices[knt][rot_one * 3 + rot_two] = u[rot_one][rot_two];
                        }
                    }
                    //printf("rmsd[%d] = %lf\n", knt, align->rmsd_mat[knt]);
                    knt++;
                }
            }
        }
    }
    free_d2t(x);
    free_d2t(y);
}

void cycle_alignment_fastclust(traj *Trajectory, alignments *align, cg_mapping *mapping) {
    /**
    * routine that computes the alignments if clustering must be fast
    * 
    * Parameters
    * ----------
    * 
    * `Trajectory` : traj object
    * 
    * `align` : alignments object
    * 
    * `mapping` : cg_mapping object
    */ 
    cycle_alignment_stride(Trajectory, align, mapping);
    printf("end of cycle_alignment_stride\n");
    int cl_id;
    int prev = Trajectory->frames -1;
    int next = Trajectory->frames - 1 - (Trajectory->frames-1)%Trajectory->stride;
    for (cl_id = Trajectory->frames - 2; cl_id > 0; cl_id--){
        if (Trajectory->strides[cl_id] == 0){
            align_two_frames(Trajectory->traj_coords[next], Trajectory->traj_coords[cl_id], next, cl_id, mapping, align);
            align_two_frames(Trajectory->traj_coords[prev], Trajectory->traj_coords[cl_id], prev, cl_id, mapping, align);
        }
        else{
            prev = next;
            next = prev - Trajectory->stride;
        }
    }
    printf("end of cycle_alignment_fastclust\n");
    //for (cl_id = 0; cl_id < Trajectory->frames*2 ; cl_id++){
    //    printf("rmsd_vector[%d] = %lf\n", cl_id,align->rmsd_vector[cl_id]);
    //}
}

void cycle_alignment(traj *Trajectory, alignments *align, cg_mapping *mapping) {
    /**
    * routine that cycles over all pairs of frames in a trajectory, calling `optimal_alignment`
    * 
    * Parameters
    * ----------
    * 
    * `Trajectory` : traj object
    * 
    * `align` : alignments object
    * 
    * `mapping` : cg_mapping object
    */ 
    int idx; // index for the CG structures 
    int n, s;
    int knt = 0;
    int i, j;
    double **x, **y;
    double u[3][3], cm1[3], cm2[3];
    int rot_one, rot_two, rot_idx;
    printf("calculating rotation matrices..\n");
    // allocating vectors
    x = d2t(mapping->n_cg, 3);
    y = d2t(mapping->n_cg, 3);
    for (n = 0; n < Trajectory->frames; n++) {
        // insert structure and RMSD computation
        zero_vec_d(cm1, 3);
        for (i = 0; i < mapping->n_at; i++) {
            if (mapping->mapping[i] == 1) {
                for (j = 0; j < 3; j++) {
                    cm1[j] += Trajectory->traj_coords[n][i * 3 + j] / mapping->n_cg;
                }
            }
        }
        // updating coms vector
        align->coms[n][0] = cm1[0];
        align->coms[n][1] = cm1[1];
        align->coms[n][2] = cm1[2];
        // filling cg structure
        idx = 0;
        for (i = 0; i < mapping->n_at; i++) {
            if (mapping->mapping[i] == 1) {
                for (j = 0; j < 3; j++) {
                    x[idx][j] = Trajectory->traj_coords[n][i * 3 + j] - cm1[j];
                }
                idx += 1;
            }
        }
        // cycling over all the other structures
        for (s = n + 1; s < Trajectory->frames; s++) {
            zero_vec_d(cm2, 3);
            for (i = 0; i < mapping->n_at; i++) {
                if (mapping->mapping[i] == 1) {
                    for (j = 0; j < 3; j++) {
                        cm2[j] += Trajectory->traj_coords[s][i * 3 + j] / mapping->n_cg;
                    }
                }
            }
            idx = 0; // once again idx is set to zero prior to the creation of the structures
            for (i = 0; i < mapping->n_at; i++) {
                if (mapping->mapping[i] == 1) {
                    for (j = 0; j < 3; j++) {
                        y[idx][j] = Trajectory->traj_coords[s][i * 3 + j] - cm2[j];
                    }
                    idx += 1;
                }
            }
            // Alignment AND RMSD computation
            align->rmsd_mat[knt] = optimal_alignment(x, y, mapping->n_cg, u);
            if (align->rsd != 0) {
                if (align->rsd == 1) {
                    align->rmsd_mat[knt] = align->rmsd_mat[knt] * sqrt(mapping->n_cg);
                }
            }
            for (rot_one = 0; rot_one < 3; rot_one++) {
                for (rot_two = 0; rot_two < 3; rot_two++) {
                    align->rotation_matrices[knt][rot_one * 3 + rot_two] = u[rot_one][rot_two];
                }
            }
            //printf("knt %d : rmsd = %lf\n", knt, align->rmsd_mat[knt]);
            //printf("knt %d : rot_mat[:3] = %lf %lf %lf\n", knt, align->rotation_matrices[knt][0], align->rotation_matrices[knt][1], align->rotation_matrices[knt][2]);
            knt++;
        }
    }
    free_d2t(x);
    free_d2t(y);
}

void correct_rmsd(alignments *new_align, traj *Trajectory, alignments *align, int cgnum, int removed, int added) {
    /**
    * routine that computes the rmsd matrix without aligning frames over frames
    * 
    * Parameters
    * ----------
    * 
    * `new_align` : trial alignments object
    * 
    * `Trajectory` : traj object
    * 
    * `align` : alignments object
    * 
    * `cgnum` : number of CG sites (useful to normalize)
    * 
    * `removed` : index of removed atom
    * 
    * `added` : index of added atom
    */
    int n, s, rot_idx, real_rot_idx;
    int knt = 0;
    int i, j;
    double str_i_j_rm, str_i_j_ad;
    double msd_removed, msd_added, v[3];
    for (n = 0; n < Trajectory->frames; n++) {
        rot_idx = n * Trajectory->frames - (n * (n + 1)) / 2 - n - 1;
        for (s = n + 1; s < Trajectory->frames; s++) {
            msd_removed = 0.0;
            msd_added = 0.0;
            real_rot_idx = rot_idx + s;
            //printf("knt %d real_rot_idx %d\n", knt, real_rot_idx);
            for (j = 0; j < 3; j++) {
                str_i_j_rm = align->coms[s][j] +
                             align->rotation_matrices[real_rot_idx][j * 3] * (Trajectory->traj_coords[n][removed * 3] - align->coms[n][0]) +
                             align->rotation_matrices[real_rot_idx][j * 3 + 1] *
                             (Trajectory->traj_coords[n][removed * 3 + 1] - align->coms[n][1]) +
                             align->rotation_matrices[real_rot_idx][j * 3 + 2] * (Trajectory->traj_coords[n][removed * 3 + 2] - align->coms[n][2]);
                str_i_j_ad =
                        align->coms[s][j] + align->rotation_matrices[real_rot_idx][j * 3] * (Trajectory->traj_coords[n][added * 3] - align->coms[n][0]) +
                        align->rotation_matrices[real_rot_idx][j * 3 + 1] * (Trajectory->traj_coords[n][added * 3 + 1] - align->coms[n][1]) +
                        align->rotation_matrices[real_rot_idx][j * 3 + 2] * (Trajectory->traj_coords[n][added * 3 + 2] - align->coms[n][2]);
                v[j] = str_i_j_rm - Trajectory->traj_coords[s][removed * 3 + j];
                msd_removed += v[j] * v[j];
                v[j] = str_i_j_ad - Trajectory->traj_coords[s][added * 3 + j];
                msd_added += v[j] * v[j];
            }
            if (align->rsd == 0) {// no rsd, standard rmsd
                new_align->rmsd_mat[knt] = sqrt(align->rmsd_mat[knt] * align->rmsd_mat[knt] - msd_removed / cgnum + msd_added / cgnum);
            } else if (align->rsd == 1) {
                new_align->rmsd_mat[knt] = sqrt(align->rmsd_mat[knt] * align->rmsd_mat[knt] - msd_removed + msd_added);
            }
            knt += 1;
        }
    }
}

double correct_rmsd_two_frames(traj *Trajectory, double u[9], double com_ref[3], double com_other[3], int cgnum, int removed, int added, int ref_id, int other_id, double prev_rmsd){
    /**
    * routine that corrects the rmsd between two frames
    * 
    * Parameters
    * ----------
    * 
    * `Trajectory` : traj object
    * 
    * `u` : rotation matrix
    * 
    * `com_ref` : reference center of mass
    * 
    * `com_other` : other center of mass
    * 
    * `removed` : index of removed atom
    * 
    * `added` : index of added atom
    * 
    * `ref_id` : index of reference frame
    * 
    * `other_id` : index of other frame
    * 
    * `prev_rmsd` : previous rmsd
    * 
    */
    double msd_removed, msd_added, str_i_j_rm, str_i_j_ad, rmsd, v[3];
    int j;
    msd_removed = 0.0;
    msd_added = 0.0;
    //printf("com[%d](%lf %lf %lf) com[%d](%lf %lf %lf)\n", other_id, prev_align->coms[other_id][0],prev_align->coms[other_id][1], prev_align->coms[other_id][2], ref_id, prev_align->coms[ref_id][0],prev_align->coms[ref_id][1],prev_align->coms[ref_id][2] );
    for (j = 0; j < 3; j++) {
        str_i_j_rm = com_other[j] + u[j * 3] * (Trajectory->traj_coords[ref_id][removed * 3] - com_ref[0]) +
                    u[j * 3 + 1] * (Trajectory->traj_coords[ref_id][removed * 3 + 1] - com_ref[1]) +
                    u[j * 3 + 2] * (Trajectory->traj_coords[ref_id][removed * 3 + 2] - com_ref[2]);
        str_i_j_ad = com_other[j] + u[j * 3] * (Trajectory->traj_coords[ref_id][added * 3] - com_ref[0]) +
                    u[j * 3 + 1] * (Trajectory->traj_coords[ref_id][added * 3 + 1] - com_ref[1]) +
                    u[j * 3 + 2] * (Trajectory->traj_coords[ref_id][added * 3 + 2] - com_ref[2]);
        v[j] = str_i_j_rm - Trajectory->traj_coords[other_id][removed * 3 + j];
        msd_removed += v[j] * v[j];
        v[j] = str_i_j_ad - Trajectory->traj_coords[other_id][added * 3 + j];
        msd_added += v[j] * v[j];
    }
    // no rsd, standard rmsd
    rmsd = sqrt(prev_rmsd * prev_rmsd - msd_removed / cgnum + msd_added / cgnum);
    //} else if (prev_align->rsd == 1) {
    //    rmsd = sqrt(prev_rmsd * prev_rmsd - msd_removed + msd_added);
    //}
    //printf("previous rmsd %d vs %d = %lf , new one %lf, delta = %lf\n", ref_id, other_id, prev_rmsd, rmsd, fabs(prev_rmsd-rmsd));
    return rmsd;
}

void correct_rmsd_fastclust(alignments *new_align, traj *Trajectory, alignments *prev_align, int cgnum, int removed, int added) {
    /**
    * routine that computes the rmsd matrix without aligning frames over frames
    * 
    * Parameters
    * ----------
    * 
    * `new_rmsd_mat` : new condensed pairwise RMSD matrix
    * 
    * `Trajectory` : traj object
    * 
    * `align` : alignments object
    * 
    * `cgnum` : number of CG sites (useful to normalize)
    * 
    * `removed` : index of removed atom
    * 
    * `added` : index of added atom
    */
    int n, s, rot_idx, real_rot_idx;
    int n_stride = 0; // current stride
    int s_stride; // second stride
    int knt = 0;
    //printf("\nnew_rmsd_mat\n");
    for (n = 0; n < Trajectory->frames; n++) {
        if (Trajectory->strides[n] == 1){
            s_stride = n_stride + 1;
            // rot_idx = n * Trajectory->frames - (n * (n + 1)) / 2 - n - 1;
            rot_idx = n_stride * (Trajectory->eff_frames) - (n_stride * (n_stride + 1)) / 2 - n_stride - 1;
            for (s = n + 1; s < Trajectory->frames; s++) {
                if (Trajectory->strides[s] == 1){
                    real_rot_idx = rot_idx + s_stride;
                    new_align->rmsd_mat[real_rot_idx] = correct_rmsd_two_frames(Trajectory, prev_align->rotation_matrices[real_rot_idx], prev_align->coms[n],prev_align->coms[s],  cgnum, removed, added, n, s, prev_align->rmsd_mat[real_rot_idx]);
                    s_stride += 1;
                }
            }
            n_stride += 1;
        }
    }
    //printf("1d vector");
    int cl_id;
    int prev = Trajectory->frames -1;
    int next = Trajectory->frames - 1 - (Trajectory->frames-1)%Trajectory->stride;
    //printf("starting with prev - next = %d %d\n", prev, next);
    for (cl_id = Trajectory->frames -2; cl_id > 0; cl_id--){
        if (Trajectory->strides[cl_id] == 0){
            //printf("passing rmsds %d, %d %lf and %lf\n", cl_id*2, cl_id*2+1,prev_align->rmsd_vector[cl_id*2], prev_align->rmsd_vector[cl_id*2+1]);
            new_align->rmsd_vector[cl_id*2+1] = correct_rmsd_two_frames(Trajectory, prev_align->rotation_matrices_vector[cl_id*2+1], prev_align->coms[prev], prev_align->coms[cl_id], cgnum, removed, added, prev, cl_id, prev_align->rmsd_vector[cl_id*2+1]);
            new_align->rmsd_vector[cl_id*2] = correct_rmsd_two_frames(Trajectory, prev_align->rotation_matrices_vector[cl_id*2], prev_align->coms[next], prev_align->coms[cl_id], cgnum, removed, added, next, cl_id, prev_align->rmsd_vector[cl_id*2]);
        }
        else{
            prev = next;
            next = prev - Trajectory->stride;
        }
    }
}

void align_traj_to_reference(traj *Trajectory,int ref_id){
    /**
    * routine that aligns the trajectory to a reference frame
    * 
    * Parameters
    * ---------- 
    * `Trajectory` : traj object
    * 
    * `ref_id` : reference frame
    */
    int i, j, f;
    double **x, **y;
    double u[3][3], cm_ref[3], cm2[3];
    zero_vec_d(cm_ref, 3);
    double rmsd, nx, ny, nz;
    for (i = 0; i < Trajectory->n_at; i++) {
            for (j = 0; j < 3; j++) {
            cm_ref[j] += Trajectory->traj_coords[ref_id][i * 3 + j] / Trajectory->n_at;
        }
    }
    //printf("cmref = %lf %lf %lf\n",cm_ref[0],cm_ref[1],cm_ref[2]);
    x = d2t(Trajectory->n_at, 3);
    y = d2t(Trajectory->n_at, 3);
    for (i = 0; i < Trajectory->n_at; i++) {
        for (j = 0; j < 3; j++) {
            x[i][j] = Trajectory->traj_coords[ref_id][i * 3 + j] - cm_ref[j];
            Trajectory->traj_coords[ref_id][i * 3 + j] = Trajectory->traj_coords[ref_id][i * 3 + j] - cm_ref[j];
        }
    }
    // cycling over the other frames
    for (f = 0; f<Trajectory->frames; f++){
        if (f != ref_id){
            zero_vec_d(cm2, 3);
            for (i = 0; i < Trajectory->n_at; i++) {
                for (j = 0; j < 3; j++) {
                    cm2[j] += Trajectory->traj_coords[f][i * 3 + j] / Trajectory->n_at;
                }
            }
            //printf("cm[%d] = %lf %lf %lf \n", f, cm2[0], cm2[1], cm2[2]);
            //printf("coms[%d] = %lf %lf %lf \n", ref_id, cm1[0], cm1[1], cm1[2]);
            for (i = 0; i < Trajectory->n_at; i++) {
                    for (j = 0; j < 3; j++) {
                        y[i][j] = Trajectory->traj_coords[f][i * 3 + j] - cm2[j];
                    }
            }
            rmsd = optimal_alignment(y,x, Trajectory->n_at, u);
            //printf("rmsd = %lf\n",rmsd);
            // write new trajectory
            // first: remove com
            for (i = 0; i < Trajectory->n_at; i++) {
                for (j = 0; j < 3; j++) {
                    Trajectory->traj_coords[f][i * 3 + j] = Trajectory->traj_coords[f][i * 3 + j] - cm2[j];
                }
            }
            // second: align to selected frame
            for (i = 0; i < Trajectory->n_at; i++) {
                nx = Trajectory->traj_coords[f][i*3]*u[0][0] + Trajectory->traj_coords[f][i*3 + 1]*u[0][1] + Trajectory->traj_coords[f][i*3 + 2]*u[0][2];
                ny = Trajectory->traj_coords[f][i*3]*u[1][0] + Trajectory->traj_coords[f][i*3 + 1]*u[1][1] + Trajectory->traj_coords[f][i*3 + 2]*u[1][2];
                nz = Trajectory->traj_coords[f][i*3]*u[2][0] + Trajectory->traj_coords[f][i*3 + 1]*u[2][1] + Trajectory->traj_coords[f][i*3 + 2]*u[2][2];
                Trajectory->traj_coords[f][i*3] = nx;
                Trajectory->traj_coords[f][i*3 + 1] = ny;
                Trajectory->traj_coords[f][i*3 + 2] = nz;
            }
        }
    }
    //printf("aligned trajectory\n");
    //for (f = 0; f<Trajectory->frames; f++){
    //    printf("frame %d\n", f);
    //    for (i = 0; i < Trajectory->n_at; i++){
    //        printf("%lf %lf %lf\n",Trajectory->traj_coords[f][i*3],Trajectory->traj_coords[f][i*3+1],Trajectory->traj_coords[f][i*3+2]);
    //    }
    //}
    free_d2t(x);
    free_d2t(y);
}
