#ifndef HDR_GEOM
#define HDR_GEOM

// geometry.c
void vecprod_d(double *a, double *b, double *c); /* c = axb (in 3d). */
double scal_d(double *a, double *b, int dim);    /* prod. scal a.b */
double coseno(double *vec1, double *vec2, int dim);

double norm_d(double *a, int dim);               /* norma vettore */
void normalize_d(double *a, int dim);            /* normalizza vettore */
double dist_d(double *a, double *b, int dim);    /* norma a-b */

double det(double a1, double a2, double a3,double b1, double b2, double b3,double c1, double c2, double c3); /*determinant */
void vec_sum_d(double *a, double *b, double *c, double d, int dim);
void print_vec_d(double *a, int dim);
void zero_vec_d(double *a, int dim); /* zeroing vector of double */
void zero_vec_i(int *a, int dim); /* zeroing vector of int */
void zero_matrix_d(double **a, int dim1, int dim2);

void myjacobi(double a[][3], int n, double *d, double v[][3], int *nrot);

#endif