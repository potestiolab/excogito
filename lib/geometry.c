/**
* \class geometry
* \brief library of functions that perform simple geometrical calculations
*/

#include <stdio.h>
#include <math.h>
#include <geometry.h>

/*******************************/
void vecprod_d(double *a, double *b, double *c) {

    c[0] = a[1] * b[2] - a[2] * b[1];
    c[1] = a[2] * b[0] - a[0] * b[2];
    c[2] = a[0] * b[1] - a[1] * b[0];
}

/*******************************/

double scal_d(double *a, double *b, int dim) {

    int i;
    double temp;

    temp = 0.0;
    for (i = 0; i < dim; i++) {
        temp += a[i] * b[i];
    }
    return (temp);
}

/*******************************/


double coseno(double *vec1, double *vec2, int dim) {

    double temp;

    temp = scal_d(vec1, vec2, dim) / (norm_d(vec1, dim) * norm_d(vec2, dim));

    return (temp);

}

/*******************************/

double norm_d(double *a, int dim) {

    return (sqrt(scal_d(a, a, dim)));
}

/*******************************/
void normalize_d(double *a, int dim) {
    int i;
    double temp;

    temp = norm_d(a, dim);
    for (i = 0; i < dim; i++) {
        a[i] = a[i] / temp;
    }
}


/*******************************/

double dist_d(double *a, double *b, int dim) {

    int i;
    double temp;

    temp = 0.0;
    for (i = 0; i < dim; i++) {
        temp += (a[i] - b[i]) * (a[i] - b[i]);
    }

    temp = sqrt(temp);
    return (temp);
}

/*******************************/

double det(double a1, double a2, double a3,
           double b1, double b2, double b3,
           double c1, double c2, double c3) {

    double temp;

    temp = a1 * (b2 * c3 - b3 * c2);
    temp += a2 * (b3 * c1 - b1 * c3);
    temp += a3 * (b1 * c2 - b2 * c1);
    return (temp);
}

/*******************************/

void vec_sum_d(double *a, double *b, double *c, double d, int dim) {


    int i;
    for (i = 0; i < dim; i++) {
        c[i] = a[i] + d * b[i];
    }
}

/*******************************/

void print_vec_d(double *a, int dim) {

    int i;

    for (i = 0; i < dim; i++) {
        printf("%4d %10.5lf\n", i, a[i]);
    }
}

/*******************************/

void zero_vec_d(double *a, int dim) {

    int i;
    for (i = 0; i < dim; i++) {
        a[i] = 0;
    }
}

/*******************************/

void zero_vec_i(int *a, int dim) {

    int i;
    for (i = 0; i < dim; i++) {
        a[i] = 0;
    }
}

void myjacobi(double a[][3], int n, double *d, double v[][3], int *nrot) {

    int j, iq, ip, i;
    double tresh, theta, tau, t, sm, s, h, g, c;
    double b[3], z[3];

#define ROTATE(a, i, j, k, l) g=a[i][j];h=a[k][l];a[i][j]=g-s*(h+g*tau);a[k][l]=h+s*(g-h*tau);

    for (ip = 0; ip <= (n - 1); ip++) {
        for (iq = 0; iq <= (n - 1); iq++)
            v[ip][iq] = 0.0;
        v[ip][ip] = 1.0;
    }
    for (ip = 0; ip <= (n - 1); ip++) {
        b[ip] = d[ip] = a[ip][ip];
        z[ip] = 0.0;
    }
    *nrot = 0;
    for (i = 1; i <= 500; i++) {
        sm = 0.0;
        for (ip = 0; ip <= n - 2; ip++) {
            for (iq = ip + 1; iq <= (n - 1); iq++)
                sm += fabs(a[ip][iq]);
        }
        if (sm == 0.0) {
            return;
        }
        if (i < 4)
            tresh = 0.2 * sm / (n * n);
        else
            tresh = 0.0;
        for (ip = 0; ip <= n - 2; ip++) {
            for (iq = ip + 1; iq <= (n - 1); iq++) {
                g = 100.0 * fabs(a[ip][iq]);
                if (i > 4 && (fabs((fabs(d[ip]) + g) - fabs(d[ip])) < 1.0e-6)
                    && (fabs((fabs(d[iq]) + g) - fabs(d[iq])) < 1.0e-6))
                    a[ip][iq] = 0.0;
                else if (fabs(a[ip][iq]) > tresh) {
                    h = d[iq] - d[ip];
                    if (fabs((fabs(h) + g) - fabs(h)) < 1.0e-6)
                        t = (a[ip][iq]) / h;
                    else {
                        theta = 0.5 * h / (a[ip][iq]);
                        t = 1.0 / (fabs(theta) + sqrt(1.0 + theta * theta));
                        if (theta < 0.0)
                            t = -t;
                    }
                    c = 1.0 / sqrt(1 + t * t);
                    s = t * c;
                    tau = s / (1.0 + c);
                    h = t * a[ip][iq];
                    z[ip] -= h;
                    z[iq] += h;
                    d[ip] -= h;
                    d[iq] += h;
                    a[ip][iq] = 0.0;
                    for (j = 0; j <= ip - 1; j++) {
                        ROTATE (a, j, ip, j, iq)
                    }
                    for (j = ip + 1; j <= iq - 1; j++) {
                        ROTATE (a, ip, j, j, iq)
                    }
                    for (j = iq + 1; j <= (n - 1); j++) {
                        ROTATE (a, ip, j, iq, j)
                    }
                    for (j = 0; j <= (n - 1); j++) {
                        ROTATE (v, j, ip, j, iq)
                    }
                    ++(*nrot);
                }
            }
        }
        for (ip = 0; ip <= (n - 1); ip++) {
            b[ip] += z[ip];
            d[ip] = b[ip];
            z[ip] = 0.0;
        }
    }
    printf("Too many iterations in routine JACOBI %lf\n", sm);
    /*  exit (1); */
#undef ROTATE
}

void zero_matrix_d(double **a, int dim1, int dim2){
    //Assign 0 values to all the "a" matrix elements
    int i,j;
    for(i=0;i<dim1;i++){
        for(j=0;j<dim2;j++){
            a[i][j]=0;
        }
    }
}