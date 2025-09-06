#include <stdio.h>
#include <math.h>
#ifdef PARALLEL
#include <mpi.h>
#endif
//#include <crtdbg.h>

//#include "memory_debug.h"
#include "macros.h"
#include "pluto.h"

#include "threeDiagonal.h"

void sequentialThreeDiagonalSolverX(double**** x, double**** rightPart, double**** a, double**** b, double**** c, int Nx, int Ny, int Nz, int Nmomentum) {
    //double**** d = create4Darray(Nx, Ny, Nz, Nmomentum);

    for (int j = 0; j < Ny; ++j) {
        for (int k = 0; k < Nz; ++k) {
            for (int l = 0; l < Nmomentum; ++l) {

                double normRightPart = 0;

                for (int i = 0; i < Nx; ++i) {
                    if ((rightPart[k][j][i][l] != rightPart[k][j][i][l]) || (0 * rightPart[k][j][i][l] != 0 * rightPart[k][j][i][l])) {
                        printf("rightPart = NaN in solver X 1, k = %d , j = %d, i = %d, l = %d\n", k, j, i, l);
                        exit(0);
                    }
                    normRightPart = normRightPart + rightPart[k][j][i][l] * rightPart[k][j][i][l];
                }

                if (normRightPart <= 0) {
                    for (int i = 0; i < Nx; ++i) {
                        x[k][j][i][l] = 0;
                    }

                    continue;
                }

                double u = a[k][j][0][l]/b[k][j][0][l];
                double v = c[k][j][Nx - 1][l]/b[k][j][Nx - 1][l];
                for (int i = 0; i < Nx; ++i) {
                    a[k][j][i][l] /= b[k][j][i][l];
                    c[k][j][i][l] /= b[k][j][i][l];
                    rightPart[k][j][i][l] /= b[k][j][i][l];
                    b[k][j][i][l] = 1.0;
                }
                //double* d = new double[Nx];
                //d[0] = u;
                //d[1] = 0;

                for (int i = 2; i < Nx; ++i) {
                    double r = 1.0 / (1 - a[k][j][i][l] * c[k][j][i - 1][l]);
                    //d[i] = -r * a[k][j][i][l] * d[i - 1];
                    rightPart[k][j][i][l] = r * (rightPart[k][j][i][l] - a[k][j][i][l] * rightPart[k][j][i - 1][l]);
                    a[k][j][i][l] = -r * a[k][j][i][l] * a[k][j][i - 1][l];
                    if (i == Nx - 1) {
                        a[k][j][i][l] += v * r;
                    }
                    c[k][j][i][l] = r * c[k][j][i][l];

                }

                for (int i = Nx - 3; i >= 1; i = i - 1) {
                    rightPart[k][j][i][l] = rightPart[k][j][i][l] - rightPart[k][j][i + 1][l] * c[k][j][i][l];
                    a[k][j][i][l] = a[k][j][i][l] - c[k][j][i][l] * a[k][j][i + 1][l];
                    c[k][j][i][l] = -c[k][j][i][l] * c[k][j][i + 1][l];
                }

                double r = 1.0 / (1.0 - a[k][j][1][l] * c[k][j][0][l]);
                rightPart[k][j][0][l] = r * (rightPart[k][j][0][l] - rightPart[k][j][1][l] * c[k][j][0][l]);
                c[k][j][0][l] = r * (u - c[k][j][0][l] * c[k][j][1][l]);
                a[k][j][0][l] = r*a[k][j][0][l];

                double a1 = 1.0;
                double c1 = c[k][j][0][l];
                double d1 = rightPart[k][j][0][l];

                double a2 = a[k][j][Nx - 1][l];
                double c2 = 1.0;
                double d2 = rightPart[k][j][Nx - 1][l];

                double y2 = (d2 - d1 * a2) / (c2 - c1 * a2);
                double y1 = d1 - c1 * y2;

                x[k][j][0][l] = y1;
                x[k][j][Nx - 1][l] = y2;

                if ((x[k][j][0][l] != x[k][j][0][l]) || (0 * x[k][j][0][l] != 0 * x[k][j][0][l])) {
                    printf("x = NaN in solver X, k = %d , j = %d, i = %d, l = %d\n", k, j, 0, l);
                    exit(0);
                }
                if ((x[k][j][Nx-1][l] != x[k][j][Nx-1][l]) || (0 * x[k][j][Nx-1][l] != 0 * x[k][j][Nx-1][l])) {
                    printf("x = NaN in solver X, k = %d , j = %d, i = %d, l = %d\n", k, j, Nx-1, l);
                    exit(0);
                }

                for (int i = 1; i < Nx - 1; ++i) {
                    x[k][j][i][l] = rightPart[k][j][i][l] - a[k][j][i][l] * y1 - c[k][j][i][l] * y2;
                    if ((x[k][j][i][l] != x[k][j][i][l]) || (0 * x[k][j][i][l] != 0 * x[k][j][i][l])) {
                        printf("x = NaN in solver X, k = %d , j = %d, i = %d, l = %d\n", k, j, i, l);
                        exit(0);
                    }
                }

                //delete[] d;
            }
        }
    }
}

void sequentialThreeDiagonalSolverY(double**** x, double**** rightPart, double**** a, double**** b, double**** c, int Nx, int Ny, int Nz, int Nmomentum) {
    //double**** d = create4Darray(Nx, Ny, Nz, Nmomentum);

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    for (int i = 0; i < Nx; ++i) {
        for (int k = 0; k < Nz; ++k) {
            for (int l = 0; l < Nmomentum; ++l) {

                double normRightPart = 0;

                for (int j = 0; j < Ny; ++j) {
                    if ((rightPart[k][j][i][l] != rightPart[k][j][i][l]) || (0 * rightPart[k][j][i][l] != 0 * rightPart[k][j][i][l])) {
                        printf("rightPart = NaN in sequential solver Y 1, k = %d , j = %d, i = %d, l = %d, rank = %d\n", k, j, i, l, rank);
                        exit(0);
                    }
                    normRightPart = normRightPart + rightPart[k][j][i][l] * rightPart[k][j][i][l];
                }

                if (normRightPart <= 0) {
                    for (int j = 0; j < Ny; ++j) {
                        x[k][j][i][l] = 0;
                    }

                    continue;
                }

                double u = a[k][0][i][l]/b[k][0][i][l];
                double v = c[k][Ny - 1][i][l]/b[k][Ny - 1][i][l];
                for (int j = 0; j < Ny; ++j) {
                    a[k][j][i][l] /= b[k][j][i][l];
                    c[k][j][i][l] /= b[k][j][i][l];
                    rightPart[k][j][i][l] /= b[k][j][i][l];
                    b[k][j][i][l] = 1.0;
                }
                //double* d = new double[Ny];
                //d[0] = u;
                //d[1] = 0;

                for (int j = 2; j < Ny; ++j) {
                    double r = 1.0 / (1 - a[k][j][i][l] * c[k][j - 1][i][l]);
                    //d[j] = -r * a[k][j][i][l] * d[j - 1];
                    rightPart[k][j][i][l] = r * (rightPart[k][j][i][l] - a[k][j][i][l] * rightPart[k][j - 1][i][l]);
                    a[k][j][i][l] = -r * a[k][j][i][l] * a[k][j- 1][i][l];
                    if (j == Ny - 1) {
                        a[k][j][i][l] += v * r;
                    }
                    c[k][j][i][l] = r * c[k][j][i][l];

                }

                for (int j = Ny - 3; j >= 1; j = j - 1) {
                    rightPart[k][j][i][l] = rightPart[k][j][i][l] - rightPart[k][j + 1][i][l] * c[k][j][i][l];
                    a[k][j][i][l] = a[k][j][i][l] - c[k][j][i][l] * a[k][j + 1][i][l];
                    c[k][j][i][l] = - c[k][j][i][l] * c[k][j + 1][i][l];
                }

                double r = 1.0 / (1.0 - a[k][1][i][l] * c[k][0][i][l]);
                rightPart[k][0][i][l] = r * (rightPart[k][0][i][l] - rightPart[k][1][i][l] * c[k][0][i][l]);
                c[k][0][i][l] = r * (u - c[k][0][i][l] * c[k][1][i][l]);
                a[k][0][i][l] = r*a[k][0][i][l];

                double a1 = 1.0;
                double c1 = c[k][0][i][l];
                double d1 = rightPart[k][0][i][l];

                double a2 = a[k][Ny - 1][i][l];
                double c2 = 1.0;
                double d2 = rightPart[k][Ny - 1][i][l];

                double y2 = (d2 - d1 * a2) / (c2 - c1 * a2);
                double y1 = d1 - c1 * y2;

                x[k][0][i][l] = y1;
                x[k][Ny - 1][i][l] = y2;

                if ((x[k][0][i][l] != x[k][0][i][l]) || (0 * x[k][0][i][l] != 0 * x[k][0][i][l])) {
                    printf("x = NaN in solver Y, k = %d , j = %d, i = %d, l = %d\n", k, 0, i, l);
                    exit(0);
                }
                if ((x[k][Ny-1][i][l] != x[k][Ny-1][i][l]) || (0 * x[k][Ny-1][i][l] != 0 * x[k][Ny-1][i][l])) {
                    printf("x = NaN in solver Y, k = %d , j = %d, i = %d, l = %d\n", k, Ny-1, i, l);
                    exit(0);
                }

                for (int j = 1; j < Ny - 1; ++j) {
                    x[k][j][i][l] = rightPart[k][j][i][l] - a[k][j][i][l] * y1 - c[k][j][i][l] * y2;
                    if ((x[k][j][i][l] != x[k][j][i][l]) || (0 * x[k][j][i][l] != 0 * x[k][j][i][l])) {
                        printf("x = NaN in solver Y, k = %d , j = %d, i = %d, l = %d\n", k, j, i, l);
                        exit(0);
                    }
                }

                //delete[] d;
            }
        }
    }
}

void sequentialThreeDiagonalSolverZ(double**** x, double**** rightPart, double**** a, double**** b, double**** c, int Nx, int Ny, int Nz, int Nmomentum) {
    //double**** d = create4Darray(Nx, Ny, Nz, Nmomentum);

    for (int i = 0; i < Nx; ++i) {
        for (int j = 0; j < Ny; ++j) {
            for (int l = 0; l < Nmomentum; ++l) {

                double normRightPart = 0;

                for (int k = 0; k < Nz; ++k) {
                    if ((rightPart[k][j][i][l] != rightPart[k][j][i][l]) || (0 * rightPart[k][j][i][l] != 0 * rightPart[k][j][i][l])) {
                        printf("rightPart = NaN in solver Z 1, k = %d , j = %d, i = %d, l = %d\n", k, j, i, l);
                        exit(0);
                    }
                    normRightPart = normRightPart + rightPart[k][j][i][l] * rightPart[k][j][i][l];
                }

                if (normRightPart <= 0) {
                    for (int k = 0; k < Nz; ++k) {
                        x[k][j][i][l] = 0;
                    }

                    continue;
                }

                double u = a[0][j][i][l]/b[0][j][i][l];
                double v = c[Nz - 1][j][i][l]/b[Nz - 1][j][i][l];
                for (int k = 0; k < Nz; ++k) {
                    a[k][j][i][l] /= b[k][j][i][l];
                    c[k][j][i][l] /= b[k][j][i][l];
                    rightPart[k][j][i][l] /= b[k][j][i][l];
                    b[k][j][i][l] = 1.0;
                }
                //double* d = new double[Nz];
                //d[0] = u;
                //d[1] = 0;

                for (int k = 2; k < Nz; ++k) {
                    double r = 1.0 / (1 - a[k][j][i][l] * c[k - 1][j][i][l]);
                    //d[k] = -r * a[k][j][i][l] * d[k - 1];
                    rightPart[k][j][i][l] = r * (rightPart[k][j][i][l] - a[k][j][i][l] * rightPart[k - 1][j][i][l]);
                    a[k][j][i][l] = -r * a[k][j][i][l] * a[k - 1][j][i][l];
                    if (k == Nz - 1) {
                        a[k][j][i][l] += v * r;
                    }
                    c[k][j][i][l] = r * c[k][j][i][l];

                }

                for (int k = Nz - 3; k >= 1; k = k - 1) {
                    rightPart[k][j][i][l] = rightPart[k][j][i][l] - rightPart[k + 1][j][i][l] * c[k][j][i][l];
                    a[k][j][i][l] = a[k][j][i][l] - c[k][j][i][l] * a[k + 1][j][i][l];
                    c[k][j][i][l] =  - c[k][j][i][l] * c[k + 1][j][i][l];
                }

                double r = 1.0 / (1.0 - a[1][j][i][l] * c[0][j][i][l]);
                rightPart[0][j][i][l] = r * (rightPart[0][j][i][l] - rightPart[1][j][i][l] * c[0][j][i][l]);
                c[0][j][i][l] = r * (u - c[0][j][i][l] * c[1][j][i][l]);
                a[0][j][i][l] = r*a[0][j][i][l];

                double a1 = 1.0;
                double c1 = c[0][j][i][l];
                double d1 = rightPart[0][j][i][l];

                double a2 = a[Nz - 1][j][i][l];
                double c2 = 1.0;
                double d2 = rightPart[Nz - 1][j][i][l];

                double y2 = (d2 - d1 * a2) / (c2 - c1 * a2);
                double y1 = d1 - c1 * y2;

                x[0][j][i][l] = y1;
                x[Nz - 1][j][i][l] = y2;

                if ((x[0][j][i][l] != x[0][j][i][l]) || (0 * x[0][j][i][l] != 0 * x[0][j][i][l])) {
                    printf("x = NaN in solver Z, k = %d , j = %d, i = %d, l = %d\n", 0, j, i, l);
                    exit(0);
                }

                if ((x[Nz-1][j][i][l] != x[Nz-1][j][i][l]) || (0 * x[Nz-1][j][i][l] != 0 * x[Nz-1][j][i][l])) {
                    printf("x = NaN in solver Z, k = %d , j = %d, i = %d, l = %d\n", Nz-1, j, i, l);
                    exit(0);
                }

                for (int k = 1; k < Nz - 1; ++k) {
                    x[k][j][i][l] = rightPart[k][j][i][l] - a[k][j][i][l] * y1 - c[k][j][i][l] * y2;
                    if ((x[k][j][i][l] != x[k][j][i][l]) || (0 * x[k][j][i][l] != 0 * x[k][j][i][l])) {
                        printf("x = NaN in solver Z, k = %d , j = %d, i = %d, l = %d\n", k, j, i, l);
                        exit(0);
                    }
                }

                //delete[] d;
            }
        }
    }
}

void sequentialThreeDiagonalSolverP(double**** x, double**** rightPart, double**** a, double**** b, double**** c, int Nx, int Ny, int Nz, int Nmomentum) {
    //double**** d = create4Darray(Nx, Ny, Nz, Nmomentum);

    for (int j = 0; j < Ny; ++j) {
        for (int k = 0; k < Nz; ++k) {
            for (int i = 0; i < Nx; ++i) {

                double normRightPart = 0;

                for (int l = 0; l < Nmomentum; ++l) {
                    normRightPart = normRightPart + rightPart[k][j][i][l] * rightPart[k][j][i][l];
                }

                if (normRightPart <= 0) {
                    for (int l = 0; l < Nmomentum; ++l) {
                        x[k][j][i][l] = 0;
                    }

                    continue;
                }

                double u = a[k][j][i][0] / b[k][j][i][0];
                double v = c[k][j][i][Nmomentum - 1] / b[k][j][i][Nmomentum - 1];
                for (int l = 0; l < Nmomentum; ++l) {
                    a[k][j][i][l] /= b[k][j][i][l];
                    c[k][j][i][l] /= b[k][j][i][l];
                    rightPart[k][j][i][l] /= b[k][j][i][l];
                    b[k][j][i][l] = 1.0;

                    if ((a[k][j][i][l] != a[k][j][i][l]) || (0 * a[k][j][i][l] != 0 * a[k][j][i][l])) {
                        printf("a = NaN k = %d, j = %d, i = %d, l = %d\n", k, j, i, l);
                        exit(0);
                    }

                    if ((c[k][j][i][l] != c[k][j][i][l]) || (0 * c[k][j][i][l] != 0 * c[k][j][i][l])) {
                        printf("c = NaN k = %d, j = %d, i = %d, l = %d\n", k, j, i, l);
                        exit(0);
                    }

                    if ((rightPart[k][j][i][l] != rightPart[k][j][i][l]) || (0 * rightPart[k][j][i][l] != 0 * rightPart[k][j][i][l])){
                        printf("rightPart = NaN k = %d, j = %d, i = %d, l = %d\n", k, j, i, l);
                        exit(0);
                    }
                }
                //double* d = new double[Nx];
                //d[0] = u;
                //d[1] = 0;

                for (int l = 2; l < Nmomentum; ++l) {
                    double r = 1.0 / (1 - a[k][j][i][l] * c[k][j][i][l-1]);
                    //d[i] = -r * a[k][j][i][l] * d[i - 1];
                    rightPart[k][j][i][l] = r * (rightPart[k][j][i][l] - a[k][j][i][l] * rightPart[k][j][i][l-1]);
                    a[k][j][i][l] = -r * a[k][j][i][l] * a[k][j][i][l-1];
                    if (l == Nmomentum - 1) {
                        a[k][j][i][l] += v * r;
                    }
                    c[k][j][i][l] = r * c[k][j][i][l];

                }

                for (int l = Nmomentum - 3; l >= 1; l = l - 1) {
                    rightPart[k][j][i][l] = rightPart[k][j][i][l] - rightPart[k][j][i][l+1] * c[k][j][i][l];
                    a[k][j][i][l] = a[k][j][i][l] - c[k][j][i][l] * a[k][j][i][l+1];
                    c[k][j][i][l] = -c[k][j][i][l] * c[k][j][i][l+1];
                }

                double y1;
                double y2;

                if (a[k][j][i][1] * c[k][j][i][0] == 1.0) {
                    rightPart[k][j][i][0] = rightPart[k][j][i][0] - rightPart[k][j][i][1] * c[k][j][i][0];
                    c[k][j][i][0] = u - c[k][j][i][0] * c[k][j][i][1];
                    double a1 = 0.0;
                    double c1 = c[k][j][i][0];
                    double d1 = rightPart[k][j][i][0];

                    double a2 = a[k][j][i][Nmomentum - 1];
                    double c2 = 1.0;
                    double d2 = rightPart[k][j][i][Nmomentum - 1];

                    y2 = (d2 - d1 * a2) / (c2 - c1 * a2);
                    y1 = d1 - c1 * y2;

                    x[k][j][i][0] = y1;
                    x[k][j][i][Nmomentum - 1] = y2;
                    if ((x[k][j][i][0] != x[k][j][i][0]) || (0 * x[k][j][i][0] != 0 * x[k][j][i][0])) {
                        printf("x = NaN in solver P, k = %d , j = %d, i = %d, l = %d\n", k, j, i, 0);
                        exit(0);
                    }

                    if ((x[k][j][i][Nmomentum-1] != x[k][j][i][Nmomentum-1]) || (0 * x[k][j][i][Nmomentum-1] != 0 * x[k][j][i][Nmomentum-1])) {
                        printf("x = NaN in solver P, k = %d , j = %d, i = %d, l = %d\n", k, j, i, 0);
                        exit(0);
                    }
                }
                else {

                    double r = 1.0 / (1.0 - a[k][j][i][1] * c[k][j][i][0]);
                    rightPart[k][j][i][0] = r * (rightPart[k][j][i][0] - rightPart[k][j][i][1] * c[k][j][i][0]);
                    c[k][j][i][0] = r * (u - c[k][j][i][0] * c[k][j][i][1]);
                    a[k][j][i][0] = r * a[k][j][i][0];

                    double a1 = 1.0;
                    double c1 = c[k][j][i][0];
                    double d1 = rightPart[k][j][i][0];

                    double a2 = a[k][j][i][Nmomentum - 1];
                    double c2 = 1.0;
                    double d2 = rightPart[k][j][i][Nmomentum - 1];

                    y2 = (d2 - d1 * a2) / (c2 - c1 * a2);
                    y1 = d1 - c1 * y2;

                    x[k][j][i][0] = y1;
                    x[k][j][i][Nmomentum - 1] = y2;

                    if ((x[k][j][i][0] != x[k][j][i][0]) || (0 * x[k][j][i][0] != 0 * x[k][j][i][0])) {
                        printf("x = NaN in solver P, k = %d , j = %d, i = %d, l = %d\n", k, j, i, 0);
                        exit(0);
                    }

                    if ((x[k][j][i][Nmomentum-1] != x[k][j][i][Nmomentum-1]) || (0 * x[k][j][i][Nmomentum-1] != 0 * x[k][j][i][Nmomentum-1])) {
                        printf("x = NaN in solver P, k = %d , j = %d, i = %d, l = %d\n", k, j, i, 0);
                        exit(0);
                    }
                }


                for (int l = 1; l < Nmomentum-1; ++l) {
                    x[k][j][i][l] = rightPart[k][j][i][l] - a[k][j][i][l] * y1 - c[k][j][i][l] * y2;
                    if ((x[k][j][i][l] != x[k][j][i][l]) || (0 * x[k][j][i][l] != 0 * x[k][j][i][l])) {
                        printf("x = NaN in solver P, k = %d , j = %d, i = %d, l = %d\n", k, j, i, l);
                        exit(0);
                    }
                }

                //delete[] d;
            }
        }
    }
}

void noparallelThreeDiagonalSolverX(double**** x, double**** rightPart, double**** a, double**** b, double**** c, int Nmomentum){
    double**** noghostRightPart = (double****) Array4D(NX3, NX2, NX1, Nmomentum, sizeof(double));
    double**** noghostA = (double****) Array4D(NX3, NX2, NX1, Nmomentum, sizeof(double));
    double**** noghostB = (double****) Array4D(NX3, NX2, NX1, Nmomentum, sizeof(double));
    double**** noghostC = (double****) Array4D(NX3, NX2, NX1, Nmomentum, sizeof(double));
    double**** noghostX = (double****) Array4D(NX3, NX2, NX1, Nmomentum, sizeof(double));

    int k, j, i;

    DOM_LOOP(k,j,i){
        for(int l = 0; l < Nmomentum; ++l){
            noghostRightPart[k-KBEG][j-JBEG][i-IBEG][l] = rightPart[k][j][i][l];
            noghostA[k-KBEG][j-JBEG][i-IBEG][l] = a[k][j][i][l];
            noghostB[k-KBEG][j-JBEG][i-IBEG][l] = b[k][j][i][l];
            noghostC[k-KBEG][j-JBEG][i-IBEG][l] = c[k][j][i][l];
        }
    }

    sequentialThreeDiagonalSolverX(noghostX, noghostRightPart, noghostA, noghostB, noghostC, NX1, NX2, NX3, Nmomentum);

    DOM_LOOP(k,j,i){
        for(int l = 0; l < Nmomentum; ++l){
            x[k][j][i][l] = noghostX[k-KBEG][j-JBEG][i-IBEG][l];
        }
    }

    FreeArray4D(noghostRightPart);
    FreeArray4D(noghostX);
    FreeArray4D(noghostA);
    FreeArray4D(noghostB);
    FreeArray4D(noghostC);
}

void noparallelThreeDiagonalSolverY(double**** x, double**** rightPart, double**** a, double**** b, double**** c, int Nmomentum){
    double**** noghostRightPart = (double****) Array4D(NX3, NX2, NX1, Nmomentum, sizeof(double));
    double**** noghostA = (double****) Array4D(NX3, NX2, NX1, Nmomentum, sizeof(double));
    double**** noghostB = (double****) Array4D(NX3, NX2, NX1, Nmomentum, sizeof(double));
    double**** noghostC = (double****) Array4D(NX3, NX2, NX1, Nmomentum, sizeof(double));
    double**** noghostX = (double****) Array4D(NX3, NX2, NX1, Nmomentum, sizeof(double));

    int k, j, i;

    DOM_LOOP(k,j,i){
        for(int l = 0; l < Nmomentum; ++l){
            noghostRightPart[k-KBEG][j-JBEG][i-IBEG][l] = rightPart[k][j][i][l];
            noghostA[k-KBEG][j-JBEG][i-IBEG][l] = a[k][j][i][l];
            noghostB[k-KBEG][j-JBEG][i-IBEG][l] = b[k][j][i][l];
            noghostC[k-KBEG][j-JBEG][i-IBEG][l] = c[k][j][i][l];
        }
    }

    sequentialThreeDiagonalSolverY(noghostX, noghostRightPart, noghostA, noghostB, noghostC, NX1, NX2, NX3, Nmomentum);

    DOM_LOOP(k,j,i){
        for(int l = 0; l < Nmomentum; ++l){
            x[k][j][i][l] = noghostX[k-KBEG][j-JBEG][i-IBEG][l];
        }
    }

    FreeArray4D(noghostRightPart);
    FreeArray4D(noghostX);
    FreeArray4D(noghostA);
    FreeArray4D(noghostB);
    FreeArray4D(noghostC);
}

void noparallelThreeDiagonalSolverZ(double**** x, double**** rightPart, double**** a, double**** b, double**** c, int Nmomentum){
    double**** noghostRightPart = (double****) Array4D(NX3, NX2, NX1, Nmomentum, sizeof(double));
    double**** noghostA = (double****) Array4D(NX3, NX2, NX1, Nmomentum, sizeof(double));
    double**** noghostB = (double****) Array4D(NX3, NX2, NX1, Nmomentum, sizeof(double));
    double**** noghostC = (double****) Array4D(NX3, NX2, NX1, Nmomentum, sizeof(double));
    double**** noghostX = (double****) Array4D(NX3, NX2, NX1, Nmomentum, sizeof(double));

    int k, j, i;

    DOM_LOOP(k,j,i){
        for(int l = 0; l < Nmomentum; ++l){
            noghostRightPart[k-KBEG][j-JBEG][i-IBEG][l] = rightPart[k][j][i][l];
            noghostA[k-KBEG][j-JBEG][i-IBEG][l] = a[k][j][i][l];
            noghostB[k-KBEG][j-JBEG][i-IBEG][l] = b[k][j][i][l];
            noghostC[k-KBEG][j-JBEG][i-IBEG][l] = c[k][j][i][l];
        }
    }

    sequentialThreeDiagonalSolverZ(noghostX, noghostRightPart, noghostA, noghostB, noghostC, NX1, NX2, NX3, Nmomentum);

    DOM_LOOP(k,j,i){
        for(int l = 0; l < Nmomentum; ++l){
            x[k][j][i][l] = noghostX[k-KBEG][j-JBEG][i-IBEG][l];
        }
    }

    FreeArray4D(noghostRightPart);
    FreeArray4D(noghostX);
    FreeArray4D(noghostA);
    FreeArray4D(noghostB);
    FreeArray4D(noghostC);
}

void noparallelThreeDiagonalSolverP(double**** x, double**** rightPart, double**** a, double**** b, double**** c, int Nmomentum){
    double**** noghostRightPart = (double****) Array4D(NX3, NX2, NX1, Nmomentum, sizeof(double));
    double**** noghostA = (double****) Array4D(NX3, NX2, NX1, Nmomentum, sizeof(double));
    double**** noghostB = (double****) Array4D(NX3, NX2, NX1, Nmomentum, sizeof(double));
    double**** noghostC = (double****) Array4D(NX3, NX2, NX1, Nmomentum, sizeof(double));
    double**** noghostX = (double****) Array4D(NX3, NX2, NX1, Nmomentum, sizeof(double));

    int k, j, i;

    DOM_LOOP(k,j,i){
        for(int l = 0; l < Nmomentum; ++l){
            noghostRightPart[k-KBEG][j-JBEG][i-IBEG][l] = rightPart[k][j][i][l];
            noghostA[k-KBEG][j-JBEG][i-IBEG][l] = a[k][j][i][l];
            noghostB[k-KBEG][j-JBEG][i-IBEG][l] = b[k][j][i][l];
            noghostC[k-KBEG][j-JBEG][i-IBEG][l] = c[k][j][i][l];
        }
    }

    sequentialThreeDiagonalSolverP(noghostX, noghostRightPart, noghostA, noghostB, noghostC, NX1, NX2, NX3, Nmomentum);

    DOM_LOOP(k,j,i){
        for(int l = 0; l < Nmomentum; ++l){
            x[k][j][i][l] = noghostX[k-KBEG][j-JBEG][i-IBEG][l];
        }
    }

    FreeArray4D(noghostRightPart);
    FreeArray4D(noghostX);
    FreeArray4D(noghostA);
    FreeArray4D(noghostB);
    FreeArray4D(noghostC);
}

#ifdef PARALLEL

void parallelThreeDiagonalSolverX(double**** x, double**** rightPart, double**** a, double**** b, double**** c, int Nmomentum, int Nprocs, int rank, MPI_Comm comm) {
    double**** parallelRightPart = (double****) Array4D(NX3, NX2, 2*Nprocs, Nmomentum, sizeof(double));
    double**** parallelA = (double****) Array4D(NX3, NX2, 2*Nprocs, Nmomentum, sizeof(double));
    double**** parallelB = (double****) Array4D(NX3, NX2, 2*Nprocs, Nmomentum, sizeof(double));
    double**** parallelC = (double****) Array4D(NX3, NX2, 2*Nprocs, Nmomentum, sizeof(double));
    double**** parallelX = (double****) Array4D(NX3, NX2, 2*Nprocs, Nmomentum, sizeof(double));

    double* outcoef = (double*) malloc(NX2*NX3*Nmomentum*6*sizeof(double));
    double* incoef = (double*) malloc(NX2*NX3*Nmomentum*6*Nprocs*sizeof(double));

    int i, j, k;
    int outCount = 0;
    int inCount = 0;

    JDOM_LOOP(j) {
        KDOM_LOOP(k) {
            for (int l = 0; l < Nmomentum; ++l) {

                double normRightPart = 0;

                IDOM_LOOP(i) {
                    normRightPart = normRightPart + rightPart[k][j][i][l] * rightPart[k][j][i][l];
                }
                double norm[1];
                double temp[1];
                norm[0] = normRightPart;
                MPI_Allreduce(norm, temp, 1, MPI_DOUBLE, MPI_SUM, comm);
                normRightPart = temp[0];

                if (normRightPart <= 0) {
                    IDOM_LOOP(i) {
                        x[k][j][i][l] = 0;
                    }
                    outcoef[outCount] = 0;
                    outCount++;
                    outcoef[outCount] = 0;
                    outCount++;
                    outcoef[outCount] = 0;
                    outCount++;
                    outcoef[outCount] = 0;
                    outCount++;
                    outcoef[outCount] = 0;
                    outCount++;
                    outcoef[outCount] = 0;
                    outCount++;
                    continue;
                }

                //double u = a[k][j][0][l]/b[k][j][0][l];
                //double v = c[k][j][Nx - 1][l]/b[k][j][Nx - 1][l];
                IDOM_LOOP(i) {
                    a[k][j][i][l] /= b[k][j][i][l];
                    c[k][j][i][l] /= b[k][j][i][l];
                    rightPart[k][j][i][l] /= b[k][j][i][l];
                    b[k][j][i][l] = 1.0;
                }
                //double* d = new double[Nx];
                //d[0] = u;
                //d[1] = 0;

                for (int i = IBEG + 2; i <= IEND; ++i) {
                    double r = 1.0 / (1 - a[k][j][i][l] * c[k][j][i - 1][l]);
                    //d[i] = -r * a[k][j][i][l] * d[i - 1];
                    rightPart[k][j][i][l] = r * (rightPart[k][j][i][l] - a[k][j][i][l] * rightPart[k][j][i - 1][l]);
                    if((rightPart[k][j][i][l] != rightPart[k][j][i][l]) || (0*rightPart[k][j][i][l] != 0*rightPart[k][j][i][l])){
                        printf("rightPart = NaN in parallel solver X 1, k = %d, j = %d, i = %d, l = %d\n", k,j,i,l);
                        exit(0);
                    }
                    a[k][j][i][l] = -r * a[k][j][i][l] * a[k][j][i - 1][l];

                    c[k][j][i][l] = r * c[k][j][i][l];

                }

                for (int i = IEND - 2; i >= IBEG + 1; i = i - 1) {
                    rightPart[k][j][i][l] = rightPart[k][j][i][l] - rightPart[k][j][i + 1][l] * c[k][j][i][l];
                    if((rightPart[k][j][i][l] != rightPart[k][j][i][l]) || (0*rightPart[k][j][i][l] != 0*rightPart[k][j][i][l])){
                        printf("rightPart = NaN in parallel solver X 2, k = %d, j = %d, i = %d, l = %d\n", k,j,i,l);
                        exit(0);
                    }
                    a[k][j][i][l] = a[k][j][i][l] - c[k][j][i][l] * a[k][j][i + 1][l];
                    c[k][j][i][l] = -c[k][j][i][l] * c[k][j][i + 1][l];
                }

                double r = 1.0 / (1.0 - a[k][j][IBEG+1][l] * c[k][j][IBEG][l]);
                rightPart[k][j][IBEG][l] = r * (rightPart[k][j][IBEG][l] - rightPart[k][j][IBEG+1][l] * c[k][j][IBEG][l]);
                if((rightPart[k][j][IBEG][l] != rightPart[k][j][IBEG][l]) || (0*rightPart[k][j][IBEG][l] != 0*rightPart[k][j][IBEG][l])){
                    printf("rightPart = NaN in parallel solver X, 3 k = %d, j = %d, i = %d, l = %d\n", k,j,IBEG,l);
                    exit(0);
                }
                //c[k][j][0][l] = r * (u - c[k][j][0][l] * c[k][j][1][l]);
                c[k][j][IBEG][l] = - r * c[k][j][IBEG][l] * c[k][j][IBEG+1][l];
                a[k][j][IBEG][l] = r*a[k][j][IBEG][l];


                //double* outcoef = (double*) malloc(6*sizeof(double));
                outcoef[outCount] = a[k][j][IBEG][l];
                if((outcoef[outCount] != outcoef[outCount]) || (0*outcoef[outCount] != 0*outcoef[outCount])){
                    printf("outcoef = NaN 1, outCount = %d, rank = %d\n", outCount, rank);
                    exit(0);
                }
                outCount++;
                outcoef[outCount] = c[k][j][IBEG][l];
                if((outcoef[outCount] != outcoef[outCount]) || (0*outcoef[outCount] != 0*outcoef[outCount])){
                    printf("outcoef = NaN 1, outCount = %d, rank = %d\n", outCount, rank);
                    exit(0);
                }
                outCount++;
                outcoef[outCount] = rightPart[k][j][IBEG][l];
                if((outcoef[outCount] != outcoef[outCount]) || (0*outcoef[outCount] != 0*outcoef[outCount])){
                    printf("outcoef = NaN 1, outCount = %d, rank = %d\n", outCount, rank);
                    exit(0);
                }
                outCount++;
                if((rightPart[k][j][IBEG][l] != rightPart[k][j][IBEG][l]) || (0*rightPart[k][j][IBEG][l] != 0*rightPart[k][j][IBEG][l])){
                    printf("rightPart = NaN in parallel solver X 4, k = %d, j = %d, i = %d, l = %d\n", k,j,IBEG,l);
                    exit(0);
                }
                outcoef[outCount] = a[k][j][IEND][l];
                if((outcoef[outCount] != outcoef[outCount]) || (0*outcoef[outCount] != 0*outcoef[outCount])){
                    printf("outcoef = NaN 1, outCount = %d, rank = %d\n", outCount, rank);
                    exit(0);
                }
                outCount++;
                outcoef[outCount] = c[k][j][IEND][l];
                if((outcoef[outCount] != outcoef[outCount]) || (0*outcoef[outCount] != 0*outcoef[outCount])){
                    printf("outcoef = NaN 1, outCount = %d, rank = %d\n", outCount, rank);
                    exit(0);
                }
                outCount++;
                outcoef[outCount] = rightPart[k][j][IEND][l];
                if((outcoef[outCount] != outcoef[outCount]) || (0*outcoef[outCount] != 0*outcoef[outCount])){
                    printf("outcoef = NaN 1, outCount = %d, rank = %d\n", outCount, rank);
                    exit(0);
                }
                outCount++;
                if((rightPart[k][j][IEND][l] != rightPart[k][j][IEND][l]) || (0*rightPart[k][j][IEND][l] != 0*rightPart[k][j][IEND][l])){
                    printf("rightPart = NaN in parallel solver X 5, k = %d, j = %d, i = %d, l = %d\n", k,j,IEND,l);
                    exit(0);
                }

                //double* incoef = (double*) malloc(6*Nprocs*sizeof(double));
            }
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);

    MPI_Gather(outcoef, 6*NX2*NX3*Nmomentum, MPI_DOUBLE, incoef, 6*NX2*NX3*Nmomentum, MPI_DOUBLE, 0, comm);

    for(int m = 0; m < Nprocs; ++m){
    JDOM_LOOP(j) {
        KDOM_LOOP(k) {
            for (int l = 0; l < Nmomentum; ++l) {
                if(rank == 0){

                        /*parallelA[k-KBEG][j-JBEG][2*m][l] = incoef[6*m];
                        parallelC[k-KBEG][j-JBEG][2*m][l] = incoef[6*m+1];
                        parallelRightPart[k-KBEG][j-JBEG][2*m][l] = incoef[6*m+2];
                        parallelB[k-KBEG][j-JBEG][2*m][l] = 1.0;
                        parallelA[k-KBEG][j-JBEG][2*m+1][l] = incoef[6*m + 3];
                        parallelC[k-KBEG][j-JBEG][2*m+1][l] = incoef[6*m + 4];
                        parallelRightPart[k-KBEG][j-JBEG][2*m+1][l] = incoef[6*m + 5];
                        parallelB[k-KBEG][j-JBEG][2*m+1][l] = 1.0;*/

                    parallelA[k-KBEG][j-JBEG][2*m][l] = incoef[inCount];
                    if((incoef[inCount] != incoef[inCount]) || (0*incoef[inCount] != 0*incoef[inCount])){
                        printf("incoef = NaN, inCount = %d, rank = %d, m = %d\n", inCount, rank, m);
                        exit(0);
                    }
                    inCount++;
                    parallelC[k-KBEG][j-JBEG][2*m][l] = incoef[inCount];
                    if((incoef[inCount] != incoef[inCount]) || (0*incoef[inCount] != 0*incoef[inCount])){
                        printf("incoef = NaN, inCount = %d, rank = %d, m = %d\n", inCount, rank, m);
                        exit(0);
                    }
                    inCount++;
                    parallelRightPart[k-KBEG][j-JBEG][2*m][l] = incoef[inCount];
                    if((incoef[inCount] != incoef[inCount]) || (0*incoef[inCount] != 0*incoef[inCount])){
                        printf("incoef = NaN, inCount = %d, rank = %d, m = %d\n", inCount, rank, m);
                        exit(0);
                    }
                    inCount++;
                    if((parallelRightPart[k-KBEG][j-JBEG][2*m][l] != parallelRightPart[k-KBEG][j-JBEG][2*m][l]) || (0*parallelRightPart[k-KBEG][j-JBEG][2*m][l] != 0*parallelRightPart[k-KBEG][j-JBEG][2*m][l])){
                        printf("parallel right part = NaN k = %d, j = %d, i = %d, l = %d, rank = %d\n ", k-KBEG, j-JBEG, 2*m, l, rank);
                        exit(0);
                    }
                    parallelB[k-KBEG][j-JBEG][2*m][l] = 1.0;
                    parallelA[k-KBEG][j-JBEG][2*m+1][l] = incoef[inCount];
                    if((incoef[inCount] != incoef[inCount]) || (0*incoef[inCount] != 0*incoef[inCount])){
                        printf("incoef = NaN, inCount = %d, rank = %d, m = %d\n", inCount, rank, m);
                        exit(0);
                    }
                    inCount++;
                    parallelC[k-KBEG][j-JBEG][2*m+1][l] = incoef[inCount];
                    if((incoef[inCount] != incoef[inCount]) || (0*incoef[inCount] != 0*incoef[inCount])){
                        printf("incoef = NaN, inCount = %d, rank = %d, m = %d\n", inCount, rank, m);
                        exit(0);
                    }
                    inCount++;
                    parallelRightPart[k-KBEG][j-JBEG][2*m+1][l] = incoef[inCount];
                    if((incoef[inCount] != incoef[inCount]) || (0*incoef[inCount] != 0*incoef[inCount])){
                        printf("incoef = NaN, inCount = %d, rank = %d, m = %d\n", inCount, rank, m);
                        exit(0);
                    }
                    inCount++;
                    if((parallelRightPart[k-KBEG][j-JBEG][2*m+1][l] != parallelRightPart[k-KBEG][j-JBEG][2*m+1][l]) || (0*parallelRightPart[k-KBEG][j-JBEG][2*m+1][l] != 0*parallelRightPart[k-KBEG][j-JBEG][2*m+1][l])){
                        printf("parallel right part = NaN k = %d, j = %d, i = %d, l = %d, rank = %d\n", k-KBEG, j-JBEG, 2*m+1, l, rank);
                        exit(0);
                    }
                    parallelB[k-KBEG][j-JBEG][2*m+1][l] = 1.0;
                    }
                }


            }
        }
    }

    free(incoef);
    free(outcoef);

    MPI_Barrier(comm);

    if(rank == 0){
        sequentialThreeDiagonalSolverX(parallelX, parallelRightPart, parallelA, parallelB, parallelC, 2*Nprocs, NX2, NX3, Nmomentum);
        double* send = (double*) malloc(2*NX2*NX3*Nmomentum*sizeof(double));
        for(int m = 1; m < Nprocs; ++m){
            int count = 0;
            KDOM_LOOP(k){
                JDOM_LOOP(j){
                    for(int l = 0; l < Nmomentum; ++l){
                        send[count] = parallelX[k-KBEG][j-JBEG][2*m][l];
                        count = count + 1;
                        send[count] = parallelX[k-KBEG][j-JBEG][2*m+1][l];
                        count = count + 1;
                    }
                }
            }
            MPI_Send(send, 2*NX2*NX3*Nmomentum, MPI_DOUBLE, m, 0, comm);
        }
        free(send);
        KDOM_LOOP(k){
            JDOM_LOOP(j){
                for(int l = 0; l < Nmomentum; ++l){
                    x[k][j][IBEG][l] = parallelX[k-KBEG][j-JBEG][0][l];
                    if ((x[k][j][IBEG][l] != x[k][j][IBEG][l]) || (0 * x[k][j][IBEG][l] != 0 * x[k][j][IBEG][l])) {
                        printf("x = NaN in solver X 1, k = %d , j = %d, i = %d, l = %d, rank = %d\n", k, j, IBEG, l, rank);
                        exit(0);
                    }
                    x[k][j][IEND][l] = parallelX[k-KBEG][j-JBEG][1][l];
                    if ((x[k][j][IEND][l] != x[k][j][IEND][l]) || (0 * x[k][j][IEND][l] != 0 * x[k][j][IEND][l])) {
                        printf("x = NaN in solver X 2, k = %d , j = %d, i = %d, l = %d, rank = %d\n", k, j, IEND, l, rank);
                        exit(0);
                    }
                }
            }
        }
    } else {
        double* recv = (double*) malloc(2*NX2*NX3*Nmomentum*sizeof(double));
        MPI_Status status;
        MPI_Recv(recv, 2*NX2*NX3*Nmomentum, MPI_DOUBLE, 0, 0, comm, &status);
        int count = 0;
        KDOM_LOOP(k){
            JDOM_LOOP(j){
                for(int l = 0; l < Nmomentum; ++l){
                    x[k][j][IBEG][l] = recv[count];
                    if ((x[k][j][IBEG][l] != x[k][j][IBEG][l]) || (0 * x[k][j][IBEG][l] != 0 * x[k][j][IBEG][l])) {
                        printf("x = NaN in solver X 3, k = %d , j = %d, i = %d, l = %d, rank = %d\n", k, j, IBEG, l, rank);
                        exit(0);
                    }
                    count = count + 1;
                    x[k][j][IEND][l] = recv[count];
                    if ((x[k][j][IEND][l] != x[k][j][IEND][l]) || (0 * x[k][j][IEND][l] != 0 * x[k][j][IEND][l])) {
                        printf("x = NaN in solver X 4, k = %d , j = %d, i = %d, l = %d, rank = %d\n", k, j, IEND, l, rank);
                        exit(0);
                    }
                    count = count + 1;
                }
            }
        }
        free(recv);
    }

    KDOM_LOOP(k){
        JDOM_LOOP(j){
            for(int l = 0; l < Nmomentum; ++l){
                for (int i = IBEG + 1; i < IEND; ++i) {
                    x[k][j][i][l] = rightPart[k][j][i][l] - a[k][j][i][l] * x[k][j][IBEG][l] - c[k][j][i][l] * x[k][j][IEND][l];
                    if ((x[k][j][i][l] != x[k][j][i][l]) || (0 * x[k][j][i][l] != 0 * x[k][j][i][l])) {
                        printf("x = NaN in solver X 5, k = %d , j = %d, i = %d, l = %d, rank = %d\n", k, j, i, l, rank);
                        exit(0);
                    }
                }

                //delete[] d;
            }
        }
    }

    FreeArray4D(parallelA);
    FreeArray4D(parallelB);
    FreeArray4D(parallelC);
    FreeArray4D(parallelRightPart);
    FreeArray4D(parallelX);
}

void parallelThreeDiagonalSolverY(double**** x, double**** rightPart, double**** a, double**** b, double**** c, int Nmomentum, int Nprocs, int rank, MPI_Comm comm) {
    //double**** d = create4Darray(Nx, Ny, Nz, Nmomentum);
    double**** parallelRightPart = (double****) Array4D(NX3, 2*Nprocs, NX1, Nmomentum, sizeof(double));
    double**** parallelA = (double****) Array4D(NX3, 2*Nprocs, NX1, Nmomentum, sizeof(double));
    double**** parallelB = (double****) Array4D(NX3, 2*Nprocs, NX1, Nmomentum, sizeof(double));
    double**** parallelC = (double****) Array4D(NX3, 2*Nprocs, NX1, Nmomentum, sizeof(double));
    double**** parallelX = (double****) Array4D(NX3, 2*Nprocs, NX1, Nmomentum, sizeof(double));

    double* outcoef = (double*) malloc(NX1*NX3*Nmomentum*6*sizeof(double));
    double* incoef = (double*) malloc(NX1*NX3*Nmomentum*6*Nprocs*sizeof(double));

    int i, j, k;
    int outCount = 0;
    int inCount = 0;

    IDOM_LOOP(i){
        KDOM_LOOP(k){
            for (int l = 0; l < Nmomentum; ++l) {

                double normRightPart = 0;

                JDOM_LOOP(j){
                    normRightPart = normRightPart + rightPart[k][j][i][l] * rightPart[k][j][i][l];
                }
                double norm[1];
                double temp[1];
                norm[0] = normRightPart;
                MPI_Allreduce(norm, temp, 1, MPI_DOUBLE, MPI_SUM, comm);
                normRightPart = temp[0];

                if (normRightPart <= 0) {
                    JDOM_LOOP(j) {
                        x[k][j][i][l] = 0;
                    }
                    outcoef[outCount] = 0;
                    outCount++;
                    outcoef[outCount] = 0;
                    outCount++;
                    outcoef[outCount] = 0;
                    outCount++;
                    outcoef[outCount] = 0;
                    outCount++;
                    outcoef[outCount] = 0;
                    outCount++;
                    outcoef[outCount] = 0;
                    outCount++;
                    continue;
                }

                //double u = a[k][j][0][l]/b[k][j][0][l];
                //double v = c[k][j][Nx - 1][l]/b[k][j][Nx - 1][l];
                JDOM_LOOP(j){
                    a[k][j][i][l] /= b[k][j][i][l];
                    c[k][j][i][l] /= b[k][j][i][l];
                    rightPart[k][j][i][l] /= b[k][j][i][l];
                    b[k][j][i][l] = 1.0;
                }
                //double* d = new double[Nx];
                //d[0] = u;
                //d[1] = 0;

                for (j = JBEG + 2; j <= JEND; ++j) {
                    double r = 1.0 / (1 - a[k][j][i][l] * c[k][j-1][i][l]);
                    //d[i] = -r * a[k][j][i][l] * d[i - 1];
                    rightPart[k][j][i][l] = r * (rightPart[k][j][i][l] - a[k][j][i][l] * rightPart[k][j-1][i][l]);
                    a[k][j][i][l] = -r * a[k][j][i][l] * a[k][j-1][i][l];

                    c[k][j][i][l] = r * c[k][j][i][l];

                }

                if((c[k][JEND][i][l] != c[k][JEND][i][l])||(0*c[k][JEND][i][l] != 0*c[k][JEND][i][l])){
                    printf("c[JEND] = NaN\n");
                    exit(0);
                }

                for (int j = JEND - 2; j >= JBEG + 1; j = j - 1) {
                    rightPart[k][j][i][l] = rightPart[k][j][i][l] - rightPart[k][j+1][i][l] * c[k][j][i][l];
                    a[k][j][i][l] = a[k][j][i][l] - c[k][j][i][l] * a[k][j+1][i][l];
                    c[k][j][i][l] = -c[k][j][i][l] * c[k][j+1][i][l];
                }

                double r = 1.0 / (1.0 - a[k][JBEG + 1][i][l] * c[k][JBEG][i][l]);
                rightPart[k][JBEG][i][l] = r * (rightPart[k][JBEG][i][l] - rightPart[k][JBEG + 1][i][l] * c[k][JBEG][i][l]);

                c[k][JBEG][i][l] = - r * c[k][JBEG][i][l] * c[k][JBEG+1][i][l];
                a[k][JBEG][i][l] = r*a[k][JBEG][i][l];


                //double* outcoef = (double*) malloc(6*sizeof(double));
                outcoef[outCount] = a[k][JBEG][i][l];
                if((outcoef[outCount] != outcoef[outCount]) || (0*outcoef[outCount] != 0*outcoef[outCount])){
                    printf("outcoef = NaN 1y, outCount = %d, rank = %d\n", outCount, rank);
                    exit(0);
                }
                outCount++;
                outcoef[outCount] = c[k][JBEG][i][l];
                if((outcoef[outCount] != outcoef[outCount]) || (0*outcoef[outCount] != 0*outcoef[outCount])){
                    printf("outcoef = NaN 2y, outCount = %d, rank = %d\n", outCount, rank);
                    exit(0);
                }
                outCount++;
                outcoef[outCount] = rightPart[k][JBEG][i][l];
                if((outcoef[outCount] != outcoef[outCount]) || (0*outcoef[outCount] != 0*outcoef[outCount])){
                    printf("outcoef = NaN 3y, outCount = %d, rank = %d\n", outCount, rank);
                    exit(0);
                }
                outCount++;
                outcoef[outCount] = a[k][JEND][i][l];
                if((outcoef[outCount] != outcoef[outCount]) || (0*outcoef[outCount] != 0*outcoef[outCount])){
                    printf("outcoef = NaN 4y, outCount = %d, rank = %d\n", outCount, rank);
                    exit(0);
                }
                outCount++;
                outcoef[outCount] = c[k][JEND][i][l];
                if((outcoef[outCount] != outcoef[outCount]) || (0*outcoef[outCount] != 0*outcoef[outCount])){
                    printf("outcoef = NaN 5y, outCount = %d, rank = %d\n", outCount, rank);
                    exit(0);
                }
                outCount++;
                outcoef[outCount] = rightPart[k][JEND][i][l];
                if((outcoef[outCount] != outcoef[outCount]) || (0*outcoef[outCount] != 0*outcoef[outCount])){
                    printf("outcoef = NaN 6y, outCount = %d, rank = %d\n", outCount, rank);
                    exit(0);
                }
                outCount++;

                //double* incoef = (double*) malloc(6*Nprocs*sizeof(double));
            }
        }
    }

                //MPI_Gather(outcoef, 6, MPI_DOUBLE, incoef, 6*Nprocs, MPI_DOUBLE, 0, comm);
    MPI_Gather(outcoef, 6*NX1*NX3*Nmomentum, MPI_DOUBLE, incoef, 6*NX1*NX3*Nmomentum, MPI_DOUBLE, 0, comm);
    for(int m = 0; m < Nprocs; ++m){
    IDOM_LOOP(i){
        KDOM_LOOP(k){
            for (int l = 0; l < Nmomentum; ++l) {
                if(rank == 0){
                        parallelA[k-KBEG][2*m][i-IBEG][l] = incoef[inCount];
                        if((incoef[inCount] != incoef[inCount]) || (0*incoef[inCount] != 0*incoef[inCount])){
                            printf("incoef = NaN 1y, inCount = %d, rank = %d, m = %d\n", inCount, rank, m);
                            exit(0);
                        }
                        inCount++;
                        parallelC[k-KBEG][2*m][i-IBEG][l] = incoef[inCount];
                        if((incoef[inCount] != incoef[inCount]) || (0*incoef[inCount] != 0*incoef[inCount])){
                            printf("incoef = NaN 2y, inCount = %d, rank = %d, m = %d\n", inCount, rank, m);
                            exit(0);
                        }
                        inCount++;
                        parallelRightPart[k-KBEG][2*m][i-IBEG][l] = incoef[inCount];
                        if((incoef[inCount] != incoef[inCount]) || (0*incoef[inCount] != 0*incoef[inCount])){
                            printf("incoef = NaN 3y, inCount = %d, rank = %d, m = %d\n", inCount, rank, m);
                            exit(0);
                        }
                        inCount++;
                        parallelB[k-KBEG][2*m][i-IBEG][l] = 1.0;
                        parallelA[k-KBEG][2*m+1][i-IBEG][l] = incoef[inCount];
                        if((incoef[inCount] != incoef[inCount]) || (0*incoef[inCount] != 0*incoef[inCount])){
                            printf("incoef = NaN 4y, inCount = %d, rank = %d, m = %d\n", inCount, rank, m);
                            exit(0);
                        }
                        inCount++;
                        parallelC[k-KBEG][2*m+1][i-IBEG][l] = incoef[inCount];
                        if((incoef[inCount] != incoef[inCount]) || (0*incoef[inCount] != 0*incoef[inCount])){
                            printf("incoef = NaN 5y, inCount = %d, rank = %d, m = %d, i = %d, k = %d, l = %d\n", inCount, rank, m, i, k, l);
                            exit(0);
                        }
                        inCount++;
                        parallelRightPart[k-KBEG][2*m+1][i-IBEG][l] = incoef[inCount];
                        if((incoef[inCount] != incoef[inCount]) || (0*incoef[inCount] != 0*incoef[inCount])){
                            printf("incoef = NaN 6y, inCount = %d, rank = %d, m = %d\n", inCount, rank, m);
                            exit(0);
                        }
                        inCount++;
                        parallelB[k-KBEG][2*m+1][i-IBEG][l] = 1.0;
                    }
                }


            }
        }
    }

    free(incoef);
    free(outcoef);

    MPI_Barrier(comm);

    if(rank == 0){
        sequentialThreeDiagonalSolverY(parallelX, parallelRightPart, parallelA, parallelB, parallelC, NX1, 2*Nprocs, NX3, Nmomentum);
        double* send = (double*) malloc(2*NX1*NX3*Nmomentum*sizeof(double));
        for(int m = 1; m < Nprocs; ++m){
            int count = 0;
            KDOM_LOOP(k){
                IDOM_LOOP(i){
                    for(int l = 0; l < Nmomentum; ++l){
                        send[count] = parallelX[k-KBEG][2*m][i-IBEG][l];
                        count = count + 1;
                        send[count] = parallelX[k-KBEG][2*m+1][i-IBEG][l];
                        count = count + 1;
                    }
                }
            }
            MPI_Send(send, 2*NX1*NX3*Nmomentum, MPI_DOUBLE, m, 0, comm);
        }
        free(send);
        KDOM_LOOP(k){
            IDOM_LOOP(i){
                for(int l = 0; l < Nmomentum; ++l){
                    x[k][JBEG][i][l] = parallelX[k-KBEG][0][i-IBEG][l];
                    if ((x[k][JBEG][i][l] != x[k][JBEG][i][l]) || (0 * x[k][JBEG][i][l] != 0 * x[k][JBEG][i][l])) {
                        printf("x = NaN in solver Y 1, k = %d , j = %d, i = %d, l = %d, rank = %d\n", k, JBEG, i, l, rank);
                        exit(0);
                    }
                    x[k][JEND][i][l] = parallelX[k-KBEG][1][i-IBEG][l];
                    if ((x[k][JEND][i][l] != x[k][JEND][i][l]) || (0 * x[k][JEND][i][l] != 0 * x[k][JEND][i][l])) {
                        printf("x = NaN in solver Y 2, k = %d , j = %d, i = %d, l = %d, rank = %d\n", k, JEND, i, l, rank);
                        exit(0);
                    }
                }
            }
        }
    } else {
        double* recv = (double*) malloc(2*NX1*NX3*Nmomentum*sizeof(double));
        MPI_Status status;
        MPI_Recv(recv, 2*NX1*NX3*Nmomentum, MPI_DOUBLE, 0, 0, comm, &status);
        int count = 0;
        KDOM_LOOP(k){
            IDOM_LOOP(i){
                for(int l = 0; l < Nmomentum; ++l){
                    x[k][JBEG][i][l] = recv[count];
                    if ((x[k][JBEG][i][l] != x[k][JBEG][i][l]) || (0 * x[k][JBEG][i][l] != 0 * x[k][JBEG][i][l])) {
                        printf("x = NaN in solver Y 3, k = %d , j = %d, i = %d, l = %d, rank = %d\n", k, JBEG, i, l, rank);
                        exit(0);
                    }
                    count = count + 1;
                    x[k][JEND][i][l] = recv[count];
                    if ((x[k][JEND][i][l] != x[k][JEND][i][l]) || (0 * x[k][JEND][i][l] != 0 * x[k][JEND][i][l])) {
                        printf("x = NaN in solver Y 4, k = %d , j = %d, i = %d, l = %d, rank = %d\n", k, JEND, i, l, rank);
                        exit(0);
                    }
                    count = count + 1;
                }
            }
        }
        free(recv);
    }

   KDOM_LOOP(k){
        IDOM_LOOP(i){
            for(int l = 0; l < Nmomentum; ++l){
                for (j = JBEG + 1; j < JEND; ++j) {
                    x[k][j][i][l] = rightPart[k][j][i][l] - a[k][j][i][l] * x[k][JBEG][i][l] - c[k][j][i][l] * x[k][JEND][i][l];
                    if ((x[k][j][i][l] != x[k][j][i][l]) || (0 * x[k][j][i][l] != 0 * x[k][j][i][l])) {
                        printf("x = NaN in solver Y, k = %d , j = %d, i = %d, l = %d, rank = %d\n", k, j, i, l, rank);
                        exit(0);
                    }
                }

                //delete[] d;
            }
        }
    }

    FreeArray4D(parallelA);
    FreeArray4D(parallelB);
    FreeArray4D(parallelC);
    FreeArray4D(parallelRightPart);
    FreeArray4D(parallelX);
}

void parallelThreeDiagonalSolverZ(double**** x, double**** rightPart, double**** a, double**** b, double**** c, int Nmomentum, int Nprocs, int rank, MPI_Comm comm) {
    double**** parallelRightPart = (double****) Array4D(2*Nprocs, NX2, NX1, Nmomentum, sizeof(double));
    double**** parallelA = (double****) Array4D(2*Nprocs, NX2, NX1, Nmomentum, sizeof(double));
    double**** parallelB = (double****) Array4D(2*Nprocs, NX2, NX1, Nmomentum, sizeof(double));
    double**** parallelC = (double****) Array4D(2*Nprocs, NX2, NX1, Nmomentum, sizeof(double));
    double**** parallelX = (double****) Array4D(2*Nprocs, NX2, NX1, Nmomentum, sizeof(double));

    double* outcoef = (double*) malloc(NX2*NX1*Nmomentum*6*sizeof(double));
    double* incoef = (double*) malloc(NX2*NX1*Nmomentum*6*Nprocs*sizeof(double));

    int i, j, k;
    int outCount = 0;
    int inCount = 0;

    JDOM_LOOP(j){
        IDOM_LOOP(i){
            for (int l = 0; l < Nmomentum; ++l) {

                double normRightPart = 0;

                KDOM_LOOP(k){
                    normRightPart = normRightPart + rightPart[k][j][i][l] * rightPart[k][j][i][l];
                }
                double norm[1];
                double temp[1];
                norm[0] = normRightPart;
                MPI_Allreduce(norm, temp, 1, MPI_DOUBLE, MPI_SUM, comm);
                normRightPart = temp[0];

                if (normRightPart <= 0) {
                    KDOM_LOOP(k){
                        x[k][j][i][l] = 0;
                    }
                    outcoef[outCount] = 0;
                    outCount++;
                    outcoef[outCount] = 0;
                    outCount++;
                    outcoef[outCount] = 0;
                    outCount++;
                    outcoef[outCount] = 0;
                    outCount++;
                    outcoef[outCount] = 0;
                    outCount++;
                    outcoef[outCount] = 0;
                    outCount++;
                    continue;
                }

                //double u = a[k][j][0][l]/b[k][j][0][l];
                //double v = c[k][j][Nx - 1][l]/b[k][j][Nx - 1][l];
                KDOM_LOOP(k){
                    a[k][j][i][l] /= b[k][j][i][l];
                    c[k][j][i][l] /= b[k][j][i][l];
                    rightPart[k][j][i][l] /= b[k][j][i][l];
                    b[k][j][i][l] = 1.0;
                }
                //double* d = new double[Nx];
                //d[0] = u;
                //d[1] = 0;

                for (k = KBEG + 2; k <= KEND; ++k) {
                    double r = 1.0 / (1 - a[k][j][i][l] * c[k-1][j][i][l]);
                    //d[i] = -r * a[k][j][i][l] * d[i - 1];
                    rightPart[k][j][i][l] = r * (rightPart[k][j][i][l] - a[k][j][i][l] * rightPart[k-1][j][i][l]);
                    a[k][j][i][l] = -r * a[k][j][i][l] * a[k-1][j][i][l];

                    c[k][j][i][l] = r * c[k][j][i][l];

                }

                for (k = KEND - 2; k >= KBEG + 1; k = k - 1) {
                    rightPart[k][j][i][l] = rightPart[k][j][i][l] - rightPart[k+1][j][i][l] * c[k][j][i][l];
                    a[k][j][i][l] = a[k][j][i][l] - c[k][j][i][l] * a[k+1][j][i][l];
                    c[k][j][i][l] = -c[k][j][i][l] * c[k+1][j][i][l];
                }

                double r = 1.0 / (1.0 - a[KBEG + 1][j][i][l] * c[KBEG][j][i][l]);
                rightPart[KBEG][j][i][l] = r * (rightPart[KBEG][j][i][l] - rightPart[KBEG + 1][j][i][l] * c[KBEG][j][i][l]);

                c[KBEG][j][i][l] = - r * c[KBEG][j][i][l] * c[KBEG+1][j][i][l];
                a[KBEG][j][i][l] = r*a[KBEG][j][i][l];


                //double* outcoef = (double*) malloc(6*sizeof(double));
                outcoef[outCount] = a[KBEG][j][i][l];
                outCount++;
                outcoef[outCount] = c[KBEG][j][i][l];
                outCount++;
                outcoef[outCount] = rightPart[KBEG][j][i][l];
                outCount++;
                outcoef[outCount] = a[KEND][j][i][l];
                outCount++;
                outcoef[outCount] = c[KEND][j][i][l];
                outCount++;
                outcoef[outCount] = rightPart[KEND][j][i][l];
                outCount++;

                //double* incoef = (double*) malloc(6*Nprocs*sizeof(double));
            }
        }
    }

                //MPI_Gather(outcoef, 6, MPI_DOUBLE, incoef, 6*Nprocs, MPI_DOUBLE, 0, comm);
    MPI_Gather(outcoef, 6*NX1*NX2*Nmomentum, MPI_DOUBLE, incoef, 6*NX1*NX2*Nmomentum, MPI_DOUBLE, 0, comm);
    for(int m = 0; m < Nprocs; ++m){
    JDOM_LOOP(j){
        IDOM_LOOP(i){
            for (int l = 0; l < Nmomentum; ++l) {
                if(rank == 0){
                        parallelA[2*m][j-JBEG][i-IBEG][l] = incoef[inCount];
                        inCount++;
                        parallelC[2*m][j-JBEG][i-IBEG][l] = incoef[inCount];
                        inCount++;
                        parallelRightPart[2*m][j-JBEG][i-IBEG][l] = incoef[inCount];
                        inCount++;
                        parallelB[2*m][j-JBEG][i-IBEG][l] = 1.0;
                        parallelA[2*m+1][j-JBEG][i-IBEG][l] = incoef[inCount];
                        inCount++;
                        parallelC[2*m+1][j-JBEG][i-IBEG][l] = incoef[inCount];
                        inCount++;
                        parallelRightPart[2*m+1][j-JBEG][i-IBEG][l] = incoef[inCount];
                        inCount++;
                        parallelB[2*m+1][j-JBEG][i-IBEG][l] = 1.0;
                    }
                }
            }
        }
    }

    free(incoef);
    free(outcoef);

    MPI_Barrier(comm);

    if(rank == 0){
        sequentialThreeDiagonalSolverZ(parallelX, parallelRightPart, parallelA, parallelB, parallelC, NX1, NX2, 2*Nprocs, Nmomentum);
        double* send = (double*) malloc(2*NX2*NX1*Nmomentum*sizeof(double));
        for(int m = 1; m < Nprocs; ++m){
            int count = 0;
            IDOM_LOOP(i){
                JDOM_LOOP(j){
                    for(int l = 0; l < Nmomentum; ++l){
                        send[count] = parallelX[2*m][j-JBEG][i-IBEG][l];
                        count = count + 1;
                        send[count] = parallelX[2*m + 1][j-JBEG][i-IBEG][l];
                        count = count + 1;
                    }
                }
            }
            MPI_Send(send, 2*NX2*NX1*Nmomentum, MPI_DOUBLE, m, 0, comm);
        }
        free(send);
        IDOM_LOOP(i){
            JDOM_LOOP(j){
                for(int l = 0; l < Nmomentum; ++l){
                    x[KBEG][j][i][l] = parallelX[0][j-JBEG][i-IBEG][l];
                    x[KEND][j][i][l] = parallelX[1][j-JBEG][i-IBEG][l];
                }
            }
        }
    } else {
        double* recv = (double*) malloc(2*NX2*NX1*Nmomentum*sizeof(double));
        MPI_Status status;
        MPI_Recv(recv, 2*NX2*NX1*Nmomentum, MPI_DOUBLE, 0, 0, comm, &status);
        int count = 0;
        IDOM_LOOP(i){
            JDOM_LOOP(j){
                for(int l = 0; l < Nmomentum; ++l){
                    x[KBEG][j][i][l] = recv[count];
                    count = count + 1;
                    x[KEND][j][i][l] = recv[count];
                    count = count + 1;
                }
            }
        }
        free(recv);
    }

    IDOM_LOOP(i){
        JDOM_LOOP(j){
            for(int l = 0; l < Nmomentum; ++l){
                for (k = KBEG + 1; k < KEND; ++k) {
                    x[k][j][i][l] = rightPart[k][j][i][l] - a[k][j][i][l] * x[KBEG][j][i][l] - c[k][j][i][l] * x[KEND][j][i][l];
                    if ((x[k][j][i][l] != x[k][j][i][l]) || (0 * x[k][j][i][l] != 0 * x[k][j][i][l])) {
                        printf("x = NaN in solver Z, k = %d , j = %d, i = %d, l = %d, rank %d\n", k, j, i, l, rank);
                        exit(0);
                    }
                }

                //delete[] d;
            }
        }
    }

    FreeArray4D(parallelA);
    FreeArray4D(parallelB);
    FreeArray4D(parallelC);
    FreeArray4D(parallelRightPart);
    FreeArray4D(parallelX);
}

#endif


