#include <stdio.h>
#include <math.h>
#include <stdbool.h>
#ifdef PARALLEL
#include <mpi.h>
#endif
//#include <crtdbg.h>

//#include "memory_debug.h"

#include "macros.h"
#include "pluto.h"
#include "specialmath.h"
#include "matrixElement.h"

double scalarMultiplyLargeVectors(double ****a, double ****b, int lnumber){
    int i,j,k;
    double result[1];
    result[0] = 0;
    double global_result[1];
    global_result[0] = 0;
    DOM_LOOP(k,j,i){
        for(int l = 0; l < lnumber; ++l){
            result[0] += a[k][j][i][l]*b[k][j][i][l];
            if(result[0] != result[0]){
                printf("scalar mult = NaN\n");
                QUIT_PLUTO(1);
            }
        }
    }

#ifdef PARALLEL
    MPI_Allreduce(result, global_result,1,MPI_DOUBLE,MPI_SUM, MPI_COMM_WORLD);
    result[0] = global_result[0];
#endif

    return result[0];
}

double evaluateError(double** hessenbergMatrix, double* vector, double beta, int n) {
    double* resVector = (double*) malloc((n+1)*sizeof(double));

	for (int i = 0; i < n + 1; ++i) {
		resVector[i] = 0;
		for (int j = 0; j < n; ++j) {
			resVector[i] += hessenbergMatrix[i][j] * vector[j];
		}
		if (i == 0) {
			resVector[i] -= beta;
		}
	}

	double norm = 0;
	for (int i = 0; i < n + 1; ++i) {
		norm += resVector[i] * resVector[i];
	}

    free(resVector);

	return sqrt(norm);
}

/*double**** multiplySpecialMatrixVector(MatrixElementNode***** matrix, double**** vector, int lnumber) {
    int i, j, k, l;

    double**** result = ARRAY_4D(NX3_TOT, NX2_TOT, NX1_TOT, lnumber, double);

    DOM_LOOP(k,j,i){
        for(l = 0; l < lnumber; ++l){
            result[k][j][i][l] = 0;
        }
    }


    DOM_LOOP(k,j,i){
        for (int l = 0; l < lnumber; ++l) {
            result[k][j][i][l] = 0;
            struct MatrixElementNode_* node = matrix[k][j][i][l];
            while(node != NULL){
                struct MatrixElement_ element = node->element;

                result[k][j][i][l] += element.value * vector[element.k][element.j][element.i][element.l];
                node = node->next;
            }
        }
    }

	return result;
}*/

void multiplySpecialMatrixVector(double**** result, MatrixElementNode***** matrix, double**** vector, int lnumber, int* par_dim, bool periodicX, bool periodicY, bool periodicZ) {
    int i,j,k;
    TOT_LOOP(k,j,i){
        for(int l = 0; l < lnumber; ++l){
            result[k][j][i][l] = 0;
        }
    }
    DOM_LOOP(k,j,i){
				for (int l = 0; l < lnumber; ++l) {
                    result[k][j][i][l] = 0;
                    struct MatrixElementNode_* node = matrix[k][j][i][l];
                    while(node != NULL){
                        struct MatrixElement_ element = node->element;

                        result[k][j][i][l] += element.value * vector[element.k][element.j][element.i][element.l];
                        node = node->next;
                        if(result[k][j][i][l] != result[k][j][i][l]){
                            printf("result[%d][%d][%d][%d] = NaN\n", k, j, i, l);
                            QUIT_PLUTO(1);
                        }
                    }
                }
	}

    exchangeLargeVector(result, lnumber, par_dim, SZ_stagx, periodicX, periodicY, periodicZ);
}

void arnoldiIterations(MatrixElementNode***** matrix, double** outHessenbergMatrix, int n,
                       LargeVectorBasis* gmresBasis, double** prevHessenbergMatrix, int lnumber, int* par_dim, bool periodicX, bool periodicY, bool periodicZ) {

    int i,j,k;
#ifdef PARALLEL
    MPI_Comm cartComm = MPI_COMM_WORLD;
	MPI_Barrier(cartComm);
#endif

	for (int i = 0; i < n; ++i) {
		if (i < n - 1) {
			for (int j = 0; j < n - 2; ++j) {
				outHessenbergMatrix[i][j] = prevHessenbergMatrix[i][j];
			}
            free(prevHessenbergMatrix[i]);
		} else {
			for (int j = 0; j < n - 2; ++j) {
				outHessenbergMatrix[i][j] = 0;
			}
		}
	}
    free(prevHessenbergMatrix);
	//printf("update hessenberg\n");
	if (n >= gmresBasis->capacity) {
        resize(gmresBasis, 2 * n);
	}
    multiplySpecialMatrixVector(gmresBasis->array[n - 1], matrix, gmresBasis->array[n - 2], lnumber, par_dim, periodicX, periodicY, periodicZ);
	gmresBasis->size += 1;
	//printf("mult special matrix");

	//printf("start exchange\n");
#ifdef PARALLEL
	MPI_Barrier(cartComm);
#endif

    //exchangeLargeVector(gmresBasis->array[n-1], xnumberAdded, ynumberAdded, znumberAdded, lnumber, additionalBinNumber, periodicX, periodicY, periodicZ, cartComm, cartCoord, cartDim, leftOutGmresBuffer, rightOutGmresBuffer, leftInGmresBuffer, rightInGmresBuffer, frontOutGmresBuffer, backOutGmresBuffer, frontInGmresBuffer, backInGmresBuffer, bottomOutGmresBuffer, topOutGmresBuffer, bottomInGmresBuffer, topInGmresBuffer);


	for (int m = 0; m < n - 1; ++m) {
		//double a = scalarMultiplyLargeVectors(resultBasis[m], tempVector, xnumber, ynumber, znumber, lnumber);
        outHessenbergMatrix[m][n - 2] = scalarMultiplyLargeVectors(gmresBasis->array[m], gmresBasis->array[n - 1], lnumber);
		//printf("outHessenbergMatrix[%d][%d] = %g\n", m, n-2, outHessenbergMatrix[m][n - 2]);

		//for (int i = 0; i < xnumber+1; ++i) {

        DOM_LOOP(k,j,i){
					for (int l = 0; l < lnumber; ++l) {
                            gmresBasis->array[n - 1][k][j][i][l] -= outHessenbergMatrix[m][n - 2] * gmresBasis->array[m][k][j][i][l];
					}
        }
	}
#ifdef PARALLEL
    MPI_Barrier(cartComm);
#endif
    exchangeLargeVector(gmresBasis->array[n - 1], lnumber, par_dim, SZ_stagx, periodicX, periodicY, periodicZ);
	//printf("finish orthogonalisation\n");
	outHessenbergMatrix[n - 1][n - 2] = sqrt(
        scalarMultiplyLargeVectors(gmresBasis->array[n - 1], gmresBasis->array[n - 1], lnumber));
	//printf("outHessenbergMatrix[%d][%d] = %g\n", n-1, n-2, outHessenbergMatrix[n - 1][n - 2]);
	if (outHessenbergMatrix[n - 1][n - 2] > 0) {
        DOM_LOOP(k,j,i){
					for (int l = 0; l < lnumber; ++l) {
                        gmresBasis->array[n - 1][k][j][i][l] /= outHessenbergMatrix[n - 1][n - 2];
					}
		}
#ifdef PARALLEL
		MPI_Barrier(cartComm);
#endif
	} else {
		printf("outHessenbergMatrix[n-1][n-2] == 0\n");
        DOM_LOOP(k,j,i){
					for (int l = 0; l < lnumber; ++l) {
                        gmresBasis->array[n - 1][k][j][i][l] = 0;
					}
		}
	}

    exchangeLargeVector(gmresBasis->array[n - 1], lnumber, par_dim, SZ_stagx, periodicX, periodicY, periodicZ);
}

void generalizedMinimalResidualMethod1(Grid* grid, MatrixElementNode***** matrix, double ****rightPart, double ****outvector, LargeVectorBasis* gmresBasis, int lnumber, double precision,
                                      int maxIteration, int verbosity, bool periodicX, bool periodicY, bool periodicZ) {
    double**** initialVector = (double****) malloc(NX3_TOT*sizeof(double***));
    for(int k = 0; k < NX3_TOT; ++k){
        initialVector[k] = (double***) malloc(NX2_TOT*sizeof(double**));
        for(int j = 0; j < NX2_TOT; ++j){
            initialVector[k][j] = (double**) malloc(NX1_TOT*sizeof(double*));
            for(int i = 0; i < NX1_TOT; ++i){
                initialVector[k][j][i] = (double*) malloc(lnumber*sizeof(double));
                for(int l = 0; l < lnumber; ++l){
                    initialVector[k][j][i][l] = 0;
                }
            }
        }
    }

    generalizedMinimalResidualMethod(grid, matrix, rightPart, outvector, initialVector, gmresBasis, lnumber, precision, maxIteration, verbosity, periodicX, periodicY, periodicZ);

    for(int k = 0; k < NX3_TOT; ++k){
        for(int j = 0; j < NX2_TOT; ++j){
            for(int i = 0; i < NX1_TOT; ++i){
                free(initialVector[k][j][i]);
            }
            free(initialVector[k][j]);
        }
        free(initialVector[k]);
    }
    free(initialVector);
}

void generalizedMinimalResidualMethod(Grid* grid, MatrixElementNode***** matrix, double ****rightPart, double ****outvector, double**** initialVector, LargeVectorBasis* gmresBasis, int lnumber, double precision,
                                      int maxIteration, int verbosity, bool periodicX, bool periodicY, bool periodicZ) {
    int rank = 0;
	int nprocs;
    int i,j,k;
#ifdef PARALLEL
    MPI_Comm cartComm = MPI_COMM_WORLD;
	MPI_Comm_size(cartComm, &nprocs);
	MPI_Comm_rank(cartComm, &rank);

	MPI_Barrier(cartComm);
#endif
	if ((rank == 0) && (verbosity > 0)) printf("start GMRES\n");

    clear(gmresBasis);

    int  par_dim[3] = {0, 0, 0};
    DIM_EXPAND(par_dim[0] = grid->nproc[IDIR] > 1;  ,
               par_dim[1] = grid->nproc[JDIR] > 1;  ,
               par_dim[2] = grid->nproc[KDIR] > 1;)

#ifdef PARALLEL
	MPI_Barrier(cartComm);
#endif
    /*exchangeLargeVector(rightPart, xnumberAdded, ynumberAdded, znumberAdded, lnumber, additionalBinNumber, periodicX,
	                    periodicY, periodicZ, cartComm, cartCoord, cartDim, leftOutGmresBuffer, rightOutGmresBuffer,
	                    leftInGmresBuffer, rightInGmresBuffer, frontOutGmresBuffer, backOutGmresBuffer, frontInGmresBuffer,
                        backInGmresBuffer, bottomOutGmresBuffer, topOutGmresBuffer, bottomInGmresBuffer, topInGmresBuffer);*/
    exchangeLargeVector(rightPart, lnumber, par_dim, SZ_stagx, periodicX, periodicY, periodicZ);
    double norm = sqrt(scalarMultiplyLargeVectors(rightPart, rightPart, lnumber));
    //printf("right part norm = %g\n", norm);
	//alertNaNOrInfinity(norm, "right partnorm = NaN in gmres\n");

	if (norm == 0) {
        DOM_LOOP(k,j,i){
					for (int l = 0; l < lnumber; ++l) {
                        outvector[k][j][i][l] = 0;
					}
		}
		return;
	}
    exchangeLargeVector(outvector, lnumber, par_dim, SZ_stagx, periodicX, periodicY, periodicZ);

	//#pragma omp parallel for
    DOM_LOOP(k,j,i){
				for (int l = 0; l < lnumber; ++l) {
                    rightPart[k][j][i][l] /= norm;
				}
	}
    exchangeLargeVector(rightPart, lnumber, par_dim, SZ_stagx, periodicX, periodicY, periodicZ);

    double normRightPart = sqrt(scalarMultiplyLargeVectors(rightPart, rightPart, lnumber));

    int matrixDimension = lnumber * grid->np_int_glob[0] * grid->np_int_glob[1]*grid->np_int_glob[2];

	double** hessenbergMatrix;
	double** newHessenbergMatrix;
    hessenbergMatrix = (double**) malloc(1*sizeof(double*));
    hessenbergMatrix[0] = (double*) malloc(1*sizeof(double));
	hessenbergMatrix[0][0] = 0;

    double** Qmatrix = (double**) malloc(2*sizeof(double*));
    double** Rmatrix = (double**) malloc(2*sizeof(double*));
    double** oldQmatrix = (double**) malloc(2*sizeof(double*));
    double** oldRmatrix = (double**) malloc(2*sizeof(double*));

	for (int i = 0; i < 2; ++i) {
        Qmatrix[i] = (double*) malloc(2*sizeof(double));
        oldQmatrix[i] = (double*) malloc(2*sizeof(double));
	}

    Rmatrix[0] = (double*) malloc(1*sizeof(double));
    Rmatrix[1] = (double*) malloc(1*sizeof(double));
    oldRmatrix[0] = (double*) malloc(1*sizeof(double));
    oldRmatrix[1] = (double*) malloc(1*sizeof(double));

    TOT_LOOP(k,j,i){
        for (int l = 0; l < lnumber; ++l) {
            initialVector[k][j][i][l] = initialVector[k][j][i][l] / norm;
        }
    }

	if (gmresBasis->capacity <= 0) {
        resize(gmresBasis, 5);
	}

    multiplySpecialMatrixVector(gmresBasis->array[0], matrix, initialVector, lnumber, par_dim, periodicX, periodicY, periodicZ);

    DOM_LOOP(k,j,i){
				for (int l = 0; l < lnumber; ++l) {
                    gmresBasis->array[0][k][j][i][l] = rightPart[k][j][i][l] - gmresBasis->array[0][k][j][i][l];
				}
	}
    exchangeLargeVector(gmresBasis->array[0], lnumber, par_dim, SZ_stagx, periodicX, periodicY, periodicZ);
	gmresBasis->size = 1;

    double norm1 = sqrt(scalarMultiplyLargeVectors(gmresBasis->array[0], gmresBasis->array[0], lnumber));

    TOT_LOOP(k, j, i){
        for (int l = 0; l < lnumber; ++l) {
            gmresBasis->array[0][k][j][i][l] = gmresBasis->array[0][k][j][i][l]/norm1;
        }
    }

	int n = 2;
    //double beta = 1.0;
    double beta = norm1;
	double error = beta;
    double* y = (double*) malloc(1*sizeof(double));

	double rho;
	double sigma;
	double cosn;
	double sinn;
	double module;

	double relativeError = 1;
	//double maxRelativeError = maxErrorLevel / (matrixDimension);
	double maxRelativeError = precision / (matrixDimension);
	//double maxRelativeError = precision;

    while ((relativeError > MAX(maxRelativeError, 1E-15) && (n < MIN(maxIteration, matrixDimension + 3)))) {
		if ((rank == 0) && (verbosity > 1)) printf("GMRES iteration %d\n", n);
        newHessenbergMatrix = (double**) malloc(n*sizeof(double*));
		for (int i = 0; i < n; ++i) {
            newHessenbergMatrix[i] = (double*) malloc((n-1)*sizeof(double));
		}
        arnoldiIterations(matrix, newHessenbergMatrix, n, gmresBasis, hessenbergMatrix, lnumber, par_dim, periodicX, periodicY, periodicZ);

		hessenbergMatrix = newHessenbergMatrix;

		if (n == 2) {
			rho = hessenbergMatrix[0][0];
			sigma = hessenbergMatrix[1][0];

			module = sqrt(rho * rho + sigma * sigma);

			cosn = rho / module;
			sinn = sigma / module;

			Qmatrix[0][0] = cosn;
			Qmatrix[0][1] = sinn;
			Qmatrix[1][0] = -sinn;
			Qmatrix[1][1] = cosn;

			oldQmatrix[0][0] = Qmatrix[0][0];
			oldQmatrix[0][1] = Qmatrix[0][1];
			oldQmatrix[1][0] = Qmatrix[1][0];
			oldQmatrix[1][1] = Qmatrix[1][1];

			Rmatrix[0][0] = module;
			Rmatrix[1][0] = 0;

			oldRmatrix[0][0] = Rmatrix[0][0];
			oldRmatrix[1][0] = Rmatrix[1][0];

		} else {
            Rmatrix = (double**) malloc(n*sizeof(double*));
			for (int i = 0; i < n; ++i) {
                Rmatrix[i] = (double*) malloc((n-1)*sizeof(double));
				if (i < n - 1) {
					for (int j = 0; j < n - 2; ++j) {
						Rmatrix[i][j] = oldRmatrix[i][j];
					}
				} else {
					for (int j = 0; j < n - 2; ++j) {
						Rmatrix[i][j] = 0;
					}
				}
			}

            Qmatrix = (double**) malloc(n*sizeof(double*));
			for (int i = 0; i < n; ++i) {
                Qmatrix[i] = (double*) malloc(n*sizeof(double));
				if (i < n - 1) {
					for (int j = 0; j < n - 1; ++j) {
						Qmatrix[i][j] = oldQmatrix[i][j];
					}
					Qmatrix[i][n - 1] = 0;
				} else {
					for (int j = 0; j < n - 1; ++j) {
						Qmatrix[i][j] = 0;
					}
					Qmatrix[n - 1][n - 1] = 1;
				}
			}

			for (int i = 0; i < n; ++i) {
				Rmatrix[i][n - 2] = 0;
				for (int j = 0; j < n; ++j) {
					Rmatrix[i][n - 2] += Qmatrix[i][j] * hessenbergMatrix[j][n - 2];
				}
			}
			rho = Rmatrix[n - 2][n - 2];
			sigma = Rmatrix[n - 1][n - 2];

			module = sqrt(rho * rho + sigma * sigma);

			cosn = rho / module;
			sinn = sigma / module;

			Rmatrix[n - 2][n - 2] = module;
			Rmatrix[n - 1][n - 2] = 0;

			for (int j = 0; j < n - 1; ++j) {
				Qmatrix[n - 2][j] = cosn * oldQmatrix[n - 2][j];
				Qmatrix[n - 1][j] = -sinn * oldQmatrix[n - 2][j];
			}
			Qmatrix[n - 2][n - 1] = sinn;
			Qmatrix[n - 1][n - 1] = cosn;
		}

        free(y);
        y = (double*) malloc((n-1)*sizeof(double));
		//printf("n = %d\n", n);

		for (int i = n - 2; i >= 0; --i) {
			y[i] = beta * Qmatrix[i][0];
			//printf("y[%d] = %g\n", i, y[i]);
			for (int j = n - 2; j > i; --j) {
				y[i] -= Rmatrix[i][j] * y[j];
			}
			if (Rmatrix[i][i] > 0) {
				y[i] /= Rmatrix[i][i];
			} else {
				y[i] = 0;
				printf("Rmatrix[%d][%d] = 0\n", i, i);
			}
			//printf("y[%d] = %g\n", i, y[i]);
			//alertNaNOrInfinity(y[i], "y = NaN\n");
		}

		error = fabs(beta * Qmatrix[n - 1][0]);
		//double error1 = evaluateError(hessenbergMatrix, y, beta, n);

		relativeError = error / normRightPart;

		for (int i = 0; i < n - 1; ++i) {
            free(oldQmatrix[i]);
            free(oldRmatrix[i]);
		}
		if (n == 2) {
            free(oldQmatrix[1]);
            free(oldRmatrix[1]);
		}
        free(oldQmatrix);
        free(oldRmatrix);

		oldQmatrix = Qmatrix;
		oldRmatrix = Rmatrix;

		n++;
	}

	n = n - 1;
	if((rank == 0) && (verbosity > 0)) printf("total GMRES iteration = %d\n", n);
    if (rank == 0) printf("relative GMRES error = %g\n", relativeError);
    if (rank == 0) printf("relative GMRES error*Nx*Ny*Nz*Np = %g\n", relativeError*matrixDimension);

	//out result


    DOM_LOOP(k,j,i){
				for (int l = 0; l < lnumber; ++l) {
                    outvector[k][j][i][l] = initialVector[k][j][i][l]*norm;
					for (int m = 0; m < n - 1; ++m) {
                        outvector[k][j][i][l] += gmresBasis->array[m][k][j][i][l] * y[m] * norm;
						//outvector[i][l] += basis[m][i][l] * y[m];
					}
				}
	}

    clear(gmresBasis);
#ifdef PARALLEL
	MPI_Barrier(cartComm);
#endif
    exchangeLargeVector(outvector, lnumber, par_dim, SZ_stagx, periodicX, periodicY, periodicZ);
    //exchangeLargeVector(outvector, xnumberAdded, ynumberAdded, znumberAdded, lnumber, additionalBinNumber, periodicX, periodicY, periodicZ, cartComm, cartCoord, cartDim, leftOutGmresBuffer, rightOutGmresBuffer, leftInGmresBuffer, rightInGmresBuffer, frontOutGmresBuffer, backOutGmresBuffer, frontInGmresBuffer, backInGmresBuffer, bottomOutGmresBuffer, topOutGmresBuffer, bottomInGmresBuffer, topInGmresBuffer);
#ifdef PARALLEL
    MPI_Barrier(cartComm);
#endif
	

	for (int i = 0; i < n; ++i) {
        free(Qmatrix[i]);
        free(Rmatrix[i]);
        free(hessenbergMatrix[i]);
	}
    free(Qmatrix);
    free(Rmatrix);
    free(hessenbergMatrix);
    free(y);

    DOM_LOOP(k,j,i){
				for (int l = 0; l < lnumber; ++l) {
                    rightPart[k][j][i][l] *= norm;
					/*for (int m = 0; m < matrix[i][j][k][l].size(); ++m) {
						double value = matrix[i][j][k][l][m].value;
						//matrix[i][l][m].value *= norm;
						value = matrix[i][j][k][l][m].value;
					}*/
				}
	}
    exchangeLargeVector(rightPart, lnumber, par_dim, SZ_stagx, periodicX, periodicY, periodicZ);
}


/*void conjugateGradientMethod(std::vector < MatrixElement >**** matrix, double**** rightPart, double**** outVector,
                             int xnumberAdded, int ynumberAdded, int znumberAdded, int additionalBinNumber, int lnumber,
                             double precision, int maxIteration, bool periodicX, bool
                             periodicY, bool periodicZ, int verbosity, MPI_Comm& cartComm, int* cartCoord,
                             int* cartDim) {
	int rank;
	int nprocs;
	MPI_Comm_size(cartComm, &nprocs);
	MPI_Comm_rank(cartComm, &rank);
	if (rank == 0 && verbosity > 0) printf("start conjugate gradient\n");

	double**** residual = new double***[xnumberAdded];
	double**** prevResidual = new double***[xnumberAdded];
	double**** z = new double***[xnumberAdded];
	double**** tempVector = new double***[xnumberAdded];

	for (int i = 0; i < xnumberAdded; ++i) {
		residual[i] = new double**[ynumberAdded];
		prevResidual[i] = new double**[ynumberAdded];
		z[i] = new double**[ynumberAdded];
		tempVector[i] = new double**[ynumberAdded];
		for (int j = 0; j < ynumberAdded; ++j) {
			residual[i][j] = new double*[znumberAdded];
			prevResidual[i][j] = new double*[znumberAdded];
			z[i][j] = new double*[znumberAdded];
			tempVector[i][j] = new double*[znumberAdded];
			for (int k = 0; k < znumberAdded; ++k) {
				residual[i][j][k] = new double[lnumber];
				prevResidual[i][j][k] = new double[lnumber];
				z[i][j][k] = new double[lnumber];
				tempVector[i][j][k] = new double[lnumber];
				for (int l = 0; l < lnumber; ++l) {
					outVector[i][j][k][l] = 0;
					prevResidual[i][j][k][l] = rightPart[i][j][k][l];
					z[i][j][k][l] = rightPart[i][j][k][l];
					tempVector[i][j][k][l] = 0;
				}
			}
		}
	}


	int iteration = 0;


	double prevResidualNorm2 = scalarMultiplyLargeVectors(prevResidual, prevResidual, xnumberAdded, ynumberAdded,
	                                                      znumberAdded,
	                                                      additionalBinNumber, lnumber, periodicX, periodicY, periodicZ,
	                                                      rank, nprocs, cartComm, cartCoord, cartDim);
	double residualNorm2 = prevResidualNorm2;
	double rightPartNorm2 = scalarMultiplyLargeVectors(rightPart, rightPart, xnumberAdded, ynumberAdded, znumberAdded,
	                                                   additionalBinNumber, lnumber, periodicX, periodicY, periodicZ, rank,
	                                                   nprocs, cartComm, cartCoord, cartDim);

	double relativeError = sqrt(residualNorm2 / rightPartNorm2);

	while ((iteration < maxIteration) && (iteration < xnumberAdded * ynumberAdded * znumberAdded * lnumber) && (
		relativeError > (precision / (xnumberAdded * ynumberAdded * znumberAdded * lnumber)))) {
		if (rank == 0 && verbosity > 1) printf("conjugate gradient iteration %d\n", iteration);

		multiplySpecialMatrixVector(tempVector, matrix, z, xnumberAdded, ynumberAdded, znumberAdded, additionalBinNumber,
		                            lnumber, periodicX, periodicY, periodicZ, rank, nprocs, cartComm, cartCoord, cartDim);

		double alpha = prevResidualNorm2 / scalarMultiplyLargeVectors(tempVector, z, xnumberAdded, ynumberAdded, znumberAdded,
		                                                              additionalBinNumber,
		                                                              lnumber, periodicX, periodicY, periodicZ, rank, nprocs,
		                                                              cartComm, cartCoord, cartDim);

		for (int i = 0; i < xnumberAdded + 1; ++i) {
			for (int j = 0; j < ynumberAdded; ++j) {
				for (int k = 0; k < znumberAdded; ++k) {
					for (int l = 0; l < lnumber; ++l) {
						outVector[i][j][k][l] += alpha * z[i][j][k][l];
						residual[i][j][k][l] = prevResidual[i][j][k][l] - alpha * tempVector[i][j][k][l];
					}
				}
			}
		}

		residualNorm2 = scalarMultiplyLargeVectors(residual, residual, xnumberAdded, ynumberAdded, znumberAdded,
		                                           additionalBinNumber, lnumber, periodicX, periodicY, periodicZ, rank,
		                                           nprocs, cartComm, cartCoord, cartDim);

		double beta = residualNorm2 / prevResidualNorm2;

		for (int i = 0; i < xnumberAdded; ++i) {
			for (int j = 0; j < ynumberAdded; ++j) {
				for (int k = 0; k < znumberAdded; ++k) {
					for (int l = 0; l < lnumber; ++l) {
						z[i][j][k][l] = residual[i][j][k][l] + beta * z[i][j][k][l];
					}
				}
			}
		}

		prevResidualNorm2 = residualNorm2;

		relativeError = sqrt(residualNorm2 / rightPartNorm2);
		iteration++;
	}

	for (int i = 0; i < xnumberAdded; ++i) {
		for (int j = 0; j < ynumberAdded; ++j) {
			for (int k = 0; k < znumberAdded; ++k) {
				delete[] residual[i][j][k];
				delete[] prevResidual[i][j][k];
				delete[] z[i][j][k];
				delete[] tempVector[i][j][k];
			}
			delete[] residual[i][j];
			delete[] prevResidual[i][j];
			delete[] z[i][j];
			delete[] tempVector[i][j];
		}
		delete[] residual[i];
		delete[] prevResidual[i];
		delete[] z[i];
		delete[] tempVector[i];
	}
	delete[] residual;
	delete[] prevResidual;
	delete[] z;
	delete[] tempVector;
}*/

void biconjugateGradientMethod(Grid* grid, MatrixElementNode***** matrix, double ****rightPart, double ****outvector, int lnumber, double precision,
                               int maxIteration, int verbosity, bool periodicX, bool periodicY, bool periodicZ) {
    int rank = 0;
	int nprocs;
    int i, j, k;
#ifdef PARALLEL
    MPI_Comm cartComm = MPI_COMM_WORLD;
    MPI_Comm_size(cartComm, &nprocs);
    MPI_Comm_rank(cartComm, &rank);

    MPI_Barrier(cartComm);
#endif
	if (rank == 0 && verbosity > 0) printf("start biconjugate gradient\n");

    double rightPartNorm2 = scalarMultiplyLargeVectors(rightPart, rightPart, lnumber);

    if(rightPartNorm2 <= 0){
        DOM_LOOP(k,j,i){
                    for (int l = 0; l < lnumber; ++l) {
                        outvector[k][j][i][l] = 0;
                    }
        }
        return;
    }

    int  par_dim[3] = {0, 0, 0};
    DIM_EXPAND(par_dim[0] = grid->nproc[IDIR] > 1;  ,
               par_dim[1] = grid->nproc[JDIR] > 1;  ,
               par_dim[2] = grid->nproc[KDIR] > 1;)

    int matrixDimension = lnumber * grid->np_int_glob[0] * grid->np_int_glob[1]*grid->np_int_glob[2];

    double**** residual = (double****) malloc(NX3_TOT*sizeof(double***));
    double**** prevResidual = (double****) malloc(NX3_TOT*sizeof(double***));
    double**** z = (double****) malloc(NX3_TOT*sizeof(double***));
    double**** p = (double****) malloc(NX3_TOT*sizeof(double***));
    double**** s = (double****) malloc(NX3_TOT*sizeof(double***));
    double**** tempVector = (double****) malloc(NX3_TOT*sizeof(double***));
    double**** tempVector2 = (double****) malloc(NX3_TOT*sizeof(double***));
    MatrixElementNode***** transposedMatrix = (MatrixElementNode*****) malloc(NX3_TOT*sizeof(MatrixElementNode****));

    for (int i = 0; i < NX3_TOT; ++i) {
        residual[i] = (double***) malloc(NX2_TOT*sizeof(double**));
        prevResidual[i] = (double***) malloc(NX2_TOT*sizeof(double**));
        z[i] = (double***) malloc(NX2_TOT*sizeof(double**));
        p[i] = (double***) malloc(NX2_TOT*sizeof(double**));
        s[i] = (double***) malloc(NX2_TOT*sizeof(double**));
        tempVector[i] = (double***) malloc(NX2_TOT*sizeof(double**));
        tempVector2[i] = (double***) malloc(NX2_TOT*sizeof(double**));
        transposedMatrix[i] = (MatrixElementNode****) malloc(NX2_TOT*sizeof(MatrixElementNode***));
        for (int j = 0; j < NX2_TOT; ++j) {
            residual[i][j] = (double**) malloc(NX1_TOT*sizeof(double*));
            prevResidual[i][j] = (double**) malloc(NX1_TOT*sizeof(double*));
            z[i][j] = (double**) malloc(NX1_TOT*sizeof(double*));
            p[i][j] = (double**) malloc(NX1_TOT*sizeof(double*));
            s[i][j] = (double**) malloc(NX1_TOT*sizeof(double*));
            tempVector[i][j] = (double**) malloc(NX1_TOT*sizeof(double*));
            tempVector2[i][j] = (double**) malloc(NX1_TOT*sizeof(double*));
            transposedMatrix[i][j] = (MatrixElementNode***) malloc(NX1_TOT*sizeof(MatrixElementNode**));
            for (int k = 0; k < NX1_TOT; ++k) {
                residual[i][j][k] = (double*) malloc(lnumber*sizeof(double));
                prevResidual[i][j][k] = (double*) malloc(lnumber*sizeof(double));
                z[i][j][k] = (double*) malloc(lnumber*sizeof(double));
                p[i][j][k] = (double*) malloc(lnumber*sizeof(double));
                s[i][j][k] = (double*) malloc(lnumber*sizeof(double));
                tempVector[i][j][k] = (double*) malloc(lnumber*sizeof(double));
                tempVector2[i][j][k] = (double*) malloc(lnumber*sizeof(double));
                transposedMatrix[i][j][k] = (MatrixElementNode**) malloc(lnumber*sizeof(MatrixElementNode*));
				for (int l = 0; l < lnumber; ++l) {
                    outvector[i][j][k][l] = 0;
					prevResidual[i][j][k][l] = rightPart[i][j][k][l];
					z[i][j][k][l] = rightPart[i][j][k][l];
					p[i][j][k][l] = rightPart[i][j][k][l];
					s[i][j][k][l] = rightPart[i][j][k][l];
					tempVector[i][j][k][l] = 0;
					tempVector2[i][j][k][l] = 0;
                    transposedMatrix[i][j][k][l] = NULL;
				}
			}
		}
	}

    transposeSpecialMatrix(transposedMatrix, matrix, lnumber);


	int iteration = 0;


    double prevResidualNorm2 = scalarMultiplyLargeVectors(p, prevResidual, lnumber);
	double residualNorm2 = prevResidualNorm2;

	double relativeError = sqrt(residualNorm2 / rightPartNorm2);

    while ((iteration < maxIteration) && (iteration < matrixDimension) && (
        relativeError > (precision / matrixDimension))) {
		if (rank == 0 && verbosity > 1) printf("biconjugate gradient iteration %d\n", iteration);

        multiplySpecialMatrixVector(tempVector, matrix, z, lnumber, par_dim, periodicX, periodicY, periodicZ);
        multiplySpecialMatrixVector(tempVector2, transposedMatrix, s, lnumber, par_dim, periodicX, periodicY, periodicZ);

        double alpha = prevResidualNorm2 / scalarMultiplyLargeVectors(tempVector, s, lnumber);

        for (int i = 0; i < NX3_TOT; ++i) {
            for (int j = 0; j < NX2_TOT; ++j) {
                for (int k = 0; k < NX1_TOT; ++k) {
					for (int l = 0; l < lnumber; ++l) {
                        outvector[i][j][k][l] += alpha * z[i][j][k][l];
						residual[i][j][k][l] = prevResidual[i][j][k][l] - alpha * tempVector[i][j][k][l];
						p[i][j][k][l] = p[i][j][k][l] - alpha * tempVector2[i][j][k][l];
					}
				}
			}
		}

        exchangeLargeVector(outvector, lnumber, par_dim, SZ_stagx, periodicX, periodicY, periodicZ);
        exchangeLargeVector(residual, lnumber, par_dim, SZ_stagx, periodicX, periodicY, periodicZ);
        exchangeLargeVector(p, lnumber, par_dim, SZ_stagx, periodicX, periodicY, periodicZ);

        residualNorm2 = scalarMultiplyLargeVectors(p, residual, lnumber);

		double beta = residualNorm2 / prevResidualNorm2;

        for (int i = 0; i < NX3_TOT; ++i) {
            for (int j = 0; j < NX2_TOT; ++j) {
                for (int k = 0; k < NX1_TOT; ++k) {
					for (int l = 0; l < lnumber; ++l) {
						z[i][j][k][l] = residual[i][j][k][l] + beta * z[i][j][k][l];
						s[i][j][k][l] = p[i][j][k][l] + beta * s[i][j][k][l];
					}
				}
			}
		}

        exchangeLargeVector(z, lnumber, par_dim, SZ_stagx, periodicX, periodicY, periodicZ);
        exchangeLargeVector(s, lnumber, par_dim, SZ_stagx, periodicX, periodicY, periodicZ);

		prevResidualNorm2 = residualNorm2;
        //todo???
		relativeError = sqrt(residualNorm2 / rightPartNorm2);
		iteration++;
	}

    if((rank == 0) && (verbosity > 0)) printf("total GMRES iteration = %d\n", iteration);
    if (rank == 0) printf("relative biconjugate gradient error = %g\n", relativeError);
    if (rank == 0) printf("relative biconjugate gradient error*Nx*Ny*Nz*Np = %g\n", relativeError*matrixDimension);

    for (int i = 0; i < NX3_TOT; ++i) {
        for (int j = 0; j < NX2_TOT; ++j) {
            for (int k = 0; k < NX1_TOT; ++k) {
                for(int l = 0; l < lnumber; ++l){
                    MatrixElementNode* curNode = transposedMatrix[i][j][k][l];
                    while(curNode != NULL){
                        MatrixElementNode* tempNode = curNode;
                        curNode = curNode->next;
                        free(tempNode);
                    }
                }
                free(residual[i][j][k]);
                free(prevResidual[i][j][k]);
                free(z[i][j][k]);
                free(p[i][j][k]);
                free(s[i][j][k]);
                free(tempVector[i][j][k]);
                free(tempVector2[i][j][k]);
                free(transposedMatrix[i][j][k]);
			}
            free(residual[i][j]);
            free(prevResidual[i][j]);
            free(z[i][j]);
            free(p[i][j]);
            free(s[i][j]);
            free(tempVector[i][j]);
            free(tempVector2[i][j]);
            free(transposedMatrix[i][j]);
		}
        free(residual[i]);
        free(prevResidual[i]);
        free(z[i]);
        free(p[i]);
        free(s[i]);
        free(tempVector[i]);
        free(tempVector2[i]);
        free(transposedMatrix[i]);
	}
    free(residual);
    free(prevResidual);
    free(z);
    free(p);
    free(s);
    free(tempVector);
    free(tempVector2);
    free(transposedMatrix);
}

void transposeSpecialMatrix(MatrixElementNode***** result, MatrixElementNode***** matrix, int lnumber){
    for(int k = 0; k < NX3_TOT; ++k){
        for(int j = 0; j < NX2_TOT; ++j){
            for(int i = 0; i < NX1_TOT; ++i){
                for(int l = 0; l < lnumber; ++l){
                    MatrixElementNode* curNode = result[k][j][i][l];
                    while(curNode != NULL){
                        MatrixElementNode* tempNode = curNode;
                        curNode = curNode->next;
                        free(tempNode);
                    }
                }
            }
        }
    }

    for(int k = 0; k < NX3_TOT; ++k){
        for(int j = 0; j < NX2_TOT; ++j){
            for(int i = 0; i < NX1_TOT; ++i){
                for(int l = 0; l < lnumber; ++l){
                    MatrixElementNode* curNode = matrix[k][j][i][l];
                    while(curNode != NULL){
                        double value = curNode->element.value;
                        int ei = curNode->element.i;
                        int ej = curNode->element.j;
                        int ek = curNode->element.k;
                        int el = curNode->element.l;

                        MatrixElementNode* tempNode = result[ek][ej][ei][el];
                        result[ek][ej][ei][el] = (MatrixElementNode*) malloc(sizeof(MatrixElementNode));
                        result[ek][ej][ei][el]->element = createMatrixElement(value, k, j, i, l);
                        result[ek][ej][ei][el]->prev = NULL;
                        result[ek][ej][ei][el]->next = tempNode;

                        curNode = curNode->next;
                    }
                }
            }
        }
    }
}


void biconjugateStabilizedGradientMethod(Grid* grid, MatrixElementNode***** matrix, double**** rightPart,
                                         double**** outVector, int lnumber, double precision, int maxIteration, int verbosity, bool periodicX, bool periodicY, bool periodicZ) {
    int rank = 0;
    int nprocs;
    int i, j, k;
#ifdef PARALLEL
    MPI_Comm cartComm = MPI_COMM_WORLD;
    MPI_Comm_size(cartComm, &nprocs);
    MPI_Comm_rank(cartComm, &rank);

    MPI_Barrier(cartComm);
#endif
    if (rank == 0 && verbosity > 0) printf("start biconjugate stabilized gradient\n");

    double rightPartNorm2 = scalarMultiplyLargeVectors(rightPart, rightPart, lnumber);

    if(rightPartNorm2 <= 0){
        DOM_LOOP(k,j,i){
                    for (int l = 0; l < lnumber; ++l) {
                        outVector[k][j][i][l] = 0;
                    }
        }
        return;
    }

    int  par_dim[3] = {0, 0, 0};
    DIM_EXPAND(par_dim[0] = grid->nproc[IDIR] > 1;  ,
               par_dim[1] = grid->nproc[JDIR] > 1;  ,
               par_dim[2] = grid->nproc[KDIR] > 1;)

    bool converges = false;

    int matrixDimension = lnumber * grid->np_int_glob[0] * grid->np_int_glob[1]*grid->np_int_glob[2];

	double alpha = 1;
	double rho = 1;
	double omega = 1;
    exchangeLargeVector(rightPart, lnumber, par_dim, SZ_stagx, periodicX, periodicY, periodicZ);
    rightPartNorm2 = scalarMultiplyLargeVectors(rightPart, rightPart, lnumber);

    double**** residual = (double****) malloc(NX3_TOT*sizeof(double***));
    double**** firstResidual = (double****) malloc(NX3_TOT*sizeof(double***));
    double**** v = (double****) malloc(NX3_TOT*sizeof(double***));
    double**** p = (double****) malloc(NX3_TOT*sizeof(double***));
    double**** s = (double****) malloc(NX3_TOT*sizeof(double***));
    double**** t = (double****) malloc(NX3_TOT*sizeof(double***));

    for (int i = 0; i < NX3_TOT; ++i) {
        residual[i] = (double***) malloc(NX2_TOT*sizeof(double**));
        firstResidual[i] = (double***) malloc(NX2_TOT*sizeof(double**));
        v[i] = (double***) malloc(NX2_TOT*sizeof(double**));
        p[i] = (double***) malloc(NX2_TOT*sizeof(double**));
        s[i] = (double***) malloc(NX2_TOT*sizeof(double**));
        t[i] = (double***) malloc(NX2_TOT*sizeof(double**));
        for (int j = 0; j < NX2_TOT; ++j) {
            residual[i][j] = (double**) malloc(NX1_TOT*sizeof(double*));
            firstResidual[i][j] = (double**) malloc(NX1_TOT*sizeof(double*));
            v[i][j] = (double**) malloc(NX1_TOT*sizeof(double*));
            p[i][j] = (double**) malloc(NX1_TOT*sizeof(double*));
            s[i][j] = (double**) malloc(NX1_TOT*sizeof(double*));
            t[i][j] = (double**) malloc(NX1_TOT*sizeof(double*));
            for (int k = 0; k < NX1_TOT; ++k) {
                residual[i][j][k] = (double*) malloc(lnumber*sizeof(double));
                firstResidual[i][j][k] = (double*) malloc(lnumber*sizeof(double));
                v[i][j][k] = (double*) malloc(lnumber*sizeof(double));
                p[i][j][k] = (double*) malloc(lnumber*sizeof(double));
                s[i][j][k] = (double*) malloc(lnumber*sizeof(double));
                t[i][j][k] = (double*) malloc(lnumber*sizeof(double));
                for (int l = 0; l < lnumber; ++l) {
                    outVector[i][j][k][l] = 0;
                    firstResidual[i][j][k][l] = rightPart[i][j][k][l];
                    v[i][j][k][l] = rightPart[i][j][k][l];
                    p[i][j][k][l] = rightPart[i][j][k][l];
                    s[i][j][k][l] = rightPart[i][j][k][l];
                    t[i][j][k][l] = 0;
                }
            }
        }
    }

    multiplySpecialMatrixVector(t, matrix, outVector, lnumber, par_dim, periodicX, periodicY, periodicZ);


    TOT_LOOP(i,j,k){
				for (int l = 0; l < lnumber; ++l) {
					outVector[i][j][k][l] = 0;
					firstResidual[i][j][k][l] = (rightPart[i][j][k][l]) - t[i][j][k][l];
					residual[i][j][k][l] = (rightPart[i][j][k][l]) - t[i][j][k][l];
					v[i][j][k][l] = 0;
					p[i][j][k][l] = 0;
					s[i][j][k][l] = 0;
					t[i][j][k][l] = 0;
				}
	}


	int iteration = 0;


    double prevResidualNorm2 = scalarMultiplyLargeVectors(residual, residual, lnumber);
	double residualNorm2 = prevResidualNorm2;

	double relativeError = sqrt(residualNorm2 / prevResidualNorm2);

    /*if (fabs(rightPartNorm2) < 1E-300) {
		return;
    }*/

    while ((iteration < maxIteration) && (iteration < matrixDimension) && (relativeError > (precision / matrixDimension))) {
        if ((rank == 0) && (verbosity > 1)) printf("biconjugate stabilized gradient iteration %d\n", iteration);

        double newRho = scalarMultiplyLargeVectors(firstResidual, residual, lnumber);
		if (fabs(rho) < 1E-100) {
			if (rank == 0) printf("rho = 0 in biconjugate\n");
			return;
		}
		if (fabs(omega) < 1E-100) {
			if (rank == 0) printf("omega = 0 in biconjugate\n");
			return;
		}
		double beta = (newRho / rho) * (alpha / omega);
		rho = newRho;

        TOT_LOOP(i,j,k){
					for (int l = 0; l < lnumber; ++l) {
						s[i][j][k][l] = 0;
					}
		}

		//for (int i = 0; i < xnumberAdded; ++i) {
		//for (int j = 0; j < ynumberAdded; ++j) {
		//for (int k = 0; k < znumberAdded; ++k) {
        DOM_LOOP(i,j,k){
					for (int l = 0; l < lnumber; ++l) {
						p[i][j][k][l] = residual[i][j][k][l] + beta * (p[i][j][k][l] - omega * v[i][j][k][l]);
					}
		}
        exchangeLargeVector(p, lnumber, par_dim, SZ_stagx, periodicX, periodicY, periodicZ);
        multiplySpecialMatrixVector(v, matrix, p, lnumber, par_dim, periodicX, periodicY, periodicZ);
        exchangeLargeVector(v, lnumber, par_dim , SZ_stagx, periodicX, periodicY, periodicZ);
        double firstRscalarV = scalarMultiplyLargeVectors(firstResidual, v, lnumber);
		if (fabs(firstRscalarV) < 1E-100) {
			if (rank == 0) printf("firstRscalarV = 0 in biconjugate\n");
			return;
		}
		alpha = rho / firstRscalarV;

		//for (int i = 0; i < xnumberAdded; ++i) {
		//for (int j = 0; j < ynumberAdded; ++j) {
		//for (int k = 0; k < znumberAdded; ++k) {
        DOM_LOOP(i , j, k){
					for (int l = 0; l < lnumber; ++l) {
						s[i][j][k][l] = residual[i][j][k][l] - alpha * v[i][j][k][l];
					}
		}

        exchangeLargeVector(s, lnumber, par_dim, SZ_stagx, periodicX, periodicY, periodicZ);
        multiplySpecialMatrixVector(t, matrix, s, lnumber, par_dim, periodicX, periodicY, periodicZ);
        exchangeLargeVector(t, lnumber, par_dim, SZ_stagx, periodicX, periodicY, periodicZ);
        double tnorm2 = scalarMultiplyLargeVectors(t, t, lnumber);
		if (fabs(tnorm2) < 1E-100) {
			if (rank == 0) printf("tnorm2 = 0 in biconjugate\n");
			return;
		}
        omega = scalarMultiplyLargeVectors(t, s, lnumber) / tnorm2;

		//for (int i = 0; i < xnumberAdded; ++i) {
		//for (int j = 0; j < ynumberAdded; ++j) {
		//for (int k = 0; k < znumberAdded; ++k) {
        DOM_LOOP(i,j,k){
					for (int l = 0; l < lnumber; ++l) {
						outVector[i][j][k][l] = outVector[i][j][k][l] + omega * s[i][j][k][l] + alpha * p[i][j][k][l];
						residual[i][j][k][l] = s[i][j][k][l] - omega * t[i][j][k][l];
					}
		}
        exchangeLargeVector(outVector, lnumber, par_dim, SZ_stagx, periodicX, periodicY, periodicZ);
        exchangeLargeVector(residual, lnumber, par_dim, SZ_stagx, periodicX, periodicY, periodicZ);

        residualNorm2 = scalarMultiplyLargeVectors(residual, residual, lnumber);

		prevResidualNorm2 = residualNorm2;

		relativeError = sqrt(residualNorm2 / rightPartNorm2);
		iteration++;
	}

    if((rank == 0) && (verbosity > 0)) printf("total biconjugate stabilized gradientiteration = %d\n", iteration);
    if (rank == 0) printf("relative biconjugate stabilized gradient error = %g\n", relativeError);
    if (rank == 0) printf("relative biconjugate stabilized gradient error*Nx*Ny*Nz*Np = %g\n", relativeError*matrixDimension);

	converges = true;

	if (relativeError > 0.1) {
		converges = false;
	}

    /*for (int i = 0; i < xnumberAdded; ++i) {
		for (int j = 0; j < ynumberAdded; ++j) {
			for (int k = 0; k < znumberAdded; ++k) {
				for (int l = 0; l < lnumber; ++l) {
					//outVector[i][j][k][l] *= tempNorm;
				}
			}
		}
    }*/

    for (int i = 0; i < NX3_TOT; ++i) {
        for (int j = 0; j < NX2_TOT; ++j) {
            for (int k = 0; k < NX1_TOT; ++k) {
                free(residual[i][j][k]);
                free(firstResidual[i][j][k]);
                free(v[i][j][k]);
                free(p[i][j][k]);
                free(s[i][j][k]);
                free(t[i][j][k]);
            }
            free(residual[i][j]);
            free(firstResidual[i][j]);
            free(v[i][j]);
            free(p[i][j]);
            free(s[i][j]);
            free(t[i][j]);
        }
        free(residual[i]);
        free(firstResidual[i]);
        free(v[i]);
        free(p[i]);
        free(s[i]);
        free(t[i]);
    }
    free(residual);
    free(firstResidual);
    free(v);
    free(p);
    free(s);
    free(t);
}

/*void gaussSeidelMethod(std::vector < MatrixElement >**** matrix, double**** rightPart, double**** outVector,
                       int xnumberAdded,
                       int ynumberAdded, int znumberAdded, int additionalBinNumber, int lnumber, int xnumberGeneral,
                       int znumberGeneral, int ynumberGeneral, double precision, int maxIteration, bool periodicX,
                       bool periodicY, bool
                       periodicZ, bool startFromZero, int verbocity, MPI_Comm& cartComm, int* cartCoord, int* cartDim) {
	int rank;
	int nprocs;
	MPI_Comm_size(cartComm, &nprocs);
	MPI_Comm_rank(cartComm, &rank);
	double* rightOutBuffer = new double[ynumberAdded * znumberAdded * lnumber];
	double* rightInBuffer = new double[ynumberAdded * znumberAdded * lnumber];
	double* leftOutBuffer = new double[ynumberAdded * znumberAdded * lnumber];
	double* leftInBuffer = new double[ynumberAdded * znumberAdded * lnumber];

	double* backOutBuffer = new double[xnumberAdded * znumberAdded * lnumber];
	double* backInBuffer = new double[xnumberAdded * znumberAdded * lnumber];
	double* frontOutBuffer = new double[xnumberAdded * znumberAdded * lnumber];
	double* frontInBuffer = new double[xnumberAdded * znumberAdded * lnumber];

	double* topOutBuffer = new double[ynumberAdded * xnumberAdded * lnumber];
	double* topInBuffer = new double[ynumberAdded * xnumberAdded * lnumber];
	double* bottomOutBuffer = new double[ynumberAdded * xnumberAdded * lnumber];
	double* bottomInBuffer = new double[ynumberAdded * xnumberAdded * lnumber];

	if ((rank == 0) && (verbocity > 0)) printf("start gauss-seidel\n");
	double normRightPart = scalarMultiplyLargeVectors(rightPart, rightPart, xnumberAdded, ynumberAdded, znumberAdded,
	                                                  additionalBinNumber, lnumber, periodicX, periodicY, periodicZ, rank,
	                                                  nprocs, cartComm, cartCoord, cartDim) / (xnumberGeneral *
		ynumberGeneral * znumberGeneral * lnumber);
	for (int i = 0; i < xnumberAdded + 1; ++i) {
		for (int j = 0; j < ynumberAdded; ++j) {
			for (int k = 0; k < znumberAdded; ++k) {
				for (int l = 0; l < lnumber; ++l) {
					//outVector[i][j][k][l] = uniformDistribution()*normRightPart/matrix[0][0][0][0][0].vfalue;
					if (startFromZero) {
						outVector[i][j][k][l] = 0;
					}
				}
			}
		}
	}


	int curIteration = 0;
	while (curIteration < maxIteration) {
		if ((rank == 0) && (verbocity > 1)) printf("Gauss-Seidel iteration %d\n", curIteration);
		if (rank > 0)
			receiveLargeVectorFromLeft(rightPart, leftInBuffer, xnumberAdded, ynumberAdded, znumberAdded, lnumber,
			                           additionalBinNumber, cartComm);
		for (int i = 0; i < xnumberAdded + 1; ++i) {
			for (int j = 0; j < ynumberAdded; ++j) {
				for (int k = 0; k < znumberAdded; ++k) {
					for (int l = 0; l < lnumber; ++l) {
						double sum = rightPart[i][j][k][l];
						double a = 1;
						for (int m = 0; m < matrix[i][j][k][l].size(); ++m) {
							MatrixElement element = matrix[i][j][k][l][m];
							if (!indexEqual(element, i, j, k, l)) {
								sum -= element.value * outVector[element.i][element.j][element.k][element.l];
							} else {
								a = element.value;
							}
						}
						outVector[i][j][k][l] = sum / a;
					}
				}
			}
		}
		exchangeLargeVector(outVector, xnumberAdded, ynumberAdded, znumberAdded, lnumber, additionalBinNumber, periodicX, periodicY, periodicZ, cartComm, cartCoord, cartDim, leftOutBuffer, rightOutBuffer, leftInBuffer, rightInBuffer, frontOutBuffer, backOutBuffer, frontInBuffer, backInBuffer, bottomOutBuffer, topOutBuffer, bottomInBuffer, topInBuffer);

		curIteration++;
	}

	delete[] leftOutBuffer;
	delete[] rightInBuffer;
	delete[] rightOutBuffer;
	delete[] leftInBuffer;

	delete[] frontOutBuffer;
	delete[] frontInBuffer;
	delete[] backOutBuffer;
	delete[] backInBuffer;

	delete[] bottomOutBuffer;
	delete[] bottomInBuffer;
	delete[] topOutBuffer;
	delete[] topInBuffer;
}
*/

bool indexLower(MatrixElement* element, int i, int j, int k, int l) {
    if (element->i < i) {
		return true;
	}
    if (element->i > i) {
		return false;
	}
    if (element->j < j) {
		return true;
	}
    if (element->j > j) {
		return false;
	}
    if (element->k < k) {
		return true;
	}
    if (element->k > k) {
		return false;
	}
    if (element->l < l) {
		return true;
	}
    if (element->l > l) {
		return false;
	}

	return false;
}

bool indexEqual(MatrixElement* element, int i, int j, int k, int l) {
    return (element->i == i) && (element->j == j) && (element->k == k) && (element->l == l);
}

bool indexUpper(MatrixElement* element, int i, int j, int k, int l) {
	return ((!indexEqual(element, i, j, k, l)) && (!indexLower(element, i, j, k, l)));
}


void conjugateGradientMethod(Grid* grid, MatrixElementNode***** matrix, double**** rightPart, double**** outVector, int lnumber,
                             double precision, int maxIteration, int verbosity, bool periodicX, bool periodicY, bool periodicZ) {
    int rank = 0;
    int nprocs;
    int i,j,k;
#ifdef PARALLEL
    MPI_Comm cartComm = MPI_COMM_WORLD;
    MPI_Comm_size(cartComm, &nprocs);
    MPI_Comm_rank(cartComm, &rank);

    MPI_Barrier(cartComm);
#endif
    int  par_dim[3] = {0, 0, 0};
    DIM_EXPAND(par_dim[0] = grid->nproc[IDIR] > 1;  ,
               par_dim[1] = grid->nproc[JDIR] > 1;  ,
               par_dim[2] = grid->nproc[KDIR] > 1;)

    if ((rank == 0) && (verbosity > 0)) printf("start conjugate gradient\n");

    double norm = sqrt(scalarMultiplyLargeVectors(rightPart, rightPart, lnumber));
    int matrixDimension = lnumber * grid->np_int_glob[0] * grid->np_int_glob[1]*grid->np_int_glob[2];
    //printf("right part norm = %g\n", norm);
    //alertNaNOrInfinity(norm, "right partnorm = NaN in gmres\n");

    if (norm == 0) {
        DOM_LOOP(k,j,i){
                    for (int l = 0; l < lnumber; ++l) {
                        outVector[k][j][i][l] = 0;
                    }
        }
        return;
    }

        double**** residual = (double****) malloc(NX3_TOT*sizeof(double***));
        double**** prevResidual = (double****) malloc(NX3_TOT*sizeof(double***));
        double**** z = (double****) malloc(NX3_TOT*sizeof(double***));
        double**** tempVector = (double****) malloc(NX3_TOT*sizeof(double***));
        for (int i = 0; i < NX3_TOT; ++i) {
            residual[i] = (double***) malloc(NX2_TOT*sizeof(double**));
            prevResidual[i] = (double***) malloc(NX2_TOT*sizeof(double**));
            z[i] = (double***) malloc(NX2_TOT*sizeof(double**));
            tempVector[i] = (double***) malloc(NX2_TOT*sizeof(double**));
            for (int j = 0; j < NX2_TOT; ++j) {
                residual[i][j] = (double**) malloc(NX1_TOT*sizeof(double*));
                prevResidual[i][j] = (double**) malloc(NX1_TOT*sizeof(double*));
                z[i][j] = (double**) malloc(NX1_TOT*sizeof(double*));
                tempVector[i][j] = (double**) malloc(NX1_TOT*sizeof(double*));
                for(int k = 0; k < NX1_TOT; ++k){
                    residual[i][j][k] = (double*) malloc(lnumber*sizeof(double));
                    prevResidual[i][j][k] = (double*) malloc(lnumber*sizeof(double));
                    z[i][j][k] = (double*) malloc(lnumber*sizeof(double));
                    tempVector[i][j][k] = (double*) malloc(lnumber*sizeof(double));
                    for(int l = 0; l < lnumber; ++l){
                        residual[i][j][k][l] = rightPart[i][j][k][l];
                        outVector[i][j][k][l] = 0;
                        prevResidual[i][j][k][l] = rightPart[i][j][k][l];
                        z[i][j][k][l] = rightPart[i][j][k][l];
                        tempVector[i][j][k][l] = 0;
                    }
                }
            }
        }


    /*for (int i = 0; i < xnumber; ++i) {
        for (int j = 0; j < ynumber; ++j) {
            for (int k = 0; k < znumber; ++k) {
                for (int l = 0; l < lnumber; ++l) {
                    for (int m = 0; m < matrix[i][j][k][l].size(); ++m) {
                        MatrixElement element = matrix[i][j][k][l][m];
                        prevResidual[i][j][k][l] += element.value * outVector[element.i][element.j][element.k][element.l];
                    }
                }
            }
        }
    }*/

    int iteration = 0;


    double prevResidualNorm2 = scalarMultiplyLargeVectors(prevResidual, prevResidual, lnumber);
    double residualNorm2 = prevResidualNorm2;
    double rightPartNorm2 = scalarMultiplyLargeVectors(rightPart, rightPart, lnumber);

    double relativeError = sqrt(residualNorm2 / rightPartNorm2);

    while ((iteration < maxIteration) && (iteration < matrixDimension) && (
        relativeError > (precision / matrixDimension))) {
        if ((rank == 0) && (verbosity > 1)) printf("conjugate gradient iteration %d\n", iteration);

        multiplySpecialMatrixVector(tempVector, matrix, z, lnumber, par_dim, periodicX, periodicY, periodicZ);

        double alpha = prevResidualNorm2 / scalarMultiplyLargeVectors(tempVector, z, lnumber);

        int i, j, k;

        DOM_LOOP(i,j,k){
                    for (int l = 0; l < lnumber; ++l) {
                        outVector[i][j][k][l] += alpha * z[i][j][k][l];
                        residual[i][j][k][l] = prevResidual[i][j][k][l] - alpha * tempVector[i][j][k][l];
                    }
        }
        exchangeLargeVector(outVector, lnumber, par_dim, SZ_stagx, periodicX, periodicY, periodicZ);
        exchangeLargeVector(residual, lnumber, par_dim, SZ_stagx, periodicX, periodicY, periodicZ);

        residualNorm2 = scalarMultiplyLargeVectors(residual, residual, lnumber);

        double beta = residualNorm2 / prevResidualNorm2;

        DOM_LOOP(i,j,k){
                    for (int l = 0; l < lnumber; ++l) {
                        z[i][j][k][l] = residual[i][j][k][l] + beta * z[i][j][k][l];
                    }
        }
        exchangeLargeVector(z, lnumber, par_dim, SZ_stagx, periodicX, periodicY, periodicZ);

        prevResidualNorm2 = residualNorm2;

        relativeError = sqrt(residualNorm2 / rightPartNorm2);
        iteration++;
    }
    if((rank == 0) && (verbosity > 0)) printf("total conjugate gradient iteration = %d\n", iteration);
    if (rank == 0) printf("relative conjugate gradient error = %g\n", relativeError);
    if (rank == 0) printf("relative conjugate gradient error*Nx*Ny*Nz*Np = %g\n", relativeError*matrixDimension);

    for (int i = 0; i < NX3_TOT; ++i) {
        for (int j = 0; j < NX2_TOT; ++j) {
            for (int k = 0; k < NX1_TOT; ++k) {
                free(residual[i][j][k]);
                free(prevResidual[i][j][k]);
                free(z[i][j][k]);
                free(tempVector[i][j][k]);
            }
            free(residual[i][j]);
            free(prevResidual[i][j]);
            free(z[i][j]);
            free(tempVector[i][j]);
        }
        free(residual[i]);
        free(prevResidual[i]);
        free(z[i]);
        free(tempVector[i]);
    }
    free(residual);
    free(prevResidual);
    free(z);
    free(tempVector);
}



