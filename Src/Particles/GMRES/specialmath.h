#ifndef _SPECIAL_MATH_H_
#define _SPECIAL_MATH_H_

#include "stdlib.h"
#include "stdio.h"
#include "stdbool.h"
#include "macros.h"
//#include "structs.h"

struct Grid;

#include "matrixElement.h"
#include "largeVectorBasis.h"

void generalizedMinimalResidualMethod1(Grid* grid, MatrixElementNode***** matrix, double ****rightPart, double ****outvector, LargeVectorBasis* gmresBasis, int lnumber, double precision,
                                      int maxIteration, int verbosity, bool periodicX, bool periodicY, bool periodicZ);
void generalizedMinimalResidualMethod(Grid* grid, MatrixElementNode***** matrix, double ****rightPart, double ****outvector, double**** initialVector, LargeVectorBasis* gmresBasis, int lnumber, double precision,
                                      int maxIteration, int verbosity, bool periodicX, bool periodicY, bool periodicZ);
void arnoldiIterations(MatrixElementNode*****matrix, double **outHessenbergMatrix, int n,
                       LargeVectorBasis* gmresBasis, double **prevHessenbergMatrix, int lnumber, int* par_dim, bool periodicX, bool periodicY, bool periodicZ);
//double**** multiplySpecialMatrixVector(MatrixElementNode***** matrix, double**** vector, int lnumber);

void multiplySpecialMatrixVector(double**** result, MatrixElementNode***** matrix, double**** vector, int lnumber, int* par_dims, bool periodicX, bool periodicY, bool periodicZ);
double evaluateError(double** hessenbergMatrix, double* vector, double beta, int n);

double scalarMultiplyLargeVectors(double ****a, double ****b, int lnumber);


double normDifferenceLargeVectors(double**** a, double**** b, int lnumber);

void conjugateGradientMethod(Grid* grid, MatrixElementNode***** matrix, double**** rightPart, double**** outVector, int lnumber,
                             double precision, int maxIteration, int verbosity, bool periodicX, bool periodicY, bool periodicZ);


void transposeSpecialMatrix(MatrixElementNode***** result, MatrixElementNode***** matrix, int lnumber);

void biconjugateGradientMethod(Grid* grid, MatrixElementNode***** matrix, double ****rightPart, double ****outvector, int lnumber, double precision,
                               int maxIteration, int verbosity, bool periodicX, bool periodicY, bool periodicZ);

void biconjugateStabilizedGradientMethod(Grid* grid, MatrixElementNode***** matrix, double**** rightPart,
                                         double**** outVector, int lnumber, double precision, int maxIteration, int verbosity, bool periodicX, bool periodicY, bool periodicZ);

bool indexLower(MatrixElement* element, int i, int j, int k, int l);
bool indexEqual(MatrixElement* element, int i, int j, int k, int l);
bool indexUpper(MatrixElement* element, int i, int j, int k, int l);

#endif
