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

void generalizedMinimalResidualMethod(Grid* grid, MatrixElementNode***** matrix, double ****rightPart, double ****outvector, LargeVectorBasis* gmresBasis, int lnumber, double precision,
                                      int maxIteration, int verbosity);
void arnoldiIterations(MatrixElementNode*****matrix, double **outHessenbergMatrix, int n,
                       LargeVectorBasis* gmresBasis, double **prevHessenbergMatrix, int lnumber, int* par_dim);
//double**** multiplySpecialMatrixVector(MatrixElementNode***** matrix, double**** vector, int lnumber);

void multiplySpecialMatrixVector(double**** result, MatrixElementNode***** matrix, double**** vector, int lnumber, int* par_dims);
double evaluateError(double** hessenbergMatrix, double* vector, double beta, int n);

double scalarMultiplyLargeVectors(double ****a, double ****b, int lnumber);


double normDifferenceLargeVectors(double**** a, double**** b, int lnumber);

void conjugateGradientMethod(Grid* grid, MatrixElementNode***** matrix, double**** rightPart, double**** outVector, int lnumber,
                             double precision, int maxIteration, int verbosity);


//void transposeSpecialMatrix(MatrixElementNode***** result, MatrixElementNode***** matrix);



bool indexLower(MatrixElement* element, int i, int j, int k, int l);
bool indexEqual(MatrixElement* element, int i, int j, int k, int l);
bool indexUpper(MatrixElement* element, int i, int j, int k, int l);

#endif
