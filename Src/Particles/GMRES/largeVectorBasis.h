#ifndef LARGEVECTORBASiS_H
#define LARGEVECTORBASIS_H

//#include "pluto.h"
#include <stdbool.h>

typedef struct LargeVectorBasis_ {
	int size;
    int znumber;
	int ynumber;
    int xnumber;
    int lnumber;
	int capacity;
    double***** array;
} LargeVectorBasis;

LargeVectorBasis createLargeVectorBasis(int sizev, int znumberv, int ynumberv, int xnumberv, int lnumberv);
void createLargeVectorBasis1(LargeVectorBasis* basis, int sizev, int znumberv, int ynumberv, int xnumberv, int lnumberv);
void resize(LargeVectorBasis* basis, int capacityv);
void clear(LargeVectorBasis* basis);

void exchangeLargeVector(double**** vector, int lnumber, int *dims, int sz_ptr, bool periodicX, bool preiodicY, bool periodicZ);

#endif
