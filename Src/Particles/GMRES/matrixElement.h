#ifndef _MATRIX_ELEMENT_H_
#define _MATRIX_ELEMENT_H_

#include "stdbool.h"

typedef struct MatrixElement_
{
	double value;
    int k;
    int j;
    int i;
    int l;
} MatrixElement;

MatrixElement createMatrixElement(double v, int kv, int jv, int iv, int lv);

bool equalsIndex(MatrixElement element1, MatrixElement element2);

typedef struct MatrixElementNode_{
    struct MatrixElement_ element;
    struct MatrixElementNode_* next;
    struct MatrixElementNode_* prev;
} MatrixElementNode;

MatrixElementNode* addElement(MatrixElementNode* curNode, double value, int k, int j, int i, int l);

#endif
