//#include <crtdbg.h>

//#include "memory_debug.h"
#include <stdio.h>
#include "matrixElement.h"


MatrixElement createMatrixElement(double v, int kv, int jv, int iv, int lv) {
    MatrixElement m;
    m.value = v;
    m.i = iv;
    m.j = jv;
    m.k = kv;
    m.l = lv;
    return m;
}

bool equalsIndex(MatrixElement element1, MatrixElement element2) {
    if (element1.i != element2.i) return false;
    if (element1.j != element2.j) return false;
    if(element1.k != element2.k) return false;
    if(element1.l != element2.l) return false;

	return true;
}

MatrixElementNode* addElement(MatrixElementNode* curNode, double value, int k, int j, int i, int l){
    MatrixElementNode* tempNode = (MatrixElementNode*) malloc(sizeof(MatrixElementNode));
    tempNode->element = createMatrixElement(value, k,j,i,l);
    tempNode->next = curNode->next;
    if(curNode->next != NULL){
        curNode->prev = tempNode;
    }
    tempNode->prev = curNode;
    curNode->next = tempNode;
    return tempNode;
}
