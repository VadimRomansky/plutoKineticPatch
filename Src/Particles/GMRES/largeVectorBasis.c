#include "stdlib.h"
#include "stdio.h"
//#include <crtdbg.h>
#include <stdbool.h>

//#include "memory_debug.h"
#include "largeVectorBasis.h"

#ifdef PARALLEL
#include "al_hidden.h"
#endif

#include "macros.h"
#include "pluto.h"

#ifdef PARALLEL
extern struct SZ *sz_stack[AL_MAX_ARRAYS];
extern int stack_ptr[AL_MAX_ARRAYS];
#endif

LargeVectorBasis createLargeVectorBasis(int capacityv, int znumberv, int ynumberv, int xnumberv, int lnumberv) {
    LargeVectorBasis basis;
    basis.size = 0;
    basis.znumber = znumberv;
    basis.ynumber = ynumberv;
    basis.xnumber = xnumberv;
    basis.lnumber = lnumberv;
    basis.capacity = capacityv;

    basis.array = (double *****) malloc(capacityv*sizeof(double****));
    for (int m = 0; m < capacityv; ++m) {
        basis.array[m] = (double****) malloc(znumberv*sizeof(double***));
        for (int i = 0; i < znumberv; ++i) {
            basis.array[m][i] = (double***) malloc(ynumberv*sizeof(double**));
            for (int j = 0; j < ynumberv; ++j) {
                basis.array[m][i][j] = (double**) malloc(xnumberv*sizeof(double*));
                for(int k = 0; k < xnumberv; ++k){
                    basis.array[m][i][j][k] = (double*) malloc(lnumberv*sizeof(double));
                    for(int l = 0; l < lnumberv; ++l){
                        basis.array[m][i][j][k][l] = 0;
                    }
                }
            }
        }
    }
    return basis;
}

void createLargeVectorBasis1(LargeVectorBasis* basis, int capacityv, int znumberv, int ynumberv, int xnumberv, int lnumberv){
    basis->size = 0;
    basis->znumber = znumberv;
    basis->ynumber = ynumberv;
    basis->xnumber = xnumberv;
    basis->lnumber = lnumberv;
    basis->capacity = capacityv;

    basis->array = (double *****) malloc(capacityv*sizeof(double****));
    for (int m = 0; m < capacityv; ++m) {
        basis->array[m] = (double****) malloc(znumberv*sizeof(double***));
        for (int i = 0; i < znumberv; ++i) {
            basis->array[m][i] = (double***) malloc(ynumberv*sizeof(double**));
            for (int j = 0; j < ynumberv; ++j) {
                basis->array[m][i][j] = (double**) malloc(xnumberv*sizeof(double*));
                for(int k = 0; k < xnumberv; ++k){
                    basis->array[m][i][j][k] = (double*) malloc(lnumberv*sizeof(double));
                    for(int l = 0; l < lnumberv; ++l){
                        basis->array[m][i][j][k][l] = 0;
                    }
                }
            }
        }
    }
    return basis;
}

void resize(LargeVectorBasis* basis, int capacityv) {
    if (capacityv > basis->capacity) {
        double***** newArray = (double *****) malloc(capacityv*sizeof(double****));
        for (int m = 0; m < basis->capacity; ++m) {
            newArray[m] = basis->array[m];
        }
        for (int m = basis->capacity; m < capacityv; ++m) {
            newArray[m] = (double****) malloc(basis->znumber*sizeof(double***));
            for (int i = 0; i < basis->znumber; ++i) {
                newArray[m][i] = (double***) malloc(basis->ynumber*sizeof(double**));
                for (int j = 0; j < basis->ynumber; ++j) {
                    newArray[m][i][j] = (double**) malloc(basis->xnumber*sizeof(double*));
                    for(int k = 0; k < basis->xnumber; ++k){
                        newArray[m][i][j][k] = (double*) malloc(basis->lnumber*sizeof(double));
                        for(int l = 0; l < basis->lnumber; ++l){
                            newArray[m][i][j][k][l] = 0.0;
                        }
                    }
                }
            }
        }

        basis->capacity = capacityv;
        free(basis->array);
        basis->array = newArray;
        return;
    }
    if (capacityv >= basis->size) {
        double***** newArray = (double*****) malloc(capacityv*sizeof(double****));
        for (int m = 0; m < capacityv; ++m) {
            newArray[m] = basis->array[m];
        }
        for (int m = capacityv; m < basis->capacity; ++m) {
            for (int i = 0; i < basis->znumber; ++i) {
                for(int j = 0; j < basis->ynumber; ++j){
                    for(int k = 0; k < basis->xnumber; ++k){
                        free(basis->array[m][i][j][k]);
                    }
                    free(basis->array[m][i][j]);
                }
                free(basis->array[m][i]);
            }
            free(basis->array[m]);
        }
        free(basis->array);
        basis->array = newArray;
        basis->capacity = capacityv;
        return;
    }
    printf("capacity < size\n");
    exit(0);
}

void clear(LargeVectorBasis* basis) {
    basis->size = 0;
    for (int m = 0; m < basis->capacity; ++m) {
        for (int i = 0; i < basis->znumber; ++i) {
            for (int j = 0; j < basis->ynumber; ++j) {
                for(int k = 0; k < basis->xnumber; ++k){
                    for(int l = 0; l < basis->lnumber; ++l){
                        basis->array[m][i][j][k][l] = 0;
                    }
                    //free(basis->array[m][i][j][k]);
                }
                //free(basis->array[m][i][j]);
            }
            //free(basis->array[m][i]);
        }
        //free(basis->array[m]);
    }
    //free(basis->array);
    //basis->capacity = 0;
}

void exchangeLargeVector(double**** vector, int lnumber, int *dims, int sz_ptr, bool periodicX, bool periodicY, bool periodicZ){
    int i,j,k;
#ifdef PARALLEL
    register int nd;
    int myrank, nproc;
    int ndim, gp, nleft, nright, tag1, tag2;
    int sendb, recvb;
    MPI_Datatype itype;
    MPI_Comm comm;
    MPI_Status status;
    struct szz *s;

    s = sz_stack[sz_ptr];

    myrank = s->rank;
    nproc = s->size;
    comm = s->comm;
    ndim = s->ndim;

    for(nd=0;nd<ndim;nd++){
        gp = s->bg[nd];

        /* If gp=0, do nothing */
        if( gp > 0){
            if(dims[nd] != 0 ){
                nleft = s->left[nd];
                nright = s->right[nd];
                itype = s->type_rl[nd];
                tag1 = s->tag1[nd];

                sendb = s->sendb1[nd];
                recvb = s->recvb1[nd];

                double* outBuffer;
                double* inBuffer;

                int send_count;

                if(nd == 0){
                    send_count = gp*NX2_TOT*NX3_TOT*lnumber;
                    outBuffer = (double*) malloc(sizeof(double)*send_count);
                    inBuffer = (double*) malloc(sizeof(double)*send_count);
                    int count = 0;
                    KTOT_LOOP(k){
                        JTOT_LOOP(j){
                            for(i = 0; i < IBEG; ++i){
                                for(int l = 0; l < lnumber; ++l){
                                    outBuffer[count] = vector[k][j][IBEG + i][l];
                                    count++;
                                }
                            }
                        }
                    }
                } else if (nd == 1){
                    send_count = gp*NX1_TOT*NX3_TOT*lnumber;
                    outBuffer = (double*) malloc(sizeof(double)*send_count);
                    inBuffer = (double*) malloc(sizeof(double)*send_count);
                    int count = 0;
                    KTOT_LOOP(k){
                        for(j = 0; j < JBEG; ++j){
                            ITOT_LOOP(i){
                                for(int l = 0; l < lnumber; ++l){
                                    outBuffer[count] = vector[k][JBEG + j][i][l];
                                    count++;
                                }
                            }
                        }
                    }
                } else if (nd == 2){
                    send_count = gp*NX1_TOT*NX2_TOT*lnumber;
                    outBuffer = (double*) malloc(sizeof(double)*send_count);
                    inBuffer = (double*) malloc(sizeof(double)*send_count);
                    int count = 0;
                    for(k = 0; k < KBEG; ++k){
                        JTOT_LOOP(j){
                            ITOT_LOOP(i){
                                for(int l = 0; l < lnumber; ++l){
                                    outBuffer[count] = vector[KBEG + k][j][i][l];
                                    count++;
                                }
                            }
                        }
                    }
                }


                MPI_Sendrecv(outBuffer, send_count, MPI_DOUBLE, nleft, tag1,
                             inBuffer, send_count, MPI_DOUBLE, nright,tag1,
                             comm, &status);

                if(nd == 0){
                    int count = 0;
                    KTOT_LOOP(k){
                        JTOT_LOOP(j){
                            for(i = 0; i < IBEG; ++i){
                                for(int l = 0; l < lnumber; ++l){
                                    vector[k][j][IEND + i + 1][l] = inBuffer[count];
                                    outBuffer[count] = vector[k][j][IEND - IBEG + i + 1][l];
                                    count++;
                                }
                            }
                        }
                    }
                } else if(nd == 1){
                    int count = 0;
                    KTOT_LOOP(k){
                        for(j = 0; j < JBEG; ++j){
                            ITOT_LOOP(i){
                                for(int l = 0; l < lnumber; ++l){
                                    vector[k][JEND + j + 1][i][l] = inBuffer[count];
                                    outBuffer[count] = vector[k][JEND - JBEG + j + 1][i][l];
                                    count++;
                                }
                            }
                        }
                    }
                } else if(nd == 2){
                    int count = 0;
                    for(k = 0; k < KBEG; ++k){
                        JTOT_LOOP(j){
                            ITOT_LOOP(i){
                                for(int l = 0; l < lnumber; ++l){
                                    vector[KEND + k + 1][j][i][l] = inBuffer[count];
                                    outBuffer[count] = vector[KEND - KBEG + k + 1][j][i][l];
                                    count++;
                                }
                            }
                        }
                    }
                }

                nleft = s->left[nd];
                nright = s->right[nd];
                itype = s->type_lr[nd];
                tag2 = s->tag2[nd];


                sendb = s->sendb2[nd];
                recvb = s->recvb2[nd];

                MPI_Sendrecv(outBuffer, send_count, MPI_DOUBLE, nright, tag2,
                             inBuffer, send_count, MPI_DOUBLE, nleft,tag2,
                             comm, &status);

                if(nd == 0){
                    int count = 0;
                    KTOT_LOOP(k){
                        JTOT_LOOP(j){
                            for(i = 0; i < IBEG; ++i){
                                for(int l = 0; l < lnumber; ++l){
                                    vector[k][j][i][l] = inBuffer[count];
                                    count++;
                                }
                            }
                        }
                    }
                } else if(nd == 1){
                    int count = 0;
                    KTOT_LOOP(k){
                        for(j = 0; j < JBEG; ++j){
                            ITOT_LOOP(i){
                                for(int l = 0; l < lnumber; ++l){
                                    vector[k][j][i][l] = inBuffer[count];
                                    count++;
                                }
                            }
                        }
                    }
                } else if(nd == 2){
                    int count = 0;
                    for(k = 0; k < KBEG; ++k){
                        JTOT_LOOP(j){
                            ITOT_LOOP(i){
                                for(int l = 0; l < lnumber; ++l){
                                    vector[k][j][i][l] = inBuffer[count];
                                    count++;
                                }
                            }
                        }
                    }
                }

                free(outBuffer);
                free(inBuffer);
            } else {
                if(nd == 0){
                    KTOT_LOOP(k){
                        JTOT_LOOP(j){
                            for(i = 0; i < IBEG; ++i){
                                for(int l = 0; l < lnumber; ++l){
                                    vector[k][j][IEND + i + 1][l] = vector[k][j][IBEG + i][l];
                                    vector[k][j][i][l] = vector[k][j][IEND - IBEG + i + 1][l];
                                }
                            }
                        }
                    }
                } else if(nd == 1){
                    KTOT_LOOP(k){
                        for(j = 0; j < JBEG; ++j){
                            ITOT_LOOP(i){
                                for(int l = 0; l < lnumber; ++l){
                                    vector[k][JEND + j + 1][i][l] = vector[k][JBEG + j][i][l];
                                    vector[k][j][i][l] = vector[k][JEND - JBEG + j + 1][i][l];
                                }
                            }
                        }
                    }
                } else if(nd == 2){
                    for(k = 0; k < KBEG; ++k){
                        JTOT_LOOP(j){
                            ITOT_LOOP(i){
                                for(int l = 0; l < lnumber; ++l){
                                    vector[KEND + k + 1][j][i][l] = vector[KBEG + k][j][i][l];
                                    vector[k][j][i][l] = vector[KEND - KBEG + k + 1][j][i][l];
                                }
                            }
                        }
                    }
                }
            }
        }
    }
#else
#if INCLUDE_IDIR
    if(periodicX){
    KTOT_LOOP(k){
        JTOT_LOOP(j){
            for(i = 0; i < IBEG; ++i){
                for(int l = 0; l < lnumber; ++l){
                    vector[k][j][IEND + i + 1][l] = vector[k][j][IBEG + i][l];
                    vector[k][j][i][l] = vector[k][j][IEND - IBEG + i + 1][l];
                }
            }
        }
    }
    }
#endif

#if INCLUDE_JDIR
    if(periodicY){
    KTOT_LOOP(k){
        for(j = 0; j < JBEG; ++j){
            ITOT_LOOP(i){
                for(int l = 0; l < lnumber; ++l){
                    vector[k][JEND + j + 1][i][l] = vector[k][JBEG + j][i][l];
                    vector[k][j][i][l] = vector[k][JEND - JBEG + j + 1][i][l];
                }
            }
        }
    }
    }
#endif

#if INCLUDE_KDIR
    if(periodicZ){
    for(k = 0; k < KBEG; ++k){
        JTOT_LOOP(j){
            ITOT_LOOP(i){
                for(int l = 0; l < lnumber; ++l){
                    vector[KEND + k + 1][j][i][l] = vector[KBEG + k][j][i][l];
                    vector[k][j][i][l] = vector[KEND - KBEG + k + 1][j][i][l];
                }
            }
        }
    }
    }
#endif

#endif
}

