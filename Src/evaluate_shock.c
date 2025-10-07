#include <stdio.h>
#include <stdbool.h>

#include "definitions.h"

#ifdef PARALLEL
#include "al_hidden.h"
#endif

#ifdef PARALLEL
extern struct SZ *sz_stack[AL_MAX_ARRAYS];
extern int stack_ptr[AL_MAX_ARRAYS];
#endif

#include "pluto.h"

typedef struct CellTracerNode_{
    int rank0;
    int i0;
    int j0;
    int k0;
    double x1;
    double x2;
    double x3;
    int i;
    int j;
    int k;
    double v1;
    double v2;
    double v3;
    double rho;
    double pressure;
    double prevgradx;
    double prevgrady;
    double prevgradz;
    int N;
    struct CellTracerNode_* next;
    struct CellTracerNode_* prev;
} CellTracerNode;

CellTracerNode* createTracer(int vrank0, int vi0, int vj0, int vk0, double x1v, double x2v, double x3v, int vi, int vj, int vk, double vv1, double vv2, double vv3, double vrho, double vpressure, double gradx, double grady, double gradz, int vN){
    CellTracerNode* tempNode = (CellTracerNode*) malloc(sizeof(CellTracerNode));
    CheckNanOrInfinity(x1v, "x1 = NaN\n");
    CheckNanOrInfinity(x2v, "x2 = NaN\n");
    CheckNanOrInfinity(x3v, "x3 = NaN\n");
    CheckNanOrInfinity(vv1, "v1 = NaN\n");
    CheckNanOrInfinity(vv2, "v2 = NaN\n");
    CheckNanOrInfinity(vv3, "v3 = NaN\n");
    CheckNanOrInfinity(vrho, "rho = Nan\n");
    tempNode->rank0 = vrank0;
    tempNode->i0 = vi0;
    tempNode->j0 = vj0;
    tempNode->k0 = vk0;
    tempNode->x1 = x1v;
    tempNode->x2 = x2v;
    tempNode->x3 = x3v;
    tempNode->i = vi;
    tempNode->j = vj;
    tempNode->k = vk;
    tempNode->v1 = vv1;
    tempNode->v2 = vv2;
    tempNode->v3 = vv3;
    tempNode->rho = vrho;
    tempNode->pressure = vpressure;
    tempNode->prevgradx = gradx;
    tempNode->prevgrady = grady;
    tempNode->prevgradz = gradz;
    tempNode->N = vN;
    tempNode->next = NULL;
    tempNode->prev = NULL;
    return tempNode;
}

CellTracerNode* addElementAfter(CellTracerNode* curNode, int vrank0, int vi0, int vj0, int vk0, double x1v, double x2v, double x3v, int vi, int vj, int vk, double vv1, double vv2, double vv3, double vrho, double vpressure, double gradx, double grady, double gradz, int N){
    CellTracerNode* tempNode = createTracer(vrank0, vi0, vj0, vk0, x1v, x2v, x3v, vi, vj, vk, vv1, vv2, vv3, vrho, vpressure, gradx, grady, gradz, N);
    tempNode->next = curNode->next;
    if(curNode->next != NULL){
        curNode->prev = tempNode;
    }
    tempNode->prev = curNode;
    curNode->next = tempNode;
    return tempNode;
}

void putTracerListToArray(CellTracerNode* tracersList, int* outbuf, double* outbufd){
    int count = 0;
    int countd = 0;
    while(tracersList != NULL){
        outbuf[count] = tracersList->rank0;
        count++;
        outbuf[count] = tracersList->i0;
        count++;
        outbuf[count] = tracersList->j0;
        count++;
        outbuf[count] = tracersList->k0;
        count++;
        outbuf[count] = tracersList->i;
        count++;
        outbuf[count] = tracersList->j;
        count++;
        outbuf[count] = tracersList->k;
        count++;
        outbuf[count] = tracersList->N;
        count++;

        outbufd[countd] = tracersList->x1;
        countd++;
        outbufd[countd] = tracersList->x2;
        countd++;
        outbufd[countd] = tracersList->x3;
        countd++;
        outbufd[countd] = tracersList->v1;
        countd++;
        outbufd[countd] = tracersList->v2;
        countd++;
        outbufd[countd] = tracersList->v3;
        countd++;
        outbufd[countd] = tracersList->rho;
        countd++;
        outbufd[countd] = tracersList->pressure;
        countd++;
        outbufd[countd] = tracersList->prevgradx;
        countd++;
        outbufd[countd] = tracersList->prevgrady;
        countd++;
        outbufd[countd] = tracersList->prevgradz;
        countd++;

        CellTracerNode* temp = tracersList;
        tracersList = tracersList->next;
        if(tracersList != NULL){
            tracersList->prev = NULL;
        }
        free(temp);
    }
}

CellTracerNode* putArrayToTracerList(int* inbuf, double* inbufd, int Nin){
	int count = 0;
	int countd = 0;
	CellTracerNode* list = NULL;
	for(int l = 0; l < Nin; ++l){
		int rank0 = inbuf[count];
		count++;
		int i0 = inbuf[count];
		count++;
		int j0 = inbuf[count];
		count++;
		int k0 = inbuf[count];
		count++;
		int i = inbuf[count];
		count++;
		int j = inbuf[count];
		count++;
		int k = inbuf[count];
		count++;
        int N = inbuf[count];
        count++;
		double x = inbufd[countd];
		countd++;
		double y = inbufd[countd];
		countd++;
		double z = inbufd[countd];
		countd++;
		double vx = inbufd[countd];
		countd++;
		double vy = inbufd[countd];
		countd++;
		double vz = inbufd[countd];
		countd++;
		double rho = inbufd[countd];
		countd++;
        double pressure = inbufd[countd];
        countd++;
        double gradx = inbufd[countd];
        countd++;
        double grady = inbufd[countd];
        countd++;
        double gradz = inbufd[countd];
        countd++;
		
        CellTracerNode* temp = createTracer(rank0, i0, j0, k0, x, y, z, i, j, k, vx, vy, vz, rho, pressure, gradx, grady, gradz, N);
		
		temp->next = list;
		if(list != NULL){
			list->prev = temp;
		}
		list = temp;
	}
	return list;
}

#ifdef PARALLEL

void traceShockParallel(Data* d, Grid* grid, int direction, double*** x1, double*** x2, double*** x3, double*** v1, double*** v2, double*** v3, double*** rho, double*** pressure){
    register int nd;
    int i, j, k;
    int myrank, nprocs;
    int globrank;
    int ndim, gp, nleft, nright, tag1, tag2;
    int sendb, recvb;
    MPI_Datatype itype;
    MPI_Comm comm;
    MPI_Status status;
    struct szz *s;
    int sz_ptr = SZ_stagx;

    s = sz_stack[sz_ptr];

    ndim = s->ndim;
    globrank = s->rank;
    comm = s->comm;

    CellTracerNode* tracers = NULL;
    CellTracerNode* stoppedTracers = NULL;
    int NactiveTracers[1];
#if INCLUDE_IDIR
    int NtoLeft = 0;
    int NtoRight = 0;
    CellTracerNode* tracersToLeft = NULL;
    CellTracerNode* tracersToRight = NULL;
#endif
#if INCLUDE_JDIR
    int NtoUp = 0;
    int NtoDown = 0;
    CellTracerNode* tracersToUp = NULL;
    CellTracerNode* tracersToDown = NULL;
#endif
#if INCLUDE_KDIR
    int NtoBack = 0;
    int NtoFront = 0;
    CellTracerNode* tracersToFront = NULL;
    CellTracerNode* tracersToBack = NULL;
#endif
    NactiveTracers[0] = 0;
    DOM_LOOP(k,j,i){
        if(!(d->flag[k][j][i] & FLAG_ENTROPY)){
            CellTracerNode* tempTracer = createTracer(globrank, i,j,k, grid->x[0][i], grid->x[1][j], grid->x[2][k], i,j,k, 0.0, 0.0, 0.0, d->Vc[RHO][k][j][i], d->Vc[RHO][k][j][i], 0, 0, 0, 0);
            if(tracers == NULL){
                tracers = tempTracer;
            } else {
                tempTracer->next = tracers;
                tracers->prev = tempTracer;
                tracers = tempTracer;
            }
            NactiveTracers[0] = NactiveTracers[0] + 1;
        }
    }
    int a[1];
    a[0] = NactiveTracers[0];
    MPI_Allreduce(a, NactiveTracers, 1, MPI_INT, MPI_SUM, comm);

    while(NactiveTracers[0] > 0){
#if INCLUDE_IDIR
        NtoLeft = 0;
        NtoRight = 0;
#endif
#if INCLUDE_JDIR
        NtoDown = 0;
        NtoUp = 0;
#endif
#if INCLUDE_KDIR
        NtoBack = 0;
        NtoFront = 0;
#endif
        while(tracers != NULL){
            int currenti = tracers->i;
            int currentj = tracers->j;
            int currentk = tracers->k;
            double x = tracers->x1;
            double y = tracers->x2;
            double z = tracers->x3;

            bool stopped = true;

            int countn = 0;

            while(!(d->flag[currentk][currentj][currenti] & FLAG_ENTROPY)){
                countn++;
                if(tracers->N > 10){
                    printf("n traced cells = %d, x = %g, y = %g, z = %g, current i = %d, current j = %d, current k = %d, rank = %d, prevgradx = %g, prevgrady = %g, prevgradz = %g\n", tracers->N, tracers->x1, tracers->x2, tracers->x3, tracers->i, tracers->j, tracers->k, globrank, tracers->prevgradx, tracers->prevgrady, tracers->prevgradz);
                }
                if(tracers->N > 1000){
                    printf("n traced cells > 1000\n");
                    printLog("n traced cells > 1000\n");
                    QUIT_PLUTO(0);
                }
                tracers->N = tracers->N + 1;
#if INCLUDE_IDIR
                if(currenti < IBEG){
                    CellTracerNode* temp = tracersToLeft;
                    tracersToLeft = tracers;
                    tracers->i = currenti - IBEG + 1;
                    tracers = tracers->next;
                    if(tracers != NULL){
                        tracers->prev = NULL;
                    }
                    tracersToLeft->next = temp;
                    tracersToLeft->prev = NULL;
                    if(temp != NULL){
                        temp->prev = tracersToLeft;
                    }
                    NtoLeft = NtoLeft + 1;
                    stopped = false;
                    if(s->lrank[0] == 0){
                        if((s->isperiodic[0] == AL_TRUE)){
                            double L = (grid->xend_glob[0] - grid->xbeg_glob[0]);
                            tracersToLeft->x1 = tracersToLeft->x1 + L;
                        } else {
                            stopped = true;
                            NtoLeft = NtoLeft - 1;
                        }
                    }
                    break;
                }
                if(currenti > IEND){
                    CellTracerNode* temp = tracersToRight;
                    tracersToRight = tracers;
                    tracers->i = currenti - IEND - 1;
                    tracers = tracers->next;
                    if(tracers != NULL){
                        tracers->prev = NULL;
                    }
                    tracersToRight->next = temp;
                    tracersToRight->prev = NULL;
                    if(temp != NULL){
                        temp->prev = tracersToRight;
                    }
                    NtoRight = NtoRight + 1;
                    stopped = false;
                    if(s->lrank[0] == s->lsize[0] - 1){
                        if((s->isperiodic[0] == AL_TRUE)){
                            double L = (grid->xend_glob[0] - grid->xbeg_glob[0]);
                            tracersToRight->x1 = tracersToRight->x1 - L;
                        } else {
                            stopped = true;
                            NtoRight = NtoRight - 1;
                        }
                    }
                    break;
                }
#endif
#if INCLUDE_JDIR
                if(currentj < JBEG){
                    CellTracerNode* temp = tracersToDown;
                    tracersToDown = tracers;
                    tracers->j = currentj - JBEG + 1;
                    tracers = tracers->next;
                    if(tracers != NULL){
                        tracers->prev = NULL;
                    }
                    tracersToDown->next = temp;
                    tracersToDown->prev = NULL;
                    if(temp != NULL){
                        temp->prev = tracersToDown;
                    }
                    NtoDown = NtoDown + 1;
                    stopped = false;
                    if(s->lrank[1] == 0){
                        if((s->isperiodic[1] == AL_TRUE)){
                            double L = (grid->xend_glob[1] - grid->xbeg_glob[1]);
                            tracersToDown->x2 = tracersToDown->x2 + L;
                        } else {
                            stopped = true;
                            NtoDown = NtoDown - 1;
                        }
                    }
                    break;
                }
                if(currentj > JEND){
                    CellTracerNode* temp = tracersToUp;
                    tracersToUp = tracers;
                    tracers->j = currentj - JEND - 1;
                    tracers = tracers->next;
                    if(tracers != NULL){
                        tracers->prev = NULL;
                    }
                    tracersToUp->next = temp;
                    tracersToUp->prev = NULL;
                    if(temp != NULL){
                        temp->prev = tracersToUp;
                    }
                    NtoUp = NtoUp + 1;
                    stopped = false;
                    if(s->lrank[1] == s->lsize[1] - 1){
                        if((s->isperiodic[1] == AL_TRUE)){
                            double L = (grid->xend_glob[1] - grid->xbeg_glob[1]);
                            tracersToUp->x2 = tracersToUp->x2 - L;
                        } else {
                            stopped = true;
                            NtoUp = NtoUp - 1;
                        }
                    }
                    break;
                }
#endif
#if INCLUDE_KDIR
                if(upstreamk < KBEG){
                    CellTracerNode* temp = tracersToBack;
                    tracersToBack = tracers;
                    tracers->k = currentk - KBEG + 1;
                    tracers = tracers->next;
                    if(tracers != NULL){
                        tracers->prev = NULL;
                    }
                    tracersToBack->next = temp;
                    tracersToBack->prev = NULL;
                    if(temp != NULL){
                        temp->prev = tracersToBack;
                    }
                    NtoBack = NtoBack + 1;
                    stopped = false;
                    if(s->lrank[2] == 0){
                        if((s->isperiodic[2] == AL_TRUE)){
                            double L = (grid->xend_glob[2] - grid->xbeg_glob[2]);
                            tracersToBack->x3 = tracersToBack->x3 + L;
                        } else {
                            stopped = true;
                            NtoBack = NtoBack - 1;
                        }
                    }
                    break;
                }
                if(upstreamk > KEND){
                    CellTracerNode* temp = tracersToFront;
                    tracersToFront = tracers;
                    tracers->k = currentk - KEND - 1;
                    tracers = tracers->next;
                    if(tracers != NULL){
                        tracers->prev = NULL;
                    }
                    tracersToFront->next = temp;
                    tracersToFront->prev = NULL;
                    if(temp != NULL){
                        temp->prev = tracersToFront;
                    }
                    NtoFront = NtoFront + 1;
                    stopped = false;
                    if(s->lrank[2] == s->lsize[2] - 1){
                        if((s->isperiodic[2] == AL_TRUE)){
                            double L = (grid->xend_glob[2] - grid->xbeg_glob[2]);
                            tracersToFront->x3 = tracersToFront->x3 - L;
                        } else {
                            stopped = true;
                            NtoFront = NtoFront - 1;
                        }
                    }
                    break;
                }
#endif
                double pgradx = 0;
                double pgrady = 0;
                double pgradz = 0;

#if INCLUDE_IDIR
                pgradx = grid->dx_dl[IDIR][currentj][currenti]*(d->Vc[PRS][currentk][currentj][currenti+1] - d->Vc[PRS][currentk][currentj][currenti-1])/(grid->x[0][currenti+1] - grid->x[0][currenti-1]);
#endif
#if INCLUDE_JDIR
                pgrady = grid->dx_dl[JDIR][currentj][currenti]*(d->Vc[PRS][currentk][currentj+1][currenti] - d->Vc[PRS][currentk][currentj-1][currenti])/(grid->x[1][currentj+1] - grid->x[1][currentj-1]);
#endif
#if INCLUDE_KDIR
                pgradz = grid->dx_dl[KDIR][currentj][currenti]*(d->Vc[PRS][currentk+1][currentj][currenti] - d->Vc[PRS][currentk-1][currentj][currenti])/(grid->x[2][currentk+1] - grid->x[2][currentk-1]);
#endif
                double gradnorm = sqrt(pgradx*pgradx + pgrady*pgrady + pgradz*pgradz);
                double scalarmult = (tracers->prevgradx*pgradx + tracers->prevgrady*pgrady + tracers->prevgradz*pgradz);
                if(scalarmult < 0){
                    break;
                }
                tracers->prevgradx = pgradx;
                tracers->prevgrady = pgrady;
                tracers->prevgradz = pgradz;

                tracers->v1 = d->Vc[VX1][currentk][currentj][currenti];
                tracers->v2 = d->Vc[VX2][currentk][currentj][currenti];
                tracers->v3 = d->Vc[VX3][currentk][currentj][currenti];
                tracers->rho = d->Vc[RHO][currentk][currentj][currenti];



                traceNextCell(grid, &x, &y, &z, direction*pgradx/gradnorm, direction*pgrady/gradnorm, direction*pgradz/gradnorm, &currenti, &currentj, &currentk);

                tracers->i = currenti;
                tracers->j = currentj;
                tracers->k = currentk;
                tracers->x1 = x;
                tracers->x2 = y;
                tracers->x3 = z;
                /*tracers->v1 = d->Vc[VX1][currentk][currentj][currenti];
                tracers->v2 = d->Vc[VX2][currentk][currentj][currenti];
                tracers->v3 = d->Vc[VX3][currentk][currentj][currenti];
                tracers->rho = d->Vc[RHO][currentk][currentj][currenti];*/
            }
            if(stopped){
                CellTracerNode* temp = stoppedTracers;
                stoppedTracers = tracers;
                tracers = tracers->next;
                if(tracers != NULL){
                    tracers->prev = NULL;
                }
                stoppedTracers->next = temp;
                stoppedTracers->prev = NULL;
            }
        }
        NactiveTracers[0] = 0;

        int Nout[1];
        int Nin[1];
        int* outbuf;
        double* outbufd;
        int* inbuf;
        double* inbufd;

        int intBufCount = 8;
        int doubleBufCount = 11;

#if INCLUDE_IDIR
        Nout[0] = NtoLeft;
        Nin[0] = 0;

        nleft = s->left[0];
        nright = s->right[0];
        tag1 = s->tag1[0];

        MPI_Sendrecv(Nout, 1, MPI_INT, nleft, tag1,
                     Nin, 1, MPI_INT, nright, tag1,
                     comm, &status);

        NactiveTracers[0] = NactiveTracers[0] + Nin[0];

        if(Nout[0] != 0){
            outbuf = (int*) malloc(intBufCount*Nout[0]*sizeof(int));
            outbufd = (double*) malloc(doubleBufCount*Nout[0]*sizeof(double));
            putTracerListToArray(tracersToLeft, outbuf, outbufd);
            tracersToLeft = NULL;
            NtoLeft = 0;
        }

        if(Nin[0] != 0){
            inbuf = (int*) malloc(intBufCount*Nin[0]*sizeof(int));
            inbufd = (double*) malloc(doubleBufCount*Nin[0]*sizeof(double));
        }

        MPI_Sendrecv(outbuf, intBufCount*Nout[0], MPI_INT, nleft, tag1,
                     inbuf, intBufCount*Nin[0], MPI_INT, nright, tag1,
                     comm, &status);
        MPI_Sendrecv(outbufd, doubleBufCount*Nout[0], MPI_DOUBLE, nleft, tag1,
                     inbufd, doubleBufCount*Nin[0], MPI_DOUBLE, nright, tag1,
                     comm, &status);

        if(Nin[0] != 0){
            CellTracerNode* tracersFrom = putArrayToTracerList(inbuf, inbufd, Nin[0]);
            while(tracersFrom != NULL){
                tracersFrom->i = tracersFrom->i + IEND;
                //tracersFrom->x1 = grid->xr[0][IEND];
                CellTracerNode* temp = tracersFrom;
                tracersFrom = tracersFrom->next;
                if(tracersFrom != NULL){
                    tracersFrom->prev = NULL;
                }
                temp->next = tracers;
                if(tracers != NULL){
                    tracers->prev = temp;
                }
                temp->prev = NULL;
                tracers = temp;
            }
        }

        if(Nout[0] > 0){
            free(outbuf);
            free(outbufd);
        }
        if(Nin[0] > 0){
            free(inbuf);
            free(inbufd);
        }

        Nout[0] = NtoRight;
        Nin[0] = 0;

        MPI_Sendrecv(Nout, 1, MPI_INT, nright, tag1,
                     Nin, 1, MPI_INT, nleft, tag1,
                     comm, &status);

        NactiveTracers[0] = NactiveTracers[0] + Nin[0];

        if(Nout[0] != 0){
            outbuf = (int*) malloc(intBufCount*Nout[0]*sizeof(int));
            outbufd = (double*) malloc(doubleBufCount*Nout[0]*sizeof(double));
            putTracerListToArray(tracersToRight, outbuf, outbufd);
            tracersToRight = NULL;
            NtoRight = 0;
        }

        if(Nin[0] != 0){
            inbuf = (int*) malloc(intBufCount*Nin[0]*sizeof(int));
            inbufd = (double*) malloc(doubleBufCount*Nin[0]*sizeof(double));
        }

        MPI_Sendrecv(outbuf, intBufCount*Nout[0], MPI_INT, nright, tag1,
                     inbuf, intBufCount*Nin[0], MPI_INT, nleft, tag1,
                     comm, &status);
        MPI_Sendrecv(outbufd, doubleBufCount*Nout[0], MPI_DOUBLE, nright, tag1,
                     inbufd, doubleBufCount*Nin[0], MPI_DOUBLE, nleft, tag1,
                     comm, &status);

        if(Nin[0] != 0){
            CellTracerNode* tracersFrom = putArrayToTracerList(inbuf, inbufd, Nin[0]);
            while(tracersFrom != NULL){
                tracersFrom->i = tracersFrom->i + IBEG;
                //tracersFrom->x1 = grid->xl[0][IBEG];
                CellTracerNode* temp = tracersFrom;
                tracersFrom = tracersFrom->next;
                if(tracersFrom != NULL){
                    tracersFrom->prev = NULL;
                }
                temp->next = tracers;
                if(tracers != NULL){
                    tracers->prev = temp;
                }
                temp->prev = NULL;
                tracers = temp;
            }
        }

        if(Nout[0] > 0){
            free(outbuf);
            free(outbufd);
        }
        if(Nin[0] > 0){
            free(inbuf);
            free(inbufd);
        }
#endif
#if INCLUDE_JDIR
        Nout[0] = NtoDown;
        Nin[0] = 0;

        nleft = s->left[1];
        nright = s->right[1];
        tag1 = s->tag1[1];

        MPI_Sendrecv(Nout, 1, MPI_INT, nleft, tag1,
                     Nin, 1, MPI_INT, nright, tag1,
                     comm, &status);

        NactiveTracers[0] = NactiveTracers[0] + Nin[0];

        if(Nout[0] != 0){
            outbuf = (int*) malloc(intBufCount*Nout[0]*sizeof(int));
            outbufd = (double*) malloc(doubleBufCount*Nout[0]*sizeof(double));
            putTracerListToArray(tracersToDown, outbuf, outbufd);
            tracersToDown = NULL;
            NtoDown = 0;
        }

        if(Nin[0] != 0){
            inbuf = (int*) malloc(intBufCount*Nin[0]*sizeof(int));
            inbufd = (double*) malloc(doubleBufCount*Nin[0]*sizeof(double));
        }

        MPI_Sendrecv(outbuf, intBufCount*Nout[0], MPI_INT, nleft, tag1,
                     inbuf, intBufCount*Nin[0], MPI_INT, nright, tag1,
                     comm, &status);
        MPI_Sendrecv(outbufd, doubleBufCount*Nout[0], MPI_DOUBLE, nleft, tag1,
                     inbufd, doubleBufCount*Nin[0], MPI_DOUBLE, nright, tag1,
                     comm, &status);

        if(Nin[0] != 0){
            CellTracerNode* tracersFrom = putArrayToTracerList(inbuf, inbufd, Nin[0]);
            while(tracersFrom != NULL){
                tracersFrom->j = tracersFrom->j + JEND;
                //tracersFrom->x2 = grid->xr[1][JEND];
                CellTracerNode* temp = tracersFrom;
                tracersFrom = tracersFrom->next;
                if(tracersFrom != NULL){
                    tracersFrom->prev = NULL;
                }
                temp->next = tracers;
                if(tracers != NULL){
                    tracers->prev = temp;
                }
                temp->prev = NULL;
                tracers = temp;
            }
        }

        if(Nout[0] > 0){
            free(outbuf);
            free(outbufd);
        }
        if(Nin[0] > 0){
            free(inbuf);
            free(inbufd);
        }

        Nout[0] = NtoUp;
        Nin[0] = 0;

        MPI_Sendrecv(Nout, 1, MPI_INT, nright, tag1,
                     Nin, 1, MPI_INT, nleft, tag1,
                     comm, &status);

        NactiveTracers[0] = NactiveTracers[0] + Nin[0];

        if(Nout[0] != 0){
            outbuf = (int*) malloc(intBufCount*Nout[0]*sizeof(int));
            outbufd = (double*) malloc(doubleBufCount*Nout[0]*sizeof(double));
            putTracerListToArray(tracersToUp, outbuf, outbufd);
            tracersToUp = NULL;
            NtoUp = 0;
        }

        if(Nin[0] != 0){
            inbuf = (int*) malloc(intBufCount*Nin[0]*sizeof(int));
            inbufd = (double*) malloc(doubleBufCount*Nin[0]*sizeof(double));
        }

        MPI_Sendrecv(outbuf, intBufCount*Nout[0], MPI_INT, nright, tag1,
                     inbuf, intBufCount*Nin[0], MPI_INT, nleft, tag1,
                     comm, &status);
        MPI_Sendrecv(outbufd, doubleBufCount*Nout[0], MPI_DOUBLE, nright, tag1,
                     inbufd, doubleBufCount*Nin[0], MPI_DOUBLE, nleft, tag1,
                     comm, &status);

        if(Nin[0] != 0){
            CellTracerNode* tracersFrom = putArrayToTracerList(inbuf, inbufd, Nin[0]);
            while(tracersFrom != NULL){
                tracersFrom->j = tracersFrom->j + JBEG;
                //tracersFrom->x2 = grid->xl[1][JBEG];
                CellTracerNode* temp = tracersFrom;
                tracersFrom = tracersFrom->next;
                if(tracersFrom != NULL){
                    tracersFrom->prev = NULL;
                }
                temp->next = tracers;
                if(tracers != NULL){
                    tracers->prev = temp;
                }
                temp->prev = NULL;
                tracers = temp;
            }
        }

        if(Nout[0] > 0){
            free(outbuf);
            free(outbufd);
        }
        if(Nin[0] > 0){
            free(inbuf);
            free(inbufd);
        }
#endif
#if INCLUDE_KDIR
        Nout[0] = NtoBack;
        Nin[0] = 0;

        nleft = s->left[2];
        nright = s->right[2];
        tag1 = s->tag1[2];

        MPI_Sendrecv(Nout, 1, MPI_INT, nleft, tag1,
                     Nin, 1, MPI_INT, nright, tag1,
                     comm, &status);

        NactiveTracers[0] = NactiveTracers[0] + Nin[0];

        if(Nout[0] != 0){
            outbuf = (int*) malloc(intBufCount*Nout[0]*sizeof(int));
            outbufd = (double*) malloc(doubleBufCount*Nout[0]*sizeof(double));
            putTracerListToArray(tracersToBack, outbuf, outbufd);
            tracersToBack = NULL;
            NtoBack = 0;
        }

        if(Nin[0] != 0){
            inbuf = (int*) malloc(intBufCount*Nin[0]*sizeof(int));
            inbufd = (double*) malloc(doubleBufCount*Nin[0]*sizeof(double));
        }

        MPI_Sendrecv(outbuf, intBufCount*Nout[0], MPI_INT, nleft, tag1,
                     inbuf, intBufCount*Nin[0], MPI_INT, nright, tag1,
                     comm, &status);
        MPI_Sendrecv(outbufd, doubleBufCount*Nout[0], MPI_DOUBLE, nleft, tag1,
                     inbufd, doubleBufCount*Nin[0], MPI_DOUBLE, nright, tag1,
                     comm, &status);

        if(Nin[0] != 0){
            CellTracerNode* tracersFrom = putArrayToTracerList(inbuf, inbufd, Nin[0]);
            while(tracersFrom != NULL){
                tracersFrom->k = tracersFrom->k + KEND;
                //tracersFrom->x3 = grid->xr[2][KEND];
                CellTracerNode* temp = tracersFrom;
                tracersFrom = tracersFrom->next;
                if(tracersFrom != NULL){
                    tracersFrom->prev = NULL;
                }
                temp->next = tracers;
                if(tracers != NULL){
                    tracers->prev = temp;
                }
                temp->prev = NULL;
                tracers = temp;
            }
        }

        if(Nout[0] > 0){
            free(outbuf);
            free(outbufd);
        }
        if(Nin[0] > 0){
            free(inbuf);
            free(inbufd);
        }

        Nout[0] = NtoFront;
        Nin[0] = 0;

        MPI_Sendrecv(Nout, 1, MPI_INT, nright, tag1,
                     Nin, 1, MPI_INT, nleft, tag1,
                     comm, &status);

        NactiveTracers[0] = NactiveTracers[0] + Nin[0];

        if(Nout[0] != 0){
            outbuf = (int*) malloc(intBufCount*Nout[0]*sizeof(int));
            outbufd = (double*) malloc(doubleBufCount*Nout[0]*sizeof(double));
            putTracerListToArray(tracersToFront, outbuf, outbufd);
            tracersToFront = NULL;
            NtoFront = 0;
        }

        if(Nin[0] != 0){
            inbuf = (int*) malloc(intBufCount*Nin[0]*sizeof(int));
            inbufd = (double*) malloc(doubleBufCount*Nin[0]*sizeof(double));
        }

        MPI_Sendrecv(outbuf, intBufCount*Nout[0], MPI_INT, nright tag1,
                     inbuf, intBufCount*Nin[0], MPI_INT, nleft, tag1,
                     comm, &status);
        MPI_Sendrecv(outbufd, doubleBufCount*Nout[0], MPI_DOUBLE, nright, tag1,
                     inbufd, doubleBufCount*Nin[0], MPI_DOUBLE, nleft, tag1,
                     comm, &status);

        if(Nin[0] != 0){
            CellTracerNode* tracersFrom = putArrayToTracerList(inbuf, inbufd, Nin[0]);
            while(tracersFrom != NULL){
                tracersFrom->k = tracersFrom->k + KBEG;
                //tracersFrom->x3 = grid->xl[2][KBEG];
                CellTracerNode* temp = tracersFrom;
                tracersFrom = tracersFrom->next;
                if(tracersFrom != NULL){
                    tracersFrom->prev = NULL;
                }
                temp->next = tracers;
                if(tracers != NULL){
                    tracers->prev = temp;
                }
                temp->prev = NULL;
                tracers = temp;
            }
        }

        if(Nout[0] > 0){
            free(outbuf);
            free(outbufd);
        }
        if(Nin[0] > 0){
            free(inbuf);
            free(inbufd);
        }
#endif

        a[0] = NactiveTracers[0];
        MPI_Allreduce(a, NactiveTracers, 1, MPI_INT, MPI_SUM, comm);
    }

    nprocs = s->size;

    int* sendcounts = (int*) malloc(nprocs*sizeof(int));
    int* sendcountsd = (int*) malloc(nprocs*sizeof(int));
    int* sdispls = (int*) malloc(nprocs*sizeof(int));
    int* sdisplsd = (int*) malloc(nprocs*sizeof(int));
    int* srelposition = (int*) malloc(nprocs*sizeof(int));
    int* recvcounts = (int*) malloc(nprocs*sizeof(int));
    int* recvcountsd = (int*) malloc(nprocs*sizeof(int));
    int* rdispls = (int*) malloc(nprocs*sizeof(int));
    int* rdisplsd = (int*) malloc(nprocs*sizeof(int));

    int intDataN = 3;
    int doubleDataN = 7;

    CellTracerNode* tempTracer = stoppedTracers;

    for(int i = 0; i < nprocs; ++i){
        sendcounts[i] = 0;
        sendcountsd[i] = 0;
        sdispls[i] = 0;
        sdisplsd[i] = 0;
        srelposition[i] = 0;
        recvcounts[i] = 0;
        recvcountsd[i] = 0;
        rdispls[i] = 0;
        rdisplsd[i] = 0;
    }

    while(tempTracer != NULL){
        int trank = tempTracer->rank0;
        sendcounts[trank] = sendcounts[trank] + intDataN;
        sendcountsd[trank] = sendcountsd[trank] + doubleDataN;
        tempTracer = tempTracer->next;
    }

    MPI_Alltoall(sendcounts, 1, MPI_INT, recvcounts, 1, MPI_INT, comm);
    MPI_Alltoall(sendcountsd, 1, MPI_INT, recvcountsd, 1, MPI_INT, comm);

    for(int i = 1; i < nprocs; ++i){
        sdispls[i] = sdispls[i-1] + sendcounts[i-1];
        sdisplsd[i] = sdisplsd[i-1] + sendcountsd[i-1];
        rdispls[i] = rdispls[i-1] + recvcounts[i-1];
        rdisplsd[i] = rdisplsd[i-1] + recvcountsd[i-1];
    }

    int totalSend = 0;
    int totalRecv = 0;
    for(int i = 0; i < nprocs; ++i){
        //todo division
        totalSend = totalSend + sendcounts[i]/intDataN;
        totalRecv = totalRecv + recvcounts[i]/intDataN;
    }

    int* outbuf;
    double* outbufd;
    int* inbuf;
    double* inbufd;

    if(totalSend > 0){
        outbuf = (int*) malloc(intDataN*totalSend*sizeof(int));
        outbufd = (double*) malloc(doubleDataN*totalSend*sizeof(double));
    }
    if(totalRecv > 0){
        inbuf = (int*) malloc(intDataN*totalRecv*sizeof(int));
        inbufd = (double*) malloc(doubleDataN*totalRecv*sizeof(double));
    }

    tempTracer = stoppedTracers;

    while(tempTracer != NULL){
        int rankt = tempTracer->rank0;

        outbuf[sdispls[rankt] + intDataN*srelposition[rankt]] = tempTracer->i0;
        outbuf[sdispls[rankt] + intDataN*srelposition[rankt]+1] = tempTracer->j0;
        outbuf[sdispls[rankt] + intDataN*srelposition[rankt]+2] = tempTracer->k0;

        outbufd[sdisplsd[rankt] + doubleDataN*srelposition[rankt]] = tempTracer->x1;
        outbufd[sdisplsd[rankt] + doubleDataN*srelposition[rankt]+1] = tempTracer->x2;
        outbufd[sdisplsd[rankt] + doubleDataN*srelposition[rankt]+2] = tempTracer->x3;
        outbufd[sdisplsd[rankt] + doubleDataN*srelposition[rankt]+3] = tempTracer->v1;
        outbufd[sdisplsd[rankt] + doubleDataN*srelposition[rankt]+4] = tempTracer->v2;
        outbufd[sdisplsd[rankt] + doubleDataN*srelposition[rankt]+5] = tempTracer->v3;
        outbufd[sdisplsd[rankt] + doubleDataN*srelposition[rankt]+6] = tempTracer->rho;

        srelposition[rankt] = srelposition[rankt] + 1;

        tempTracer = tempTracer->next;
    }

    MPI_Alltoallv(outbuf, sendcounts, sdispls, MPI_INT, inbuf, recvcounts, rdispls, MPI_INT, comm);
    MPI_Alltoallv(outbufd, sendcountsd, sdisplsd, MPI_DOUBLE, inbufd, recvcountsd, rdisplsd, MPI_DOUBLE, comm);

    int count = 0;
    int countd = 0;
    if(totalRecv > 0){
        for(int l = 0; l < totalRecv; ++l){
            int i = inbuf[count];
            count++;
            int j = inbuf[count];
            count++;
            int k = inbuf[count];
            count++;

            x1[k][j][i] = inbufd[countd];
            countd++;
            x2[k][j][i] = inbufd[countd];
            countd++;
            x3[k][j][i] = inbufd[countd];
            countd++;
            v1[k][j][i] = inbufd[countd];
            countd++;
            v2[k][j][i]= inbufd[countd];
            countd++;
            v3[k][j][i] = inbufd[countd];
            countd++;
            rho[k][j][i] = inbufd[countd];
            countd++;
        }
    }

    while(stoppedTracers != NULL){
        CellTracerNode* tempTracer = stoppedTracers->next;
        free(stoppedTracers);
        stoppedTracers = tempTracer;
    }

    if(totalSend > 0){
        free(outbuf);
        free(outbufd);
    }
    if(totalRecv > 0){
        free(inbuf);
        free(inbufd);
    }


    free(sendcounts);
    free(sendcountsd);
    free(sdispls);
    free(sdisplsd);
    free(srelposition);
    free(recvcounts);
    free(recvcountsd);
    free(rdispls);
    free(rdisplsd);
}

#endif

void traceShock(Data* d, Grid* grid, int direction, double*** x1, double*** x2, double*** x3, double*** v1, double*** v2, double*** v3, double*** rho, double*** pressure){
    if((direction != 1) && (direction != -1)){
        printf("direction for tracing must be 1 or -1\n");
        printLog("direction for tracing must be 1 or -1\n");
        QUIT_PLUTO(0);
    }

    int i, j, k;

#ifdef PARALLEL
    traceShockParallel(d, grid, direction, x1, x2, x3, v1, v2, v3, rho, pressure);
#else
    DOM_LOOP(k,j,i){
        if(!(d->flag[k][j][i] & FLAG_ENTROPY)){

            int currenti = i;
            int currentj = j;
            int currentk = k;

            double x = grid->x[0][i];
            double y = grid->x[1][j];
            double z = grid->x[2][k];

            double prevgradx = 0;
            double prevgrady = 0;
            double prevgradz = 0;

            int count = 0;

            while(!(d->flag[currentk][currentj][currenti] & FLAG_ENTROPY)){
#if INCLUDE_IDIR
                if(currenti < IBEG){
                    if(grid->lbound[0] == PERIODIC){
                        double L = (grid->xend_glob[0] - grid->xbeg_glob[0]);
                        x = x + L;
                        currenti = IEND + currenti - IBEG + 1;
                    } else {
                        break;
                    }
                }
                if(currenti > IEND){
                    if(grid->rbound[0] == PERIODIC){
                        double L = (grid->xend_glob[0] - grid->xbeg_glob[0]);
                        x = x - L;
                        currenti = IBEG + currenti - IEND - 1;
                    } else {
                        break;
                    }
                }
#endif
#if INCLUDE_JDIR
                if(currentj < JBEG){
                    if(grid->lbound[1] == PERIODIC){
                        double L = (grid->xend_glob[1] - grid->xbeg_glob[1]);
                        y = y + L;
                        currentj = JEND + currentj - JBEG + 1;
                    } else {
                        break;
                    }
                }
                if(currentj > JEND){
                    if(grid->rbound[1] == PERIODIC){
                        double L = (grid->xend_glob[1] - grid->xbeg_glob[1]);
                        y = y - L;
                        currentj = JBEG + currentj - JEND - 1;
                    } else {
                        break;
                    }
                }
#endif
#if INCLUDE_KDIR
                if(upstreamk < KBEG){
                    if(grid->lbound[2] == PERIODIC){
                        double L = (grid->xend_glob[2] - grid->xbeg_glob[2]);
                        z = z + L;
                        currentk = KEND + currentk - KBEG + 1;
                    } else {
                        break;
                    }
                }
                if(upstreamk > KEND){
                    if(grid->rbound[2] == PERIODIC){
                        double L = (grid->xend_glob[2] - grid->xbeg_glob[2]);
                        z = z - L;
                        currentk = KBEG + currentk - KEND - 1;
                    } else {
                        break;
                    }
                }
#endif
                count++;
                //printf("n traced cell = %d\n", count);
                if(count > 1000){
                    printf("n traced cells > 1000\n");
                    printLog("n traced cells > 1000\n");
                    QUIT_PLUTO(0);
                }
                double pgradx = 0;
                double pgrady = 0;
                double pgradz = 0;

#if INCLUDE_IDIR
                pgradx = grid->dx_dl[IDIR][currentj][currenti]*(d->Vc[PRS][currentk][currentj][currenti+1] - d->Vc[PRS][currentk][currentj][currenti-1])/(grid->x[0][currenti+1] - grid->x[0][currenti-1]);
#endif
#if INCLUDE_JDIR
                pgrady = grid->dx_dl[JDIR][currentj][currenti]*(d->Vc[PRS][currentk][currentj+1][currenti] - d->Vc[PRS][currentk][currentj-1][currenti])/(grid->x[1][currentj+1] - grid->x[1][currentj-1]);
#endif
#if INCLUDE_KDIR
                pgradz = grid->dx_dl[KDIR][upstreamj][upstreami]*(d->Vc[PRS][upstreamk+1][upstreamj][upstreami] - d->Vc[PRS][upstreamk-1][upstreamj][upstreami])/(grid->x[2][upstreamk+1] - grid->x[2][upstreamk-1]);
#endif
                double gradnorm = sqrt(pgradx*pgradx+ pgrady*pgrady + pgradz*pgradz);
                double scalarmult = (prevgradx*pgradx + prevgrady*pgrady + prevgradz*pgradz);
                if(scalarmult < 0){
                    break;
                }
                prevgradx = pgradx;
                prevgrady = pgrady;
                prevgradz = pgradz;
                v1[k][j][i] = d->Vc[VX1][currentk][currentj][currenti];
                v2[k][j][i] = d->Vc[VX2][currentk][currentj][currenti];
                v3[k][j][i] = d->Vc[VX3][currentk][currentj][currenti];

                rho[k][j][i] = d->Vc[RHO][currentk][currentj][currenti];
                pressure[k][j][i] = d->Vc[PRS][currentk][currentj][currenti];

                traceNextCell(grid, &x, &y, &z, direction*pgradx/gradnorm, direction*pgrady/gradnorm, direction*pgradz/gradnorm, &currenti, &currentj, &currentk);
            }

            x1[k][j][i] = grid->x[0][currenti];
            x2[k][j][i] = grid->x[1][currentj];
            x3[k][j][i] = grid->x[2][currentk];

            /*v1[k][j][i] = d->Vc[VX1][currentk][currentj][currenti];
            v2[k][j][i] = d->Vc[VX2][currentk][currentj][currenti];
            v3[k][j][i] = d->Vc[VX3][currentk][currentj][currenti];

            rho[k][j][i] = d->Vc[RHO][currentk][currentj][currenti];*/
        }
    }
#endif
}

void updateShockFront(Data* d, Grid* grid){
    int i,j,k;
    //printf("evaluating shock\n");
    FlagShock(d, grid);
    TOT_LOOP(k,j,i){
        d->shockWidth[k][j][i] = grid->dx[0][i];
        d->velocityJump[k][j][i] = 0.0;
        d->upstreamDensity[k][j][i] = d->Vc[RHO][k][j][i];
        d->downstreamDensity[k][j][i] = d->Vc[RHO][k][j][i];
        d->upstreamPressure[k][j][i] = d->Vc[PRS][k][j][i];
        d->downstreamPressure[k][j][i] = d->Vc[PRS][k][j][i];
        d->downstreamV1[k][j][i] = 0.0;
        d->downstreamV2[k][j][i] = 0.0;
        d->downstreamV3[k][j][i] = 0.0;
        d->upstreamV1[k][j][i] = 0.0;
        d->upstreamV2[k][j][i] = 0.0;
        d->upstreamV3[k][j][i] = 0.0;
    }


    traceShock(d, grid, -1, d->upstreamx1, d->upstreamx2, d->upstreamx3, d->upstreamV1, d->upstreamV2, d->upstreamV3, d->upstreamDensity, d->upstreamPressure);
    traceShock(d, grid, 1, d->downstreamx1, d->downstreamx2, d->downstreamx3, d->downstreamV1, d->downstreamV2, d->downstreamV3, d->downstreamDensity, d->downstreamPressure);

    DOM_LOOP(k,j,i){
        if(!(d->flag[k][j][i] & FLAG_ENTROPY)){
            double xd, yd, zd;
            double xu, yu, zu;
            double Vux, Vuy, Vuz;
            double Vdx, Vdy, Vdz;

#if GEOMETRY == CARTESIAN
            xd = d->downstreamx1[k][j][i];
            yd = d->downstreamx2[k][j][i];
            zd = d->downstreamx3[k][j][i];

            xu = d->upstreamx1[k][j][i];
            yu = d->upstreamx2[k][j][i];
            zu = d->upstreamx3[k][j][i];

            Vdx = d->downstreamV1[k][j][i];
            Vdy = d->downstreamV2[k][j][i];
            Vdz = d->downstreamV3[k][j][i];

            Vux = d->upstreamV1[k][j][i];
            Vuy = d->upstreamV2[k][j][i];
            Vuz = d->upstreamV3[k][j][i];
            //todo other velocity
#elif GEOMETRY == CYLINDRICAL
            xd = d->downstreamx1[k][j][i]*cos(d->downstreamx3[k][j][i]);
            yd = d->downstreamx1[k][j][i]*sin(d->downstreamx3[k][j][i]);
            zd = d->downstreamx2[k][j][i];

            xu = d->upstreamx1[k][j][i]*cos(d->upstreamx3[k][j][i]);
            yu = d->upstreamx1[k][j][i]*sin(d->upstreamx3[k][j][i]);
            zu = d->upstreamx2[k][j][i];

            Vdx = d->downstreamV1[k][j][i]*cos(d->downstreamx3[k][j][i]) - d->downstreamV3[k][j][i]*sin(d->downstreamx3[k][j][i]);
            Vdy = d->downstreamV1[k][j][i]*sin(d->downstreamx3[k][j][i]) + d->downstreamV3[k][j][i]*cos(d->downstreamx3[k][j][i]);
            Vdz = d->downstreamV2[k][j][i];

            Vux = d->upstreamV1[k][j][i]*cos(d->upstreamx3[k][j][i]) - d->upstreamV3[k][j][i]*sin(d->upstreamx3[k][j][i]);
            Vuy = d->upstreamV1[k][j][i]*sin(d->upstreamx3[k][j][i]) + d->upstreamV3[k][j][i]*cos(d->upstreamx3[k][j][i]);
            Vuz = d->upstreamV2[k][j][i];
#elif GEOMETRY == POLAR
            xd = d->downstreamx1[k][j][i]*cos(d->downstreamx2[k][j][i]);
            yd = d->downstreamx1[k][j][i]*sin(d->downstreamx2[k][j][i]);
            zd = d->downstreamx3[k][j][i];

            xu = d->upstreamx1[k][j][i]*cos(d->upstreamx2[k][j][i]);
            yu = d->upstreamx1[k][j][i]*sin(d->upstreamx2[k][j][i]);
            zu = d->upstreamx3[k][j][i];

            Vdx = d->downstreamV1[k][j][i]*cos(d->downstreamx2[k][j][i]) - d->downstreamV2[k][j][i]*sin(d->downstreamx2[k][j][i]);
            Vdy = d->downstreamV1[k][j][i]*sin(d->downstreamx2[k][j][i]) + d->downstreamV2[k][j][i]*cos(d->downstreamx2[k][j][i]);
            Vdz = d->downstreamV3[k][j][i];

            Vux = d->upstreamV1[k][j][i]*cos(d->upstreamx2[k][j][i]) - d->upstreamV2[k][j][i]*sin(d->upstreamx2[k][j][i]);
            Vuy = d->upstreamV1[k][j][i]*sin(d->upstreamx2[k][j][i]) + d->upstreamV2[k][j][i]*cos(d->upstreamx2[k][j][i]);
            Vuz = d->upstreamV3[k][j][i];
#elif GEOMETRY == SPHERICAL
            xd = d->downstreamx1[k][j][i]*sin(d->downstreamx2[k][j][i])*cos(d->downstreamx3[k][j][i]);
            yd = d->downstreamx1[k][j][i]*sin(d->downstreamx2[k][j][i])*sin(d->downstreamx3[k][j][i]);
            zd = d->downstreamx1[k][j][i]*cos(d->downstreamx2[k][j][i]);

            xu = d->upstreamx1[k][j][i]*sin(grid->x[1][upstreamj])*cos(d->upstreamx3[k][j][i]);
            yu = d->upstreamx1[k][j][i]*sin(d->upstreamx2[k][j][i])*sin(d->upstreamx3[k][j][i]);
            zu = d->upstreamx1[k][j][i]*cos(d->upstreamx2[k][j][i]);

            Vdx = d->downstreamV1[k][j][i]*sin(d->downstreamx2[k][j][i])*cos(d->downstreamx3[k][j][i]) + d->downstreamV2[k][j][i]*cos(d->downstreamx2[k][j][i])*cos(d->downstreamx3[k][j][i]) - d->downstreamV3[k][j][i]*sin(d->downstreamx3[k][j][i]);
            Vdy = d->downstreamV1[k][j][i]*sin(d->downstreamx2[k][j][i])*sin(d->downstreamx3[k][j][i]) + d->downstreamV2[k][j][i]*cos(d->downstreamx2[k][j][i])*sin(d->downstreamx3[k][j][i]) + d->downstreamV3[k][j][i]*cos(d->downstreamx3[k][j][i]);
            Vdz = d->downstreamV1[k][j][i]*cos(d->downstreamx2[k][j][i]) - d->downstreamV2[k][j][i]*cos(d->downstreamx2[k][j][i]);

            Vux = d->upstreamV1[k][j][i]*sin(d->upstreamx2[k][j][i])*cos(d->upstreamx3[k][j][i]) + d->upstreamV2[k][j][i]*cos(d->upstreamx2[k][j][i])*cos(d->upstreamx3[k][j][i]) - d->upstreamV3[k][j][i]*sin(d->upstreamx3[k][j][i]);
            Vuy = d->upstreamV1[k][j][i]*sin(d->upstreamx2[k][j][i])*sin(d->upstreamx3[k][j][i]) + d->upstreamV2[k][j][i]*cos(d->upstreamx2[k][j][i])*sin(d->upstreamx3[k][j][i]) + d->upstreamV3[k][j][i]*cos(d->upstreamx3[k][j][i]);
            Vuz = d->upstreamV1[k][j][i]*cos(d->upstreamx2[k][j][i]) - d->upstreamV2[k][j][i]*cos(d->upstreamx2[k][j][i]);
#else
#endif
            double width = sqrt((xd-xu)*(xd-xu) + (yd-yu)*(yd-yu) + (zd-zu)*(zd-zu));
            //double V = sqrt((Vdx - Vux)*(Vdx - Vux) + (Vdy - Vuy)*(Vdy - Vuy) + (Vdz - Vuz)*(Vdz - Vuz));
            double V = fabs((Vdx - Vux)*(xd - xu) + (Vdy - Vuy)*(yd - yu) + (Vdz - Vuz)*(zd - zu))/width;

            d->shockWidth[k][j][i] = width;
            d->velocityJump[k][j][i] = V;
        }
    }
}

void traceNextCell(Grid* grid, double* x1, double* x2, double* x3, double v1, double v2, double v3, int* i, int* j, int* k){
    //todo proper stright lines for other geometries
    double vsqr = v1*v1 + v2*v2 + v3*v3;
    if(vsqr <= 0){
        printf("v = 0 in traceNextCell\n");
        printLog("v = 0 in traceNextCell\n");
        QUIT_PLUTO(1);
    }
#if INCLUDE_KDIR
    double lx = grid->xl[0][*i];
    double rx = grid->xr[0][*i];
    double ly = grid->xl[1][*j];
    double ry = grid->xr[1][*j];
    double lz = grid->xl[2][*k];
    double rz = grid->xr[2][*k];

    double dx;
    double dy;
    double dz;
    if(vx > 0){
        dx = (rx - *x1)/grid->dx_dl[IDIR][*j][*i];
        if(dx == 0){
            (*i) = (*i) + 1;
            traceNextCell(grid, x1, x2, x3, v1, v2, v3, i, j, k);
            return;
        }
    } else {
        dx = (*x1 - lx)/grid->dx_dl[IDIR][*j][*i];
        if(dx == 0){
            (*i) = (*i) - 1;
            traceNextCell(grid, x1, x2, x3, v1, v2, v3, i, j, k);
            return;
        }
    }
    if(vy > 0){
        dy = (ry - *x2)/grid->dx_dl[JDIR][*j][*i];
        if(dy == 0){
            (*j) = (*j) + 1;
            traceNextCell(grid, x1, x2, x3, v1, v2, v3, i, j, k);
            return;
        }
    } else {
        dy = (*x2 - ly)/grid->dx_dl[JDIR][*j][*i];
        if(dy == 0){
            (*j) = (*j) - 1;
            traceNextCell(grid, x1, x2, x3, v1, v2, v3, i, j, k);
            return;
        }
    }
    if(vz > 0){
        dz = (rz - *x3)/grid->dx_dl[KDIR][*j][*i];
        if(dz == 0){
            (*k) = (*k) + 1;
            traceNextCell(grid, x1, x2, x3, v1, v2, v3, i, j, k);
            return;
        }
    } else {
        dz = (*x3 - lz)/grid->dx_dl[KDIR][*j][*i];
        if(dz == 0){
            (*k) = (*k) - 1;
            traceNextCell(grid, x1, x2, x3, v1, v2, v3, i, j, k);
            return;
        }
    }
#if GEOMETRY == CARTESIAN

    if(fabs(dx*vy) > fabs(dy*vx)){
        if(fabs(dz*vy) > fabs(dy*vz)){
            double dt = fabs(dy/vy);
            *x1 = *x1 + dt*v1;
            *x3 = *x3 + dt*v3;
            if(vy > 0){
                *x2 = ry;
                *j = (*j)+1;
            } else {
                *x2 = ly;
                *j = (*j)-1;
            }
        } else {
            double dt = fabs(dz/vz);
            *x1 = *x1 + dt*v1;
            *x2 = *x2 + dt*v2;
            if(vz > 0){
                *x3 = rz;
                *k = (*k)+1;
            } else {
                *x3 = lz;
                *k = (*k)-1;
            }
        }
    } else {
        if(fabs(dz*vx) > fabs(dx*vz)){
            double dt = fabs(dx/vx);
            *x2 = *x2 + dt*v2;
            *x3 = *x3 + dt*v3;
            if(vx > 0){
                *x1 = rx;
                *i = (*i)+1;
            } else {
                *x1 = lx;
                *i = (*i)-1;
            }
        } else {
            double dt = fabs(dz/vz);
            *x1 = *x1 + dt*v1;
            *x2 = *x2 + dt*v2;
            if(vz > 0){
                *x3 = rz;
                *k = (*k)+1;
            } else {
                *x3 = lz;
                *k = (*k)-1;
            }
        }
    }
    CheckNanOrInfinity(*x1, "x1 = NaN\n");
    CheckNanOrInfinity(*x2, "x2 = NaN\n");
    CheckNanOrInfinity(*x3, "x3 = NaN\n");
#elif GEOMETRY == CYLINDRICAL
    double invdtz = fabs(v2/dy);
    double invdtr = 0;
    double invdtphi = 0;

    double vx = v1*cos(grid->x[2][k]) - v3*sin(grid->x[2][k]);
    double vy = v1*sin(grid->x[2][k]) + v3*cos(grid->x[2][k]);

    double x0 = x1*cos(x3);
    double y0 = x1*sin(x3);

    double vxy2 = vx*vx + vy*vy;


    if(v1 > 0){
        double D = vxy2*xr*xr - sqr(vx*y0 + vy*x0);
        invdtr = vxy2/(sqrt(D) - vx*x0 - vy*y0);
        if(invdtr < 0){
            printf("invdtr < 0\n");
            QUIT_PLUTO(0);
        }
    } else {
        double D = vxy2*xl*xl - sqr(vx*y0 + vy*x0);
        if(D >= 0){
            invdtr = vxy2/(-sqrt(D) - vx*x0 - vy*y0);
            if(invdtr < 0){
                printf("invdtr < 0\n");
                QUIT_PLUTO(0);
            }
        }
    }

    if(v3 > 0){
        invdtphi = (vy - tan(rz)*vx)/(tan(rz)*x0 - y0);
    } else {
        invdtphi = (vy - tan(lz)*vx)/(tan(lz)*x0 - y0);
    }
    if(invdtphi < 0){
        printf("invdtphi < 0\n");
        QUIT_PLUTO(0);
    }

    if(invdtz > invdtr){
        if(invdtz > invdtphi){
            double dt = 1.0/invdtz;
            (*x1) = sqrt(sqr(x0 + vx*dt) + sqr(y0 + vy*dt));
            (*x3) = atan2(y0 + vy*dt, x0 + vx*dt);
            if(v2 > 0){
                (*x2) = ry;
                (*j) = (*j) + 1;
            } else {
                (*x2) = ly;
                (*j) = (*j) - 1;
            }
        } else {
            double dt = 1.0/invdtphi;
            (*x1) = sqrt(sqr(x0 + vx*dt) + sqr(y0 + vy*dt));
            (*x2) = (*x2) + v2*dt;
            if(v3 > 0){
                (*x3) = rz;
                (*k) = (*k) + 1;
            } else {
                (*x3) = lz;
                (*k) = (*k) - 1;
            }
        }
    } else {
        if(invdtr > invdtphi){
            double dt = 1.0/invdtr;
            (*x2) = (*x2) + v2*dt;
            (*x3) = atan2(y0 + vy*dt, x0 + vx*dt);
            if(v1 > 0){
                (*x1) = rx;
                (*i) = (*i) + 1;
            } else {
                (*x1) = lx;
                (*i) = (*i) - 1;
            }
        } else {
            double dt = 1.0/invdtphi;
            (*x1) = sqrt(sqr(x0 + vx*dt) + sqr(y0 + vy*dt));
            (*x2) = (*x2) + v2*dt;
            if(v3 > 0){
                (*x3) = rz;
                (*k) = (*k) + 1;
            } else {
                (*x3) = lz;
                (*k) = (*k) - 1;
            }
        }
    }
    CheckNanOrInfinity(*x1, "x1 = NaN\n");
    CheckNanOrInfinity(*x2, "x2 = NaN\n");
    CheckNanOrInfinity(*x3, "x3 = NaN\n");

#elif GEOMETRY == POLAR
    double invdtz = fabs(v3/dz);
    double invdtr = 0;
    double invdtphi = 0;

    double vx = v1*cos(grid->x[1][j]) - v2*sin(grid->x[1][j]);
    double vy = v1*sin(grid->x[1][j]) + v2*cos(grid->x[1][j]);

    double x0 = x1*cos(x2);
    double y0 = x1*sin(x2);

    double vxy2 = vx*vx + vy*vy;


    if(v1 > 0){
        double D = vxy2*xr*xr - sqr(vx*y0 + vy*x0);
        invdtr = vxy2/(sqrt(D) - vx*x0 - vy*y0);
        if(invdtr < 0){
            printf("invdtr < 0\n");
            QUIT_PLUTO(0);
        }
    } else {
        double D = vxy2*xl*xl - sqr(vx*y0 + vy*x0);
        if(D >= 0){
            invdtr = vxy2/(-sqrt(D) - vx*x0 - vy*y0);
            if(invdtr < 0){
                printf("invdtr < 0\n");
                QUIT_PLUTO(0);
            }
        }
    }

    if(v2 > 0){
        invdtphi = (vy - tan(ry)*vx)/(tan(ry)*x0 - y0);
    } else {
        invdtphi = (vy - tan(ly)*vx)/(tan(ly)*x0 - y0);
    }
    if(invdtphi < 0){
        printf("invdtphi < 0\n");
        QUIT_PLUTO(0);
    }

    if(invdtz > invdtr){
        if(invdtz > invdtphi){
            double dt = 1.0/invdtz;
            (*x1) = sqrt(sqr(x0 + vx*dt) + sqr(y0 + vy*dt));
            (*x2) = atan2(y0 + vy*dt, x0 + vx*dt);
            if(v2 > 0){
                (*x3) = rz;
                (*k) = (*k) + 1;
            } else {
                (*x3) = lz;
                (*k) = (*k) - 1;
            }
        } else {
            double dt = 1.0/invdtphi;
            (*x1) = sqrt(sqr(x0 + vx*dt) + sqr(y0 + vy*dt));
            (*x3) = (*x3) + v3*dt;
            if(v2 > 0){
                (*x2) = ry;
                (*j) = (*j) + 1;
            } else {
                (*x2) = ly;
                (*j) = (*j) - 1;
            }
        }
    } else {
        if(invdtr > invdtphi){
            double dt = 1.0/invdtr;
            (*x3) = (*x3) + v3*dt;
            (*x2) = atan2(y0 + vy*dt, x0 + vx*dt);
            if(v1 > 0){
                (*x1) = rx;
                (*i) = (*i) + 1;
            } else {
                (*x1) = lx;
                (*i) = (*i) - 1;
            }
        } else {
            double dt = 1.0/invdtphi;
            (*x1) = sqrt(sqr(x0 + vx*dt) + sqr(y0 + vy*dt));
            (*x3) = (*x3) + v3*dt;
            if(v2 > 0){
                (*x2) = ry;
                (*j) = (*j) + 1;
            } else {
                (*x2) = ly;
                (*j) = (*j) - 1;
            }
        }
    }
    CheckNanOrInfinity(*x1, "x1 = NaN\n");
    CheckNanOrInfinity(*x2, "x2 = NaN\n");
    CheckNanOrInfinity(*x3, "x3 = NaN\n");
#elif GEOMETRY == SPHERICAL
    double invdttheta = 0;
    double invdtr = 0;
    double invdtphi = 0;

    double vx = v1*sin(grid->x[1][j])*cos(grid->x[2][k]) + v2*cos(grid->x[1][j])*cos(grid->x[2][k]) - v3*sin(grid->x[2][k]);
    double vy = v1*sin(grid->x[1][j])*sin(grid->x[2][k]) + v2*cos(grid->x[1][j])*sin(grid->x[2][k]) + v3*cos(grid->x[2][k]);
    double vz = v1*cos(grid->x[1][j]) - v2*sin(grid->x[1][j]);

    double x0 = x1*sin(x2)*cos(x3);
    double y0 = x1*sin(x2)*sin(x3);
    double z0 = x1*cos(x2);


    if(v1 > 0){
        double D = sqr(x0*vx + y0*vy + z0*vz) - vsqr*(x1*x1 - rx*rx);
        invdtr = vsqr/(sqrt(D) - vx*x0 - vy*y0 - vz*z0);
        if(invdtr < 0){
            printf("invdtr < 0\n");
            QUIT_PLUTO(0);
        }
    } else {
        double D = sqr(x0*vx + y0*vy + z0*vz) - vsqr*(x1*x1 - lx*lx);
        if(D >= 0){
            invdtr = vsqr/(-sqrt(D) - vx*x0 - vy*y0);
            if(invdtr < 0){
                printf("invdtr < 0\n");
                QUIT_PLUTO(0);
            }
        }
    }

    if(v2 > 0){
        double cosry = cos(ry);
        double a = vz*vz - vsqr*cosry*cosry;
        double b = z0*vz - cosry*cosry*(vx*x0 + vy*y0 + vz*z0);
        double c = z0*z0 - x1*x1*cosry*cosry;
        double D = b*b - a*c;
        if(D > 0){
            invdttheta = a/(-b + sqrt(D));
            if(invdttheta < 0){
                printf("invdttheta < 0\n");
                QUIT_PLUTO(0);
            }
        }
    } else {
        double cosly = cos(ly);
        double a = vz*vz - vsqr*cosly*cosly;
        double b = z0*vz - cosly*cosly*(vx*x0 + vy*y0 + vz*z0);
        double c = z0*z0 - x1*x1*cosly*cosly;
        double D = b*b - a*c;
        if(D > 0){
            invdttheta = a/(-b + sqrt(D));
            if(invdttheta < 0){
                printf("invdttheta < 0\n");
                QUIT_PLUTO(0);
            }
        }
    }

    if(v3 > 0){
        invdtphi = (vy - tan(rz)*vx)/(tan(rz)*x0 - y0);
    } else {
        invdtphi = (vy - tan(lz)*vx)/(tan(lz)*x0 - y0);
    }
    if(invdtphi < 0){
        printf("invdtphi < 0\n");
        QUIT_PLUTO(0);
    }

    if(invdtphi > invdtr){
        if(invdtphi > invdttheta){
            double dt = 1.0/invdtphi;
            (*x1) = sqrt(sqr(x0 + vx*dt) + sqr(y0 + vy*dt) + sqr(z0 + vz*dt));
            (*x2) = acos(z0 + vz*dt, (*x1));
            if(v3 > 0){
                (*x3) = rz;
                (*k) = (*k) + 1;
            } else {
                (*x3) = lz;
                (*k) = (*k) - 1;
            }
        } else {
            double dt = 1.0/invdttheta;
            (*x1) = sqrt(sqr(x0 + vx*dt) + sqr(y0 + vy*dt) + sqr(z0 + vz*dt));
            (*x3) = atan2(y0 + vy*dt, x0 + vx*dt);
            if(v2 > 0){
                (*x2) = ry;
                (*j) = (*j) + 1;
            } else {
                (*x2) = ly;
                (*j) = (*j) - 1;
            }
        }
    } else {
        if(invdtr > invdttheta){
            double dt = 1.0/invdtr;
            if(v1 > 0){
                (*x1) = rx;
                (*i) = (*i) + 1;
            } else {
                (*x1) = lx;
                (*i) = (*i) - 1;
            }
            (*x2) = acos(z0 + vz*dt, (*x1));
            (*x3) = atan2(y0 + vy*dt, x0 + vx*dt);
        } else {
            double dt = 1.0/invdttheta;
            (*x1) = sqrt(sqr(x0 + vx*dt) + sqr(y0 + vy*dt) + sqr(z0 + vz*dt));
            (*x3) = atan2(y0 + vy*dt, x0 + vx*dt);
            if(v2 > 0){
                (*x2) = ry;
                (*j) = (*j) + 1;
            } else {
                (*x2) = ly;
                (*j) = (*j) - 1;
            }
        }
    }
    CheckNanOrInfinity(*x1, "x1 = NaN\n");
    CheckNanOrInfinity(*x2, "x2 = NaN\n");
    CheckNanOrInfinity(*x3, "x3 = NaN\n");
#else
#endif
#elif INCLUDE_JDIR
    double lx = grid->xl[0][*i];
    double rx = grid->xr[0][*i];
    double ly = grid->xl[1][*j];
    double ry = grid->xr[1][*j];

    double dx;
    double dy;
    if(v1 > 0){
        dx = (rx - (*x1))/grid->dx_dl[IDIR][*j][*i];
        if(dx == 0){
            (*i) = (*i) + 1;
            traceNextCell(grid, x1, x2, x3, v1, v2, v3, i, j, k);
            return;
        }
    } else {
        dx = ((*x1) - lx)/grid->dx_dl[IDIR][*j][*i];
        if(dx == 0){
            (*i) = (*i) - 1;
            traceNextCell(grid, x1, x2, x3, v1, v2, v3, i, j, k);
            return;
        }
    }
    if(v2 > 0){
        dy = (ry - (*x2))/grid->dx_dl[JDIR][*j][*i];
        if(dy == 0){
            (*j) = (*j) + 1;
            traceNextCell(grid, x1, x2, x3, v1, v2, v3, i, j, k);
            return;
        }
    } else {
        dy = ((*x2) - ly)/grid->dx_dl[JDIR][*j][*i];
        if(dy == 0){
            (*j) = (*j) - 1;
            traceNextCell(grid, x1, x2, x3, v1, v2, v3, i, j, k);
            return;
        }
    }
#if ((GEOMETRY == CARTESIAN) || (GEOMETRY == CYLINDRICAL))
    int a = *i;

    if(fabs(dx*v2) > fabs(dy*v1)){
        double dt = fabs(dy/v2);
        (*x1) = (*x1) + dt*v1;
        if(v2 > 0){
            (*x2) = ry;
            (*j) = (*j)+1;
        } else {
            (*x2) = ly;
            (*j) = (*j)-1;
        }
    } else {
        double dt = fabs(dx/v1);
        (*x2) = (*x2) + dt*v2;
        if(v1 > 0){
            (*x1) = rx;
            (*i) = (*i)+1;
        } else {
            (*x1 )= lx;
            (*i) = (*i)-1;
        }
    }
    CheckNanOrInfinity(*x1, "x1 = NaN\n");
    CheckNanOrInfinity(*x2, "x2 = NaN\n");
#elif GEOMETRY == POLAR
    double invdtr = 0;
    double invdtphi = 0;

    double vx = v1*cos(grid->x[1][j]) - v2*sin(grid->x[1][j]);
    double vy = v1*sin(grid->x[1][j]) + v2*cos(grid->x[1][j]);

    double x0 = x1*cos(x2);
    double y0 = x1*sin(x2);

    double vxy2 = vx*vx + vy*vy;


    if(v1 > 0){
        double D = vxy2*xr*xr - sqr(vx*y0 + vy*x0);
        invdtr = vxy2/(sqrt(D) - vx*x0 - vy*y0);
        if(invdtr < 0){
            printf("invdtr < 0\n");
            QUIT_PLUTO(0);
        }
    } else {
        double D = vxy2*xl*xl - sqr(vx*y0 + vy*x0);
        if(D >= 0){
            invdtr = vxy2/(-sqrt(D) - vx*x0 - vy*y0);
            if(invdtr < 0){
                printf("invdtr < 0\n");
                QUIT_PLUTO(0);
            }
        }
    }

    if(v2 > 0){
        invdtphi = (vy - tan(ry)*vx)/(tan(ry)*x0 - y0);
    } else {
        invdtphi = (vy - tan(ly)*vx)/(tan(ly)*x0 - y0);
    }
    if(invdtphi < 0){
        printf("invdtphi < 0\n");
        QUIT_PLUTO(0);
    }

    if(invdtr > invdtphi){
        double dt = 1.0/invdtr;
        (*x2) = atan2(y0 + vy*dt, x0 + vx*dt);
        if(v1 > 0){
            (*x1) = rx;
            (*i) = (*i) + 1;
        } else {
            (*x1) = lx;
            (*i) = (*i) - 1;
        }
    } else {
        double dt = 1.0/invdtphi;
        (*x1) = sqrt(sqr(x0 + vx*dt) + sqr(y0 + vy*dt));
        if(v2 > 0){
            (*x2) = ry;
            (*j) = (*j) + 1;
        } else {
            (*x2) = ly;
            (*j) = (*j) - 1;
        }
    }
    CheckNanOrInfinity(*x1, "x1 = NaN\n");
    CheckNanOrInfinity(*x2, "x2 = NaN\n");
#elif GEOMETRY == SPHERICAL
    double invdttheta = 0;
    double invdtr = 0;

    double vx = v1*sin(grid->x[1][j])*cos(grid->x[2][k]) + v2*cos(grid->x[1][j])*cos(grid->x[2][k]) - v3*sin(grid->x[2][k]);
    double vy = v1*sin(grid->x[1][j])*sin(grid->x[2][k]) + v2*cos(grid->x[1][j])*sin(grid->x[2][k]) + v3*cos(grid->x[2][k]);
    double vz = v1*cos(grid->x[1][j]) - v2*sin(grid->x[1][j]);

    double x0 = x1*sin(x2)*cos(x3);
    double y0 = x1*sin(x2)*sin(x3);
    double z0 = x1*cos(x2);


    if(v1 > 0){
        double D = sqr(x0*vx + y0*vy + z0*vz) - vsqr*(x1*x1 - rx*rx);
        invdtr = vsqr/(sqrt(D) - vx*x0 - vy*y0 - vz*z0);
        if(invdtr < 0){
            printf("invdtr < 0\n");
            QUIT_PLUTO(0);
        }
    } else {
        double D = sqr(x0*vx + y0*vy + z0*vz) - vsqr*(x1*x1 - lx*lx);
        if(D >= 0){
            invdtr = vsqr/(-sqrt(D) - vx*x0 - vy*y0);
            if(invdtr < 0){
                printf("invdtr < 0\n");
                QUIT_PLUTO(0);
            }
        }
    }

    if(v2 > 0){
        double cosry = cos(ry);
        double a = vz*vz - vsqr*cosry*cosry;
        double b = z0*vz - cosry*cosry*(vx*x0 + vy*y0 + vz*z0);
        double c = z0*z0 - x1*x1*cosry*cosry;
        double D = b*b - a*c;
        if(D > 0){
            invdttheta = a/(-b + sqrt(D));
            if(invdttheta < 0){
                printf("invdttheta < 0\n");
                QUIT_PLUTO(0);
            }
        }
    } else {
        double cosly = cos(ly);
        double a = vz*vz - vsqr*cosly*cosly;
        double b = z0*vz - cosly*cosly*(vx*x0 + vy*y0 + vz*z0);
        double c = z0*z0 - x1*x1*cosly*cosly;
        double D = b*b - a*c;
        if(D > 0){
            invdttheta = a/(-b + sqrt(D));
            if(invdttheta < 0){
                printf("invdttheta < 0\n");
                QUIT_PLUTO(0);
            }
        }
    }


    if(invdtr > invdttheta){
        double dt = 1.0/invdtr;
        if(v1 > 0){
            (*x1) = rx;
            (*i) = (*i) + 1;
        } else {
            (*x1) = lx;
            (*i) = (*i) - 1;
        }
        (*x2) = acos(z0 + vz*dt, (*x1));
    } else {
       double dt = 1.0/invdttheta;
       (*x1) = sqrt(sqr(x0 + vx*dt) + sqr(y0 + vy*dt) + sqr(z0 + vz*dt));
       if(v2 > 0){
            (*x2) = ry;
            (*j) = (*j) + 1;
       } else {
            (*x2) = ly;
            (*j) = (*j) - 1;
        }
    }

    CheckNanOrInfinity(*x1, "x1 = NaN\n");
    CheckNanOrInfinity(*x2, "x2 = NaN\n");
#else
#endif
#else
    if(v1 > 0){
        *x1 = grid->xr[0][*i];
        *i = (*i) + 1;
    } else if (v1 < 0){
        *x1 = grid->xl[0][*i];
        *i = (*i) - 1;
    } else {
        printLog("vx = 0 in trace cell\n");
        exit(0);
    }
    CheckNanOrInfinity(*x1, "x1 = NaN\n");
#endif
}
