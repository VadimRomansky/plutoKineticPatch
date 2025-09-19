#include <stdio.h>
#include <stdbool.h>

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
    struct CellTracerNode_* next;
    struct CellTracerNode_* prev;
} CellTracerNode;

CellTracerNode* createTracer(int vrank0, int vi0, int vj0, int vk0, double x1v, double x2v, double x3v, int vi, int vj, int vk, double vv1, double vv2, double vv3, double vrho){
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
    tempNode->next = NULL;
    tempNode->prev = NULL;
    return tempNode;
}

CellTracerNode* addElementAfter(CellTracerNode* curNode, int vrank0, int vi0, int vj0, int vk0, double x1v, double x2v, double x3v, int vi, int vj, int vk, double vv1, double vv2, double vv3, double vrho){
    CellTracerNode* tempNode = createTracer(vrank0, vi0, vj0, vk0, x1v, x2v, x3v, vi, vj, vk, vv1, vv2, vv3, vrho);
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
		
		CellTracerNode* temp = createTracer(rank0, i0, j0, k0, x, y, z, i, j, k, vx, vy, vz, rho);
		
		temp->next = list;
		if(list != NULL){
			list->prev = temp;
		}
		list = temp;
	}
	return list;
}

void traceShock(Data* d, Grid* grid, int direction, double*** x1, double*** x2, double*** x3, double*** v1, double*** v2, double*** v3, double*** rho){
    if((direction != 1) && (direction != -1)){
        printf("direction for tracing must be 1 or -1\n");
        printLog("direction for tracing must be 1 or -1\n");
        QUIT_PLUTO(0);
    }

    int i, j, k;

#ifdef PARALLEL

    register int nd;
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
            CellTracerNode* tempTracer = createTracer(globrank, i,j,k, grid->x[0][i], grid->x[1][j], grid->x[2][k], i,j,k, 0.0, 0.0, 0.0, d->Vc[RHO][k][j][i]);
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
    MPI_Allreduce(a, NactiveTracers, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

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
            int upstreami = tracers->i;
            int upstreamj = tracers->j;
            int upstreamk = tracers->k;
            double x = tracers->x1;
            double y = tracers->x2;
            double z = tracers->x3;

            bool stopped = true;

            while(!(d->flag[upstreamk][upstreamj][upstreami] & FLAG_ENTROPY)){
#if INCLUDE_IDIR
                if(upstreami < IBEG){
                    CellTracerNode* temp = tracersToLeft;
                    tracersToLeft = tracers;
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
                    break;
                }
                if(upstreami > IEND){
                    CellTracerNode* temp = tracersToRight;
                    tracersToRight = tracers;
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
                    break;
                }
#endif
#if INCLUDE_JDIR
                if(upstreamj < JBEG){
                    CellTracerNode* temp = tracersToDown;
                    tracersToDown = tracers;
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
                    break;
                }
                if(upstreamj > JEND){
                    CellTracerNode* temp = tracersToUp;
                    tracersToUp = tracers;
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
                    break;
                }
#endif
#if INCLUDE_KDIR
                if(upstreamk < KBEG){
                    CellTracerNode* temp = tracersToBack;
                    tracersToBack = tracers;
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
                    break;
                }
                if(upstreamk > KEND){
                    CellTracerNode* temp = tracersToFront;
                    tracersToFront = tracers;
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
                    break;
                }
#endif
                double pgradx = 0;
                double pgrady = 0;
                double pgradz = 0;

#if INCLUDE_IDIR
                pgradx = grid->dx_dl[IDIR][upstreamj][upstreami]*(d->Vc[PRS][upstreamk][upstreamj][upstreami+1] - d->Vc[PRS][upstreamk][upstreamj][upstreami-1])/(grid->x[0][upstreami+1] - grid->x[0][upstreami-1]);
#endif
#if INCLUDE_JDIR
                pgrady = grid->dx_dl[JDIR][upstreamj][upstreami]*(d->Vc[PRS][upstreamk][upstreamj+1][upstreami] - d->Vc[PRS][upstreamk][upstreamj-1][upstreami])/(grid->x[1][upstreamj+1] - grid->x[1][upstreamj-1]);
#endif
#if INCLUDE_KDIR
                pgradz = grid->dx_dl[KDIR][upstreamj][upstreami]*(d->Vc[PRS][upstreamk+1][upstreamj][upstreami] - d->Vc[PRS][upstreamk-1][upstreamj][upstreami])/(grid->x[2][upstreamk+1] - grid->x[2][upstreamk-1]);
#endif
                traceNextCell(grid, &x, &y, &z, direction*pgradx, direction*pgrady, direction*pgradz, &upstreami, &upstreamj, &upstreamk);

                tracers->i = upstreami;
                tracers->j = upstreamj;
                tracers->k = upstreamk;
                tracers->x1 = x;
                tracers->x2 = y;
                tracers->x3 = z;
                tracers->v1 = d->Vc[VX1][k][j][i];
                tracers->v2 = d->Vc[VX2][k][j][i];
                tracers->v3 = d->Vc[VX3][k][j][i];
                tracers->rho = d->Vc[RHO][k][j][i];
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
            outbuf = (int*) malloc(7*Nout[0]*sizeof(int));
            outbufd = (double*) malloc(7*Nout[0]*sizeof(double));
            putTracerListToArray(tracersToLeft, outbuf, outbufd);
            tracersToLeft = NULL;
            NtoLeft = 0;
        }
        if(Nin[0] != 0){
            inbuf = (int*) malloc(7*Nin[0]*sizeof(int));
            inbufd = (double*) malloc(7*Nin[0]*sizeof(double));
            CellTracerNode* tracersFrom = putArrayToTracerList(inbuf, inbufd, Nin[0]);
            while(tracersFrom != NULL){
                tracersFrom->i = IEND;
                tracersFrom->x1 = grid->xr[0][IEND];
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
            outbuf = (int*) malloc(7*Nout[0]*sizeof(int));
            outbufd = (double*) malloc(7*Nout[0]*sizeof(double));
            putTracerListToArray(tracersToRight, outbuf, outbufd);
            tracersToRight = NULL;
            NtoRight = 0;
        }

        if(Nin[0] != 0){
            inbuf = (int*) malloc(7*Nin[0]*sizeof(int));
            inbufd = (double*) malloc(7*Nin[0]*sizeof(double));
            CellTracerNode* tracersFrom = putArrayToTracerList(inbuf, inbufd, Nin[0]);
            while(tracersFrom != NULL){
                tracersFrom->i = IBEG;
                tracersFrom->x1 = grid->xl[0][IBEG];
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
            outbuf = (int*) malloc(7*Nout[0]*sizeof(int));
            outbufd = (double*) malloc(7*Nout[0]*sizeof(double));
            putTracerListToArray(tracersToDown, outbuf, outbufd);
            tracersToDown = NULL;
            NtoLeft = 0;
        }

        if(Nin[0] != 0){
            inbuf = (int*) malloc(7*Nin[0]*sizeof(int));
            inbufd = (double*) malloc(7*Nin[0]*sizeof(double));
            CellTracerNode* tracersFrom = putArrayToTracerList(inbuf, inbufd, Nin[0]);
            while(tracersFrom != NULL){
                tracersFrom->j = JEND;
                tracersFrom->x2 = grid->xr[1][JEND];
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
            outbuf = (int*) malloc(7*Nout[0]*sizeof(int));
            outbufd = (double*) malloc(7*Nout[0]*sizeof(double));
            putTracerListToArray(tracersToUp, outbuf, outbufd);
            tracersToUp = NULL;
            NtoRight = 0;
        }

        if(Nin[0] != 0){
            inbuf = (int*) malloc(7*Nin[0]*sizeof(int));
            inbufd = (double*) malloc(7*Nin[0]*sizeof(double));
            CellTracerNode* tracersFrom = putArrayToTracerList(inbuf, inbufd, Nin[0]);
            while(tracersFrom != NULL){
                tracersFrom->j = JBEG;
                tracersFrom->x2 = grid->xl[1][JBEG];
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
            outbuf = (int*) malloc(7*Nout[0]*sizeof(int));
            outbufd = (double*) malloc(7*Nout[0]*sizeof(double));
            putTracerListToArray(tracersToBack, outbuf, outbufd);
            tracersToBack = NULL;
            NtoLeft = 0;
        }

        if(Nin[0] != 0){
            inbuf = (int*) malloc(7*Nin[0]*sizeof(int));
            inbufd = (double*) malloc(7*Nin[0]*sizeof(double));
            CellTracerNode* tracersFrom = putArrayToTracerList(inbuf, inbufd, Nin[0]);
            while(tracersFrom != NULL){
                tracersFrom-k = KEND;
                tracersFrom->x3 = grid->xr[2][KEND];
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
            outbuf = (int*) malloc(7*Nout[0]*sizeof(int));
            outbufd = (double*) malloc(7*Nout[0]*sizeof(double));
            putTracerListToArray(tracersToFront, outbuf, outbufd);
            tracersToFront = NULL;
            NtoRight = 0;
        }

        if(Nin[0] != 0){
            inbuf = (int*) malloc(7*Nin[0]*sizeof(int));
            inbufd = (double*) malloc(7*Nin[0]*sizeof(double));
            CellTracerNode* tracersFrom = putArrayToTracerList(inbuf, inbufd, Nin[0]);
            while(tracersFrom != NULL){
                tracersFrom->k = KBEG;
                tracersFrom->x3 = grid->xl[2][KBEG];
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

    for(int i = 1; i < nprocs; ++i){
        sdispls[i] = sdispls[i-1] + sendcounts[i];
        sdisplsd[i] = sdisplsd[i-1] + sendcounts[i];
        rdispls[i] = rdispls[i-1] + recvcounts[i];
        rdisplsd[i] = rdisplsd[i-1] + recvcounts[i];
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
        inbufd = (double*) malloc(doubleDataN*totalRecv*sizeof(int));
    }

    tempTracer = stoppedTracers;

    while(tempTracer != NULL){
        int rankt = tempTracer->rank0;

        outbuf[sdispls[rankt] + 3*srelposition[rankt]] = tempTracer->i0;
        outbuf[sdispls[rankt] + 3*srelposition[rankt]+1] = tempTracer->j0;
        outbuf[sdispls[rankt] + 3*srelposition[rankt]+2] = tempTracer->k0;

        outbufd[sdisplsd[rankt] + 7*srelposition[rankt]] = tempTracer->x1;
        outbufd[sdisplsd[rankt] + 7*srelposition[rankt]+1] = tempTracer->x2;
        outbufd[sdisplsd[rankt] + 7*srelposition[rankt]+2] = tempTracer->x3;
        outbufd[sdisplsd[rankt] + 7*srelposition[rankt]+3] = tempTracer->v1;
        outbufd[sdisplsd[rankt] + 7*srelposition[rankt]+4] = tempTracer->v2;
        outbufd[sdisplsd[rankt] + 7*srelposition[rankt]+5] = tempTracer->v3;
        outbufd[sdisplsd[rankt] + 7*srelposition[rankt]+6] = tempTracer->rho;

        srelposition[rankt] = srelposition[rankt] + 1;

        tempTracer = tempTracer->next;
    }

    MPI_Alltoallv(outbuf, sendcounts, sdispls, MPI_INT, inbuf, recvcounts, rdispls, MPI_INT, comm);
    MPI_Alltoallv(outbufd, sendcountsd, sdisplsd, MPI_DOUBLE, inbufd, recvcountsd, rdisplsd, MPI_DOUBLE, comm);

    if(totalSend > 0){
        free(outbuf);
        free(outbufd);
    }
    if(totalRecv > 0){
        free(inbuf);
        free(inbufd);
    }

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


    free(sendcounts);
    free(sendcountsd);
    free(sdispls);
    free(sdisplsd);
    free(srelposition);
    free(recvcounts);
    free(recvcountsd);
    free(rdispls);
    free(rdisplsd);
#else
#endif
    DOM_LOOP(k,j,i){
        if(!(d->flag[k][j][i] & FLAG_ENTROPY)){

            int upstreami = i;
            int upstreamj = j;
            int upstreamk = k;

            double x = grid->x[0][i];
            double y = grid->x[1][j];
            double z = grid->x[2][k];

            while(!(d->flag[upstreamk][upstreamj][upstreami] & FLAG_ENTROPY)){
                double pgradx = 0;
                double pgrady = 0;
                double pgradz = 0;

#if INCLUDE_IDIR
                pgradx = grid->dx_dl[IDIR][upstreamj][upstreami]*(d->Vc[PRS][upstreamk][upstreamj][upstreami+1] - d->Vc[PRS][upstreamk][upstreamj][upstreami-1])/(grid->x[0][upstreami+1] - grid->x[0][upstreami-1]);
#endif
#if INCLUDE_JDIR
                pgrady = grid->dx_dl[JDIR][upstreamj][upstreami]*(d->Vc[PRS][upstreamk][upstreamj+1][upstreami] - d->Vc[PRS][upstreamk][upstreamj-1][upstreami])/(grid->x[1][upstreamj+1] - grid->x[1][upstreamj-1]);
#endif
#if INCLUDE_KDIR
                pgradz = grid->dx_dl[KDIR][upstreamj][upstreami]*(d->Vc[PRS][upstreamk+1][upstreamj][upstreami] - d->Vc[PRS][upstreamk-1][upstreamj][upstreami])/(grid->x[2][upstreamk+1] - grid->x[2][upstreamk-1]);
#endif
                traceNextCell(grid, &x, &y, &z, direction*pgradx, direction*pgrady, direction*pgradz, &upstreami, &upstreamj, &upstreamk);
            }

            x1[k][j][i] = grid->x[0][upstreami];
            x2[k][j][i] = grid->x[1][upstreamj];
            x3[k][j][i] = grid->x[2][upstreamk];

            v1[k][j][i] = d->Vc[VX1][upstreamk][upstreamj][upstreami];
            v2[k][j][i] = d->Vc[VX2][upstreamk][upstreamj][upstreami];
            v3[k][j][i] = d->Vc[VX3][upstreamk][upstreamj][upstreami];

            rho[k][j][i] = d->Vc[RHO][upstreamk][upstreamj][upstreami];
        }
    }
}

void updateShockFront(Data* d, Grid* grid){
    int i,j,k;
    printf("evaluating shock\n");
    FlagShock(d, grid);
    TOT_LOOP(k,j,i){
        d->shockWidth[k][j][i] = grid->dx[0][i];
        d->velocityJump[k][j][i] = 0.0;
        d->upstreamDensity[k][j][i] = d->Vc[RHO][k][j][i];
        d->downstreamDensity[k][j][i] = d->Vc[RHO][k][j][i];
    }


    traceShock(d, grid, -1, d->upstreamx1, d->upstreamx2, d->upstreamx3, d->upstreamV1, d->upstreamV2, d->upstreamV3, d->upstreamDensity);
    traceShock(d, grid, 1, d->downstreamx1, d->downstreamx2, d->downstreamx3, d->downstreamV1, d->downstreamV2, d->downstreamV3, d->downstreamDensity);

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
#elif GEOMETRY == POLAR
            xd = d->downstreamx1[k][j][i]*cos(d->downstreamx2[k][j][i]);
            yd = d->downstreamx1[k][j][i]*sin(d->downstreamx2[k][j][i]);
            zd = d->downstreamx3[k][j][i];

            xu = d->upstreamx1[k][j][i]*cos(d->upstreamx2[k][j][i]);
            yu = d->upstreamx1[k][j][i]*sin(d->upstreamx2[k][j][i]);
            zu = d->upstreamx3[k][j][i];
#elif GEOMETRY == SPHERICAL
            xd = d->downstreamx1[k][j][i]*sin(d->downstreamx2[k][j][i])*cos(d->downstreamx3[k][j][i]);
            yd = d->downstreamx1[k][j][i]*sin(d->downstreamx2[k][j][i])*sin(d->downstreamx3[k][j][i]);
            zd = d->downstreamx1[k][j][i]*cos(d->downstreamx2[k][j][i]);

            xu = d->upstreamx1[k][j][i]*sin(grid->x[1][upstreamj])*cos(d->upstreamx3[k][j][i]);
            yu = d->upstreamx1[k][j][i]*sin(d->upstreamx2[k][j][i])*sin(d->upstreamx3[k][j][i]);
            zu = d->upstreamx1[k][j][i]*cos(d->upstreamx2[k][j][i]);
#else
#endif
            double width = sqrt((xd-xu)*(xd-xu) + (yd-yu)*(yd-yu) + (zd-zu)*(zd-zu));
            double V = sqrt((Vdx - Vux)*(Vdx - Vux) + (Vdy - Vuy)*(Vdy - Vuy) + (Vdz - Vuz)*(Vdz - Vuz));

            d->shockWidth[k][j][i] = width;
            d->velocityJump[k][j][i] = V;
        }
    }
}

void traceNextCell(Grid* grid, double* x1, double* x2, double* x3, double vx, double vy, double vz, int* i, int* j, int* k){
    //todo proper stright lines for other geometries
    double v2 = vx*vx + vy*vy + vz*vz;
    if(v2 <= 0){
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
    } else {
        dx = (*x1 - lx)/grid->dx_dl[IDIR][*j][*i];
    }
    if(vy > 0){
        dy = (ry - *x2)/grid->dx_dl[JDIR][*j][*i];
    } else {
        dy = (*x2 - ly)/grid->dx_dl[JDIR][*j][*i];
    }
    if(vz > 0){
        dz = (rz - *x3)/grid->dx_dl[KDIR][*j][*i];
    } else {
        dz = (*x3 - lz)/grid->dx_dl[KDIR][*j][*i];
    }

    if(fabs(dx*vy) > fabs(dy*vx)){
        if(fabs(dz*vy) > fabs(dy*vz)){
            double dt = fabs(dy/vy);
            *x1 = *x1 + dt*vx;
            *x3 = *x3 + dt*vz;
            if(vy > 0){
                *x2 = ry;
                *j = (*j)+1;
            } else {
                *x2 = ly;
                *j = (*j)-1;
            }
        } else {
            double dt = fabs(dz/vz);
            *x1 = *x1 + dt*vx;
            *x2 = *x2 + dt*vy;
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
            *x2 = *x2 + dt*vy;
            *x3 = *x3 + dt*vz;
            if(vx > 0){
                *x1 = rx;
                *i = (*i)+1;
            } else {
                *x1 = lx;
                *i = (*i)-1;
            }
        } else {
            double dt = fabs(dz/vz);
            *x1 = *x1 + dt*vx;
            *x2 = *x2 + dt*vy;
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
#elif INCLUDE_JDIR
    int a = *i;
    double lx = grid->xl[0][*i];
    double rx = grid->xr[0][*i];
    double ly = grid->xl[1][*j];
    double ry = grid->xr[1][*j];

    double dx;
    double dy;
    if(vx > 0){
        dx = (rx - *x1)/grid->dx_dl[IDIR][*j][*i];
    } else {
        dx = (*x1 - lx)/grid->dx_dl[IDIR][*j][*i];
    }
    if(vy > 0){
        dy = (ry - *x2)/grid->dx_dl[JDIR][*j][*i];
    } else {
        dy = (*x2 - ly)/grid->dx_dl[JDIR][*j][*i];
    }

    if(fabs(dx*vy) > fabs(dy*vx)){
        double dt = fabs(dy/vy);
        *x1 = *x1 + dt*vx;
        if(vy > 0){
            *x2 = ry;
            *j = (*j)+1;
        } else {
            *x2 = ly;
            *j = (*j)-1;
        }
    } else {
        double dt = fabs(dx/vx);
        *x2 = *x2 + dt*vy;
        if(vx > 0){
            *x1 = rx;
            *i = (*i)+1;
        } else {
            *x1 = lx;
            *i = (*i)-1;
        }
    }
    CheckNanOrInfinity(*x1, "x1 = NaN\n");
    CheckNanOrInfinity(*x2, "x2 = NaN\n");
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
