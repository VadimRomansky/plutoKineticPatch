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

void updateShockFront(Data* d, Grid* grid){
    int i,j,k;
    FlagShock(d, grid);
    TOT_LOOP(k,j,i){
        d->shockWidth[k][j][i] = grid->dx[0][i];
        d->velocityJump[k][j][i] = 0.0;
        d->upstreamDensity[k][j][i] = d->Vc[RHO][k][j][i];
        d->downstreamDensity[k][j][i] = d->Vc[RHO][k][j][i];
    }

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
                pgradx = grid->dx_dl[IDIR][j][i]*(d->Vc[PRS][k][j][i+1] - d->Vc[PRS][k][j][i-1])/(grid->x[0][i+1] - grid->x[0][i-1]);
#endif
#if INCLUDE_JDIR
                pgrady = grid->dx_dl[JDIR][j][i]*(d->Vc[PRS][k][j+1][i] - d->Vc[PRS][k][j-1][i])/(grid->x[1][j+1] - grid->x[1][j-1]);
#endif
#if INCLUDE_KDIR
                pgradz = grid->dx_dl[KDIR][j][i]*(d->Vc[PRS][k+1][j][i] - d->Vc[PRS][k-1][j][i])/(grid->x[2][k+1] - grid->x[2][k-1]);
#endif
                traceNextCell(grid, &x, &y, &z, -pgradx, -pgrady, -pgradz, &upstreami, &upstreamj, &upstreamk);

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

#if INCLUDE_IDIR
        int Nout[1];
        int Nin[1];
        Nout[0] = NtoLeft;

        nleft = s->left[0];
        nright = s->right[0];
        tag1 = s->tag1[0];

        MPI_Sendrecv(Nout, 1, MPI_INT, nleft, tag1,
                     Nin, 1, MPI_INT, nright, tag1,
                     comm, &status);
                     
        NactiveTracers[0] = NactiveTracers[0] + Nin[0];

        int* outbuf;
        double* outbufd;
        if(Nout[0] != 0){
            outbuf = (int*) malloc(7*Nout[0]*sizeof(int));
            outbufd = (double*) malloc(7*Nout[0]*sizeof(double));
            putTracerListToArray(tracersToLeft, outbuf, outbufd);
            tracersToLeft = NULL;
            NtoLeft = 0;
        }
        int* inbuf;
        double* inbufd;
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

        nleft = s->left[1];
        nright = s->right[1];
        tag1 = s->tag1[1];

        MPI_Sendrecv(Nout, 1, MPI_INT, nleft, tag1,
                     Nin, 1, MPI_INT, nright, tag1,
                     comm, &status);
                     
        NactiveTracers[0] = NactiveTracers[0] + Nin[0];

        int* outbuf;
        double* outbufd;
        if(Nout[0] != 0){
            outbuf = (int*) malloc(7*Nout[0]*sizeof(int));
            outbufd = (double*) malloc(7*Nout[0]*sizeof(double));
            putTracerListToArray(tracersToDown, outbuf, outbufd);
            tracersToDown = NULL;
            NtoLeft = 0;
        }
        int* inbuf;
        double* inbufd;
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

        nleft = s->left[2];
        nright = s->right[2];
        tag1 = s->tag1[2];

        MPI_Sendrecv(Nout, 1, MPI_INT, nleft, tag1,
                     Nin, 1, MPI_INT, nright, tag1,
                     comm, &status);
                     
        NactiveTracers[0] = NactiveTracers[0] + Nin[0];

        int* outbuf;
        double* outbufd;
        if(Nout[0] != 0){
            outbuf = (int*) malloc(7*Nout[0]*sizeof(int));
            outbufd = (double*) malloc(7*Nout[0]*sizeof(double));
            putTracerListToArray(tracersToBack, outbuf, outbufd);
            tracersToBack = NULL;
            NtoLeft = 0;
        }
        int* inbuf;
        double* inbufd;
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
#else

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
                pgradx = grid->dx_dl[IDIR][j][i]*(d->Vc[PRS][k][j][i+1] - d->Vc[PRS][k][j][i-1])/(grid->x[0][i+1] - grid->x[0][i-1]);
#endif
#if INCLUDE_JDIR
                pgrady = grid->dx_dl[JDIR][j][i]*(d->Vc[PRS][k][j+1][i] - d->Vc[PRS][k][j-1][i])/(grid->x[1][j+1] - grid->x[1][j-1]);
#endif
#if INCLUDE_KDIR
                pgradz = grid->dx_dl[KDIR][j][i]*(d->Vc[PRS][k+1][j][i] - d->Vc[PRS][k-1][j][i])/(grid->x[2][k+1] - grid->x[2][k-1]);
#endif
                traceNextCell(grid, &x, &y, &z, -pgradx, -pgrady, -pgradz, &upstreami, &upstreamj, &upstreamk);
            }

            int downstreami = i;
            int downstreamj = j;
            int downstreamk = k;

            x = grid->x[0][i];
            y = grid->x[1][j];
            z = grid->x[2][k];

            while(!(d->flag[downstreamk][downstreamj][downstreami] & FLAG_ENTROPY)){
                double pgradx = 0;
                double pgrady = 0;
                double pgradz = 0;

#if INCLUDE_IDIR
                pgradx = grid->dx_dl[IDIR][j][i]*(d->Vc[PRS][k][j][i+1] - d->Vc[PRS][k][j][i-1])/(grid->x[0][i+1] - grid->x[0][i-1]);
#endif
#if INCLUDE_JDIR
                pgrady = grid->dx_dl[JDIR][j][i]*(d->Vc[PRS][k][j+1][i] - d->Vc[PRS][k][j-1][i])/(grid->x[1][j+1] - grid->x[1][j-1]);
#endif
#if INCLUDE_KDIR
                pgradz = grid->dx_dl[KDIR][j][i]*(d->Vc[PRS][k+1][j][i] - d->Vc[PRS][k-1][j][i])/(grid->x[2][k+1] - grid->x[2][k-1]);
#endif
                traceNextCell(grid, &x, &y, &z, pgradx, pgrady, pgradz, &downstreami, &downstreamj, &downstreamk);
            }

            double Vd1 = d->Vc[VX1][downstreamk][downstreamj][downstreami];
            double Vd2 = d->Vc[VX2][downstreamk][downstreamj][downstreami];
            double Vd3 = d->Vc[VX3][downstreamk][downstreamj][downstreami];

            double Vu1 = d->Vc[VX1][upstreamk][upstreamj][upstreami];
            double Vu2 = d->Vc[VX2][upstreamk][upstreamj][upstreami];
            double Vu3 = d->Vc[VX3][upstreamk][upstreamj][upstreami];

            double rhod = d->Vc[RHO][downstreamk][downstreamj][downstreami];
            double rhou = d->Vc[RHO][upstreamk][upstreamj][upstreami];

            double xd, yd, zd;
            double xu, yu, zu;
            double Vdx, Vdy, Vdz;
            double Vux, Vuy, Vuz;

            double width = 1E100;

#if GEOMETRY == CARTESIAN
            xd = grid->x[0][downstreami];
            yd = grid->x[1][downstreamj];
            zd = grid->x[2][downstreamk];

            xu = grid->x[0][upstreami];
            yu = grid->x[1][upstreamj];
            zu = grid->x[2][upstreamk];

            Vdx = Vd1;
            Vdy = Vd2;
            Vdz = Vd3;

            Vux = Vu1;
            Vuy = Vu2;
            Vuz = Vu3;
#elif GEOMETRY == CYLINDRICAL
            xd = grid->x[0][downstreami]*cos(grid->x[2][downstreamk]);
            yd = grid->x[0][downstreami]*sin(grid->x[2][downstreamk]);
            zd = grid->x[1][downstreamj];

            xu = grid->x[0][upstreami]*cos(grid->x[2][upstreamk]);
            yu = grid->x[0][upstreami]*sin(grid->x[2][upstreamk]);
            zu = grid->x[1][upstreamj];
#elif GEOMETRY == POLAR
            xd = grid->x[0][downstreami]*cos(grid->x[1][downstreamj]);
            yd = grid->x[0][downstreami]*sin(grid->x[1][downstreamj]);
            zd = grid->x[2][downstreamk];

            xu = grid->x[0][upstreami]*cos(grid->x[1][upstreamj]);
            yu = grid->x[0][upstreami]*sin(grid->x[1][upstreamj]);
            zu = grid->x[2][upstreamk];
#elif GEOMETRY == SPHERICAL
            xd = grid->x[0][downstreami]*sin(grid->x[1][downstreamj])*cos(grid->x[2][downstreamk]);
            yd = grid->x[0][downstreami]*sin(grid->x[1][downstreamj])*sin(grid->x[2][downstreamk]);
            zd = grid->x[0][downstreami]*cos(grid->x[1][downstreamj]);

            xu = grid->x[0][upstreami]*sin(grid->x[1][upstreamj])*cos(grid->x[2][upstreamk]);
            yu = grid->x[0][upstreami]*sin(grid->x[1][upstreamj])*sin(grid->x[2][upstreamk]);
            zu = grid->x[0][upstreami]*cos(grid->x[1][upstreamj]);
#else
#endif
            width = sqrt((xd-xu)*(xd-xu) + (yd-yu)*(yd-yu) + (zd-zu)*(zd-zu));
            double V = sqrt((Vdx - Vux)*(Vdx - Vux) + (Vdy - Vuy)*(Vdy - Vuy) + (Vdz - Vuz)*(Vdz - Vuz));

            d->shockWidth[k][j][i] = width;
            d->velocityJump[k][j][i] = V;
            d->upstreamDensity[k][j][i] = rhou;
            d->downstreamDensity[k][j][i] = rhod;
        }
    }
#endif
}

void traceNextCell(Grid* grid, double* x1, double* x2, double* x3, double vx, double vy, double vz, int* i, int* j, int* k){
    //todo proper stright lines for other geometries
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
#endif
}
