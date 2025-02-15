#include "pluto.h"

#include "matrixElement.h"
#include "specialmath.h"

#if TURBULENT_FIELD == YES
void AdvanceTurbulentField(Data *d, timeStep *Dts, double dt, Grid *grid){
    int k,j,i;
    //printf("trubulent field\n");
    //return;
    TOT_LOOP(k,j,i){
        for(int l = 0; l < NTURB; ++l){
                MatrixElementNode* curNode = d->turbulent_matrix[k][j][i][l];
                while(curNode != NULL){
                    MatrixElementNode* tempNode = curNode;
                    curNode = curNode->next;
                    free(tempNode);
                }
                curNode = (MatrixElementNode*) malloc(sizeof(MatrixElementNode));
                curNode->element = createMatrixElement(1.0, k,j,i,l);
                curNode->next = NULL;
                curNode->prev = NULL;
                d->turbulent_rightPart[k][j][i][l] = 0.0;
                d->turbulent_matrix[k][j][i][l] = curNode;
        }
    }

    DOM_LOOP(k, j, i){
#if INCLUDE_IDIR
        if(grid->lbound[0] != 0){
            if(grid->lbound[0] != PERIODIC){
            if(i == IBEG){
                for(int l = 0; l < NTURB; ++l){
                    MatrixElementNode* curNode = d->turbulent_matrix[k][j][i][l];
                    curNode = addElement(curNode, -1.0, k,j,i+1,l);
                }
                continue;
            }
            }
        }

        if(grid->rbound[0] != 0){
            if(grid->rbound[0] != PERIODIC){
            if(i == IEND){
                for(int l = 0; l < NTURB; ++l){
                    MatrixElementNode* curNode = d->turbulent_matrix[k][j][i][l];
                    curNode = addElement(curNode, -1.0, k,j,i-1,l);
                }
                continue;
            }
            }
        }
#endif

#if INCLUDE_JDIR
        if(grid->lbound[1] != 0){
            if(grid->lbound[1] != PERIODIC){
                if(j == JBEG){
                    for(int l = 0; l < NTURB; ++l){
                        MatrixElementNode* curNode = d->turbulent_matrix[k][j][i][l];
                        curNode = addElement(curNode, -1.0, k,j + 1,i,l);
                    }
                    continue;
                }
            }
        }

        if(grid->rbound[1] != 0){
            if(grid->rbound[1] != PERIODIC){
                if(j == JEND){
                    for(int l = 0; l < NTURB; ++l){
                        MatrixElementNode* curNode = d->turbulent_matrix[k][j][i][l];
                        curNode = addElement(curNode, -1.0, k,j - 1,i,l);
                    }
                    continue;
                }
            }
        }
#endif

#if INCLUDE_KDIR
        if(grid->lbound[2] != 0){
            if(grid->lbound[2] != PERIODIC){
                if(k == KBEG){
                    for(int l = 0; l < NTURB; ++l){
                        MatrixElementNode* curNode = d->turbulent_matrix[k][j][i][l];
                        curNode = addElement(curNode, -1.0, k + 1,j,i,l);
                    }
                    continue;
                }
            }
        }

        if(grid->rbound[2] != 0){
            if(grid->rbound[2] != PERIODIC){
                if(k == KEND){
                    for(int l = 0; l < NTURB; ++l){
                        MatrixElementNode* curNode = d->turbulent_matrix[k][j][i][l];
                        curNode = addElement(curNode, -1.0, k - 1,j,i,l);
                    }
                    continue;
                }
            }
        }

#endif
        double value;
        for(int l = 0; l < NTURB; ++l){
            d->turbulent_rightPart[k][j][i][l] = d->Wt[k][j][i][l];
            MatrixElementNode* curNode = d->turbulent_matrix[k][j][i][l];

            double u1 = d->Vc[VX1][k][j][i];
            double u2 = d->Vc[VX2][k][j][i];
            double u3 = d->Vc[VX3][k][j][i];

            double G = 0;

            curNode = addElement(curNode, -dt*G, k, j, i, l);

#if GEOMETRY == CARTESIAN
#if INCLUDE_IDIR
            // 1.5 div u W - 0.5 u grad W
            if(u1 > 0){
                value = 1.5 *dt*u1/(grid->x[0][i] - grid->x[0][i-1]) - 0.5*dt*u1/(grid->x[0][i] - grid->x[0][i-1]);
                curNode = addElement(curNode, value, k, j, i, l);
                if((value != value) || (value*0 != value*0)){
                    printf("1 turbulent value = NaN for (%d %d %d %d) (%d %d %d %d)\n", k,j,i,l, curNode->element.k, curNode->element.j, curNode->element.i, curNode->element.l);
                    QUIT_PLUTO(1);
                }

                value = -1.5*dt*d->Vc[VX1][k][j][i-1]/(grid->x[0][i] - grid->x[0][i-1]) + 0.5*dt*u1/(grid->x[0][i] - grid->x[0][i-1]);
                curNode = addElement(curNode, value, k, j, i-1, l);
                if((value != value) || (value*0 != value*0)){
                    printf("2 turbulent value = NaN for (%d %d %d %d) (%d %d %d %d)\n", k,j,i,l, curNode->element.k, curNode->element.j, curNode->element.i, curNode->element.l);
                    QUIT_PLUTO(1);
                }

            } else {
                value = 1.5 *dt*u1/(grid->x[0][i+1] - grid->x[0][i]) - 0.5*dt*u1/(grid->x[0][i+1] - grid->x[0][i]);
                curNode = addElement(curNode, value, k, j, i+1, l);
                if((value != value) || (value*0 != value*0)){
                    printf("3 turbulent value = NaN for (%d %d %d %d) (%d %d %d %d)\n", k,j,i,l, curNode->element.k, curNode->element.j, curNode->element.i, curNode->element.l);
                    QUIT_PLUTO(1);
                }

                value = -1.5*dt*d->Vc[VX1][k][j][i]/(grid->x[0][i+1] - grid->x[0][i]) + 0.5*dt*u1/(grid->x[0][i+1] - grid->x[0][i]);
                curNode = addElement(curNode, value, k, j, i, l);
                if((value != value) || (value*0 != value*0)){
                    printf("4 turbulent value = NaN for (%d %d %d %d) (%d %d %d %d)\n", k,j,i,l, curNode->element.k, curNode->element.j, curNode->element.i, curNode->element.l);
                    QUIT_PLUTO(1);
                }

            }
#endif
#if INCLUDE_JDIR
            if(u2 > 0){
                value = 1.5 *dt*u2/(grid->x[1][j] - grid->x[1][j-1]) - 0.5*dt*u2/(grid->x[1][j] - grid->x[1][j-1]);
                curNode = addElement(curNode, value, k, j, i, l);
                if((value != value) || (value*0 != value*0)){
                    printf("5 turbulent value = NaN for (%d %d %d %d) (%d %d %d %d)\n", k,j,i,l, curNode->element.k, curNode->element.j, curNode->element.i, curNode->element.l);
                    QUIT_PLUTO(1);
                }

                value = -1.5*dt*d->Vc[VX2][k][j-1][i]/(grid->x[1][j] - grid->x[1][j-1]) + 0.5*dt*u2/(grid->x[1][j] - grid->x[1][j-1]);
                curNode = addElement(curNode, value, k, j-1, i, l);
                if((value != value) || (value*0 != value*0)){
                    printf("6 turbulent value = NaN for (%d %d %d %d) (%d %d %d %d)\n", k,j,i,l, curNode->element.k, curNode->element.j, curNode->element.i, curNode->element.l);
                    QUIT_PLUTO(1);
                }

            } else {
                value = 1.5 *dt*d->Vc[VX2][k][j+1][i]/(grid->x[1][j+1] - grid->x[1][j]) - 0.5*dt*u2/(grid->x[1][j+1] - grid->x[1][j]);
                curNode = addElement(curNode, value, k, j+1, i, l);
                if((value != value) || (value*0 != value*0)){
                    printf("7 turbulent value = NaN for (%d %d %d %d) (%d %d %d %d)\n", k,j,i,l, curNode->element.k, curNode->element.j, curNode->element.i, curNode->element.l);
                    QUIT_PLUTO(1);
                }

                value = -1.5*dt*u2/(grid->x[1][j+1] - grid->x[1][j]) + 0.5*dt*u2/(grid->x[1][j+1] - grid->x[1][j]);
                curNode = addElement(curNode, value, k, j, i, l);
                if((value != value) || (value*0 != value*0)){
                    printf("8 turbulent value = NaN for (%d %d %d %d) (%d %d %d %d)\n", k,j,i,l, curNode->element.k, curNode->element.j, curNode->element.i, curNode->element.l);
                    QUIT_PLUTO(1);
                }

            }
#endif
#if INCLUDE_KDIR
            if(u3 > 0){
                value = 1.5 *dt*u3/(grid->x[2][k] - grid->x[2][k-1]) - 0.5*dt*u3/(grid->x[2][k] - grid->x[2][k-1]);
                curNode = addElement(curNode, value, k, j, i, l);
                if((value != value) || (value*0 != value*0)){
                    printf("9 turbulent value = NaN for (%d %d %d %d) (%d %d %d %d)\n", k,j,i,l, curNode->element.k, curNode->element.j, curNode->element.i, curNode->element.l);
                    QUIT_PLUTO(1);
                }

                value = -1.5*dt*d->Vc[VX3][k-1][j][i]/(grid->x[2][k] - grid->x[2][k-1]) + 0.5*dt*u3/(grid->x[2][k] - grid->x[2][k-1]);
                curNode = addElement(curNode, value, k-1, j, i, l);
                if((value != value) || (value*0 != value*0)){
                    printf("10 turbulent value = NaN for (%d %d %d %d) (%d %d %d %d)\n", k,j,i,l, curNode->element.k, curNode->element.j, curNode->element.i, curNode->element.l);
                    QUIT_PLUTO(1);
                }

            } else {
                value = 1.5 *dt*u3/(grid->x[2][k+1] - grid->x[2][k]) - 0.5*dt*u3/(grid->x[2][k+1] - grid->x[2][k]);
                curNode = addElement(curNode, value, k+1, j, i, l);
                if((value != value) || (value*0 != value*0)){
                    printf("11 turbulent value = NaN for (%d %d %d %d) (%d %d %d %d)\n", k,j,i,l, curNode->element.k, curNode->element.j, curNode->element.i, curNode->element.l);
                    QUIT_PLUTO(1);
                }

                value = -1.5*dt*d->Vc[VX3][k][j][i]/(grid->x[2][k+1] - grid->x[2][k]) + 0.5*dt*u3/(grid->x[2][k+1] - grid->x[2][k]);
                curNode = addElement(curNode, value, k, j, i, l);
                if((value != value) || (value*0 != value*0)){
                    printf("12 turbulent value = NaN for (%d %d %d %d) (%d %d %d %d)\n", k,j,i,l, curNode->element.k, curNode->element.j, curNode->element.i, curNode->element.l);
                    QUIT_PLUTO(1);
                }

            }
#endif
#elif GEOMETRY == CYLINDRICAL
#if INCLUDE_IDIR
            if(u1 > 0){
                value = 1.5 *dt*u1/(grid->x[0][i] - grid->x[0][i-1]) - 0.5*dt*u1/(grid->x[0][i] - grid->x[0][i-1]);
                curNode = addElement(curNode, value, k, j, i, l);
                if((value != value) || (value*0 != value*0)){
                    printf("13 turbulent value = NaN for (%d %d %d %d) (%d %d %d %d)\n", k,j,i,l, curNode->element.k, curNode->element.j, curNode->element.i, curNode->element.l);
                    QUIT_PLUTO(1);
                }

                value = -1.5*dt*grid->x[i-1][0]*d->Vc[VX1][k][j][i-1]/(grid->x[0][i]*(grid->x[0][i] - grid->x[0][i-1])) + 0.5*dt*u1/(grid->x[0][i] - grid->x[0][i-1]);
                curNode = addElement(curNode, value, k, j, i-1, l);
                if((value != value) || (value*0 != value*0)){
                    printf("14 turbulent value = NaN for (%d %d %d %d) (%d %d %d %d)\n", k,j,i,l, curNode->element.k, curNode->element.j, curNode->element.i, curNode->element.l);
                    QUIT_PLUTO(1);
                }

            } else {
                value = 1.5 *dt*grid->x[0][i+1]*u1/(grid->x[0][i]*(grid->x[0][i+1] - grid->x[0][i])) - 0.5*dt*u1/(grid->x[0][i+1] - grid->x[0][i]);
                curNode = addElement(curNode, value, k, j, i+1, l);
                if((value != value) || (value*0 != value*0)){
                    printf("15 turbulent value = NaN for (%d %d %d %d) (%d %d %d %d)\n", k,j,i,l, curNode->element.k, curNode->element.j, curNode->element.i, curNode->element.l);
                    QUIT_PLUTO(1);
                }

                value = -1.5*dt*d->Vc[VX1][k][j][i]/(grid->x[0][i+1] - grid->x[0][i]) + 0.5*dt*u1/(grid->x[0][i+1] - grid->x[0][i]);
                curNode = addElement(curNode, value, k, j, i, l);
                if((value != value) || (value*0 != value*0)){
                    printf("16 turbulent value = NaN for (%d %d %d %d) (%d %d %d %d)\n", k,j,i,l, curNode->element.k, curNode->element.j, curNode->element.i, curNode->element.l);
                    QUIT_PLUTO(1);
                }

            }
#endif
#if INCLUDE_JDIR
            if(u2 > 0){
                value = 1.5 *dt*u2/(grid->x[1][j] - grid->x[1][j-1]) - 0.5*dt*u2/(grid->x[1][j] - grid->x[1][j-1]);
                curNode = addElement(curNode, value, k, j, i, l);
                if((value != value) || (value*0 != value*0)){
                    printf("17 turbulent value = NaN for (%d %d %d %d) (%d %d %d %d)\n", k,j,i,l, curNode->element.k, curNode->element.j, curNode->element.i, curNode->element.l);
                    QUIT_PLUTO(1);
                }

                value = -1.5*dt*d->Vc[VX2][k][j-1][i]/(grid->x[1][j] - grid->x[1][j-1]) + 0.5*dt*u2/(grid->x[1][j] - grid->x[1][j-1]);
                curNode = addElement(curNode, value, k, j-1, i, l);
                if((value != value) || (value*0 != value*0)){
                    printf("18 turbulent value = NaN for (%d %d %d %d) (%d %d %d %d)\n", k,j,i,l, curNode->element.k, curNode->element.j, curNode->element.i, curNode->element.l);
                    QUIT_PLUTO(1);
                }

            } else {
                value = 1.5 *dt*d->Vc[VX2][k][j+1][i]/(grid->x[1][j+1] - grid->x[1][j]) - 0.5*dt*u2/(grid->x[1][j+1] - grid->x[1][j]);
                curNode = addElement(curNode, value, k, j+1, i, l);
                if((value != value) || (value*0 != value*0)){
                    printf("19 turbulent value = NaN for (%d %d %d %d) (%d %d %d %d)\n", k,j,i,l, curNode->element.k, curNode->element.j, curNode->element.i, curNode->element.l);
                    QUIT_PLUTO(1);
                }

                value = -1.5*dt*u2/(grid->x[1][j+1] - grid->x[1][j]) + 0.5*dt*u2/(grid->x[1][j+1] - grid->x[1][j]);
                curNode = addElement(curNode, value, k, j, i, l);
                if((value != value) || (value*0 != value*0)){
                    printf("20 turbulent value = NaN for (%d %d %d %d) (%d %d %d %d)\n", k,j,i,l, curNode->element.k, curNode->element.j, curNode->element.i, curNode->element.l);
                    QUIT_PLUTO(1);
                }

            }
#endif
#if INCLUDE_KDIR
            if(u3 > 0){
                value = 1.5 *dt*u3/(grid->x[0][i]*(grid->x[2][k] - grid->x[2][k-1])) - 0.5*dt*u3/(grid->x[0][i]*(grid->x[2][k] - grid->x[2][k-1]));
                curNode = addElement(curNode, value, k, j, i, l);
                if((value != value) || (value*0 != value*0)){
                    printf("21 turbulent value = NaN for (%d %d %d %d) (%d %d %d %d)\n", k,j,i,l, curNode->element.k, curNode->element.j, curNode->element.i, curNode->element.l);
                    QUIT_PLUTO(1);
                }

                value = -1.5*dt*d->Vc[VX3][k-1][j][i]/(grid->x[0][i]*(grid->x[2][k] - grid->x[2][k-1])) + 0.5*dt*u3/(grid->x[0][i]*(grid->x[2][k] - grid->x[2][k-1]));
                curNode = addElement(curNode, value, k-1, j, i, l);
                if((value != value) || (value*0 != value*0)){
                    printf("22 turbulent value = NaN for (%d %d %d %d) (%d %d %d %d)\n", k,j,i,l, curNode->element.k, curNode->element.j, curNode->element.i, curNode->element.l);
                    QUIT_PLUTO(1);
                }

            } else {
                value = 1.5 *dt*u3/(grid->x[0][i]*(grid->x[2][k+1] - grid->x[2][k])) - 0.5*dt*u3/(grid->x[0][i]*(grid->x[2][k+1] - grid->x[2][k]));
                curNode = addElement(curNode, value, k+1, j, i, l);
                if((value != value) || (value*0 != value*0)){
                    printf("23 turbulent value = NaN for (%d %d %d %d) (%d %d %d %d)\n", k,j,i,l, curNode->element.k, curNode->element.j, curNode->element.i, curNode->element.l);
                    QUIT_PLUTO(1);
                }

                value = -1.5*dt*d->Vc[VX3][k][j][i]/(grid->x[0][i]*(grid->x[2][k+1] - grid->x[2][k])) + 0.5*dt*u3/(grid->x[0][i]*(grid->x[2][k+1] - grid->x[2][k]));
                curNode = addElement(curNode, value, k, j, i, l);
                if((value != value) || (value*0 != value*0)){
                    printf("24 turbulent value = NaN for (%d %d %d %d) (%d %d %d %d)\n", k,j,i,l, curNode->element.k, curNode->element.j, curNode->element.i, curNode->element.l);
                    QUIT_PLUTO(1);
                }

            }
#endif
#elif GEOMETRY == POLAR
#if INCLUDE_IDIR
            if(u1 > 0){
                value = 1.5 *dt*u1/(grid->x[0][i] - grid->x[0][i-1]) - 0.5*dt*u1/(grid->x[0][i] - grid->x[0][i-1]);
                curNode = addElement(curNode, value, k, j, i, l);
                if((value != value) || (value*0 != value*0)){
                    printf("25 turbulent value = NaN for (%d %d %d %d) (%d %d %d %d)\n", k,j,i,l, curNode->element.k, curNode->element.j, curNode->element.i, curNode->element.l);
                    QUIT_PLUTO(1);
                }

                value = -1.5*dt*grid->x[i-1][0]*d->Vc[VX1][k][j][i-1]/(grid->x[0][i]*(grid->x[0][i] - grid->x[0][i-1])) + 0.5*dt*u1/(grid->x[0][i] - grid->x[0][i-1]);
                curNode = addElement(curNode, value, k, j, i-1, l);
                if((value != value) || (value*0 != value*0)){
                    printf("26 turbulent value = NaN for (%d %d %d %d) (%d %d %d %d)\n", k,j,i,l, curNode->element.k, curNode->element.j, curNode->element.i, curNode->element.l);
                    QUIT_PLUTO(1);
                }

            } else {
                value = 1.5 *dt*grid->x[0][i+1]*u1/(grid->x[0][i]*(grid->x[0][i+1] - grid->x[0][i])) - 0.5*dt*u1/(grid->x[0][i+1] - grid->x[0][i]);
                curNode = addElement(curNode, value, k, j, i+1, l);
                if((value != value) || (value*0 != value*0)){
                    printf("27 turbulent value = NaN for (%d %d %d %d) (%d %d %d %d)\n", k,j,i,l, curNode->element.k, curNode->element.j, curNode->element.i, curNode->element.l);
                    QUIT_PLUTO(1);
                }

                value = -1.5*dt*d->Vc[VX1][k][j][i]/(grid->x[0][i+1] - grid->x[0][i]) + 0.5*dt*u1/(grid->x[0][i+1] - grid->x[0][i]);
                curNode = addElement(curNode, value, k, j, i, l);
                if((value != value) || (value*0 != value*0)){
                    printf("28 turbulent value = NaN for (%d %d %d %d) (%d %d %d %d)\n", k,j,i,l, curNode->element.k, curNode->element.j, curNode->element.i, curNode->element.l);
                    QUIT_PLUTO(1);
                }

            }
#endif
#if INCLUDE_JDIR
            if(u2 > 0){
                value = 1.5 *dt*u2/(grid->x[0][i]*(grid->x[1][j] - grid->x[1][j-1])) - 0.5*dt*u2/(grid->x[0][i]*(grid->x[1][j] - grid->x[1][j-1]));
                curNode = addElement(curNode, value, k, j, i, l);
                if((value != value) || (value*0 != value*0)){
                    printf("29 turbulent value = NaN for (%d %d %d %d) (%d %d %d %d)\n", k,j,i,l, curNode->element.k, curNode->element.j, curNode->element.i, curNode->element.l);
                    QUIT_PLUTO(1);
                }

                value = -1.5*dt*d->Vc[VX2][k][j-1][i]/(grid->x[0][i]*(grid->x[1][j] - grid->x[1][j-1])) + 0.5*dt*u2/(grid->x[0][i]*(grid->x[1][j] - grid->x[1][j-1]));
                curNode = addElement(curNode, value, k, j-1, i, l);
                if((value != value) || (value*0 != value*0)){
                    printf("30 turbulent value = NaN for (%d %d %d %d) (%d %d %d %d)\n", k,j,i,l, curNode->element.k, curNode->element.j, curNode->element.i, curNode->element.l);
                    QUIT_PLUTO(1);
                }

            } else {
                value = 1.5 *dt*d->Vc[VX2][k][j+1][i]/(grid->x[0][i]*(grid->x[1][j+1] - grid->x[1][j])) - 0.5*dt*u2/(grid->x[0][i]*(grid->x[1][j+1] - grid->x[1][j]));
                curNode = addElement(curNode, value, k, j+1, i, l);
                if((value != value) || (value*0 != value*0)){
                    printf("31 turbulent value = NaN for (%d %d %d %d) (%d %d %d %d)\n", k,j,i,l, curNode->element.k, curNode->element.j, curNode->element.i, curNode->element.l);
                    QUIT_PLUTO(1);
                }

                value = -1.5*dt*u2/(grid->x[0][i]*(grid->x[1][j+1] - grid->x[1][j])) + 0.5*dt*u2/(grid->x[0][i]*(grid->x[1][j+1] - grid->x[1][j]));
                curNode = addElement(curNode, value, k, j, i, l);
                if((value != value) || (value*0 != value*0)){
                    printf("32 turbulent value = NaN for (%d %d %d %d) (%d %d %d %d)\n", k,j,i,l, curNode->element.k, curNode->element.j, curNode->element.i, curNode->element.l);
                    QUIT_PLUTO(1);
                }

            }
#endif
#if INCLUDE_KDIR
            if(u3 > 0){
                value = 1.5 *dt*u3/(grid->x[2][k] - grid->x[2][k-1]) - 0.5*dt*u3/(grid->x[2][k] - grid->x[2][k-1]);
                curNode = addElement(curNode, value, k, j, i, l);
                if((value != value) || (value*0 != value*0)){
                    printf("33 turbulent value = NaN for (%d %d %d %d) (%d %d %d %d)\n", k,j,i,l, curNode->element.k, curNode->element.j, curNode->element.i, curNode->element.l);
                    QUIT_PLUTO(1);
                }

                value = -1.5*dt*d->Vc[VX3][k-1][j][i]/(grid->x[2][k] - grid->x[2][k-1]) + 0.5*dt*u3/(grid->x[2][k] - grid->x[2][k-1]);
                curNode = addElement(curNode, value, k-1, j, i, l);
                if((value != value) || (value*0 != value*0)){
                    printf("34 turbulent value = NaN for (%d %d %d %d) (%d %d %d %d)\n", k,j,i,l, curNode->element.k, curNode->element.j, curNode->element.i, curNode->element.l);
                    QUIT_PLUTO(1);
                }

            } else {
                value = 1.5 *dt*u3/(grid->x[2][k+1] - grid->x[2][k]) - 0.5*dt*u3/(grid->x[2][k+1] - grid->x[2][k]);
                curNode = addElement(curNode, value, k+1, j, i, l);
                if((value != value) || (value*0 != value*0)){
                    printf("35 turbulent value = NaN for (%d %d %d %d) (%d %d %d %d)\n", k,j,i,l, curNode->element.k, curNode->element.j, curNode->element.i, curNode->element.l);
                    QUIT_PLUTO(1);
                }

                value = -1.5*dt*d->Vc[VX3][k][j][i]/(grid->x[2][k+1] - grid->x[2][k]) + 0.5*dt*u3/(grid->x[2][k+1] - grid->x[2][k]);
                curNode = addElement(curNode, value, k, j, i, l);
                if((value != value) || (value*0 != value*0)){
                    printf("36 turbulent value = NaN for (%d %d %d %d) (%d %d %d %d)\n", k,j,i,l, curNode->element.k, curNode->element.j, curNode->element.i, curNode->element.l);
                    QUIT_PLUTO(1);
                }

            }
#endif
#elif GEOMETRY == SPHERICAL
#if INCLUDE_IDIR
            if(u1 > 0){
                value = 1.5 *dt*u1/(grid->x[0][i] - grid->x[0][i-1]) - 0.5*dt*u1/(grid->x[0][i] - grid->x[0][i-1]);
                curNode = addElement(curNode, value, k, j, i, l);
                if((value != value) || (value*0 != value*0)){
                    printf("37 turbulent value = NaN for (%d %d %d %d) (%d %d %d %d)\n", k,j,i,l, curNode->element.k, curNode->element.j, curNode->element.i, curNode->element.l);
                    QUIT_PLUTO(1);
                }

                value = -1.5*dt*grid->x[0][i-1]*grid->x[0][i-1]*d->Vc[VX1][k][j][i-1]/(grid->x[0][i]*grid->x[0][i]*(grid->x[0][i] - grid->x[0][i-1])) + 0.5*dt*u1/(grid->x[0][i] - grid->x[0][i-1]);
                curNode = addElement(curNode, value, k, j, i-1, l);
                if((value != value) || (value*0 != value*0)){
                    printf("38 turbulent value = NaN for (%d %d %d %d) (%d %d %d %d)\n", k,j,i,l, curNode->element.k, curNode->element.j, curNode->element.i, curNode->element.l);
                    QUIT_PLUTO(1);
                }

            } else {
                value = 1.5 *dt*grid->x[0][i+1]*grid->x[0][i+1]**u1/(grid->x[0][i]*grid->x[0][i]*(grid->x[0][i+1] - grid->x[0][i])) - 0.5*dt*u1/(grid->x[0][i+1] - grid->x[0][i]);
                curNode = addElement(curNode, value, k, j, i+1, l);
                if((value != value) || (value*0 != value*0)){
                    printf("39 turbulent value = NaN for (%d %d %d %d) (%d %d %d %d)\n", k,j,i,l, curNode->element.k, curNode->element.j, curNode->element.i, curNode->element.l);
                    QUIT_PLUTO(1);
                }

                value = -1.5*dt*d->Vc[VX1][k][j][i]/(grid->x[0][i+1] - grid->x[0][i]) + 0.5*dt*u1/(grid->x[0][i+1] - grid->x[0][i]);
                curNode = addElement(curNode, value, k, j, i, l);
                if((value != value) || (value*0 != value*0)){
                    printf("40 turbulent value = NaN for (%d %d %d %d) (%d %d %d %d)\n", k,j,i,l, curNode->element.k, curNode->element.j, curNode->element.i, curNode->element.l);
                    QUIT_PLUTO(1);
                }

            }
#endif
#if INCLUDE_JDIR
            if(u2 > 0){
                value = 1.5 *dt*u2/(grid->x[0][i]*(grid->x[1][j] - grid->x[1][j-1])) - 0.5*dt*u2/(grid->x[0][i]*(grid->x[1][j] - grid->x[1][j-1]));
                curNode = addElement(curNode, value, k, j, i, l);
                if((value != value) || (value*0 != value*0)){
                    printf("41 turbulent value = NaN for (%d %d %d %d) (%d %d %d %d)\n", k,j,i,l, curNode->element.k, curNode->element.j, curNode->element.i, curNode->element.l);
                    QUIT_PLUTO(1);
                }

                value = -1.5*dt*sin(grid->x[1][j-1])*d->Vc[VX2][k][j-1][i]/(grid->x[0][i]*sin(grid->x[1][j])*(grid->x[1][j] - grid->x[1][j-1])) + 0.5*dt*u2/(grid->x[0][i]*(grid->x[1][j] - grid->x[1][j-1]));
                curNode = addElement(curNode, value, k, j-1, i, l);
                if((value != value) || (value*0 != value*0)){
                    printf("42 turbulent value = NaN for (%d %d %d %d) (%d %d %d %d)\n", k,j,i,l, curNode->element.k, curNode->element.j, curNode->element.i, curNode->element.l);
                    QUIT_PLUTO(1);
                }

            } else {
                value = 1.5 *dt*sin(grid->x[1][j+1])*d->Vc[VX2][k][j+1][i]/(grid->x[0][i]*sin(grid->x[1][j])*(grid->x[1][j+1] - grid->x[1][j])) - 0.5*dt*u2/(grid->x[0][i]*(grid->x[1][j+1] - grid->x[1][j]));
                curNode = addElement(curNode, value, k, j+1, i, l);
                if((value != value) || (value*0 != value*0)){
                    printf("43 turbulent value = NaN for (%d %d %d %d) (%d %d %d %d)\n", k,j,i,l, curNode->element.k, curNode->element.j, curNode->element.i, curNode->element.l);
                    QUIT_PLUTO(1);
                }

                value = -1.5*dt*u2/(grid->x[0][i]*(grid->x[1][j+1] - grid->x[1][j])) + 0.5*dt*u2/(grid->x[0][i]*(grid->x[1][j+1] - grid->x[1][j]));
                curNode = addElement(curNode, value, k, j, i, l);
                if((value != value) || (value*0 != value*0)){
                    printf("44 turbulent value = NaN for (%d %d %d %d) (%d %d %d %d)\n", k,j,i,l, curNode->element.k, curNode->element.j, curNode->element.i, curNode->element.l);
                    QUIT_PLUTO(1);
                }

            }
#endif
#if INCLUDE_KDIR
            printf("todo kdir\n");
            exit(0);
#endif
#endif
        }
    }

//printf("finish creating turbulent matrix\n");

    generalizedMinimalResidualMethod(grid, d->turbulent_matrix, d->turbulent_rightPart, d->Wt, d->turbulentBasis, NTURB, 1E-5, MAX_GMRES_ITERATIONS, 1);

    DOM_LOOP(k,j,i){
        for(int l = 0; l < NTURB; ++l){
            if(d->Wt[k][j][i][l] < 0){
                //printf("Wt < 0\n");
                d->Wt[k][j][i][l] = 0;
            }
            if(d->Wt[k][j][i][l] != d->Wt[k][j][i][l]){
                print("turbulent field = NaN\n");
                QUIT_PLUTO(1);
            }
        }
    }
//printf("finish turbulent field\n");
}
#endif
