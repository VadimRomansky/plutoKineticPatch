#include "pluto.h"
#include <stdbool.h>

#include "matrixElement.h"
#include "specialmath.h"

#if TURBULENT_FIELD == YES

void complexsqrt(double a, double b, double* c, double* d){
    double rho = sqrt(a*a + b*b);
    double phi = atan2(b, a);
    if(phi < 0){
        phi = phi + 2*CONST_PI;
    }
    rho = sqrt(rho);
    phi = phi/2;
    (*c) = rho*cos(phi);
    (*d) = rho*sin(phi);
}

void evaluateGrowthRate(Data *d, Grid *grid){
    int i, j, k;
    TOT_LOOP(k,j,i){
        for(int l = 0; l < NTURB; ++l){
            d->growthRate[k][j][i][l] = 0;
        }
    }

    DOM_LOOP(k,j,i){
        double B1 = d->Vc[BX1][k][j][i];
        double B2 = d->Vc[BX2][k][j][i];
        double B3 = d->Vc[BX3][k][j][i];
        double Bls = sqrt(B1*B1 + B2*B2 + B3*B3);

        double J1 = 0;
        double J2 = 0;
        double J3 = 0;

        J1 = J1 + d->Jkin1[k][j][i][0]*(d->k_turb[1] - d->k_turb[0]);
        J2 = J2 + d->Jkin2[k][j][i][0]*(d->k_turb[1] - d->k_turb[0]);
        J3 = J3 + d->Jkin3[k][j][i][0]*(d->k_turb[1] - d->k_turb[0]);
        for(int m = 1; m < NMOMENTUM; ++m){
            J1 = J1 + d->Jkin1[k][j][i][m]*(d->k_turb[m] - d->k_turb[m-1]);
            J2 = J2 + d->Jkin2[k][j][i][m]*(d->k_turb[m] - d->k_turb[m-1]);
            J3 = J3 + d->Jkin3[k][j][i][m]*(d->k_turb[m] - d->k_turb[m-1]);
        }

        double J = sqrt(J1*J1 + J2*J2 + J3*J3)*(PARTICLES_KIN_E_MC*PARTICLES_KIN_MASS*PARTICLES_KIN_C);

        for(int l = 1; l < NTURB; ++l){
            Bls = sqrt(4*CONST_PI*d->Wt[k][j][i][l]*(d->k_turb[l]-d->k_turb[l-1]) + Bls*Bls);
            CheckNanOrInfinity(Bls, "Bls = NaN");
            //todo check!
            double kc = fabs(4*CONST_PI*J/((CONST_c/UNIT_VELOCITY)*Bls));
            double Va = Bls/sqrt(4*CONST_PI*d->Vc[RHO][k][j][i]);

            double A1re = 0;
            double A1im = 0;
            double A2re = 0;
            double A2im = 0;
            for(int m = 0; m < NMOMENTUM; ++m){
                double z = d->k_turb[l]*d->p_grid[m]/(PARTICLES_KIN_E_MC*Bls);
                CheckNanOrInfinity(z, "z = NaN");
                double sigma1re = 0;
                double sigma1im = 0;
                double sigma2re = 0;
                double sigma2im = 0;
                if( fabs(z - 1.0) < 0.00001){
                    sigma1re = 3.0/2.0;
                } else if(z > 1) {
                    sigma1re = (1.5/(z*z)) + 0.75*(1 - 1/((z*z)))*log(fabs((z+1)/(z-1)))/z;
                } else if(0.01 < z) {
                    sigma1re = (1.5/(z*z)) + 0.75*(1 - 1/((z*z)))*log(fabs((z+1)/(z-1)))/z;
                }  else {
                    sigma1re = 1 + 0.2*z*z;
                }
                //sigma1 = 1;

                if(z > 1){
                    sigma1im = sigma1im -(3*CONST_PI/(4*z))*(1 - 1/(z*z));
                }
                CheckNanOrInfinity(sigma1re, "sigma = NaN");
                CheckNanOrInfinity(sigma1im, "sigma = NaN");
                sigma2re = sigma1re;
                sigma2im = -sigma1im;

                double dp = d->p_grid[1] - d->p_grid[0];
                if(m > 0){
                    dp = d->p_grid[m] - d->p_grid[m-1];
                }

                double crflux = sqrt(d->Jkin1[k][j][i][m]*d->Jkin1[k][j][i][m]
                                     + d->Jkin2[k][j][i][m]*d->Jkin2[k][j][i][m]
                                     + d->Jkin3[k][j][i][m]*d->Jkin3[k][j][i][m])*dp*(PARTICLES_KIN_E_MC*PARTICLES_KIN_MASS*PARTICLES_KIN_C);

                A1re = A1re + sigma1re*crflux;
                A1im = A1im + sigma1im*crflux;

                CheckNanOrInfinity(A1re, "A = NaN");
                CheckNanOrInfinity(A1im, "A = NaN");
                A2re = A2re + sigma2re*crflux;
                A2im = A2im + sigma2im*crflux;
            }

            //Complex complex1 = Complex(1);
            //Complex delta = complex1 - (A1/J - 1)*(kc/kgrid[k]);
            double deltare = 1.0 - (A1re/J - 1)*(kc/d->k_turb[l]);
            double deltaim = 0.0 - (A1im/J)*(kc/d->k_turb[l]);

            //Complex b1 = (complex1 - (A1/J - 1)*(kc/d->k_turb[l]))*POW2(d->k_turb[l]*Va);
            //Complex b2 = (complex1 + (A2/J - 1)*(kc/kgrid[k]))*POW2(kgrid[k]*Va);

            double b1re = (1.0 - (A1re/J - 1)*(kc/d->k_turb[l]))*POW2(d->k_turb[l]*Va);
            double b1im = (0.0 - (A1im/J)*(kc/d->k_turb[l]))*POW2(d->k_turb[l]*Va);

            double b2re = (1.0 - (A2re/J - 1)*(kc/d->k_turb[l]))*POW2(d->k_turb[l]*Va);
            double b2im = (0.0 - (A2im/J)*(kc/d->k_turb[l]))*POW2(d->k_turb[l]*Va);
            double alpha = 1.5;
            //double alpha = 0;
            //Complex d1 = Complex(0, -1)*((A1*0.5/J) + 1.5)*(kgrid[k]*kc)*alpha/(4*pi*middleDensity[i]);
            //Complex d2 = Complex(0, 1)*((A2*0.5/J) + 1.5)*(kgrid[k]*kc)*alpha/(4*pi*middleDensity[i]);

            double d1im = -((A1re*0.5/J) + 1.5)*(d->k_turb[l]*kc)*alpha/(4*CONST_PI*d->Vc[RHO][k][j][i]);
            double d1re = (A1im*0.5/J)*(d->k_turb[l]*kc)*alpha/(4*CONST_PI*d->Vc[RHO][k][j][i]);

            double d2im = ((A2re*0.5/J) + 1.5)*(d->k_turb[l]*kc)*alpha/(4*CONST_PI*d->Vc[RHO][k][j][i]);
            double d2re = -(A2im*0.5/J)*(d->k_turb[l]*kc)*alpha/(4*CONST_PI*d->Vc[RHO][k][j][i]);

            //Complex G1p = (csqrt(d1*d1 +b1*4) - d1)/2;
            //Complex G1m = (csqrt(d1*d1 +b1*4) + d1)/(-2);

            double D1re, D1im;
            complexsqrt(d1re*d1re-d1im*d1im +4*b1re, 2*d1re*d1im + 4*b1im, &D1re, &D1im);

            //Complex G2p = (csqrt(d2*d2 +b2*4) - d2)/2;
            //Complex G2m = (csqrt(d2*d2 +b2*4) + d2)/(-2);

            double D2re, D2im;
            complexsqrt(d2re*d2re-d2im*d2im +4*b2re, 2*d2re*d2im + 4*b2im, &D2re, &D2im);

            double rate1 = (D1im - d1im)/2.0;
            double rate3 = (-D1im - d1im)/2.0;
            double rate2 = (D2im - d2im)/2.0;
            double rate4 = (-D2im - d2im)/2.0;


            double rate = MAX(MAX(rate1, rate2), MAX(rate3, rate4));
            if(rate > 0){
                rate *= 2;
            } else {
                rate = 0;
            }

            if((rate != rate) || (0*rate != 0*rate)){
                printf("rate = NaN\n");
                exit(0);
            }

            d->growthRate[k][j][i][l] = rate;
            CheckNanOrInfinity(rate, "growth_rate = NaN");
            /*if(rate > 0){
                printf("%g\n", rate);
            }*/
        }
    }
}

void AdvanceTurbulentField(Data *d, timeStep *Dts, double dt, Grid *grid){
    int k,j,i;
    //printf("trubulent field\n");
    //return;

    bool periodicX = (grid->lbound[0] == PERIODIC);
    bool periodicY = (grid->lbound[1] == PERIODIC);
    bool periodicZ = (grid->lbound[2] == PERIODIC);

    double inv_dt, inv_dt_new;
    inv_dt = 1.e-18;

    Dts->invDt_magnetic = inv_dt;

    evaluateGrowthRate(d, grid);

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

            double G = d->growthRate[k][j][i][l];

            inv_dt_new = G;

            Dts->invDt_magnetic = MAX(Dts->invDt_magnetic, inv_dt_new);

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

    generalizedMinimalResidualMethod1(grid, d->turbulent_matrix, d->turbulent_rightPart, d->Wt, d->turbulentBasis, NTURB, 1E-5, MAX_GMRES_ITERATIONS, 1, periodicX, periodicY, periodicZ);

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
