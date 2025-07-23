/* ///////////////////////////////////////////////////////////////////// */
/*! 
  \file  
  \brief Update CR particles (with feedback) using Boris scheme.
 
  Boris pusher for updating Cosmic ray particles.
  \code
   dt = dt_g/Nsub
    for n = 0,Nsub-1 {
     Compute Fcr
     Loop on particles{
       x(n+1/2) = x(n)     + (dt/2)*v(n)     (drift)
       v(n+1)   = v(n)     +  dt*a(n+1/2)    (kick+rotate+kick)
       x(n+1)   = x(n+1/2) + (dt/2)*v(n+1)   (drift)
        save dL += dt*(v(n) + v(n+1))/2
       Compute time step:
         dL < Nz*dx --> dt < Nz*dx/[v(0)/2 + v(1) + .. V(Nsub-1)/2)]
     }
   }
  \endcode
 
  Time step restriction is computed by requiring that no particle
  travels more than \c Nmax = \c PARTICLES_CR_NCELL_MAX zones and that
  the Larmor scale is resolved with more than 1 cycle:
  \f[
   \left\{
     \begin{array}{lcl}
       \Delta s_d &=& \DS\Delta t\max_{n,p}\left(
                      \frac{|v_{p,d}^n + v_{p,d}^{n+1}|}{2}\right)
                   < \epsilon_s N_{\max}\Delta x_d \\ \noalign{\medskip}
       \Delta t  &<& \epsilon_L \Omega_L^{-1}
     \end{array}
   \right.
  \f]
  where the maximum extends to all sub-steps (\c n) and particles
  (p), \f$ \Omega_L = qB_\perp/(\gamma m c) \f$ is the Larmor frequency while
  \f$\epsilon_s (\sim 0.9)\f$ and \f$\epsilon_L\sim 0.3\f$ are safety factors
  while 
  \f[
   \vec{B}^2_\perp = \vec{B}^2 - \frac{(\vec{v}\cdot\vec{B})^2}
                                      {\vec{v}\cdot\vec{v}}                                
  \f]
  
  \authors A. Mignone (mignone@to.infn.it)\n

  \b References
    - "MAGNETOHYDRODYNAMIC-PARTICLE-IN-CELL METHOD FOR COUPLING COSMIC RAYS 
       WITH A THERMAL PLASMA: APPLICATION TO NON-RELATIVISTIC SHOCKS"\n
       Bai et al., ApJ (2015) 809, 55

  \date   May 09, 2020
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"
#include <stdbool.h>
#include "specialmath.h"
#include "matrixElement.h"

#ifdef PARALLEL
#include "al_hidden.h"
#endif

#ifdef PARALLEL
extern struct SZ *sz_stack[AL_MAX_ARRAYS];
extern int stack_ptr[AL_MAX_ARRAYS];
#endif

#if TURBULENT_FIELD == YES
double evaluateTurbulentField(Data* data, int i, int j, int k, double* B0, double momentum_q){
    double B20 =DOT_PRODUCT(B0, B0);
    int kindex = 0;

    double B2turb = 0;
    double rg = momentum_q * PARTICLES_KIN_C/sqrt(B20 + B2turb);

    for(int l = 0; l < NTURB; ++l){
        if(rg > 2*CONST_PI/data->k_turb[l]){
            return sqrt(B20 + B2turb);
        }
        double dk;
        if(l == 0){
            dk = data->k_turb[l+1] - data->k_turb[l];
        } else {
            dk = data->k_turb[l] - data->k_turb[l-1];
        }
        B2turb += 8*CONST_PI*data->Wt[k][j][i][l]*dk;

        rg = momentum_q*PARTICLES_KIN_C/sqrt(B20 + B2turb);
    }
    return sqrt(B20 + B2turb);
}
#endif

double evaluateDiffusionCoefficient(Data* data, int i, int j, int k, double u){
    const double coef = (PARTICLES_KIN_C / (PARTICLES_KIN_E_MC))*(CONST_c/UNIT_VELOCITY)/3.0;
    double B[3];
    B[IDIR] = data->Vc[BX1][k][j][i];
    B[JDIR] = data->Vc[BX2][k][j][i];
    B[KDIR] = data->Vc[BX3][k][j][i];
#if TURBULENT_FIELD == YES
    double Bmag = evaluateTurbulentField(data,i,j,k, B, u/PARTICLES_KIN_E_MC);
#else
    double Bmag = sqrt(B[IDIR]*B[IDIR] + B[JDIR]*B[JDIR] + B[KDIR]*B[KDIR]);
#endif
    double D = coef*u/Bmag;

    if ((D != D) || (D*0 != D*0)){
        printf("D = NaN or Infinity\n");
        printf("u = %g, B = %g\n", u, Bmag);
        QUIT_PLUTO(1);
    }

    //D = 1E26/(UNIT_LENGTH*UNIT_VELOCITY);
    //D = coef/Bmag;
    D = D;

    return D;
}

/* ********************************************************************* */
void Particles_KIN_Update(Data *data, timeStep *Dts, double dt, Grid *grid)
/*!
 * Advance particle by a step dt.
 * 
 * \param [in,out]  d     Data structure (contains paeticles list)
 * \param [in,out]  Dts   timeStep structure
 * \param [in]      dt    time time increment
 * \param [in]      grid  pointer to Grid structure
 *********************************************************************** */
{
    int i, j, k, dir;
    int kcycle;
    long int dsize = NX1_TOT * NX2_TOT * NX3_TOT;
    double inv_dt, inv_dt_new;
    inv_dt = 1.e-18;

    Dts->invDt_advection = inv_dt;
    Dts->invDt_acceleration = inv_dt;
    Dts->invDt_diffusion = inv_dt;
    double Dr_uD = 1.e-18;
    Dts->Dr_uD = Dr_uD;


    //crank-nickelson
    double factor = 0.5;
    //factor = 1.0;

    FlagShock(data, grid);
    //double err = ConsToPrim3D()

    bool periodicX = (grid->lbound[0] == PERIODIC);
    bool periodicY = (grid->lbound[1] == PERIODIC);
    bool periodicZ = (grid->lbound[2] == PERIODIC);

    TOT_LOOP(k,j,i){
        for(int l = 0; l < NMOMENTUM; ++l){
                MatrixElementNode* curNode = data->matrix[k][j][i][l];
                while(curNode != NULL){
                    MatrixElementNode* tempNode = curNode;
                    curNode = curNode->next;
                    free(tempNode);
                }
                curNode = (MatrixElementNode*) malloc(sizeof(MatrixElementNode));
                curNode->element = createMatrixElement(1.0, k,j,i,l);
                curNode->next = NULL;
                curNode->prev = NULL;
                data->rightPart[k][j][i][l] = 0.0;
                data->matrix[k][j][i][l] = curNode;

                curNode = data->rightPartMatrix[k][j][i][l];
                while(curNode != NULL){
                    MatrixElementNode* tempNode = curNode;
                    curNode = curNode->next;
                    free(tempNode);
                }
                curNode = (MatrixElementNode*) malloc(sizeof(MatrixElementNode));
                curNode->element = createMatrixElement(1.0, k,j,i,l);
                curNode->next = NULL;
                curNode->prev = NULL;
                data->rightPartMatrix[k][j][i][l] = curNode;

#if INCLUDE_IDIR
                data->ax[k][j][i][l] = 0;
                data->bx[k][j][i][l] = 1.0;
                data->cx[k][j][i][l] = 0;
#endif
#if INCLUDE_JDIR
                data->ay[k][j][i][l] = 0;
                data->by[k][j][i][l] = 1.0;
                data->cy[k][j][i][l] = 0;
#endif
#if INCLUDE_KDIR
                data->az[k][j][i][l] = 0;
                data->bz[k][j][i][l] = 1.0;
                data->cz[k][j][i][l] = 0;
#endif
                data->ap[k][j][i][l] = 0;
                data->bp[k][j][i][l] = 1.0;
                data->cp[k][j][i][l] = 0;
        }
    }

    DOM_LOOP(k,j,i){
        double divu = 0;

        int maxNU = NMOMENTUM - 1;

#if INCLUDE_IDIR
        if(grid->lbound[0] != 0){
            if(grid->lbound[0] != PERIODIC){
            if(i == IBEG){
                if(grid->lbound[0] != OUTFLOW){
                    for(int l = 0; l < maxNU; ++l){
                        MatrixElementNode* curNode = data->matrix[k][j][i][l];
                        curNode = addElement(curNode, -1.0, k,j,i+1,l);
                        data->cx[k][j][i][l] = -1.0;
                    }
                }
                continue;
            }
            }
        }

        if(grid->rbound[0] != 0){
            if(grid->rbound[0] != PERIODIC){
            if(i == IEND){
                if(grid->rbound[0] != OUTFLOW){
                    for(int l = 0; l < maxNU; ++l){
                        MatrixElementNode* curNode = data->matrix[k][j][i][l];
                        curNode = addElement(curNode, -1.0, k,j,i-1,l);
                        data->ax[k][j][i][l] = -1.0;
                    }
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
                    if(grid->lbound[1] != OUTFLOW){
                        for(int l = 0; l < maxNU; ++l){
                            MatrixElementNode* curNode = data->matrix[k][j][i][l];
                            curNode = addElement(curNode, -1.0, k,j + 1,i,l);
                            data->cy[k][j][i][l] = -1.0;
                        }
                    }
                    continue;
                }
            }
        }

        if(grid->rbound[1] != 0){
            if(grid->rbound[1] != PERIODIC){
                if(j == JEND){
                    if(grid->rbound[1] != OUTFLOW){
                        for(int l = 0; l < maxNU; ++l){
                            MatrixElementNode* curNode = data->matrix[k][j][i][l];
                            curNode = addElement(curNode, -1.0, k,j - 1,i,l);
                            data->ay[k][j][i][l] = -1.0;
                        }
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
                    if(grid->lbound[2] != OUTFLOW){
                        for(int l = 0; l < maxNU; ++l){
                            MatrixElementNode* curNode = data->matrix[k][j][i][l];
                            curNode = addElement(curNode, -1.0, k + 1,j,i,l);
                            data->cz[k][j][i][l] = -1.0;
                        }
                    }
                    continue;
                }
            }
        }

        if(grid->rbound[2] != 0){
            if(grid->rbound[2] != PERIODIC){
                if(k == KEND){
                    if(grid->rbound[2] != OUTFLOW){
                        for(int l = 0; l < maxNU; ++l){
                            MatrixElementNode* curNode = data->matrix[k][j][i][l];
                            curNode = addElement(curNode, -1.0, k - 1,j,i,l);
                            data->az[k][j][i][l] = -1.0;
                        }
                    }
                    continue;
                }
            }
        }

#endif

#if GEOMETRY == CARTESIAN
#if INCLUDE_IDIR
        //divu += (data->Vc[VX1][k][j][i+1] - data->Vc[VX1][k][j][i-1])/(grid->x[0][i+1] - grid->x[0][i-1]);
        if(data->Vc[VX1][k][j][i] > 0){
            divu += (data->Vc[VX1][k][j][i] - data->Vc[VX1][k][j][i-1])/(grid->x[0][i] - grid->x[0][i-1]);
        } else {
            divu += (data->Vc[VX1][k][j][i+1] - data->Vc[VX1][k][j][i])/(grid->x[0][i+1] - grid->x[0][i]);
        }
#endif
#if INCLUDE_JDIR
        divu += (data->Vc[VX2][k][j+1][i] - data->Vc[VX2][k][j-1][i])/(grid->x[1][j+1] - grid->x[1][j-1]);
#endif
#if INCLUDE_KDIR
        divu += (data->Vc[VX3][k+1][j][i] - data->Vc[VX3][k-1][j][i])/(grid->x[2][k+1] - grid->x[2][k-1]);
#endif
#elif GEOMETRY == CYLINDRICAL
#if INCLUDE_IDIR
        divu += (grid->x[0][i+1]*data->Vc[VX1][k][j][i+1] - grid->x[0][i-1]*data->Vc[VX1][k][j][i-1])/(grid->x[0][i]*(grid->x[0][i+1] - grid->x[0][i-1]));
        if((divu != divu) || (0*divu != 0*divu)){
            printf("divu = NaN\n");
            QUIT_PLUTO(1);
        }
#endif
#if INCLUDE_JDIR
        divu += (data->Vc[VX2][k][j+1][i] - data->Vc[VX2][k][j-1][i])/(grid->x[1][j+1] - grid->x[1][j-1]);
        if((divu != divu) || (0*divu != 0*divu)){
            printf("divu = NaN\n");
            QUIT_PLUTO(1);
        }
#endif
#if INCLUDE_KDIR
        divu += (data->Vc[VX3][k+1][j][i] - data->Vc[VX3][k-1][j][i])/(grid->x[0][i]*(grid->x[2][k+1] - grid->x[2][k-1]));
        if((divu != divu) || (0*divu != 0*divu)){
            printf("divu = NaN\n");
            QUIT_PLUTO(1);
        }
#endif
#elif GEOMETRY == POLAR
#if INCLUDE_IDIR
        divu += (grid->x[0][i+1]*data->Vc[VX1][k][j][i+1] - grid->x[0][i-1]*data->Vc[VX1][k][j][i-1])/(grid->x[0][i]*(grid->x[0][i+1] - grid[0][i-1]));
        if((divu != divu) || (0*divu != 0*divu)){
            printf("divu = NaN\n");
            QUIT_PLUTO(1);
        }

#endif
#if INCLUDE_JDIR
        divu += (data->Vc[VX2][k][j+1][i] - data->Vc[VX2][k][j-1][i])/(grid->x[1][j+1] - grid->x[1][j-1]);
        if((divu != divu) || (0*divu != 0*divu)){
            printf("divu = NaN\n");
            QUIT_PLUTO(1);
        }

#endif
#if INCLUDE_KDIR
        divu += (data->Vc[VX3][k+1][j][i] - data->Vc[VX3][k-1][j][i])/(grid->x[0][i]*(grid->x[2][k+1] - grid->x[2][k-1]));
        if((divu != divu) || (0*divu != 0*divu)){
            printf("divu = NaN\n");
            QUIT_PLUTO(1);
        }

#endif
#elif GEOMETRY == SPHERICAL
#if INCLUDE_IDIR
        divu += (grid->x[0][i+1]*grid->x[0][i+1]*data->Vc[VX1][k][j][i+1] - grid->x[0][i-1]*grid->x[0][i-1]*data->Vc[VX1][k][j][i-1])/(grid->x[0][i]*grid->x[0][i]*(grid->x[0][i+1] - grid[0][i-1]));
        if((divu != divu) || (0*divu != 0*divu)){
            printf("divu = NaN\n");
            QUIT_PLUTO(1);
        }

#endif
#if INCLUDE_JDIR
        double sintheta = sin(grid->x[1][j]);
        divu += (sin(grid->x[1][j+1])*data->Vc[VX2][k][j+1][i] - sin(grid->x[1][j-1])*data->Vc[VX2][k][j-1][i])/(grid->x[0][i]*sintheta*(grid->x[1][j+1] - grid->x[1][j-1]));
        if((divu != divu) || (0*divu != 0*divu)){
            printf("divu = NaN\n");
            QUIT_PLUTO(1);
        }

#endif
#if INCLUDE_KDIR
        divu += (data->Vc[VX3][k+1][j][i] - data->Vc[VX3][k-1][j][i])/(grid->x[0][i]*sintheta*(grid->x[2][k+1] - grid->x[2][k-1]));
        if((divu != divu) || (0*divu != 0*divu)){
            printf("divu = NaN\n");
            QUIT_PLUTO(1);
        }

#endif
#else
        printf("wrong geometry\n");
        QUIT_PLUTO(1);
#endif

        ////why 100?
        inv_dt_new = fabs(3*divu*data->p_grid[0]/(data->p_grid[1] - data->p_grid[0]))/PARTICLES_KIN_EPS;
        Dts->invDt_acceleration = MAX(Dts->invDt_acceleration, inv_dt_new);
       
        double dx1 = 1E100;
        double dx2 = 1E100;
        double dx3 = 1E100;
#if GEOMETRY == CARTESIAN
        dx1 = grid->x[0][i+1] - grid->x[0][i];
#if INCLUDE_JDIR
        dx2 = grid->x[1][j+1] - grid->x[1][j];
#endif
#if INCLUDE_KDIR
        dx3 = grid->x[2][k+1] - grid->x[2][k];
#endif
#elif GEOMETRY == CYLINDRICAL
        dx1 = grid->x[0][i+1] - grid->x[0][i];
#if INCLUDE_JDIR
        dx2 = grid->x[1][j+1] - grid->x[1][j];
#endif
#if INCLUDE_KDIR
        dx3 = grid->x[0][i]*(grid->x[2][k+1] - grid->x[2][k]);
#endif

#elif GEOMETRY == POLAR
        dx1 = grid->x[0][i+1] - grid->x[0][i];
#if INCLUDE_JDIR
        dx2 = grid->x[0][i]*(grid->x[1][j+1] - grid->x[1][j]);
#endif
#if INCLUDE_KDIR
        dx3 = grid->x[2][k+1] - grid->x[2][k];
#endif

#elif GEOMETRY == SPHERICAL
        dx1 = grid->x[0][i+1] - grid->x[0][i];
#if INCLUDE_JDIR
        dx2 = grid->x[0][i]*(grid->x[1][j+1] - grid->x[1][j]);
#endif
#if INCLUDE_KDIR
        dx3 = grid->x[0][i]*sin(grid->x[1][j])*(grid->x[2][k+1] - grid->x[2][k]);
#endif

#endif

        inv_dt_new = fabs(2*data->Vc[VX1][k][j][i]/dx1)/PARTICLES_KIN_EPS;
        Dts->invDt_advection = MAX(Dts->invDt_advection, inv_dt_new);
#if INCLUDE_JDIR
        inv_dt_new = fabs(2*data->Vc[VX2][k][j][i]/dx2)/PARTICLES_KIN_EPS;
        Dts->invDt_advection = MAX(Dts->invDt_advection, inv_dt_new);
#endif
#if INCLUDE_KDIR
        inv_dt_new = fabs(2*data->Vc[VX3][k][j][i]/dx3)/PARTICLES_KIN_EPS;
        Dts->invDt_advection = MAX(Dts->invDt_advection, inv_dt_new);
#endif


        for(int l = 0; l < maxNU; ++l){
            if(l > 0){
                inv_dt_new = fabs(3*divu*data->p_grid[l]/(data->p_grid[l] - data->p_grid[l-1]))/PARTICLES_KIN_EPS;
                Dts->invDt_acceleration = MAX(Dts->invDt_acceleration, inv_dt_new);
            }

            data->rightPart[k][j][i][l] = data->Fkin[k][j][i][l];

            MatrixElementNode* curNode = data->matrix[k][j][i][l];
            MatrixElementNode* curNode2 = data->rightPartMatrix[k][j][i][l];
            /*while(curNode != NULL){
                MatrixElementNode* tempNode = curNode;
                curNode = curNode->next;
                free(tempNode);
            }*/

            double D = evaluateDiffusionCoefficient(data, i,j,k, data->p_grid[l]);
            double Dr_uD_new = 0.0;
#if INCLUDE_IDIR
            double Dleft =  evaluateDiffusionCoefficient(data, i-1,j,k, data->p_grid[l]);

            double Dright =  evaluateDiffusionCoefficient(data, i+1,j,k, data->p_grid[l]);
            Dr_uD_new = Dr_uD_new + fabs(dx1*data->Vc[VX1][k][j][i])/D;
#endif

#if INCLUDE_JDIR
            double Dbottom =  evaluateDiffusionCoefficient(data, i,j-1,k, data->p_grid[l]);

            double Dtop = evaluateDiffusionCoefficient(data, i,j+1,k, data->p_grid[l]);
            Dr_uD_new = Dr_uD_new + fabs(dx2*data->Vc[VX2][k][j][i])/D;
#endif

#if INCLUDE_KDIR
            double Dback = evaluateDiffusionCoefficient(data, i,j,k-1, data->p_grid[l]);

            double Dfront = evaluateDiffusionCoefficient(data, i,j,k+1, data->p_grid[l]);
            Dr_uD_new = Dr_uD_new + fabs(dx3*data->Vc[VX3][k][j][i])/D;
#endif
            Dts->Dr_uD = MAX(Dts->Dr_uD, Dr_uD_new);
            inv_dt_new = fabs(2*4*D/(dx1*dx1))/PARTICLES_KIN_EPS;
            Dts->invDt_diffusion = MAX(Dts->invDt_diffusion, inv_dt_new);
#if INCLUDE_JDIR
            inv_dt_new = fabs(2*4*D/(dx2*dx2))/PARTICLES_KIN_EPS;
            Dts->invDt_diffusion= MAX(Dts->invDt_diffusion, inv_dt_new);
#endif
#if INCLUDE_KDIR
            inv_dt_new = fabs(2*4*D/(dx3*dx3))/PARTICLES_KIN_EPS;
            Dts->invDt_diffusion = MAX(Dts->invDt_diffusion, inv_dt_new);
#endif

            double dp = 0;
            if(divu > 0){
                if(l == NMOMENTUM - 1){
                    dp = data->p_grid[NMOMENTUM - 1] - data->p_grid[NMOMENTUM - 2];
                } else {
                    dp = data->p_grid[l+1] - data->p_grid[l];
                }
            } else {
                if(l == 0){
                    dp = data->p_grid[1] - data->p_grid[0];
                } else {
                    dp = data->p_grid[l] - data->p_grid[l-1];
                }
            }

            double value = 0;

            if(divu >= 0){
                value = factor*dt*divu*data->p_grid[l]/(3.0*dp);
                curNode = addElement(curNode, value, k,j,i,l);
                data->bp[k][j][i][l] += value;
                curNode2 = addElement(curNode2, -value, k, j, i, l);
                if((value != value) || (value*0 != value*0)){
                    printf("1 value = NaN for (%d %d %d %d) (%d %d %d %d)\n", k,j,i,l, curNode->element.k, curNode->element.j, curNode->element.i, curNode->element.l);
                    QUIT_PLUTO(1);
                }
                if(l < NMOMENTUM - 1){
                    value = -factor*dt*divu*data->p_grid[l]/(3.0*dp);
                    curNode = addElement(curNode, value, k,j,i,l+1);
                    data->cp[k][j][i][l] += value;
                    curNode2 = addElement(curNode2, -value, k, j, i, l+1);
                    if((value != value) || (value*0 != value*0)){
                        printf("2 value = NaN for (%d %d %d %d) (%d %d %d %d)\n", k,j,i,l, curNode->element.k, curNode->element.j, curNode->element.i, curNode->element.l);
                        QUIT_PLUTO(1);
                    }
                }
            } else {
                value = -factor*dt*divu*data->p_grid[l]/(3.0*dp);
                curNode = addElement(curNode, value, k,j,i,l);
                data->bp[k][j][i][l] += value;
                curNode2 = addElement(curNode2, -value, k, j, i, l);
                if((value != value) || (value*0 != value*0)){
                    printf("3 value = NaN for (%d %d %d %d) (%d %d %d %d)\n", k,j,i,l, curNode->element.k, curNode->element.j, curNode->element.i, curNode->element.l);
                    QUIT_PLUTO(1);
                }
                if( l > 0){
                    value = factor*dt*divu*data->p_grid[l]/(3.0*dp);
                    curNode = addElement(curNode, value, k,j,i,l-1);
                    data->ap[k][j][i][l] += value;
                    curNode2 = addElement(curNode2, -value, k, j, i, l-1);
                    if((value != value) || (value*0 != value*0)){
                        printf("4 value = NaN for (%d %d %d %d) (%d %d %d %d)\n", k,j,i,l, curNode->element.k, curNode->element.j, curNode->element.i, curNode->element.l);
                        QUIT_PLUTO(1);
                    }
                }
            }

#if GEOMETRY == CARTESIAN
            value = factor*dt*(((Dright + D)/(grid->x[0][i+1] - grid->x[0][i]))+((D+Dleft)/(grid->x[0][i] - grid->x[0][i-1])))/(grid->x[0][i+1] - grid->x[0][i-1]);
            curNode = addElement(curNode, value, k,j,i,l);
            data->bx[k][j][i][l] += value;
            curNode2 = addElement(curNode2, -value, k,j,i,l);
                if((value != value) || (value*0 != value*0)){
                    printf("5 value = NaN for (%d %d %d %d) (%d %d %d %d)\n", k,j,i,l, curNode->element.k, curNode->element.j, curNode->element.i, curNode->element.l);
                    QUIT_PLUTO(1);
                }


            value = -factor*dt*((Dright + D)/(grid->x[0][i+1] - grid->x[0][i]))/(grid->x[0][i+1] - grid->x[0][i-1]);
            curNode = addElement(curNode, value, k, j, i+1, l);
            data->cx[k][j][i][l] += value;
            curNode2 = addElement(curNode2, -value, k, j, i+1, l);
                if((value != value) || (value*0 != value*0)){
                    printf("6 value = NaN for (%d %d %d %d) (%d %d %d %d)\n", k,j,i,l, curNode->element.k, curNode->element.j, curNode->element.i, curNode->element.l);
                    QUIT_PLUTO(1);
                }


            value = -factor*dt*((D+Dleft)/(grid->x[0][i] - grid->x[0][i-1]))/(grid->x[0][i+1] - grid->x[0][i-1]);
            curNode = addElement(curNode, value, k, j, i-1, l);
            data->ax[k][j][i][l] += value;
            curNode2 = addElement(curNode2, -value, k, j, i-1, l);
                if((value != value) || (value*0 != value*0)){
                    printf("7 value = NaN for (%d %d %d %d) (%d %d %d %d)\n", k,j,i,l, curNode->element.k, curNode->element.j, curNode->element.i, curNode->element.l);
                    QUIT_PLUTO(1);
                }


            if(data->Vc[VX1][k][j][i] > 0){
                value = factor*dt*data->Vc[VX1][k][j][i]/(grid->x[0][i] - grid->x[0][i-1]);
                curNode = addElement(curNode, value, k,j,i,l);
                data->bx[k][j][i][l] += value;
                curNode2 = addElement(curNode2, -value, k,j,i,l);
                if((value != value) || (value*0 != value*0)){
                    printf("8 value = NaN for (%d %d %d %d) (%d %d %d %d)\n", k,j,i,l, curNode->element.k, curNode->element.j, curNode->element.i, curNode->element.l);
                    QUIT_PLUTO(1);
                }

                value = -factor*dt*data->Vc[VX1][k][j][i-1]/(grid->x[0][i] - grid->x[0][i-1]);
                curNode = addElement(curNode, value, k,j,i-1,l);
                data->ax[k][j][i][l] += value;
                curNode2 = addElement(curNode2, -value, k,j,i-1,l);
                if((value != value) || (value*0 != value*0)){
                    printf("9 value = NaN for (%d %d %d %d) (%d %d %d %d)\n", k,j,i,l, curNode->element.k, curNode->element.j, curNode->element.i, curNode->element.l);
                    QUIT_PLUTO(1);
                }

            } else {
                value = -factor*dt*data->Vc[VX1][k][j][i]/(grid->x[0][i+1] - grid->x[0][i]);
                curNode = addElement(curNode, value, k,j,i,l);
                data->bx[k][j][i][l] += value;
                curNode2 = addElement(curNode2, -value, k,j,i,l);
                if((value != value) || (value*0 != value*0)){
                    printf("10 value = NaN for (%d %d %d %d) (%d %d %d %d)\n", k,j,i,l, curNode->element.k, curNode->element.j, curNode->element.i, curNode->element.l);
                    QUIT_PLUTO(1);
                }

                value = factor*dt*data->Vc[VX1][k][j][i+1]/(grid->x[0][i+1] - grid->x[0][i]);
                curNode = addElement(curNode, value, k,j,i+1,l);
                data->cx[k][j][i][l] += value;
                curNode2 = addElement(curNode2, -value, k,j,i+1,l);
                if((value != value) || (value*0 != value*0)){
                    printf("11 value = NaN for (%d %d %d %d) (%d %d %d %d)\n", k,j,i,l, curNode->element.k, curNode->element.j, curNode->element.i, curNode->element.l);
                    QUIT_PLUTO(1);
                }

            }
#if INCLUDE_JDIR
            value = factor*dt*(((Dtop + D)/(grid->x[1][j+1] - grid->x[1][j]))+((D+Dbottom)/(grid->x[1][j] - grid->x[1][j-1])))/(grid->x[1][j+1] - grid->x[1][j-1]);
            curNode = addElement(curNode, value, k,j,i,l);
            data->by[k][j][i][l] += value;
            curNode2 = addElement(curNode2, -value, k,j,i,l);
                if((value != value) || (value*0 != value*0)){
                    printf("12 value = NaN for (%d %d %d %d) (%d %d %d %d)\n", k,j,i,l, curNode->element.k, curNode->element.j, curNode->element.i, curNode->element.l);
                    QUIT_PLUTO(1);
                }


            value = -factor*dt*((Dtop + D)/(grid->x[1][j+1] - grid->x[1][j]))/(grid->x[1][j+1] - grid->x[1][j-1]);
            curNode = addElement(curNode, value, k, j+1, i, l);
            data->cy[k][j][i][l] += value;
            curNode2 = addElement(curNode2, -value, k, j+1, i, l);
                if((value != value) || (value*0 != value*0)){
                    printf("13 value = NaN for (%d %d %d %d) (%d %d %d %d)\n", k,j,i,l, curNode->element.k, curNode->element.j, curNode->element.i, curNode->element.l);
                    QUIT_PLUTO(1);
                }


            value = -factor*dt*((D+Dbottom)/(grid->x[1][j] - grid->x[1][j-1]))/(grid->x[1][j+1] - grid->x[1][j-1]);
            curNode = addElement(curNode, value, k, j-1, i, l);
            data->ay[k][j][i][l] += value;
            curNode2 = addElement(curNode2, -value, k, j-1, i, l);
                if((value != value) || (value*0 != value*0)){
                    printf("14 value = NaN for (%d %d %d %d) (%d %d %d %d)\n", k,j,i,l, curNode->element.k, curNode->element.j, curNode->element.i, curNode->element.l);
                    QUIT_PLUTO(1);
                }


            if(data->Vc[VX2][k][j][i] > 0){
                value = factor*dt*data->Vc[VX2][k][j][i]/(grid->x[1][j] - grid->x[1][j-1]);
                curNode = addElement(curNode, value, k,j,i,l);
                data->by[k][j][i][l] += value;
                curNode2 = addElement(curNode2, -value, k,j,i,l);
                if((value != value) || (value*0 != value*0)){
                    printf("15 value = NaN for (%d %d %d %d) (%d %d %d %d)\n", k,j,i,l, curNode->element.k, curNode->element.j, curNode->element.i, curNode->element.l);
                    QUIT_PLUTO(1);
                }

                value = -factor*dt*data->Vc[VX2][k][j-1][i]/(grid->x[1][j] - grid->x[1][j-1]);
                curNode = addElement(curNode, value, k,j-1,i,l);
                data->ay[k][j][i][l] += value;
                curNode2 = addElement(curNode2, -value, k,j-1,i,l);
                if((value != value) || (value*0 != value*0)){
                    printf("16 value = NaN for (%d %d %d %d) (%d %d %d %d)\n", k,j,i,l, curNode->element.k, curNode->element.j, curNode->element.i, curNode->element.l);
                    QUIT_PLUTO(1);
                }

            } else {
                value = -factor*dt*data->Vc[VX2][k][j][i]/(grid->x[1][j+1] - grid->x[1][j]);
                curNode = addElement(curNode, value, k,j,i,l);
                data->by[k][j][i][l] += value;
                curNode2 = addElement(curNode2, -value, k,j,i,l);
                if((value != value) || (value*0 != value*0)){
                    printf("17 value = NaN for (%d %d %d %d) (%d %d %d %d)\n", k,j,i,l, curNode->element.k, curNode->element.j, curNode->element.i, curNode->element.l);
                    QUIT_PLUTO(1);
                }

                value = factor*dt*data->Vc[VX2][k][j+1][i]/(grid->x[1][j+1] - grid->x[1][j]);
                curNode = addElement(curNode, value, k,j+1,i,l);
                data->cy[k][j][i][l] += value;
                curNode2 = addElement(curNode2, -value, k,j+1,i,l);
                if((value != value) || (value*0 != value*0)){
                    printf("18 value = NaN for (%d %d %d %d) (%d %d %d %d)\n", k,j,i,l, curNode->element.k, curNode->element.j, curNode->element.i, curNode->element.l);
                    QUIT_PLUTO(1);
                }

            }
#endif

#if INCLUDE_KDIR
            value = factor*dt*(((Dfront + D)/(grid->x[2][k+1] - grid->x[2][k]))+((D+Dback)/(grid->x[2][k] - grid->x[2][k-1])))/(grid->x[2][k+1] - grid->x[2][k-1]);
            curNode = addElement(curNode, value, k,j,i,l);
            data->bz[k][j][i][l] += value;
            curNode2 = addElement(curNode2, -value, k,j,i,l);
                if((value != value) || (value*0 != value*0)){
                    printf("19 value = NaN for (%d %d %d %d) (%d %d %d %d)\n", k,j,i,l, curNode->element.k, curNode->element.j, curNode->element.i, curNode->element.l);
                    QUIT_PLUTO(1);
                }


            value = -factor*dt*((Dfront + D)/(grid->x[2][k+1] - grid->x[2][k]))/(grid->x[2][k+1] - grid->x[2][k-1]);
            curNode = addElement(curNode, value, k+1, j, i, l);
            data->cz[k][j][i][l] += value;
            curNode2 = addElement(curNode2, -value, k+1, j, i, l);
                if((value != value) || (value*0 != value*0)){
                    printf("20 value = NaN for (%d %d %d %d) (%d %d %d %d)\n", k,j,i,l, curNode->element.k, curNode->element.j, curNode->element.i, curNode->element.l);
                    QUIT_PLUTO(1);
                }


            value = -factor*dt*((D+Dback)/(grid->x[2][k] - grid->x[2][k-1]))/(grid->x[2][k+1] - grid->x[2][k-1]);
            curNode = addElement(curNode, value, k-1, j, i, l);
            data->az[k][j][i][l] += value;
            curNode2 = addElement(curNode2, -value, k-1, j, i, l);
                if((value != value) || (value*0 != value*0)){
                    printf("21 value = NaN for (%d %d %d %d) (%d %d %d %d)\n", k,j,i,l, curNode->element.k, curNode->element.j, curNode->element.i, curNode->element.l);
                    QUIT_PLUTO(1);
                }


            if(data->Vc[VX3][k][j][i] > 0){
                value = factor*dt*data->Vc[VX3][k][j][i]/(grid->x[2][k] - grid->x[2][k-1]);
                curNode = addElement(curNode, value, k,j,i,l);
                data->bz[k][j][i][l] += value;
                curNode2 = addElement(curNode2, -value, k,j,i,l);
                if((value != value) || (value*0 != value*0)){
                    printf("22 value = NaN for (%d %d %d %d) (%d %d %d %d)\n", k,j,i,l, curNode->element.k, curNode->element.j, curNode->element.i, curNode->element.l);
                    QUIT_PLUTO(1);
                }

                value = -factor*dt*data->Vc[VX3][k-1][j][i]/(grid->x[2][k] - grid->x[2][k-1]);
                curNode = addElement(curNode, value, k-1,j,i,l);
                data->az[k][j][i][l] += value;
                curNode2 = addElement(curNode2, -value, k-1,j,i,l);
                if((value != value) || (value*0 != value*0)){
                    printf("23 value = NaN for (%d %d %d %d) (%d %d %d %d)\n", k,j,i,l, curNode->element.k, curNode->element.j, curNode->element.i, curNode->element.l);
                    QUIT_PLUTO(1);
                }

            } else {
                value = -factor*dt*data->Vc[VX3][k][j][i]/(grid->x[2][k+1] - grid->x[2][k]);
                curNode = addElement(curNode, value, k,j,i,l);
                data->bz[k][j][i][l] += value;
                curNode2 = addElement(curNode2, -value, k,j,i,l);
                if((value != value) || (value*0 != value*0)){
                    printf("24 value = NaN for (%d %d %d %d) (%d %d %d %d)\n", k,j,i,l, curNode->element.k, curNode->element.j, curNode->element.i, curNode->element.l);
                    QUIT_PLUTO(1);
                }

                value = factor*dt*data->Vc[VX3][k+1][j][i]/(grid->x[2][k+1] - grid->x[2][k]);
                curNode = addElement(curNode, value, k+1,j,i,l);
                data->cz[k][j][i][l] += value;
                curNode2 = addElement(curNode2, -value, k+1,j,i,l);
                if((value != value) || (value*0 != value*0)){
                    printf("25 value = NaN for (%d %d %d %d) (%d %d %d %d)\n", k,j,i,l, curNode->element.k, curNode->element.j, curNode->element.i, curNode->element.l);
                    QUIT_PLUTO(1);
                }

            }
#endif

#elif GEOMETRY == CYLINDRICAL
            value = factor*dt*(((grid->x[0][i+1]*Dright + grid->x[0][i]*D)/(grid->x[0][i+1] - grid->x[0][i]))+((grid->x[0][i]*D + grid->x[0][i-1]*Dleft)/(grid->x[0][i] - grid->x[0][i-1])))/(grid->x[0][i]*(grid->x[0][i+1] - grid->x[0][i-1]));
            curNode = addElement(curNode, value, k,j,i,l);
            data->bx[k][j][i][l] += value;
            curNode2 = addElement(curNode2, -value, k,j,i,l);
            if((value != value) || (value*0 != value*0)){
                printf("26 value = NaN for (%d %d %d %d) (%d %d %d %d)\n", k,j,i,l, curNode->element.k, curNode->element.j, curNode->element.i, curNode->element.l);
                QUIT_PLUTO(1);
            }

            value = -factor*dt*((grid->x[0][i+1]*Dright + grid->x[0][i]*D)/(grid->x[0][i+1] - grid->x[0][i]))/(grid->x[0][i]*(grid->x[0][i+1] - grid->x[0][i-1]));
            curNode = addElement(curNode, value, k, j, i+1, l);
            data->cx[k][j][i][l] += value;
            curNode2 = addElement(curNode2, -value, k, j, i+1, l);
            if((value != value) || (value*0 != value*0)){
                printf("27 value = NaN for (%d %d %d %d) (%d %d %d %d)\n", k,j,i,l, curNode->element.k, curNode->element.j, curNode->element.i, curNode->element.l);
                QUIT_PLUTO(1);
            }

            value = -factor*dt*((grid->x[0][i]*D+grid->x[0][i-1]*Dleft)/(grid->x[0][i] - grid->x[0][i-1]))/(grid->x[0][i]*(grid->x[0][i+1] - grid->x[0][i-1]));
            curNode = addElement(curNode, value, k, j, i-1, l);
            data->ax[k][j][i][l] += value;
            curNode2 = addElement(curNode2, -value, k, j, i-1, l);
            if((value != value) || (value*0 != value*0)){
                printf("28 value = NaN for (%d %d %d %d) (%d %d %d %d)\n", k,j,i,l, curNode->element.k, curNode->element.j, curNode->element.i, curNode->element.l);
                QUIT_PLUTO(1);
            }

            if(data->Vc[VX1][k][j][i] > 0){
                value = dt*grid->x[0][i]*data->Vc[VX1][k][j][i]/(grid->x[0][i]*(grid->x[0][i] - grid->x[0][i-1]));
                curNode = addElement(curNode, value, k,j,i,l);
                data->bx[k][j][i][l] += value;
                curNode2 = addElement(curNode2, -value, k,j,i,l);
                if((value != value) || (value*0 != value*0)){
                    printf("29 value = NaN for (%d %d %d %d) (%d %d %d %d)\n", k,j,i,l, curNode->element.k, curNode->element.j, curNode->element.i, curNode->element.l);
                    QUIT_PLUTO(1);
                }
                value = -factor*dt*grid->x[0][i-1]*data->Vc[VX1][k][j][i-1]/(grid->x[0][i]*(grid->x[0][i] - grid->x[0][i-1]));
                curNode = addElement(curNode, value, k,j,i-1,l);
                data->ax[k][j][i][l] += value;
                curNode2 = addElement(curNode2, -value, k,j,i-1,l);
                if((value != value) || (value*0 != value*0)){
                    printf("30 value = NaN for (%d %d %d %d) (%d %d %d %d)\n", k,j,i,l, curNode->element.k, curNode->element.j, curNode->element.i, curNode->element.l);
                    QUIT_PLUTO(1);
                }
            } else {
                value = -factor*dt*grid->x[0][i]*data->Vc[VX1][k][j][i]/(grid->x[0][i]*(grid->x[0][i+1] - grid->x[0][i]));
                curNode = addElement(curNode, value, k,j,i,l);
                data->bx[k][j][i][l] += value;
                curNode2 = addElement(curNode2, -value, k,j,i,l);
                if((value != value) || (value*0 != value*0)){
                    printf("31 value = NaN for (%d %d %d %d) (%d %d %d %d)\n", k,j,i,l, curNode->element.k, curNode->element.j, curNode->element.i, curNode->element.l);
                    QUIT_PLUTO(1);
                }
                value = factor*dt*grid->x[0][i+1]*data->Vc[VX1][k][j][i+1]/(grid->x[0][i]*(grid->x[0][i+1] - grid->x[0][i]));
                curNode = addElement(curNode, value, k,j,i+1,l);
                data->cx[k][j][i][l] += value;
                curNode2 = addElement(curNode2, -value, k,j,i+1,l);
                if((value != value) || (value*0 != value*0)){
                    printf("32 value = NaN for (%d %d %d %d) (%d %d %d %d)\n", k,j,i,l, curNode->element.k, curNode->element.j, curNode->element.i, curNode->element.l);
                    QUIT_PLUTO(1);
                }
            }

#if INCLUDE_JDIR
            value = factor*dt*(((Dtop + D)/(grid->x[1][j+1] - grid->x[1][j]))+((D+Dbottom)/(grid->x[1][j] - grid->x[1][j-1])))/(grid->x[1][j+1] - grid->x[1][j-1]);
            curNode = addElement(curNode, value, k,j,i,l);
            data->by[k][j][i][l] += value;
            curNode2 = addElement(curNode2, -value, k,j,i,l);
            if((value != value) || (value*0 != value*0)){
                printf("33 value = NaN for (%d %d %d %d) (%d %d %d %d)\n", k,j,i,l, curNode->element.k, curNode->element.j, curNode->element.i, curNode->element.l);
                QUIT_PLUTO(1);
            }

            value = -factor*dt*((Dtop + D)/(grid->x[1][j+1] - grid->x[1][j]))/(grid->x[1][j+1] - grid->x[1][j-1]);
            curNode = addElement(curNode, value, k, j+1, i, l);
            data->cy[k][j][i][l] += value;
            curNode2 = addElement(curNode2, -value, k, j+1, i, l);
            if((value != value) || (value*0 != value*0)){
                printf("34 value = NaN for (%d %d %d %d) (%d %d %d %d)\n", k,j,i,l, curNode->element.k, curNode->element.j, curNode->element.i, curNode->element.l);
                QUIT_PLUTO(1);
            }

            value = -factor*dt*((D+Dbottom)/(grid->x[1][j] - grid->x[1][j-1]))/(grid->x[1][j+1] - grid->x[1][j-1]);
            curNode = addElement(curNode, value, k, j-1, i, l);
            data->ay[k][j][i][l] += value;
            curNode2 = addElement(curNode2, -value, k, j-1, i, l);
            if((value != value) || (value*0 != value*0)){
                printf("35 value = NaN for (%d %d %d %d) (%d %d %d %d)\n", k,j,i,l, curNode->element.k, curNode->element.j, curNode->element.i, curNode->element.l);
                QUIT_PLUTO(1);
            }

            if(data->Vc[VX2][k][j][i] > 0){
                value = factor*dt*data->Vc[VX2][k][j][i]/(grid->x[1][j] - grid->x[1][j-1]);
                curNode = addElement(curNode, value, k,j,i,l);
                data->by[k][j][i][l] += value;
                curNode2 = addElement(curNode2, -value, k,j,i,l);
                if((value != value) || (value*0 != value*0)){
                    printf("36 value = NaN for (%d %d %d %d) (%d %d %d %d)\n", k,j,i,l, curNode->element.k, curNode->element.j, curNode->element.i, curNode->element.l);
                    QUIT_PLUTO(1);
                }
                value = -factor*dt*data->Vc[VX2][k][j-1][i]/(grid->x[1][j] - grid->x[1][j-1]);
                curNode = addElement(curNode, value, k,j-1,i,l);
                data->ay[k][j][i][l] += value;
                curNode2 = addElement(curNode2, -value, k,j-1,i,l);
                if((value != value) || (value*0 != value*0)){
                    printf("37 value = NaN for (%d %d %d %d) (%d %d %d %d)\n", k,j,i,l, curNode->element.k, curNode->element.j, curNode->element.i, curNode->element.l);
                    QUIT_PLUTO(1);
                }
            } else {
                value = -factor*dt*data->Vc[VX2][k][j][i]/(grid->x[1][j+1] - grid->x[1][j]);
                curNode = addElement(curNode, value, k,j,i,l);
                data->by[k][j][i][l] += value;
                curNode2 = addElement(curNode2, -value, k,j,i,l);
                if((value != value) || (value*0 != value*0)){
                    printf("38 value = NaN for (%d %d %d %d) (%d %d %d %d)\n", k,j,i,l, curNode->element.k, curNode->element.j, curNode->element.i, curNode->element.l);
                    QUIT_PLUTO(1);
                }
                value = factor*dt*data->Vc[VX2][k][j+1][i]/(grid->x[1][j+1] - grid->x[1][j]);
                curNode = addElement(curNode, value, k,j+1,i,l);
                data->cy[k][j][i][l] += value;
                curNode2 = addElement(curNode2, -value, k,j+1,i,l);
                if((value != value) || (value*0 != value*0)){
                    printf("39 value = NaN for (%d %d %d %d) (%d %d %d %d)\n", k,j,i,l, curNode->element.k, curNode->element.j, curNode->element.i, curNode->element.l);
                    QUIT_PLUTO(1);
                }
            }
#endif

#if INCLUDE_KDIR
            value = factor*dt*(((Dfront + D)/(grid->x[2][k+1] - grid->x[2][k]))+((D+Dback)/(grid->x[2][k] - grid->x[2][k-1])))/(grid->x[0][i]*grid->x[0][i]*(grid->x[2][k+1] - grid->x[2][k-1]));
            curNode = addElement(curNode, value, k,j,i,l);
            data->bz[k][j][i][l] += value;
            curNode2 = addElement(curNode2, -value, k,j,i,l);
                if((value != value) || (value*0 != value*0)){
                    printf("40 value = NaN for (%d %d %d %d) (%d %d %d %d)\n", k,j,i,l, curNode->element.k, curNode->element.j, curNode->element.i, curNode->element.l);
                    QUIT_PLUTO(1);
                }


            value = -factor*dt*((Dfront + D)/(grid->x[2][k+1] - grid->x[2][k]))/(grid->x[0][i]*grid->x[0][i]*(grid->x[2][k+1] - grid->x[2][k-1]));
            curNode = addElement(curNode, value, k+1, j, i, l);
            data->cz[k][j][i][l] += value;
            curNode2 = addElement(curNode2, -value, k+1, j, i, l);
                if((value != value) || (value*0 != value*0)){
                    printf("41 value = NaN for (%d %d %d %d) (%d %d %d %d)\n", k,j,i,l, curNode->element.k, curNode->element.j, curNode->element.i, curNode->element.l);
                    QUIT_PLUTO(1);
                }


            value = -factor*dt*((D+Dback)/(grid->x[2][k] - grid->x[2][k-1]))/(grid->x[0][i]*grid->x[0][i]*(grid->x[2][k+1] - grid->x[2][k-1]));
            curNode = addElement(curNode, value, k-1, j, i, l);
            data->az[k][j][i][l] += value;
            curNode2 = addElement(curNode2, -value, k-1, j, i, l);
                if((value != value) || (value*0 != value*0)){
                    printf("42 value = NaN for (%d %d %d %d) (%d %d %d %d)\n", k,j,i,l, curNode->element.k, curNode->element.j, curNode->element.i, curNode->element.l);
                    QUIT_PLUTO(1);
                }


            if(data->Vc[VX3][k][j][i] > 0){
                value = factor*dt*data->Vc[VX3][k][j][i]/(grid->x[0][i]*(grid->x[2][k] - grid->x[2][k-1]));
                curNode = addElement(curNode, value, k,j,i,l);
                data->bz[k][j][i][l] += value;
                curNode2 = addElement(curNode2, -value, k,j,i,l);
                if((value != value) || (value*0 != value*0)){
                    printf("43 value = NaN for (%d %d %d %d) (%d %d %d %d)\n", k,j,i,l, curNode->element.k, curNode->element.j, curNode->element.i, curNode->element.l);
                    QUIT_PLUTO(1);
                }

                value = -factor*dt*data->Vc[VX3][k-1][j][i]/(grid->x[0][i]*(grid->x[2][k] - grid->x[2][k-1]));
                curNode = addElement(curNode, value, k-1,j,i,l);
                data->az[k][j][i][l] += value;
                curNode2 = addElement(curNode2, -value, k-1,j,i,l);
                if((value != value) || (value*0 != value*0)){
                    printf("44 value = NaN for (%d %d %d %d) (%d %d %d %d)\n", k,j,i,l, curNode->element.k, curNode->element.j, curNode->element.i, curNode->element.l);
                    QUIT_PLUTO(1);
                }

            } else {
                value = -factor*dt*data->Vc[VX3][k][j][i]/(grid->x[0][i]*(grid->x[2][k+1] - grid->x[2][k]));
                curNode = addElement(curNode, value, k,j,i,l);
                data->bz[k][j][i][l] += value;
                curNode2 = addElement(curNode2, -value, k,j,i,l);
                if((value != value) || (value*0 != value*0)){
                    printf("45 value = NaN for (%d %d %d %d) (%d %d %d %d)\n", k,j,i,l, curNode->element.k, curNode->element.j, curNode->element.i, curNode->element.l);
                    QUIT_PLUTO(1);
                }

                value = factor*dt*data->Vc[VX3][k+1][j][i]/(grid->x[0][i]*(grid->x[2][k+1] - grid->x[2][k]));
                curNode = addElement(curNode, value, k+1,j,i,l);
                data->cz[k][j][i][l] += value;
                curNode2 = addElement(curNode2, -value, k+1,j,i,l);
                if((value != value) || (value*0 != value*0)){
                    printf("46 value = NaN for (%d %d %d %d) (%d %d %d %d)\n", k,j,i,l, curNode->element.k, curNode->element.j, curNode->element.i, curNode->element.l);
                    QUIT_PLUTO(1);
                }

            }
#endif

#elif GEOMETRY == POLAR
            value = factor*dt*(((grid->x[0][i+1]*Dright + grid->x[0][i]*D)/(grid->x[0][i+1] - grid->x[0][i]))+((grid->x[0][i]*D + grid->x[0][i-1]*Dleft)/(grid->x[0][i] - grid->x[0][i-1])))/(grid->x[0][i]*(grid->x[0][i+1] - grid->x[0][i-1]));
            curNode = addElement(curNode, value, k,j,i,l);
            data->bx[k][j][i][l] += value;
            curNode2 = addElement(curNode2, -value, k,j,i,l);
                if((value != value) || (value*0 != value*0)){
                    printf("47 value = NaN for (%d %d %d %d) (%d %d %d %d)\n", k,j,i,l, curNode->element.k, curNode->element.j, curNode->element.i, curNode->element.l);
                    QUIT_PLUTO(1);
                }


            value = -factor*dt*((grid->x[0][i+1]*Dright + grid->x[0][i]*D)/(grid->x[0][i+1] - grid->x[0][i]))/(grid->x[0][i]*(grid->x[0][i+1] - grid->x[0][i-1]));
            curNode = addElement(curNode, value, k, j, i+1, l);
            data->cx[k][j][i][l] += value;
            curNode2 = addElement(curNode2, -value, k, j, i+1, l);
                if((value != value) || (value*0 != value*0)){
                    printf("48 value = NaN for (%d %d %d %d) (%d %d %d %d)\n", k,j,i,l, curNode->element.k, curNode->element.j, curNode->element.i, curNode->element.l);
                    QUIT_PLUTO(1);
                }


            value = -factor*dt*((grid->x[0][i]*D+grid->x[0][i-1]*Dleft)/(grid->x[0][i] - grid->x[0][i-1]))/(grid->x[0][i]*(grid->x[0][i+1] - grid->x[0][i-1]));
            curNode = addElement(curNode, value, k, j, i-1, l);
            data->ax[k][j][i][l] += value;
            curNode2 = addElement(curNode2, -value, k, j, i-1, l);
                if((value != value) || (value*0 != value*0)){
                    printf("49 value = NaN for (%d %d %d %d) (%d %d %d %d)\n", k,j,i,l, curNode->element.k, curNode->element.j, curNode->element.i, curNode->element.l);
                    QUIT_PLUTO(1);
                }


            if(data->Vc[VX1][k][j][i] > 0){
                value = factor*dt*grid->x[0][i]*data->Vc[VX1][k][j][i]/(grid->x[0][i]*(grid->x[0][i] - grid->x[0][i-1]));
                curNode = addElement(curNode, value, k,j,i,l);
                data->bx[k][j][i][l] += value;
                curNode2 = addElement(curNode2, -value, k,j,i,l);
                if((value != value) || (value*0 != value*0)){
                    printf("50 value = NaN for (%d %d %d %d) (%d %d %d %d)\n", k,j,i,l, curNode->element.k, curNode->element.j, curNode->element.i, curNode->element.l);
                    QUIT_PLUTO(1);
                }

                value = -factor*dt*grid->x[0][i-1]*data->Vc[VX1][k][j][i-1]/(grid->x[0][i]*(grid->x[0][i] - grid->x[0][i-1]));
                curNode = addElement(curNode, value, k,j,i-1,l);
                data->ax[k][j][i][l] += value;
                curNode2 = addElement(curNode2, -value, k,j,i-1,l);
                if((value != value) || (value*0 != value*0)){
                    printf("51 value = NaN for (%d %d %d %d) (%d %d %d %d)\n", k,j,i,l, curNode->element.k, curNode->element.j, curNode->element.i, curNode->element.l);
                    QUIT_PLUTO(1);
                }

            } else {
                value = -factor*dt*grid->x[0][i]*data->Vc[VX1][k][j][i]/(grid->x[0][i]*(grid->x[0][i+1] - grid->x[0][i]));
                curNode = addElement(curNode, value, k,j,i,l);
                data->bx[k][j][i][l] += value;
                curNode2 = addElement(curNode2, -value, k,j,i,l);
                if((value != value) || (value*0 != value*0)){
                    printf("52 value = NaN for (%d %d %d %d) (%d %d %d %d)\n", k,j,i,l, curNode->element.k, curNode->element.j, curNode->element.i, curNode->element.l);
                    QUIT_PLUTO(1);
                }

                value = factor*dt*grid->x[0][i+1]*data->Vc[VX1][k][j][i+1]/(grid->x[0][i]*(grid->x[0][i+1] - grid->x[0][i]));
                curNode = addElement(curNode, value, k,j,i+1,l);
                data->cx[k][j][i][l] += value;
                curNode2 = addElement(curNode2, -value, k,j,i+1,l);
                if((value != value) || (value*0 != value*0)){
                    printf("53 value = NaN for (%d %d %d %d) (%d %d %d %d)\n", k,j,i,l, curNode->element.k, curNode->element.j, curNode->element.i, curNode->element.l);
                    QUIT_PLUTO(1);
                }

            }

#if INCLUDE_JDIR
            value = factor*dt*(((Dtop + D)/(grid->x[1][j+1] - grid->x[1][j]))+((D+Dbottom)/(grid->x[1][j] - grid->x[1][j-1])))/(grid->x[0][i]*grid->x[0][i]*(grid->x[1][j+1] - grid->x[1][j-1]));
            curNode = addElement(curNode, value, k,j,i,l);
            data->by[k][j][i][l] += value;
            curNode2 = addElement(curNode2, -value, k,j,i,l);
                if((value != value) || (value*0 != value*0)){
                    printf("54 value = NaN for (%d %d %d %d) (%d %d %d %d)\n", k,j,i,l, curNode->element.k, curNode->element.j, curNode->element.i, curNode->element.l);
                    QUIT_PLUTO(1);
                }


            value = -factor*dt*((Dtop + D)/(grid->x[1][j+1] - grid->x[1][j]))/(grid->x[0][i]*grid->x[0][i]*(grid->x[1][j+1] - grid->x[1][j-1]));
            curNode = addElement(curNode, value, k, j+1, i, l);
            data->cy[k][j][i][l] += value;
            curNode2 = addElement(curNode2, -value, k, j+1, i, l);
                if((value != value) || (value*0 != value*0)){
                    printf("55 value = NaN for (%d %d %d %d) (%d %d %d %d)\n", k,j,i,l, curNode->element.k, curNode->element.j, curNode->element.i, curNode->element.l);
                    QUIT_PLUTO(1);
                }


            value = -factor*dt*((D+Dbottom)/(grid->x[1][j] - grid->x[1][j-1]))/(grid->x[0][i]*grid->x[0][i]*(grid->x[1][j+1] - grid->x[1][j-1]));
            curNode = addElement(curNode, value, k, j-1, i, l);
            data->ay[k][j][i][l] += value;
            curNode2 = addElement(curNode2, -value, k, j-1, i, l);
                if((value != value) || (value*0 != value*0)){
                    printf("56 value = NaN for (%d %d %d %d) (%d %d %d %d)\n", k,j,i,l, curNode->element.k, curNode->element.j, curNode->element.i, curNode->element.l);
                    QUIT_PLUTO(1);
                }


            if(data->Vc[VX2][k][j][i] > 0){
                value = factor*dt*data->Vc[VX2][k][j][i]/(grid->x[0][i]*(grid->x[1][j] - grid->x[1][j-1]));
                curNode = addElement(curNode, value, k,j,i,l);
                data->by[k][j][i][l] += value;
                curNode2 = addElement(curNode2, -value, k,j,i,l);
                if((value != value) || (value*0 != value*0)){
                    printf("57 value = NaN for (%d %d %d %d) (%d %d %d %d)\n", k,j,i,l, curNode->element.k, curNode->element.j, curNode->element.i, curNode->element.l);
                    QUIT_PLUTO(1);
                }

                value = -factor*dt*data->Vc[VX2][k][j-1][i]/(grid->x[0][i]*(grid->x[1][j] - grid->x[1][j-1]));
                curNode = addElement(curNode, value, k,j-1,i,l);
                data->ay[k][j][i][l] += value;
                curNode2 = addElement(curNode2, -value, k,j-1,i,l);
                if((value != value) || (value*0 != value*0)){
                    printf("58 value = NaN for (%d %d %d %d) (%d %d %d %d)\n", k,j,i,l, curNode->element.k, curNode->element.j, curNode->element.i, curNode->element.l);
                    QUIT_PLUTO(1);
                }

            } else {
                value = -factor*dt*data->Vc[VX2][k][j][i]/(grid->x[0][i]*(grid->x[1][j+1] - grid->x[1][j]));
                curNode = addElement(curNode, value, k,j,i,l);
                data->by[k][j][i][l] += value;
                curNode2 = addElement(curNode2, -value, k,j,i,l);
                if((value != value) || (value*0 != value*0)){
                    printf("59 value = NaN for (%d %d %d %d) (%d %d %d %d)\n", k,j,i,l, curNode->element.k, curNode->element.j, curNode->element.i, curNode->element.l);
                    QUIT_PLUTO(1);
                }

                value = factor*dt*data->Vc[VX2][k][j+1][i]/(grid->x[0][i]*(grid->x[1][j+1] - grid->x[1][j]));
                curNode = addElement(curNode, value, k,j+1,i,l);
                data->cy[k][j][i][l] += value;
                curNode2 = addElement(curNode2, -value, k,j+1,i,l);
                if((value != value) || (value*0 != value*0)){
                    printf("60 value = NaN for (%d %d %d %d) (%d %d %d %d)\n", k,j,i,l, curNode->element.k, curNode->element.j, curNode->element.i, curNode->element.l);
                    QUIT_PLUTO(1);
                }

            }
#endif

#if INCLUDE_KDIR
            value = factor*dt*(((Dfront + D)/(grid->x[2][k+1] - grid->x[2][k]))+((D+Dback)/(grid->x[2][k] - grid->x[2][k-1])))/(grid->x[2][k+1] - grid->x[2][k-1]);
            curNode = addElement(curNode, value, k,j,i,l);
            data->bz[k][j][i][l] += value;
            curNode2 = addElement(curNode2, -value, k,j,i,l);
                if((value != value) || (value*0 != value*0)){
                    printf("61 value = NaN for (%d %d %d %d) (%d %d %d %d)\n", k,j,i,l, curNode->element.k, curNode->element.j, curNode->element.i, curNode->element.l);
                    QUIT_PLUTO(1);
                }


            value = -factor*dt*((Dfront + D)/(grid->x[2][k+1] - grid->x[2][k]))/(grid->x[2][k+1] - grid->x[2][k-1]);
            curNode = addElement(curNode, value, k+1, j, i, l);
            data->cz[k][j][i][l] += value;
            curNode2 = addElement(curNode2, -value, k+1, j, i, l);
                if((value != value) || (value*0 != value*0)){
                    printf("62 value = NaN for (%d %d %d %d) (%d %d %d %d)\n", k,j,i,l, curNode->element.k, curNode->element.j, curNode->element.i, curNode->element.l);
                    QUIT_PLUTO(1);
                }


            value = -factor*dt*((D+Dback)/(grid->x[2][k] - grid->x[2][k-1]))/(grid->x[2][k+1] - grid->x[2][k-1]);
            curNode = addElement(curNode, value, k-1, j, i, l);
            data->az[k][j][i][l] += value;
            curNode2 = addElement(curNode2, -value, k-1, j, i, l);
                if((value != value) || (value*0 != value*0)){
                    printf("63 value = NaN for (%d %d %d %d) (%d %d %d %d)\n", k,j,i,l, curNode->element.k, curNode->element.j, curNode->element.i, curNode->element.l);
                    QUIT_PLUTO(1);
                }


            if(data->Vc[VX3][k][j][i] > 0){
                value = factor*dt*data->Vc[VX3][k][j][i]/(grid->x[2][k] - grid->x[2][k-1]);
                curNode = addElement(curNode, value, k,j,i,l);
                data->bz[k][j][i][l] += value;
                curNode2 = addElement(curNode2, -value, k,j,i,l);
                if((value != value) || (value*0 != value*0)){
                    printf("64 value = NaN for (%d %d %d %d) (%d %d %d %d)\n", k,j,i,l, curNode->element.k, curNode->element.j, curNode->element.i, curNode->element.l);
                    QUIT_PLUTO(1);
                }

                value = -factor*dt*data->Vc[VX3][k-1][j][i]/(grid->x[2][k] - grid->x[2][k-1]);
                curNode = addElement(curNode, value, k-1,j,i,l);
                data->az[k][j][i][l] += value;
                curNode2 = addElement(curNode2, -value, k-1,j,i,l);
                if((value != value) || (value*0 != value*0)){
                    printf("65 value = NaN for (%d %d %d %d) (%d %d %d %d)\n", k,j,i,l, curNode->element.k, curNode->element.j, curNode->element.i, curNode->element.l);
                    QUIT_PLUTO(1);
                }

            } else {
                value = -factor*dt*data->Vc[VX3][k][j][i]/(grid->x[2][k+1] - grid->x[2][k]);
                curNode = addElement(curNode, value, k,j,i,l);
                data->bz[k][j][i][l] += value;
                curNode2 = addElement(curNode2, -value, k,j,i,l);
                if((value != value) || (value*0 != value*0)){
                    printf("66 value = NaN for (%d %d %d %d) (%d %d %d %d)\n", k,j,i,l, curNode->element.k, curNode->element.j, curNode->element.i, curNode->element.l);
                    QUIT_PLUTO(1);
                }

                value = factor*dt*data->Vc[VX3][k+1][j][i]/(grid->x[2][k+1] - grid->x[2][k]);
                curNode = addElement(curNode, value, k+1,j,i,l);
                data->cz[k][j][i][l] += value;
                curNode2 = addElement(curNode2, -value, k+1,j,i,l);
                if((value != value) || (value*0 != value*0)){
                    printf("67 value = NaN for (%d %d %d %d) (%d %d %d %d)\n", k,j,i,l, curNode->element.k, curNode->element.j, curNode->element.i, curNode->element.l);
                    QUIT_PLUTO(1);
                }
            }
#endif

#elif GEOMETRY == SPHERICAL
            value = factor*dt*(((grid->x[0][i+1]*grid->x[0][i+1]*Dright + grid->x[0][i]*grid->x[0][i]*D)/(grid->x[0][i+1] - grid->x[0][i]))+((grid->x[0][i]*grid->x[0][i]*D+grid->x[0][i-1]*grid->x[0][i-1]*Dleft)/(grid->x[0][i] - grid->x[0][i-1])))/(grid->x[0][i]*grid->x[0][i]*(grid->x[0][i+1] - grid->x[0][i-1]));
            curNode = addElement(curNode, value, k,j,i,l);
            data->bx[k][j][i][l] += value;
            curNode2 = addElement(curNode2, -value, k,j,i,l);
                if((value != value) || (value*0 != value*0)){
                    printf("68 value = NaN for (%d %d %d %d) (%d %d %d %d)\n", k,j,i,l, curNode->element.k, curNode->element.j, curNode->element.i, curNode->element.l);
                    QUIT_PLUTO(1);
                }


            value = -factor*dt*((grid->x[0][i+1]*grid->x[0][i+1]*Dright + grid->x[0][i]*grid->x[0][i]*D)/(grid->x[0][i+1] - grid->x[0][i]))/(grid->x[0][i]*grid->x[0][i]*(grid->x[0][i+1] - grid->x[0][i-1]));
            curNode = addElement(curNode, value, k, j, i+1, l);
            data->cx[k][j][i][l] += value;
            curNode2 = addElement(curNode2, -value, k, j, i+1, l);
                if((value != value) || (value*0 != value*0)){
                    printf("69 value = NaN for (%d %d %d %d) (%d %d %d %d)\n", k,j,i,l, curNode->element.k, curNode->element.j, curNode->element.i, curNode->element.l);
                    QUIT_PLUTO(1);
                }


            value = -factor*dt*((grid->x[0][i]*grid->x[0][i]*D+grid->x[0][i-1]*grid->x[0][i-1]*Dleft)/(grid->x[0][i] - grid->x[0][i-1]))/(grid->x[0][i]*grid->x[0][i]*(grid->x[0][i+1] - grid->x[0][i-1]));
            curNode = addElement(curNode, value, k, j, i-1, l);
            data->ax[k][j][i][l] += value;
            curNode2 = addElement(curNode2, -value, k, j, i-1, l);
                if((value != value) || (value*0 != value*0)){
                    printf("70 value = NaN for (%d %d %d %d) (%d %d %d %d)\n", k,j,i,l, curNode->element.k, curNode->element.j, curNode->element.i, curNode->element.l);
                    QUIT_PLUTO(1);
                }


            if(data->Vc[VX1][k][j][i] > 0){
                value = factor*dt*grid->x[0][i]*grid->x[0][i]*data->Vc[VX1][k][j][i]/(grid->x[0][i]*grid->x[0][i]*(grid->x[0][i] - grid->x[0][i-1]));
                curNode = addElement(curNode, value, k,j,i,l);
                data->bx[k][j][i][l] += value;
                curNode2 = addElement(curNode2, -value, k,j,i,l);
                if((value != value) || (value*0 != value*0)){
                    printf("71 value = NaN for (%d %d %d %d) (%d %d %d %d)\n", k,j,i,l, curNode->element.k, curNode->element.j, curNode->element.i, curNode->element.l);
                    QUIT_PLUTO(1);
                }

                value = -factor*dt*grid->x[0][i-1]*grid->x[0][i-1]*data->Vc[VX1][k][j][i-1]/(grid->x[0][i]*grid->x[0][i]*(grid->x[0][i] - grid->x[0][i-1]));
                curNode = addElement(curNode, value, k,j,i-1,l);
                data->ax[k][j][i][l] += value;
                curNode2 = addElement(curNode2, -value, k,j,i-1,l);
                if((value != value) || (value*0 != value*0)){
                    printf("72 value = NaN for (%d %d %d %d) (%d %d %d %d)\n", k,j,i,l, curNode->element.k, curNode->element.j, curNode->element.i, curNode->element.l);
                    QUIT_PLUTO(1);
                }

            } else {
                value = -factor*dt*grid->x[0][i]*grid->x[0][i]*data->Vc[VX1][k][j][i]/(grid->x[0][i]*grid->x[0][i]*(grid->x[0][i+1] - grid->x[0][i]));
                curNode = addElement(curNode, value, k,j,l);
                data->bx[k][j][i][l] += value;
                curNode2 = addElement(curNode2, -value, k,j,i,l);
                if((value != value) || (value*0 != value*0)){
                    printf("73 value = NaN for (%d %d %d %d) (%d %d %d %d)\n", k,j,i,l, curNode->element.k, curNode->element.j, curNode->element.i, curNode->element.l);
                    QUIT_PLUTO(1);
                }

                value = factor*dt*grid->x[0][i+1]*grid->x[0][i+1]*data->Vc[VX1][k][j][i+1]/(grid->x[0][i]*grid->x[0][i]*(grid->x[0][i+1] - grid->x[0][i]));
                curNode = addElement(curNode, value, k,j,i+1,l);
                data->cx[k][j][i][l] += value;
                curNode2 = addElement(curNode2, -value, k,j,i+1,l);
                if((value != value) || (value*0 != value*0)){
                    printf("74 value = NaN for (%d %d %d %d) (%d %d %d %d)\n", k,j,i,l, curNode->element.k, curNode->element.j, curNode->element.i, curNode->element.l);
                    QUIT_PLUTO(1);
                }

            }
#if INCLUDE_JDIR
            value = factor*dt*(((Dtop*sin(grid->x[1][j+1]) + D*sin(grid->x[1][j]))/(grid->x[1][j+1] - grid->x[1][j]))+((D*sin(grid->x[1][j])+Dbottom*sin(grid->x[1][j-1]))/(grid->x[1][j] - grid->x[1][j-1])))/(grid->x[0][i]*grid->x[0][i]*sin(grid->x[1][j])*(grid->x[1][j+1] - grid->x[1][j-1]));
            curNode = addElement(curNode, value, k,j,i,l);
            data->by[k][j][i][l] += value;
            curNode2 = addElement(curNode2, -value, k,j,i,l);
                if((value != value) || (value*0 != value*0)){
                    printf("75 value = NaN for (%d %d %d %d) (%d %d %d %d)\n", k,j,i,l, curNode->element.k, curNode->element.j, curNode->element.i, curNode->element.l);
                    QUIT_PLUTO(1);
                }


            value = -factor*dt*((Dtop*sin(grid->x[1][j+1]) + D*sin(x->grid->x[1][j]))/(grid->x[1][j+1] - grid->x[1][j]))/(grid->x[0][i]*grid->x[0][i]*sin(grid->x[1][j])*(grid->x[1][j+1] - grid->x[1][j-1]));
            curNode = addElement(curNode, value, k, j+1, i, l);
            data->cy[k][j][i][l] += value;
            curNode2 = addElement(curNode2, -value, k, j+1, i, l);
                if((value != value) || (value*0 != value*0)){
                    printf("76 value = NaN for (%d %d %d %d) (%d %d %d %d)\n", k,j,i,l, curNode->element.k, curNode->element.j, curNode->element.i, curNode->element.l);
                    QUIT_PLUTO(1);
                }


            value = -factor*dt*((D*sin(grid->x[1][j])+Dbottom*sin(grid->x[1][j-1]))/(grid->x[1][j] - grid->x[1][j-1]))/(grid->x[0][i]*grid->x[0][i]*sin(grid->x[1][j])*(grid->x[1][j+1] - grid->x[1][j-1]));
            curNode = addElement(curNode, value, k, j-1, i, l);
            data->ay[k][j][i][l] += value;
            curNode2 = addElement(curNode2, -value, k, j-1, i, l);
                if((value != value) || (value*0 != value*0)){
                    printf("77 value = NaN for (%d %d %d %d) (%d %d %d %d)\n", k,j,i,l, curNode->element.k, curNode->element.j, curNode->element.i, curNode->element.l);
                    QUIT_PLUTO(1);
                }


            if(data->Vc[VX2][k][j][i] > 0){
                value = factor*dt*data->Vc[VX2][k][j][i]/(grid->x[0][i]*grid->x[0][i]*(grid->x[1][j] - grid->x[1][j-1]));
                curNode = addElement(curNode, value, k,j,i,l);
                data->by[k][j][i][l] += value;
                curNode2 = addElement(curNode2, -value, k,j,i,l);
                if((value != value) || (value*0 != value*0)){
                    printf("78 value = NaN for (%d %d %d %d) (%d %d %d %d)\n", k,j,i,l, curNode->element.k, curNode->element.j, curNode->element.i, curNode->element.l);
                    QUIT_PLUTO(1);
                }

                value = -factor*dt*sin(grid->x[1][j-1])*data->Vc[VX2][k][j-1][i]/(grid->x[0][i]*grid->x[0][i]*sin(grid->x[1][j])*(grid->x[1][j] - grid->x[1][j-1]));
                curNode = addElement(curNode, value, k,j-1,i,l);
                data->ay[k][j][i][l] += value;
                curNode2 = addElement(curNode2, -value, k,j-1,i,l);
                if((value != value) || (value*0 != value*0)){
                    printf("79 value = NaN for (%d %d %d %d) (%d %d %d %d)\n", k,j,i,l, curNode->element.k, curNode->element.j, curNode->element.i, curNode->element.l);
                    QUIT_PLUTO(1);
                }

            } else {
                value = -factor*dt*data->Vc[VX2][k][j][i]/(grid->x[0][i]*grid->x[0][i]*(grid->x[1][j+1] - grid->x[1][j]));
                curNode = addElement(curNode, value, k,j,i,l);
                data->by[k][j][i][l] += value;
                curNode2 = addElement(curNode2, -value, k,j,i,l);
                if((value != value) || (value*0 != value*0)){
                    printf("80 value = NaN for (%d %d %d %d) (%d %d %d %d)\n", k,j,i,l, curNode->element.k, curNode->element.j, curNode->element.i, curNode->element.l);
                    QUIT_PLUTO(1);
                }

                value = factor*dt*sin(grid->x[1][j+1])*data->Vc[VX2][k][j+1][i]/(grid->x[0][i]*grid->x[0][i]*sin(grid->x[1][j])*(grid->x[1][j+1] - grid->x[1][j]));
                curNode = addElement(curNode, value, k,j+1,i,l);
                data->cy[k][j][i][l] += value;
                curNode2 = addElement(curNode2, -value, k,j+1,i,l);
                if((value != value) || (value*0 != value*0)){
                    printf("81 value = NaN for (%d %d %d %d) (%d %d %d %d)\n", k,j,i,l, curNode->element.k, curNode->element.j, curNode->element.i, curNode->element.l);
                    QUIT_PLUTO(1);
                }

            }
#endif

#if INCLUDE_KDIR
            double r2sin2theta = grid->x[0][i]*grid->x[0][i]*sin(grid->x[1][j])*sin(grid->x[1][j]);
            value = factor*dt*(((Dfront + D)/(grid->x[2][k+1] - grid->x[2][k]))+((D+Dback)/(grid->x[2][k] - grid->x[2][k-1])))/(r2sin2theta*(grid->x[2][k+1] - grid->x[2][k-1]));
            curNode = addElement(curNode, value, k,j,i,l);
            data->bz[k][j][i][l] += value;
            curNode2 = addElement(curNode2, -value, k,j,i,l);
                if((value != value) || (value*0 != value*0)){
                    printf("82 value = NaN for (%d %d %d %d) (%d %d %d %d)\n", k,j,i,l, curNode->element.k, curNode->element.j, curNode->element.i, curNode->element.l);
                    QUIT_PLUTO(1);
                }


            value = -factor*dt*((Dfront + D)/(grid->x[2][k+1] - grid->x[2][k]))/(r2sin2theta*(grid->x[2][k+1] - grid->x[2][k-1]));
            curNode = addElement(curNode, value, k+1, j, i, l);
            data->cz[k][j][i][l] += value;
            curNode2 = addElement(curNode2, -value, k+1, j, i, l);
                if((value != value) || (value*0 != value*0)){
                    printf("83 value = NaN for (%d %d %d %d) (%d %d %d %d)\n", k,j,i,l, curNode->element.k, curNode->element.j, curNode->element.i, curNode->element.l);
                    QUIT_PLUTO(1);
                }


            value = -factor*dt*((D+Dback)/(grid->x[2][k] - grid->x[2][k-1]))/(r2sin2theta*(grid->x[2][k+1] - grid->x[2][k-1]));
            curNode = addElement(curNode, value, k-1, j, i, l);
            data->az[k][j][i][l] += value;
            curNode2 = addElement(curNode2, -value, k-1, j, i, l);
                if((value != value) || (value*0 != value*0)){
                    printf("84 value = NaN for (%d %d %d %d) (%d %d %d %d)\n", k,j,i,l, curNode->element.k, curNode->element.j, curNode->element.i, curNode->element.l);
                    QUIT_PLUTO(1);
                }


            if(data->Vc[VX3][k][j][i] > 0){
                value = factor*dt*data->Vc[VX3][k][j][i]/(grid->x[0][i]*sin(grid->x[1][j])*(grid->x[2][k] - grid->x[2][k-1]));
                curNode = addElement(curNode, value, k,j,i,l);
                data->bz[k][j][i][l] += value;
                curNode2 = addElement(curNode2, -value, k,j,i,l);
                if((value != value) || (value*0 != value*0)){
                    printf("85 value = NaN for (%d %d %d %d) (%d %d %d %d)\n", k,j,i,l, curNode->element.k, curNode->element.j, curNode->element.i, curNode->element.l);
                    QUIT_PLUTO(1);
                }

                value = -factor*dt*data->Vc[VX3][k-1][j][i]/(grid->x[0][i]*sin(grid->x[1][j])*(grid->x[2][k] - grid->x[2][k-1]));
                curNode = addElement(curNode, value, k-1,j,i,l);
                data->az[k][j][i][l] += value;
                curNode2 = addElement(curNode2, -value, k-1,j,i,l);
                if((value != value) || (value*0 != value*0)){
                    printf("86 value = NaN for (%d %d %d %d) (%d %d %d %d)\n", k,j,i,l, curNode->element.k, curNode->element.j, curNode->element.i, curNode->element.l);
                    QUIT_PLUTO(1);
                }

            } else {
                value = -factor*dt*data->Vc[VX3][k][j][i]/(grid->x[0][i]*sin(grid->x[1][j])*(grid->x[2][k+1] - grid->x[2][k]));
                curNode = addElement(curNode, value, k,j,i,l);
                data->bz[k][j][i][l] += value;
                curNode2 = addElement(curNode2, -value, k,j,i,l);
                if((value != value) || (value*0 != value*0)){
                    printf("87 value = NaN for (%d %d %d %d) (%d %d %d %d)\n", k,j,i,l, curNode->element.k, curNode->element.j, curNode->element.i, curNode->element.l);
                    QUIT_PLUTO(1);
                }

                value = factor*dt*data->Vc[VX3][k+1][j][i]/(grid->x[0][i]*sin(grid->x[1][j])*(grid->x[2][k+1] - grid->x[2][k]));
                curNode = addElement(curNode, value, k+1,j,i,l);
                data->cz[k][j][i][l] += value;
                curNode2 = addElement(curNode2, -value, k+1,j,i,l);
                if((value != value) || (value*0 != value*0)){
                    printf("88 value = NaN for (%d %d %d %d) (%d %d %d %d)\n", k,j,i,l, curNode->element.k, curNode->element.j, curNode->element.i, curNode->element.l);
                    QUIT_PLUTO(1);
                }

            }
#endif
#else
            printf("wrong geometry\n");
            QUIT_PLUTO(1);
#endif
            double diagonal = 0;
            double notdiagonal = 0;
            curNode = data->matrix[k][j][i][l];

            while(curNode != NULL){
                if((curNode->element.i == i) && (curNode->element.j == j) && (curNode->element.k == k) && (curNode->element.l == l)){
                    diagonal = diagonal + curNode->element.value;
                } else {
                    notdiagonal = notdiagonal + fabs(curNode->element.value);
                }
                curNode = curNode->next;
            }

            if(diagonal < notdiagonal){
                //printf("diagonal = %g < notdiagonal = %g, %d %d %d %d\n", diagonal, notdiagonal, i, j, k, l);
                //exit(0);
            }
        }
    }
    //crank-nickelson
    setBoundaryRightPartToZero(data, grid);
    int  par_dim[3] = {0, 0, 0};
    DIM_EXPAND(par_dim[0] = grid->nproc[IDIR] > 1;  ,
               par_dim[1] = grid->nproc[JDIR] > 1;  ,
               par_dim[2] = grid->nproc[KDIR] > 1;)
    exchangeLargeVector(data->Fkin, NMOMENTUM, par_dim, SZ_stagx, periodicX, periodicY, periodicZ);
    multiplySpecialMatrixVector(data->rightPart, data->rightPartMatrix, data->Fkin, NMOMENTUM, par_dim);
    /*TOT_LOOP(k,j,i){
       for(int l = 0; l < NMOMENTUM; ++l){
           data->rightPart[k][j][i][l] = data->Fkin[k][j][i][l];
       }
    }*/

    setBoundaryRightPartToZero(data, grid);

#if PARTICLES_KIN_SOLVER == THREE_DIAGONAL
#ifdef PARALLEL
    register int nd;
    int myrank, nproc;
    int ndim, gp, nleft, nright, tag1, tag2;
    int sendb, recvb;
    MPI_Datatype itype;
    MPI_Comm comm;
    MPI_Status status;
    struct szz *s;

    s = sz_stack[SZ_stagx];


#if INCLUDE_IDIR
    comm = s->oned_comm[0];
    myrank = s->lrank[0];
    nproc = s->lsize[0];
    ndim = s->ndim;
    //printf("parallel solver x\n");
    if(par_dim[0] != 0){
    parallelThreeDiagonalSolverX(data->Fkin, data->rightPart, data->ax, data->bx, data->cx, NMOMENTUM, nproc, myrank, comm);
    } else {
    noparallelThreeDiagonalSolverX(data->Fkin, data->rightPart, data->ax, data->bx, data->cx, NMOMENTUM);
    }
    exchangeLargeVector(data->Fkin, NMOMENTUM, par_dim, SZ_stagx);
    TOT_LOOP(k,j,i){
        for(int l = 0; l < NMOMENTUM; ++l){
            data->rightPart[k][j][i][l] = data->Fkin[k][j][i][l];
        }
    }
    setBoundaryRightPartToZero(data, grid);
#endif
#if INCLUDE_JDIR
    comm = s->oned_comm[1];
    myrank = s->lrank[1];
    nproc = s->lsize[1];
    ndim = s->ndim;
    //printf("parallel solver y\n");
    /*if(par_dim[1] != 0){
    parallelThreeDiagonalSolverY(data->Fkin, data->rightPart, data->ay, data->by, data->cy, NMOMENTUM, nproc, myrank, comm);
    } else {
    noparallelThreeDiagonalSolverY(data->Fkin, data->rightPart, data->ay, data->by, data->cy, NMOMENTUM);
    }*/
    exchangeLargeVector(data->Fkin, NMOMENTUM, par_dim, SZ_stagx);
    TOT_LOOP(k,j,i){
        for(int l = 0; l < NMOMENTUM; ++l){
            data->rightPart[k][j][i][l] = data->Fkin[k][j][i][l];
        }
    }
    setBoundaryRightPartToZero(data, grid);
#endif
#if INCLUDE_KDIR
    comm = s->oned_comm[2];
    myrank = s->lrank[2];
    nproc = s->lsize[2];
    ndim = s->ndim;
    //printf("parallel solver z\n");
    if(par_dim[2] != 0){
    parallelThreeDiagonalSolverZ(data->Fkin, data->rightPart, data->az, data->bz, data->cz, NMOMENTUM, nproc, myrank, comm);
    } else {
    noparallelThreeDiagonalSolverZ(data->Fkin, data->rightPart, data->az, data->bz, data->cz, NMOMENTUM);
    }
    exchangeLargeVector(data->Fkin, NMOMENTUM, par_dim, SZ_stagx);
    TOT_LOOP(k,j,i){
        for(int l = 0; l < NMOMENTUM; ++l){
            data->rightPart[k][j][i][l] = data->Fkin[k][j][i][l];
        }
    }
    setBoundaryRightPartToZero(data, grid);
#endif
    //printf("noparallel solver p\n");
    noparallelThreeDiagonalSolverP(data->Fkin, data->rightPart, data->ap, data->bp, data->cp, NMOMENTUM);
    exchangeLargeVector(data->Fkin, NMOMENTUM, par_dim, SZ_stagx);
    TOT_LOOP(k,j,i){
        for(int l = 0; l < NMOMENTUM; ++l){
            data->rightPart[k][j][i][l] = data->Fkin[k][j][i][l];
        }
    }


#else
#if INCLUDE_IDIR
    //printf("noparallel solver x\n");
    noparallelThreeDiagonalSolverX(data->Fkin, data->rightPart, data->ax, data->bx, data->cx, NMOMENTUM);
    exchangeLargeVector(data->Fkin, NMOMENTUM, par_dim, SZ_stagx, periodicX, periodicY, periodicZ);
    DOM_LOOP(k,j,i){
        for(int l = 0; l < NMOMENTUM; ++l){
            data->rightPart[k][j][i][l] = data->Fkin[k][j][i][l];
        }
    }
    setBoundaryRightPartToZero(data, grid);
#endif
#if INCLUDE_JDIR
    //printf("noparallel solver y\n");
    //noparallelThreeDiagonalSolverY(data->Fkin, data->rightPart, data->ay, data->by, data->cy, NMOMENTUM);
    exchangeLargeVector(data->Fkin, NMOMENTUM, par_dim, SZ_stagx, periodicX, periodicY, periodicZ);
    DOM_LOOP(k,j,i){
        for(int l = 0; l < NMOMENTUM; ++l){
            data->rightPart[k][j][i][l] = data->Fkin[k][j][i][l];
        }
    }
    setBoundaryRightPartToZero(data, grid);
#endif
#if INCLUDE_KDIR
    //printf("noparallel solver z\n");
    noparallelThreeDiagonalSolverZ(data->Fkin, data->rightPart, data->az, data->bz, data->cz, NMOMENTUM);
    exchangeLargeVector(data->Fkin, NMOMENTUM, par_dim, SZ_stagx);
    DOM_LOOP(k,j,i){
        for(int l = 0; l < NMOMENTUM; ++l){
            data->rightPart[k][j][i][l] = data->Fkin[k][j][i][l];
        }
    }
    setBoundaryRightPartToZero(data, grid);
#endif
    //printf("noparallel solver p\n");
    noparallelThreeDiagonalSolverP(data->Fkin, data->rightPart, data->ap, data->bp, data->cp, NMOMENTUM);
    exchangeLargeVector(data->Fkin, NMOMENTUM, par_dim, SZ_stagx, periodicX, periodicY, periodicZ);
    DOM_LOOP(k,j,i){
        for(int l = 0; l < NMOMENTUM; ++l){
            data->rightPart[k][j][i][l] = data->Fkin[k][j][i][l];
        }
    }

#endif
#else

    double**** initialVector = (double****) malloc(NX3_TOT*sizeof(double***));
    for(int k = 0; k < NX3_TOT; ++k){
        initialVector[k] = (double***) malloc(NX2_TOT*sizeof(double**));
        for(int j = 0; j < NX2_TOT; ++j){
            initialVector[k][j] = (double**) malloc(NX1_TOT*sizeof(double*));
            for(int i = 0; i < NX1_TOT; ++i){
                initialVector[k][j][i] = (double*) malloc(NMOMENTUM*sizeof(double));
                for(int l = 0; l < NMOMENTUM; ++l){
                    initialVector[k][j][i][l] = data->Fkin[k][j][i][l];
                }
            }
        }
    }

    double precision = 0.1/data->p_grid[NMOMENTUM - 1];
    precision = 1E-5;
    generalizedMinimalResidualMethod(grid, data->matrix, data->rightPart, data->Fkin, initialVector, data->gmresBasis, NMOMENTUM, precision, MAX_GMRES_ITERATIONS, 1);
    //generalizedMinimalResidualMethod1(grid, data->matrix, data->rightPart, data->Fkin, data->gmresBasis, NMOMENTUM, precision, MAX_GMRES_ITERATIONS, 1);
    //biconjugateStabilizedGradientMethod(grid, data->matrix, data->rightPart, data->Fkin, NMOMENTUM, precision, MAX_GMRES_ITERATIONS, 1);

    for(int k = 0; k < NX3_TOT; ++k){
        for(int j = 0; j < NX2_TOT; ++j){
            for(int i = 0; i < NX1_TOT; ++i){
                free(initialVector[k][j][i]);
            }
            free(initialVector[k][j]);
        }
        free(initialVector[k]);
    }
    free(initialVector);
#endif

    TOT_LOOP(k,j,i){
        for(int l = 0; l < NMOMENTUM; ++l){
            if(data->Fkin[k][j][i][l] < 0){
                data->Fkin[k][j][i][l] = 0;
            }
        }
    }

    inv_dt = MAX(inv_dt, Dts->invDt_advection);
    inv_dt = MAX(inv_dt, Dts->invDt_acceleration);

    Dts->invDt_particles = inv_dt;
    Dts->omega_particles = inv_dt;

    DOM_LOOP(k,j,i){
        data->Pkin[k][j][i] = 0.0;
        for(int l = 0; l < NMOMENTUM; ++l){
            double dp = 0;
            if(l == 0){
                dp = data->p_grid[1] - data->p_grid[0];
            } else {
                dp = data->p_grid[l] - data->p_grid[l-1];
            }

            double beta = data->p_grid[l]/sqrt(1.0 + data->p_grid[l]*data->p_grid[l]);

            data->Pkin[k][j][i] = data->Pkin[k][j][i] + (1.0/3.0)*dp*data->p_grid[l]*PARTICLES_KIN_MASS*PARTICLES_KIN_C*beta*PARTICLES_KIN_C*data->Fkin[k][j][i][l]/(data->p_grid[l]);

            data->Jkin1[k][j][i][l] = 0.0;
            data->Jkin2[k][j][i][l] = 0.0;
            data->Jkin3[k][j][i][l] = 0.0;
            double D = evaluateDiffusionCoefficient(data, i, j, k, data->p_grid[l]);
#if GEOMETRY == CARTESIAN
#if INCLUDE_IDIR
            data->Jkin1[k][j][i][l] = - (D*(data->Fkin[k][j][i+1][l] - data->Fkin[k][j][i-1][l])/(grid->x[0][i+1] - grid->x[0][i-1]))/data->p_grid[l];
#endif
#if INCLUDE_JDIR
            data->Jkin2[k][j][i][l] = - (D*(data->Fkin[k][j+1][i][l] - data->Fkin[k][j-1][i][l])/(grid->x[1][j+1] - grid->x[1][j-1]))/data->p_grid[l];
#endif
#if INCLUDE_KDIR
            data->Jkin3[k][j][i][l] = - (D*(data->Fkin[k+1][j][i][l] - data->Fkin[k-1][j][i][l])/(grid->x[2][k+1] - grid->x[2][k-1]))/data->p_grid[l];
#endif
#elif GEOMETRY == CYLINDRICAL
#if INCLUDE_IDIR
            data->Jkin1[k][j][i][l] = - (D*(data->Fkin[k][j][i+1][l] - data->Fkin[k][j][i-1][l])/(grid->x[0][i+1] - grid->x[0][i-1]))/data->p_grid[l];
#endif
#if INCLUDE_JDIR
            data->Jkin2[k][j][i][l] = - (D*(data->Fkin[k][j+1][i][l] - data->Fkin[k][j-1][i][l])/(grid->x[1][j+1] - grid->x[1][j-1]))/data->p_grid[l];
#endif
#if INCLUDE_KDIR
            data->Jkin3[k][j][i][l] = - (D*(data->Fkin[k+1][j][i][l] - data->Fkin[k-1][j][i][l])/(grid->x[0][i]*(grid->x[2][k+1] - grid->x[2][k-1])))/data->p_grid[l];
#endif
#elif GEOMETRY == POLAR
#if INCLUDE_IDIR
            data->Jkin1[k][j][i][l] = - (D*(data->Fkin[k][j][i+1][l] - data->Fkin[k][j][i-1][l])/(grid->x[0][i+1] - grid->x[0][i-1]))/data->p_grid[l];
#endif
#if INCLUDE_JDIR
            data->Jkin2[k][j][i][l] = - (D*(data->Fkin[k][j+1][i][l] - data->Fkin[k][j-1][i][l])/(grid->x[0][i]*(grid->x[1][j+1] - grid->x[1][j-1])))/data->p_grid[l];
#endif
#if INCLUDE_KDIR
            data->Jkin3[k][j][i][l] = - (D*(data->Fkin[k+1][j][i][l] - data->Fkin[k-1][j][i][l])/(grid->x[2][k+1] - grid->x[2][k-1]))/data->p_grid[l];
#endif
#elif GEOMETRY == SPHERICAL
#if INCLUDE_IDIR
            data->Jkin1[k][j][i][l] = - (D*(data->Fkin[k][j][i+1][l] - data->Fkin[k][j][i-1][l])/(grid->x[0][i+1] - grid->x[0][i-1]))/data->p_grid[l];
#endif
#if INCLUDE_JDIR
            data->Jkin2[k][j][i][l] = - (D*(data->Fkin[k][j+1][i][l] - data->Fkin[k][j-1][i][l])/(grid->x[0][i]*(grid->x[1][j+1] - grid->x[1][j-1])))/data->p_grid[l];
#endif
#if INCLUDE_KDIR
            data->Jkin3[k][j][i][l] = - (D*(data->Fkin[k+1][j][i][l] - data->Fkin[k-1][j][i][l])/(grid->x[0][i]*sin(grid->x[1][j])*(grid->x[2][k+1] - grid->x[2][k-1])))/data->p_grid[l];
#endif
#endif
        }
    }

    DEBUG_FUNC_END ("KIN_Update");
}

void setBoundaryRightPartToZero(Data* data, Grid* grid){
        int k,j,i;
        TOT_LOOP(k,j,i){
            int maxNU = NMOMENTUM;

    #if INCLUDE_IDIR
            if(grid->lbound[0] != 0){
                if(grid->lbound[0] != PERIODIC){
                if(i == IBEG){
                        for(int l = 0; l < maxNU; ++l){
                            data->rightPart[k][j][i][l] = 0;
                        }
                }
                }
            }

            if(grid->rbound[0] != 0){
                if(grid->rbound[0] != PERIODIC){
                if(i == IEND){
                        for(int l = 0; l < maxNU; ++l){
                            data->rightPart[k][j][i][l] = 0;
                        }
                }
                }
            }
    #endif

    #if INCLUDE_JDIR
            if(grid->lbound[1] != 0){
                if(grid->lbound[1] != PERIODIC){
                    if(j == JBEG){
                            for(int l = 0; l < maxNU; ++l){
                                data->rightPart[k][j][i][l] = 0;
                            }
                    }
                }
            }

            if(grid->rbound[1] != 0){
                if(grid->rbound[1] != PERIODIC){
                    if(j == JEND){
                            for(int l = 0; l < maxNU; ++l){
                                data->rightPart[k][j][i][l] = 0;
                            }
                    }
                }
            }
    #endif
    }
}

void Particles_KIN_LorentzTransformation(double *u, const double *vg) {
    double v = sqrt(vg[0] * vg[0] + vg[1] * vg[1] + vg[2] * vg[2]);
    double gamma = 1.0 / sqrt(1 - (vg[0] * vg[0] + vg[1] * vg[1] + vg[2] * vg[2]) / (PARTICLES_MC_C * PARTICLES_MC_C));
    double srch = 1.0 / gamma;

    if (v <= 0) {
        return;
    }

    double u1[3];
    double u0 = sqrt(PARTICLES_MC_C * PARTICLES_MC_C + u[0] * u[0] + u[1] * u[1] + u[2] * u[2]);

    double crossProduct1[3];
    crossProduct1[0] = (u[1] * vg[2] - u[2] * vg[1])/v;
    crossProduct1[1] = (u[2] * vg[0] - u[0] * vg[2])/v;
    crossProduct1[2] = (u[0] * vg[1] - u[1] * vg[0])/v;

    double crossProduct2[3];
    crossProduct2[0] = (crossProduct1[1] * vg[2] - crossProduct1[2] * vg[1])/v;
    crossProduct2[1] = (crossProduct1[2] * vg[0] - crossProduct1[0] * vg[2])/v;
    crossProduct2[2] = (crossProduct1[0] * vg[1] - crossProduct1[1] * vg[0])/v;

    u1[0] = gamma * (u[0] - vg[0] * u0) + (gamma - 1.0) * crossProduct2[0];
    u1[1] = gamma * (u[1] - vg[1] * u0) + (gamma - 1.0) * crossProduct2[1];
    u1[2] = gamma * (u[2] - vg[2] * u0) + (gamma - 1.0) * crossProduct2[2];

    if(u1[0] != u1[0]){
        printLog("u1[0] = nan\n");
        QUIT_PLUTO(1);
    }

    u[0] = u1[0];
    u[1] = u1[1];
    u[2] = u1[2];
}
