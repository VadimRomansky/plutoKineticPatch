/* ///////////////////////////////////////////////////////////////////// */
#include "math.h"
#include "string.h"
#include <stdbool.h>

#include "pluto.h"
#include "time.h"

double g_usersTemporal;
const double SNRradius = 0.3E19/UNIT_LENGTH;

/* ********************************************************************* */
void Init (double *us, double x1, double x2, double x3)
/*
 *
 *********************************************************************** */
{
    g_maxCoolingRate = 0.99;
    double mu;
    mu = MeanMolecularWeight(us);
    g_usersTemporal = MeanMolecularWeight(us);

    double UNIT_MASS = UNIT_DENSITY*pow(UNIT_LENGTH, 3.0);
    double UNIT_TIME = UNIT_LENGTH/UNIT_VELOCITY;
    double UNIT_ENERGY = UNIT_MASS*UNIT_VELOCITY*UNIT_VELOCITY;
    double UNIT_MFIELD = (UNIT_VELOCITY*sqrt(4.0*CONST_PI*UNIT_DENSITY));

    double B_amb   = g_inputParam[B_AMB]/UNIT_MFIELD;
    double rho = g_inputParam[RHO];
    double T = g_inputParam[T_AMB];
    double T1 = g_inputParam[T_COR];

    double p = (rho*1.38E-16*T/mu)/UNIT_ENERGY;


    double rmax = 1.0;

    us[RHO] = rho;
    us[PRS] = p;
    us[VX1] = 0;
    us[VX2] = 0;
    us[VX3] = 0;
    us[BX1] = 0;
    us[BX2] = 0;
    us[BX3] = B_amb;


    us[TRC] = 0.0;        // Zero valued scalar tracer
    g_minCoolingTemp = 100.0; // Low temperature floor
    g_smallPressure = 1.0e-6; // Low pressure floor
    g_usersTemporal = MeanMolecularWeight(us);

}

/* ********************************************************************* */
void InitDomain (Data *d, Grid *grid)
/*!
 * Assign initial condition by looping over the computational domain.
 * Called after the usual Init() function to assign initial conditions
 * on primitive variables.
 * Value assigned here will overwrite those prescribed during Init().
 *
 *
 *********************************************************************** */
{
    double UNIT_MASS = UNIT_DENSITY*pow(UNIT_LENGTH, 3.0);
    double UNIT_TIME = UNIT_LENGTH/UNIT_VELOCITY;
    double UNIT_ENERGY = UNIT_MASS*UNIT_VELOCITY*UNIT_VELOCITY;
    double UNIT_MFIELD = (UNIT_VELOCITY*sqrt(4.0*CONST_PI*UNIT_DENSITY));

    double B_amb   = g_inputParam[B_AMB]/UNIT_MFIELD;
    double rho = g_inputParam[RHO];
    double T = g_inputParam[T_AMB];
    double T1 = g_inputParam[T_COR];


    printf("NX1 = %d\n", NX1);
    printf("NX1_TOT = %d\n", NX1_TOT);
    printf("IBEG = %d\n", IBEG);
    printf("IEND = %d\n", IEND);

    printf("gbeg = %d\n", grid->gbeg[0]);
    printf("beg = %d\n", grid->beg[0]);
    printf("lbeg = %d\n", grid->lbeg[0]);

    int rank = 0;
    int nprocs;
#ifdef PARALLEL
    MPI_Comm cartComm = MPI_COMM_WORLD;
    MPI_Comm_size(cartComm, &nprocs);
    MPI_Comm_rank(cartComm, &rank);
    printf("rank = %d\n", rank);

    MPI_Barrier(cartComm);
#endif

    int i,j,k;

    int Nstar = 5;

    double xc[5];
    double yc[5];
    xc[0] = 0.4*grid->xend_glob[0] + 0.6*grid->xbeg_glob[0];
    yc[0] = 0.45*grid->xend_glob[1] + 0.55*grid->xbeg_glob[1];

    xc[1] = 0.55*grid->xend_glob[0] + 0.45*grid->xbeg_glob[0];
    yc[1] = 0.4*grid->xend_glob[1] + 0.5*grid->xbeg_glob[1];

    xc[2] = 0.2*grid->xend_glob[0] + 0.8*grid->xbeg_glob[0];
    yc[2] = 0.7*grid->xend_glob[1] + 0.3*grid->xbeg_glob[1];

    xc[3] = 0.7*grid->xend_glob[0] + 0.3*grid->xbeg_glob[0];
    yc[3] = 0.55*grid->xend_glob[1] + 0.45*grid->xbeg_glob[1];

    xc[4] = 0.35*grid->xend_glob[0] + 0.65*grid->xbeg_glob[0];
    yc[4] = 0.2*grid->xend_glob[1] + 0.8*grid->xbeg_glob[1];

    double rmax = 50.0;

    printf("prs = %d\n", PRS);

    /*for(int i = 0; i < 5000000; ++i){
        printf("%d\n",i);
    }*/

    TOT_LOOP(k,j,i){
        //if(grid->x[0][i] < (grid->xend_glob[0] + grid->xbeg_glob[0])/2.0){

            double mu = MeanMolecularWeight(d->Uc[k][j][i]);

            double p = (rho*1.38E-16*T/((CONST_mp/UNIT_MASS)*mu))/UNIT_ENERGY;

            d->Vc[RHO][k][j][i] = rho;
            d->Vc[PRS][k][j][i] = p;
            d->Vc[VX1][k][j][i] = 0;
            d->Vc[VX2][k][j][i] = 0;
            d->Vc[VX3][k][j][i] = 0;
            d->Vc[BX1][k][j][i] = 0;
            d->Vc[BX2][k][j][i] = 0;
            d->Vc[BX3][k][j][i] = B_amb;

            for(int l = 0; l < Nstar; ++l){
            double r = sqrt((grid->x[0][i] - xc[l])*(grid->x[0][i] - xc[l]) + (grid->x[1][j] - yc[l])*(grid->x[1][j] - yc[l]));
            if(r < rmax){
                d->Vc[PRS][k][j][i] = p*T1/T;
                //d->Vc[VX1][k][j][i] = 0.1;
            }
            }

    }

#if TURBULENT_FIELD == YES
    k_turb_min = 2*CONST_PI/grid->dl_min[0];
    k_turb_max = 100*k_turb_min;
    double factor = pow(k_turb_max/k_turb_min, 1.0/(NTURB - 1.0));
    d->k_turb[0] = k_turb_min;
    for(int i = 1; i < NTURB; ++i){
        d->k_turb[i] = d->k_turb[i-1]*factor;
    }

    double A = B_amb*B_amb/(12*CONST_PI*(pow(k_turb_min, -2.0/3.0) - pow(k_turb_max, -2.0/3.0)));

    for(int l = 0; l < NTURB; ++l){
        double W = A*pow(d->k_turb[l], -5.0/3.0);
        if(W < 0){
            printLog("W turbuelnt < 0, %g\n", W);
            QUIT_PLUTO(0);
        }
        TOT_LOOP(k,j,i){
            d->Wt[k][j][i][l] = W;
            d->turbulent_rightPart[k][j][i][l] = W;
        }
    }
#endif

#if (PARTICLES == PARTICLES_KIN) || (TURBULENT_FIELD == YES)
    double emc = PARTICLES_MC_E_MC;
    double c = PARTICLES_MC_C;
    double temp_p_min = 0.01*(grid->dl_min[0]*UNIT_LENGTH)*g_inputParam[B_AMB]*CONST_e/(CONST_mp*CONST_c*CONST_c);
    p_grid_min = P_GRID_MIN;
    double ratio = temp_p_min/p_grid_min;
    p_grid_max = P_GRID_MAX;
    double factor1 = pow(p_grid_max/p_grid_min, 1.0/(NMOMENTUM - 1.0));
    d->p_grid[0] = p_grid_min;
    for(i = 1; i < NMOMENTUM; ++i){
        d->p_grid[i] = d->p_grid[i-1]*factor1;
    }

    for(int l = 0; l < NMOMENTUM; ++l){
        TOT_LOOP(k,j,i){
            d->Jkin1[k][j][i][l] = 0;
            d->Jkin2[k][j][i][l] = 0;
            d->Jkin3[k][j][i][l] = 0;
        }
    }
#endif

#if PARTICLES == PARTICLES_KIN
    for(int l = 0; l < NMOMENTUM; ++l){
        TOT_LOOP(k,j,i){
            d->Fkin[k][j][i][l] = 0.0;
            d->rightPart[k][j][i][l] = 0.0;
            d->Pkin[k][j][i] = 0.0;
            d->injectedEnergy[k][j][i] = 0.0;
        }
    }
#endif

    TOT_LOOP(k,j,i){
        d->shockWidth[k][j][i] = grid->dx[0][i];
        d->velocityJump[k][j][i] = 0.0;
        d->upstreamDensity[k][j][i] = d->Vc[RHO][k][j][i];
        d->downstreamDensity[k][j][i] = d->Vc[RHO][k][j][i];
        d->upstreamPressure[k][j][i] = d->Vc[PRS][k][j][i];
        d->downstreamPressure[k][j][i] = d->Vc[PRS][k][j][i];
        d->upstreamx1[k][j][i] = grid->x[0][i];
        d->upstreamx2[k][j][i] = grid->x[1][j];
        d->upstreamx3[k][j][i] = grid->x[2][k];
        d->downstreamx1[k][j][i] = grid->x[0][i];
        d->downstreamx2[k][j][i] = grid->x[1][j];
        d->downstreamx3[k][j][i] = grid->x[2][k];
        d->upstreamV1[k][j][i] = 0;
        d->upstreamV2[k][j][i] = 0;
        d->upstreamV3[k][j][i] = 0;
        d->downstreamV1[k][j][i] = 0;
        d->downstreamV2[k][j][i] = 0;
        d->downstreamV3[k][j][i] = 0;
    }


    //initTurbulence(d, grid, 0.1);
}

/* **************************************************************** */
void Analysis (const Data *d, Grid *grid)
/*
 *
 ****************************************************************** */
{
}

/* ********************************************************************* */
void UserDefBoundary (const Data *d, RBox * box, int side, Grid *grid)
/*
 * Sets inflow boundary condition at the top boundary (side == X2_END)
 * and the stellar wind region when side == 0.
 *
 *********************************************************************** */
{



}
