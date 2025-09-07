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
    double T_amb, cs_amb, rho_amb, V_amb, B_amb, mu;
    mu = MeanMolecularWeight(us);
    g_usersTemporal = MeanMolecularWeight(us);
    double UNIT_MASS = UNIT_DENSITY*pow(UNIT_LENGTH, 3.0);
    double UNIT_TIME = UNIT_LENGTH/UNIT_VELOCITY;
    double UNIT_ENERGY = UNIT_MASS*UNIT_VELOCITY*UNIT_VELOCITY;
    double UNIT_MFIELD = (UNIT_VELOCITY*sqrt(4.0*CONST_PI*UNIT_DENSITY));

    cs_amb  = 0.1*g_inputParam[V_d] / UNIT_VELOCITY;
    double rho_ism   = g_inputParam[RHO_d]; //interstellar medium density
    B_amb   = g_inputParam[B_AMB]/UNIT_MFIELD;


    double r = x1;


    us[RHO] = rho_ism;
    us[PRS] = us[RHO]*cs_amb*cs_amb/g_gamma;
    us[VX1] = 0;
    us[VX2] = 0;
    us[VX3] = 0;

#if GEOMETRY == SPHERICAL
    us[BX1] = B_amb*cos(x2);      // Uniform MF in z direction
    us[BX2] = B_amb*sin(x2);
    us[BX3] = 0.0;
#elif GEOMETRY == CYLINDRICAL
    us[BX1] = 0;    // Uniform MF in z direction
    us[BX2] = B_amb;
    us[BX3] = 0.0;
#elif GEOMETRY == CARTESIAN
    us[BX1] = 0;      // Uniform MF in z direction
    us[BX2] = B_amb;
    us[BX3] = 0;
#endif


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
    double rho_2 = g_inputParam[RHO_d];
    double V_2 = g_inputParam[V_d]/UNIT_VELOCITY;
    double sigma = g_inputParam[SIGMA];

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

    double rho_1 = rho_2/sigma;
    double V_1 = V_2*sigma;

    double flux = rho_2*V_2;

    double p_1 = ((g_gamma - 1.0)/g_gamma)*flux*flux*((g_gamma + 1.0)/(2*rho_2*(g_gamma - 1.0)) - 1.0/(2*rho_1));
    double p_2 = p_1 + flux*flux*(1.0/rho_1 - 1.0/rho_2);

    int i,j,k;

    TOT_LOOP(k,j,i){
        bool flag = (grid->x[0][i] < (grid->xend_glob[0] + grid->xbeg_glob[0])/2.0);
        //bool flag = (grid->x[1][j] < (grid->xend_glob[1] + grid->xbeg_glob[1])/2.0);
        if(flag){
            d->Vc[RHO][k][j][i] = rho_1;
            d->Vc[PRS][k][j][i] = p_1;
            d->Vc[VX1][k][j][i] = V_1;
            d->Vc[VX2][k][j][i] = 0;
            d->Vc[VX3][k][j][i] = 0;
            d->Vc[BX1][k][j][i] = B_amb;
            d->Vc[BX2][k][j][i] = 0;
            d->Vc[BX3][k][j][i] = 0;
            //printf("1\n");
        } else {
            d->Vc[RHO][k][j][i] = rho_2;
            d->Vc[PRS][k][j][i] = p_2;
            d->Vc[VX1][k][j][i] = V_2;
            d->Vc[VX2][k][j][i] = 0;
            d->Vc[VX3][k][j][i] = 0;
            d->Vc[BX1][k][j][i] = B_amb;
            d->Vc[BX2][k][j][i] = 0;
            d->Vc[BX3][k][j][i] = 0;
            //print("0\n");
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
    p_grid_min = 100;
    double ratio = temp_p_min/p_grid_min;
    p_grid_max = 100000*p_grid_min;
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
