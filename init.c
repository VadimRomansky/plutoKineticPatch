/* ///////////////////////////////////////////////////////////////////// */
#include "math.h"
#include "string.h"
#include "pluto.h"
#include "time.h"

double g_usersTemporal;
const double SNRradius = 1.0E19/UNIT_LENGTH;

void BodyForceVector(double *v, double *g,
                     double x1, double x2, double x3, int i, int j, int k, const Grid* grid, const Data* data){
    g[IDIR] = 0;
    g[JDIR] = 0;
    g[KDIR] = 0;
    //printf("body force\n");
#if PARTICLES == PARTICLES_KIN
#if GEOMETRY == CARTESIAN
#if INCLUDE_IDIR
    g[IDIR] = -(data->Pkin[k][j][i+1] - data->Pkin[k][j][i-1])/(grid->x[IDIR][i+1] - grid->x[IDIR][i-1]);
#endif
#if INCLUDE_JDIR
    g[JDIR] = -(data->Pkin[k][j+1][i] - data->Pkin[k][j-1][i])/(grid->x[JDIR][j+1] - grid->x[JDIR][j-1]);
#endif
#if INCLUDE_KDIR
    g[KDIR] = -(data->Pkin[k+1][j][i] - data->Pkin[k-1][j][i])/(grid->x[KDIR][k+1] - grid->x[KDIR][k-1]);
#endif
#elif GEOMETRY == CYLINDRICAL
#if INCLUDE_IDIR
    g[IDIR] = -(data->Pkin[k][j][i+1] - data->Pkin[k][j][i-1])/(grid->x[IDIR][i+1] - grid->x[IDIR][i-1]);
#endif
#if INCLUDE_JDIR
    g[JDIR] = -(data->Pkin[k][j+1][i] - data->Pkin[k][j-1][i])/(grid->x[JDIR][j+1] - grid->x[JDIR][j-1]);
#endif
#if INCLUDE_KDIR
    g[KDIR] = -(data->Pkin[k+1][j][i] - data->Pkin[k-1][j][i])/(grid->x[IDIR][i]*(grid->x[KDIR][k+1] - grid->x[KDIR][k-1]));
#endif
#elif GEOMETRY == POLAR
#if INCLUDE_IDIR
    g[IDIR] = -(data->Pkin[k][j][i+1] - data->Pkin[k][j][i-1])/(grid->x[IDIR][i+1] - grid->x[IDIR][i-1]);
#endif
#if INCLUDE_JDIR
    g[JDIR] = -(data->Pkin[k][j+1][i] - data->Pkin[k][j-1][i])/(grid->x[IDIR][i]*(grid->x[JDIR][j+1] - grid->x[JDIR][j-1]));
#endif
#if INCLUDE_KDIR
    g[KDIR] = -(data->Pkin[k+1][j][i] - data->Pkin[k-1][j][i])/(grid->x[KDIR][k+1] - grid->x[KDIR][k-1]);
#endif
#elif GEOMETRY == SPHERICAL
#if INCLUDE_IDIR
    g[IDIR] = -(data->Pkin[k][j][i+1] - data->Pkin[k][j][i-1])/(grid->x[IDIR][i+1] - grid->x[IDIR][i-1]);
#endif
#if INCLUDE_JDIR
    g[JDIR] = -(data->Pkin[k][j+1][i] - data->Pkin[k][j-1][i])/(grid->x[IDIR][i]*(grid->x[JDIR][j+1] - grid->x[JDIR][j-1]));
#endif
#if INCLUDE_KDIR
    g[KDIR] = -(data->Pkin[k+1][j][i] - data->Pkin[k-1][j][i])/(grid->x[IDIR][i]*sin(grid->x[JDIR][j])*(grid->x[KDIR][k+1] - grid->x[KDIR][k-1]));
#endif
#endif
#endif
    g[IDIR] = g[IDIR]/v[RHO];
    g[JDIR] = g[JDIR]/v[RHO];
    g[KDIR] = g[KDIR]/v[RHO];
}

void convertIntToString(char* result, int a) {
    /*for(int i = 0; i < 4; ++i){
        result[i] = '0';
    }*/
    if (a == 0) {
        return;
    }
    if(a < 0){
        printLog("int to str a < 0\n");
        QUIT_PLUTO(1);
    }
    if (a > 0) {
        int n = 0;
        while (a > 0) {
            int last = a % 10;
            a = a / 10;
            char c = last + '0';
            result[3-n] = c;
            n++;
        }
        return;
    }
}

const double lazzatiMaxGB = 39.6357;
double getLazzatiDistribution(double beta){
//lazzati table
//what is this table?
/*0.108017	0.145412
0.119277	0.137471
0.14071	        0.122866
0.160599	0.111887
0.179303	0.10189
0.206914	0.0877181
0.233572	0.0783987
0.287959	0.0603237
0.351119	0.0455551
0.400747	0.0422684
0.45739	        0.0392189
0.510663	0.0377777
0.622671	0.0350522
0.750928	0.0285289
0.936035	0.0215443*/

/*0.102228	0.152885
0.209206	0.0888524
0.284804	0.0610803
0.383471	0.0435947
0.636543	0.034846
0.967489	0.0226488
1.51991	0.0136568
2.25978	0.00838971
5.1632	0.00419689
8.20117	0.00170591
11.4134	6.802582e-004
17.3474	3.876769e-004
24.68	2.049093e-004
29.1149	1.123946e-004
39.6357	1.016398e-005*/

const int Nlazzati = 15;
double lazzatiGB[Nlazzati+1];
lazzatiGB[0] = 0;
lazzatiGB[1] = 0.102228;
lazzatiGB[2] = 0.209206;
lazzatiGB[3] = 0.284804;
lazzatiGB[4] = 0.383471;
lazzatiGB[5] = 0.636543;
lazzatiGB[6] = 0.967489;
lazzatiGB[7] = 1.51991;
lazzatiGB[8] = 2.25978;
lazzatiGB[9] = 5.1632;
lazzatiGB[10] = 8.20117;
lazzatiGB[11] = 11.4134;
lazzatiGB[12] = 17.3474;
lazzatiGB[13] = 24.68;
lazzatiGB[14] = 29.1149;
lazzatiGB[15] = 39.6357;
double lazzatiF[Nlazzati+1];
lazzatiF[0] = 1.0;
lazzatiF[1] = 0.152885;
lazzatiF[2] = 0.0888524;
lazzatiF[3] = 0.0610803;
lazzatiF[4] = 0.0435947;
lazzatiF[5] = 0.034846;
lazzatiF[6] = 0.0226488;
lazzatiF[7] = 0.0136568;
lazzatiF[8] = 0.00838971;
lazzatiF[9] = 0.00419689;
lazzatiF[10] = 0.00170591;
lazzatiF[11] = 6.802582e-4;
lazzatiF[12] = 3.876769e-4;
lazzatiF[13] = 2.049093e-4;
lazzatiF[14] = 1.123946e-4;
lazzatiF[15] = 1.016398e-5;
double lazzatiBeta[Nlazzati+1];
double lazzatidEdBeta[Nlazzati+1];

for(int i = 0; i < Nlazzati+1; ++i){
    lazzatiBeta[i] = lazzatiGB[i]/sqrt(1 + lazzatiGB[i]*lazzatiGB[i]);
}
for(int i = 0; i < Nlazzati; ++i){
    lazzatidEdBeta[i] = -(lazzatiF[i+1] - lazzatiF[i])/(lazzatiBeta[i+1] - lazzatiBeta[i]);
}

double int0 = lazzatiF[0] - lazzatiF[1];
double zeroKoef = 3*int0/(lazzatiBeta[1]*lazzatiBeta[1]*lazzatiBeta[1]);
lazzatidEdBeta[Nlazzati] = 0;
if(beta > lazzatiBeta[Nlazzati+1]){
return 0;
}
if(beta < lazzatiBeta[1]){
    return zeroKoef*beta*beta;
}
for(int i = 0; i < Nlazzati + 1; ++i){
    if(lazzatiBeta[i] > beta){
        return (lazzatidEdBeta[i-1]*(lazzatiBeta[i] - beta) + lazzatidEdBeta[i]*(beta - lazzatiBeta[i-1]))/(lazzatiBeta[i] - lazzatiBeta[i-1]);
}
}
}


/* ********************************************************************* */
void Init (double *us, double x1, double x2, double x3)
/*
 *********************************************************************** */
{
    //InitLazzati(us, x1, x2, x3);
    InitPower(us, x1, x2, x3);
}

/* ********************************************************************* */
void InitLazzati (double *us, double x1, double x2, double x3)
/*
 *********************************************************************** */
{
    double T_amb, cs_amb, rho_amb, V_amb, mu;
    mu = MeanMolecularWeight(us);
    g_usersTemporal = MeanMolecularWeight(us);
    double UNIT_MASS = UNIT_DENSITY*pow(UNIT_LENGTH, 3.0);
    double UNIT_TIME = UNIT_LENGTH/UNIT_VELOCITY;
    double UNIT_ENERGY = UNIT_MASS*UNIT_VELOCITY*UNIT_VELOCITY;

    T_amb   = g_inputParam[T_AMB];
    double m_dot = ((g_inputParam[M_dot] * CONST_Msun / CONST_YrSec) / UNIT_MASS) * UNIT_TIME;
    V_amb   = g_inputParam[V_AMB]   / UNIT_VELOCITY;
    cs_amb  = sqrt(g_gamma*CONST_kB*T_amb/(mu*CONST_mp)) / UNIT_VELOCITY;
    double rho_ism   = g_inputParam[RHO_ISM]; //interstellar medium density

    /* 'n' and 's1' - SNR type parameters. 'r_max' and 'v_max' - SNR blastwave initial radius and top speed respectively */
    //double UNIT_FC = UNIT_DENSITY*pow(UNIT_LENGTH/UNIT_VELOCITY,3.0);

    double Mej    = g_inputParam[M_COR] * CONST_Msun /UNIT_MASS;       /* Units [g*s3/cm3], Defines SNR's core density at time 't' */
    double Eej = g_inputParam[E_COR] / UNIT_ENERGY; /* SNR core's velocity at its outer boundary */

    double lazzatiMaxBeta = lazzatiMaxGB/sqrt(1 + lazzatiMaxGB*lazzatiMaxGB);
    double v_max = lazzatiMaxBeta*CONST_c/UNIT_VELOCITY;

    double Eej_k, Eej_t, r_max = 0.03;

    double t = r_max / v_max; // Time since SN explosion

    rho_amb = m_dot/(4*CONST_PI*V_amb*r_max*r_max);

    Eej_k = Eej;

    //printf("Energy = %g\n", Eej_k*UNIT_VELOCITY*UNIT_VELOCITY*UNIT_MASS);
    Eej_t = 0.01*Eej_k/0.99;

    //new
    //double rmax = 0.003;
    //double rhocor = Mej/(4*CONST_PI*rmax*rmax*rmax);
    //only if unit_velocity = c
    //double gam = 1.0/sqrt(1.0 - v_cor*v_cor);
    //double Eej = Mej*(gam - 1);
    //double Eej_t = 0.01*Eej/0.99;

    double M_Star_O   = 30.0;
    double M_Star_WR  = 18.0;
    double M_Star_YSG = 29.0;
    double M_Star_RSG = 22.0;

    double T_effv_O   = 36.0;
    double T_effv_WR  = 50.0;
    double T_effv_YSG = 7.50;
    double T_effv_RSG = 3.90;

    double L_Star_O   = 220000.0;
    double L_Star_WR  = 510000.0;
    double L_Star_YSG = 720000.0;
    double L_Star_RSG = 260000.0;

    double M_Loss_O   = 5.35e-6;
    double M_Loss_WR  = 6.50e-5;
    double M_Loss_YSG = 2.30e-4;
    double M_Loss_RSG = 1.75e-4;

    double OmegaS_O   = 2.75e-5;
    double OmegaS_WR  = 6.60e-7;
    double OmegaS_YSG = 1.00e-10;
    double OmegaS_RSG = 4.50e-11;

    double Gedd_O   = 0.1800;
    double Gedd_WR  = 0.4780;
    double Gedd_YSG = 0.5480;
    double Gedd_RSG = 0.3500;
    // todo how to get dr
    double dr = 2.5/100000;

    double lazzatiMassKoef = 225;
    double lazzatiMass03 = 0.09;

#if DIMENSIONS == 1 /* Suitable for 1D spherical geometry */
    double r = x1;
    if (r <= r_max){
    us[VX1] = r/t;
    double beta = us[VX1]*UNIT_VELOCITY/CONST_c;
    double gam = 1.0/sqrt(1.0 - beta*beta);
    double energy = Eej_k*getLazzatiDistribution(beta)/t;

        us[RHO] = energy/(4*CONST_PI*r*r*(gam-1.0)*(CONST_c/UNIT_VELOCITY)*(CONST_c/UNIT_VELOCITY));
        us[PRS] = us[RHO]*(g_gamma-1.0)*Eej_t/Mej;
    }else{ /* Initialization of the progenitor's CSM */
        us[RHO] = m_dot/(4*CONST_PI*V_amb*r*r);
        us[PRS] = us[RHO]*cs_amb*cs_amb/g_gamma;
        us[VX1] = V_amb;
        //SetWind(us, r, 4*M_Star_WR, T_effv_WR, L_Star_WR, 5*M_Loss_WR, OmegaS_WR, Gedd_WR);
    }
#endif
}

/* ********************************************************************* */
void InitPower (double *us, double x1, double x2, double x3)
/*
 *********************************************************************** */
{
    double T_amb, cs_amb, rho_amb, V_amb, mu, B_amb;
    mu = MeanMolecularWeight(us);
    g_usersTemporal = MeanMolecularWeight(us);
    double UNIT_MASS = UNIT_DENSITY*pow(UNIT_LENGTH, 3.0);
    double UNIT_TIME = UNIT_LENGTH/UNIT_VELOCITY;
    double UNIT_ENERGY = UNIT_MASS*UNIT_VELOCITY*UNIT_VELOCITY;
    double UNIT_MFIELD = UNIT_VELOCITY*sqrt(4.0*CONST_PI*UNIT_DENSITY);

    T_amb   = g_inputParam[T_AMB];
    double m_dot = ((g_inputParam[M_dot] * CONST_Msun / CONST_YrSec) / UNIT_MASS) * UNIT_TIME;
    V_amb   = g_inputParam[V_AMB]   / UNIT_VELOCITY;
    cs_amb  = sqrt(g_gamma*CONST_kB*T_amb/(mu*CONST_mp)) / UNIT_VELOCITY;
    double rho_ism   = g_inputParam[RHO_ISM]; //interstellar medium density
    B_amb   = g_inputParam[B_AMB]/UNIT_MFIELD;

    /* 'n' and 's1' - SNR type parameters. 'r_max' and 'v_max' - SNR blastwave initial radius and top speed respectively */
    //double UNIT_FC = UNIT_DENSITY*pow(UNIT_LENGTH/UNIT_VELOCITY,3.0);

    double Mej    = g_inputParam[M_COR] * CONST_Msun /UNIT_MASS;       /* Units [g*s3/cm3], Defines SNR's core density at time 't' */
    double Eej = g_inputParam[E_COR] / UNIT_ENERGY; /* SNR core's velocity at its outer boundary */

    double powerMaxBeta = 0.99;
    double v_max = powerMaxBeta*CONST_c/UNIT_VELOCITY;

    double Eej_k, Eej_t, r_max = 0.03;

    double t = r_max / v_max; // Time since SN explosion

    rho_amb = m_dot/(4*CONST_PI*V_amb*r_max*r_max);

    Eej_k = Eej;

    //printf("Energy = %g\n", Eej_k*UNIT_VELOCITY*UNIT_VELOCITY*UNIT_MASS);
    Eej_t = 0.01*Eej_k/0.99;

    //new
    //double rmax = 0.003;
    //double rhocor = Mej/(4*CONST_PI*rmax*rmax*rmax);
    //only if unit_velocity = c
    //double gam = 1.0/sqrt(1.0 - v_cor*v_cor);
    //double Eej = Mej*(gam - 1);
    //double Eej_t = 0.01*Eej/0.99;

    double M_Star_O   = 30.0;
    double M_Star_WR  = 18.0;
    double M_Star_YSG = 29.0;
    double M_Star_RSG = 22.0;

    double T_effv_O   = 36.0;
    double T_effv_WR  = 50.0;
    double T_effv_YSG = 7.50;
    double T_effv_RSG = 3.90;

    double L_Star_O   = 220000.0;
    double L_Star_WR  = 510000.0;
    double L_Star_YSG = 720000.0;
    double L_Star_RSG = 260000.0;

    double M_Loss_O   = 5.35e-6;
    double M_Loss_WR  = 6.50e-5;
    double M_Loss_YSG = 2.30e-4;
    double M_Loss_RSG = 1.75e-4;

    double OmegaS_O   = 2.75e-5;
    double OmegaS_WR  = 6.60e-7;
    double OmegaS_YSG = 1.00e-10;
    double OmegaS_RSG = 4.50e-11;

    double Gedd_O   = 0.1800;
    double Gedd_WR  = 0.4780;
    double Gedd_YSG = 0.5480;
    double Gedd_RSG = 0.3500;
    // todo how to get dr
    double dr = 2.5/100000;

    double powerIndex = 1.5;

#if DIMENSIONS == 1 /* Suitable for 1D spherical geometry */
    double gambeta0 = 0.1;
    double A = Eej_k*(powerIndex - 1.0)/(powerIndex*gambeta0);
    double r = x1;
    if (r <= r_max){
    us[VX1] = r/t;
    double beta = us[VX1]*UNIT_VELOCITY/CONST_c;
    double gam = 1.0/sqrt(1.0 - beta*beta);
    double energy;
    if (gam*beta < gambeta0){
        energy = A*(1 + beta/pow(1-beta*beta,1.5))/t;
    } else {
        energy = A*((1 + beta/pow(1-beta*beta, 1.5))/(pow(gam*beta/gambeta0, powerIndex)))/t;
    }

        us[RHO] = energy/(4*CONST_PI*r*r*(gam-1.0)*(CONST_c/UNIT_VELOCITY)*(CONST_c/UNIT_VELOCITY));
        us[PRS] = us[RHO]*(g_gamma-1.0)*Eej_t/Mej;
    }else{ /* Initialization of the progenitor's CSM */
        us[RHO] = m_dot/(4*CONST_PI*V_amb*r*r);
        us[PRS] = us[RHO]*cs_amb*cs_amb/g_gamma;
        us[VX1] = V_amb;
        //SetWind(us, r, 4*M_Star_WR, T_effv_WR, L_Star_WR, 5*M_Loss_WR, OmegaS_WR, Gedd_WR);
    }
#endif
    us[BX1] = B_amb;
    us[BX2] = 0.0;
    us[BX3] = 0.0;
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
    int i,j, k;

    double UNIT_MFIELD = UNIT_VELOCITY*sqrt(4.0*CONST_PI*UNIT_DENSITY);
    double UNIT_PRS    = UNIT_DENSITY*UNIT_VELOCITY*UNIT_VELOCITY;
    double UNIT_TIME = UNIT_LENGTH/UNIT_VELOCITY;
    double UNIT_MASS = UNIT_DENSITY*UNIT_LENGTH*UNIT_LENGTH*UNIT_LENGTH;
    double UNIT_ENERGY = UNIT_MASS*UNIT_VELOCITY*UNIT_VELOCITY;


double B_amb   = g_inputParam[B_AMB]/UNIT_MFIELD;

#if TURBULENT_FIELD == YES
    k_turb_min = 2000*CONST_PI/grid->dl_min[0];
    k_turb_max = 1E8*k_turb_min;
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
