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

static void Particles_CR_GetElectricField(Data *, Data_Arr, Data_Arr, Grid *);

#if TURBULENT_FIELD == YES
double evaluateTurbulentField(Data* data, int i, int j, int k, double* B0, double momentum_q){
    double B20 =DOT_PRODUCT(B0, B0);
    int kindex = 0;

    double B2turb = 0;
    double rg = momentum_q * PARTICLES_MC_C/sqrt(B20 + B2turb);

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

        rg = momentum_q*PARTICLES_MC_C/sqrt(B20 + B2turb);
    }
    return sqrt(B20 + B2turb);
}
#endif

/* ********************************************************************* */
void Particles_MC_Update(Data *data, timeStep *Dts, double dt, Grid *grid)
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

    double pcoord_old[3], pspeed_old[3];
    double vg[3], B[3], E[3], vfluid[256];

    double um[3], up[3], u1[3], u[3], b[3], u_old[3], u_pf[3], v_pf[3];
    double u2_pf, gamma_pf, scrh_pf, vg2, preobr_dt, gamma_vg, vp_vg;
    double gamma, gamma_old, Bperp;
    double b2, Bmag2, u2, u2_old, v2; /* Squares of vector */
    double qd[4], scrh, omL;
    double dt0, dt_half, h_half, inv_dt, omegaL, qg;
    double inv_dtL = 0.0, inv_dts = 0.0;
    double wF, wM;
    const double c2 = PARTICLES_MC_C * PARTICLES_MC_C;
    static double ****dM_tot;  /* Cumulative momentum-energy deposition */
    static double ****dM;      /* Single sub-cycle mom-en varation      */
    static double ****Fcr0;    /* CR Force at t^n (used for extrapol. with RK2) */
    static double ****emf0;    /* Convective electric field (-v X B) */
    static double ****emf_res; /* Resistive electric field */
    static double ***w;
#if PARTICLES_DEPOSIT == INTEGER
    double Cnorm[4];  /* Used only when PARTICLES_DEPOSIT == INTEGER */
    double Fcr_max[4] = {0.0};
#endif

    particleNode *CurNode;
    Particle *p;

    DEBUG_FUNC_BEG ("MC_Update");

    if (g_time < Dts->particles_tstart) return;

#if SHOW_TIMING
    clock_t clock_beg = clock(), clock0;
#endif

#if PARTICLES_MC_FEEDBACK == YES
    Particles_MC_ResetForce(data);
#endif

#if (GEOMETRY != CYLINDRICAL) && (GEOMETRY != CARTESIAN)
#error MC Particles work in CYLINDRICAL and CARTESIAN geometry only
#endif

    Boundary(data, ALL_DIR, grid);

/* --------------------------------------------------------
   0. Allocate memory
   -------------------------------------------------------- */

    if (w == NULL) {
        w = ARRAY_BOX (-1, 1, -1, 1, -1, 1, double);
        emf0 = ARRAY_4D(3, NX3_TOT, NX2_TOT, NX1_TOT, double);

#if   ((PHYSICS == MHD) && (RESISTIVITY != NO)) \
 || (PHYSICS == ResRMHD)
        emf_res = ARRAY_4D(3, NX3_TOT, NX2_TOT, NX1_TOT, double);
#endif
    }


    inv_dt = 1.e-18;  /* Used to compute inverse time step due to */
    /* gyration + zone crossing                 */
    omegaL = 1.e-18;  /* Used to compute inverse time step due to */
    /* gyration only                            */
    dt0 = dt;  /* Save initial dt */
    dt = dt / (double) Dts->Nsub_particles;
    dt_half = 0.5 * dt;
    h_half = 0.5 * dt * PARTICLES_MC_E_MC;

/* --------------------------------------------------------
   3. Start sub-cycling
   -------------------------------------------------------- */

    Dts->invDt_particles = 1.e-38;
    Dts->omega_particles = 1.e-38;

    // reset flux
    for(int l = 0; l < NMOMENTUM; ++l){
        TOT_LOOP(k,j,i){
            data->Jkin1[k][j][i][l] = 0;
            data->Jkin2[k][j][i][l] = 0;
            data->Jkin3[k][j][i][l] = 0;
        }
    }

    PARTICLES_LOOP(CurNode, data->PHead) {
        p = &(CurNode->p);
        p->scattering_time = 0;
#if GEOMETRY == CYLINDRICAL
        if(p->coord[0] <= 0){
            print("Ooops, particle r[0] = %g\n",p->coord[0]);
        }
#endif
    }

    bool propagationCompleted = false;

    while(! propagationCompleted){
        PARTICLES_LOOP(CurNode, data->PHead) {

            p = &(CurNode->p);
            bool dt_over = (p->scattering_time >= dt);
            double dt_passed = p->scattering_time;
            double dt_left = dt-p->scattering_time;
            //double fraction = 0;


            //be carefull! while
            int iterations = 0;
            while (!dt_over) {
                /*iterations++;
                if(iterations > 10000) {
                    printf("iterations = %d\n", iterations);
                }*/
                double v_old[3], v[3], r[3];

                /* -- A. Get particle 4-velocity and compute \gamma^n -- */


                double u2 = DOT_PRODUCT(p->speed, p->speed);
                gamma = sqrt(1.0 + u2 / c2);
                scrh = 1.0 / gamma;

                for (dir = 0; dir < 3; dir++) {
                    r[dir] = p->coord[dir];
                    pcoord_old[dir] = p->coord[dir];
                    u_old[dir] = p->speed[dir]; /* Compute 4-vel  */
                    v_old[dir] = u_old[dir] * scrh;
                    u[dir] = p->speed[dir];
                    u_pf[dir] = p->speed[dir];
                }
                gamma_old = gamma;
                u2_old = u2;
                v[IDIR] = u[IDIR] * scrh;
                v[JDIR] = u[JDIR] * scrh;
                v[KDIR] = u[KDIR] * scrh;


                Particles_GetWeights(p, p->cell, w, grid);
                i = p->cell[IDIR];
                j = p->cell[JDIR];
                k = p->cell[KDIR];


#if PARTICLES_SHAPE >= 0

                //todo
#if GEOMETRY == CARTESIAN
                B[IDIR] = Particles_Interpolate(data->Vc[BX1], w, p->cell);
                B[JDIR] = Particles_Interpolate(data->Vc[BX2], w, p->cell);
                B[KDIR] = Particles_Interpolate(data->Vc[BX3], w, p->cell);
#elif GEOMETRY == CYLINDRICAL
                B[IDIR] = Particles_Interpolate(data->Vc[iBR], w, p->cell);
                B[JDIR] = Particles_Interpolate(data->Vc[iBZ], w, p->cell);
                B[KDIR] = Particles_Interpolate(data->Vc[iBPHI], w, p->cell);
#else
            #error MC Particles work in CYLINDRICAL and CARTESIAN geometry only
#endif

#if GEOMETRY == CARTESIAN
                vg[IDIR] = Particles_Interpolate(data->Vc[VX1], w, p->cell);
                vg[JDIR] = Particles_Interpolate(data->Vc[VX2], w, p->cell);
                vg[KDIR] = Particles_Interpolate(data->Vc[VX3], w, p->cell);
#elif GEOMETRY == CYLINDRICAL
                vg[IDIR] = Particles_Interpolate(data->Vc[iVR], w, p->cell);
                vg[JDIR] = Particles_Interpolate(data->Vc[iVZ], w, p->cell);
                vg[KDIR] = Particles_Interpolate(data->Vc[iVPHI], w, p->cell);
#endif

                E[IDIR] = -(vg[JDIR] * B[KDIR] - vg[KDIR] * B[JDIR]);
                E[JDIR] = -(vg[KDIR] * B[IDIR] - vg[IDIR] * B[KDIR]);
                E[KDIR] = -(vg[IDIR] * B[JDIR] - vg[JDIR] * B[IDIR]);


                //i don't understand what is it
#elif PARTICLES_SHAPE == -1

                Init(vfluid, p->coord[IDIR], p->coord[JDIR], p->coord[KDIR]);
                B[IDIR] = vfluid[BX1]; B[JDIR] = vfluid[BX2]; B[KDIR] = vfluid[BX3];
                E[IDIR] = vfluid[EX1]; E[JDIR] = vfluid[EX2]; E[KDIR] = vfluid[EX3];

#endif
                Particles_MC_LorentzTransformation(u_pf, vg);
                u2_pf=u_pf[0]*u_pf[0]+u_pf[1]*u_pf[1]+u_pf[2]*u_pf[2];
                gamma_pf = sqrt(1.0 + u2_pf / c2);
                scrh_pf = 1.0 / gamma_pf;
                v_pf[IDIR] = u_pf[IDIR] * scrh_pf;
                v_pf[JDIR] = u_pf[JDIR] * scrh_pf;
                v_pf[KDIR] = u_pf[KDIR] * scrh_pf;
                vg2 = vg[0] * vg[0] + vg[1] * vg[1] + vg[2] * vg[2];
                if (vg2 == 0) {
                    gamma_vg=1;
                    preobr_dt=1;
                }  else {
                    gamma_vg=1/sqrt(1.0 - vg2 / c2);
                    vp_vg=v[0]*vg[0]+v[1]*vg[1]+v[2]*vg[2];
                    preobr_dt=gamma_vg*(1 - vp_vg / c2);
                }

                Bmag2 = DOT_PRODUCT(B, B);
                //todo initialization
                Bmag2 = 1.0E-22;

#if TURBULENT_FIELD == YES
                double momentum_q = sqrt(u2_pf)/PARTICLES_MC_E_MC;
                //todo B into pf
                double Bmag = evaluateTurbulentField(data,i,j,k, B, momentum_q);
                Bmag2 = Bmag*Bmag;
#endif
                //Bohm diffusion length
                double lambda = sqrt(u2_pf) * PARTICLES_MC_C / (PARTICLES_MC_E_MC * sqrt(Bmag2));
                double tau = lambda / sqrt(v_pf[0] * v_pf[0] + v_pf[1] * v_pf[1] + v_pf[2] * v_pf[2]);
                /* monte-carlo pusher */
                double t_step1 = 0;
                double t_cross = 0;
#if GEOMETRY == CYLINDRICAL
                double rnext = grid->xr[0][i];
                double rprev = grid->xl[0][i];

                double znext = grid->xr[1][j];
                double zprev = grid->xl[1][j];
                double vrphi = sqrt(v[0] * v[0] + v[2] * v[2]);
                double cosalpha = v[0] / vrphi;
                double lr;
                double lz;
                double invdt0 = 0, invdt1 = 0;

                if(r[0] <= 0){
                    print("Ooops, particle r[0] = %g\n",r[0]);
                }

                if (v[0] > 0) {
                    if(r[0] >= rnext){
                        double rnext = grid->xr[0][i+1];
                        double rprev = grid->xl[0][i+1];
                    }
                    lr = sqrt(r[0] * r[0] * cosalpha * cosalpha - r[0] * r[0] + rnext * rnext) - r[0] * cosalpha;
                    //if(lr <= 0){
                    //    rnext = grid->xr[0][i+1];
                    //    rprev = grid->xl[0][i+1];
                    //    lr = sqrt(r[0] * r[0] * cosalpha * cosalpha - r[0] * r[0] + rnext * rnext) - r[0] * cosalpha;
                    //}
                    if(lr/r[0] < 1E-10){
                        lr = r[0]*1E-10;
                    }
                    invdt0 = vrphi/lr;
                } else {
                    //printf("%16.20lf ",r[0]);
                    double a = r[0] - 0.6;
                    if(r[0] <= rprev){
                        rnext = grid->xr[0][i-1];
                        rprev = grid->xl[0][i-1];
                    }
                    double D = r[0] * r[0] * cosalpha * cosalpha - r[0] * r[0] + rprev * rprev;
                    if(D > 0){
                        lr = -sqrt(D) - r[0] * cosalpha;
                        if(lr/r[0] < 1E-10){
                            lr = r[0]*1E-10;
                        }
                        //if(lr <= 0){
                        //    rnext = grid->xr[0][i-1];
                        //    rprev = grid->xl[0][i-1];
                        //    lr = r[0] - rprev;
                        //}
                    } else {
                        lr = sqrt(r[0] * r[0] * cosalpha * cosalpha - r[0] * r[0] + rnext * rnext) - r[0] * cosalpha;
                    }
                    invdt0 = vrphi/lr;
                }

#if INCLUDE_JDIR
                if (v[1] > 0) {
                    lz = znext - r[1];
                    if(lz <= 0){
                        znext = grid->xr[1][j+1];
                        zprev = grid->xl[1][j+1];
                        lz = znext - r[1];
                    }
                    invdt1 = v[1]/lz;
                } else {
                    lz = r[1] - zprev;
                    if(lz <= 0){
                        znext = grid->xr[1][j-1];
                        zprev = grid->xl[1][j-1];
                        lz = r[1] - zprev;
                    }
                    invdt1 = fabs(v[1])/lz;
                }
#endif

                t_cross = 1.0/MAX(invdt0, invdt1);

#elif GEOMETRY == CARTESIAN
                double xnext = grid->xr[0][i];
                double xprev = grid->xl[0][i];

                double ynext = grid->xr[1][j];
                double yprev = grid->xl[1][j];

                double znext = grid->xr[2][k];
                double zprev = grid->xl[2][k];
                double lx, ly, lz;
                double invdt0 = 0, invdt1 = 0, invdt2 = 0;


                if (v[0] > 0) {
                    lx = xnext - r[0];
                    if(lx <= 0){
                        xnext = grid->xr[0][i+1];
                        xprev = grid->xl[0][i+1];
                        lx = xnext - r[0];
                    }
                    invdt0 = v[0]/lx;
                } else {
                    lx = r[0] - xprev;
                    if(lx <= 0){
                        xnext = grid->xr[0][i-1];
                        xprev = grid->xl[0][i-1];
                        lx = r[0] - xprev;
                    }
                    invdt0 = fabs(v[0])/lx;
                }
#if INCLUDE_JDIR
                if (v[1] > 0) {
                    ly = ynext - r[1];
                    if(ly <= 0){
                        ynext = grid->xr[1][j+1];
                        yprev = grid->xl[1][j+1];
                        ly = ynext - r[1];
                    }
                    invdt1 = v[1]/ly;
                } else{
                    ly = r[1] - yprev;
                    if(ly <= 0){
                        ynext = grid->xr[1][j-1];
                        yprev = grid->xl[1][j-1];
                        ly = r[1] - yprev;
                    }
                    invdt1 = fabs(v[1])/ly;
                }
#endif

#if INCLUDE_KDIR
                if (v[2] > 0) {
                    lz = znext - r[2];
                    if(lz <= 0){
                        znext = grid->xr[2][k+1];
                        zprev = grid->xl[2][k+1];
                        lz = znext - r[2];
                    }
                    invdt2 = v[2]/lz;
                } else {
                    lz = r[2] - zprev;
                    if(lz <= 0){
                        znext = grid->xr[2][k-1];
                        zprev = grid->xl[2][k-1];
                        lz = r[2] - zprev;
                    }
                    invdt2 = fabs(v[2])/lz;
                }
#endif
                t_cross = 1.0/MAX(invdt0, MAX(invdt1, invdt2));
                /*if(ly*fabs(v[0]) > lx*fabs(v[1])) {
                    if (lz * fabs(v[0]) > lx * fabs(v[2])) {
                        t_cross = lx / fabs(v[0]);
                    } else {
                        t_cross = lz / fabs(v[2]);
                    }
                } else {
                    if(lz*fabs(v[1]) > ly*fabs(v[2])){
                        t_cross = ly/fabs(v[1]);
                    } else {
                        t_cross = lz/fabs(v[2]);
                    }
                }*/


#endif
                //todo check, mke more clever method
                t_cross = t_cross * 1.01;
                double tau1 = tau/preobr_dt;


                if (t_cross > dt_left) {
                    if (tau1 > dt_left) {
                        dt_passed += dt_left;
                        p->scattering_time += dt_left;
                        t_step1 = dt_left;
                        dt_left = 0;
                        dt_over = true;
                    } else {
                        dt_passed += tau1;
                        p->scattering_time += tau1;
                        t_step1 = tau1;
                        dt_left -= tau1;
                    }
                } else {
                    if (tau1 > t_cross) {
                        dt_passed += t_cross;
                        p->scattering_time += t_cross;
                        t_step1 = t_cross;
                        dt_left -= t_cross;
                    } else {
                        dt_passed += tau1;
                        p->scattering_time += tau1;
                        t_step1 = tau1;
                        dt_left -= tau1;
                    }
                }
                Particles_MC_BallisticMove(p, t_step1, grid);
#if TURBULENT_FIELD == YES
                Particles_MC_addFlux(p, t_step1, dt, data, grid, w);
#endif
                double prob = 1.0 - exp(-t_step1 / tau1);
                double xi = RandomNumber(0, 1);
                if (xi < prob) {
                    Particles_MC_LargeScaleScattering(p->speed, vg);
#if PARTICLES_MC_FEEDBACK == YES
                    Particles_MC_addForce(data, grid, p);
#endif
                }
                bool leave = CheckParticleLeft(data, grid, p);
                if(leave){
                    break;
                }
            }
        }


    /* ----------------------------------------------------
       3e. Set boundary condition after deposition at
           x^{n+1/2} has been done
       ---------------------------------------------------- */

        Particles_Boundary(data, grid);

        int numberToExchange = FindNumberToExchange(data, grid);

        propagationCompleted = (numberToExchange == 0);

        Particles_BoundaryExchange(data, grid);
    }
#if SHOW_TIMING
    clock0 = clock();
#endif
    /* End loop on sub-cycles */

    Dts->invDt_particles = inv_dt;
    Dts->omega_particles = omegaL;

/* ----------------------------------------------------------
   4. Compute feedback array Fcr at t(n+1/2) needed in the
      corrector step.
   ---------------------------------------------------------- */

#if PARTICLES_MC_FEEDBACK == YES
    Particles_MC_ComputeForce4(data)
#endif


#if SHOW_TIMING
    {
      double dclock_tot;

      Dts->clock_particles = (double)(clock() - clock_beg)/CLOCKS_PER_SEC;
      dclock_tot = (double)(clock() - clock_beg)/CLOCKS_PER_SEC;

      printLog ("  Total: %f, [%8.3e per particle]\n",dclock_tot,
                                                   dclock_tot/p_nparticles);
    }
#endif
    DEBUG_FUNC_END ("CR_Update");
}

void Particles_MC_BallisticMove(Particle* p, double t_step, Grid* grid){
    const double c2 = PARTICLES_MC_C * PARTICLES_MC_C;
    double u2 = DOT_PRODUCT(p->speed, p->speed);
    double gamma = sqrt(1.0 + u2 / c2);
    double scrh = 1.0 / gamma;

    double pcoord_old[3], u[3], v[3];

    for (int dir = 0; dir < 3; dir++) {
        pcoord_old[dir] = p->coord[dir];
        u[dir] = p->speed[dir];
    }
    v[IDIR] = u[IDIR] * scrh;
    v[JDIR] = u[JDIR] * scrh;
    v[KDIR] = u[KDIR] * scrh;
#if GEOMETRY == CYLINDRICAL
    double vrphi = sqrt(v[0] * v[0] + v[2] * v[2]);
    if(vrphi != vrphi){
        printLog("vrphi = nan \n");
        QUIT_PLUTO(1);
    }
    double cosalpha = v[0] / vrphi;
    double l = vrphi*t_step;
    if(l/pcoord_old[0] < 1E-10){
        p->coord[0] = pcoord_old[0] + l*cosalpha;
    } else {
        p->coord[0] = sqrt(pcoord_old[0]*pcoord_old[0] + l*l + 2*pcoord_old[0]*l*cosalpha);
    }
    //double a = p->coord[0] - 0.6;
    if(p->coord[0] <= 0){
        print("Ooops, particle r[0] = %g\n",p->coord[0]);
    }
    p->coord[1] = pcoord_old[1] + v[1]*t_step;
    double arccos = (pcoord_old[0]*pcoord_old[0] + p->coord[0]*p->coord[0] -l*l)/(2*pcoord_old[0]*p->coord[0]);
    if((arccos > 1.0) && (arccos < 1.00000001)){
        arccos = 1.0;
    }
    if(arccos >= 1.00000001){
        printLog("acos > 1 \n");
        QUIT_PLUTO(1);
    }

    if((arccos < -1.0) && (arccos > -1.00000001)){
        arccos = -1.0;
    }
    if(arccos <= -1.00000001){
        printLog("acos < -1 \n");
        QUIT_PLUTO(1);
    }
    double deltaPhi = acos(arccos);

    double urphi = sqrt(p->speed[0]*p->speed[0] + p->speed[2]*p->speed[2]);
    double newalpha = acos(cosalpha) - deltaPhi;
    p->speed[0] = urphi*cos(newalpha);
    if(p->speed[0] != p->speed[0]){
        printLog("p->speed[0] = nan \n");
        QUIT_PLUTO(1);
    }
    p->speed[2] = urphi*sin(newalpha);

    if(v[2] < 0){
        deltaPhi = -deltaPhi;
        p->speed[2] = -p->speed[2];
    }
    p->coord[2] = pcoord_old[2] + deltaPhi;

#elif GEOMETRY == CARTESIAN
    p->coord[0] = p->coord[0] + v[0]*t_step;
    p->coord[1] = p->coord[1] + v[1]*t_step;
    p->coord[2] = p->coord[2] + v[2]*t_step;
#elif
    printLog("monte-carlo pusher works only with cylindrical and cartesian geometry \n");
    QUIT_PLUTO(1);
#endif
    int dir;
    DIM_LOOP(dir) {
      double ngh = grid->nghost[dir];
      double xL  = grid->xbeg[dir] - grid->dx[dir][0]*ngh;
      int i   = (int)((p->coord[dir] - xL)*grid->inv_dx[dir][IBEG]);

      if (i < 0 || i >= grid->np_tot[dir]){
        printLog ("! Particles_LocateCell(): particle outside the ");
        printLog ("i = %d dir = %d xp[dir] = %lf xL = %lf", i, dir, p->coord[dir], xL);
        printLog ("computational domain\n");
        QUIT_PLUTO(1);
      }

    }

    for (int dir = 0; dir < 3; dir++) {
        if (p->coord[dir] != p->coord[dir]){
            printLog ("! Particles_MC_BallisticMove(): particle (# %d) coord is NaN\n", p->id);
            Particles_Display(p);
            QUIT_PLUTO(1);
        }
    }
}


/* ********************************************************************* */
void Particles_CR_GetElectricField(Data *data, Data_Arr emfc, Data_Arr emfr,
                                   Grid *grid)
/*
 *
 * \param [in]   d       pointer to PLUTO data structure
 * \param [out]  emfc    convective electric field [only with feedback]
 * \param [out]  emfr    resistive electric field
 * \param [in]   grid    pointer to Grid structure
 *
 *********************************************************************** */
{
    int i, j, k;
    double qg, E[3], vg[3], B[3];

#if PARTICLES_MC_FEEDBACK == YES
    TOT_LOOP(k,j,i){
      vg[IDIR] = data->Vc[VX1][k][j][i]; B[IDIR] = data->Vc[BX1][k][j][i];
      vg[JDIR] = data->Vc[VX2][k][j][i]; B[JDIR] = data->Vc[BX2][k][j][i];
      vg[KDIR] = data->Vc[VX3][k][j][i]; B[KDIR] = data->Vc[BX3][k][j][i];

      emfc[IDIR][k][j][i] = -(vg[JDIR]*B[KDIR] - vg[KDIR]*B[JDIR]);
      emfc[JDIR][k][j][i] = -(vg[KDIR]*B[IDIR] - vg[IDIR]*B[KDIR]);
      emfc[KDIR][k][j][i] = -(vg[IDIR]*B[JDIR] - vg[JDIR]*B[IDIR]);

      qg = PARTICLES_MC_E_MC_GAS*data->Vc[RHO][k][j][i];
      data->Ex1[k][j][i] = emfc[IDIR][k][j][i] - data->Fcr[IDIR][k][j][i]/qg;
      data->Ex2[k][j][i] = emfc[JDIR][k][j][i] - data->Fcr[JDIR][k][j][i]/qg;
      data->Ex3[k][j][i] = emfc[KDIR][k][j][i] - data->Fcr[KDIR][k][j][i]/qg;
    }
#endif


#if PARTICLES_MC_FEEDBACK == NO

#if (PHYSICS == MHD) && (RESISTIVITY != NO)
    double ***Bx1 = data->Vc[BX1];
    double ***Bx2 = data->Vc[BX2];
    double ***Bx3 = data->Vc[BX3];
    double vgas[NVAR];

    for (k = INCLUDE_KDIR; k < NX3_TOT-INCLUDE_KDIR; k++){
    for (j = INCLUDE_JDIR; j < NX2_TOT-INCLUDE_JDIR; j++){
    for (i = INCLUDE_IDIR; i < NX1_TOT-INCLUDE_IDIR; i++){

      int nv;
      double J[3], eta[3];
      double *x1 = grid->x[IDIR], *dx1 = grid->dx[IDIR];
      double *x2 = grid->x[JDIR], *dx2 = grid->dx[JDIR];
      double *x3 = grid->x[KDIR], *dx3 = grid->dx[KDIR];

      J[IDIR] =  CDIFF_X2(Bx3,k,j,i)/dx2[j] - CDIFF_X3(Bx2,k,j,i)/dx3[k];
      J[JDIR] =  CDIFF_X3(Bx1,k,j,i)/dx3[k] - CDIFF_X1(Bx3,k,j,i)/dx1[i];
      J[KDIR] =  CDIFF_X1(Bx2,k,j,i)/dx1[i] - CDIFF_X2(Bx1,k,j,i)/dx2[j];

      NVAR_LOOP(nv) vgas[nv] = data->Vc[nv][k][j][i];

      /* -- Compute current and resistivity at cell center -- */

      Resistive_eta (vgas, x1[i], x2[j], x3[k], J, eta);
      emfr[IDIR][k][j][i] = eta[IDIR]*J[IDIR];
      emfr[JDIR][k][j][i] = eta[JDIR]*J[JDIR];
      emfr[KDIR][k][j][i] = eta[KDIR]*J[KDIR];
    }}}

  /* ------------------------------------
     2b. In parallel, fill ghost zones
         for resistive electric field
         since it is computed using
         central differences
     ------------------------------------ */

#ifdef PARALLEL
    MPI_Barrier (MPI_COMM_WORLD);
    AL_Exchange ((char *)emfr[IDIR][0][0], SZ1);
    AL_Exchange ((char *)emfr[JDIR][0][0], SZ1);
    AL_Exchange ((char *)emfr[KDIR][0][0], SZ1);
#endif

#endif  /* (PHYSICS == MHD) && (RESISTIVITY != NO) */

#if PHYSICS == ResRMHD
    TOT_LOOP(k,j,i){
      vg[IDIR] = data->Vc[VX1][k][j][i]; B[IDIR] = data->Vc[BX1][k][j][i];
      vg[JDIR] = data->Vc[VX2][k][j][i]; B[JDIR] = data->Vc[BX2][k][j][i];
      vg[KDIR] = data->Vc[VX3][k][j][i]; B[KDIR] = data->Vc[BX3][k][j][i];

      E[IDIR] = -(vg[JDIR]*B[KDIR] - vg[KDIR]*B[JDIR]);
      E[JDIR] = -(vg[KDIR]*B[IDIR] - vg[IDIR]*B[KDIR]);
      E[KDIR] = -(vg[IDIR]*B[JDIR] - vg[JDIR]*B[IDIR]);

      emfr[IDIR][k][j][i] = (data->Vc[EX1][k][j][i] - E[IDIR]);
      emfr[JDIR][k][j][i] = (data->Vc[EX2][k][j][i] - E[JDIR]);
      emfr[KDIR][k][j][i] = (data->Vc[EX3][k][j][i] - E[KDIR]);
    }
#endif

#endif  /* PARTICLES_CR_FEEDBACK == NO */

}

void Particles_MC_LargeScaleScattering(double *u, const double *vg) {
    Particles_MC_LorentzTransformation(u, vg);
    double unorm = sqrt(u[0] * u[0] + u[1] * u[1] + u[2] * u[2]);
    double phi = RandomNumber(0, 8 * atan(1.0));
    double mu = RandomNumber(-1.0, 1.0);
    double sintheta = sqrt(1.0 - mu * mu);
    u[0] = unorm * mu;
    u[1] = unorm * sintheta * cos(phi);
    u[2] = unorm * sintheta * sin(phi);
    double v1[3];
    v1[0] = -vg[0];
    v1[1] = -vg[1];
    v1[2] = -vg[2];
    Particles_MC_LorentzTransformation(u, v1);
}

void Particles_MC_LorentzTransformation(double *u, const double *vg) {
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

void    Particles_MC_addFlux(Particle* p, double t_step1, double dt, Data* data, Grid* grid, double*** w){
    int i, j, k;
    double vg[3], u_pf[3];
    Particles_GetWeights(p, p->cell, w, grid);
    i = p->cell[IDIR];
    j = p->cell[JDIR];
    k = p->cell[KDIR];

    u_pf[0] = p->speed[0];
    u_pf[1] = p->speed[1];
    u_pf[2] = p->speed[2];

#if GEOMETRY == CARTESIAN
                vg[IDIR] = Particles_Interpolate(data->Vc[VX1], w, p->cell);
                vg[JDIR] = Particles_Interpolate(data->Vc[VX2], w, p->cell);
                vg[KDIR] = Particles_Interpolate(data->Vc[VX3], w, p->cell);
#elif GEOMETRY == CYLINDRICAL
                vg[IDIR] = Particles_Interpolate(data->Vc[iVR], w, p->cell);
                vg[JDIR] = Particles_Interpolate(data->Vc[iVZ], w, p->cell);
                vg[KDIR] = Particles_Interpolate(data->Vc[iVPHI], w, p->cell);
#endif
    Particles_MC_LorentzTransformation(u_pf, vg);
    double u = sqrt(u_pf[0]*u_pf[0] + u_pf[1]*u_pf[1] + u_pf[2]*u_pf[2]);
    int pindex = 0;
    if(u < p_grid_min){
        return;
    } else if(u > p_grid_max){
        return;
    } else {
        double dlog = (log(p_grid_max) - log(p_grid_min))/NMOMENTUM;
        pindex = (log(u) - log(p_grid_min))/dlog;
    }

    double v1 = (p->speed[0]/sqrt(1+u*u))*(CONST_c/UNIT_VELOCITY);
    double v2 = (p->speed[1]/sqrt(1 +u*u))*(CONST_c/UNIT_VELOCITY);
    double v3 = (p->speed[2]/sqrt(1+u*u))*(CONST_c/UNIT_VELOCITY);

    double q = p->mass*PARTICLES_MC_E_MC*PARTICLES_MC_C/grid->dV[k][j][i];


    for (int kk = -INCLUDE_KDIR; kk <= INCLUDE_KDIR; kk++) {
    for (int jj = -INCLUDE_JDIR; jj <= INCLUDE_JDIR; jj++) {
    for (int ii = -INCLUDE_IDIR; ii <= INCLUDE_IDIR; ii++) {
        data->Jkin1[k+kk][j+jj][i+ii][pindex] += q*v1*w[kk][jj][ii]*t_step1/dt;
        data->Jkin2[k+kk][j+jj][i+ii][pindex] += q*v2*w[kk][jj][ii]*t_step1/dt;
        data->Jkin3[k+kk][j+jj][i+ii][pindex] += q*v3*w[kk][jj][ii]*t_step1/dt;
    }}}
}
