/* ///////////////////////////////////////////////////////////////////// */
/*!
 \file
 \brief Initialize particle distrbution, set B.C. and injection.
 
 This file contains routines to initialize particles on the grid,
 assign user-defined boundary conditions and inject particles.
 Particles attributes that can be set are: position, velocity, color.
 In case of evolution with spectra, the initial spectral profile is also
 prescribed here.
 
 \authors A. Mignone (andrea.mignone@unito.it)\n
          B. Vaidya (bvaidya@unito.it)\n
          D. Mukherjee
          
 \date    March 21, 2021
 */
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

/* ********************************************************************* */
void Particles_Init(Data *d, Grid *grid)
/*!
 *  Sets initial conditions on particles.
 *
 *  \param [in]    d       Pointer to the PLUTO data structure.
 *  \param [in]    grid    Pointer to the PLUTO grid structure.
 *
 *********************************************************************** */
{
    double UNIT_MASS = UNIT_DENSITY*pow(UNIT_LENGTH, 3.0);
#if PARTICLES !=  PARTICLES_KIN
  int i,j,k, np, dir, nc;
  int np_glob = RuntimeGet()->Nparticles_glob;
  int np_cell = RuntimeGet()->Nparticles_cell;
  static int first_call = 1;
  double xbeg[3], xend[3];
  static double ***w;
  Particle p;

  if (first_call){
    RandomSeed(time(NULL),0);
    first_call = 0;
  }

  if (w == NULL) {
      w  = ARRAY_BOX (-1, 1, -1, 1, -1, 1, double);
  }

/* --------------------------------------------------------------
   1. Global initialization
   -------------------------------------------------------------- */
  printf("particle init\n");
  if (np_glob > 0){
    printf("start particle global initialization np_glob = %d\n", np_glob);
    for (dir = 0; dir < 3; dir++){
      xbeg[dir] = grid->xbeg_glob[dir];
      xend[dir] = grid->xend_glob[dir];
    }

    for (np = 0; np < np_glob; np++){
      Particles_LoadUniform(np, np_glob, xbeg, xend, p.coord);
#if GEOMETRY == CYLINDRICAL || GEOMETRY == SPHERICAL
            if(p.coord[0] <= 0){
                print("Ooops, particle r[0] = %g\n",p.coord[0]);
            }
#endif
      #if (PARTICLES == PARTICLES_LP) && (PARTICLES_LP_SPECTRA == YES)
      Particles_LP_InitSpectra(&p);  
      for (nc = 0; nc < PARTICLES_LP_NCOLORS; nc++) p.color[nc] = 0.0;
      #endif
#if PARTICLES == PARTICLES_CR
            double u = 1E4;
            double theta = RandomNumber(0, CONST_PI);
            double phi = RandomNumber(0, 2*CONST_PI);
            p.speed[IDIR] = u*sin(theta)*cos(phi);
            p.speed[JDIR] = u*sin(theta)*sin(phi);
            p.speed[KDIR] = u*cos(theta);
            Particles_GetWeights(&p, p.cell, w, grid);
            double u2 = sqrt(p.speed[IDIR]*p.speed[IDIR] + p.speed[JDIR]*p.speed[JDIR] + p.speed[KDIR]*p.speed[KDIR]);

            i = p.cell[IDIR];
            j = p.cell[JDIR];
            k = p.cell[KDIR];
            //printf("k=%d, j = %d, i = %d\n");
            p.mass        = 1.e-3*d->Vc[RHO][k][j][i]*grid->dV[k][j][i];
            //p.mass        = 1E23*1.67E-24;
            p.color       = 0.0;
            //p.coord[0] = 0.0;
            //p.coord[1] = 0.0;
            //printf("3\n");
#endif
#if PARTICLES == PARTICLES_MC
            double u = 1E4;
            double theta = RandomNumber(0, CONST_PI);
            double phi = RandomNumber(0, 2*CONST_PI);
            p.speed[IDIR] = u*sin(theta)*cos(phi);
            p.speed[JDIR] = u*sin(theta)*sin(phi);
            p.speed[KDIR] = u*cos(theta);
            Particles_GetWeights(&p, p.cell, w, grid);

            i = p.cell[IDIR];
            j = p.cell[JDIR];
            k = p.cell[KDIR];
            //printf("k=%d, j = %d, i = %d\n");
            p.mass        = 1.e-3*d->Vc[RHO][k][j][i]*grid->dV[k][j][i];
            p.color       = 0.0;
            p.scattering_time = 0;
            //printf("3\n");
#endif
      Particles_Insert (&p, d, PARTICLES_CREATE, grid);
    }
    printf("finish particle global initialization\n");
  }

/* ------------------------------------------------------------------
   2. Cell by cell initialization.
      Note: You may use Particles_LoadRandom() to initialize
            velocity components but not spatial coordinates.
   ------------------------------------------------------------------ */

  if (np_cell > 0){

    DOM_LOOP(k,j,i){
      xbeg[IDIR] = grid->xl[IDIR][i]; xend[IDIR] = grid->xr[IDIR][i];
      xbeg[JDIR] = grid->xl[JDIR][j]; xend[JDIR] = grid->xr[JDIR][j];
      xbeg[KDIR] = grid->xl[KDIR][k]; xend[KDIR] = grid->xr[KDIR][k];

      for (np = 0; np < np_cell; np++){

      /* -- Spatial distribution -- */

        Particles_LoadUniform(np, np_cell, xbeg, xend, p.coord);

        #if (PARTICLES == PARTICLES_LP) && (PARTICLES_LP_SPECTRA == YES)
        Particles_LP_InitSpectra(&p);
        for (nc = 0; nc < PARTICLES_LP_NCOLORS; nc++) p.color[nc] = 0.0;
        #endif

#if PARTICLES == PARTICLES_CR
                double u = 1E4;
                double theta = RandomNumber(0, CONST_PI);
                double phi = RandomNumber(0, 2*CONST_PI);
                p.speed[IDIR] = u*sin(theta)*cos(phi);
                p.speed[JDIR] = u*sin(theta)*sin(phi);
                p.speed[KDIR] = u*cos(theta);
                Particles_GetWeights(&p, p.cell, w, grid);
                double u2 = sqrt(p.speed[IDIR]*p.speed[IDIR] + p.speed[JDIR]*p.speed[JDIR] + p.speed[KDIR]*p.speed[KDIR]);

                i = p.cell[IDIR];
                j = p.cell[JDIR];
                k = p.cell[KDIR];
                //printf("k=%d, j = %d, i = %d\n");
                p.mass        = 1.e-3*d->Vc[RHO][k][j][i]*grid->dV[k][j][i];
                //p.mass        = 1E23*1.67E-24;
                p.color       = 0.0;
                //p.coord[0] = 0.0;
                //p.coord[1] = 0.0;
                //printf("3\n");
#endif
#if PARTICLES == PARTICLES_MC
                double u = 1E4;
                double theta = RandomNumber(0, CONST_PI);
                double phi = RandomNumber(0, 2*CONST_PI);
                p.speed[IDIR] = u*sin(theta)*cos(phi);
                p.speed[JDIR] = u*sin(theta)*sin(phi);
                p.speed[KDIR] = u*cos(theta);
                Particles_GetWeights(&p, p.cell, w, grid);

                i = p.cell[IDIR];
                j = p.cell[JDIR];
                k = p.cell[KDIR];
                //printf("k=%d, j = %d, i = %d\n");
                p.mass        = 1.e-3*d->Vc[RHO][k][j][i]*grid->dV[k][j][i];
                //p.mass        = 1E23*1.67E-24;
                p.color       = 0.0;
                //p.coord[0] = 0.0;
                //p.coord[1] = 0.0;
                //printf("3\n");
                p.scattering_time = 0;
#endif
        Particles_Insert (&p, d, PARTICLES_CREATE, grid);
      }
    }
    printf("finish particle local initialization\n");
  }  
  Particles_SetID(d->PHead);
#else
    int i,j,k,l;
    TOT_LOOP(k,j,i){
        for(l = 0; l < NMOMENTUM; ++l){
        }
    }

    DOM_LOOP(k,j,i){
            d->injectedEnergy[k][j][i] = 0;
            for(int l = 0; l < NMOMENTUM; ++l){
                double E = PARTICLES_KIN_MASS*PARTICLES_KIN_C*PARTICLES_KIN_C*sqrt(1.0 + d->p_grid[l]*d->p_grid[l]);
                double dp = 0;
                if( l == 0){
                    dp = d->p_grid[1] - d->p_grid[0];
                } else {
                    dp = d->p_grid[l] - d->p_grid[l-1];
                }
                d->injectedEnergy[k][j][i] += E*d->Fkin[k][j][i][l]*grid->dV[k][j][i]*dp/d->p_grid[l];
            }
    }

#endif
    printf("finish particle initialization\n");
}

#if PARTICLES_LP_SPECTRA == YES
/* ********************************************************************** */
void Particles_LP_InitSpectra(Particle* pl)
/*!
 *  Initialize spectra for each particle (only for LAGRANGIAN).
 *  Specify here the initial distribution of N(E) with E for each particle
 *
 *  \param [in]      pl      Pointer to the Particle structure.
 * 
 ********************************************************************** */
{
  int i;
  double Emin, Emax, DeltaE, N_0, alpha1, alpha2, Eb, N1, N2;
  double lnEmin, lnEmax, dlnE, scrh;
  Emin = 1.0e-2;
  Emax = 1.0e4;
  lnEmin = log10(Emin);
  lnEmax = log10(Emax);
  dlnE = (lnEmax - lnEmin)/((double)PARTICLES_LP_NEBINS);
    
  pl->nmicro = 0.001; /* The number density of micro-particles */

/* --------------------------------------------------------
   0. Initialize structure members
   -------------------------------------------------------- */

  pl->cmp_ratio = 1.0;
  pl->shkflag = 0;
  pl->Vshk_upst[RHO] = -1.0;
  pl->Vshk_dnst[RHO] = -1.0;
  pl->cr = 0.0;

/* --------------------------------------------------------
   1. Initialize energy grid.
      Note that pl->eng[i] gives the location of the
      *left* interface of an energy bins.
      There're therefore NBINS+1 nodes in the grid.
   -------------------------------------------------------- */

  for (i = 0; i <= PARTICLES_LP_NEBINS; i++){
    scrh       = lnEmin + i*dlnE;
    pl->eng[i] = pow(10.0, scrh);

    /* Relativistic Maxwellian Distribution*/
    /* double mu = 2.0; // mu = m_0c^2/kT_e 
     double k2_mu = BesselKn(2,mu); 
     pl->chi[i] = (pl->nmicro)*(mu/k2_mu)*(pl->eng[i] + 1.0)
                  *sqrt(pl->eng[i]*(pl->eng[i] + 2.0))*exp(-mu*(pl->eng[i] + 1.0)); 
     pl->cmp_ratio = 1.0; 
     pl->shkflag = 0;
     pl->rho_min = -1.0;
     pl->rho_max = -1.0;
     pl->cr = 0.0; */

  }

/* ----------------------------------------------------------------------
   2. Initialize spectral distribution to power law.
      Chi in each bin is assigned as an average over the energy interval.
      chi_i = \int^Eh_El E^(-p) dE /(Eh - El)
   ---------------------------------------------------------------------- */
  alpha1 = 3.0;
  scrh = (pl->nmicro)/(pow(Emax,1.0-alpha1)-pow(Emin,1.0-alpha1));
  for (i = 0; i < PARTICLES_LP_NEBINS; i++){
    pl->chi[i]  = scrh*(  pow(pl->eng[i+1],1.0-alpha1) 
                        - pow(pl->eng[i],1.0-alpha1))/(pl->eng[i+1] - pl->eng[i]);
  }

}
#endif

/* ********************************************************************* */
void Particles_Inject(Data *data, Grid *grid)
/*!
 *  Inject particles as you wish.
 *
 *  \param [in]  data    Pointer to the PLUTO data structure.
 *  \param [in]  grid    Pointer to the PLUTO grid structure.
 *********************************************************************** */
{
    static int first_call = 1;
    double UNIT_MASS = UNIT_DENSITY*pow(UNIT_LENGTH, 3.0);

#if PARTICLES != PARTICLES_KIN
    Particle p;
    int i,j,k, np;
    double t, dt, v, E, gamma;
    int num_part;
    double xbeg[3], xend[3];
    double *x1 = grid->xgc[IDIR]; /* -- array pointer to x1 coordinate -- */
    double *x2 = grid->xgc[JDIR]; /* -- array pointer to x2 coordinate -- */
    double *x3 = grid->xgc[KDIR]; /* -- array pointer to x3 coordinate -- */
    int xx,yy,ii,jj;

    //printf("Particles_Inject\n");
    t = g_time;
    dt = g_dt;

    double t1, t2;
    t1 = 1000000000000.0;
    t2 = t1 + dt;

    num_part = 100;

    if ((t >= t1) && (t <= t2)) {
        print (">> Injection\n");

        if (first_call) RandomSeed(time(NULL),0);
        first_call = 0;
        /* Seed random number */

        for (xx = 0; xx < 100; xx++) {
            for (yy = 0; yy < 100; yy++) {

                p.coord[0] =  2.0 + xx*0.005;
                p.coord[1] = -1.0 + yy*0.005;
                p.coord[2] = 0.0;

                gamma = 2.4e4;
                //v = sqrt(1 - 1/gamma/gamma);

                double phi = 2*3.1415*RandomNumber(0, 1);

                p.speed[IDIR] = gamma * cos(phi);
                p.speed[JDIR] = gamma * sin(phi);
                p.speed[KDIR] = 0;

                //p.mass = 1.0*grid->dV[k][j][i];
                //p.mass        = 1.0E-26;
                //p.color = prank;

#if PARTICLES == PARTICLES_MC
                p.scattering_time = 0;
#endif

                //status = Particles_Insert(&p, data, PARTICLES_CREATE, grid);
                Particles_Insert(&p, data, PARTICLES_CREATE, grid);
                Particles_SetID(data->PHead);
            }
        }
    }
#else
    updateShockFront(data, grid);
    FlagShock(data, grid);
    int i,j,k;
    DOM_LOOP(k,j,i){
        if(!(data->flag[k][j][i] & FLAG_ENTROPY)){
            double us[256];
            int nv;
            NVAR_LOOP(nv) us[nv]=data->Vc[nv][k][j][i];
            us[RHO] = data->downstreamDensity[k][j][i];
            us[PRS] = data->downstreamPressure[k][j][i];
            us[VX1] = data->downstreamV1[k][j][i];
            us[VX2] = data->downstreamV2[k][j][i];
            us[VX3] = data->downstreamV3[k][j][i];
            double mu = MeanMolecularWeight(us)*CONST_mp/UNIT_MASS;
            double T = (data->downstreamPressure[k][j][i]*mu/data->downstreamDensity[k][j][i]);

            double compression = data->downstreamDensity[k][j][i]/data->upstreamDensity[k][j][i];
            double injection = 0;
            double xd, yd, zd;
            double xu, yu, zu;
            double nx, ny, nz;

#if GEOMETRY == CARTESIAN
            xd = data->downstreamx1[k][j][i];
            yd = data->downstreamx2[k][j][i];
            zd = data->downstreamx3[k][j][i];

            xu = data->upstreamx1[k][j][i];
            yu = data->upstreamx2[k][j][i];
            zu = data->upstreamx3[k][j][i];
#elif GEOMETRY == CYLINDRICAL
            xd = data->downstreamx1[k][j][i]*cos(d->downstreamx3[k][j][i]);
            yd = data->downstreamx1[k][j][i]*sin(d->downstreamx3[k][j][i]);
            zd = data->downstreamx2[k][j][i];

            xu = data->upstreamx1[k][j][i]*cos(d->upstreamx3[k][j][i]);
            yu = data->upstreamx1[k][j][i]*sin(d->upstreamx3[k][j][i]);
            zu = data->upstreamx2[k][j][i];

#elif GEOMETRY == POLAR
            xd = data->downstreamx1[k][j][i]*cos(d->downstreamx2[k][j][i]);
            yd = data->downstreamx1[k][j][i]*sin(d->downstreamx2[k][j][i]);
            zd = data->downstreamx3[k][j][i];

            xu = data->upstreamx1[k][j][i]*cos(d->upstreamx2[k][j][i]);
            yu = data->upstreamx1[k][j][i]*sin(d->upstreamx2[k][j][i]);
            zu = data->upstreamx3[k][j][i];
#elif GEOMETRY == SPHERICAL
            xd = data->downstreamx1[k][j][i]*sin(d->downstreamx2[k][j][i])*cos(d->downstreamx3[k][j][i]);
            yd = data->downstreamx1[k][j][i]*sin(d->downstreamx2[k][j][i])*sin(d->downstreamx3[k][j][i]);
            zd = data->downstreamx1[k][j][i]*cos(d->downstreamx2[k][j][i]);

            xu = data->upstreamx1[k][j][i]*sin(grid->x[1][upstreamj])*cos(d->upstreamx3[k][j][i]);
            yu = data->upstreamx1[k][j][i]*sin(d->upstreamx2[k][j][i])*sin(d->upstreamx3[k][j][i]);
            zu = data->upstreamx1[k][j][i]*cos(d->upstreamx2[k][j][i]);
#else
#endif
            nx = (xu - xd)/data->shockWidth[k][j][i];
            ny = (yu - yd)/data->shockWidth[k][j][i];
            nz = (zu - zd)/data->shockWidth[k][j][i];
            double n1, n2, n3;
#if GEOMETRY == CARTESIAN
            n1 = nx;
            n2 = ny;
            n3 = nz;
#elif GEOMETRY == CYLINDRICAL
            n1 = nx*cos(grid->x[2][k]) + ny*sin(grid->x[2][k]);
            n3 = ny*cos(grid->x[2][k]) - nx*sin(grid->x[2][k]);
            n2 = nz;
#elif GEOMETRY == POLAR
            n1 = nx*cos(grid->x[1][j]) + ny*sin(grid->x[1][j]);
            n2 = ny*cos(grid->x[1][j]) - nx*sin(grid->x[1][j]);
            n3 = nz;
#elif GEOMETRY == SPHERICAL
            n1 = nx*sin(grid->x[1][j])*cos(grid->x[2][k]) + ny*sin(grid->x[1][j])*sin(grid->x[2][k]) + nz*cos(grid->x[1][j]);
            n2 = nx*cos(grid->x[1][j])*cos(grid->x[2][k]) + ny*cos(grid->x[1][j])*sin(grid->x[2][k]) - nz*sin(grid->x[1][j]);
            n3 = -nx*sin(grid->x[2][k]) + ny*cos(grid->x[2][k]);
#else
            printLog("wrong geometry/n");
            QUIT_PLUTO(1);
#endif
            //to compensate that there are several cells along shock
            double localFlux = fabs(n1*grid->A[0][k][j][i]) + fabs(n2*grid->A[1][k][j][i]) + fabs(n3*grid->A[2][k][j][i]);
            if(compression > COMPRESSION_TRESHOLD){
                //injection = INJECTION_PARAMETER*((data->downstreamDensity[k][j][i]/mu)*data->velocityJump[k][j][i]/(compression - 1.0))*(grid->dV[k][j][i]/(localFlux*data->shockWidth[k][j][i]))*data->p_grid[0]/(data->p_grid[1] - data->p_grid[0]);
                injection = INJECTION_PARAMETER*((data->downstreamDensity[k][j][i]/mu)*data->velocityJump[k][j][i]/(compression - 1.0))*(1.0/data->shockWidth[k][j][i])*data->p_grid[0]/(data->p_grid[1] - data->p_grid[0]);
            }
            CheckNanOrInfinity(injection, "Injection = NaN\n");
            data->Fkin[k][j][i][0] += injection*g_dt;
            double E = PARTICLES_KIN_MASS*PARTICLES_KIN_C*PARTICLES_KIN_C*sqrt(1.0 + data->p_grid[0]*data->p_grid[0]);
            data->injectedEnergy[k][j][i] += E*g_dt*injection*grid->dV[k][j][i]*(data->p_grid[1] - data->p_grid[0])/data->p_grid[0];
        }
    }
#endif
}


/* ********************************************************************* */
void Particles_UserDefBoundary(Data *d, int side, Grid *grid)
/*
 *
 *********************************************************************** */
{
#if PARTICLES != PARTICLES_KIN
  int    dir;
  double xbeg[3], xend[3];
  particleNode *curr = d->PHead, *next;
  Particle *p;
  
  for (dir = 0; dir < 3; dir++) {
    xbeg[dir] = grid->xbeg_glob[dir];
    xend[dir] = grid->xend_glob[dir];
  }

  if (side == X1_BEG){
    dir = IDIR;
    while (curr != NULL){
      p    = &(curr->p);  
      next = curr->next; 
      if (p->coord[dir] < xbeg[dir]) Particles_Destroy (curr,d);
      curr = next;
    }
  }

  if (side == X1_END){
    dir = IDIR;
    while (curr != NULL){
      p    = &(curr->p);  
      next = curr->next; 
      if (p->coord[dir] > xend[dir]) Particles_Destroy (curr,d);
      curr = next;
    }
  }

  if (side == X2_BEG){
    dir = JDIR;
    while (curr != NULL){
      p    = &(curr->p);  
      next = curr->next; 
      if (p->coord[dir] < xbeg[dir]) Particles_Destroy (curr,d);
      curr = next;
    }
  }

  if (side == X2_END){
    dir = JDIR;
    while (curr != NULL){
      p    = &(curr->p);  
      next = curr->next; 
      if (p->coord[dir] > xend[dir]) Particles_Destroy (curr,d);
      curr = next;
    }
  }

  if (side == X3_BEG){
    dir = KDIR;
    while (curr != NULL){
      p    = &(curr->p);  
      next = curr->next; 
      if (p->coord[dir] < xbeg[dir]) Particles_Destroy (curr,d);
      curr = next;
    }
  }

  if (side == X3_END){
    dir = KDIR;
    while (curr != NULL){
      p    = &(curr->p);  
      next = curr->next; 
      if (p->coord[dir] > xend[dir]) Particles_Destroy (curr,d);
      curr = next;
    }
  }
  #endif
}
