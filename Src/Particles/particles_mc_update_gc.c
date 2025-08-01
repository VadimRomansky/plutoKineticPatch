/* ///////////////////////////////////////////////////////////////////// */
/*! 
  \file  
  \brief Update CR particles (without feedback) using guiding center
         approximation.
         
  GCA algorithm for updating charged particles.
  \code
   Loop on particles if first call {
    Interpolate EM fields at particle position
    Convert velocities to match GCA variables
    p->speed[0] = \gamma v  (4-velocity)
    p->speed[1] = \gamma    (Lorentz factor)
    p->speed[2] = \mu       (magnetic moment in lab frame)
   }
   Loop on particles{
    Interpolate EM fields at particle position
    Compute GC velocity and parallel acceleration dRdt[4]
    x_i(n+1/2)          = x_i(n)          + (dt/2)*dRdt_i(n)                    (position)
    u_\parallel(n+1/2)  = u_\parallel(n)  + (dt/2)*dRdt_3(n)                    (parallel velocity)
    \gamma(n+1/2)       = sqrt{1 + 2 \mu \omega (n+1/2) + u_\parallel(n+1/2)}   (Lorentz factor)
    Compute time step:
     dt = min(dt_{old}, Ncell_{max}*dx/(v))
    }
  \endcode
  
  Time step restriction is computed by requiring that no particle
  travels more than \c Nmax = \c PARTICLES_CR_NCELL_MAX zones and that
  the Larmor scale is resolved with more than 1 cycle (see
  particles\_cr\_update.c documentation for more information)
   
  \authors H. Haudemand (herveh96@hotmail.it),
           A. Mignone (mignone@to.infn.it),
           E. Puzzoni\n

  \b References

  \date   Mar 18, 2021
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

/*Use to visualize variable with "PRINT_VAR(name of variable)"*/
#define  PRINT_VAR(VAR_NAME)  PrintValue(#VAR_NAME, (VAR_NAME))
static void PrintValue(char *var_name, double var){  printLog("  %s\t=\t%8.10e\n", var_name, var);  }

#if PARTICLES_CR_GC == YES   /* Compile only if guiding center is enabled. */

#define CROSS_X1(a,b)   (a[1]*b[2] - a[2]*b[1])
#define CROSS_X2(a,b)   (a[2]*b[0] - a[0]*b[2])
#define CROSS_X3(a,b)   (a[0]*b[1] - a[1]*b[0])
#define EPS_B                   1.e-40

#ifndef PARTICLES_CR_GC_TIME_DER
  #define PARTICLES_CR_GC_TIME_DER  NO
#endif

#define gc_cERR_INVALID          1   /* typically arises when GCA fails, e.g. */
                                    /*  sqrt(<0)                             */
#define gc_cERR_DOMAIN_OVERSTEP  2   /* particle outside allowed domain */
#define gc_cERR_dRdt             3   /* Right hand side becomes superluminal */
#define gc_cERR_LARMOR_OVERSTEP  4   /* Larmor radius exceeds dx */
#define gc_cERR_FAST_VARYING     5

static void Particles_CR_setGC (Data *d, Grid *grid);
static int  Particles_CR_getGC (Particle *, Data *, double *, double, Grid *);
static void Particles_CR_GCinfo(Particle *, double *);
static void Particles_CR_destroyGCParticle (Particle *, particleNode *, Data *);

static int    gc_rkStage;
static double ***gc_cE[3];
static double ***gc_cEres[3];
static double ***gc_b[3];
static double ***gc_b_old[3];
static double ***gc_vE[3];
static double ***gc_vE_old[3];
static double ***gc_bGrad_b[3];
static double ***gc_vEGrad_b[3];
static double ***gc_bGrad_vE[3];  
static double ***gc_vEGrad_vE[3];
static double ***gc_Grad_Bf1[3];
static double ***gc_Grad_B[3];
static double ***gc_Bmf1;
static double ***gc_Bmf1_old;

/* ********************************************************************* */
void Particles_CR_Update(Data *data, timeStep *Dts, double dt, Grid *grid)
/*!
 * Update particles position and velocity by a step dt considering
 * the guiding center approximation.
 * 
 * \param [in,out]  d     Data structure (contains particles list)
 * \param [in,out]  Dts   timeStep structure
 * \param [in]      dt    time increment
 * \param [in]      grid  pointer to Grid structure
 *********************************************************************** */
{ 
  static int first_call = 1;
  static int gc_convert = 1;
  static int restart = 0;
  int dir, i, j, k, nfields = 6, collected_invalid[2] = {0};
  int err, ndestroy[8] = {0};
  double ***Fields[nfields], Interpolated_fields[nfields];
  double pcoord0[3], pspeed0[3], dRdt[6], v[3], current_v[3];
  double B[3], b[3], vE[3], Bmag, Bmag_inv, Bmag_inv2, vEmag2, vEb;
  double u[3], upar, uperp_mag2, umag2, gamma;
  double mu, scrh, k1[4], k2[4], k3[4], k4[4];
  double c2 = (PARTICLES_CR_C*PARTICLES_CR_C), inv_c2 = 1./c2;
  static double ***w;
  particleNode *CurNode;
  Particle *p;

/*
  if (first_call && g_time <= Dts->particles_tstart) gc_convert = 1;
  if (first_call && g_time >  Dts->particles_tstart) gc_convert = 0;
  first_call = 0;
*/

  if (g_time < Dts->particles_tstart) return;

/* --------------------------------------------------------
   0. Allocate memory and clear FLAG_GCA_FAILURE
   -------------------------------------------------------- */

  if (w == NULL){
    w = ARRAY_BOX (-1, 1, -1, 1, -1, 1, double);
  }

/* --------------------------------------------------------
   1a. Define and compute 3D global arrays for
      GC integration
   -------------------------------------------------------- */

  Boundary (data, ALL_DIR, grid);
  Particles_CR_setGC (data, grid);

  Fields[0] = data->Vc[BX1];
  Fields[1] = data->Vc[BX2];
  Fields[2] = data->Vc[BX3];
  Fields[3] = gc_vE[IDIR];
  Fields[4] = gc_vE[JDIR];
  Fields[5] = gc_vE[KDIR];
 
/* --------------------------------------------------------
   1b. At the very first step, we convert the particle
       velocity into parallel velocity and magnetic moment.
       These will be used onwards.
   -------------------------------------------------------- */

  double max_inv_dt  = 1.e-18;
  /* -- first_call controls inv_dt, if inv_dt < 0 then
        temporal derivatives are not computed by 
        Particles_CR_getGC -- */
  double inv_dt = 1./dt;


inv_dt = -1.0;
  if (gc_convert){
    print ("> Particles_CR_Update(): converting particle speed to momentum\n");

    PARTICLES_LOOP(CurNode, data->PHead){
      p = &(CurNode->p);
  
    /* -- A.  Compute weights and indices at time n, check if
              particle is in the very last ghost zone -- */

      if (!Particles_CheckSingle(p, 0, grid)){ /* This is unlikely */
        printLog("! Particles_CR_Update(): particle (id = %d) is in last ghost zone.\n", p->id);
        Particles_CR_destroyGCParticle(p, CurNode, data);
        continue;
      }

      Particles_GetWeights(p, p->cell, w, grid); 
      i = p->cell[IDIR];
      j = p->cell[JDIR];
      k = p->cell[KDIR];
      if (data->flag[k][j][i] & (FLAG_GCA_FAILURE)){
        printLog("! Particles_CR_Update(): particle (id = %d) close to ", p->id);
        printLog (" unphysical singularity (E(perp) > B).\n");
        Particles_CR_destroyGCParticle(p, CurNode, data);
        continue;
      }
      
    /* -- u is velocity only if PARTICLES_CR_GC_4VEL == NO,
          else is 4-velocity -- */
      u[IDIR] = p->speed[IDIR];
      u[JDIR] = p->speed[JDIR];
      u[KDIR] = p->speed[KDIR];
      umag2 = DOT_PRODUCT(u,u);
      
    /* -- B. Interpolate electromagnetic fields at time n -- */

      Particles_InterpolateArr(Fields, nfields, w, p->cell, Interpolated_fields);   
      for (dir = 0; dir < 3; dir++) B[dir] = Interpolated_fields[dir];

    /* -- C.  Compute important quantities in the lab frame -- */
    
      Bmag = sqrt(DOT_PRODUCT(B,B));

    /* -- Initial magnetic moment can't be 0 -- */
      if(Bmag == 0) { 
        printLog("! Particles_CR_Update(): initial B is zero for particle %d.\n", p->id);
        Particles_CR_destroyGCParticle(p, CurNode, data);
        continue;
      }
      
      Bmag_inv = 1./Bmag;
      b[IDIR] = B[IDIR]*Bmag_inv;  
      b[JDIR] = B[JDIR]*Bmag_inv;
      b[KDIR] = B[KDIR]*Bmag_inv;    
          
    /* -- D.  Copy upar, gamma and mu at time n into the
              p->speed array -- */

      #if PARTICLES_CR_GC_4VEL == NO
          
      gamma = 1.0/sqrt(1.0 - umag2*inv_c2);
      upar   = DOT_PRODUCT(u,b);
      uperp_mag2 = umag2 - upar*upar;
      
      mu = 0.5*gamma*gamma*uperp_mag2*Bmag_inv;
        
      p->speed[IDIR] = gamma*upar;   
      p->speed[JDIR] = gamma;           
      p->speed[KDIR] = mu;

      #elif PARTICLES_CR_GC_4VEL == YES
                  
      gamma = sqrt(1. + umag2*inv_c2);

      upar   = DOT_PRODUCT(u,b);
      uperp_mag2 = umag2 - upar*upar;
      
      mu = 0.5*uperp_mag2*Bmag_inv;
            
      p->speed[IDIR] = upar;
      p->speed[JDIR] = gamma;
      p->speed[KDIR] = mu;
      
      #endif
    }
  }
  gc_convert = 0;

/* --------------------------------------------------------
   2. Push particles.
      NOTE: we now have at disposal the particles position,
      parallel velocity and magnetic moment.
   -------------------------------------------------------- */

  PARTICLES_LOOP(CurNode, data->PHead){
    err = 0; /* 0 => no error; 1  => R_L > dx;
                          2 => \epsilon \simeq 1  */
    p = &(CurNode->p);

  /* --------------------------------------------
     2a.  Save particles position and velocity at time n,
          check if particle is near unphysical singularity
     -------------------------------------------- */
    
    for (dir = 0; dir < 3; dir++) {
      pcoord0[dir] = p->coord[dir];  
      pspeed0[dir] = p->speed[dir];       
    }

    Particles_GetWeights(p, p->cell, w, grid); 
    i = p->cell[IDIR];
    j = p->cell[JDIR];
    k = p->cell[KDIR];
    if (data->flag[k][j][i] & (FLAG_GCA_FAILURE)){
      printLog("! Particles_CR_Update(): particle (id = %d) near unphysical",p->id);
      printLog(" singularity (E_perp > B, incomplete stencil error).\n");
      Particles_CR_destroyGCParticle(p, CurNode, data);
      ndestroy[0]++;
      continue;
    }

    #if PARTICLES_CR_GC_INTEGRATOR == RK4
    
    /* --------------------------------------------
     2b.  RK4 predictor(s)
     -------------------------------------------- */
  
    gc_rkStage = 1;
    err = Particles_CR_getGC(p, data, k1, inv_dt, grid);
    if(err == gc_cERR_INVALID){
      printLog("! Particles_CR_Update(): (p->id = %d) E(perp) > B\n", p->id);
      Particles_CR_destroyGCParticle(p, CurNode, data);
      ndestroy[err]++;
      continue;
    }
    p->coord[IDIR] += 0.5*dt*k1[IDIR];
    p->coord[JDIR] += 0.5*dt*k1[JDIR];
    p->coord[KDIR] += 0.5*dt*k1[KDIR];
    p->speed[IDIR] += 0.5*dt*k1[3];

    //p->speed[JDIR] = sqrt(1. + dRdt[4] + p->speed[IDIR]*p->speed[IDIR]);    
  
    gc_rkStage = 2;
    err = Particles_CR_getGC(p, data, k2, inv_dt, grid);
    if(err == gc_cERR_INVALID){
      printLog("! Particles_CR_Update(): (p->id = %d) E(perp) > B\n", p->id);
      Particles_CR_destroyGCParticle(p, CurNode, data);
      ndestroy[err]++;
      continue;
    }
    p->coord[IDIR] = pcoord0[IDIR] + 0.5*dt*k2[IDIR];
    p->coord[JDIR] = pcoord0[JDIR] + 0.5*dt*k2[JDIR];
    p->coord[KDIR] = pcoord0[KDIR] + 0.5*dt*k2[KDIR];
    p->speed[IDIR] = pspeed0[IDIR] + 0.5*dt*k2[3];

    gc_rkStage = 3;
    err = Particles_CR_getGC(p, data, k3, inv_dt, grid);
    if(err == gc_cERR_INVALID){
      printLog("! Particles_CR_Update(): (p->id = %d) E(perp) > B\n", p->id);
      Particles_CR_destroyGCParticle(p, CurNode, data);
      ndestroy[err]++;
      continue;
    }
    p->coord[IDIR] = pcoord0[IDIR] + dt*k3[IDIR];
    p->coord[JDIR] = pcoord0[JDIR] + dt*k3[JDIR];
    p->coord[KDIR] = pcoord0[KDIR] + dt*k3[KDIR];
    p->speed[IDIR] = pspeed0[IDIR] + dt*k3[3];

    gc_rkStage = 4;
    err = Particles_CR_getGC(p, data, k4, inv_dt, grid);

  /* --------------------------------------------
     2b.  RK4 corrector
     -------------------------------------------- */
        
    p->coord[IDIR] = pcoord0[IDIR] + dt*(k1[IDIR] + 2.*k2[IDIR] + 2.*k3[IDIR] + k4[IDIR])/6.;
    p->coord[JDIR] = pcoord0[JDIR] + dt*(k1[JDIR] + 2.*k2[JDIR] + 2.*k3[JDIR] + k4[JDIR])/6.;
    p->coord[KDIR] = pcoord0[KDIR] + dt*(k1[KDIR] + 2.*k2[KDIR] + 2.*k3[KDIR] + k4[KDIR])/6.;
    p->speed[IDIR] = pspeed0[IDIR] + dt*(k1[3]    + 2.*k2[3]    + 2.*k3[3]    + k4[3])/6.;
    
    #elif PARTICLES_CR_GC_INTEGRATOR == RK2
 
  /* --------------------------------------------
     2c.  RK2 predictor
     -------------------------------------------- */
        
    gc_rkStage = 1;
    err = Particles_CR_getGC(p, data, k1, inv_dt, grid);
    if(err == gc_cERR_INVALID){
      printLog("! Particles_CR_Update(): (p->id = %d) Eperp > B.\n", p->id);
      Particles_CR_destroyGCParticle(p, CurNode, data);
      ndestroy[gc_cERR_INVALID]++;
      continue;
    }
    p->coord[IDIR] += dt*k1[IDIR];
    p->coord[JDIR] += dt*k1[JDIR];
    p->coord[KDIR] += dt*k1[KDIR];
    p->speed[IDIR] += dt*k1[3];
    
  /* --------------------------------------------
     2c.  RK2 corrector
     -------------------------------------------- */
        
    gc_rkStage = 2;
    err = Particles_CR_getGC(p, data, k2, inv_dt, grid);
    if(err == gc_cERR_INVALID){
      printLog("! Particles_CR_Update(): (p->id = %d) Eperp > B.\n", p->id);
      Particles_CR_destroyGCParticle(p, CurNode, data);
      ndestroy[gc_cERR_INVALID]++;
      continue;
    }
    p->coord[IDIR] = pcoord0[IDIR] + 0.5*dt*(k1[IDIR] + k2[IDIR]);
    p->coord[JDIR] = pcoord0[JDIR] + 0.5*dt*(k1[JDIR] + k2[JDIR]);
    p->coord[KDIR] = pcoord0[KDIR] + 0.5*dt*(k1[KDIR] + k2[KDIR]);
    p->speed[IDIR] = pspeed0[IDIR] + 0.5*dt*(k1[3]    + k2[3]);
    
    #elif PARTICLES_CR_GC_INTEGRATOR == RK_MIDPOINT
    
    /* --------------------------------------------
     2d.  RK2 midpoint predictor
     -------------------------------------------- */
     
    gc_rkStage = 1;
    err = Particles_CR_getGC(p, data, k1, inv_dt, grid);
    if (   err == gc_cERR_DOMAIN_OVERSTEP
        || err == gc_cERR_INVALID
        || err == gc_cERR_dRdt){
      Particles_CR_destroyGCParticle(p, CurNode, data);
      ndestroy[err]++;
      continue;
    }
    p->coord[IDIR] += 0.5*dt*k1[IDIR];
    p->coord[JDIR] += 0.5*dt*k1[JDIR];
    p->coord[KDIR] += 0.5*dt*k1[KDIR];
    p->speed[IDIR] += 0.5*dt*k1[3];  
 
  /* --------------------------------------------
     2d.  RK2 midpoint corrector
     -------------------------------------------- */
        
    gc_rkStage = 2;
    err = Particles_CR_getGC(p, data, k2, inv_dt, grid);
    if (   err == gc_cERR_DOMAIN_OVERSTEP
        || err == gc_cERR_INVALID
        || err == gc_cERR_dRdt){
      Particles_CR_destroyGCParticle(p, CurNode, data);
      ndestroy[err]++;
      continue;
    }

    p->coord[IDIR] = pcoord0[IDIR] + dt*k2[IDIR];
    p->coord[JDIR] = pcoord0[JDIR] + dt*k2[JDIR];
    p->coord[KDIR] = pcoord0[KDIR] + dt*k2[KDIR];
    p->speed[IDIR] = pspeed0[IDIR] + dt*k2[3];
    
    #endif

  /* --------------------------------------------
     2f. Compute time step restriction based on
         the maximum allowed distance that a
         particle can travel at its current speed:
         1/dt_1 = v^{n+1/2}/(eps * dx)
         where eps = PARTICLES_CR_NCELL_EPS
     -------------------------------------------- */
        
    for (dir = 0; dir < DIMENSIONS; dir++) {
      scrh   = fabs((p->coord[dir] - pcoord0[dir])/dt);
      scrh  /= PARTICLES_CR_NCELL_MAX*grid->dx[dir][p->cell[dir]];
      max_inv_dt = MAX(max_inv_dt,scrh);
    }

    if      (err == gc_cERR_LARMOR_OVERSTEP) collected_invalid[0] += 1;
    else if (err == gc_cERR_FAST_VARYING)    collected_invalid[1] += 1;
  }  /* End loop on particles */

  WARNING(
  /* -- Print warning for collected invalidities -- */
  if (collected_invalid[0] != 0){      
    printLog ("! Particles_CR_Update(): WARNING,"
    " gyroradiuses of %d particles are bigger than cell dimension."
    " GCA may lose validity.\n", collected_invalid[0]);
  }
  if (collected_invalid[1] != 0){      
    printLog ("! Particles_CR_Update(): WARNING,"
    " %d particles do not respect the slow varying field condition,"
    " GCA may lose validity.\n", collected_invalid[1]);
  }
  );
  
  data->particles_GC_InvalidCount = MAX(collected_invalid[0],collected_invalid[1]);
  Dts->invDt_particles = max_inv_dt;
  Dts->omega_particles = 1.e-18;  /* Does not enter in the time step computation */

/* --------------------------------------------------------
   3. Update old grid terms needed to compute temporal
   derivatives
   -------------------------------------------------------- */
  
#if PARTICLES_CR_GC_TIME_DER == YES
  TOT_LOOP(k, j, i){
    gc_b_old[IDIR][k][j][i] = gc_b[IDIR][k][j][i];
    gc_b_old[JDIR][k][j][i] = gc_b[JDIR][k][j][i];
    gc_b_old[KDIR][k][j][i] = gc_b[KDIR][k][j][i];
    
    gc_vE_old[IDIR][k][j][i] = gc_vE[IDIR][k][j][i];
    gc_vE_old[JDIR][k][j][i] = gc_vE[JDIR][k][j][i];
    gc_vE_old[KDIR][k][j][i] = gc_vE[KDIR][k][j][i];
    
    gc_Bmf1_old[k][j][i] = gc_Bmf1[k][j][i];
  }
#endif

/* --------------------------------------------------------
   4. Set boundary condition
   -------------------------------------------------------- */
  
  Particles_Boundary(data, grid);
  Particles_BoundaryExchange(data, grid);
  
/* --------------------------------------------------------
   5. Error diagnostic
   -------------------------------------------------------- */

  #ifdef PARALLEL
  int ndestroy_glob[8];
  MPI_Allreduce (ndestroy, ndestroy_glob, 8, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  for (i = 0; i < 8; i++) ndestroy[i] = ndestroy_glob[i];
  #endif
  if (ndestroy[gc_cERR_INVALID] > 0) {
    print ("! %d particles have been removed [gc_cERR_INVALID]\n",
              ndestroy[gc_cERR_INVALID]);
  }
  if (ndestroy[gc_cERR_DOMAIN_OVERSTEP] > 0) {
    print ("! %d particles have been removed [gc_cERR_DOMAIN_OVERSTEP]\n",
              ndestroy[gc_cERR_DOMAIN_OVERSTEP]);
  }
  if (ndestroy[gc_cERR_dRdt] > 0) {
    print ("! %d particles have been removed [gc_cERR_dRdt]\n",
              ndestroy[gc_cERR_dRdt]);
  }

/* --------------------------------------------------------
   6.  Compute gamma_{n+1} to output p->speed[JDIR] AFTER
        boundary check. 
   -------------------------------------------------------- */

  PARTICLES_LOOP(CurNode, data->PHead){
    p = &(CurNode->p);

    Particles_GetWeights(p, p->cell, w, grid);
    Particles_InterpolateArr(Fields, nfields, w, p->cell, Interpolated_fields);   
    for(dir = 0; dir < 3; dir++){
      B[dir]  = Interpolated_fields[dir];
      vE[dir] = Interpolated_fields[dir + 3];
    }
    
    /* -- Cleaning step to ensure vE is perpendicular to B -- */

    Bmag_inv2 = 1.0/(DOT_PRODUCT(B,B));
    vEb = DOT_PRODUCT(vE,B);  
    vE[IDIR] = vE[IDIR] - vEb*B[IDIR]*Bmag_inv2;
    vE[JDIR] = vE[JDIR] - vEb*B[JDIR]*Bmag_inv2;
    vE[KDIR] = vE[KDIR] - vEb*B[KDIR]*Bmag_inv2;
    vEmag2 = DOT_PRODUCT(vE,vE); 
    
    p->speed[JDIR] = sqrt((1. 
                   + 2.*p->speed[KDIR]*sqrt(DOT_PRODUCT(B,B)*(1 - vEmag2*inv_c2)*inv_c2) 
                   + p->speed[IDIR]*p->speed[IDIR]*inv_c2)/(1. - vEmag2*inv_c2));
    
    /* -- Destroy particle if vE > 1 -- */
    if(isnan(p->speed[JDIR])){
      printLog("! Particles_CR_Update(): (p->id = %d) vE > 1\n", p->id);
      Particles_CR_destroyGCParticle(p, CurNode, data);
      continue;
    }
  }
  first_call = 0;
}

/* ********************************************************************* */
void Particles_CR_destroyGCParticle(Particle *p, particleNode *CurNode, Data *data)
/*!
 * Destroys particle at CurNode and updates CurNode to the next node
 *
 * \param  [in]      p          Pointer to PLUTO particle data structure.
 * \param  [in,out]  CurNode    Node of the particle that will be destroyed
 * \param  [in]      data       Pointer to PLUTO data structure
 *********************************************************************** */
{
  particleNode *NextNode;
  
//  printLog ("! Particles_CR_destroyGCParticle(): p->id = %d has been destroyed.\n",
//             p->id);
  NextNode = CurNode->next;
  Particles_Destroy(CurNode, data);
  CurNode = NextNode;
}

/* ********************************************************************* */
void Particles_CR_setGC (Data *data, Grid *grid)
/*!
 * Define and compute 3D arrays on the grid 
 * (later needed for particle interpolation)
 * NOTE: Vector fields are saved as static 
 * field[x/y/z][k cell number][j cell number][i cell number]
 * b is magnetic versor field, u_E is the fluid
 * velocity field, gc_Bmf1 is a scalar field
 *
 * \param [in]      d         Data structure (contains particles list)
 * \param [in]      grid      pointer to Grid structure
 * \param [out]     bGradb    (b \cdot \nabla) b vector field
 * \param [out]     gc_vEGrad_b  (vE \cdot \nabla) b vector field  
 * \param [out]     gc_bGrad_vE  (b \cdot \nabla vE) vector field
 * \param [out]     gc_vEGrad_vE (vE \cdot \nabla vE) vector field
 * \param [out]     gc_Grad_Bf1  \nabla gc_Bmf1 vector field
 * \param [out]     gc_Grad_B    \nabla B vector field, needed to
 *                            check validity of GCA
 *********************************************************************** */
{
  int i,j,k;
  double B[3], cE[3], vg[3], b[3], vE[3];
  double Eperp_mag2, Bmag, BE, Bmag2, cEmag2, Bmag_inv, f1, vEmag2;
  double inv_c2 = 1./(PARTICLES_CR_C*PARTICLES_CR_C);
  
  double *dx = grid->dx[IDIR];
  double *dy = grid->dx[JDIR];
  double *dz = grid->dx[KDIR];
  double dx_inv, dy_inv, dz_inv;

/* --------------------------------------------------------
   0. Allocate memory
   -------------------------------------------------------- */

  if (gc_bGrad_b[0] == NULL){
    gc_b[0] = ARRAY_3D(NX3_TOT, NX2_TOT, NX1_TOT, double);
    gc_b[1] = ARRAY_3D(NX3_TOT, NX2_TOT, NX1_TOT, double);
    gc_b[2] = ARRAY_3D(NX3_TOT, NX2_TOT, NX1_TOT, double);

    gc_vE[0] = ARRAY_3D(NX3_TOT, NX2_TOT, NX1_TOT, double);
    gc_vE[1] = ARRAY_3D(NX3_TOT, NX2_TOT, NX1_TOT, double);
    gc_vE[2] = ARRAY_3D(NX3_TOT, NX2_TOT, NX1_TOT, double);

    gc_bGrad_b[0] = ARRAY_3D(NX3_TOT, NX2_TOT, NX1_TOT, double);
    gc_bGrad_b[1] = ARRAY_3D(NX3_TOT, NX2_TOT, NX1_TOT, double);
    gc_bGrad_b[2] = ARRAY_3D(NX3_TOT, NX2_TOT, NX1_TOT, double);

    gc_vEGrad_b[0] = ARRAY_3D(NX3_TOT, NX2_TOT, NX1_TOT, double);
    gc_vEGrad_b[1] = ARRAY_3D(NX3_TOT, NX2_TOT, NX1_TOT, double);
    gc_vEGrad_b[2] = ARRAY_3D(NX3_TOT, NX2_TOT, NX1_TOT, double);  

    gc_bGrad_vE[0] = ARRAY_3D(NX3_TOT, NX2_TOT, NX1_TOT, double);
    gc_bGrad_vE[1] = ARRAY_3D(NX3_TOT, NX2_TOT, NX1_TOT, double);
    gc_bGrad_vE[2] = ARRAY_3D(NX3_TOT, NX2_TOT, NX1_TOT, double);

    gc_vEGrad_vE[0] = ARRAY_3D(NX3_TOT, NX2_TOT, NX1_TOT, double);
    gc_vEGrad_vE[1] = ARRAY_3D(NX3_TOT, NX2_TOT, NX1_TOT, double);
    gc_vEGrad_vE[2] = ARRAY_3D(NX3_TOT, NX2_TOT, NX1_TOT, double);

    gc_Grad_Bf1[0] = ARRAY_3D(NX3_TOT, NX2_TOT, NX1_TOT, double);
    gc_Grad_Bf1[1] = ARRAY_3D(NX3_TOT, NX2_TOT, NX1_TOT, double);
    gc_Grad_Bf1[2] = ARRAY_3D(NX3_TOT, NX2_TOT, NX1_TOT, double);
  
  /* -- Grad_B is needed to check GCA validity -- */
    gc_Grad_B[0] = ARRAY_3D(NX3_TOT, NX2_TOT, NX1_TOT, double);
    gc_Grad_B[1] = ARRAY_3D(NX3_TOT, NX2_TOT, NX1_TOT, double);
    gc_Grad_B[2] = ARRAY_3D(NX3_TOT, NX2_TOT, NX1_TOT, double);

    gc_cE[0] = ARRAY_3D(NX3_TOT, NX2_TOT, NX1_TOT, double);
    gc_cE[1] = ARRAY_3D(NX3_TOT, NX2_TOT, NX1_TOT, double);
    gc_cE[2] = ARRAY_3D(NX3_TOT, NX2_TOT, NX1_TOT, double);

    gc_cEres[0] = ARRAY_3D(NX3_TOT, NX2_TOT, NX1_TOT, double);
    gc_cEres[1] = ARRAY_3D(NX3_TOT, NX2_TOT, NX1_TOT, double);
    gc_cEres[2] = ARRAY_3D(NX3_TOT, NX2_TOT, NX1_TOT, double);

    gc_Bmf1 = ARRAY_3D(NX3_TOT, NX2_TOT, NX1_TOT, double);
  
  /* -- Step (n-1) arrays related to time-dependent terms -- */
  
    gc_vE_old[0] = ARRAY_3D(NX3_TOT, NX2_TOT, NX1_TOT, double);
    gc_vE_old[1] = ARRAY_3D(NX3_TOT, NX2_TOT, NX1_TOT, double);
    gc_vE_old[2] = ARRAY_3D(NX3_TOT, NX2_TOT, NX1_TOT, double);
    
    gc_b_old[0] = ARRAY_3D(NX3_TOT, NX2_TOT, NX1_TOT, double);
    gc_b_old[1] = ARRAY_3D(NX3_TOT, NX2_TOT, NX1_TOT, double);
    gc_b_old[2] = ARRAY_3D(NX3_TOT, NX2_TOT, NX1_TOT, double);
    
    gc_Bmf1_old = ARRAY_3D(NX3_TOT, NX2_TOT, NX1_TOT, double);
  }

/* --------------------------------------------------------
   1a. Compute time-independent plasma quantities not
       containing derivatives
   -------------------------------------------------------- */

  TOT_LOOP(k,j,i){
    
    B[IDIR]  = data->Vc[BX1][k][j][i];
    B[JDIR]  = data->Vc[BX2][k][j][i];
    B[KDIR]  = data->Vc[BX3][k][j][i];
    vg[IDIR] = data->Vc[VX1][k][j][i];
    vg[JDIR] = data->Vc[VX2][k][j][i];
    vg[KDIR] = data->Vc[VX3][k][j][i];

  /* -- Here gc_cE contains the full electric field  -- */

    #if PHYSICS == ResRMHD
    gc_cE[IDIR][k][j][i] = data->Vc[EX1][k][j][i];
    gc_cE[JDIR][k][j][i] = data->Vc[EX2][k][j][i];
    gc_cE[KDIR][k][j][i] = data->Vc[EX3][k][j][i];
    #elif (PHYSICS == MHD && RESISTIVITY != NO)
      #error GC not working yet with non-zero resistivity
    #else
    gc_cE[IDIR][k][j][i] = CROSS_X1(B,vg);
    gc_cE[JDIR][k][j][i] = CROSS_X2(B,vg);
    gc_cE[KDIR][k][j][i] = CROSS_X3(B,vg);
    #endif

  /* -- Here gc_cEres contains only the resistive terms -- */

    #if PHYSICS == ResRMHD
    gc_cEres[IDIR][k][j][i] = data->Vc[EX1][k][j][i] - CROSS_X1(B,vg);
    gc_cEres[JDIR][k][j][i] = data->Vc[EX2][k][j][i] - CROSS_X2(B,vg);
    gc_cEres[KDIR][k][j][i] = data->Vc[EX3][k][j][i] - CROSS_X3(B,vg);
    #endif

    cE[IDIR] = gc_cE[IDIR][k][j][i];
    cE[JDIR] = gc_cE[JDIR][k][j][i];
    cE[KDIR] = gc_cE[KDIR][k][j][i];
    
    cEmag2 = DOT_PRODUCT(cE,cE);
    BE    = DOT_PRODUCT(B,cE);
    Bmag  = sqrt(DOT_PRODUCT(B,B)) + EPS_B;
    Bmag_inv = 1./Bmag;
    
    gc_b[IDIR][k][j][i] = B[IDIR]*Bmag_inv;  
    gc_b[JDIR][k][j][i] = B[JDIR]*Bmag_inv;
    gc_b[KDIR][k][j][i] = B[KDIR]*Bmag_inv;
    
    b[IDIR] = gc_b[IDIR][k][j][i];
    b[JDIR] = gc_b[JDIR][k][j][i];
    b[KDIR] = gc_b[KDIR][k][j][i];
    
    vE[IDIR] = CROSS_X1(cE,b)*Bmag_inv;
    vE[JDIR] = CROSS_X2(cE,b)*Bmag_inv;
    vE[KDIR] = CROSS_X3(cE,b)*Bmag_inv;
    vEmag2   = DOT_PRODUCT(vE,vE);

    gc_vE[IDIR][k][j][i] = vE[IDIR];
    gc_vE[JDIR][k][j][i] = vE[JDIR];
    gc_vE[KDIR][k][j][i] = vE[KDIR];
    
    Eperp_mag2 = cEmag2 - BE*BE*Bmag_inv*Bmag_inv;

    f1 = sqrt(1.0 - vEmag2*inv_c2);
    gc_Bmf1[k][j][i] = Bmag*f1;
    
    if (isnan(gc_Bmf1[k][j][i])){
      int iL = (i == 0         ? 0:1);  /* Avoid extending flag outside */
      int iR = (i == NX1_TOT-1 ? 0:1);  /* total computational dom.     */
      int jL = (j == 0         ? 0:1);
      int jR = (j == NX2_TOT-1 ? 0:1);
      int kL = (k == 0         ? 0:1);
      int kR = (k == NX3_TOT-1 ? 0:1);

      DIM_EXPAND(
          data->flag[k][j][i]    |= FLAG_GCA_FAILURE;
          data->flag[k][j][i+iR] |= FLAG_GCA_FAILURE;
          data->flag[k][j][i-iL] |= FLAG_GCA_FAILURE;  ,
          data->flag[k][j-jL][i] |= FLAG_GCA_FAILURE;  
          data->flag[k][j+jR][i] |= FLAG_GCA_FAILURE;  ,
          data->flag[k-kL][j][i] |= FLAG_GCA_FAILURE;
          data->flag[k+kR][j][i] |= FLAG_GCA_FAILURE;)
      printLog ("! Particles_CR_setGC(): |vE| = %12.6e > 1\n", sqrt(vEmag2));
/*
      printLog ("                        Bmag = %12.6e\n",Bmag);
      printLog ("                        Emag = %12.6e\n",sqrt(Emag2));      
      printLog ("                        vE   = ");ShowVector(vE,3);
      double vEmag = sqrt(vEmag2);
      vE[IDIR] *= 0.999/vEmag;
      vE[JDIR] *= 0.999/vEmag;
      vE[KDIR] *= 0.999/vEmag;
*/
    }
  }

/* --------------------------------------------------------
   1b. Compute plasma quantities containing derivatives
   -------------------------------------------------------- */

  for (k = INCLUDE_KDIR; k < NX3_TOT-INCLUDE_KDIR; k++){
  for (j = INCLUDE_JDIR; j < NX2_TOT-INCLUDE_JDIR; j++){
  for (i = INCLUDE_IDIR; i < NX1_TOT-INCLUDE_IDIR; i++){
  
    dx_inv = 1./dx[i];
    dy_inv = 1./dy[j];
    dz_inv = 1./dz[k];
        
    double dBmf1_dx = CDIFF_X1(gc_Bmf1, k, j, i)*dx_inv;
    double dBmf1_dy = CDIFF_X2(gc_Bmf1, k, j, i)*dy_inv;
    double dBmf1_dz = CDIFF_X3(gc_Bmf1, k, j, i)*dz_inv;

    double dbx_dx = CDIFF_X1(gc_b[IDIR], k, j, i)*dx_inv;
    double dbx_dy = CDIFF_X2(gc_b[IDIR], k, j, i)*dy_inv;
    double dbx_dz = CDIFF_X3(gc_b[IDIR], k, j, i)*dz_inv;

    double dby_dx = CDIFF_X1(gc_b[JDIR], k, j, i)*dx_inv;
    double dby_dy = CDIFF_X2(gc_b[JDIR], k, j, i)*dy_inv;
    double dby_dz = CDIFF_X3(gc_b[JDIR], k, j, i)*dz_inv;

    double dbz_dx = CDIFF_X1(gc_b[KDIR], k, j, i)*dx_inv;
    double dbz_dy = CDIFF_X2(gc_b[KDIR], k, j, i)*dy_inv;
    double dbz_dz = CDIFF_X3(gc_b[KDIR], k, j, i)*dz_inv;

    double dux_dx = CDIFF_X1(gc_vE[IDIR], k, j, i)*dx_inv;   
    double dux_dy = CDIFF_X2(gc_vE[IDIR], k, j, i)*dy_inv;
    double dux_dz = CDIFF_X3(gc_vE[IDIR], k, j, i)*dz_inv;

    double duy_dx = CDIFF_X1(gc_vE[JDIR], k, j, i)*dx_inv; 
    double duy_dy = CDIFF_X2(gc_vE[JDIR], k, j, i)*dy_inv;
    double duy_dz = CDIFF_X3(gc_vE[JDIR], k, j, i)*dz_inv;

    double duz_dx = CDIFF_X1(gc_vE[KDIR], k, j, i)*dx_inv;
    double duz_dy = CDIFF_X2(gc_vE[KDIR], k, j, i)*dy_inv;
    double duz_dz = CDIFF_X3(gc_vE[KDIR], k, j, i)*dz_inv;
    

    gc_bGrad_b[IDIR][k][j][i] = DIM_EXPAND(  gc_b[IDIR][k][j][i]*dbx_dx,
                                           + gc_b[JDIR][k][j][i]*dbx_dy,
                                           + gc_b[KDIR][k][j][i]*dbx_dz);
                     
    gc_bGrad_b[JDIR][k][j][i] = DIM_EXPAND(  gc_b[IDIR][k][j][i]*dby_dx,
                                           + gc_b[JDIR][k][j][i]*dby_dy,
                                           + gc_b[KDIR][k][j][i]*dby_dz);

    gc_bGrad_b[KDIR][k][j][i] = DIM_EXPAND(  gc_b[IDIR][k][j][i]*dbz_dx,
                                           + gc_b[JDIR][k][j][i]*dbz_dy,
                                           + gc_b[KDIR][k][j][i]*dbz_dz);
    
    #if PARTICLES_CR_GC_DEBUG
    if (   isnan (gc_bGrad_b[IDIR][k][j][i]) 
        || isnan (gc_bGrad_b[JDIR][k][j][i]) 
        || isnan (gc_bGrad_b[KDIR][k][j][i]) ){
      printLog ("! Particles_CR_setGC(): nan in bgradB\n");
      QUIT_PLUTO(1);
    }
    #endif
        
    gc_vEGrad_b[IDIR][k][j][i] = DIM_EXPAND(  gc_vE[IDIR][k][j][i]*dbx_dx,
                                            + gc_vE[JDIR][k][j][i]*dbx_dy,
                                            + gc_vE[KDIR][k][j][i]*dbx_dz);
                     
    gc_vEGrad_b[JDIR][k][j][i] = DIM_EXPAND(  gc_vE[IDIR][k][j][i]*dby_dx,
                                            + gc_vE[JDIR][k][j][i]*dby_dy,
                                            + gc_vE[KDIR][k][j][i]*dby_dz);

    gc_vEGrad_b[KDIR][k][j][i] = DIM_EXPAND(  gc_vE[IDIR][k][j][i]*dbz_dx,
                                            + gc_vE[JDIR][k][j][i]*dbz_dy,
                                            + gc_vE[KDIR][k][j][i]*dbz_dz);

    gc_bGrad_vE[IDIR][k][j][i] = DIM_EXPAND(  gc_b[IDIR][k][j][i]*dux_dx,
                                            + gc_b[JDIR][k][j][i]*dux_dy,
                                            + gc_b[KDIR][k][j][i]*dux_dz);
                     
    gc_bGrad_vE[JDIR][k][j][i] = DIM_EXPAND(  gc_b[IDIR][k][j][i]*duy_dx,
                                         + gc_b[JDIR][k][j][i]*duy_dy,
                                         + gc_b[KDIR][k][j][i]*duy_dz);

    gc_bGrad_vE[KDIR][k][j][i] = DIM_EXPAND(  gc_b[IDIR][k][j][i]*duz_dx,
                                            + gc_b[JDIR][k][j][i]*duz_dy,
                                            + gc_b[KDIR][k][j][i]*duz_dz);

    gc_vEGrad_vE[IDIR][k][j][i] = DIM_EXPAND( gc_vE[IDIR][k][j][i]*dux_dx,
                                            + gc_vE[JDIR][k][j][i]*dux_dy,
                                            + gc_vE[KDIR][k][j][i]*dux_dz);
                     
    gc_vEGrad_vE[JDIR][k][j][i] = DIM_EXPAND( gc_vE[IDIR][k][j][i]*duy_dx,
                                            + gc_vE[JDIR][k][j][i]*duy_dy,
                                            + gc_vE[KDIR][k][j][i]*duy_dz);

    gc_vEGrad_vE[KDIR][k][j][i] = DIM_EXPAND( gc_vE[IDIR][k][j][i]*duz_dx,
                                            + gc_vE[JDIR][k][j][i]*duz_dy,
                                            + gc_vE[KDIR][k][j][i]*duz_dz);

    gc_Grad_Bf1[IDIR][k][j][i] = dBmf1_dx; 
    gc_Grad_Bf1[JDIR][k][j][i] = dBmf1_dy;
    gc_Grad_Bf1[KDIR][k][j][i] = dBmf1_dz;
        
    gc_Grad_B[IDIR][k][j][i] = CDIFF_X1(data->Vc[BX1], k, j, i)*dx_inv;
    gc_Grad_B[JDIR][k][j][i] = CDIFF_X2(data->Vc[BX2], k, j, i)*dy_inv;
    gc_Grad_B[KDIR][k][j][i] = CDIFF_X3(data->Vc[BX3], k, j, i)*dz_inv;
  }}}

}

/* ********************************************************************* */
int Particles_CR_getGC(Particle *p, Data *data, double *dRdt, double inv_dt, Grid *grid)
/*!
 * Compute the right hand side of particle equation in
 * the guiding center approximation.
 * Interpolate fields provided by Particles_CR_setGC
 * at cell location.
 * Conditions:
 * 1) Larmor radius is required to be smaller than cell dimension
 * 2) Larmor radius is required to be small compared to magnetic
 *    gradient scale
 *
 * \param [in]     p                 Pointer to PLUTO particle data structure.
 * \param [in]     data              Pointer to PLUTO data structure.
 * \param [in,out] dRdt              double array containing GC speed and parallel
 *                                      4-velocity.
 * \param [in]     inv_dt        Inverse of time step, if <0 time derivatives 
 *                               are not computed.
 * \param [in]     grid          Pointer to Grid PLUTO structure.
 *
 * \return  Return 0 on success, otherwise -1 for math errors (vE > 1, nan),
 *          1 if gyroradius > dx, 2 if weakly varying field appproximation is
 *          not met.
 *
 *********************************************************************** */
{
  int i, j, k, dir, nfields;
  int ngh = grid->nghost[IDIR];
  int err = 0;
  double B[3], cE[3], cEres[3], b[3], vE[3], sum[3], timeterms[4];
  double bdb[3], vEdb[3], bdvE[3], vEdvE[3], dBf1[3], dB[3], dvEdt[3];
  double b_old[3], vE_old[3], Bmagf1_old, BtoB_old[3];
  double Bmag, Bmag2, Bmag_inv, Bstar, cEmag2, cEpar, Bmagf1, vEmag2, dBmag;
  double omega, omega_inv, vpar, dRdtmag2, R_L;
  double gammaE, gammaE2, gamma, gamma_inv;
  double scrh;
  double upar = p->speed[IDIR];
  double mu   = p->speed[KDIR];
  double e_mc = PARTICLES_CR_E_MC, inv_c = 1./(PARTICLES_CR_C), c2 = PARTICLES_CR_C*PARTICLES_CR_C;
  double mc_e = 1./(PARTICLES_CR_E_MC), inv_c2 = 1./c2;
  static double ***w;
/*!CAMBIARE c*/

/* --------------------------------------------------------
   0. Allocate memory
   -------------------------------------------------------- */

  if (w == NULL) w = ARRAY_BOX (-1, 1, -1, 1, -1, 1, double);

/* --------------------------------------------------------
   1. Compute weights and indices, return if incomplete stencils
   -------------------------------------------------------- */

  if (!Particles_CheckSingle(p, ngh-2, grid)){ 
//    printLog("! Particles_CR_getGC() [stage = %d]: p->id = %d outside allowed domain.\n",
//                gc_rkStage, p->id);
    err = gc_cERR_DOMAIN_OVERSTEP;
    return err;
  }

  Particles_GetWeights(p, p->cell, w, grid); 
  i = p->cell[IDIR];
  j = p->cell[JDIR];
  k = p->cell[KDIR];
   
  if ( DIM_EXPAND(   data->flag[k][j][i]   & (FLAG_GCA_FAILURE)
                  || data->flag[k][j][i-1] & (FLAG_GCA_FAILURE)
                  || data->flag[k][j][i+1] & (FLAG_GCA_FAILURE),
                  || data->flag[k][j+1][i] & (FLAG_GCA_FAILURE)
                  || data->flag[k][j-1][i] & (FLAG_GCA_FAILURE),
                  || data->flag[k+1][j][i] & (FLAG_GCA_FAILURE)
                  || data->flag[k-1][j][i] & (FLAG_GCA_FAILURE))) {
    return  gc_cERR_INVALID;
  }

/* --------------------------------------------------------
   2. Interpolate all fields
   -------------------------------------------------------- */
  
  nfields = 36;
  /* -- WARNING: initializing these two arrays before GetWeights
        leads to a error -- */
  double ***Fields[nfields];
  double Interpolated_fields[nfields];

  Fields[0] = data->Vc[BX1];
  Fields[1] = data->Vc[BX2];
  Fields[2] = data->Vc[BX3];
  
  for(dir = 0; dir < 3; dir++){
    Fields[dir +  3] = gc_cE[dir];
    Fields[dir +  6] = gc_bGrad_b[dir];
    Fields[dir +  9] = gc_vEGrad_b[dir];
    Fields[dir + 12] = gc_bGrad_vE[dir];
    Fields[dir + 15] = gc_vEGrad_vE[dir];
    Fields[dir + 18] = gc_Grad_Bf1[dir];
    Fields[dir + 21] = gc_Grad_B[dir];
    Fields[dir + 24] = gc_b_old[dir];
    Fields[dir + 27] = gc_vE_old[dir];
    Fields[dir + 30] = gc_vE[dir];
    Fields[dir + 33] = gc_cEres[dir];
  }
  Particles_InterpolateArr(Fields, nfields, w, p->cell, Interpolated_fields);
  
  for (dir = 0; dir < 3; dir++){
    B[dir]      = Interpolated_fields[dir];
    cE[dir]     = Interpolated_fields[dir +  3];
    bdb[dir]    = Interpolated_fields[dir +  6];
    vEdb[dir]   = Interpolated_fields[dir +  9];
    bdvE[dir]   = Interpolated_fields[dir + 12];
    vEdvE[dir]  = Interpolated_fields[dir + 15];
    dBf1[dir]   = Interpolated_fields[dir + 18];
    dB[dir]     = Interpolated_fields[dir + 21];
    b_old[dir]  = Interpolated_fields[dir + 24];
    vE_old[dir] = Interpolated_fields[dir + 27];
    vE[dir]     = Interpolated_fields[dir + 30];
    cEres[dir]   = Interpolated_fields[dir + 33];
  }
  Bmagf1_old    = Particles_Interpolate(gc_Bmf1_old, w, p->cell);

/* --------------------------------------------------------
   3a. Compute Electric field separately so that ideal
       term will be orthogonal to B.
   -------------------------------------------------------- */
/*
  double vg[3];
  vg[IDIR] = Particles_Interpolate(data->Vc[VX1], w, p->cell);
  vg[JDIR] = Particles_Interpolate(data->Vc[VX2], w, p->cell);
  vg[KDIR] = Particles_Interpolate(data->Vc[VX3], w, p->cell);
  cE[IDIR] = -(vg[JDIR]*B[KDIR] - vg[KDIR]*B[JDIR]);
  cE[JDIR] = -(vg[KDIR]*B[IDIR] - vg[IDIR]*B[KDIR]);
  cE[KDIR] = -(vg[IDIR]*B[JDIR] - vg[JDIR]*B[IDIR]);
  #if PHYSICS == ResRMHD
  cE[IDIR] += Eres[IDIR];
  cE[JDIR] += Eres[JDIR];
  cE[KDIR] += Eres[KDIR];
  #endif
*/
// Use Analytic
//double vfluid[256];
//Init (vfluid, p->coord[IDIR], p->coord[JDIR], p->coord[KDIR]);
//B[IDIR] = vfluid[BX1]; B[JDIR] = vfluid[BX2]; B[KDIR] = vfluid[BX3];
//cE[IDIR] = vfluid[EX1]; cE[JDIR] = vfluid[EX2]; cE[KDIR] = vfluid[EX3];      

/* --------------------------------------------------------
   3a. Compute important quantities in the lab frame
   -------------------------------------------------------- */
   
  Bmag2  = DOT_PRODUCT(B,B);
  cEmag2 = DOT_PRODUCT(cE,cE);
  dBmag  = sqrt(DOT_PRODUCT(dB,dB));

  Bmag     = sqrt(Bmag2);
  Bmag_inv = 1.0/(Bmag + EPS_B);

  b[IDIR] = B[IDIR]*Bmag_inv;
  b[JDIR] = B[JDIR]*Bmag_inv;
  b[KDIR] = B[KDIR]*Bmag_inv;
  cEpar = DOT_PRODUCT(cE,b);

  vE[IDIR] = CROSS_X1(cE,b)*Bmag_inv;  /* vE is orthogonal to both E and B */
  vE[JDIR] = CROSS_X2(cE,b)*Bmag_inv;  /* by construction                  */
  vE[KDIR] = CROSS_X3(cE,b)*Bmag_inv;
  vEmag2   = DOT_PRODUCT(vE,vE); 

/* --------------------------------------------------------
   3b. Cleaning step to ensure that:
       - vE  is perpendicular to b (done already)
       - bdb is perpendicular to b
   -------------------------------------------------------- */
/*
  scrh = DOT_PRODUCT(vE,b);  
  vE[IDIR] = vE[IDIR] - scrh*b[IDIR];
  vE[JDIR] = vE[JDIR] - scrh*b[JDIR];
  vE[KDIR] = vE[KDIR] - scrh*b[KDIR];
  vEmag2   = DOT_PRODUCT(vE,vE); 
*/

  scrh = DOT_PRODUCT(bdb,b);  
  bdb[IDIR] -= scrh*b[IDIR];
  bdb[JDIR] -= scrh*b[JDIR];
  bdb[KDIR] -= scrh*b[KDIR];

  gammaE  = 1./sqrt(1. - vEmag2*inv_c2);
  gammaE2 = gammaE*gammaE;
  Bmagf1 = Bmag/gammaE;

  /*omega = e_mc*sqrt(0.5*(Bmag2 - Emag2) //complete omega expression
          + 0.5*sqrt( (Bmag2 - Emag2)*(Bmag2 - Emag2) 
                      + 4.*(DOT_PRODUCT(E,B)*DOT_PRODUCT(E,B))));*/
/*  omega  = e_mc*Bmag*sqrt(1.0 - vEmag2);  */
  omega     = e_mc*Bmag;
  omega_inv = 1./omega;
  gamma     = sqrt((1.0 + upar*upar*inv_c2 + 2.0*mc_e*mu*omega*inv_c2)/(1.0 - vEmag2*inv_c2));
  gamma_inv = 1.0/gamma;
  vpar      = upar*gamma_inv;

/* --------------------------------------------------------
   4. Multiple checks for GC conditions
   -------------------------------------------------------- */

  /* -- A.  Compute gyroradius R_L and check if it is larger than
            cell dimension -- */
            
  R_L = mc_e*sqrt(2.*mu*Bmag)*Bmag_inv;
  for (dir = 0; dir < DIMENSIONS; dir++) {
    if ( R_L > grid->dx[dir][p->cell[dir]]){
      err = gc_cERR_LARMOR_OVERSTEP;
      #if PARTICLES_CR_GC_DEBUG == YES
      //PRINT_VAR(dir);
      //PRINT_VAR(R_L);
      //PRINT_VAR(grid->dx[dir][p->cell[dir]]);
      #endif
    }
  }
  
  /* -- B.  Check if \epsilon = |u/(\omega L)| << 1. \epsilon
            is the Taylor series parameter for the GCA, u is
            the 4 velocity, omega the gyration frequency and
            L is the lenght scale for which 
            \Delta B is comparable to B-- */
  
  if (   fabs(p->speed[IDIR])*dBmag*omega_inv        > 0.1*Bmag
      || fabs(gammaE*sqrt(vEmag2))*dBmag*omega_inv   > 0.1*Bmag){
    err = gc_cERR_FAST_VARYING;
  }
    
  /* -- Time -- */ /*time is not necessary
  for (dir = 0; dir < 3; dir++){
    BtoB_old[dir] = Bmag*(b[dir] - b_old[dir]);
  }
  if (gamma*sqrt(DOT_PRODUCT(BtoB_old,BtoB_old))*inv_dt/omega > 0.1*Bmag){
    *err = 0;
  }*/
  
  /* -- C.  Check if E_perp > B, would lead to analytical
            singulartity for the GCA -- */
  if( ((cEmag2 - cEpar*cEpar)*inv_c2) > (Bmag2) || (vEmag2*inv_c2) > 1.) {
    int i1,j1,k1;
    printLog ("! Particles_CR_getGC(): |Eperp| > |B|\n");
    #if PHYSICS == ResRMHD
    printLog ("Emag(i,j,k)    = %12.6e  [%12.6e]\n", sqrt(cEmag2),
            sqrt(  data->Vc[EX1][k][j][i]*data->Vc[EX1][k][j][i]
                 + data->Vc[EX2][k][j][i]*data->Vc[EX2][k][j][i]
                 + data->Vc[EX3][k][j][i]*data->Vc[EX3][k][j][i]));
    #endif
    printLog ("Bmag(i,j,k)    = %12.6e  [%12.6e]\n", sqrt(Bmag2),
            sqrt( data->Vc[BX1][k][j][i]*data->Vc[BX1][k][j][i]
                + data->Vc[BX2][k][j][i]*data->Vc[BX2][k][j][i]
                + data->Vc[BX3][k][j][i]*data->Vc[BX3][k][j][i]));
    printLog ("B(p) = "); ShowVector(B,3);
    printLog ("E(p) = "); ShowVector(cE,3);
    printLog ("vEmag2 = %12.6e\n",vEmag2);
/*
    for (k1 = -INCLUDE_KDIR; k1 <= INCLUDE_KDIR; k1++){
    for (j1 = -INCLUDE_JDIR; j1 <= INCLUDE_JDIR; j1++){
    for (i1 = -INCLUDE_IDIR; i1 <= INCLUDE_IDIR; i1++){
      printLog ("    B(%d,%d,%d)  = %12.6e  %12.6e  %12.6e\n",
                     i1,j1,k1,
                     data->Vc[BX1][k+k1][j+j1][i+i1],
                     data->Vc[BX2][k+k1][j+j1][i+i1],
                     data->Vc[BX3][k+k1][j+j1][i+i1]);
    }}}
*/

    return gc_cERR_INVALID;
  }
  /* -- D.  Check if particle is located in last ghost zone,
            bdb would be NaN -- */
            
  if ( isnan (bdb[IDIR]) || isnan (bdb[JDIR]) || isnan (bdb[KDIR]) ){
    int i1, j1, k1;
    err = gc_cERR_INVALID;
    printLog ("! Particles_CR_getGC() [stage = %d]:\n", gc_rkStage);
    printLog ("    bdB = "); ShowVector(bdb,3);
    printLog ("    particle is in last ghost zone ?\n");
    printLog ("    i = %d  (x-dom size = %d, %d)\n",i,0, NX1_TOT-1);
    printLog ("    j = %d  (y-dom size = %d, %d)\n",j,0, NX2_TOT-1);
    printLog ("    k = %d  (k-dom size = %d, %d)\n",k,0, NX3_TOT-1);

    for (k1 = -INCLUDE_KDIR; k1 <= INCLUDE_KDIR; k1++){
    for (j1 = -INCLUDE_JDIR; j1 <= INCLUDE_JDIR; j1++){
    for (i1 = -INCLUDE_IDIR; i1 <= INCLUDE_IDIR; i1++){
      printLog ("    bgradb(%d,%d,%d)  = %12.6e  %12.6e  %12.6e\n",
                     i1,j1,k1,
                     gc_bGrad_b[IDIR][k+k1][j+j1][i+i1],
                     gc_bGrad_b[JDIR][k+k1][j+j1][i+i1],
                     gc_bGrad_b[KDIR][k+k1][j+j1][i+i1]);
    }}}
    for (k1 = -INCLUDE_KDIR; k1 <= INCLUDE_KDIR; k1++){
    for (j1 = -INCLUDE_JDIR; j1 <= INCLUDE_JDIR; j1++){
    for (i1 = -INCLUDE_IDIR; i1 <= INCLUDE_IDIR; i1++){
      printLog ("    B(%d,%d,%d)  = %12.6e  %12.6e  %12.6e\n",
                     i1,j1,k1,
                     data->Vc[BX1][k+k1][j+j1][i+i1],
                     data->Vc[BX2][k+k1][j+j1][i+i1],
                     data->Vc[BX3][k+k1][j+j1][i+i1]);
    }}}

    for (k1 = -INCLUDE_KDIR; k1 <= INCLUDE_KDIR; k1++){
    for (j1 = -INCLUDE_JDIR; j1 <= INCLUDE_JDIR; j1++){
    for (i1 = -INCLUDE_IDIR; i1 <= INCLUDE_IDIR; i1++){
      printLog ("w[%d][%d][%d]   = %12.6e\n",i1,j1,k1, w[k1][j1][i1]);
    }}}
    printLog ("    Lower PARTICLES_CR_NCELL_MAX is required. \n");
  }
  
/* --------------------------------------------------------
   5. Compute dRdt
   -------------------------------------------------------- */
  
  sum[IDIR] =   vpar*upar*bdb[IDIR]
              + upar*(vEdb[IDIR] + bdvE[IDIR])
              + gamma*vEdvE[IDIR]
              + e_mc*inv_c2*vpar*cEpar*vE[IDIR]
              + mu*gamma_inv*dBf1[IDIR];
  sum[JDIR] =   vpar*upar*bdb[JDIR]
              + upar*(vEdb[JDIR] + bdvE[JDIR])
              + gamma*vEdvE[JDIR]
              + e_mc*inv_c2*vpar*cEpar*vE[JDIR]
              + mu*gamma_inv*dBf1[JDIR];
  sum[KDIR] =   vpar*upar*bdb[KDIR]
              + upar*(vEdb[KDIR] + bdvE[KDIR])
              + gamma*vEdvE[KDIR]
              + e_mc*inv_c2*vpar*cEpar*vE[KDIR]
              + mu*gamma_inv*dBf1[KDIR];
  
  /* -- Compute and add time terms, should not be added
        in first call because of missing (n-1) term -- */
#if PARTICLES_CR_GC_TIME_DER == YES
  if(inv_dt > 0.){
    timeterms[IDIR] = inv_dt*(
                        upar*(b[IDIR] - b_old[IDIR])
                      + gamma*(vE[IDIR] - vE_old[IDIR]) 
                      + mu*gamma_inv*inv_c*vE[IDIR]*(Bmagf1 - Bmagf1_old));
    timeterms[JDIR] = inv_dt*(
                        upar*(b[JDIR] - b_old[JDIR])
                      + gamma*(vE[JDIR] - vE_old[JDIR]) 
                      + mu*gamma_inv*inv_c*vE[JDIR]*(Bmagf1 - Bmagf1_old));
    timeterms[KDIR] = inv_dt*(
                        upar*(b[KDIR] - b_old[KDIR])
                      + gamma*(vE[KDIR] - vE_old[KDIR]) 
                      + mu*gamma_inv*inv_c*vE[KDIR]*(Bmagf1 - Bmagf1_old));
  
    sum[IDIR] += timeterms[IDIR];
    sum[JDIR] += timeterms[JDIR];
    sum[KDIR] += timeterms[KDIR];
  }
#endif

  scrh = mc_e*gammaE2*Bmag_inv;
  dRdt[IDIR] = vpar*b[IDIR] + vE[IDIR] + scrh*CROSS_X1(b,sum);
  dRdt[JDIR] = vpar*b[JDIR] + vE[JDIR] + scrh*CROSS_X2(b,sum);
  dRdt[KDIR] = vpar*b[KDIR] + vE[KDIR] + scrh*CROSS_X3(b,sum);              

  dRdt[3]    = - upar*DOT_PRODUCT(b,bdvE)
               - gamma*DOT_PRODUCT(b,vEdvE) 
               + e_mc*cEpar
               - mu*gamma_inv*DOT_PRODUCT(b,dBf1);

  /* -- Compute and add time terms, should not be added
        in first call because of missing (n-1) term -- */
#if PARTICLES_CR_GC_TIME_DER == YES
  if(inv_dt > 0.){
    dvEdt[IDIR] = inv_dt*(vE[IDIR] - vE_old[IDIR]);
    dvEdt[JDIR] = inv_dt*(vE[JDIR] - vE_old[JDIR]);
    dvEdt[KDIR] = inv_dt*(vE[KDIR] - vE_old[KDIR]);
    dRdt[3] -= gamma*DOT_PRODUCT(b,dvEdt);
  }
#endif

  #if PARTICLES_CR_GC_DEBUG == YES
  if(b[JDIR] - b_old[JDIR]!=0 && inv_dt != -1. && -1==1){
    print("\n");
    PRINT_VAR(gammaE2);
    PRINT_VAR(Bmag_inv);
    PRINT_VAR(inv_dt);
    PRINT_VAR(mu);
    PRINT_VAR(gamma_inv);
    PRINT_VAR(DOT_PRODUCT(b,dBf1));
    PRINT_VAR(b[IDIR]);
    PRINT_VAR(b[JDIR]);
    PRINT_VAR(b[KDIR]);
    printLog ("  timeterms\t= ");ShowVector(timeterms,3);
    PRINT_VAR(CROSS_X1(b,sum));
    PRINT_VAR(CROSS_X2(b,sum));
    PRINT_VAR(CROSS_X3(b,sum));
    printLog ("  dRdt\t= ");ShowVector(dRdt,4);
    printLog ("  dvEdt\t= ");ShowVector(dvEdt,3);
    PRINT_VAR(DOT_PRODUCT(b,dvEdt));
    PRINT_VAR(b[IDIR] - b_old[IDIR]);
    PRINT_VAR(b[JDIR] - b_old[JDIR]);
    PRINT_VAR(b[KDIR] - b_old[KDIR]);
    printLog ("  b_old\t= ");ShowVector(b_old,3);
    printLog ("  vE\t= ");ShowVector(vE,3);
    PRINT_VAR(vpar*b[IDIR]);
    PRINT_VAR(vpar*b[JDIR]);
    PRINT_VAR(vpar*b[KDIR]);
  }
  #endif

  dRdtmag2 = DOT_PRODUCT(dRdt, dRdt);
  if(sqrt(dRdtmag2) > 2.0){   // Check if dR/dt exceeds speed of light. 
//    printLog("! Particles_CR_getGC(): |dR/dt| = %12.6e > c\n", sqrt(dRdtmag2));
    return gc_cERR_dRdt;
  }
  
  #if PARTICLES_CR_GC_DEBUG == YES
  /* -- Use gamma as a debug tool -- */
  if(gamma < 1. || isnan(gamma)){
    printLog("\n");
    PRINT_VAR(p->id);
    PRINT_VAR(gamma);
    PRINT_VAR(mu);
    PRINT_VAR(omega);
    PRINT_VAR(1 - 2*mu*omega);
    PRINT_VAR(dRdtmag2);
    printLog("\n");
  }
  #endif
  
  /*******************--DEBUG--*******************/
  if (isnan(dRdtmag2)) {
    printLog ("! Particles_CR_getGC(): nan found in dRdtmag\n");
    #if PARTICLES_CR_GC_DEBUG == YES
    printLog ("  p(id)\t= %d; pcell = %d %d %d\n",p->id, p->cell[IDIR], p->cell[JDIR], p->cell[KDIR]);
    PRINT_VAR(gamma);
    PRINT_VAR(gammaE);
    PRINT_VAR(mu);
    PRINT_VAR(omega);
    PRINT_VAR((1. - 2*mu*omega)/(1. - dRdtmag2));
    PRINT_VAR(vEmag2);
    PRINT_VAR(dRdtmag2);
    PRINT_VAR(cEpar);
    PRINT_VAR(Bmagf1);
    PRINT_VAR(Bmagf1_old);
    printLog ("  p->speed\t= ");ShowVector(p->speed,3);
    printLog ("  dRdt\t= ");ShowVector(dRdt,5);
    printLog ("  b\t= "); ShowVector(b,3);
    printLog ("  b_old\t= "); ShowVector(b_old,3);
    printLog ("  bdb\t= "); ShowVector(bdb,3);
    printLog ("  vEdb\t= "); ShowVector(vEdb,3);
    printLog ("  vE\t= "); ShowVector(vE,3);
    printLog ("  vE_old\t= "); ShowVector(vE_old,3);
    printLog ("  E\t= "); ShowVector(cE,3);
    printLog ("  B\t= "); ShowVector(B,3);
    printLog ("  bdvE \t= "); ShowVector(bdvE,3);
    printLog ("  vEdvE\t= "); ShowVector(vEdvE,3);
    printLog ("  dBf1\t= "); ShowVector(dBf1,3);
    #endif
    err = gc_cERR_INVALID;
    //QUIT_PLUTO(1);
  }
  /*******************--DEBUG--*******************/

  return err;
}


/* ********************************************************************* */
void Particles_CR_GCinfo(Particle *p, double *pcoord0)
/*
 *
 *********************************************************************** */
{
  printLog ("! Particles_CR_GCInfo(): [rk stage = %d]\n",gc_rkStage);
  printLog ("!   pcoord0  = "); ShowVector(pcoord0,3);
  printLog ("!   p->coord = "); ShowVector(p->coord,3);
}

#endif /* PARTICLES_CR_GC == YES */
