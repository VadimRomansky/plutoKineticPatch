/* ///////////////////////////////////////////////////////////////////// */
/*! 
  \file  
  \brief  Solves the linear transport step of the FARGO-MHD algorithm. 
 
  Implements a uniform shift of the conservative fluid variables in 
  the direction of orbital motion which is assumed to be 
  \f$ x_2=y,\phi\f$ (for Cartesian or polar coordinates) or 
  \f$ x_3 = \phi\f$  (for spherical coordinates).
  The shift is converted into an integer part "m" and a fractional
  remainder "eps".\n
 
  Parallelization is performed by adding extra data layers
  above and below the current data set so as to enlarge the required
  computational stencil in the orbital direction.
  The additional buffers are received from the processors lying,
  respectively, below and above and their size changes dynamically 
  from one step to the next.
 
  \b Reference
    - "A conservative orbital advection scheme for simulations 
       of magnetized shear flows with the PLUTO code"\n
       Mignone et al, A&A (2012) 545, A152  (--> Section 2.4)

  \author A. Mignone (andrea.mignone@unito.it)\n
          G. Muscianisi (g.muscianisi@cineca.it)
  \date   July 10, 2019
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

#define FARGO_MOD(s) ((s) < SBEG ? (s)+NS:((s) > SEND ? (s)-NS:(s)))
#define MAX_BUF_SIZE 64

static void FARGO_Flux (double *, double *, double);

static void FARGO_ParallelExchange (Data_Arr U, Data_Arr Us,
                             double ***rbl[], RBox *,
                             double ***rbu[], RBox *, Grid *grid);

/* ********************************************************************* */
void FARGO_ShiftSolution(Data_Arr U, Data_Arr Us, Grid *grid)
/*!
 * Shifts conserved variables in the orbital direction.
 *
 * \param [in,out]  U    a 3D array of conserved, zone-centered values
 * \param [in,out]  Us   a 3D array of staggered magnetic fields
 * \param [in]      grid  pointer to Grid structure;
 *
 * \return  This function has no return value.
 * \todo    Optimization: avoid using  too many if statements like
 *          on nproc_s > 1 or == 1
 *********************************************************************** */
#if GEOMETRY != SPHERICAL
 #define s j   /* -- orbital direction index -- */
#else
 #define s k
#endif
{
  int    i, j, k, nv, ngh, nvar_c;
  int    sm, m;
  static int mmax, mmin;  /* Used to determine buffer size from time to time */
  int    nbuf_lower, nbuf_upper, nproc_s = 1;
  double dphi, w, Lphi, dL, eps;
  double **wA, ***Uc[NVAR];
  double *q, *flux;
  double *x1  = grid->x[IDIR],   *x2 = grid->x[JDIR],  *x3 = grid->x[KDIR];
  double *x1p = grid->xr[IDIR], *x2p = grid->xr[JDIR], *x3p = grid->xr[KDIR];
  double *x1m = grid->xl[IDIR], *x2m = grid->xl[JDIR], *x3m = grid->xl[KDIR];
  double *dx1 = grid->dx[IDIR], *dx2 = grid->dx[JDIR], *dx3 = grid->dx[KDIR];

  static double *q00, *flux00;

  static double ***Ez, ***Ex;
  #ifdef PARALLEL
   double ***upper_buf[NVAR];  /* upper layer buffer  */
   double ***lower_buf[NVAR];  /* lower layer buffer  */
   RBox    lbox, ubox;              /* lower and upper box structures  */
  #endif

#if GEOMETRY == CYLINDRICAL
  printLog ("! FARGO cannot be used in this geometry\n");
  QUIT_PLUTO(1);
#endif

/* --------------------------------------------------------
   1. Add source term (if any)
   -------------------------------------------------------- */

#ifdef SHEARINGBOX
  FARGO_Source(U, Us, g_dt, grid);  
#endif

#if (GEOMETRY != SPHERICAL) && !INCLUDE_JDIR
  return;   /* No shift needed in axisymmetric geometry */
#endif

/* --------------------------------------------------------
   2. Allocate memory for static arrays.
      In parallel mode, one-dimensional arrays must be
      large enough to contain data values coming from
      neighbour processors.
      For this reason we augment them with an extra zone
      buffer containing MAX_BUF_SIZE cells on both sides.
      In this way, the array indices for q can be negative.
   
    
      <..MAX_BUF_SIZE..>|---|---| ... |---|----|<..MAX_BUF_SIZE..>
                          0   1            NX2-1
   -------------------------------------------------------- */

  if (q00 == NULL){
    #if GEOMETRY == SPHERICAL && DIMENSIONS == 3
    q00    = ARRAY_1D(NX3_TOT + 2*MAX_BUF_SIZE,double);
    flux00 = ARRAY_1D(NX3_TOT + 2*MAX_BUF_SIZE,double);
    #else
    q00    = ARRAY_1D(NX2_TOT + 2*MAX_BUF_SIZE,double);
    flux00 = ARRAY_1D(NX2_TOT + 2*MAX_BUF_SIZE,double);
    #endif
    #ifdef STAGGERED_MHD
    Ez = ARRAY_3D(NX3_TOT, NX2_TOT, NX1_TOT, double);
    Ex = ARRAY_3D(NX3_TOT, NX2_TOT, NX1_TOT, double);
    #endif
 
  /* ------------------------------------------------------
      The value of mmax and mmin is kept from call to call
      because they are used to determine the size of the
      send/receive buffers each time step.
      At the very beginning, we initialize them to the
      largest possible value. 
     ------------------------------------------------------ */

    mmax = MIN(NS, MAX_BUF_SIZE-1) - grid->nghost[SDIR] - 2;
    mmin = mmax;
  }

  q    = q00    + MAX_BUF_SIZE;
  flux = flux00 + MAX_BUF_SIZE;

/* --------------------------------------------------------
   3. Get the pointer to the orbital speed array
   -------------------------------------------------------- */

  wA = FARGO_Velocity();

  nproc_s = grid->nproc[SDIR];   /* -- # of processors along shift dir -- */
  ngh     = grid->nghost[SDIR]; 
  Lphi    = g_domEnd[SDIR] - g_domBeg[SDIR];

  dphi = grid->dx[SDIR][SBEG];

/* --------------------------------------------------------
   4. Fargo parallelization is performed by adding extra
      data layers above and below the current data set so
      as to enlarge the required computational stencil in
      the orbital direction:

           ^
           |    |recv_lower|
           |    +----------+
           |    |          |
           |    |     U    |
           |    |          |
           |    +----------+
           |    |recv_upper|

      The buffers recv_upper and recv_lower are received
      from the processors lying, respectively, below and
      above and their size changes dynamically from
      one step to the next.
   -------------------------------------------------------- */
    
#ifdef PARALLEL
  if (nproc_s > 1){ 
    nbuf_lower = abs(mmin) + ngh + 1;
    nbuf_upper =     mmax  + ngh + 1;

  /* -- check that buffer is not larger than processor size -- */

    if (nbuf_upper > NS || nbuf_lower > NS){
      printLog ("! FARGO_ShiftSolution: buffer size exceeds processsor size\n");
      QUIT_PLUTO(1);
    }

  /* -- check that buffer is less than MAX_BUF_SIZE -- */
 
    if (nbuf_upper > MAX_BUF_SIZE || nbuf_lower > MAX_BUF_SIZE){
      printLog ("! FARGO_ShiftSolution: buffer size exceeds maximum length\n");
      QUIT_PLUTO(1);
    }

  /* -- set grid indices of the lower & upper boxes -- */

    #if GEOMETRY == CARTESIAN || GEOMETRY == POLAR

    lbox.ibeg = IBEG; lbox.jbeg =            0; lbox.kbeg = KBEG; 
    lbox.iend = IEND; lbox.jend = nbuf_lower-1; lbox.kend = KEND;
    lbox.di   = lbox.dj = lbox.dk = 1;

    ubox.ibeg = IBEG; ubox.jbeg =            0; ubox.kbeg = KBEG; 
    ubox.iend = IEND; ubox.jend = nbuf_upper-1; ubox.kend = KEND;
    ubox.di = ubox.dj = ubox.dk = 1;

    #ifdef STAGGERED_MHD  /* -- enlarge boxes to fit staggered fields --*/
    DIM_EXPAND(lbox.ibeg--; ubox.ibeg--;  ,
                                     ;  ,
             lbox.kbeg--; ubox.kbeg--;)
    #endif 

    #elif GEOMETRY == SPHERICAL

    lbox.ibeg = IBEG; lbox.jbeg = JBEG; lbox.kbeg =            0; 
    lbox.iend = IEND; lbox.jend = JEND; lbox.kend = nbuf_lower-1; 
    lbox.di = lbox.dj = lbox.dk = 1;

    ubox.ibeg = IBEG; ubox.jbeg = JBEG; ubox.kbeg =            0; 
    ubox.iend = IEND; ubox.jend = JEND; ubox.kend = nbuf_upper-1; 
    ubox.di = ubox.dj = ubox.dk = 1;

    #ifdef STAGGERED_MHD  /* -- enlarge boxes to fit staggered fields --*/
    DIM_EXPAND(lbox.ibeg--; ubox.ibeg--;  ,
             lbox.jbeg--; ubox.jbeg--;  ,
                                     ; )
    #endif 

    #endif  /* GEOMETRY */

  /* -- allocate memory -- */

    for (nv = 0; nv < NVAR; nv++){
      upper_buf[nv] = ARRAY_BOX(ubox.kbeg, ubox.kend,
                                ubox.jbeg, ubox.jend,
                                ubox.ibeg, ubox.iend, double);
      lower_buf[nv] = ARRAY_BOX(lbox.kbeg, lbox.kend,
                                lbox.jbeg, lbox.jend,
                                lbox.ibeg, lbox.iend, double);
    }

  /* -- exchange data values between neighbour processes -- */
 
     FARGO_ParallelExchange(U, Us, lower_buf, &lbox, upper_buf, &ubox, grid);
   }
#endif  /* PARALLEL */

/* --------------------------------------------------------
   5. Shift cell-centered quantities
   -------------------------------------------------------- */

  mmax = mmin = 0;
#if GEOMETRY == CARTESIAN || GEOMETRY == POLAR
  KDOM_LOOP(k) IDOM_LOOP(i) {
#else
  JDOM_LOOP(j) IDOM_LOOP(i) {
#endif

    #if GEOMETRY == CARTESIAN
    w = wA[k][i];
    #elif GEOMETRY == POLAR
    w = wA[k][i]/x1[i];
    #elif GEOMETRY == SPHERICAL
    w = wA[j][i]/(x1[i]*sin(x2[j]));
    #endif

  /* --------------------------------------------
     5a. Obtain shift in terms of an integer
         number of cells (m) plus a remainder eps. 
         The integer part is obtained by rounding
         dL/dphi  to the nearest integer.
         Examples:

         dL/dphi = 3.4   -->   m =  3, eps =  0.4
         dL/dphi = 4.7   -->   m =  5, eps = -0.3
         dL/dphi = -0.8  -->   m = -1, eps =  0.2
     -------------------------------------------- */

    dL  = w*g_dt;                 /* spatial shift */
    dL  = fmod(dL, Lphi);
    m   = (int)floor(dL/dphi + 0.5); /* nearest integer part */
    eps = dL/dphi - m;               /* fractional remainder */

  /* -- check if shift is compatible with estimated buffer size -- */

    if (nproc_s > 1){
      if (w > 0.0 && (m + ngh) > nbuf_upper){
        printLog ("! FARGO_ShiftSolution: positive shift too large\n");
        QUIT_PLUTO(1);
      }
      if (w < 0.0 && (-m + ngh) > nbuf_lower){
        printLog ("! FARGO_ShiftSolution: negative shift too large\n");
        QUIT_PLUTO(1);
      }
    }

    mmax = MAX(m, mmax);
    mmin = MIN(m, mmin);

  /* --------------------------------------------
     5b. Start main loop on variables 
     -------------------------------------------- */
   
    for (nv = NVAR; nv--;   ){

      #ifdef STAGGERED_MHD
      DIM_EXPAND(if (nv == BX1) continue;  ,
               if (nv == BX2) continue;  ,
               if (nv == BX3) continue;)
      #endif

  /* -- copy 3D data into 1D array -- */

      SDOM_LOOP(s) q[s] = U[k][j][i][nv];

      if (nproc_s > 1){ /* -- copy values from lower and upper buffers -- */
        #ifdef PARALLEL
        for (s = SBEG-1; s >= SBEG - nbuf_upper; s--){
          q[s] = FARGO_ARRAY_INDEX(upper_buf[nv], SBEG-1-s, k, j, i);
        }
        for (s = SEND+1; s <= SEND + nbuf_lower; s++){
          q[s] = FARGO_ARRAY_INDEX(lower_buf[nv], s-SEND-1, k, j, i);
        } 
        FARGO_Flux (q-m, flux-m, eps); 
        #endif
      }else{               /* -- impose periodic b.c. -- */
        for (s = 1; s <= ngh; s++) {
          q[SBEG - s] = q[SBEG - s + NS]; 
          q[SEND + s] = q[SEND + s - NS]; 
        }    
        FARGO_Flux (q, flux, eps); 
      }

  /* -- update main solution array -- */

      if (nproc_s > 1){
        #ifdef PARALLEL
        for (s = SBEG; s <= SEND; s++){
          sm = s - m;
          U[k][j][i][nv] = q[sm] - eps*(flux[sm] - flux[sm-1]);
        }
        #endif
      }else{
        for (s = SBEG; s <= SEND; s++){
          sm = FARGO_MOD(s - m);
          U[k][j][i][nv] = q[sm] - eps*(flux[sm] - flux[sm-1]);
        }
      }
    }
  }  /* -- end main spatial loop -- */

/* --------------------------------------------------------
   6. Shift staggered magnetic fields
   -------------------------------------------------------- */

#ifdef STAGGERED_MHD 
{
  int    ss;
  double A1p, A1m, A2p, A2m;
  double ***bx1, ***bx2, ***bx3;

/* -- pointer shortcuts -- */

  DIM_EXPAND(bx1 = Us[BX1s];  ,
             bx2 = Us[BX2s];  ,
             bx3 = Us[BX3s];)

/* ----------------------------------------------
   6a. Compute:
        Ez   at x(i+1/2), y(j+1/2), z(k)  (Cartesian or Polar) 
       -Eth  at x(i+1/2), y(j), z(k+1/2)  (Spherical) 
   ---------------------------------------------- */

  #if GEOMETRY == CARTESIAN || GEOMETRY == POLAR
  KDOM_LOOP(k) for (i = IBEG-1; i <= IEND; i++){
  #else
  JDOM_LOOP(j) for (i = IBEG-1; i <= IEND; i++){
  #endif

    #if GEOMETRY == CARTESIAN
    w = 0.5*(wA[k][i] + wA[k][i+1]);
    #elif GEOMETRY == POLAR
    w = 0.5*(wA[k][i] + wA[k][i+1])/x1p[i];
    #elif GEOMETRY == SPHERICAL
    w = 0.5*(wA[j][i] + wA[j][i+1])/(x1p[i]*sin(x2[j]));
    #endif

  /* --------------------------------------------
     6b. Obtain shift in terms of an integer
         number of cells (m) plus a remainder eps. 
     -------------------------------------------- */

    dL  = w*g_dt;                 /* spatial shift */
    dL  = fmod(dL, Lphi);
    m   = (int)floor(dL/dphi + 0.5); /* nearest integer part */
    eps = dL/dphi - m;               /* remainder in units of dphi */

  /* -- copy 3D array into 1D buffer -- */

    SDOM_LOOP(s) q[s] = bx1[k][j][i];
    if (nproc_s > 1){  /* -- copy values from lower and upper buffers -- */
      #ifdef PARALLEL
      for (s = SBEG-1; s >= SBEG - nbuf_upper; s--){
        q[s] = FARGO_ARRAY_INDEX(upper_buf[BX1], SBEG-1-s, k, j, i); 
      }
      for (s = SEND+1; s <= SEND + nbuf_lower; s++){
        q[s] = FARGO_ARRAY_INDEX(lower_buf[BX1], s-SEND-1, k, j, i);
      }
      FARGO_Flux (q-m, flux-m, eps);
      #endif
    }else{                  /* -- assign periodic b.c. -- */
      for (s = 1; s <= ngh; s++) {
        q[SBEG - s] = q[SBEG - s + NS]; 
        q[SEND + s] = q[SEND + s - NS];  
      }    
      FARGO_Flux (q, flux, eps);
    }

  /* -- update main solution array -- */

    if (nproc_s > 1){ 
      #ifdef PARALLEL
      if (m >= 0){
        for (s = SBEG - 1; s <= SEND; s++) {
          sm = s - m;
          Ez[k][j][i] = eps*flux[sm];
          for (ss = s; ss > (s-m); ss--) Ez[k][j][i] += q[ss];
          Ez[k][j][i] *= dphi;
        }
      }else{
        for (s = SBEG-1; s <= SEND; s++) {
          sm = s - m;
          Ez[k][j][i] = eps*flux[sm];
          for (ss = s+1; ss < (s-m+1); ss++) Ez[k][j][i] -= q[ss];
          Ez[k][j][i] *= dphi;
        }
      } 
      #endif
    }else{
      if (m >= 0){
        for (s = SBEG - 1; s <= SEND; s++) {
          sm = FARGO_MOD(s - m);
          Ez[k][j][i] = eps*flux[sm];
          for (ss = s; ss > (s-m); ss--) Ez[k][j][i] += q[FARGO_MOD(ss)];
          Ez[k][j][i] *= dphi;
        }
      }else {
        for (s = SBEG-1; s <= SEND; s++) {
          sm = FARGO_MOD(s - m);
          Ez[k][j][i] = eps*flux[sm];
          for (ss = s+1; ss < (s-m+1); ss++) Ez[k][j][i] -= q[FARGO_MOD(ss)];
          Ez[k][j][i] *= dphi;
        }
      }
    }
  }

#if INCLUDE_KDIR

/* ----------------------------------------------
   6c. Compute
        Ex  at x(i), y(j+1/2), z(k+1/2)  (Cartesian or Polar) 
       -Er  at x(i), y(j+1/2), z(k+1/2)  (Spherical) 
   ---------------------------------------------- */

#if GEOMETRY == CARTESIAN || GEOMETRY == POLAR
  for (k = KBEG-1; k <= KEND; k++) IDOM_LOOP(i){
    int BS = BX3;
    double ***bs = Us[BX3s];
#else
  for (j = JBEG-1; j <= JEND; j++) IDOM_LOOP(i){
    int BS = BX2;
    double ***bs = Us[BX2s];
#endif
  
    #if GEOMETRY == CARTESIAN
    w = 0.5*(wA[k][i] + wA[k+1][i]);
    #elif GEOMETRY == POLAR
    w = 0.5*(wA[k][i] + wA[k+1][i])/x1[i];
    #elif GEOMETRY == SPHERICAL
    w = 0.5*(wA[j][i] + wA[j+1][i])/(x1[i]*sin(x2p[j]));
    #endif

  /* --------------------------------------------
     6d. Obtain shift in terms of an integer
         number of cells (m) plus a remainder eps. 
     -------------------------------------------- */

    dL  = w*g_dt;                 /* spatial shift */
    dL  = fmod(dL, Lphi);
    m   = (int)floor(dL/dphi + 0.5); /* nearest integer part */
    eps = dL/dphi - m;               /* remainder in units of dphi */

  /* -- copy 3D array into 1D buffer -- */

    SDOM_LOOP(s) q[s] = -bs[k][j][i];
    if (nproc_s > 1){    /* -- copy values from lower and upper buffers -- */
      #ifdef PARALLEL
      for (s = SBEG-1; s >= SBEG - nbuf_upper; s--){
        q[s] = -FARGO_ARRAY_INDEX(upper_buf[BS], SBEG-1-s, k, j, i); 
      }
      for (s = SEND+1; s <= SEND + nbuf_lower; s++){
        q[s] = -FARGO_ARRAY_INDEX(lower_buf[BS], s-SEND-1, k, j, i);
      }
      FARGO_Flux (q-m, flux-m, eps);
      #endif
    }else{               /* -- assign periodic b.c. -- */
      for (s = 1; s <= ngh; s++) {                        
        q[SBEG - s] = q[SBEG - s + NS]; 
        q[SEND + s] = q[SEND + s - NS]; 
      }    
      FARGO_Flux (q, flux, eps);
    }

  /* -- update main solution array -- */

    if (nproc_s > 1){
      #ifdef PARALLEL
      if (m >= 0){
        for (s = SBEG-1; s <= SEND; s++) {
          sm = s - m;
          Ex[k][j][i] = eps*flux[sm];
          for (ss = s; ss > (s-m); ss--) Ex[k][j][i] += q[ss];
          Ex[k][j][i] *= dphi;
        }
      }else {
        for (s = SBEG-1; j <= SEND; s++) {
          sm = s - m;
          Ex[k][j][i] = eps*flux[sm];
          for (ss = s+1; ss < (s-m+1); ss++) Ex[k][j][i] -= q[ss];
          Ex[k][j][i] *= dphi;
        }
      }
      #endif
    }else{
      if (m >= 0){
        for (s = SBEG-1; s <= SEND; s++) {
          sm = FARGO_MOD(s - m);
          Ex[k][j][i] = eps*flux[sm];
          for (ss = s; ss > (s-m); ss--) Ex[k][j][i] += q[FARGO_MOD(ss)];
          Ex[k][j][i] *= dphi;
        }
      }else {
        for (s = SBEG-1; j <= SEND; s++) {
          sm = FARGO_MOD(s - m);
          Ex[k][j][i] = eps*flux[sm];
          for (ss = s+1; ss < (s-m+1); ss++) Ex[k][j][i] -= q[FARGO_MOD(ss)];
          Ex[k][j][i] *= dphi;
        }
      }
    }
  }
#endif /*  DIMENSIONS == 3 */

/* ----------------------------------------------
   6e. Update bx1 staggered magnetic field
   ---------------------------------------------- */

  KDOM_LOOP(k) JDOM_LOOP(j) for (i = IBEG-1; i <= IEND; i++){
    #if GEOMETRY == CARTESIAN || GEOMETRY == POLAR
    bx1[k][j][i] -= (Ez[k][j][i] - Ez[k][j-1][i])/dphi;
    #elif GEOMETRY == SPHERICAL /* in spherical coordinates Ez = -Eth */
//    bx1[k][j][i] -= (Ez[k][j][i] - Ez[k-1][j][i])/dphi;
    double dmu = cos(x2m[j]) - cos(x2p[j]);
    bx1[k][j][i] -= dx2[j]*sin(x2[j])/(dmu*dphi)*(Ez[k][j][i] - Ez[k-1][j][i]);
    #endif
  }

/* ----------------------------------------------
   6f. Update bx2 staggered magnetic field
   ---------------------------------------------- */

  KDOM_LOOP(k) for (j = JBEG-1; j <= JEND; j++) IDOM_LOOP(i){
    #if GEOMETRY == CARTESIAN
    bx2[k][j][i] += DIM_EXPAND(                                    ,
                           (Ez[k][j][i] - Ez[k][j][i-1])/dx1[i]  ,
                         - (Ex[k][j][i] - Ex[k-1][j][i])/dx3[k]);
    #elif GEOMETRY == POLAR
     bx2[k][j][i] += DIM_EXPAND(                                             ,
                       (x1p[i]*Ez[k][j][i] -  x1p[i-1]*Ez[k][j][i-1])/dx1[i]  ,
                       - x1[i]*(Ex[k][j][i] - Ex[k-1][j][i])/dx3[k]); 
    #elif GEOMETRY == SPHERICAL /* in spherical coordinates Ex = -Er */
     bx2[k][j][i] += (Ex[k][j][i] - Ex[k-1][j][i])/dphi;
    #endif
  }

/* ----------------------------------------------
   6g. Update bx3 staggered magnetic field
   ---------------------------------------------- */

#if INCLUDE_KDIR
  for (k = KBEG-1; k <= KEND; k++) JDOM_LOOP(j) IDOM_LOOP(i){
    #if GEOMETRY == CARTESIAN || GEOMETRY == POLAR
    bx3[k][j][i] += (Ex[k][j][i] - Ex[k][j-1][i])/dphi;
    #elif GEOMETRY == SPHERICAL /* Ez = -Eth, Ez = -Er */
    A1p = x1p[i]*x1p[i];
    A1m = x1m[i]*x1m[i];
    A2p = fabs(sin(x2p[j]));
    A2m = fabs(sin(x2m[j]));

    bx3[k][j][i] += sin(x2[j])*(   A1p*Ez[k][j][i] 
                                -  A1m*Ez[k][j][i-1])/(x1[i]*dx1[i])
                    - (A2p*Ex[k][j][i] - A2m*Ex[k][j-1][i])/dx2[j];
    #endif
  }
#endif

/* -- redefine cell-centered field -- */
  
  DOM_LOOP(k,j,i){
    DIM_EXPAND(
      U[k][j][i][BX1] = 0.5*(Us[BX1s][k][j][i] + Us[BX1s][k][j][i-1]);  ,
      U[k][j][i][BX2] = 0.5*(Us[BX2s][k][j][i] + Us[BX2s][k][j-1][i]);  ,
      U[k][j][i][BX3] = 0.5*(Us[BX3s][k][j][i] + Us[BX3s][k-1][j][i]);
    )
  }
}
#endif /* STAGGERED_MHD */

  #ifdef PARALLEL
   if (nproc_s > 1){
     for (nv = 0; nv < NVAR; nv++){
       FreeArrayBox(upper_buf[nv], ubox.kbeg, ubox.jbeg, ubox.ibeg);
       FreeArrayBox(lower_buf[nv], lbox.kbeg, lbox.jbeg, lbox.ibeg);
     }
   }
  #endif
}

#ifdef PARALLEL
/* ********************************************************************* */
void FARGO_ParallelExchange (Data_Arr U, Data_Arr Us,
                             double ***recv_buf_lower[], RBox *lbox,
                             double ***recv_buf_upper[], RBox *ubox,
                             Grid *grid)
/*!
 * Exchange data values between  adjacent processors lying in the
 * the same direction (x2 for CARTESIAN/POLAR or x3 for SPHERICAL)
 * Values will be stored inside recv_buf_lower & recv_buf_upper.
 *
 * \param [in,out]  U    a 3D array of conserved, zone-centered values
 * \param [in,out]  Us   a 3D array of staggered magnetic fields
 * \param [out]     recv_buf_lower buffer array received from the 
 *                                 processor ahead
 * \param [out]     recv_buf_upper buffer array received from the 
 *                                 processor behind
 * \param [in]      ubox           pointer to RBox structure defining 
 *                                 the index range of the upper data layer
 * \param [in]      lbox           pointer to RBox structure defining 
 *                                 the index range of the lower data layer
 * \param [in]      grid  pointer to Grid structure;
 *
 * \return  This function has no return value.
 *********************************************************************** */
{
  int   nv, i, j, k;
  int   ib, jb, kb;
  int   ie, je, ke;
  int   coords[3], dst, src;

  long int  ntot;
  double ***send_buf[NVAR], *sbuf, *rbuf;
  MPI_Comm cartcomm;
  MPI_Status status;

/* -------------------------------------------------------------------
    get ranks of the processors lying above (upper) and below (lower)
   ------------------------------------------------------------------ */

  AL_Get_cart_comm(SZ1, &cartcomm);
  for (i = 0; i < DIMENSIONS; i++) coords[i] = grid->rank_coord[i];
  coords[SDIR] += 1;
  MPI_Cart_rank(cartcomm, coords, &dst);

  for (i = 0; i < DIMENSIONS; i++) coords[i] = grid->rank_coord[i];
  coords[SDIR] -= 1;
  MPI_Cart_rank(cartcomm, coords, &src);

/* -------------------------------------------------------------
                Send/Recv of the upper layer
   ------------------------------------------------------------- */

  /* -- determine buffer size & dimensions -- */

  ib = ubox->ibeg; jb = ubox->jbeg; kb = ubox->kbeg;
  ie = ubox->iend; je = ubox->jend; ke = ubox->kend;

  ntot  = (ie - ib + 1)*(je - jb + 1)*(ke - kb + 1);
  for (nv = 0; nv < NVAR; nv++){
    send_buf[nv] = ARRAY_BOX(kb, ke, jb, je, ib, ie, double);
  }

  /* -- copy data -- */

  BOX_LOOP(ubox, k, j, i){
    #if GEOMETRY == CARTESIAN || GEOMETRY == POLAR
/* send_buf[nv][k][j][i] = FARGO_ARRAY_INDEX(U, SEND-s, k,j,i); */

     for (nv = NVAR; nv--; ) send_buf[nv][k][j][i] = U[k][JEND-j][i][nv]; 
     #ifdef STAGGERED_MHD
      DIM_EXPAND(send_buf[BX1][k][j][i] = Us[BX1s][k][JEND-j][i];  ,
                                                            ;  ,
               send_buf[BX3][k][j][i] = Us[BX3s][k][JEND-j][i];)
     #endif

    #elif GEOMETRY == SPHERICAL

     for (nv = NVAR; nv--; ) send_buf[nv][k][j][i] = U[KEND-k][j][i][nv]; 
     #ifdef STAGGERED_MHD
      DIM_EXPAND(send_buf[BX1][k][j][i] = Us[BX1s][KEND-k][j][i];  ,
               send_buf[BX2][k][j][i] = Us[BX2s][KEND-k][j][i];  ,
                                                            ;)
     #endif

    #endif
  }

/* -- data communication -- */

  for (nv = 0; nv < NVAR; nv++){
    sbuf = &(send_buf[nv][kb][jb][ib]);
    rbuf = &(recv_buf_upper[nv][kb][jb][ib]);
    MPI_Sendrecv(sbuf, ntot, MPI_DOUBLE, dst, 1,
                 rbuf, ntot, MPI_DOUBLE, src, 1, cartcomm, &status);
  }

/* -- free send buffer memory -- */

  for (nv = 0; nv < NVAR; nv++) FreeArrayBox(send_buf[nv], kb, jb, ib);

/* -------------------------------------------------------------
                  Send/Recv of the lower layer
   ------------------------------------------------------------- */

  /* -- determine buffer size & dimensions -- */

  ib = lbox->ibeg; jb = lbox->jbeg; kb = lbox->kbeg;
  ie = lbox->iend; je = lbox->jend; ke = lbox->kend;

  ntot  = (ie - ib + 1)*(je - jb + 1)*(ke - kb + 1);

  for (nv = 0; nv < NVAR; nv++){
    send_buf[nv] = ARRAY_BOX(kb, ke, jb, je, ib, ie, double);
  }

  BOX_LOOP(lbox, k, j, i) {
    #if GEOMETRY == CARTESIAN || GEOMETRY == POLAR

     for (nv = NVAR; nv--;  ) send_buf[nv][k][j][i] = U[k][JBEG+j][i][nv]; 
     #ifdef STAGGERED_MHD
      DIM_EXPAND(send_buf[BX1][k][j][i] = Us[BX1s][k][JBEG+j][i];  ,
                                                                ;  ,
                 send_buf[BX3][k][j][i] = Us[BX3s][k][JBEG+j][i];)
     #endif

    #elif GEOMETRY == SPHERICAL

     for (nv = NVAR; nv--;  ) send_buf[nv][k][j][i] = U[KBEG+k][j][i][nv]; 
     #ifdef STAGGERED_MHD
      DIM_EXPAND(send_buf[BX1][k][j][i] = Us[BX1s][KBEG+k][j][i];  ,
                 send_buf[BX2][k][j][i] = Us[BX2s][KBEG+k][j][i];  ,
                                                                ;)
     #endif

    #endif
  }

  /* -- data communication -- */

  for (nv = 0; nv < NVAR; nv++){
    sbuf = &(send_buf[nv][kb][jb][ib]);
    rbuf = &(recv_buf_lower[nv][kb][jb][ib]);
    MPI_Sendrecv(sbuf, ntot, MPI_DOUBLE, src, 100,
                 rbuf, ntot, MPI_DOUBLE, dst, 100, cartcomm, &status);
  }

  /* -- free send buffer memory -- */

  for (nv = 0; nv < NVAR; nv++) FreeArrayBox(send_buf[nv], kb, jb, ib);
}
#endif /* PARALLEL */

/* ********************************************************************* */
void FARGO_Flux (double *q, double *flx, double eps)
/*!
 *  Compute the transport flux corresponding to a fractional 
 *  increment eps.
 *
 * \param [in]  q   1D array of a conserved quantity
 * \param [out] flx the conservative flux
 * \param [in]  eps the fractional increment
 *
 * The following MAPLE script may be used to check the 
 * correctness of implementation
 * \code
   restart;
   P0 := 3/2*q[i] - (q[p]+q[m])/4;
   P1 := q[p] - q[m];
   P2 := 3*(q[p] - 2*q[i] + q[m]);
 
   FS[p] := P0 + P1/2*(1 - epsilon) + P2*(1/4 - epsilon/2 + epsilon^2/3);
   FS[m] := P0 - P1/2*(1 - epsilon) + P2*(1/4 - epsilon/2 + epsilon^2/3);
 
   dq  := q[p] - q[m];
   d2q := q[p] - 2*q[i] + q[m];
   FF[p] := q[p] - epsilon/2*(dq + d2q*(3 - 2*epsilon));
   FF[m] := q[m] + epsilon/2*(dq - d2q*(3 - 2*epsilon));
 * \endcode
 *
 ************************************************************************* */
{
  int    s;
  double dqc, d2q, dqp, dqm;
  static double *dq, *dql, *qp, *qm;
  
  if (dq == NULL){
    dq  = ARRAY_1D(NMAX_POINT, double);
    dql = ARRAY_1D(NMAX_POINT, double);
  }

/* -- compute limited slopes -- */

  dq[0] = q[1] - q[0];
  for (s = 1; s < NS_TOT-1; s++){
    dq[s] = q[s+1] - q[s];
    if (dq[s]*dq[s-1] > 0.0){
      dqc = 0.5*(dq[s] + dq[s-1]);
      d2q = 2.0*(fabs(dq[s]) < fabs(dq[s-1]) ? dq[s]:dq[s-1]);
      dql[s] = fabs(d2q) < fabs(dqc) ? d2q: dqc;
    }else dql[s] = 0.0;
  } 

/* vanleer_lim(dq,dq-1,dql,JBEG-1,JEND+1,NULL);  */

  #if FARGO_ORDER == 2  /* MUSCL-Hancock, 2nd order */
   if (eps > 0.0){
     for (s=SBEG-1; s<=SEND; s++) flx[s] = q[s] + 0.5*dql[s]*(1.0 - eps);
   }else { /* -- remember: eps > 0 by definition -- */
     for (s=SBEG-1; s<=SEND; s++) flx[s] = q[s+1] - 0.5*dql[s+1]*(1.0 + eps);
   }
  #elif FARGO_ORDER == 3  /* PPM, 3rd order */
   if (qp == NULL){
     qp = ARRAY_1D(NMAX_POINT, double);
     qm = ARRAY_1D(NMAX_POINT, double);
   }
   for (s = SBEG-1; s <= SEND+1; s++){
     dqp =  0.5*dq[s]   - (dql[s+1] - dql[s])/6.0;
     dqm = -0.5*dq[s-1] + (dql[s-1] - dql[s])/6.0;
     if (dqp*dqm > 0.0) dqp = dqm = 0.0;
     else{
       if (fabs(dqp) >= 2.0*fabs(dqm)) dqp = -2.0*dqm;
       if (fabs(dqm) >= 2.0*fabs(dqp)) dqm = -2.0*dqp;
     }
     qp[s] = q[s] + dqp;
     qm[s] = q[s] + dqm;
   }

   if (eps > 0.0){
     for (s = SBEG-1; s <= SEND; s++){
       dqc = qp[s] - qm[s];
       d2q = qp[s] - 2.0*q[s] + qm[s];
       flx[s] = qp[s] - 0.5*eps*(dqc + d2q*(3.0 - 2.0*eps));
     }
   }else{  
     for (s = SBEG-1; s <= SEND; s++){
       dqc = qp[s+1] - qm[s+1];
       d2q = qp[s+1] - 2.0*q[s+1] + qm[s+1];
       flx[s] = qm[s+1] - 0.5*eps*(dqc - d2q*(3.0 + 2.0*eps));
     }
   }
  #else
   printLog ("! FARGO_ORDER should be either 2 or 3\n");
   QUIT_PLUTO(1);
  #endif

}

#if PARTICLES
/* ********************************************************************* */
void FARGO_ShiftParticles(Data *d, Grid *grid, double dt)
/*!
 *
 *
 *
 *********************************************************************** */
{
  int    i,  j,  k;
  int    i1, j1, k1;
  double w, v, **wA;
  static double ***wp;
  particleNode *curNode;
  Particle *p;

  #if !INCLUDE_JDIR
  return;
  #endif
  
/* -----------------------------------------------------
   0. Allocate memory
   ----------------------------------------------------- */

  if (wp == NULL){
    wp = ARRAY_BOX (-1, 1, -1, 1, -1, 1, double);
  }

/* -----------------------------------------------------
   1. Shift particles in the orbital direction
   ----------------------------------------------------- */
  
  wA = FARGO_Velocity();
  PARTICLES_LOOP(curNode, d->PHead){
    p = &(curNode->p);

  /* -- Interpolate fargo shift velocity at particle position -- */

    Particles_GetWeights(p, p->cell, wp, grid);  
    i = p->cell[IDIR];
    j = p->cell[JDIR];
    k = p->cell[KDIR];

    v = 0.0;
    for (k1 = -INCLUDE_KDIR; k1 <= INCLUDE_KDIR; k1++){
    for (j1 = -INCLUDE_JDIR; j1 <= INCLUDE_JDIR; j1++){
    for (i1 = -INCLUDE_IDIR; i1 <= INCLUDE_IDIR; i1++){
      v += wp[k1][j1][i1]*wA[k + k1][i + i1];
    }}}
// v = -SB_Q*SB_OMEGA*p->coord[IDIR];
    p->coord[JDIR] += v*dt;
  }

/* -----------------------------------------------------
   2. Shift particles in the orbital direction
   ----------------------------------------------------- */

  Particles_Boundary(d, grid);
  Particles_BoundaryExchange(d, grid);
}
#endif
#undef MAX_BUF_SIZE
