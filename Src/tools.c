/* ///////////////////////////////////////////////////////////////////// */
/*!
  \file
  \brief Collection of general-purpose functions.

  \author A. Mignone (andrea.mignone@unito.it)
  \date   Jan 27, 2020
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"
   
/* ********************************************************************* */
void GetNeighbourRanks (Grid *grid, int **nranks)
/*! 
 *  Find the ranks of neighbour processors direction by direction.
 *  Store them in the array \c nrank[dir][s] where \c dir is the
 *  direction and <tt> s = 0,1 </tt> stands for the left (0) or right
 *  (1) processor.
 *  If a processor touches the physical boundary that is not periodic,
 *  the corresponding neighbour rank will be set to -1.
 *
 * \param [in]   grid    pointer to an array of grid structures.
 * \param [out]  nranks  an array of integers containing the ranks of
 *                       the neighbouring processors
 *********************************************************************** */
{
#ifdef PARALLEL
  int dir, coords[3], rnk;
  int lbound, rbound;
  MPI_Comm cartcomm;
  

  AL_Get_cart_comm(SZ1, &cartcomm);
    
/* ------------------------------------------------------------
    Neighbour exists when there's no physical boundary, or
    when PERIODIC or SHEARING boundary conditions are imposed
    along that direction.
   ------------------------------------------------------------ */

  for (dir = 0; dir < 3;dir++){

  /* -- Obtain local processor coordinates -- */

    DIM_EXPAND(coords[IDIR] = grid->rank_coord[IDIR];  ,
              coords[JDIR] = grid->rank_coord[JDIR];  ,
              coords[KDIR] = grid->rank_coord[KDIR];)

    nranks[dir][0] = nranks[dir][1] = -1;

    lbound = grid->lbound[dir];
    if (lbound == 0 || lbound == PERIODIC || lbound == SHEARING){
      coords[dir] = grid->rank_coord[dir] - 1;
      MPI_Cart_rank(cartcomm, coords, &rnk);
      nranks[dir][0] = rnk;
    }

    rbound = grid->rbound[dir];

    if (rbound == 0 || rbound == PERIODIC || rbound == SHEARING){
      coords[dir] = grid->rank_coord[dir] + 1;
      MPI_Cart_rank(cartcomm, coords, &rnk);
      nranks[dir][1] = rnk;
    }
  }
  
  MPI_Barrier(MPI_COMM_WORLD);

  return;
#endif  /* PARALLEL */
}

#if !(defined CHOMBO) && (defined GNUPLOT_HEADER)
/* ********************************************************************* */
void GnuplotSetting(Runtime *runtime, Grid *grid)
/*
 * Set-up a gnuplot script containing grid info, variable names, etc...
 *********************************************************************** */
{
  int  k, nv;
  int  dbl_out = 0;
  int  flt_out = 0;
  char *var_name[256];
  Output *output;
  time_t  now;
  FILE *fp;

  if (prank != 0) return;

/* --------------------------------------------------------
   0. Allocate memory for an array of strings containing
      the names of the variables.
      Get the variable names being written to disk.
   -------------------------------------------------------- */

  for (nv = 0; nv < 255; nv++) var_name[nv] = ARRAY_1D(32,char);

/* --------------------------------------------------------
   1. Determine which output datatypes have been enabled
   -------------------------------------------------------- */

  for (k = 0; k < MAX_OUTPUT_TYPES; k++){
    output = runtime->output + k;
    if (output->type == DBL_OUTPUT && (output->dt > 0 || output->dn > 0)){
      dbl_out = 1;
    }
    if (output->type == FLT_OUTPUT && (output->dt > 0 || output->dn > 0)){
      flt_out = 1;
    }
  }

/* --------------------------------------------------------
   2. Open file and write relevant info.
   -------------------------------------------------------- */
  
  time(&now);
  fp = fopen("pluto.gp","w");

  fprintf (fp, "# ***********************************************************\n");
  fprintf (fp, "#\n");
  fprintf (fp, "# Initialize simulation quantities for reading binary data   \n");
  fprintf (fp, "# with Gnuplot. \n");
  fprintf (fp, "# Generated by PLUTO on  %s",asctime(localtime(&now)));
  fprintf (fp, "#\n");
  fprintf (fp, "# ***********************************************************\n");

/* --------------------------------------------------------
   2a. Select default data type and index name
   -------------------------------------------------------- */

  if (dbl_out){
    fprintf (fp, "if (!exists('dtype')) {\n");
    fprintf (fp,"  dtype = 'dbl' # Default data type\n");
    fprintf (fp,"}\n");
  }else if (flt_out){
    fprintf (fp, "if (!exists('dtype')) {\n");
    fprintf (fp,"  dtype = 'flt' # Default data type\n");
    fprintf (fp,"}\n");
  }
  fprintf (fp,"print '> dtype set to ',dtype\n");
 
/* --------------------------------------------------------
   2b. Set grid size
   -------------------------------------------------------- */

  fprintf (fp, "\n");
  fprintf (fp, "# --------------------------------------- \n");
  fprintf (fp, "# Grid setting                            \n");
  fprintf (fp, "# --------------------------------------- \n");

  fprintf (fp, "nx1 = %d\n", grid->np_int_glob[IDIR]);
  fprintf (fp, "nx2 = %d\n", grid->np_int_glob[JDIR]);
  fprintf (fp, "nx3 = %d\n", grid->np_int_glob[KDIR]);
  
  fprintf (fp, "x1beg = %f; x1end = %f\n", grid->xbeg_glob[IDIR],
                                           grid->xend_glob[IDIR]);

  fprintf (fp, "x2beg = %f; x2end = %f\n", grid->xbeg_glob[JDIR],
                                           grid->xend_glob[JDIR]);

  fprintf (fp, "x3beg = %f; x3end = %f\n", grid->xbeg_glob[KDIR],
                                           grid->xend_glob[KDIR]);
  fprintf (fp, "dx1   = %f\n", grid->dx[IDIR][IBEG]);
  fprintf (fp, "dx2   = %f\n", grid->dx[JDIR][JBEG]);
  fprintf (fp, "dx3   = %f\n", grid->dx[KDIR][KBEG]);

  fprintf (fp, "\n");

/* --------------------------------------------------------
   2c. Variable names
   -------------------------------------------------------- */

  fprintf (fp, "# --------------------------------------- \n");
  fprintf (fp, "# Variable indices                        \n");
  fprintf (fp, "# --------------------------------------- \n");

  if (flt_out){
    k = GetOutputVarNames(FLT_OUTPUT, var_name);
    fprintf (fp, "if (dtype eq 'flt' || dtype eq 'tab') {\n");
    for (nv = 0; nv < k; nv++) {
      fprintf (fp, "  %s = %d\n",var_name[nv], nv);
    }
    char all_var[512] = {"> Var Names: "};
    for (nv = 0; nv < k; nv++) {
      strcat (all_var, var_name[nv]);
      strcat (all_var, "  ");
    }
    fprintf (fp, "  print '%s'\n",all_var);
    fprintf (fp, "}\n");
  }

  if (dbl_out){
    k = GetOutputVarNames(DBL_OUTPUT, var_name);
    fprintf (fp, "if (dtype eq 'dbl') {\n");
    for (nv = 0; nv < k; nv++) {
      fprintf (fp, "  %s = %d\n",var_name[nv], nv);
    }
    char all_var[512] = {"> Var Names: "};
    for (nv = 0; nv < k; nv++) {
      strcat (all_var, var_name[nv]);
      strcat (all_var, "  ");
    }
    fprintf (fp, "  print '%s'\n",all_var);
    fprintf (fp, "}\n");
  }
  
  fprintf (fp, "if (!exists('nvar')) {\n");
  fprintf (fp, "  nvar = 0\n");
  fprintf (fp, "  print '> nvar = ',nvar\n");
  fprintf (fp, "}\n");

/* --------------------------------------------------------
   3. Make gnuplot count the number of lines in "flt.out"
      and "dbl.out"
   -------------------------------------------------------- */

  fprintf (fp, "\n");
  fprintf (fp, "# --------------------------------------- \n");
  fprintf (fp, "# Output datatype                         \n");
  fprintf (fp, "# --------------------------------------- \n");
  fprintf (fp,"#set term push       # remember current terminal\n");
  fprintf (fp,"#set term unknown\n");

  if (dbl_out){
    fprintf (fp,"#plot 'dbl.out' u 1\n");
    fprintf (fp,"#ndbl = int(GPVAL_DATA_Y_MAX)\n");
  }
  if (flt_out){
    fprintf (fp,"#plot 'flt.out' u 1\n");
    fprintf (fp,"#nflt = int(GPVAL_DATA_Y_MAX)\n");
  }
  fprintf (fp,"#set term pop         # restore terminal\n");

/* --------------------------------------------------------
   4. Load some other gnuplot scripts
   -------------------------------------------------------- */

  fprintf (fp, "\n");
  fprintf (fp, "# --------------------------------------- \n");
  fprintf (fp, "# Load gnuplot scripts                    \n");
  fprintf (fp, "# --------------------------------------- \n");
  fprintf (fp, "load 'macros.gp'\n");
  fprintf (fp, "if (nx2 > 1) load 'pm3d_setting.gp'\n");

  fclose(fp);
}
#endif

/* ********************************************************************* */
int IsLittleEndian (void) 
/*!
 * Return 1 if the current architecture has little endian order
 *
 *********************************************************************** */
{
  int TestEndian = 1;
  return *(char*)&TestEndian;
}

/* ********************************************************************* */
void MakeState (Sweep *sweep)
/*!
 * Allocate memory areas for arrays inside the sweep
 * structure.
 *********************************************************************** */
{
  State *stateC = &(sweep->stateC);
  State *stateL = &(sweep->stateL);
  State *stateR = &(sweep->stateR);

/* --------------------------------------------------------
   0. Allocate memory for sweep structure members
   -------------------------------------------------------- */

  sweep->vn      = ARRAY_2D(NMAX_POINT, NVAR, double);
  sweep->flux    = ARRAY_2D(NMAX_POINT, NVAR, double);
  sweep->src     = ARRAY_2D(NMAX_POINT, NVAR, double);

  sweep->tc_flux = ARRAY_2D(NMAX_POINT, NVAR, double);    

  sweep->rhs     = ARRAY_2D(NMAX_POINT, NVAR, double);
  sweep->press   = ARRAY_1D(NMAX_POINT, double);
  sweep->Bn      = ARRAY_1D(NMAX_POINT, double);
  sweep->En      = ARRAY_1D(NMAX_POINT, double);
  sweep->SL      = ARRAY_1D(NMAX_POINT, double);
  sweep->SR      = ARRAY_1D(NMAX_POINT, double);
  #if RADIATION
  sweep->SrL     = ARRAY_1D(NMAX_POINT, double);
  sweep->SrR     = ARRAY_1D(NMAX_POINT, double);
  #endif

  sweep->pnt_flux = ARRAY_2D(NMAX_POINT, NVAR, double);
  sweep->dff_flux = ARRAY_2D(NMAX_POINT, NVAR, double);
  sweep->SaL      = ARRAY_1D(NMAX_POINT, double);
  sweep->SaR      = ARRAY_1D(NMAX_POINT, double);
  sweep->Sc       = ARRAY_1D(NMAX_POINT, double);
  sweep->dL       = ARRAY_1D(NMAX_POINT, double);
  sweep->dR       = ARRAY_1D(NMAX_POINT, double);
  sweep->aL       = ARRAY_1D(NMAX_POINT, double);
  sweep->aR       = ARRAY_1D(NMAX_POINT, double);

/* -- eigenvectors -- */

  sweep->lmax    = ARRAY_1D(NVAR, double);
  sweep->flag    = ARRAY_1D(NMAX_POINT, uint16_t);

/* --------------------------------------------------------
   1. Allocate memory for state structure members.
      C/L stand, respectively, for cell center and
      left interfaces (i+1/2, L).
   -------------------------------------------------------- */
  
  StateStructAllocate (stateC);
  StateStructAllocate (stateL);

/* --------------------------------------------------------
   2. Allocate memory for the right state structure.
      Note that we add an offset +1 in order to access
      stateR->v[i-1].
   -------------------------------------------------------- */
  
  stateR->v      = ARRAY_2D(NMAX_POINT, NVAR, double)+1;
  stateR->u      = ARRAY_2D(NMAX_POINT, NVAR, double)+1;
  stateR->flux   = ARRAY_2D(NMAX_POINT, NVAR, double)+1;
  stateR->lambda = ARRAY_2D(NMAX_POINT, NVAR, double)+1;
  stateR->prs    = ARRAY_1D(NMAX_POINT, double)+1;
  stateR->a2     = ARRAY_1D(NMAX_POINT, double)+1;
  stateR->cw     = ARRAY_1D(NMAX_POINT, double)+1;
  stateR->h      = ARRAY_1D(NMAX_POINT, double)+1;
  stateR->Lp     = ARRAY_3D(NMAX_POINT, NFLX, NFLX, double)+1;
  stateR->Rp     = ARRAY_3D(NMAX_POINT , NFLX, NFLX, double)+1;
  stateR->J      = stateL->J;
  stateR->cCR    = stateL->cCR;
  stateR->Fcr    = ARRAY_2D(NMAX_POINT, 4, double) + 1; 
  stateR->fluxCR = ARRAY_2D(NMAX_POINT, NVAR, double)+1;
  stateR->Bbck   = stateL->Bbck;
}

/* ********************************************************************* */
void PlutoError (int condition, char *str)
/*!
 * If condition is true, issue an error and quit the code.
 *
 *********************************************************************** */
{
  char *str_err="! Error: ";

  if (condition) {
    printLog (str_err);
    printLog (str);
    printLog ("\n");
    QUIT_PLUTO(1);
  }
}

/* ********************************************************************* */
void PrintColumnLegend(char *legend[], int nfields, FILE *fp)
/*!
 *  Print nicely formatted legend to an ASCII datfile header.
 *
 * \param [in] legend    an array of strings describing the fields
 * \param [in] nfields   the number of fields (also the dimensions of
 *                       legend)
 * \param [in] fp        a valid file descriptor.
 *********************************************************************** */
{
  int n;
  char num_string[15];

  for (n = 0; n < nfields; n++){
    sprintf (num_string,"[%d]",n);
    if (n == 0) fprintf (fp,"# %-13s ",num_string);
    else        fprintf (fp,"%-13s  ",num_string);
  }
  fprintf (fp,"\n");

  for (n = 0; n < nfields; n++){
    if (n == 0) fprintf (fp,"# %-13s ",legend[n]);
    else        fprintf (fp,"%-13s  ",legend[n]);
  }
  fprintf (fp,"\n");
}


/* ********************************************************************* */
void StateStructAllocate (State *p)
/*!
 * Allocate memory areas for arrays inside the State structure.
 *********************************************************************** */
{
  p->v      = ARRAY_2D(NMAX_POINT, NVAR, double);
  p->u      = ARRAY_2D(NMAX_POINT, NVAR, double);
  p->flux   = ARRAY_2D(NMAX_POINT, NVAR, double);
  p->lambda = ARRAY_2D(NMAX_POINT, NVAR, double);
  p->prs    = ARRAY_1D(NMAX_POINT, double);
  p->a2     = ARRAY_1D(NMAX_POINT, double);
  p->cw     = ARRAY_1D(NMAX_POINT, double);
  p->h      = ARRAY_1D(NMAX_POINT, double);
  p->Lp     = ARRAY_3D(NMAX_POINT, NFLX, NFLX, double);
  p->Rp     = ARRAY_3D(NMAX_POINT, NFLX, NFLX, double);
  p->J      = ARRAY_2D(NMAX_POINT, 3, double);
  p->cCR    = ARRAY_2D(NMAX_POINT, 3, double);
  p->Fcr    = ARRAY_2D(NMAX_POINT, 4, double); 
  p->fluxCR = ARRAY_2D(NMAX_POINT, NVAR, double);
  p->Bbck   = ARRAY_2D(NMAX_POINT, 3, double);
}

/* ********************************************************************* */
int StringArrayIndex (char *str_arr[], int size, char *str)
/*!
 *  Find the index of an array of strings whose value matches the
 *  input string \c *str.
 *
 * \param [in]  str_arr  the array of strings
 * \param [in]  size     the size of the array
 * \param [in]  str      the string to search for
 *********************************************************************** */
{
  int i;

  for (i = 0; i < size; i++){
    if (!strcmp(str_arr[i],str)) return i;
  }
  return -1;
}

/* ********************************************************************* */
void SymmetryCheck (Data_Arr V, int where, RBox *box)
/*!
 * Check if vectors preserve symmetry / antisymmetry properties
 * at a symmetry axis (r = 0 for cylindrical or theta = 0 for spherical)
 * Vectors can only be cell-centered (meaning that *all* components must
 * be cell-centered) or face-centered (meaning that *all* components must
 * be face-centered).
 * Staggered fields are not allowed.
 *********************************************************************** */
{
#if GEOMETRY == CARTESIAN
  return;
#else
  int i,j,k;
  int ip, jp;
  double ***Vx1 = V[0];
  double ***Vx2 = V[1];
  double ***Vx3 = V[2];
  double scrh1 = 0.0, scrh2 = 0.0, scrh3 = 0.0;

  BOX_LOOP(box,k,j,i){

    #if GEOMETRY == CYLINDRICAL
    if (i >= IBEG) continue;

    if (where == X1FACE){   /* Vector is defined at radial i+1/2 interfaces */
      ip = 2*IBEG-i-2;
      jp = j;

    /* At IBEG-1/2 interface, check that there're no Inf or Nan values */
      if (i == IBEG-1){
        if (   fabs(Vx1[k][j][i]) > 1.e8 || Vx1[k][j][i] != Vx1[k][j][i]
            || fabs(Vx2[k][j][i]) > 1.e8 || Vx2[k][j][i] != Vx2[k][j][i]
            || fabs(Vx3[k][j][i]) > 1.e8 || Vx3[k][j][i] != Vx3[k][j][i]){
          printLog ("! SymmetryCheck(): invalid value\n");
          QUIT_PLUTO(1);
        }
        continue; /* Skip this point since i = ip would be the same */
      }
    }else {                /* Vector is defined at center or j+1/2 interface */
      ip = 2*IBEG-i-1;
      jp = j;
    }
    scrh1 = fabs(Vx1[k][j][i] + Vx1[k][jp][ip]); /* r-component:   antisymmetric */
    scrh2 = fabs(Vx2[k][j][i] - Vx2[k][jp][ip]); /* z-component:   symmetric */
    scrh3 = fabs(Vx3[k][j][i] + Vx3[k][jp][ip]); /* phi-component: antisymmetric */
    #elif GEOMETRY == SPHERICAL
    if (j >= JBEG) continue;

    if (where == X2FACE){  /* Vector is define at meridional j+1/2 interfaces */
      ip = i;
      jp = 2*JBEG-j-2;
    /* At JBEG-1/2 interface, check that there're no Inf or Nan values */
      if (j == JBEG-1){
        if (   fabs(Vx1[k][j][i]) > 1.e8 || Vx1[k][j][i] != Vx1[k][j][i]
            || fabs(Vx2[k][j][i]) > 1.e8 || Vx2[k][j][i] != Vx2[k][j][i]
            || fabs(Vx3[k][j][i]) > 1.e8 || Vx3[k][j][i] != Vx3[k][j][i]){
          printLog ("! SymmetryCheck(): invalid value\n");
          printLog ("! V = (%12.6e, %12.6e, %12.6e\n",
                    Vx1[k][j][i],Vx2[k][j][i],Vx3[k][j][i]);
          QUIT_PLUTO(1);
        }
        continue; /* Skip this point since j = jp would be the same */
      }
    }else {             /* Vector is define at center or radial interfaces */     
      ip = i;
      jp = 2*JBEG-j-1;
    }
    scrh1 = fabs(Vx1[k][j][i] - Vx1[k][jp][ip]); /* r-component:   symmetric */
    scrh2 = fabs(Vx2[k][j][i] + Vx2[k][jp][ip]); /* th-component:  antisymmetric */
    scrh3 = fabs(Vx3[k][j][i] + Vx3[k][jp][ip]); /* phi-component: antisymmetric */
    #endif


    if (scrh1 > 1.e-14 || scrh2 > 1.e-14 || scrh3 > 1.e-14){
      #if GEOMETRY == CYLINDRICAL
      if (where == X1FACE){
        printLog ("! SymmetryCheck(): Vector not symmetric at i+1/2,j = %d+1/2, %d\n",i,j);
      }else if (where == X2FACE){
        printLog ("! SymmetryCheck(): Vector not symmetric at i,j+1/2 = %d, %d+1/2\n",i,j);
      }else {
        printLog ("! SymmetryCheck(): Vector not symmetric at i,j = %d, %d\n",i,j);
      }
      #elif GEOMETRY == SPHERICAL
      if (where == X1FACE){
        printLog ("! SymmetryCheck(): Vector not symmetric at i+1/2,j = %d+1/2, %d\n",i,j);
      }else if (where == X2FACE){
        printLog ("! SymmetryCheck(): Vector not symmetric at i,j+1/2 = %d, %d+1/2\n",i,j);
      }else{
        printLog ("! SymmetryCheck(): Vector not symmetric at i,j = %d, %d\n",i,j);
      }
      #endif
      printLog ("! i, j, k = %d, %d, %d --> V = (%12.6e, %12.6e, %12.6e) \n",
                   i,j,k,Vx1[k][j][i],Vx2[k][j][i], Vx3[k][j][i]);
      printLog ("! ip,jp,k = %d, %d, %d --> V = (%12.6e, %12.6e, %12.6e)\n",
                 ip,jp,k, Vx1[k][jp][ip],Vx2[k][jp][ip], Vx3[k][jp][ip]);

      QUIT_PLUTO(1);
    } 
  }
#endif /* GEOMETRY != CARTESIAN */
}

/* ********************************************************************* */
void SwapEndian (void *x, const int nbytes) 
/*!
 * Swap the byte order of x.
 * 
 * \param [in] x      pointer to the variable being swapped
 * \param [in] nbytes data type size
 * \return This function has no return value.
 *********************************************************************** */
{
  int k;
  static char Swapped[16];
  char *c;

  c = (char *) x;

  for (k = nbytes; k--; ) Swapped[k] = *(c + nbytes - 1 - k);
  for (k = nbytes; k--; ) c[k] = Swapped[k];
}

/* ********************************************************************* */
const char *VarName (int nv, int mode) 
/*!
 * Return the name of the variable with index nv.
 * Mode 1/2 can be used to choose prim / conservative notations.
 * 
 * \return A string
 *********************************************************************** */
{
  if (nv == RHO) return "RHO";
  if (mode == 1) {
    if (nv == VX1) return "VX1";
    if (nv == VX2) return "VX2";
    if (nv == VX3) return "VX3";
  }else{
    if (nv == MX1) return "MX1";
    if (nv == MX2) return "MX2";
    if (nv == MX3) return "MX3";
  }

  #if HAVE_ENERGY
  if (mode == 1) {
    if (nv == PRS) return "PRS";
  }else{
    if (nv == ENG) return "ENG";
  }
  #endif

  #if (PHYSICS == MHD) || (PHYSICS == RMHD) || (PHYSICS == ResRMHD)
  if (nv == BX1) return "BX1";
  if (nv == BX2) return "BX2";
  if (nv == BX3) return "BX3";
  #endif

  #if (PHYSICS == ResRMHD)
  if (nv == EX1) return "EX1";
  if (nv == EX2) return "EX2";
  if (nv == EX3) return "EX3";
  #endif
  
  #if NSCL > 0
  if (nv == TRC) return "TRC";
  #endif

  return "(null)";
}

/* ********************************************************************* */
void WriteAsciiFile (char *fname, double *q, int nvar)
/*! 
 *  Write one row a multi-column ascii file 
 *
 *********************************************************************** */
{
  int n;
  static char old_file[64] = "\0";
  FILE *fp;
  
/* --------------------------------------------------------
    If the old file name matches the new one, then open 
    in append mode. Otherwise, open in write mode.
   -------------------------------------------------------- */
   
  if (strcmp (fname,old_file) == 0){  
    fp = fopen(fname,"a");
  }else{
    fp = fopen(fname,"w");
    strcpy (old_file, fname);
  }
  
  for (n = 0; n < nvar; n++) fprintf (fp,"%12.6e  ",q[n]);
  fprintf (fp,"\n");
  fclose(fp);
}

