/* ///////////////////////////////////////////////////////////////////// */
/*!
 \file
 \brief Set turbulence output data attributes.
 
 \authors A. Mignone (andrea.mignone@unito.it)\n
 
 \date    Aug 27, 2020
 */
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

#if TURBULENT_FIELD == YES

static Output *all_outputs;

/* ********************************************************************* */
void Turbulence_SetOutput (Data *d, Runtime *runtime)
/*!
 *  Set default attributes (variable names, pointers to data structures, 
 *  filename extensions, etc...) of the particle output structures.
 *
 * \param [in] d        pointer to Data structure
 * \param [in] runtime  pointer to Runtime structure
 *
 *********************************************************************** */
{
  int nv, i, k;
  char vname[64];
  Output *output;

  all_outputs = runtime->output;

/* --------------------------------------------------------
   1. Loop on output types 
   -------------------------------------------------------- */

  for (k = 0; k < MAX_OUTPUT_TYPES; k++){

    output = runtime->output + k;

    //for (i = 0; i < MAX_OUTPUT_VARS; i++) output->field_dim[i] = 1;

  /* -- 1a. Exclude fluid output types -- */

    if (! (output->type == TURBULENCE_DBL_OUTPUT ||
           output->type == TURBULENCE_FLT_OUTPUT ||
           output->type == TURBULENCE_VTK_OUTPUT ||
           output->type == TURBULENCE_TAB_OUTPUT ||
           output->type == TURBULENCE_HDF5_OUTPUT)) continue;

    for (i = 0; i < MAX_OUTPUT_VARS; i++) output->field_dim[i] = 1;

  /* ------------------------------------------------------
     1b. Allocate memory for field names.
         Note that, for LP particles with spectra we must
         use a large number in order to fit spectral
         information.
     ------------------------------------------------------ */


    output->var_name = ARRAY_2D(16, 32, char);

    strcpy(output->dir, runtime->output_dir); /* output directory is the same     */
                                              /* for all outputs (easy to change) */
    output->nfile    = -1;

  /* -- 1c. Set particles filename extensions -- */

    if (output->type == TURBULENCE_DBL_OUTPUT) strcpy (output->ext,"dbl");
    if (output->type == TURBULENCE_FLT_OUTPUT) strcpy (output->ext,"flt");
    if (output->type == TURBULENCE_VTK_OUTPUT) strcpy (output->ext,"vtk");
    if (output->type == TURBULENCE_TAB_OUTPUT) strcpy (output->ext,"tab");

  /* ------------------------------------------------------
     1d. Set default field names (all of them).
         Important: the order must be preserved when
         mapping strucure members to array!
     ------------------------------------------------------ */

    i = 0;

    strcpy(output->var_name[i++], "id");
    strcpy(output->var_name[i++], "x1");
    strcpy(output->var_name[i++], "x2");
    strcpy(output->var_name[i++], "x3");
    strcpy(output->var_name[i++], "i");
    strcpy(output->var_name[i++], "j");
    strcpy(output->var_name[i++], "k");
    strcpy(output->var_name[i++], "kturb");
    strcpy(output->var_name[i++], "W");
    strcpy(output->var_name[i++], "dV");
    output->field_dim[7] = NTURB;
    output->field_dim[8] = NTURB;

    output->nvar = i;
    
  /* -- Initialize dump_var to YES for all fields -- */  
    for (nv = output->nvar; nv--; ) output->dump_var[nv] = YES;
  }

}

#endif
