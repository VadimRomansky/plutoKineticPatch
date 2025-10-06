#include "turbulence.h"

/* ///////////////////////////////////////////////////////////////////// */
/*!
 \file
 \brief Main driver of turbulence output data called during the CheckForOutput
        stage.

 \authors B. Vaidya (bvaidya@unito.it)\n
          A. Mignone (andrea.mignone@unito.it)\n

 \date   Aug 20, 2020
 */
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"
#if TURBULENT_FIELD == YES
/* ************************************************************************* */
void Turbulence_WriteData(Data *d, Output *output, Grid *grid)
/*!
 *  Main driver of particles output data.
 *
 *  \param[in]    d         Pointer to the PLUTO data structure.
 *  \param[in]    output    Pointer to the PLUTO Output strcuture.
 *
 **************************************************************************** */
{
  char filename[512];
  time_t tbeg, tend;




  output->nfile++;  /* Increment file output number */

  sprintf (filename,"%s/turbulence.%04d.%s", output->dir, output->nfile,
                                            output->ext);
  print ("> Writing turbulence file #%d (%s) to disk...\n",
             output->nfile, output->ext);

#ifdef PARALLEL
  MPI_Barrier (MPI_COMM_WORLD);
  if (prank == 0) time(&tbeg);
#endif

  if (   output->type == PARTICLES_DBL_OUTPUT
      || output->type == PARTICLES_FLT_OUTPUT) {

    Turbulence_WriteBinary(d, grid, 1.0/d->Dts->invDt_magnetic,
                          output, filename);

  }else if (output->type == PARTICLES_VTK_OUTPUT) {

    Turbulence_WriteVTK(d, grid, output, filename);

  }else if (output->type == PARTICLES_TAB_OUTPUT) {

    Turbulence_WriteTab(d, grid, filename);

  }

#ifdef PARALLEL
  MPI_Barrier (MPI_COMM_WORLD);
  if (prank == 0){
    time(&tend);
    print ("  [%5.2f sec ]\n",difftime(tend,tbeg));
  }
#endif

}

#endif
