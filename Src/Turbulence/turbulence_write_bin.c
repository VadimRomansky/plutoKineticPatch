#include "turbulence.h"

/* ///////////////////////////////////////////////////////////////////// */
/*!
 \file
 \brief Hack for writing turbulence as particles binary data in .dbl, .flt and
        ASCII in .tab (only for serial version) format .

 \authors A. Mignone (andrea.mignone@unito.it)\n
          B. Vaidya (bvaidya@unito.it)\n
          D. Mukherjee

 \date    Oct 31, 2020
 */
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

#if PARTICLES != PARTICLES_LP
/* ********************************************************************* */
void Turbulence_WriteBinary(Data* data, Grid* grid, double dt_magnetic,
                           Output *output, char *filename)
/*!
 *  Write particle data in single or double precision binary format.
 *  The binary file structure is:
 *
 *  <Header (ASCII) section>
 *   .
 *   .
 *   .
 *  {field1, field2, ...}_p
 *   .
 *   .
 *   .
 *
 * All fields are mapped from the particle structure into an array
 * and thus converted to either single or double precision, depending
 * on the function that is being called.
 * Fields can be scalars or array: for each particle nelem indicates
 * the number of elements written for each particle.
 * The field order does not have to match the structure member order but it
 * must be consistent with the sequence declared in Particles_SetOutput()
 * function.
 * Individual fields may be excluded by calling SetOutputVar()
 * from ChangeOutputVar ().
 *
 *  \param [in]  PheadRef      Pointer to the Head Node of Particle List.
 *  \param [in]  dt_particles  Particle time step
 *  \param [in]  output        Pointer to output structure
 *  \param [in]  filename      File name of particle data: particles.nnnn.flt,
 *                             where nnnn is the data file number.
 ************************************************************************* */
{
  char     fheader[1024];
  size_t   size;
  int      dir, nv, nfields;
  int     *dump_var = output->dump_var;
  int      nvar     = output->nvar;
  long int nelem;
  long int i,j,k, offset, nparticles_glob, proc_npart[g_nprocs+1];
  void    *arr;
  float   *farr;
  double  *darr, gamma;

#if PARTICLES != PARTICLES_KIN
  particleNode * CurNode;
#endif

/* --------------------------------------------------------
   0. Allocate memory for required fields
   -------------------------------------------------------- */

  nfields = 0; /* Count how many fields are written to disk */
  nelem   = 0; /* The number of elements to be written (some fields may
                  be arrays). */
  for (nv = 0; nv < nvar; nv++) {
    if (dump_var[nv]) {
      nfields++;
      nelem += output->field_dim[nv];
    }
  }
#if PARTICLES == PARTICLES_KIN
  long out_particles = 0;
  DOM_LOOP(k,j,i){
          int ig = i + grid->beg[0] - grid->lbeg[0] - IBEG;
          int jg = j + grid->beg[1] - grid->lbeg[1] - JBEG;
          int kg = k + grid->beg[2] - grid->lbeg[2] - KBEG;
          int write_p = 1;
#if INCLUDE_IDIR
          if( ig % PARTICLES_KIN_OUTPUT_STEP != 0){
              write_p = 0;
          }
#endif

#if INCLUDE_JDIR
          if(jg % PARTICLES_KIN_OUTPUT_STEP != 0){
              write_p = 0;
          }
#endif
#if INCLUDE_KDIR
          if(kg % PARTICLES_KIN_OUTPUT_STEP != 0){
              write_p = 0;
          }
#endif
          if(write_p){
              out_particles = out_particles + 1;
          }
}
#ifdef PARALLEL
  //long nparl[1];
  //nparl[0] = out_particles;
  //long nparg[1];
  //MPI_Allreduce(nparl, nparg, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
  //out_particles = nparg[0];
#endif
#else
  long out_particles = p_nparticles;
#endif
  darr = ARRAY_1D(nelem*out_particles, double);
  farr = ARRAY_1D(nelem*out_particles, float);


/* --------------------------------------------------------
   1. Loop over particles and map struct. members to array
   -------------------------------------------------------- */

  i  = 0;  /* Array index */
  int iarr = 0;

  #if PARTICLES != PARTICLES_KIN
  PARTICLES_LOOP(CurNode, data->PHead){

  /* ------------------------------------------------------
     1b. Map structure members to array.
         Important: field order should match the order
         given in Particles_SetOutput().
         Here nv scan *all* fields (nv <= nfields)
     ------------------------------------------------------ */

    nv = 0;
    #if (PARTICLES == PARTICLES_CR) || (PARTICLES == PARTICLES_MC)
    if (dump_var[nv++]) darr[i++] = CurNode->p.id;
    if (dump_var[nv++]) darr[i++] = CurNode->p.coord[IDIR];
    if (dump_var[nv++]) darr[i++] = CurNode->p.coord[JDIR];
    if (dump_var[nv++]) darr[i++] = CurNode->p.coord[KDIR];
    if (dump_var[nv++]) darr[i++] = CurNode->p.speed[IDIR];
    if (dump_var[nv++]) darr[i++] = CurNode->p.speed[JDIR];
    if (dump_var[nv++]) darr[i++] = CurNode->p.speed[KDIR];
    if (dump_var[nv++]) darr[i++] = CurNode->p.mass;
    if (dump_var[nv++]) darr[i++] = CurNode->p.tinj;
    if (dump_var[nv++]) darr[i++] = CurNode->p.color;
    #endif

    #if PARTICLES == PARTICLES_DUST
    if (dump_var[nv++]) darr[i++] = CurNode->p.id;
    if (dump_var[nv++]) darr[i++] = CurNode->p.coord[IDIR];
    if (dump_var[nv++]) darr[i++] = CurNode->p.coord[JDIR];
    if (dump_var[nv++]) darr[i++] = CurNode->p.coord[KDIR];
    if (dump_var[nv++]) darr[i++] = CurNode->p.speed[IDIR];
    if (dump_var[nv++]) darr[i++] = CurNode->p.speed[JDIR];
    if (dump_var[nv++]) darr[i++] = CurNode->p.speed[KDIR];
    if (dump_var[nv++]) darr[i++] = CurNode->p.mass;
    if (dump_var[nv++]) darr[i++] = CurNode->p.tau_s;
    if (dump_var[nv++]) darr[i++] = CurNode->p.tinj;
    if (dump_var[nv++]) darr[i++] = CurNode->p.color;
    #endif

  } /* End PARTICLES_LOOP() */
#else
  DOM_LOOP(k,j,i){
          int ig = i + grid->beg[0] - grid->lbeg[0] - IBEG;
          int jg = j + grid->beg[1] - grid->lbeg[1] - JBEG;
          int kg = k + grid->beg[2] - grid->lbeg[2] - KBEG;
          int write_p = 1;
#if INCLUDE_IDIR
          if( ig % PARTICLES_KIN_OUTPUT_STEP != 0){
              write_p = 0;
          }
#endif

#if INCLUDE_JDIR
          if(jg % PARTICLES_KIN_OUTPUT_STEP != 0){
              write_p = 0;
          }
#endif

#if INCLUDE_KDIR
          if(kg % PARTICLES_KIN_OUTPUT_STEP != 0){
              write_p = 0;
          }
#endif
          if(write_p){
          nv = 0;
          if (dump_var[nv++]) darr[iarr++] = 0.0;
          if (dump_var[nv++]) darr[iarr++] = grid->x[0][i];
          if (dump_var[nv++]) darr[iarr++] = grid->x[1][j];
          if (dump_var[nv++]) darr[iarr++] = grid->x[2][k];
          if (dump_var[nv++]) darr[iarr++] = (i + grid->beg[0] - grid->lbeg[0] - IBEG);
          //printf("i = %d\n",i + grid->beg[0] - grid->lbeg[0] - IBEG);
          if (dump_var[nv++]) darr[iarr++] = (j + grid->beg[1] - grid->lbeg[1] - JBEG);
          if (dump_var[nv++]) darr[iarr++] = (k + grid->beg[2] - grid->lbeg[2] - KBEG);

          if (dump_var[nv++]) {
              for(int l = 0; l < NMOMENTUM; ++l){
                  darr[iarr++] = data->p_grid[l];
                  if(darr[iarr-1] != darr[iarr-1]){
                      printf("darr = NaN\n");
                  }
              }
          }
          if (dump_var[nv++]) {
              for(int l = 0; l < NMOMENTUM; ++l){
                  darr[iarr++] = data->Fkin[k][j][i][l]/(data->p_grid[l]*data->p_grid[l]*data->p_grid[l]);
                  if(darr[iarr-1] != darr[iarr-1]){
                      printf("darr = NaN\n");
                  }
              }
          }
          if (dump_var[nv++]) {
              darr[iarr++] = grid->dV[k][j][i];
              if(darr[iarr-1] != darr[iarr-1]){
                  printf("darr = NaN\n");
              }
          }
          if(dump_var[nv++]){
              darr[iarr++] = data->injectedEnergy[k][j][i];
              if(darr[iarr-1] != darr[iarr-1]){
                  printf("darr = NaN\n");
              }
          }
      }
  }
#endif

/* --------------------------------------------------------
   2. Compute the total number of particles and gather
      each proc. particle number into the proc_npart[1:]
      array. This serves to compute the file offset.
   -------------------------------------------------------- */

  offset = 0L;
#ifdef PARALLEL
  MPI_Allreduce(&out_particles, &nparticles_glob, 1,
                MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allgather(&out_particles, 1, MPI_LONG, &(proc_npart[1]), 1,
                MPI_LONG, MPI_COMM_WORLD);

  /* Compute individual processor offset (in particle units) */

  for(i = 0; i < prank; i++) offset += proc_npart[i+1];

  MPI_File fhw;
  MPI_File_open(MPI_COMM_WORLD, filename,
                MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fhw);
  MPI_File_close(&fhw);
#else
  nparticles_glob = out_particles;
#endif

/* --------------------------------------------------------
   3. Write file header section.
   -------------------------------------------------------- */

  sprintf(fheader,"# PLUTO %s binary particle data file\n", PLUTO_VERSION);
  sprintf(fheader+strlen(fheader),"# dimensions     %d\n",DIMENSIONS);
  sprintf(fheader+strlen(fheader),"# dt_particles   %12.6e\n",dt_magnetic);
  if (IsLittleEndian())
      sprintf(fheader+strlen(fheader),"# endianity      little\n");
  else
      sprintf(fheader+strlen(fheader),"# endianity      big\n");

  sprintf(fheader+strlen(fheader),"# nparticles     %ld\n",nparticles_glob);
  sprintf(fheader+strlen(fheader),"# idCounter      %ld\n",p_idCounter);
  sprintf(fheader+strlen(fheader),"# particletype   %d\n",PARTICLES);
  if (output->type == PARTICLES_FLT_OUTPUT){
    sprintf(fheader+strlen(fheader),"# precision      float\n");
  }else if (output->type == PARTICLES_DBL_OUTPUT){
    sprintf(fheader+strlen(fheader),"# precision      double\n");
  }

  sprintf(fheader+strlen(fheader),"# time           %12.6e\n",g_time);
  sprintf(fheader+strlen(fheader),"# stepNumber     %ld\n",g_stepNumber);
  sprintf(fheader+strlen(fheader),"# nfields        %d\n",nfields);
  sprintf(fheader+strlen(fheader),"# field_names    ");
  for (i = 0; i < nvar; i++){
    if (dump_var[i]) {
      sprintf(fheader+strlen(fheader),"%s", output->var_name[i]);
      sprintf(fheader+strlen(fheader),"  ");
    }
  }
  sprintf(fheader+strlen(fheader),"\n");

  sprintf(fheader+strlen(fheader),"# field_dim      ");
  for (i = 0; i < nvar; i++){
    if (dump_var[i]) {
      sprintf(fheader+strlen(fheader),"%d",output->field_dim[i]);
      sprintf(fheader+strlen(fheader),"  ");
    }
  }
  sprintf(fheader+strlen(fheader),"\n");

  FileWriteHeader(fheader, filename, -1);

/* --------------------------------------------------------
   4. Write data
   -------------------------------------------------------- */

  if (output->type == PARTICLES_DBL_OUTPUT) size = sizeof(double);
  if (output->type == PARTICLES_FLT_OUTPUT) size = sizeof(float);

  offset *= nelem*size; /* Convert offset to bytes */
  if (output->type == PARTICLES_FLT_OUTPUT){
    for (i = 0; i < nelem*out_particles; i++) farr[i] = (float)darr[i];
    arr = (void *) farr;
  }else if (output->type == PARTICLES_DBL_OUTPUT){
    arr = (void *) darr;
  }
/*
if (output->type == PARTICLES_DBL_OUTPUT){
  i = 0;
  for (k = 0; k < p_nparticles; k++){
    for (nv = 0; nv < nfields; nv++){
      if (darr[i] != darr[i]){
        printLog ("nan found\n");
        QUIT_PLUTO(1);
      }
      printLog ("np = %d, darr = %f\n",k,darr[i]);
      i++;
    }
  }
}
*/
  FileWriteArray(arr, offset, nelem*out_particles, size, filename);

/* --------------------------------------------------------
   5. Write data
   -------------------------------------------------------- */

  FreeArray1D(darr);
  FreeArray1D(farr);
}
#endif /* PARTICLES != PARTICLES_LP */

/* ********************************************************************* */
void Particles_WriteTab(Data* data, Grid* grid, char filename[128])
/*
 * Write particle coordinates, ids and speeds in Tab ASCII format
 * only for *serial* version.
 *
 *  \param [in]  PheadRef    Pointer to the Head Node of Particle List.
 *  \param [in]  filename    File name of particle data: particles.nnnn.tab,
 *                           where nnnn is the data file number.
 *********************************************************************** */
{
#ifdef PARALLEL
  print("! WARNING: Particle Data in Tabulated ASCII format is only written \
            with serial version");
#else
  int i,j,k;
  FILE *stream;

#if PARTICLES == PARTICLES_KIN
  long out_particles = NX1*NX2*NX3*NMOMENTUM;
#else
  long out_particles = p_nparticles;
#endif

  stream = fopen(filename, "w");
  fprintf(stream, "# Nparticles: %ld\n", out_particles);
  fprintf(stream, "# Step:       %ld\n", g_stepNumber);
  fprintf(stream, "# Time:       %f\n",  g_time);

#if PARTICLES != PARTICLES_KIN
  particleNode* CurNode;

  CurNode = data->PHead;
  while(CurNode != NULL) {
    fprintf(stream, "%d", CurNode->p.id);

    for (i = 0; i < 3; ++i) {
      fprintf(stream, "  %lf", CurNode->p.coord[i]);
    }
    for (i = 0; i < 3; ++i) {
      fprintf(stream, "  %lf", CurNode->p.speed[i]);
    }
    fprintf(stream, "\n");

    CurNode = CurNode->next;
  }
#else
  DOM_LOOP(k,j,i){
      for(int l = 0; l < NMOMENTUM; ++l){
          fprintf(stream, "%d", 0);


          fprintf(stream, "  %lf %lf %lf", grid->x[0][i], grid->x[1][j], grid->x[2][k]);

          fprintf(stream, "  %d %d %d", i + grid->beg[0] - grid->lbeg[0], j + grid->beg[1] - grid->lbeg[1], k + grid->beg[2] - grid->lbeg[2]);


          fprintf(stream, "  %lf %lf %lf", data->p_grid[l], 0.0, 0.0);

          fprintf(stream, " %lf %lf %lf", data->Fkin[k][j][i], grid->dV[k][j][i], 0.0);

          fprintf(stream, "\n");
      }
  }
#endif
  fclose(stream);

#endif
}

