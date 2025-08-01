/* ///////////////////////////////////////////////////////////////////// */
/*!
 \file
 \brief Writer for particles binary data in .vtk format.
 
 VTK particle files are written using the standard legacy
 format as unstructured grid dataset.
 The file consists of binary data (both float and int)
 interlevead by ASCII headers.
 Only the following fields may be written:
 
 <coord>, <id>, <tinj>, <color>, <speed>

 If \c <n> is the number of particles, the structure of a
 VTK particle file is the following:

 
     PLUTO <xxx> VTK Data
     BINARY
     DATASET UNSTRUCTURED_GRID
     POINTS <n> float
       <x1,y1,z1> <x2,y2,z2> ...  <xn,yn,zn>
     CELL_TYPE <n>
       1 ... 1
     POINT_DATA <n>
     SCALARS Identity int 1
     LOOKUP_TABLE default
       <id1> <id2> ... <idn>
     VECTORS Velocity float 3
       <vx1, vy1, vz1> <vx2, vy2, vz2> ... <vxn, vyn, vzn>
     SCALARS Color float 1
     LOOKUP_TABLE defaults
       <color0> <color1> ... <colorn>
     SCALAR t_inj float 1
     LOOKUP_TABLE default
       <tinj1> <tinj2> ... <tinjn>

 
 \authors B. Vaidya (bvaidya@unito.it)\n
          A. Mignone (andrea.mignone@unito.it)\n
          D. Mukherjee
 
 \date    Aug 20, 2020
 */
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

/* ********************************************************************* */
void Particles_WriteVTK(Data* data, Grid* grid, Output *output, char filename[128])
 /*! 
 *  Write particle Data in Float binary format.
 *
 *  \param [in]  PheadRef    Pointer to the Head Node of Particle List.
 *  \param [in]  filename    File name of particle data: particles.nnnn.vtk,
 *                           where nnnn is the data file number.
 *
 *  NOTE: Can be visualized with VisIt, Paraview.
 *        Particles has to visualised using Pseudocolor plot 
 *        in Visit and Glyps in Paraview.
 *********************************************************************** */
{
  size_t s;
  int   i,j,k,np, dir;
  int   indx, indx1, indx2, indx3;
  int   *pids, *typarr;
  long int offset, tot_np, proc_npart[g_nprocs+1];
  char  vtkheader[512];
  float *pos, *vel, *color, *tinj;
  double gamma, u[3], rcoord[3];

#if PARTICLES == PARTICLES_KIN
  long out_particles = NX1*NX2*NX3*NMOMENTUM;
#else
  particleNode * CurNode;
  long out_particles = p_nparticles;
#endif
 
  pids   = ARRAY_1D(out_particles, int);
  typarr = ARRAY_1D(out_particles, int);
  pos    = ARRAY_1D(3*out_particles, float);
  vel    = ARRAY_1D(3*out_particles, float);
  color  = ARRAY_1D(out_particles,float);
  tinj   = ARRAY_1D(out_particles,float);
  
/* --------------------------------------------------------------
   0. Create data arrays to be written (coordinates, velocities, 
      etc...). Swap byte ordering to big endian if necessary.
   -------------------------------------------------------------- */   
   
#if PARTICLES != PARTICLES_KIN
  i = 0;
  PARTICLES_LOOP(CurNode, data->PHead){

  /* -- Compute four-velocity if necessary -- */

    for (dir = 0; dir < 3; dir++) u[dir] = CurNode->p.speed[dir];

    #if GEOMETRY == POLAR
    {
      double r   = CurNode->p.coord[IDIR];
      double phi = (DIMENSIONS >= 2 ? CurNode->p.coord[JDIR]:0.0);

      rcoord[IDIR] = r*cos(phi);
      rcoord[JDIR] = r*sin(phi);
      rcoord[KDIR] = CurNode->p.coord[KDIR];
    }
    #elif GEOMETRY == SPHERICAL
    {
      double r   = CurNode->p.coord[IDIR];
      double th  = (DIMENSIONS >= 2 ? CurNode->p.coord[JDIR]:0.0);
      double phi = (DIMENSIONS == 3 ? CurNode->p.coord[KDIR]:0.0);

      rcoord[IDIR] = r*sin(th)*cos(phi);
      rcoord[JDIR] = r*sin(th)*sin(phi);
      rcoord[KDIR] = r*cos(th);
    }
    #else
    rcoord[IDIR] = CurNode->p.coord[IDIR];  /* Standard Cartesian  */
    rcoord[JDIR] = CurNode->p.coord[JDIR];  /* coordinates         */
    rcoord[KDIR] = CurNode->p.coord[KDIR];
    #endif

    for (dir = 0; dir < 3; dir++) {
      pos[3*i + dir] = rcoord[dir];
      vel[3*i + dir] = u[dir];
    }
    #if PARTICLES == PARTICLES_LP
    color[i]  = CurNode->p.color[0];
    #else
    color[i]  = CurNode->p.color;
    #endif

    pids[i]   = CurNode->p.id;
    tinj[i]   = CurNode->p.tinj;
    typarr[i] = 1;
 
    if (IsLittleEndian()) {
      for (dir = 0; dir < 3; dir++) {
        SWAP_VAR(pos[3*i+ dir]);
        SWAP_VAR(vel[3*i+ dir]);
      }
      SWAP_VAR(color[i]);
      SWAP_VAR(pids[i]);
      SWAP_VAR(tinj[i]);
      SWAP_VAR(typarr[i]);
    }   
    i++;
  }
#else
  np = 0;
  DOM_LOOP(k,j,i){
      for(int l = 0; l < NMOMENTUM; ++l){
          u[0] = data->p_grid[l];
          u[1] = 0;
          u[2] = 0;
          //todo other coordinates?

#if GEOMETRY == POLAR
          double r   = grid->x[IDIR][i];
          double phi = (DIMENSIONS >= 2 ? grid->x[JDIR][j]:0.0);

          rcoord[IDIR] = r*cos(phi);
          rcoord[JDIR] = r*sin(phi);
          rcoord[KDIR] = grid->x[KDIR][k];
#elif GEOMETRY == SPHERICAL
          double r   = grid->x[IDIR][i];
          double th  = (DIMENSIONS >= 2 ? grid->x[JDIR][j]:0.0);
          double phi = (DIMENSIONS == 3 ? grid->x[KDIR][k]:0.0);

          rcoord[IDIR] = r*sin(th)*cos(phi);
          rcoord[JDIR] = r*sin(th)*sin(phi);
          rcoord[JDIR] = r*cos(th);
#else
          rcoord[IDIR] = grid->x[IDIR][i];  /* Standard Cartesian  */
          rcoord[JDIR] = grid->x[JDIR][j];  /* coordinates         */
          rcoord[KDIR] = grid->x[KDIR][k];
#endif

          for (dir = 0; dir < 3; dir++) {
            pos[3*np + dir] = rcoord[dir];
            vel[3*np + dir] = u[dir];
          }
          color[np] = 0.0;
          pids[np] = 0.0;
          tinj[np] = 0.0;
          typarr[np] = 1;
          if (IsLittleEndian()) {
            for (dir = 0; dir < 3; dir++) {
              SWAP_VAR(pos[3*i+ dir]);
              SWAP_VAR(vel[3*i+ dir]);
            }
            SWAP_VAR(color[i]);
            SWAP_VAR(pids[i]);
            SWAP_VAR(tinj[i]);
            SWAP_VAR(typarr[i]);
          }
          np++;
      }
  }
#endif
    
/* ------------------------------------------------------------------
   1. Compute the total number of particles and individual processor
      offsets required to write file.
      The offsets are computed by summing the number of particles
      of all  ranks < prank and stored into the header[1:] array. 
   ------------------------------------------------------------------ */   

#ifdef PARALLEL
  MPI_Allreduce(&p_nparticles, &tot_np, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allgather(&p_nparticles, 1, MPI_LONG, &(proc_npart[1]), 1,
               MPI_LONG, MPI_COMM_WORLD);

/* Compute individual processor offset (in particle units) */

  offset = 0;
  for(i = 0; i < prank; i++) offset += proc_npart[i+1];

#else
  offset = 0;
  tot_np = proc_npart[1] = p_nparticles;
#endif  
  
/* ----------------------------------------------------------------
   2. Proc. #0 write initial header using UNSTRUCTURED_GRID data
      and particles coordinates.
   ---------------------------------------------------------------- */

  sprintf(vtkheader,"# vtk DataFile Version 2.0\n");
  sprintf(vtkheader+strlen(vtkheader),"PLUTO %s VTK Data\n",PLUTO_VERSION);
  sprintf(vtkheader+strlen(vtkheader),"BINARY\n");
  sprintf(vtkheader+strlen(vtkheader),"DATASET UNSTRUCTURED_GRID\n");
  sprintf(vtkheader+strlen(vtkheader),"POINTS %ld float\n",tot_np);
  FileWriteHeader(vtkheader, filename, -1);  

  s = sizeof(float);
  FileWriteArray(pos, 3*offset*s, 3*p_nparticles, s, filename);
  
/* ----------------------------------------------------------
   3. Write particle CELL Type Header and data.
   ---------------------------------------------------------- */
   
  sprintf(vtkheader, "\nCELL_TYPES %ld\n",tot_np);
  FileWriteHeader(vtkheader, filename, 0);      

  s = sizeof(int);
  FileWriteArray(typarr, offset*s, p_nparticles, s, filename);

  sprintf(vtkheader, "\nPOINT_DATA %ld\n",tot_np); 
  FileWriteHeader(vtkheader, filename, 0);      

/* ----------------------------------------------------------
   4. Write particle ID header and data.
   ---------------------------------------------------------- */

  indx = StringArrayIndex(output->var_name, output->nvar, "id");
  if (output->dump_var[indx]){
    sprintf(vtkheader, "\nSCALARS Identity int 1\n"); 
    sprintf(vtkheader+strlen(vtkheader), "LOOKUP_TABLE default\n");    
    FileWriteHeader(vtkheader, filename, 0);  
  
    s = sizeof(int);
    FileWriteArray(pids, offset*s, p_nparticles, s, filename);
  }

/* ----------------------------------------------------------
   5. Write particle injection time header and data.
   ---------------------------------------------------------- */

  indx = StringArrayIndex(output->var_name, output->nvar, "tinj");
  if (output->dump_var[indx]){
    sprintf(vtkheader, "\nSCALARS tinj float 1\n"); 
    sprintf(vtkheader+strlen(vtkheader), "LOOKUP_TABLE default\n"); 
    FileWriteHeader(vtkheader, filename, 0);  

    s = sizeof(float);
    FileWriteArray(tinj, offset*s, p_nparticles, s, filename);
  }

/* ----------------------------------------------------------
   6. Write particle color header and data.
   ---------------------------------------------------------- */
   
  indx = StringArrayIndex(output->var_name, output->nvar, "color");
  if (output->dump_var[indx]){
    sprintf(vtkheader, "\nSCALARS Color float 1\n"); 
    sprintf(vtkheader+strlen(vtkheader), "LOOKUP_TABLE default\n"); 
    FileWriteHeader(vtkheader, filename, 0);  

    s = sizeof(float);
    FileWriteArray(color, offset*s, p_nparticles, s, filename);
  }

/* ----------------------------------------------------------
   7. Write particle velocity header and (vector) data
   ---------------------------------------------------------- */

  indx1 = StringArrayIndex(output->var_name, output->nvar, "vx1");
  indx2 = StringArrayIndex(output->var_name, output->nvar, "vx2");
  indx3 = StringArrayIndex(output->var_name, output->nvar, "vx3");
  if (   output->dump_var[indx1]
      && output->dump_var[indx2]
      && output->dump_var[indx3]){
  
    #if (PARTICLES == PARTICLES_CR) || (PARTICLES == PARTICLES_MC)
    sprintf(vtkheader, "\nVECTORS Four-Velocity float 3\n"); 
    FileWriteHeader(vtkheader, filename, 0);  
    #else  
    sprintf(vtkheader, "\nVECTORS Velocity float 3\n"); 
    FileWriteHeader(vtkheader, filename, 0);  
    #endif

    s = sizeof(float);
    FileWriteArray(vel, 3*offset*s, 3*p_nparticles, s,filename);
  }

/* ----------------------------------------------------------
   8. Free memory.
   ---------------------------------------------------------- */

  FreeArray1D(pids);
  FreeArray1D(typarr);
  FreeArray1D(pos);
  FreeArray1D(vel);
  FreeArray1D(color);
  FreeArray1D(tinj);
}
