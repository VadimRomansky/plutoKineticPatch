#ifndef THREEDIAGONAL_H
#define THREEDIAGONAL_H

#include <stdio.h>
#include <math.h>
#ifdef PARALLEL
#include <mpi.h>
#endif

void sequentialThreeDiagonalSolverX(double**** x, double**** rightPart, double**** a, double**** b, double**** c, int Nx, int Ny, int Nz, int Nmomentum);

void sequentialThreeDiagonalSolverY(double**** x, double**** rightPart, double**** a, double**** b, double**** c, int Nx, int Ny, int Nz, int Nmomentum);

void sequentialThreeDiagonalSolverZ(double**** x, double**** rightPart, double**** a, double**** b, double**** c, int Nx, int Ny, int Nz, int Nmomentum);

void sequentialThreeDiagonalSolverP(double**** x, double**** rightPart, double**** a, double**** b, double**** c, int Nx, int Ny, int Nz, int Nmomentum);

void noparallelThreeDiagonalSolverX(double**** x, double**** rightPart, double**** a, double**** b, double**** c, int Nmomentum);

void noparallelThreeDiagonalSolverY(double**** x, double**** rightPart, double**** a, double**** b, double**** c, int Nmomentum);

void noparallelThreeDiagonalSolverZ(double**** x, double**** rightPart, double**** a, double**** b, double**** c, int Nmomentum);

void noparallelThreeDiagonalSolverP(double**** x, double**** rightPart, double**** a, double**** b, double**** c, int Nmomentum);

#ifdef PARALLEL
void parallelThreeDiagonalSolverX(double**** x, double**** rightPart, double**** a, double**** b, double**** c, int Nmomentum, int Nprocs, int rank, MPI_Comm comm);

void parallelThreeDiagonalSolverY(double**** x, double**** rightPart, double**** a, double**** b, double**** c, int Nmomentum, int Nprocs, int rank, MPI_Comm comm);

void parallelThreeDiagonalSolverZ(double**** x, double**** rightPart, double**** a, double**** b, double**** c, int Nmomentum, int Nprocs, int rank, MPI_Comm comm);
#endif

#endif // THREEDIAGONAL_H
