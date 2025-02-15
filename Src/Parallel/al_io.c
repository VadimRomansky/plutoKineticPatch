/* ///////////////////////////////////////////////////////////////////// */
/*!
  \file
  \brief Miscellaneous functions for IO operations
  
  \author A. Malagoli (University of Chicago)
  \author A. Mignone (mignone@ph.unito.it)
  \author G. Muscianisi (g.muscianisi@cineca.it)

  \date Aug 26, 2012
*/ 
/* ///////////////////////////////////////////////////////////////////// */
#include "al_hidden.h"  /*I "al_hidden.h" I*/

#include <string.h>

#include "io.h"

extern long int g_usedMemory;
extern char ***Array3D (int nx, int ny, int nz, size_t dsize);

/* 
   The SZ structure stack is defined and maintained
   in al_szptr_.c
   Here we include an external reference to it in
   order to be able to make internal references to it.
*/
extern SZ *sz_stack[AL_MAX_ARRAYS];
extern int stack_ptr[AL_MAX_ARRAYS];

/* ********************************************************************* */
int AL_File_open(char *filename, int sz_ptr )
/*!
 * Open a file associated with a distributed array
 *
 * \param [in] filename    name of the file
 * \param [in] sz_ptr      integer pointer to the distributed array descriptor
 *
 * \return  ifp   pointer to integer file pointer
 *********************************************************************** */
{
    int myrank, nproc;
    int errcode;
    MPI_File ifp;
    MPI_Comm comm;
    SZ *s;

    /* DIAGNOSTICS
    Check that sz_ptr points to an allocated SZ
  */
    if( stack_ptr[sz_ptr] == AL_STACK_FREE){
        printf("AL_Decompose: wrong SZ pointer\n");
    }

    s = sz_stack[sz_ptr];

    myrank = s->rank;
    nproc = s->size;

    comm = s->comm;

    MPI_Barrier(comm);


    errcode = MPI_File_open(comm, filename,
                            MPI_MODE_CREATE | MPI_MODE_RDWR | MPI_MODE_UNIQUE_OPEN,
                            MPI_INFO_NULL, &ifp);

    s->io_offset = 0;
    s->ifp = ifp;

    /* DIAGNOSTICS */
#ifdef DEBUG
    int myid, len;
    char es[128];
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    if( errcode ){
        MPI_Error_string(errcode, es, &len);
        printf("Errcode from MPI_File_open: %d | %s\n", errcode, es);
    }
    printf("myid %d, AL_File_open: Opened file %s\n",myid,filename);
#endif

    return (int) AL_SUCCESS;
}


/* ********************************************************************* */
int AL_File_close(int sz_ptr)
/*!
 * Close a file associate with and array
 *
 * \param [in] sz_ptr   integer pointer to array descriptor
 *********************************************************************** */
{
    int myrank;
    int errcode;
    SZ *s;

    s = sz_stack[sz_ptr];
    myrank = s->rank;

    errcode =  MPI_File_close(&(s->ifp));

    /* DIAGNOSTICS */
#ifdef DEBUG
    int myid, len;
    char es[128];
    if( errcode ){
        MPI_Comm_rank(MPI_COMM_WORLD, &myid);
        MPI_Error_string(errcode, es, &len);
        printf("Errcode from MPI_File_close: %d | %s\n", errcode, es);
    }
    printf("myid %d, AL_File_close: Closed file, errcode %d \n",myid, errcode);
#endif

    return (int) AL_SUCCESS;
}



#ifdef NSAMPLING

double*** AL_Create_sampled_array(double*** V, int sz_ptr)
/*!
 * Creates sampled to parallel IO
 *
 * \param [in] V  poinoutput sampled array
 * \param [in] sz_ptr  integer pointer to the distributed array descriptor
 *********************************************************************** */
{
    SZ *s;
    s = sz_stack[sz_ptr];
    int myrank = s->rank;
    size_t dsize = sizeof(double);

    int koffset = s->samplinggstart[2]%NSAMPLING;
    if(koffset > 0){
        koffset = NSAMPLING - koffset;
    }
    int joffset = s->samplinggstart[1]%NSAMPLING;
    if(joffset > 0){
        joffset = NSAMPLING - joffset;
    }
    int ioffset = s->samplinggstart[0]%NSAMPLING;
    if(ioffset > 0){
        ioffset = NSAMPLING - ioffset;
    }


    double*** a = (double***) Array3D(s->samplinglarrdim[2], s->samplinglarrdim[1], s->samplinglarrdim[0], dsize);
    for(int k =0; k < s->samplinglarrdim[2]; ++k){
        for(int j = 0; j < s->samplinglarrdim[1]; ++j){
            for(int i = 0; i < s->samplinglarrdim[0]; ++i){
                //printf("rank = %d, k = %d, j = %d, i = %d\n", myrank, k, j, i);
                //printf("rank = %d, koffset = %d, joffset = %d, ioffset = %d\n", myrank, koffset, joffset, ioffset);
                //printf("V = %g\n", V[koffset + k*NSAMPLING][joffset + j*NSAMPLING][ioffset + i*NSAMPLING]);
                a[k][j][i] = V[s->lbeg[2] + koffset + k*NSAMPLING][s->lbeg[1] + joffset + j*NSAMPLING][s->lbeg[0] + ioffset + i*NSAMPLING];
            }
        }
    }
    return a;
}

double*** AL_Get_sampled_array(double*** V, int sz_ptr)
/*!
 * Creates sampled to parallel IO
 *
 * \param [in] V  poinoutput sampled array
 * \param [in] sz_ptr  integer pointer to the distributed array descriptor
 *********************************************************************** */
{
    SZ *s;
    s = sz_stack[sz_ptr];
    int myrank = s->rank;
    size_t dsize = sizeof(double);

    int koffset = s->samplinggstart[2]%NSAMPLING;
    if(koffset > 0){
        koffset = NSAMPLING - koffset;
    }
    int joffset = s->samplinggstart[1]%NSAMPLING;
    if(joffset > 0){
        joffset = NSAMPLING - joffset;
    }
    int ioffset = s->samplinggstart[0]%NSAMPLING;
    if(ioffset > 0){
        ioffset = NSAMPLING - ioffset;
    }
    for(int k =0; k < s->samplinglarrdim[2]; ++k){
        for(int j = 0; j < s->samplinglarrdim[1]; ++j){
            for(int i = 0; i < s->samplinglarrdim[0]; ++i){
                //printf("rank = %d, k = %d, j = %d, i = %d\n", myrank, k, j, i);
                //printf("rank = %d, koffset = %d, joffset = %d, ioffset = %d\n", myrank, koffset, joffset, ioffset);
                //printf("V = %g\n", V[koffset + k*NSAMPLING][joffset + j*NSAMPLING][ioffset + i*NSAMPLING]);
                s->sampledArray[k][j][i] = V[s->lbeg[2] + koffset + k*NSAMPLING][s->lbeg[1] + joffset + j*NSAMPLING][s->lbeg[0] + ioffset + i*NSAMPLING];
            }
        }
    }
    return s->sampledArray;
}

#endif

/* ********************************************************************* */
int AL_Write_array(void *va, int sz_ptr, int istag)
/*!
 * Write a distributed array to a file in parallel using synchronous and
 * collective IO operations
 *
 * \param [in] buffer  pointer to the buffer to write
 * \param [in] sz_ptr  integer pointer to the distributed array descriptor
 * \param [in] istag   set it to -1 for centred variables,
 *                     set to 0,1,2 for staggered field in
 *                     the x,y,z direction
 *********************************************************************** */
{
    char *a;
    register int i;
    int myrank;
    long long nelem;
    int errcode;
    SZ *s;

    MPI_Status status;
    AL_Datatype gsub_arr, lsub_arr;

    MPI_File ifp;
    MPI_Offset offset;
    int  size;


    a = (char *) va;

    s = sz_stack[sz_ptr];
    myrank = s->rank;
    ifp = s->ifp;

    offset   = s->io_offset;
    //printf("offset = %d, rank = %d\n", offset, myrank);
#ifdef NSAMPLING
    if (istag == -1){
        gsub_arr = s->gsamplingsubarr;
        lsub_arr = s->lsamplingsubarr;
    }else{
        //todo
        printf("istg\n");
        gsub_arr = s->gsubarr_stag[istag];
        lsub_arr = s->lsubarr_stag[istag];
    }

#else

    if (istag == -1){
        gsub_arr = s->gsubarr;
        lsub_arr = s->lsubarr;
    }else{
        gsub_arr = s->gsubarr_stag[istag];
        lsub_arr = s->lsubarr_stag[istag];
    }

#endif

    MPI_Barrier(s->comm);

    //printf("rank = %d, offset = %d\n", myrank, offset);
    errcode = MPI_File_set_view(ifp, offset, MPI_BYTE, gsub_arr,
                                "native", MPI_INFO_NULL);
    //printf("rank = %d, offset = %d\n", myrank, offset);
#ifdef DEBUG
    int myid, len;
    char es[256];
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    if( errcode ){
        MPI_Error_string(errcode, es, &len);
        printf("Errcode from MPI_File_set_view: %d | %s\n", errcode, es);
    }
    printf("myid %d, offset file_set_view %lld\n", myid,offset);
#endif

    errcode = MPI_File_write_all(ifp, a, 1, lsub_arr, &status);

#ifdef DEBUG
    if( errcode ){
        MPI_Error_string(errcode, es, &len);
        printf("Errcode from MPI_File_write_all: %d | %s\n", errcode, es);
    }
#endif
    
    MPI_Type_size( s->type, &size);

    nelem = 1;
    for(i = 0; i < s->ndim; i++) {
        /* -------------------------------------------------------------------- */
        /*!  Bugs fixed on Aug 26, 2012:
  \li nelem has been declared long long in order to manage dbl and flt
     output 'single_file' in which each PLUTO variable is >= 4GB
  \li the updating of nelem after the writing of a variable has been
     modified from  ''nelem *= (long long)(s->arrdim[i] +(istag == i));''
     to ''nelem *= (long long)(s->arrdim[i]);'' because when a staggered
     variable is written, in the sz_ptr descriptor there is the
     information about the really number of element of this variable.
     This appens because the function AL_Set_staggered_dim is called
     for a staggered array. */
        /* -------------------------------------------------------------------- */
        /*    nelem *= (long long)(s->arrdim[i] +(istag == i)); */

#ifdef NSAMPLING
        nelem *= (long long)(s->samplingarrdim[i]);
#else
        nelem *= (long long)(s->arrdim[i]);
#endif
    }

    s->io_offset += (long long)(size)*nelem;


    return (int) AL_SUCCESS;
}


/* ********************************************************************* */
int AL_Read_array(void *va, int sz_ptr, int istag)
/*!
 * Read a distributed array to a file in parallel using synchronous and
 * collective IO operations
 *
 * \param [in] buffer  pointer to the buffer to write
 * \param [in] sz_ptr  integer pointer to a distributed array descriptor
 * \param [in] istag   set it to -1 for centred variables,
 *                     set to 0,1,2 for staggered field in
 *                     the x,y,z direction
 *********************************************************************** */
{
    char *a;
    register int i;
    int myrank;
    long long nelem;
    SZ *s;
    MPI_Status status;
    AL_Datatype gsub_arr, lsub_arr;


    MPI_File ifp;
    MPI_Offset offset;
    int size;

    a = (void *)va;

    s = sz_stack[sz_ptr];
    myrank = s->rank;

    ifp = s->ifp;

    offset   = s->io_offset;
    if (istag == -1){
        gsub_arr = s->gsubarr;
        lsub_arr = s->lsubarr;
    }else{
        gsub_arr = s->gsubarr_stag[istag];
        lsub_arr = s->lsubarr_stag[istag];
    }

    MPI_Barrier(s->comm);

    MPI_File_set_view(ifp, offset, MPI_BYTE, gsub_arr,
                      "native", MPI_INFO_NULL);
    MPI_File_read_all(ifp, a, 1, lsub_arr, &status);

    MPI_Type_size( s->type, &size);

    nelem = 1;
    for(i = 0; i < s->ndim; i++) {
        /* -------------------------------------------------------------------- */
        /*!  Bugs fixed on Aug 26, 2012:
  \li nelem has been declared long long in order to manage dbl and flt
     output 'single_file' in which each PLUTO variable is >= 4GB
  \li the updating of nelem after the writing of a variable has been
     modified from  ''nelem *= (long long)(s->arrdim[i] +(istag == i));''
     to ''nelem *= (long long)(s->arrdim[i]);'' because when a staggered
     variable is written, in the sz_ptr descriptor there is the
     information about the really number of element of this variable.
     This appens because the function AL_Set_staggered_dim is called
     for a staggered array. */
        /* -------------------------------------------------------------------- */
        /*    nelem *= (long long)(s->arrdim[i] + (istag == i)); */
        nelem *= (long long)(s->arrdim[i]);
    }
    s->io_offset += (long long)(size)*nelem;

    return (int) AL_SUCCESS;
}

/* -------------------------------------------------------------------
    New ArrayLib additions:
   ------------------------------------------------------------------- */

/* ******************************************************************* */
int AL_Write_header(void *vbuffer, int nelem, AL_Datatype type, int sz_ptr)
/*
 *
 *
 * Replace AL_Write_common which does not work.
 *
 *
 ********************************************************************* */
{
    char *buffer;
    int myrank;
    int nbytes, size;
    SZ *s;
    MPI_Status status;
    MPI_File ifp;


    MPI_Offset offset;
    MPI_Type_size(type, &size);
    buffer = (char *) vbuffer;

    nbytes = nelem*size;
    s      = sz_stack[sz_ptr];
    myrank = s->rank;

    ifp    = s->ifp;
    offset = s->io_offset;
    MPI_Barrier(s->comm);
    MPI_File_set_view(ifp, offset, MPI_BYTE, MPI_CHAR,
                      "native", MPI_INFO_NULL);
    if( myrank == 0 ){
        /*    MPI_File_write(ifp, buffer, strlen(buffer), MPI_CHAR, &status);  */
        MPI_File_write(ifp, buffer, nelem, type, &status);
    }

    s->io_offset += nbytes;


    return (int) AL_SUCCESS;
}

/* ************************************************************ */
long long AL_Get_offset(int sz_ptr)
/*
 *
 *
 ************************************************************** */
{
    SZ *s;

    if( stack_ptr[sz_ptr] == AL_STACK_FREE){
        printf("AL_Get_offset: wrong SZ pointer\n");
    }
    s  = sz_stack[sz_ptr];
    return s->io_offset;
}
/* ************************************************************ */
int AL_Set_offset(int sz_ptr, long long offset)
/*
 *
 *
 ************************************************************** */
{
    SZ *s;

    if( stack_ptr[sz_ptr] == AL_STACK_FREE){
        printf("AL_Get_offset: wrong SZ pointer\n");
    }
    s  = sz_stack[sz_ptr];
    s->io_offset = offset;
    return (int)AL_SUCCESS;
}
