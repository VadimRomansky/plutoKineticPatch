/* ///////////////////////////////////////////////////////////////////// */
/*!
  \file
  \brief Particle module header file

  Contains basic macro definitions, structure definitions and global
  variable declarations used by the Particle module.

  \author A. Mignone (andrea.mignone@unito.it)\n
          B. Vaidya
	        D. Mukherjee

  \date   Sep 05, 2022
*/
/* ///////////////////////////////////////////////////////////////////// */

#define INTEGER       1  /* Used for particle deposition */
#define REAL          2

#define PARTICLES_CREATE           1
#define PARTICLES_CREATE_WITH_ID   2
#define PARTICLES_TRANSFER         3
#define PARTICLES_RESTART          4

#ifndef PARTICLES_LP_SPECTRA
  #define PARTICLES_LP_SPECTRA   NO
#endif

/* *********************************************************************
    Define shapes in different directions

    Weight functions: 1 = Nearestg Grid Point (NGP)
                      2 = Cloud-In-Cell (CIC)
                      3 = Triangular Shape Cloud (TSC)
                      4 = Volume coordinate
   ********************************************************************* */

#ifndef PARTICLES_SHAPE
  #define PARTICLES_SHAPE  3
#endif

#ifndef PARTICLES_ISHAPE
  #define PARTICLES_ISHAPE  PARTICLES_SHAPE
#endif
#ifndef PARTICLES_JSHAPE
  #define PARTICLES_JSHAPE  PARTICLES_SHAPE
#endif
#ifndef PARTICLES_KSHAPE
  #define PARTICLES_KSHAPE  PARTICLES_SHAPE
#endif

/* *********************************************************************
    Test/Debug directive
   ********************************************************************* */

#ifndef PARTICLES_DEPOSIT
  #define PARTICLES_DEPOSIT   REAL
#endif

/* *********************************************************************
    MPI Data type
    Set PARTICLES_USE_MPI_DATATYPE to YES to define an appropriate
    datatype for MPI send/recv of structure.
    Set PARTICLES_USE_MPI_DATATYPE to NO to pass the structure as a
    sequence of BYTEs (requires no specific datatype)
   ********************************************************************* */

#ifndef PARTICLES_USE_MPI_DATATYPE
  #if PARTICLES == PARTICLES_LP
    #define PARTICLES_USE_MPI_DATATYPE  YES
  #else
    #define PARTICLES_USE_MPI_DATATYPE  NO
  #endif
#endif

/* *********************************************************************
    Particle type: COSMIC RAYS
   ********************************************************************* */

#if PARTICLES == PARTICLES_CR

/*! Set Guiding center approx to NO by default. */
  #ifndef PARTICLES_CR_GC
    #define PARTICLES_CR_GC  NO
  #endif

/*! Guiding Center flags. */
  #if PARTICLES_CR_GC == YES

    #ifndef  PARTICLES_CR_GC_TIME_STEPPING
      #define  PARTICLES_CR_GC_TIME_STEPPING   RK_MIDPOINT
    #endif
    #ifndef  PARTICLES_CR_GC_DEBUG
      #define PARTICLES_CR_GC_DEBUG   NO
    #endif

    #ifndef PARTICLES_CR_GC_TIME_DER
     #define PARTICLES_CR_GC_TIME_DER  NO
    #endif

    #ifndef PARTICLES_CR_GC_FULL_OMEGA
     #define PARTICLES_CR_GC_FULL_OMEGA   YES   /* Use the full expression  */
                                               /* for gyrofrequency */
    #endif

  /* Error codes */
    #define GC_ERR_INVALID          1   /* typically arises when GCA fails, e.g. */
                                        /*  sqrt(<0)                             */
    #define GC_ERR_ZERO_B           2   /* Magnetic null */
    #define GC_ERR_DOMAIN_OVERSTEP  3   /* particle outside allowed domain */
    #define GC_ERR_dRdt             4   /* Right hand side becomes superluminal */
    #define GC_ERR_LARMOR_OVERSTEP  5   /* Larmor radius exceeds dx */
    #define GC_ERR_FAST_VARYING     6

  #endif /* PARTICLES_CR_GC == YES */

/*! Sets the artificial value of the speed of light */
  #ifndef PARTICLES_CR_C
    #if PHYSICS == RMHD
      #define PARTICLES_CR_C       1.0
    #else
      #define PARTICLES_CR_C       10000.0
    #endif
  #endif

/*! Check and print the Larmor radius (debug only) */
  #ifndef PARTICLES_CR_CHECK_RL
    #define PARTICLES_CR_CHECK_RL    NO
  #endif

/*! Charge to mass ratio for CR particles [CR particles]. */
  #ifndef PARTICLES_CR_E_MC
    #define PARTICLES_CR_E_MC      (CONST_e/(CONST_mp*CONST_c))*UNIT_LENGTH*sqrt(4*CONST_PI*UNIT_DENSITY)
  #endif

/*! Charge to mass ratio for the fluid [CR Particles]. */
  #ifndef PARTICLES_CR_E_MC_GAS
    #define PARTICLES_CR_E_MC_GAS  (CONST_e/(CONST_mp*CONST_c))*UNIT_LENGTH*sqrt(4*CONST_PI*UNIT_DENSITY)
  #endif

/*! Enable / disable CR particle feedback onto the gas [CR particles]. */
  #ifndef PARTICLES_CR_FEEDBACK
    #define PARTICLES_CR_FEEDBACK  YES
  #endif

  #if PARTICLES_CR_FEEDBACK == NO
    #undef PARTICLES_CR_UPWIND_FLUX
    #define PARTICLES_CR_UPWIND_FLUX  NO

  /* Define the electric field indices for the test-particle approach
     when the PHYSICS module is different from ResRMHD
     These indices are directly taken from the Resistive RMHD module,
     where the electric field is an independent variable. */
    #if !(defined EX1) && (PARTICLES_SHAPE < 0)
      #define  EX1  8
      #define  EX2  9
      #define  EX3 10
    #endif
  #endif

/*! Fraction of Larmor time (1/OmegaL) covered during one hydro time step
    [CR particles] */
  #ifndef PARTICLES_CR_LARMOR_EPS
    #define PARTICLES_CR_LARMOR_EPS  0.3
  #endif

/*! Maximum number of zones crossed during one sub-cycle [CR particles].  */
  #ifndef PARTICLES_CR_NCELL_MAX
    #define PARTICLES_CR_NCELL_MAX  1.8
  #endif

/*! Max number of particles time steps per hydro step. When set to a
    negative integer, force sub-stepping to be  the same [CR particles]. */
  #ifndef PARTICLES_CR_NSUB
    #if PARTICLES_CR_GC == YES
      #define PARTICLES_CR_NSUB      1
    #else
      #define PARTICLES_CR_NSUB      4
    #endif
  #endif

/*! Predictor step in mover. This makes the scheme 2nd order [CR particles]. */
  #ifndef PARTICLES_CR_PREDICTOR
    #define PARTICLES_CR_PREDICTOR  2
  #endif

  #ifndef PARTICLES_CR_UPWIND_FLUX
    #define PARTICLES_CR_UPWIND_FLUX  NO
  #endif

/*! Default size of particle binary files */
  #ifndef PARTICLES_CR_MAX_FILE_SIZE
		#define PARTICLES_CR_MAX_FILE_SIZE 1.9 /* In GB */
  #endif

#endif

/* *********************************************************************
    Particle type: DUST
   ********************************************************************* */

#if PARTICLES == PARTICLES_DUST

/*! Enable / disable dust particle feedback onto the gas [Dust particles]. */
  #ifndef PARTICLES_DUST_FEEDBACK
    #define PARTICLES_DUST_FEEDBACK   YES
  #endif

  #ifndef PARTICLES_DUST_SB_ETA_VK
    #define PARTICLES_DUST_SB_ETA_VK    0.0
  #endif

/*! Particles pusher (\c EXP_MIDPOINT or \c SEMI_IMPLICIT) [Dust particles]. */
  #ifndef PARTICLES_DUST_TIME_STEPPING
    #define PARTICLES_DUST_TIME_STEPPING  EXP_MIDPOINT
  #endif

/*! Set the particle stopping time to be consant (CONSTANT) or
    user-defined (USERDEF). */
  #ifndef PARTICLES_DUST_STOPPING_TIME
    #define PARTICLES_DUST_STOPPING_TIME     CONSTANT
  #endif

/*! Default size of particle binary files */
	#ifndef PARTICLES_DUST_MAX_FILE_SIZE
		 #define PARTICLES_DUST_MAX_FILE_SIZE 1.9 /* In GB */
	#endif

#endif

/* *********************************************************************
    Particle type: Monte-Carlo
   ********************************************************************* */

#if PARTICLES == PARTICLES_MC
/*! Sets the artificial value of the speed of light */
  #ifndef PARTICLES_MC_C
    #if PHYSICS == RMHD
      #define PARTICLES_MC_C       1.0
    #else
      #define PARTICLES_MC_C       10000.0
    #endif
  #endif

/*! Charge to mass ratio for MC particles [MC particles]. */
  #ifndef PARTICLES_MC_E_MC
    #define PARTICLES_MC_E_MC      1.0
  #endif

/*! Charge to mass ratio for the fluid [MC Particles]. */
  #ifndef PARTICLES_MC_E_MC_GAS
    #define PARTICLES_MC_E_MC_GAS  1.0
  #endif

/*! Enable / disable MC particle feedback onto the gas [MC particles]. */
  #ifndef PARTICLES_MC_FEEDBACK
    #define PARTICLES_MC_FEEDBACK  NO
  #endif

  #if PARTICLES_MC_FEEDBACK == NO

  #endif

#endif

/* *********************************************************************
   Partile type: kinetic
   ********************************************************************* */

#if PARTICLES == PARTICLES_KIN

#ifndef PARTICLES_KIN_SOLVER
  #define PARTICLES_KIN_SOLVER = THREE_DIAGONAL
#endif

#ifndef PARTICLES_KIN_OUTPUT_STEP
  #define PARTICLES_KIN_OUTPUT_STEP 1
#endif

/*! Sets the artificial value of the speed of light */
  #ifndef PARTICLES_KIN_C
    #if PHYSICS == RMHD
      #define PARTICLES_KIN_C       1.0
    #else
      #define PARTICLES_KIN_C       10000.0
    #endif
  #endif

/*! Charge to mass ratio for kinetic particles [KIN particles]. */
  #ifndef PARTICLES_KIN_E_MC
    #define PARTICLES_KIN_E_MC      1.0
  #endif

/*! Charge to mass ratio for the fluid [KIN Particles]. */
  #ifndef PARTICLES_KIN_E_MC_GAS
    #define PARTICLES_KIN_E_MC_GAS  1.0
  #endif

/*! Enable / disable KIN particle feedback onto the gas [KIN particles]. */
  #ifndef PARTICLES_KIN_FEEDBACK
    #define PARTICLES_KIN_FEEDBACK  NO
  #endif

/*! Epsilon for time step (advection, acceleration and diffusion)
    [KIN particles] */
  #ifndef PARTICLES_KIN_EPS
    #define PARTICLES_KIN_EPS  0.99
  #endif

#endif

/* *********************************************************************
    Particle type: LAGRANGIAN
   ********************************************************************* */

#if PARTICLES == PARTICLES_LP

  #ifndef PARTICLES_LP_NCOLORS
    #define PARTICLES_LP_NCOLORS  4
  #endif

/*! Default size of particle binary files */
  #ifndef PARTICLES_LP_MAX_FILE_SIZE
		 #define PARTICLES_LP_MAX_FILE_SIZE 1.9 /* In GB */
	#endif
#endif


#if PARTICLES == PARTICLES_LP && PARTICLES_LP_SPECTRA == YES
  #ifndef PARTICLES_LP_NEBINS        /* --> PARTICLES_LP_NEBINS */
    #define PARTICLES_LP_NEBINS 100
  #endif

  #ifndef PARTICLES_LP_SHK_THRESHOLD
    #define PARTICLES_LP_SHK_THRESHOLD   50
  #endif

  #ifndef PARTICLES_LP_SPEC_ENERGY
    #define PARTICLES_LP_SPEC_ENERGY    (CONST_me*CONST_c*CONST_c)   /* Unit energy [ergs] to scale spectra */
  #endif

  #ifndef PARTICLES_LP_SHK_GRADP
    #define PARTICLES_LP_SHK_GRADP   0.1 /* To ensure the particles dont cross the same shock twice */
  #endif

  #ifndef PARTICLES_LP_NONTH_FRACN
    #define  PARTICLES_LP_NONTH_FRACN  0.1
  #endif

  #ifndef PARTICLES_LP_NONTH_FRACE
    #define  PARTICLES_LP_NONTH_FRACE  0.1
  #endif

  #ifndef PARTICLES_LP_MICROETA
    #define PARTICLES_LP_MICROETA      1.4142135623730951
  #endif

  #ifndef PARTICLES_LP_ICCMBZ
   #define PARTICLES_LP_ICCMBZ		0.0
  #endif

#ifndef PARTICLES_LP_CONV_SPECTRA
    #define PARTICLES_LP_CONV_SPECTRA    YES
#endif

  #define UNIT_MAGFIELD       (UNIT_VELOCITY*sqrt(4.0*CONST_PI*UNIT_DENSITY))
  #define UNIT_TIME           (UNIT_LENGTH/UNIT_VELOCITY)
  #define SYNCHROTRON_CONST   (0.0015829*UNIT_MAGFIELD*UNIT_MAGFIELD           \
                               *PARTICLES_LP_SPEC_ENERGY*UNIT_TIME)
#endif

/* *********************************************************************
    Global variables (see globals.h for doc)
   ********************************************************************* */

extern long int p_nparticles;
extern long int p_idCounter;
extern int p_nrestart;

/* *********************************************************************
    Structure definitions
   ********************************************************************* */
#if PARTICLES != PARTICLES_KIN
#if PARTICLES == PARTICLES_CR
#if PARTICLES_CR_GC == NO
typedef struct Particle_{
  double   coord[3];     /**< Particle coordinates */
  double   speed[3];     /**< Particle velocity */
  double   coord_old[3]; /**< Particle coordinates at previous time level */
  double   speed_old[3]; /**< Particles velocity at previous time level */
  double   mass;         /**< Mass density of a single particle         */
  float    tinj;         /**< Particle injection time */
  float    color;        /**< User-supplied real number to distinguish particles */
  int      cell[3];      /**< Indices (i,j,k) of the cell hosting the particle */
  uint32_t id;           /**< Particle id. */
} Particle;
#endif

#if PARTICLES_CR_GC == YES
typedef struct Particle_{
  double   coord[3];     /**< Particle coordinates */
  double   speed[3];     /**< Particle velocity */
  double   coord_old[3]; /**< Particle coordinates at previous time level */
  double   speed_old[3]; /**< Particles velocity at previous time level */
  double   rhs_m1[4];    /**< Right hand side in the past at n-1 */
  double   rhs_m2[4];    /**< Right hand side in the past at n-2 */
  double   mass;         /**< Mass density of a single particle         */
  float    tinj;         /**< Particle injection time */
  float    color;        /**< User-supplied real number to distinguish particles */
  int      cell[3];      /**< Indices (i,j,k) of the cell hosting the particle */
  uint32_t id;           /**< Particle id. */
} Particle;

typedef struct GCaux_{
  double b[4];
  double vE[4];
  double dBf1[4];  /* grad (B/gammaE)  */
  double Bmag;
  double Bmag_inv;
  double vEmag2;
  double gammaE;
  double gammaE2;
  double cEpar;
  double mu;
  double gamma;
  double omega;
  double vpar;
  double upar;
} GCaux;

#endif

#endif

#if PARTICLES == PARTICLES_MC
typedef struct Particle_{
  double   coord[3];     /**< Particle coordinates */
  double   speed[3];     /**< Particle velocity */ //4-velocity?
  double   coord_old[3]; /**< Particle coordinates at previous time level */
  double   speed_old[3]; /**< Particles velocity at previous time level */
  double   mass;         /**< Mass density of a single particle         */
  float    tinj;         /**< Particle injection time */
  float    color;        /**< User-supplied real number to distinguish particles */
  int      cell[3];      /**< Indices (i,j,k) of the cell hosting the particle */
  uint32_t id;           /**< Particle id. */
  double   scattering_time; /**< shows how much time of current time step particle already develops */
} Particle;
#endif

#if PARTICLES == PARTICLES_DUST
typedef struct Particle_{
  double   coord[3];  /**< Particle coordinates */
  double   speed[3];  /**< Particle velocity    */
  double   mass;      /**< Particle mass        */
  double   tau_s;     /**< Particle stopping time */
  float    tinj;      /**< Particle injection time */
  float    color;     /**< User-supplied real number to distinguish particles */
  int      cell[3];   /**< Indices (i,j,k) of the cell hosting the particle */
  uint32_t id;        /**< Particle id. */
} Particle;
#endif

#if (PARTICLES == PARTICLES_LP) && (PARTICLES_LP_SPECTRA == NO)
typedef struct Particle_{
  double coord[3];     /**< Particle coordinates */
  double speed[3];     /**< Particle velocity */
  double coord_old[3]; /**< Particle coordinates at previous time level */
  double speed_old[3]; /**< Particles velocity at previous time level */
  double density;
  float  tinj;         /**< Particle injection time */
  float  color[PARTICLES_LP_NCOLORS]; /**< User-supplied real number to distinguish particles */
  int    cell[3];                  /**< Indices (i,j,k) of the cell hosting the particle */
  uint32_t id;                     /**< Particle id. Numbering is sequential for each rank */
} Particle;
#endif

#if (PARTICLES == PARTICLES_LP) &&  (PARTICLES_LP_SPECTRA == YES)
typedef struct Particle_{
  double coord[3];     /**< Particle coordinates */
  double speed[3];     /**< Particle velocity */
  double coord_old[3]; /**< Particle coordinates at previous time level */
  double speed_old[3]; /**< Particles velocity at previous time level */
  double density;

  double pressure;
  double shk_gradp;   /**< Maximum of grad(p) as the particle travels inside shock */
  double cr;
  double cmp_ratio;   /* can we remove it ???  */
  double lorG;
  double nmicro;

  double Vshk_upst[NVAR];
  double Vshk_dnst[NVAR];
  double eng[PARTICLES_LP_NEBINS+1];
  double chi[PARTICLES_LP_NEBINS];

  double mag[3];      /* can we remove it ???  */
  float  tinj;                      /**< Particle injection time */
  float  color[PARTICLES_LP_NCOLORS];  /**< User-supplied real number to distinguish particles */
  int    cell[3];                   /**< Indices (i,j,k) of the cell hosting the particle */
  uint32_t id;                      /**< Particle id. Numbering is sequential for each rank */

  char shkflag;
} Particle;
#endif
typedef struct particleNode_{
  struct Particle_     p;
  struct particleNode_ *next;
  struct particleNode_ *prev;
} particleNode;
#endif

typedef struct particleProbe_{
  double vg[NVAR];
  double lorentz;
} particleProbe;


/* ***************************************************************
    Macro definitions
   *************************************************************** */

#ifndef PARTICLES_USE_ARRAY
  #define PARTICLES_USE_ARRAY   NO
#endif

/* Useful macro to loop over particles (ATTENTION: do *NOT* use this
  macro inside loops involving creation / destruction of particles) */
#define PARTICLES_LOOP(a,b)   for (a = b; a != NULL; a = a->next)

/* ***************************************************************
    Function Prototypes
   *************************************************************** */

void    Particles_Boundary(Data *, Grid *);
#if PARTICLES != PARTICLES_KIN
bool    CheckParticleLeft(Data *d, Grid *grid, Particle* p);
#endif
void    Particles_BoundaryExchange(Data *, Grid *);
#if PARTICLES != PARTICLES_KIN
int     Particles_BoundaryCheck(Particle *p, Grid *grid);
long    Particles_CheckAll   (particleNode *, int, Grid *);
int     Particles_CheckSingle(Particle *, int, Grid *);
#endif
#if PARTICLES == PARTICLES_CR
void    Particles_CR_ComputeCurrent(const Data *, Grid *);
void    Particles_CR_ComputeForce(Data_Arr, const Data *, Grid *);
void    Particles_CR_ConservativeFeedback(Data_Arr, Data_Arr, double, RBox *);
void    Particles_CR_Flux(const State *, int, int);
void    Particles_CR_Predictor(Data *, double, Grid *);
void    Particles_CR_States1DCopy(const Data *, const Sweep *, int, int);
void    Particles_CR_StatesSource(const Sweep *, double, int, int, Grid *);
void    Particles_CR_Update(Data *, timeStep *, double, Grid *);

#if PARTICLES_CR_GC == YES
void Particles_CR_GC_Convert(Data *, Grid *);
void Particles_CR_GC_Lorentz(Particle *, Data *, Grid *);
int  Particles_CR_GC_RK2(Particle *, Data *, double, Grid *);

void Particles_CR_setGC (Data *, double, Grid *);
int  Particles_CR_getGC (Particle *, Data *, double *, double, Grid *);
int  Particles_CR_GC_RHS (Particle *, Data *, double *, Grid *grid);

void Particles_CR_GCinfo(Particle *, double *);
void Particles_CR_destroyGCParticle (Particle *, particleNode *, Data *);
#endif
#if PARTICLES == PARTICLES_MC
void    Particles_MC_ComputeCurrent(const Data *, Grid *);
void    Particles_MC_ComputeForce(Data_Arr, const Data *, Grid *);
void    Particles_MC_ConservativeFeedback(Data_Arr, Data_Arr, double, RBox *);
void    Particles_MC_EMFields(double *, double *, double *);
void    Particles_MC_Flux(const State *, int, int);
void    Particles_MC_States1DCopy(const Data *, const Sweep *, int, int);
void    Particles_MC_StatesSource(const Sweep *, double, int, int, Grid *);
void    Particles_MC_Update(Data *, timeStep *, double, Grid *);
void    Particles_MC_BallisticMove(Particle* p, double t_step, Grid* grid);
void    Particles_MC_LargeScaleScattering(double* u, const double* vg);
void    Particles_MC_LorentzTransformation(double* u, const double* vg);
void    Particles_MC_AddForce(const Data* data, Grid* grid, Particle* p);
void    Particles_MC_ResetForce(const Data* data);
void    Particles_MC_ComputeForce4(const Data* data);

#if TURBULENT_FIELD == YES
void    Particles_MC_addFlux(Particle* p, double t_step1, double dt, Data* data, Grid* grid, double*** w);
#endif
#endif
#if PARTICLES == PARTICLES_KIN
void    setBoundaryRightPartToZero(Data* data, Grid* grid);
void    Particles_KIN_ComputeCurrent(const Data *, Grid *);
void    Particles_KIN_ComputeForce(Data_Arr, const Data *, Grid *);
void    Particles_KIN_ConservativeFeedback(Data_Arr, Data_Arr, double, RBox *);
void    Particles_KIN_EMFields(double *, double *, double *);
void    Particles_KIN_Flux(const State *, int, int);
void    Particles_KIN_States1DCopy(const Data *, const Sweep *, int, int);
void    Particles_KIN_StatesSource(const Sweep *, double, int, int, Grid *);
void    Particles_KIN_Update(Data *, timeStep *, double, Grid *);
void    Particles_KIN_LorentzTransformation(double* u, const double* vg);
void    Particles_KIN_ResetForce(const Data* data);
void    Particles_KIN_ComputeForce4(const Data* data);
double evaluateDiffusionCoefficient(Data* data, int i, int j, int k, double u);
#endif

#endif

#if PARTICLES != PARTICLES_KIN
void    Particles_Density(Particle *, double *);
void    Particles_Deposit(particleNode *, void (*Func)(Particle *, double *),
                          Data_Arr, int, Grid *);
void    Particles_DepositBoundaryExchange(Data_Arr, int, Grid *);
void    Particles_Destroy(particleNode *, Data *);
void    Particles_Display(Particle*);
#endif

#if PARTICLES == PARTICLES_DUST
void    Particles_Dust_ComputeForce (Data_Arr, Data *, Grid *);
void    Particles_Dust_ConservativeFeedback(Data_Arr, Data_Arr, double, RBox *);
void    Particles_Dust_StatesSource(Sweep *, Data_Arr, double, int, int, Grid *);
double  Particles_Dust_StoppingTime(double *, Particle *);
void    Particles_Dust_Update(Data *, timeStep *, double, Grid *);
#endif

void    Particles_Init(Data *, Grid *);
void    Particles_Inject(Data *d, Grid *grid);

#if PARTICLES != PARTICLES_KIN
int     Particles_Insert(Particle *, Data *, char, Grid *);
void    Particles_GetWeights (Particle *, int *, double ***, Grid *);
#endif

double  Particles_Interpolate(double ***, double ***, int *);
void    Particles_InterpolateArr(Data_Arr, int, double ***, int *, double *);

void    Particles_ListToArray (Data *d);

void    Particles_LoadRandom(double *xbeg, double *xend,
                          double (*DistribFunc)(double, double, double),
                          double coor[]);
void    Particles_LoadUniform(int, int, double *, double *, double *);

int     Particles_LocateCell(double *, int *, Grid *);

#if PARTICLES == PARTICLES_LP
void    Particles_LP_Update (Data *, timeStep *, double, Grid *);
#if PARTICLES_LP_SPECTRA == YES

void    Particles_LP_ComputeShockNormalSpeed(Particle *, double *, double *);
void    Particles_LP_FixValue(Particle *, particleProbe *, Data *, Grid *);
void    Particles_LP_FlagShock(Data *, char ***, Grid *);
void    Particles_LP_GradP(double*,int, double ***u[], Grid *, int indici[]);
double  Particles_LP_GetEmaxAtShock(Particle *, particleProbe *, double, Grid *);
#if PARTICLES_LP_CONV_SPECTRA == NO
double  Particles_LP_GetEminAtShock(Particle *, particleProbe *, double, double);
#endif

void    Particles_LP_Get_4vel(double *, double *, double *);
void    Particles_LP_Get_3vel(double *, double *, double *);
void    Particles_LP_Get_4mag(double *, double, double *);
void    Particles_LP_Get_3mag(double *, double, double *);

void    Particles_LP_IC_Emissivity(Particle *, double, double, double *);
void    Particles_LP_InitSpectra(Particle*);
void    Particles_LP_IntegrateSpectra(Particle *, double, double *);
void    Particles_LP_SampleShock (Particle *, particleProbe *, Grid *);
void    Particles_LP_Spectra(Data *, float ***, Particle*, Grid *, double);  // Useless

#if PARTICLES_LP_CONV_SPECTRA
void    Particles_LP_UpdateUpstreamSpectra(Particle*, Grid *, double, double,
                                           long int, double, double);
#else
void    Particles_LP_UpdateUpstreamSpectra(Particle* , double , double ,
                                           double );
#endif

void    Particles_LP_UpdateSpectra(Data *, double, Grid *);
void    Particles_LP_Sync_Emissivity(Particle*, double, double, double *, double *, double *);
#endif
#endif

void    Particles_Restart(Runtime *, Data *, int, Grid *);
void    Particles_Set(Data *, Grid *);
void    Particles_SetOutput (Data *, Runtime *);

#if PARTICLES != PARTICLES_KIN

Particle *Particles_Select(particleNode *, int);
void    Particles_SetID(particleNode *);
int     Particles_Number(particleNode *);
void    Particles_ShowList(particleNode *, int);
void    Particles_WriteTrajectory (Particle *, char);
#endif
void    Particles_UserDefBoundary(Data *d, int, Grid *);

void    Particles_WriteBinary(Data* d, Grid* grid, double, Output *, char *);
void    Particles_WriteData(Data *d, Output *, Grid *);
void    Particles_WriteTab   (Data* d, Grid* grid, char filename[]);
void    Particles_WriteVTK   (Data* d, Grid* grid, Output *, char filename[]);

#ifdef PARALLEL
 extern MPI_Datatype MPI_PARTICLE;
 extern MPI_Datatype PartOutputType;
 void Particles_StructDatatype();
#endif
