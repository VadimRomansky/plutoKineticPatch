#define  PHYSICS                        RMHD
#define  DIMENSIONS                     1
#define  GEOMETRY                       SPHERICAL
#define  BODY_FORCE                     VECTOR
#define  COOLING                        NO
#define  RECONSTRUCTION                 LINEAR
#define  TIME_STEPPING                  RK2
#define  NTRACER                        1
#define  PARTICLES                      PARTICLES_KIN
#define  USER_DEF_PARAMETERS            7

/* -- physics dependent declarations -- */

#define  EOS                            IDEAL
#define  ENTROPY_SWITCH                 SELECTIVE
#define  RADIATION                      NO
#define  DIVB_CONTROL                   DIV_CLEANING

/* -- user-defined parameters (labels) -- */

#define  M_COR                          0
#define  E_COR                          1
#define  Rf                             2
#define  T_AMB                          3
#define  Vf                             4
#define  RHO_ISM                        5
#define  B_AMB                          6

/* [Beg] user-defined constants (do not change this line) */

#define  INITIAL_SMOOTHING              YES
#define  INTERNAL_BOUNDARY              YES
#define  FAILSAFE                       YES
#define  SHOW_TIME_STEPS                YES
#define  SHOCK_FLATTENING               MULTID
#define  LIMITER                        VANLEER_LIM
#define  UNIT_DENSITY                   CONST_mp
#define  UNIT_LENGTH                    3.086e17
#define  UNIT_VELOCITY                  CONST_c
#define  PARTICLES_CR_E_MC              (CONST_e/(CONST_mp*CONST_c))*UNIT_LENGTH*sqrt(4*CONST_PI*UNIT_DENSITY)
#define  PARTICLES_CR_C                 CONST_c/UNIT_VELOCITY
#define  PARTICLES_CR_FEEDBACK          NO
#define  PARTICLES_MC_E_MC              (CONST_e/(CONST_mp*CONST_c))*UNIT_LENGTH*sqrt(4*CONST_PI*UNIT_DENSITY)
#define  PARTICLES_MC_C                 CONST_c/UNIT_VELOCITY
#define  PARTICLES_MC_FEEDBACK          NO
#define  PARTICLES_CR_NSUB              1
#define  PARTICLES_KIN_E_MC             ((CONST_e/(CONST_mp*CONST_c))*UNIT_LENGTH*sqrt(4*CONST_PI*UNIT_DENSITY))
#define  PARTICLES_KIN_C                (CONST_c/UNIT_VELOCITY)
#define  PARTICLES_KIN_MASS             (CONST_mp/(UNIT_DENSITY*UNIT_LENGTH*UNIT_LENGTH*UNIT_LENGTH))
#define  PARTICLES_KIN_FEEDBACK         NO
#define  PARTICLES_KIN_OUTPUT_STEP      8
#define  PARTICLES_KIN_SOLVER           THREE_DIAGONAL
#define  PARTICLES_KIN_EPS              0.2
#define  P_GRID_MIN                     1000
#define  P_GRID_MAX                     1E5
#define  INJECTION_PARAMETER            1E-5
#define  COMPRESSION_TRESHOLD           1.4
#define  TURBULENT_FIELD                NO
#define  TURBULENCE_OUTPUT_STEP         8
#define  TURBULENCE_SOLVER              THREE_DIAGONAL
#define  NTURB                          10
#define  NMOMENTUM                      20
#define  MAX_GMRES_ITERATIONS           20
#define  WARNING_MESSAGES               NO
#define  MULTIPLE_LOG_FILES             NO
#define  CONST_Year                     3.15E7

/* [End] user-defined constants (do not change this line) */
