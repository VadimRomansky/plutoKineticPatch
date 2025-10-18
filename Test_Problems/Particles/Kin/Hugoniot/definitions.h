#define  PHYSICS                        RMHD
#define  DIMENSIONS                     2
#define  GEOMETRY                       CARTESIAN
#define  BODY_FORCE                     NO
#define  COOLING                        NO
#define  RECONSTRUCTION                 LINEAR
#define  TIME_STEPPING                  RK2
#define  NTRACER                        1
#define  PARTICLES                      PARTICLES_KIN
#define  USER_DEF_PARAMETERS            4

/* -- physics dependent declarations -- */

#define  EOS                            IDEAL
#define  ENTROPY_SWITCH                 SELECTIVE
#define  RADIATION                      NO
#define  DIVB_CONTROL                   DIV_CLEANING

/* -- user-defined parameters (labels) -- */

#define  RHO_d                          0
#define  V_d                            1
#define  B_AMB                          2
#define  SIGMA                          3

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
#define  PARTICLES_CR_E_MC              ((CONST_e/(CONST_mp*CONST_c))*UNIT_LENGTH*sqrt(4*CONST_PI*UNIT_DENSITY))
#define  PARTICLES_CR_C                 (CONST_c/UNIT_VELOCITY)
#define  PARTICLES_CR_FEEDBACK          NO
#define  PARTICLES_MC_E_MC              ((CONST_e/(CONST_mp*CONST_c))*UNIT_LENGTH*sqrt(4*CONST_PI*UNIT_DENSITY))
#define  PARTICLES_MC_C                 (CONST_c/UNIT_VELOCITY)
#define  PARTICLES_MC_FEEDBACK          NO
#define  PARTICLES_CR_NSUB              1
#define  PARTICLES_KIN_E_MC             ((CONST_e/(CONST_mp*CONST_c))*UNIT_LENGTH*sqrt(4*CONST_PI*UNIT_DENSITY))
#define  PARTICLES_KIN_C                (CONST_c/UNIT_VELOCITY)
#define  PARTICLES_KIN_MASS             (CONST_mp/(UNIT_DENSITY*UNIT_LENGTH*UNIT_LENGTH*UNIT_LENGTH))
#define  PARTICLES_KIN_FEEDBACK         NO
#define  PARTICLES_KIN_OUTPUT_STEP      1
#define  PARTICLES_KIN_SOLVER           THREE_DIAGONAL
#define  PARTICLES_KIN_EPS              0.1
#define  P_GRID_MIN                     100
#define  P_GRID_MAX                     1E7
#define  INJECTION_PARAMETER            1E-7
#define  COMPRESSION_TRESHOLD           1.4
#define  TURBULENT_FIELD                YES
#define  TURBULENCE_OUTPUT_STEP         1
#define  NTURB                          50
#define  NMOMENTUM                      100
#define  MAX_GMRES_ITERATIONS           20
#define  WARNING_MESSAGES               NO
#define  MULTIPLE_LOG_FILES             NO
#define  CONST_Year                     3.15E7

/* [End] user-defined constants (do not change this line) */
