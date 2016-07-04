!> list of aliases used in maphys
!! 
!! Here we give aliases to a few array components, namely :
!!   - ICNTL, RCNTL   : user controls 
!!     ( see doc_icntl.inc and doc_rcntl.inc )
!!   - IKEEP, RKEEP   : internal controls 
!!     ( see doc_ikeep.inc and doc_rkeep.inc )
!!   - IINFO, RINFO   : informations provided to the user (local) 
!!     ( see doc_iinfo.inc and doc_rinfo.inc )
!!   - IINFOG, RINFOG : informations provided to the user (global) 
!!     ( see doc_iinfog.inc and doc_rinfog.inc )
!!   .
!! And also to several enumerations :
!!   - SDS_*  
!!
!! @author Yohan Lee-tin-yien
!!
Module FTS_ftsolver_enum
  Implicit none

  !-----------------------------------------------------------------------------
  ! [1] Integer controls aliases
  !-----------------------------------------------------------------------------

  !-------------------------------------------------------------------
  ! [1.0] Partitionner Related
  !-------------------------------------------------------------------

  ! Icntl
  Integer, Parameter :: ICNTL_PART_Strategy    = 7

  ! Ikeep
  Integer, Parameter :: IKEEP_PART_Strategy  = 1

  ! Significant values
  Integer, Parameter :: PART_STRATEGY_isUnset         = 0
  Integer, Parameter :: PART_STRATEGY_isMetisNodeND   = 1
  Integer, Parameter :: PART_STRATEGY_isMetisEdgeND   = 2
  Integer, Parameter :: PART_STRATEGY_isMetisNodeWND  = 3
  Integer, Parameter :: PART_STRATEGY_isScotchNodeND  = 4

  !-------------------------------------------------------------------
  ! [1.1] Sparse direct Solver related
  !-------------------------------------------------------------------
  
  ! Icntl
  Integer, Parameter :: ICNTL_SDS_Default  = 13
  Integer, Parameter :: ICNTL_SDS_Facto    = 32
  Integer, Parameter :: ICNTL_SDS_Precond  = 15

  ! Ikeep
  Integer, Parameter :: IKEEP_SDS_Facto    = 32
  Integer, Parameter :: IKEEP_SDS_Precond  = 15

  ! Significant values
  Integer, Parameter :: SDS_DEFAULT_Mumps      = 1
  Integer, Parameter :: SDS_DEFAULT_Pastix     = 2
  Integer, Parameter :: SDS_DEFAULT_UsrDefined = 3
  
  ! 2nd level of parallelism

  ! Icntl
  Integer, Parameter :: ICNTL_2LVLS_Bind          = 36
  Integer, Parameter :: ICNTL_2LVLS_NNodes        = 37
  Integer, Parameter :: ICNTL_2LVLS_NcoresPerNode = 38
  Integer, Parameter :: ICNTL_2LVLS_NThrdsPerProc = 39
  Integer, Parameter :: ICNTL_2LVLS_NProcs        = 40
  Integer, Parameter :: ICNTL_FTSOLVERBIND          = 42

  ! Ikeep
  Integer, Parameter :: IKEEP_NODES    = 5
  Integer, Parameter :: IKEEP_NODEID   = 6
  Integer, Parameter :: IKEEP_NODECOMM = 7

  !-------------------------------------------------------------------
  ! [1.2] ILU Related
  !-------------------------------------------------------------------

  Integer, Parameter :: ICNTL_ILU_LUFillPct    = 10
  Integer, Parameter :: ICNTL_ILU_SCHURFillPct = 11

  Integer, Parameter :: ICNTL_ILU_LUFill        = 8
  Integer, Parameter :: ICNTL_ILU_SCHURFill     = 9
  Integer, Parameter :: RCNTL_ILU_LUThreshold    = 8
  Integer, Parameter :: RCNTL_ILU_SCHURThreshold = 9

  Integer, Parameter :: IKEEP_ILU_SCHURFill     = 8
  Integer, Parameter :: IKEEP_ILU_LUFill        = 9
  Integer, Parameter :: RKEEP_ILU_SCHURThreshold = 8
  Integer, Parameter :: RKEEP_ILU_LUThreshold    = 9


  !-------------------------------------------------------------------
  ! [1.3] Schur Related
  !-------------------------------------------------------------------

  ! General
  Integer, Parameter :: ICNTL_SCHUR_Strategy      = 30
  Integer, Parameter :: IKEEP_SCHUR_Strategy      = 30
  Integer, Parameter :: IKEEP_SCHUR_UselessAfter  = 17

  ! When Exact
  Integer, Parameter :: ICNTL_EXSCHUR_BlocSize  = 31
  Integer, Parameter :: IKEEP_EXSCHUR_BlocSize  = 31

  ! Significant values
  Integer, Parameter :: SCHUR_STRATEGY_Exact          = 0
  Integer, Parameter :: SCHUR_STRATEGY_ExactFromFacto = 1
  Integer, Parameter :: SCHUR_STRATEGY_ApprxWithIluT  = 2

  !-------------------------------------------------------------------
  ! [1.4] Preconditioner related
  !-------------------------------------------------------------------

  ! General 
  Integer, Parameter :: ICNTL_PCD_Strategy  = 21
  Integer, Parameter :: IKEEP_PCD_Strategy  = 21
  Integer, Parameter :: RCNTL_PCD_SparsifyThreshold = 11
  Integer, Parameter :: RKEEP_PCD_SparsifyThreshold = 11

  

  ! Significant values
  Integer, Parameter :: PCD_STRATEGY_Unset           =-1 
  Integer, Parameter :: PCD_STRATEGY_isLocalExact    = 1 
  Integer, Parameter :: PCD_STRATEGY_isLocalApprox   = 2
  Integer, Parameter :: PCD_STRATEGY_isForcedByILUT  = 3
  Integer, Parameter :: PCD_STRATEGY_isNone       = 4
  Integer, Parameter :: PCD_STRATEGY_isAutodetected  = 5
  Integer, Parameter :: PCD_STRATEGY_isGlobalApprox  = 6


  !-------------------------------------------------------------------
  ! [1.5] Iterative Solver related (GMRES,CG)
  !-------------------------------------------------------------------
  
  ! Icntl
  Integer, Parameter :: ICNTL_ITS_Solver      = 20
  Integer, Parameter :: ICNTL_ITS_OrtStrategy = 22
  Integer, Parameter :: ICNTL_ITS_InitGuess   = 23
  Integer, Parameter :: ICNTL_ITS_MaxIters    = 24
  Integer, Parameter :: ICNTL_ITS_ResStrategy = 25
  Integer, Parameter :: ICNTL_ITS_Restart     = 26
  Integer, Parameter :: ICNTL_ITS_MatVect     = 27

  ! Ikeep
  Integer, Parameter :: IKEEP_ITS_Solver      = 20
  Integer, Parameter :: IKEEP_ITS_OrtStrategy = 22
  Integer, Parameter :: IKEEP_ITS_InitGuess   = 23
  Integer, Parameter :: IKEEP_ITS_MaxIters    = 24
  Integer, Parameter :: IKEEP_ITS_ResStrategy = 25
  Integer, Parameter :: IKEEP_ITS_Restart     = 26
  Integer, Parameter :: IKEEP_ITS_MatVect     = 27

  ! Rcntl
  Integer, Parameter :: RCNTL_ITS_Tolerance = 21

  ! Rkeep
  Integer, Parameter :: RKEEP_ITS_Tolerance = 1

  ! Significant Values
  Integer, Parameter :: ITS_Solver_Unset       = 0
  Integer, Parameter :: ITS_Solver_isPackGMRES = 1
  Integer, Parameter :: ITS_Solver_isPackCG    = 2
  Integer, Parameter :: ITS_Solver_isAuto      = 3
  Integer, Parameter :: ITS_Solver_isPackFGMRES = 4

  Integer, Parameter :: ORT_STRATEGY_Unset      = -1
  Integer, Parameter :: ORT_STRATEGY_IsMGS      =  0
  Integer, Parameter :: ORT_STRATEGY_IsIMGS     =  1
  Integer, Parameter :: ORT_STRATEGY_IsCGS      =  2
  Integer, Parameter :: ORT_STRATEGY_IsICGS     =  3

  Integer, Parameter :: INIT_GUESS_Unset   = -1
  Integer, Parameter :: INIT_GUESS_isSet   =  1
  Integer, Parameter :: INIT_GUESS_isUnset =  0

  Integer, Parameter :: RES_STRATEGY_Unset            = -1
  Integer, Parameter :: RES_STRATEGY_isRecurrence     =  0
  Integer, Parameter :: RES_STRATEGY_isMatVectProduct =  1

  Integer, Parameter :: MAT_VECT_Unset          = -1
  Integer, Parameter :: MAT_VECT_isAutodetected = 0
  Integer, Parameter :: MAT_VECT_isExplicit     = 1
  Integer, Parameter :: MAT_VECT_isImplicit     = 2

  !-------------------------------------------------------------------
  ! [1.6] IINFO
  !-------------------------------------------------------------------

  Integer, Parameter :: IINFO_STATUS           =  1
  Integer, Parameter :: IINFO_STATUS_DETAIL    =  2
  Integer, Parameter :: IINFO_STRAT_PART       =  3
  Integer, Parameter :: IINFO_LCMAT_ORDER      =  4
  Integer, Parameter :: IINFO_LCMAT_SIZEOF    =  5
  Integer, Parameter :: IINFO_LCMAT_NBENTRIES  =  6
  Integer, Parameter :: IINFO_STRAT_FACTO      =  7
  Integer, Parameter :: IINFO_STRAT_SCHUR      =  8   
  Integer, Parameter :: IINFO_LCFACTO_SIZEOF  =  9
  Integer, Parameter :: IINFO_LCFACTO_NBPIVOTS = 10
  Integer, Parameter :: IINFO_SCHUR_ORDER      = 11
  Integer, Parameter :: IINFO_SCHUR_SIZEOF    = 12
  Integer, Parameter :: IINFO_SCHUR_NBENTRIES  = 13
  Integer, Parameter :: IINFO_STRAT_PCD        = 14
  Integer, Parameter :: IINFO_STRAT_PCDSDS     = 15
  Integer, Parameter :: IINFO_PCD_ORDER        = 16
  Integer, Parameter :: IINFO_PCD_SIZEOF      = 17
  Integer, Parameter :: IINFO_PCD_NBPIVOTS     = 18
  Integer, Parameter :: IINFO_PCD_PERKEPT      = 19
  Integer, Parameter :: IINFO_STRAT_ITS        = 20
  Integer, Parameter :: IINFO_STRAT_ITSSMV     = 21
  Integer, Parameter :: IINFO_FAC_SDSMEMPEAK   = 22
  Integer, Parameter :: IINFO_PCD_SDSMEMPEAK   = 23
  Integer, Parameter :: IINFO_FTSOLVER_MEMPEAK   = 24
  Integer, Parameter :: IINFO_FTSOLVER_MEMUSED   = 25
  Integer, Parameter :: IINFO_ITS_MEMUSED      = 26
  Integer, Parameter :: IINFO_ITS_MEMPEAK      = 27
  Integer, Parameter :: IINFO_BUFF_MEMUSED     = 28
  Integer, Parameter :: IINFO_BUFF_MEMPEAK     = 29
  Integer, Parameter :: IINFO_LIBSDS_MEMUSED   = 30
  Integer, Parameter :: IINFO_LIBSDS_MEMPEAK   = 31
  Integer, Parameter :: IINFO_NODEID           = 32
  Integer, Parameter :: IINFO_NODEDOMAINS      = 33
  Integer, Parameter :: IINFO_NODE_MEMUSED     = 34
  Integer, Parameter :: IINFO_NODE_MEMPEAK     = 35



  !-------------------------------------------------------------------
  ! [1.7] RINFO
  !-------------------------------------------------------------------

  Integer, Parameter :: RINFO_TIMING_Fault        = 1
  Integer, Parameter :: RINFO_TIMING_Solve        = 2
  Integer, Parameter :: RINFO_TIMING_Interp       = 3
  Integer, Parameter :: RINFO_TIMING_InterpRhs    = 4
  Integer, Parameter :: RINFO_TIMING_InterpSlv    = 5
  Integer, Parameter :: RINFO_TIMING_Restart      = 6
  Integer, Parameter :: RINFO_FLOP_LI             = 7
  Integer, Parameter :: RINFO_FLOP_LSI            = 8
  Integer, Parameter :: RINFO_M_QR                = 9
  Integer, Parameter :: RINFO_N_QR                = 10
  !-------------------------------------------------------------------
  ! [1.8] RINFOG
  !-------------------------------------------------------------------

  Integer, Parameter :: RINFOG_SOL_NORM22         = 1
  Integer, Parameter :: RINFOG_BACKWARD_ERROR     = 2
  Integer, Parameter :: RINFOG_ITS_BckErr         = 3
  Integer, Parameter :: RINFOG_BckwrdErrorDist    = 4

  !-------------------------------------------------------------------
  ! [1.9] IINFOG
  !-------------------------------------------------------------------

  Integer, Parameter :: IINFOG_STATUS        = 1
  Integer, Parameter :: IINFOG_STATUS_DETAIL = 2
  Integer, Parameter :: IINFOG_MAT_ORDER     = 3
  Integer, Parameter :: IINFOG_MAT_NBENTRIES = 4
  Integer, Parameter :: IINFOG_ITS_NbIters   = 5

  !-------------------------------------------------------------------
  ! [1.*] Others
  !-------------------------------------------------------------------

  Integer, Parameter :: ICNTL_OUTPUT_ErrUnit   = 1
  Integer, Parameter :: ICNTL_OUTPUT_WarnUnit  = 2
  Integer, Parameter :: ICNTL_OUTPUT_StdUnit   = 3
  Integer, Parameter :: ICNTL_OUTPUT_Verbosity = 4

  Integer, Parameter :: ICNTL_PRINT_Cntls     = 5
  Integer, Parameter :: ICNTL_PRINT_Infos     = 6
  Integer, Parameter :: IKEEP_PRINT_Cntls     = 13
  Integer, Parameter :: IKEEP_PRINT_Infos     = 14
  Integer, Parameter :: PRINT_Not             = 0
  Integer, Parameter :: PRINT_Once            = 1
  Integer, Parameter :: PRINT_EveryStep       = 2

  Integer, Parameter :: ICNTL_SHIFT_Do         = 12

  Integer, Parameter :: IKEEP_HOSTRANK = 4
  Integer, Parameter :: IKEEP_MPIRANK  = 3
  Integer, Parameter :: IKEEP_CURJOB   = 2
  Integer, Parameter :: IKEEP_NBDOMAINS   = 16
  Integer, Parameter :: CURJOB_IsInit     = -1
  Integer, Parameter :: CURJOB_IsAnalysis =  1
  Integer, Parameter :: CURJOB_IsFacto    =  2
  Integer, Parameter :: CURJOB_IsPrecond  =  3
  Integer, Parameter :: CURJOB_IsSolve    =  4
  Integer, Parameter :: CURJOB_IsFinalize = -2

  Integer, Parameter :: ICNTL_INSYSTEM = 43
  Integer, Parameter :: IKEEP_INSYSTEM = 12
  Integer, Parameter :: INSYSTEM_IsUnset       = 0
  Integer, Parameter :: INSYSTEM_IsCentralized = 1
  Integer, Parameter :: INSYSTEM_IsDistributed = 2
  Integer, Parameter :: IKEEP_SYMMETRY = 18
  Integer, Parameter :: FAULT_SIMPLE = 1
  Integer, Parameter :: FAULT_DOUBLE = 2

  

End Module FTS_ftsolver_enum
