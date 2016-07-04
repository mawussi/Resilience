! Warning: XFTS_GENFILE_COMMENT
#include "fts_defs_f.h"
#include "fts_macros_f.h"

! [+] module : XFTS_ftsolver_aux_mod ---------------------------------------------
!
!> Auxialiary module for ftsolver
!!
!! Basically contains useful routines for ftsolver.
!!
Module XFTS_ftsolver_mod

      Use fts_ftsolver_enum
      Use XFTS_ftsolver_type
      Use XFTS_sls_mod
      Use FTS_error_mod
      Implicit None

  !* Private constants *!
  Character(len=FTSOLVER_STRL), Private, Parameter :: FLNAME= &
       "XFTS_ARITHfts_ftsolver_aux_mod.F90"

  ! List of routines
  Public :: XFTS_FTSOLVER_Nullify
  Public :: XFTS_FTSOLVER_Free
  Public :: XFTS_FTSOLVER_init    
  Public :: XFTS_FTSOLVER_PretreatInputMatrix
  Public :: XFTS_FTSOLVER_Set_default_icntl  
Contains

    ! [+] routine : XFTS_FTSOLVER_Init  -----------------------------------------------
    !
    !> initialize a ftsl_instance
    !! 
    !! @param[in,out ] ftsl     the ftsl instance to initialize
    !!

    !!
    Subroutine XFTS_ftsolver_Init(ftsl) ! inout

      !* Modules & co *!

      Use fts_env_mod
      Implicit None
      Include "mpif.h"

      !* Arguments *!

      Type(XFTS_ftsolver_t), Intent(inout) :: ftsl

      !* Local variables *!

      ! Scalars
      Integer                      :: iinfo
      Integer                      :: rank,np

      !- End of header ---------------------------------------------------------

      !-------------------------------------------------------------------------
      ! [-] Initialize local variables
      !-------------------------------------------------------------------------

      iinfo = FTS_SUCCESS
      
      ! temporaly set the logging error to standard output.


      Call MPI_Comm_size(ftsl%comm, np  , iinfo )
      ASSRT(iinfo == MPI_SUCCESS)
      Call MPI_Comm_rank(ftsl%comm, rank, iinfo )
      ASSRT(iinfo == MPI_SUCCESS)

      !-------------------------------------------------------------------------
      ! [-] Nullify all components & Give default values to parameters
      !-------------------------------------------------------------------------

      Call XFTS_FTSOLVER_Nullify(ftsl,iinfo)
      FTS_ONFAILURE_GOTO9999(iinfo)
      Call XFTS_FTSOLVER_Set_default_icntl(ftsl%icntl,iinfo)
      FTS_ONFAILURE_GOTO9999(iinfo)

      !-------------------------------------------------------------------------
      ! [-] Set other components
      !-------------------------------------------------------------------------


      Call fts_env_init(ftsl%env)
      Call fts_env_setMPIComm(ftsl%env,ftsl%comm)

      ftsl%ikeep(IKEEP_MPIRANK)  = rank
      ftsl%ikeep(IKEEP_HOSTRANK) = 0
      ftsl%ikeep(IKEEP_NBDOMAINS) = np
      ftsl%ikeep(IKEEP_NODES   ) = fts_env_getNbNodes (ftsl%env)
      ftsl%ikeep(IKEEP_NODEID  ) = fts_env_getNodeId  (ftsl%env)
      ftsl%ikeep(IKEEP_NODECOMM) = fts_env_getNodeComm(ftsl%env)

      ftsl%iinfo(IINFO_NODEID)   = ftsl%ikeep(IKEEP_NODEID)
      Call MPI_Comm_size(&
           fts_env_getNodeComm(ftsl%env), &
           ftsl%iinfo(IINFO_NODEDOMAINS),&
           iinfo)
      ASSRT( iinfo == MPI_SUCCESS )

      !-------------------------------------------------------------------------
      ! [*] Finish
      !-------------------------------------------------------------------------
      
9999  Continue
      ftsl%iinfo( IINFO_STATUS ) =  iinfo

    end Subroutine XFTS_ftsolver_Init





  !routine : XFTS_FTSOLVER_Nullify  ----------------------------------------------
  !
  !> nullify a ftsolver_instance
  !! 
  !! @param[in,out ] ftsl     the ftsolver instance to nullify
  !!
  !! @author Yohan Lee-tin-yien
  !! @todo   nullify derived type : ftsl%dls_precond_schur%dds
  !!
  Subroutine XFTS_ftsolver_Nullify & ! intents
       (ftsl ,           & ! inout
       info )              ! out
    
    !* Modules *!

    Use XFTS_ftsolver_type
    Use XFTS_sparse_matrix_mod, Only : &
         XFTS_sm_nullify ! routine
    Use XFTS_dense_matrix_mod, Only : &
         XFTS_dm_nullify ! routine

  
    implicit none

    !* Subroutine arguments *!
    Type(XFTS_ftsolver_t), intent(inout) :: ftsl
    integer       , intent(  out) :: info

    !* Local variables *!

    ! contants
    integer       , parameter :: DEFAULT_INT   = -9999
    real(kind=8)  , parameter :: DEFAULT_REAL  = -9999.d0

    ! others
    integer                       :: iinfo

    !- End of header -----------------------------------------------------------

    !---------------------------------------------------------------------------
    ! [1] USER PARAMETERS
    !---------------------------------------------------------------------------

    !-------------------------------------------------------------------
    ! [1.1] MPI communicator (do not touch)
    !-------------------------------------------------------------------
    ! ftsl%comm = 

    !-------------------------------------------------------------------
    ! [1.2] Input matrix (in coordinate format)
    !-------------------------------------------------------------------
    ftsl%sym = DEFAULT_INT
    ftsl%n   = DEFAULT_INT
    ftsl%nnz = DEFAULT_INT
    nullify( ftsl%rows )
    nullify( ftsl%cols )
    nullify( ftsl%values )
    ftsl%write_matrix = ""

    !-------------------------------------------------------------------
    ! [1.3] Right-hand-side & Solution (in dense format ordered by columns)
    !-------------------------------------------------------------------

    ftsl%nrhs  = DEFAULT_INT
    nullify( ftsl%rhs ) 
    nullify( ftsl%sol ) 
    nullify( ftsl%lc_rhs ) 

    !-------------------------------------------------------------------
    ! [1.4] Controls
    !-------------------------------------------------------------------

    ftsl%job    = DEFAULT_INT


    !-------------------------------------------------------------------
    ! [1.5] Statistics
    !-------------------------------------------------------------------

    ! on this process (MPI)
    ftsl%iinfo = -1
    ftsl%rinfo = -1

    ! on all process (MPI)
    ftsl%iinfog = -1
    ftsl%rinfog = -1

    ftsl%iinfomin = -1
    ftsl%iinfomax = -1
    ftsl%iinfoavg = -1
    ftsl%iinfosig = -1
    ftsl%rinfomin = -1
    ftsl%rinfomax = -1
    ftsl%rinfoavg = -1
    ftsl%rinfosig = -1

    !---------------------------------------------------------------------------
    ! [2] Internal working data
    !---------------------------------------------------------------------------


    !-------------------------------------------------------------------
    ! [2.1] internal controls
    !-------------------------------------------------------------------

    ftsl%ikeep  = DEFAULT_INT
    ftsl%rkeep  = DEFAULT_REAL

    !-------------------------------------------------------------------
    ! [2.2]  Description of the Non-overlapping domain Decomposition 
    !-------------------------------------------------------------------

    !---------------------------------------------------------
    ! [2.2.1] Local domain description 
    !---------------------------------------------------------
    ! Scalars
      ftsl%lc_domain%size_domain         = DEFAULT_INT
      ftsl%lc_domain%size_domain_delta1  = DEFAULT_INT
      ftsl%lc_domain%size_domain_delta2  = DEFAULT_INT
      ftsl%lc_domain%size_vertperneig    = DEFAULT_INT
      ftsl%lc_domain%size_neig_edge      = DEFAULT_INT
      ftsl%lc_domain%size_sec_neig_edge  = DEFAULT_INT
      ftsl%lc_domain%size_secin_neig_edge = DEFAULT_INT


      ftsl%lc_domain%nb_neig             = DEFAULT_INT
      ftsl%lc_domain%nb_sec_neig         = DEFAULT_INT
      ftsl%lc_domain%nb_secin_neig       = DEFAULT_INT


    ! Arrays

      Nullify( ftsl%lc_domain%index_local2global      )
      Nullify( ftsl%lc_domain%index_neig      )
      Nullify( ftsl%lc_domain%ptr_direct_neig_edge   )
      Nullify( ftsl%lc_domain%ptr_vertperneig         )
      Nullify( ftsl%lc_domain%index_vertperneig       )
      Nullify( ftsl%lc_domain%ptr_neig_edge           )
      Nullify( ftsl%lc_domain%index_neig_edge         )
      Nullify( ftsl%lc_domain%index_sec_neig          )
      Nullify( ftsl%lc_domain%ptr_sec_neig_edge       )
      Nullify( ftsl%lc_domain%index_sec_neig_edge     )
      Nullify( ftsl%lc_domain%index_secin_neig        )
      Nullify( ftsl%lc_domain%ptr_secin_neig_edge     )
      Nullify( ftsl%lc_domain%index_secin_neig_edge   )

    

    !---------------------------------------------------------
    ! [2.2.3] Blocs on the matrix of the domain 
    !---------------------------------------------------------

    ! [  Aii     Aib ] 
    ! - i: interior  related (nodes inside the domain)
    ! - b: border    related (nodes on its interface )
    Call XFTS_sm_nullify(ftsl%sm_blockrow_delta1,iinfo)  
    Call XFTS_sm_nullify(ftsl%sm_blockrow,iinfo)
    Call XFTS_sls_nullify(ftsl%sm_precond,iinfo)
    Call XFTS_sm_nullify(ftsl%sm_lsi,iinfo)
    Call XFTS_sm_nullify(ftsl%sm_li,iinfo)


    !---------------------------------------------------------
    ! [2.2.5] saved scaling
    !---------------------------------------------------------
    nullify( ftsl%row_scaling )  
    nullify( ftsl%col_scaling )  


    !-------------------------------------------------------------------
    ! [3] Finish
    !-------------------------------------------------------------------
    
    ! force the success of the routine
    info          = 0
    ftsl%iinfo(1) = 0


  end Subroutine XFTS_ftsolver_Nullify

! [+] routine : XFTS_ftsolver_Free  ---------------------------------------------------
!
!> Free a ftsolver_instance
!! 
!! @param[in,out ] ftsl     the ftsolver instance to nullify
!!
!! @author Yohan Lee-tin-yien
!! @todo   Implement
Subroutine XFTS_ftsolver_Free    & ! intents
     (ftsl ,           & ! inout
     info )              ! out

  Use XFTS_ftsolver_type
  Use FTS_ftsolver_enum
  Use XFTS_sparse_matrix_mod, Only : &
       XFTS_sm_nullify,   &
       XFTS_sm_free
  Use XFTS_dense_matrix_mod, Only : &
       XFTS_dm_nullify,   &
       XFTS_dm_free




  Implicit None
  Include 'mpif.h'

  !* Subroutine arguments *!
  Type(XFTS_ftsolver_t), Intent(inout) :: ftsl
  Integer       , Intent(  out) :: info

  !* Local variables *!

  ! others
  Integer                       :: iinfo

  !- End of header -------------------------------------------------------------


  !-----------------------------------------------------------------------------
  ! [-] USER PARAMETERS (nothing to do)
  !-----------------------------------------------------------------------------

  !-----------------------------------------------------------------------------
  ! [-] Internal working data
  !-----------------------------------------------------------------------------

  !-------------------------------------------------------------------
  ! [--]  Description of the Non-overlapping domain Decomposition 
  !-------------------------------------------------------------------

  !---------------------------------------------------------
  ! [---] Local domain description (interface + interior)
  !---------------------------------------------------------
  !@ Type(XFTS_ftsolver_domain_t)     :: lc_domain 
  

  If (Associated(ftsl%lc_domain%index_local2global))&
       Deallocate( ftsl%lc_domain%index_local2global)

  If (Associated(ftsl%lc_domain%index_neig))&
       Deallocate( ftsl%lc_domain%index_neig)

  If (Associated(ftsl%lc_domain%ptr_direct_neig_edge))&
       Deallocate( ftsl%lc_domain%ptr_direct_neig_edge)

  If (Associated(ftsl%lc_domain%ptr_vertperneig))&
       Deallocate( ftsl%lc_domain%ptr_vertperneig)

  If (Associated(ftsl%lc_domain%index_vertperneig))&
      Deallocate( ftsl%lc_domain%index_vertperneig)
      
  If (Associated(ftsl%lc_domain%ptr_neig_edge))&
       Deallocate( ftsl%lc_domain%ptr_neig_edge )

  If (Associated(ftsl%lc_domain%index_neig_edge))&
       Deallocate( ftsl%lc_domain%index_neig_edge)

  If (Associated(ftsl%lc_domain%index_sec_neig))&
       Deallocate( ftsl%lc_domain%index_sec_neig )

  If (Associated(ftsl%lc_domain%ptr_sec_neig_edge))&
       Deallocate( ftsl%lc_domain%ptr_sec_neig_edge)

  If (Associated(ftsl%lc_domain%index_sec_neig_edge))&
       Deallocate( ftsl%lc_domain%index_sec_neig_edge)

  If (Associated(ftsl%lc_domain%index_secin_neig))&
       Deallocate( ftsl%lc_domain%index_secin_neig )

  If (Associated(ftsl%lc_domain%ptr_secin_neig_edge))&
       Deallocate( ftsl%lc_domain%ptr_secin_neig_edge )

  If (Associated(ftsl%lc_domain%index_secin_neig_edge))&
       Deallocate( ftsl%lc_domain%index_secin_neig_edge )


  !---------------------------------------------------------
  ! [---] Blocs on the matrix of the domain 
  !---------------------------------------------------------
  Call XFTS_sm_free(ftsl%sm_blockrow_delta1, iinfo)
  Call XFTS_sm_free(ftsl%sm_blockrow, iinfo)
  Call XFTS_sls_free(ftsl%sm_precond,iinfo)
  !Call XFTS_sm_free(ftsl%sm_lsi,iinfo)
  Call XFTS_sm_free(ftsl%sm_li,iinfo)
  !coucou
  !---------------------------------------------------------
  ! [---] saved scaling
  !---------------------------------------------------------


  If (Associated( ftsl%row_scaling )) Deallocate( ftsl%row_scaling )  
  If (Associated( ftsl%col_scaling )) Deallocate( ftsl%col_scaling )  

  !------------------------------------------------------------------
  ! [-] Finish
  !------------------------------------------------------------------

  ! This routine always return SUCCESS ( Do not check for errors )
  info = 0

End Subroutine XFTS_ftsolver_free




    ! [+] routine : XFTS_ftsolver_PretreatInputMatrix ---------------------------------
    !
    !> Pretreat the input matrix
    !!
    !!  Pretreat the matrix given by the user.
    !!  
    !!
    !!  @param[in,out] ftsl   the ftsolver instance
    !!  @param[   in ] master the MPI rank where gb_A is defined
    !!  @param[   in ] rank   the MPI rank
    !!  @param[   out] gb_A   the pretreated input matrix 
    !!                        only relevant on "master"
    !!  @param[   out] info the routine status
    !!
    !!  @author Yohan Lee-tin-yien
    !!
    Subroutine XFTS_ftsolver_PretreatInputMatrix &
         (ftsl, master, rank, gb_A, info )

      !* Modules *!
      Use FTS_ftsolver_enum
      Use FTS_mem_mod
      Use XFTS_ftsolver_type
      Use XFTS_sparse_matrix_mod
      Implicit None
      Include "mpif.h"

      !* Arguments *!
      Type(XFTS_ftsolver_t)       , Intent(inout) :: ftsl
      Integer              , Intent(  in) :: master
      Integer              , Intent(  in) :: rank
      Type(XFTS_sparse_matrix_t), Intent(  out) :: gb_A
      Integer              , Intent(  out) :: info

      !* Local variables *!
      Real(kind=8) :: time_pretreatInputMatrix
      Real(kind=8) :: time_copyAssemble
      Real(kind=8) :: time_symStruct

      !- End of header ---------------------------------------------------------

      ! [-] Init

      info = 0
      Call XFTS_sm_nullify(gb_A, info)
      FTS_ONFAILURE_RETURN(info)

      ! [--] Exit early
      If (rank /= master) Return

      time_pretreatInputMatrix = MPI_Wtime()

      ! [-] Copy inputs into gb_A
      !Write(*,*) "Warning:", Trim(FLNAME), __LINE__
      Call XFTS_sm_CreateFromData(&
           gb_A,SM_FMT_IJV,ftsl%sym,&
           ftsl%n, ftsl%n, ftsl%nnz, &
           ftsl%rows, ftsl%cols, ftsl%values,&
           info )
      FTS_ONFAILURE_RETURN(info)

      ! [-] On symmetric matrices transpose the upper triangle
      !Write(*,*) "Warning:", Trim(FLNAME), __LINE__
      If ( (gb_A%sym == SM_SYM_IsSPD) .or. &
           (gb_A%sym == SM_SYM_IsSymmetric))Then
         Call XFTS_sm_transposeUpper(gb_A)
      End If
      
      ! [-] Assemble matrix
      !Write(*,*) "Warning:", Trim(FLNAME), __LINE__
      time_copyAssemble = MPI_Wtime()
      Call XFTS_sm_Assemble(gb_A, gb_A%n, info)
      CHCKASSRT( info >= 0, info )
      If (info < 0) Return
      time_copyAssemble = MPI_Wtime() - time_copyAssemble
      
      ! [-] Symetrize matrix structure on General matrices
      !Write(*,*) "Warning:", Trim(FLNAME), __LINE__
      If ( gb_A%sym == SM_SYM_IsGeneral )Then
         
         time_symStruct = MPI_Wtime()
         Call XFTS_sm_symStruct(gb_A, info)
         CHCKASSRT( info >= 0, info )
         If (info < 0) Return
         time_symStruct = MPI_Wtime() - time_symStruct
         
      End If
      
      ! [-] Convert the matrix into IJV format
      !Write(*,*) "Warning:", Trim(FLNAME), __LINE__
      Call XFTS_sm_convert( gb_A, SM_FMT_IJV, info )
      CHCKASSRT( info >= 0, info )
      If (info < 0) Return

      ! [-] write preprocessed matrix
      !Write(*,*) "Warning:", Trim(FLNAME), __LINE__
      If ( Trim(ftsl%write_matrix) /= "" )Then
         
         Open(UNIT=111,FILE=ftsl%write_matrix,&
              ACTION="WRITE",STATUS="REPLACE", IOSTAT=info)
         CHCKASSRT( info >= 0, info )
         If (info < 0) Return
         
         Call XFTS_sm_mmwrite(gb_A, 111, info )
         CHCKASSRT( info >= 0, info )
         If (info < 0) Return
         
      End If

      time_pretreatInputMatrix = MPI_Wtime() - time_pretreatInputMatrix

      ! Save statistics

      ftsl%iinfog(IINFOG_MAT_ORDER) = gb_A%m
      ftsl%iinfog(IINFOG_MAT_NBENTRIES) = gb_A%nnz
      Call FTS_mem_add2usage(ftsl%mem,XFTS_sm_sizeof( gb_A ))

    End Subroutine XFTS_ftsolver_PretreatInputMatrix


  ! [+] routine : XFTS_ftsolver_Set_default_icntl  -------------------------------
  !
  !> initialize controls (integers)
  !!
  !! give default values to ICNTLs
  !!
  !! @param[in,out ] icntl   integer controls
  !! @param[   out ] info    routine status
  !!
  !! @author Yohan Lee-tin-yien
  !! @todo Complete 
  Subroutine XFTS_ftsolver_Set_default_icntl(icntl,info)

    Use FTS_ftsolver_enum ! only ICNTL_*

    Implicit None

    !* Subroutine arguments *!

    Integer, Intent(inout) :: icntl( FTSOLVER_ICNTL_SIZE )
    Integer, Intent(  out) :: info

    !* Local Variables *!

    ! Constants
    Integer, Parameter :: ZERO = 0

    ! Scalars
    Integer :: iinfo
    Integer :: i


    !- End of header -----------------------------------------------------------

    !---------------------------------------------------------------------------
    ! [0] Initialize data
    !---------------------------------------------------------------------------
    !

    ! Initialize by default to ZERO
    Do i = 1, FTSOLVER_ICNTL_SIZE
       icntl(i) = ZERO
    End Do

    !---------------------------------------------------------------------------
    ! [1] Set static values
    !---------------------------------------------------------------------------
    !

    ! outputs
    icntl( ICNTL_OUTPUT_ErrUnit   ) = 0
    icntl( ICNTL_OUTPUT_WarnUnit  ) = 0
    icntl( ICNTL_OUTPUT_StdUnit   ) = 6
    icntl( ICNTL_OUTPUT_Verbosity ) = 3

    icntl( ICNTL_PRINT_Cntls      ) = PRINT_Not
    icntl( ICNTL_PRINT_Infos      ) = PRINT_Not

    ! shift
    icntl( ICNTL_SHIFT_Do ) = FTS_FALSE



    ! Iterative Solver
    icntl( ICNTL_ITS_OrtStrategy ) = 3
    icntl( ICNTL_ITS_InitGuess   ) = 0
    icntl( ICNTL_ITS_MaxIters    ) = 0
    icntl( ICNTL_ITS_ResStrategy ) = 0 
    icntl( ICNTL_ITS_Restart     ) = 0


    ! 2nd lvl of parallelism
    icntl( ICNTL_2LVLS_Bind          ) = -1
    icntl( ICNTL_2LVLS_NNodes        ) = 0
    icntl( ICNTL_2LVLS_NcoresPerNode ) = 0
    icntl( ICNTL_2LVLS_NThrdsPerProc ) = 0
    icntl( ICNTL_2LVLS_NProcs        ) = 0


    ! input system
    icntl( ICNTL_INSYSTEM) = INSYSTEM_isCentralized

    !---------------------------------------------------------------------------

    !---------------------------------------------------------------------------
    ! [3] End routine
    !---------------------------------------------------------------------------

    ! this routine do not return errors
    info = 0

  End Subroutine XFTS_ftsolver_Set_default_icntl



  End Module XFTS_ftsolver_mod

