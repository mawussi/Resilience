  ! Warning: XFTS_GENFILE_COMMENT
#include "fts_defs_f.h"
#include "fts_macros_f.h"

! [+] module : XFTS_ftsolver_aux_mod ---------------------------------------------
!
!> Auxialiary module for ftsolver
!!
!! Basically contains useful routines for ftsolver.
!!
Module XFTS_solve_mod

  Use fts_ftsolver_enum
      Use FTS_error_mod
      Implicit None

      Public  :: XFTS_InitTimers
      Public  :: XFTS_VectorScalarProduct
      Public  :: XFTS_MatrixVectorProduct
      Public  :: XFTS_PcdVectorProduct
      Public  :: XFTS_Exchange_contrib
      Public  :: XFTS_PcdUpdate
      Public  :: XFTS_FGMRES_solve
      Public  :: XFTS_Solve
    contains


  ! [+] routine : XFTS_InitTimers ----------------------------------
  !
  !> Initialize the timers of the iterative solver
  !!
  !! @param [in,out] ftsl structure containing the timers
  !!
  Subroutine XFTS_InitTimers(ftsl)

    !* Modules *! 

    Use XFTS_ftsolver_type
    Use FTS_ftsolver_enum
    Implicit None

    !* Arguments *!

    Type(XFTS_ftsolver_t)       , Intent(inout) :: ftsl

    !- End of header---------------------------------------------------------

    ftsl%rinfo(RINFO_TIMING_Solve)        = 0.d0
    ftsl%rinfo(RINFO_TIMING_Interp)       = 0.d0
    ftsl%rinfo(RINFO_TIMING_InterpRhs)    = 0.d0
    ftsl%rinfo(RINFO_TIMING_InterpSlv)    = 0.d0
    ftsl%rinfo(RINFO_TIMING_Restart)      = 0.d0
    ftsl%rinfo(RINFO_FLOP_LI)             = 0.d0
    ftsl%rinfo(RINFO_FLOP_LSI)            = 0.d0
    ftsl%rinfo(RINFO_M_QR)                = 0.d0
    ftsl%rinfo(RINFO_N_QR)                = 0.d0
  End Subroutine XFTS_InitTimers


    ! [+] routine : FTSOLVER_Solve -----------------------------------------------
    !
    !> Performs the solve step
    !!
    !! @param[in,out ] ftsl          the ftsolver instance
    !!
    !!
    Subroutine XFTS_Solve (ftsl)

      !* Module(s) used *!

      Use XFTS_ftsolver_type
      Use fts_log_mod
      Use XFTS_dense_matrix_mod
      Use XFTS_sparse_matrix_mod
      Use FTS_time_mod
      Implicit None

      !* Arguments *!
      Type(XFTS_ftsolver_t), Intent(inout) :: ftsl

      !* Local variables *!

      ! Scalars
      Integer                        :: iinfo  
      Integer                        :: master, rank

      ! Derived types
      Type(XFTS_dense_matrix_t)           :: bound_rhs
      Type(XFTS_dense_matrix_t)           :: bound_sol

      !- End of header----------------------------------------------------------

      !-------------------------------------------------------------------------
      ! [-] Init
      !-------------------------------------------------------------------------

      
      rank  =ftsl%ikeep(IKEEP_MPIRANK)
      master=ftsl%ikeep(IKEEP_HOSTRANK)

      Call XFTS_dm_Nullify( bound_rhs , iinfo)
      FTS_ONFAILURE_GOTO9999(iinfo)

      Call XFTS_dm_Nullify( bound_sol , iinfo)
      FTS_ONFAILURE_GOTO9999(iinfo)
      
      Call XFTS_InitTimers(ftsl)
      !-------------------------------------------------------------------------
      ! [-] Generate the right-hand-sides & the init guesses
      !-------------------------------------------------------------------------

      !---------------------------------------------------------------
      ! [--] Generate the domain rhs 
      !---------------------------------------------------------------

         Call XFTS_DM_Create &
              (bound_rhs, ftsl%lc_domain%size_domain, &
               1,ftsl%lc_domain%size_domain,iinfo)
         FTS_ONFAILURE_GOTO9999(iinfo)
         Call XFTS_xcopy(bound_rhs%v,ftsl%lc_rhs, bound_rhs%m)

         Call XFTS_DM_Create &
              (bound_sol, ftsl%lc_domain%size_domain, &
               1,ftsl%lc_domain%size_domain,iinfo)
         FTS_ONFAILURE_GOTO9999(iinfo)

         Call FTS_time_start(ftsl%rinfo (RINFO_TIMING_Solve))
         Call XFTS_FGMRES_Solve(ftsl,bound_rhs,bound_sol, iinfo)
         Call FTS_time_stop(ftsl%rinfo (RINFO_TIMING_Solve))
         FTS_ONFAILURE_GOTO9999(iinfo)


      !-------------------------------------------------------------------------
      ! [-] Finish
      !-------------------------------------------------------------------------
         
9999  Continue

      ! Free memory
      Call XFTS_dm_Free(bound_rhs, iinfo)
      Call XFTS_dm_Free(bound_sol, iinfo)

      ! report status
      ftsl%iinfo( IINFO_STATUS ) = iinfo

    end Subroutine XFTS_Solve





  ! [+] routine :  XFTS_FGMRES_solve -----------------------------------------
  !
  !> solve the interface linear system using FGMRES
  !!
  !! Solve the linear system (distributed) with FGMRES.
  !! where the system is :
  !!           S. x = rhs
  !!
  !!-----
  !!
  !! @param[in,out] ftsl  
  !!
  !!       the ftsolve instance which contains :
  !!       - the parameters of the Iterative solver
  !!       - Its Preconditioner (which have multiple representation)
  !!       - necessary data to communicate between processors.
  !!       - the informations : iinfo/rinfo [in,out]
  !!
  !! @param[in    ] bound_rhs  
  !!
  !!       it holds the local part of the right-hand-side "rhs"
  !!
  !! @param[in,out] bound_sol
  !!
  !!       - On input , it may hold the local part of the initial guess x0
  !!       - On output, it holds the local part of the solution x
  !!
  !! @param[   out] info 
  !!       
  !!       Integer describing the routine status
  !!
  !!-----
  !!
  !! @note
  !! See also the documentation of PackFGMRES at :
  !!   http://www.cerfacs.fr/algor/reports/2003/TR_PA_03_03.pdf
  !!
  !! @author Azzam Haidar
  !! @author Luc   Giraud
  !! @author Yohan Lee-tin-yien
  !! @author Mawussi Zounon 
  Subroutine XFTS_FGMRES_solve(ftsl, bound_rhs, bound_sol, info)

    !* Module used *!

    Use XFTS_ftsolver_type
    Use FTS_ftsolver_enum
    Use XFTS_dense_matrix_mod
    Use XFTS_fault_mod
    Use fts_log_mod
    Use fts_error_mod
    Use FTS_time_mod


    Use FTS_mem_mod
    Implicit None

    !* Arguments *!

    Type(XFTS_ftsolver_t)        , Intent(inout) :: ftsl
    Type(XFTS_dense_matrix_t)  , Intent(in   ) :: bound_rhs 
    Type(XFTS_dense_matrix_t)  , Intent(inout) :: bound_sol 
    Integer               , Intent(  out) :: info

    !* Local Variables *!

    ! Scalars
    Integer    :: restrt
    FTS_INT :: lwork
    FTS_INT :: global_ndof
    FTS_INT :: bound_ndof
    Integer    :: i
    Integer    :: fault_info, last_failed_it
    Integer    :: rank, master
    Integer    :: verb
    Logical    :: HasInitGuess
    Integer(kind=8) :: mem

    ! Arrays
    Integer, Target :: fgmirc   (9)
    Integer         :: fgmicntl (6)
    Real(Kind=XFTS_FLOATKIND)         :: fgmrcntl (3)

    Integer         :: fgmiinfo (8)
    Real(Kind=XFTS_FLOATKIND)         :: fgmrinfo     !lty: why 4 and not 1 ?
    
    XFTS_FLOAT, Pointer :: work     (:)

    ! Aliases (Pointers)

    ! to fgmirc (reverse communication)
    Integer, Pointer :: revcom
    Integer, Pointer :: colx
    Integer, Pointer :: coly
    Integer, Pointer :: colz
    Integer, Pointer :: nbscal

    ! to fgmres workspace "work"
    Type(XFTS_dense_matrix_t) :: x 
    Type(XFTS_dense_matrix_t) :: y 
    Type(XFTS_dense_matrix_t) :: z

    !- End of header -------------------------------------------------

    !-----------------------------------------------------------------
    ! [1] Initialize local variables / FGMRES package
    !-----------------------------------------------------------------

    !---------------------------------------------------------
    ! [1.1] Get differents sizes & options
    !---------------------------------------------------------

    global_ndof    = ftsl%n
    bound_ndof    = ftsl%lc_domain%size_domain
    HasInitGuess = .false. !(ftsl%IKEEP(IKEEP_ITS_InitGuess) == 1)

    !---------------------------------------------------------
    ! [1.2] Set FGMRES options
    !---------------------------------------------------------


    ! Get default options
    Call init_XFTS_ARITHfgmres(fgmicntl,fgmrcntl)

    ! Set specific options
    fgmrcntl(1) = Real( ftsl%rcntl(RCNTL_ITS_Tolerance),KIND=XFTS_FLOATKIND)

    fgmicntl(2) = 0
    fgmicntl(4) = ftsl%ICNTL(ICNTL_ITS_OrtStrategy )
    fgmicntl(5) = ftsl%ICNTL(ICNTL_ITS_InitGuess   )
    fgmicntl(6) = ftsl%ICNTL(ICNTL_ITS_MaxIters    )

    restrt     =  ftsl%ICNTL(ICNTL_ITS_Restart )

    ! Outputs

    rank= ftsl%ikeep(IKEEP_MPIRANK)
    master= ftsl%ikeep(IKEEP_HOSTRANK)
    Call fts_log_GetVerbosity(verb)
    If ( (verb >= MSG_STD).and.(rank == master) )Then
       ! print warning messages
       Call fts_log_GetUnit(MSG_WARNING, fgmicntl(2))
       ! print convergence history to stdout
       Call fts_log_GetUnit(MSG_STD    , fgmicntl(3))
    End If

    ! Create solution if no init guess
    If (HasInitGuess .eqv. .False.) Then
       Call XFTS_dm_Nullify(bound_sol,info) 
       FTS_ONFAILURE_GOTO9999(info)
       Call XFTS_dm_Create(bound_sol,bound_ndof,1,bound_ndof,info  ) 
       FTS_ONFAILURE_GOTO9999(info)
       Do i=1,bound_sol%m
          bound_sol%v(i) = XFTS_FLOATZERO
       End Do
    End If

    !---------------------------------------------------------
    ! [1.3] Allocate FGMRES workspace "work"
    !---------------------------------------------------------

    ! Set lwork value
    lwork = restrt**2 + restrt*(2*bound_ndof+5) + 5* bound_ndof + 1
     ! the workspace should be large enough to store the m dot-products 
    ! Allocate
    If ((fgmicntl(4).Eq.2).Or.(fgmicntl(4).Eq.3)) Then
       lwork = lwork + restrt
     Else
        lwork = lwork +1
     End If
    
    Allocate(work(lwork), STAT = info )
    CHCKALLOC(info)
    FTS_ONFAILURE_GOTO9999(info)
    
    !---------------------------------------------------------
    ! [1.4] Setup the data
    !---------------------------------------------------------
    ! Setup the rhs
    ! Work(1:bound_ndof)             = local initial guess
    ! Work(bound_ndof+1:2*bound_ndof) = local right-hand-side
    If (HasInitGuess) Call XFTS_ARITHcopy &
         (bound_ndof, bound_sol%v(1),1, work(1), 1)

    Call XFTS_ARITHcopy&
         (bound_ndof,bound_rhs%v(1),1,work(bound_ndof+1), 1)

    ! initialize vectors
    info = 0
    Call XFTS_dm_nullify(x,info)
    Call XFTS_dm_nullify(y,info)
    Call XFTS_dm_nullify(z,info)
    FTS_ONFAILURE_GOTO9999(info)

    x%ld = bound_ndof
    x%m  = bound_ndof
    x%n  = 1

    y%ld = bound_ndof
    y%m  = bound_ndof
    y%n  = 1

    z%ld = bound_ndof
    z%m  = bound_ndof
    z%n  = 1

    ! Associate aliases
    revcom => fgmirc(1)
    colx   => fgmirc(2)
    coly   => fgmirc(3)
    colz   => fgmirc(4)
    nbscal => fgmirc(5)

    !-----------------------------------------------------------------
    ! [2] Call the solver, Iterate until solution is found
    !-----------------------------------------------------------------
    fault_info = 0
!    fgmiinfo(2) = -1
    last_failed_it = -1
    fgmirc(8)=0
     fgmirc(9) = 0
    revcom = 1 

    Do While ((revcom > 0).And.(info >= 0))

       !----------------------------------------------------
       ! [2.1] Call the driver
       !----------------------------------------------------


10       Call Drive_XFTS_ARITHfgmres(         &
            global_ndof,bound_ndof,  &
            restrt,lwork,work,    &
            fgmirc,fgmicntl,fgmrcntl, &
            fgmiinfo,fgmrinfo)

       !set  fgmirc(9)  to 1 continue increasing the iterations
       fgmirc(9) = 1
       !----------------------------------------------------
       ! [2.2] Associate arrays
       !----------------------------------------------------

       y%v => work(coly : (coly+bound_ndof-1) )
       z%v => work(colz : (colz+bound_ndof-1) )

       If ( revcom /= 4)Then
          x%n   = 1
          x%v => work(colx : (colx+bound_ndof-1) )
       Else
          x%n   =  nbscal
          x%v => work(colx : (colx + bound_ndof*nbscal -1) )
       End If




       !----------------------------------------------------
       ! [2.2.1] Fault simulation
       !----------------------------------------------------
       !check if the iteration will be faulty
       !and set flag to restart fgmres 

       If((fgmirc(8) .Eq. 0) .And. (fault_info < 1 )) Then

          Call XFTS_IsFault(ftsl, fgmiinfo(2),fault_info)
          If(fault_info > 0) Then !multiple fault
             Call FTS_time_start(ftsl%rinfo(RINFO_TIMING_Restart))
             fgmirc(8) = 1
          End If
       End If

       ! ! Check if the fault can be taken into account, i.e., the current iterate has been                                                                                                              
       ! !computed and the internal variables have been reset to allow for a new call           

       If ((fgmirc(1).Eq.0).And.(fgmirc(8).Eq.1))Then ! 
          fgmicntl(5) = 1
          fgmirc(8) = 0
             Call FTS_time_stop(ftsl%rinfo (RINFO_TIMING_Restart))
          Call XFTS_Check_Fault(ftsl, work) 
          Goto 10
       End If



       !----------------------------------------------------
       ! [2.3] Do what the driver asked to perform.
       !----------------------------------------------------
       !
       ! 1     - Perform Matrix Vector Product : z <-- Matrix . x
       ! 2     - Apply left preconditioning    : z <-- Precond . x
       ! 3     - Apply right precondiotioning  : z <-- Precond . x
       ! 4     - Perform scalar product        : z <-- x * y
       ! 0     - Exit                          : Solving step ended
       ! other - Throw Error : unknown revcom
       !

       Select Case (revcom)
       Case (0); Continue                      ! Exit    
       Case (1); Call XFTS_MatrixVectorProduct (ftsl, x, z, info)
       Case (2); CHCKASSRT(.False., info )     ! Ftsolver does not support left preconditioning
       Case (3); Call XFTS_PcdVectorProduct   (ftsl, x,    z, info)
       Case (4); Call XFTS_VectorScalarProduct(ftsl, x, y, z, info)
       Case Default; CHCKASSRT(.False., info ) ! Undocumented revcom
       End Select

    End Do

    !---------------------------------------------------------------------------
    ! [-] Checking the outputs 
    !---------------------------------------------------------------------------    

    !-----------------------------------------------
    ! [--] Checking revcom, print a message on error
    !-----------------------------------------------

    Select Case (revcom)
    Case (0); Continue ! Gmres exited correctly
    Case (1); Call FTS_Log(MSG_ERROR,"Error occured in matrix/vector Product")
    Case (2); Call FTS_Log(MSG_ERROR,"Ftsolver does not support left preconditioning")
    Case (3); Call FTS_Log(MSG_ERROR,"Error occured in right preconditioning")
    Case (4); Call FTS_Log(MSG_ERROR,"Error occured in scalar product")
    Case Default; Call FTS_LogWithInfo(MSG_ERROR,revcom,&
            "From PackFGMRES: undocumented revcom code. revcom =")
    End Select
    FTS_ONFAILURE_GOTO9999(info)
    
    !----------------------------------------------------
    ! [--] Checking fgmiinfo(1), the return status of PackGMRES
    !----------------------------------------------------

    Select Case (fgmiinfo(1))
    Case ( 0) ! convergence achieved
       info = FTS_SUCCESS
       Continue
    Case (-4) ! non convergence. Consider it as a warning and not an error
       info = FTS_SUCCESS + 1
       !Call FTS_LogWithInfo(MSG_WARNING,fgmicntl(7),&
        !    "From PackFGMRES: convergence not achieved. Nb iterations:")
    Case Default
       Call FTS_LogWithInfo(MSG_ERROR,fgmiinfo(1),"PackFGRMES exited with error code =")
       CHCKASSRT(.False., info )
    End Select
    FTS_ONFAILURE_GOTO9999(info)

    !-----------------------------------------------------------------
    ! [-] Exit routine
    !-----------------------------------------------------------------

    ! Save data
    Call XFTS_ARITHcopy(bound_ndof, work(1),1, bound_sol%v(1),1 )

    ftsl%IINFOG( IINFOG_ITS_NbIters ) = fgmiinfo(2)
    ftsl%RINFOG( RINFOG_ITS_BckErr  ) = fgmrinfo

    mem = INT(lwork*XFTS_FLOATBYTESIZE,8)
    Call FTS_mem_add2usage(ftsl%mem, mem )
    ftsl%IINFO ( IINFO_ITS_MEMPEAK  ) = byte2Mbyte(mem)
    ftsl%IINFO ( IINFO_ITS_MEMUSED  ) = byte2Mbyte(mem)

    ! Free memory
9999 Continue
    If (Associated(work)) Deallocate( work )

  End Subroutine XFTS_FGMRES_Solve



  ! [+] routine :  XFTS_VectorScalarProduct ----------------------------------
  !
  !> Perform the operation  : z <-- x . y 
  !!
  !! Perform the Matrix Vector Product, with :
  !!  - x       : the first  vector          (distributed)
  !!  - y       : the second vector          (distributed)
  !!  - z       : the result                 (distributed)
  !!
  !! Where all vectors x,y,z are vectors defined on the interface.
  !!
  !! @param[in,out] ftsl  
  !!
  !!       the structure containing : 
  !!       - the data necessary to propagate the local results.
  !!       - the temporary buffer on the interface "intrfbuff"
  !!       - the weight to apply on "lc_y"
  !!
  !! @param[in,out] lc_x      
  !!
  !!       the local part of "x" in the matrix product.
  !!
  !! @param[in,out] lc_y      
  !!
  !!       the local part of "y" in the matrix product.
  !!
  !! @param[in,out] lc_z      
  !!
  !!       the local part of "z" the matrix product.
  !!
  !! @param[in,out] info   
  !!         
  !!       the routine status.
  !!
  !!
  Subroutine XFTS_VectorScalarProduct(ftsl, lc_x, lc_y, lc_z, info)

    !* Module used *!

    Use XFTS_ftsolver_type     , Only : &
         XFTS_ftsolver_t                   ! type
    Use XFTS_dense_matrix_mod, Only : &
         XFTS_dense_matrix_t
    Use FTS_ftsolver_enum
    Implicit None
    Include "mpif.h"

    !* Subroutine Arguments *!

    Type(XFTS_ftsolver_t)        , intent(inout) :: ftsl
    Type(XFTS_dense_matrix_t)  , intent(in   ) :: lc_x 
    Type(XFTS_dense_matrix_t)  , intent(in   ) :: lc_y
    Type(XFTS_dense_matrix_t)  , intent(inout) :: lc_z 
    Integer               , intent(  out) :: info

    !* Local Variables *!

    ! scalar
    Integer    :: iinfo
    Integer    :: comm
    Integer    :: nb_scalar_products
    FTS_INT :: i
    FTS_INT :: nrows
    Real(kind=8) :: timing

    ! Strings
    Character*1                  :: trans 

    ! Arrays
    XFTS_FLOAT, Pointer :: orthvect (:)

    !- End of header -------------------------------------------------

    !-----------------------------------------------------------------
    ! [1] Initialize local variables
    !-----------------------------------------------------------------

    info = 0
    timing = MPI_WTime()

    ! Associate input data
    comm               =  ftsl%comm
    nb_scalar_products =  lc_x%n
    nrows              =  lc_y%m

    ! Allocate workspace
    Nullify( orthvect) 
    Allocate(&
         orthvect(nb_scalar_products),&
          STAT= info )
    CHCKASSRT( info == 0, info )
    If ( info < 0) Goto 9999


    !-----------------------------------------------------------------
    ! [3] Perform the scalar product (result = orthvect)
    !-----------------------------------------------------------------
    !
    ! Here there is multiple vector/vector product to do (multiple "x")
    ! So instead of using the BLAS "xdot", we use "xmv" 
    ! with "x^T" as the matrix.
    !
    ! Temporary result is saved in orthvect
    !

    trans = 'C'
    Call XFTS_ARITHgemv  (   & !                              
         trans,lc_x%m,       & !   orthvect(i) = x(:,i)^T . y 
         nb_scalar_products, & !                              
         XFTS_FLOATONE,      & ! 
         lc_x%v(1),lc_x%ld,  & !                              
         lc_y%v(1) ,1,   & ! With                         
         XFTS_FLOATZERO,     & !   i = 1 .. nb_scalar_products
         orthvect(1), 1      & !
         )

    !-----------------------------------------------------------------
    ! [4] Compute "lc_z" , sum of "orthvect" on all processor
    !-----------------------------------------------------------------

    Call MPI_AllReduce(      & ! 
         orthvect,lc_z%v,    & ! lc_z(i) = SUM( orthvect(i) ) 
         nb_scalar_products, & !            
         XFTS_FLOATMPI,      & ! With 
         MPI_SUM,comm,iinfo  & !  i = 1 .. nb_scalar_products
         )
    ASSRT( iinfo == MPI_SUCCESS )

    !-----------------------------------------------------------------
    ! [5] Exit routine
    !-----------------------------------------------------------------

9999 Continue

    ! Free buffers
    If (Associated(orthvect)) Deallocate(orthvect)

  End Subroutine XFTS_VectorScalarProduct



  ! [+] routine :  XFTS_MatrixVectorProduct ------------------------------------
  !
  !> Perform Matrix Vector Product : z <-- Matrix . x
  !!
  !! Perform Matrix Vector Product : z <-- Matrix . x
  !!  
  !! @param[in,out] ftsl  

  !!
  !! @param[in,out] x      
  !!
  !!       the local part of the vector in the matrix product.
  !!
  !! @param[in,out] z      
  !!
  !!       the local part of result of the matrix product.
  !!
  !! @param[in,out] info   
  !!         
  !!       the routine status.
  !!
  !! @todo Maybe replace direct calls to blas, etc. to generic ones.
  !! @todo other strategies
  !!
  !! @author Yohan Lee-tin-yien
  !!
  Subroutine XFTS_MatrixVectorProduct(ftsl, x, z, info)

    !* Module used *!

    Use FTS_ftsolver_enum
    Use XFTS_ftsolver_type
    Use XFTS_dense_matrix_mod

    Implicit None
    include 'mpif.h'    
    !* Subroutine Arguments *!

    Type(XFTS_ftsolver_t)        , intent(inout) :: ftsl
    Type(XFTS_dense_matrix_t)  , intent(in   ) :: x 
    Type(XFTS_dense_matrix_t)  , intent(inout) :: z 
    Integer               , intent(  out) :: info

    !* Local variables * !
    Real(kind=8) :: timing
!    Real(kind=8), External :: MPI_Wtime

    !- End of header -------------------------------------------------



    ! Perform local matrix vector product

    ! Scalars
    Integer :: tagg      !< MPI tagg
    Integer :: comm, rank      !< MPI Communicator, rank
    Integer :: iinfo     !< status of the called routines
    Integer :: i,k     !< dummy iterators
    Integer :: nb_neig      !< field of lc_intrf. see lc_intrf%mynbvi
    Integer :: neigh     !< neighbor
    Integer :: size_vertperneig         !> size of  index_vertperneig
    Integer :: size_neig_edge, neig_col
    Integer :: nloc, col , row, iter_col, iter_row 
    Integer :: mpistatus(MPI_STATUS_SIZE) 
    FTS_INT :: sbuffsize !< size of sbuff(:)
    FTS_INT :: rbuffsize !< size of rbuff(:)
    FTS_INT :: start 
    FTS_INT :: end
 
    ! Arrays

    Integer   , Pointer :: mpirreq(:)               !< MPI requests receivers
    Integer   , Pointer :: mpisreq(:)               !< MPI requests senders
    XFTS_FLOAT, Pointer :: sbuff (:)             !< MPI send Buffer
    XFTS_FLOAT, Pointer :: rbuff (:)             !< MPI receive Buffer
    Integer, Pointer    :: index_neig          (:)  !< field of boundary
    Integer, pointer :: ptr_direct_neig_edge (:)        !> direct pointer on index_neig_edge
    Integer, pointer :: ptr_vertperneig (:)
    Integer, pointer :: index_vertperneig (:)
    Integer, pointer :: ptr_neig_edge      (:) 
    Integer, pointer :: index_neig_edge    (:)

    


    !- End of header -----------------------------------------------------------

    !---------------------------------------------------------------------------
    ! [1] Init
    !---------------------------------------------------------------------------
    
    ! status
    info = 0
    timing = MPI_WTime()
    rank   = ftsl%ikeep(IKEEP_MPIRANK)
    ! mpi
    tagg             =  44
    comm             =  ftsl%comm


    
    ! interface description
    nb_neig             =  ftsl%lc_domain%nb_neig
    size_vertperneig    =  ftsl%lc_domain%size_vertperneig
    nloc                =  ftsl%lc_domain%size_domain
    index_neig          => ftsl%lc_domain%index_neig
    ptr_direct_neig_edge => ftsl%lc_domain%ptr_direct_neig_edge
    ptr_vertperneig     => ftsl%lc_domain%ptr_vertperneig
    index_vertperneig   => ftsl%lc_domain%index_vertperneig
    ptr_neig_edge       => ftsl%lc_domain%ptr_neig_edge
    index_neig_edge     => ftsl%lc_domain%index_neig_edge
    


    ! values


    !---------------------------------------------------------------------------
    ! [1] Allocate the buffers
    !---------------------------------------------------------------------------
    
    Nullify( sbuff, rbuff )
    Nullify( mpirreq, mpisreq )

    Allocate( mpirreq(nb_neig), mpisreq(nb_neig), STAT= info )
    CHCKASSRT( info == 0, info )
    If (info < 0) Goto 9999

    sbuffsize = ptr_vertperneig(nb_neig+1)- ptr_vertperneig(1)
    rbuffsize = ptr_direct_neig_edge(nb_neig+1)- ptr_direct_neig_edge(1)
    Allocate( sbuff(sbuffsize), rbuff(rbuffsize), STAT= info )
    CHCKASSRT( info == 0, info )
    If (info < 0) Goto 9999

    Do i=1,sbuffsize 
       sbuff(i) = XFTS_FLOATZERO
       sbuff(i) = -999999999
    End Do

    Do i=1,rbuffsize 
       rbuff(i) = XFTS_FLOATZERO
    End Do


    !---------------------------------------------------------------------------
    ! [2] Send the computed value on the interface
    !---------------------------------------------------------------------------

    Do i =1,nb_neig

       ! Pack 
       start = ptr_vertperneig(i)
       end   = ptr_vertperneig(i+1) - 1
       Do k = start, end 
          sbuff(k) = x%v( index_vertperneig(k) ) 
       End Do

       ! Send
       Call MPI_Isend (sbuff(start),(end-start+1),&
            XFTS_FLOATMPI,index_neig(i),tagg,comm,mpisreq(i),iinfo)
       ASSRT( iinfo  == MPI_SUCCESS )
       if (rank .EQ. 3) Then
       End if
    End Do



    !---------------------------------------------------------------------------
    ! [3] Start receiving 
    !---------------------------------------------------------------------------

    Do i=1,nb_neig
       
       start = ptr_direct_neig_edge(i)
       end   = ptr_direct_neig_edge(i+1) - 1
       Call MPI_Irecv( rbuff(start),(end-start+1), &
            XFTS_FLOATMPI,index_neig(i), tagg,comm,mpirreq(i),iinfo)
       ASSRT( iinfo  == MPI_SUCCESS )

    End Do



    !---------------------------------------------------
    ! compute local product
    !--------------------------------------------------
    Do i=1, nloc
       z%v(i) = XFTS_FLOATZERO
    End Do
    

    Do iter_col=1,nloc
       DO k= ftsl%sm_blockrow%csc(iter_col), ftsl%sm_blockrow%csc(iter_col+1)-1
          row = ftsl%sm_blockrow%i(k)
          z%v(row) = z%v(row) + ftsl%sm_blockrow%v(k)*x%v(iter_col)
       End Do
    End Do


    !---------------------------------------------------------------------------
    ! [3] Unpack values from neigbhors
    !---------------------------------------------------------------------------

    Do i=1,nb_neig

       Call MPI_WaitAny(nb_neig, mpirreq, neigh, MPI_STATUS_IGNORE, iinfo )
       ASSRT( iinfo  == MPI_SUCCESS )
       start = ptr_direct_neig_edge(neigh)
       end   = ptr_direct_neig_edge(neigh+1) - 1
       Do neig_col= start, end
          col = neig_col+ nloc 
          Do k= ftsl%sm_blockrow%csc(col), ftsl%sm_blockrow%csc(col+1)-1
             row = ftsl%sm_blockrow%i(k)
             z%v(row) = z%v(row) + ftsl%sm_blockrow%v(k)*rbuff(neig_col)
          End Do
       End Do
       
    End Do
    Call MPI_WaitAll(nb_neig, mpisreq, MPI_STATUSES_IGNORE, iinfo )
    ASSRT( iinfo  == MPI_SUCCESS )

    ! --------------------------------------------------------------------------
    ! [3] Exit routine
    ! --------------------------------------------------------------------------

9999 Continue
    If (Associated(sbuff)) Deallocate(sbuff)
    If (Associated(rbuff)) Deallocate(rbuff)
    If (Associated(mpirreq)) Deallocate(mpirreq)
    If (Associated(mpisreq)) Deallocate(mpisreq)


  End Subroutine XFTS_MatrixVectorProduct


  ! [+] routine :  XFTS_PcdVectorProduct ------------------------------------
  !
  !> Apply the preconditioning  : z <-- Precond . x
  !!
  !! Perform the Matrix Vector Product, with :
  !!  - Precond : the preconditioning matrix (usually distributed)
  !!  - x       : the vector                 (usually distributed)
  !!  - z       : the result                 (usually distributed)
  !!
  !! On error this routine abort the system, to prevent MPI locking.
  !! 
  !! @param[in,out] ftsl  
  !!
  !!       the structure containing : 
  !!       - the Preconditioner in different formats
  !!       - the data necessary to propagate the local results.
  !!
  !! @param[in,out] lc_x      
  !!
  !!       the (local part of the) vector in the matrix product.
  !!
  !! @param[in,out] lc_z      
  !!
  !!       the (local part of the) result of the matrix product.
  !!
  !! @param[in,out] info   
  !!         
  !!       the routine status.
  !!
  !!
  Subroutine XFTS_PcdVectorProduct(ftsl, lc_x, lc_z, info)

    !* Module used *!

    Use XFTS_ftsolver_type, Only :         &
         XFTS_ftsolver_t                     ! type(s)
    Use XFTS_dense_matrix_mod
    Use FTS_ftsolver_enum
    Use XFTS_mumps_mod
    Use fts_log_mod
    Implicit None


    !* Arguments *!

    Type(XFTS_ftsolver_t)        , intent(inout) :: ftsl
    Type(XFTS_dense_matrix_t)  , intent(in   ) :: lc_x 
    Type(XFTS_dense_matrix_t)  , intent(inout) :: lc_z 
    Integer               , intent(  out) :: info

    !* Local Variables *!
    Type(XFTS_dense_matrix_t)  :: AS_z


    ! scalar
    Integer :: Pcd_Strategy, size_Pcd
    Real(kind=8) :: timing
    Real(kind=8), External :: MPI_Wtime

    !- End of header -------------------------------------------------

    !-----------------------------------------------------------------
    ! [1] Initialize local variables
    !-----------------------------------------------------------------

    
    size_Pcd = ftsl%lc_domain%size_domain_delta1
    timing = MPI_WTime()
    info = 0


    Call XFTS_DM_Create (AS_z, size_Pcd, 1,size_Pcd,info)
    
    ! Copy lc_x into AS_z
    Call XFTS_ARITHcopy(lc_x%m, lc_x%v, 1,AS_z%v,1)

    !-----------------------------------------------!
    !  Exchange contribution. Each process          !
    !  send and receive additional entries (detla1) !
    !  required to apply the AS preconditioner      !
    !-----------------------------------------------!

    Call XFTS_Exchange_contrib(ftsl, lc_x, AS_z, info) 
    

    ! Perform local product : lc_z <-- Pcd. lc_z
     Call XFTS_mumps_solve_RHS (&
     ftsl%sm_precond%sds,&
          AS_z%v,AS_z%n,AS_z%ld,&
          info)
     CHCKRINFO(info)
     If (info < 0) Goto 9999

     ! Copy lc_x into AS_z
     Call XFTS_ARITHcopy(lc_z%m, AS_z%v, 1,lc_z%v,1)


        
9999 Continue

    timing = MPI_WTime() - timing
    ASSRT( info >= 0 )
    Call XFTS_DM_Free(AS_z, info)

  End Subroutine XFTS_PcdVectorProduct



  ! [+] routine :  XFTS_Exchange_contrib ------------------------------------
  !
  !> 
  !!
  !!  Receive contribution from neighbor to apply 
  !!  the AS preconditioner
  !! @param[in,out] ftsl  

  !!
  !! @param[in,out] x      
  !!
  !!       the local part of the vector 
  !!
  !! @param[in,out] z      
  !!
  !!       the local part of x + entries of neighbors
  !!
  !! @param[in,out] info   
  !!         
  !!       the routine status.
  !!
  !! @todo Maybe replace direct calls to blas, etc. to generic ones.
  !! @todo other strategies
  !!
  !! @author Mawussi Zounon
  !!
  Subroutine XFTS_Exchange_contrib(ftsl, x, z, info)

    !* Module used *!

    Use FTS_ftsolver_enum
    Use XFTS_ftsolver_type
    Use XFTS_dense_matrix_mod

    Implicit None
    include 'mpif.h'    
    !* Subroutine Arguments *!

    Type(XFTS_ftsolver_t)        , intent(inout) :: ftsl
    Type(XFTS_dense_matrix_t)  , intent(in   ) :: x 
    Type(XFTS_dense_matrix_t)  , intent(inout) :: z 
    Integer               , intent(  out) :: info

    !* Local variables * !
    Real(kind=8) :: timing
!    Real(kind=8), External :: MPI_Wtime

    !- End of header -------------------------------------------------



    ! Perform local matrix vector product

    ! Scalars
    Integer :: tagg      !< MPI tagg
    Integer :: comm, rank      !< MPI Communicator, rank
    Integer :: iinfo     !< status of the called routines
    Integer :: i,k     !< dummy iterators
    Integer :: nb_neig      !< field of lc_intrf. see lc_intrf%mynbvi
    Integer :: neigh     !< neighbor
    Integer :: size_vertperneig         !> size of  index_vertperneig
    Integer :: size_neig_edge, neig_col
    Integer :: nloc, col , row, iter_col, iter_row 
    Integer :: mpistatus(MPI_STATUS_SIZE) 
    FTS_INT :: sbuffsize !< size of sbuff(:)
    FTS_INT :: rbuffsize !< size of rbuff(:)
    FTS_INT :: start 
    FTS_INT :: end
 
    ! Arrays

    Integer   , Pointer :: mpirreq(:)               !< MPI requests receivers
    Integer   , Pointer :: mpisreq(:)               !< MPI requests senders
    XFTS_FLOAT, Pointer :: sbuff (:)             !< MPI send Buffer
    Integer, Pointer    :: index_neig          (:)  !< field of boundary
    Integer, pointer :: ptr_direct_neig_edge (:)        !> direct pointer on index_neig_edge
    Integer, pointer :: ptr_vertperneig (:)
    Integer, pointer :: index_vertperneig (:)
    Integer, pointer :: ptr_neig_edge      (:) 
    Integer, pointer :: index_neig_edge    (:)

    


    !- End of header -----------------------------------------------------------

    !---------------------------------------------------------------------------
    ! [1] Init
    !---------------------------------------------------------------------------
    
    ! status
    info = 0
    timing = MPI_WTime()
    rank   = ftsl%ikeep(IKEEP_MPIRANK)
    ! mpi
    tagg             =  44
    comm             =  ftsl%comm


    
    ! interface description
    nb_neig             =  ftsl%lc_domain%nb_neig
    size_vertperneig    =  ftsl%lc_domain%size_vertperneig
    nloc                =  ftsl%lc_domain%size_domain
    index_neig          => ftsl%lc_domain%index_neig
    ptr_direct_neig_edge => ftsl%lc_domain%ptr_direct_neig_edge
    ptr_vertperneig     => ftsl%lc_domain%ptr_vertperneig
    index_vertperneig   => ftsl%lc_domain%index_vertperneig

    

    ! values


    !---------------------------------------------------------------------------
    ! [1] Allocate the buffers
    !---------------------------------------------------------------------------
    
    Nullify( sbuff)
    Nullify( mpirreq, mpisreq )

    Allocate( mpirreq(nb_neig), mpisreq(nb_neig), STAT= info )
    CHCKASSRT( info == 0, info )
    If (info < 0) Goto 9999

    sbuffsize = ptr_vertperneig(nb_neig+1)- ptr_vertperneig(1)

    Allocate( sbuff(sbuffsize), STAT= info )
    CHCKASSRT( info == 0, info )
    If (info < 0) Goto 9999

    Do i=1,sbuffsize 
       sbuff(i) = XFTS_FLOATZERO
    End Do



    !---------------------------------------------------------------------------
    ! [2] Send the computed value on the interface
    !---------------------------------------------------------------------------

    Do i =1,nb_neig

       ! Pack 
       start = ptr_vertperneig(i)
       end   = ptr_vertperneig(i+1) - 1
       Do k = start, end 
          sbuff(k) = x%v( index_vertperneig(k) ) 
       End Do

       ! Send
       Call MPI_Isend (sbuff(start),(end-start+1),&
            XFTS_FLOATMPI,index_neig(i),tagg,comm,mpisreq(i),iinfo)
       ASSRT( iinfo  == MPI_SUCCESS )
       if (rank .EQ. 3) Then
       End if
    End Do



    !---------------------------------------------------------------------------
    ! [3] Start receiving 
    !---------------------------------------------------------------------------

    Do i=1,nb_neig
       
       start = nloc + ptr_direct_neig_edge(i)
       end   = nloc + ptr_direct_neig_edge(i+1) - 1
       Call MPI_Irecv( z%v(start),(end-start+1), &
            XFTS_FLOATMPI,index_neig(i), tagg,comm,mpirreq(i),iinfo)
       ASSRT( iinfo  == MPI_SUCCESS )

    End Do


    ! Wait
    Call MPI_WaitAll(nb_neig, mpirreq, MPI_STATUS_IGNORE, iinfo )
    ASSRT( iinfo  == MPI_SUCCESS )

    Call MPI_WaitAll(nb_neig, mpisreq, MPI_STATUSES_IGNORE, iinfo )
    ASSRT( iinfo  == MPI_SUCCESS )

    ! --------------------------------------------------------------------------
    ! [3] Exit routine
    ! --------------------------------------------------------------------------

9999 Continue
    If (Associated(sbuff)) Deallocate(sbuff)
    If (Associated(mpirreq)) Deallocate(mpirreq)
    If (Associated(mpisreq)) Deallocate(mpisreq)

    timing = MPI_WTime() - timing

  End Subroutine XFTS_Exchange_contrib


  ! [+] routine :  XFTS_PcdUpdate ------------------------------------
  !
  !> Update z after the application of the preconditionner
  !!
  !! Perform Matrix Vector Product : z_loc <--z_loc + contrib_of_Neighbors
  !!  
  !! @param[in,out] ftsl  

  !!
  !! @param[in,out] z      
  !!
  !!       the local part of the vector  extended with entries of delta1
  !!
  !!
  !!  !!
  !! @param[in,out] info   
  !!         
  !!       the routine status.
  !!
  !!Mawussi Zounon
  !!
  Subroutine XFTS_PcdUpdate(ftsl, z, info)

    !* Module used *!

    Use FTS_ftsolver_enum
    Use XFTS_ftsolver_type
    Use XFTS_dense_matrix_mod

    Implicit None
    include 'mpif.h'    
    !* Subroutine Arguments *!

    Type(XFTS_ftsolver_t)        , intent(inout) :: ftsl
    Type(XFTS_dense_matrix_t)  , intent(inout) :: z 
    Integer               , intent(  out) :: info

    !* Local variables * !
    Real(kind=8) :: timing
!    Real(kind=8), External :: MPI_Wtime

    !- End of header -------------------------------------------------



    ! Perform local matrix vector product

    ! Scalars
    Integer :: tagg      !< MPI tagg
    Integer :: comm, rank      !< MPI Communicator, rank
    Integer :: iinfo     !< status of the called routines
    Integer :: i,k     !< dummy iterators
    Integer :: nb_neig      !< field of lc_intrf. see lc_intrf%mynbvi
    Integer :: neigh     !< neighbor
    Integer :: size_vertperneig         !> size of  index_vertperneig
    Integer :: size_neig_edge, neig_col
    Integer :: nloc, col , row, iter_col, iter_row 
    Integer :: mpistatus(MPI_STATUS_SIZE) 
    FTS_INT :: sbuffsize !< size of sbuff(:)
    FTS_INT :: rbuffsize !< size of rbuff(:)
    FTS_INT :: start 
    FTS_INT :: end
 
    ! Arrays

    Integer   , Pointer :: mpirreq(:)               !< MPI requests receivers
    Integer   , Pointer :: mpisreq(:)               !< MPI requests senders
    XFTS_FLOAT, Pointer :: rbuff (:)             !< MPI receive Buffer
    Integer, Pointer    :: index_neig          (:)  !< field of boundary
    Integer, pointer :: ptr_direct_neig_edge (:)        !> direct pointer on index_neig_edge
    Integer, pointer :: ptr_vertperneig (:)
    Integer, pointer :: index_vertperneig (:)
    Integer, pointer :: ptr_neig_edge      (:) 
    Integer, pointer :: index_neig_edge    (:)

    


    !- End of header -----------------------------------------------------------

    !---------------------------------------------------------------------------
    ! [1] Init
    !---------------------------------------------------------------------------
    
    ! status
    info = 0
    timing = MPI_WTime()
    rank   = ftsl%ikeep(IKEEP_MPIRANK)
    ! mpi
    tagg             =  44
    comm             =  ftsl%comm


    
    ! interface description
    nb_neig             =  ftsl%lc_domain%nb_neig
    nloc                =  ftsl%lc_domain%size_domain
    index_neig          => ftsl%lc_domain%index_neig
    ptr_direct_neig_edge => ftsl%lc_domain%ptr_direct_neig_edge
    ptr_vertperneig     => ftsl%lc_domain%ptr_vertperneig
    index_vertperneig   => ftsl%lc_domain%index_vertperneig
    
    ! values


    !---------------------------------------------------------------------------
    ! [1] Allocate the buffers
    !---------------------------------------------------------------------------
    
    Nullify( rbuff )
    Nullify( mpirreq, mpisreq )

    Allocate( mpirreq(nb_neig), mpisreq(nb_neig), STAT= info )
    CHCKASSRT( info == 0, info )
    If (info < 0) Goto 9999

    rbuffsize = ptr_vertperneig(nb_neig+1)- ptr_vertperneig(1)
    Allocate(rbuff(rbuffsize), STAT= info )
    CHCKASSRT( info == 0, info )
    If (info < 0) Goto 9999

    Do i=1,rbuffsize 
       rbuff(i) = XFTS_FLOATZERO
    End Do


    !---------------------------------------------------------------------------
    ! [2] Send the computed value on the interface
    !---------------------------------------------------------------------------

    Do i =1,nb_neig

       start = nloc+ ptr_direct_neig_edge(i)
       end   = nloc+ ptr_direct_neig_edge(i+1) - 1

       ! Send
       Call MPI_Isend (z%v(start),(end-start+1),&
            XFTS_FLOATMPI,index_neig(i),tagg,comm,mpisreq(i),iinfo)
       ASSRT( iinfo  == MPI_SUCCESS )
    End Do



    !---------------------------------------------------------------------------
    ! [3] Start receiving 
    !---------------------------------------------------------------------------

    Do i=1,nb_neig
       
        start = ptr_vertperneig(i)
        end   = ptr_vertperneig(i+1) - 1

       Call MPI_Irecv( rbuff(start),(end-start+1), &
            XFTS_FLOATMPI,index_neig(i), tagg,comm,mpirreq(i),iinfo)
       ASSRT( iinfo  == MPI_SUCCESS )

    End Do


    !---------------------------------------------------------------------------
    ! [3] Unpack values from neigbhors and add contribution
    !---------------------------------------------------------------------------

    Do i=1,nb_neig

       Call MPI_WaitAny(nb_neig, mpirreq, neigh, MPI_STATUS_IGNORE, iinfo )
       ASSRT( iinfo  == MPI_SUCCESS )
       
       start = ptr_vertperneig(i)
       end   = ptr_vertperneig(i+1) - 1
       
       Do  k = start, end
          z%v( index_vertperneig(k) ) =0.5* z%v( index_vertperneig(k) ) +0.5* rbuff(k)  
       End Do
       
    End Do
    Call MPI_WaitAll(nb_neig, mpisreq, MPI_STATUSES_IGNORE, iinfo )
    ASSRT( iinfo  == MPI_SUCCESS )

    ! --------------------------------------------------------------------------
    ! [3] Exit routine
    ! --------------------------------------------------------------------------

9999 Continue
    If (Associated(rbuff)) Deallocate(rbuff)
    If (Associated(mpirreq)) Deallocate(mpirreq)
    If (Associated(mpisreq)) Deallocate(mpisreq)

    timing = MPI_WTime() - timing


  End Subroutine XFTS_PcdUpdate






End Module XFTS_solve_mod
