! Warning: XFTS_GENFILE_COMMENT
#include "fts_defs_f.h"
! [+] module : XFTS_sls_mod ----------------------------------------------------
!
!> Module for sparse linear system 
!!
!!
Module XFTS_sls_mod


  Use XFTS_sparse_matrix_type
  Use XFTS_sparse_matrix_mod
  Use XFTS_dense_matrix_type
  Use XFTS_dense_matrix_mod
  Use XFTS_mumps_mod
  
  Implicit None

  !> We want to solve the sparse linear system 
  !! \f[  A . sol = rhs \f]
  !! where A is a sparse matrix
  !!
  !! @par [the sparse matrix]
  !! - sm_A   : the sparse matrix A 
  !!
  !! @par [the right-hand-side]
  !! The right hand side (rhs) can be given in dense or in coordinate format.
  !! - rhs_is_sparse : controls the choosen format of the rhs
  !!          - FTS_TRUE  : the rhs is sparse
  !!          - FTS_FALSE : the rhs is dense
  !!          - Default  = FTS_FALSE
  !!          .
  !! - dm_rhs : pointer to the right-hand-side given in dense format
  !!          It is modified during the solving phase.
  !! - sm_rhs : pointer to the right-hand-side given in coordinate format
  !! . 
  !! @par [the solution]
  !! - dm_sol : a pointer to the solution.
  !!          - with a dense rhs , it points to the dm_rhs 
  !!          - with a sparse rhs, the user initialize/allocate dm_sol
  !!          - before the calling the solving phase
  !!          .
  !! @par [the solver]
  !! - sds    : the sparse direct solver used to solve the system
  !! .
  Type XFTS_sls_t; sequence

     Type(XFTS_sparse_matrix_t)      :: sm_A
     Type(XFTS_ARITHmumps_struc)     :: sds
     Type(XFTS_dense_matrix_t)       :: dm_rhs
     Type(XFTS_dense_matrix_t)       :: dm_sol

  End Type XFTS_sls_t

  ! list of routines
  Public :: XFTS_SLS_Nullify
  Public :: XFTS_SLS_Free

Contains
  ! [+] routine : XFTS_sls_nullify  --------------------------------------------
  !
  !> nullify a sparse linear system instance
  !! 
  !! @param[in,out ] sls    the sparse linear system instance to nullify
  !!
  !! @author Yohan Lee-tin-yien
  !!
  Subroutine XFTS_SLS_Nullify( sls, info )

    Implicit None

    !* Subroutine arguments *!
    Type(XFTS_sls_t), Intent(inout) :: sls
    Integer, Intent(  out) :: info

    !- End of header -----------------------------------------------------------


    Call XFTS_sm_Nullify(sls%sm_A,info)
    Call XFTS_dm_Nullify(sls%dm_rhs,info)
    CALL XFTS_dm_Nullify(sls%dm_sol,info)

    info = 0

  End Subroutine XFTS_SLS_Nullify

  ! [+] routine : XFTS_SLS_Free  -----------------------------------------------
  !
  !> Free a sparse linear system instance
  !! 
  !! Liberate the memory used by a sparse linear system instance.
  !!
  !!-----
  !!
  !! @param[in,out ] sls    the sparse linear system instance to nullify
  !! @param[   out ] info   the routine status
  !!
  !!-----
  !!
  !! @author Yohan Lee-tin-yien
  !!
  !! @verbatim
  !!  - Date     : Version : Comments
  !!  - 10/01/11 : 0.1a    : Create routine
  !!  - 14/02/11 : 0.1b    : update error/log handling.
  !! @endverbatim
  !!
  Subroutine XFTS_SLS_Free( sls, info )

    !* Modules *!
      Use XFTS_sparse_matrix_mod, Only : &
      XFTS_sm_free              ! routine
      Use XFTS_dense_matrix_mod, Only : &
      XFTS_dm_free, &           ! routine
      XFTS_dm_Nullify
      Use XFTS_mumps_mod, Only : &
      XFTS_mumps_Finalize       ! routine
      Use fts_log_mod
    Implicit None

    !* Arguments *!
    Type(XFTS_sls_t), Intent(inout) :: sls
    Integer , Intent(  out) :: info

    !* Local variables *!

    ! scalars
    Integer :: iinfo
    Integer :: istep
    Integer :: msgclass

    ! strings
    Character(len=FTSOLVER_STRL), Parameter :: rname = "SLS_Free"

    !- End of header -----------------------------------------------------------

    !---------------------------------------------------------------------------
    ! [1] Init
    !---------------------------------------------------------------------------
    iinfo = 0

    !---------------------------------------------------------------------------
    ! [2] Free the matrix
    !---------------------------------------------------------------------------
    
       Call XFTS_SM_Free(sls%sm_A  , iinfo )


    !---------------------------------------------------------------------------
    ! [3] Free the solver
    !---------------------------------------------------------------------------

       Call XFTS_mumps_Finalize (sls%sds  , iinfo )


    !---------------------------------------------------------------------------
    ! [4] Free the second member 
    !---------------------------------------------------------------------------


          Call XFTS_DM_Free (sls%dm_sol, iinfo )

          Call XFTS_DM_Free    (sls%dm_rhs, iinfo )

          Call XFTS_DM_Nullify (sls%dm_sol, iinfo )



    !---------------------------------------------------------------------------
    ! [5] Free the components
    !---------------------------------------------------------------------------



    !---------------------------------------------------------------------------
    ! [6] Finish
    !---------------------------------------------------------------------------

9999 Continue
  ! Print error/warning messages
  If ( iinfo /=  0 ) Then
     
     If ( iinfo > 0) msgclass = MSG_WARNING
     If ( iinfo < 0) msgclass = MSG_ERROR
     
     Select Case(istep) 
     Case Default
        Call fts_logWithInfo (msgclass,istep, Trim(rname)//&
             " at internal step =")
     End Select
     
  End If

  ! Set return code
  If ( iinfo == 0 ) info =  0
  If ( iinfo <  0 ) info = -istep
  If ( iinfo >  0 ) info = +istep


End Subroutine XFTS_SLS_FREE




End Module XFTS_sls_mod
