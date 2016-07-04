! Warning: XFTS_GENFILE_COMMENT
#include "fts_defs_f.h"

!> set of utilitary routines to generate data
Module XFTS_toolkit_gen_mod

  !* No Implicit Typing *!

  Implicit None
  
  !* Private constants *!
  
  Integer, Private, Parameter :: stderr = 6

  !* Access specifiers *!

  Public :: XFTS_gen_sol_ones
  Public :: XFTS_gen_sol_prand
  Public :: XFTS_gen_rhs

  !* Implementations *!

Contains


  ! [+] routine : XFTS_gen_sol_ones ----------------------------------------------
  !
  !> gen a SOL with solution of size n, with equals to vector ones
  !!
  Subroutine XFTS_gen_sol_ones ( n, sol )

    Implicit None

    ! Arguments
    FTS_INT           , Intent(in ) :: n
    XFTS_FLOAT, Pointer, Intent(out) :: sol (:)

    ! Local variables
    FTS_INT           :: i

    ! End of header ------------------------------------------------------------
    
    Allocate (sol( n ))

    Do i=1,n
      sol(i) = XFTS_FLOATONE
    End Do

  End Subroutine XFTS_gen_sol_ones


  ! [+] routine : XFTS_gen_sol_prand ----------------------------------------------
  !
  !> gen a pseudo random solution of size n
  !!
  Subroutine XFTS_gen_sol_prand ( n , sol )

    Implicit None

    ! Arguments
    FTS_INT           , Intent(in ) :: n
    XFTS_FLOAT, Pointer, Intent(out) :: sol (:)

    ! Interfaces to lapack routines
    XFTS_FLOAT, External :: XFTS_ARITHlarnd

    ! Local variables
    FTS_INT           :: i
    Integer              :: idist
    Integer              :: iseed(4)

    ! End of header ------------------------------------------------------------

    idist = 1
    iseed = (/0,0,0,1/)
    Allocate (sol( n ))

    Do i=1,n
       sol(i) = XFTS_ARITHlarnd( IDIST, ISEED )
    End Do

  End Subroutine XFTS_gen_sol_prand


  ! [+] routine : XFTS_gen_rhs ---------------------------------------------------
  !
  !> gen a second member RHS from the matrix SMATRIX and the solution SOL
  !!
  Subroutine XFTS_gen_rhs(smatrix, sol, rhs )
    Use XFTS_sparse_matrix_mod
    Use XFTS_dense_matrix_mod
    Implicit None

    ! Arguments
    Type(XFTS_sparse_matrix_t), Intent(in ) :: smatrix
    XFTS_FLOAT, Pointer, Intent(in ) :: sol (:)
    XFTS_FLOAT, Pointer, Intent(out) :: rhs (:)

    ! Local variables
    Integer              :: iinfo
    FTS_INT           :: n
    FTS_INT           :: i
    Type(XFTS_dense_matrix_t) :: dm_sol
    Type(XFTS_dense_matrix_t) :: dm_rhs

    ! End of header ------------------------------------------------------------
    
    n = smatrix%n
    ! init
    Allocate (rhs( n ))

    Call XFTS_dm_Create( dm_sol, n, 1, n, iinfo)
    Call XFTS_dm_Create( dm_rhs, n, 1, n, iinfo)
    
    Do i=1,n
       dm_sol%v(i) = sol(i)
    End Do
    ! perform matrix vector product
    Call XFTS_sm_VectorProduct(smatrix,0,0, dm_sol, dm_rhs,iinfo)
    Do i=1,n
       rhs(i) = dm_rhs%v(i)
    End Do

    ! Finish
    Call XFTS_dm_Free(dm_sol,iinfo)
    Call XFTS_dm_Free(dm_rhs,iinfo)


  End Subroutine XFTS_gen_rhs


End Module XFTS_toolkit_gen_mod

