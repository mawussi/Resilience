! Warning: XFTS_GENFILE_COMMENT
!> ftsolver dense matrix representation
!!
!! It is a simple encapsulation of blas matrices.
!!

Module XFTS_dense_matrix_type
  Implicit None

  !-----------------------------------------------------------------------------
  ! [1.1] Derived types 
  !-----------------------------------------------------------------------------

  !> structure representing a dense matrix
  Type XFTS_dense_matrix_t; sequence

     !> The number of rows.
     !!
     !! Specifies the number of rows of the matrix.
     !! "m" must be at least zero.
     Integer               :: m

     !> The number of colums.
     !!
     !! Specifies the number of columns of the matrix.
     !! "n" must be at least zero.
     Integer               :: n

     !> The leading dimension.
     !!
     !! The leading dimension of the array "val".  
     !! ld >= max(1,m).
     Integer               :: ld
     
     !> The matrix coefficients.
     !!
     !! array of DIMENSION ( ld, n ).
     !! the leading "m" by "n" part of the array "val" must
     !! contain the matrix of coefficients.
     XFTS_FLOAT, Pointer :: v (:)

     !> The arithmetic of the matrix. 
     !!
     !! It must be one referenced by the enumerations
     !! DM_ARITH_is*
     Integer               :: arith

     !> The symmetry of the matrix. (Optional Attribute)
     !!
     !! It must be one referenced by the enumerations
     !! DM_SYM_is*
     Integer               :: sym

     !> Which entries do you have. (Optional attribute)
     !!
     !! Specifies what entries you have in the matrix
     !! for example when the matrix is symmetric, you may only store
     !! the lower or the upper part.
     !! It must be one referenced by the enumerations
     !! DM_STORED_*
     Integer               :: stored

  End Type XFTS_dense_matrix_t

  !-----------------------------------------------------------------------------
  ! [1.2] Enumerations
  !-----------------------------------------------------------------------------

  !> @warning
  !!
  !! Only DM_TYPE_isGeneral is supported for the moment.
  !! 
  !! @todo support SPD
  !!
  !! @note DM_SYM_IsSPD means (widespread convention) :
  !!   - in real   , symmetric positive definite
  !!   - in complex, hermitian positive definite
  !!
  Integer, Parameter :: DM_SYM_Unset                    = -1
  Integer, Parameter :: DM_SYM_MaximumIndex             = 2

  Integer, Parameter :: DM_SYM_IsGeneral                = 0 ! DEFAULT
  Integer, Parameter :: DM_SYM_IsSPD                    = 1
  Integer, Parameter :: DM_SYM_IsSymmetric              = 2

  Integer, Parameter :: DM_ARITH_IsD                     = 0 
  Integer, Parameter :: DM_ARITH_IsS                     = 1 
  Integer, Parameter :: DM_ARITH_IsC                     = 2 
  Integer, Parameter :: DM_ARITH_IsZ                     = 3 

  Integer, Parameter :: DM_STORED_Unset = -1
  Integer, Parameter :: DM_STORED_ALL   = 0 ! Default
  Integer, Parameter :: DM_STORED_UPPER = 1
  Integer, Parameter :: DM_STORED_LOWER = 2

End Module XFTS_dense_matrix_type
