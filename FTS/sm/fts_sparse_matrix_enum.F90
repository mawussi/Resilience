!> Enumerations for the sparse matrix module
!!
!!----
!! @author Yohan Lee-tin-yien
!! @version v 0.2a 
!!
!! @Date History
!!   Date      : version : Description 
!! - 2011/06/21: v0.2a   : Import from the "type" (Yohan Lee-tin-yien)
!!

Module FTS_sparse_matrix_enum
  Implicit None

  !-----------------------------------------------------------------------------
  ! [1.1] available format storages
  !-----------------------------------------------------------------------------
  ! Default format is the coordinate format.
  ! As third party libraries, or specific routines, requires compressed format
  ! The compressed vector can be appended to the classic coordinate format.
  
  !> i,j,value (coordinate) format (Default)
  Integer, Parameter :: SM_FMT_IJV = 0   
  !> ijv and compressed row    format
  Integer, Parameter :: SM_FMT_CSR = 1   
  !> ijv and compressed column format
  Integer, Parameter :: SM_FMT_CSC = 2   
  !> unset format
  Integer, Parameter :: SM_FMT_UNS = -1
  !> specific temporary format
  Integer, Parameter :: SM_FMT_SPE = 4

  !-----------------------------------------------------------------------------
  ! [1.2] available symmetries 
  !-----------------------------------------------------------------------------

  ! On symmetric and SPD matrices only the entries beneath and at the diagonal
  ! are stored.
  Integer, Parameter :: SM_SYM_IsGeneral   = 0 !< general (Default)
  Integer, Parameter :: SM_SYM_IsSPD       = 1 !< symetric positive definite
  Integer, Parameter :: SM_SYM_IsSymmetric = 2 !< symmetric

End Module FTS_sparse_matrix_enum
