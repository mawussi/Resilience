! Warning: XFTS_GENFILE_COMMENT
!> FTSOLVER sparse matrix derived type 
!!
!!----
!! @author Yohan Lee-tin-yien
!! @version v 0.2a 
!!
!! @Date History
!!   Date    : version : Description 
!! - 21/06/11: v0.2b   : Extract FTS_sparse_matrix_enum (Y.Lee-tin-yien)
!! - 13/01/11: v0.2a   : Append comments (Yohan Lee-tin-yien)
!!                       - specifies symmetric cases
!     !                       - Append other comments
!     ! - 2010    : v0.2a   : first implementation (Yohan Lee-tin-yien)
#include "fts_defs_f.h"
      
      Module XFTS_sparse_matrix_type
      
!     Dependencies
      Use FTS_sparse_matrix_enum
      Implicit None
      
!     
      Type XFTS_sparse_matrix_t; Sequence
      
!     > storage format  (see FMT_*      )
      Integer                             :: fmt      
!     > symmetry        (see MATRIX_IS_*)
      Integer                             :: sym      
      
!     > number of rows
      FTS_INT                          :: m        
!     > number of columns
      FTS_INT                          :: n        
!     > number of non zero
      FTS_INT                          :: nnz      
      
     !     > rows               [nnz]   
      FTS_INT, pointer                 :: i   (:)  
!     > columns            [nnz]   
      FTS_INT, pointer                 :: j   (:)  
!     > values             [nnz]  
      XFTS_FLOAT, pointer               :: v   (:)  
      
     !     > compressed vector  [(ni|nj)+1]
      !integer, pointer                  :: cs  (:)  
      !> row    compressed vector
      FTS_INT  , dimension(:), pointer :: csr  
      !> column compressed vector
      FTS_INT  , dimension(:), pointer :: csc  
      
      End Type XFTS_sparse_matrix_t
      
      
      End Module XFTS_sparse_matrix_type
      

