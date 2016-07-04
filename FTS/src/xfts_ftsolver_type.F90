! Warning: XFTS_GENFILE_COMMENT
#include "fts_defs_f.h"

!> Module defining the derived type "XFTS_ftsolver_t"
  Module XFTS_ftsolver_type
    
    !* Modules *!
    Use XFTS_sparse_matrix_type
    Use XFTS_sls_mod
    Use FTS_part_type, Only : &
           local_domain_t
    Use FTS_mem_mod, Only : &
         FTS_mem_t
         

   Implicit none
   Include "fts_env_type_f.inc"
      


    Type XFTS_ftsolver_t; Sequence
     ! ===============
     ! USER PARAMETERS
     ! ===============
     !   
     !> MPI communicator
     !!   
     !! The MPI Communicator used by ftsolver.\n
     !! It must be set with a valid MPI Communicator before any call
     !! to ftsolver.
     !!    
     Integer :: comm 

     ! Input matrix (in coordinate format)
     ! ---------------------------------------
     
     !> matrix symmetry
     !! - 0 = unsymmetric
     !! - 1 = symmetric
     !! - 2 = SPD
     !! .
     Integer                :: sym 
     !> order of the matrix
     Integer                :: n
     !> number of entries in the sparse matrix.
     !! If the matrix is symmetric or SPD,
     !! user only gives half of the entries.
     Integer                :: nnz
     Integer      , pointer :: rows (:)
     Integer      , pointer :: cols (:)
     XFTS_FLOAT , pointer :: values    (:)

     !> ask to write the input matrix to a file
     Character(len = FTSOLVER_STRL) :: write_matrix

     ! Right-hand-side (in dense format ordered by columns)
     ! ----------------------------------------------------
     Integer  :: nrhs  
     XFTS_FLOAT , pointer :: rhs (:) 

     ! Solution (in dense format ordered by columns)
     ! ----------------------------------------------------
     XFTS_FLOAT , pointer :: sol (:) 

     ! Controls
     ! --------
     Integer      :: job 
     Integer(kind=4) :: icntl  ( FTSOLVER_ICNTL_SIZE  ) 
     Real   (kind=8) :: rcntl  ( FTSOLVER_RCNTL_SIZE  ) 
     ! Statistics
     ! -----------
     ! version
     Character(len=FTSOLVER_STRL) :: version

     ! on this process (MPI)

     Integer(kind=4) :: iinfo ( FTSOLVER_IINFO_SIZE ) 
     Real   (kind=8) :: rinfo ( FTSOLVER_RINFO_SIZE ) 
     
     ! on all process (MPI)

     Integer(kind=4) :: iinfomin ( FTSOLVER_IINFO_SIZE ) 
     Integer(kind=4) :: iinfomax ( FTSOLVER_IINFO_SIZE ) 
     Real   (kind=8) :: iinfoavg ( FTSOLVER_IINFO_SIZE ) 
     Real   (kind=8) :: iinfosig ( FTSOLVER_IINFO_SIZE ) 

     Real   (kind=8) :: rinfomin ( FTSOLVER_RINFO_SIZE ) 
     Real   (kind=8) :: rinfomax ( FTSOLVER_RINFO_SIZE ) 
     Real   (kind=8) :: rinfoavg ( FTSOLVER_RINFO_SIZE )
     Real   (kind=8) :: rinfosig ( FTSOLVER_RINFO_SIZE ) 

     Integer(kind=4) :: iinfog   ( FTSOLVER_IINFOG_SIZE ) 
     Real   (kind=8) :: rinfog   ( FTSOLVER_RINFOG_SIZE ) 

     ! =====================
     ! Internal working data
     ! =====================
     
     ! internal controls
     ! -----------------
     Integer(kind=4) :: ikeep ( FTSOLVER_IKEEP_SIZE )    
     Real   (kind=8) :: rkeep ( FTSOLVER_RKEEP_SIZE )    

     ! memory
     Type(FTS_mem_t) :: mem

     ! environement
     Type(fts_env_t) :: env


     ! -------------------------------------------------------
     ! Local domain description 
     type(local_domain_t)     :: lc_domain 


     ! Local right hand side
     XFTS_FLOAT, pointer :: lc_rhs (:)

     ! bloc row delta1 matrix
     Type(XFTS_sparse_matrix_t) :: sm_blockrow_delta1

      ! local block row 
      Type(XFTS_sparse_matrix_t) :: sm_blockrow

     ! saved scalings on the local matrix
     XFTS_FLOAT, pointer :: row_scaling (:)  
     XFTS_FLOAT, pointer :: col_scaling (:)  
      


     ! Additive Schwartz  preconditioner
      Type(XFTS_sls_t)  :: sm_precond
     
     ! Rectangular matrix to solve LSI     
      Type(XFTS_sparse_matrix_t)  :: sm_lsi

      ! squarre matrix to solve LI     
      Type(XFTS_sparse_matrix_t)  :: sm_li
      
  End Type XFTS_ftsolver_t


End Module XFTS_ftsolver_type
