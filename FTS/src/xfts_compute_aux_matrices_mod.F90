#include "fts_defs_f.h"
#include "fts_macros_f.h"



Module XFTS_compute_aux_matrices_mod


      Implicit None 
      Public :: XFTS_compute_aux_matrices

      
      Character(len=FTSOLVER_STRL), Private, Parameter :: &
           FLNAME= "compute_aux_matrices_mod.90"
      contains


    ! [+] routine : XFTS_compute_aux_matrices---------------------------------------------
    !
    !>  
    !!
    !! Compute block row, precond and lsi matrices
    !! 
    !!
    !! @param[in,out ] ftsl     the ftsolver instance 
    !!
    !!
    Subroutine XFTS_compute_aux_matrices(ftsl)


      !* Module *!
      Use XFTS_ftsolver_type
      Use FTS_ftsolver_enum
      Use XFTS_dense_matrix_mod
      Use XFTS_sparse_matrix_mod
      Implicit None
      Include 'mpif.h'
            

      !* Arguments *!
      Type(XFTS_ftsolver_t), Intent(inout) :: ftsl

      !* Local variables *!
      
      ! Scalars
      Integer      :: rank
      Integer      :: info, i, k
      Integer      :: n, m, nnz
      Integer      :: iter_row, iter_col
      Integer      :: row, col
      !pointer
      Integer,pointer      :: csr(:),  bloc_row_csr(:), tmp_csr(:), csc(:)
      Integer,pointer      :: j(:), bloc_row_j(:), ai(:)
      XFTS_FLOAT, pointer  :: v(:), bloc_row_v(:), v_lsi(:)

      ! Derived types
      Type(XFTS_sparse_matrix_t)    :: sm_A
      ! End of header ----------------------------------------------------------

      !------------------------------------------------------------------
      ! Compute the local additive Schwarz preconditioner to factorize !
      !------------------------------------------------------------------
      rank   = ftsl%ikeep(IKEEP_MPIRANK)
      m =      ftsl%lc_domain%size_domain_delta1
      nnz =    ftsl%sm_blockrow_delta1%csr(m+1) -1


      Nullify(csr, tmp_csr, csc)
      Nullify(j, ai)
      Nullify(v, v_lsi)
      Allocate(csr(nnz))
      Allocate(tmp_csr(nnz))
      Allocate(j(nnz))
      Allocate(v(nnz))


      Do i=1,m+1
         csr(i) =ftsl%sm_blockrow_delta1%csr(i)
      End Do

      Do i=1,nnz
         j(i)   =ftsl%sm_blockrow_delta1%j(i)
         v(i)   =ftsl%sm_blockrow_delta1%v(i)
      End Do

      n = maxval(j)

      !delete extra columns
      tmp_csr(1) = 1
      iter_col = 0
      Do iter_row=1,m
         Do K=csr(iter_row),csr(iter_row+1)-1
           If( j(k) .LE. m) Then
              iter_col = iter_col+1
              j(iter_col) = j(k)
              v(iter_col) = v(k)
           End IF
        End Do
        tmp_csr(iter_row+1) = iter_col+1
      End Do

      nnz = iter_col
      Call XFTS_sm_CreateFromData(ftsl%sm_precond%sm_A,&
           SM_FMT_CSR,&
           ftsl%sm_blockrow_delta1%sym,&
           m,&
           m,&
           nnz,&
           tmp_csr,&
           j,&
           v,&
           info )
      Call XFTS_sm_convert(ftsl%sm_precond%sm_A,&
           SM_FMT_IJV, info )
      CHCKASSRT( info >= 0, info )

      !-------------------------------------
      !Create the local block row
      !------------------------------------

      m =      ftsl%lc_domain%size_domain
      nnz =    ftsl%sm_blockrow_delta1%csr(m+1) -1
      
      Nullify(bloc_row_csr)
      Nullify(bloc_row_j)
      Nullify(bloc_row_v)
      Allocate(bloc_row_csr(nnz))
      Allocate(bloc_row_j(nnz))
      Allocate(bloc_row_v(nnz))
      
      Do i=1,m+1
         bloc_row_csr(i) =ftsl%sm_blockrow_delta1%csr(i)
      End Do
      
      Do i=1,nnz
         bloc_row_j(i)   =ftsl%sm_blockrow_delta1%j(i)
         bloc_row_v(i)   =ftsl%sm_blockrow_delta1%v(i)
      End Do
      n = maxval( bloc_row_j)
      
      Call XFTS_sm_CreateFromData(&
           ftsl%sm_blockrow,&
           SM_FMT_CSR,&
           ftsl%sm_blockrow_delta1%sym,&
           m,&
           n,&
           nnz,&
            bloc_row_csr,&
            bloc_row_j,&
            bloc_row_v,&
           info )
      Call XFTS_sm_convert(ftsl%sm_blockrow, SM_FMT_CSC, info )
      
      !---------------------------------------------
      ! delete column out of local domain          !
      ! to get squarre matrix used for the LI      !
      !--------------------------------------------!

      tmp_csr(1) = 1
      iter_col = 0
      Do iter_row=1,m
         Do K= bloc_row_csr(iter_row), bloc_row_csr(iter_row+1)-1
           If( j(k) .LE. m) Then
              iter_col = iter_col+1
               bloc_row_j(iter_col) =  bloc_row_j(k)
               bloc_row_v(iter_col) =  bloc_row_v(k)
           End IF
        End Do
         tmp_csr(iter_row+1) = iter_col+1
      End Do

      nnz = iter_col

      Call XFTS_sm_CreateFromData(ftsl%sm_li,&
           SM_FMT_CSR,&
           ftsl%sm_blockrow_delta1%sym,&
           m,&
           m,&
           nnz,&
           tmp_csr,&
           bloc_row_j,&
           bloc_row_v,&
           info )


      Call XFTS_sm_convert(ftsl%sm_li,&
           SM_FMT_IJV, info ) 
      CHCKASSRT( info >= 0, info )

      !--------------------------------------------------
      ! Convert ftsl%sm_blockrow_delta1 into CSC format !
      ! extrat the rectangular matrix for LSI           !
      !--------------------------------------------------


      Call XFTS_sm_convert(ftsl%sm_blockrow_delta1, SM_FMT_CSC, info )
      CHCKASSRT( info >= 0, info )
      
      n =      ftsl%lc_domain%size_domain
      nnz =    ftsl%sm_blockrow_delta1%csc(n+1) -1

      Allocate(csc(nnz))
      Allocate(ai(nnz))
      Allocate(v_lsi(nnz))
      
      Do i=1,n+1
         csc(i) =ftsl%sm_blockrow_delta1%csc(i)
      End Do
      
      Do i=1,nnz
         ai(i)  =ftsl%sm_blockrow_delta1%i(i)
         v_lsi(i)   =ftsl%sm_blockrow_delta1%v(i)
      End Do
      m = maxval(ai)

      Call XFTS_sm_CreateFromData(&
           ftsl%sm_lsi,&
           SM_FMT_CSC,&
           ftsl%sm_blockrow_delta1%sym,&
           m,&
           n,&
           nnz,&
           ai,&
           csc,&
           v_lsi,&
           info )
      Call XFTS_sm_convert(ftsl%sm_lsi, SM_FMT_IJV, info )


      IF (Associated(csr)) Deallocate(csr)
      IF (Associated(bloc_row_csr)) Deallocate(bloc_row_csr)
      IF (Associated(tmp_csr)) Deallocate(tmp_csr)
      IF (Associated(j)) Deallocate(j)
      IF (Associated(bloc_row_j)) Deallocate(bloc_row_j)
      IF (Associated(v)) Deallocate(v)
      IF (Associated(bloc_row_v)) Deallocate(bloc_row_v)
      IF (Associated(csc)) Deallocate(csc)
      IF (Associated(ai)) Deallocate(ai)
      IF (Associated(v_lsi)) Deallocate(v_lsi)      
    End Subroutine XFTS_compute_aux_matrices


End Module XFTS_compute_aux_matrices_mod
