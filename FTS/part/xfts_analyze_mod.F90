
#include"fts_defs_f.h"
#include "fts_macros_f.h"

      Module XFTS_analyze_mod
      Use FTS_get_neighbor_mod
      Use FTS_part_type
      Use XFTS_ftsolver_mod 
      Use FTS_ftsolver_enum
      Use fts_error_mod
      Implicit None
      Include "scotchf.h"




      Type scotch_graph_t
      sequence
     ! scalars
      Integer(kind=SCOTCH_NUMKIND) :: baseval
      Integer(kind=SCOTCH_NUMKIND) :: vertnbr
      Integer(kind=SCOTCH_NUMKIND) :: edgenbr
     
     ! array
      Integer(kind=SCOTCH_NUMKIND), Pointer :: verttab(:)
      Integer(kind=SCOTCH_NUMKIND), Pointer :: vendtab(:)
      Integer(kind=SCOTCH_NUMKIND), Pointer :: velotab(:)
      Integer(kind=SCOTCH_NUMKIND), Pointer :: edgetab(:)
     
     ! opaque scotch structure 
     Real(kind=8) :: data(SCOTCH_GRAPHDIM)
  End type scotch_graph_t


  Character(len=FTSOLVER_STRL), Private, Parameter :: &
       FLNAME= "XMPH_ARITHfts_analyze_mod.90"

      Public    :: XFTS_analyze
      Public    :: XFTS_distsystem_compute
      Public    :: XFTS_CreateGraph
      Public    :: FTS_scotch_partgraph
      Private   :: FTS_sm_graph2scotch_graph
      Private   :: FTS_init_describe_local
      Private   :: FTS_Free_describe_local 
      Private   :: FTS_print_describe_local
      Private   :: FTS_get_local_index
      contains 


  !routine : FTS_init_describe_local  ----------------------
  !
  !> nullify pointers in describe_local type
  !! 
  !! @param[in,out ] describe_lc   
  !!


subroutine FTS_init_describe_local(describe_lc, info)

      
      Implicit None
      !subroutine argument
      Type(describe_local_t), Intent(inout) :: describe_lc
      Integer              , Intent(out)    :: info

      integer       , parameter :: DEFAULT_INT   = -9999
      
    ! Scalars
      describe_lc%size_domain         = DEFAULT_INT
      describe_lc%size_vertperneig    = DEFAULT_INT
      describe_lc%size_neig_edge      = DEFAULT_INT
      describe_lc%size_sec_neig_edge  = DEFAULT_INT
      
      describe_lc%nb_neig             = DEFAULT_INT
      describe_lc%nb_sec_neig         = DEFAULT_INT

! Arrays
      Nullify( describe_lc%index_domain        )
      Nullify( describe_lc%index_neig          )
      Nullify( describe_lc%ptr_direct_neig_edge       )
      Nullify( describe_lc%ptr_vertperneig         )
      Nullify( describe_lc%index_vertperneig       )
      Nullify( describe_lc%index_neig_edge         )
      Nullify( describe_lc%index_sec_neig          )
      Nullify( describe_lc%ptr_sec_neig_edge       )
      Nullify( describe_lc%index_sec_neig_edge     )
      info = 0
End subroutine  FTS_init_describe_local



  !routine : FTS_Free_describe_local  ----------------------
  !
  !> Free pointers in describe_local type
  !! 
  !! @param[in,out ] describe_lc   
  !! @param[out ]    info   


      Subroutine FTS_Free_describe_local(describe_lc, info)


      Implicit None      
      !subroutine argument
      Type(describe_local_t), Intent(inout) :: describe_lc
      Integer              , Intent(out)    :: info

      

      If (Associated(describe_lc%index_domain))&
      Deallocate( describe_lc%index_domain)
      
      If (Associated(describe_lc%index_neig))&
      Deallocate( describe_lc%index_neig)
      
      If (Associated(describe_lc%ptr_direct_neig_edge))&
      Deallocate( describe_lc%ptr_direct_neig_edge)
      
      If (Associated(describe_lc%ptr_vertperneig))&
      Deallocate( describe_lc%ptr_vertperneig)
      
      If (Associated(describe_lc%index_vertperneig))&
      Deallocate( describe_lc%index_vertperneig)
      
      If (Associated(describe_lc%index_neig_edge))&
      Deallocate( describe_lc%index_neig_edge)
      
      If (Associated(describe_lc%index_sec_neig))&
      Deallocate( describe_lc%index_sec_neig )
      
      If (Associated(describe_lc%ptr_sec_neig_edge))&
      Deallocate( describe_lc%ptr_sec_neig_edge)

      If (Associated(describe_lc%index_sec_neig_edge))&
      Deallocate( describe_lc%index_sec_neig_edge)
      
      info =0
      
      End Subroutine FTS_Free_describe_local

    ! [+] routine : XFTS_Analyze ---------------------------------------------
    !> perform the analysis phase
    !!
    !!  distribute the global linear system into local linear systems (one
    !!  per MPI process).
    !!
    !!  @param[in,out] ftsl the FTsolver instance
    !!
    !!

      Subroutine XFTS_analyze(ftsl)

      !* Modules & co. *!
      Use XFTS_ftsolver_type
      Use XFTS_sparse_matrix_mod

      Implicit None
      Include 'mpif.h'

      !* Arguments *!
      Type(XFTS_ftsolver_t), Intent(inout) :: ftsl

      
      Integer :: iinfo, ignore
      Integer    :: rank, master
      ! Derived types
      Type(XFTS_sparse_matrix_t) :: gb_A


      ![--] Get MPI information
      rank   = ftsl%ikeep(IKEEP_MPIRANK)
      master = ftsl%ikeep(IKEEP_HOSTRANK)


      !gb_A will contain the sm in ijval format
      Call XFTS_ftsolver_pretreatInputMatrix(ftsl,master,rank,gb_A, iinfo)
      FTS_ONFAILURE_GOTO9999(iinfo)


      Call XFTS_distsystem_compute(ftsl,master,gb_A,ftsl%comm,iinfo)
      FTS_ONFAILURE_GOTO9999(iinfo)

      !   CHCKASSRT( .False.,ignore)
      Call XFTS_sm_free(gb_A, iinfo)
      FTS_ONFAILURE_GOTO9999(iinfo)


      !-------------------------------------------------------------------------
      ! [-] Exit
      !-------------------------------------------------------------------------

9999  Continue

      ! [--] set return code
      ftsl%iinfo( IINFO_STATUS ) = iinfo 

      End   Subroutine XFTS_analyze

    ! [+] routine : XFTS_distsystem_compute ---------------------------------------------
    !> convert the sparse matrix into a scotch graph format
    !! partition the graph and compute neighbors 
    !!  distribute the global linear system into local linear systems (one
    !!  per MPI process).
    !!

    !!  @param[in,out] ftss the FTsolver instance
    !!  @param[in,out] gb_A the global matrix t
    !!
    !!


 Subroutine XFTS_distsystem_compute( ftsl,master,sm_global_A, comm,iinfo)

     
      Use XFTS_ftsolver_type
      Use XFTS_sparse_matrix_mod

      Implicit None
      Include 'mpif.h'

      !* Subroutine arguments *!
      Type(XFTS_ftsolver_t)       , Intent(inout) :: ftsl
      Integer                   , Intent(in)    :: master
      Type(XFTS_sparse_matrix_t), Intent(inout) :: sm_global_A
      Integer                   , Intent(inout) :: comm
      Integer                   , Intent(  out) :: iinfo


      !* Local Variables *!

      ! Scalars

      Integer :: i, k, m, n
      Integer :: nnz_a, size_aj, col, ptr_row
      Integer :: current_row, current_col 
      Integer(8) :: max_size
      Integer :: info, rank, maxrank, tagg, domain
      Integer :: nbdom, partstrat,  current_part
      Integer :: vertnbr, edgenbr, current_index 
      Integer ::  iter_index, iter_begin, iter_end
      Integer :: lc_index, gl_index, current_begin
      Integer :: size_local2global, size_buffer_neig
      Integer :: size_buffer_send, size_buffer_recv
      Integer :: size_buffer_delta2, size_delta1
      Integer :: size_buffer_sec_neig
      Integer :: sec_neig, iter_sec_neig, size_sec_neig_edge
      Integer :: nb_secin_neig, secin_neig, max_secin
      Integer :: status(MPI_STATUS_SIZE)
      Integer :: mpistatus(MPI_STATUS_SIZE) 

      ! pointers 
      integer,  pointer :: parttab  (:) , verttab(:), edgetab(:)
      Integer, pointer  :: index_local2global (:) , index_buffer_neig(:)
      Integer, pointer  :: index_buffer_delta2(:)
      Integer, pointer  :: index_buffer_send(:)
      Integer, pointer  :: index_buffer_recv(:)
      Integer, pointer  :: tmp_lc2gl(:),  mpiReq(:)
      Integer, pointer  :: csr(:), aj(:)
      Integer, pointer  :: local_sum(:)
      Integer, pointer  :: global_sum(:)
      Integer, pointer  :: index_buffer_sec_neig(:)

      XFTS_FLOAT, pointer  :: val(:) , index_buffer_rhs(:)
      
      
      ! Derived types
      Type(matrix_graph_t) :: graph_global_A
      Type(describe_local_t)  :: describe_lc
      
      Call MPI_Comm_rank(ftsl%comm,rank,info)
      ASSRT( info == MPI_SUCCESS )
      
      Call MPI_Comm_size(ftsl%comm,nbdom,info)
      ASSRT( info == MPI_SUCCESS )

      !Initialize describe type
      Call  FTS_init_describe_local(describe_lc, info)
      Nullify(index_local2global)
      Nullify(index_buffer_neig)
      Nullify(index_buffer_sec_neig)
      Nullify(index_buffer_delta2)
      Nullify(index_buffer_recv)
      Nullify(index_buffer_send)
      Nullify(tmp_lc2gl)
      Nullify(mpiReq)
      Nullify(parttab)
      Nullify(verttab)
      Nullify(edgetab)
      Nullify(local_sum)
      Nullify(global_sum)
      Nullify(csr)
      Nullify(aj)
      Nullify(val)
      Nullify(index_buffer_rhs)

      If (master == rank) Then

         ! transform the sparse matrix into a adjacency graph
         Call XFTS_CreateGraph(sm_global_A, graph_global_A , iinfo)
         If (iinfo < 0) Goto 9999      
           
          ! patition the graph in nbdom with scotch
         Call  FTS_scotch_partgraph (nbdom, graph_global_A, parttab, info )        
         
         ! initialize variables to call get_distance 
         vertnbr = graph_global_A%ndof
         edgenbr = graph_global_A%xadj(graph_global_A%ndof+1)-1

         Nullify(verttab)
         Nullify(edgetab) 
         Allocate(verttab(vertnbr+1))
         Allocate(edgetab(edgenbr))

         Do i=1, vertnbr+1
            verttab(i) = INT(graph_global_A%xadj(i),SCOTCH_NUMKIND)
         End Do

         Do i=1, edgenbr
            edgetab(i) = INT(graph_global_A%adjncy(i),SCOTCH_NUMKIND)
         End Do
         
         maxrank = nbdom -1

         ! master extracts it neighbor information

         CALL FTS_distance(vertnbr, edgenbr, verttab, edgetab,&
              parttab,  0, describe_lc,  info) 


          ftsl%lc_domain%size_domain = describe_lc%size_domain
 
          ftsl%lc_domain%size_domain_delta1 = &
          describe_lc%size_domain + describe_lc%size_neig_edge
          
          ftsl%lc_domain%size_domain_delta2 =&
          ftsl%lc_domain%size_domain_delta1+describe_lc%size_sec_neig_edge
         
         !Allocate memory for index_local2global
          Allocate(ftsl%lc_domain%index_local2global(&
          ftsl%lc_domain%size_domain_delta2))
          
          !Allocate memory local rhs
          Allocate(ftsl%lc_rhs(ftsl%lc_domain%size_domain))
          
          !fill index_local2global with interior indices
          ! and save local rhs
          iter_begin =1
          iter_end   =  describe_lc%size_domain 
          
          
          Do iter_index=iter_begin, iter_end
             current_index = describe_lc%index_domain(iter_index)
             !save index
             ftsl%lc_domain%index_local2global(iter_index) = current_index
             ! save local rhs
             ftsl%lc_rhs(iter_index) = ftsl%rhs(current_index)
          End Do
          
          !fill index_local2global with distance one indices
          iter_begin = iter_end +1
          iter_end  = iter_end + describe_lc%size_neig_edge
 
          i= 1
          Do iter_index= iter_begin, iter_end
             ftsl%lc_domain%index_local2global(iter_index) =&
             describe_lc%index_neig_edge(i)
             i = i+1
          End Do
             
          !fill index_local2global with distance two indices
          iter_begin = iter_end +1
          iter_end   = iter_end + describe_lc%size_sec_neig_edge 
          
          i=1
          Do iter_index= iter_begin, iter_end
             ftsl%lc_domain%index_local2global(iter_index) =&
             describe_lc%index_sec_neig_edge(i)
             i = i+1
          End Do

                      
          ! number of neighbors
          ftsl%lc_domain%nb_neig  = describe_lc%nb_neig
         
          Allocate(ftsl%lc_domain%index_neig(ftsl%lc_domain%nb_neig))
             ! neighbors indices
             iter_begin = 1
             iter_end   = describe_lc%nb_neig
             
             i =1
             Do iter_index=iter_begin, iter_end
                ftsl%lc_domain%index_neig(iter_index) = describe_lc%index_neig(i)
                i = i+1
             End Do
             
             ! direct pointer on neighbors edge
             Allocate(ftsl%lc_domain%ptr_direct_neig_edge(ftsl%lc_domain%nb_neig+1))
 
             iter_begin = 1
             iter_end =  describe_lc%nb_neig+1
             
             Do iter_index=iter_begin,iter_end
                ftsl%lc_domain%ptr_direct_neig_edge(iter_index) =&
                describe_lc%ptr_direct_neig_edge(iter_index)
             End Do
 
             ! pointer vert connected to neighbors
             Allocate(ftsl%lc_domain%ptr_vertperneig(ftsl%lc_domain%nb_neig+1))
             iter_begin = 1
             iter_end =  describe_lc%nb_neig+1
             
             i =1
             Do iter_index =iter_begin, iter_end
                ftsl%lc_domain%ptr_vertperneig(iter_index) =&
                describe_lc%ptr_vertperneig(i)
                i= i+1
             End Do
 
             ! size of index_vertperneig
              ftsl%lc_domain%size_vertperneig = describe_lc%size_vertperneig
 
 
             ! vertices associated with neighbors
             ! need to be converted to local index
             Allocate(ftsl%lc_domain%index_vertperneig(describe_lc%size_vertperneig))
             iter_begin = 1
             iter_end =  describe_lc%size_vertperneig
             
             i=1
             Do iter_index=iter_begin, iter_end
                gl_index = describe_lc%index_vertperneig(i)
                lc_index = FTS_get_local_index(&
                ftsl%lc_domain%index_local2global, gl_index)
                ftsl%lc_domain%index_vertperneig(iter_index) = lc_index
                i = i+1
             End Do
             
 
             ! size index_neig_edge
              ftsl%lc_domain%size_neig_edge = describe_lc%size_neig_edge
             
             ! vertices on neighbors associated with 
             ! connection edges
             ! need to be converted to local index
              Allocate(ftsl%lc_domain%index_neig_edge(describe_lc%size_neig_edge))
              iter_begin = 1
              iter_end =  describe_lc%size_neig_edge
              
              i=1
              Do iter_index=iter_begin, iter_end
                 gl_index = describe_lc%index_neig_edge(i)
                 lc_index = FTS_get_local_index(&
                 ftsl%lc_domain%index_local2global, gl_index)
                 ftsl%lc_domain%index_neig_edge(iter_index) = lc_index
                 i = i+1
              End Do
              
 
             !-----------------------------------------------------
             ! extract  distance two neighbors description !
             !-----------------------------------------------------
             
             
             ! number of neighbors
              ftsl%lc_domain%nb_sec_neig = describe_lc%nb_sec_neig
             
             ! neighbors indices
              Allocate(ftsl%lc_domain%index_sec_neig(describe_lc%nb_sec_neig))
              iter_begin = 1
              iter_end   = describe_lc%nb_sec_neig 
             
             i =1
             Do iter_index=iter_begin, iter_end
               ftsl%lc_domain%index_sec_neig(iter_index) = describe_lc%index_sec_neig(i)
                i = i+1
             End Do
             
 
             !  pointer on distance two neighbors edges
              Allocate(ftsl%lc_domain%ptr_sec_neig_edge(describe_lc%nb_sec_neig+1))
             iter_begin = 1
             iter_end   = describe_lc%nb_sec_neig+1
 
             i = 1
             Do iter_index=iter_begin,iter_end
                ftsl%lc_domain%ptr_sec_neig_edge(iter_index) = describe_lc%ptr_sec_neig_edge(i)
                i = i+1
             End Do
 
             
             !size_sec_neig_edge
             ftsl%lc_domain%size_sec_neig_edge = describe_lc%size_sec_neig_edge
             
             !index_sec_neig_edge
             Allocate(ftsl%lc_domain%index_sec_neig_edge(describe_lc%size_sec_neig_edge))
             iter_begin = 1
             iter_end =  describe_lc%size_sec_neig_edge
             
             i = 1
             Do iter_index=iter_begin,iter_end
                gl_index = describe_lc%index_sec_neig_edge(i)
                lc_index = FTS_get_local_index(&
                ftsl%lc_domain%index_local2global, gl_index)
                ftsl%lc_domain%index_sec_neig_edge(iter_index) = lc_index
                i = i+1
             End Do


             ! Allocate memory for master and other procs
             !TODO : optimize allocation with mpi_max
             size_local2global = ftsl%lc_domain%size_domain_delta2
             max_size = size_local2global*1000

             !Allocate max memory for ai, aj, val
             Allocate(csr(max_size))
             Allocate(aj(max_size))
             Allocate(val(max_size))
             Allocate(index_buffer_rhs(max_size))

             !Allocate memory for index_local2global
             Allocate(index_local2global(max_size))
             
             !Allocate memory for index_buffer_neig
             Allocate(index_buffer_neig(max_size))

             Allocate(index_buffer_sec_neig(max_size))

             !Allocation memory index_buffer_delta2
             Allocate(index_buffer_delta2(max_size))


             !--------------------------------
             ! Compute matrix block row delta1 entries !
             !--------------------------------

             ! [-] Convert the matrix into CSR format

             Call XFTS_sm_convert( sm_global_A, SM_FMT_CSR, info )
             CHCKASSRT( info >= 0, info )
             nnz_a   = 0
             size_aj = 0
             ptr_row = 1
             csr(ptr_row) = 1
             
             Do iter_index=1, ftsl%lc_domain%size_domain_delta1
                current_row = ftsl%lc_domain%index_local2global(iter_index)
                
                Do i = 1, size_local2global
                   current_col = ftsl%lc_domain%index_local2global(i)
                   
                   Do k=sm_global_A%csr(current_row),sm_global_A%csr(current_row+1)-1
                      col=sm_global_A%j(k)
                      If ( col .EQ. current_col) Then
                         
                         lc_index = FTS_get_local_index(&
                         ftsl%lc_domain%index_local2global, col)
                         nnz_a = nnz_a +1
                         aj(nnz_a) = lc_index
                         val(nnz_a) = sm_global_A%v(k)
                         Exit
                      End If
                   End Do
                End Do
                ptr_row = ptr_row+1
                csr(ptr_row) = nnz_a+1
             End Do

            !--------------------------------------
            ! create  the block row delta1 matrix !
            !--------------------------------------
             m = ftsl%lc_domain%size_domain_delta1
             n = size_local2global
             Call XFTS_sm_CreateFromData(&
             ftsl%sm_blockrow_delta1, SM_FMT_CSR, SM_SYM_IsGeneral,&
             m, n, nnz_a, &
             csr, aj, val,&
             info )
             FTS_ONFAILURE_RETURN(info)
             !------------------------------------!
             !send neighbor information to other  ! 
             !------------------------------------!

             Do domain=1,maxrank

                CALL FTS_distance(vertnbr, edgenbr, verttab, edgetab,&
                parttab,  domain, describe_lc,  info)
                !compute size of matrix delata2
                size_local2global = 4 + describe_lc%size_domain+&
                     describe_lc%size_neig_edge + describe_lc%size_sec_neig_edge

                size_delta1 = describe_lc%size_domain+&
                     describe_lc%size_neig_edge
                
                            
                ! fill index_local2global 
                index_local2global(1) = describe_lc%size_domain
                index_local2global(2) = describe_lc%size_neig_edge
                index_local2global(3) = describe_lc%size_sec_neig_edge
                index_local2global(4) = ftsl%n

                !fill index_local2global with interior indices
                iter_begin =5
                iter_end = iter_begin + describe_lc%size_domain -1
            
                i=1
                Do iter_index=iter_begin, iter_end
                   current_index  = describe_lc%index_domain(i)
                   !save index local2global
                   index_local2global(iter_index) = current_index
                   !save rhs
                   index_buffer_rhs(i) = ftsl%rhs(current_index)
                   i = i+1
                End Do

              !fill index_local2global with distance one indices
                iter_begin = iter_end +1
                iter_end  = iter_end + describe_lc%size_neig_edge
                
                i= 1
                Do iter_index= iter_begin, iter_end
                   index_local2global(iter_index) = describe_lc%index_neig_edge(i)
                   i = i+1
                End Do
            
               !fill index_local2global with distance two indices
                iter_begin = iter_end +1
                iter_end   = iter_end + describe_lc%size_sec_neig_edge 
                
                i=1
                Do iter_index= iter_begin, iter_end
                   index_local2global(iter_index) = describe_lc%index_sec_neig_edge(i)
                   i = i+1
                End Do
               ! send index_local2global
                tagg = domain+1
                CALL MPI_send(index_local2global(1), size_local2global,&
                MPI_INTEGER, domain, tagg, ftsl%comm, info )
                ASSRT( info == MPI_SUCCESS )

               !-----------------------------------------------------
               ! extract and send distance one neighbor description !
               !-----------------------------------------------------

                !compute size of the buffer that contains neighbor information
                size_buffer_neig = 5 + 3*describe_lc%nb_neig &
                + describe_lc%size_vertperneig + describe_lc%size_neig_edge
            
                ! number of neighbors
                index_buffer_neig(1) = describe_lc%nb_neig
                
                ! neighbors indices
                iter_begin = 2
                iter_end   = describe_lc%nb_neig + 1
            
                i =1
                Do iter_index=iter_begin, iter_end
                   index_buffer_neig(iter_index) = describe_lc%index_neig(i)
                   i = i+1
                End Do
            
                ! direct pointer on neighbors edge
                iter_begin = iter_end+1
                iter_end = iter_end + describe_lc%nb_neig+1
                
                i = 1
                Do iter_index=iter_begin,iter_end
                   index_buffer_neig(iter_index) = describe_lc%ptr_direct_neig_edge(i)
                   i = i+1
                End Do

                ! pointer vert connected to neighbors
                iter_begin = iter_end+1
                iter_end = iter_end+  describe_lc%nb_neig+1
                
                i =1
                Do iter_index =iter_begin, iter_end
                   index_buffer_neig(iter_index) = describe_lc%ptr_vertperneig(i)
                   i= i+1
                End Do
                
                ! size of index_vertperneig
                iter_end = iter_end+1
                index_buffer_neig(iter_end) = describe_lc%size_vertperneig
                
                ! vertices associated with neighbors
                ! need to be converted to local index
                iter_begin = iter_end +1
                iter_end = iter_end + describe_lc%size_vertperneig
                
                tmp_lc2gl => index_local2global(5: size_local2global)
                i=1
                Do iter_index=iter_begin, iter_end
                   gl_index = describe_lc%index_vertperneig(i)
                   lc_index = FTS_get_local_index(tmp_lc2gl, gl_index)
                   index_buffer_neig(iter_index) = lc_index
                   i = i+1
                End Do

                ! size index_neig_edge
                iter_end = iter_end+1
                index_buffer_neig(iter_end) = describe_lc%size_neig_edge
                
                ! vertices on neighbors associated with 
                ! connection edges
                ! need to be converted to local index
                iter_begin = iter_end +1
                iter_end = iter_end + describe_lc%size_neig_edge
            
                i=1
                Do iter_index=iter_begin, iter_end
                   gl_index = describe_lc%index_neig_edge(i)
                   lc_index = FTS_get_local_index(tmp_lc2gl, gl_index)
                   index_buffer_neig(iter_index) = lc_index
                   i = i+1
                End Do

                ! send index_buffer_neig
                tagg = domain+2
                CALL MPI_send(index_buffer_neig(1) , size_buffer_neig,&
                     MPI_INTEGER, domain, tagg, ftsl%comm, info )
                ASSRT( info == MPI_SUCCESS )

                !-----------------------------------------------------
                ! extract and send distance two neighbors description !
                !-----------------------------------------------------
                
                !compute send buffer size
                size_buffer_sec_neig = 3+ 2*describe_lc%nb_sec_neig &
                     + describe_lc%size_sec_neig_edge
                
                ! number of neighbors
                index_buffer_sec_neig(1) = describe_lc%nb_sec_neig
            
                ! neighbors indices
                iter_begin = 2
                iter_end   = describe_lc%nb_sec_neig + 1
                
                i =1
                Do iter_index=iter_begin, iter_end
                   index_buffer_sec_neig(iter_index) = describe_lc%index_sec_neig(i)
                   i = i+1
                End Do
                
                !  pointer on distance two neighbors edges
                iter_begin = iter_end+1
                iter_end = iter_end + describe_lc%nb_sec_neig+1
                
                i = 1
                Do iter_index=iter_begin,iter_end
                   index_buffer_sec_neig(iter_index) = describe_lc%ptr_sec_neig_edge(i)
                   i = i+1
                End Do
                
                !size_sec_neig_edge
                iter_end = iter_end +1
                index_buffer_sec_neig(iter_end) = describe_lc%size_sec_neig_edge
                
                !index_sec_neig_edge
                iter_begin = iter_end+1
                iter_end = iter_end + describe_lc%size_sec_neig_edge
                
                i = 1
                Do iter_index=iter_begin,iter_end
                   gl_index = describe_lc%index_sec_neig_edge(i)
                   lc_index = FTS_get_local_index(tmp_lc2gl, gl_index)
                   index_buffer_sec_neig(iter_index) = lc_index
                   i = i+1
                End Do
                
                !send index_buffer_sec_neig

                tagg = domain+3
                CALL MPI_send(index_buffer_sec_neig(1) , size_buffer_sec_neig,&
                MPI_INTEGER, domain, tagg, ftsl%comm, info )
                ASSRT( info == MPI_SUCCESS )
                
                ! compute block row matrix delta1 for each neighbor 
                ! convert the global sparse matrix in csr
                ! format. Then extract delta2 submatrix
                ! ATTENTION : convert indices into local 
                
                !decrement size_local2global to map
                ! the size of block row matrix delta2
                size_local2global = size_local2global - 4
                nnz_a   = 0
                size_aj = 0
                ptr_row = 1
                csr(ptr_row) = 1
                
                Do iter_index=1, size_delta1
                   current_row = tmp_lc2gl(iter_index)
                   
                   Do i = 1, size_local2global
                      current_col = tmp_lc2gl(i)
                  
                      Do k=sm_global_A%csr(current_row),sm_global_A%csr(current_row+1)-1
                         col=sm_global_A%j(k)
                         If ( col .EQ. current_col) Then
                            
                            lc_index = FTS_get_local_index(tmp_lc2gl, col)
                            nnz_a = nnz_a +1
                            aj(nnz_a) = lc_index
                            val(nnz_a) = sm_global_A%v(k)
                            Exit
                     End If
                  End Do
               End Do
               ptr_row = ptr_row+1
               csr(ptr_row) = nnz_a+1
            End Do
            !--------------------------------------------------------------------
            ! Prepare buffer end send the block row delta1 matrix  and local rhs!
            !--------------------------------------------------------------------

            size_buffer_delta2 = 1 & !  nnz
            + ptr_row &         ! size csr
            + nnz_a             ! size indices col 
            
            
            index_buffer_delta2(1) = nnz_a
            
            !CSR
            Do iter_index=1, ptr_row
               index_buffer_delta2(1+iter_index) = csr(iter_index)
            End Do
            
            ! j and val
            iter_begin = 1+ ptr_row+1
            iter_end   = iter_begin + nnz_a -1


            i = 1
            Do iter_index= iter_begin, iter_end
               index_buffer_delta2(iter_index) = aj(i)
               i= i+1
            End Do

            !send index_buffer_delta2
            tagg = domain+5
            CALL MPI_send(index_buffer_delta2(1) , size_buffer_delta2,&
            MPI_INTEGER, domain, tagg, ftsl%comm, info )
            ASSRT( info == MPI_SUCCESS )

            ! send val
            tagg = domain+6
            CALL MPI_send(val(1) , nnz_a,&
            XFTS_FLOATMPI, domain, tagg, ftsl%comm, info )
            ASSRT( info == MPI_SUCCESS )
            tagg = domain+7
            CALL MPI_send(index_buffer_rhs(1) ,describe_lc%size_domain,&
            XFTS_FLOATMPI, domain, tagg, ftsl%comm, info )
            ASSRT( info == MPI_SUCCESS )
         End Do

      Else

         !recv index_local2global
         tagg = rank+1
         Call MPI_Probe(master,tagg,ftsl%comm,mpiStatus,info)
         Call MPI_Get_count(mpiStatus,MPI_INTEGER,size_local2global,info)
         ASSRT( info == MPI_SUCCESS )


            Allocate(index_local2global(size_local2global))
            
            Call MPI_recv(index_local2global(1),size_local2global ,&
            MPI_INTEGER, master, tagg, ftsl%comm, status,info)
            ASSRT( info == MPI_SUCCESS )
         
         
            !store received data 
            ftsl%lc_domain%size_domain = index_local2global(1)
            ftsl%lc_domain%size_domain_delta1 = index_local2global(1)+&
            index_local2global(2)
            ftsl%lc_domain%size_domain_delta2 = ftsl%lc_domain%size_domain_delta1 +&
            index_local2global(3)
            ftsl%n = index_local2global(4)
            
            iter_begin =5
            iter_end = size_local2global
         
            !Allocate memory for ftsl%lc_domain%index_local2global
            Allocate(ftsl%lc_domain%index_local2global(size_local2global-4))
            i=1
            Do iter_index=iter_begin, iter_end
               ftsl%lc_domain%index_local2global(i) = index_local2global(iter_index) 
               i = i+1
            End Do

            !-------------------------------------------
            ! Receive distance one neighbor description!
            !------------------------------------------

            tagg = rank+2
            Call MPI_Probe(master,tagg,ftsl%comm,mpiStatus,info)
            Call MPI_Get_count(mpiStatus,MPI_INTEGER,size_buffer_neig,info)
            ASSRT( info == MPI_SUCCESS )
            
            Allocate(index_buffer_neig(size_buffer_neig))
            
            Call MPI_recv(index_buffer_neig(1),size_buffer_neig ,&
            MPI_INTEGER, master, tagg, ftsl%comm, status,info)
            ASSRT( info == MPI_SUCCESS )
         

            ! store neighbor information
            ftsl%lc_domain%nb_neig = index_buffer_neig(1)

            ! neighbors indices
            iter_begin = 2
            iter_end   = ftsl%lc_domain%nb_neig + 1


            ! Allocate memory for index_neig
            Allocate(ftsl%lc_domain%index_neig(ftsl%lc_domain%nb_neig))
            i =1
            Do iter_index=iter_begin, iter_end
              ftsl%lc_domain%index_neig(i) = index_buffer_neig(iter_index) 
               i = i+1
            End Do


            ! direct pointer
            iter_begin = iter_end +1
            iter_end= iter_end + ftsl%lc_domain%nb_neig + 1
            
            !Allocate memory for ptr_direct_neig_edge
            Allocate(ftsl%lc_domain%ptr_direct_neig_edge(ftsl%lc_domain%nb_neig+1))
            i=1
            Do iter_index=iter_begin, iter_end
               ftsl%lc_domain%ptr_direct_neig_edge(i) = index_buffer_neig(iter_index)
               i = i+1
            End Do

            ! vertperneig pointer
            iter_begin = iter_end +1
            iter_end= iter_end + ftsl%lc_domain%nb_neig + 1
            
            !Allocate memory for ptr_vertperneig
            Allocate(ftsl%lc_domain%ptr_vertperneig(ftsl%lc_domain%nb_neig+1))
            i=1
            Do iter_index=iter_begin, iter_end
               ftsl%lc_domain%ptr_vertperneig(i) = index_buffer_neig(iter_index)
               i = i+1
            End Do

            ! size vertperneig
            iter_end = iter_end+1
            ftsl%lc_domain%size_vertperneig = index_buffer_neig(iter_end)
            
            ! index vertperneig
            iter_begin = iter_end+1
            iter_end = iter_end + ftsl%lc_domain%size_vertperneig

            ! allocate memory for ftsl%lc_domain%index_vertperneig
            Allocate(ftsl%lc_domain%index_vertperneig(ftsl%lc_domain%size_vertperneig))
            i=1
            Do iter_index=iter_begin, iter_end
               ftsl%lc_domain%index_vertperneig(i) = index_buffer_neig(iter_index)
               i=i+1
            End Do

            ! size of index_neig_edge
            iter_end = iter_end + 1
            ftsl%lc_domain%size_neig_edge = index_buffer_neig(iter_end)
            
            
            !vertices on neighbors associated with 
            ! connection edges
            iter_begin = iter_end +1
            iter_end = iter_end + ftsl%lc_domain%size_neig_edge 
            

            !Allocate memory for ftsl%lc_domain%index_neig_edge
            Allocate(ftsl%lc_domain%index_neig_edge(ftsl%lc_domain%size_neig_edge))
            i = 1
            Do iter_index= iter_begin, iter_end
               ftsl%lc_domain%index_neig_edge(i) = index_buffer_neig(iter_index)
               i = i+1
            End Do

            !-------------------------------------------
            ! Receive distance two neighbor description!
            !------------------------------------------
            tagg = rank+3
            Call MPI_Probe(master,tagg,ftsl%comm,mpiStatus,info)
            Call MPI_Get_count(mpiStatus,MPI_INTEGER,size_buffer_sec_neig,info)
            ASSRT( info == MPI_SUCCESS )
            
            Allocate(index_buffer_sec_neig(size_buffer_sec_neig))
            
            Call MPI_recv(index_buffer_sec_neig(1),size_buffer_sec_neig ,&
            MPI_INTEGER, master, tagg, ftsl%comm, status,info)
            ASSRT( info == MPI_SUCCESS )
         

            ! store neighbor information
            ftsl%lc_domain%nb_sec_neig = index_buffer_sec_neig(1)

            ! neighbors indices
            iter_begin = 2
            iter_end   = ftsl%lc_domain%nb_sec_neig + 1


            ! Allocate memory for index_sec_neig
            Allocate(ftsl%lc_domain%index_sec_neig(ftsl%lc_domain%nb_sec_neig))
            i =1
            Do iter_index=iter_begin, iter_end
              ftsl%lc_domain%index_sec_neig(i) = index_buffer_sec_neig(iter_index) 
               i = i+1
            End Do

            ! direct pointerg
            iter_begin = iter_end +1
            iter_end= iter_end + ftsl%lc_domain%nb_sec_neig + 1
            
            !Allocate memory for ptr_sec_neig_edge
            Allocate(ftsl%lc_domain%ptr_sec_neig_edge(ftsl%lc_domain%nb_sec_neig+1))
            i=1
            Do iter_index=iter_begin, iter_end
               ftsl%lc_domain%ptr_sec_neig_edge(i) = index_buffer_sec_neig(iter_index)
               i = i+1
            End Do

             ! size of index_sec_neig_edge
            iter_end = iter_end + 1
            ftsl%lc_domain%size_sec_neig_edge = index_buffer_sec_neig(iter_end)
            
            !vertices on neighbors associated with 
            ! contributing to distance two relationship

            iter_begin = iter_end +1
            iter_end = iter_end + ftsl%lc_domain%size_sec_neig_edge 
            

            !Allocate memory for ftsl%lc_domain%index_neig_edge
            Allocate(ftsl%lc_domain%index_sec_neig_edge(ftsl%lc_domain%size_sec_neig_edge))
            i = 1
            Do iter_index= iter_begin, iter_end
               ftsl%lc_domain%index_sec_neig_edge(i) = index_buffer_sec_neig(iter_index)
               i = i+1
            End Do

            !-------------------------------------------
            ! Receive block row delta1 matrix data!
            !------------------------------------------
            tagg = rank+5
            Call MPI_Probe(master,tagg,ftsl%comm,mpiStatus,info)
            Call MPI_Get_count(mpiStatus,MPI_INTEGER,size_buffer_delta2,info)
            ASSRT( info == MPI_SUCCESS )
                        
            Allocate(index_buffer_delta2(size_buffer_delta2))

            Call MPI_recv(index_buffer_delta2(1),size_buffer_delta2 ,&
            MPI_INTEGER, master, tagg, ftsl%comm, status,info)
            ASSRT( info == MPI_SUCCESS )
            
            nnz_a = index_buffer_delta2(1)
            m     =    ftsl%lc_domain%size_domain_delta1
            n     =    ftsl%lc_domain%size_domain_delta2
            iter_begin = 2
            iter_end = m+2
            
            Allocate(csr(m+1))
            
            i= 1
            Do iter_index = iter_begin, iter_end
               csr(i) = index_buffer_delta2(iter_index)
               i = i+1
            End Do

            iter_begin = iter_end+1
            iter_end   = iter_end + nnz_a 

            Allocate(aj(nnz_a))
            i = 1
            Do iter_index = iter_begin, iter_end
               aj(i)  =  index_buffer_delta2(iter_index)
               i = i+1
            End Do
            ! receive val
            Allocate(val(nnz_a))

            tagg = rank+6
            Call MPI_recv(val(1),nnz_a ,&
            XFTS_FLOATMPI, master, tagg, ftsl%comm, status,info)
            ASSRT( info == MPI_SUCCESS )
            
           
            Call XFTS_sm_CreateFromData(&
            ftsl%sm_blockrow_delta1 ,SM_FMT_CSR, SM_SYM_IsGeneral,&
            m, n, nnz_a, &
            csr, aj, val,&
            info )
            FTS_ONFAILURE_RETURN(info)


            ! receive local rhs
            Allocate(ftsl%lc_rhs(ftsl%lc_domain%size_domain))
            tagg = rank+7
            Call MPI_recv(ftsl%lc_rhs(1),ftsl%lc_domain%size_domain ,&
            XFTS_FLOATMPI, master, tagg, ftsl%comm, status,info)
            ASSRT( info == MPI_SUCCESS )

         End If

         ! Deallocate memory 
         Call FTS_Free_describe_local(describe_lc, info)
         Call XFTS_sm_free(sm_global_A, info)
         If (Associated (parttab)) Deallocate(parttab)
         If (Associated (verttab)) Deallocate(verttab)
         If (Associated (edgetab)) Deallocate(edgetab)
         If (Associated (index_local2global)) Deallocate(index_local2global)
         If (Associated (index_buffer_neig)) Deallocate(index_buffer_neig)
         If (Associated (index_buffer_delta2)) Deallocate(index_buffer_delta2) 
         If (Associated (index_buffer_sec_neig)) Deallocate(index_buffer_sec_neig)
         If (Associated (csr)) Deallocate(csr)
         If (Associated (aj)) Deallocate(aj)      
         If (Associated (val)) Deallocate(val)      
         If (Associated (index_buffer_rhs)) Deallocate(index_buffer_rhs)      



         !-------------------------------------------------!
         ! Each domain will receive information from       !
         ! domains that see it as distance two neighbors.  !
         ! The information mainly consists of indices      ! 
         ! quired for LSI rhs computation contribution     !
         !-------------------------------------------------!

         Allocate(mpiReq(nbdom))
         Allocate(local_sum(nbdom))
         Allocate(global_sum(nbdom))
         Allocate(index_buffer_send(&
         ftsl%lc_domain%size_sec_neig_edge))
         tagg = 777

         ! informe each proc of the number of procs
         ! that seen it as distance two neighbors
         ! This helps to know how many time a proc
         ! sould set a receive fonction 
         Do i=1, nbdom
            local_sum(i) = 0
            global_sum(i) = 0
         End Do

         Do iter_sec_neig=1,ftsl%lc_domain%nb_sec_neig
            sec_neig = ftsl%lc_domain%index_sec_neig(iter_sec_neig)
            local_sum(sec_neig +1) = 1
         End Do

         Call MPI_Allreduce(local_sum, global_sum, nbdom , MPI_INTEGER,&
         MPI_SUM, ftsl%comm, info)
         nb_secin_neig = global_sum(rank+1)

         ftsl%lc_domain%nb_secin_neig = nb_secin_neig


         Call MPI_Allreduce(ftsl%lc_domain%size_sec_neig_edge,&
         max_secin, 1 , MPI_INTEGER, MPI_MAX, ftsl%comm, info)

         Allocate(index_buffer_recv(max_secin*nb_secin_neig))
         Allocate(ftsl%lc_domain%index_secin_neig(nb_secin_neig))
         Allocate(ftsl%lc_domain%ptr_secin_neig_edge(nb_secin_neig+1))
         Allocate(ftsl%lc_domain%index_secin_neig_edge(max_secin*nb_secin_neig))
         
         
         ! Send information 
         current_begin = 1
         Do iter_sec_neig=1, ftsl%lc_domain%nb_sec_neig
            
            sec_neig = ftsl%lc_domain%index_sec_neig(iter_sec_neig)
            iter_begin =  ftsl%lc_domain%ptr_sec_neig_edge( iter_sec_neig)
            iter_end   =  ftsl%lc_domain%ptr_sec_neig_edge(iter_sec_neig+1) -1
            size_sec_neig_edge = iter_end - iter_begin +1
            
            !convert local indices to global
            i = 1
            Do iter_index=iter_begin, iter_end
               lc_index = ftsl%lc_domain%index_sec_neig_edge(iter_index)
               gl_index =  ftsl%lc_domain%index_local2global(lc_index)
               index_buffer_send(iter_index) = gl_index
               i = i+1
            End Do

            CALL MPI_Isend(index_buffer_send(current_begin), &
            size_sec_neig_edge, MPI_INTEGER, sec_neig,&
            tagg, ftsl%comm, mpiReq(iter_sec_neig), info)
            ASSRT( info == MPI_SUCCESS )
            current_begin = current_begin + size_sec_neig_edge
          End Do

         
         !--------------------!
         !       receive      !        
         !--------------------!
         current_begin = 0
         ftsl%lc_domain%ptr_secin_neig_edge(1) = 1
         Do secin_neig=1, nb_secin_neig
            
            Call MPI_Probe(MPI_ANY_SOURCE,tagg,ftsl%comm,mpiStatus,info)
            ASSRT( info == MPI_SUCCESS )
            
            Call MPI_Get_count(mpiStatus,MPI_INTEGER,size_buffer_recv,info)
            ASSRT( info == MPI_SUCCESS )
            
            Call MPI_Recv(index_buffer_recv(1), size_buffer_recv,&
            MPI_INTEGER, mpiStatus(MPI_SOURCE),&
            tagg,ftsl%comm, MPI_STATUS_IGNORE,info)
            ASSRT( info == MPI_SUCCESS )

            ftsl%lc_domain%index_secin_neig(secin_neig) = mpiStatus(MPI_SOURCE)
            ftsl%lc_domain%ptr_secin_neig_edge(secin_neig+1) = &
            ftsl%lc_domain%ptr_secin_neig_edge(secin_neig)+ size_buffer_recv
            
            ! set edge, need convertion global to local
            Do iter_index=1, size_buffer_recv
               gl_index = index_buffer_recv(iter_index)
               lc_index = FTS_get_local_index(&
               ftsl%lc_domain%index_local2global, gl_index)
               ftsl%lc_domain%index_secin_neig_edge(&
               current_begin+iter_index) = lc_index
            End Do
            current_begin = current_begin+ size_buffer_recv
         End Do 

         Call MPI_Waitall(ftsl%lc_domain%nb_sec_neig,&
         mpiReq, MPI_STATUSES_IGNORE,info)
         ASSRT( info == MPI_SUCCESS )
         ftsl%lc_domain%size_secin_neig_edge = current_begin

9999  Continue
  
         If (Associated (local_sum)) Deallocate(local_sum)      
         If (Associated (global_sum)) Deallocate(global_sum)
         If (Associated (mpiReq)) Deallocate(mpiReq)          
         If (Associated (index_buffer_send)) Deallocate(index_buffer_send)
         If (Associated (index_buffer_recv)) Deallocate(index_buffer_recv)
End Subroutine XFTS_distsystem_compute


 ! [+] routine : FTS_scotch_partgraph ---------------------------------
    !
    !> partition the matrix graph with scotch
    !!
    !!----
    !!
    !! @param nbdom       wanted number of domain
    !! @param parttab     the permutation
    !! @param info        status of the subroutine (0=SUCCESS)
    !!
    !!----
    !!
    !!
    !!
      subroutine FTS_scotch_partgraph   & ! intents
      (nbdom, graph_global_A,            & ! in
      parttab, info )                    ! out

      !* Modules *!

      Implicit None
      !* Arguments *!

      integer   , intent(in)                  :: nbdom
      type(matrix_graph_t), intent(in) :: graph_global_A

      integer, intent(inout), pointer :: parttab  (:)
      integer, intent(out)          :: info
      


      ! scalars
      Integer(kind=SCOTCH_SYSKIND) :: iinfo 
      Integer(kind=SCOTCH_NUMKIND) :: vertnbr
      Integer                      :: i
      ! arrays     
      Integer(kind=SCOTCH_NUMKIND), Pointer :: permtab (:)



      ! types
      Type(scotch_graph_t) :: scotch_graph
      Real(kind=8) :: scotch_strat(SCOTCH_STRATDIM)
     


      !-------------------------------------------------------------------------
      ! [1.0] Init
      !-------------------------------------------------------------------------
      
      ! scalars
      vertnbr = graph_global_A%ndof

      ! arrays
      Nullify(permtab)
      Allocate(permtab(vertnbr  ), stat=iinfo)
      Allocate(parttab(vertnbr  ), stat=iinfo)
      If (iinfo /= 0 ) Goto 9999

      
      !-------------------------------------------------------------------------
      ! [2.0] Init SCOTCH structures
      !-------------------------------------------------------------------------
      
      ! convert graph
      Call FTS_sm_graph2SCOTCH_Graph(graph_global_A, scotch_graph, iinfo)
      If (iinfo < 0 ) Goto 9999

      ! init the strategy
      Call SCOTCHFStratInit(scotch_strat,iinfo)
      If (iinfo /= 0 ) Goto 9999

      !-------------------------------------------------------------------------
      ! [3.0] Partition the SCOTCH GRAPH
      !-------------------------------------------------------------------------

      ! call scotch

      Call SCOTCHFGraphPart(&
           scotch_graph%data, nbdom, scotch_strat, &
           permtab, iinfo )

      Do i=1, vertnbr
         parttab(i) = INT(permtab(i), kind=SCOTCH_NUMKIND)
      End Do

      !-------------------------------------------------------------------------
      ! [5.0] End
      !-------------------------------------------------------------------------

9999  Continue
            ! Free Memory
      Call SCOTCHfGraphExit(scotch_graph%data, iinfo)
      Call SCOTCHFStratExit(scotch_strat,iinfo)

      If( Associated(permtab)) DeAllocate(permtab)

      If( Associated(scotch_graph%vendtab)) Nullify   (scotch_graph%vendtab)
      If( Associated(scotch_graph%verttab)) DeAllocate(scotch_graph%verttab)
      If( Associated(scotch_graph%edgetab)) DeAllocate(scotch_graph%edgetab)
      If( Associated(scotch_graph%velotab)) DeAllocate(scotch_graph%velotab)
      
      End Subroutine FTS_scotch_partgraph




    ! [+] routine : FTS_sm_graph2scotch_graph ----------------------------------------
    !
    !> Convert a sparse matrix graph into a scotch graph
    !!
    !!----
    !!
    !! @param [in ] sm_graph   the sparse matrix graph
    !! @param [out] scotch_graph   the scotch matrix graph
    !! @param [out] info           the status
    !!
    !!----
    !!
    !! @author Yohan Lee-tin-yien
    Subroutine FTS_sm_graph2scotch_graph(sm_graph, scotch_graph,info)

      !* Modules *!
      Implicit None

      !* Arguments *!
      Type(matrix_graph_t) , Intent(in)  :: sm_graph
      Type(scotch_graph_t) , Intent(out) :: scotch_graph
      Integer                     , Intent(out) :: info

      !* Local variables *!


      Integer(kind=SCOTCH_SYSKIND) :: iinfo
      Integer(kind=SCOTCH_NUMKIND) :: baseval, vertnbr, edgenbr
      Integer(kind=SCOTCH_NUMKIND) :: i

      !- End of header ---------------------------------------------------------

      !-------------------------------------------------------------------------
      ! [1.0] Init
      !-------------------------------------------------------------------------

      ! warning : 
      ! - we need to copy data for 32/64 bits integer conversions
      ! - vendtab is not defined (see SCOTCH documentation about optional arg.)
      ! - edlotab is not defined
      ! - vlbltab is not defined
      ! - velotab is only defined when the vertices are weighted
      

      ! scalars

      scotch_graph%baseval =  1
      scotch_graph%vertnbr =  sm_graph%ndof
      scotch_graph%edgenbr =  sm_graph%xadj(sm_graph%ndof+1)-1
      scotch_graph%data(:) =  0.d0 ! unsignificant init

      baseval = scotch_graph%baseval
      vertnbr = scotch_graph%vertnbr
      edgenbr = scotch_graph%edgenbr

      Nullify(scotch_graph%verttab)
      Nullify(scotch_graph%vendtab)
      Nullify(scotch_graph%edgetab)
      Nullify(scotch_graph%velotab)

      Allocate(scotch_graph%verttab(vertnbr+1), stat=iinfo)
      Allocate(scotch_graph%edgetab(edgenbr), stat=iinfo)


      Do i=1, vertnbr+1
         scotch_graph%verttab(i) = INT(sm_graph%xadj(i),SCOTCH_NUMKIND)
      End Do

      Do i=1, edgenbr
         scotch_graph%edgetab(i) = INT(sm_graph%adjncy(i),SCOTCH_NUMKIND)
      End Do

      scotch_graph%vendtab => scotch_graph%verttab(2:vertnbr+1)


      Call SCOTCHfGraphInit(scotch_graph%data, iinfo)
      If(iinfo /= 0) iinfo = -iinfo
      If(iinfo /= 0) Goto 9999
      
      !-------------------------------------------------------------------------
      ! [2.0] Build
      !-------------------------------------------------------------------------


         Call SCOTCHfGraphBuild(scotch_graph%data, &
              scotch_graph%baseval    , scotch_graph%vertnbr, &
              scotch_graph%verttab(1) , scotch_graph%verttab(2), &
              scotch_graph%verttab(1), scotch_graph%verttab(1) , &
              scotch_graph%edgenbr, scotch_graph%edgetab(1) , &
              scotch_graph%edgetab(1), iinfo)

      If (iinfo /= 0) iinfo = -iinfo
      If (iinfo /= 0) Goto 9999

      !-------------------------------------------------------------------------
      ! [3.0] Check
      !-------------------------------------------------------------------------

   
      Call SCOTCHfGraphCheck(scotch_graph%data, iinfo )
      If (iinfo /= 0) iinfo = -iinfo
      If (iinfo /= 0) Goto 9999

      !-------------------------------------------------------------------------
      ! [4.0] End
      !-------------------------------------------------------------------------

9999  Continue

      ! free memory on error

    End Subroutine FTS_sm_graph2scotch_graph


! [+] routine : CreateGraph -----------------------------------------------
!
!> Create the graph corresponding to the matrix rows and columns 
!!
!! @param sm_global_A    the sparse matrix of the global linear system
!! @param graph_global_A the associated graph
!! @param info           status of the Subroutine (0=SUCCESS)
!!
!! @todo update algo comments
Subroutine XFTS_CreateGraph        &
     (sm_global_A,   & ! in
     & graph_global_A, info)     ! out

  !* Module(s) & co. *!
  Use XFTS_sparse_matrix_type
  Use fts_log_mod
  Implicit none
         
  !* Arguments *!

  Type(XFTS_sparse_matrix_t), intent(in) :: sm_global_A
  Type(matrix_graph_t), intent(out) :: graph_global_A
  Integer                    , intent(out) :: info
  
  !* Local variables *!

  ! Constants
  Integer, Parameter :: UNSET= -1

  ! Scalars
  Integer               :: iinfo
  Integer               :: istep
  Integer               :: msg_class
  Integer               :: symtype
  FTS_INT :: i, j, k
  FTS_INT :: ind
  FTS_INT            :: nnza, ndof
  FTS_INT            :: maxadj
  Logical               :: TO_BE_SYMETRIZED

  ! Arrays
  FTS_INT, Pointer :: ia(:), ja(:)
  FTS_INT, Pointer :: xadj  (:)   ! xadjacency
  FTS_INT, Pointer :: adjncy(:)   ! adjacency 

  ! Strings
  Character(len=FTSOLVER_STRL) :: rname="PART_CreateGraph"

  !- End of header ------------------------------------------------------------

  !----------------------------------------------------------------------------
  ! [1] Initialize and associate local variables
  !----------------------------------------------------------------------------

  ! init graph
  graph_global_A%ndof   = UNSET
  graph_global_A%nnz    = UNSET
  graph_global_A%maxadj = UNSET
  graph_global_A%algo   = UNSET

  Nullify(graph_global_A%xadj)
  Nullify(graph_global_A%adjncy)
  Nullify(graph_global_A%vwgt)

  ! init local scalars
  
  iinfo = 0
  symtype = sm_global_A%sym
  nnza = sm_global_A%nnz
  ndof = sm_global_A%m
  ia   => sm_global_A%i
  ja   => sm_global_A%j
  

  Select Case(symtype)
  Case(SM_SYM_IsGeneral  ) 
     TO_BE_SYMETRIZED = .False. ! structure was already symmetrized
  Case(SM_SYM_IsSymmetric,SM_SYM_IsSPD) 
     TO_BE_SYMETRIZED = .True.  
  Case Default            
     TO_BE_SYMETRIZED = .True.  
  End Select

  !----------------------------------------------------------------------------
  ! [2] Compute the graph according to its symmetry
  !----------------------------------------------------------------------------
  If( TO_BE_SYMETRIZED )Then

     !-------------------------------------------------------------------------
     ! [2.1] Transform (ia,ja) into graph data structure
     ! [ - ] (for symmetric --> triangular entries)
     !------------------------------------------------------------------------- 

     !-------------------------------------------------------------------------
     ! [2.1.1] Allocate xadj,adjcncy, vwgt
     !------------------------------------------------------------------------- 


     Allocate(xadj(ndof+1),STAT=iinfo)
     IF( iinfo /= 0) iinfo = -iinfo
     If (iinfo /= 0) Goto 9999

     Allocate(adjncy(2*(nnza)),STAT=iinfo)
     IF( iinfo /= 0) iinfo = -iinfo
     If (iinfo /= 0) Goto 9999

     !-------------------------------------------------------------------------
     ! [2.1.2] Construct xadj 
     !------------------------------------------------------------------------- 

     ! init
     Do i=1,ndof+1
        xadj(i) = 0
     End Do

     ! count the number of edges associated to each row 
     xadj(1)=1
     do k=1,nnza
        i=ia(k)
        j=ja(k)
        if(i /= j)then
           xadj(i+1)= xadj(i+1) + 1
           xadj(j+1)= xadj(j+1) + 1
        endif
     enddo

     ! sum the counts to form xadj
     xadj(1)=1
     do i=2,ndof+1
        xadj(i)= xadj(i) + xadj(i-1)
     enddo

     !-------------------------------------------------------------------------
     ! [2.1.3] Form the adjncy
     !------------------------------------------------------------------------- 
     
     !> @warning - about xadj
     !! xadj is complete at this step.
     !! but momentarily modify it to save memory. 
     !! value is retrieve after shifting.

     ! construct adjcncy (& modify xadj)
     Do k=1,nnza
        i=ia(k)
        j=ja(k)
        If(i /= j)Then

           ind          = xadj(i)
           adjncy(ind)  = j
           xadj(i)      = xadj(i) + 1 

           ind          = xadj(j)
           adjncy(ind)  = i
           xadj(j)      = xadj(j) + 1 

        End If
     End Do

     ! retrieve xadj by shifting it by 1
     Do k=1,ndof
        xadj(ndof-k+2) = xadj(ndof-k+1)
     End Do
     xadj(1)=1


  Else ! NOT TO_BE_SYMMETRIZED

     !-------------------------------------------------------------------------
     ! [2.2] transform (ia,ja) into graph data structure
     ! [ - ] (for symmetric pattern matrices --> all entries)
     !-------------------------------------------------------------------------
     ! See also section [2.1] for full comments

     ! allocate xadj, adjcny and vwgt 

     Allocate(xadj(ndof+1),STAT=iinfo)
     IF( iinfo /= 0) iinfo = -iinfo
     If (iinfo /= 0) Goto 9999

     Allocate(adjncy(nnza),STAT=iinfo)
     IF( iinfo /= 0) iinfo = -iinfo
     If (iinfo /= 0) Goto 9999

     ! form xadj
     Do i=1,ndof+1
        xadj(i) = 0
     End Do

     xadj(1)=1
     do k=1,nnza
        i=ia(k)
        j=ja(k)
        if(i /= j)then
           xadj(i+1)= xadj(i+1) + 1
        endif
     enddo

     xadj(1)=1
     do i=2,ndof+1
        xadj(i)= xadj(i) + xadj(i-1)
     enddo

     ! form the adjncy
     do k=1,nnza
        i=ia(k)
        j=ja(k)
        if(i /= j)then
           ind          = xadj(i)
           adjncy(ind)  = j
           xadj(i)    = xadj(i) + 1 
        endif
     enddo

     do i=1,ndof
        xadj(ndof-i+2) = xadj(ndof-i+1)
     enddo
     xadj(1)=1

  end IF

  !-----------------------------------------------------------------------------
  ! [3] Compute statistics
  !-----------------------------------------------------------------------------

  maxadj = 0
  Do i=1,ndof
     k= xadj(i+1)-xadj(i) 
     maxadj= Max(maxadj,k)
  End Do
  
  !-----------------------------------------------------------------------------
  ! [4] Finish
  !-----------------------------------------------------------------------------

  ! Save computed data
  graph_global_A%ndof  = ndof
  graph_global_A%nnz   = nnza
  graph_global_A%maxadj   = maxadj

  graph_global_A%xadj    => xadj
  graph_global_A%adjncy  => adjncy


  ! Print error messages
9999 Continue
  If ( iinfo /=  0 ) Then

     If ( iinfo > 0) msg_class = MSG_WARNING
     If ( iinfo < 0) msg_class = MSG_ERROR     
     
  End If

  ! report success, failure
  If ( iinfo == 0 ) info =  0
  If ( iinfo <  0 ) info = - istep
  If ( iinfo >  0 ) info =   istep

      End Subroutine XFTS_CreateGraph


      subroutine FTS_print_describe_local(describe_lc, rank, info )

      Implicit None 

      ! subroutine arguments 
      Type(describe_local_t), Intent(in) :: describe_lc
      Integer, Intent(in)                :: rank
      Integer, Intent(inout)            :: info
      
      !local variables
      Integer  iter_neig, iter_sec_neig,i
      Integer nb_connected_vertices, nb_sec_connected_vertices
      Integer start_vertperneig, end_vertperneig
      Integer current_neig, current_sec_neig
      Integer start_edge
      Integer end_edge
      

      write(*,*) "Domain", rank
      write(*,*) "It has ", describe_lc%size_domain, "vertices: "
      write(*,*)"local vertices:",  describe_lc%index_domain
      write(*,*) "It has", describe_lc%nb_neig, "distance one neighbors"      
      write(*,*) "It has", describe_lc%nb_sec_neig, &
      "distance two neighbors"      
      write(*,*)
      write(*,*) "---------Distance one neighgors-----------"

      write(*,*)
      write(*,*) "---------Distance two neighgors-----------"
      
      Do iter_sec_neig=1,describe_lc%nb_sec_neig
         current_sec_neig = describe_lc%index_sec_neig(iter_sec_neig)  
         
         start_edge = describe_lc%ptr_sec_neig_edge(iter_sec_neig)
         end_edge = describe_lc%ptr_sec_neig_edge(iter_sec_neig+1)-1
         nb_sec_connected_vertices = end_edge - start_edge+1
         write(*,*) current_sec_neig, "contribute with",  nb_sec_connected_vertices, "vertices"
         write(*,*)  "distance two vertices:",&
         describe_lc%index_sec_neig_edge(start_edge:end_edge)
         write(*,*)
      End Do

      End subroutine FTS_print_describe_local


      Function FTS_get_local_index(this_array, gl_index)&
      result(lc_index)

      Implicit None 
      
      !Arguments
      Integer, pointer, Intent(in) :: this_array(:)
      Integer,          Intent(in) :: gl_index
      Integer                      :: lc_index


      !local variables 
      Integer  i, size_array
      
      lc_index = -1
      size_array = size(this_array)

      Do i=1,size_array
         If (this_array(i) .EQ. gl_index) Then
            lc_index = i
            exit 
         End If
      End Do
      
      End  Function FTS_get_local_index



      End Module XFTS_analyze_mod



      










