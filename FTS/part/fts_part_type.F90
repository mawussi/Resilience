#include "fts_defs_f.h"

!> derived types for FTSOLVER partitioner.
!!
Module FTS_part_type


  !* No Implicit typing *!
      Implicit None

  !lty> begin unsure comments

  ! [+] Type : ftsolver_rhs_partition_t ---------------------------------------
  !
  !> data used to split the global system into local system 
    Type ftsolver_rhs_partition_t
       sequence

       !* global data *! (only on master = 0)

       !> number of domains
       Integer*FTS_INTKIND :: nbdom 

       !> order of the global matrix 
       Integer*FTS_INTKIND :: gballndof          
       
       !> order of the schur complement
       !! @Note  it is a copy of ftsolver_domain_t%gballintrf.
       Integer*FTS_INTKIND :: gballintrf  

       !> ? vertex ? size
       INTEGER :: combivtxpcsz                 

       !> processus displacement for domain's interior 
       INTEGER, pointer :: procintdisp (:)      
       !> processus displacement for domain's interface
       INTEGER, pointer :: procintrfdisp  (:)   
       !> domain langrage  
       INTEGER, pointer :: domLg (:)            
       !> degree of freedom of each domain's interior 
       INTEGER, pointer :: domintdof (:)        

       !> degree of freedom of each domain's interface
       INTEGER, pointer :: domintrfdof (:)     
       
       !> permutation done on the global matrix
       INTEGER, pointer :: metperm (:)          

       !* Strategy *!
       INTEGER :: rhsway ! selected strategy (1..2) , forced to 2
       ! 1 rst strategy 
       INTEGER, pointer :: scatindices (:)      ! allow to scatter data ?
       ! 2 nd strategy   
       INTEGER, pointer :: scatlogicindices (:) ! allow to scatter data ?
       INTEGER, pointer :: domlogicintrfdof(:)  ! degree of freedom 

       !> global to local numbering in the interface (size = gballintrf)
       !! Note: (Yohan Lee-tin-yien 2011-10-03)
       !! gbtoloc is only useful with a global preconditioner.
       !! It is unused in the current version,
       !! and not exported in the read/write routine,
       !! but are left there for future use.
       Integer*FTS_INTKIND, pointer :: gbtoloc      (:)     

       
    end Type ftsolver_rhs_partition_t

    ! [+] Type : ftsolver_domain_t -----------------------------------------------
    !
    !> defines the data on each domain
    Type local_domain_t
       sequence
       !--------------------------------------------------------------------
       !> @par [ local matrix A_delta0, A_delta1, and A_delta2  description ]
       !--------------------------------------------------------------------

       !> number of rows (columns) in the local matrix A_delta0
       Integer           :: size_domain
       


       !> number of rows (columns) in the local matrix  A_delta1 
       !(only vertices at distance one)
       Integer           :: size_domain_delta1


       !> number of rows (columns) in the local matrix A_delta2
       !(only vertices at distance two)
       Integer           :: size_domain_delta2


       !> convertion local index to global 
       Integer, pointer  :: index_local2global (:)

       !------------------------------------------------
       !> @par [ distance one' neighbors  description  ]
       !------------------------------------------------
       
       !> number of distance one' neighbors 
       Integer          :: nb_neig

       !> index of distance one' neighbors
       Integer, pointer :: index_neig  (:)
       

       !> direct pointer on index_neig_edge
       Integer, pointer :: ptr_direct_neig_edge (:)

       !> pointer on  vertices of index_vertperneig (size: nb_neig +1)
       !> this data structure contains the same vertices stored 
       !> index_periph_vert. However, index_vertperneig is different
       !> because for each neighbor(domain) , it store successively the vertices 
       !> it is linked to. So if a vertex is linked to many neighbors
       !> it replicated 
       Integer, pointer :: ptr_vertperneig (:)

       !> size of  index_vertperneig
       Integer          :: size_vertperneig

       !> contains for each
       ! neighbor domain, the vertices of the periphery 
       ! of A_delta0 they are linked to. 
       Integer, pointer :: index_vertperneig (:)

       !> Pointer on edges at distance one
       !> useful if some vertices in index_vertperneig
       !> have many edges with a given domain 
       Integer, pointer :: ptr_neig_edge      (:) 
       
       !> size of  index_d1_edge
       Integer          :: size_neig_edge

       !> contains end vertices of edges at distance one
       Integer, pointer :: index_neig_edge    (:)

       !------------------------------------------------
       !> @par [ distance two' neighbors  description  ]
       !------------------------------------------------

       !> number of distance two' neighbors(domain)
       Integer          :: nb_sec_neig

       !> index of distance two' neighbors
       Integer, pointer :: index_sec_neig     (:)

       !> pointer on index_sec_neig_edge
       Integer, pointer :: ptr_sec_neig_edge       (:)

       !> size of index_sec_neig_edge
       Integer          :: size_sec_neig_edge 

       !> contains vertices at distance two
       Integer, pointer :: index_sec_neig_edge     (:)

       !------------------------------------------------
       !> @par  Description of neighbors that see 
       ! this domain as distanc two neighbor           
       !------------------------------------------------

       !> number incoming of distance two' neighbors(domain)
       Integer          :: nb_secin_neig

       !> index incoming of distance two' neighbors
       Integer, pointer :: index_secin_neig     (:)

       !> pointer on index_secin_neig_edge
       Integer, pointer :: ptr_secin_neig_edge       (:)

       !> size of index_secin_neig_edge
       Integer          :: size_secin_neig_edge 

       !> contains incoming vertices at distance two
       Integer, pointer :: index_secin_neig_edge     (:)

    end Type local_domain_t

    !< end unsure comments

    ! [+] Type : ftsolver_matrix_graph_t -----------------------------------------
    !
    !> the corresponding graph to a sparse square matrix
    !! matrix graph (optionally weighted) is represented by an adjacency list.
    !!
    !! @see metis manual section 5.1
    Type matrix_graph_t
       sequence
       !> order of the sparse matrix
       Integer*FTS_INTKIND          :: ndof       
       !> number of entries of sparse matrix (number of non-zeros)
       FTS_INT       :: nnz 

       !> see metis manual section 5.1 (size ndof+1)
       Integer*FTS_INTKIND, Pointer :: xadj   (:) 
       !> see metis manual section 5.1 (size nnz   )
       Integer*FTS_INTKIND, Pointer :: adjncy (:) 
       !> specifies the weight on vertices (size ndof)
       Integer*FTS_INTKIND, Pointer :: vwgt   (:) 

       !> algorithm used to split the graph
       Integer*FTS_INTKIND          :: algo       

       !* statistics *!

       !> maximum number of adjacent vertices in the graph
       Integer*FTS_INTKIND :: maxadj    
    end Type matrix_graph_t

    ! [+] Type : ftsolver_binary_node_t -------------------------------------------
    !
    !> a node in the binary tree
    !!

       ! [+] Type : ftsolver_domains_t ----------------------------------------
       !
       !> data used in the analysis phase to describes
       !! all the domains.
       !!
       Type ftsolver_domains_t
         sequence

         !> number of domains
         Integer*FTS_INTKIND :: nbdom 

         !-------------------------------
         !> @par [ interior description ] 
         !-------------------------------

         !> nnz in interiors (size nbdom)
         FTS_INT, Pointer :: domintnnz  (:) 
         !> dof of each interior (size nbdom)
         Integer*FTS_INTKIND , Pointer :: domintdof (:)
         !> Pointer to partitioning tree's leafs
         Integer*FTS_INTKIND, Pointer :: domstptr (:)  

         !-------------------------------
         !> @par [ interface description ] 
         !-------------------------------

         !> nnz in interfaces (size nbdom)
         FTS_INT, Pointer :: domintrfnnz(:)
         !> nnz on all interfaces
         FTS_INT :: nnzallinterface     
         !> degree of freedom of each domain's interface
         Integer*FTS_INTKIND, Pointer  :: domintrfdof (:)
         !> weight (lty : ? which one ?) on vertex(size totinterface)
         Integer*FTS_INTKIND, Pointer :: vtxweight (:) 
         !> lty :? vertex position ? (size totinterface+1)
         Integer*FTS_INTKIND, Pointer :: vtxpos    (:) 

         ! To construct each domain, we need to identify 
         ! the interfaces by analyzing the binary tree.
         ! To each node "i" on an interface, we associate 
         ! a unique identifier and its domain id,
         ! which are respectively stored into 
         ! 'intrindices(i)' and 'intrfproc(i)'
         ! (range for i is 1..combivtxpcsz) 
         
         !> lty: ? sum on interfaces of "myintrfndof" ?
         Integer*FTS_INTKIND :: combivtxpcsz               
         !> interface's vertex index  (size >= combivtxpcsz) 
         Integer*FTS_INTKIND, Pointer  :: intrfindices (:) 
         !> interface's vertex domain (size >= combivtxpcsz) 
         Integer*FTS_INTKIND, Pointer  :: intrfproc    (:) 

         !-------------------------------
         !> @par [ statistics ] 
         !-------------------------------

         !> maximum of interfaces' degrees of freedom 
         Integer*FTS_INTKIND :: maxprocintrf               
         !> minimum of interfaces' degrees of freedom 
         Integer*FTS_INTKIND :: minprocintrf  
         !> see ftsolver_binary_tree%totinterface             
         Integer*FTS_INTKIND :: totinterface  

         !-------------------------------
         !> @par [ other ] 
         !-------------------------------


         !> permutation on rows/columns
         Integer*FTS_INTKIND, Pointer :: metperm  (:) 

      End Type ftsolver_domains_t





    !> defines the data on each domain
    Type describe_local_t
       sequence
       !--------------------------------------------------------------------
       !> @par [ local matrix A_delta0, A_delta1, and A_delta2  description ]
       !--------------------------------------------------------------------

       !> number of rows (columns) in the local matrix A_delta0
       Integer           :: size_domain

       !> number of rows (columns) in the local matrix A_delta0
       Integer           :: size_periph
       
       !> vertices (rows) in the local matrix A_delta0
       Integer,  pointer :: index_domain     (:)


       !------------------------------------------------
       !> @par [ distance one' neighbors  description  ]
       !------------------------------------------------
       
       !> number of distance one' neighbors 
       Integer          :: nb_neig

       !> index of distance one' neighbors
       Integer, pointer :: index_neig  (:)
       
       !> direct pointer on index_neig_edge
       Integer, pointer :: ptr_direct_neig_edge (:)
       
       !> pointer on  vertices of index_vertperneig (size: nb_d1_neig +1)
       !> this data structure contains the same vertices stored 
       !> index_periph_vert. However, index_vertperneig is different
       !> because for each neighbor(domain) , it store successively the vertices 
       !> it is linked to. So if a vertex is linked to many neighbors
       !> it replicated 
       Integer, pointer :: ptr_vertperneig (:)

       !> size of  index_vertperneig
       Integer          :: size_vertperneig

       !> contains for each
       ! neighbor domain, the vertices of the periphery 
       ! of A_delta0 they are linked to. 
       Integer, pointer :: index_vertperneig (:)

       !> Pointer on edges at distance one
       !> useful if some vertices in index_vertperneig
       !> have many edges with a given domain 
       Integer, pointer :: ptr_neig_edge      (:) 

       !> size of  index_d1_edge
       Integer          :: size_neig_edge

       !> contains end vertices of edges at distance one
       Integer, pointer :: index_neig_edge    (:)

       !------------------------------------------------
       !> @par [ distance two' neighbors  description  ]
       !------------------------------------------------

       !> number of distance two' neighbors(domain)
       Integer          :: nb_sec_neig

       !> index of distance two' neighbors
       Integer, pointer :: index_sec_neig     (:)

       !> pointer on index_d2_edge
       Integer, pointer :: ptr_sec_neig_edge       (:)

       !> size of index_d2_edge
       Integer          :: size_sec_neig_edge 

       !> contains vertices at distance two
       Integer, pointer :: index_sec_neig_edge     (:)

    end Type describe_local_t

    End Module FTS_part_Type
