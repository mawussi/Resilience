!This Module compute vertices at distance 1 and 2 
! for a given domain/partition
#include"fts_defs_f.h"
      
Module FTS_get_neighbor_mod

      Use FTS_part_type
      Implicit None 


 
  ! [+] type : FTS_queue_t ----------------
  !!
  !!
  !>  queuetab  array to simulate the queue
  !! 
  !>  queuebeg   index of the beginning of the queue
  !!
  !>  queueend   index of the end of the queue
  !!-----------------
  !! @author Mawussi Zounon
  !!
      Type  FTS_queue_t; Sequence

      Integer, dimension(:), ALLOCATABLE  :: queuetab
      Integer                             :: queuebeg 
      Integer                             :: queueend
      End Type FTS_queue_t


      Public  :: FTS_distance  
      Private :: FTS_queue_init
      Private :: FTS_queue_insert
      Private :: FTS_queue_top_and_pop
      Private :: FTS_notin
      Private :: FTS_quick_sort
      Contains



  ! [+] routine : FTS_queue_init ----------------
  !
  !> Initialize a queue
  !  
  !!
  !!-------
  !! @param[in,out ]  queue_t           queue
  !! @param[ in]      max_size_queue    maximal size of the queue
  !!
  !!-------
  !! @author Mawussi Zounon
  !!
      Subroutine FTS_queue_init(& ! intents
      queue_t,                  & ! inout
      max_size_queue            & ! in
      )

      Type(FTS_queue_t) :: queue_t      
      Integer         :: vertnbr
      Integer         :: max_size_queue
      
      !allocate memory for queue_t
      allocate(queue_t%queuetab(max_size_queue) )
  
      queue_t%queuebeg = 1
      queue_t%queueend = 1
      
      End Subroutine FTS_queue_init

  ! [+] routine : FTS_queue_insert ----------------
  !
  !> Insert new element in the queue 
  !  
  !!
  !!-------
  !! @param[in,out ]  queue_t           queue
  !! @param[ in]     elem              new element
  !!
  !!-------
  !! @author Mawussi Zounon
  !!

      Subroutine FTS_queue_insert(& ! intents
      queue_t,                    & ! inout
      elem                        & ! in
      )
      
      Type(FTS_queue_t) :: queue_t
      Integer         :: elem
      
      queue_t%queuetab(queue_t%queueend) = elem
      queue_t%queueend = (queue_t%queueend+1)

      End Subroutine FTS_queue_insert




  ! [+] function : FTS_notin ----------------
  !
  !> Check if an array contains a given value
  ! 
  !!
  !!-------
  !! @param[in]     myarray               target array
  !! @param[in]     arraysize             size of the array
  !! @param[in]     maxarraysize           max  size of the array
  !! @param[in ]    elem             the value to search in the array
  !! @param[out ]   itcontains
  !!-------
  !! @author Mawussi Zounon

      Function FTS_notin(&
      myarray,&
      arraysize,&
      elem &
      )&
      result(itcontains)
      
      Integer, Intent(in)  :: myarray(:)
      Integer, Intent(in)  :: arraysize
      Integer, Intent(in)  :: elem
      Logical                :: itcontains

      !Local variables 
      Integer  :: iter_array

      itcontains = .true.
      
      If (arraysize .EQ. 0) return
      
      Do iter_array=1,arraysize,1
         If(myarray(iter_array) .EQ. elem) Then
            itcontains = .false.
            Exit
         End IF
      End Do

      End Function FTS_notin
      


  ! [+] routine : FTS_queue_top_and_pop ----------------
  !
  !> Removes the next element in the queue, 
  !  effectively reducing its size by one.
  !!
  !!-------
  !! @param[in,out ]  queue_t           queue
  !! @param[ out]     elem              next element
  !!
  !!-------
  !! @author Mawussi Zounon
  !!
      Subroutine FTS_queue_top_and_pop (& ! intents 
      queue_t,                          & ! inout
      elem                              & ! out
      )
      Implicit None 

      Type(FTS_queue_t) :: queue_t
      Integer         :: elem


      elem = queue_t%queuetab(queue_t%queuebeg)

      queue_t%queuebeg = ( queue_t%queuebeg +1)

      End subroutine FTS_queue_top_and_pop


  ! [+] routine : FTS_distance ----------------
  !
  !> Search all vertices at distance 1 or 2 
  ! of a given domain/partition 
  !  
  !!
  !!-------
  !! @param[in ]      vertnbr        number of vertices
  !! @param[in ]      edgenbr        number of edges               
  !! @param[in ]      verttab        array of start indice in edgetab of  vertices
  !! @param[in ]      edgetab        array of edges
  !! @param[in ]      parttab        array of partition
  !! @param[in ]      idpart         id of the target domain
  !! @param[inout ]   disttab        array of distance 
  !! @param[inout ]   vertinterftab     array of interface vertices
  !! @param[inout ]   vertdist1nbr   number interface vertices
  !! @param[inout ]   edgedist1nbr   number vertices at distance one
  !! @param[inout ]   vertdist1tab    array of start indice in ngedgetab of distance one vertices 
  !! @param[inout ]   edgedist1tab    array of egde at distance on
  !!-------
  !! @author Mawussi Zounon
  !!


      Subroutine FTS_distance (&
      vertnbr,&                 ! number of vertex
      edgenbr,&                 ! number of edges
      verttab,&                 ! array of start indice in edgetab of  vertices
      edgetab,&                 ! array of edges
      parttab,&                 ! partitioning array 
      idpart,&                  ! id of the target  domain 
      describe_lc,&             ! neighbor information
      info)                     ! info


      !routine argumentes
      Integer, Intent(in)          :: vertnbr
      Integer, Intent(in)          :: edgenbr
      Integer, Intent(in), pointer :: verttab(:)
      Integer, Intent(in), pointer :: edgetab(:)
      Integer, Intent(in), pointer :: parttab (:)
      Integer, Intent(in)          :: idpart
      Type(describe_local_t), Intent(inout)  :: describe_lc
      Integer, Intent(inout)       :: info
      
      
      !Parameters 
      Integer, parameter :: Max_distance = 2


      !Local variable
      Logical         :: have_neighbor
      Integer         :: i, j ,k
      Integer         :: vertnum, edgenum, vertend
      Integer         :: qsize, edgebeg, edgeend, part 
      Integer         :: iter_vertnum,  iter_edge, iter_vertend, iter_part
      Integer         :: vertinterfnbr, tmp_vertinterfnbr
      Integer         :: current_size_domain, size_domain
      Integer         :: nb_neig, maxnb_neig, size_neig_edge
      Integer         :: nb_sec_neig, tmp_size_sec_neig_edge, size_sec_neig_edge
      Integer         :: neigproc, neig2proc

      Integer, dimension(:), allocatable   :: index_neig, index_sec_neig
      Integer, pointer  :: disttab(:), order(:),  index_domain(:)
      Integer, pointer  :: vertinterftab(:), ptr_vertinterftab(:), tmp_vertinterftab(:)
      Integer, pointer  :: index_neig_edge(:), ptr_neig(:), tmp_index_neig_edge(:), ptr_direct_neig_edge(:)
      Integer, pointer  :: index_sec_neig_edge(:), ptr_sec_neig_edge(:), tmp_index_sec_neig_edge(:)

      !Initialization of queue_t
      Type(FTS_queue_t) :: queue_t
      CALL FTS_queue_init(queue_t, vertnbr)
      
      !Initialization 
      maxnb_neig                     = maxval(parttab) 
      part                           = maxnb_neig
      size_domain                    = 0
      nb_neig                        = 0
      tmp_vertinterfnbr              = 0
      nb_sec_neig                    = 0
      tmp_size_sec_neig_edge         = 0
      size_sec_neig_edge             = 0
      
      Nullify(disttab, order, index_domain)
      Nullify(vertinterftab, tmp_vertinterftab) 
      Nullify(ptr_vertinterftab, ptr_direct_neig_edge)
      Nullify(index_neig_edge, tmp_index_neig_edge) 
      Nullify(index_sec_neig_edge, tmp_index_sec_neig_edge)
      Nullify(ptr_sec_neig_edge)
      ! Allocation 
      allocate(disttab(vertnbr))
      allocate(order(edgenbr))
      allocate(index_domain(vertnbr))
      allocate(vertinterftab(vertnbr))
      allocate(tmp_vertinterftab(vertnbr))
      allocate(ptr_vertinterftab(vertnbr))
      allocate(index_neig(maxnb_neig))
      allocate(ptr_direct_neig_edge(vertnbr))
      allocate(index_neig_edge(edgenbr))
      allocate(tmp_index_neig_edge(edgenbr))
      allocate(index_sec_neig(maxnb_neig+1))
      allocate(index_sec_neig_edge(edgenbr))
      allocate(tmp_index_sec_neig_edge(edgenbr))
      allocate(ptr_sec_neig_edge(vertnbr))

      ptr_vertinterftab(1)           = 1
      ptr_sec_neig_edge(1)           = 1

      ! Initialization of distance array
      Do vertnum=1, vertnbr, 1
         If (parttab(vertnum) .eq. idpart) Then
          disttab (vertnum) = 0
          !increment size_domain
          size_domain = size_domain+1
          index_domain(size_domain) = vertnum
          CALL FTS_queue_insert(queue_t, vertnum)
         Else
            disttab(vertnum) = -1
         End If
      End Do

      Do While( modulo((vertnbr + queue_t%queueend - queue_t%queuebeg),&
         vertnbr) > 0)
         qsize = modulo((vertnbr + queue_t%queueend - queue_t%queuebeg),&
         vertnbr)
         CALL FTS_queue_top_and_pop (queue_t, vertnum)
         If (disttab(vertnum) .Ge. Max_distance) Exit

         !Iterate through neigbours of vertnum! 
         have_neighbor = .false.

         edgebeg = verttab(vertnum)
         edgeend = verttab(vertnum+1)-1
         Do edgenum=edgebeg,edgeend,1
            vertend = edgetab(edgenum)

            If ((disttab(vertnum) .EQ. 0) .And. disttab(vertend) .Ne. 0) Then
               have_neighbor = .true.
            End If
            If (disttab(vertend) .EQ. -1) Then
               disttab(vertend) = disttab(vertnum)+1
               CALL FTS_queue_insert(queue_t, vertend)
               
               ! Add vertend to find vertex at distance 1
               If (disttab(vertend) .EQ. 1) Then
               !If the edge belongs to a new neighbor proc then 
               ! increment number of neighbor and add the proc(partition)
               ! to the array of neighbor proc
                  neigproc = parttab(vertend)
                  If ((nb_neig .EQ. 0)  .or. (FTS_notin(index_neig, nb_neig, neigproc) )) Then
                     nb_neig = nb_neig+1
                     index_neig(nb_neig) = neigproc
                  End If
               
                  !Add neighbor at distance 2
               ElseIf(disttab(vertend) .EQ. 2) Then
                  neig2proc = parttab(vertend)
                  If ((nb_sec_neig .EQ. 0) .or. (FTS_notin(index_sec_neig, nb_sec_neig, neig2proc))) Then 
                     nb_sec_neig = nb_sec_neig+1
                     index_sec_neig(nb_sec_neig) = neig2proc
                  End If
                  tmp_size_sec_neig_edge = tmp_size_sec_neig_edge +1
                  tmp_index_sec_neig_edge(tmp_size_sec_neig_edge) = vertend
               End If
            End If
         End Do
         If ( have_neighbor) Then
            tmp_vertinterfnbr               = tmp_vertinterfnbr +1
            tmp_vertinterftab(tmp_vertinterfnbr)   = vertnum
         End If
      End Do
      
      ! Distribute neighbor vertices dist1 per neighbor proc (partition) 
      vertinterfnbr = 0
      size_neig_edge  = 0
      ptr_direct_neig_edge(1) = 1
      Do neigproc =1,nb_neig,1
         Do iter_vertnum =1,tmp_vertinterfnbr
            vertnum = tmp_vertinterftab(iter_vertnum)
            edgebeg = verttab(vertnum)
            edgeend = verttab(vertnum+1)-1
            have_neighbor = .false.
            Do edgenum=edgebeg,edgeend,1
               vertend = edgetab(edgenum)
               If (parttab(vertend) .EQ. index_neig(neigproc)) Then
                  have_neighbor = .true.
                  If ( FTS_notin(index_neig_edge,size_neig_edge, vertend))Then
                     size_neig_edge = size_neig_edge + 1
                     index_neig_edge(size_neig_edge) = vertend
                  End If
               End If
            End DO
            If (have_neighbor) Then
               vertinterfnbr = vertinterfnbr +1               
               vertinterftab(vertinterfnbr) = vertnum
            End If
            
         End Do
         ptr_vertinterftab(neigproc+1) = vertinterfnbr+1
         ptr_direct_neig_edge(neigproc+1) = size_neig_edge +1
         Nullify(ptr_neig)
         ptr_neig => index_neig_edge(ptr_direct_neig_edge(neigproc):&
              ptr_direct_neig_edge(neigproc+1)-1)
       Call FTS_quick_sort(ptr_neig, order)
      End Do


      ! Distribute neighbor vertices dist2 per neighbor proc (partition) 
      Do neig2proc =1,nb_sec_neig
         
         Do iter_edge=1,tmp_size_sec_neig_edge
            edgenum = tmp_index_sec_neig_edge(iter_edge)

            If (parttab(edgenum) .EQ. index_sec_neig(neig2proc)) Then
               size_sec_neig_edge = size_sec_neig_edge +1              
               index_sec_neig_edge(size_sec_neig_edge) = edgenum
            End If
            ptr_sec_neig_edge(neig2proc+1) = size_sec_neig_edge +1
         End Do
      End Do


      ! fill describe_lc
      describe_lc%size_domain = size_domain
      Allocate(describe_lc%index_domain(size_domain))

      describe_lc%size_periph = tmp_vertinterfnbr
      ! fill index_domain first with periphery indices 

      !fill index_domain with other indices of  idpart
      current_size_domain= 0
      Do i=1,size_domain 
         vertnum = index_domain(i)
         If (FTS_notin(describe_lc%index_domain, current_size_domain, vertnum))Then
            current_size_domain = current_size_domain+1
            describe_lc%index_domain(current_size_domain) = vertnum
         End If
      End Do

      !save neighboring describtion
       describe_lc%nb_neig = nb_neig
       Allocate(describe_lc%index_neig(nb_neig))
       
       Do i=1,nb_neig
          describe_lc%index_neig(i) = index_neig(i)
       End Do
  
       Allocate(describe_lc%ptr_direct_neig_edge(nb_neig+1))
       Do i=1,nb_neig
          describe_lc%ptr_direct_neig_edge(i) = ptr_direct_neig_edge(i)
       End Do
       describe_lc%ptr_direct_neig_edge(nb_neig+1) = ptr_direct_neig_edge(nb_neig+1)
       
       Allocate(describe_lc%ptr_vertperneig(nb_neig+1))

      Do i=1, nb_neig
         describe_lc%ptr_vertperneig(i) = ptr_vertinterftab(i)
      End Do
      describe_lc%ptr_vertperneig(nb_neig+1) = ptr_vertinterftab(nb_neig+1)

      describe_lc%size_vertperneig = vertinterfnbr
      Allocate(describe_lc%index_vertperneig(vertinterfnbr))
      
      Do i=1,vertinterfnbr
         describe_lc%index_vertperneig(i) = vertinterftab(i)
      End Do

      describe_lc%size_neig_edge = size_neig_edge
      Allocate(describe_lc%index_neig_edge(size_neig_edge))
      
      Do i=1,size_neig_edge
         describe_lc%index_neig_edge(i) = index_neig_edge(i)
      End Do
       !------------------------------------------------
       ! [ distance two' neighbors  description        ]
       !------------------------------------------------


      !save distance two neighboring  indices


      describe_lc%nb_sec_neig = nb_sec_neig
      Allocate(describe_lc%index_sec_neig(nb_sec_neig))
      
      Do i=1,nb_sec_neig
         describe_lc%index_sec_neig(i) = index_sec_neig(i)
      End Do
      
      Allocate(describe_lc%ptr_sec_neig_edge(nb_sec_neig+1))
      Do i=1,nb_sec_neig
         describe_lc%ptr_sec_neig_edge(i) = ptr_sec_neig_edge(i)
      End Do
      describe_lc%ptr_sec_neig_edge(nb_sec_neig+1) =&
      ptr_sec_neig_edge(nb_sec_neig+1)

      describe_lc%size_sec_neig_edge = size_sec_neig_edge
      Allocate(describe_lc%index_sec_neig_edge(size_sec_neig_edge))
      
      
      Do i=1,size_sec_neig_edge
         describe_lc%index_sec_neig_edge(i) = index_sec_neig_edge(i)
      End Do

      ! End deallocate 
      If (Associated (disttab))  Deallocate(disttab)
      If (Associated (order)) Deallocate(order)
      If (Associated (vertinterftab)) Deallocate (vertinterftab)
      If (Associated (tmp_vertinterftab)) Deallocate( tmp_vertinterftab)
      If (Associated (index_domain)) Deallocate(index_domain)
      If (Associated (ptr_vertinterftab)) Deallocate (ptr_vertinterftab)

      If (Associated ( index_neig_edge))  Deallocate (index_neig_edge) 
      If (Associated (tmp_index_neig_edge)) Deallocate(tmp_index_neig_edge)
      If (Associated (ptr_direct_neig_edge)) Deallocate (ptr_direct_neig_edge)
      If (Associated ( index_sec_neig_edge))  Deallocate (index_sec_neig_edge) 
      If (Associated (tmp_index_sec_neig_edge)) Deallocate (tmp_index_sec_neig_edge)
      If (Associated (ptr_sec_neig_edge))   Deallocate(ptr_sec_neig_edge)

      End Subroutine FTS_distance 



      RECURSIVE SUBROUTINE FTS_quick_sort(list, order)
        
        IMPLICIT NONE
        INTEGER, DIMENSION (:), INTENT(IN OUT)  :: list
        INTEGER, DIMENSION (:), INTENT(OUT)  :: order
  
        ! Local variable
        INTEGER :: i
        
        DO i = 1, SIZE(list)
           order(i) = i
        END DO
  
        CALL quick_sort_1(1, SIZE(list))
        
      CONTAINS
        
        RECURSIVE SUBROUTINE quick_sort_1(left_end, right_end)
          
          INTEGER, INTENT(IN) :: left_end, right_end
          
          !     Local variables
          INTEGER             :: i, j, itemp
          INTEGER                :: reference, temp
          INTEGER, PARAMETER  :: max_simple_sort_size = 6
          
          IF (right_end < left_end + max_simple_sort_size) THEN
             ! Use interchange sort for small lists
             CALL interchange_sort(left_end, right_end)
             
          ELSE
             ! Use partition ("quick") sort
             reference = list((left_end + right_end)/2)
             i = left_end - 1; j = right_end + 1
             
             DO
                ! Scan list from left end until element >= reference is found
                DO
                   i = i + 1
                   IF (list(i) >= reference) EXIT
                END DO
                ! Scan list from right end until element <= reference is found
                DO
                   j = j - 1
                   IF (list(j) <= reference) EXIT
                END DO
                
                
                IF (i < j) THEN
                   ! Swap two out-of-order elements
                   temp = list(i); list(i) = list(j); list(j) = temp
                   itemp = order(i); order(i) = order(j); order(j) = itemp
                ELSE IF (i == j) THEN
                   i = i + 1
                   EXIT
                ELSE
                   EXIT
                END IF
             END DO
             
             IF (left_end < j) CALL quick_sort_1(left_end, j)
             IF (i < right_end) CALL quick_sort_1(i, right_end)
          END IF
          
        END SUBROUTINE quick_sort_1
        
        
        SUBROUTINE interchange_sort(left_end, right_end)
    
          INTEGER, INTENT(IN) :: left_end, right_end
          
          !     Local variables
          INTEGER             :: i, j, itemp
          INTEGER                :: temp
          
          DO i = left_end, right_end - 1
             DO j = i+1, right_end
                IF (list(i) > list(j)) THEN
                   temp = list(i); list(i) = list(j); list(j) = temp
                   itemp = order(i); order(i) = order(j); order(j) = itemp
                END IF
             END DO
          END DO
          
        END SUBROUTINE interchange_sort
        
      END SUBROUTINE FTS_quick_sort

End Module FTS_get_neighbor_mod
