! Warning: XFTS_GENFILE_COMMENT
#include "fts_defs_f.h"
#include "fts_macros_f.h"
!! module for fault pluging.
!!
!! 
Module XFTS_fault_mod

  !* Modules * !
  Use fts_error_mod
  Implicit None
  !* Private constants *!
  Character(len=FTSOLVER_STRL), Private, Parameter :: &
       FLNAME = "xfts_fault_mod.F90"

  !* List of routines *!
  Public   :: XFTS_CHECK_FAULT
  Public   :: XFTS_ISFAULT
Contains
  
  Subroutine XFTS_Check_Fault( &! intents
       ftsl,                   & ! in
       work                   & ! in
       )
    
    !*modules *!
    Use XFTS_ftsolver_type
    Use FTS_ftsolver_enum
    Use XFTS_dense_matrix_mod  
    Use FTS_time_mod
!!DBG
    Implicit None 
    Include 'mpif.h'
    !* Arguments *!
    Type(XFTS_ftsolver_t)   , Intent(inout) :: ftsl
    XFTS_FLOAT            :: work(*)

    !* Local variables *!
    Integer      :: i, FT_strategy 
    Integer      :: Diff
    Integer      :: myrank
    Integer      :: bound_ndof
    Integer      :: restrt ! restart parameter
    Integer      :: info
    Integer      :: ArrayOfFaultyPorc(3)
    Integer      :: FailedNeighbs(3)
    Integer      :: intrf_ndof
    Integer      :: NbOfFailedNeighb
    Integer      :: NbOfFaultyProc
    Integer      :: Iamfailed
    XFTS_FLOAT   :: test

    Type(XFTS_dense_matrix_t)  :: bound_rhs 
    Type(XFTS_dense_matrix_t)  :: x
    !* initialize variables*!
    Diff = 0
    Iamfailed = 0
    FT_strategy = ftsl%icntl(29)
    myrank= ftsl%ikeep(IKEEP_MPIRANK)
    bound_ndof    = ftsl%lc_domain%size_domain

    
    !---------------------------------------------------------
    ! [1] Fault detection
    !---------------------------------------------------------
    ! [1.1] Timers
    !sett initialize timer
    Call FTS_time_start(ftsl%rinfo (RINFO_TIMING_Fault))
    Call XFTS_GenerateFault(ftsl,ArrayOfFaultyPorc,NbOfFaultyProc,Iamfailed)
    Call FTS_time_stop(ftsl%rinfo (RINFO_TIMING_Fault))
    Call FTS_time_start(ftsl%rinfo (RINFO_TIMING_Interp))

    !---------------------------------------------------------
    ! [2] Fault recovery
    !---------------------------------------------------------

    If (NbOfFaultyProc .Eq. 1)  Then 
       !---------------------------------------------------------
       ! [2.1] single fault case
       !---------------------------------------------------------
       If (Iamfailed .Eq. 1) Then 
          !---------------------------------------------------------
          ! [2.1.a] set  workspace to zero
          !---------------------------------------------------------

          If(FT_strategy .Le. 2) Then
             Do i =1, bound_ndof
                work(i) = XFTS_FLOATZERO
             End Do
          End If
          !---------------------------------------------------------
          ! [2.1.b] Receive contribution  from neighbors and update rhs
          !---------------------------------------------------------   
          
          If (FT_strategy .EQ. 1) Then ! LI
             
             Call FTS_time_start(ftsl%rinfo (RINFO_TIMING_InterpRhs))
             Call XFTS_FRecv_LI(ftsl, work, info)
             Call FTS_time_stop(ftsl%rinfo (RINFO_TIMING_InterpRhs))

             Call FTS_time_start(ftsl%rinfo (RINFO_TIMING_InterpSlv))
             CAll XFTS_InterpLI(ftsl, work, info)
             Call FTS_time_stop(ftsl%rinfo (RINFO_TIMING_InterpSlv))

          Else If (FT_strategy .EQ. 2) Then  !LSI
             Call FTS_time_start(ftsl%rinfo (RINFO_TIMING_InterpRhs))
             Call XFTS_FRecv_LSI(ftsl, work, info)
             Call FTS_time_stop(ftsl%rinfo (RINFO_TIMING_InterpRhs))

             Call FTS_time_start(ftsl%rinfo (RINFO_TIMING_InterpSlv))
             CAll XFTS_InterpLSI(ftsl, work, info)
             Call FTS_time_stop(ftsl%rinfo (RINFO_TIMING_InterpSlv))
             
          End If
       Else
          !---------------------------------------------------------
          ! [2.1.c] Send data to failed neighbor
          !---------------------------------------------------------           
          If (FT_strategy .EQ. 1) Then ! LI 
             
             Call XFTS_IsNeighbor(ftsl, ArrayOfFaultyPorc, NbOfFaultyProc,&
                  FailedNeighbs, NbOfFailedNeighb,info)
             If(info == 1 ) Then  ! The faulty processor is my neighbor
                Call XFTS_FSend_LI(ftsl, work, FailedNeighbs,NbOfFailedNeighb, info)           
             End If
          End If
          
          If (FT_strategy .EQ. 2) Then  !LSI
             
             Call XFTS_IsNeighbor(ftsl, ArrayOfFaultyPorc, NbOfFaultyProc,&
                  FailedNeighbs, NbOfFailedNeighb,info)
             If(info == 1 ) Then  ! The faulty processor is my neighbor
                Call XFTS_FSend_Dist1(ftsl, work, FailedNeighbs,NbOfFailedNeighb, info)           
             End If

             ! Distance two neighbors
             Call XFTS_IsSecNeighbor(ftsl, ArrayOfFaultyPorc, NbOfFaultyProc,&
                  FailedNeighbs, NbOfFailedNeighb,info)
             If(info == 1 ) Then  ! The faulty processor is my neighbor
                Call XFTS_FSend_Dist2(ftsl, work, FailedNeighbs,NbOfFailedNeighb, info)           
             End If
          End If
       End If
    End if
    
    ![--] timers                                                                                                                   
    Call FTS_time_stop(ftsl%rinfo (RINFO_TIMING_Interp))
    !----------------------------------------------------------------------
    ! [5] Exit routine
    !----------------------------------------------------------------------
    
9999 Continue
       
  END SUBROUTINE XFTS_Check_Fault
  
     !XFTS_ISFAULT return
     ! 0 if there is no fault
     !else the number of faults
     SUBROUTINE XFTS_ISFAULT(ftsl, iter,info)
       
       Use XFTS_ftsolver_type
       Use FTS_ftsolver_enum
       Implicit None
       !*Arguments*!
       Type(XFTS_ftsolver_t) ,Intent(in ) :: ftsl
       Integer             ,Intent(in ) :: iter
       Integer             ,Intent(inout) :: info

       Integer FT_taux
       Integer FT_MaxFault
       Integer nbfault
       FT_taux = ftsl%icntl(28)
       FT_MaxFault = ftsl%icntl(33)
       nbfault = ftsl%icntl(34)
       info = 0
       If((iter .eq. FT_taux) .and. (ftsl%icntl(29) < 4)) info = nbfault
       
     END SUBROUTINE XFTS_ISFAULT
     

     SUBROUTINE XFTS_GenerateFault(ftsl,ArrayOfFaultyPorc, NbOfFaultyProc,Iamfailed)

       
       Use XFTS_ftsolver_type
       Use FTS_ftsolver_enum
       
       Implicit None
       Include'mpif.h'
       
       !Arguments *!
       Type(XFTS_ftsolver_t) ,Intent(in ) :: ftsl
       Integer             ,Intent(inout) :: ArrayOfFaultyPorc(3)
       Integer             ,Intent(out )  :: nbOfFaultyProc
       Integer             ,Intent(inout ):: Iamfailed
       !local variables*!
       Integer i, j, nNeighb
       Integer myrank, root
       Integer nbproc
       Integer minCriterion
       Integer candidate
       Integer info
       Integer Faulty(2)
       Integer nbfault
       Integer k
       Integer Min_k
       Integer Best_candidate
       Real    rand_number
       Integer, Pointer      :: failedproc(:)
       Integer, Pointer      :: failedcandidate(:)
       Integer, Pointer :: indexVi       (:)
       
       nbproc= ftsl%ikeep(IKEEP_NBDOMAINS)
       myrank= ftsl%ikeep(IKEEP_MPIRANK)
       nbfault = ftsl%icntl(34)
       nNeighb = ftsl%lc_domain%nb_neig
       indexVi => ftsl%lc_domain%index_neig
       root = 0
       candidate = 0
       Min_k=10000000
       !* memory allocation
       Allocate(failedproc(nbproc), STAT = info)
       CHCKASSRT( info == 0, info )
       If( info < 0 ) Goto 9999
       
       !* memory allocation
       Allocate(failedcandidate(nbproc), STAT = info)
       CHCKASSRT( info == 0, info )
       If( info < 0 ) Goto 9999

       
       NbOfFaultyProc = 0
       !* All processor determine the one with less neighbors 
       Call MPI_Allreduce(nNeighb, minCriterion, 1, MPI_INT, MPI_MIN, ftsl%comm, info)
       
       If (minCriterion .Eq. nNeighb) Then
          candidate = 1
       End If
       
       Call MPI_Gather(candidate, 1, MPI_INT, &
            failedcandidate, 1, MPI_INT, root, ftsl%comm, info); 
       ASSRT( info == MPI_SUCCESS )

       !* The proc of rank 0 choose the first faulty 
       !* Proc ramdomly
       If (myrank .Eq. 0) Then
          ! Select the processor with smallest rank
          Do i = 1, nbproc
             If (failedcandidate(i) .Eq. 1) Then
                Faulty(1) = i-1
                !Faulty(1) = nbproc/2
             End If
          End Do
       End If
       !* Broadcast the id of the faulty proc 
       Call MPI_Bcast(Faulty(1), 1, MPI_INT, root, ftsl%comm, info)
       ASSRT( info == MPI_SUCCESS )
       !*Each proc If he is the failed proc
       If( Faulty(1) .eq. myrank) Then
          Iamfailed = 1
       End If
       
       ! Each processor send Iamfailed
       ! to others 
       Call MPI_Allgather(Iamfailed, 1, MPI_INT, &
            failedproc, 1, MPI_INT, ftsl%comm, info); 
       ASSRT( info == MPI_SUCCESS )
       
       
       ! count number of fault
       Do i = 1, nbproc
          If (failedproc(i) .Eq. 1) Then
             NbOfFaultyProc = NbOfFaultyProc+1
             ArrayOfFaultyPorc(NbOfFaultyProc) = i-1
          End If
       End Do
       
9999 continue
    !Free memory
    If (Associated(failedproc)) Deallocate(failedproc)

  End SUBROUTINE XFTS_GENERATEFAULT


  !* FTS_IsNeighbor check if the faulty proc
  !* is my neighbor
  SUBROUTINE XFTS_IsNeighbor( & ! intents
       ftsl,                 & ! in 
       ArrayOfFaultyPorc,    & ! in
       NbOfFaultyProc,       & ! in
       FailedNeighbs,        & ! in
       NbOfFailedNeighb,     & ! in
       info                  & ! out 
       )
       
    !* Module(s) *!
    Use XFTS_ftsolver_type
    Use fts_error_mod        
    Use FTS_ftsolver_enum
    
    Implicit None
    Include 'mpif.h'
    
    !* Arguments *!
    Type(XFTS_ftsolver_t)   , Intent(in) :: ftsl
    Integer               , Intent(in  )  :: ArrayOfFaultyPorc(3)
    Integer               , Intent(in  )  :: NbOfFaultyProc
    Integer               , Intent(out  )  :: FailedNeighbs(3)
    Integer               , Intent(out )  :: NbOfFailedNeighb
    Integer               , Intent(out ) :: info    
    !* Local variables *!
    
    ! Scalars
    Integer :: nNeighb 
    Integer :: neighb  
    Integer :: failedrank
    Integer :: i,myrank
    ! Arrays 
    Integer     , Pointer :: indexVi       (:)
    
    !initialization
    info = 0
    myrank= ftsl%ikeep(IKEEP_MPIRANK)
    nNeighb      = ftsl%lc_domain%nb_neig
    indexVi      => ftsl%lc_domain%index_neig
    NbOfFailedNeighb = 0
    Do i =1,NbOfFaultyProc
       failedrank = ArrayOfFaultyPorc(i) 
       Do neighb=1, nNeighb
          If ((failedrank .Eq. indexVi(neighb)) .And. (failedrank .Ne. myrank)) Then
             info = 1
             NbOfFailedNeighb = NbOfFailedNeighb +1
             FailedNeighbs(NbOfFailedNeighb) = failedrank
          End If
       End Do
    End DO
9999 Continue
      
  END SUBROUTINE XFTS_IsNeighbor
  

  !* FTS_IsSecNeighbor check if the faulty proc
  !* is my distance 2 neighbor
  SUBROUTINE XFTS_IsSecNeighbor( & ! intents
       ftsl,                 & ! in 
       ArrayOfFaultyPorc,    & ! in
       NbOfFaultyProc,       & ! in
       FailedNeighbs,        & ! in
       NbOfFailedNeighb,     & ! in
       info                  & ! out 
       )
       
    !* Module(s) *!
    Use XFTS_ftsolver_type
    Use fts_error_mod        
    Use FTS_ftsolver_enum
    
    Implicit None
    Include 'mpif.h'
    
    !* Arguments *!
    Type(XFTS_ftsolver_t)   , Intent(in) :: ftsl
    Integer               , Intent(in  )  :: ArrayOfFaultyPorc(3)
    Integer               , Intent(in  )  :: NbOfFaultyProc
    Integer               , Intent(out  )  :: FailedNeighbs(3)
    Integer               , Intent(out )  :: NbOfFailedNeighb
    Integer               , Intent(out ) :: info    
    !* Local variables *!
    
    ! Scalars
    Integer :: nb_secin_neig 
    Integer :: secin_neig  
    Integer :: failedrank
    Integer :: i,myrank
    ! Arrays 
    Integer     , Pointer :: index_secin_neig       (:)
    
    !initialization
    info = 0
    myrank= ftsl%ikeep(IKEEP_MPIRANK)
    nb_secin_neig      = ftsl%lc_domain%nb_secin_neig
    index_secin_neig      => ftsl%lc_domain%index_secin_neig
    NbOfFailedNeighb = 0
    Do i =1,NbOfFaultyProc
       failedrank = ArrayOfFaultyPorc(i) 
       Do secin_neig=1, nb_secin_neig
          If ((failedrank .Eq. index_secin_neig(secin_neig)) .And. (failedrank .Ne. myrank)) Then
             info = 1
             NbOfFailedNeighb = NbOfFailedNeighb +1
             FailedNeighbs(NbOfFailedNeighb) = failedrank
          End If
       End Do
    End DO
9999 Continue
      
  END SUBROUTINE XFTS_IsSecNeighbor


  SUBROUTINE   XFTS_FSend_LI( & ! intents
       ftsl,                 & ! inout
       work,                 & ! in
       FailedNeighbs,        & ! in
       NbOfFailedNeighb,     & ! in
       info                  & ! out
       )       
       
    !* Module(s) *!
    Use XFTS_ftsolver_type
    Use fts_error_mod        
    Use FTS_ftsolver_enum
    Use FTS_time_mod
    Implicit None
    Include 'mpif.h'
    
    !* Arguments *!
    Type(XFTS_ftsolver_t)    , Intent(inout) :: ftsl
    XFTS_FLOAT                             :: work(*)
    Integer                ,Intent(in  )   :: FailedNeighbs(3)
    Integer                ,Intent(in  )   :: NbOfFailedNeighb
    Integer                , Intent(out )  :: info
    
    !* Local variables *!
    ! Scalars
    Integer :: tagg      !< MPI tagg
    Integer :: comm, rank, failedrank      !< MPI Communicator, rank
    Integer :: iinfo     !< status of the called routines
    Integer :: i,k     !< dummy iterators
    Integer :: nb_neig      !< field of lc_intrf. see lc_intrf%mynbvi
    Integer :: neigh, indexFailed_Neighbor     !< neighbor
    Integer :: size_vertperneig         !> size of  index_vertperneig
    Integer :: size_neig_edge, neig_col
    Integer :: nloc, col , row, iter_col, iter_row 
    Integer :: mpistatus(MPI_STATUS_SIZE) 
    FTS_INT :: sbuffsize !< size of sbuff(:)
    FTS_INT :: rbuffsize !< size of rbuff(:)
    FTS_INT :: start 
    FTS_INT :: end
 
    ! Arrays


    Integer   , Pointer :: mpisreq(:)               !< MPI requests senders
    XFTS_FLOAT, Pointer :: sbuff (:)             !< MPI send Buffer
    Integer, Pointer    :: index_neig          (:)  !< field of boundary
    Integer, pointer :: ptr_vertperneig (:)
    Integer, pointer :: index_vertperneig (:)



    


    !- End of header -----------------------------------------------------------

    !---------------------------------------------------------------------------
    ! [1] Init
    !---------------------------------------------------------------------------
    
    ! status
    info = 0
    rank   = ftsl%ikeep(IKEEP_MPIRANK)
    ! mpi
    tagg             =  44
    comm             =  ftsl%comm

    ! interface description
    nb_neig             =  ftsl%lc_domain%nb_neig
    size_vertperneig    =  ftsl%lc_domain%size_vertperneig
    nloc                =  ftsl%lc_domain%size_domain
    index_neig          => ftsl%lc_domain%index_neig
    ptr_vertperneig     => ftsl%lc_domain%ptr_vertperneig
    index_vertperneig   => ftsl%lc_domain%index_vertperneig
    


    ! values


    !---------------------------------------------------------------------------
    ! [1] Allocate the buffers
    !---------------------------------------------------------------------------
    
    Nullify( sbuff)
    Nullify( mpisreq )

    Allocate(mpisreq(nb_neig), STAT= info )
    CHCKASSRT( info == 0, info )
    If (info < 0) Goto 9999



    Do i=1,sbuffsize 
       sbuff(i) = -999999999
    End Do


    !---------------------------------------------------------------------------
    ! [2] Send the computed value on the interface
    !---------------------------------------------------------------------------

    failedrank = FailedNeighbs(1)

    ! search the local index of
    ! the failed neighbor 
       Do neigh=1, nb_neig
          If (failedrank .Eq. index_neig(neigh)) Then
             indexFailed_Neighbor = neigh
             Exit  
          End If
       End Do

       sbuffsize = ptr_vertperneig(indexFailed_Neighbor+1)-&
            ptr_vertperneig(indexFailed_Neighbor)
       Allocate( sbuff(sbuffsize), STAT= info )
       CHCKASSRT( info == 0, info )
       If (info < 0) Goto 9999

          ! Pack 
       start = ptr_vertperneig(indexFailed_Neighbor)
       end= ptr_vertperneig(indexFailed_Neighbor+1) - 1

       i=1
       Do k = start, end 
          sbuff(i) = work( index_vertperneig(k) ) 
          i=i+1
       End Do
       
       ! Send
       Call MPI_Isend (sbuff(1),(end-start+1),&
            XFTS_FLOATMPI, failedrank, tagg,comm,mpisreq(1),iinfo)
       ASSRT( iinfo  == MPI_SUCCESS )

       
       Call MPI_WaitAll(1, mpisreq, MPI_STATUSES_IGNORE, iinfo )
       ASSRT( iinfo  == MPI_SUCCESS )

    ! --------------------------------------------------------------------------
    !  Exit routine
    ! --------------------------------------------------------------------------

9999 Continue
    If (Associated(sbuff)) Deallocate(sbuff)
    If (Associated(mpisreq)) Deallocate(mpisreq)

End SUBROUTINE XFTS_FSend_LI



  SUBROUTINE   XFTS_FSend_Dist1( & ! intents
       ftsl,                 & ! inout
       work,                 & ! in
       FailedNeighbs,        & ! in
       NbOfFailedNeighb,     & ! in
       info                  & ! out
       )       
       
    !* Module(s) *!
    Use XFTS_ftsolver_type
    Use fts_error_mod        
    Use FTS_ftsolver_enum
    Use FTS_time_mod
    Implicit None
    Include 'mpif.h'
    
    !* Arguments *!
    Type(XFTS_ftsolver_t)    , Intent(inout) :: ftsl
    XFTS_FLOAT                             :: work(*)
    Integer                ,Intent(in  )   :: FailedNeighbs(3)
    Integer                ,Intent(in  )   :: NbOfFailedNeighb
    Integer                , Intent(out )  :: info
    
    !* Local variables *!
    ! Scalars
    Integer :: tagg      !< MPI tagg
    Integer :: comm, rank, failedrank      !< MPI Communicator, rank
    Integer :: iinfo     !< status of the called routines
    Integer :: i,k     !< dummy iterators
    Integer :: nb_neig      !< field of lc_intrf. see lc_intrf%mynbvi
    Integer :: neigh, indexFailed_Neighbor     !< neighbor
    Integer :: size_vertperneig         !> size of  index_vertperneig
    Integer :: size_neig_edge, neig_col
    Integer :: nloc, col , row, iter_col, iter_row 
    Integer :: mpistatus(MPI_STATUS_SIZE) 
    FTS_INT :: sbuffsize !< size of sbuff(:)
    FTS_INT :: rbuffsize !< size of rbuff(:)
    FTS_INT :: start 
    FTS_INT :: end
 
    ! Arrays

    Integer   , Pointer :: mpisreq(:)               !< MPI requests senders
    XFTS_FLOAT, Pointer :: sbuff (:)             !< MPI send Buffer
    Integer, Pointer    :: index_neig          (:)  !< field of boundary
    Integer, pointer :: ptr_vertperneig (:)
    Integer, pointer :: index_vertperneig (:)

    !- End of header -----------------------------------------------------------

    !---------------------------------------------------------------------------
    ! [1] Init
    !---------------------------------------------------------------------------
    
    ! status
    info = 0
    rank   = ftsl%ikeep(IKEEP_MPIRANK)
    ! mpi
    tagg             =  44
    comm             =  ftsl%comm

    ! interface description
    nb_neig             =  ftsl%lc_domain%nb_neig
    size_vertperneig    =  ftsl%lc_domain%size_vertperneig
    nloc                =  ftsl%lc_domain%size_domain
    index_neig          => ftsl%lc_domain%index_neig
    ptr_vertperneig     => ftsl%lc_domain%ptr_vertperneig
    index_vertperneig   => ftsl%lc_domain%index_vertperneig

    !---------------------------------------------------------------------------
    ! [1] Allocate the buffers
    !---------------------------------------------------------------------------
    
    Nullify( sbuff)
    Nullify( mpisreq )

    Allocate(mpisreq(nb_neig), STAT= info )
    CHCKASSRT( info == 0, info )
    If (info < 0) Goto 9999


    !---------------------------------------------------------------------------
    ! [2] Send x and rhs
    !---------------------------------------------------------------------------

    failedrank = FailedNeighbs(1)

    ! search the local index of
    ! the failed neighbor 
       Do neigh=1, nb_neig
          If (failedrank .Eq. index_neig(neigh)) Then
             indexFailed_Neighbor = neigh
             Exit  
          End If
       End Do

       sbuffsize = 2*(ptr_vertperneig(indexFailed_Neighbor+1)-&
            ptr_vertperneig(indexFailed_Neighbor))
       Allocate( sbuff(sbuffsize), STAT= info )
       CHCKASSRT( info == 0, info )
       If (info < 0) Goto 9999

       Do i=1,sbuffsize 
          sbuff(i) = -999999999
       End Do
       
       ! Pack 
       start = ptr_vertperneig(indexFailed_Neighbor)
       end= ptr_vertperneig(indexFailed_Neighbor+1) - 1
       

       i=1
       Do k = start, end 
          sbuff(i) = work( index_vertperneig(k) ) 
          sbuff((end-start+1)+i) = work(nloc+ index_vertperneig(k) ) 
          i=i+1
       End Do
       
       ! Send
       Call MPI_Isend (sbuff(1), sbuffsize,&
            XFTS_FLOATMPI, failedrank, tagg,comm,mpisreq(1),iinfo)
       ASSRT( iinfo  == MPI_SUCCESS )

       
       Call MPI_WaitAll(1, mpisreq, MPI_STATUSES_IGNORE, iinfo )
       ASSRT( iinfo  == MPI_SUCCESS )

    ! --------------------------------------------------------------------------
    !  Exit routine
    ! --------------------------------------------------------------------------

9999 Continue
    If (Associated(sbuff)) Deallocate(sbuff)
    If (Associated(mpisreq)) Deallocate(mpisreq)

End SUBROUTINE XFTS_FSend_Dist1


  SUBROUTINE   XFTS_FSend_Dist2( & ! intents
       ftsl,                 & ! inout
       work,                 & ! in
       FailedNeighbs,        & ! in
       NbOfFailedNeighb,     & ! in
       info                  & ! out
       )       
       
    !* Module(s) *!
    Use XFTS_ftsolver_type
    Use fts_error_mod        
    Use FTS_ftsolver_enum
    Use FTS_time_mod
    Implicit None
    Include 'mpif.h'
    
    !* Arguments *!
    Type(XFTS_ftsolver_t)    , Intent(inout) :: ftsl
    XFTS_FLOAT                             :: work(*)
    Integer                ,Intent(in  )   :: FailedNeighbs(3)
    Integer                ,Intent(in  )   :: NbOfFailedNeighb
    Integer                , Intent(out )  :: info
    
    !* Local variables *!
    ! Scalars
    Integer :: tagg      !< MPI tagg
    Integer :: comm, rank, failedrank      !< MPI Communicator, rank
    Integer :: iinfo     !< status of the called routines
    Integer :: i,k     !< dummy iterators
    Integer :: nb_secin_neig      !< field of lc_intrf. see lc_intrf%mynbvi
    Integer :: neigh, indexFailed_Neighbor     !< neighbor
    Integer :: size_vertperneig         !> size of  index_vertperneig
    Integer :: size_neig_edge, neig_col
    Integer :: nloc, col , row, iter_col, iter_row 
    Integer :: mpistatus(MPI_STATUS_SIZE) 
    FTS_INT :: sbuffsize !< size of sbuff(:)
    FTS_INT :: rbuffsize !< size of rbuff(:)
    FTS_INT :: start 
    FTS_INT :: end
 
    ! Arrays
    Integer   , Pointer :: mpisreq(:)               !< MPI requests senders
    XFTS_FLOAT, Pointer :: sbuff (:)             !< MPI send Buffer
    Integer, Pointer    :: index_secin_neig          (:)  !< field of boundary
    Integer, pointer :: ptr_secin_neig_edge (:)
    Integer, pointer :: index_secin_neig_edge (:)

    !- End of header -----------------------------------------------------------

    !---------------------------------------------------------------------------
    ! [1] Init
    !---------------------------------------------------------------------------
    
    ! status
    info = 0
    rank   = ftsl%ikeep(IKEEP_MPIRANK)
    ! mpi
    tagg             =  45
    comm             =  ftsl%comm
    ! interface description
    nb_secin_neig           =  ftsl%lc_domain%nb_secin_neig
    nloc                    =  ftsl%lc_domain%size_domain
    index_secin_neig        => ftsl%lc_domain%index_secin_neig
    ptr_secin_neig_edge     => ftsl%lc_domain%ptr_secin_neig_edge
    index_secin_neig_edge   => ftsl%lc_domain%index_secin_neig_edge
    
    ! values

    !---------------------------------------------------------------------------
    ! [1] Allocate the buffers
    !---------------------------------------------------------------------------
    
    Nullify( sbuff)
    Nullify( mpisreq )

    Allocate(mpisreq(nb_secin_neig), STAT= info )
    CHCKASSRT( info == 0, info )
    If (info < 0) Goto 9999

    !---------------------------------------------------------------------------
    ! [2] Send x and rhs
    !---------------------------------------------------------------------------

    failedrank = FailedNeighbs(1)

    ! search the local index of
    ! the failed neighbor 
       Do neigh=1, nb_secin_neig
          If (failedrank .Eq. index_secin_neig(neigh)) Then
             indexFailed_Neighbor = neigh
             Exit  
          End If
       End Do

       sbuffsize = ptr_secin_neig_edge(indexFailed_Neighbor+1)-&
            ptr_secin_neig_edge(indexFailed_Neighbor)
       Allocate( sbuff(sbuffsize), STAT= info )
       CHCKASSRT( info == 0, info )
       If (info < 0) Goto 9999

       Do i=1,sbuffsize 
          sbuff(i) = -999999999
       End Do
       ! Pack 
       start = ptr_secin_neig_edge(indexFailed_Neighbor)
       end= ptr_secin_neig_edge(indexFailed_Neighbor+1) - 1

       i=1
       Do k = start, end 
          !x
          sbuff(i) = work( index_secin_neig_edge(k) ) 
          i=i+1
       End Do
       
       ! Send
       Call MPI_Isend (sbuff(1), sbuffsize,&
            XFTS_FLOATMPI, failedrank, tagg,comm,mpisreq(1),iinfo)
       ASSRT( iinfo  == MPI_SUCCESS )
       
       Call MPI_WaitAll(1, mpisreq, MPI_STATUSES_IGNORE, iinfo )
       ASSRT( iinfo  == MPI_SUCCESS )

    ! --------------------------------------------------------------------------
    !  Exit routine
    ! --------------------------------------------------------------------------

9999 Continue
    If (Associated(sbuff)) Deallocate(sbuff)
    If (Associated(mpisreq)) Deallocate(mpisreq)

End SUBROUTINE XFTS_FSend_Dist2




SUBROUTINE XFTS_FRecv_LI(& ! intents
     ftsl,                 & ! in
     work,                 & ! inout
     info                  & ! out
     )
  
  !* Module(s) *!
  Use XFTS_ftsolver_type
  Use fts_error_mod        
  Use FTS_ftsolver_enum
  Use FTS_time_mod
  Implicit None
  Include 'mpif.h'
  
  !* Arguments *!
  Type(XFTS_ftsolver_t)   , Intent(inout) :: ftsl
  XFTS_FLOAT            , Intent(inout) :: work(*)
  Integer               , Intent(out )  :: info

    ! Scalars
    Integer :: tagg      !< MPI tagg
    Integer :: comm, rank      !< MPI Communicator, rank
    Integer :: iinfo     !< status of the called routines
    Integer :: i,k     !< dummy iterators
    Integer :: nb_neig      !< field of lc_intrf. see lc_intrf%mynbvi
    Integer :: neigh     !< neighbor
    Integer :: size_vertperneig         !> size of  index_vertperneig
    Integer :: size_neig_edge, neig_col
    Integer :: nloc, col , row, iter_col, iter_row 
    Integer :: mpistatus(MPI_STATUS_SIZE) 
    FTS_INT :: rbuffsize !< size of rbuff(:)
    FTS_INT :: start 
    FTS_INT :: end
 
    ! Arrays

    Integer   , Pointer :: mpirreq(:)               !< MPI requests receivers
    XFTS_FLOAT, Pointer :: rbuff (:)             !< MPI receive Buffer
    Integer, Pointer    :: index_neig          (:)  !< field of boundary
    Integer, pointer :: ptr_direct_neig_edge (:)        !> direct pointer on index_neig_edge

    !- End of header -----------------------------------------------------------

    !---------------------------------------------------------------------------
    ! [1] Init
    !---------------------------------------------------------------------------
    
    ! status
    info = 0
    rank   = ftsl%ikeep(IKEEP_MPIRANK)
    ! mpi
    tagg             =  44
    comm             =  ftsl%comm
    
    ! interface description
    nb_neig             =  ftsl%lc_domain%nb_neig
    size_vertperneig    =  ftsl%lc_domain%size_vertperneig
    nloc                =  ftsl%lc_domain%size_domain
    index_neig          => ftsl%lc_domain%index_neig
    ptr_direct_neig_edge => ftsl%lc_domain%ptr_direct_neig_edge

    !---------------------------------------------------------------------------
    ! [1] Allocate the buffers
    !---------------------------------------------------------------------------
    
    Nullify(rbuff )
    Nullify( mpirreq )

    Allocate( mpirreq(nb_neig), STAT= info )
    CHCKASSRT( info == 0, info )
    If (info < 0) Goto 9999

    rbuffsize = ptr_direct_neig_edge(nb_neig+1)- ptr_direct_neig_edge(1)
    Allocate(rbuff(rbuffsize), STAT= info )
    CHCKASSRT( info == 0, info )
    If (info < 0) Goto 9999


    Do i=1,rbuffsize 
       rbuff(i) = XFTS_FLOATZERO
    End Do


    !---------------------------------------------------------------------------
    ! [2] Start receiving 
    !---------------------------------------------------------------------------

    Do i=1,nb_neig
       
       start = ptr_direct_neig_edge(i)
       end   = ptr_direct_neig_edge(i+1) - 1
       Call MPI_Irecv( rbuff(start),(end-start+1), &
            XFTS_FLOATMPI,index_neig(i), tagg,comm,mpirreq(i),iinfo)
       ASSRT( iinfo  == MPI_SUCCESS )

    End Do

    !---------------------------------------------------------------------------
    !  [3] Unpack values from neigbhors
    !---------------------------------------------------------------------------

    ! save contribution in work(2*nloc+1:3*nloc) to avoid overwriting  rhs
    Do i=1,nloc
       work(2*nloc+i) =  work(nloc+i) 
    End Do

    Do i=1,nb_neig

       Call MPI_WaitAny(nb_neig, mpirreq, neigh, MPI_STATUS_IGNORE, iinfo )
       ASSRT( iinfo  == MPI_SUCCESS )
       start = ptr_direct_neig_edge(neigh)
       end   = ptr_direct_neig_edge(neigh+1) - 1
       Do neig_col= start, end
          col = neig_col+ nloc 
          Do k= ftsl%sm_blockrow%csc(col), ftsl%sm_blockrow%csc(col+1)-1
             row = ftsl%sm_blockrow%i(k)
             ! save contribution in work(2*nloc+1:3*nloc) to avoid overwriting  rhs
             work(2*nloc+row) = work(2*nloc+row) - ftsl%sm_blockrow%v(k)*rbuff(neig_col)
          End Do
       End Do
       
    End Do

    ! --------------------------------------------------------------------------
    !  Exit routine
    ! --------------------------------------------------------------------------

9999 Continue
    If (Associated(rbuff)) Deallocate(rbuff)
    If (Associated(mpirreq)) Deallocate(mpirreq)
  End SUBROUTINE  XFTS_FRecv_LI



SUBROUTINE XFTS_FRecv_LSI(& ! intents
     ftsl,                 & ! in
     work,                 & ! inout
     info                  & ! out
     )
  
  !* Module(s) *!
  Use XFTS_ftsolver_type
  Use fts_error_mod        
  Use FTS_ftsolver_enum
  Use FTS_time_mod
  Implicit None
  Include 'mpif.h'
  
  !* Arguments *!
  Type(XFTS_ftsolver_t)   , Intent(inout) :: ftsl
  XFTS_FLOAT            , Intent(inout) :: work(*)
  Integer               , Intent(out )  :: info

    ! Scalars
    Integer :: tagg, tagg2      !< MPI tagg
    Integer :: comm, rank      !< MPI Communicator, rank
    Integer :: iinfo     !< status of the called routines
    Integer :: i, k, size_delta1
    Integer :: nb_neig, nb_sec_neig     
    Integer :: neigh     !< neighbor
    Integer :: size_vertperneig         !> size of  index_vertperneig
    Integer :: size_neig_edge, neig_col
    Integer :: nloc, col , row, iter_neig, iter_row 
    Integer :: mpistatus(MPI_STATUS_SIZE) 
    FTS_INT :: rbuffsize, rbuffsize2 !< size of rbuff(:)
    FTS_INT :: start, shift_start 
    FTS_INT :: end, interf_size
 
    ! Arrays
    Integer   , Pointer :: mpirreq(:)               !< MPI requests receivers
    XFTS_FLOAT, Pointer :: rbuff (:)             !< MPI receive Buffer
    Integer, Pointer    :: index_neig (:), index_sec_neig(:)
    Integer, pointer :: ptr_direct_neig_edge (:)        !> direct pointer on index_neig_edge

    Integer   , Pointer :: mpirreq2(:)             
    XFTS_FLOAT, Pointer :: rbuff2 (:)             
    Integer, Pointer    :: index_sec_neig_edge (:) 
    Integer, pointer :: ptr_sec_neig_edge (:)     

    !- End of header -----------------------------------------------------------

    !---------------------------------------------------------------------------
    ! [1] Init
    !---------------------------------------------------------------------------
    
    ! status
    info = 0
    rank   = ftsl%ikeep(IKEEP_MPIRANK)
    ! mpi
    tagg             =  44
    tagg2            =  45
    comm             =  ftsl%comm
    
    ! interface description
    size_delta1         = ftsl%lc_domain%size_domain_delta1
    nb_neig             =  ftsl%lc_domain%nb_neig
    nloc                =  ftsl%lc_domain%size_domain
    index_neig          => ftsl%lc_domain%index_neig
    ptr_direct_neig_edge => ftsl%lc_domain%ptr_direct_neig_edge

    nb_sec_neig         =  ftsl%lc_domain%nb_sec_neig
    index_sec_neig      => ftsl%lc_domain%index_sec_neig
    ptr_sec_neig_edge   => ftsl%lc_domain%ptr_sec_neig_edge
    !---------------------------------------------------------------------------
    ! [1] Allocate the buffers
    !---------------------------------------------------------------------------
    
    Nullify(rbuff, rbuff2 )
    Nullify( mpirreq, mpirreq2 )


    Allocate( mpirreq(nb_neig), STAT= info )
    CHCKASSRT( info == 0, info )
    If (info < 0) Goto 9999


    rbuffsize = 2*(ptr_direct_neig_edge(nb_neig+1)- ptr_direct_neig_edge(1))
    Allocate(rbuff(rbuffsize), STAT= info )
    CHCKASSRT( info == 0, info )
    If (info < 0) Goto 9999

    Allocate( mpirreq2(nb_sec_neig), STAT= info )
    CHCKASSRT( info == 0, info )
    If (info < 0) Goto 9999


    rbuffsize2 = ptr_sec_neig_edge(nb_sec_neig+1)- ptr_sec_neig_edge(1)
    Allocate(rbuff2(rbuffsize2), STAT= info )
    CHCKASSRT( info == 0, info )
    If (info < 0) Goto 9999    

    Do i=1,rbuffsize 
       rbuff(i) = XFTS_FLOATZERO
    End Do

    Do i=1,rbuffsize2 
       rbuff2(i) = XFTS_FLOATZERO
    End Do

    ! save rhs of LSI of size delta1  in work(2*delta1+1: 3*delta1) 

    ! set entries to zero
    start = 2*size_delta1+1
    end = 3*size_delta1

    Do i= start,end
       work(i) = XFTS_FLOATZERO
    End Do
    
    ! set the nloc values of rhs
    Do i=1,nloc
       work(2*size_delta1+i) =  work(nloc+i) 
    End Do
    !---------------------------------------------------------------------------
    ! [2] Start receiving neig distance one
    !---------------------------------------------------------------------------

    Do i=1,nb_neig
       
       start = ptr_direct_neig_edge(i)
       end   = ptr_direct_neig_edge(i+1) - 1
       shift_start = 2*(start-1) + 1
       Call MPI_Irecv( rbuff(shift_start), 2*(end-start+1), &
            XFTS_FLOATMPI,index_neig(i), tagg,comm,mpirreq(i),iinfo)
       ASSRT( iinfo  == MPI_SUCCESS )

    End Do

    !---------------------------------------------------------------------------
    ! [3] Start receiving neig distance two
    !---------------------------------------------------------------------------

    Do i=1,nb_sec_neig
       
       start = ptr_sec_neig_edge(i)
       end   = ptr_sec_neig_edge(i+1) - 1
       Call MPI_Irecv( rbuff2(start), (end-start+1), &
            XFTS_FLOATMPI,index_sec_neig(i), tagg2,comm,mpirreq2(i),iinfo)
       ASSRT( iinfo  == MPI_SUCCESS )

    End Do

    !---------------------------------------------------------------------------
    !  [4] Unpack values from neigbhors neig distance one
    !---------------------------------------------------------------------------

    Do iter_neig=1,nb_neig

       Call MPI_WaitAny(nb_neig, mpirreq, neigh, MPI_STATUS_IGNORE, iinfo )
       ASSRT( iinfo  == MPI_SUCCESS )
       start = ptr_direct_neig_edge(neigh)
       end   = ptr_direct_neig_edge(neigh+1) - 1
       interf_size = end-start+1
       shift_start = 2*(start-1) + 1
       i=1
       Do neig_col= start, end
          col = neig_col+ nloc 
          !save rhs
          work(2*size_delta1+col) = work(2*size_delta1+col) + rbuff(shift_start -1 + interf_size + i)
          Do k= ftsl%sm_blockrow_delta1%csc(col), ftsl%sm_blockrow_delta1%csc(col+1)-1
             row = ftsl%sm_blockrow_delta1%i(k)
             ! save contribution in work(2*delta1+1:3*delta1) 
             work(2*size_delta1+row) = work(2*size_delta1+row) - ftsl%sm_blockrow_delta1%v(k)*&
                  rbuff(shift_start-1 +i)
          End Do
          i=i+1
       End Do
    End Do

    !---------------------------------------------------------------------------
    !  [5] Unpack values from neigbhors neig distance two
    !---------------------------------------------------------------------------


    Do iter_neig=1,nb_sec_neig
       Call MPI_WaitAny(nb_sec_neig, mpirreq2, neigh, MPI_STATUS_IGNORE, iinfo )
       
       ASSRT( iinfo  == MPI_SUCCESS )
       start = ptr_sec_neig_edge(neigh)
       end   = ptr_sec_neig_edge(neigh+1) - 1
       Do neig_col= start, end
          col = neig_col+ ftsl%lc_domain%size_domain_delta1 
          Do k= ftsl%sm_blockrow_delta1%csc(col), ftsl%sm_blockrow_delta1%csc(col+1)-1
             row = ftsl%sm_blockrow_delta1%i(k)
             ! save contribution in work(2*delta1+1:3*delta1) 
             work(2*size_delta1+row) = work(2*size_delta1+row) - ftsl%sm_blockrow_delta1%v(k)*&
                  rbuff2(neig_col)
          End Do
       End Do
    End Do
!    write(*,*) "rank", "rhs", work(2*size_delta1+1:3*size_delta1) 
    ! --------------------------------------------------------------------------
    !  Exit routine
    ! --------------------------------------------------------------------------

9999 Continue
    If (Associated(rbuff)) Deallocate(rbuff)
    If (Associated(mpirreq)) Deallocate(mpirreq)
    If (Associated(rbuff2)) Deallocate(rbuff2)
    If (Associated(mpirreq2)) Deallocate(mpirreq2)
  End SUBROUTINE  XFTS_FRecv_LSI




  SUBROUTINE XFTS_InterpLI(& ! intents
       ftsl,                 & ! in
       work,                 & ! inout
       info                  & ! out
       )
  
    
    !* Module *!
    Use XFTS_dense_matrix_mod
    Use FTS_mem_mod
    Use FTS_ftsolver_enum
    Use XFTS_ftsolver_type
    Use XFTS_mumps_mod
    Implicit None
    Include 'mpif.h'
  
    !* Arguments *!
    Type(XFTS_ftsolver_t)   , Intent(inout) :: ftsl
    XFTS_FLOAT            , Intent(inout) :: work(*)
    Integer               , Intent(out )  :: info
    

    
    !Local variables 

    Integer :: nloc, i, unit
    Type(XFTS_ARITHmumps_struc)     :: id_mumps
    Type(XFTS_dense_matrix_t)  :: rhs

     
    nloc = ftsl%lc_domain%size_domain
    
    Call XFTS_DM_Create (rhs, nloc, 1,nloc,info)

     ! Copy work(nloc+1:2*nloc) to into rhs
    Call XFTS_ARITHcopy(nloc, work(2*nloc+1), 1,rhs%v,1)
    

    Call XFTS_mumps_Set_MPICommunicator(id_mumps,MPI_COMM_SELF,info)
    CHCKASSRT( info >= 0, info )
    If (info < 0 ) Return
  
    Call XFTS_mumps_Set_matrix (id_mumps,ftsl%sm_li,info)
    CHCKASSRT( info >= 0, info )
    If (info < 0 ) Return
    
  
    Call XFTS_mumps_Analyze   (id_mumps,info)
    CHCKASSRT( info >= 0, info )
    If (info < 0 ) Return
  
    Call XFTS_mumps_Factorize (id_mumps,info)
    CHCKASSRT( info >= 0, info )
    If (info < 0 ) Return

    Call XFTS_mumps_solve_RHS (&
         id_mumps,&
         rhs%v, rhs%n,rhs%ld,&
         info)
    CHCKRINFO(info)
    If (info < 0) Goto 9999

    Call XFTS_ARITHcopy(nloc, rhs%v, 1, work(1),1)

    ftsl%rinfo(RINFO_FLOP_LI) = id_mumps%rinfog(3)
9999 Continue

End SUBROUTINE XFTS_InterpLI



SUBROUTINE XFTS_InterpLSI(& ! intents
       ftsl,                 & ! in
       work,                 & ! inout
       info                  & ! out
       )
  
    !* Module *!
    Use XFTS_sparse_matrix_mod
    Use FTS_ftsolver_enum
    Use XFTS_ftsolver_type
    Use XFTS_ARITHqrm_mod
    Implicit None
    Include 'mpif.h'

    !* Arguments *!
    Type(XFTS_ftsolver_t)   , Intent(inout) :: ftsl
    XFTS_FLOAT            , Intent(inout) :: work(*)
    Integer               , Intent(out )  :: info
    !Local variables 

    Integer :: nnz, m, n, i, unit
    type(XFTS_ARITHqrm_spmat_type)          :: qrm_mat
    XFTS_FLOAT, allocatable  :: b(:), x(:), r(:)
    Real(kind=XFTS_FLOATKIND)      :: rnrm, onrm
    character                       :: matfile*30='', transp 

    ! nnz = ftsl%sm_lsi%nnz
    ! m   = ftsl%sm_lsi%m
    ! n   = ftsl%sm_lsi%n

    nnz = ftsl%sm_li%nnz
    m   = ftsl%sm_li%m
    n   = ftsl%sm_li%n

    ftsl%rinfo(RINFO_M_QR)                = m
    ftsl%rinfo(RINFO_N_QR)                = n
    ! initialize the control data structure. 
    call qrm_spmat_init(qrm_mat)
    call qrm_set('qrm_eunit', 6)
    
    ! allocate arrays for the input matrix

    call qrm_palloc(qrm_mat%irn, nnz)
    call qrm_palloc(qrm_mat%jcn, nnz)
    call qrm_palloc(qrm_mat%val, nnz)

    
    !     ! initialize the input matrix
    !Do i=1,nnz
    ! qrm_mat%irn => ftsl%sm_lsi%i
    ! qrm_mat%jcn => ftsl%sm_lsi%j
    ! qrm_mat%val => ftsl%sm_lsi%v

    qrm_mat%irn => ftsl%sm_li%i
    qrm_mat%jcn => ftsl%sm_li%j
    qrm_mat%val => ftsl%sm_li%v
    !End Do
    
    qrm_mat%m   = m
    qrm_mat%n   = n
    qrm_mat%nz  = nnz

    !set some control parameters
    call qrm_set(qrm_mat, 'qrm_ib',120)
    call qrm_set(qrm_mat, 'qrm_nb',120)
    call qrm_set(qrm_mat, 'qrm_nthreads',1)
    call qrm_set(qrm_mat, 'qrm_ordering',qrm_metis_)

    call qrm_analyse(qrm_mat, transp)
    call qrm_factorize(qrm_mat, transp)

    call qrm_aalloc(b, qrm_mat%m)
    call qrm_aalloc(r, qrm_mat%m)
    call qrm_aalloc(x, qrm_mat%n)
    
    
    !Copy the updated rhs into b
    Do i=1,m
       b(i) = work(2*m+i)
    End Do
    
    
    call qrm_apply(qrm_mat, 't', b)
    call qrm_solve(qrm_mat, 'n', b, x)

    ! Copy solution to work(1:nloc)
    Do i=1,n
       work(i) = x(i)
    End Do

    ftsl%rinfo(RINFO_FLOP_LSI) = qrm_mat%gstats(qrm_facto_flops_)
  ! compute the residual
  !call qrm_residual_norm(qrm_mat, r, x, rnrm)
  !call qrm_residual_orth(qrm_mat, r, onrm)   
  !write(*,'("||r||/||A||    = ",e10.2)')rnrm
  !write(*,'("||A^tr||/||r|| = ",e10.2)')onrm

9999 Continue
    open (UNIT=unit,file="qr.mtx",action="write",status="replace")
    Call  XFTS_sm_mmwrite(ftsl%sm_lsi,unit,info)
    close(unit)
    call qrm_adealloc(b)
    call qrm_adealloc(r)
    call qrm_adealloc(x)
    !call qrm_spmat_destroy(qrm_mat, all=.true.)
  
End SUBROUTINE XFTS_InterpLSI




END MODULE XFTS_fault_mod


