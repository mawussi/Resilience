! Warning: XFTS_GENFILE_COMMENT
#include "fts_defs_f.h"
#include "fts_macros_f.h"  
  !> Update the state of ftsolver instance
  Module XFTS_state_update_mod

    ! List of routines 

    Public :: XFTS_state_updateStatus

    Private :: XFTS_state_rsummarize

    ! Private constants
    Character(len=FTSOLVER_STRL), Parameter, Private :: &
         FLNAME= "XFTS_ARITHfts_ftsolver_read_mod.F90"
    
    Contains


    ! [+] routine : XFTS_state_UpdateStatus -----------------------------------------
    !
    !> update the status of the ftsolver instance.
    !!
    !! Namely update the informations (iinfo*/rinfo*) of the instance.
    !!  
    !!----
    !!
    !! @param [in,out] ftsl
    !!        the ftsolver instance, it uses the fields:
    !!        - comm        [in ]
    !!        - rinfo*      [in,out]
    !!
    !!----
    !!

    !!
    Subroutine XFTS_state_UpdateStatus(ftsl)

      !* Modules & co. *!

      Use FTS_mem_mod
      Use fts_log_mod
      Use fts_error_mod
      Use FTS_ftsolver_enum
      Use XFTS_ftsolver_type

      Implicit None

      !* Subroutine arguments *!
      Type(XFTS_ftsolver_t), Intent(inout) :: ftsl

      !* Local Variables *!

      ! Scalars
      Integer    :: iinfo
      Integer    :: rank, master
      Integer    :: sverb,sunit

      !-------------------------------------------------------------------------
      ! [1] Init
      !-------------------------------------------------------------------------

      rank   = ftsl%ikeep(IKEEP_MPIRANK)
      master = ftsl%ikeep(IKEEP_HOSTRANK)

      Call XFTS_state_rsummarize(&
           ftsl%comm, FTSOLVER_RINFO_SIZE, ftsl%rinfo,&
           ftsl%rinfomin,ftsl%rinfomax,&
           ftsl%rinfoavg,ftsl%rinfosig,&
           iinfo)
      CHCKASSRT(iinfo >= 0, iinfo)
      If ( iinfo < 0 ) Goto 9999

      !-------------------------------------------------------------------------
      ! Finish
      !-------------------------------------------------------------------------

9999  Continue
      ! set return code
    End Subroutine XFTS_state_UpdateStatus



  ! [+] routine : XFTS_state_rsummarize ----------------------------------------
  !
  !> Get the summarised statistics from an float statistic.
  !!
  !! @param[in ] comm  The MPI Communicator
  !! @param[in ] nstat The number of statistics
  !! @param[in ] rstat  The list of statistics
  !! @param[out] rmin   The minimum value of the statistic 
  !! @param[out] rmax   The maximum value of the statistic 
  !! @param[out] ravg   The average value of the statistic 
  !! @param[out] rsig   The standard deviation of the statistic 
  !! @param[out] info  the routine status
  !!
  !! @author Yohan Lee-tin-yien
  !!
  Subroutine XFTS_state_rsummarize&
       ( comm, nstat, rstat, rmin, rmax, ravg, rsig, info )
    
    !* Module    *!

    Use fts_error_mod
    Implicit None
    Include "mpif.h"

    !* Arguments *!

    Integer, Intent(in ) ::  comm  
    Integer, Intent(in ) ::  nstat 
    Real(kind=8), Intent(in ) ::  rstat (nstat)
    Real(kind=8), Intent(out) ::  rmin  (nstat) 
    Real(kind=8), Intent(out) ::  rmax  (nstat) 
    Real(kind=8), Intent(out) ::  ravg  (nstat) 
    Real(kind=8), Intent(out) ::  rsig  (nstat)
    Integer, Intent(out) ::  info  

    !* Local variables *!

    ! scalars
    Integer              :: np
    Integer              :: p ! counter of processes
    Integer              :: s ! counter on statistics
    Real(kind=8) :: rnp
    Real(kind=8) :: rval

    ! arrays
    Real(kind=8), Allocatable :: statlist(:)

    !- End of header -----------------------------------------------------------
    !---------------------------------------------------------------------------
    ! [1] Init 
    !---------------------------------------------------------------------------

    Call MPI_Comm_size( comm, np, info )
    ASSRT( info == MPI_SUCCESS )

    Allocate(statlist(nstat*np),STAT=info )
    CHCKASSRT( info == 0, info )
    If (info < 0) Return
    
    rnp = REAL(np,KIND=8)

    !---------------------------------------------------------------------------
    ! [2] Gather the values
    !---------------------------------------------------------------------------

    Call MPI_AllGather( &
         rstat(1)   , nstat, MPI_DOUBLE_PRECISION, &
         statlist(1), nstat, MPI_DOUBLE_PRECISION, &
         comm, info )
    ASSRT( info == MPI_SUCCESS )



    !---------------------------------------------------------------------------
    ! [3] Compute min/max/avg 
    !---------------------------------------------------------------------------
    
    Do s=1, nstat

       rmin(s) = statlist(s)
       rmax(s) = statlist(s)
       ravg(s) = statlist(s)
       Do p=2,np
          rval  = statlist(s+(p-1)*nstat)
          If (rmin(s) > rval) rmin(s) = rval
          If (rmax(s) < rval) rmax(s) = rval
          ravg(s) = ravg(s) + rval
       End Do
       ravg(s) = ravg(s)/rnp 
       
    End Do

    !---------------------------------------------------------------------------
    ! [4] Compute sig the standard deviation
    !---------------------------------------------------------------------------
    
    Do s=1, nstat

       rsig(s) = 0.d0
       Do p=1,np
          rval     = statlist(s+(p-1)*nstat)
          rsig(s)  = rsig(s) + ( rval - ravg(s) )**2
       End Do
       rsig(s) = SQRT( rsig(s)/rnp )

    End Do

    !---------------------------------------------------------------------------
    ! [5] Finish
    !---------------------------------------------------------------------------

    Deallocate( statlist )

  End Subroutine XFTS_state_rsummarize


 
  End Module XFTS_state_update_mod
