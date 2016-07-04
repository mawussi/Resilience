#include "fts_defs_f.h"
#include "fts_macros_f.h"

! [+] module : fts_env_mod ---------------------------------------------------
!
!> Ftsolver Environement module
!! 
!! Provides environment detection/manipulation for ftsolver.
!!
Module fts_env_mod
  
  !* No implicit typing *!
  
  Implicit None

  !* Defined type(s) *!
  
  Include "fts_env_type_f.inc"

  !* Private constants *!
  Character(len=FTSOLVER_STRL), Private, Parameter :: FLNAME = &
       "fts_env_mod.F90"

  !* Access specifiers *!

  Public  :: fts_env_init
  Public  :: fts_env_exit
  Public  :: fts_env_setMPIComm
  Public  :: fts_env_getNodeId
  Public  :: fts_env_getNodeComm
  Public  :: fts_env_getNbNodes

  Private :: comm_splitbynode

  !* Routines *!

  Contains
    !
    ! constructors
    ! 

    Subroutine fts_env_init(menv)
      Implicit None
      Type(fts_env_t), Intent(inout) :: menv
      menv%comm     = -1
      menv%nodecomm = -1
      menv%nodeid   = -1
    End Subroutine fts_env_init

    Subroutine fts_env_exit(menv)
      Use fts_error_mod
      Implicit None
      Include "mpif.h"
      Type(fts_env_t), Intent(inout) :: menv

      Logical :: MPICommIsSet 
      Integer :: info
      
      MPICommIsSet = ( menv%nodeid /= -1 )
      If ( MPICommIsSet )Then
         Call MPI_Comm_free(menv%comm,info)
         ASSRT( info == MPI_SUCCESS )

         Call MPI_Comm_free(menv%nodecomm,info)
         ASSRT( info == MPI_SUCCESS )
      End If

      menv%comm     = -1
      menv%nodecomm = -1
      menv%nodeid   = -1
    End Subroutine fts_env_exit
    
    !
    ! setters
    !

    Subroutine fts_env_setMPIComm(menv, comm)

      !* Module(s) *!

      Use fts_error_mod
      Implicit None
      Include "mpif.h"
      
      !* Argument(s) *!

      Type(fts_env_t), Intent(inout) :: menv
      Integer           , Intent(in   ) :: comm

      !* Local variable(s) *!

      Integer :: rank 
      Integer :: info
      Integer :: np

      !- End of header -----------------------------------------------------------

      menv%comm = comm

      Call MPI_Comm_size(menv%comm,np  , info)
      ASSRT( info == MPI_SUCCESS )
      Call MPI_Comm_rank(menv%comm,rank, info)
      ASSRT( info == MPI_SUCCESS )
      Call comm_splitbynode(menv%comm,np,menv%nodecomm,menv%nodeid, info)
      ASSRT( info == MPI_SUCCESS )

#if FTS_ENV_DEBUG
      Call MPI_Comm_size( menv%nodecomm, np, info )
      ASSRT( info == MPI_SUCCESS )

      print *,"NodeId", menv%nodeid, &
           "have ", np, "elements ", &
           "(from Process rank ",rank,")"
#endif


    End Subroutine fts_env_setMPIComm


    !
    ! getters
    !

    Integer Function fts_env_getNodeId(menv)
      Type(fts_env_t), Intent(in) :: menv
      fts_env_getnodeid=menv%nodeid
    End Function fts_env_getnodeid

    Integer Function fts_env_getNodeComm(menv)
      Type(fts_env_t), Intent(in) :: menv
      fts_env_getnodecomm=menv%nodecomm
    End Function fts_env_getnodecomm

    Integer Function fts_env_getNbNodes(menv)
      Use fts_error_mod
      Implicit None
      Include "mpif.h"
      Type(fts_env_t), Intent(in) :: menv
      Integer :: info, res
      Call MPI_AllReduce(&
           menv%nodeid,res,1,MPI_INTEGER,&
           MPI_MAX,menv%comm,info)
      ASSRT( info == MPI_SUCCESS )
      fts_env_getnbnodes = res
    End Function fts_env_getnbnodes

    
    !
    ! Methods
    !
    Subroutine comm_splitbynode(comm, commnp, subcomm, subcommid, info )

      !* Module(s) *!

      Use fts_error_mod
      Implicit None
      Include "mpif.h"
      
      !* Argument(s) *!

      Integer, Intent(in ) :: comm
      Integer, Intent(in ) :: commnp
      Integer, Intent(out) :: subcomm
      Integer, Intent(out) :: subcommid
      Integer, Intent(out) :: info
      
      !* Local variable(s) *!

      ! Scalars
      Integer :: i, st, ed
      Integer :: rank
      Integer :: hkey

      ! Strings
      Character*(MPI_MAX_PROCESSOR_NAME)     :: hname,hname2
      Character*(commnp*MPI_MAX_PROCESSOR_NAME) :: hnames

      !- End of header -----------------------------------------------------------

      hname=""
      hnames=""
      subcomm   = -1
      subcommid = -1
      info      =  MPI_SUCCESS

      Call MPI_Comm_rank(comm,rank, info)
      ASSRT( info == MPI_SUCCESS )

      Call MPI_Get_Processor_name(hname, i ,info)
      ASSRT( info == MPI_SUCCESS )

      Call MPI_AllGather( &
           hname  ,MPI_MAX_PROCESSOR_NAME  ,MPI_CHARACTER, &
           hnames ,MPI_MAX_PROCESSOR_NAME  ,MPI_CHARACTER, &
           comm,info)
      ASSRT( info == MPI_SUCCESS )

      ! get the color
      hname2=hnames(1:MPI_MAX_PROCESSOR_NAME)
      subcommid = 1
      Do i=1,commnp
         st=(i-1)*MPI_MAX_PROCESSOR_NAME + 1
         ed=i *MPI_MAX_PROCESSOR_NAME
         If ( hname2 /= hnames(st:ed)) subcommid = subcommid+1
         If ( hname  == hnames(st:ed)) exit
      End Do

      ! get the key
      hkey=0
      Do i=1,rank+1
         st=(i-1)*MPI_MAX_PROCESSOR_NAME + 1
         ed=i *MPI_MAX_PROCESSOR_NAME
         If ( hnames(st:ed) == hname ) hkey=hkey+1
      End Do

      Call MPI_COMM_SPLIT( comm, subcommid, hkey, subcomm, info)
      ASSRT( info == MPI_SUCCESS )

    End Subroutine comm_splitbynode

End Module fts_env_mod
