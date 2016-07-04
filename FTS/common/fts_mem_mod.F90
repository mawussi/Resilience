! ftsolver memory module.
#include "fts_defs_f.h"
#include "fts_macros_f.h"

! [+] module : FTS_mem_mod ---------------------------------------------------
!
!> Ftsolver memory module
!! 
Module FTS_mem_mod
  
  !* No implicit typing *!
  Implicit None

  !* Private constants *!
  Character(len=FTSOLVER_STRL), Private, Parameter ::&
       FLNAME = "FTS_mem_mod.F90"

  !* Type definition *! 
  Type FTS_mem_t; sequence

     ! memory usage related to the instance (in Bytes)
     Integer(Kind=8) :: memusage 

     ! memory peak  related to the instance (in Bytes)
     Integer(Kind=8) :: mempeak

     ! memory usage of a package related to the instance (in Bytes).
     ! It is process dependent. 
     Integer(Kind=8) :: pidmemusage 

     ! memory peak of a package related to the instance (in Bytes).
     ! It is process dependent. 
     Integer(Kind=8) :: pidmempeak

  End type FTS_mem_t

  !* Access specifiers *!
  Public :: byte2Mbyte

  Public :: FTS_mem_init
  Public :: FTS_mem_exit
  Public :: FTS_mem_build

  Public :: FTS_mem_setusage
  Public :: FTS_mem_setpeak
  Public :: FTS_mem_setpidusage
  Public :: FTS_mem_setpidpeak

  Public :: FTS_mem_getpeak
  Public :: FTS_mem_getusage
  Public :: FTS_mem_getpidpeak
  Public :: FTS_mem_getpidusage
  Public :: FTS_mem_getnodepeak
  Public :: FTS_mem_getnodeusage

  Public :: FTS_mem_update
  Public :: FTS_mem_add2mem
  Public :: FTS_mem_add2usage
  Public :: FTS_mem_add2peak
  Public :: FTS_mem_add2pidpeak
  Public :: FTS_mem_sub2usage
  
  !* routines *!

  Contains


    ! [+] function : Byte2MByte ------------------------------------------------
    ! 
    !> convert a size in bytes into a size in Mega bytes
    Integer(kind=4) Function Byte2MByte(bytes)
      Implicit None
      Integer(kind=8) :: bytes
      Real   (kind=8) :: conv
      
      conv   = 1_8/Real(1048576,KIND=8)
      Byte2MByte = Ceiling( Real(bytes,KIND=8) * conv )
    End Function Byte2MByte
    
    !
    ! constructors & destructors
    !

  Subroutine FTS_mem_init( mem )
    Type(FTS_mem_t), Intent(out) :: mem
    mem%memusage     = 0
    mem%mempeak      = 0
    mem%pidmempeak   = 0
    mem%pidmemusage  = 0
  End Subroutine FTS_mem_init

  Subroutine FTS_mem_exit( mem )
    Type(FTS_mem_t), Intent(out) :: mem
    mem%memusage     = 0
    mem%mempeak      = 0
    mem%pidmempeak   = 0
    mem%pidmemusage  = 0
  End Subroutine FTS_mem_exit

  Subroutine FTS_mem_build( mem, memusage, mempeak, pidmemusage, pidmempeak )
    Type(FTS_mem_t), Intent(out) :: mem
    Integer(Kind=8) , Intent(in) :: memusage 
    Integer(Kind=8) , Intent(in) :: mempeak
    Integer(Kind=8) , Intent(in) :: pidmemusage 
    Integer(Kind=8) , Intent(in) :: pidmempeak
    mem%memusage     = memusage
    mem%mempeak      = mempeak
    mem%pidmempeak   = pidmemusage
    mem%pidmemusage  = pidmempeak
  End Subroutine FTS_mem_build

  !
  ! setters
  !

  Subroutine FTS_mem_setusage( mem, memset )
    Type(FTS_mem_t), Intent(inout) :: mem
    Integer(Kind=8) , Intent(in) :: memset
    mem%memusage  = memset
    Call FTS_mem_update(mem)
  End Subroutine FTS_mem_setusage

  Subroutine FTS_mem_setpeak( mem, memset )
    Type(FTS_mem_t), Intent(inout) :: mem
    Integer(Kind=8) , Intent(in) :: memset
    mem%mempeak    = memset
    Call FTS_mem_update(mem)
  End Subroutine FTS_mem_setpeak

  Subroutine FTS_mem_setpidusage( mem, memset )
    Type(FTS_mem_t), Intent(inout) :: mem
    Integer(Kind=8) , Intent(in) :: memset
    mem%pidmemusage  = memset
    Call FTS_mem_update(mem)
  End Subroutine FTS_mem_setpidusage

  Subroutine FTS_mem_setpidpeak( mem, memset )
    Type(FTS_mem_t), Intent(inout) :: mem
    Integer(Kind=8) , Intent(in) :: memset
    mem%pidmempeak    = memset
    Call FTS_mem_update(mem)
  End Subroutine FTS_mem_setpidpeak

  !
  ! getters
  !

  Integer(Kind=8) Function FTS_mem_getpidusage(mem)
    Type(FTS_mem_t), Intent(in) :: mem
    FTS_mem_getpidusage = mem%pidmemusage
  End Function  FTS_mem_getpidusage

  Integer(Kind=8) Function FTS_mem_getpidpeak(mem)
    Type(FTS_mem_t), Intent(in) :: mem
    FTS_mem_getpidpeak = mem%pidmempeak
  End Function  FTS_mem_getpidpeak

  Integer(Kind=8) Function FTS_mem_getusage(mem)
    Type(FTS_mem_t), Intent(in) :: mem
    FTS_mem_getusage = mem%memusage
  End Function  FTS_mem_getusage

  Integer(Kind=8) Function FTS_mem_getpeak(mem)
    Type(FTS_mem_t), Intent(in) :: mem
    FTS_mem_getpeak = mem%mempeak
  End Function  FTS_mem_getpeak

  Integer(Kind=8) Function FTS_mem_getallpeak(mem)
    Type(FTS_mem_t), Intent(in) :: mem
    FTS_mem_getallpeak = mem%mempeak + mem%pidmempeak
  End Function  FTS_mem_getallpeak

  Integer(Kind=8) Function FTS_mem_getallusage(mem)
    Type(FTS_mem_t), Intent(in) :: mem
    FTS_mem_getallusage = mem%memusage + mem%pidmemusage
  End Function  FTS_mem_getallusage


  Integer(Kind=8) Function FTS_mem_getnodeusage(mem, env)

    Use fts_error_mod
    Use fts_env_mod
    Implicit None
    Include "mpif.h"

    ! Return the memory usage of the node

    Type(FTS_mem_t), Intent(in) :: mem
    Type(fts_env_t), Intent(in) :: env
    
    Integer :: nodcomm, ierr
    Integer(Kind=8) :: pidmem, nodmem

    nodcomm = fts_env_getNodeComm(env)
    pidmem=mem%memusage + mem%pidmemusage
    Call MPI_AllReduce(pidmem,nodmem,1,&
         MPI_INTEGER8,MPI_SUM,&
         nodcomm, ierr)
    FTS_mem_getnodeusage = nodmem

  End Function  FTS_mem_getnodeusage

  Integer(Kind=8) Function FTS_mem_getnodepeak(mem, env)

    Use fts_env_mod
    Use fts_error_mod
    Implicit None
    Include "mpif.h"

    Type(FTS_mem_t), Intent(in) :: mem
    Type(fts_env_t), Intent(in) :: env
    
    Integer :: nodcomm, ierr
    Integer(Kind=8) :: pidmem, nodmem

    ! Return the memory peak of the node

    nodcomm = fts_env_getNodeComm(env)
    pidmem=mem%mempeak + mem%pidmempeak
    Call MPI_AllReduce(pidmem,nodmem,1,&
         MPI_INTEGER8,MPI_SUM,&
         nodcomm, ierr)
    FTS_mem_getnodepeak = nodmem
  End Function  FTS_mem_getnodepeak

  !
  ! methods
  !

  Subroutine FTS_mem_update(mem)
    Type(FTS_mem_t), Intent(inout) :: mem
    mem%pidmempeak = Max(mem%pidmempeak, mem%pidmemusage)
    mem%mempeak = Max(mem%mempeak, mem%memusage)
  End Subroutine FTS_mem_update

  Subroutine FTS_mem_add2mem( mem, memadded )
    Type(FTS_mem_t), Intent(inout) :: mem
    Type(FTS_mem_t), Intent(in   ) :: memadded

    mem%mempeak     = Max &
         (mem%memusage + memadded%mempeak,mem%mempeak)
    mem%pidmempeak  = Max &
         (mem%pidmemusage + memadded%pidmempeak,mem%pidmempeak)

    mem%memusage    = mem%memusage      + memadded%memusage
    mem%pidmemusage = mem%pidmemusage   + memadded%pidmemusage

    Call FTS_mem_update(mem)
  End Subroutine FTS_mem_add2mem

  Subroutine FTS_mem_add2usage( mem, memadded )
    Type(FTS_mem_t), Intent(inout) :: mem
    Integer(Kind=8) , Intent(in)   :: memadded 
    mem%memusage     = mem%memusage + memadded
    Call FTS_mem_update(mem)
  End Subroutine FTS_mem_add2usage

  Subroutine FTS_mem_add2peak( mem, memadded )
    Type(FTS_mem_t), Intent(inout) :: mem
    Integer(Kind=8) , Intent(in) :: memadded
    Call FTS_mem_update(mem)
    mem%mempeak     = mem%mempeak + memadded
  End Subroutine FTS_mem_add2peak

  Subroutine FTS_mem_add2pidpeak( mem, memadded )
    Type(FTS_mem_t), Intent(inout) :: mem
    Integer(Kind=8) , Intent(in) :: memadded
    Call FTS_mem_update(mem)
    mem%pidmempeak     = mem%pidmempeak + memadded
  End Subroutine FTS_mem_add2pidpeak

  Subroutine FTS_mem_sub2usage( mem, memsubstracted )
    Type(FTS_mem_t), Intent(inout) :: mem
    Integer(Kind=8) , Intent(in)   :: memsubstracted
    mem%memusage     = mem%memusage - memsubstracted
  End Subroutine FTS_mem_sub2usage

  End Module FTS_mem_mod
