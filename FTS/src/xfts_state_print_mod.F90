! Warning: XFTS_GENFILE_COMMENT
#include "fts_defs_f.h"
#include "fts_macros_f.h"

! [+] module : XFTS_state_print_mod --------------------------------------------------
!
!> module to handle statistics in ftsolver
!!
Module XFTS_state_print_mod
  !* Modules *!
  Use FTS_log_mod
  Use FTS_error_mod
  Use FTS_ftsolver_enum
  Use XFTS_ftsolver_type
  
  !* No Implicit Typing *!
  Implicit None

  !* Access specifiers *!

  Public :: XFTS_state_print

  Private :: PrintSH  ! section header
  Private :: PrintTH1 ! table header (1 value)
  Private :: PrintTH4 ! table header (4 values)
  Private :: PrintTF1 ! table footer (1 value)
  Private :: PrintTF4 ! table footer (4 values)
  Private :: PrintTL  ! table line 
  Private :: PrintINTERP  ! table line   

  !* Private constants *!

  Character(len=FTSOLVER_STRL), Private, Parameter ::&
       FLNAME = "XFTS_state_print_mod.F90"


  !* routines *!

Contains

  ! [+] routine : XFTS_state_print ---------------------------------------------
  !
  !> Print the ftsolver statistics
  !!
  !! @param[in ] unit  The unit where to print the statistics
  !! @param[in ] ftsl  The ftsolver instance containing the stats.
  !!          - comm
  !!          - {rinfo}{,min,max,avg}
  !! @param[out] info  the routine status
  !!
  !!
  Subroutine XFTS_state_print( ftsl, info )
    Use FTS_log_mod
    Use XFTS_ftsolver_type
    Use FTS_ftsolver_enum
    Implicit None

    !* Arguments *!

    Type(XFTS_ftsolver_t) , Intent(in ) :: ftsl
    Integer             , Intent(out) :: info

    !* Local variables *!

    ! scalars
    Integer :: verb
    Integer :: unit
    Integer :: lend            ! lenght of the description

    !- End of header -----------------------------------------------------------

    !---------------------------------------------------------------------------
    ! [0] Init
    !---------------------------------------------------------------------------

    info = 0
    If (ftsl%ikeep(IKEEP_MPIRANK) /= ftsl%ikeep(IKEEP_HOSTRANK) ) Return 

#define PRINFO(x,y) Call PrintTL(x,lend,ftsl,unit,y)
#define PRINTERP(x,y) Call PrintInterp(x,lend,ftsl,unit,y)


    !---------------------------------------------------------------------------
    ! [] Print the timings
    !---------------------------------------------------------------------------

    lend=1 + Len_Trim("Timing -> Analyze -> Pretreat Input Mtx.") 

    Call PrintSH(unit,'Timings')

    PRINFO(RINFO_TIMING_Solve,"Timing -> Solve")
    
    Call PrintSH(unit,'Interpolation')

       PRINTERP(RINFO_TIMING_Fault,"Timing ->  Choice of faulty proc")       
       PRINTERP(RINFO_TIMING_Restart,"Timing -> Restart")
       PRINTERP(RINFO_TIMING_Interp,"Timing ->  Intepolation")
       PRINTERP(RINFO_TIMING_InterpRhs,"Timing ->  Intepolation -> Rhs ")
       PRINTERP(RINFO_TIMING_InterpSlv,"Timing ->  Intepolation -> Solve")

    Call PrintSH(unit,'More details')
       PRINTERP(RINFO_FLOP_LI,"Flop -> LI")
       PRINTERP(RINFO_FLOP_LSI,"Flop -> LSI")
       PRINTERP(RINFO_M_QR,"Size -> QR -> LINE")
       PRINTERP(RINFO_N_QR,"Size -> QR -> COL")

    Call PrintTF4(unit,lend)


#undef PRINFO
#undef PRINTERP


    info = 0

  End Subroutine XFTS_state_print

  !> Print a section header
  Subroutine PrintSH(unit,msg)
    Integer, Intent(in) :: unit
    Character(len=*) :: msg
    Integer :: i

    Write(unit,*)
    Write(unit,*) '     * ',Trim(msg)
    Write(unit,*)
    ! Write(unit,*) '* ',Trim(msg)
    ! Write(unit,*) ('=',i=1,Len_Trim(msg)+2)

  End Subroutine PrintSH

  !> Print a table's header 
  !! where the table's line consisting have only 1 value
  Subroutine PrintTH1(unit,lend)
    Integer, Intent(in) :: unit
    Integer, Intent(in) :: lend
    
    Character(len=FTSOLVER_STRL) :: stmp
    Integer :: i
    !-

    stmp = "Description"
    Write(unit,*)
    Write(unit,FMT='(A,A9,A5,A,A11)') &
         stmp(1:lend),"  FIELD (","INDEX",")","VALUE"
    Write(unit,*)  ('-', i=1,lend+15+11)


  End Subroutine PrintTH1

  !> Print a table's header 
  !! where the table's line consisting have
  !! 4 values (min,max,avg,std)
  Subroutine PrintTH4(unit,lend)
    Integer, Intent(in) :: unit
    Integer, Intent(in) :: lend

    Character(len=FTSOLVER_STRL) :: stmp
    Integer :: i
    !-

    Write(unit,*)
    stmp = "Description"
    Write(unit,FMT='(A,A9,A5,A1,4A11)') &
         stmp(1:lend),"  FIELD (","INDEX",")",&
         'MIN','MAX','AVERAGE','STD.DEV.'
    Write(unit,*) ('-', i=1,lend+15+4*11)


  End Subroutine PrintTH4



  !> Print a table's footer
  !! where the table's line consisting have 1 value
  !!  values (min,max,avg,std)
  Subroutine PrintTF1(unit,lend)
    Integer, Intent(in) :: unit
    Integer, Intent(in) :: lend

    Write(unit,*)
    ! Integer :: i
    ! Write(unit,*) ('-', i=1,lend+15+11)


  End Subroutine PrintTF1

  !> Print a table's footer
  !! where the table's line consisting have
  !! 4 values (min,max,avg,std)
  Subroutine PrintTF4(unit,lend)
    Integer, Intent(in) :: unit
    Integer, Intent(in) :: lend

    Write(unit,*)
    
    ! Integer :: i
    ! Write(unit,*) ('-', i=1,lend+15+4*11)

  End Subroutine PrintTF4

  ! Print a table's line 
  Subroutine PrintTL(ind,lend,ftsl,unit,msg)
    Use XFTS_ftsolver_type
    Implicit None

    Integer, Intent(in) :: ind
    Integer, Intent(in) :: lend
    Type(XFTS_ftsolver_t), Intent(in) :: ftsl
    Integer, Intent(in) :: unit
    Character(len=*), Intent(in) :: msg

    Character(len=FTSOLVER_STRL) :: stmp

    !-
    stmp = Trim(msg) 
    Write(unit,'(A,A9,I5,A,4(1PE11.3))') &
         stmp(1:lend),"  RINFO (",ind,")",&
         ftsl%rinfomin(ind), ftsl%rinfomax(ind),&
         ftsl%rinfoavg(ind), ftsl%rinfosig(ind)
    
  End Subroutine PrintTL


  Subroutine PrintInterp(ind,lend,ftsl,unit,msg)
    Use XFTS_ftsolver_type
    Implicit None

    Integer, Intent(in) :: ind
    Integer, Intent(in) :: lend
    Type(XFTS_ftsolver_t), Intent(in) :: ftsl
    Integer, Intent(in) :: unit
    Character(len=*), Intent(in) :: msg

    Character(len=FTSOLVER_STRL) :: stmp


    stmp = Trim(msg) 
    Write(unit,'(A,A9,I5,A,1PE11.3)') &
         stmp(1:lend),"  RINFO (",ind,")",&
         ftsl%rinfomax(ind)

    
  End Subroutine PrintInterp
End Module XFTS_state_print_mod

