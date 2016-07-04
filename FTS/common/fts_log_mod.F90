!> Module to log ftsolver messages
Module fts_log_mod

  !* No Implicit Typing *!
  Implicit None

  !* Private variables *!
  Integer, Private, Save :: unit_err  = 0
  Integer, Private, Save :: unit_warn = 0
  Integer, Private, Save :: unit_std  = 6
  Integer, Private, Save :: verbosity = 3
  Integer, Private, Save :: rank      = -1

  !* Enumerations *!
  Integer, Public, Parameter :: MSG_ERROR    = 1
  Integer, Public, Parameter :: MSG_WARNING  = 2
  Integer, Public, Parameter :: MSG_STD      = 3
  Integer, Public, Parameter :: MSG_VERBOSE  = 4
  Integer, Public, Parameter :: MSG_VERBOSE2 = 5
  Integer, Public, Parameter :: MSG_DEBUG    = 6

  !* Access specifiers to routines *!

  Public :: fts_log_init
  Public :: fts_log_exit
  Public :: fts_log
  Public :: fts_logWithInfo

  Public :: fts_log_GetVerbosity
  Public :: fts_log_GetUnit

  !* Routines *!

Contains

  !> Initialize the module.
  !!
  !! Set the different unit used to log usual, error and warning messages.
  !! To deactivate an output, give a negative value to the unit.
  !!
  !! @param[in] stdunit   Fortran Unit used to log ftsolver messages 
  !! @param[in] errunit   Fortran Unit used to log ftsolver errors   
  !! @param[in] warnunit  Fortran Unit used to log ftsolver warnings
  !! @param[in] verblvl   wanted level of verbosity 
  !! @param[in] rank      MPI rank (prefixing all messages)
  !! 
  !! @note 
  !! This routine would be called at the begining of ftsolver_driver.
  !! (ie, each time the user call ftsolver)
  Subroutine fts_log_init ( errunit, warnunit, stdunit, verblvl, logrank )

    !* Arguments *!
    
    Integer, Intent(in) ::  errunit
    Integer, Intent(in) ::  warnunit
    Integer, Intent(in) ::  stdunit
    Integer, Intent(in) ::  verblvl
    Integer, Intent(in) ::  logrank

    !- End of header -----------------------------------------------------------

    unit_err  = errunit
    unit_warn = warnunit
    unit_std  = stdunit
    verbosity = verblvl
    rank      = logrank

  End Subroutine fts_log_init

  !> Exit the module.
  !!
  !! Reset the unit used to log usual, error and warning messages
  !! to their default values.
  !! 
  Subroutine fts_log_exit

    unit_err  = 0
    unit_warn = 0
    unit_std  = 6
    verbosity = -1

  End Subroutine fts_log_exit


  !> Log a message.
  !!
  !! @param[in] msg_class  The class of the message
  !! @param[in] msg        The message to be logged
  !! 
  Subroutine fts_log(msg_class, msg)

    !* Arguments *!

    Integer, Intent(in) ::  msg_class
    Character(len=*), Intent(in) :: msg
    
    !- End of header -----------------------------------------------------------

    If ( verbosity >= msg_class ) Then

       Select Case(msg_class)
       Case (MSG_ERROR  ); If(unit_err  >= 0) &
            Write(unit_err ,'(A,I5.5,2A)') "[",rank,"] ERROR: ",Trim(msg)
       Case (MSG_WARNING); If(unit_warn >= 0) &
            Write(unit_warn,'(A,I5.5,2A)') "[",rank,"] WARNING: ",Trim(msg)
       Case (MSG_STD    ); If(unit_std  >= 0) &
            Write(unit_std ,'(A,I5.5,2A)') "[",rank,"] :",Trim(msg)
       Case (MSG_VERBOSE,MSG_VERBOSE2); If(unit_std  >= 0) &
            Write(unit_std ,'(A,I5.5,2A)') "[",rank,"] :",Trim(msg)
       Case (MSG_DEBUG  ); If(unit_err  >= 0) &
            Write(unit_err ,'(A,I5.5,2A)') "[",rank,"] DEBUG: ",Trim(msg)
       End Select

    End If

  End Subroutine fts_log

  !> Log a message with a status as additional information.
  !!
  !! @param[in] msg_class  The class of the message
  !! @param[in] msg_stat   Integer specifying the status of the message
  !! @param[in] msg        The message to be logged
  !! 
  Subroutine fts_logWithInfo(msg_class,msg_stat, msg)

    !* Arguments *!

    Integer, Intent(in) ::  msg_class
    Integer, Intent(in) ::  msg_stat
    Character(len=*), Intent(in) :: msg
    
    !- End of header -----------------------------------------------------------

    If ( verbosity >= msg_class ) Then

       Select Case(msg_class)
       Case (MSG_ERROR  ); If(unit_err  >= 0) &
            Write(unit_err ,'(A,I5.5,2A,I15)')&
            "[",rank,"] ERROR: ",Trim(msg), msg_stat

       Case (MSG_WARNING); If(unit_warn >= 0) &
            Write(unit_warn,'(A,I5.5,2A,I15)')&
            "[",rank,"] WARNING: ",Trim(msg), msg_stat

       Case (MSG_STD    ); If(unit_std  >= 0) &
            Write(unit_std ,'(A,I5.5,2A,I15)')&
            "[",rank,"] :",Trim(msg), msg_stat

       Case (MSG_VERBOSE,MSG_VERBOSE2); If(unit_std  >= 0) &
            Write(unit_std ,'(A,I5.5,2A,I15)')&
            "[",rank,"] :",Trim(msg), msg_stat

       Case (MSG_DEBUG  ); If(unit_err  >= 0) &
            Write(unit_err ,'(A,I5.5,2A,I15)')&
            "[",rank,"] DEBUG: ",Trim(msg), msg_stat
       End Select

    End If

  End Subroutine fts_logWithInfo

  ! [+] routine : fts_log_GetVerbosity ---------------------------------------
  !
  !> Provide to the user the verbosity level
  !!
  !! @param[out] verblvl   The verbosity
  !! 
  Subroutine fts_log_GetVerbosity(verblvl)

    !* Arguments *!

    Integer, Intent(out) ::  verblvl

    !- End of header -----------------------------------------------------------
    
    verblvl = verbosity

  End Subroutine fts_log_GetVerbosity


  ! [+] routine : fts_log_GetUnit -------------------------------------------
  !
  !> Provide to the user the unit associated to a class of messages.
  !!
  !! @param[in ] msg_class the class of messages
  !! @param[out] unit      the unit associated to the class.
  !! 
  Subroutine fts_log_GetUnit(msg_class, unit )

    !* Arguments *!

    Integer, Intent(in ) ::  msg_class
    Integer, Intent(out) ::  unit

    !- End of header -----------------------------------------------------------

    Select Case(msg_class)
    Case (MSG_ERROR  ); unit = unit_err
    Case (MSG_WARNING); unit = unit_warn
    Case (MSG_STD    ); unit = unit_std
    Case (MSG_VERBOSE); unit = unit_std
    Case (MSG_VERBOSE2); unit = unit_std
    Case (MSG_DEBUG  ); unit = unit_err
    Case Default      ; unit = unit_std 
    End Select

  End Subroutine fts_log_GetUnit


End Module fts_log_mod
