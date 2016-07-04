  ! handle timing 
  Module FTS_time_mod
    Implicit None
    
    Public :: FTS_time_start
    Public :: FTS_time_stop
    Public :: FTS_time_add
    Public :: FTS_time_reset

    Contains

      Subroutine FTS_time_start(time)
        Implicit None
        Real(Kind=8), Intent(out) :: time
        Real(Kind=8), External    :: MPI_Wtime
        time = MPI_Wtime()
      End Subroutine FTS_time_start

      Subroutine FTS_time_stop(time)
        Implicit None
        Real(Kind=8), Intent(inout) :: time
        Real(Kind=8), External    :: MPI_Wtime
        time = MPI_Wtime() - time
      End Subroutine FTS_time_stop

      Subroutine FTS_time_add(time,timetoadd)
        Implicit None
        Real(Kind=8), Intent(inout) :: time
        Real(Kind=8), Intent(in   ) :: timetoadd
        time = time + timetoadd
      End Subroutine FTS_time_add

      Subroutine FTS_time_reset(time)
        Implicit None
        Real(Kind=8), Intent(out) :: time
        time = 0.d0
      End Subroutine FTS_time_reset

  End Module FTS_time_mod
