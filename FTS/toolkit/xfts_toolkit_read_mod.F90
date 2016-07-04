! Warning: XFTS_GENFILE_COMMENT
#include "fts_defs_f.h"
! [+] module : XFTS_toolkit_read_mod -------------------------------------------
!
!> This module contains useful routines
!! to read files.
!!
Module XFTS_toolkit_read_mod

      
!* No Implicit Typing *!

      Implicit None
  
     !* Access specifiers *!
      
      Public :: XFTS_read_matrix
      Public :: XFTS_read_rhs
      Public :: XFTS_read_param_fixedformat
      Public :: XFTS_read_param_freeformat
      Public :: XFTS_check_assertions      
      !* Implementations
      
      Contains

  ! [+] routine : read_matrix --------------------------------------------------
  !
  !> Get the Input Matrix according to its file extension.
  !!
  !! @param[in,out] smatrix      the input matrix
  !! @param[   out] info         the routine status
  !!
  !! @author Yohan Lee-tin-yien
  !!
      Subroutine XFTS_read_matrix ( smatrix , filename )

    !* Modules used *!

      Use XFTS_sparse_matrix_mod  
      Implicit none
      Include 'mpif.h'

    !* Subroutine Arguments *!
      Type(XFTS_sparse_matrix_t), Intent(inout) :: smatrix
      Character(len=*)     , Intent(in   ) :: filename
      

    !* local variables *!

    ! Constants
      Integer, Parameter :: stderr = 6

    ! Scalars
      Integer :: info
      Integer :: unit 
      Logical :: have_mpi

    ! Strings
      Character(len=FTSOLVER_STRL ) :: fext

    !- End of header ------------------------------------------------------

    ! read according to file extension

      info = 0
      unit = 1111
      fext = filename(  Index(filename,'.',BACK=.True.) : Len_trim(filename) )
      Open(UNIT=unit, FILE= filename, STATUS='OLD', ACTION='READ')

      
      Select Case(Trim(fext)) 
      Case ('.ijv')             ! Coordinate       
         Call XFTS_sm_ijvread (smatrix, unit, info ) 
      Case ('.mm','.mtx')       ! Matrix Market 
         Call XFTS_sm_mmread(smatrix, unit, info ) 
      Case ('.rsa','.rua','.rse','.rue','.csa','.cua','.cse','.cue') ! Harwell boeing 
         Call XFTS_sm_hbread(smatrix, unit, info ) 
      Case ('.psa','.pua','.pse','.pue') ! Harwell boeing, unsupported
         Write(stderr,*) "Error : unsupported Harwell boeing file"
      Case Default
         Write(stderr,*) "Error : unknown file extension"
      End Select
      Close (unit)

      If (info <0) Then
         Write(stderr,*) "Error : while reading the matrix, aborting." 
         Call MPI_Initialized(have_mpi,info)
         If(have_mpi) Call MPI_Abort(info,info)
         Stop "error while reading inputs."
      End If

      End Subroutine XFTS_read_matrix
      

  ! [+] routine : read_rhs - --------------------------------------------------
  !
  !> Get the Input RHS, in a file in ijv format .
  !!
  !! @param[in,out] rhs          the input rhs
  !! @param[   out] info         the routine status
  !! @param[in    ] filename     the rhs filename
  !!
  !! @author Yohan Lee-tin-yien
  !!
  Subroutine XFTS_read_rhs ( rhs, info , filename )

    !* Modules used *!
      Use XFTS_sparse_matrix_mod  
    Implicit none

    !* Subroutine Arguments *!
    XFTS_FLOAT, Pointer, Intent(  out) :: rhs (:)
    Integer              , Intent(  out) :: info
    Character(len=*)     , Intent(in   ) :: filename

    !* local variables *!

    ! Scalars
    Integer :: iinfo
    Integer :: unit 
    FTS_INT :: k, n, row, col
    XFTS_FLOAT :: val
    ! Strings
    Type(XFTS_sparse_matrix_t) :: sm_rhs

    !- End of header ------------------------------------------------------

    ! [1] init
    unit = 1111
    Nullify( rhs )
    Call XFTS_sm_nullify (sm_rhs, iinfo)
    If (iinfo /= 0) Goto 9999

    ! [2] read data
    Open(UNIT=unit, FILE=filename,STATUS='OLD', ACTION='READ', IOSTAT = iinfo )
    If (iinfo /= 0) Goto 9999
    Call XFTS_sm_ijvread (sm_rhs, unit, iinfo ) 
    If (iinfo /= 0) Goto 9999
    Close( unit )
    
    ! [4] RHS = the first column of sm_rhs
    n = sm_rhs%m
    Allocate( rhs( n ), STAT = iinfo )
    If (iinfo /= 0) Goto 9999
    Do k=1,n
       rhs(k)= XFTS_FLOATZERO
    End Do

    Do k=1,sm_rhs%nnz
       row = sm_rhs%i(k)
       col = sm_rhs%j(k)
       val = sm_rhs%v(k)

       If ( col  == 1 ) rhs(row) = val
    End Do

#if FTSOLVER_DEBUG
    Do k=1,n
      Write(11,*) k, rhs(k)
    End Do
#endif    

    ! [5] Finish
9999 Continue
    If ( iinfo == 0 ) info =  0
    If ( iinfo /= 0 ) info = -1

    If (info < 0 )Then
       If (Associated( rhs )) Deallocate( rhs )
    End If

    Call XFTS_sm_free( sm_rhs , iinfo )

  End Subroutine XFTS_read_rhs

  !
  !-----------------------------------------------------------------------------
  !
  !> Write a vector in a file in ijv format.
  !! 
  !! do nothing on empty filename.
  !!
  !! @param[in    ] n            the size of the vector
  !! @param[in    ] vect         the vector
  !! @param[   out] info         the routine status
  !! @param[in    ] filename     the filename
  !!
  !! @author Yohan Lee-tin-yien
  !!
  Subroutine XFTS_write_vect ( vect, n, info , filename )

    !* Modules used *!
    Implicit none

    !* Arguments *!
    FTS_INT, Intent(  in) :: n
    XFTS_FLOAT, Intent(  in) :: vect (n)
    Integer              , Intent(  out) :: info
    Character(len=*)     , Intent(in   ) :: filename

    !* local variables *!

    Integer :: unit 
    FTS_INT :: k

    !- End of header ------------------------------------------------------

    ! [1] init
    info = FTS_SUCCESS

    ! exit early 
    If (Len_Trim(filename) == 0) Return
    
    unit = 1111
    Open(UNIT=unit, FILE=filename, ACTION='WRITE', IOSTAT = info )
    If (info /= 0)Then
       Write(0,*) "Cannot open file:", Trim(filename)
       info = - __LINE__
       Return
    End If

    ! [2] write 
    ! header
    Write(unit,*) n, 1, n

    ! data
    Do k=1,n
      Write(unit,*) k, 1, vect(k)
    End Do

    ! [5] Finish
    Close( unit )

  End Subroutine XFTS_write_vect


  ! [+] routine : read_param_fixedformat --------------------------------------------
  !
  !> Read the parameters of a simulation from a file in fixed format.
  !! Format description :
  !! 
  !!   DESCRIPTION OF VALUE 1
  !!   VALUE 1
  !!   DESCRIPTION OF VALUE 2
  !!   VALUE 2
  !!   ...
  !!
  !!----
  !!
  !! @author Yohan Lee-tin-yien
  !! @author Azzam Haidar
  !! @author Luc   Giraud
  !!
  !!----
  !!
  !! History
  !! - 14/01/11 : v0.1b : import from stable version. (Yohan Lee-tin-yien )
  !!                      
  Subroutine XFTS_read_param_fixedformat(unit,icntl,rcntl,sym,&
       matrixfile,rhsfile,initguessfile,topofile)

    Implicit None

    !* Arguments *!
    Integer, Intent(in   ) :: unit
    Integer, Intent(inout) :: icntl(FTSOLVER_ICNTL_SIZE)
    Real(kind=8), Intent(inout) :: rcntl(FTSOLVER_RCNTL_SIZE)
    Integer, Intent(  out) :: sym
    Character(len=FTSOLVER_STRL), Intent(out) :: matrixfile
    Character(len=FTSOLVER_STRL), Intent(out) :: rhsfile
    Character(len=FTSOLVER_STRL), Intent(out) :: initguessfile
    Character(len=FTSOLVER_STRL), Intent(out) :: topofile

    !- End of header ------------------------------------------------------
    matrixfile=""
    rhsfile=""
    initguessfile=""
    topofile=""

    Read(unit,*) 
    Read(unit,'(a)') matrixfile
    Read(unit,*) 
    Read(unit,'(a)') rhsfile
    Read(unit,*) 
    Read(unit,'(a)') initguessfile
    Read(unit,*) 
    Read(unit,*) sym
    Read(unit,*) 
    Read(unit,*) icntl(7)
    Read(unit,*) 
    Read(unit,*) icntl(13)
    Read(unit,*) 
    Read(unit,*) icntl(21)
    Read(unit,*) 
    Read(unit,*) rcntl(11)
    Read(unit,*) 
    Read(unit,*) rcntl(21)
    Read(unit,*) 
    Read(unit,*) icntl(20)
    Read(unit,*) 
    Read(unit,*) icntl(24)
    Read(unit,*) 
    Read(unit,*) icntl(27)
    Read(unit,*) 
    Read(unit,*) icntl(26)
    Read(unit,*) 
    Read(unit,*) icntl(22)
    Read(unit,*) 
    Read(unit,*) icntl(28)
  End Subroutine XFTS_read_param_fixedformat


  ! [+] routine : read_param_freeformat --------------------------------------------
  !
  !> Read the parameters of a simulation from a file in free format.
  !! Format description :
  !! 
  !!   # a COMMENT 
  !!   INCTL(XX) = YY
  !!   ...
  !!
  !!----
  !!
  !! @author Yohan Lee-tin-yien
  !!
  !!----
  !!
  !! History
  !! - 25/02/11 : v0.1a : write routine
  !! - 03/03/11 : v0.1a : add init guess file
  !!
  Subroutine XFTS_read_param_freeformat(&
       unit,icntl,rcntl,sym,job,&
       matrixfile,rhsfile,initguessfile,&
       outrhsfile,outsolfile)

    Use XFTS_ftsolver_type
    Use FTS_toolkit_parser_mod
    Implicit None

    !* Arguments *!
    Integer, Intent(in   ) :: unit
    Integer, Intent(inout) :: icntl(FTSOLVER_ICNTL_SIZE)
    Real(kind=8), Intent(inout) :: rcntl(FTSOLVER_RCNTL_SIZE)
    Integer, Intent(  out) :: sym
    Integer, Intent(  out) :: job
    Character(len=FTSOLVER_STRL), Intent(out) :: matrixfile
    Character(len=FTSOLVER_STRL), Intent(out) :: rhsfile
    Character(len=FTSOLVER_STRL), Intent(out) :: initguessfile
    Character(len=FTSOLVER_STRL), Intent(out) :: outrhsfile
    Character(len=FTSOLVER_STRL), Intent(out) :: outsolfile

    !* Local variables *!

    Character(len=charlen)  :: keystr
    Character(len=charlen)  :: string
    Character(len=charlen)  :: ikey
    Character(len=charlen)  :: rkey


    !- End of header ------------------------------------------------------

    outrhsfile=""
    outsolfile=""
    matrixfile=""
    initguessfile=""
    rhsfile=""

    ! Read strings
    keystr="MATFILE"
    Call FTS_parser_read_string (unit, keystr, string )
    matrixfile = Trim(string)

    keystr="RHSFILE"
    Call FTS_parser_read_string (unit, keystr, string )
    rhsfile = Trim(string)

    keystr="INITGUESSFILE"
    Call FTS_parser_read_string (unit, keystr, string )
    initguessfile = Trim(string)

    keystr="OUTRHSFILE"
    Call FTS_parser_read_string (unit, keystr, string )
    outrhsfile = Trim(string)

    keystr="OUTSOLFILE"
    Call FTS_parser_read_string (unit, keystr, string )
    outsolfile = Trim(string)

    ! Read symmetry
    keystr="SYM"
    Call FTS_parser_read_string (unit, keystr, string )
    If( Len_Trim(string) == 0 )Then
       sym = -1
    Else
       Read(string,*) sym
    End If

    ! Read job
    keystr="JOB"
    Call FTS_parser_read_string (unit, keystr, string )
    If( Len_Trim(string) == 0 )Then
       job = -1
    Else
       Read(string,*) job
    End If


    ! Read parameters
    ikey="ICNTL"
    rkey="RCNTL"
    Call FTS_read_array(FTSOLVER_ICNTL_SIZE, FTSOLVER_RCNTL_SIZE, &
         icntl, rcntl, unit, ikey, rkey)

  End Subroutine XFTS_read_param_freeformat


  ! [+] routine : check_assertions ---------------------------------------------------
  !
  !> Apply all assertions stored in the stream "UNIT"
  !!
  !! 
  !! @param ftsl the maphys instance to check 
  !! @param unit the unit of file containings the assertion list
  !! @param info routine status 
  !!             - 0     = all assertion are true, 
  !!             - other = at least one  assertion are false
  !!
  Subroutine XFTS_check_assertions(ftsl,unit, info) 

    !* Modules *!
    Use XFTS_ftsolver_mod
    Use XFTS_ftsolver_type
    Use FTS_toolkit_parser_mod
    Implicit None
    Include 'mpif.h'

    !* Arguments
    Type(XFTS_ftsolver_t), Intent(in) :: ftsl 
    Integer       , Intent(in) :: unit 
    Integer       , Intent(out) :: info

    !* Local variables *!

    ! constants
    Integer, Parameter :: master = 0

    ! scalars
    Integer                :: comm
    Integer                :: linfo
    Integer                :: nassert
    Integer                :: i,j
    Real(kind=8)           :: relem, r1,r2,r3
    Integer                :: have_more_or_less
    Logical                :: test

    ! Strings
    Character(len=charlen) :: elem1,elem2,elem3,op
    Character(len=charlen) :: elem
    Character(len=charlen) :: keystr
    Character(len=charlen) :: msg

    Character(len=charlen) :: desc

    ! Arrays
    Character(len=charlen), Pointer :: assert (:)

    !- End of header -------------------------------------------------------------

    !-----------------------------------------------------------------------------
    ! [1.0] Init
    !-----------------------------------------------------------------------------
    comm = ftsl%comm
    keystr="DESCRIPTION"
    Call FTS_parser_read_string (unit, keystr, desc )


    !-----------------------------------------------------------------------------
    ! [2.0] Get assertions list
    !-----------------------------------------------------------------------------

    ! get
    Call FTS_Parser_getAssertions(unit, nassert, assert, linfo)

    ! exit early if no assert asked
    If ( nassert < 1 ) Then
       info = 0
       Return
    End If

    !-----------------------------------------------------------------------------
    ! [3.0] Treat each assertion
    !-----------------------------------------------------------------------------
    Do i=1, nassert

#if VERBOSE
       Write(*,*) "Asserting :", Trim(assert(i))
#endif 

       !----------------------------------------------------------------
       ! [2.1] Split the assertion string
       !----------------------------------------------------------------
       Call FTS_Parser_splitAssertion (assert(i),op,elem1,elem2,elem3)

       !----------------------------------------------------------------
       ! [2.2] Evaluate the elements
       !----------------------------------------------------------------
       Do j=1,3
          If( j == 1) elem = elem1
          If( j == 2) elem = elem2
          If( j == 3) elem = elem3

          relem=FTS_parser_Reval&
               (ftsl%comm, 0,&
               FTSOLVER_ICNTL_SIZE, ftsl%icntl, &
               FTSOLVER_RCNTL_SIZE, ftsl%rcntl, &
               FTSOLVER_IKEEP_SIZE, ftsl%ikeep, &
               FTSOLVER_RKEEP_SIZE, ftsl%rkeep, &
               FTSOLVER_IINFO_SIZE, ftsl%iinfo, &
               FTSOLVER_RINFO_SIZE, ftsl%rinfo, &
               FTSOLVER_IINFOG_SIZE, ftsl%iinfog, & 
               FTSOLVER_RINFOG_SIZE, ftsl%rinfog, &
               elem)

          If( j == 1) r1   = relem
          If( j == 2) r2   = relem
          If( j == 3) r3   = relem
       End Do

       !----------------------------------------------------------------
       ! [2.3] Evaluate the assertion
       !----------------------------------------------------------------

       If ( Trim(elem3) == "" ) have_more_or_less = 0
       If ( Trim(elem3) /= "" ) have_more_or_less = 1

       linfo = FTS_parser_eval_operator(r1,r2,op, have_more_or_less, r3 )

       If ( linfo == 0) test = .True.
       If ( linfo /= 0) test = .False.

       If ( .not. test ) Goto 9999

    End Do

    !----------------------------------------------------------------
    ! [*.*] Finish
    !----------------------------------------------------------------


9999 Continue

    ! Handle Failed test  
    If ( have_more_or_less == 0 ) Then
       write(msg,'(A,A,G10.3,A,G10.3,A)') , &
            Trim(elem1)//Trim(op)//Trim(elem2), &
            " (",r1,Trim(op),r2,')'
    Else
       write(msg,'(A,A,G10.3,A,G10.3,A,G10.3,A)') , &
            Trim(elem1)//Trim(op)//Trim(elem2)//'+-'//Trim(elem3), &
            " (",r1,Trim(op),r2,"+-",r3,')'
    End If

    ! print message
    If ( test .eqv. .False. ) write(*,'(4A)') "[FAIL] ", desc(1:60), "::", Trim(msg)
    If ( test .eqv. .True.  ) write(*,'(4A)') "[PASS] ", desc(1:60), "::", Trim(msg)
    ! update routine status
    If ( test .eqv. .False. ) info = -1
    If ( test .eqv. .True.  ) info =  0

    ! Free memory
    Deallocate(assert)

  End Subroutine XFTS_check_assertions




End Module XFTS_toolkit_read_mod
