#include "fts_defs_f.h"

Module FTS_toolkit_parser_mod
  Implicit None

  ! Constants
  Integer, Parameter :: charlen = FTSOLVER_STRL
  Integer, Parameter :: noplist = 5 
  Character(len=charlen), Parameter :: oplist (noplist)= &
       (/ "==",">=","<=","< ","> " /)

  Public  :: FTS_read_array
  Public  :: FTS_parser_read_string
  Private :: FTS_getIndex
  Private :: FTS_getRValue
  Public  :: FTS_test_getIndex

  Public  :: FTS_Parser_getAssertions
  Public  :: FTS_Parser_splitAssertion
  Public  :: FTS_Parser_Reval
  Public :: FTS_Parser_eval_operator
  
Contains 

  ! [+] routine : FTS_test_print_status ------------------------------------------
  !
  !> Print the status of a routine.
  !!
  Subroutine FTS_test_print_status (rname,isSuccess)
    Implicit None

    Character(len=charlen) :: rname
    Logical                :: isSuccess

    If (isSuccess) Then
       Write(*,*) rname, "[PASS]"
    Else
       Write(*,*) rname, "[FAIL]"
    End If
  End Subroutine FTS_test_print_status

  ! [+] routine : FTS_parser_read_array ------------------------------------------
  !
  !> Parse the data from stream "UNIT" to update the values 
  !! of "iarray" or "rarray".
  !!
  !! Read from the stream "UNIT", several arrays values.
  !! Where an array value is in the form :
  !! 
  !! RKEY(INDEX) = VALUE ! For rarray 
  !! IKEY(INDEX) = VALUE ! For iarray
  !!
  !!----
  !!
  !! @param  input  the INPUT string
  !! @param  key    the key
  !! @return INDEX the index
  !!
  !!----
  !! @author Yohan Lee-tin-yien
  !!
  Subroutine FTS_read_array(niarray, nrarray, iarray, rarray, unit, ikey, rkey)

    !* Args *!
    Integer :: niarray
    Integer :: nrarray
    Integer :: iarray(niarray)
    Real*8 :: rarray(nrarray)
    Integer :: unit
    Character(len=charlen) :: ikey
    Character(len=charlen) :: rkey

    !* Locals *!
    Character(len=charlen) :: input 
    Integer :: istat
    Logical :: isIARRAY
    Logical :: isRARRAY
    !- End of header ---------------------------------------------------------

    istat=0

    Do While (istat == 0)
       !--------------------------------------------------
       ! [1.0] Read
       !--------------------------------------------------
       Read(unit,'(A)',IOSTAT= istat ) input

       !--------------------------------------------------
       ! [2.0] Jump irrelevant entries
       !--------------------------------------------------

       ! jump comments
       If ( input(1:1) == '#' ) Cycle

       ! test if iarray or rarray
       isIARRAY = (Index(input, Trim(ikey)) == 1 )
       isRARRAY = (Index(input, Trim(rkey)) == 1 )

       ! jump if not
       If ((.not. isIARRAY) .and. (.not. isRARRAY )) Cycle

       !--------------------------------------------------
       ! [3.0] Append the entry
       !--------------------------------------------------

       If (isIARRAY) iarray(FTS_getIndex(input)) = INT(FTS_getRValue(input))
       If (isRARRAY) rarray(FTS_getIndex(input)) = FTS_getRValue(input)

    End Do

    !--------------------------------------------------
    ! [4.0] Rewind the file
    !--------------------------------------------------
    Rewind(unit)

  End Subroutine FTS_read_array


  ! [+] routine : FTS_parser_read_string ------------------------------------------
  !
  !> Parse the stream "UNIT" to  get the string "STRING"
  !! with keyword 'KEYSTRING'
  !! 
  !! KEYSTRING   = STRING 
  !!
  !! If multiple KEYSTRING = STRING is found,
  !! It returns the last one.
  !!
  !!----
  !!
  !! @param  unit   the stream unit
  !! @param  keystr the key of the string
  !! @param  string the string 
  !! @param  info   routine status 0 if found, 1 if not
  !!
  !!----
  !! @author Yohan Lee-tin-yien
  !!
  Subroutine FTS_parser_read_string &
       (unit, keystr, string)

    !* Args *!
    Integer               , Intent(in)  :: unit
    Character(len=charlen), Intent(in)  :: keystr
    Character(len=charlen), Intent(out) :: string

    !* Locals *!
    Character(len=charlen) :: line 
    Integer :: istat
    Logical :: isSTRING
    !- End of header ---------------------------------------------------------

    istat=0
    string= ""
    
    Do While (istat == 0)
       !--------------------------------------------------
       ! [1.0] Read
       !--------------------------------------------------
       Read(unit,'(A)',IOSTAT= istat ) line

       !--------------------------------------------------
       ! [2.0] Jump irrelevant entries
       !--------------------------------------------------

       ! jump comments
       If ( line(1:1) == '#' ) Cycle

       ! test if 
       isSTRING = ( Index(line, Trim(keystr)) == 1 )
#if DEBUG
       write(*,*) "line  = ", Trim(line)
       write(*,*) "keystr= ", Trim(keystr)
       write(*,*) "index = ", Index(line, Trim(keystr))
       write(*,*) "isSTRING =", isSTRING 
#endif

       ! jump if not
       If (.not. isSTRING) Cycle

       ! save if found
       If ( isSTRING) string = line(Index(line,"=")+1:)

    End Do

    ! remove the empty spaces before the string
    If ( Trim(string) /= "" ) string = Trim(AdjustL(string))

    !--------------------------------------------------
    ! [4.0] Rewind the file
    !--------------------------------------------------
    Rewind(unit)

  End Subroutine FTS_parser_read_string


  ! [+] routine : FTS_getIndex ---------------------------------------------------
  !
  !> Get the INPUT from INPUT where input is :
  !! INPUT = ".*(INDEX).*"
  !!
  !!----
  !!
  !! @param  input  the INPUT string
  !! @param  key    the key
  !! @return INDEX the index
  !!
  !!----
  !! @author Yohan Lee-tin-yien
  !!
  Integer Function FTS_getIndex(input)
    Implicit None
    Character(len=charlen) :: input
    Integer :: start
    Integer :: end

    start= Scan(input,"(") + 1
    If( start == 1) FTS_getIndex = -1
    If( start /= 1) end = Scan(input,")") - 1

#if DEBUG
    If( (start /= 1) .and. (end /= -1) ) &
         write(*,*) input(start:end)
#endif

    If( (start /= 1) .and. (end /= -1) ) &
         read(input(start:end),*) FTS_getIndex

  End Function FTS_getIndex

  ! its tester
  Subroutine FTS_test_getIndex
    Implicit None
    Character(len=charlen) :: input = "INCTL(205)"
    Integer, Parameter     :: expected = 205
    Character(len=charlen) :: rname = "test_FTS_getIndex"
    Logical                :: isSuccess

    isSuccess=( FTS_getIndex(input) == expected )
    Call FTS_test_print_status(rname, isSuccess)
  End Subroutine FTS_test_getIndex


  ! [+] routine : FTS_getRValue ---------------------------------------------------
  !
  !> Get the VALUE from INPUT where input is :
  !! INPUT = "* = VALUE"
  !!
  !!----
  !!
  !! @param  input  the INPUT string
  !! @return Value  the value (real*8)
  !!
  !!----
  !! @author Yohan Lee-tin-yien
  !!
  Real(kind=8) Function FTS_getRValue(input)
    Implicit None
    Character(len=charlen) :: input
    Integer :: pos
    pos = scan(input,"=") + 1
    read(input(pos:),*) FTS_getRValue 
  End Function FTS_getRValue

  ! [+] routine : FTS_Parser_getAssertions ----------------------------------------
  !
  !> Get the list of assertions from stream "unit".
  !!
  !! ASSERT assertion1
  !! ASSERT assertion2
  !!
  !!----
  !!
  !! @param[in ] unit    the unit stream
  !! @param[out] nassert number of assertions found
  !! @param[out] assert  list of asssertions
  !!
  !!----
  !! @author Yohan Lee-tin-yien
  !!
  Subroutine FTS_Parser_getAssertions(unit, nassert, assert, info)
    Implicit None

    !* Arguments *!
    Integer, Intent(in)  :: unit
    Integer, Intent(out) :: nassert
    Character(len=charlen), Pointer, Intent(out) :: assert (:)
    Integer, Intent(out) :: info

    !* Local variables *!
    Integer :: istat
    Integer :: i
    Character(len=charlen) :: input

    !- End of header ---------------------------------------------------------

    !-------------------------------------------------------------------------
    ! [1.0] Get nassert
    !-------------------------------------------------------------------------
    istat=0
    nassert=0
    input=""
    Do While (istat == 0)
       If (input(1:Len('ASSERT')) == 'ASSERT') nassert = nassert + 1
       Read(unit,'(A)',IOSTAT= istat ) input
    End Do
    Rewind(unit)

    !-------------------------------------------------------------------------
    ! [2.0] Create assert(:)
    !-------------------------------------------------------------------------

    Allocate(assert(nassert),stat=info)
    If (info /= 0) Return
    assert= ""

    !-------------------------------------------------------------------------
    ! [3.0] Populate assert(:)
    !-------------------------------------------------------------------------
    input=""
    i= 0
    istat=0
    Do While (istat == 0)
       If (input(1:Len('ASSERT')) == 'ASSERT') Then
          i = i + 1
          assert(i) = input(Len('ASSERT')+1:Len_Trim(input))
       End If
       Read(unit,'(A)',IOSTAT= istat ) input
    End Do
    Rewind(unit)

  End Subroutine FTS_Parser_getAssertions

  ! [+] routine : FTS_Parser_splitAssertion --------------------------------------
  !
  !> Split an Assertion into different part 
  !!
  !! ELEM1 OP ELEM2 [+- ELEM3 ]
  !!
  !! If a field or an element is not found return empty string.
  !!----
  !!
  !! @param[in ] assertion the assertion to process
  !! @param[out] op        the operator found is 
  !! @param[out] elem1     the first element
  !! @param[out] elem2     the second element
  !! @param[out] elem3     the fird element
  !!
  !!----
  !!
  !! @author Yohan Lee-tin-yien
  !!
  Subroutine FTS_Parser_splitAssertion &
       (assertion                  &
       ,op,elem1,elem2,elem3       &
       )
    Implicit None

    !* Arguments *!
    Character(len=charlen) , Intent(in ) :: assertion
    Character(len=charlen) , Intent(out) :: op,elem1,elem2,elem3

    !* Local variables *!
    Integer :: i,j
    Integer :: op_start, op_end
    Integer :: elem1_start, elem1_end
    Integer :: elem2_start, elem2_end
    Integer :: elem3_start, elem3_end

    !- End of header ---------------------------------------------------------

    !-------------------------------------------------------------------------
    ! [1.0] Init
    !-------------------------------------------------------------------------
    op=""
    elem1=""
    elem2=""
    elem3=""

    op_start=0
    op_end=0
    elem1_start=0
    elem1_end=0
    elem2_start=0
    elem2_end=0
    elem3_start=0
    elem3_end=0

    !-------------------------------------------------------------------------
    ! [2.0] If present - Get elem3
    !-------------------------------------------------------------------------
    ! detect
    elem3_start = Index(assertion,'+-')
    have_elem3: If (elem3_start /= 0 ) Then
       ! adjust the begining
       elem3_start = elem3_start+Len('+-')
       elem3_end   = Len_Trim(assertion)
       elem3=assertion(elem3_start:elem3_end)
    End If have_elem3

    !-------------------------------------------------------------------------
    ! [3.0] Get operator "op"
    !-------------------------------------------------------------------------

    Do i=1,noplist
       j =INDEX(assertion,Trim(oplist(i)))

       have_op:If (j /= 0) Then
          op_start= j
          op_end=op_start+Len_Trim(oplist(i))
          op=assertion(op_start:op_end)
       End If have_op

    End Do

    !-------------------------------------------------------------------------
    ! [3.0] Get elem1
    !-------------------------------------------------------------------------

    If( op_start /= 0 ) elem1_start = 1
    If( op_start /= 0 ) elem1_end   = op_start - 1 

    elem1 = assertion(elem1_start:elem1_end)

    !-------------------------------------------------------------------------
    ! [4.0] Get elem2
    !-------------------------------------------------------------------------

    have_elem2:If( op_start /= 0 )Then
       elem2_start = op_end + 1
       If (elem3_start /=0) elem2_end = elem3_start - Len('+-') -1
       If (elem3_start ==0) elem2_end = Len_Trim(assertion)
    End If have_elem2

    elem2 = assertion(elem2_start:elem2_end)

  End Subroutine FTS_Parser_splitAssertion

  !> evaluate an element contained in a string
  !!
  !! @warning comm and master are unused in current version.
  !! @todo add comm and master handling.
  Real*8 Function FTS_Parser_Reval&
       (comm, master,&
       nicntl,icntl, nrcntl,rcntl, &
       nikeep,ikeep, nrkeep,rkeep,  &
       niinfo,iinfo, nrinfo,rinfo, &
       niinfog,iinfog, nrinfog,rinfog, &
       string) 

    !* Dependencies *!
    Implicit None
    Include 'mpif.h'

    !* Arguments *!
    Integer, Intent(in) :: comm,master

    Integer, Intent(in) :: nicntl, nrcntl
    Integer, Intent(in) :: nikeep, nrkeep
    Integer, Intent(in) :: niinfo, nrinfo
    Integer, Intent(in) :: niinfog, nrinfog

    Integer, Intent(in) :: icntl (nicntl )
    Integer, Intent(in) :: ikeep (nikeep ) 
    Integer, Intent(in) :: iinfo (niinfo )
    Integer, Intent(in) :: iinfog(niinfog)

    Real(kind=8), Intent(in) :: rcntl (nrcntl )
    Real(kind=8), Intent(in) :: rkeep (nrkeep )
    Real(kind=8), Intent(in) :: rinfo (nrinfo )
    Real(kind=8), Intent(in) :: rinfog(nrinfog)


    Character(len=charlen) :: string

    !* Local variables *!
    Character(len=charlen) :: s
    Integer :: i
    Integer :: istat
    !- End of header ---------------------------------------------------------

    !-------------------------------------------------------------------------
    ! [1.0] Init
    !-------------------------------------------------------------------------

    FTS_Parser_Reval=0.d0

    ! If string contains empty, return
    If (Len_Trim(string) == 0) Goto 9999

    ! If string contains a real, read it and return
    Read(string,*,iostat=istat) FTS_Parser_Reval
    If (istat == 0 ) Goto 9999

    ! 
    s=ADJUSTL(string)

    !-------------------------------------------------------------------------
    ! [2.0] Get the Index
    !-------------------------------------------------------------------------

    ! get index
    i = FTS_getIndex(s)
    ! check index
    If (i < 1) Goto 9999

    !-------------------------------------------------------------------------
    ! [2.0] Give the value
    !-------------------------------------------------------------------------

    IF((0 /= INDEX(s,'ICNTL(')).and.(i <= nicntl)) FTS_parser_Reval = icntl(i)
    IF((0 /= INDEX(s,'RCNTL(')).and.(i <= nrcntl)) FTS_parser_Reval = rcntl(i)

    IF((0 /= INDEX(s,'IKEEP(')).and.(i <= nikeep)) FTS_parser_Reval = ikeep(i)
    IF((0 /= INDEX(s,'RKEEP(')).and.(i <= nrkeep)) FTS_parser_Reval = rkeep(i)

    IF((0 /= INDEX(s,'IINFO(')).and.(i <= niinfo)) FTS_parser_Reval = iinfo(i)
    IF((0 /= INDEX(s,'RINFO(')).and.(i <= nrinfo)) FTS_parser_Reval = rinfo(i)

    IF((0 /= INDEX(s,'IINFOG(')).and.(i <= niinfog)) FTS_parser_Reval = iinfog(i)
    IF((0 /= INDEX(s,'RINFOG(')).and.(i <= nrinfog)) FTS_parser_Reval = rinfog(i)

    !-------------------------------------------------------------------------
    ! [*.0] Exit
    !-------------------------------------------------------------------------

9999 Continue
    Return

  End Function FTS_Parser_Reval

  !> evaluate the operator op
  !!

  Integer Function FTS_Parser_eval_operator&
       (v1,v2, op, have_more_or_less, v3)

    Implicit None

    !* Arguments *!
    Real(kind=8) :: v1
    Real(kind=8) :: v2
    Character(len=charlen) :: op

    Integer :: have_more_or_less
    Real(kind=8) :: v3

    !* Local variable *!
    Logical :: test

    If ( have_more_or_less == 0 ) Then

       !-----------------------------------------
       ! [1.0] no +- 
       !-----------------------------------------

       select case(Trim(op))
       case("=="); test = (v1 == v2)
       case("<="); test = (v1 <= v2)
       case(">="); test = (v1 >= v2)
       case default; test = .False.
       End select

    Else

       !-----------------------------------------
       ! [2.0] have +- 
       !-----------------------------------------

       select case(Trim(op))
       case("=="); test = ( Abs(v1 -v2) <= v3)
       case default; test = .False.
       End select

    Endif

    !-----------------------------------------
    ! [3.0] Form the return value
    !-----------------------------------------

    If(      test ) FTS_parser_eval_operator = 0
    If(.not. test ) FTS_parser_eval_operator = -1

  End Function FTS_Parser_eval_operator


End Module FTS_toolkit_parser_mod
