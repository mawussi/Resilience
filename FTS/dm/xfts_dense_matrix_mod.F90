! Warning: XFTS_GENFILE_COMMENT
#include "fts_defs_f.h"
#include "fts_macros_f.h"

!> Module for dense matrices
!!
!!  Module containing routines to manipulate dense matrices
!!
Module XFTS_dense_matrix_mod

  !* Module(s) used *!
  Use FTS_error_mod
  Use XFTS_dense_matrix_type

  !* No implicit typing *!
  Implicit none

  !* Private constants *!
  Character(len=FTSOLVER_STRL), Private, Parameter :: &
       FLNAME = "XFTS_ARITHfts_dense_matrix_mod.F90"

  !* Defined routine(s) *!
  Public :: XFTS_transpose_matrix_1D ! lty: obsolete ?
  Public :: XFTS_dm_unsymmetrize

  Public :: XFTS_dm_nullify
  Public :: XFTS_dm_build
  Public :: XFTS_dm_create
  Public :: XFTS_dm_dup
  Public :: XFTS_dm_free
  Public :: XFTS_dm_copyBloc
  Public :: XFTS_dm_mmwrite
  Public :: XFTS_dm_mmwritecoo
  Public :: XFTS_dm_sizeof
  Public :: XFTS_dm_vectorproduct


Contains

  !> Transpose a real nxn matrix in 1d representation
  !!
  !! this is an optimized version (from Azzam Haidar).
  !! 
  !! 
  subroutine XFTS_transpose_matrix_1D(A,n)
    Implicit None
    FTS_INT  ,  intent(in)    :: n
    XFTS_FLOAT,  intent(inout) :: A(n*n)

    XFTS_FLOAT :: tmp
    integer :: i,j

    do i =  1, n
       do j = 1, i-1
          tmp          = A((j-1)*n+i)
          A((j-1)*n+i) = A((i-1)*n+j)
          A((i-1)*n+j) = tmp
       end do
    end do

    ! slow implementation (with intrinsic functions)
    ! A = pack( reshape(A,(/n,n/),order=(/2,1/)), .TRUE.)
  end subroutine XFTS_transpose_matrix_1D

  ! [+] routine : XFTS_dm_unsymmetrize ------------------------------------
  !
  !> Unsymmetrize a dense matrix
  !!
  !!
  !! @param[in,out ] M             the dense matrix to unsymmetrize
  !!
  !! @author Yohan Lee-tin-yien
  !!
  !! @version 0.1 

  Subroutine XFTS_dm_unsymmetrize( dm )

    !* Arguments *!

    Type(XFTS_dense_matrix_t), Intent(inout) :: dm

    !* Local variable(s) *!
    FTS_INT :: i, j, ij, ji,ld

    !- End of header------------------------------------------------------------

    !---------------------------------------------------------------------------
    ! [1] Exit early
    !---------------------------------------------------------------------------

    If (dm%sym == DM_SYM_IsGeneral ) Return

    dm%sym = DM_SYM_IsGeneral
    If (dm%stored == DM_STORED_ALL ) Return

    !---------------------------------------------------------------------------
    ! [2] Case only have lower
    !---------------------------------------------------------------------------

    If (dm%stored == DM_STORED_LOWER )Then 
       ld=dm%ld
       Do j=1,dm%n
          Do i=j+1,dm%m
             ij = (j-1)*ld + i
             ji = (i-1)*ld + j
             dm%v(ji) = dm%v(ij) 
          End Do
       End Do
    End If

    !---------------------------------------------------------------------------
    ! [3] Case only have upper
    !---------------------------------------------------------------------------

    If (dm%stored == DM_STORED_UPPER )Then 
       ld=dm%ld
       Do j=1,dm%n
          Do i=j+1,dm%m
             ij = (j-1)*ld + i
             ji = (i-1)*ld + j
             dm%v(ij) = dm%v(ji) 
          End Do
       End Do
    End If

    !---------------------------------------------------------------------------
    ! [4] Set flags
    !---------------------------------------------------------------------------
    
    ! storage
    dm%stored = DM_STORED_ALL

  End Subroutine XFTS_dm_unsymmetrize




  ! [+] routine : XFTS_dm_nullify  ----------------------------------------
  !
  !> Nullify a dense matrix
  !!
  !!
  !! @param[in,out ] M             the dense matrix to nullify
  !! @param[out    ] info          the routine status
  !!
  !!
  !! @author Yohan Lee-tin-yien
  !!
  !! @version 0.1 

  Subroutine XFTS_dm_nullify( & ! intents
       M   ,                      & ! inout
       info                       & ! out
       )
    Implicit None

    !* Subroutine arguments *!

    Type(XFTS_dense_matrix_t), Intent(inout) :: M
    Integer             , Intent(  out) :: info

    !- End of header------------------------------------------------------------

    !---------------------------------------------------------------------------
    ! [1] Set size to ZEROS
    !---------------------------------------------------------------------------
    !
    M%m  = 0 
    M%n  = 0 
    M%ld = 0 

    !---------------------------------------------------------------------------
    ! [2] Set pointers to NULL
    !---------------------------------------------------------------------------
    !
    Nullify(M%v)


    !---------------------------------------------------------------------------
    ! [3] Set others flags
    !---------------------------------------------------------------------------

    M%sym    = DM_SYM_IsGeneral
    M%stored = DM_STORED_ALL

#if XFTS_HAVE_ARITH_D
    M%arith  = DM_ARITH_IsD

#elif XFTS_HAVE_ARITH_S
    M%arith  = DM_ARITH_IsS

#elif XFTS_HAVE_ARITH_Z
    M%arith  = DM_ARITH_IsZ

#elif XFTS_HAVE_ARITH_C
    M%arith  = DM_ARITH_IsC
#endif

    !---------------------------------------------------------------------------
    ! [*.*] Handle errors
    !---------------------------------------------------------------------------
    !
    ! Everything went well
    info = 0

    ! Returning
    Return

  End Subroutine XFTS_dm_nullify


  ! [+] routine : XFTS_DM_Build -------------------------------------------
  !
  !> Build a dense matrix.
  !!
  !! Create dense matrix :
  !!   - of size mxn 
  !!   - stored in a vector of size ldxn.
  !!
  !! @param[in,out] DM           
  !!
  !!       the dense matrix to be created.
  !!
  !! @param[in    ] m            
  !!
  !!       the number of rows in the matrix.
  !!
  !! @param[in    ] n            
  !!
  !!        the number of columns in the matrix.
  !!
  !! @param[in    ] ld            
  !!
  !!       the leading dimmension of the matrix.
  !!
  !! @param[   out] info          
  !!
  !!       the routine status
  !!
  !!
  Subroutine XFTS_DM_Build(DM, m, n, ld, sym, stored, data, info)

    !* Arguments *!
    Type(XFTS_dense_matrix_t)         , Intent(inout) :: DM
    FTS_INT                   , Intent(in   ) :: m
    FTS_INT                   , Intent(in   ) :: n
    FTS_INT                   , Intent(in   ) :: ld
    FTS_INT                   , Intent(in   ) :: sym
    FTS_INT                   , Intent(in   ) :: stored
    XFTS_FLOAT, Pointer        , Intent(in   ) :: data (:)
    Integer                      , Intent(  out) :: info


    !* Local Variables *!

    !- End of header -------------------------------------------------


    !-----------------------------------------------------------------
    ! [1] Check arguments
    !-----------------------------------------------------------------

    info = 0
#if FTSOLVER_DEBUG
    If ( m  <= 0 ) info = -1
    If ( n  <= 0 ) info = -1
    If ( ld <  m ) info = -1
    If (info < 0)  Return
#endif

    !-----------------------------------------------------------------
    ! [2] Copy attributes
    !-----------------------------------------------------------------

    dm%m   = m
    dm%n   = n
    dm%ld  = ld
    dm%sym = sym
    dm%stored = stored
    
    !-----------------------------------------------------------------
    ! [3] Link data
    !-----------------------------------------------------------------
    
    dm%v => data
     
  End Subroutine XFTS_DM_Build



  ! [+] routine : XFTS_DM_Create ------------------------------------------
  !
  !> Create a dense matrix.
  !!
  !! Create dense matrix :
  !!   - of size mxn 
  !!   - stored in a vector of size ldxn.
  !!   - with initial values = 0
  !!   - and default initial flags.
  !!
  !! @param[in,out] DM           
  !!
  !!       the dense matrix to be created.
  !!
  !! @param[in    ] m            
  !!
  !!       the number of rows in the matrix.
  !!
  !! @param[in    ] n            
  !!
  !!        the number of columns in the matrix.
  !!
  !! @param[in    ] ld            
  !!
  !!       the leading dimmension of the matrix.
  !!
  !! @param[   out] info          
  !!
  !!       the routine status
  !!
  !!
  Subroutine XFTS_DM_Create(DM, m, n, ld, info)

    !* Module(s) used *!
    Use fts_error_mod
    Implicit None

    !* Arguments *!
    Type(XFTS_dense_matrix_t)         , Intent(inout) :: DM
    FTS_INT                   , Intent(in   ) :: m
    FTS_INT                   , Intent(in   ) :: n
    FTS_INT                   , Intent(in   ) :: ld
    Integer                      , Intent(  out) :: info

    !* Local Variables *!
    ! scalars
    FTS_INT :: k

    !- End of header -------------------------------------------------

    !-----------------------------------------------------------------
    ! [1] Check arguments
    !-----------------------------------------------------------------

    info = 0
    Nullify (DM%v)

    CHCKASSRT( m > 0, info)
    FTS_ONFAILURE_RETURN( info )

    CHCKASSRT( n > 0  , info)
    FTS_ONFAILURE_RETURN( info )

    CHCKASSRT( ld >= m, info)
    FTS_ONFAILURE_RETURN( info )

    !-----------------------------------------------------------------
    ! [2] Allocate memory
    !-----------------------------------------------------------------

    Allocate(DM%v(ld*n), stat = info)
    CHCKALLOC( info )
    FTS_ONFAILURE_RETURN( info )

    !-----------------------------------------------------------------
    ! [3] Give values
    !-----------------------------------------------------------------

    ! attributes
    DM%sym   = DM_SYM_IsGeneral
    DM%stored = DM_STORED_ALL

#if defined(FTSOLVER_ARITH_d)
    DM%arith  = DM_ARITH_IsD

#elif defined(FTSOLVER_ARITH_s)
    DM%arith  = DM_ARITH_IsS

#elif defined(FTSOLVER_ARITH_z)
    DM%arith  = DM_ARITH_IsZ

#elif defined(FTSOLVER_ARITH_c)
    DM%arith  = DM_ARITH_IsC

#endif

    ! Sizes
    DM%m      = m
    DM%n      = n
    DM%ld     = ld

    ! values
    Do k = 1,ld*n
       DM%v(k) = XFTS_FLOATZERO
    End Do
     
   End Subroutine XFTS_DM_Create

   
  ! [+] routine : XFTS_dm_dup   -------------------------------------
  !
  !> Duplicate a dense matrix
  !!
  !!
  !! @param[in    ] In            the dense matrix to copy
  !! @param[   out] Out           the copy of the  dense matrix 
  !! @param[out   ] info          the routine status
  !!
  !! @author Yohan Lee-tin-yien
  !!
  !! @version 0.1 
  !!
  Subroutine XFTS_dm_dup    ( & ! intents
       in,                         & ! in
       out,                        & ! out
       info                        & ! 
       )

    !* Module(s) *!

    Use fts_error_mod
    Implicit None

    !* Arguments *!

    Type(XFTS_dense_matrix_t), Intent(in   ) :: in
    Type(XFTS_dense_matrix_t), Intent(  out) :: out
    Integer             , Intent(  out) :: info

    !* Local variables *!

    ! Scalars
    FTS_INT :: k

    !- End of header ---------------------------------------------

    !-------------------------------------------------------------
    ! [1] Copy everything
    !-------------------------------------------------------------

    ! Copy the scalars
    out%ld      = in%ld      
    out%m       = in%m        
    out%n       = in%n        

    out%stored  = in%stored   
    out%sym     = in%sym     

    ! Allocate new arrays

    Nullify( out%v )
    Allocate(out%v( out%ld * out%n), STAT = info)
    CHCKASSRT( info == 0, info )
    If (info < 0 ) Return

    ! Copy the arrays
    Do k=1,(out%ld*out%n)
       out%v(k) =  in%v(k)
    End Do

  End Subroutine XFTS_dm_dup

  ! [+] routine : XFTS_DM_Free ------------------------------------------
  !
  !> Free a dense matrix.
  !!
  !! Free dense matrix :
  !!
  !! @param[in,out] DM           
  !!
  !!       the dense matrix to be freed.
  !!
  !! @param[   out] info          
  !!
  !!       the routine status
  !!
  !!
  Subroutine XFTS_DM_Free(DM, info)

    !* Arguments *!

    Type(XFTS_dense_matrix_t) , Intent(inout) :: DM
    Integer                   , Intent(  out) :: info

    !- End of header -------------------------------------------------

    IF( Associated(DM%v) ) Deallocate(DM%v, STAT = info)
    CHCKALLOC( info )

  End Subroutine XFTS_DM_Free

  ! [+] routine : XFTS_dm_copyBloc ----------------------------------------
  !
  !> Copy into a bloc B of size mxn from a bloc A,
  !! where 
  !!
  !! @param[in   ] m        the number of rows to be copied
  !! @param[in   ] n        the number of columns to be copied
  !! @param[in   ] ldA      the leading dimmension of bloc A
  !! @param[in   ] A        the bloc A
  !! @param[in   ] ldB      the leading dimmension of bloc A
  !! @param[inout] B        the bloc B
  !! @param[  out] info     the return status
  !!
  !! @author Yohan Lee-tin-yien
  !!
  Subroutine XFTS_dm_copyBloc(m,n,ldA,A,ldB,B,info)

    !* Arguments *!

    FTS_INT          , Intent(in   )  :: m
    FTS_INT          , Intent(in   )  :: n
    FTS_INT          , Intent(in   )  :: ldA
    XFTS_FLOAT        , Intent(in   )  :: A(ldA*n)
    FTS_INT          , Intent(in   )  :: ldB
    XFTS_FLOAT        , Intent(inout)  :: B(ldB*n)
    Integer             , Intent(  out)  :: info

    !* Local variables *!

    FTS_INT :: i
    FTS_INT :: j

    !- End of header -----------------------------------------------------------

    !---------------------------------------------------------------------------
    ! [1] Checks arguments
    !---------------------------------------------------------------------------

    info = 0

    If ( m   < 0 ) info = -1
    If ( n   < 0 ) info = -2
    If ( ldA < m ) info = -3
    If ( ldB < m ) info = -5

    If ( info < 0) Return

    !---------------------------------------------------------------------------
    ! [2] Perform the copy
    !---------------------------------------------------------------------------

    Do j=1,n
       Do i=1,m
          B(i+(j-1)*ldB) = A(i+(j-1)*ldA)
       End Do
    End Do

  End Subroutine XFTS_dm_copyBloc



  ! [+] routine : XFTS_dm_mmwrite -----------------------------------------
  !
  !> export to matrix market format
  !!
  !! write a dense matrix in a file with Matrix Market Format. 
  !!
  !! @param[in   ] DM       the matrix 
  !! @param[in   ] funit     the fortran unit stream
  !! @param[in   ] info      the routine status
  !!
  !!-----
  !!
  !! @note 
  !! See the specification of the format at :
  !! http://people.sc.fsu.edu/~jburkardt/pdf/mm_format.pdf
  !!
  !! @author Yohan Lee-tin-yien
  !!
  Subroutine XFTS_dm_mmwrite(DM, funit,info )

    !* Arguments *!
    Type(XFTS_dense_matrix_t), intent(in) :: DM
    Integer             , intent(in)  :: funit
    Integer             , intent(out) :: info 

    !* Local variables *!
    
    ! scalars
    Integer :: m
    Integer :: n
    Integer :: ldA
    FTS_INT :: i
    FTS_INT :: j
    FTS_INT :: ij,ji

    ! arrays
    XFTS_FLOAT, Pointer :: A(:)

    ! string
    Character*10 :: rep
    Character*7  :: field
    Character*19 :: symm

    !- End of header -----------------------------------------------------------

    !-----------------------------------------------------------------
    ! [1] Initialize data / Check
    !-----------------------------------------------------------------

    info  =  0
    m     =  DM%m
    n     =  DM%n
    ldA   =  DM%ld
    A     => DM%v

    ! check sizes in symmetric cases
    If (dm%sym /= DM_SYM_IsGeneral) Then
       CHCKASSRT(m == n, info)
       FTS_ONFAILURE_RETURN(info)
    End If

    ! check the unit
    Inquire(UNIT=funit,WRITE=rep)
    CHCKASSRT( rep(1:1) /= 'N', info )
    FTS_ONFAILURE_RETURN(info)

    !-----------------------------------------------------------------
    ! [2] Write the header
    !-----------------------------------------------------------------

    rep = "array"
#if  XFTS_HAVE_ARITH_S || XFTS_HAVE_ARITH_D
    field = "real"
#elif XFTS_HAVE_ARITH_C || XFTS_HAVE_ARITH_Z
    field = "complex"
#else
#error "bad arithmetic"    
#endif

    Select Case (dm%sym)
    Case(DM_SYM_IsGeneral  ); symm = "general"
    Case(DM_SYM_IsSymmetric); symm = "symmetric"
    Case(DM_SYM_IsSPD      ); symm = "symmetric"
    End Select

    ! types
    Write(funit,'(3(A," "),A)') "%%MatrixMarket matrix",&
         Trim(rep), Trim(field), Trim(symm) 
    Write (funit,'(2(I10))') m, n

    !-----------------------------------------------------------------
    ! [3] Write the data
    !-----------------------------------------------------------------

    Select Case (dm%sym)

    Case(DM_SYM_IsGeneral  )
       !-----------------------------------------------------------------
       ! [3.1] General -> write all values
       !-----------------------------------------------------------------

       Do j = 1,n
          Do i = 1,m
             ij = i+(j-1)*ldA

#if  XFTS_HAVE_ARITH_S || XFTS_HAVE_ARITH_D
             Write(funit,'(1(1PE26.8))') A(ij)
#elif XFTS_HAVE_ARITH_C || XFTS_HAVE_ARITH_Z
             Write(funit,'(2(1PE26.8))') Real(A(ij)), Aimag(A(ij))
#endif

          End Do
       End Do

    Case(DM_SYM_IsSymmetric,DM_SYM_IsSPD)

       !-----------------------------------------------------------------
       ! [3.2] Symmetric -> write only lower triangular part
       !-----------------------------------------------------------------

       If( dm%stored == DM_STORED_UPPER )Then        ! Upper
          Do j = 1,n
             Do i = j,n
                ji = j+(i-1)*ldA
#if  XFTS_HAVE_ARITH_S || XFTS_HAVE_ARITH_D
                Write(funit,'(1(1PE26.8))') A(ji)
#elif XFTS_HAVE_ARITH_C || XFTS_HAVE_ARITH_Z
                Write(funit,'(2(1PE26.8))') Real(A(ji)), Aimag(A(ji))
#endif
             End Do
          End Do

       Else                        ! All or Lower

          Do j = 1,n
             Do i = j,n
                ij = i+(j-1)*ldA
#if  XFTS_HAVE_ARITH_S || XFTS_HAVE_ARITH_D
                Write(funit,'(1(1PE26.8))') A(ij)
#elif XFTS_HAVE_ARITH_C || XFTS_HAVE_ARITH_Z
                Write(funit,'(2(1PE26.8))') Real(A(ij)), Aimag(A(ij))
#endif
             End Do
          End Do

       End If

    End Select

  End Subroutine XFTS_dm_mmwrite

  ! [+] routine : XFTS_dm_mmwritecoo -------------------------------------------
  !
  !> export to matrix market coordinate format 
  !!
  !! write a dense matrix in a file with Matrix Market Format. 
  !! Write the entries in coordinate mode (not in array mode).
  !!
  !! @param[in   ] DM       the matrix 
  !! @param[in   ] funit     the fortran unit stream
  !! @param[in   ] info      the routine status
  !!
  !!-----
  !!
  !! @note 
  !! See the specification of the format at :
  !! http://people.sc.fsu.edu/~jburkardt/pdf/mm_format.pdf
  !!
  !! @author Yohan Lee-tin-yien
  !!
  Subroutine XFTS_dm_mmwritecoo(DM, funit,info )

    !* Arguments *!
    Type(XFTS_dense_matrix_t), intent(in) :: DM
    Integer             , intent(in)  :: funit
    Integer             , intent(out) :: info 

    !* Local variables *!
    
    ! scalars
    Integer :: m
    Integer :: n
    Integer :: ldA
    FTS_INT :: i
    FTS_INT :: j
    FTS_INT :: ij,ji
    FTS_INT :: cntnullval

    ! arrays
    XFTS_FLOAT, Pointer :: A(:)

    ! string
    Character*10 :: rep

    Character*19 :: symm
#if ( XFTS_HAVE_ARITH_S ) || ( XFTS_HAVE_ARITH_D )
    Character(len=FTSOLVER_STRL), Parameter :: MMfmt   = '(2I10,1PE26.18)'
    Character*7, Parameter  :: field = "real"
#elif ( XFTS_HAVE_ARITH_C ) || ( XFTS_HAVE_ARITH_Z )
    Character(len=FTSOLVER_STRL), Parameter :: MMfmt   = "(2I10,1PE26.18,1PE26.18)"
    Character*7, Parameter  :: field = "complex"
#else
#error "Bad arithmetic"
#endif

    !- End of header -----------------------------------------------------------

    !-----------------------------------------------------------------
    ! [1] Initialize data / Check
    !-----------------------------------------------------------------

    info  =  0
    m     =  DM%m
    n     =  DM%n
    ldA   =  DM%ld
    A     => DM%v

    ! check sizes in symmetric cases
    If (dm%sym /= DM_SYM_IsGeneral) Then
       CHCKASSRT(m == n, info)
       FTS_ONFAILURE_RETURN(info)
    End If

    ! check the unit
    Inquire(UNIT=funit,WRITE=rep)
    CHCKASSRT( rep(1:1) /= 'N', info )
    FTS_ONFAILURE_RETURN(info)

    !-----------------------------------------------------------------
    ! [2] Write the header
    !-----------------------------------------------------------------

    rep = "coordinate"

    Select Case (dm%sym)
    Case(DM_SYM_IsGeneral  ); symm = "general"
    Case(DM_SYM_IsSymmetric); symm = "symmetric"
    Case(DM_SYM_IsSPD      ); symm = "symmetric"
    End Select

    ! types
    Write(funit,'(3(A," "),A)') "%%MatrixMarket matrix",&
         Trim(rep), Trim(field), Trim(symm) 



    !-----------------------------------------------------------------
    ! [3] Write the data
    !-----------------------------------------------------------------

    Select Case (dm%sym)

    Case(DM_SYM_IsGeneral  )
       !-----------------------------------------------------------------
       ! [3.1] General -> write all values
       !-----------------------------------------------------------------

       ! sizes
       cntnullval = 0
       Do j = 1,n
          Do i = 1,m
             ij = i+(j-1)*ldA
             If (A(ij) == XFTS_FLOATZERO) &
                  cntnullval=cntnullval+1
          End Do
       End Do
       Write (funit,'(3(I10))') m, n, m*n - cntnullval
       
       ! values
       Do j = 1,n
          Do i = 1,m
             ij = i+(j-1)*ldA

             If (A(ij) == XFTS_FLOATZERO) Cycle
             Write(funit,FMT=MMFmt) i,j, A(ij)

          End Do
       End Do


    Case(DM_SYM_IsSymmetric,DM_SYM_IsSPD)

       !-----------------------------------------------------------------
       ! [3.2] Symmetric -> write only lower triangular part
       !-----------------------------------------------------------------

       If( dm%stored == DM_STORED_UPPER )Then        ! Upper

          cntnullval = 0
          Do j = 1,n
             Do i = j,n
                ji = j+(i-1)*ldA
                If (A(ji) == XFTS_FLOATZERO) &
                     cntnullval=cntnullval+1
             End Do
          End Do
          Write (funit,'(3(I10))') m, n, m*n - cntnullval

          Do j = 1,n
             Do i = j,n

                ji = j+(i-1)*ldA
                If (A(ji) == XFTS_FLOATZERO) Cycle
                Write(funit,FMT=MMFmt) i,j, A(ji)

             End Do
          End Do

       Else                        ! All or Lower

          cntnullval = 0
          Do j = 1,n
             Do i = j,n
                ij = i+(j-1)*ldA
                If (A(ij) == XFTS_FLOATZERO) &
                     cntnullval=cntnullval+1
             End Do
          End Do
          Write (funit,'(3(I10))') m, n, m*n - cntnullval

          Do j = 1,n
             Do i = j,n

                ij = i+(j-1)*ldA
                If (A(ij) == XFTS_FLOATZERO) Cycle
                Write(funit,FMT=MMFmt) i,j, A(ij)

             End Do
          End Do

       End If

    End Select

  End Subroutine XFTS_dm_mmwritecoo

  ! [+] function : XFTS_dm_sizeof --------------------------------------
  !
  !>  Give the memory size in bytes used to stored the dense matrix.
  !!
  !!----
  !!
  !! @param [in] sm   the sparse matrix to check
  !!
  function XFTS_dm_sizeof(dm) RESULT(mem)

    !* Arguments *!
    
    Integer(kind=8)                   :: mem
    Type(XFTS_dense_matrix_t), intent(in) :: dm

    !- End of header------------------------------------------------------------

    mem = (dm%ld * dm%n) * XFTS_FLOATBYTESIZE

  End function XFTS_dm_sizeof


  ! [+] routine :  XFTS_DM_VectorProduct ----------------------------------
  !
  !> Perform Matrix Vector Product : z <-- Mat . x
  !!
  !! Perform Matrix Vector Product : z <-- Mat . x
  !! With as Matrix a Dense matrix general, symmetric, etc.
  !! This is done by calling BLAS routines.
  !!
  !! @param[in    ] DM
  !!
  !!       The matrix in the matrix vector Product.
  !!
  !! @param[in    ] x
  !!
  !!       The vector in the matrix vector Product.
  !!       
  !! @param[in,out] z
  !!
  !!       The result of the matrix vector Product.
  !!       Its values are overwritten.
  !!
  !! @param[   out] info
  !!
  !!       The routine status.
  !!
  !! @author Yohan Lee-tin-yien
  !!
  !!
  Subroutine XFTS_DM_VectorProduct(dm, x, z, info)

    !* Arguments *!

    Type(XFTS_dense_matrix_t), Intent(in   ) :: dm
    Type(XFTS_dense_matrix_t), Intent(in   ) :: x
    Type(XFTS_dense_matrix_t), Intent(inout) :: z
    Integer             , Intent(  out) :: info

    !* Local Variables *!

    ! Strings
    Character*1 :: trans
    Character*1 :: uplo

    !- End of header -------------------------------------------------

    !-----------------------------------------------------------------
    ! [1] Call BLAS according to the MatrixFormat
    !-----------------------------------------------------------------
    
    info = 0
    Select Case (dm%sym)
    Case (DM_SYM_IsGeneral)

       trans = 'N'
       Call XFTS_ARITHgemv  (                  &
            trans,dm%m, dm%n, XFTS_FLOATONE, &
            dm%v(1), dm%ld,                   &
            x%v(1), 1,                        &
            XFTS_FLOATZERO,                  &
            z%v(1), 1                         &
            )

    Case ( DM_SYM_IsSymmetric, DM_SYM_IsSPD )

       If ( dm%stored == DM_STORED_UPPER ) uplo='U'
       If ( dm%stored == DM_STORED_LOWER ) uplo='L'
       If ( dm%stored == DM_STORED_ALL   ) uplo='U'! or 'L' either works

       Call XFTS_ARITHsymv (                & 
            uplo, dm%n, XFTS_FLOATONE,       &
            dm%v(1), dm%ld,                   &
            x%v(1), 1,                        &
            XFTS_FLOATZERO,                  &
            z%v(1), 1                         &
            )

    Case Default

       CHCKASSRT(.False., info )

    End Select

  End Subroutine XFTS_DM_VectorProduct



End Module XFTS_dense_matrix_mod
