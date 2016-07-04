! shared macros in maphys

! error messages
#ifdef MAPHYS_BACKTRACE_ERROR
! CHCKASSRT and ASSRT macros requires :
!   - the usage of "maphys_error_mod" 
!   - the definition of the string FLNAME.
!
! Warning on their usage, since F90 line are truncated to 132 characters.
! - "x" + "res" must take ? characters
! - FTS_ONFAILURE_* "x" must take less than 20 characters 
!   (with 20 spaces at the begining of the line)
!
#define CHCKASSRT( x, res ) Call fts_check( (x) , __LINE__, res, FLNAME );
#define ASSRT( x ) Call fts_assert( (x) , __LINE__, FLNAME );
#define FTS_ONFAILURE_GOTO9999(x) Call fts_check((x>=0),__LINE__,x,FLNAME);If(x<0)Goto 9999
#define FTS_ONFAILURE_ABORT(x) Call fts_assert((x>=0), __LINE__, FLNAME );
#define FTS_ONFAILURE_RETURN(x) Call fts_check((x>=0),__LINE__,x,FLNAME);If(x<0)Return
#else
#define CHCKASSRT( x, res ) If(.NOT.(x)) res = - __LINE__ ;
#define ASSRT( x ) If(.NOT.(x)) Call fts_abort ;
#define FTS_ONFAILURE_GOTO9999(x) If(x<0)Then; x=__LINE__; Goto 9999; End If
#define FTS_ONFAILURE_ABORT(x) If(x<0) Call fts_abort ;
#define FTS_ONFAILURE_RETURN(x) If(x<0) Return;
#endif

! mkl -- num threads controls
#if HAVE_MKL_SET_NUM_THREADS
#define MKL_SET_NUM_THREADS( x ) Call mkl_set_num_threads( (x) )
#else
#define MKL_SET_NUM_THREADS( x ) 
#endif

!
! specifics checks / assertions
!
! a routine info 
#define CHCKRINFO( x ) CHCKASSRT( x >= 0, x )
! an allocation
#define CHCKALLOC( x ) CHCKASSRT( x == 0, x )
! an read or write
#define CHCKIO( x ) CHCKASSRT( x == 0, x )
! an MPI communication
#define ASSRTMPI( x )  ASSRT( x == MPI_SUCCESS )

