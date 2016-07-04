! Warning: XFTS_GENFILE_COMMENT
#include "fts_defs_f.h"
#include "fts_macros_f.h"

!> FTSOLVER sparse matrix module.
!!
!! set of routines to create and manipulate sparse matrices
!!
Module XFTS_sparse_matrix_mod

  !* Module(s) *!
  Use fts_error_mod
  Use XFTS_sparse_matrix_type 

  !* No Implicit typing *!
  Implicit None

  !* Private constants *!
  Character(len=FTSOLVER_STRL), Private, Parameter :: FLNAME = &
       "XFTS_ARITHfts_sparse_matrix_mod.F90"

  !* Access specifiers *!

  ! array manipulation
  Public :: XFTS_icopy
  Public :: XFTS_xcopy
  Public :: XFTS_idup
  Public :: XFTS_xdup
  Public :: XFTS_iset
  Public :: XFTS_xset
  Public :: XFTS_irealloc
  Public :: XFTS_iperm
  Public :: XFTS_xperm
  Public :: XFTS_invperm
  Public :: XFTS_qsort
  Public :: XFTS_sm_ind2ptr
  Public :: XFTS_sm_ptr2ind


  ! creation / deletion
  Public :: XFTS_sm_nullify
  Public :: XFTS_sm_free
  Public :: XFTS_sm_ijv
  Public :: XFTS_sm_create
  Public :: XFTS_sm_createFromData
  Public :: XFTS_sm_dup
  Public :: XFTS_sm_insertEntry 
  Public :: XFTS_sm_realloc

 
  ! IO
  Public :: XFTS_sm_mmread
  Public :: XFTS_sm_ijvread
  Private :: XFTS_sm_ijvReadData
  Public :: XFTS_sm_hbread
  Public :: XFTS_sm_mmwrite

  ! others 
  Public :: XFTS_sm_sizeof
  Public :: XFTS_sm_check
  Public :: XFTS_sm_mirror
  Public :: XFTS_sm_convert
  Public :: XFTS_sm_get_submatrices   



  ! mathematic operations
  Public :: XFTS_sm_transposeUpper
  Public :: XFTS_sm_assemble
  Public :: XFTS_sm_symStruct
  Private :: HasEntry
  Public :: XFTS_sm_vectorproduct 

  ! MPI operations
  Public :: XFTS_sm_bcast
  Public :: XFTS_sm_allgather


  !* Routines *!

Contains

  ! [+] routine : XFTS_icopy -------------------------------------------------
  !
  !> copy an MPI_INT integer array
  !!
  !! @param[in,out] array    the array (size)
  !! @param[in    ] tocopy   the array to copy     (size)
  !! @param[in    ] size   the size of the array
  !!                      
  Subroutine XFTS_icopy(array, tocopy , size) 

    !* Arguments *!

    FTS_INT , pointer , intent(inout) :: array  (:)
    FTS_INT , pointer , intent(in   ) :: tocopy (:)
    FTS_INT           , intent(in   ) :: size

    !* Local variables *!

    FTS_INT           :: k

    ! End of header ----------------------------------------------------------

    Do k = 1, size
       array(k) = tocopy(k)
    End Do

  End Subroutine XFTS_icopy


  ! [+] routine : XFTS_xcopy -------------------------------------------------
  !
  !> copy an XFTS_FLOAT array
  !!
  !! @param[in,out] array    the array             (size)
  !! @param[in    ] tocopy   the array to copy     (size)
  !! @param[in    ] size   the size of the array
  !!                      
  Subroutine XFTS_xcopy(array, tocopy , size) 

    !* Arguments *!

    XFTS_FLOAT , pointer , intent(inout) :: array  (:)
    XFTS_FLOAT , pointer , intent(in   ) :: tocopy (:)
    FTS_INT           , intent(in   ) :: size

    !* Local variables *!

    FTS_INT           :: k

    ! End of header ----------------------------------------------------------

    Do k = 1, size
       array(k) = tocopy(k)
    End Do

  End Subroutine XFTS_xcopy

  ! [+] routine : XFTS_idup -------------------------------------------------
  !
  !> duplicate an FTS_INT integer array
  !!
  !! @param[in,out] array    the array             (size)
  !! @param[in    ] tocopy   the array to copy     (size)
  !! @param[in    ] size     the size of the array
  !! @param[   out] info     the return status
  !!                      
  Subroutine XFTS_idup(array, tocopy , size, info) 

    !* Arguments *!

    FTS_INT , Pointer , Intent(  out) :: array  (:)
    FTS_INT , Pointer , Intent(in   ) :: tocopy (:)
    FTS_INT           , Intent(in   ) :: size
    Integer           , Intent(  out) :: info

    ! End of header ----------------------------------------------------------

    Nullify(array)
    Allocate(array(size),STAT=info)
    CHCKALLOC(info)
    FTS_ONFAILURE_RETURN(info)
    
    Call XFTS_icopy(array,tocopy,size)

  End Subroutine XFTS_idup

  ! [+] routine : XFTS_xdup -------------------------------------------------
  !
  !> duplicate an XFTS_FLOAT integer array
  !!
  !! @param[in,out] array    the array             (size)
  !! @param[in    ] tocopy   the array to copy     (size)
  !! @param[in    ] size     the size of the array
  !! @param[   out] info     the return status
  !!                      
  Subroutine XFTS_xdup(array, tocopy , size, info) 

    !* Arguments *!

    XFTS_FLOAT , Pointer , Intent(  out) :: array  (:)
    XFTS_FLOAT , Pointer , Intent(in   ) :: tocopy  (:)
    FTS_INT           , Intent(in   ) :: size
    Integer           , Intent(  out) :: info

    ! End of header ----------------------------------------------------------

    Nullify(array)
    Allocate(array(size),STAT=info)
    CHCKALLOC(info)
    FTS_ONFAILURE_RETURN(info)
    
    Call XFTS_xcopy(array,tocopy,size)

  End Subroutine XFTS_xdup


  ! [+] routine : XFTS_iset -------------------------------------------------
  !
  !> set a FTS_INT array to a value
  !!
  !! @param[in,out] i         the array             (size)
  !! @param[in    ] istart    the index to start
  !! @param[in    ] iend      the index to end
  !! @param[in    ] ival      the value to set
  !!                      
  Subroutine XFTS_iset(i, istart , iend, ival ) 

    !* Arguments *!

    FTS_INT , Pointer , Intent(inout) :: i     (:)
    FTS_INT           , Intent(in   ) :: istart
    FTS_INT           , Intent(in   ) :: iend
    FTS_INT           , Intent(in   ) :: ival

    !* Local variable
    FTS_INT :: k

    ! End of header ----------------------------------------------------------

    Do k=istart,iend
       i(k) = ival
    End Do
    
  End Subroutine XFTS_iset

  ! [+] routine : XFTS_xset -------------------------------------------------
  !
  !> set a XFTS_FLOAT array to a value
  !!
  !! @param[in,out] x         the array             (size)
  !! @param[in    ] xstart    the index to start
  !! @param[in    ] xend      the index to end
  !! @param[in    ] xval      the value to set
  !!                      
  Subroutine XFTS_xset(x, xstart , xend, xval ) 

    !* Arguments *!

    XFTS_FLOAT , Pointer , Intent(inout) :: x     (:)
    FTS_INT           , Intent(in   ) :: xstart
    FTS_INT           , Intent(in   ) :: xend
    XFTS_FLOAT        , Intent(in   ) :: xval

    !* Local variable
    FTS_INT :: k

    ! End of header ----------------------------------------------------------

    Do k=xstart,xend
       x(k) = xval
    End Do
    
  End Subroutine XFTS_xset


  ! [+] routine : XFTS_irealloc -------------------------------------------------
  !
  !> reallocate a FTS_INT array to a new size
  !!
  !! @param[in,out] array     the array             (size)
  !! @param[in    ] newsize   the new size
  !! @param[in    ] oldsize   the old size
  !! @param[in    ] info      the return status
  !!                      
  Subroutine XFTS_irealloc( array , newsize, oldsize, info ) 

    !* Arguments *!

    FTS_INT , Pointer , Intent(inout) :: array (:)
    FTS_INT           , Intent(in   ) :: newsize
    FTS_INT           , Intent(in   ) :: oldsize
    Integer           , Intent(  out) :: info

    !* Local variable
    FTS_INT, Pointer :: tmp(:)

    ! End of header ----------------------------------------------------------

    !---------------------------------------------------------------------------
    ! [1] Avoid unecessary work
    !---------------------------------------------------------------------------

    info = FTS_SUCCESS
    If ( newsize ==  oldsize ) Return

    !---------------------------------------------------------------------------
    ! [2] Allocate the new array
    !---------------------------------------------------------------------------

    Nullify(tmp)
    Allocate(tmp(newsize), STAT=info)
    CHCKALLOC(info)
    FTS_ONFAILURE_RETURN(info)

    !---------------------------------------------------------------------------
    ! [3] Copy
    !---------------------------------------------------------------------------

    Call XFTS_icopy(tmp,array,Min(oldsize,newsize))

    !---------------------------------------------------------------------------
    ! [4] Fill with null elements
    !---------------------------------------------------------------------------

    If ( newsize > oldsize )Then
       Call XFTS_iset(tmp,oldsize+1,newsize,0)
    End If

    !---------------------------------------------------------------------------
    ! [5] Finish
    !---------------------------------------------------------------------------

    ! replace data
    Deallocate(array)
    array   => tmp
    
  End Subroutine XFTS_irealloc




  ! [+] routine : XFTS_iperm -------------------------------------------------
  !
  !> permute an integer array
  !!
  !! Perform A = A(perm)
  !! where A is a pointer
  !! and  perm a vector permutation.
  !! 
  !! This subroutine use a temporary array,
  !! to avoid issue with huge arrays.     
  !!
  !!-----
  !! 
  !! @param[in,out] array  the array (size)
  !! @param[in    ] perm   the permutation vector (size)
  !! @param[in    ] size   the size of the array
  !! @param[   out] istat  the routine status 
  !!                      
  Subroutine XFTS_iperm(array, perm, size, istat) 

    !* Arguments *!

    FTS_INT , pointer , intent(inout) :: array(:)
    FTS_INT , pointer , intent(in   ) :: perm (:)
    FTS_INT           , intent(in   ) :: size
    Integer         , intent(  out) :: istat

    !* Local variables *!

    FTS_INT , pointer :: new_array (:)
    FTS_INT           :: i

    ! End of header ----------------------------------------------------------

    istat=0 
    Allocate(new_array(size),STAT=iSTAT)
    If (istat/=0) istat=-1

    If (istat==0)Then

       Do i = 1, size
          new_array(i) = array(perm(i))
       End Do

       Do i = 1, size 
          array(i) = new_array(i)
       End Do

       Deallocate(new_array,STAT=istat)
       If (istat/=0) istat = -2

    Endif

  End Subroutine XFTS_iperm

  ! [+] routine : XFTS_xperm -------------------------------------------------
  !
  !> permute a float array
  !!
  !! These subroutines perform A = A(perm)
  !! where A is a pointer
  !! and  perm a vector permuation.
  !! 
  !! They are neccessary for huge arrays.     
  !!
  !!-----
  !! 
  !! @param[in,out] array  the array to permute   (size)
  !! @param[in    ] perm   the permutation vector (size)
  !! @param[in    ] size   the size of the array
  !! @param[   out] istat  the routine status
  !!
  Subroutine XFTS_xperm(array, perm, size, istat) 
    
    !* Arguments *!

    XFTS_FLOAT, Pointer , Intent(inout) :: array(:)
    FTS_INT  , Pointer , Intent(in   ) :: perm (:)
    FTS_INT            , Intent(in   ) :: size
    Integer          , Intent(  out) :: istat

    !* Local variables *!

    XFTS_FLOAT , pointer :: new_array(:)
    FTS_INT             :: i

    !- End of header----------------------------------------------------------

    istat=0 
    Allocate(new_array(size),stat=istat)
    If (istat/=0) istat=-1

    If (istat==0)then
       Do i = 1, size
          new_array(i) = array(perm(i))
       End do

       Do i = 1, size 
          array(i) = new_array(i)
       End Do

       Deallocate(new_array,STAT=istat)
       If (istat/=0) istat = -2

    Endif

  End Subroutine XFTS_xperm

  ! [+] routine : XFTS_invperm ---------------------------------------------èèèè
  !
  !> perform that inverse permutation of an integer array
  !!
  !! These subroutines perform A = perm(A)
  !! where A is a pointer
  !! and  perm a vector permuation.
  !! 
  !! They are neccessary for huge arrays.     
  !!
  !!-----
  !! 
  !! @param[in,out] array  the array to permute   (size)
  !! @param[in    ] perm   the permutation vector (size)
  !! @param[in    ] size   the size of the array
  !! @param[   out] istat  the routine status 
  !!                      
  Subroutine XFTS_invperm(array, perm, size, istat) 

    !* Arguments *!

    FTS_INT , pointer , intent(inout) :: array(:)
    FTS_INT , pointer , intent(in   ) :: perm (:)
    FTS_INT           , intent(in   ) :: size
    Integer         , intent(  out) :: istat

    !* Local variables *!

    FTS_INT , pointer :: new_array (:)
    FTS_INT           :: i

    ! End of header ----------------------------------------------------------

    istat=0 
    Allocate(new_array(size),STAT=iSTAT)
    If (istat/=0) istat=-1

    If (istat==0)Then

       Do i = 1, size
          new_array(i) = perm(array(i))
       End Do

       Do i = 1, size 
          array(i) = new_array(i)
       End Do

       Deallocate(new_array,STAT=istat)
       If (istat/=0) istat = -2

    Endif

  End Subroutine XFTS_invperm



  ! [+] routine : XFTS_qsort ---------------------------------------------------
  !
  !> Sort an array, and get the permutation made.
  !! This routine call the qsort 
  !! from the C standard library to sort the elements.
  !!
  !! @param [in    ] size   The size of "array"
  !! @param [in,out] array  The integer array to be sorted 
  !! @param [in,out] perm   The permutation array allocated by caller
  !! @param [   out] info   The return status
  !!
  !! @author Yohan Lee-tin-yien
  !! @warning biggest matrix tested 
  !!   - mat_scale_100x100.mtx (size = 468027)
  !! 
  Subroutine XFTS_qsort(size, array, perm, info )

    !* Modules *!
    Use fts_error_mod
    Implicit None

    !* Arguments *!
    FTS_INT           ,Intent(in   ) :: size
    FTS_INT  ,Pointer ,Intent(inout) :: array (:)
    FTS_INT  ,Pointer ,Intent(inout) :: perm (:)
    Integer              ,Intent(  out) :: info

    !* Local variables *!
    Integer, Pointer :: arrayAndperm(:) ! concatenation of array and perm.
    Integer          :: i

    !- End of header----------------------------------------------------------
    
    Allocate(arrayAndperm(2*size),STAT=info)
    CHCKASSRT(info == 0, info)
    If (info < 0) Return

    Do i=1,size
       arrayAndperm(2*(i-1)+1) = array(i)
       arrayAndperm(2*i) = i
    End Do

    Call FTS_qsortc( arrayAndperm(1), size )

    Do i=1,size
       array(i) = arrayAndperm(2*(i-1)+1) 
    End Do

    Do i=1,size
       perm(i)  = arrayAndperm(2*i)
    End Do

    Deallocate(arrayAndperm)

  End Subroutine XFTS_qsort


  ! [+] routine : XFTS_sm_ind2ptr ----------------------------------------
  !
  !> Create the pointer list associated to a list of indices.
  !! (for example : row_ind -> row_ptr in the CSR conversion). 
  !!
  !! @param [in ] indsize  The size of "ind"
  !! @param [in ] ind      The list of indices 
  !! @param [in ] ptrsize  The size of "ptr"
  !! @param [out] ptr      The pointer list
  !!
  !! @author Azzam Haidar
  !! @author Yohan Lee-tin-yien
  Subroutine XFTS_sm_ind2ptr( indsize, ind, ptrsize, ptr )

    !* Arguments *!

    FTS_INT, Intent(in ) :: indsize
    FTS_INT, Intent(in ) :: ind(:)
    FTS_INT, Intent(in ) :: ptrsize
    FTS_INT, Intent(out) :: ptr(:)

    !* Local variables *!

    FTS_INT              :: i, ivp

    ! End of header -----------------------------------------------------------

    ! Empty ptr
    Do i=1, ptrsize
       ptr(i) = 0
    End Do

    ! count the number of occurrences of an index.
    Do i=1,indsize
       ivp=ind(i)
       ptr(ivp+1)= ptr(ivp+1) + 1
    End Do

    ! add this nb of indsize by row to find the starting point of each.
    ptr(1)=1
    Do i=2,ptrsize
       ptr(i)= ptr(i) + ptr(i-1)
    End Do

  End Subroutine XFTS_sm_ind2ptr

  ! [+] routine : XFTS_sm_ptr2ind ----------------------------------------
  !
  !> Create the list of indices associated to a pointer list.
  !! (for example : row_ptr -> row_ind  in the CSR conversion). 
  !!
  !! @param [in ] ptrsize  The size of "ptr"
  !! @param [in ] ptr      The pointer list
  !! @param [out] ind      The list of indices (size >= (ptr(ptrsize)-1)) 
  !!
  !! @author Azzam Haidar
  !! @author Yohan Lee-tin-yien
  !!
  Subroutine XFTS_sm_ptr2ind(ptrsize,ptr,ind)

    !* Arguments *!

    FTS_INT, Intent(in ) :: ptrsize
    FTS_INT, Intent(in ) :: ptr(:)
    FTS_INT, Intent(out) :: ind(:)

    !* Local variables *!

    FTS_INT              :: i

    ! End of header -----------------------------------------------------------

    Do i = 1, ptrsize-1   
       ind(ptr(i):ptr(i+1)-1) = i
    End Do

  End Subroutine XFTS_sm_ptr2ind


  ! [+] routine : XFTS_sm_nullify ----------------------------------------
  !
  !> Initialize a sparse matrix
  !!
  !! Initialiaze a sparse matrix, by nullifying all its pointers
  !! and giving defaults values to its attributes.
  !!
  !!-----
  !!
  !! @param[  out] sm    the sparse matrix
  !! @param[  out] istat the routine status
  !!
  !!-----
  !!
  !! @author Yohan Lee-tin-yien
  !!
  Subroutine XFTS_sm_nullify(sm, istat)

    !* Arguments *!

    Type(XFTS_sparse_matrix_t), intent(  out) :: sm
    Integer              , intent(  out) :: istat

    !- End of header----------------------------------------------------------

    ! Set default attributes & sizes
    sm%fmt = SM_FMT_IJV    !< storage format  
    sm%sym = SM_SYM_IsGeneral  !< symmetry        

    sm%m  = 0        !< number of rows
    sm%n  = 0        !< number of columns
    sm%nnz = 0        !< number of non zero

    ! Nullify all arrays
    Nullify(sm%i)     !< rows               [nnz]   
    Nullify(sm%j)     !< columns            [nnz]   

    ! Nullify(sm%cs ) !< compressed vector  [(ni|nj)+1]
    Nullify(sm%csr)   !< compressed vector  [(ni|nj)+1]
    Nullify(sm%csc)   !< compressed vector  [(ni|nj)+1]

    Nullify(sm%v)     !< values             [nnz]  

    ! Finish 
    istat = 0

  End Subroutine XFTS_sm_nullify


  ! [+] routine : XFTS_sm_free -------------------------------------------
  !
  !> Delete a sparse matrix
  !!
  !! Delete a sparse matrix, by nullifying all its pointers
  !! and giving defaults values to its attributes.
  !!
  !!-----
  !!
  !! @param[  out] sm    the sparse matrix
  !! @param[  out] info  the routine status
  !!
  !!-----
  !!
  !! @author Yohan Lee-tin-yien
  !!
  Subroutine XFTS_sm_free(sm,info )

    !* Arguments *!

    Type(XFTS_sparse_matrix_t), intent(inout) :: sm
    Integer              , intent(  out) :: info 

    !- End of header----------------------------------------------------------
    
    !---------------------------------------------------------------------------
    ! [1] free IJV components
    !---------------------------------------------------------------------------

    If( Associated(sm%i) ) Deallocate(sm%i)
    If( Associated(sm%j) ) Deallocate(sm%j)
    If( Associated(sm%v) ) Deallocate(sm%v)

    !---------------------------------------------------------------------------
    ! [2] free CSR/CSC components
    !---------------------------------------------------------------------------

    If ((sm%fmt == SM_FMT_CSR).and.(Associated(sm%csr))) Deallocate(sm%csr)
    If ((sm%fmt == SM_FMT_CSC).and.(Associated(sm%csc))) Deallocate(sm%csc)

    !---------------------------------------------------------------------------
    ! [3] Nullify the sparse matrix
    !---------------------------------------------------------------------------

    Call XFTS_sm_nullify(sm,info)

  End Subroutine XFTS_sm_free


  ! [+] routine : XFTS_sm_ijv --------------------------------------------
  !
  !> create a sparse matrix (in coordinate format) giving its fields
  !! 
  !!
  !! This is the constructor of a sparse matrix in coordinate format,
  !! by giving each required fields.
  !!
  !! Notice that the fields i,j,v of M are pointers
  !! to iA,jA,A sub-arrays (from index 1 to index nnz) 
  !! and not their copy.
  !!
  !!-----
  !!
  !! @param[  out] sm    the sparse matrix
  !! @param[in   ] i     the list of rows
  !! @param[in   ] j     the list of columns
  !! @param[in   ] v     the list of values
  !! @param[in   ] m     the number of rows
  !! @param[in   ] n     the number of columns
  !! @param[in   ] nnz   the number of entries
  !! @param[in   ] sym   the symmetry of the sparse matrix
  !! @param[  out] istat the routine status
  !!
  !!-----
  !!
  !! @author Yohan Lee-tin-yien
  !!
  Subroutine XFTS_sm_ijv &
       (sm, i, j, v, m, n, nnz, sym, istat )

    !* Arguments *!

    Type(XFTS_sparse_matrix_t), Intent(out) :: sm

    Integer, Intent(in)         :: m
    Integer, Intent(in)         :: n
    Integer, Intent(in)         :: nnz
    Integer, Intent(in), Target :: i(:)   !< rows
    Integer, Intent(in), Target :: j(:)   !< columns
    XFTS_FLOAT, Intent(in), Target ::  v(:)   !< values

    Integer, Intent(in)         :: sym ! matrix symmetry
    Integer, Intent(  out)      :: istat

    !- End of header----------------------------------------------------------

    sm%fmt = SM_FMT_IJV

    sm%m = m
    sm%n = n
    sm%nnz = nnz
    sm%i => i
    sm%j => j
    sm%v => v
    sm%sym  = sym

    ! verify valid symmetry
    istat=-1
    if ( sym == SM_SYM_IsSPD ) istat = 0 
    if ( sym == SM_SYM_IsSymmetric ) istat = 0 
    if ( sym == SM_SYM_IsGeneral  ) istat = 0 

    Return
  End Subroutine XFTS_sm_ijv


  ! [+] routine : XFTS_sm_create -----------------------------------------
  !> Create a sparse matrix
  !!
  !! Create a sparse matrix with nnz elements.
  !! Give defaults values to its attributes.
  !! User must fill/modify the data.
  !!
  !!-----
  !!
  !! @param [out]  sm   the sparse matrix 
  !! @param [in ]  nnz  the number of entries in the sparse matrix 
  !! @param [out]  info the routine status
  !!
  !!------
  !!
  !! @author Yohan Lee-tin-yien
  !!
  Subroutine XFTS_sm_create(sm,nnz,info)

    !* Arguments *! 
    Type(XFTS_sparse_matrix_t), Intent(inout) :: sm
    FTS_INT           , Intent(in   ) :: nnz
    Integer              , Intent(  out) :: info

    ! End of header ----------------------------------------------------------

    Nullify(sm%i,sm%j,sm%v)

    ! NB : Give negative "m" and "n" to force user to edit the matrix data.
    sm%nnz = nnz
    sm%m   = -1
    sm%n   = -1

    Allocate( sm%i(nnz), sm%j(nnz), sm%v(nnz), STAT= info )
    CHCKALLOC(info)

  End Subroutine XFTS_sm_create

  ! [+] routine : XFTS_sm_createFromData -----------------------------
  !> Create a sparse matrix and set its data
  !!
  !! Create a sparse matrix by specifying all its data
  !!
  !!-----
  !!
  !! @param [out]  sm   the sparse matrix 
  !! @param [in ]  fmt  the matrix format
  !! @param [in ]  sym  the symmetry
  !! @param [in ]  m    the number of rows 
  !! @param [in ]  n    the number of columns
  !! @param [in ]  nnz  the number of entries in the sparse matrix 
  !!                  it is ignored if if fmt /= SM_FMT_isIJV
  !! @param [in ]  i    pointer to the rows    (or csr if SM_FMT_IsCSR)
  !! @param [in ]  i    pointer to the columns (or csc if SM_FMT_IsCSC)
  !! @param [in ]  v    pointer to the values  (or csc if SM_FMT_IsCSC)
  !! @param [out]  info the routine status
  !!
  !!------
  !!
  !! @author Yohan Lee-tin-yien
  !!
  Subroutine XFTS_sm_createFromData(sm,fmt,sym,&
       m,n,nnz, i,j,v, info)
    Use XFTS_sparse_matrix_type
    Implicit None

    !* Arguments *! 
    Type(XFTS_sparse_matrix_t), Intent(out) :: sm
    Integer             , Intent(in) :: fmt      
    Integer             , Intent(in) :: sym      
    FTS_INT             , Intent(in) :: m        
    FTS_INT             , Intent(in) :: n        
    FTS_INT             , Intent(in) :: nnz      
    FTS_INT, pointer    , Intent(in) :: i   (:)  
    FTS_INT, pointer    , Intent(in) :: j   (:)  
    XFTS_FLOAT, pointer , Intent(in) :: v   (:)  
    Integer              , Intent(  out) :: info
    
    ! End of header ----------------------------------------------------------

    !-------------------------------------------------------------------------
    ! [-] Init
    !-------------------------------------------------------------------------

    Call XFTS_sm_nullify(sm, info )
    FTS_ONFAILURE_RETURN(info)

    Select Case (fmt)
    Case( SM_FMT_IJV, SM_FMT_CSR, SM_FMT_CSC ); info = 0
    Case Default                              ; info = -1 ;
    End Select
    FTS_ONFAILURE_RETURN(info)

    !-------------------------------------------------------------------------
    ! [-] Allocate
    !-------------------------------------------------------------------------
    
    Select Case (fmt)
    Case( SM_FMT_IJV ); sm%nnz = nnz ;
    Case( SM_FMT_CSR ); sm%nnz = i(m+1) -1 ;
    Case( SM_FMT_CSC ); sm%nnz = j(n+1) -1 ;
    End Select
    


    Allocate(sm%i(sm%nnz),STAT=info)
    CHCKALLOC(info)
    FTS_ONFAILURE_RETURN(info);

    Allocate(sm%j(sm%nnz),STAT=info)
    CHCKALLOC(info)
    FTS_ONFAILURE_RETURN(info);

    Allocate(sm%v(sm%nnz),STAT=info)
    CHCKALLOC(info)
    FTS_ONFAILURE_RETURN(info);


    Select Case (fmt)
    Case( SM_FMT_IJV ); continue;
    Case( SM_FMT_CSR ); Allocate(sm%csr(m+1),STAT=info) ;
    Case( SM_FMT_CSC ); Allocate(sm%csc(n+1),STAT=info) ;
    End Select

    !-------------------------------------------------------------------------
    ! [-] Copy the data
    !-------------------------------------------------------------------------

    sm%fmt = fmt
    sm%sym = sym
    sm%m   = m
    sm%n   = n

!    ASSRT(fmt == SM_FMT_IJV )
    Select Case (fmt)
    Case( SM_FMT_IJV ) 
       Call XFTS_icopy(sm%i,i, sm%nnz)
       Call XFTS_icopy(sm%j,j, sm%nnz)
       Call XFTS_xcopy(sm%v,v, sm%nnz)

    Case( SM_FMT_CSR )
       Call XFTS_icopy(sm%csr,i,m+1)
       Call XFTS_sm_ptr2ind(m+1,sm%csr,sm%i)
       Call XFTS_icopy(sm%j,j, sm%nnz)
       Call XFTS_xcopy(sm%v,v, sm%nnz)
    Case( SM_FMT_CSC )
       Call XFTS_icopy(sm%i,i, sm%nnz)
       Call XFTS_icopy(sm%csc,j, n+1)
       Call XFTS_sm_ptr2ind(n+1,sm%csc,sm%j)
       Call XFTS_xcopy(sm%v,v, sm%nnz)
    End Select

  End Subroutine XFTS_sm_createFromData



  ! [+] routine : XFTS_sm_dup --------------------------------------------
  !
  !> Duplicate the sparse matrix smin into smout.
  !!
  !!-----
  !!
  !! @param[in    ] smin   the sparse matrix to be duplicated
  !! @param[   out] smout  the duplicate
  !! @param[   out] info   the return status 
  !!
  Subroutine XFTS_sm_dup(smout, smin, info )
    Implicit None

    !* Arguments *!

    Type(XFTS_sparse_matrix_t), Intent(in ) :: smin
    Type(XFTS_sparse_matrix_t), Intent(out) :: smout
    Integer                   , Intent(out) :: info 

    !- End of header----------------------------------------------------------

    Call XFTS_sm_nullify(smout,info)
    FTS_ONFAILURE_GOTO9999(info)
    
    smout%m   =  smin%m
    smout%n   =  smin%n
    smout%nnz  =  smin%nnz
    smout%fmt  =  smin%fmt
    smout%sym  =  smin%sym
    
    !
    Call XFTS_idup(smout%i,smin%i,smout%nnz,info)
    FTS_ONFAILURE_GOTO9999(info)

    Call XFTS_idup(smout%j,smin%j,smout%nnz,info)
    FTS_ONFAILURE_GOTO9999(info)

    Call XFTS_xdup(smout%v,smin%v,smout%nnz,info)
    FTS_ONFAILURE_GOTO9999(info)

    ! 
    If ( smout%fmt == SM_FMT_CSR ) Then ! CSR
       Call XFTS_idup(smout%csr,smin%csr,smout%m+1,info)
       FTS_ONFAILURE_GOTO9999(info)
    End If

    If ( smout%fmt == SM_FMT_CSC ) Then ! CSC
       Call XFTS_idup(smout%csc,smin%csc,smout%n+1,info)
       FTS_ONFAILURE_GOTO9999(info)
    End If

    ! Exit 
9999 Continue

    If (info<0) Call XFTS_sm_free(smout,info)

  End Subroutine XFTS_sm_dup



  ! [+] routine : XFTS_sm_insertEntry ------------------------------------
  !
  !> insert an entry (i,j,v) at index "index" in the sparse matrix "sm".
  !! Return -1, if "index > sm%nnz", 0 otherwise.
  !!
  !!----
  !!
  !! @param [in,out] sm      The sparse matrix
  !! @param [in    ] index   The index of the entry
  !! @param [in    ] i       The row    of the entry
  !! @param [in    ] j       The column of the entry
  !! @param [in    ] v       The value  of the entry
  !! @param [   out] info    The return status 
  !!
  !!----
  !!
  !! @author Yohan Lee-tin-yien
  !!
  Subroutine XFTS_sm_insertEntry(sm,index, i,j,v, info)

    !* Arguments *!

    Type(XFTS_sparse_matrix_t), intent(inout) :: sm
    FTS_INT           , intent(in   ) :: index
    FTS_INT           , intent(in   ) :: i
    FTS_INT           , intent(in   ) :: j
    XFTS_FLOAT         , intent(in   ) :: v
    Integer              , intent(  out) :: info

    ! End of header ------------------------------------------------------------

    !---------------------------------------------------------------------------
    ! [1] Check
    !---------------------------------------------------------------------------

    If ( index > sm%nnz )Then
       info = -1
       Return
    Else
       info =  0
    End If

    !---------------------------------------------------------------------------
    ! [2] Insert
    !---------------------------------------------------------------------------
    
    sm%i(index) = i
    sm%j(index) = j
    sm%v(index) = v

  End Subroutine XFTS_sm_insertEntry



  ! [+] routine : XFTS_sm_realloc ----------------------------------------
  !
  !> Reallocate a sparse matrice to "newnnz" entries.
  !!
  !! In "newnnz" is inferior than "sm%nnz", only the "newnnz" first elements
  !! of "sm%i,sm%j,sm%k" are kept.
  !!
  !!----
  !!
  !! @param [in,out] sm      The sparse matrix
  !! @param [in    ] newnnz  The new number of entries
  !! @param [   out] info    The return status
  !!
  !!----
  !!
  !! @author Yohan Lee-tin-yien
  !!
  Subroutine XFTS_sm_realloc(sm, newnnz, info)
    
    !* Module(s) *!
    Use fts_log_mod

    !* Arguments *!

    Type(XFTS_sparse_matrix_t), intent(inout) :: sm
    FTS_INT           , intent(in   ) :: newnnz
    Integer              , intent(  out) :: info

    !* Local variables *!

    ! scalars
    FTS_INT :: nnzcopied

    ! arrays
    FTS_INT , Pointer       :: newi(:)
    FTS_INT , Pointer       :: newj(:)
    XFTS_FLOAT , Pointer       :: newv(:)

    ! End of header ------------------------------------------------------------

    !---------------------------------------------------------------------------
    ! [1] Avoid unecessary work
    !---------------------------------------------------------------------------
    
    If ( newnnz == sm%nnz ) Return

    !---------------------------------------------------------------------------
    ! [2] Allocate the new arrays
    !---------------------------------------------------------------------------

    Nullify(newi,newj,newv)

    Allocate(newi(newnnz),newj(newnnz),newv(newnnz),STAT=info)
    CHCKALLOC(info)
    FTS_ONFAILURE_GOTO9999(info)

    !---------------------------------------------------------------------------
    ! [3] Copy
    !---------------------------------------------------------------------------
    
    nnzcopied = Min(sm%nnz,newnnz)
    Call XFTS_icopy(newi,sm%i,nnzcopied)
    Call XFTS_icopy(newj,sm%j,nnzcopied)
    Call XFTS_xcopy(newv,sm%v,nnzcopied)

    !---------------------------------------------------------------------------
    ! [4] Fill with null elements
    !---------------------------------------------------------------------------

    If ( newnnz > sm%nnz )Then
       Call XFTS_iset(newi,sm%nnz+1,newnnz,0)
       Call XFTS_iset(newj,sm%nnz+1,newnnz,0)
       Call XFTS_xset(newv,sm%nnz+1,newnnz,XFTS_FLOATZERO)
    End If

    !---------------------------------------------------------------------------
    ! [5] Finish
    !---------------------------------------------------------------------------

    ! replace data
    Deallocate(sm%i)
    Deallocate(sm%j)
    Deallocate(sm%v)
    
    sm%nnz = newnnz
    sm%i   => newi
    sm%j   => newj
    sm%v   => newv

9999 Continue
    ! handle errors
    If (info<0)Then
       If (Associated(newi)) Deallocate(newi)
       If (Associated(newj)) Deallocate(newj)
       If (Associated(newv)) Deallocate(newv)
    End If
    
  End Subroutine XFTS_sm_Realloc


  ! [+] routine : XFTS_sm_mmread -----------------------------------------
  !
  !> read a sparse matrix from matrix market file. 
  !!
  !! The matrix in the file must follow ftsolver float representation 
  !! that is to say either real or complex.
  !! This Subroutine only support symmetric and general matrixes.
  !! On symmetric matrices, only the lower half entries are given.
  !!
  !!-----
  !!
  !! @param[  out] sm    the read sparse matrix
  !! @param[in   ] iu    the unit file (must be opened/closed by the caller)
  !! @param[  out] info  the return status
  !!
  !!-----
  !!
  !! @par File description
  !!
  !! @verbatim
  !! %%MatrixMarket matrix coordinate (complex|real) (symmetric|general)
  !! % [comments 1]
  !! % [comments 2]
  !! % ...
  !!  m n  nnz
  !!  i j  v
  !!  ...
  !! @endverbatim
  Subroutine XFTS_sm_mmread(sm, iu, info )
    Implicit None

    !* Arguments *!
    
    Type(XFTS_sparse_matrix_t), intent(  out) :: sm
    Integer              , intent(in   ) :: iu
    Integer              , intent(  out) :: info 

    !* Local Variables *!
    Integer :: k
    Character(len=FTSOLVER_STRL) ::  string

    !- End of header----------------------------------------------------------

    !----------------------------------------------------------------------------- 
    ! [1] Init
    !----------------------------------------------------------------------------- 

    Call XFTS_sm_nullify(sm,info)
    FTS_ONFAILURE_RETURN(info)

    sm%fmt = SM_FMT_IJV

    !----------------------------------------------------------------------------- 
    ! [2] Read header
    !----------------------------------------------------------------------------- 

    ! Read 

    Read(iu,'(A)',IOSTAT=info) string
    string = Adjustl(string)
 

    ! Check that matrix is in coordinate format
    k=Index(string,"%%MatrixMarket matrix coordinate")
    CHCKASSRT(k /= 0 ,info)
    FTS_ONFAILURE_RETURN(info)
    
    ! Check the arithmetic
#if ( XFTS_HAVE_ARITH_S ) || ( XFTS_HAVE_ARITH_D )
    k=Index(string,"real")
#elif ( XFTS_HAVE_ARITH_C ) || ( XFTS_HAVE_ARITH_Z )
    k=Index(string,"complex")
#endif
    CHCKASSRT(k /= 0,info)
    FTS_ONFAILURE_RETURN(info)

    ! Set & Check the symmetry

    info = -1
    
    k=Index(string,"symmetric")
    If (k/=0) sm%sym = SM_SYM_IsSymmetric
    If (k/=0) info = FTS_SUCCESS

    k=Index(string,"general")
    If (k/=0) sm%sym = SM_SYM_IsGeneral
    If (k/=0) info = FTS_SUCCESS

    FTS_ONFAILURE_RETURN(info)

    !----------------------------------------------------------------------------- 
    ! [3] Jump comments
    !----------------------------------------------------------------------------- 

    string="%"
    Do While (string(1:1) == '%')
       If(info ==0) Read(iu,*,IOSTAT=info) string
    End Do
    If( info == 0) Backspace(iu,IOSTAT=info )
    CHCKIO(info)
    FTS_ONFAILURE_RETURN(info)

    !----------------------------------------------------------------------------- 
    ! [4] Read the data
    !----------------------------------------------------------------------------- 

    Call XFTS_sm_ijvReadData(sm, iu, info )
    FTS_ONFAILURE_RETURN(info)
    
  End Subroutine XFTS_sm_mmread


  ! [+] routine : XFTS_sm_ijvread ----------------------------------------
  !
  !> read a sparse matrix in coordinate format.
  !!
  !! The matrix in the file must follow ftsolver float representation 
  !! that is to say either real or complex.
  !! This file format do not differenciate general/symmetric matrixes.
  !! On symmetric matrices all the entries are given.
  !!
  !!-----
  !!
  !! @param[  out] sm    the read sparse matrix
  !! @param[in   ] iu    the unit file (must be opened/closed by the caller)
  !! @param[  out] info  the return status
  !!
  !!-----
  !!
  !! @par File description
  !!
  !! @verbatim
  !!  m n nnz
  !!  i j  v
  !!  ...
  !! @endverbatim
  !!
  Subroutine XFTS_sm_ijvread(sm, iu, info)

    !* Arguments *!

    Type(XFTS_sparse_matrix_t), intent(  out) :: sm
    Integer              , intent(in   ) :: iu
    Integer              , intent(  out) :: info

    !* Local Variables *!

    !- End of header----------------------------------------------------------

    !-------------------------------------------------------------------------
    ! [1] Set the attributes
    !-------------------------------------------------------------------------

    sm%fmt = SM_FMT_IJV
    sm%sym = SM_SYM_IsGeneral

    !-------------------------------------------------------------------------
    ! [2] Get the data
    !-------------------------------------------------------------------------

    Call XFTS_sm_ijvReadData(sm,iu,info)

  End Subroutine XFTS_sm_ijvread

  ! [+] routine : XFTS_sm_ijvReadData ------------------------------------
  !
  !> read the data of a sparse matrix from a unit (see the format below).
  !!
  !! Read the components "m,n,nnz, i,j,v" of a the sparse matrix.
  !!
  !!-----
  !!
  !! @param[in,out] sm    the sparse matrix 
  !! @param[in    ] iu    the unit file (must be opened/closed by the caller)
  !! @param[   out] info  the return status
  !!
  !!-----
  !!
  !! @par File description
  !!
  !! @verbatim
  !!  m n nnz
  !!  i j  v
  !!  ...
  !! @endverbatim
  !!
  Subroutine XFTS_sm_ijvReadData(sm, iu, info )

    !* Arguments *!
    
    Type(XFTS_sparse_matrix_t), Intent(inout) :: sm
    Integer              , Intent(in   ) :: iu
    Integer              , Intent(  out) :: info 

    !* Local Variables *!

    Integer      :: k
#if ( XFTS_HAVE_ARITH_C ) || ( XFTS_HAVE_ARITH_Z )
    Real(kind=8) :: Rpart, Ipart 
#endif

    !- End of header----------------------------------------------------------

    info = FTS_SUCCESS
    Nullify(sm%i,sm%j,sm%v)

    !----------------------------------------------------------------------------- 
    ! [1] Read sizes
    !----------------------------------------------------------------------------- 

    Read(iu,*,IOSTAT= info ) sm%m,sm%n,sm%nnz
    CHCKIO(info)
    CHCKASSRT( sm%m   >  0, info )
    CHCKASSRT( sm%n   >  0, info )
    CHCKASSRT( sm%nnz >= 0, info )
    FTS_ONFAILURE_RETURN(info)

    !----------------------------------------------------------------------------- 
    ! [2] Allocate memory
    !----------------------------------------------------------------------------- 

    Allocate(sm%i(sm%nnz),sm%j(sm%nnz),sm%v(sm%nnz),STAT=info )
    CHCKALLOC(info)
    FTS_ONFAILURE_RETURN(info)

    !----------------------------------------------------------------------------- 
    ! [3] Read matrix values
    !----------------------------------------------------------------------------- 

    k=0
    Do While( (k < sm%nnz).And.(info == 0))
       k = k +1

#if ( XFTS_HAVE_ARITH_S ) || ( XFTS_HAVE_ARITH_D )
       Read(iu,*,IOSTAT=info ) sm%i(k),sm%j(k),sm%v(k)
#elif ( XFTS_HAVE_ARITH_C ) || ( XFTS_HAVE_ARITH_Z )
       Read(iu,*,IOSTAT=info ) sm%i(k),sm%j(k),Rpart, Ipart
       sm%v(k)=CMPLX(Rpart,Ipart)
#endif

    End do
    CHCKIO(info)
    FTS_ONFAILURE_RETURN(info)

  End Subroutine XFTS_sm_ijvReadData


  ! [+] routine : XFTS_sm_hbread -----------------------------------------
  !
  !> read a sparse matrix from an harwell boeing file
  !!
  !! The matrix in the file must follow ftsolver float representation 
  !! that is to say either real or complex.
  !! This Subroutine only support symmetric and general matrixes.
  !!
  !!-----
  !!
  !! @param[  out] sm    the read sparse matrix
  !! @param[in   ] iu    the unit file (must be opened/closed by the caller)
  !! @param[  out] info  the return status
  !!
  !!-----
  !!
  !! @par File description
  !!
  !! see <http://math.nist.gov/MatrixMarket/formats.html#hb>
  !!  
  Subroutine XFTS_sm_hbread(sm, iu, info )

    Implicit None

    !* Arguments *!

    Type(XFTS_sparse_matrix_t), intent(  out) :: sm
    Integer              , intent(in   ) :: iu
    Integer              , intent(  out) :: info 

    !* Local Variables *!

    Character(len=72)  ::  title
    Character(len=8)   ::  key
    Character(len=3)   ::  mxtype
    Character(len=16)  ::  ptrfmt, indfmt
    Character(len=20)  ::  valfmt, rhsfmt

    Integer :: totcrd, ptrcrd, indcrd, valcrd, rhscrd
    Integer :: nrow  , ncol  , nnzero, neltvl

    Integer      :: k

    !- End of header----------------------------------------------------------

    Call XFTS_sm_nullify(sm,info)
    FTS_ONFAILURE_RETURN(info)

    !---------------------------------------------------------------------------
    ! [1] read & check header
    !---------------------------------------------------------------------------

    ! read
    Read(iu,'( A72, A8 / 5I14 / A3, 11X, 4I14 / 2A16, 2A20 )',IOSTAT=info) &
         title, key, &
         & totcrd, ptrcrd, indcrd, valcrd, rhscrd, &
         & mxtype, nrow  , ncol  , nnzero, neltvl, &
         & ptrfmt, indfmt, valfmt, rhsfmt
    CHCKIO(info)
    FTS_ONFAILURE_RETURN(info)

    ! check

    ! matrix arithmetic

#if ( XFTS_HAVE_ARITH_S ) || ( XFTS_HAVE_ARITH_D )
    if ((mxtype(1:1) /= 'R').and.(mxtype(1:1) /= 'r')) info = -1
#elif XFTS_HAVE_ARITH_C || XFTS_HAVE_ARITH_Z
    if ((mxtype(1:1) /= 'C').and.(mxtype(1:1) /= 'c')) info = -1
#endif
    FTS_ONFAILURE_RETURN(info)

    ! matrix symetry

    sm%sym = -1
    If (mxtype(2:2) == 'S') sm%sym = SM_SYM_IsSymmetric
    If (mxtype(2:2) == 's') sm%sym = SM_SYM_IsSymmetric
    If (mxtype(2:2) == 'U') sm%sym = SM_SYM_IsGeneral
    If (mxtype(2:2) == 'u') sm%sym = SM_SYM_IsGeneral
    If (sm%sym == -1 ) info = -1
    FTS_ONFAILURE_RETURN(info)

    ! assembled 

    If ((mxtype(3:3) /= 'A').and.(mxtype(3:3) /= 'a')) info = -1
    FTS_ONFAILURE_RETURN(info)

    ! set the attributes

    sm%fmt =  SM_FMT_CSC
    sm%m  = nrow
    sm%n  = ncol
    sm%nnz = nnzero

    !---------------------------------------------------------------------------
    ! [3] Allocate memory
    !---------------------------------------------------------------------------

    Allocate(sm%csc(sm%n+1),STAT=info)
    CHCKALLOC(info)
    FTS_ONFAILURE_GOTO9999(info)

    Allocate(sm%i(sm%nnz),sm%j(sm%nnz),sm%v(sm%nnz),STAT=info)
    CHCKALLOC(info)
    FTS_ONFAILURE_GOTO9999(info)

    !---------------------------------------------------------------------------
    ! [4] Read matrix structure
    !---------------------------------------------------------------------------

    Read(iu,PTRFMT,IOSTAT=info) (sm%csc(k),k =1,ncol+1)
    CHCKIO(info)
    FTS_ONFAILURE_GOTO9999(info)

    Read(iu,INDFMT,IOSTAT=info) (sm%i(k)  ,k =1,nnzero )
    CHCKIO(info)
    FTS_ONFAILURE_GOTO9999(info)

    !---------------------------------------------------------------------------
    ! [5] Read matrix values
    !---------------------------------------------------------------------------

    If (VALCRD .GT. 0)Then

       ! Read (iu,VALFMT,iostat=info) (sm%v(k), k=1,nnzero)
       Read (iu,VALFMT,IOSTAT=info) sm%v(1:nnzero)
       CHCKIO(info)
       FTS_ONFAILURE_GOTO9999(info)

    Endif 

    !---------------------------------------------------------------------------
    ! [6] Expand sm%csc into sm%j
    !---------------------------------------------------------------------------

    Call XFTS_sm_ptr2ind(sm%n+1, sm%csc, sm%j)

    !---------------------------------------------------------------------------
    ! [7] Finish
    !---------------------------------------------------------------------------

9999 Continue
    If (info < 0) Call XFTS_sm_free(sm,info)
    
  End Subroutine XFTS_sm_hbread


  ! [+] routine : XFTS_sm_mmwrite ----------------------------------------
  !
  !> write a sparse matrix in a file with Matrix Market Format. 
  !!
  !!-----
  !!
  !! @param[in   ] sm    the sparse matrix to write
  !! @param[in   ] funit the unit file (must be opened/closed by the caller)
  !! @param[  out] info  the return status
  !!
  !!-----
  !!
  !! @see XFTS_sm_mmread()
  !!
  Subroutine XFTS_sm_mmwrite(sm,funit,info)

    !* Module(s) *!
    Use fts_log_mod

    !* Arguments *!

    Type(XFTS_sparse_matrix_t), intent(in) :: sm
    integer, intent(in)               :: funit
    integer, intent(out)              :: info

    !* Local Variables *!

    ! Constants
    Character(len=FTSOLVER_STRL), Parameter :: MMbegin = &
         "%%MatrixMarket matrix coordinate"
    Character(len=FTSOLVER_STRL), Parameter :: st_err = &
         "Error: in XFTS_sm_mmwrite, "

#if ( XFTS_HAVE_ARITH_S ) || ( XFTS_HAVE_ARITH_D )
    Character(len=FTSOLVER_STRL), Parameter :: MMfmt   = '(2I10,1PE26.18)'
    Character(len=FTSOLVER_STRL), Parameter :: MMarith = " real"
#elif ( XFTS_HAVE_ARITH_C ) || ( XFTS_HAVE_ARITH_Z )
    Character(len=FTSOLVER_STRL), Parameter :: MMfmt   = "(2I10,1PE26.18,1PE26.18)"
    Character(len=FTSOLVER_STRL), Parameter :: MMarith = " complex"
#endif

    ! Scalars 
    FTS_INT :: k

    ! Strings
    Character(len=FTSOLVER_STRL) :: MMsym 

    !- End of header----------------------------------------------------------

    !-------------------------------------------------------------------------
    ! [1] Write header
    !-------------------------------------------------------------------------
    Select Case( sm%sym )
    Case(SM_SYM_IsGeneral  ); MMsym = " general  "
    Case(SM_SYM_IsSymmetric); MMsym = " symmetric"
    Case(SM_SYM_IsSPD      ); MMsym = " symmetric"
    End Select

    Write (funit,'(3A)', IOSTAT=info) Trim(MMbegin),Trim(MMarith), Trim(MMsym)
    CHCKIO(info)
    FTS_ONFAILURE_RETURN(info)

    !-------------------------------------------------------------------------
    ! [2] Write sizes
    !-------------------------------------------------------------------------

    Write (funit, '(3I10)',IOSTAT=info ) sm%m, sm%n, sm%nnz
    CHCKIO(info)
    FTS_ONFAILURE_RETURN(info)

    !-------------------------------------------------------------------------
    ! [3] Write entries
    !-------------------------------------------------------------------------

    Write(funit, FMT=MMfmt, IOSTAT=info ) &
         (sm%i(k),sm%j(k),sm%v(k),k=1, sm%nnz)
    CHCKIO(info)
    FTS_ONFAILURE_RETURN(info)

  End Subroutine XFTS_sm_mmwrite

  ! [+] function : XFTS_sm_sizeof --------------------------------------
  !
  !>  Give the memory size in bytes used to stored the sparse matrix.
  !!
  !!----
  !!
  !! @param [in] sm   the sparse matrix to check
  !!
  function XFTS_sm_sizeof(sm) RESULT(mem)

    !* Arguments *!
    
    Integer(kind=8)                   :: mem
    Type(XFTS_sparse_matrix_t), intent(in) :: sm
    

    !- End of header------------------------------------------------------------

    mem = 0
    mem = mem + sm%nnz*XFTS_FLOATBYTESIZE     ! v
    mem = mem + sm%nnz*FTS_INTBYTESIZE  *2 ! i,j

    If ( sm%fmt == SM_FMT_CSC ) mem = mem + (sm%n+1)*(FTS_INTBYTESIZE)
    If ( sm%fmt == SM_FMT_CSR ) mem = mem + (sm%m+1)*(FTS_INTBYTESIZE)
    
    mem = mem + 2*4                  ! flags (fmt,sym)
    mem = mem + 3*FTS_INTBYTESIZE ! sizes (m,n,nnz)

  End function XFTS_sm_sizeof

  ! [+] routine : XFTS_sm_check ------------------------------------------
  !
  !>  check a sparse matrix.
  !!
  !! check for inconsistencies of sparse matrix fields (array dimensions)
  !!
  !!----
  !!
  !! @param [in] sm   the sparse matrix to check
  !!
  Subroutine XFTS_sm_check(sm)

    !* Arguments *!

    Type(XFTS_sparse_matrix_t), intent(in) :: sm

    !* Local Variables *!

    Integer :: ierr(9)

    !- End of header----------------------------------------------------------

    ierr(:) = 0

    ! errors
    if ( sm%nnz  <= 0 ) ierr(7) = 1
    if ( sm%m   <= 0 ) ierr(8) = 1
    if ( sm%n   <= 0 ) ierr(9) = 1

    if ( size(sm%i) < sm%nnz ) ierr(2) = 2
    if ( size(sm%j) < sm%nnz ) ierr(3) = 2
    if ( size(sm%v) < sm%nnz ) ierr(4) = 2

    select case(sm%fmt)     ! check the compressed array
    case( SM_FMT_CSR )
       if ( (.not. associated(sm%csr)         ) .or. &
            (size(sm%csr) /= (sm%m +1)      ) .or. &
            ( sm%nnz /=  (sm%csr(sm%m+1)-1) )) ierr(5)=1
    case( SM_FMT_CSC )
       if ( (.not. associated(sm%csc)        )  .or. &
            (size(sm%csc) /= (sm%n +1)      ) .or. &
            ( sm%nnz /=  (sm%csc(sm%n+1)-1) )) ierr(6)=1
    End select

    ! warnings
    if ( size(sm%i) /= sm%nnz ) ierr(2) = 1
    if ( size(sm%j) /= sm%nnz ) ierr(3) = 1
    if ( size(sm%v) /= sm%nnz ) ierr(4) = 1


    ! global
    ierr(1) = maxval(ierr(2:6))

    ! inform user
    if(ierr(1) == 0 ) write(6,*) "SUCCESS : matrix checked"

    if(ierr(2) == 1 ) write(6,*) "WARNING : size(sm%i) is ", size(sm%i),"instead of", sm%nnz
    if(ierr(3) == 1 ) write(6,*) "WARNING : size(sm%j) is ", size(sm%j),"instead of", sm%nnz
    if(ierr(4) == 1 ) write(6,*) "WARNING : size(sm%v) is ", size(sm%v),"instead of", sm%nnz

    if(ierr(2) == 2 ) write(6,*) "ERROR : size(sm%i) is ", size(sm%i),"instead of", sm%nnz
    if(ierr(3) == 2 ) write(6,*) "ERROR : size(sm%j) is ", size(sm%j),"instead of", sm%nnz
    if(ierr(4) == 2 ) write(6,*) "ERROR : size(sm%v) is ", size(sm%v),"instead of", sm%nnz

    if(ierr(5) /= 0 ) write(6,*) "ERROR : mismatch beetween sm%fmt and sm%csr"
    if(ierr(6) /= 0 ) write(6,*) "ERROR : mismatch beetween sm%fmt and sm%csc"

    if(ierr(7) /= 0 ) write(6,*) "ERROR : sm%nnz  is strictly positive", sm%nnz
    if(ierr(8) /= 0 ) write(6,*) "ERROR : sm%m   is strictly positive", sm%m
    if(ierr(9) /= 0 ) write(6,*) "ERROR : sm%n   is strictly positive", sm%n

  End Subroutine XFTS_sm_check

  ! [+] routine : XFTS_sm_mirror -----------------------------------------
  !
  !> Mirror the sparse matrix entries according to the diagonal.
  !!
  !! For each non-diagonal entry (i,j,v) , add a new entry (j,i,conjugate(v))
  !! 
  !! @warning small extra memory consumption.
  !!  This routine may allocate a more memory than necessary.
  !!  Namely, it allocates "2*nnz" elements instead of counting them.
  !!
  !!----
  !!
  !! @param[in,out] sm   the sparse matrix 
  !!    - On input,  the sparse matrix to be symmetrised according to its diagonal.
  !!    - On output, the sparse matrix symmetry is General
  !!                 with the new entries appended.
  !! @param[  ,out] info the routine status
  !!
  !!----
  !!
  !! @par History
  !!
  !! @verbatim
  !!   Date     : Version : Comments
  !! - 26/01/11 :         : handle matrix with non diagonal element.
  !! - 19/11/10 : v0.1    : Create routine (Yohan Lee-tin-yien)
  !! @endverbatim
  !!
  Subroutine XFTS_sm_mirror(sm,info)

    !* Module(s) *!
    Use fts_log_mod
    Implicit None

    !* Arguments *!
    Type(XFTS_sparse_matrix_t) , Intent(inout) :: sm
    Integer               , Intent(  out) :: info

    !* Local variables *!

    ! scalar
    FTS_INT :: i_order ! prefixes : i_ = input, o_ = output
    FTS_INT :: i_nnz   
    FTS_INT :: o_nnz
    FTS_INT :: i       ! dummy integer on inputs
    FTS_INT :: o       ! dummy integer on outputs

    ! arrays
    FTS_INT  , Pointer :: i_row (:)
    FTS_INT  , Pointer :: i_col (:)
    XFTS_FLOAT, Pointer :: i_val (:)

    FTS_INT  , Pointer :: o_row (:)
    FTS_INT  , Pointer :: o_col (:)
    XFTS_FLOAT, Pointer :: o_val (:)


    !- End of header ---------------------------------------------------------

    !-------------------------------------------------------------------------
    ! [0] Check the matrix 
    !-------------------------------------------------------------------------
    info = FTS_SUCCESS

    ! matrix must be square
    CHCKASSRT( sm%m == sm%n, info)
    FTS_ONFAILURE_RETURN(info)

    !-------------------------------------------------------------------------
    ! [1] Init
    !-------------------------------------------------------------------------

    Nullify(i_row, i_col, i_val )
    Nullify(o_row, o_col, o_val )

    ! get sizes
    i_order = sm%m
    i_nnz   = sm%nnz

    ! get the output's number of entries 
    ! we estimate it here by its upper limit.
    o_nnz = 2 * i_nnz 

    !
    i_row => sm%i
    i_col => sm%j
    i_val => sm%v

    Allocate(o_row(o_nnz),o_col(o_nnz),o_val(o_nnz),STAT=info)
    CHCKALLOC(info)
    FTS_ONFAILURE_RETURN(info)

    !-------------------------------------------------------------------------
    ! [2] Constructs o_row, o_col, o_val
    !-------------------------------------------------------------------------

    ! copy previous values
    Call XFTS_icopy(o_row,i_row,i_nnz)
    Call XFTS_icopy(o_col,i_col,i_nnz)
    Call XFTS_xcopy(o_val,i_val,i_nnz)

    ! for each entry (i,j) (i/=j) add entry (j,i)
    o = i_nnz
    Do i = 1, i_nnz
       ! jump the diagonal
       If( i_row(i) == i_col(i) ) Cycle

       ! add entry
       o = o + 1
       o_row(o) = i_col(i) 
       o_col(o) = i_row(i) 

#if ( XFTS_HAVE_ARITH_S ) || ( XFTS_HAVE_ARITH_D )
       o_val(o) = i_val(i) 
#elif ( XFTS_HAVE_ARITH_C ) || ( XFTS_HAVE_ARITH_Z )
       o_val(o) = CONJG( i_val(i) )
#endif

    End Do

    !-------------------------------------------------------------------------
    ! [3] Finish
    !-------------------------------------------------------------------------

    ! Update the structure

    If ( sm%fmt == SM_FMT_CSR ) Deallocate( sm%csr )
    If ( sm%fmt == SM_FMT_CSC ) Deallocate( sm%csc )

    Deallocate(sm%i)
    Deallocate(sm%j)
    Deallocate(sm%v)

    sm%fmt =  SM_FMT_IJV
    sm%sym =  SM_SYM_IsGeneral
    sm%nnz =  o 
    sm%i   => o_row
    sm%j   => o_col
    sm%v   => o_val

  End subroutine XFTS_sm_mirror



  ! [+] routine : XFTS_sm_convert ----------------------------------------
  !
  !>  modify the format of a sparse matrix.
  !!
  !! convert the format from IJV,CSR,CSC to IJV,CSR,CSC
  !!
  !!-----
  !!
  !! @param[in,out] sm    the sparse matrix to be updated
  !! @param[in    ] fmt   the wanted format (IJV,CSR,CSC)
  !! @param[   out] info  the routine status
  !!
  !!-----
  !!
  !! @author Yohan Lee-tin-yien
  !!
  Subroutine XFTS_sm_convert( sm, fmt, info )

    !* Modules *!
    Use fts_error_mod

    !* Arguments *!

    Type(XFTS_sparse_matrix_t), intent(inout) :: sm
    Integer              , intent(in)    :: fmt
    Integer              , intent(out)   :: info

    !* Local Variables *!

    ! scalars
    FTS_INT :: cssize

    ! related to IJV to CSR or CSC conversion
    ! Integer, Parameter :: inc_flag=2;
        
    Integer, Pointer :: perm (:)
    Integer, Pointer :: cs (:)

    Integer, Pointer :: master (:)
    Integer, Pointer :: slave (:)

    !- End of header----------------------------------------------------------

    !-------------------------------------------------------------------------
    ! [1] CSR or CSC to IJV 
    !-------------------------------------------------------------------------

    info = 0
    Nullify(cs)
    Nullify(perm)
    Nullify(master)
    Nullify(slave)

    !-------------------------------------------------------------------------
    ! [1.1] Avoid unecessary conversion
    !-------------------------------------------------------------------------

    If( fmt == sm%fmt ) Return

    !-------------------------------------------------------------------------
    ! [1.2] Convert from CSR to IJV
    !-------------------------------------------------------------------------

    csr2ijv: if ( sm%fmt == SM_FMT_CSR )then
       Call XFTS_sm_ptr2ind(sm%m+1,sm%csr,sm%i)
       If ( Associated(sm%csr) ) Deallocate(sm%csr)
       sm%fmt = SM_FMT_IJV
    End if csr2ijv

    !-------------------------------------------------------------------------
    ! [1.3] Convert from CSC to IJV
    !-------------------------------------------------------------------------

    csc2ijv: if ( sm%fmt == SM_FMT_CSC )then
       Call XFTS_sm_ptr2ind(sm%n+1,sm%csc,sm%j)
       If ( Associated(sm%csc) ) Deallocate(sm%csc)
       sm%fmt = SM_FMT_IJV
    End if csc2ijv

    !-------------------------------------------------------------------------
    ! [2] Convert from IJV to CSR or CSC
    !-------------------------------------------------------------------------

    ijv2csrc:if ( fmt /= SM_FMT_IJV )then
       ! init

       Nullify(master,slave,perm,cs)

       Select Case (fmt)
       Case ( SM_FMT_CSR )
          cssize = sm%m+1
          master => sm%i
          slave  => sm%j
       Case ( SM_FMT_CSC )
          cssize = sm%n+1
          master => sm%j
          slave  => sm%i
       End Select

       ! sort ascendently

       ! allocate(perm(sm%nnz) )   ; perm(:) = 0
       ! Call KB07AI_mod(sm%i, sm%nnz, perm)  !* safer not-free alternative
       ! Call IPSORT( master(1:sm%nnz), sm%nnz, perm, inc_flag, info)

       Allocate(perm(sm%nnz), STAT= info )
       CHCKALLOC(info)
       FTS_ONFAILURE_RETURN(info)

       Call XFTS_qsort(sm%nnz, master, perm, info)
       FTS_ONFAILURE_RETURN(info)

       Call XFTS_iperm(slave,perm,sm%nnz,info)
       FTS_ONFAILURE_RETURN(info)
       
       Call XFTS_xperm(sm%v  ,perm,sm%nnz,info)
       FTS_ONFAILURE_RETURN(info)

       Deallocate(perm)

       ! construct the compressed array
       Allocate(cs(cssize), STAT= info )
       CHCKALLOC(info)
       FTS_ONFAILURE_RETURN(info)

       Call XFTS_sm_ind2ptr(sm%nnz, master, cssize, cs)

       If( fmt == SM_FMT_CSR) sm%csr => cs
       If( fmt == SM_FMT_CSC) sm%csc => cs

       ! finish
       Nullify(slave)
       Nullify(master)
       Nullify(cs)

       sm%fmt = fmt
    Endif ijv2csrc

  End Subroutine XFTS_sm_convert

  ! [+] routine : XFTS_sm_get_submatrices --------------------------------
  !
  !> Compute the 4 blocs of a square sparse matrix A.
  !!
  !! Computes the blocs Aii, Aib, Abi and Abb of the sparse matrix sm.
  !! The routine do not save the memory on symmetric matrices (Aib = Abi^T).
  !! Matrix sm is in that form : 
  !! @verbatim
  !!              [Aii Aib]     
  !!              [Abi Abb]     
  !! @endverbatim
  !!
  !!----
  !!
  !! @param[in     ] Aii_ndof      the number of rows/columns in the Aii bloc
  !! @param[in     ] A             the local matrix which contains the blocs
  !! @param[out    ] Aii           the bloc ii 
  !! @param[out    ] Aib           the bloc ib 
  !! @param[out    ] Abi           the bloc bi 
  !! @param[out    ] Abb           the bloc bb 
  !! @param[out    ] info          the routine status
  !!
  !!----
  !!
  !! @author Luc   Giraud
  !! @author Azzam Haidar
  !! @author Yohan Lee-tin-yien
  !!
  !! @version 0.1
  Subroutine XFTS_sm_get_submatrices( & ! intents
       Aii_ndof, sm,                         & ! in
       Aii, Aib,                            & ! out
       Abi, Abb,                            &
       info                                 &
       )

    Implicit None

    !* Arguments *!

    Integer              , Intent(in   ) :: Aii_ndof
    Type(XFTS_sparse_matrix_t), Intent(in   ) :: sm
    Type(XFTS_sparse_matrix_t), Intent(  out) :: Aii
    Type(XFTS_sparse_matrix_t), Intent(  out) :: Aib
    Type(XFTS_sparse_matrix_t), Intent(  out) :: Abi
    Type(XFTS_sparse_matrix_t), Intent(  out) :: Abb
    Integer              , Intent(  out) :: info

    !* Local variables *!

    ! scalars
    Integer :: info_ignore
    Integer    :: sym      !< symmetry of a matrix

    FTS_INT :: A_ndof   !< number of rows/columns in A
    FTS_INT :: Abb_ndof !< number of rows/columns in Abb 

    FTS_INT :: Aii_nnz
    FTS_INT :: Aib_nnz 
    FTS_INT :: Abi_nnz
    FTS_INT :: Abb_nnz

    FTS_INT :: k,l,m,n,p !< dummy counters

    FTS_INT   :: i       !< the row    index of an entry
    FTS_INT   :: j       !< the column index of an entry
    XFTS_FLOAT :: v       !< the value        of an entry

    !- End of header -------------------------------------------------------------

    !-----------------------------------------------------------------------------
    ! [1] Check arugments, Initialize local variables, etc.
    !-----------------------------------------------------------------------------

    ! init pointers & structures
    Call XFTS_sm_nullify(Aii,info_ignore)
    Call XFTS_sm_nullify(Abb,info_ignore)
    Call XFTS_sm_nullify(Abi,info_ignore)
    Call XFTS_sm_nullify(Aib,info_ignore)

    ! init values

    info    = FTS_SUCCESS
    sym      = sm%sym
    A_ndof   = sm%m
    Abb_ndof = A_ndof - Aii_ndof

    ! check
    CHCKASSRT(A_ndof > Aii_ndof, info )
    FTS_ONFAILURE_GOTO9999(info)

    !-----------------------------------------------------------------------------
    ! [2] Create the matrices
    !-----------------------------------------------------------------------------

    !-----------------------------------------------------------------------------
    ! [2.1] Identify the differents nnz
    !-----------------------------------------------------------------------------

    ! Reset counters
    Aii_nnz = 0
    Aib_nnz = 0            
    Abi_nnz = 0
    Abb_nnz = 0

    l  = Aii_ndof ! limit beetween the blocks

    ! Increment counters
    Do k=1,sm%nnz

       i=sm%i(k)
       j=sm%j(k)

       If ((i >  l ).and.(j >  l)) Abb_nnz = Abb_nnz + 1
       If ((i <= l ).and.(j <= l)) Aii_nnz = Aii_nnz + 1
       If ((i >  l ).and.(j <= l)) Abi_nnz = Abi_nnz + 1
       If ((i <= l ).and.(j >  l)) Aib_nnz = Aib_nnz + 1

    End Do

    ! Handle symmetric case :  Abi = Aib^T
    If ( sym /= SM_SYM_IsGeneral )Then 
       Aib_nnz = Abi_nnz + Aib_nnz
       Abi_nnz = Aib_nnz
    End If

    ! check
    If ( sym == SM_SYM_IsGeneral ) k= Aii_nnz + Aib_nnz + Abi_nnz + Abb_nnz
    If ( sym /= SM_SYM_IsGeneral ) k= Aii_nnz + Aib_nnz + Abb_nnz
    CHCKASSRT( sm%nnz == k, info ) 
    FTS_ONFAILURE_GOTO9999(info)

    !-----------------------------------------------------------------------------
    ! [2.1] Create the sparse matrices
    !-----------------------------------------------------------------------------

    ! allocate memory
    Call XFTS_sm_create(Aii,Aii_nnz,info)
    FTS_ONFAILURE_GOTO9999(info)

    Call XFTS_sm_create(Abb,Abb_nnz,info)
    FTS_ONFAILURE_GOTO9999(info)

    Call XFTS_sm_create(Abi,Abi_nnz,info)
    FTS_ONFAILURE_GOTO9999(info)

    Call XFTS_sm_create(Aib,Aib_nnz,info)
    FTS_ONFAILURE_GOTO9999(info)

    ! set attributes & sizes

    Aii%m   = Aii_ndof
    Aii%n   = Aii_ndof
    Aii%sym = sm%sym

    Abb%m   = Abb_ndof
    Abb%n   = Abb_ndof
    Abb%sym = sm%sym

    Aib%m   = Aii_ndof
    Aib%n   = Abb_ndof
    Aib%sym = SM_SYM_IsGeneral

    Abi%m   = Abb_ndof
    Abi%n   = Aii_ndof
    Abi%sym = SM_SYM_IsGeneral

    !---------------------------------------------------------------------------
    ! [3] Fill the sparse matrices data
    !---------------------------------------------------------------------------

    ! Reset verification counters
    k = 0
    l = 0
    m = 0
    n = 0

    Select Case ( sym )

    Case ( SM_SYM_IsGeneral ) 

       !------------------------------------------------------------------------
       ! [3.1] UNSYMMETRIC matrices  
       !------------------------------------------------------------------------

       ! Get the values
       do p=1,sm%nnz

          i = sm%i(p)
          j = sm%j(p)
          v = sm%v(p)

          If (i <= Aii_ndof) Then  ! [Aii Aib]
             If (j >  Aii_ndof) Then
                ! the coef is related to Aib
                k         = k + 1
                Aib%v(k)  = v
                Aib%i(k)  = i
                Aib%j(k)  = j
             Else
                ! Coefficient is related to Aii 
                l         = l + 1
                Aii%v(l)  = v
                Aii%i(l)  = i
                Aii%j(l)  = j
             Endif
          Else ! [Abi Abb]
             If (j <= Aii_ndof) Then
                ! the coef is related to Abi
                m         = m + 1
                Abi%v(m)  = v
                Abi%i(m)  = i
                Abi%j(m)  = j
             Else
                ! the coef is related to Abb
                n         = n + 1
                Abb%v(n)  = v
                Abb%i(n)  = i   
                Abb%j(n)  = j   
             Endif
          Endif
       Enddo

       ! Verify
       CHCKASSRT(k == Aib_nnz, info)
       CHCKASSRT(l == Aii_nnz, info)
       CHCKASSRT(m == Abi_nnz, info)
       CHCKASSRT(n == Abb_nnz, info)
       FTS_ONFAILURE_GOTO9999(info)

    Case ( SM_SYM_IsSymmetric, SM_SYM_IsSPD ) 

       !------------------------------------------------------------------------
       ! [3.2] SYMMETRIC matrices  
       !------------------------------------------------------------------------

       ! Get the values
       do p=1,sm%nnz

          i = sm%i(p)
          j = sm%j(p)
          v = sm%v(p)

          If (i <= Aii_ndof) Then  ! [Aii Aib]
             If (j >  Aii_ndof) Then
                k         = k + 1
                ! the coef is related to Aib & Abi
                Aib%v(k)  = v
                Aib%i(k)  = i
                Aib%j(k)  = j

                Abi%v(k)  = v
                Abi%i(k)  = j
                Abi%j(k)  = i
             Else
                l         = l + 1
                ! Coefficient is related to Aii 
                Aii%v(l)  = v
                Aii%i(l)  = i
                Aii%j(l)  = j
             Endif
          Else ! [Abi Abb]
             If (j <= Aii_ndof) Then
                k         = k + 1
                ! the coef is related to Abi & Aib
                Abi%v(k)  = v
                Abi%i(k)  = i
                Abi%j(k)  = j

                Aib%v(k)  = v
                Aib%i(k)  = j
                Aib%j(k)  = i
             Else
                n         = n + 1
                ! the coef is related to Abb
                Abb%v(n)  = v
                Abb%i(n)  = i   
                Abb%j(n)  = j   
             Endif
          Endif
       Enddo

       ! Verify

       CHCKASSRT(k == Aib_nnz, info)
       CHCKASSRT(l == Aii_nnz, info)
       CHCKASSRT(k == Abi_nnz, info)
       CHCKASSRT(n == Abb_nnz, info)
       FTS_ONFAILURE_GOTO9999(info)


    End Select

    !-----------------------------------------------------------------------------
    ! [4] Finish
    !-----------------------------------------------------------------------------
    ! 

9999 Continue

    ! On error, free the allocated memory
    If (info < 0) Then
       Call XFTS_sm_free(Aii,info_ignore)
       Call XFTS_sm_free(Abb,info_ignore)
       Call XFTS_sm_free(Aib,info_ignore)
       Call XFTS_sm_free(Abi,info_ignore)
    End If

  End Subroutine XFTS_sm_get_submatrices



  ! [+] routine :  XFTS_SM_VectorProduct ---------------------------------
  !
  !> Perform Matrix Vector Product : y <-- Mat . x
  !!
  !! Perform Matrix Vector Product : y <-- Mat . x
  !! With as Matrix a sparse matrix general, symmetric, etc.,
  !! Or a section of a sparse matrix.
  !!
  !!----
  !!
  !! @param[in    ] row_displ
  !!
  !!       displacement to apply on the row indexes of "Mat".
  !!       This means that inside the routine, the row "row" on the matrix is :
  !!                    row = sm%row + row_displ
  !!
  !!       If "Mat" is the whole matrix (not a section), 
  !!       row_displ should be equal to "0".
  !!
  !!       If "Mat" is a section of a sparse Matrix, 
  !!       row_displ should be equal to -row_start;
  !!       where row_start the starting index of the section of Mat.
  !!
  !!       If during computation, "row <= 0" the entry is ignored.
  !!
  !! @param[in    ] col_displ
  !!
  !!       displacement to apply on the columns indexes of Mat.
  !!       Same role as "row_displ" but on columns.
  !!
  !! @param[in    ] sm
  !!
  !!       The matrix or its section in the matrix vector Product.
  !!       
  !!
  !! @param[in    ] x
  !!
  !!       The vector in the matrix vector Product.
  !!       
  !! @param[in,out] y
  !!
  !!       The result of the matrix vector Product.
  !!       Its values are overwritten.
  !!
  !! @param[   out] info
  !!
  !!       The routine status.
  !!
  !!----
  !!
  !! @author Azzam Haidar
  !! @author Luc   Giraud
  !! @author Yohan Lee-tin-yien
  !!
  !!
  Subroutine XFTS_SM_VectorProduct &
       (sm, row_displ, col_displ, x, y, info)

    !* Module(s) *!
    Use XFTS_dense_matrix_mod, Only : &
         XFTS_dense_matrix_t
    Use fts_log_mod

    !* Arguments *!

    Type(XFTS_sparse_matrix_t), Intent(in   ) :: sm
    FTS_INT          , Intent(in   ) :: row_displ
    FTS_INT          , Intent(in   ) :: col_displ
    Type(XFTS_dense_matrix_t), Intent(in   ) :: x
    Type(XFTS_dense_matrix_t), Intent(inout) :: y
    Integer             , Intent(  out) :: info

    !* Local Variables *!
    ! Scalar 
    Integer    :: sym
    FTS_INT :: nnz
    FTS_INT :: n
    FTS_INT :: k, row, col

    ! Arrays
    FTS_INT  , Pointer :: Mi (:)
    FTS_INT  , Pointer :: Mj (:)
    XFTS_FLOAT, Pointer :: Mv (:)
    XFTS_FLOAT, Pointer :: Xv (:)
    XFTS_FLOAT, Pointer :: Yv (:)

    !- End of header -------------------------------------------------
    
    !-----------------------------------------------------------------
    ! [1] Initialize variable
    !-----------------------------------------------------------------

    ! Exit early if bad symmetry
    Select Case (SM%sym)
    Case ( SM_SYM_IsGeneral, SM_SYM_IsSPD, SM_SYM_IsSymmetric ) 
       Continue
    Case Default
       CHCKASSRT(.False., info)
       FTS_ONFAILURE_RETURN(info)
    End Select


    sym  =  SM%sym
    nnz  =  sm%nnz 
    Mv   => sm%v
    Mi   => sm%i
    Mj   => sm%j

    n    =  sm%m
    Yv   => y%v
    Xv   => x%v

    !-----------------------------------------------------------------
    ! [2] Perform the matrix product
    !-----------------------------------------------------------------

    ! delete previous value
    Call XFTS_xset(Yv,1,n,XFTS_FLOATZERO)

    ! Do the matrix product
    Select Case (sym)
    Case ( SM_SYM_IsGeneral ) 
       Do k=1,nnz
          
          row = Mi(k)+row_displ
          If (row <= 0) Cycle
          
          col = Mj(k)+col_displ
          If (col <= 0) Cycle
          
          Yv(row) = Yv(row) + Mv(k) * Xv(col)
       End do
       
    Case ( SM_SYM_IsSPD, SM_SYM_IsSymmetric ) 
       
       Do k=1,nnz
          
          row = Mi(k)+row_displ
          col = Mj(k)+col_displ
          
          If (row <= 0) Cycle
          If (col <= 0) Cycle
          
          If(row /= col)Then 
             Yv(row) = Yv(row) + Mv(k)*Xv(col)
             Yv(col) = Yv(col) + Mv(k)*Xv(row)
          Else
             Yv(row) = Yv(row) + Mv(k)*Xv(col)
          Endif
       End do
       
    End Select

  End Subroutine XFTS_SM_VectorProduct

  ! [+] routine : XFTS_sm_bcast ------------------------------------------
  !
  !> broadcast a sparse matrix
  !!
  !! broadcast the sparse matrix "sm"
  !! from the MPI process "root" to all processors in "comm".
  !!
  !!----
  !!
  !! @param[in,out] sm    the sparse matrix to broadcast 
  !! @param[in    ] root  the rank according to "comm" where sm is.
  !! @param[in,out] comm  the MPI Communicator
  !! @param[   out] info  the routine status
  !!
  !!----
  !!
  !! @author Yohan Lee-tin-yien
  !!
  Subroutine XFTS_sm_bcast(sm, root, comm, info  )

    !* Module(s) & co. *!
    Use fts_log_mod
    Use fts_error_mod
    Implicit None
    Include 'mpif.h'

    !* Arguments *!

    Type(XFTS_sparse_matrix_t) , Intent(inout) :: sm
    Integer               , Intent(in)    :: root
    Integer               , Intent(inout) :: comm
    Integer               , Intent(out)   :: info 

    !* Local Variables *!

    ! scalars
    Integer    :: rank
    FTS_INT :: attributes(5)


    !- End of header----------------------------------------------------------

    !-------------------------------------------------------------------------
    ! [1] Get mpi rank
    !-------------------------------------------------------------------------

    Call MPI_comm_rank(comm, rank, info)
    ASSRTMPI(info)

    !-------------------------------------------------------------------------
    ! [2] Broadcast the attributes (sizes, format, etc.)
    !-------------------------------------------------------------------------

    If (rank == root)Then
       attributes(1) = sm%m
       attributes(2) = sm%n
       attributes(3) = sm%nnz
       attributes(4) = sm%fmt
       attributes(5) = sm%sym
    End If

    Call MPI_Bcast(attributes, 5,FTS_INTMPI,root,comm,info)
    ASSRTMPI(info)

    !-------------------------------------------------------------------------
    ! [3] non root, allocate the memory & set the attributes
    !-------------------------------------------------------------------------

    If (rank /= root)Then

       Call XFTS_sm_Create(sm, attributes(3),info)
       FTS_ONFAILURE_ABORT(info)
       
       sm%m    = attributes(1)
       sm%n    = attributes(2)
       sm%nnz  = attributes(3)
       sm%fmt  = attributes(4)
       sm%sym  = attributes(5)

       Select Case(sm%fmt)
       Case(SM_FMT_CSR)
          Allocate(sm%csr(sm%m+1),STAT=info)
          FTS_ONFAILURE_ABORT(info)
       Case(SM_FMT_CSC)
          Allocate(sm%csc(sm%n+1),STAT=info)
          FTS_ONFAILURE_ABORT(info)
       Case Default
          Continue
       End Select

    End If

    !-------------------------------------------------------------------------
    ! [4] Broadcast the data
    !-------------------------------------------------------------------------

    !

    Call MPI_Bcast (sm%i, sm%nnz, FTS_INTMPI  , root, comm, info)
    ASSRTMPI(info)

    Call MPI_Bcast (sm%j, sm%nnz, FTS_INTMPI  , root, comm, info)
    ASSRTMPI(info)
    
    Call MPI_Bcast (sm%v, sm%nnz, XFTS_FLOATMPI, root, comm, info)
    ASSRTMPI(info)

    !
    Select Case(sm%fmt)
    Case(SM_FMT_CSR)
       Call MPI_Bcast(sm%csr, sm%m+1, FTS_INTMPI,root, comm, info)
       ASSRTMPI(info)

    Case(SM_FMT_CSC)
       Call MPI_Bcast(sm%csc, sm%n+1, FTS_INTMPI,root, comm, info)
       ASSRTMPI(info)

    Case Default
       Continue
    End Select

  End Subroutine XFTS_sm_bcast

  ! [+] routine : XFTS_sm_Allgather --------------------------------------
  !
  !> Gather a sparse matrice from a group of processes.
  !! all the input matrices must have the same symmetry.
  !! The output matrice is in coordinate format.
  !! 
  !!---
  !!
  !! @param[in    ] sm_send The sparse matrix values to send (with the same %sym)
  !! @param[in,out] sm_recv The received sparse matrix values 
  !! @param[in,out] comm    The MPI Communicator of the group of processes.
  !! @param[   out] info    The subroutines status (0 = success, else = error)
  !!
  !!---
  !!
  !! @warning on output, sm_recv is in coordinate format 
  !!
  Subroutine XFTS_sm_Allgather &
       ( sm_send, sm_recv, global_n , comm, info )

    !* Module(s) *!
    Use fts_log_mod
    Use fts_error_mod
    implicit none

    !* Arguments *!

    include 'mpif.h'
    Type(XFTS_sparse_matrix_t), intent(in   ):: sm_send
    Type(XFTS_sparse_matrix_t), intent(inout):: sm_recv
    Integer              , intent(in   ):: comm
    FTS_INT           , intent(in   ):: global_n
    Integer              , intent(  out):: info

    !* Local Variables *!

    ! scalars
    Integer :: np
    FTS_INT :: i

    ! arrays
    FTS_INT, pointer :: nnz_all(:), displ(:)

    !- End of header----------------------------------------------------------

    !-------------------------------------------------------------------------
    ! [1] Init
    !------------------------------------------------------------------------- 

    ! get number of processes
    Call MPI_Comm_size(comm, np, info) 
    ASSRTMPI(info)

    ! allocate local memory
    Nullify(nnz_all,displ)
    Allocate(nnz_all(np),displ(np), STAT= info)
    CHCKALLOC(info)
    FTS_ONFAILURE_ABORT(info)

    ! Check the symmetry
    Call MPI_Allgather( &
         sm_send%sym  , 1,MPI_INTEGER,   &
         nnz_all, 1,MPI_INTEGER,         &
         comm, info   )
    ASSRTMPI(info)

    info = FTS_SUCCESS
    If ( Any( nnz_all(1:np) /= sm_send%sym) ) info = -1
    FTS_ONFAILURE_ABORT(info)

    !-------------------------------------------------------------------------
    ! [2] Create sm_recv structure
    !------------------------------------------------------------------------- 

    ! gather nnz
    nnz_all(1:np) = 0
    Call MPI_AllGather(                &
         sm_send%nnz  , 1,MPI_INTEGER, &
         nnz_all , 1,MPI_INTEGER,      &
         comm, info)
    ASSRTMPI(info)

    ! allocate memory
    sm_recv%nnz = Sum( nnz_all(1:np) )
    Call XFTS_sm_create(sm_recv, sm_recv%nnz, info )
    FTS_ONFAILURE_ABORT(info)

    ! set attributes & sizes
    sm_recv%m =  global_n
    sm_recv%n = global_n
    sm_recv%fmt = SM_FMT_IJV
    sm_recv%sym = sm_send%sym

    !-------------------------------------------------------------------------
    ! [3] Gather the data
    !------------------------------------------------------------------------- 

    ! construct displacement vector for gatherv
    displ = 0
    do i=2,np
       displ(i)= displ(i-1) + nnz_all(i-1)
    End do

    ! gather the fields 
    Call MPI_Allgatherv                            &
         (sm_send%i     , sm_send%nnz,MPI_INTEGER, &
         sm_recv%i, nnz_all, displ, MPI_INTEGER,   &
         comm, info)
    ASSRTMPI(info)

    Call MPI_Allgatherv                            &
         (sm_send%j     , sm_send%nnz,MPI_INTEGER, &
         sm_recv%j, nnz_all, displ, MPI_INTEGER,   &
         comm, info)
    ASSRTMPI(info)

    Call MPI_Allgatherv                            &
         (sm_send%v     , sm_send%nnz,XFTS_FLOATMPI, &
         sm_recv%v, nnz_all, displ, XFTS_FLOATMPI,   &
         comm, info)
    ASSRTMPI(info)

    ! free temporary memory
    If (Associated(nnz_all)) Deallocate(nnz_all)
    If (Associated(displ  )) Deallocate(displ )

  End Subroutine XFTS_sm_Allgather

  ! [+] routine : XFTS_sm_transposeUpper ---------------------------------
  !
  !> Transpose the entries in the sparse matrix upper triangular part.
  !!
  !! @param[in,out ] M             the sparse matrix
  !!
  Subroutine XFTS_sm_transposeUpper(sm)
    
    !* Arguments *!

    Type(XFTS_sparse_matrix_t), intent(inout) :: sm

    !* Local variables *!

    ! scalars
    FTS_INT :: ind
    FTS_INT :: k

    !- End of header -----------------------------------------------------------

    Do k=1,sm%nnz
       If ( sm%i(k) < sm%j(k) )Then 
          ind     = sm%i(k)
          sm%i(k) = sm%j(k)
          sm%j(k) = ind
       End If
    End Do

  End Subroutine XFTS_sm_transposeUpper

  ! [+] routine : XFTS_sm_assemble ---------------------------------------
  !
  !> assemble a sparse matrix (sum duplicate entries)
  !!
  !! @param[in,out ] M             the matrix to assemble
  !! @param[in,out ] global_n      the order of the matrix
  !! @param[out    ] info          the Subroutine status
  !!
  !!----
  !!
  !! @Warning:
  !!
  !!  this Subroutine have stability issue for
  !!  very large duppicated entries on very large problems. 
  !!  This is due to IPSORT function stability issue.
  !! 
  !!  A partial solution is to use HSL sort functions : kb07ai_mod.
  !!  But the HSL Subroutine was removed due to copyrigth issues.
  !!  
  !!  Bug found by Azzam Haidar (30/04/2009)
  !!
  !!
  !! @author Azzam Haidar
  !! @author Yohan Lee-tin-yien
  !!
  !! @version 0.1
  !!
  Subroutine XFTS_sm_assemble( sm , global_n, info )

    !* Module(s) *!
    Use fts_log_mod
    Use fts_error_mod
    Implicit None
    Include "mpif.h"

    !* Arguments *!

    Type(XFTS_sparse_matrix_t), intent(inout) :: sm
    FTS_INT           , intent(in   ) :: global_n
    Integer              , intent(  out) :: info

    !* Local variables *!

    ! scalars
    FTS_INT :: i
    FTS_INT :: k
    FTS_INT :: start
    FTS_INT :: oldnnz
    FTS_INT :: newk
    FTS_INT :: st,ed

    ! arrays
    FTS_INT, Pointer :: perm (:), jptr(:)
    XFTS_FLOAT, Pointer :: vptr(:)

    !- End of header -----------------------------------------------------------

    !---------------------------------------------------------------------------
    ! [0] Initialize local variables
    !---------------------------------------------------------------------------

    Nullify(perm)
    oldnnz   = sm%nnz

    !---------------------------------------------------------------------------
    ! [2] Sort the entries on M
    !---------------------------------------------------------------------------
    !>  Warning: on IPSORT (slatec)
    !!
    !!  the XFTS_sm_convert subroutine have stability issue for
    !!  very large duppicated entries on very large problems.
    !!  This is due to IPSORT function stability issue.
    !!
    !!  A partial solution is to use HSL sort functions : kb07ai_mod.
    !!  But the HSL subroutine was removed due to copyright issues.
    !!
    !!  Bug found by Azzam Haidar (30/04/2009)
    !!

    Call XFTS_sm_convert(sm, SM_FMT_CSR, info )
    FTS_ONFAILURE_GOTO9999(info)

    ! sort the row indices
    Allocate(perm(oldnnz),STAT=info)
    CHCKALLOC(info)
    FTS_ONFAILURE_GOTO9999(info)

    ! sort each columns of each lines.
    start = 0
    Do i=1,global_n

       ! k = number of element on row i.
       k= sm%csr(i+1) - sm%csr(i)

       jptr => sm%j(start+1:start+k)
       vptr => sm%v(start+1:start+k)

       ! Call IPSORT      (jptr,k,perm,kflag,info)
       Call XFTS_qsort(k, jptr, perm, info )
       FTS_ONFAILURE_GOTO9999(info)

       Call XFTS_xperm(vptr, perm, k, info)
       FTS_ONFAILURE_GOTO9999(info)

       start = start + k 

    End Do

    !---------------------------------------------------------------------------
    ! [2] Sum  duplicate entries
    !---------------------------------------------------------------------------
    !
    ! Algorithm explained on an example.
    !
    ! initial:                  ! before sum1:  
    ! csr= [ 1 2 3 6 8 10]      ! csr= [ 1 2 3 6 8 10]
    !
    ! i= [ 1 2 3 3 3 4 4 5 5 ]  ! i= [ 1 2 3 3 3 4 4 5 5 ]  
    ! j= [ 1 2 4 4 4 4 5 6 6 ]  ! j= [ 1 2 4 4 4 4 5 6 6 ]   
    ! v= [ 0 0 1 1 1 0 0 2 2 ]  ! v= [ 0 0 1 1 1 0 0 2 2 ] 
    !      ^   . . .     . .    !        ^ ^ . . .   . .
    !      st                   !       st ed
    !     newk      ed unset    !     newk
    !                           
    ! sum1:                     ! next iteration:           
    ! csr= [ 1 2 3 6 8 10]      ! csr= [ 1 2 3 4 6 8]
    !
    ! i= [ 1 2 3 3 3 4 4 5 5 ]  ! i= [ 1 2 3 x x 4 4 5 5 ]  
    ! j= [ 1 2 4 4 4 4 5 6 6 ]  ! j= [ 1 2 4 x x 4 5 6 6 ]  
    ! v= [ 0 0 1 1 1 0 0 2 2 ]  ! v= [ 0 0 3 x x 0 0 2 2 ]  
    !          ^     ^   . .    !            ^   ^ ^ . .   
    !          st    ed         !         newk st ed     
    !        newk               !
    !                             
    ! sum2:                     ! end:                 
    ! csr= [ 1 2 3 4 6 8]       ! csr= [ 1 2 3 4 6 7]
    !
    ! i= [ 1 2 3 4 4 x x 5 5 ]  ! i= [ 1 2 3 4 4 5 x x x ]
    ! j= [ 1 2 4 4 5 x x 6 6 ]  ! j= [ 1 2 4 5 4 6 x x x ]
    ! v= [ 0 0 3 0 0 x x 2 2 ]  ! v= [ 0 0 3 0 0 4 x x x ]
    !                ^   ^   ^  !                  ^     ^
    !             newk  st   ed !               newk     st
    !                                                    ed
    st = 1
    newk = st
    Do While (st<=oldnnz)

       ! 
       ed = st + 1

       ! sum identical entries 
       If (ed <= oldnnz)Then
          Do While(  (sm%i(st) == sm%i(ed))&
               .And. (sm%j(st) == sm%j(ed)))
             sm%v(st) = sm%v(st) + sm%v(ed)
             ed = ed + 1
          End Do
       End If

       ! put entry "st" into entry "newk"
       If (newk /= st)Then

          sm%i(newk) = sm%i(st)
          sm%j(newk) = sm%j(st)
          sm%v(newk) = sm%v(st)

       End If

       ! update counters
       newk = newk + 1
       st = ed

    End Do
    sm%nnz = newk - 1

    ! update csr
    Call XFTS_sm_ind2ptr(sm%nnz, sm%i, sm%m+1, sm%csr )
    
    
    !-------------------------------------------------------------------------
    ! [5] Finish
    !------------------------------------------------------------------------- 

9999 Continue
    ! Free memory
    If (Associated(perm)) Deallocate(perm)
    
  End Subroutine XFTS_sm_assemble

  ! [+] routine : XFTS_sm_symStruct --------------------------------------
  !
  !> symmetrise a sparse matrix structure.
  !!
  !! for each non-diagonal entries, add a null symmetric entry if symmetric entry 
  !! do not exist.
  !!
  !!----
  !!
  !! @param [in,out] sm     the sparse matrix to be symmetrized
  !! @param [in,out] info   the routine status
  !!
  !!----
  !!
  !! @author Yohan Lee-tin-yien
  !! 
  Subroutine XFTS_sm_symStruct(sm, info)

    !* Module(s) *!
    Use fts_log_mod
    Implicit None

    !* Arguments *!

    Type(XFTS_sparse_matrix_t), intent(inout) :: sm
    Integer              , intent(  out) :: info

    !* Local variables *!

    ! scalars
    FTS_INT :: cnt
    FTS_INT :: nnz
    FTS_INT :: i
    FTS_INT :: j
    FTS_INT :: k
    FTS_INT :: smaddsize
    
    ! type
    Type(XFTS_sparse_matrix_t) :: smadd
    
    ! End of header ------------------------------------------------------------

    !---------------------------------------------------------------------------
    ! [1] Init
    !---------------------------------------------------------------------------

    ! init smadd
    smaddsize = INT(sm%n/4,kind=FTS_INTKIND)
    Call XFTS_sm_nullify(smadd, info )
    FTS_ONFAILURE_RETURN(info)
    Call XFTS_sm_create(smadd, smaddsize, info)
    FTS_ONFAILURE_RETURN(info)
    
    ! sort and append to a compressed vector if necessary
    If (sm%fmt == SM_FMT_IJV )Then
       Call XFTS_sm_convert(sm, SM_FMT_CSC, info)
       FTS_ONFAILURE_GOTO9999(info)
    End If

    !---------------------------------------------------------------------------
    ! [2] Compute smadd (list of elements to add)
    !---------------------------------------------------------------------------

    cnt = 0
    Do k=1,sm%nnz

       i = sm%i(k)
       j = sm%j(k)

       ! Avoid diagonal
       If ( i == j ) Cycle

       ! Do not add the symmetric entry if it already exist
       If ( HasEntry(sm,j,i) ) Cycle 

       ! add the symmetric entry
       cnt = cnt+1

       If (cnt > smadd%nnz )Then ! smadd too small
          smaddsize = INT(smadd%nnz*1.5,KIND=FTS_INTKIND )+1
          Call XFTS_sm_realloc(smadd, smaddsize,info)
          FTS_ONFAILURE_GOTO9999(info)
       End If

       Call XFTS_sm_insertEntry(smadd,cnt,j,i,XFTS_FLOATZERO,info)
       FTS_ONFAILURE_GOTO9999(info)

    End Do

    !---------------------------------------------------------------------------
    ! [3] Fuse sm and smadd 
    !---------------------------------------------------------------------------

    nnz = sm%nnz
    Call XFTS_sm_realloc(sm, nnz + cnt,info)
    FTS_ONFAILURE_GOTO9999(info)

    Do k=1,cnt
       sm%i(nnz+k) = smadd%i(k)
    End Do
    Do k=1,cnt
       sm%j(nnz+k) = smadd%j(k)
    End Do
    Do k=1,cnt
       sm%v(nnz+k) = smadd%v(k)
    End Do

    If ( cnt > 0 )Then ! appending data broke the csr/csc format.
       Call XFTS_sm_convert(sm,SM_FMT_IJV,info)
       FTS_ONFAILURE_GOTO9999(info)
    End If

    !---------------------------------------------------------------------------
    ! [4] Finish
    !---------------------------------------------------------------------------

9999 Continue
    Call XFTS_sm_free(smadd,info)


  End Subroutine XFTS_sm_symStruct

  ! [+] function : HasEntry ----------------------------------------------------
  !
  !> check if the sm has the entry (i,j).
  !!
  !! the matrix must be in CSC or CSR format.
  !!
  !!----
  !!
  !! @param [in    ] sm      The sparse matrix
  !! @param [in    ] i       The row    of the entry
  !! @param [in    ] j       The column of the entry
  !!
  !!-----
  !!
  !! @author Yohan Lee-tin-yien
  !! 
  Logical Function HasEntry(sm,i,j)
    
    !* Arguments *!
    
    Type(XFTS_sparse_matrix_t), Intent(in) :: sm
    FTS_INT, Intent(in) :: i
    FTS_INT, Intent(in) :: j

    !* Local variables 
    FTS_INT :: k
    FTS_INT :: row
    FTS_INT :: col

    ! End of header ------------------------------------------------------------

    !---------------------------------------------------------------------------
    ! [1] Init
    !---------------------------------------------------------------------------

    HasEntry = .False.


    If ( sm%fmt == SM_FMT_CSC )Then

       !------------------------------------------------------------------------
       ! [2] Search in CSC
       !------------------------------------------------------------------------

       Do k=sm%csc(j),sm%csc(j+1)-1
          row=sm%i(k)
          If (row == i)Then
             HasEntry = .True.
             Exit
          End If
       End Do

    Else If( sm%fmt == SM_FMT_CSR )Then

       !------------------------------------------------------------------------
       ! [3] Search in CSR
       !------------------------------------------------------------------------

       Do k=sm%csr(i),sm%csr(i+1)-1
          col=sm%j(k)
          If (col == j)Then
             HasEntry = .True.
             Exit
          End If
       End Do

    End If

  End Function HasEntry

End Module XFTS_sparse_matrix_mod


