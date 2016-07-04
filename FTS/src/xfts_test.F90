!This program reads input files
! and extracts configuration informations 
! for a given domain/partition

#include"fts_defs_f.h"
      
Program Test

      Use XFTS_ftsolver_mod
      Use XFTS_sparse_matrix_mod
      Use XFTS_toolkit_mod
      Use XFTS_analyze_mod
      Use XFTS_compute_aux_matrices_mod
      Use XFTS_precond_mod
      Use XFTS_Solve_mod
      Use XFTS_state_print_mod
      Use XFTS_state_update_mod
      Implicit None 
      Include 'mpif.h'

      Integer        :: istep
      Integer        :: iinfo
      Integer        :: rank
      Integer        :: unit
      Integer        :: sym     ! symmetry of the input matrix
      Integer        :: job
      Integer        :: mpi_provided
      Integer        :: nbdom
      FTS_INT        :: i
      Real(kind=XFTS_FLOATKIND) :: norm2,m,tmp

      XFTS_FLOAT, Pointer :: sol(:)
      !Strings to read input file
      Character(len=FTSOLVER_STRL) :: paramfile
      Character(len=FTSOLVER_STRL)   :: prefix
      Character(len=FTSOLVER_STRL) :: matrixfile
      Character(len=FTSOLVER_STRL) :: rhsfile
      Character(len=FTSOLVER_STRL) :: outsolfile
      Character(len=FTSOLVER_STRL) :: outrhsfile
      Character(len=FTSOLVER_STRL) :: initguessfile
      Character(len=FTSOLVER_STRL) :: nbdom_str


      !Derived types
      Type (XFTS_ftsolver_t) :: ftsl 
      Type(XFTS_sparse_matrix_t) :: smatrix

        !- End of header -------------------------------------------------------------

      !-----------------------------------------------------------------------------
      ! [1] Initialize parallelism
      !-----------------------------------------------------------------------------
  
      Call MPI_Init_thread(MPI_THREAD_MULTIPLE,mpi_provided,iinfo) !
      istep = 11
      if (iinfo /= MPI_SUCCESS) iinfo = -1
      if (iinfo < 0) goto 9999


      Call MPI_Comm_rank(MPI_COMM_WORLD, rank, iinfo)
      istep = 12
      if (iinfo /= MPI_SUCCESS) iinfo = -1
      if (iinfo < 0) goto 9999

      Nullify(sol)
      ftsl%comm = MPI_COMM_WORLD

      !initialize the instance of ftsolver
      Call XFTS_ftsolver_init(ftsl)

      !Set the number of domain on all mpi domain
      call MPI_REDUCE(rank,nbdom,1,MPI_INTEGER,MPI_MAX,0,MPI_COMM_WORLD,iinfo)


      If ( rank == 0 ) then
            ! Read argument 
         Write(nbdom_str,*) nbdom+1
         nbdom_str=Adjustl(nbdom_str)
         Call Get_command_argument(1, paramfile )
         Call Get_command_argument(2, prefix    )

            !Open parameter file and read matrix file' name
         unit = 11
         Open(UNIT=unit, FILE=paramfile, STATUS='old', ACTION='read')

         Call XFTS_read_param_freeformat &
              (unit,ftsl%icntl,ftsl%rcntl,sym,job, &
              matrixfile,rhsfile,initguessfile,&
              outrhsfile,outsolfile)
         Close(unit)
         
         Write(*,*) 'arith    is :', Trim('XFTS_FLOAT')
         Write(*,*) 'prefix   is :', Trim(prefix)
         Write(*,*) 'matrix   is :', Trim(matrixfile)
         Write(*,*) 'rhs      is :', Trim(rhsfile)
         Write(*,*) 'initguess is :', Trim(initguessfile)
         Write(*,*) 'outsolfile is :', Trim(outsolfile)
         Write(*,*) 'nb of MPI processes :', Trim(nbdom_str)


         !--------------------------------------------------------------------------
         ! [3.2.2] Read topology file
         !--------------------------------------------------------------------------
         
         Select Case( mpi_provided )
         Case( MPI_THREAD_SINGLE     ); Write(*,*) "Have MPI_THREAD_SINGLE"
         Case( MPI_THREAD_FUNNELED   ); Write(*,*) "Have MPI_THREAD_FUNNELED"
         Case( MPI_THREAD_SERIALIZED ); Write(*,*) "Have MPI_THREAD_SERIALIZED"
         Case( MPI_THREAD_MULTIPLE   ); Write(*,*) "Have MPI_THREAD_MULTIPLE"
         End Select
         
         !Read the coefficient matrix and its characteristics 

     
         Call XFTS_read_matrix( smatrix, matrixfile )

         ! Check if the matrix is square
         If( smatrix%m /= smatrix%n ) iinfo = -2 
         If (iinfo < 0) goto 9999
      
         ! when no symmetry provided the one in the file
         If ( sym < 0 ) sym = smatrix%sym
         ! verify symmetry
         If ( sym /= smatrix%sym ) Then
            If ( ( sym == SM_SYM_IsSPD ).And.( smatrix%sym == SM_SYM_isSymmetric ))Then
               Continue
            Else
               Write(*,*), "Warning : symmetry do not match :",&
               sym, smatrix%sym, "(inputfile matrixfile)"
            End If
         End If

         smatrix%sym = sym
         ftsl%sym  =  sym
         ftsl%n    =  smatrix%n
         ftsl%nnz  =  smatrix%nnz
         ftsl%nrhs =  1
         
         ftsl%rows => smatrix%i
         ftsl%cols => smatrix%j
         ftsl%values => smatrix%v

         !Set the rhs : read from file; on failure, generate one.

          Call XFTS_read_rhs( ftsl%rhs, iinfo, rhsfile )
         
         If (iinfo < 0 ) Then
            write(*,*) "Warning : failed to open rhsfile"
            write(*,*) "Warning : Generates RHS from random solution"
            Call XFTS_gen_sol_prand(smatrix%m, sol )
            Call XFTS_gen_rhs( smatrix, sol, ftsl%rhs )
         End If
         !*allocate memory for the solution
         Allocate( ftsl%sol( smatrix%m ))
      End If

      !-----------------------------------------------------------------------------
      ! broadcast the icntl
      !-----------------------------------------------------------------------------

      Call MPI_Bcast(ftsl%icntl, FTSOLVER_ICNTL_SIZE, MPI_INTEGER, 0, MPI_COMM_WORLD, iinfo)
       If (iinfo /= MPI_SUCCESS) iinfo = -1
       If (iinfo < 0) goto 9999

       Call MPI_Bcast(ftsl%rcntl, FTSOLVER_RCNTL_SIZE, MPI_REAL8, 0, MPI_COMM_WORLD, iinfo)
       If (iinfo /= MPI_SUCCESS) iinfo = -1
       If (iinfo < 0) goto 9999

       
       Call XFTS_analyze(ftsl)
       Call XFTS_compute_aux_matrices(ftsl)
       CAll XFTS_compute_precond(ftsl)

       Do i=1,4
          ftsl%icntl(29) = i
          Call MPI_Barrier(MPI_COMM_WORLD, iinfo)
          If ( i == 0 .and. rank ==0) then
             write(*,*)
             write(*,*) "-----------0:-RESET----------------"
             write(*,*)
          Else if (i == 1 .and. rank ==0) then
             write(*,*)
             write(*,*) "---------1:-LI------------------"
             write(*,*)
          Else if (i == 2 .and. rank ==0) then
             write(*,*)
             write(*,*) "---------2:-LSI------------------"
             write(*,*)
          Else if (i == 3 .and. rank ==0) then
             write(*,*)
             write(*,*) "---------3:-RESTART------------------"
             write(*,*)
          Else if (i == 4 .and. rank ==0) then
             write(*,*)
             write(*,*) "---------4:-REFERENCE------------------"
             write(*,*)
          End If



             Call XFTS_solve(ftsl)
             Call XFTS_state_UpdateStatus(ftsl)
          Call XFTS_state_print( ftsl, iinfo )
       End Do
9999 continue
      Call XFTS_sm_free(smatrix,iinfo)
      If (Associated (ftsl%rhs)) Deallocate(ftsl%rhs)
      If (Associated (sol)) Deallocate(sol)
      If (Associated (ftsl%lc_rhs)) Deallocate(ftsl%lc_rhs)

      ! destroy the instance of ftsl
      Call XFTS_ftsolver_free(ftsl,iinfo)
      Call MPI_Finalize(iinfo)
      
      End Program Test
      
