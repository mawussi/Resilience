!This program reads input files
! and extracts configuration informations 
! for a given domain/partition

#include"fts_defs_f.h"
      
Program Bench_qrm
      Use XFTS_dense_matrix_mod
      Use XFTS_sparse_matrix_mod
      Use XFTS_toolkit_mod
      Use XFTS_mumps_mod 
      Use XFTS_ARITHqrm_mod
      Implicit None 
      Include 'mpif.h'

      Integer        :: istep
      Integer        :: iinfo
      Integer        :: rank
      Integer        :: unit
      Integer        :: sym     ! symmetry of the input matrix
      Integer        :: job
      Integer        :: mpi_provided
      Integer        :: nbdom, nthreads
      FTS_INT        :: i
      Real(Kind=8) :: time, time_facto, time_analyze, flop, flop_elim, flop_ass

      !Strings to read input file
      Character(len=FTSOLVER_STRL) :: paramfile
      Character(len=FTSOLVER_STRL)   :: prefix
      Character(len=FTSOLVER_STRL) :: matrixfile
      Character(len=charlen)  :: keystr
      Character(len=charlen)  :: string
      Character(len=FTSOLVER_STRL) :: nbdom_str
      character                       :: matfile*30='', transp 

      !Derived types
      Type(XFTS_sparse_matrix_t) :: smatrix
      type(XFTS_ARITHqrm_spmat_type)          :: qrm_mat
      Type(XFTS_ARITHmumps_struc)     :: id_mumps

        !- End of header -------------------------------------------------------------

      !-----------------------------------------------------------------------------
      ! [1] Initialize parallelism
      !-----------------------------------------------------------------------------
  
      Call MPI_Init_thread(MPI_THREAD_MULTIPLE,mpi_provided,iinfo) !
      istep = 11
      if (iinfo /= MPI_SUCCESS) iinfo = -1
      if (iinfo < 0) goto 9999
      write(*,*) "provided", mpi_provided

      Call MPI_Comm_rank(MPI_COMM_WORLD, rank, iinfo)
      istep = 12
      if (iinfo /= MPI_SUCCESS) iinfo = -1
      if (iinfo < 0) goto 9999


      ! Read argument 
      Write(nbdom_str,*) nbdom+1
      nbdom_str=Adjustl(nbdom_str)
      Call Get_command_argument(1, paramfile )
      Call Get_command_argument(2, prefix    )
         
      !Open parameter file and read matrix file' name
      unit = 11
      Open(UNIT=unit, FILE=paramfile, STATUS='old', ACTION='read')
      

      keystr="MATFILE"
      Call FTS_parser_read_string (unit, keystr, string )
      matrixfile = Trim(string)
      Close(unit)
         
      Write(*,*) 'arith    is :', Trim('XFTS_FLOAT')
      Write(*,*) 'prefix   is :', Trim(prefix)
      Write(*,*) 'matrix   is :', Trim(matrixfile)
      
      
      Call XFTS_read_matrix( smatrix, matrixfile )
      
      write(*,*)"*********Qr_mumps( QR )******************"
      write(*,*)
      !* Call to qr_mumps!
      !
      !********************
      call qrm_spmat_init(qrm_mat)
      call qrm_set('qrm_eunit', 6)

      call qrm_palloc(qrm_mat%irn, smatrix%nnz)
      call qrm_palloc(qrm_mat%jcn, smatrix%nnz)
      call qrm_palloc(qrm_mat%val, smatrix%nnz)

      qrm_mat%irn => smatrix%i
      qrm_mat%jcn => smatrix%j
      qrm_mat%val => smatrix%v

      qrm_mat%m   = smatrix%m
      qrm_mat%n   = smatrix%n
      qrm_mat%nz  = smatrix%nnz
    
      call qrm_set(qrm_mat, 'qrm_ib',40)
      call qrm_set(qrm_mat, 'qrm_nb',40)
      call qrm_set(qrm_mat, 'qrm_ordering',qrm_metis_)
      
      
      call qrm_set(qrm_mat, 'qrm_nthreads',1)
      time_analyze = MPI_Wtime()
      call qrm_analyse(qrm_mat, transp)
      time_analyze = MPI_Wtime() - time_analyze

      time_facto = MPI_Wtime()
      call qrm_factorize(qrm_mat, transp)
      time_facto = MPI_Wtime() - time_facto

      flop =  qrm_mat%gstats(qrm_facto_flops_)
      
      write(*,'("  Time to do the analysis       : ",es10.3)')time_analyze 
      write(*,'("  Time to do the factorization  : ",es10.3)')time_facto
      write(*,'("  Time Analyze + factorization  : ",es10.3)')time_facto+time_analyze
      write(*,'("  Total flops at facto          : ",1PE11.3)')flop
      write(*,*)
      write(*,'("  Nonzeroes in R                : ",i11)')qrm_mat%gstats(qrm_nnz_r_)
      write(*,'("  Nonzeroes in H                : ",i11)')qrm_mat%gstats(qrm_nnz_h_)

      !* Call to mumps 
      !
      !**************
      if (smatrix%m > smatrix%n ) goto 9999
      write(*,*) "********     MUMPS( LU )   **************"
    Call XFTS_mumps_Set_MPICommunicator(id_mumps,MPI_COMM_SELF,iinfo)

    Call XFTS_mumps_Set_matrix (id_mumps,smatrix, iinfo)
    
    time_analyze = MPI_Wtime()
    Call XFTS_mumps_Analyze   (id_mumps,iinfo)
    time_analyze = MPI_Wtime() - time_analyze

    time_facto = MPI_Wtime()
    Call XFTS_mumps_Factorize (id_mumps,iinfo)
    time_facto = MPI_Wtime() - time_facto

    flop_elim = id_mumps%rinfog(3)
    flop_ass = id_mumps%rinfog(2)

      write(*,'("  Time to do the analysis       :",es10.3)')time_analyze 
      write(*,'("  Time to do the factorization  : ",es10.3)')time_facto
      write(*,'("  Time Analyze + factorization  : ",es10.3)')time_facto+time_analyze
      write(*,*)
      write(*,'("  Total flops at assembly       : ",1PE11.3)')flop_ass
      write(*,'("  Total flops at elimination    : ",1PE11.3)')flop_elim
      write(*,'("  Total entries in the factors  : ",i11)')id_mumps%infog(29)
      write(*,*)
      
    
9999 continue
      call qrm_spmat_destroy(qrm_mat, all=.true.) 
     ! Call XFTS_sm_free(smatrix,iinfo)
      Call MPI_Finalize(iinfo)
      
      End Program Bench_qrm
      
