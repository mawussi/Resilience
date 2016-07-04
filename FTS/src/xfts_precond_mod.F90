#include "fts_defs_f.h"
#include "fts_macros_f.h"



Module XFTS_precond_mod


      Implicit None 

      Public :: XFTS_Compute_precond

      
      Character(len=FTSOLVER_STRL), Private, Parameter :: &
           FLNAME= "XFTS_ARITHfts_precond_mod.90"
      contains


    ! [+] routine : XFTS_Precond ---------------------------------------------
    !
    !>  
    !!
    !! The preconditionner is a factorization of delta1 matrix
    !! 
    !!
    !! @param[in,out ] ftsl     the ftsolver instance 
    !!
    !!
    Subroutine XFTS_Compute_precond(ftsl)

      !Test 
      Use XFTS_solve_mod
      Use XFTS_dense_matrix_mod


      !* Module *!
      Use FTS_mem_mod
      Use FTS_ftsolver_enum
      Use XFTS_ftsolver_type
      Use XFTS_mumps_mod
      Implicit None
      Include 'mpif.h'
            

      !* Arguments *!
      Type(XFTS_ftsolver_t), Intent(inout) :: ftsl

      !* Local variables *!
      
      ! Scalars


      Integer(kind=8) :: piv
      Real(kind=8) :: flops_estielim
      Real(kind=8) :: flops_assemb
      Real(kind=8) :: flops_elim
      Integer      :: Not_needed
      Integer      :: rank
      Integer      :: info, i, n

      ! Derived types
      Type(FTS_mem_t)  :: sds_mem
      ! End of header ----------------------------------------------------------

      rank   = ftsl%ikeep(IKEEP_MPIRANK)
      !-------------------------------------------------------------------------
      ! [] Compute PCD <- Factorize(SpSBar)
      !-------------------------------------------------------------------------
      !

      Call XFTS_mumps_Set_MPICommunicator(ftsl%sm_precond%sds,MPI_COMM_SELF,info)
      CHCKASSRT( info >= 0, info )
      If (info < 0 ) Return

      Call XFTS_mumps_Set_matrix (ftsl%sm_precond%sds,ftsl%sm_precond%sm_A,info)
      CHCKASSRT( info >= 0, info )
      If (info < 0 ) Return

      Call XFTS_mumps_Analyze   (ftsl%sm_precond%sds,info)
      CHCKASSRT( info >= 0, info )
      If (info < 0 ) Return
      
      Call XFTS_mumps_Factorize (ftsl%sm_precond%sds,info)
      CHCKASSRT( info >= 0, info )
      If (info < 0 ) Return


      !-------------------------------------------------------------------------
      ! [5] Save statistics
      !-------------------------------------------------------------------------



      Call XFTS_mumps_get_numstats&
           (ftsl%sm_precond%sds,piv)
      Call XFTS_mumps_get_memstats&
           (ftsl%sm_precond%sds,sds_mem)
      Call FTS_mem_add2mem(ftsl%mem,sds_mem)




      ftsl%iinfo (IINFO_PCD_SIZEOF   ) = &
           Byte2MByte(FTS_mem_getallusage(sds_mem))
      ftsl%iinfo (IINFO_PCD_SDSMEMPEAK    ) = &
           Byte2MByte(FTS_mem_getallpeak   (sds_mem))

      ftsl%iinfo (IINFO_PCD_NBPIVOTS  ) = INT(piv,4) 


      
      !-------------------------------------------------------------------------
      ! [4] Finish
      !-------------------------------------------------------------------------

9999 continue
    End Subroutine XFTS_Compute_precond


End Module XFTS_precond_mod
