!! ##############################################################################################
!!
!! Copyright 2012 CNRS, INPT
!!  
!! This file is part of qr_mumps.
!!  
!! qr_mumps is free software: you can redistribute it and/or modify
!! it under the terms of the GNU Lesser General Public License as 
!! published by the Free Software Foundation, either version 3 of 
!! the License, or (at your option) any later version.
!!  
!! qr_mumps is distributed in the hope that it will be useful,
!! but WITHOUT ANY WARRANTY; without even the implied warranty of
!! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!! GNU Lesser General Public License for more details.
!!  
!! You can find a copy of the GNU Lesser General Public License
!! in the qr_mumps/doc directory.
!!
!! ##############################################################################################
!!##############################################################################################
!> @file smalltest.F90
!! Basic example
!!
!! $Date: 2012-03-20 14:20:07 +0100 (mar. 20 mars 2012) $
!! $Author: abuttari $
!! $Version 0.0.1$
!! $Revision: 345 $
!!
!!##############################################################################################


program _qrm_test_small

!  use _qrm_mod
  implicit none

  type(_qrm_spmat_type)          :: qrm_mat
  integer                        :: ierr, nargs, i, nrhs
  character                      :: matfile*30='', transp 
  _qrm_data, allocatable, target :: b(:), x(:), r(:)
  _qrm_real                      :: rnrm, onrm

  ! initialize the control data structure. 
  call qrm_spmat_init(qrm_mat)

  call qrm_set('qrm_eunit', 6)
  
  ! allocate arrays for the input matrix
  call qrm_palloc(qrm_mat%irn, 13)
  call qrm_palloc(qrm_mat%jcn, 13)
  call qrm_palloc(qrm_mat%val, 13)

  ! initialize the input matrix
  qrm_mat%jcn = (/1,1,1,2,2,3,3,3,3,4,4,5,5/)
  qrm_mat%irn = (/2,3,6,1,6,2,4,5,7,2,3,2,4/)
  qrm_mat%val = (/0.7,0.6,0.4,0.1,0.1,0.3,0.6,0.7,0.2,0.5,0.2,0.1,0.6/)
  qrm_mat%m   = 7
  qrm_mat%n   = 5
  qrm_mat%nz  = 13

  ! set some control parameters
  call qrm_set(qrm_mat, 'qrm_ib',2)
  call qrm_set(qrm_mat, 'qrm_nb',2)
  call qrm_set(qrm_mat, 'qrm_nthreads',1)
  call qrm_set(qrm_mat, 'qrm_ordering',qrm_natural_)

  write(*,'(30("="))')

  write(*,'("Starting Analysis")')
  call qrm_analyse(qrm_mat, transp)

  write(*,'("Starting Factorization")')
  call qrm_factorize(qrm_mat, transp)

  call qrm_aalloc(b, qrm_mat%m)
  call qrm_aalloc(r, qrm_mat%m)
  call qrm_aalloc(x, qrm_mat%n)
  
  b = _qrm_one
  ! as by is changed when applying Q', we save a copy in r for later use
  r = b

  call qrm_apply(qrm_mat, 't', b)
  call qrm_solve(qrm_mat, 'n', b, x)

  ! compute the residual
  call qrm_residual_norm(qrm_mat, r, x, rnrm)
  call qrm_residual_orth(qrm_mat, r, onrm)   
  write(*,'("||r||/||A||    = ",e10.2)')rnrm
  write(*,'("||A^tr||/||r|| = ",e10.2)')onrm

  call qrm_adealloc(b)
  call qrm_adealloc(r)
  call qrm_adealloc(x)
  call qrm_spmat_destroy(qrm_mat, all=.true.)
  write(*,'("Done.")')
  write(*,'(" ")')
  write(*,'(" ")')
 
  write(*,'("  Nonzeroes in R           : ",i20)')qrm_mat%gstats(qrm_nnz_r_)
  write(*,'("  Total flops at facto     : ",i20)')qrm_mat%gstats(qrm_facto_flops_)
  write(*,'("  Total unallocated memory : ",i20)')qrm_tot_mem
  write(*,'("  Memory peak              : ",f9.3," MB")') &
       &real(qrm_max_mem,kind(1.d0))/1024.d0/1024.d0

  stop

end program _qrm_test_small

