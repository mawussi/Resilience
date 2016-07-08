# FTS
**FTS** stands  for Fault Tolerant Solver. FTS is a  **_Fortran 90_**  prototype software designed to evaluate the potentials of numerical resilient algorithms in parallel distributed memory environments. Full information on the numerical resilient algorithms are detailed in [paper](https://hal.inria.fr/hal-01323192/file/final_nlaa.pdf).
In order to evaluate the numerical behaviour and assessing the computation time the new resilient algorithm, we use a fault injection strategy based on fault simulation. Process crashes are simulated by overwriting the critical data. Data overwriting is aims mainly to lead to data loss condition which is the main consequence when a fault occurs in real world applications. This project focuses on linear algebra solver, more precisely the Krylov subspace methods for large sparse linear systems. The numerical resilient strategies developed in this project have demonstrated promising results for GMRES, FGMRES, CG, BiCGStaB, but here we provide the code only for GMRES and FGMRES. 
    
The main folder **FTS** contains the following subfolders.

* **src** contains subroutines for fault simulation and data recovery. 

* **packfgmres** contains the **FGMRES** solver subroutine and the associated auxiliary subroutines. 

* **packgmres** contains the  **GMRES** solver subroutine and the associated auxiliary subroutines.

* **sm** a toolkit for sparse matrix operations. 

* **dm** a toolkit for dense matrix operations.

* **part** Contains  subroutines for graph partitioning, to subdivide the input matrix into sub-domains and compute the necessary data for the computation of the preconditioner.

* **toolkit** I/O toolkit for reading/writing and parsing files.

* **slatec** tookit for vector operations, contains files : fdump.f  i1mach.f  icopy.f  ipsort.f  isort.f  j4save.f xercnt.f  xerhlt.f  xermsg.f  xerprn.f  xersve.f  xgetua.f available at : [http://www.netlib.org/slatec/](http://www.netlib.org/slatec/).

* **common** contains auxiliary subroutines for timing, error messages printing, measuring memory consuming etc.    

* **modules** contains .mod files generated at compilation time.

* **lib** contains the libraries generated once FTS is compiled.

# Requirements

#### Graph partitioning 
* **Metis**: available at [metis-5.1.0.tar.gz](http://glaros.dtc.umn.edu/gkhome/fetch/sw/metis/metis-5.1.0.tar.gz).
* **Scotch**: available at [scotch_6.0.4.tar.gz](http://gforge.inria.fr/frs/download.php/file/34618/scotch_6.0.4.tar.gz).
* **Colamd** available at [SuiteSaprse](https://github.com/PetterS/SuiteSparse)

#### Sparse solver
* **Mumps** for sparse LU factorization available at [ http://mumps.enseeiht.fr/MUMPS_5.0.1.tar.gz](MUMPS_5.0.1.tar.gz)
* **QR_Mumps** for sparse QR factorization available at [http://buttari.perso.enseeiht.fr/qr_mumps/](http://buttari.perso.enseeiht.fr/qr_mumps/)

#### BLAS & LAPACK
* **MKL** is now free for academics (students and researchers) available at [https://software.intel.com/en-us/articles/free-mkl](https://software.intel.com/en-us/articles/free-mkl)

OR
* **Netlib BLAS** no optimized BLAS routines, available at [BLAS-3.6.0.tgz](http://www.netlib.org/blas/blas-3.6.0.tgz)
* **Netlib LAPACK** no optimized LAPACK routines, available at [LAPACK-3.6.1.tgz](http://www.netlib.org/lapack/lapack-3.6.1.tgz) 
 
#### Hardware locality
* **Hwloc** for more information on the hardware topology,  available at [hwloc-1.11.3.tar.gz](https://www.open-mpi.org/software/hwloc/v1.11/downloads/hwloc-1.11.3.tar.gz)


 

