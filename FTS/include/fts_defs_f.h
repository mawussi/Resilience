! fortran definitions for ftsolver, common to all arithmetics
#ifndef FTS_DEFS_F_H__
#define FTS_DEFS_F_H__



! Status
#define FTS_SUCCESS         0

! Logical
#define FTS_LOGICAL   INTEGER
#define FTS_TRUE            0
#define FTS_FALSE           1
#define FTS_LOGICAL_UNSET   2

! Integers
#define FTS_INT        Integer
#define FTS_INTKIND         4
#define SCOTCH_NUMKIND      4
#define SCOTCH_SYSKIND      4
#define FTS_INTBYTESIZE     4
#define FTS_INTMPI     MPI_INTEGER


! global size
#define FTSOLVER_STRL        1024
#define FTSOLVER_ICNTL_SIZE  45
#define FTSOLVER_RCNTL_SIZE  40 
#define FTSOLVER_IINFO_SIZE  40
#define FTSOLVER_RINFO_SIZE  10
#define FTSOLVER_IINFOG_SIZE 5
#define FTSOLVER_RINFOG_SIZE 20 
#define FTSOLVER_IKEEP_SIZE  41
#define FTSOLVER_RKEEP_SIZE  40   

#define MTHREAD_ICNTL_SIZE        5
#define ANA_TIMING_SIZE          10

#define SDS_MAX_INDEX      3

#endif 
