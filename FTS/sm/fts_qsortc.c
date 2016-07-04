/*  */
#include <stdio.h>
#include <stdlib.h>

/**
 * compare if "a" is superior to "b"
 *
 * @param [in] a The integer to be compared
 * @param [in] b The integer to be compared to
 * @return "a-b"
 */
int fts_qsort_comp(const void * a,const void * b)
{
  return ( (*(int *)a) - (*(int *)b) );
}

/**
 *  sort an array and get its permutation using standard quicksort.
 *  each entry is 2 ints : [ value, index ]
 *
 * @param [in,out] arrayAndperm
 *        an array of "2*size" elements
 *        Containing justaposedly [value, index].
 *        That is to say, for i=0 to size,
 *        arrayAndperm(2*i) is the array, perm(2*i+1) is the permutation.
 *        - On input  : the array is not sorted and perm(i) should be equals to "i+1"
 *        - On output : the array is sorted     and perm(i) contains the permutation made.
 * @param [in    ] size
 *        The number of elements in the array "array"
 * 
*/
void fts_qsortc(int *arrayAndperm, int* size)
{
  qsort(arrayAndperm,(*size), 2*sizeof(int), fts_qsort_comp) ;
}

/* Fortran binding to maphys_qsortc */
void FTS_QSORTC(int *arrayAndperm, int* size)
{
  fts_qsortc(arrayAndperm,size);
}

void fts_qsortc_(int *arrayAndperm, int* size)
{
  fts_qsortc(arrayAndperm,size);
}

void fts_qsortc__(int *arrayAndperm, int* size)
{
  fts_qsortc(arrayAndperm,size);
}


