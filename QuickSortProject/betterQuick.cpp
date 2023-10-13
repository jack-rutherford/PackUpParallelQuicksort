#include "betterQuick.h"
#include "standardSortAlgorithms.h"

//-----------------------------------------------------------------------------
// ignores threads
void better_quick_sort_wrapper(int *A,int n,int threads) {
	// Currently just calls the standard quick_sort algorithm.
	// Replace this with your algorithm. It will probably 
	// be best to implement methods called better_quick_sort
	// and better_partition and have this algorithm call
	// better_quick_sort and have that use better_partition.
	 quick_sort_wrapper(A,n,threads);
}