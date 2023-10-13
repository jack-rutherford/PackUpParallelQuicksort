#include "betterMerge.h"
#include "standardSortAlgorithms.h"

void better_merge_sort_wrapper(int *A,int n, int threads) { 
    // Currently just calls the standard merge_sort algorithm.
	// Replace this with your algorithm. It will probably 
	// be best to implement a method called better_merge_sort
	// that is your implementation and have this method call 
	// that method after doing any initialization of extra 
	// variables, etc. 
	// You will also need to implement better_merge.
    merge_sort_wrapper(A,n,threads);
}
