// The three standard sorting algorithms.
// There are two versions of each: one for a whole array and
// one for part of an array.  
// The ones for a whole array take an additional parameter so they
// an be called by the runSorts code.
//
void insertion_sort_wrapper(int *A,int n, int threads);
// A version that can be used in the middle of an array.
// Both indices l and r are inclusive.
void insertion_sort(int *A,int l,int r);

void merge_sort_wrapper(int *A,int n, int threads);
void merge_sort(int *A,int l,int r);
void merge(int *A,int l,int m,int r);

void quick_sort_wrapper(int *A,int n, int threads);
void quick_sort(int *A,int l,int r);
int partition(int *A, int l, int r);

// For the qsort algorithm.
void built_in_quick_sort_wrapper(int *A,int n, int threads);
