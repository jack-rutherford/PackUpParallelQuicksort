#include<omp.h>
#include "standardSortAlgorithms.h"
#include "betterQuick.h"
#include "util.h"


#include<iostream>
#include<fstream>
using namespace std;
int threadsss = 1;


//---------------------------------------------------------------------
// Partition
//
// Here are the function prototypes for my partitionP algorithm
// that does partition of an array in parallel using one up and one
// down pass.  It is a bit complicated, so feel free to ignore this
// and instead implement a packHigh algorithm and use pack and packHigh
// to implement partition.


// A node struct with a few more useful fields.
struct pnode {
     pnode* left;
     pnode* right;
     int sum; // the number of items in the range < threshold
     int sum2; // the number of items in the range > threshold
     int lo; // left endpoint, inclusive
     int hi; // right endpoint, exclusive
     int fl; // fromLeft for items < threshold
     int fl2; // fromLeft for items > threshold
     pnode() {
        left=0;
        right=0;
        lo=0;
        hi=0;
        sum=0;
        fl=0;
        fl2=0;
        sum2=0;
     }
    pnode(int low,int high) {
        left=0;
        right=0;
        sum=0;
        fl=0;
            fl2=0;
        sum2=0;
        lo=low;
        hi=high;  
    }
    void print() { // for debugging
        cout<<"["<<lo<<", "<<hi<<")  "<<sum<<"   "<<fl<<"|"<<sum2<<"   "<<fl2<<"\n";
    }
};


// A useful method to call at the end of your algorithm to delete the tree.
void deleteTreeP(pnode* root) {
	if(root!=0) {
		#pragma omp task untied
		deleteTreeP(root->left);
		deleteTreeP(root->right);
		#pragma omp taskwait
		delete root;
	}
}


//---------------------------------------------------------------------
//---------------------------------------------------------------------
//---------------------------------------------------------------------
// Pack


// This struct is useful since pack returns an array of unknown size
// and C/C++ doesn't store the array size as part of an array.
struct Array {
    int size;
    int* a;
    Array() {
        size=0;
        a=NULL;
    }
    Array(int *array, int s) {
        size = s;
        a = array;
    }
    void print() {
        printArray(a,size);
    }
};
struct Range {
	int lo;
	int hi;
};

//-----------------------------------------------------------------------------
// thresh is the pivot value.
// root contains the left and right endpoints of the subarray so we don't need
// to send them separately.
void partitionUp(pnode* root, int *in, int thresh, int n) {
	//implement me (or not)
	// implement me.
    // create a bitmap and do step one of parallel prefix step
    static const int SEQUENTIAL_CUTOFF = 50000;
    int lo=root->lo;
    int hi=root->hi;
    // Use a sequential cutoff to do the work in
    // serial when it is small enough
    if (hi-lo < SEQUENTIAL_CUTOFF) {
        root->sum=0;
    	root->sum2 = 0;
        for (int i = lo; i < hi; i++) {
            if(in[i] < thresh)
                root->sum += 1;
     		else if(in[i] > thresh)
       			root->sum2 += 1;
        }
    } else {
        int mid=(hi+lo)/2;
    	pnode* left = new pnode(lo,mid);
        pnode* right = new pnode(mid,hi);
        root->left=left;
        root->right=right;
        // Run first part on another thread
        #pragma omp task untied
        {
            partitionUp(left, in, thresh, n);
        }
        // Run second part on this thread
        partitionUp(right, in, thresh, n);
        #pragma omp taskwait
        root->sum = left->sum + right->sum;
    	root->sum2 = left->sum2 + right->sum2;
    }
}

// out is where the results of partition should go.
// l is the left endpoint of this call to partitionDown.
// l2 is the left endpoint of where the "higher" elements go.
// ****Both l and l2 should be passed as-is in recursive calls. 
// thresh is the pivot value.
//
void partitionDown(pnode* root, int *in, int *out, int l, int l2, int thresh, int n) {
	// implement me (or not)
	pnode* left = root->left;
    if(left!=0) {
        left->fl = root->fl;
    	left->fl2 = root->fl2;
        root->right->fl = root->fl + left->sum;
    	root->right->fl2 = root->fl2 + left->sum2;

	    // cout << "after root assignments" << "\n";
        // Run first part on another thread
        #pragma omp task untied
        {
        	partitionDown(left,in, out,l,l2, thresh, n);
        }
        // Run second part on this thread
        partitionDown(root->right,in, out,l,l2, thresh, n);
		#pragma omp taskwait
    } else {
        // The sequential-cutoff part of the down pass is
        // determined by the nodes that have already been
        // constructed.  Since root has no children, we
        // are at the cutoff and proceed sequentially.
        int lo=root->lo;
		int hi=root->hi;
		//cout << "range: " << lo << " - " << hi << "\n";
		int index = l + root->fl;
		int index2 = l2 + root->fl2;
	    cout << "Index: " << index << "\tIndex2: " << index2 << "\n";
	    cout << "Lo: " << lo << "\tHi: " << hi << "\n";
	    cout << "l: " << l << "\tl2: " << l2 << "\n";
        int i;
		for(i = lo; i < hi; i++){
            // cout << "Inside for loop" << endl;
            // cout << "Thresh: " << thresh << "\ti: " << i << "\t";
            // cout << "size: " << n << "\t";
            // cout<<"in[i]:\t" << in[i] << "\n";
      		if(in[i] < thresh){
				//cout << index << ": " << in[i] << "\n";
        		out[index] = in[i];
				index++;
	            // cout << "After first else" << "\n";

      		}
			else if(in[i] > thresh)
			{
				// cout << "in here" << "\n";
        		out[index2] = in[i];
				index2++;
	            // cout << "After second else" << "\n";
			}
  		}
  }
}

// Partition array in into array out from l (inclusive) to r (exclusive) using
// thresh as the pivot value.
// Return the range of pivot values (inclusive of lo, exclusive of hi).
// That is, the range of indices where the value thresh ended up (in case there
// are more than one).
Range partitionP(int *in,int *out, int l, int r, int thresh, int n) {
	// Implement me (or not)
	
	pnode *root = new pnode(l,r);
	cout << "before partition" << "\n";
	partitionUp(root, in, thresh, n);
	cout << "between partition up and down" << "\n";
	//int *out = new int[r-l+1];
	// cout << "Out array made" << "\n";

	// int out[r-l];
	partitionDown(root, in, out, l, r - root->sum2, thresh, n);
	
	cout << "finished partitionDown" << "\n";

	// #pragma omp for schedule(auto)
	for(int i = l; i < r; i++)
	{
        // cout << "out[i]: " << out[i] << "\t";
		in[i] = out[i];
		// cout << "Index: " << i << "\tin: " << in[i] << "\t" << l << " - " << r << "\n";
	}
	cout << "a" << "\n";

	Range pivots;
	pivots.lo = l + root->sum;
	pivots.hi = r - root->sum2;

	#pragma omp task untied
	{
		deleteTreeP(root);
	}
	
	cout << "b" << "\n";

	for(int i = pivots.lo; i < pivots.hi; i++)
		in[i] = thresh;
	
	cout << "c" << "\n"; 
	#pragma omp taskwait
	return pivots; // Just so it compiles.
}

//correct
void fj_quick_sort(int *A, int *out, int l, int r)
{
    int n = r;
	// cout << "in main" << "\n";
    if(l<r)
    {
        if((r - l) < 50)
        {
            insertion_sort(A, l, r);
			// cout << "insertion sort" << "\n";
        }
        else if((r - l) < 1000)
        {
			// cout << "in sequential fj_quick_sort" << "\n";
            int p = partition(A, l, r);
			// cout << "after sequential partition fj_quick_sort" << "\n";
            fj_quick_sort(A, out, l, p-1);
            fj_quick_sort(A, out, p+1, r);
        }
        else{
            // cout << "Inside else statement\n";
            Range p = partitionP(A, out, l, r, A[l], n);
            #pragma omp task untied
            {
                fj_quick_sort(A, out, l, p.lo);
            }

            fj_quick_sort(A, out, p.hi, r);
			// cout << "finished else statement" << "\n";
            
        }
       
    }


}


void fj_quick_sort_wrapper(int *A,int n,int threads) {
    // Currently just calls the standard quick_sort algorithm.
    // Replace this with your algorithm. It will probably
    // be best to implement methods called better_quick_sort
    // and better_partition and have this algorithm call
    // better_quick_sort and have that use better_partition.

    int *out = new int[n];
    omp_set_num_threads(threads);
    omp_set_nested(true);
    threadsss = threads;
    #pragma omp parallel
    {
        #pragma omp single
        {
            fj_quick_sort(A, out, 0, n);
        }
           
    }
    delete[] out;
}
