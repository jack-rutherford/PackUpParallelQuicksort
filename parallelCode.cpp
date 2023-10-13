#include<iostream>
#include<fstream>
#include<cstdlib>
#include<omp.h>
using namespace std;

#include<chrono>
#include "util.h"

//---------------------------------------------------------------------
// The sumWorker and sum code taken from SIPAC OpenMP version of notes.
// The rest of the code by Chuck Cusack, November 2015.
//

//---------------------------------------------------------------------
// Parallel sum of the elements of an array.
// Using the standard divide-and-conquer algorithm
int sumWorker(int arr[], int lo, int hi) {
	static const int SEQUENTIAL_CUTOFF = 10000;
	if (hi - lo < SEQUENTIAL_CUTOFF) {
		int result = 0;
		for (int i = lo; i < hi; i++) {
			result += arr[i];
		}
		return result;
	} else {
		int leftSum, rightSum;
		// Run half of the work on another thread
		#pragma omp task untied shared(leftSum)
		leftSum = sumWorker(arr, lo, lo+(hi-lo)/2);

		// Run half of the work on this thread
		rightSum = sumWorker(arr, lo+(hi-lo)/2, hi);
	
		// Wait until the other thread is done
		#pragma omp taskwait

		return leftSum + rightSum;
	}
}

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

int sum(int arr[], int len, int threads) {
	int result;
	omp_set_num_threads(threads);
	omp_set_nested(1);

	#pragma omp parallel
	#pragma omp single
	{
	// The initial call is in a parallel block but as a single
	// The threads will be spawned by sumWorker.
	result = sumWorker(arr, 0, len);
	}
	return result;
}

//---------------------------------------------------------------------
// Parallel Prefix Sum
// Start with a node to store the data we need to construct the tree
struct node {
     node* left;
     node* right;
     int sum;
     int lo; // left endpoint, inclusive
     int hi; // right endpoint, exclusive
     int fl; // fromLeft
     node() {
        left=0;
        right=0;
	lo=0;
	hi=0;
	sum=0;
	fl=0;
     }
	node(int low,int high) {
		left=0;
		right=0;
		sum=0;
		fl=0;
		lo=low;
		hi=high;   
	}
	void print() { // for debugging.
		cout<<"["<<lo<<", "<<hi<<")  "<<sum<<"   "<<fl<<"\n";
	}
};

// The "Up" pass of the paralle-prefix sum algorithm
void prefixSumUp(node* root, int in[]) {
	static const int SEQUENTIAL_CUTOFF = 50000;
	int lo=root->lo;
	int hi=root->hi;
	// Use a sequential cutoff to do the work in 
	// serial when it is small enough
	if (hi-lo < SEQUENTIAL_CUTOFF) {
		root->sum=0;
		for (int i = lo; i < hi; i++) {
			root->sum += in[i];
		}
	} else {
		int mid=(hi+lo)/2;
        node* left = new node(lo,mid);
  		node* right = new node(mid,hi); 
		root->left=left;
		root->right=right;
		// Run first part on another thread
		#pragma omp task untied
		{
			prefixSumUp(left,in);
		}
		// Run second part on this thread
		prefixSumUp(right,in);
		#pragma omp taskwait
		root->sum = left->sum+right->sum;
	}
}
void prefixSumDown(node* root, int in[], int out[]) {
	node* left = root->left;
	if(left!=0) {
		left->fl = root->fl;
		root->right->fl = root->fl + left->sum;

		// Run first part on another thread
		#pragma omp task untied
		{
			prefixSumDown(left,in,out);
		}
		// Run second part on this thread
		prefixSumDown(root->right,in,out);
	} else {
		// The sequential-cutoff part of the down pass is
		// determined by the nodes that have already been
		// constructed.  Since root has no children, we
		// are at the cutoff and proceed sequentially.
		int lo=root->lo;
		int hi=root->hi;
		out[lo]=root->fl+in[lo];
		for(int i=lo+1;i<hi;i++) {
			out[i]=in[i]+out[i-1];
		}
	}
}

// A useful method to call at the end of your algorithm to delete the tree.
void deleteTree(node* root) {
	if(root!=0) {
		#pragma omp task untied
		deleteTree(root->left);
		deleteTree(root->right);
		#pragma omp taskwait
		delete root;
	}
}

void deleteTreeP(pnode* root) {
	if(root!=0) {
		#pragma omp task untied
		deleteTreeP(root->left);
		deleteTreeP(root->right);
		#pragma omp taskwait
		delete root;
	}
}

// The main parallel prefix sum algorithm.
void prefixSumP(int in[],int out[], int len, int threads) {
	omp_set_num_threads(threads);
	omp_set_nested(1);
	node *root = new node(0,len);
	#pragma omp parallel
	#pragma omp single
	{
		prefixSumUp(root, in);
		prefixSumDown(root, in, out);
	}
}

// The sequential prefix sum algorithm
void prefixSum(int in[],int out[], int len) {
	out[0]=in[0];
	for(int i=1;i<len;i++) {
		out[i]=in[i]+out[i-1];
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
void packUp(node* root, int in[], int thresh) {
	// implement me.
	// create a bitmap and do step one of parallel prefix step
	static const int SEQUENTIAL_CUTOFF = 8;
	int lo=root->lo;
	int hi=root->hi;
	// Use a sequential cutoff to do the work in 
	// serial when it is small enough
	if (hi-lo < SEQUENTIAL_CUTOFF) {
		root->sum=0;
		for (int i = lo; i < hi; i++) {
			if(in[i] < thresh)
			{
				root->sum += 1;
			}
		}
	} else {
		int mid=(hi+lo)/2;
    	node* left = new node(lo,mid);
		node* right = new node(mid,hi); 
		root->left=left;
		root->right=right;
		// Run first part on another thread
		#pragma omp task untied
		{
			packUp(left, in, thresh);
		}
		// Run second part on this thread
		packUp(right, in, thresh);
		#pragma omp taskwait
		root->sum = left->sum+right->sum;
	}
	}

void packDown(node* root, int in[], int out[], int thresh) {
	node* left = root->left;
	if(left!=0) {
		left->fl = root->fl;
		root->right->fl = root->fl + left->sum;

		// Run first part on another thread
		#pragma omp task untied
		{
			packDown(left,in,out, thresh);
		}
		// Run second part on this thread
		packDown(root->right,in,out, thresh);
		#pragma omp taskwait
	} else {
		// The sequential-cutoff part of the down pass is
		// determined by the nodes that have already been
		// constructed.  Since root has no children, we
		// are at the cutoff and proceed sequentially.
		
		int lo=root->lo;
		int hi=root->hi;
		int index = root->fl;
		
		for(int i = lo; i < hi; i++){
      		if(in[i] < thresh){
        		out[index] = in[i];
				index++;
      		}
  		}
 
	}
}

Array packP(int in[],int len, int thresh, int threads) {
	// Set number of threads, set nested, create the root node,
	// call packUp, create an array of the appropriate size to
	// pass to packDown, then call packDown.  
	// Finally, return an Array object based on the out array
	// sent to packDown.
 
	omp_set_num_threads(threads);
	omp_set_nested(1);
	node *root = new node(0,len); 
	
	int length = 0;
	int *out= new int[len];
	#pragma omp parallel
	#pragma omp single
	{
		packUp(root, in, thresh);
		
		length = root->sum;
		packDown(root, in, out, thresh);

		deleteTree(root);
	}

	return Array(out, length);
	
	 // just so it compiles.
}
// The sequential pack algorithm.
// You shouldn't use this in your pack, but you SHOULD use it
// to test against your algorithm.
Array pack(int in[], int len, int thresh) {
	int size=0;
	for(int i=0;i<len;i++) {
		if(in[i]<thresh) size++;
	}
	int j=0;
	int *out = new int[size];
	for(int i=0;i<len;i++) {
		if(in[i]<thresh) {
			out[j]=in[i];
			j++;
		}
	}
	return Array(out,size);
}

// Useful to return the range of pivot values in case there
// are more than one.  My partition 
struct Range {
	int lo;
	int hi;
};



// thresh is the pivot value.
// root contains the left and right endpoints of the subarray so we don't need
// to send them separately.
void partitionUp(pnode* root, int *in, int thresh) {
	//implement me (or not)
	// implement me.
    // create a bitmap and do step one of parallel prefix step
    static const int SEQUENTIAL_CUTOFF = 8;
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
     		if(in[i] > thresh)
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
            partitionUp(left, in, thresh);
        }
        // Run second part on this thread
        partitionUp(right, in, thresh);
        #pragma omp taskwait
        root->sum = left->sum+right->sum;
    	root->sum2 = left->sum2+ right->sum2;
    }
}

// out is where the results of partition should go.
// l is the left endpoint of this call to partitionDown.
// l2 is the left endpoint of where the "higher" elements go.
// ****Both l and l2 should be passed as-is in recursive calls. 
// thresh is the pivot value.
//
void partitionDown(pnode* root, int *in, int *out, int l, int l2, int thresh) {
	// implement me (or not)
	pnode* left = root->left;
    if(left!=0) {
        left->fl = root->fl;
    	left->fl2 = root-> fl2;
        root->right->fl = root->fl + left->sum;
    	root->right->fl2 = root->fl2 + left->sum2;


        // Run first part on another thread
        #pragma omp task untied
        {
        	partitionDown(left,in,out,l,l2, thresh);
        }
        // Run second part on this thread
        partitionDown(root->right,in,out,l,l2, thresh);
		#pragma omp taskwait
    } else {
        // The sequential-cutoff part of the down pass is
        // determined by the nodes that have already been
        // constructed.  Since root has no children, we
        // are at the cutoff and proceed sequentially.
        int lo=root->lo;
		int hi=root->hi;
		cout << "range: " << lo << " - " << hi << "\n";
		int index = l + root->fl;
		int index2 = l2 + root->fl2;
		for(int i = lo; i < hi; i++){
      		if(in[i] < thresh){
				cout << index << ": " << in[i] << "\n";
        		out[index] = in[i];
				index++;
      		}
			else if(in[i] > thresh)
			{
				cout << index2 << ": " << in[i] << "\n";
        		out[index2] = in[i];
				index2++;
			}
  		}
  }
}

// Partition array in into array out from l (inclusive) to r (exclusive) using
// thresh as the pivot value.
// Return the range of pivot values (inclusive of lo, exclusive of hi).
// That is, the range of indices where the value thresh ended up (in case there
// are more than one).
Range partitionP(int *in, int *out, int l, int r, int thresh) {
	// Implement me (or not)
	pnode *root = new pnode(l,r);
	partitionUp(root, in, thresh);
	partitionDown(root, in, out, l, r - root->sum2, thresh);
	cout << "finished partitionDown" << "\n";
	#pragma omp for
	for(int i = l; i < r; i++)
	{
		in[i] = out[i];
	}

	Range pivots;
	pivots.lo = root->lo + root->sum;
	pivots.hi = root->hi - root->sum2;

	cout << "pivot low: " << pivots.lo << endl;
	cout << "pivot high: " << pivots.hi << endl;


	#pragma omp task untied
	deleteTreeP(root);

	for(int i = pivots.lo; i < pivots.hi; i++)
		in[i] = thresh;

	cout << "pivot set" << endl;
	
	return pivots; // Just so it compiles.
}

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------
void testSum () {
	typedef std::chrono::high_resolution_clock Clock;
	int size = 10000000;
	//int threads = 2;
	int *a = new int[size];
	for(int i=0;i<size;i++) {
		a[i] = rand();
	}
	for(int i=1;i<=4;i++) {
		auto startTime = Clock::now();
		int s = sum(a,size, i);
		auto endTime = Clock::now();
		auto elapsedTime =std::chrono::duration_cast<std::chrono::nanoseconds>(endTime - startTime).count();
		long CPUTime=(long)((elapsedTime*1.0)/1000);
		cout<<"sum: "<<s<<"  Threads: "<<i<<"  Time: "<<CPUTime<<"\n";
	}
}
//-------------------------------------------------------------------------
void testPrefixSum() {
	typedef std::chrono::high_resolution_clock Clock;
	int size = 100000000;
	int *a = new int[size];
	int *b = new int[size];
	int *c = new int[size];
	for(int i=0;i<size;i++) {
		a[i] = rand()%20;
	}
	prefixSum(a,c,size);

	for(int i=1;i<=4;i++) {
		clearArray(b,size);
		auto startTime = Clock::now();
		prefixSumP(a,b,size,i);
		auto endTime = Clock::now();
		if(!equals(b,c,size)) {
			cout<<"**********Algorithm did not work!*********\n";
		}

		auto elapsedTime =std::chrono::duration_cast<std::chrono::nanoseconds>(endTime - startTime).count();
		long CPUTime=(long)((elapsedTime*1.0)/1000);
		cout<<"Threads: "<<i<<"  Time: "<<CPUTime<<"\n";
		printArray(b,size);
	}
}
//-------------------------------------------------------------------------
void testPack() {
	typedef std::chrono::high_resolution_clock Clock;
	int size = 30;
	int thresh = 10;
	int max = 20;
	int *a = new int[size];
	for(int i=0;i<size;i++) {
		a[i] = rand()%max;
	}
	cout<<"A: ";
	printArray(a,size);
	Array c = pack(a, size,thresh);
	cout<<"C: ";
	c.print();

	for(int i=1;i<=4;i++) {
		auto startTime = Clock::now();
		Array b = packP(a,size,thresh,i);
		cout<<"B: ";
		b.print();
		auto endTime = Clock::now();
		if(b.size!=c.size || !equals(b.a,c.a,b.size)) {
			cout<<"**********Algorithm did not work!*********\n";
		}

		auto elapsedTime =std::chrono::duration_cast<std::chrono::nanoseconds>(endTime - startTime).count();
		long CPUTime=(long)((elapsedTime*1.0)/1000);
		cout<<"Threads: "<<i<<"  Time: "<<CPUTime<<"\n";
		//printArray(b.a,size);
	}
}
//-------------------------------------------------------------------------
void testPartition() {
	typedef std::chrono::high_resolution_clock Clock;
	int size = 30;
	int thresh = 10;
	int max = 20;
	int *a = new int[size];
	for(int i=0;i<size;i++) {
		a[i] = rand()%max;
	}
	cout<<"A: ";
	for(int i = 0; i < size; i++)
		cout << a[i] << " ";
	cout << "\n";
	
	int out[size];
	Range c = partitionP(a, out, 0, size - 1, thresh);
	cout << "B: ";
	for(int i = 0; i < size; i++)
		cout << a[i] << " ";
	cout << "\n";
}
//------------------------------------------------------------------------
int main (int argc,char* argv[]) {
	//testSum();
	//cout << "Sum worked" << endl;
	//testPrefixSum();
	//cout << "prefix sum worked" << endl;
	//testPack();	
	testPartition();
}
