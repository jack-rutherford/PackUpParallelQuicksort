#include<omp.h>
#include "standardSortAlgorithms.h"
#include "betterQuick.h"
#include "util.h"


#include<iostream>
#include<fstream>
using namespace std;
int threadss = 1;


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
void deleteTree(pnode* root) {
    if(root!=0) {
        #pragma omp task untied
        deleteTree(root->left);
        deleteTree(root->right);
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
void packUp(pnode* root, int in[], int thresh) {
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
            packUp(left, in, thresh);
        }
        // Run second part on this thread
        packUp(right, in, thresh);
        #pragma omp taskwait
        root->sum = left->sum+right->sum;
    root->sum2 = left->sum2+ right->sum2;
    }
}
void packDown(pnode* root, int in[], int out[], int thresh, int len) {
    // implement me
  pnode* left = root->left;
    if(left!=0) {
        left->fl = root->fl;
    left->fl2 = root-> fl2;
        root->right->fl = root->fl + left->sum;
    root->right->fl2 = root->fl2 + left->sum2;


        // Run first part on another thread
        #pragma omp task untied
        {
        packDown(left,in,out, thresh, len);
        }
        // Run second part on this thread
        packDown(root->right,in,out, thresh, len);
    } else {
        // The sequential-cutoff part of the down pass is
        // determined by the nodes that have already been
        // constructed.  Since root has no children, we
        // are at the cutoff and proceed sequentially.
        int lo=root->lo;
        int hi=root->hi;
    int index = 0;
        for(int i = lo; i < hi; i++){
      if(in[i] < thresh){
        index = root->sum + root->fl;
        out[index - 1] = in[i];
      }
     else if(in[i] > thresh)
     {
       index = root->sum2 + root->fl2;
       out[len - index] = in[i];
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
   int out[len];
     pnode *root = new pnode(0,len);


    #pragma omp parallel
    #pragma omp single
    {
        packUp(root, in, thresh);
    for(int i = root->sum; i < len - root->sum2; i++)
    {
      out[i] = thresh;
    }
        packDown(root, in, out, thresh, len);
        cout << "Delete tree" << endl;
        deleteTree(root);
    }
    return Array(out,root->sum + root->sum2); // just so it compiles.
}


//-----------------------------------------------------------------------------
int fj_partition(int *A, int l, int r)
{
  cout << "in partition" << endl;
  cout << r << endl;
  cout << l << endl;
  int length = r - l;
  int* B[length];
  cout << "initialized B" << endl;
  #pragma omp for
  for(int i = 0; i < length; i++){
    cout << "in loop1" << endl;
      *B[i] = A[l+i];
      cout << "in loop2" << endl;
  }
  cout << "a" << endl;
  Array p = packP(*B, length, A[l], threadss);
  //B = p.a;
  int thresh = A[l];
  int location = -1;
  cout << "b" << endl;
  #pragma omp for
    for(int i = l; i < r; i++){
      A[i] = p.a[i-l];
      if(p.a[i-l] == thresh){
        location = i;
      }
    }
    //delete[] B;
  return location;
}


void fj_quick_sort(int *A, int l, int r)
{
    if(l<r)
    {
        if((r - l) < 50)
        {
            insertion_sort(A, l, r);
        }
        else if((r - l) < 1000)
        {
            int p = fj_partition(A, l, r);
            fj_quick_sort(A,l, p-1);
            fj_quick_sort(A, p+1, r);
        }
        else{
            cout << "Inside else statement\n";
            int p = fj_partition(A, l, r);
            #pragma omp task untied
            {
                fj_quick_sort(A, l, p-1);
            }


            #pragma omp task untied
            {
                fj_quick_sort(A, p+1, r);
            }
        }
       
    }


}


void fj_quick_sort_wrapper(int *A,int n,int threads) {
    // Currently just calls the standard quick_sort algorithm.
    // Replace this with your algorithm. It will probably
    // be best to implement methods called better_quick_sort
    // and better_partition and have this algorithm call
    // better_quick_sort and have that use better_partition.


    omp_set_num_threads(threads);
    omp_set_nested(true);
    threadss = threads;
    #pragma omp parallel
    {
        #pragma omp single
        {
            fj_quick_sort(A, 0, n-1);
        }
           
    }
}
