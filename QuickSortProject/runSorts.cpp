// runSorts.cpp
// A program that runs 5 different sorting algorithms on various data and
// outputs data that allows the algorithms to be compared.
//
// Your goal is to modify the quick sort and merge sort implementations
// in yourSorts.cpp to make them as efficient as possible.
//
// YOU SHOULD NOT MOFIDY THIS FILE AT ALL except commenting/uncommenting
// some of the Algorithms below (so it only runs the ones you want).
//
// Charles A. Cusack, September, 2003
// Revised, February, 2008
// Revised, August, 2015
//-------------------------------------------------------------------------
#include<iostream>
#include<fstream>
using namespace std;
#include<chrono>
#include<math.h>
#include "standardSortAlgorithms.h"
#include "betterQuick.h"
#include "betterMerge.h"
#include "PQuick.h"
#include "FJQuick.h"
#include "listGenerator.h"
#include "util.h"

struct Algorithm {
	string name;
	void(*alg)(int*,int,int);
	int threads;
	Algorithm(string n, void(*f)(int*,int,int), int t) {
		name = n;
		alg=f;
		threads=t;
	}
};

// You may comment out the unimplemented algorithms below until you implement
// them to save time when running your tests.  Make sure to change numAlgs to
// match the correct size.
const int numAlgs= 19;

Algorithm *algs[numAlgs] = {
	// Regular merge sort
	new Algorithm("Merge     ",merge_sort_wrapper,1),
	// The improved version of merge sort.
	new Algorithm("Merge2    ",better_merge_sort_wrapper,1),
	// the built-in qsort
	new Algorithm("qsort     ",built_in_quick_sort_wrapper,1),
	// regular quick sort
	new Algorithm("Quick     ",quick_sort_wrapper,1),
	// The improved version quick sort
	new Algorithm("Quick2    ",better_quick_sort_wrapper,1),
	// The parallel version quick sort on different numbers of threads
	new Algorithm("pQuickP1  ",p_quick_sort_wrapper,1),
	new Algorithm("pQuickP2  ",p_quick_sort_wrapper,2),
	new Algorithm("pQuickP4  ",p_quick_sort_wrapper,4),
	new Algorithm("pQuickP8  ",p_quick_sort_wrapper,8),
	new Algorithm("pQuickP16 ",p_quick_sort_wrapper,16),
	new Algorithm("pQuickP24 ",p_quick_sort_wrapper,24),
	new Algorithm("pQuickP32 ",p_quick_sort_wrapper,32),
	// The fork-join version of quick sort on different numbers of threads
	new Algorithm("fQuickP1  ",fj_quick_sort_wrapper,1),
	new Algorithm("fQuickP2  ",fj_quick_sort_wrapper,2),
	new Algorithm("fQuickP4  ",fj_quick_sort_wrapper,4),
	new Algorithm("fQuickP8  ",fj_quick_sort_wrapper,8),
	new Algorithm("fQuickP16 ",fj_quick_sort_wrapper,16),
	new Algorithm("fQuickP24 ",fj_quick_sort_wrapper,24),
	new Algorithm("fQuickP32 ",fj_quick_sort_wrapper,32)
};

//-------------------------------------------------------------------------
void doSort(int k,int *A,int n);
void outputAlgorithmName(int i);
//-------------------------------------------------------------------------

int main (int argc,char* argv[]) {
typedef std::chrono::high_resolution_clock Clock;


if(argc!=3 ) {       // Check for correct number of command line arguments.
    cout<<"USAGE: "<<argv[0]<<" dataset seed\n"; 
    cout<<"where dataset is a specification file for tests to run\n";
    return(1); 
}

//------------------------------------------------------------------------
// Seed the random number generator
int seed=atoi(argv[2]);
//------------------------------------------------------------------------
// Open the test cases input file
ifstream datafile(argv[1]);
int tests;
datafile>>tests;

//------------------------------------------------------------------------
// Set up the acccounting variables

int inputTypes=5;

// Counts for the sorts
int typeCount[inputTypes+1];
double costs[numAlgs][inputTypes+1];

for(int i=0;i<=inputTypes;i++) {
   typeCount[i]=0;
   }
for(int i=0;i<numAlgs;i++) {
   for(int j=0;j<=inputTypes;j++) {
      costs[i][j]=0;
   }
}

//------------------------------------------------------------------------
// The variable used by the list generator functions
char type;   // data type
int intType;   // data type
int mod;     // modulus for algorithms with mod.
int shift;   //end of clump
int percent; // clump percentage
int start;   // start of clump
int end;     //end of clump


//------------------------------------------------------------------------
// The variables used to store the arrays, etc.
int n;            //number of elements to sort
int* A=NULL;  // The array of elements to sort
int* B=NULL;  // The array the algorithms will sort
int* C=NULL;  // The array that will store the correct sorted array

//------------------------------------------------------------------------
// Print the header for the output of each trial.
cout<<"\n\nResults for test cases from file "<<argv[1]<<"\n\n";
cout<<"Algorithm Type Size  CPUTime (microseconds)\n";
cout<<"------------------------------\n";

//-----------------------------------
// Runs all of the trials
for(int i=0;i<tests;i++) {
     datafile>>type;          // Read the line from the file to determine
     datafile>>n;             // Type and size of input
     if(type=='m' || type=='c') { // Get additional information if needed
        datafile>>mod;
        datafile>>shift;
     }
     if(type=='c') {
        datafile>>percent;
        datafile>>start;
        datafile>>end;
     }
     //------------------------------------------------------------------------
     // Generate a random list to test with
     srand(seed);  // Seed the generator with the previously saved seed.
                  // See below why this is done.
     switch(type) {
     case 'i':  
         A=inorder(n);
         intType=0;
         break;
     case 'I':  
         A=almostInorder(n);
         intType=0;
         break;
     case 'o':  
         A=outorder(n);
         intType=1;
         break;
     case 'O':  
         A=almostOutorder(n);
         intType=1;
         break;
     case 'r':  
         A=randomlist(n);
         intType=2;
         break;
     case 'm':  
         A=randomlistMod(n,mod,shift);
         intType=3;
         break;
     case 'c':  
         A=clump(n,mod,shift,percent,start,end);
         intType=4;
         break;
     default:  
         A=randomlist(n);
         intType=2;
     }
     seed=rand(); // Store the seed so we can re-seed the generator with the
                  // same seeds.  This is used in case someone uses rand() in
                  // their sorting algorithm, which would make the lists 
                  // generated for them different than those who do not use rand.
     typeCount[intType]++;
     typeCount[inputTypes]++;

     //------------------------------------------------------------------------
     // Sort the list once with known correct algorithm.
     // Use it later to check that the others are correct.
     if(C!=NULL) {
         delete []C;
     }
     C=new int[n];
     for(int l=0;l<n;l++) {
         C[l]=A[l];
     }
     // I am using the UNIX qsort to generate the correct answer for each array.
     // You are certainly NOT allowed to call this or any other utility sort
     // algorithm.  I am using this because, unlike your algorithms, I know
     // that it actually works.
     qsort((void *)C,n,sizeof(int),compare);

     for(int i=0;i<numAlgs;i++) {
          //--------------------------------------------------------------------
          // Copy the array so it can be sorted multiple times
          if(B!=NULL) delete []B;
          B=new int[n];       // Copy data so it can be used later.
          for(int k=0;k<n;k++)  
             B[k]=A[k];

          //--------------------------------------------------------------------
          // Start the timers, run the sort, stop the timers
          auto startTime = Clock::now();

          doSort(i,B,n); // Perform the ith sort.

          auto endTime =Clock::now();

          //--------------------------------------------------------------------
          // Check that the algorithm actually sorted the array correctly
          int k=0;
          while(k<n && B[k]==C[k])  {
               k++;
          }
          if(k!=n)  {
             cout<<"Terminating.  This algorithm did not work: "<<(*algs[i]).name<<"\n";
             return(1);
          }

          //--------------------------------------------------------------------
          // Compute and store the running time of the algorithm.

          auto elapsedTime =std::chrono::duration_cast<std::chrono::nanoseconds>(endTime - startTime).count();
          long CPUTime=(long)((elapsedTime*1.0)/1000);

          costs[i][intType]+=CPUTime;
          costs[i][inputTypes]+=CPUTime;

          //--------------------------------------------------------------------
          // Output the details of this trial

          outputAlgorithmName(i);
          cout.width(2);
          cout<<type;
          cout.width(8);
          cout<<n;
          cout.width(10);
          cout<<fixed<<CPUTime;
          cout<<"\n"<<flush;

          if(B!=NULL) {
              delete []B;                      // Get rid of the sorted array.
              B=NULL;
              }
     }
     cout<<"------------------------------\n";
     if(A!=NULL) {
        delete []A;                      // Get rid of the sorted array.
        A=NULL;
     }
     if(C!=NULL) {
        delete []C;                      // Get rid of the sorted array.
        C=NULL;
     }
   }

    //-------------------------------------------------------------------------
    // Print out the results
    cout<<"\n\n";
    cout<<"           ";
    cout<<"In Order      Reversed        Random      Repeated";
    cout<<"       Clumped       Average\n";
    for(int i=0;i<numAlgs;i++) {
       outputAlgorithmName(i);
       for(int j=0;j<=inputTypes;j++) {
       cout.width(10);
       if(costs[i][j]!=0)
            cout<<(long)(costs[i][j]/typeCount[j])<<" ("<<typeCount[j]<<")";
       else
         cout<<"No Data";
    }
    cout<<"\n";
    //-------------------------------------------------------------------------
}

datafile.close();

return 0;
}
//-------------------------------------------------------------------------
void doSort(int k,int *A,int n) {
	Algorithm *alg = algs[k];
	(*alg->alg)(A,n,alg->threads);
}
//-------------------------------------------------------------------------
void outputAlgorithmName(int i) {
    if(i>=0 && i<numAlgs) {
       cout<<(*algs[i]).name;
    } else {
       cout<<"Invalid  ";
    }
}
//-------------------------------------------------------------------------
