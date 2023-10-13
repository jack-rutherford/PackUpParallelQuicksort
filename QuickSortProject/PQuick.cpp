#include<omp.h>

#include "PQuick.h"
#include "standardSortAlgorithms.h"
#include "betterQuick.h"
#include "util.h"
#include<iostream>
#include<fstream>
using namespace std;
int threadss = 0;


void myP_quick_sort(int *A, int l, int r)
{
	if(l<r)
	{
		if((r - l) < 50)
		{
			insertion_sort(A, l, r);
		}
		else if((r - l) < 1000)
		{
			int p = partition(A, l, r);
			myP_quick_sort(A,l, p-1);
			myP_quick_sort(A, p+1, r);
		}
		else{
			// cout << "Inside else statement\n";
			int p = partition(A, l, r);
			#pragma omp task untied
			{
				myP_quick_sort(A, l, p-1);
			}

			#pragma omp task untied
			{
				myP_quick_sort(A, p+1, r);
			}
		}
		
	}

}

//-------------------------------------------------------------------------
void p_quick_sort_wrapper(int *A,int n, int threads) {
	// Get the OMP stuff started and make an initial call
	// to your sorting algorithm, which I recommend you call
	// p_quick_sort
	
	// For now it just calls the better version of quick_sort
	
	omp_set_num_threads(threads);
	omp_set_nested(true);
	threadss = threads;
	#pragma omp parallel 
	{
		#pragma omp single
		{
			myP_quick_sort(A, 0, n-1);
		}
			
	}
}
//-------------------------------------------------------------------------

//------------------------------------------------------------------------