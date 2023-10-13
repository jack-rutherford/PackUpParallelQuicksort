#include<iostream>
#include<fstream>
using namespace std;
#include "util.h"

void Swap(int &x,int &y) {
   int temp = x;
   x = y;
   y = temp;
}
// For use with qsort
int compare(const void *x, const void *y) {
    return (*(int *)x - *(int *) y);
}

// Misc useful stuff
void printArray(int *a,int n) {
	if(n<100) {
		for(int i=0;i<n;i++) {
			cout<<a[i]<<", ";
		}
		cout<<"\n";
	} else {
		cout<<"too long to print\n";
	}
}
void clearArray(int *a,int n) {
	for(int i=0;i<n;i++) {
		a[i]=0;
	}
	cout<<"\n";
}
int equals(int *a,int *b,int n) {
	for(int i=0;i<n;i++) {
		if(a[i]!=b[i]) return 0;
	}
	return 1;
}
