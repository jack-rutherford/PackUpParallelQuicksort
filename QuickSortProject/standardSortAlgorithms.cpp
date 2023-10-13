//-------------------------------------------------------------------------
#include <cstdlib>
#include "standardSortAlgorithms.h"
#include "util.h"
//-------------------------------------------------------------------------
//Written by Chuck Cusack, 2003.
//Modified by Chuck Cusack, 2008
//-------------------------------------------------------------------------
void insertion_sort(int *A,int l,int r) {
        int temp;
        int i;
        for (int j=l+1;j<=r;j++) {
            temp=A[j];
            i=j-1;
            while (i>=l && temp < A[i]) {
                   A[i+1]=A[i];
                   i--;
                   }
            A[i+1]=temp;
            }
}
//-------------------------------------------------------------------------
int partition(int *A, int l, int r) {
   // Easiest way to avoid worst-case behavior
    Swap(A[l],A[(l+r)/2]);
    int p;
    p = A[l];
    int i = l+1;
    int j = r;
    while (1) {
        while (A[i] <= p && i < r) ++i;
        while (A[j] >= p && j > l) --j;
        if (i >= j) {
           Swap(A[j],A[l]);
           return j;
           }
        else {
             Swap(A[i],A[j]); 
             }
        }
}
//-------------------------------------------------------------------------
void quick_sort(int *A,int l,int r) {
    if(l<r) {
       int p=partition(A,l,r);
       quick_sort(A,l,p-1);
       quick_sort(A,p+1,r);
   }
}
//-------------------------------------------------------------------------
void merge(int *A,int l,int m,int r) {
    int size=r-l+1; 
    int mid=m-l+1;
    int *B=new int[size];
    for(int i=0;i<mid;i++)
       B[i]=A[l+i];
    int blah=r+mid;
    for(int j=mid;j<size;j++)
       B[j]=A[blah-j];
    int i=0;
    int j=size-1;
    for(int k=l;k<=r;k++) {
       if(B[i]<B[j]) {
         A[k]=B[i];
         i++;
       }
       else {
         A[k]=B[j];
         j--;
       }
    }
    delete []B;
}
//-------------------------------------------------------------------------
void merge_sort(int *A,int l,int r) {
    if(l<r) {
        int mid = (l+r)/2;
        merge_sort(A, l, mid);
        merge_sort(A, mid + 1, r);
        merge(A, l, mid, r);
    }
}
//-------------------------------------------------------------------------
//
// threads is ignored for these algorithms
void insertion_sort_wrapper(int *A,int n, int threads) {
    insertion_sort(A,0,n-1);
}
void quick_sort_wrapper(int *A,int n, int threads) {
    quick_sort(A,0,n-1);
}
void merge_sort_wrapper(int *A,int n, int threads) {
    merge_sort(A,0,n-1);
}
void built_in_quick_sort_wrapper(int *A,int n, int threads) {
    qsort((void *)A,n,sizeof(int),compare);
}
//-------------------------------------------------------------------------
