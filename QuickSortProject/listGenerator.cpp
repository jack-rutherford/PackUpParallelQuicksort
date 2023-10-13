//--------------------------------------------------------------------
// listGenerator.cpp
// Generate various types of lists
// Chuck Cusack, September 2003
//--------------------------------------------------------------------
#include "listGenerator.h"
//--------------------------------------------------------------------
int random(int mod) {
//if(mod==0)
    //return rand()+rand()*32768;
//else
    //return (rand()+rand()*32768)%mod;
if(mod==0)
    return rand();
else
    return rand()%mod;
}
//--------------------------------------------------------------------
// Generate an array of size elements consisting of 0 through size-1,
//  in order, but then mess a few up.
int *almostInorder(int size) {
    int *A=new int[size];
    for(int i=0;i<size;i++) 
       A[i]=i;
    // Mess up a couple of elements.
    for(int j=0;j<10;j++) {
	    int i1=random(size);
	    int i2=random(size);
            int temp = A[i1];
            A[i1]=A[i2];
            A[i2]=temp;
    }
    return A;
}
//--------------------------------------------------------------------
// Generate an array of size elements consisting of 0 through size-1,
//  out of order, but then mess a few up.
int *almostOutorder(int size) {
    int *A=new int[size];
    for(int i=0;i<size;i++) 
       A[i]=size-i-1;
    // Mess up a couple of elements.
    for(int j=0;j<10;j++) {
	    int i1=random(size);
	    int i2=random(size);
            int temp = A[i1];
            A[i1]=A[i2];
            A[i2]=temp;
    }
    return A;
}
//--------------------------------------------------------------------
// Generate an array of size elements consisting of 0 through size-1,
//  in order
int *inorder(int size) {
    int *A=new int[size];
    for(int i=0;i<size;i++) 
       A[i]=i;
    return A;
}
//--------------------------------------------------------------------
// Generate an array of size elements consisting of 0 through size-1,
//  in reverse order
int *outorder(int size) {
    int *A=new int[size];
    for(int i=0;i<size;i++) 
       A[i]=size-i-1;
    return A;
}
//--------------------------------------------------------------------
// Generate a random list of numbers from 0 to 2^30.
int *randomlist(int size) {
    int *A=new int[size];
    for(int i=size-1;i>-1;i--) 
       A[i]=random(0);
    return A;
}
//--------------------------------------------------------------------
// Generate a random list of numbers from 0 to mod, shifted up by shift.
int *randomlistMod(int size,int mod,int shift) {
    int *A=new int[size];
    for(int i=size-1;i>-1;i--) 
       A[i]=shift + random(mod);
    return A;
}
//--------------------------------------------------------------------
// Create a list of size size where each entry is between 
// shift-1 and mod-shift-1,
// percent percent of the elements are between start and end.
//
int *clump(int size,int mod,int shift,int percent,int start,int end) {
    if(start<shift || start>shift+mod-1) 
          start=shift;               // make sure start is in range
    if(end>mod+shift-1 || end < shift)
          end=mod+shift-1;    // Make sure end is in range

    int *A=new int[size];

    for(int i=0;i<size;i++)  {
       int number=shift + random(mod);
       if(rand()%100<percent)
         number=number%(end-start)+start;
       A[i]=number;
    }
    return A;
}
//--------------------------------------------------------------------
