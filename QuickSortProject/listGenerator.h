//--------------------------------------------------------------------
// listGenerator.h
// Generate various types of lists
// Chuck Cusack, September 2003
//--------------------------------------------------------------------
#include<cstdlib>
//--------------------------------------------------------------------
int *almostInorder(int size);
int *almostOutorder(int size);
int *inorder(int size);
int *outorder(int size);
int *randomlist(int size);
int *randomlistMod(int size,int mod,int shift);
int *clump(int size,int mod,int shift,int percent,int start,int end);
//--------------------------------------------------------------------

