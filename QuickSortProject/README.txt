This is code that is set up to make it easy to compare different
versions of several popular sorting algorithms.
It is currently set up to run on Ubuntu.  Because of how the times
are computed it may or may not run on other OSes.

It runs two versions of mergesort and 3 versions of quicksort.
The idea is to make improvements to the second version of quicksort
and mergesort (quicksort2 and mergesort2, both available in the
yourSorts.cpp file) and compare them with the unmodified versions.
The third version of quicksort is qsort, the version of quicksort
available in the standard C++ library.

Use "make" to compile.
Run it using something like:

runSorts file 11

where "file" is the name of the file with test cases and
"11" is any number--it
is used as a seed in the random number generator.

There are three sample test case files.  The "fullRange" tests will
take too long to run with the current versions of quicksort, so to run
that both quicksort and quicksort2 need to be improved.  Or fake it
and call qsort instead of quicksort in runSorts.cpp.

// Chuck Cusack, 8/19/15