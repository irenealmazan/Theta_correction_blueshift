#!/bin/bash

FILEROOT=$1
MATLABROOT="/Applications/MATLAB_R2008a"
PPHASEFOLDER="/Users/b66494/PythonPhasing"
SP4UTILFOLDER="$PPHASEFOLDER/Sp4utils"
FFTWHEADER="/Users/b66494/local/include"
FFTWLIB="/Users/b66494/local/lib"

rm $FILEROOT.mexmaci

gcc $FILEROOT.c Sp4utils/Sp4*.c -I$MATLABROOT/extern/include -I$SP4UTILFOLDER -I$PPHASEFOLDER -I$FFTWHEADER -c

$MATLABROOT/bin/mex CC=gcc LD=gcc -L $FFTWLIB/libfftw3.a -lm -output $FILEROOT $FILEROOT.o Sp4*.o

rm *.o