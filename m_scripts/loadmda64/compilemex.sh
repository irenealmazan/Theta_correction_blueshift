#!/bin/bash

FILEROOT=mda2matlab
MATLABROOT='/Applications/MATLAB_R2015b.app'
MDAUTILSFOLDER='./mdautils_64bit'

gcc $FILEROOT.c -I$MATLABROOT/extern/include -I$MDAUTILSFOLDER -c
$MATLABROOT/bin/mex -largeArrayDims CC=gcc LD=gcc $FILEROOT.o -L$MDAUTILSFOLDER -lm -lmda-load -output $FILEROOT 
