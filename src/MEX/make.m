mex -O 'CXXOPTIMFLAGS=-DNDEBUG -O3' -largeArrayDims ExactLinesearchFull.c COMPFLAGS="$COMPFLAGS -openmp" LINKFALGS="$LINKFALGS -openmp"

mex -O 'CXXOPTIMFLAGS=-DNDEBUG -O3' -largeArrayDims mult_sparse.c bcdic.c COMPFLAGS="$COMPFLAGS -openmp" LINKFALGS="$LINKFALGS -openmp"

mex -O 'CXXOPTIMFLAGS=-DNDEBUG -O3' -largeArrayDims softshrink.c bcdic.c COMPFLAGS="$COMPFLAGS -openmp" LINKFALGS="$LINKFALGS -openmp"

mex -O 'CXXOPTIMFLAGS=-DNDEBUG -O3' -largeArrayDims mult_GinvAi.c bcdic.c COMPFLAGS="$COMPFLAGS -openmp" LINKFALGS="$LINKFALGS -openmp"

mex -O 'CXXOPTIMFLAGS=-DNDEBUG -O3' -largeArrayDims calcParamsForQuadLinesearch_mex.c bcdic.c COMPFLAGS="$COMPFLAGS -openmp" LINKFALGS="$LINKFALGS -openmp"

mex -O 'CXXOPTIMFLAGS=-DNDEBUG -O3' -largeArrayDims mult_SparseTranspose_Dense.c bcdic.c COMPFLAGS="$COMPFLAGS -openmp" LINKFALGS="$LINKFALGS -openmp"

mex -O 'CXXOPTIMFLAGS=-DNDEBUG -O3' -largeArrayDims calcGradAndActiveSet.c libmwblas.lib COMPFLAGS="$COMPFLAGS -openmp" LINKFALGS="$LINKFALGS -openmp"

mex -O 'CXXOPTIMFLAGS=-DNDEBUG -O3' -largeArrayDims mulXXt_sparse.c libmwblas.lib

mex -O 'CXXOPTIMFLAGS=-DNDEBUG -O3' -largeArrayDims PCG_block_mex_BCDIC.c COMPFLAGS="$COMPFLAGS -openmp" LINKFALGS="$LINKFALGS -openmp"

mex -O 'CXXOPTIMFLAGS=-DNDEBUG -O3' -largeArrayDims PCG_block_mex_BCDIC_Hybrid.c  MAT_hybrid.c COMPFLAGS="$COMPFLAGS -openmp" LINKFALGS="$LINKFALGS -openmp"
