%MAKE_TETGEN_MEX Makes the tetgen interface
mex -DPTR_EXCEEDS_LONG -DTETLIBRARY -output tetgen_mex ...
    -I../../ThirdParty/tetgen1.4.3/source...
    ../../ThirdParty/tetgen1.4.3/source/predicates.cxx ...
    ../../ThirdParty/tetgen1.4.3/source/tetgen.cxx ...
    tetgen_mex.cpp 