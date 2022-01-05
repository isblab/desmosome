Place the public header files in this directory. They will be
available to your code (and other modules) with

     #include <IMP/desmosome/myheader.h>

All headers should include `IMP/desmosome/desmosome_config.h` as their
first include and surround all code with `IMPDESMOSOME_BEGIN_NAMESPACE`
and `IMPDESMOSOME_END_NAMESPACE` to put it in the IMP::desmosome namespace
and manage compiler warnings.

Headers should also be exposed to SWIG in the `pyext/swig.i-in` file.
