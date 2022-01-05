Place the private header files in this directory. They will be
available to your code with

     #include <IMP/desmosome/internal/myheader.h>

All headers should include `IMP/desmosome/desmosome_config.h` as their
first include and surround all code with `IMPDESMOSOME_BEGIN_INTERNAL_NAMESPACE`
and `IMPDESMOSOME_END_INTERNAL_NAMESPACE` to put it in the
IMP::desmosome::internal namespace and manage compiler warnings.
