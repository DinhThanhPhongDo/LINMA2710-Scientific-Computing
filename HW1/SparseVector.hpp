#ifndef SPARSEVECTORHEADERDEF
#define SPARSEVECTORHEADERDEF
#include "AbstractVector.hpp"
#include "Vector.hpp"

class SparseMatrix;

class SparseVector : public AbstractVector
{
    private:
        int nnz;
        int *rowidx;
        double *nzval;
    public:
        SparseVector(int nnz, int const *rowidx, double const *nzval, int size);
        SparseVector();
        ~SparseVector();
        // read-only zero-based indexing
        double Read(int i) const;
        // assignment
        SparseVector& operator=(const SparseVector& otherVector);
        // binary +
        SparseVector operator+(const SparseVector& v1) const;
};
#endif
