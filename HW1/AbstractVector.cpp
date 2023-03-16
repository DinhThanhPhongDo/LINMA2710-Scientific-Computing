#include <cassert>
#include "AbstractVector.hpp"

AbstractVector::AbstractVector(int size)
{
    assert(size > 0);
    mSize = size;
}
// Method to get the size of a vector
int AbstractVector::GetSize() const
{
    return mSize;
}
// MATLAB style friend to get the size of a vector
int length(const AbstractVector& v)
{
    return v.mSize;
}
AbstractVector::~AbstractVector()
{
}
