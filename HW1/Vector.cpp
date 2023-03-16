#include <cmath>
#include <iostream>
#include <cassert>
#include "Vector.hpp"

// Overridden copy constructor
// Allocates memory for new vector, and copies
// entries of other vector into it
Vector::Vector(const Vector& otherVector) : AbstractVector(otherVector.GetSize())
{
    mData = new double [GetSize()];
    for (int i=0; i<GetSize(); i++)
    {
        mData[i] = otherVector.mData[i];
    }
}
// Constructor for vector of a given size
// Allocates memory, and initialises entries
// to zero
Vector::Vector(int size) : AbstractVector(size)
{
    mData = new double [GetSize()];
    for (int i=0; i<GetSize(); i++)
    {
        mData[i] = 0.0;
    }
}
// Overridden destructor to correctly free memory
Vector::~Vector()
{
    delete[] mData;
}
// Overloading square brackets
// Note that this uses ‘zero-based’ indexing,
// and a check on the validity of the index
double& Vector::operator[](int i)
{
    assert(i > -1);
    assert(i < GetSize());
    return mData[i];
}
// Read-only variant of []
// Note that this uses ‘zero-based’ indexing,
// and a check on the validity of the index
double Vector::Read(int i) const
{
    assert(i > -1);
    assert(i < GetSize());
    return mData[i];
}
// Overloading round brackets
// Note that this uses ‘zero-based’ indexing,
// and a check on the validity of the index
double& Vector::operator()(int i)
{
    assert(i > -1);
    assert(i < GetSize());
    return mData[i];
}
// Overloading the assignment operator
Vector& Vector::operator=(const Vector& otherVector)
{
    if (GetSize() == otherVector.GetSize())
    {
        for (int i=0; i<GetSize(); i++)
            mData[i] = otherVector.mData[i];
    }else{
        delete[] mData;
        mSize = otherVector.GetSize();
        mData = new double[GetSize()];
        for (int i=0; i<GetSize(); i++)
            mData[i] = otherVector.mData[i];
    }

    return * this;
}
// Overloading the unary + operator
Vector Vector::operator+() const
{
    Vector v(GetSize());
    for (int i=0; i<GetSize(); i++)
    {
        v[i] = mData[i];
    }
    return v;
}
// Overloading the unary - operator
Vector Vector::operator-() const
{
    Vector v(GetSize());
    for (int i=0; i<GetSize(); i++)
    {
        v[i] = -mData[i];
    }
    return v;
}
// Overloading the binary + operator
Vector Vector::operator+(const Vector& v1) const
{
    assert(GetSize() == v1.GetSize());
    Vector v(GetSize());
    for (int i=0; i<GetSize(); i++)
    {
        v[i] = mData[i] + v1.mData[i];
    }
    return v;
}
// Overloading the binary - operator
Vector Vector::operator-(const Vector& v1) const
{
    assert(GetSize() == v1.GetSize());
    Vector v(GetSize());
    for (int i=0; i<GetSize(); i++)
    {
        v[i] = mData[i] - v1.mData[i];
    }
    return v;
}
// Overloading scalar multiplication
Vector Vector::operator *(double a) const
{
    Vector v(GetSize());
    for (int i=0; i<GetSize(); i++)
    {
        v[i] = a*mData[i];
    }
    return v;
}
// Method to calculate norm (with default value p=2)
// corresponding to the Euclidean norm
double Vector::CalculateNorm(int p) const
{
    double norm_val, sum = 0.0;
    for (int i=0; i<GetSize(); i++)
    {
        sum += pow(fabs(mData[i]), p);
    }
    norm_val = pow(sum, 1.0/((double)(p)));
    return norm_val;
}
