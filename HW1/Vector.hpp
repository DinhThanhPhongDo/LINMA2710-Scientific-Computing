#ifndef VECTORHEADERDEF
#define VECTORHEADERDEF
#include "AbstractVector.hpp"
class Vector : public AbstractVector
{
    private:
        double * mData; // data stored in vector
    public:
        Vector(const Vector& otherVector);
        Vector(int size);
        ~Vector();
        double& operator[](int i); // zero-based indexing
        // read-only zero-based indexing
        double Read(int i) const;
        double& operator()(int i); // zeros-based indexing
        // assignment
        Vector& operator=(const Vector& otherVector);
        Vector operator+() const; // unary +
        Vector operator-() const; // unary -
        Vector operator+(const Vector& v1) const; // binary +
        Vector operator-(const Vector& v1) const; // binary -
        // scalar multiplication
        Vector operator *(double a) const;
        // p-norm method
        double CalculateNorm(int p=2) const;
};
#endif
