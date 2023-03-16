#ifndef ABSTRACTVECTORHEADERDEF
#define ABSTRACTVECTORHEADERDEF
class AbstractVector
{
    protected:
        int mSize; // size of vector
    public:
        AbstractVector(int size);
        virtual ~AbstractVector();
        int GetSize() const;
        // read-only zero-based indexing
        virtual double Read(int i) const = 0;
        // declare length function as a friend
        friend int length(const AbstractVector& v);
};
// Prototype signature of length() friend function
int length(const AbstractVector& v);
#endif
