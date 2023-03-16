#include <cstdlib>
#include <cstdio>
#include <cmath>
#include "AbstractVector.hpp"
#include "SparseVector.hpp"
#include "Vector.hpp"
#include "SparseMatrix.hpp"

#define IWANTCOLOR // comment this line to remove color
#ifdef IWANTCOLOR
#define ANSI_COLOR_RED     "\x1b[31m"
#define ANSI_COLOR_GREEN   "\x1b[32m"
#define ANSI_COLOR_YELLOW  "\x1b[33m"
#define ANSI_COLOR_BLUE    "\x1b[34m"
#define ANSI_COLOR_MAGENTA "\x1b[35m"
#define ANSI_COLOR_CYAN    "\x1b[36m"
#define ANSI_COLOR_RESET   "\x1b[0m"
#else
#define ANSI_COLOR_RED     ""
#define ANSI_COLOR_GREEN   ""
#define ANSI_COLOR_YELLOW  ""
#define ANSI_COLOR_BLUE    ""
#define ANSI_COLOR_MAGENTA ""
#define ANSI_COLOR_CYAN    ""
#define ANSI_COLOR_RESET   ""
#endif

#include <iostream>
using namespace std;

bool isapprox(double x, double y, double rtol=1e-8, double atol=1e-8)
{

    return fabs(x - y) <= atol + rtol*fmax(fabs(x), fabs(y));
}

bool isapprox(const AbstractVector& x, const AbstractVector& y, double rtol=1e-8, double atol=1e-8)
{
    int n = length(x);
    
    if (length(y) != n) {
        /*cout<<"false returned"<<endl;
        printf("lenght(y)= %d \n",length(y));
        printf("n= %d \n",n);*/
        return false;
    }
    
    
    
    for (int i = 0; i < n; i++) {
        double a = x.Read(i);
        double b = y.Read(i);
        /*printf("a = %lf\n",a);
        printf("b = %lf\n",b);*/
        if (!isapprox(a, b, rtol, atol))
        {
            printf(ANSI_COLOR_RED " x[%d] = %lf != %lf = y[%d]\n" ANSI_COLOR_RESET, i, a, b, i);
            return false;
        }
       
    }
    return true;
}

void testisapprox(AbstractVector const &x, AbstractVector const &y, char const *s)
{
    //printf("inside testisapprox");
    printf(ANSI_COLOR_BLUE "%s" ANSI_COLOR_RESET ":", s);
    if (!isapprox(x, y)) {
        //printf("!isapprox(x, y); return exit 1 \n");
        exit(EXIT_FAILURE);
    } else {
        printf(ANSI_COLOR_GREEN " OK\n" ANSI_COLOR_RESET);
    }
}

void testisapprox(double x, double y, char const *s)
{
    printf(ANSI_COLOR_BLUE "%s" ANSI_COLOR_RESET ":", s);
    if (!isapprox(x, y)) {
        exit(EXIT_FAILURE);
    } else {
        printf(ANSI_COLOR_GREEN " OK\n" ANSI_COLOR_RESET);
    }
}

int main()
{
    /*
        SparseVector tests
    */
    // instantiation empty SparseVector
    {
    SparseVector sv(0,NULL,NULL,2);
    Vector ans(2);
    testisapprox(ans, sv, "Instantiate empty sparse vector ");
    }

    // instantiation full SparseVector
    {
    double data[3] = {1.,2.,3.};
    int ridx[3] = {0,1,2};
    SparseVector sv(3,ridx,data,3);
    Vector ans(3);
    ans(0) = 1.; ans(1) = 2.; ans(2) = 3.;
    testisapprox(ans, sv, "Instantiate full sparse vector ");
    }

    // instantiation partial SparseVector
    {
    double data[3] = {1.,2.,3.};
    int ridx[3] = {1,2,4};
    SparseVector sv(3,ridx,data,5);
    Vector ans(5);
    ans(1) = 1.; ans(2) = 2.; ans(4) = 3.;
    testisapprox(ans, sv, "Instantiate partial sparse vector ");
    }

    // checks data is well copied and delete operation (with valgrind)
    {
    double* data = new double[3];
    data[0] = 1.; data[1] = 2.; data[2] = 3.;
    int* ridx = new int[3]; ridx[0] = 0; ridx[1] = 1; ridx[2] = 2;
    SparseVector *sv = new SparseVector(3,ridx,data,3);
    data[0] = 2.;
    ridx[1] = 2;
    Vector ans(3);
    ans(0) = 1.; ans(1) = 2.; ans(2) = 3.;
    testisapprox(ans, *sv, "data copied at construction of sparse vector ");
    delete[] data;
    delete[] ridx;
    delete sv;
    }

    // assignation = for SparseVector test (with valgrind for memory copy)
    {
    double data[3] = {1.,2.,3.};
    int ridx[3] = {1,2,4};
    SparseVector sv(3,ridx,data,5);
    Vector ans(5);
    ans(1) = 1.; ans(2) = 2.; ans(4) = 3.;

    SparseVector sv_;
    sv_ = sv;
    testisapprox(ans, sv_, "assignement for sparse vectors ");
    }

    // sum of SparseVectors
    {
    Vector v1(5), v2(5);
    v1(0) = 1.0;
    v1(1) = 2.;
    v1(3) = -3.;
    
    v2(1) = 1.0;
    v2(2) = 2.0;
    v2(4) = 4.5;

    int rowvals1[] = {0,1,3};
    int rowvals2[] = {1,2,4};
    double nzvals1[] = {1.,2.,-3.};
    double nzvals2[] = {1.,2.,4.5};
    SparseVector sv1(3,rowvals1,nzvals1, 5);
    SparseVector sv2(3,rowvals2,nzvals2, 5);

    testisapprox(v1+v2, sv1+sv2, "sum SparseVector");
    }
    
    /*
        SparseMatrix tests note that they rely on a correct matrix,vector product
    */
    
    Vector e0(4); e0(0) = 1.;
    Vector e1(4); e1(1) = 1.;
    Vector e2(4); e2(2) = 1.;
    Vector e3(4); e3(3) = 1.;
    Vector null(4);
    // SparseMatrix initialization
    {
    SparseMatrix sm0(4,4);
    testisapprox(null,sm0*e0, "instanciation null sparse matrix 0");
    testisapprox(null,sm0*e3, "instanciation null sparse matrix 4");
    int nnz = 7;
    int ridx[nnz] = {0,2,0,1,3,2,0};
    int cidx[nnz] = {0,0,1,1,1,2,3};
    double nzval[nnz] = {2.1,1.,3.5,2.2,2.0,2.3,1.9};
    SparseMatrix sm(nnz, ridx, cidx, nzval, 4, 4);
    Vector ans0(4); ans0(0) = 2.1; ans0(2) = 1.;
    Vector ans1(4); ans1(0) = 3.5; ans1(1) = 2.2; ans1(3) = 2.;
    Vector ans2(4); ans2(2) = 2.3;
    Vector ans3(4); ans3(0) = 1.9;
    testisapprox(ans0,sm*e0, "instanciation sparse matrix 0");
    testisapprox(ans1,sm*e1, "instanciation sparse matrix 1");
    testisapprox(ans2,sm*e2, "instanciation sparse matrix 2");
    testisapprox(ans3,sm*e3, "instanciation sparse matrix 3");

    //with disorder in the argumentq
    int ri[7] = {2,3,2,0,1,0,0};
    int ci[7] = {0,1,2,3,1,1,0};
    double nzv[7] = {1.,2.,2.3,1.9,2.2,3.5,2.1};
    sm = SparseMatrix(nnz, ri, ci, nzv, 4, 4);
    testisapprox(ans0,sm*e0, "instanciation sparse matrix with disorder 0");
    testisapprox(ans1,sm*e1, "instanciation sparse matrix with disorder 1");
    testisapprox(ans2,sm*e2, "instanciation sparse matrix with disorder 2");
    testisapprox(ans3,sm*e3, "instanciation sparse matrix with disorder 3");
    }

    // = and delete operation (need valgrind)
    {
    
    int* ridx = new int[2]; ridx[0]=0;ridx[1]=1;
    int* cidx = new int[2]; cidx[0]=0; cidx[1]=1;
    double* nzval = new double[2]; nzval[0]=0.; nzval[1]=0.;
    
    SparseMatrix *sm = new SparseMatrix(2,ridx,cidx,nzval,3,3);
    SparseMatrix sm_(3,3);
    sm_ = *sm;
    
    delete[] ridx;
    delete[] cidx;
    delete[] nzval;
    delete sm;
    }

    // unary +,- and scalar multilication
    {
    int nnz = 7;
    int ridx[nnz] = {0,2,0,1,3,2,0};
    int cidx[nnz] = {0,0,1,1,1,2,3};
    double nzval[nnz] = {2.1,1.,3.5,2.2,2.0,2.3,1.9};
    SparseMatrix sm(nnz, ridx, cidx, nzval, 4, 4);
    
    Vector ans0(4); ans0(0) = 2.1; ans0(2) = 1.;
    Vector ans1(4); ans1(0) = 3.5; ans1(1) = 2.2; ans1(3) = 2.;
    Vector ans2(4); ans2(2) = 2.3;
    Vector ans3(4); ans3(0) = 1.9;
    testisapprox(ans0,(+sm)*e0, "unary + sparse matrix 0");
    testisapprox(ans1,(+sm)*e1, "unary + sparse matrix 1");
    testisapprox(ans2,(+sm)*e2, "unary + sparse matrix 2");
    testisapprox(ans3,(+sm)*e3, "unary + sparse matrix 3");

    testisapprox(-ans0,(-sm)*e0, "unary - sparse matrix 0");
    testisapprox(-ans1,(-sm)*e1, "unary - sparse matrix 1");
    testisapprox(-ans2,(-sm)*e2, "unary - sparse matrix 2");
    testisapprox(-ans3,(-sm)*e3, "unary - sparse matrix 3");

    testisapprox(ans0*2.,(sm*2.)*e0, "scalar product sparse matrix 0");
    testisapprox(ans1*2.,(sm*2.)*e1, "scalar product sparse matrix 1");
    testisapprox(ans2*2.,(sm*2.)*e2, "scalar product sparse matrix 2");
    testisapprox(ans3*2.,(sm*2.)*e3, "scalar product sparse matrix 3");
    }

    // matrix vector product
    {
    int nnz = 7;
    int ridx[nnz] = {0,2,0,1,3,2,0};
    int cidx[nnz] = {0,0,1,1,1,2,3};
    double nzval[nnz] = {2.1,1.,3.5,2.2,2.0,2.3,1.9};
    SparseMatrix sm(nnz, ridx, cidx, nzval, 4, 4);

    Vector v(4); v(0)=0.7; v(1)=-3.4; v(2)=-1.; v(3)=1.;
    Vector ans(4); ans(0)=-8.53;ans(1)=-7.48;ans(2)=-1.6;ans(3)=-6.8;
    testisapprox(ans,sm*v, "sparse matrix * vector");
    }
    
    // rectangular matrix - vector product
    {
    int nnz = 6;
    int ridx[nnz] = {0,2,0,1,2,0};
    int cidx[nnz] = {0,0,1,1,2,3};
    double nzval[nnz] = {2.1,1.,3.5,2.2,2.3,1.9};
    SparseMatrix sm(nnz, ridx, cidx, nzval, 3, 4);
    
    Vector v(4); v(0)=0.7; v(1)=-3.4; v(2)=-1.; v(3)=1.;
    Vector ans(3); ans(0)=-8.53;ans(1)=-7.48;ans(2)=-1.6;
    testisapprox(ans,sm*v, "rectangular sparse matrix * vector");
    
    }
    


    return 0;
}
