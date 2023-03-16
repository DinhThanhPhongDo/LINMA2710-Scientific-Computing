#include <cmath>
#include <iostream>
#include <cassert>
#include <array>
#include "SparseVector.hpp"
#include "Vector.hpp"
#include "SparseMatrix.hpp"
#include "Vector.hpp"

using namespace std;

SparseMatrix::SparseMatrix(int numRows, int numCols){
    // not too sure it is exact...
    
    this->m=numRows;        // Number of rows (size1?)
    this->n=numCols;        // Number of columns (size2?)
    this->M=0;              // the maximal number of non-zero entries on one column in the matrix
    this->nnz=0;            // number of non-zero values
    this->rowidx=NULL;      // row indices
    this->nzval=NULL;       // Stored values
}
SparseMatrix::SparseMatrix(int m, int n, int nnz, int M, int* rowidx, double* nzval){
    //logiquement, il est bon car c'est un deep copy
    this->m=m;
    this->n=n;
    this->M=M;
    this->nnz = nnz;

    int len = this->n * this->M; //c'est normalement bon. meme si on a une colonne de 0, n repprésente toujours nombre de col. Et on aura une ligne de -1

    this->rowidx = new int [len];
    for (int i = 0; i < len; i++)
    {
        this->rowidx[i]=rowidx[i];
    }   

    this->nzval=new double [len];
    for (int i = 0; i < len; i++)
    {
        this->nzval[i]=nzval[i];
    }
    
    
}
SparseMatrix::SparseMatrix(int nnz, int const *ridx, int const *cidx, double const *nzval, int size1, int size2){
    this->m=size1;
    this->n=size2;
    this->nnz = nnz;
    this->M=0;// compter les columns et trouver le max ( A FAIRE)

    for (int j=0; j<this->n ; j++)
    {
        int count =0;

        for (int i =0; i< this->nnz; i++){

            if (cidx[i]==j){

                count ++;
            }

        }

        if (count > this->M){

            this->M = count;
        }

    }


    
    int len = this->n * this->M;

    this->rowidx = new int [len];
    this->nzval=new double [len];

    for (int i = 0; i < len; i++)
    {
        this->rowidx[i]=-1;
        this->nzval[i]=0.;
    }   

    for (int j=0; j<this->n;j++){
        int count=0;
        for (int i =0; i< this->nnz; i++){

            if (cidx[i]==j){
                /*
                cout<< "ridx[i]=";
                cout<< ridx[i]<< endl;  
                cout<< "cidx[i]=";
                cout<< cidx[i]<< endl;   

                cout<<"i=";
                cout<<this->M*j+count<<endl;

                */
                this->rowidx[this->M*j+count]=ridx[i];
                this->nzval[this->M*j+count]=nzval[i];
                /*
                cout<< "rowidx=";
                cout<< rowidx[this->M*j+count] << endl;
                cout<< "val=";
                cout<< nzval[this->M*j+count]<< endl;
                cout<<"-------------------------"<< endl;
                */
                count=count+1;
                
            }

        }
    }
    /*
    cout<< this->nnz << endl;
    cout<< this->m << endl;
    cout<< this->n << endl;
    cout<< this->M << endl;
    */



}
SparseMatrix::~SparseMatrix(){
    delete[] rowidx;
    delete[] nzval;
}
int SparseMatrix::GetSize(int i) const{
    if (i==1){
        return this-> m;
    }
    else{
        return this->n;
    }
}
SparseMatrix& SparseMatrix::operator=(const SparseMatrix& otherSparseMatrix){
   
    delete[] rowidx;
    delete[] nzval;
    
    this->m=otherSparseMatrix.m;        // Number of rows (size1?)
    this->n=otherSparseMatrix.n;        // Number of columns (size2?)
    this->M=otherSparseMatrix.M;              // the maximal number of non-zero entries on one column in the matrix
    this->nnz=otherSparseMatrix.nnz;
    int len = this->n * this->M;

    this->rowidx = new int [len];
    this->nzval = new double[len];
    

    for (int i = 0; i < len; i++)
    {
        this->rowidx[i]=otherSparseMatrix.rowidx[i];
        this->nzval[i]=otherSparseMatrix.nzval[i];
    }   
    
    return *this;
}
SparseMatrix SparseMatrix::operator+() const{
    int len = this->n * this->M;

    int *new_rowidx = new int [len];
    double *new_nzval = new double[len];
    
    for (int i = 0; i < len; i++)
    {
        new_rowidx[i] = this->rowidx[i];
        new_nzval[i] = this->nzval[i];
    }   
    
    SparseMatrix new_spm(this->m, this->n, this->nnz, this->M, new_rowidx, new_nzval);
    
    delete [] new_rowidx;
    delete [] new_nzval;
    return new_spm;
}
SparseMatrix SparseMatrix::operator-() const{
    int len = this->n * this->M;

    int *new_rowidx = new int [len];
    double *new_nzval = new double[len];
    
    for (int i = 0; i < len; i++)
    {
        new_rowidx[i] = this->rowidx[i];
        new_nzval[i] = this->nzval[i]*(-1.0);
    }
    
    SparseMatrix new_spm(this->m, this->n, this->nnz, this->M, new_rowidx, new_nzval);
    
    delete [] new_rowidx;
    delete [] new_nzval;
    return new_spm;
}
SparseMatrix SparseMatrix::operator*(double a) const{
    int len = this->n * this->M;

    int *new_rowidx = new int [len];
    double *new_nzval = new double[len];
    
    for (int i = 0; i < len; i++)
    {
        new_rowidx[i] = this->rowidx[i];
        new_nzval[i] = this->nzval[i]*a;
    }   
    
    SparseMatrix new_spm(this->m, this->n, this->nnz, this->M, new_rowidx, new_nzval);
    
    delete [] new_rowidx;
    delete [] new_nzval;
    return new_spm;
}

Vector operator*(const SparseMatrix& m, const Vector& v){

    
    int len = m.M * m.n;

    Vector ans(m.m);
    for (int i =0; i< len; i++){

        if (m.rowidx[i] != -1){
            int index = (int) i/m.M ;
            ans(m.rowidx[i]) += v.Read(index)* m.nzval[i];
        }
    }
    /*
    for (int i =0; i< m.m; i++){
        cout<< ans.Read(i) <<endl;
    }
    */
    return ans;
    
}