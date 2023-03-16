
#include <cmath>
#include <iostream>
#include <cassert>
#include "SparseVector.hpp"

using namespace std;


SparseVector::SparseVector() : AbstractVector(1)
{
    nnz=0;
    nzval=NULL;
    rowidx=NULL;
}

SparseVector::SparseVector(int nnz, int const *rowidx, double const *nzval, int size) : AbstractVector(size)
{
    this->nnz=nnz;
    this->nzval=new double[nnz];
    this->rowidx=new int[nnz];
  
    for(int j=0; j<nnz; j++){
        this->rowidx[j]=rowidx[j]+1;
        this->nzval[j]=nzval[j];

    }
}

// Overridden destructor to correctly free memory
SparseVector::~SparseVector()
{
    delete[] nzval;
    delete[] rowidx;
}

double SparseVector::Read(int i) const
{
    assert(i > -1);
    //printf("GetSize() = %d \n",GetSize());
    assert(i < GetSize());
    for (int j=0; j<this->nnz; j++) {
        if (i==this->rowidx[j]-1) {
            return this->nzval[j];
        }
        
    }
    return 0.0;
}


// Overloading the assignment operator

SparseVector& SparseVector::operator=(const SparseVector& otherVector)
{
    delete [] rowidx;
    delete [] nzval;

    this->mSize = otherVector.mSize;

    this->nnz = otherVector.nnz;
    //printf("nnz = %d \n",this->nnz);
    this->rowidx = new int [otherVector.GetSize()];
    this->nzval = new double [otherVector.GetSize()];

    for (int i=0;i< this->nnz ;i++)
    {
        this->nzval[i] = otherVector.nzval[i];
        this->rowidx[i] = otherVector.rowidx[i]; 
        /*
        cout<<"rowidx[i]=";
        cout<<this->rowidx[i]<<endl;
        cout<<"nzval[i]=";
        cout<<this->nzval[i]<<endl;
        cout<<"--------------"<<endl;*/
    }

    return *this;
}


// source : https://www.geeksforgeeks.org/operations-sparse-matrices/

SparseVector SparseVector::operator+(const SparseVector& v1) const
{
    int max=nnz+v1.nnz; //sum of both nonzero
    int count=0;
    int count1=0;
    int temp;
    int bv=0;
    int bv1=0;
    // find the number of non zero (given by max)
    for(int i=0;(count<nnz && count1<v1.nnz);i++)
    {
        // les if else sont ils utiles?? Oui_ //cout et //cout1 servent de suport pour avancer dans les vecteurs 
        if(this->rowidx[count]>v1.rowidx[count1])
        //
        {
            count1++;
            
        }
        else if (this->rowidx[count]<v1.rowidx[count1])
            count++;
        else //if they are equal
        {
            count1++;
            count++;
            max--;
        }
            
           
    }

    double *values = new double [max];
    int *rows = new int [max];
    count=0;
    count1=0;

    
    for (int i = 0; i <max; i++)
    {
        //cout<< "i=";
        //cout<< i << endl;
        
        if(rowidx[count]==v1.rowidx[count1] )
        //adding togheter
        {
            rows[i]=rowidx[count]-1;
            values[i] = nzval[count]+ v1.nzval[count1];
            //cout<< values[i] << endl;
            //cout<< rows[i] << endl;
            count1++;
            count++;
            
        }
        else if (rowidx[count]>v1.rowidx[count1])
        // adding this row
        {
            rows[i]=v1.rowidx[count1]-1;
            values[i] = v1.nzval[count1];
            //cout<< values[i] << endl;
            //cout<< rows[i] << endl;
            count1++;
        }
        else
        //adding other row
        {
            rows[i]=rowidx[count]-1;
            values[i] = nzval[count];
            //cout<< values[i] << endl;
            //cout<< rows[i] << endl;
            count++;
        }
        
        if(count==nnz){
            bv=1;
            temp=i+1;
            break;
        }
        if(count1==v1.nnz){
            bv1=1;
            temp=i+1;
            break;
        }
        
    
    }
    if(bv==1){
        for(int j=temp;j<max;j++){
            rows[j]=v1.rowidx[count1]-1;
            values[j] = v1.nzval[count1];
            //cout<< values[j] << endl;
            //cout<< rows[j] << endl;
            count1++;
        }
        
        
    }
    if(bv1==1){
        for(int j=temp;j<max;j++){
            rows[j]=rowidx[count]-1;
            values[j] = nzval[count];
            //cout<< values[j] << endl;
            //cout<< rows[j] << endl;
            count++;
        }
        
        
    }


    
    SparseVector res(max,rows,values,GetSize());
    delete [] rows;
    delete [] values;
    return res;
}


