//
//  util.cpp
//  uncertainGraph
//
//  Created by dongqingxiao on 1/18/17.
//  Copyright © 2017 dongqingxiao. All rights reserved.
//

#include "util.hpp"
using namespace arma;
using namespace std;



int buildSimilarityMatrix(fmat data, fmat& S)
{
    int Nc = data.n_rows;
    
    for(int i = 0;i<Nc;i++)
    {
        S(i,i)=0;
        for(int j = i+1;j<Nc;j++)
        {
            vec diff = conv_to<vec>::from(data.row(i)-data.row(j));
            vec powvec = pow(diff,2);
            float s = sum(powvec);
            S(i,j)=s;
            S(j,i)=s;
        }
    }
    
    
    float stdev = as_scalar(stddev(stddev(S,0,1),0,0));
    S = exp(-S/stdev);
    return 0;
}


int outputMatrix(fmat matrix, string matrix_filename)
{
    
    ofstream file;
    file.open(matrix_filename.c_str());
    for(int i = 0; i < matrix.n_rows; i++)
    {
        for(int j = 0; j < matrix.n_cols; j++)
        {
            
            file.width(10);
            file.setf(ios::left);
            file << fixed;
            file << setprecision(5) << matrix(i,j);
        }
        file << "\n";
    }
    file.close();
    return 0;
}


