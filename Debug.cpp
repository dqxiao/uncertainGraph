//
//  Debug.cpp
//  uncertainGraph
//
//  Created by dongqingxiao on 10/6/16.
//  Copyright © 2016 dongqingxiao. All rights reserved.
//

#include "Debug.hpp"
#include <iostream>
#include <vector>
#include <math.h>
#include <fstream>
#include <boost/math/distributions/normal.hpp>
using boost::math::normal;
using boost::math::iround;

std::ostream &operator<<(std::ostream &os, std::vector<double> &v)
{
    for(auto &i: v)
        os << i<< std::endl;
    return os;
}


DistPrameter::DistPrameter(double m, double std){
    mean=m;
    stddev=std;
}

DistPrameter::DistPrameter(){
    DistPrameter(0,0);
}


DistPrameter::DistPrameter(vector<double> eps){
    double pSum=0;
    double stdSum=0;
    for(double item: eps){
        pSum+=item;
        stdSum+=item*(1-item);
    }
    mean=pSum;
    stddev=sqrt(stdSum); //done
}


int DistPrameter::probMaxVal(){
    return ceil(mean+3*stddev);
}


int DistPrameter::probMinVal(){
    double minVal=mean-3*stddev;
    if (minVal>0){
        return floor(minVal);
    }
    return 0; //observation
}





double DistPrameter::CDF(int val){
    if(val>=mean+3*stddev){
        return 1;
    }
    if(val<mean-3*stddev){
        return 0;
    }
    
    if(stddev==0){
        if(abs(val-mean)<1){
            return 1;
        }
        return 0;
    }
    normal ns(mean,stddev);
    
    double result=boost::math::cdf(ns, val-0.5);
    
    return result;
}
double DistPrameter::PDF(int val){
    return CDF(val+1)-CDF(val);  // OK
}





void DistPrameter::degreeDistribution(vector<double> &result){
    int l=probMinVal();
    int h=probMaxVal();
    for(int val=l;val<=h;val++){
        result.push_back(PDF(val));
    }
    //done
}

// done
double DistPrameter::distance(DistPrameter o){
    double result=0;
    double rv=pow(stddev,2)/pow(o.stddev,2);
    result+=log(0.25*(rv+1.0/rv+2));
    result+=pow(abs(mean-o.mean),2)/(pow(stddev,2)+pow(o.stddev,2));
    result*=0.25;
    return result;
}


double DistPrameter::distance(DistPrameter o, string func){
    if(func=="EL"){
        return sqrt(pow(o.mean-mean,2)+pow(o.stddev-stddev, 2));
    }else{
        return distance(o);
    }
    //done 
}


double DistPrameter::similiarity(DistPrameter o){
    int th=probMaxVal();
    int tl=probMinVal();
    
    int oh=o.probMaxVal();
    int ol=o.probMinVal();
    
    int l=max(tl,ol);
    int h=min(th,oh);
    double result=0;
    if(h>=l){
        vector<double> tR;
        vector<double> oR;
        for(int i=l;i<=h;i++){
            result+=o.PDF(i)*PDF(i);
        }
    }
    return result;
}

double DistPrameter::similiarityNormal(DistPrameter o){
    double n_mean=mean-o.mean;
    double n_std=sqrt(pow(stddev,2)+pow(o.stddev,2));
    DistPrameter n(n_mean,n_std);
    return n.PDF(0);
}



double DistPrameter::getVariance(){
    return pow(stddev,2);
}

ExactDist::ExactDist(vector<double> probs){
    ps=probs;
}


int ExactDist::probMinVal(){
    return 0;
}

int ExactDist::probMaxVal(){
    return (int)ps.size();
}

// you need know more about boost library
void ExactDist::degreeDistribution(vector<double> & result){
    
    int l=probMinVal();
    int h=probMaxVal();
    
    vector<double> r (h-l+1,0);
    r[0]=1;
    
    for(int i=0;i<ps.size();i++){
        double p=ps[i];
        for(int j=i;j>=0;j--){
            r[j+1]=r[j+1]*(1-p)+r[j]*p;
        }
        r[0]=r[0]*(1-p);
    }
    
    result=r;
}


double entropy(double p){
    if(p==0){
        return 0;
    }
    
    return -p*log2(p);
}

// read dist parameter from file --just for debugging
void readDistParameter(string filepath, vector<DistPrameter> &vs){
    
    fstream distFile(filepath);
    string line;
    string item;
    vector<double> v;
    if(distFile.is_open()){
        while (getline(distFile,line)) {
            stringstream ss(line); // for split line
            int i=0;
            double mean;
            double std;
            while(getline(ss,item,'\t')){
                if(i==0){
                    mean=stof(item);
                }
                if(i==1){
                    std=stof(item);
                }
                i+=1;
            }
            vs.push_back(DistPrameter(mean, std));
        }
        
        
        distFile.close();
    }else{
        cout<<"can't find/open"<<filepath<<endl;
    }
    
}

//
void test_fast_discreteDistribution(const std::vector<double>& weights, const size_t num_samples){
    std::default_random_engine generator;
    fast_discrete_distribution<int> distribution(weights);
    distribution.PrintBuckets();
    
    std::vector<size_t> counts(weights.size(), 0);
    for (size_t i = 0; i < num_samples; ++i) {
        const int number = distribution(generator);
        assert(number >= 0);
        assert(number < static_cast<int>(weights.size()));
        ++counts[number];
    }
    
    std::cout << "counts:" << std::endl;
    for (size_t i = 0; i < weights.size(); ++i)
        cout << i << " (" << weights[i] << ") : "
        << std::string(counts[i], '*') << endl;
    
    cout << endl;

}

void writeUniquess(string filePath,vector<double> vs_uniq){
    ofstream myfile(filePath);
    for(double u: vs_uniq){
        myfile<<u<<endl;
    }
    myfile.close();
}


