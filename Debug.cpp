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


//done
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


