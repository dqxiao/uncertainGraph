//
//  Distribution.cpp
//  uncertainGraph
//
//  Created by dongqingxiao on 10/6/16.
//  Copyright © 2016 dongqingxiao. All rights reserved.
//

#include "Distribution.hpp"
#include <iostream>
#include <math.h>
#include "Debug.cpp"
using namespace std;


void Distribution::setModel(vector<double> probs){
    edgeProbs=probs;
}






void DiscreteDistribution::setModel(vector<double> probs){
    edgeProbs=probs;
    
    // generated obProbs for each possible value
    vector<double> obProbs;
    
    obProbs.push_back(1);
    int count;
    for(double item: probs){
       
        double preVal=0;
        int i;
        for(i=0;i<obProbs.size();i++){
          preVal=obProbs[i];
            if(preVal==0){
                obProbs[i]=obProbs[i]*(1-item);
            }else{
                obProbs[i]=preVal*item+ obProbs[i]*(1-item);
            }
            
        }
        
        obProbs[i]=preVal*item;
        //done
        
    }
    
    cout<<obProbs<<endl;
}


long DiscreteDistribution::domainRange(){
    return edgeProbs.size();
}



DiscreteDistribution::DiscreteDistribution(vector<double> probs){
    setModel(probs);
}

double DiscreteDistribution::beliefProb(DiscreteDistribution other){
//  return logProb
    long selfDomain=edgeProbs.size();
    if(other.domainRange()>=selfDomain){
        double sum=0;
        vector<double> bProbs=other.probs;
        for(int i=0;i<selfDomain;i++){
            sum+=probs[i]*log(bProbs[i]);
        }
        
        return sum;
    }
    
    return 0;
}

void ContinousDistribution::setModel(vector<double> probs){
//  normal distribution
    edgeProbs=probs;
    double probSum=0;
    double varSum=0;
    
    for(double item: probs){
        probSum+=item;
        varSum+=item*(1-item);
    }
    
    mean=probSum/probs.size();
    std=sqrt(varSum);
}



