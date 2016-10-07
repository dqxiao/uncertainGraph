//
//  Distribution.hpp
//  uncertainGraph
//
//  Created by dongqingxiao on 10/6/16.
//  Copyright © 2016 dongqingxiao. All rights reserved.
//

#ifndef Distribution_hpp
#define Distribution_hpp

#include <stdio.h>
#include <vector>

using namespace std;



class Distribution{

public:
//    Distribution(vector<double> probs);
    void setModel(vector<double> probs);
    double beliefProb(Distribution other);
protected:
    vector<double> edgeProbs;
};


class DiscreteDistribution:Distribution {
    
public:
    DiscreteDistribution(vector<double> probs);
    void setModel(vector<double> probs);
    double beliefProb(Distribution other);
    long domainRange();
protected:
    vector<double> probs;
private:
    double beliefProb(DiscreteDistribution other);
};


class ContinousDistribution:Distribution{
//    here, let us consider N(mean, std^2) distribution class
public:
    ContinousDistribution(vector<double> probs);
    void setModel(vector<double> probs);
    double beliefProb(Distribution other);
protected:
    double mean;
    double std;
};
#endif /* Distribution_hpp */
