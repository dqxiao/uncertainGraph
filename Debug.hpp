//
//  Debug.hpp
//  uncertainGraph
//
//  Created by dongqingxiao on 10/6/16.
//  Copyright © 2016 dongqingxiao. All rights reserved.
//

#ifndef Debug_hpp
#define Debug_hpp

#include <stdio.h>
#include <vector>
#include <math.h>
using namespace std;

class DistPrameter{
public:
    double mean;
    double stddev;
    DistPrameter(double mean, double stddev);
    DistPrameter(vector<double> edgeProbs);
    DistPrameter();
    int probMaxVal();
    int probMinVal();
    double PDF(int val); // done
    double CDF(int val); // done
    void degreeDistribution(vector<double> & result);
    double distance(DistPrameter o);
    double distance(DistPrameter o,string func);
    double getVariance();
};

class ExactDist{
public:
    vector<double> ps;
    ExactDist(vector<double> edgeProbs);
    int probMaxVal();
    int probMinVal();
    void degreeDistribution(vector<double> & result);
};

#endif /* Debug_hpp */
