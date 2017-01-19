//
//  UncertainGraph.hpp
//  uncertainGraph
//
//  Created by dongqingxiao on 10/6/16.
//  Copyright © 2016 dongqingxiao. All rights reserved.
//

#ifndef UncertainGraph_hpp
#define UncertainGraph_hpp

#include <stdio.h>
#include <igraph.h>
#include <iostream>
#include "Debug.hpp"
#include <vector>
using namespace std;

class UncertainGraph{
public:
    igraph_t graph;
    long nv;
    long ne;
    igraph_vector_t pe;
    UncertainGraph(long nv);      // empty constructor
    UncertainGraph(const UncertainGraph & obj); // copy const constructor
    UncertainGraph(UncertainGraph & obj); // copy constructor
    UncertainGraph& operator=(const UncertainGraph & obj); // equal copy constructor
    ~UncertainGraph(); // destructor function
    void setEdges(igraph_vector_t *edges); // set edges
    void setEdgeProbs(igraph_vector_t *probs); // set edges Probs
    vector<double> edgeProbs(igraph_integer_t v);  // get its neighbood ...
    vector<DistPrameter> degreeDistributions(); // get approximate normal distribution
    vector<DistPrameter> degreeDistributions(vector<int> vs); // get the approximate normal distribution of the target vertex
    void globalDegreeDist(vector<double> &result); // calculation
    void approxGlobalDegreeDist(vector<double> &result); // cal the degree dist of the overall uncertain graph 
    int getMaxProbDegree(); //get max Prob Degree
    void approxVertexSim(vector<double>& vs,vector<double>&result); // using sampling to calculate in the sample way
    vector<DistPrameter> degreeDistributions(double lb);
    void anonymityCheck(UncertainGraph obGraph, int k, double * r); // check whether the obfuscated graph provide enough privacy level

private:
    vector<double> gdd;

    
};



#endif /* UncertainGraph_hpp */

