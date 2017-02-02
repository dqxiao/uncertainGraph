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
#include <armadillo>
#include <string>
using namespace std;
using namespace arma;

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
    void vertexCommon(vector<double> & vs_comm, int upper,double & vc_sum); //Need implement one simple idea to estimate the node's density
    void approximateUniquness(vector<double> & vs_uniq, int k, int scaleUpper); // done
    
    void obfuscation(string obfMethod);
    void setPrivacyConstraint(int k, float tolerance);  // set the desirable privacy level
    void configuration(bool g_sc, bool opt);            // configuration
    void setColdStart(float sizeMutipler, float sigma, int nAttempt); // done
    void genObfuscation(double & robf); // genObfuscation based on your configuration
    void excludingNodes(vector<double> & vProbs); // function used to excluded the largests values
    
private:
    vector<double> gdd; // global degree distribution
    int k;  // obfuscation level
    float tolerance; // tolerance level
    bool g_sc; // used graph size constraint clustering
    bool opt;  // used randomized method for searching the approximate optimal edge modification
    float c;   // sizeMultipler
    float sigma; // sigma
    int t;        // the number of attempt
    vector<double> vs_uniq;
    vector<double> vProbs;
};



#endif /* UncertainGraph_hpp */

