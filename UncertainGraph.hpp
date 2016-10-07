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
#include <vector>
#include <Distribution.hpp>
using namespace std;

class UncertainGraph{
public:
    igraph_t graph;
    long nv;
    long ne;
    igraph_vector_t pe;

public:
    UncertainGraph(long nv);      // empty constructor
    UncertainGraph(const UncertainGraph & obj); // copy const constructor
    UncertainGraph(UncertainGraph & obj); // copy constructor
    UncertainGraph& operator=(const UncertainGraph & obj); //
    ~UncertainGraph(); // destructor function
    
    void setEdges(igraph_vector_t *edges); // set edges
    void setEdgeProbs(igraph_vector_t *probs); // set edges Probs
    vector<double> edgeProbs(igraph_integer_t v);  // get its neighbood ...
    

//    just for I/O debugging
    vector<ContinousDistribution> degreeDistributions();
    
};



#endif /* UncertainGraph_hpp */

