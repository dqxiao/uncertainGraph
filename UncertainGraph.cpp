//
//  UncertainGraph.cpp
//  uncertainGraph
//
//  Created by dongqingxiao on 10/6/16.
//  Copyright © 2016 dongqingxiao. All rights reserved.
//

#include "UncertainGraph.hpp"
#include "UncertainGraph.hpp"
using namespace std;



UncertainGraph::UncertainGraph(long n){
//    create the empty graph with n vertices
    
    igraph_empty(&graph, (igraph_real_t)n, IGRAPH_DIRECTED);
    nv=n;
    ne=0;
}

UncertainGraph::UncertainGraph(const UncertainGraph & obj){
//  copy the structure from other object
    igraph_copy(&graph, &obj.graph);
    nv=obj.nv;
    ne=obj.ne;
    igraph_vector_copy(&pe, &obj.pe);
}

UncertainGraph::UncertainGraph(UncertainGraph &obj){
//
    igraph_copy(&graph, &obj.graph);
    nv=obj.nv;
    ne=obj.ne;
    //igraph_vector_init(&pe,ne);
    igraph_vector_copy(&pe, &obj.pe);
    
}

UncertainGraph& UncertainGraph::operator=(const UncertainGraph & obj){
//    
    igraph_copy(&graph, &obj.graph);
    nv=obj.nv;
    ne=obj.ne;
    // igraph_vector_init(&pe,ne);
    igraph_vector_copy(&pe, &obj.pe);
    return *this;
}




UncertainGraph::~UncertainGraph(){
    // cout<<"call uncertain deconstructor function"<<endl;
    igraph_destroy(&graph);
    if(ne!=0){
        igraph_vector_destroy(&pe);
    }
    
}

void UncertainGraph::setEdges(igraph_vector_t *edges){
    igraph_add_edges(&graph, edges, 0);
    ne=igraph_vector_size(edges);
    ne/=2;
    if(ne!=0){
        igraph_vector_init(&pe,ne);
    }
}

void UncertainGraph::setEdgeProbs(igraph_vector_t *probs){
    igraph_vector_copy(&pe, probs);
    
}


vector<double> UncertainGraph::edgeProbs(igraph_integer_t v){
    //resize
    igraph_vector_t eids;
    vector<double> result;

    igraph_incident(&graph, &eids, v, IGRAPH_ALL);
    
    
    for(int i=0;i<igraph_vector_size(&eids);i++){
        long   eid=VECTOR(eids)[i];
        double p=VECTOR(pe)[eid];
        result.push_back(p);
    }
    igraph_vector_destroy(&eids);
    
    return result;
}

vector<ContinousDistribution> UncertainGraph::degreeDistributions(){
    vector<ContinousDistribution> result;
    
    for(int i=0;i<nv;i++){
        vector<double> eps;
        eps=edgeProbs(igraph_integer_t (i));
        ContinousDistribution dd= ContinousDistribution(eps);
        result.push_back(dd);
    }
    return result;
}