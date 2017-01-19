//
//  UncertainGraph.cpp
//  uncertainGraph
//
//  Created by dongqingxiao on 10/6/16.
//  Copyright © 2016 dongqingxiao. All rights reserved.
//

#include "UncertainGraph.hpp"
#include <boost/math/distributions/normal.hpp>
#include <random>
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
    igraph_copy(&graph, &obj.graph);
    nv=obj.nv;
    ne=obj.ne;
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
    
    igraph_vector_init(&eids, nv);
    igraph_incident(&graph, &eids, v, IGRAPH_ALL);
    
    for(int i=0;i<igraph_vector_size(&eids);i++){
        long   eid=VECTOR(eids)[i];
        double p=VECTOR(pe)[eid];
        result.push_back(p);
    }
    igraph_vector_destroy(&eids);
    
    return result;
}

vector<DistPrameter> UncertainGraph::degreeDistributions(){
    vector<DistPrameter> result;
    for(int i=0;i<nv;i++){
        vector<double> eps;
        eps=edgeProbs(igraph_integer_t (i));
        if(eps.size()==0){
            continue;
        }
        DistPrameter dp(eps);
        result.push_back(dp);
    }
    return result;
}


vector<DistPrameter> UncertainGraph::degreeDistributions(double lb){
    vector<DistPrameter> result;
    for(int i=0;i<nv;i++){
        vector<double> eps;
        eps=edgeProbs(igraph_integer_t (i));
        if(eps.size()==0){
            continue;
        }
        
        DistPrameter dp(eps);
        if(dp.mean>lb and dp.mean<100){
            result.push_back(dp);
        }
    }
    return result;
}


vector<DistPrameter> UncertainGraph::degreeDistributions(vector<int> vs){
    vector<DistPrameter> result;
    for(auto v: vs){
        vector<double> eps;
        eps=edgeProbs(igraph_integer_t (v));
        DistPrameter dp(eps);
        result.push_back(dp);
    }
    return result;
}

int UncertainGraph::getMaxProbDegree(){
    int result=0;
    for(int i=0;i<nv;i++){
        vector<double> eps;
        eps=edgeProbs(igraph_integer_t(i));
        if(eps.size()>=30){
            DistPrameter dp(eps);
            int h=dp.probMaxVal();
            result=std::max(result,h);
        }else{
            result=std::max(result,(int)eps.size());
        }
    }
    
    return result;
}


void UncertainGraph::approxGlobalDegreeDist(vector<double> &result){
    
    int nsample=1000;
    
    double p=1.0/nsample;
    random_device rand_dev;
    mt19937     generator(rand_dev());
    uniform_real_distribution<double> unDist(0.0,1.0);
    igraph_vector_t edges;
    igraph_vector_t vertexDegree;
    igraph_vector_init(&edges, 2*ne);
    igraph_get_edgelist(&graph, &edges, false);
    
    igraph_vector_init(&vertexDegree, nv);
    
    igraph_degree(&graph,&vertexDegree,igraph_vss_all(), IGRAPH_ALL, true);
    
    int maxDegree=igraph_vector_max(&vertexDegree);
    int realMax=0;
    vector<double> r(maxDegree+1,0); //
    
    for(int i=0;i<nsample;i++){
        //generate sample graph by
        vector<int> vds(int(nv),0);
        for(int eid=0;eid<ne;eid++){
            igraph_real_t eProb=VECTOR(pe)[eid];
            float x=unDist(generator);
            int source=(int)VECTOR(edges)[2*eid];
            int end=(int)VECTOR(edges)[2*eid+1];
            
            if(eProb>=x){
                vds[source]+=1;
                vds[end]+=1;
            }
        }
        
        for(int v=0;v<nv;v++){
            int d=vds[v];
            r[d]+=p;
            realMax=std::max(realMax,d);
        }

    }
    
    vector<double> rf(&r[0],&r[realMax]);
    
    result=rf;
    
}

/**
 * vs: vertex similiarity estiamtion 
 * gdd: global degree distribution estimation
 */

void UncertainGraph::approxVertexSim(vector<double> &vs,vector<double> &gdd){
    
    
    
    int nsample=1000;
    
    double p=1.0/nsample;
    random_device rand_dev;
    mt19937     generator(rand_dev());
    uniform_real_distribution<double> unDist(0.0,1.0);
    igraph_vector_t edges;
    igraph_vector_init(&edges, 2*ne);
    igraph_get_edgelist(&graph, &edges, false);
    
    int md=(int)gdd.size()-1;
    
    for(int i=0;i<nsample;i++){
        //generate sample graph by
        vector<int> vds(int(nv),0);
        for(int eid=0;eid<ne;eid++){
            igraph_real_t eProb=VECTOR(pe)[eid];
            float x=unDist(generator);
            int source=(int)VECTOR(edges)[2*eid];
            int end=(int)VECTOR(edges)[2*eid+1];
            
            if(eProb>=x){
                vds[source]+=1;
                vds[end]+=1;
            }
        }
        
        for(int i=0;i<nv;i++){
            int vd=min((int)vds[i],md);
            vs[i]+=p*gdd[vd];
        }
        
    }
    
    //done
    
}




void UncertainGraph::globalDegreeDist(vector<double> &result){
    
    int md=getMaxProbDegree();
    vector<double> r(md+1,0);
    
    for(int i=0;i<nv;i++){
        vector<double> eps;
        vector<double> vResult;
        eps=edgeProbs(igraph_integer_t(i));
        
        if(eps.size()>=30){
            
            DistPrameter dp(eps);
            
            int l=dp.probMinVal();
            int h=dp.probMaxVal();
            dp.degreeDistribution(vResult);
            
            for(int i=l;i<=h;i++){
                r[i]=r[i]+vResult[i-l];
            }
            
        }else{
            ExactDist ed(eps);
            int l=ed.probMinVal();
            int h=ed.probMaxVal();
            ed.degreeDistribution(vResult);
    
            for(int i=l;i<=h;i++){
                r[i]+=vResult[i-l];
            }
           
            
        }
        
    }
    //done
    result=r;
}


void UncertainGraph::anonymityCheck(UncertainGraph obGraph, int k, double * r){
    //inititation
    if(gdd.size()==0){
        approxGlobalDegreeDist(gdd);
    }
    
    vector<int> remaingVertex;
    vector<double> cumSum;
    vector<DistPrameter> obDist=obGraph.degreeDistributions(); // get the degree distribution
    
    for(int i=0;i<nv;i++){
        
        DistPrameter ed=obDist[i];
        vector<double> result=ed.
        int l=ed.probMinVal();
        int h=ed.probMaxVal();
        
        
    
    }
    
}






