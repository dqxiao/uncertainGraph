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
#include <algorithm>
#include <unordered_map>
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
    cout<<"ending the etimation"<<endl;
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
    
    vector<int> uniqueVS;  // storing the unique vertex in the obfuscated graph
    vector<double> vsComm;  // storing the sum of ....
    vector<DistPrameter> obDist=obGraph.degreeDistributions(); // get the degree distribution
    int md=(int) gdd.size()-1;
    
    for(int i=0;i<nv;i++){
        
        DistPrameter ed=obDist[i];
        vector<double> ddist;
        ed.degreeDistribution(ddist);
        int l=ed.probMinVal();
        int h=ed.probMaxVal();
        double result=0;
        for(int d=l;d<=h;d++){
            result+=ddist[d-l]*gdd[min(d,md)];
            if(result>k){
                break;
            }
        }
        
        if(result<k){
            uniqueVS.push_back(i);
            vsComm.push_back(result);
        }
    }
    
    int N_UOBF=(int)uniqueVS.size();
    cout<<"Raw number of under-obfuscated vertices:"<<N_UOBF<<endl;
    cout<<"Torelence level:"<<(double)N_UOBF/nv<<endl;
    
    double target=log2((double)k);
    vector<DistPrameter> dist=degreeDistributions();
    int N_OBF=0;
    for(int i=0; i<N_UOBF;i++){
        DistPrameter uvDP=obDist[uniqueVS[i]];
        double ce=0;
        double ctarget=target-log2(vsComm[i]);
        ctarget*=vsComm[i];
        for(int v=0;v<nv;v++){
            DistPrameter vDP=obDist[uniqueVS[i]];
            ce+=entropy(uvDP.similiarityNormal(vDP));
            if(ce>ctarget){
                N_OBF+=1;
                break;
            }
        }
        //I want to update this issue as quick way
    }
    
    
    //just for debugging
    //exact computing is necessary or not needed??
    N_UOBF-=N_OBF;
    cout<<"The number of under-obfuscated vertices:"<<N_UOBF<<endl;
    cout<<"Torelence level:"<<(double)N_UOBF/nv<<endl;
    *r=(double)N_UOBF/nv;
}

//excluding common vertex 
void UncertainGraph::vertexCommon(vector<double> & vs_comm, int upper,double &vc_sum){
    if(gdd.size()==0){
        approxGlobalDegreeDist(gdd);
    }
    
    vector<DistPrameter> dist=degreeDistributions();
    int md=(int) gdd.size()-1;
    double apprSum=0;     //I thought this will work and work well
    for(int i=0;i<nv;i++){
        
        DistPrameter ed=dist[i];
        vector<double> ddist;
        ed.degreeDistribution(ddist);
        int l=ed.probMinVal();
        int h=ed.probMaxVal();
        double result=0;
        for(int d=l;d<=h;d++){
            result+=ddist[d-l]*gdd[min(d,md)];
            if(result>upper){
                break;
            }
        }

        vs_comm[i]=result; // huge value *k
        apprSum+=result;
        
    }
    
    vc_sum=apprSum;

}

void UncertainGraph::approximateUniquness(vector<double> &vs_uniq, int k, int scaleUpper){

    int upper=k*scaleUpper;
    double sum=0;
    vertexCommon(vs_uniq,upper,sum); // such very common one nodes
    
//    for(int i=0;i<vs_uniq.size();i++){
//        vs_uniq[i]=sum/vs_uniq[i];
//    }
    
//    write
    
    
//    sum=1/sum;                                 // normalization
//    transform(vs_uniq.begin(), vs_uniq.end(), vs_uniq.begin(), std::bind2nd(std::multiplies<double>(), sum));
    writeUniquess("/Users/dongqingxiao/Documents/uncertainGraphDis/uniqunessLog/dblp_comm.txt",vs_uniq);
    

}

// configuration
void UncertainGraph::configuration(bool g_sc, bool opt){
    g_sc=g_sc;
    opt=opt;
}

// set privacy constraint
void UncertainGraph::setPrivacyConstraint(int pk, float t){
    k=pk;
    tolerance=t;
}

void UncertainGraph::setColdStart(float sizeMultiplier, float sigma, int nAttempt){
    c=sizeMultiplier;
    sigma=sigma;
    t=nAttempt; 
}



void UncertainGraph::excludingNodes(vector<double> & vProbs){
    vector<Node_UN> nodes;
    int skip=(int)lround(tolerance*nv/2);
    
    cout<<"nv"<<nv<<endl;
    cout<<"tor"<<tolerance<<endl;
    cout<<"skip"<<skip<<endl;
    
    //cacluate the key difference
    for(int i=0; i<nv;i++){
        float v_uniq_val=vs_uniq[i];
        nodes.push_back(Node_UN(i,v_uniq_val));
        vProbs[i]=v_uniq_val;
    }
    
    sort(nodes.begin(),nodes.end(),[] (const Node_UN & lhs, const Node_UN &rhs ){
        return lhs.uniq> rhs.uniq;
    });
    
    //excluding 
    for(int i=0;i<skip;i++){
        Node_UN nu=nodes[i];
        vProbs[nu.vid]=0;
    }
    
}

//obfuscation graph with different method
void UncertainGraph::obfuscation(string obfMethod){
    //done 
    
}


bool exist(unordered_map<Edge, double> adjMatrix,Edge e){
    
    auto search=adjMatrix.find(e);
    if(search!=adjMatrix.end()){
        return true;
    }
    
    return false;
}

double existVal(unordered_map<Edge, double> adjMatrix,Edge e){
    
    if(not exist(adjMatrix,e)){
        return -1;
    }else{
        return adjMatrix[e];
    }
}


//done
void UncertainGraph::genObfuscation(double & robf){
    //randomized method: try t times
    //compute the uniquness of vertex
    if(vs_uniq.size()==0){
        vs_uniq=vector<double>(nv,0);
        approximateUniquness(vs_uniq, k, 1000);
        cout<<"<--cal vs_uniq end-->"<<endl;
        return; // just for this time
    }else{
        cout<<"<--local vs_uniq end-->"<<endl;
    }
    
    cout<<"Excluding nodes with highest uniquness score"<<endl;
    if(vProbs.size()==0){
        vProbs=vector<double>(nv,0.0);
        excludingNodes(vProbs);
    }
    
    fast_discrete_distribution<int> nodeSampler(vProbs);
    uniform_real_distribution<double> unDist(0.0,1.0);
 
    igraph_vector_t eids;
    igraph_vector_init(&eids, 0);
    igraph_get_edgelist(&graph, &eids, false);
//    sp_mat adjMatrix(nv,nv);
    unordered_map<Edge, double> adjMatrix;
    
    
    
    int c_edge=0;
    for(int i=0;i<ne;i++){
        int from=VECTOR(eids)[2*i];
        int to=VECTOR(eids)[2*i+1];
        double edgeProb=VECTOR(pe)[i];
        
        if(from>to){
            swap(from,to);
        }
        if(from==to){
            cout<<"error exist for self loop"<<endl;
            
        }
        if (edgeProb!=0){
            adjMatrix[Edge(from,to)]=edgeProb;
//            adjMatrix(from,to)=edgeProb;
        }else{
            cout<<"error by precision"<<endl;
        }
        
    }
    
    

    //randomized generation of obfucation candidate
    //t: attempt number, default setting: 5
    for(int i=0;i<t;i++){
        
        cout<<"iter " <<i<<endl;
        
        
    
        std::default_random_engine gen;
        double cobf=1;
        //the key idea is the inc and remove edge is relative smaller than adjMatrix
        unordered_map<Edge, double> rmMatrix;
        unordered_map<Edge, double> incMatrix;
        
        
        vector<int> sampledEdges;
        vector<int> existedEdges;
        int incEdge=0;
        
        
       //done
        
        
        
        int edgeCount=(int)ne;
        int ce=c*edgeCount;
        
        while(edgeCount!=ce){
            
            int from=nodeSampler(gen);
            int to=nodeSampler(gen);
            
            if(from==to){
                continue;
            }
        
            if(from>to){
                swap(from,to);
            }
            
           
            
            //transfer prob from existing prob to non-existing edges
            //I need collect the edges
            
            //case analysis
            
            // edge not in adjMatrix: can be added into obfMatrix only once
            // edge in adjMatrix: can be removed only once???
            
            Edge e(from,to);
            
            igraph_bool_t res;
            
            
            igraph_are_connected(&graph,from,to,&res);
            if(res){
                
                double p=unDist(gen);
                double pval=adjMatrix[e];
                
                if(pval>=p){
                    // decided to remove such edge
                    
                    if(!exist(rmMatrix,e)){
                        rmMatrix[e]=-1;
                        edgeCount-=1;
                    }
                
                }
                
              
            
            }else{
                // not exist in the original graph
                if(not exist(incMatrix,e))
                {
                    incMatrix[e]=2;
                    edgeCount+=1;
                    sampledEdges.push_back(from);
                    sampledEdges.push_back(to);
                    incEdge+=1;
                    
                }
            
            }
           
            
            
        }
        
        
        cout<<"check "<<sampledEdges.size()<<endl;
        cout<<"check "<<incEdge<<endl;
        
        
        for(int e=0;e<ne;e++){
            int from=VECTOR(eids)[2*e];
            int to=VECTOR(eids)[2*e+1];
            
            if(from>to){
                swap(from,to);
            }
            
            
            Edge te(from,to);
            
            if(not exist(rmMatrix,te)){
                sampledEdges.push_back(from);
                sampledEdges.push_back(to);
            }
            
        }
        
        
        cout<<"get "<<ce<<" edges by sampling"<<endl;
        cout<<"get "<<edgeCount <<"edges as required"<<endl;
        cout<<"check "<<sampledEdges.size()<<endl;
        
    
        
        
        
//        if(cobf<robf){
//            robf=cobf;
//            
//        }
        
    }

}

