//
//  GraphIO.cpp
//  uncertainGraph
//
//  Created by dongqingxiao on 10/6/16.
//  Copyright © 2016 dongqingxiao. All rights reserved.
//

#include "GraphIO.hpp"
#include <vector>
#include <fstream>
#include <sstream>
#include <iostream>
#include <stdio.h>
//#include "Debug.hpp"
using namespace std;



UncertainGraph readUncertainGraph(string filepath,char sep){
//    File format
//    Each line: u,s,t
    fstream graphFile(filepath);
    string line;
    string item;
    vector<double> v;
    igraph_vector_t v_graph,e_probs;
    vector<double> edge_prob;
    if(graphFile.is_open()){
        while (getline(graphFile,line)) {
            stringstream ss(line); // for split line
            int i=0;
            while(getline(ss,item,sep)){
                if(i<2){
                    v.push_back(stol(item));
                }else{
                    edge_prob.push_back(stof(item));
                }
                
                i+=1;
            }
        }
        
        
        graphFile.close();
    }else{
        cout<<"can't find/open"<<filepath<<endl;
    }
    

    
    
    igraph_vector_init(&v_graph, v.size());
    
    for(long int i=0;i<v.size();i++){
        igraph_vector_set(&v_graph,i,v[i]);
    }
    igraph_real_t nv=igraph_vector_max(&v_graph);
    UncertainGraph pg(nv+1);
    pg.setEdges(&v_graph);
    
    cout<<"init edge probs:"<<edge_prob.size()<<endl;
    //igraph_vector_view(&e_probs, parray, edge_prob.size());
    igraph_vector_init(&e_probs,edge_prob.size());
    
    for(long int i=0;i<edge_prob.size();i++){
        igraph_vector_set(&e_probs,i,edge_prob[i]);
    }
    
    pg.setEdgeProbs(&e_probs);
    
    
    cout<<"init uncertain graph from"<<filepath<<" in edges, p format"<<endl;
    
    igraph_vector_destroy(&v_graph);
    igraph_vector_destroy(&e_probs);
    return pg;
}

UncertainGraph readUncertainGraphAdjacency(string filepath,long nv){
    ifstream graphFile(filepath);
    string line;
    string item;
    vector<double> v;
    igraph_vector_t v_graph,e_probs;
    vector<double> edge_prob;
    if(graphFile.is_open()){
        while (getline(graphFile,line)) {
            stringstream ss(line); // for split line
            int i=0;
            while(getline(ss,item,'\t')){
                if(i<2){
                    v.push_back(stol(item));
                }else{
                    edge_prob.push_back(stof(item));
                }
                
                i+=1;
            }
        }
        
        
        graphFile.close();
    }else{
        cout<<"can't open"<<filepath<<endl;
    }
    
    double * array=v.data();
    igraph_vector_view(&v_graph, array, v.size());
    UncertainGraph pg((igraph_real_t) nv);
    pg.setEdges(&v_graph);
    
    double * parray=edge_prob.data();
    cout<<"init edge probs:"<<edge_prob.size()<<endl;
    igraph_vector_view(&e_probs, parray, edge_prob.size());
    
    pg.setEdgeProbs(&e_probs);
    
    igraph_vector_destroy(&v_graph);
    igraph_vector_destroy(&e_probs);
    
    return pg;
}



void writeDistribution(UncertainGraph g,string filePath){
//  write mean,std for each vertices   
    vector<DistPrameter> result;
    ofstream myfile(filePath);
    result=g.degreeDistributions();
    for(DistPrameter cd: result){
        myfile<<cd.mean<<","<<cd.stddev<<endl;
    }
    
    myfile.close();
}


void writeDistribution(UncertainGraph g,string filePath,double lb){
    
    vector<DistPrameter> result;
    ofstream myfile(filePath);
    result=g.degreeDistributions();
    for(DistPrameter cd: result){
        if(cd.mean>=lb){
            myfile<<cd.mean<<","<<cd.stddev<<endl;
        }
    }
    
    myfile.close();
}

void writeDistribution(UncertainGraph g,string filePath,vector<int> vs){
    //  write mean,std for each vertices
    vector<DistPrameter> result;
    ofstream myfile(filePath);
    result=g.degreeDistributions(vs);
    for(DistPrameter cd: result){
        myfile<<cd.mean<<","<<cd.stddev<<endl;
    }
    myfile.close();
    //done 
}

// done



