//
//  SCCluster.cpp
//  uncertainGraph
//
//  Created by dongqingxiao on 1/18/17.
//  Copyright © 2017 dongqingxiao. All rights reserved.
//

#include "SCCluster.hpp"
#include <sstream>
#include <fstream>


SCCluster::SCCluster(int N_MIN, int N_MAX, int N_OVERLAP){
    minSize=N_MIN; // the min-size per cluster
    maxSize=N_MAX; // the max-size per cluster
    overlap=N_OVERLAP; // the number of overlapping point per cluster
}



void SCCluster::clustering(fmat simMatrix){
    clusterMatrix = dominantSetsClustering(simMatrix,maxSize,minSize,overlap);  //done

}


void SCCluster::clustering(vector<DistPrameter> &vs){
    //build similiarityMatrix
    fmat s=buildSimiliary(vs, "EL");
    //clustering
    clustering(s);
}


fmat SCCluster::buildSimiliary(vector<DistPrameter> &vs, string funcName){
    
    cout<<"build similiarty Matrix"<<endl;
    int NC=(int)vs.size();
    
    fmat simMatrix(NC,NC);
    
    for(int i=0;i<NC;i++){
        simMatrix(i,i)=0;
        DistPrameter vi=vs[i];
        
        for(int j=i+1;j<NC;j++){
            DistPrameter vj=vs[j];
            
            float d=vi.distance(vj,funcName);
        
            simMatrix(i,j)=d;
            simMatrix(j,i)=d;
        }
    
    }
    
    //re-scale
    float stdev = as_scalar(stddev(stddev(simMatrix,0,1),0,0));
    simMatrix= exp(-simMatrix/stdev);
    
    return simMatrix;
}



// writeToFile
void writeToFile(string filePath,vector<DistPrameter>&vs, SCCluster sc){

    ofstream myfile(filePath);
    
    if(!myfile.is_open()){
        cout<<"can not open this file"<<endl;
        return;
    }
    
    umat ol=sum(sc.clusterMatrix,0); // take sum by each column
    int clusterNum=(int)sc.clusterMatrix.n_rows;
    for(int i=0;i<vs.size();i++){
        DistPrameter v=vs[i];
        if(ol(0,i)!=1){
            myfile<<v.mean<<"\t"<<v.stddev<<"\t"<<clusterNum<<endl;
        }else{
            for(int j=0;j<clusterNum;j++){
                //row,column
                if(sc.clusterMatrix(j,i)==1){
                     myfile<<v.mean<<"\t"<<v.stddev<<"\t"<<j<<endl;
                    break;
                }
            }
        
        }
    }
    
    cout<<"finished the output"<<endl;
    
}

