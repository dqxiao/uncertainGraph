//
//  main.cpp
//  uncertainGraph
//
//  Created by dongqingxiao on 10/4/16.
//  Copyright © 2016 dongqingxiao. All rights reserved.
//

#include <iostream>
#include <string>
#include <UncertainGraph.hpp>
#include <GraphIO.hpp>
#include "SCCluster.hpp"
#include <iomanip>  //set precsion
using namespace std;
using namespace arma;


void test(string dataset){
    string homeFolder="/Users/dongqingxiao/Documents/uncertainGraphProject/allDataSet/input/";
    string filePath=homeFolder+dataset+".txt";
    string targetFolder="/Users/dongqingxiao/Documents/uncertainGraphDis/disDataset/";
    string targetFilePath=targetFolder+dataset+"BC.txt";
    string targetFilterFilePath=targetFolder+dataset+"_filter"+".txt";
    //done

    UncertainGraph ug=readUncertainGraph(filePath,'\t');
    
    ug.setColdStart(1.01, 0.1, 1);
    ug.setPrivacyConstraint(60, 0.001);
    
    double robf=1;
    
    ug.genObfuscation(robf);
    
    
}


void testSim(){
    DistPrameter x(10,3);
    DistPrameter y(11,2);
    cout<<x.similiarity(y)<<endl;
    cout<<x.similiarityNormal(y)<<endl;
}


void testMin_k_cluster(){
    
    string filePath="/Users/dongqingxiao/Documents/uncertainGraphDis/disDataset/testConverge.txt";
    string folder="/Users/dongqingxiao/Documents/uncertainGraphDis/disDataset/";
    vector<DistPrameter> vs;
    readDistParameter(filePath, vs);
    int k=40;
    
    SCCluster sc(k,2*k,0);
    
    sc.clustering(vs);
    
    string targetFilePath=folder+"test_sc_"+to_string(k)+".txt";
    writeToFile(targetFilePath, vs, sc);  //done
}






int main(int argc, const char * argv[]) {
    
    test("dblp");
//    unordered_map<int, double> sd;
//    sd[1]=2;
//    
//    
//    cout<<sd[1]<<endl;
//    cout<<sd[2]<<endl;
    
//    cout<<
    
    
    return 0;
}
