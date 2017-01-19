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
using namespace std;

void test(string dataset){
    string homeFolder="/Users/dongqingxiao/Documents/uncertainGraphProject/allDataSet/input/";
    string filePath=homeFolder+dataset+".txt";
    string targetFolder="/Users/dongqingxiao/Documents/uncertainGraphDis/disDataset/";
    string targetFilePath=targetFolder+dataset+"BC.txt";
    string targetFilterFilePath=targetFolder+dataset+"_filter"+".txt";
    //done

    UncertainGraph ug=readUncertainGraph(filePath,'\t');
    
    vector<DistPrameter> dps=ug.degreeDistributions(3);
    
//    cout<<dps.size()<<endl;
    SCCluster sc(100,200,0);
    
    sc.clustering(dps);
    
    writeToFile(targetFilePath, dps, sc);
    
    
    
}


int main(int argc, const char * argv[]) {
    
    
    test("hepth");
    
    return 0;
}
