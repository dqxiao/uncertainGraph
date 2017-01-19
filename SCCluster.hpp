//
//  SCCluster.hpp
//  uncertainGraph
//
//  Created by dongqingxiao on 1/18/17.
//  Copyright © 2017 dongqingxiao. All rights reserved.
//

#ifndef SCCluster_hpp
#define SCCluster_hpp

#include <stdio.h>
#include <Debug.hpp>
#include <vector>
#include "util.hpp"
#include "domainSets.hpp"



class SCCluster{
public:
    int minSize;
    int maxSize;
    int overlap;
    umat clusterMatrix;
    SCCluster(int N_MIN, int N_MAX, int N_OVERLAP);
    void clustering(vector<DistPrameter> &vs); //done
private:
    void clustering(fmat s); // clustering by similiarity matrix
    fmat buildSimiliary(vector<DistPrameter> &vs, string funcName); // pass
};


void writeToFile(string filePath,vector<DistPrameter>&vs, SCCluster sc);


#endif /* SCCluster_hpp */
