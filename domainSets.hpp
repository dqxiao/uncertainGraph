//
//  domainSets.hpp
//  uncertainGraph
//
//  Created by dongqingxiao on 1/18/17.
//  Copyright © 2017 dongqingxiao. All rights reserved.
//

#ifndef domainSets_hpp
#define domainSets_hpp

#include "util.hpp"
using namespace arma;

umat dominantSetsClustering(fmat S,int N_MAX_CLUSTER_SIZE = 50 , int N_MIN_CLUSTER_SIZE = 3, int N_OVERLAP_SIZE = 2);

int dominantSetExtraction(fmat S,uvec& dominant, uvec& non_dominant, int N_MAX_CLUSTER_SIZE, int N_MIN_CLUSTER_SIZE);

#endif /* domainSets_hpp */
