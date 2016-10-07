//
//  GraphIO.hpp
//  uncertainGraph
//
//  Created by dongqingxiao on 10/6/16.
//  Copyright © 2016 dongqingxiao. All rights reserved.
//

#ifndef GraphIO_hpp
#define GraphIO_hpp

#include <stdio.h>
#include <igraph.h>
#include <Graph.hpp>
#include <UncertainGraph.hpp>


UncertainGraph readUncertainGraph(string filepath, char sep);
UncertainGraph readUncertainGraphAdjacency(string filepath,long nv);

void writeDistribution(UncertainGraph g,string filepath);


#endif /* GraphIO_hpp */
