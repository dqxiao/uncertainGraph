//
//  Graph.hpp
//  uncertainGraph
//
//  Created by dongqingxiao on 10/4/16.
//  Copyright © 2016 dongqingxiao. All rights reserved.
//

#ifndef Graph_hpp
#define Graph_hpp

#include <stdio.h>
#include <igraph.h>
#include <string>
using namespace std;


class Graph{
protected:
    igraph_t graph;
public:
    long nv;
    long ne;
    Graph(long nv);      // empty constructor
    Graph(const Graph & obj); // copy const constructor
    Graph(Graph & obj); // copy constructor
    ~Graph(); // destructor function
    void set_edges(igraph_vector_t * edges); // add edges to graph
    
};

// Graph I/O
Graph init_from_file(string filepath);
Graph init_from_Adj_File(string filepath);




#endif /* Graph_hpp */

