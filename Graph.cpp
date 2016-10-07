#include "Graph.hpp"


Graph::Graph(long n){
    igraph_empty(&graph, (igraph_real_t)n, IGRAPH_DIRECTED);// create empty graph
    nv=n;
    //cout<<"init nv"<<nv<<endl;
}

Graph::Graph(const Graph & obj){
    igraph_copy(&graph, &obj.graph);
    nv=obj.nv;
    ne=obj.ne;
}

Graph::Graph(Graph &obj){
    igraph_copy(&graph, &obj.graph);
    nv=obj.nv;
    ne=obj.ne;
}



