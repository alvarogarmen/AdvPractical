//
// Created by alvar on 25/10/2023.
//
#ifndef GRADER_H
#define GRADER_H

#include "Edge.hpp"
#include "Graph.hpp"

template<typename SizeType>
bool edgeCross(Edge<SizeType>& edge1, Edge<SizeType>& edge2){
    if (edge1.source<edge2.source && edge1.target>edge2.target || edge1.source>edge2.source && edge1.target<edge2.source){
        return true;
    }
    return false;
}
template<typename SizeType>
int crossGrader(Graph<SizeType>& myGraph){
    int crossings = 0;
    for (int i = 0; i<myGraph.edges.size(); i++){
        for (int j = i; j<myGraph.edges.size(); j++){
            if(edgeCross(myGraph.edges[i], myGraph.edges[j])){
                crossings++;
            }
        }
    }
    return crossings;
}

#endif