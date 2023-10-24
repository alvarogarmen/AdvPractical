//
// Created by alvar on 23/10/2023.
//

#ifndef ADVPRACTICAL_GRADER_H
#define ADVPRACTICAL_GRADER_H

#include <vector>
#include "Edge.hh"
#include "Graph.hh"
bool edgeCross(Edge& edge1, Edge& edge2){
    if (edge1.source<edge2.source && edge1.target>edge2.target || edge1.source>edge2.source && edge1.target<edge2.source){
        return true;
    }
    return false;
}

int crossGrader(Graph& myGraph){
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

#endif //ADVPRACTICAL_GRADER_H
