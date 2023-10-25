//
// Created by alvar on 26/04/2023.
//

#include <string>
#include "argtable3.h"
#include "Graph.hh"
#include "Grader.hh"
#include "InputGraphManually.hh"
#include <chrono>

#include <iostream>

void testFunction(){
    Graph myGraph = Graph(0, 0, 0);

    Node node1 = Node(1, 2);
    Node node2 = Node(2, 1);
    Node node3 = Node(3, 1);
    Node node4 = Node(4, 1);

    myGraph.insertLeftNode(node1);
    myGraph.insertLeftNode(node2);
    myGraph.insertLeftNode(node3);
    myGraph.insertLeftNode(node4);

    myGraph.insertEdge(1, 1);
    myGraph.insertEdge(1, 2);
    myGraph.insertEdge(2, 2);
    myGraph.insertEdge(3, 1);
    myGraph.insertEdge(4, 3);

    std::cout<<crossGrader(myGraph)<<std::endl;


}

int main(int argc, char* argv[]) {
    //Argtable 3
    int leftSize;
    int rightSize;
    int edgeSize;


    struct arg_int* leftSizeArg = arg_int0(NULL, "leftSize", "<integer>", "number of nodes on left Graph");
    struct arg_int* rightSizeArg = arg_int0(NULL, "rightSize", "<integer>", "number of nodes on right Graph");
    struct arg_int* edgeSizeArg = arg_int0(NULL, "edgeSize", "<integer>", "number of edges");


    struct arg_end* end = arg_end(20); // Define the end marker for the argtable array

    void* argtable[] = { leftSizeArg, rightSizeArg, edgeSizeArg, end };

    const char* progname = "myprogram";
    int nerrors = arg_parse(argc, argv, argtable);

    if (nerrors > 0) {
        arg_print_errors(stdout, end, progname);
        arg_print_syntax(stdout, argtable, "\n");
        return 1; // Handle parsing errors
    }

    if (leftSizeArg->count > 0) {
        leftSize = leftSizeArg->ival[0];
    }

    if (rightSizeArg->count > 0) {
        rightSize = rightSizeArg->ival[0];
    }
    if (edgeSizeArg->count > 0){
        edgeSize = edgeSizeArg->ival[0];
    }

    // Call out functions

    testFunction();
    Graph myGraph = inputGraphManually();
    std::cout<<"Crossings: "<< crossGrader(myGraph)<<std::endl;
    arg_freetable(argtable, sizeof(argtable) / sizeof(argtable[0]));

    return 0;
}