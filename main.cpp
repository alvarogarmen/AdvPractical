//
// Created by alvar on 26/04/2023.
//

#include <string>
#include "argtable3.h"
#include <chrono>

#include <iostream>



int main(int argc, char* argv[]) {
    //Argtable 3
    double sourceNode;
    double targetNode;
    int numLandmarks;
    int newLandmarks;
    int newPoints;
    int numPoints;
    std::string graph;

    struct arg_dbl* sourceArg = arg_dbl0(NULL, "source", "<double>", "source node");
    struct arg_dbl* targetArg = arg_dbl0(NULL, "target", "<double>", "target node");
    struct arg_int* landmarkArg = arg_int0(NULL, "landmarks", "<integer>", "number of landmarks");
    struct arg_int* newLandmarkArg = arg_int0(NULL, "newLandmarks?", "<integer>", "0 for old landmarks");
    struct arg_str* graphArg = arg_str0(NULL, "graph", "<string>", "file with the graph, keep the .graph!");
    struct arg_int* newPointsArg = arg_int0(NULL, "newPoints?", "<integer>", "generate new Points?");
    struct arg_int* numPointsArg = arg_int0(NULL, "numNewPoints?", "<integer>", "how many new Points?");


    struct arg_end* end = arg_end(20); // Define the end marker for the argtable array

    void* argtable[] = { sourceArg, targetArg, landmarkArg, newLandmarkArg, graphArg, newPointsArg, numPointsArg, end };

    const char* progname = "myprogram";
    int nerrors = arg_parse(argc, argv, argtable);

    if (nerrors > 0) {
        arg_print_errors(stdout, end, progname);
        arg_print_syntax(stdout, argtable, "\n");
        return 1; // Handle parsing errors
    }

    if (sourceArg->count > 0) {
        sourceNode = sourceArg->dval[0];
    }

    if (targetArg->count > 0) {
        targetNode = targetArg->dval[0];
    }
    if (landmarkArg->count > 0){
        numLandmarks = landmarkArg->ival[0];
    }
    if(newLandmarkArg->count >0){
        newLandmarks = newLandmarkArg->ival[0];
    }
    if (graphArg->count > 0) {
        graph = graphArg->sval[0];
    }
    if (numPointsArg->count > 0) {
        numPoints = numPointsArg->ival[0];
    }
    if (newPointsArg->count > 0) {
        newPoints = newPointsArg->ival[0];
    }
    // Call out functions

    arg_freetable(argtable, sizeof(argtable) / sizeof(argtable[0]));

    return 0;
}