//
// Created by alvar on 24/10/2023.
//

#ifndef ADVPRACTICAL_EDGE_H
#define ADVPRACTICAL_EDGE_H
struct Edge{
    Edge(int source, int target);
    int source;
    int target;
};
Edge::Edge(int source, int target) {
    this->source=source;
    this->target=target;
}
#endif //ADVPRACTICAL_EDGE_H
