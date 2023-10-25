//
// Created by alvar on 25/10/2023.
//
#ifndef EDGE_H
#define EDGE_H
template<typename SizeType>
struct Edge{
    Edge(SizeType source, SizeType target);
    SizeType source;
    SizeType target;
};

template<typename SizeType>
Edge<SizeType>::Edge(SizeType source, SizeType target) {
    this->source=source;
    this->target=target;
}

#endif