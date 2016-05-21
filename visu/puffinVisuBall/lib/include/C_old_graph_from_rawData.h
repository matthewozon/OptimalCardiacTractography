#ifndef C_OLD_GRAPH_FROM_RAWDATA_H
#define C_OLD_GRAPH_FROM_RAWDATA_H

#include <rawData.h>
#include <iostream>
#include <stdio.h>

using namespace std;

class C_old_graph_from_rawData //assume that all data are unmosa:
{
    public:
        /** Default constructor */
        C_old_graph_from_rawData();
        /** Default destructor */
        virtual ~C_old_graph_from_rawData();

        unsigned long countNbPointInMask(rawData<double>* mask/**or a segmented rawData*/);
        Point* createPointArrayFromMask(rawData<double>* mask/**or a segmented rawData*/);
        unsigned long countNbEdgesInMask(Point* p, unsigned long Np);
        Edge* createEdgeArray(unsigned long Ne);

        bool bindPointsAndEdges(Point* p, unsigned long Np, Edge* e, unsigned long Ne);

        bool savePoints(char * filename, PointWithNeighbors* m_pointsWithNeighbors, unsigned long Np);
        PointWithNeighbors* getPointWithNeighborsFromPointAndEdge(Point* p, unsigned long Np, Edge* e, unsigned long Ne);
};

#endif // C_OLD_GRAPH_FROM_RAWDATA_H
