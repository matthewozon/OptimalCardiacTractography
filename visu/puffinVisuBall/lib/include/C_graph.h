#ifndef C_GRAPH_H
#define C_GRAPH_H
#include <settings_def.h>

#include <iostream>
#include <string.h>
#include <struct.h>

#include <math.h>
//#define SMALL_NUM 1e-14


#include <stdio.h>
#include <stdlib.h>
//#include <C_toolbox_eigen_sym.h>
#include <C_toolbox_rand2.h>
//#include <C_toolbox_SVD.h>
#include <rawData.h>


using namespace std;


class C_graph
{
    public:
        /** Default constructor: make sure that there are no points with no neighbors*/
        C_graph(rawData<double>* mask /*rawData+diffInfo, Points*+Edge*+Np+Ne*/);
        C_graph(C_graph &G);
        /** Default destructor */
        virtual ~C_graph();


        //methods to read, write and show elements of the graph the graph
        bool saveEdges(char * filename/*="edges.cformat"*/); //no need to save points, cause they are always the same now... but later they might move :-)
        unsigned long getLenPoints(void) const {return nb_points;};
        unsigned long getLenEdges(void) const {return nb_edges;};
        Point* getPoints(void) const {return m_points;};
        Edge* getEdges(void) const {return m_edges;};

        //tools
        Edge* findEdge(Point* p1, Point* p2);
        Point* findOtherPoint(Edge* e, Point* p);
        Point* findCommunPoint(Edge* e1, Edge* e2);
        unsigned long getNbConnectedEdge(void);
        void randInit(double alphaConnected);

        long getIndexOfPair(Point* p, long idx1, long idx2);
        double getPixDX(){return pix_dx;};
        double getPixDY(){return pix_dy;};
        double getPixDZ(){return pix_dz;};
        double getXlen(){return DX;};
        double getYlen(){return DY;};
        double getZlen(){return DZ;};
        rawData<double>* getMask(){return m_mask;};

    private:
        //relevant attributes
        unsigned long nb_points;
        unsigned long nb_edges;
        Point *m_points; //so that there is always a back up of the first address
        Edge *m_edges; //so that there is always a back up of the first address
        rawData<double>* m_mask;

        double pix_dx, pix_dy, pix_dz; ///lattice
        double DX, DY, DZ; //embeding box size

        unsigned long countNbPointInMask(/**rawData<double>* maskor a segmented rawData*/);
        Point* createPointArrayFromMask(/**rawData<double>* maskor a segmented rawData*/);
        unsigned long countNbEdgesInMask(void);
        Edge* createEdgeArray(void);

        bool bindPointsAndEdges(void);

        bool copyGraph(Point** points, Edge** edges, Point* m_points_source, Edge* m_edges_source, unsigned long Np, unsigned long Ne);

};

#endif // C_GRAPH_H
