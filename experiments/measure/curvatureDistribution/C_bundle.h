#ifndef C_BUNDLE_H
#define C_BUNDLE_H

#include <settings_def.h>
#include <iostream>
#include <C_graph.h>
#include <vector>
#include <fstream>
#include <string>

using namespace std;

class C_bundle
{
    public:
        /** Default constructor */
        C_bundle(C_graph* graph, bool cutAcuteEdges=true);
        /** Default destructor */
        virtual ~C_bundle();

        vector< vector<unsigned long> > FIBERS;

    private:

        //attributs
        bool extractFibers(bool cutAcuteEdges); //
        C_graph* m_graph;
        Point* m_point;
        unsigned long nb_points;
        Edge* m_edge;
        unsigned long nb_edges;

        //linear fibre extraction
        bool extractLinearFibers(void); //
        bool extractLinearFiber(unsigned long idx_start); //

        //loop-shaped extraction
        bool extractLoopFibers(void); //
        bool extractLoopFiber(unsigned long idx); //

        //graph correction
        bool cutAcuteEdge(); //
        bool cleanGraph(); //
        bool deleteSingle();

        //tools
        unsigned long findFollowingPoint(Point* p, Edge* e); //
        vector<unsigned long> getStartingIdx(); //
        void disconnectEdge(Edge* e); //




        ///under progress
        //fiber correction
        bool fiberProc(void);
        vector< vector<unsigned long> > fiberProc(vector<unsigned long> fib);
};

#endif // C_BUNDLE_H
