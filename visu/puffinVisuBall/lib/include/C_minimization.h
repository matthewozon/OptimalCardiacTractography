#ifndef C_MINIMIZATION_H
#define C_MINIMIZATION_H
#include <settings_def.h>
#include <iostream>
#include <C_graph.h>
#include <C_cooling_schedule.h>
#include <C_energy.h>
#include <C_move.h>
#include <C_toolbox_rand2.h>
#include <vector>
#include <C_toolbox_integrate.h>
#include <C_toolbox_interpolate.h>

#include <fstream>
#include <string>

using namespace std;

extern "C"
{

class C_minimization
{
    public:
        /** Default constructor */
        //C_minimization();
        /** Default destructor */
        virtual ~C_minimization();
        /** Copy constructor
         *  \param other Object to copy from
         */
        /** main ctor */
        //C_minimization(std::string graph_filename);
        C_minimization(C_graph* graph);

        virtual int run(void)=0;// const;
        virtual int iteration(void)=0;

        C_graph* getGraph(void) const {return m_graph;};
        bool saveFiber(std::string fileNameFib, vector< vector<unsigned long> > fibers);

        void randInit(double alphaConnected);
        //void circleInit(double xc, double yc, double zc, double R, bool setToNotConnectedGraph);
        //double sqrDistToSphere(double xc, double yc, double zc, double R, double x, double y, double z);
        //double distToStraightLine(double xu, double yu, double zu, double xc, double yc, double zc, double xa, double ya, double za);
        //void initBeams(double xc, double yc, double zc, double xu, double yu, double zu, bool setToNotConnectedGraph);


    protected:
        C_graph *m_graph;
        std::string m_graph_filename;
        unsigned long* sweep_idx_edges(void);
        unsigned long* sweep_idx_points(void);
        unsigned long* sweep_idx(unsigned long N);
        void flip_edge(Edge * edge,double psum_p1,double psum_p2);
    private:
        void disconnectEdge(Edge* e);
        void shuffle(unsigned long *array, size_t n);

        //bool mustDeleteGraph;


};

}

#endif // C_MINIMIZATION_H
