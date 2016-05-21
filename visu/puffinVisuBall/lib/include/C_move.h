/**this classe draw random edges belonging to a graph, according a given mode.
it can draw only one edge (move1=MOVE_T1)
or two connected edges (move2=MOVE_T2),
or four edge (move4=MOVE_T4)
It also makes the vertices move with a new method: moveVertex
this class shall be used for Simulated Annealing and Stochastic Continuation*/

#ifndef C_MOVE_H
#define C_MOVE_H
#include <iostream>
#include <C_graph.h>
#include <C_toolbox_rand2.h>
#include <vector>


using namespace std;


class C_move
{
    public:
        /** Default constructor */
        C_move(C_graph *graph);
        /** Default destructor */
        virtual ~C_move();
        bool move(long kind); //, Mixture* mix=NULL //here kind can be MOVE_T1, MOVE_T2 or MOVE_T4
        bool moveICM(long kind); //, Mixture* mix=NULL //here kind can be MOVE_T1, MOVE_T2 or MOVE_T4
        bool initMove(long kind); //here kind can be MOVE_T1, MOVE_T2, MOVE_T1+MOVE_T2 or MOVE_T4
        bool initMoveICM(long kind); //here kind can be MOVE_T1, MOVE_T2 or MOVE_T4
        bool endMove(void);
        bool endMoveICM(void);

        //array of Edge address
        MovingEdge* aukMove;
        Mixture *m_mix;
        double T;     ///current temperature
        double Tinit; ///initial temperature
        double Tend;  ///final temperature
        double m_varSigma; ///std of gaussian distribution: vertex move
        bool m_voxelConstraint;
        bool m_moveFromCenter;

    private:

        C_toolbox_rand2 *t;
        C_graph *m_graph;
        bool belongToROI(double x, double y, double z);
        bool move1(void);
        bool move1ICM(void);
        bool move2(void);
        bool move2ICM(void);
        bool move12(bool sa);
        bool move124(void);
        bool move4(void);

        bool moveVertex(void); //how to handle cases for which the proposed position is out of the original pixel? How to know the original pixel?
        bool movev1(bool sa);
        bool movev2(bool sa);

        LandScape myLandScape;
        bool getLandScape2(void);
        bool deleteLandScape(void);

        unsigned long n;
};

#endif // C_MOVE_H
