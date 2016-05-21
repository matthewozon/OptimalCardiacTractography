#ifndef C_TEMPERATURE_H
#define C_TEMPERATURE_H

#include <C_graph.h>
#include <C_energy.h>
#include <C_toolbox_rand2.h>
#include <C_move.h>

class C_temperature
{
    public:
        /** Default constructor */
        C_temperature(C_graph* graph, C_energy* U, unsigned short energy_kind);
        /** Default destructor */
        virtual ~C_temperature();

        //seed for random number generation (may be modified during a process, but I never did that)
        long idum;

        //
        double getTinit(void) const {return Tinit;}
        double getTend(void) const {return Tend;}

        //compute the estimation of the temperatures
        bool runTemperature(unsigned long kind_move_init/** = MOVE_T1*/, unsigned long kind_move_end/** = MOVE_T1*/, double tauInit/**=0.8*/, double tauEnd/**=0.0005*/, unsigned long M/**=1000*/);

        //the parameters that characterize the communication kernel
        double m_alphaMixture;
        double m_varSigma;
        double m_varSigma_end;
        bool m_voxelConstraint;
        bool m_moveFromCenter;

    protected:
        //the variables that contain the initial and final temperatures
        double Tinit;
        double Tend;

    private:

        unsigned short m_energy_kind;
        bool startingMove; //only for moving vertex with varying std: SC_MOVE_V1
        C_energy* m_U; //this is a pointer to the actual energy (that is already initialized)
        C_toolbox_rand2* m_rand;
        C_graph* m_graph;


        //root finder: Raphson-Newton
        double rootFinder(double* deltaU /*an array of positive value*/, unsigned long M /*number of elements in deltaU*/, double tau /*acceptance rate*/, double epsilon=0.00001 /*precision on beta*/);
        double fxtau(double beta, double* deltaU, double M, double tau);
        double fxtau_deriv(double beta, double* deltaU, double M);

        //generation of uphill moves
        double* uphill_energy(unsigned long kind_move, bool init, unsigned long* N);
        double uphill_energyM1(Edge** edges, bool init);
        double uphill_energyM2(Edge** edges, bool init);
        double uphill_energyM4(Edge** edges, bool init);
        double uphill_energyM1M2(MovingEdge* edges, bool init);
        double uphill_energyP1(long idxPoint, Point* newV);
        void flip_edge(Edge * edge,double psum_p1,double psum_p2);
};

#endif // C_TEMPERATURE_H
