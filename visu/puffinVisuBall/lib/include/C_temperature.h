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

        long idum;

        //
        double getTinit(void) const {return Tinit;}
        double getTend(void) const {return Tend;}

        //
        bool runTemperature(unsigned long kind_move_init/** = MOVE_T1*/, unsigned long kind_move_end/** = MOVE_T1*/, double tauInit/**=0.8*/, double tauEnd/**=0.0005*/, unsigned long M/**=1000*/);
        //bool runTemperature2(unsigned long kind_move_init/** = MOVE_T1*/, unsigned long kind_move_end/** = MOVE_T1*/, double tauInit/**=0.8*/, double tauEnd/**=0.0005*/, unsigned long M/**=1000*/); //should not be used
        double m_alphaMixture;
        double m_varSigma;
        double m_varSigma_end;
        bool m_voxelConstraint;
        bool m_moveFromCenter;

        //
//        double upBound;
//        double lowBound;
//        double tauScale;

    protected:
        double Tinit;
        double Tend;

    private:

        unsigned short m_energy_kind;
        bool startingMove; //only for moving vertex with varying std: SC_MOVE_V1
        //unsigned long CONCAVE;
        //double globU;
        //double m_alpha;
        C_energy* m_U;
        C_toolbox_rand2* m_rand;
        C_graph* m_graph;

        //
        double new_temp_by_secant(Edge * edges, int n_edges, double accept_max, double T_0);
        double annealing_sampling(Edge * edges, int n_edges, /*double alpha_topo, */double temperature,int tries);

        //new
        double rootFinder(double* deltaU /*an array of positive value*/, unsigned long M /*number of elements in deltaU*/, double tau /*acceptance rate*/, double epsilon=0.00001 /*precision on beta*/);
        double fxtau(double beta, double* deltaU, double M, double tau);
        double fxtau_deriv(double beta, double* deltaU, double M);
        double* uphill_energy(unsigned long kind_move, bool init, unsigned long* N);
        double uphill_energyM1(Edge** edges, bool init);
        double uphill_energyM2(Edge** edges, bool init);
        double uphill_energyM4(Edge** edges, bool init);
        double uphill_energyM1M2(MovingEdge* edges, bool init);
        double uphill_energyP1(long idxPoint, Point* newV);
        void flip_edge(Edge * edge,double psum_p1,double psum_p2);

//        ///new assessment methods
//        double temperatureAssessement(unsigned long kind_move, unsigned long M, double tau, double* tauSup);
//        double* uphill(unsigned long kind_move, double tauSup, unsigned long* N, unsigned long* NICM);
//        double uphillM1(Edge** edges, double tauSup);
//        double uphillM2(Edge** edges, double tauSup);
//        double uphillM4(Edge** edges, double tauSup);
//        double uphillM1M2(MovingEdge* edges, double tauSup);
//        void initGraph(unsigned long kind_move, double tauSup);
};

#endif // C_TEMPERATURE_H
