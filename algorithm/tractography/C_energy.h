/**C_energy class compute energies on a graph.
It can calculate energies for each vertex of a graph,
or it can compute the energy differnce while flipping a given edge.
Energies are divided into two terms:
    - data fidelity
    - topological constrain
Several kind of energy will be available... but only one now.
Next: energy will depend on a parameter T (temperature)*/


/**
now this class must implement methods that:
    - fill the data field of graph for energy U1 and U15
    - return the energy as it did before.
*/

#ifndef C_ENERGY_H
#define C_ENERGY_H
#include <iostream>
#include <C_graph.h>
#include <rawData.h>
#include <C_point_array_data.h>

#include <C_toolbox_integrate.h>
#include <C_toolbox_interpolate.h>
#include <C_tensorMaker.h>

#include <C_thread_init_energy.h>
#include <C_thread_deltaU_move_vertex.h>

#define F0 0.5L
#define F1 1.0L

#define PHI0 0.5L
#define PHI1 0.5L
#define PHI2 0.0L


using namespace std;

//class C_energy;
class C_thread_deltaU_move_vertex;

#define MY_TEST_MULTITHREAD

#define TEST_FOR_VERTEX_MOVE

class C_SA;
class C_SC;

class C_energy
{
    friend class C_SA;
    friend class C_SC;
    friend class C_thread_deltaU_move_vertex;
    //friend class C_thread_init_energy;
    public:
        /** Default constructor */
        C_energy(C_graph *graph, rawData<double>* _rawData, rawData<double>* _rawDataMask);
        /** Default destructor */
        virtual ~C_energy();

        /** calling function*/
        double deltaEdgeU(Edge* e, double* partial_sum_p1, double* partial_sum_p2, long kind); //ok
        double deltaVertexMoveU(Point* V /**point of the graph*/, Point* newV /**not in graph*/, long kind, vector<long>* idxPoint, vector<vector<long> >* idxPair, vector<vector<double> >* dataPair);
        bool initThreads();
        void deleteThreads();
        //void normalizeTensors(void);

        double getGlobalEnergy(long kind, double* Edata, double* Etopo); //ok
        #ifdef MY_TEST_MULTITHREAD
        void initEnergy(long kind, unsigned long nb_thread);
        #else
        void initEnergy(long kind); //ok
        #endif
        /**to be used for "fast" initialization: this initializes a part of the data array*/
        double initEnergy(unsigned long kind, unsigned long idxStart, unsigned long idxEnd); //ok
        bool initPointData(long kind); ///should be called after initEnergy //ok

        ///attributs
        double m_alpha;
        double m_beta;
        double m_beta0;
        void computeBetaIncrease(double T, double Tinit, double Tend);
        void computeBetaDecrease(double T, double Tinit, double Tend);
        //double m_gamma;
        double m_t0;
        double m_sigma;
        double m_phi0;
        double m_phi1;
        double m_R; //"radius" of the box used for interpolation: physical length
        unsigned long m_nbInterpolatingPoint;
        C_graph *m_graph;
        rawData<double>* m_rawData; ///must be 4D: either DTI if DimT==6 or DWI if DimT>6 (must provide Hpsi and b)
        rawData<double>* m_rawDataMask; ///must be 3D
        double** Hpsi; ///pseudo inverse for DTI computation: matrix 6 rows by DimT-1 columns
        double b;/// b-value for DTI computation
        unsigned long METHOD_INTERPOLATION, OBJECT_INTERPOLATED, METHOD_INTEGRATION; //, unsigned long METHOD_INTERPOLATION/**=GAUSSIAN_PHI*/, unsigned long OBJECT_INTERPOLATED/**=TENSOR*/, unsigned long METHOD_INTEGRATION/**=SIMPSON_RULE*/

    protected:

        C_thread_deltaU_move_vertex** THREAD_DELTA;

        //absolut method
        double U1data(Point * point, Vector * edge1, Vector * edge2);//ok
        double U15data(Edge* e1, Edge* e2);//ok
        double U15topo(Point* p);//ok
        double U16data(Edge* e1, Edge* e2);
        double U16topo(Point* p);

        //variation methods
        double deltaEdgeU1(Edge* e, double* partial_sum_p1, double* partial_sum_p2);//ok
        double deltaU1(Point *point, Vector * edge, double * partial_sum);//ok
        double deltaU1data(Point *point, Vector * edge, double * partial_sum);//ok
        double deltaU1topo(Point *point, Vector * edge);//ok

        double deltaEdgeU15(Edge* e, double* partial_sum_p1, double* partial_sum_p2);//ok
        double deltaU15(Point *point, Vector * edge, double * partial_sum);//ok
        double deltaU15data(Point *point, Vector * edge, double * partial_sum);//ok
        double deltaU15topo(Point *point, Vector * edge);//ok

        double deltaEdgeU16(Edge* e, double* partial_sum_p1, double* partial_sum_p2);
        double deltaU16(Point *point, Vector * edge, double * partial_sum);
        double deltaU16data(Point *point, Vector * edge, double * partial_sum);
        double deltaU16topo(Point *point, Vector * edge);

        //private tools
        double U1dataComputation(Point * point, Vector * edge1, Vector * edge2);
        double U15dataComputation(Point* p1, Point* p2);//ok
        double U16dataComputation(Point* p1, Point* p2);//ok

        double f(double t, double alpha); //ok
        double f2(double t, double t0, double sigma); //ok
        double f3(double t, double t0, double sigma); //ok
        double phi(int n); //ok
        double phi_2(int n); //ok
        double angleFunc(Edge* v1, Edge* v2); ///get cosine of angle between two adjacent edges //ok
        double angleFunc(Point* p1, Point* p2, Point* p12); ///get cosine of angle betwee
        double angleFunc2(Edge* e1, Edge* e2); ///get angle between two adjacent edges //ok

        #ifdef MY_TEST_MULTITHREAD
        double initPairEdgeData(unsigned long kind, unsigned long nb_thread);
        #else
        double initPairEdgeData(unsigned long kind);
        #endif
        double initPairEdgeDataMultiThread(unsigned long kind, unsigned long idxStart, unsigned long idxEnd); //, unsigned long METHOD_INTERPOLATION=GAUSSIAN_PHI, unsigned long OBJECT_INTERPOLATED=TENSOR, unsigned long METHOD_INTEGRATION=SIMPSON_RULE
};

#endif // C_ENERGY_H
