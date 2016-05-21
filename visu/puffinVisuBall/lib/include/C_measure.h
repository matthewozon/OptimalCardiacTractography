#ifndef C_MEASURE_H
#define C_MEASURE_H

#include <vector>
#include <fstream>
#include <C_graph.h>
#include <C_toolbox_integrate.h>
//#include <C_interpolate.h>
#include <C_toolbox_interpolate.h>
#include <C_toolbox_eigen_sym.h>
#include <C_toolbox_spline.h>
#include <C_toolbox_PCA_3D.h>
#include <C_bundle.h>

#include <C_load_fiber.h>

//#include <C_measure_spline.h>
//#include <C_measure_synthetic.h>

#include <C_point_array_data.h>
#include <rawData.h>


using namespace std;

class C_measure
{
    public:
        /** Default constructor */
        C_measure();
        /** Default destructor */
        virtual ~C_measure();

        void splinMeasure(string fileNameFib, double* mean, double* sigma, double* Etild2, double* E2);//rawData<double>* rawTensor,
        void splinMeasure(CtlCtlStruct* FIBERS, double* mean, double* sigma, double* Etild2, double* E2);//rawData<double>* rawTensor,
        double error_fiber(string fileNameFib); //rawData<double>* rawTensor,


        long getFibersStat(string fileNameEdge, double* L, double* S);
        void getLengthAndBoxDimension(string fileNameFib, vector<double>* L, vector<double>* dU, vector<double>* dV, vector<double>* dW, bool rawCond);
        void getLengthAndBoxDimension(CtlCtlStruct* FIBERS, vector<double>* L, vector<double>* dU, vector<double>* dV, vector<double>* dW);

        //CtlCtlStruct* readFiber(string fileNameFib, bool rawData); //return an approximated fiber: B-spline
        rawData<double>* m_rawTensor;
        rawData<double>* m_rawMask;

    protected:

        double error_fiber(CtlCtlStruct* FIBERS); //rawData<double>* rawTensor,
        double error_fiber(CtlStruct* phi); //rawData<double>* rawTensor,
        double len_fib(CtlStruct* phi);



        //compute fiber length (as define in (5.23))
        double LSplin(CtlStruct* phi);
        double lengthError(CtlStruct* phi, double* E2); //rawData<double>* rawTensor,

        C_toolbox_integrate* tool;
        void fillArray(double x, double y, double z, double** a);
        double TphiPrimeNorm2(double x, double y, double z, double dx, double dy, double dz);
        double FA(double* eig);

        /*interpolatePoints**/ C_point_array_data* fillInterpolatePoints(CtlStruct* phi); //rawData<double>* rawTensor,

        double meanLengthSplin(CtlCtlStruct* spline);
        double stdLengthSplin(CtlCtlStruct* spline, double mean);
        double meanLengthErrorSplinE2(CtlCtlStruct* spline, double* E2);///rawData<double>* rawTensor,  /// the same as meanLengthError computed over approximation spline

        void getLengthAndBoxDimension(CtlStruct* phi, double* L, double* dU, double* dV, double* dW);
};

#endif // C_MEASURE_H
