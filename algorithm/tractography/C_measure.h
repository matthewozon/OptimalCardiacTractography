#ifndef C_MEASURE_H
#define C_MEASURE_H

#include <vector>
#include <fstream>
#include <C_graph.h>
#include <C_toolbox_integrate.h>
#include <C_toolbox_interpolate.h>
#include <C_toolbox_eigen_sym.h>
#include <C_toolbox_spline.h>
#include <C_toolbox_PCA_3D.h>
#include <C_bundle.h>

#include <C_load_fiber.h>

#include <C_point_array_data.h>
#include <rawData.h>


using namespace std;


class facette
{
public:
    facette()
    {
        p1=-1;
        p2=-1;
        p3=-1;
    }
    facette(int _p1, int _p2, int _p3)
    {
        p1=_p1;
        p2=_p2;
        p3=_p3;
    }
    facette(const facette & f)
    {
        p1 = f.p1;
        p2 = f.p2;
        p3 = f.p3;
    }


    facette(facette & f)
    {
        p1 = f.p1;
        p2 = f.p2;
        p3 = f.p3;
    }

    ~facette(){}

    //point indices
    int p1;
    int p2;
    int p3;
};

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
        double error_fiber_tilde(string fileNameFib);
        double error_fiber(CtlCtlStruct* FIBERS); //rawData<double>* rawTensor,
        double error_fiber_tilde(CtlCtlStruct* FIBERS);


        long getFibersStat(string fileNameEdge, double* L, double* S);
        void getLengthAndBoxDimension(string fileNameFib, vector<double>* L, vector<double>* dU, vector<double>* dV, vector<double>* dW, bool rawCond);
        void getLengthAndBoxDimension(CtlCtlStruct* FIBERS, vector<double>* L, vector<double>* dU, vector<double>* dV, vector<double>* dW);
        void getLengthAndNormalizedTotalPerimeter(CtlCtlStruct* FIBERS, vector<double>* L, vector<double>* totPermim);

        //(dx, dy, dz): short axis, (dfx, dfy, dfz): fiber direction, (xc, yc, zc): origin of heart frame
        void getHelixAndTransversAngle(CtlCtlStruct* FIBERS, double xc, double yc, double zc, double dx, double dy, double dz,  vector<double>* helixAngle, vector<double>* transversAngle, vector<double>* r, vector<double>* theta, vector<double>* z);
        double getHelixAngle(double dx, double dy, double dz, double dfx, double dfy, double dfz);
        //(dux, duy, duz): position vector in the heart frame (translation of (xc, yc, zc))
        double getTransversAngle(double dx, double dy, double dz, double dux, double duy, double duz, double dfx, double dfy, double dfz);

        //CtlCtlStruct* readFiber(string fileNameFib, bool rawData); //return an approximated fiber: B-spline
        rawData<double>* m_rawTensor;
        rawData<double>* m_rawMask;

    protected:
        double error_fiber(CtlStruct* phi); //rawData<double>* rawTensor,
        double error_fiber_tilde(CtlStruct* phi);
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
        void getLengthAndNormalizedTotalPerimeter(CtlStruct* phi, double* L, double* totPerim);
};

#endif // C_MEASURE_H
