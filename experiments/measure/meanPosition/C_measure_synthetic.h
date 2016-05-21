#ifndef C_MEASURE_SYNTHETIC_H
#define C_MEASURE_SYNTHETIC_H

/**this class need to be checked, but how?
    vectors from tensor can be extracted and compare to the wished directions  //OK
    //FIBERS are made of 30 control points
    create a fiber of length 1 and calculate its length  //OK for numerical method, but analyticals are not as nice as expected and much more time consuming (but multithreading is possible)
    create 30 fibers of length 1 and check std == 0
    creare 30 fibers, 15 of length 1, 15 of length 2 and check std==0.5
    create a circle fiber, radius 10, check if returned value is 10  //
    create a straight fiber, check if returned value tend to 0
    compute the to previous on a fiber population of 15 circle fiber (R=10) and 15 straight line (R=+infty), meanCurvature should be 7.5
    check the compare curvature
    generate a uniform tensor field and two fibers one // one T to the diffusion direction and check the error values
*/

#include <C_spline.h>
#include <C_toolbox_integrate.h>
#include <math.h>
//#include <fstream>
#include <settings_def.h>
//#include <structs.h>
#include <struct.h>

#define TEST

class C_measure_synthetic
{
    public:
        /** Default constructor */
        C_measure_synthetic();
        /** Default destructor */
        virtual ~C_measure_synthetic();

        ///attribut about tensor field
        double helix_angle; ///helix angle
        double lambda1;     ///first eigen value
        double lambda2;     ///second eigen value
        double lambda3;     ///third eigen value

        ///extrinsic feature
        double meanFiberError(CtlCtlStruct* FIBERS);

        ///intrinsic featrue
            ///numerical methods
        double meanFiberLength(CtlCtlStruct* FIBERS);
        double stdFiberLength(CtlCtlStruct* FIBERS, double muL);
        double meanFiberCurvature(CtlCtlStruct* FIBERS);
        double stdMeanFiberCurvature(CtlCtlStruct* FIBERS, double muK);
        double meanSTDFiberCurvature(CtlCtlStruct* FIBERS);
        double meanFiberCurvatureCompare(CtlCtlStruct* FIBERS); ///to be implemented
        double stdMeanFiberCurvatureCompare(CtlCtlStruct* FIBERS); ///to be implemented

            ///analytical methods
        double meanFiberLengthS(CtlCtlStruct* DPHI);
        double stdFiberLengthS(CtlCtlStruct* DPHI, double muL);
        double meanFiberCurvatureS(CtlCtlStruct* FIBERS, CtlCtlStruct* DPHI, CtlCtlStruct* DDPHI);
        double stdMeanFiberCurvatureS(CtlCtlStruct* FIBERS, CtlCtlStruct* DPHI, CtlCtlStruct* DDPHI, double muK);
        //double meanSTDFiberCurvatureS(CtlCtlStruct* FIBERS, CtlCtlStruct* DPHI, CtlCtlStruct* DDPHI);
        double meanFiberCurvatureCompareS(CtlCtlStruct* FIBERS, CtlCtlStruct* DPHI, CtlCtlStruct* DDPHI);
        double stdMeanFiberCurvatureCompareS(CtlCtlStruct* FIBERS, CtlCtlStruct* DPHI, CtlCtlStruct* DDPHI, double muK);
        //double meanSTDFiberCurvatureCompareS(CtlCtlStruct* FIBERS, CtlCtlStruct* DPHI, CtlCtlStruct* DDPHI);

    #ifndef TEST
    private: ///reset to private after test
    #endif
        ///extrinsic feature
        void syntheticTensor(double x, double y, double z, double* T);
        double syntheticCurvature(double x, double y, double z);
        double fiberError(CtlStruct* phi);

        ///intrinsic featrue
            ///numerical methods
        double fiberLength(CtlStruct* phi);
        double fiberCurvature(CtlStruct* phi);
        double stdFiberCurvature(CtlStruct* phi, double muK);
        double fiberCurvatureCompare(CtlStruct* phi);
            ///analytical methods
        double fiberLengthS(CtlStruct* dphi);
        double fiberCurvatureS(CtlStruct* phi, CtlStruct* dphi, CtlStruct* ddphi);
        double fiberCurvatureCompareS(CtlStruct* phi, CtlStruct* dphi, CtlStruct* ddphi);
};

#endif // C_MEASURE_SYNTHETIC_H
