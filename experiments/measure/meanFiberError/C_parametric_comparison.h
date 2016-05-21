#ifndef C_PARAMETRIC_COMPARISON_H
#define C_PARAMETRIC_COMPARISON_H
#include <C_spline.h>
#include <math.h>
#include <C_load_fiber.h>
#include <C_toolbox_integrate.h>

/**
    test:
    test all private method using a single fiber
        1) synthetic fiber
        2) real fiber
    output in files so that you can use it easily with octave to show results
*/

#define TEST
class C_parametric_comparison
{
    public:
        /** Default constructor */
        C_parametric_comparison();
        /** Default destructor */
        virtual ~C_parametric_comparison();

        ///vector field parameter
        double alpha; //angle between horizontal plan (z=0) and the main eigenvector
        double getMeanDist(string fileNameFib, bool spl); ///OK
        double getMeanDist2(string fileNameFib, bool spl); ///OK
        double getMeanDist3(string fileNameFib, bool spl);
        double getMeanDist4(string fileNameFib, bool spl);


        double getMeanDist(CtlCtlStruct* FIBER); ///OK
        double getMeanDist2(CtlCtlStruct* FIBER); ///OK
        double getMeanDist3(CtlCtlStruct* FIBER);//?
        double getMeanDist4(CtlCtlStruct* FIBER);//?

    #ifndef TEST
    private:
    #endif
        double getDist(CtlStruct* FIBER);  ///OK used in getMeanDist and getMeanDist2
        double getDist3(CtlStruct* FIBER); //OK
        double getDist4(CtlStruct* FIBER);
        void getStreamPoint(double t, double r0, double z0, double* x, double* y, double* z);  ///OK
        void getMainVector(double x, double y, double z, double* vx, double* vy, double* vz);
        bool getTimeAndParam(CtlStruct* FIBER, double *r0, double *z0, double *t);  ///OK
        bool timeEstimation(CtlStruct* FIBER, double r0, double z0, double *t);  ///OK
        double radiusEstimation(CtlStruct* FIBER);  ///OK //, double* t, double r0);
        double altitudeEstimation(CtlStruct* FIBER, double* t);  ///OK

        void diffeomorphism(double *t); ///not yet
        //double gfun(double* p, double t, double r0, double z0);  ///OK
        //double hfun(CtlStruct* FIBER, double *t, double r0);  ///OK
        //double hderive(CtlStruct* FIBER, double *t, double r0);  ///OK
        double fiberLength(CtlStruct* phi); ///OK
        double fiberAngle(CtlStruct* phi, double r0, double z0, double* t, double *A);  ///must return the mean angle if everything went fine or NaN if not. A is an array of length equal to the number of element in the fiber and contains angles

};

#endif // C_PARAMETRIC_COMPARISON_H
