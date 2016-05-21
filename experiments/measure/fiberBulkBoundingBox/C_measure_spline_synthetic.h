#ifndef C_MEASURE_SPLINE_SYNTHETIC_H
#define C_MEASURE_SPLINE_SYNTHETIC_H

#include <C_spline.h>
#include <C_graph.h>
#include <C_toolbox_integrate.h>
#include <C_parametric_comparison.h>
#include <stdio.h>

class C_measure_spline_synthetic
{
    public:
        /** Default constructor */
        C_measure_spline_synthetic();
        /** Default destructor */
        virtual ~C_measure_spline_synthetic();

        void getLenAndErrPop(CtlCtlStruct* FIBERS, bool red); ///the fibers are the raw fibers! interpolation is not yet computed?


        ///attribut about tensor field
        double helix_angle; ///helix angle
        double lambda1;     ///first eigen value
        double lambda2;     ///second eigen value
        double lambda3;     ///third eigen value

    private:
        double fiberError(CtlStruct* phi);
        double fiberLength(CtlStruct* phi);
        void syntheticTensor(double x, double y, double z, double* T);
};

#endif // C_MEASURE_SPLINE_SYNTHETIC_H
