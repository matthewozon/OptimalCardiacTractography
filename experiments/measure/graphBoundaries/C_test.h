#ifndef C_TEST_H
#define C_TEST_H

#include <C_measure_synthetic.h>
#include <stdio.h>
#include <C_toolbox_eigen_sym.h>
#include <C_spline.h>



#include <C_parametric_comparison.h>
#include <C_measure_spline_synthetic.h>

#include <C_crop.h>

class C_test
{
    public:
        /** Default constructor */
        C_test();
        /** Default destructor */
        virtual ~C_test();

        void testSynthetic(void);
        void testSpline(void);
        void testFiberLength(void);
        void testFiberCurvature(void);
        void testFiberCurvatureCompare(void);
        void testFiberError(void);

        void testComparisonStreamPoint(void);
        void testComparisonRadius(void);
        void testComparisonTime(void);
        void testComparisonTimeAndParam(void);
        void testComparisonDist(void);
        void testComparisonTimeAndParamReal(void);
        void testComparisonDistReal(void);

        void testHisto(void);

        #ifdef CROP
        void testCrop(void);
        #endif
};

#endif // C_TEST_H
