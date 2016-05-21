#ifndef C_TOOLBOX_INTEGRATE_H
#define C_TOOLBOX_INTEGRATE_H

#include <math.h>
#include <C_toolbox.h>

class C_toolbox_integrate : public C_toolbox
{
    public:
        /** Default constructor */
        C_toolbox_integrate();
        /** Default destructor */
        virtual ~C_toolbox_integrate();
        double run(double* ctrlPoints, unsigned short nbCtrlPoints, double intervalLength, unsigned long METHOD, double* errMax=NULL);
        double runAbs(double** ctrlPoints, unsigned short nbCtrlPoints, double dX, double dY, double dZ, unsigned long METHOD);
    protected:
    private:
        double runMiddleSum(double* ctrlPoints, unsigned short nbCtrlPoints, double subIntervalLength);
        double runMiddleSumAbs(double** ctrlPoints, unsigned short nbCtrlPoints, double dx, double dy, double dz); //only for 3D array
        double runTrapezoidalRule(double* ctrlPoints, unsigned short nbCtrlPoints, double subIntervalLength);
        double runTrapezoidalRuleAbs(double** ctrlPoints, unsigned short nbCtrlPoints, double dx, double dy, double dz); //only for 3D array
        double runSimpson38Rule(double* ctrlPoints, unsigned short nbCtrlPoints, double subIntervalLength);
        double runSimpson38RuleAbs(double** ctrlPoints, unsigned short nbCtrlPoints, double dx, double dy, double dz); //only for 3D array

        //maximum error computation
        double errMaxMiddleSum(double* ctrlPoints, unsigned short nbCtrlPoints, double intervalLength, double subIntervalLength);
        double errTrapezoidalRule(double* ctrlPoints, unsigned short nbCtrlPoints, double intervalLength, double subIntervalLength);
        double errSimpson38Rule(double* ctrlPoints, unsigned short nbCtrlPoints, double intervalLength, double subIntervalLength);

        //maximum of derivative
        double maxF2(double* ctrlPoints, unsigned short nbCtrlPoints, double subIntervalLength);
        double maxF4(double* ctrlPoints, unsigned short nbCtrlPoints, double subIntervalLength);
};

#endif // C_INTEGRATE_H
