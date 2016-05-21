#ifndef C_TOOLBOX_H
#define C_TOOLBOX_H

#include <settings_def.h>
#include <vector>
#include <math.h>
#include <stdlib.h>
#include <MACRO_DEF.h>

#include <iostream> ///must be inclded
using namespace std; ///must stay there!!!

class C_toolbox
{
    public:
        /** Default constructor */
        C_toolbox();
        /** Default destructor */
        virtual ~C_toolbox()=0;


        //tool
        void reshapeSymMatrix(double* A, double** a); //top level
        double FA(double eig1, double eig2, double eig3); //top level
        double getAngle(double vx1, double vy1, double vz1, double vx2, double vy2, double vz2);
        double sqrDistToSphere(double xc, double yc, double zc, double R, double x, double y, double z); //top level
        double distToStraightLine(double xu, double yu, double zu, double xc, double yc, double zc, double xa, double ya, double za); //top level
};

#endif // C_TOOLBOX_H
