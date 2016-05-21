#include "C_toolbox.h"

C_toolbox::C_toolbox()
{
    //ctor
}

C_toolbox::~C_toolbox()
{
    //dtor
}



//tool
void C_toolbox::reshapeSymMatrix(double* A, double** a) //3D
{
    a[0][0] = A[0];
    a[0][1] = A[1];
    a[0][2] = A[2];
    a[1][0] = A[1];
    a[1][1] = A[3];
    a[1][2] = A[4];
    a[2][0] = A[2];
    a[2][1] = A[4];
    a[2][2] = A[5];
    /* symmetric matrix format:
	 * ( 0 1 2 )
	 * ( 1 3 4 )
	 * ( 2 4 5 ) */
}



double C_toolbox::FA(double eig1, double eig2, double eig3)
{
    return sqrt( 0.5 * ( SQR( eig1 - eig2 ) + SQR( eig2 - eig3 ) + SQR( eig3 - eig1 ) ) / (  SQR(eig1) + SQR(eig2) + SQR(eig3) ) );
}


double C_toolbox::getAngle(double vx1, double vy1, double vz1, double vx2, double vy2, double vz2)
{
    double psi = acos((vx1*vx2 + vy1*vy2 + vz1*vz2)/(sqrt(SQR(vx1) + SQR(vy1) + SQR(vz1))*sqrt(SQR(vx2) + SQR(vy2) + SQR(vz2))));
    if(psi>Pi/2)
    {
        return Pi-psi;
    }
    else
    {
        return psi;
    }
}


double C_toolbox::sqrDistToSphere(double xc, double yc, double zc, double R, double x, double y, double z)
{
    return SQR( ABS(1 - ( R / (sqrt( SQR(xc-x) + SQR(yc-y) + SQR(zc-z)))))) * ( SQR(xc-x) + SQR(yc-y) + SQR(zc-z) );
}
double C_toolbox::distToStraightLine(double xu, double yu, double zu, double xc, double yc, double zc, double xa, double ya, double za)
{
    return ( sqrt( SQR((ya-yc)*zu - (za-zc)*yu) + SQR((za-zc)*xu - (xa-xc)*zu) + SQR((xa-xc)*yu - (ya-yc)*xu) ) ) / sqrt(xu*xu + yu*yu + zu*zu);
}


