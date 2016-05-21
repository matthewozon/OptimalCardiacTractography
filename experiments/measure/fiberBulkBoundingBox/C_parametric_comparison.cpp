#include <C_parametric_comparison.h>

C_parametric_comparison::C_parametric_comparison()
{
    //ctor
}

C_parametric_comparison::~C_parametric_comparison()
{
    //dtor
}

double C_parametric_comparison::getMeanDist(string fileNameFib, bool spl)
{
    C_load_fiber* F = new C_load_fiber();
    CtlCtlStruct* FIBER = F->readFiber(fileNameFib);
    CtlCtlStruct* FIBER_SPLINE;
    if(spl)
    {
        C_spline* SPLINE = new C_spline();
        SPLINE->computeFiberSpline(FIBER, &FIBER_SPLINE, 5, CURVE);
        delete SPLINE;//->~C_spline();
    }
    else
    {
        FIBER_SPLINE = FIBER;
    }

    double m = getMeanDist(FIBER_SPLINE);
    delete FIBER;//->~CtlCtlStruct();
    if(spl)
    {
        delete FIBER_SPLINE;//->~CtlCtlStruct();
    }
    delete F;//->~C_load_fiber();
    return m;
}

double C_parametric_comparison::getMeanDist2(string fileNameFib, bool spl)
{
    C_load_fiber* F = new C_load_fiber();
    CtlCtlStruct* FIBER = F->readFiber(fileNameFib);
    CtlCtlStruct* FIBER_SPLINE;
    if(spl)
    {
        C_spline* SPLINE = new C_spline();
        SPLINE->computeFiberSpline(FIBER, &FIBER_SPLINE, 5, CURVE);
        delete SPLINE;//->~C_spline();
    }
    else
    {
        FIBER_SPLINE = FIBER;
    }

    double m = getMeanDist2(FIBER_SPLINE);
    delete FIBER;//->~CtlCtlStruct();
    if(spl)
    {
        delete FIBER_SPLINE;//->~CtlCtlStruct();
    }
    delete F;//->~C_load_fiber();
    return m;
}

double C_parametric_comparison::getMeanDist(CtlCtlStruct* FIBER)
{
    double ERR = 0.0;
    double temp = 0.0;
    unsigned long b = 0;
    for(unsigned long i=0 ; i<FIBER->N ; i++)
    {
        if(FIBER->fibers[i].N_element>1)
        {
            temp = getDist(&(FIBER->fibers[i]));
        }
        else
        {
            temp = -1.0;
        }

        if(temp>=0.0)
        {
            ERR += temp;
        }
        else
        {
            b++;
        }
    }
    if(b<FIBER->N)
    {
        return ERR/( ((double) (FIBER->N) ) - ((double) b) );
    }
    else
    {
        cout << "F YOU we have: " << FIBER->N << " fibers and " << b << "are fucking wierd" << endl;
        return -1.0;
    }

}

double C_parametric_comparison::getMeanDist2(CtlCtlStruct* FIBER)
{
    double ERR = 0.0;
    double temp = 0.0;
    unsigned long b = 0;
    double L;
    for(unsigned long i=0 ; i<FIBER->N ; i++)
    {
        if(FIBER->fibers[i].N_element>1)
        {
            temp = getDist(&(FIBER->fibers[i]));
            L = fiberLength(&(FIBER->fibers[i]));
        }
        else
        {
            temp = -1.0;
            L = 1.0;
        }

        temp = temp/L;
        if(temp>=0.0)
        {
            ERR += temp;
        }
        else
        {
            b++;
        }
    }
    if(b<FIBER->N)
    {
        return ERR/( ((double) (FIBER->N) ) - ((double) b) );
    }
    else
    {
        return -1.0;
    }

}

double C_parametric_comparison::getMeanDist3(CtlCtlStruct* FIBER)
{
    double ERR = 0.0;
    double temp = 0.0;
    unsigned long b = 0;
    //double L;
    for(unsigned long i=0 ; i<FIBER->N ; i++)
    {
        if(FIBER->fibers[i].N_element>1)
        {
            temp = getDist3(&(FIBER->fibers[i]));
        }
        else
        {
            temp=-1.0;
        }
        //L = fiberLength(&(FIBER->fibers[i]));
        //temp = temp/L;
        if(temp>=0.0)
        {
            ERR += temp;
        }
        else
        {
            b++;
        }
    }
    if(b<FIBER->N)
    {
        return ERR/( ((double) (FIBER->N) ) - ((double) b) );
    }
    else
    {
        return -1.0;
    }

}

double C_parametric_comparison::getMeanDist4(CtlCtlStruct* FIBER)
{
    double ERR = 0.0;
    double temp = 0.0;
    unsigned long b = 0;
    //double L;
    for(unsigned long i=0 ; i<FIBER->N ; i++)
    {
        if(FIBER->fibers[i].N_element>1)
        {
            temp = getDist4(&(FIBER->fibers[i]));
        }
        else
        {
            temp = -1.0;
        }

        //L = fiberLength(&(FIBER->fibers[i]));
        //temp = temp/L;
        if(temp>=0.0)
        {
            ERR += temp;
        }
        else
        {
            b++;
        }
    }
    if(b<FIBER->N)
    {
        return ERR/( ((double) (FIBER->N) ) - ((double) b) );
    }
    else
    {
        return -1.0;
    }

}

double C_parametric_comparison::getDist(CtlStruct* FIBER)
{
    try
    {
        double r0, z0;
        double *t = new double[FIBER->N_element];
        if(getTimeAndParam(FIBER, &r0, &z0, t))
        {
            ///
            double x, y, z;
            double m = 0.0;
            for(unsigned long k=0 ; k<FIBER->N_element ; k++)
            {
                getStreamPoint(t[k], r0, z0, &x, &y, &z);
                m +=  SQR(FIBER->elts[3*k]-x) + SQR(FIBER->elts[3*k+1]-y) + SQR(FIBER->elts[3*k+2]-z);
            }
            delete t;
            return m/((double) (FIBER->N_element));
        }
        else
        {
            //cout << -1.0 << endl;
            delete t;
            return -1.0;
        }
    }
    catch(...)
    {
        //cout << "-2" << endl;
        return -2.0;
    }
}

double C_parametric_comparison::getDist3(CtlStruct* FIBER)
{
    try
    {
        double r0, z0;
        double *t = new double[FIBER->N_element];
        if(getTimeAndParam(FIBER, &r0, &z0, t))
        {
            ///
            double x, y, z;
            double m = 0.0;
            //double* A = new double[FIBER->N_element];
            //double Amean = 0.1 + (9.9*2.0/pi)*fiberAngle(FIBER, r0, z0, t, NULL/**A*/); ///to be modified
            double Amean = fiberAngle(FIBER, r0, z0, t, NULL/**A*/);
            if(Amean!=Amean)
            {
                //delete A;
                delete t;
                throw -1.0;
            }
            for(unsigned long k=0 ; k<FIBER->N_element ; k++)
            {
                getStreamPoint(t[k], r0, z0, &x, &y, &z);
                m += (SQR(FIBER->elts[3*k]-x) + SQR(FIBER->elts[3*k+1]-y) + SQR(FIBER->elts[3*k+2]-z));
            }
            m = m/((double) (FIBER->N_element));
            delete t;
            //delete A;

            ///test
            m = sin(Amean)*m;
            return m/fiberLength(FIBER);///sin(E[angle])*E[dist]/L
        }
        else
        {
            //cout << -1.0 << endl;
            delete t;
            return -1.0;
        }
    }
    catch(...)
    {
        //cout << "-2" << endl;
        return -2.0;
    }
}

double C_parametric_comparison::getDist4(CtlStruct* FIBER)
{
    try
    {
        double r0, z0;
        double *t = new double[FIBER->N_element];
        if(getTimeAndParam(FIBER, &r0, &z0, t))
        {
            double Amean = fiberAngle(FIBER, r0, z0, t, NULL);
            if(Amean!=Amean)
            {
                //delete A;
                delete t;
                throw -1.0;
            }

            ///test
            delete t;
            return sin(Amean);//Amean
        }
        else
        {
            //cout << -1.0 << endl;
            delete t;
            return -1.0;
        }
    }
    catch(...)
    {
        //cout << "-2" << endl;
        return -2.0;
    }
}

double C_parametric_comparison::fiberAngle(CtlStruct* phi, double r0, double z0, double* t, double *A)
{
    if(phi->N_element<=1)
    {
        return -1.0;
    }
    double a = 0.0, temp;
    double x, y, z;
    double vx, vy, vz;
    double dx, dy, dz, dn;
    ///

    for(unsigned long k=0 ; k<phi->N_element ; k++)
    {
        getStreamPoint(t[k], r0, z0, &x, &y, &z);
        getMainVector(x, y, z, &vx, &vy, &vz);
        if( (ABS(vx) + ABS(vy) + ABS(vz))>0.1)
        {
            if(k==0)
            {
                dx = phi->elts[3*(k+1)] - phi->elts[3*k];
                dy = phi->elts[3*(k+1)+1] - phi->elts[3*k+1];
                dz = phi->elts[3*(k+1)+2] - phi->elts[3*k+2];
            }
            else if(k==phi->N_element-1)
            {
                dx = phi->elts[3*k] - phi->elts[3*(k-1)];
                dy = phi->elts[3*k+1] - phi->elts[3*(k-1)+1];
                dz = phi->elts[3*k+2] - phi->elts[3*(k-1)+2];
            }
            else
            {
                dx = phi->elts[3*(k+1)] - phi->elts[3*(k-1)];
                dy = phi->elts[3*(k+1)+1] - phi->elts[3*(k-1)+1];
                dz = phi->elts[3*(k+1)+2] - phi->elts[3*(k-1)+2];
            }
            dn = sqrt(SQR(dx) + SQR(dy) + SQR(dz));
            dx = dx/dn;
            dy = dy/dn;
            dz = dz/dn;

            if(A!=NULL)
            {
                A[k] = acos(dx*vx + dy*vy + dz*vz);
                if(A[k]>0.5*Pi)
                {
                    A[k] = Pi - A[k];
                }
                a += A[k];
            }
            else
            {
                temp = acos(dx*vx + dy*vy + dz*vz);
                if(temp>0.5*Pi)
                {
                    temp = Pi - temp;
                }
                a += temp;
            }

        }
        else
        {
            A[k] = NAN;
            a = NAN;
        }
    }
    return a/((double) phi->N_element);
}

void C_parametric_comparison::getStreamPoint(double t, double r0, double z0, double* x, double* y, double* z)
{
    if(r0!=0.0)
    {
        *x = r0*cos(t*cos(alpha)/r0);
        *y = -r0*sin(t*cos(alpha)/r0);
    }
    else
    {
        *x = 0.0;
        *y = 0.0;
    }
    *z = t*sin(alpha) + z0;
    return;
}

void C_parametric_comparison::getMainVector(double x, double y, double z, double* vx, double* vy, double* vz)
{
    if(sqrt(SQR(x)+SQR(y))!=0.0)
    {
        *vx = cos(alpha)*y/sqrt(SQR(x)+SQR(y));
        *vy = -cos(alpha)*x/sqrt(SQR(x)+SQR(y));
        *vz = sin(alpha);
    }
    else
    {
        *vx = 0.0;
        *vy = 0.0;
        *vz = 0.0;
    }
    return;
}

bool C_parametric_comparison::getTimeAndParam(CtlStruct* FIBER, double *r0, double *z0, double *t)
{
    if(FIBER->N_element==0)
    {
        return false;
    }
    try
    {
        ///init

        *z0 = 0;

        ///radius estimation does not need the knwledge of time
        *r0 = radiusEstimation(FIBER);

        ///iter
        for(unsigned short k=0 ; k<50 ; k++)
        {
            //time estimation first
            if(!timeEstimation(FIBER, *r0, *z0, t))
            {
                throw 987;
            }

            //parameters estimation
            *z0 = altitudeEstimation(FIBER, t);
        }
    }
    catch(...)
    {
        ///what else?
        return false;
    }

    return true;
}
bool C_parametric_comparison::timeEstimation(CtlStruct* FIBER, double r0, double z0, double *t)
{
    if(FIBER->N_element==0)
    {
        return false;
    }

    if(false)
    {
//        try
//        {
//            double* P = new double[3];
//            for(unsigned short i=0 ; i<FIBER->N_element ; i++)
//            {
//                ///for ith point, init time
//                P[0] = ((double) (FIBER->elts[3*i]));
//                P[1] = ((double) (FIBER->elts[3*i+1]));
//                P[2] = ((double) (FIBER->elts[3*i+2]));
//                if(cos(alpha)==0.0)
//                {
//                    if(i==0)
//                    {
//                        t[i] = 0.0;
//                    }
//                    else
//                    {
//                        t[i] = t[i-1];
//                    }
//                }
//                else
//                {
//                    t[i] = -(r0/cos(alpha)) * atan2(P[1], P[0]);
//                }
//
//                ///compute algorithm
//                for(unsigned short k=0 ; k<500 ; k++)
//                {
//                    t[i] = gfun(P, t[i], r0, z0);
//                }
//            }
//            delete P;
//        }
//        catch(...)
//        {
//            return false;
//        }
    }
    else
    {
        try
        {
            if(cos(alpha)==0.0)
            {
                for(unsigned short i=0 ; i<FIBER->N_element ; i++)
                {
                    t[i] = ((double) (FIBER->elts[3*i+2]))-z0;
                }
            }
            else if(sin(alpha)==0.0)
            {
                for(unsigned short i=0 ; i<FIBER->N_element ; i++)
                {
                    ///for ith point, init time
                    t[i] = -r0*atan2(((double) (FIBER->elts[3*i+1])), ((double) (FIBER->elts[3*i])));
                }
            }
            else
            {
                double* P = new double[3];
                double k;
                for(unsigned short i=0 ; i<FIBER->N_element ; i++)
                {
                    ///for ith point, init time
                    P[0] = ((double) (FIBER->elts[3*i]));
                    P[1] = ((double) (FIBER->elts[3*i+1]));
                    P[2] = ((double) (FIBER->elts[3*i+2]));
                    k = (0.5/Pi)*(atan2(P[1], P[0]) + (P[2]-z0)/(r0*tan(alpha)));
                    k = floor(k+0.5);
                    t[i] = (r0/cos(alpha)) * ( 2.0*k*Pi - atan2(P[1], P[0]) );
                }
                delete P;
            }

        }
        catch(...)
        {
            return false;
        }
    }


    return true;
}
double C_parametric_comparison::radiusEstimation(CtlStruct* FIBER)//, double* t, double r0)
{
    //double r = r0;
    //double H, dH;

    double r;
    try
    {
        r = 0.0;
        for(unsigned long i=0 ; i<FIBER->N_element ; i++)
        {
            r += sqrt(SQR(FIBER->elts[3*i]) + SQR(FIBER->elts[3*i+1]));
        }
        r = r/((double) FIBER->N_element);
//        for(unsigned long i=0 ; i<50 ; i++)
//        {
//            H = hfun(FIBER, t, r);
//            if(H!=H)
//            {
//                //cout << "H is NaN" << endl;
//                throw 741;
//            }
//            if(H==0.0)
//            {
//                i=50;
//            }
//            else
//            {
//                dH = hderive(FIBER, t, r);
//                if(dH!=dH)
//                {
//                    //cout << "dH is NaN" << endl;
//                    throw 742;
//                }
//                if(dH==0.0)
//                {
//                    //cout << "dH is NULL" << endl;
//                    throw 654; //devide by 0!!!
//                }
//                r = r - (H/dH);
//            }
//        }
    }
    catch(...)
    {
        r = -1.0;
    }
    return r;
}
double C_parametric_comparison::altitudeEstimation(CtlStruct* FIBER, double* t)
{
    double z0 = 0.0;
    try
    {
        for(unsigned long i=0 ; i<FIBER->N_element ; i++)
        {
            z0 += (FIBER->elts[3*i+2] - t[i]*sin(alpha));
        }
    }
    catch(...)
    {
        z0 = 0.0;
    }

    return z0/((double) FIBER->N_element);
}

void C_parametric_comparison::diffeomorphism(double *t) ///nothing yet, leave it like it is!
{
    return;
}
//double C_parametric_comparison::gfun(double* P, double t, double r0, double z0)
//{
//    double R = sqrt(SQR(P[0]) + SQR(P[1]));
//    double xi = t*cos(alpha)/r0;
//    double PSI = atan2(P[1],P[0]);
//    return (1.0/sin(alpha)) * (P[2] - z0 -  (R/tan(alpha))*sin(xi + PSI));
//}
//double C_parametric_comparison::hfun(CtlStruct* FIBER, double *t, double r0)
//{
//    double Rf = 0.0;
//    double PSI = 0.0;
//    double xi = 0.0;
//    double h = 0.0;
//    for(unsigned long i=0 ; i<FIBER->N_element ; i++)
//    {
//        Rf = sqrt(SQR(FIBER->elts[3*i]) + SQR(FIBER->elts[3*i+1]));
//        PSI = atan2(FIBER->elts[3*i+1],FIBER->elts[3*i]);
//        xi = t[i]*cos(alpha)/r0;
//        h += (r0 - Rf*(cos(xi+PSI) + xi*sin(xi+PSI)));
//    }
//    return h;
//}
//double C_parametric_comparison::hderive(CtlStruct* FIBER, double *t, double r0)
//{
//    double Rf = 0.0;
//    double PSI = 0.0;
//    double xi = 0.0;
//    double h = 0.0;
//    for(unsigned long i=0 ; i<FIBER->N_element ; i++)
//    {
//        Rf = sqrt(SQR(FIBER->elts[3*i]) + SQR(FIBER->elts[3*i+1]));
//        PSI = atan2(FIBER->elts[3*i+1],FIBER->elts[3*i]);
//        xi = t[i]*cos(alpha)/r0;
//        h += (1 - (Rf/r0)*(xi*sin(xi+PSI) - SQR(xi)*cos(xi+PSI)));
//    }
//    return h;
//}

double C_parametric_comparison::fiberLength(CtlStruct* phi)
{
    if(phi->N_element<=1)
    {
        return 0.0;
    }
    if(phi->N_element==2)
    {
        return sqrt(SQR(phi->elts[0]-phi->elts[3]) + SQR(phi->elts[0+1]-phi->elts[3+1]) + SQR(phi->elts[0+2]-phi->elts[3+2]));//3*i+
    }

    double *DPHI = new double[phi->N_element];
    if(DPHI==NULL)
    {
        return -1.0;
    }
    double ds = 2.0/((double) (phi->N_element-1));
    DPHI[0] = sqrt(SQR(phi->elts[3*1]-phi->elts[3*0]) + SQR(phi->elts[3*1+1]-phi->elts[3*0+1]) + SQR(phi->elts[3*1+2]-phi->elts[3*0+2]))/(ds/2.0);
    for(unsigned long i=1 ; i<phi->N_element-1 ; i++)
    {
        DPHI[i] = sqrt(SQR(phi->elts[3*(i+1)]-phi->elts[3*(i-1)]) + SQR(phi->elts[3*(i+1)+1]-phi->elts[3*(i-1)+1]) + SQR(phi->elts[3*(i+1)+2]-phi->elts[3*(i-1)+2]))/(ds);
    }
    DPHI[phi->N_element-1] = sqrt(SQR(phi->elts[3*(phi->N_element-1)]-phi->elts[3*(phi->N_element-2)]) + SQR(phi->elts[3*(phi->N_element-1)+1]-phi->elts[3*(phi->N_element-2)+1]) + SQR(phi->elts[3*(phi->N_element-1)+2]-phi->elts[3*(phi->N_element-2)+2]) )/(ds/2.0);

    ///integrate the norm of derivative parametric curve over the whole arc
    C_toolbox_integrate* t = new C_toolbox_integrate();

    double L = t->run(DPHI, phi->N_element, 1.0, SIMPSON_RULE, NULL);

    delete t;//->~C_integrate();
    delete DPHI;
    return L;
}
