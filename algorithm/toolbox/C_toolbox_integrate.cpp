#include <C_toolbox_integrate.h>

C_toolbox_integrate::C_toolbox_integrate():C_toolbox()
{
    //ctor
}

C_toolbox_integrate::~C_toolbox_integrate()
{
    //dtor
}
double C_toolbox_integrate::run(double* ctrlPoints, unsigned short nbCtrlPoints, double intervalLength, unsigned long METHOD, double* errMax)
{
    if(ctrlPoints!=0 && nbCtrlPoints>1 && intervalLength>0)
    {

        if(METHOD==MIDDLE_SUM)
        {
            //control points must be in the middle of each interval. Means there are as much points as intervals
            double subIntervalLength = intervalLength/((double) nbCtrlPoints);
            if(errMax!=NULL)
            {
                *errMax = errMaxMiddleSum(ctrlPoints, nbCtrlPoints, intervalLength, subIntervalLength);
            }
            return runMiddleSum(ctrlPoints,nbCtrlPoints,subIntervalLength);
        }
        else if(METHOD == TRAPEZOIDAL_RULE)
        {
            double subIntervalLength = intervalLength/((double) nbCtrlPoints-1.0);
            if(errMax!=NULL)
            {
                *errMax = errTrapezoidalRule(ctrlPoints, nbCtrlPoints, intervalLength, subIntervalLength);
            }
            return runTrapezoidalRule(ctrlPoints,nbCtrlPoints,subIntervalLength);
        }
        else if(METHOD == SIMPSON_RULE)
        {
            double subIntervalLength = intervalLength/((double) nbCtrlPoints-1.0);
            if(errMax!=NULL)
            {
                *errMax = errSimpson38Rule(ctrlPoints, nbCtrlPoints, intervalLength, subIntervalLength);
            }
            return runSimpson38Rule(ctrlPoints,nbCtrlPoints,subIntervalLength);
        }
        else
        {
            if(errMax!=NULL)
            {
                *errMax = NAN;
            }
            return 0; // maybe add boolean tag to say that it failed?
        }
    }
    else
    {
        return 0;
    }

}
double C_toolbox_integrate::runAbs(double** ctrlPoints, unsigned short nbCtrlPoints, double dX, double dY, double dZ, unsigned long METHOD)
{
    if(ctrlPoints!=0 && nbCtrlPoints>1 && (ABS(dX)>0 || ABS(dY)>0 || ABS(dZ)>0) )
    {
        //std::cout << "computing integral" << std::endl;
        if(METHOD==MIDDLE_SUM)
        {
            //control points must be in the middle of each interval. Means there are as much points as intervals
            return runMiddleSumAbs(ctrlPoints, nbCtrlPoints, dX/((double) nbCtrlPoints), dY/((double) nbCtrlPoints), dZ/((double) nbCtrlPoints));
        }
        else if(METHOD == TRAPEZOIDAL_RULE)
        {
            return runTrapezoidalRuleAbs(ctrlPoints, nbCtrlPoints, dX/((double) nbCtrlPoints-1.0), dY/((double) nbCtrlPoints-1.0), dZ/((double) nbCtrlPoints-1.0));
        }
        else if(METHOD == SIMPSON_RULE)
        {
            return runSimpson38RuleAbs(ctrlPoints, nbCtrlPoints, dX/((double) nbCtrlPoints-1.0), dY/((double) nbCtrlPoints-1.0), dZ/((double) nbCtrlPoints-1.0));
        }
        else
        {
            return 0; // maybe add boolean tag to say that it failed?
        }
    }
    else
    {
        throw 0.0;
        return 0;
    }

}
double C_toolbox_integrate::runMiddleSum(double* ctrlPoints, unsigned short nbCtrlPoints, double subIntervalLength)
{
    double sum = 0.0;
    for(unsigned short i=0 ; i<nbCtrlPoints ; i++)
        sum += ctrlPoints[i];
    return subIntervalLength*sum;
}
double C_toolbox_integrate::runTrapezoidalRule(double* ctrlPoints, unsigned short nbCtrlPoints, double subIntervalLength)
{
    double sum = ctrlPoints[0] + ctrlPoints[nbCtrlPoints-1];
    for(unsigned short i=1 ; i<nbCtrlPoints-1 ; i++)
        sum += 2*ctrlPoints[i];
    return subIntervalLength*sum/2.0;
}
double C_toolbox_integrate::runSimpson38Rule(double* ctrlPoints, unsigned short nbCtrlPoints, double subIntervalLength)
{
    //check nb points: must be multiple of 3 (if not, there are two cases to do +1 trapezoidal rule or +2 another)
    unsigned short m = (unsigned short) (floor( ((double) (nbCtrlPoints-1))/3.0));
    unsigned short m_rest = nbCtrlPoints - (3*m+1);
    double sum = 0.0;

    //compute 3/8 simpson rule on the 3m first points
    for(unsigned short i=0 ; i<m ; i++)
        sum += ctrlPoints[3*i] + 3.0*ctrlPoints[3*i+1] + 3.0*ctrlPoints[3*i+2] + ctrlPoints[3*i+3];

    sum = (3.0*subIntervalLength/8.0) * sum;

    //so that it bind the last values to the integral
    if(m_rest==1)
    {
        //trapezoidal case
        sum += (subIntervalLength/2.0) * (ctrlPoints[nbCtrlPoints-2] + ctrlPoints[nbCtrlPoints-1]);
    }
    else if(m_rest==2) //means s==2
    {
        //trapezoidal case 2
        sum += (subIntervalLength/2.0) * (ctrlPoints[nbCtrlPoints-3] + 2*ctrlPoints[nbCtrlPoints-2] + ctrlPoints[nbCtrlPoints-1]);
    }

    return sum;
}

double C_toolbox_integrate::runMiddleSumAbs(double** ctrlPoints, unsigned short nbCtrlPoints, double dx, double dy, double dz)
{
    double sum = 0.0;
    //std::cout << "dx: " << dx << "dy: " << dy << "dz: " << dz << std::endl;
    for(unsigned short i=0 ; i<nbCtrlPoints ; i++)
        sum += ABS(dx*ctrlPoints[0][i] + dy*ctrlPoints[1][i] + dz*ctrlPoints[2][i]);
    return sum;
}
double C_toolbox_integrate::runTrapezoidalRuleAbs(double** ctrlPoints, unsigned short nbCtrlPoints, double dx, double dy, double dz)
{
    double sum = ABS( dx*(ctrlPoints[0][0] + ctrlPoints[0][nbCtrlPoints-1]) + dy*(ctrlPoints[1][0] + ctrlPoints[1][nbCtrlPoints-1]) + dz*(ctrlPoints[2][0] + ctrlPoints[2][nbCtrlPoints-1]) ) / 2.0;
    for(unsigned short i=1 ; i<nbCtrlPoints-1 ; i++)
        sum += ABS( dx*ctrlPoints[0][i] + dy*ctrlPoints[1][i] + dz*ctrlPoints[2][i] );
    return sum;
}
double C_toolbox_integrate::runSimpson38RuleAbs(double** ctrlPoints, unsigned short nbCtrlPoints, double dx, double dy, double dz)
{
    //check nb points: must be multiple of 3 (if not, there are two cases to do +1 trapezoidal rule or +2 another)
    unsigned short m = (unsigned short) (floor( ((double) (nbCtrlPoints-1))/3.0));
    unsigned short m_rest = nbCtrlPoints - (3*m+1);
    double sum = 0.0;

    //compute 3/8 simpson rule on the 3m first points
    for(unsigned short i=0 ; i<m ; i++)
        sum += ABS(\
        dx*(ctrlPoints[0][3*i] + 3.0*ctrlPoints[0][3*i+1] + 3.0*ctrlPoints[0][3*i+2] + ctrlPoints[0][3*i+3]) +\
        dy*(ctrlPoints[1][3*i] + 3.0*ctrlPoints[1][3*i+1] + 3.0*ctrlPoints[1][3*i+2] + ctrlPoints[1][3*i+3]) +\
        dz*(ctrlPoints[2][3*i] + 3.0*ctrlPoints[2][3*i+1] + 3.0*ctrlPoints[2][3*i+2] + ctrlPoints[2][3*i+3]));

    sum = (3.0/8.0) * sum;
//    std::cout << "sum: " << sum << std::endl;
//    if(sum==0.0)
//    {
//        for(unsigned short i=0 ; i<nbCtrlPoints ; i++)
//            std::cout << "F = (" << ctrlPoints[0][i] << " " << ctrlPoints[1][i] << " " << ctrlPoints[2][i] << std::endl;
//        for(unsigned short i=0 ; i<m ; i++)
/*            std::cout << "small element: " << "dx, dy, dz:" << dx << "," << dy << "," << dz << " " << \
//            dx*(ctrlPoints[0][3*i] + 3.0*ctrlPoints[0][3*i+1] + 3.0*ctrlPoints[0][3*i+2] + ctrlPoints[0][3*i+3]) << " " <<\
//            dy*(ctrlPoints[1][3*i] + 3.0*ctrlPoints[1][3*i+1] + 3.0*ctrlPoints[1][3*i+2] + ctrlPoints[1][3*i+3]) << " " <<\
//            dz*(ctrlPoints[2][3*i] + 3.0*ctrlPoints[2][3*i+1] + 3.0*ctrlPoints[2][3*i+2] + ctrlPoints[2][3*i+3]) << std::endl;*/
//    }

    //so that it bind the last values to the integral
    if(m_rest==1)
    {
        //trapezoidal case
        sum += ABS( (dx/2.0) * (ctrlPoints[0][nbCtrlPoints-2] + ctrlPoints[0][nbCtrlPoints-1]) + (dy/2.0) * (ctrlPoints[1][nbCtrlPoints-2] + ctrlPoints[1][nbCtrlPoints-1]) + (dz/2.0) * (ctrlPoints[2][nbCtrlPoints-2] + ctrlPoints[2][nbCtrlPoints-1]) );
    }
    else if(m_rest==2) //means s==2
    {
        //trapezoidal case 2
        sum += ABS( (dx/2.0) * (ctrlPoints[0][nbCtrlPoints-3] + 2*ctrlPoints[0][nbCtrlPoints-2] + ctrlPoints[0][nbCtrlPoints-1]) + (dy/2.0) * (ctrlPoints[1][nbCtrlPoints-3] + 2*ctrlPoints[1][nbCtrlPoints-2] + ctrlPoints[1][nbCtrlPoints-1]) + (dz/2.0) * (ctrlPoints[2][nbCtrlPoints-3] + 2*ctrlPoints[2][nbCtrlPoints-2] + ctrlPoints[2][nbCtrlPoints-1]));
    }

    return sum;
}


//error computation
double C_toolbox_integrate::maxF2(double* ctrlPoints, unsigned short nbCtrlPoints, double subIntervalLength)
{
    double err = 0;
    for(unsigned short i=1 ; i<nbCtrlPoints-1 ; i++)
    {
        err = MAX(err, ABS(ctrlPoints[i-1] -2*ctrlPoints[i] + ctrlPoints[i]));
        //abs((float) 1.0);
    }
    return err/(subIntervalLength*subIntervalLength);
}
double C_toolbox_integrate::maxF4(double* ctrlPoints, unsigned short nbCtrlPoints, double subIntervalLength)
{
    double err = 0;
    for(unsigned short i=2 ; i<nbCtrlPoints-2 ; i++)
    {
        err = MAX(err, ABS(ctrlPoints[i-2] - 4*ctrlPoints[i-1] + 6*ctrlPoints[i] - 4*ctrlPoints[i+1] + ctrlPoints[i+2]));
    }
    return err/(subIntervalLength*subIntervalLength*subIntervalLength*subIntervalLength);
}


double C_toolbox_integrate::errMaxMiddleSum(double* ctrlPoints, unsigned short nbCtrlPoints, double intervalLength, double subIntervalLength)
{
    //(b-a)^3 max(f'')/(24 n^2)
    return (intervalLength*intervalLength*intervalLength)*maxF2(ctrlPoints, nbCtrlPoints, subIntervalLength)/(24.0 * ((double) nbCtrlPoints) * ((double) nbCtrlPoints));
}

double C_toolbox_integrate::errTrapezoidalRule(double* ctrlPoints, unsigned short nbCtrlPoints, double intervalLength, double subIntervalLength)
{
    //(b-a)^3 max(f'')/(12 n^2)
    return (intervalLength*intervalLength*intervalLength)*maxF2(ctrlPoints, nbCtrlPoints, subIntervalLength)/(12.0 * ((double) nbCtrlPoints) * ((double) nbCtrlPoints));
}

double C_toolbox_integrate::errSimpson38Rule(double* ctrlPoints, unsigned short nbCtrlPoints, double intervalLength, double subIntervalLength)
{
    //(b-a)h^4 max(f'''')/80
    return intervalLength*subIntervalLength*subIntervalLength*subIntervalLength*subIntervalLength*maxF4(ctrlPoints, nbCtrlPoints, subIntervalLength)/80.0;
}
