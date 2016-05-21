#include <C_measure_synthetic.h>

C_measure_synthetic::C_measure_synthetic()
{
    //ctor
}

C_measure_synthetic::~C_measure_synthetic()
{
    //dtor
}


double C_measure_synthetic::meanFiberError(CtlCtlStruct* FIBERS)
{
    if(FIBERS->N==0)
    {
        return 0.0;
    }
    double eps = 0.0;
    for(unsigned long i=0 ; i<FIBERS->N ; i++ )
    {
        eps += fiberError(&(FIBERS->fibers[i]));
    }
    return eps/((double) FIBERS->N);
}
double C_measure_synthetic::meanFiberLength(CtlCtlStruct* FIBERS)
{
    if(FIBERS->N==0)
    {
        return 0.0;
    }
    double L = 0.0;
    for(unsigned long i=0 ; i<FIBERS->N ; i++ )
    {
        L += fiberLength(&(FIBERS->fibers[i]));
    }
    return L/((double) FIBERS->N);
}
double C_measure_synthetic::stdFiberLength(CtlCtlStruct* FIBERS, double muL)
{
    if(FIBERS->N==0 || FIBERS->N==1)
    {
        return 0.0;
    }

    double L = muL;///meanFiberLength(FIBERS);
    double S = 0, Si;
    for(unsigned long i=0 ; i<FIBERS->N ; i++ )
    {
        Si = fiberLength(&(FIBERS->fibers[i]));
        Si = (Si-L)*(Si-L);
        S += Si;
    }
    return sqrt(S/((double) FIBERS->N -1));
}


double C_measure_synthetic::meanFiberCurvatureCompare(CtlCtlStruct* FIBERS)
{
    if(FIBERS->N==0)
    {
        return 0.0;
    }

    double S = 0, Si;
    for(unsigned long i=0 ; i<FIBERS->N ; i++ )
    {
        Si = fiberCurvatureCompare(&(FIBERS->fibers[i]));
        S += Si;
    }
    return S/((double) FIBERS->N);
}

double C_measure_synthetic::meanFiberCurvature(CtlCtlStruct* FIBERS)
{
    if(FIBERS->N==0)
    {
        return 0.0;
    }

    double S = 0, Si;
    for(unsigned long i=0 ; i<FIBERS->N ; i++ )
    {
        Si = fiberCurvature(&(FIBERS->fibers[i]));
        S += Si;
    }
    return S/((double) FIBERS->N);
}
double C_measure_synthetic::stdMeanFiberCurvature(CtlCtlStruct* FIBERS, double muK)
{
    if(FIBERS->N==0 || FIBERS->N==1)
    {
        return 0.0;
    }

    double S = 0, Si;
    for(unsigned long i=0 ; i<FIBERS->N ; i++ )
    {
        Si = fiberCurvature(&(FIBERS->fibers[i]));
        S += (Si-muK)*(Si-muK);
    }
    return sqrt(S/((double) (FIBERS->N-1)));
}
double C_measure_synthetic::meanSTDFiberCurvature(CtlCtlStruct* FIBERS)
{
    if(FIBERS->N==0)
    {
        return 0.0;
    }

    double S = 0, Si;
    for(unsigned long i=0 ; i<FIBERS->N ; i++ )
    {
        Si = fiberCurvature(&(FIBERS->fibers[i]));
        Si = stdFiberCurvature(&(FIBERS->fibers[i]), Si);
        S += Si;
    }
    return S/((double) FIBERS->N);
}







double C_measure_synthetic::meanFiberLengthS(CtlCtlStruct* DPHI)
{
    double L = 0.0;
    for(unsigned long i=0 ; i<DPHI->N ; i++ )
    {
        L += fiberLengthS(&(DPHI->fibers[i]));
    }
    return L/((double) DPHI->N);
}
double C_measure_synthetic::stdFiberLengthS(CtlCtlStruct* DPHI, double muL)
{
    if(DPHI->N==0 || DPHI->N==1)
    {
        return 0.0;
    }

    double S = 0, Si;
    for(unsigned long i=0 ; i<DPHI->N ; i++ )
    {
        Si = fiberLengthS(&(DPHI->fibers[i]));
        Si = (Si-muL)*(Si-muL);
        S += Si;
    }
    return sqrt(S/((double) DPHI->N -1));
}
double C_measure_synthetic::meanFiberCurvatureS(CtlCtlStruct* FIBERS, CtlCtlStruct* DPHI, CtlCtlStruct* DDPHI)
{
    if(FIBERS->N==0)
    {
        return 0.0;
    }

    double S = 0, Si;
    for(unsigned long i=0 ; i<FIBERS->N ; i++ )
    {
        Si = fiberCurvatureS(&(FIBERS->fibers[i]), &(DPHI->fibers[i]), &(DDPHI->fibers[i]));
        S += Si;
    }
    return S/((double) FIBERS->N);
}
double C_measure_synthetic::meanFiberCurvatureCompareS(CtlCtlStruct* FIBERS, CtlCtlStruct* DPHI, CtlCtlStruct* DDPHI)
{
    if(FIBERS->N==0)
    {
        return 0.0;
    }

    double S = 0, Si;
    for(unsigned long i=0 ; i<FIBERS->N ; i++ )
    {
        Si = fiberCurvatureCompareS(&(FIBERS->fibers[i]), &(DPHI->fibers[i]), &(DDPHI->fibers[i]));
        S += Si;
    }
    return S/((double) FIBERS->N);
}

double C_measure_synthetic::stdMeanFiberCurvatureS(CtlCtlStruct* FIBERS, CtlCtlStruct* DPHI, CtlCtlStruct* DDPHI, double muK)
{
    if(FIBERS->N==0 || FIBERS->N==1)
    {
        return 0.0;
    }

    double S = 0, Si;
    for(unsigned long i=0 ; i<FIBERS->N ; i++ )
    {
        Si = fiberCurvatureS(&(FIBERS->fibers[i]), &(DPHI->fibers[i]), &(DDPHI->fibers[i]));
        S += (Si-muK)*(Si-muK);
    }
    return sqrt(S/((double) (FIBERS->N-1)));
}
double C_measure_synthetic::stdMeanFiberCurvatureCompareS(CtlCtlStruct* FIBERS, CtlCtlStruct* DPHI, CtlCtlStruct* DDPHI, double muK)
{
    if(FIBERS->N==0 || FIBERS->N==1)
    {
        return 0.0;
    }

    double S = 0, Si;
    for(unsigned long i=0 ; i<FIBERS->N ; i++ )
    {
        Si = fiberCurvatureCompareS(&(FIBERS->fibers[i]), &(DPHI->fibers[i]), &(DDPHI->fibers[i]));
        S += (Si-muK)*(Si-muK);
    }
    return sqrt(S/((double) (FIBERS->N-1)));
}




void C_measure_synthetic::syntheticTensor(double x, double y, double z, double* T)
{
    double r = sqrt(SQR(x) + SQR(y));
    double v1x, v1y, v1z;
    double v2x, v2y, v2z;
    double v3x, v3y, v3z;

    v1x = cos(helix_angle)*y/r;
    v1y = -cos(helix_angle)*x/r;
    v1z = sin(helix_angle);

    if( ABS(helix_angle - 0.5*Pi) <0.000000001 )
    {
        v2x = 1.0;
        v2y = 0.0;
        v2z = 0.0;
    }
    else
    {
        v2x = x/r;
        v2y = y/r;
        v2z = 0.0;
    }

    v3x = v1y*v2z - v1z*v2y;
    v3y = v1z*v2x - v1x*v2z;
    v3z = v1x*v2y - v1y*v2x;

    T[0] = lambda1*SQR(v1x) + lambda2*SQR(v2x) + lambda3*SQR(v3x);
    T[1] = lambda1*v1x*v1y + lambda2*v2x*v2y + lambda3*v3x*v3y;
    T[2] = lambda1*v1x*v1z + lambda2*v2x*v2z + lambda3*v3x*v3z;
    T[3] = lambda1*SQR(v1y) + lambda2*SQR(v2y) + lambda3*SQR(v3y);
    T[4] = lambda1*v1y*v1z + lambda2*v2y*v2z + lambda3*v3y*v3z;
    T[5] = lambda1*SQR(v1z) + lambda2*SQR(v2z) + lambda3*SQR(v3z);

    return;
}

double C_measure_synthetic::syntheticCurvature(double x, double y, double z)
{
    return ABS(SQR(cos(helix_angle))/sqrt(SQR(x) + SQR(y)));
}

double C_measure_synthetic::fiberLength(CtlStruct* phi)
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

double C_measure_synthetic::fiberError(CtlStruct* phi)
{
    /**this function return the mean fiber error, second definition*/
    double L = fiberLength(phi);
    ///compute DPHI as it's done in fiberLength
    if(phi->N_element<=1)
    {
        return 0.0;
    }

    if(phi->N_element==2)
    {
        ///be careful with this case
        double ds = 1.0;
        double dx = (phi->elts[3]-phi->elts[0])/ds;
        double dy = (phi->elts[4]-phi->elts[1])/ds;
        double dz = (phi->elts[5]-phi->elts[2])/ds;
        double dphi = sqrt(SQR(dx) + SQR(dy) + SQR(dz));
        double* T = new double[6];
        double Tx, Ty, Tz;
        syntheticTensor(0.5*(phi->elts[3]+phi->elts[0]), 0.5*(phi->elts[4]+phi->elts[1]), 0.5*(phi->elts[5]+phi->elts[2]), T);
        Tx = T[0]*dx + T[1]*dy + T[2]*dz;
        Ty = T[1]*dx + T[3]*dy + T[4]*dz;
        Tz = T[2]*dx + T[4]*dy + T[5]*dz;
        delete T;
        return (1.0/dphi)*(1.0-(sqrt(SQR(Tx) + SQR(Ty) + SQR(Tz) )/(lambda1*dphi)));
    }

    double *DPHIX = new double[phi->N_element];
    if(DPHIX==NULL)
    {
        return -1.0;
    }
    double *DPHIY = new double[phi->N_element];
    if(DPHIY==NULL)
    {
        delete DPHIX;
        return -1.0;
    }
    double *DPHIZ = new double[phi->N_element];
    if(DPHIZ==NULL)
    {
        delete DPHIY;
        delete DPHIX;
        return -1.0;
    }
    double ds = 2.0/((double) (phi->N_element-1));
    DPHIX[0] = (phi->elts[3*1]-phi->elts[3*0])/(ds/2.0);
    DPHIY[0] = (phi->elts[3*1+1]-phi->elts[3*0+1])/(ds/2.0);
    DPHIZ[0] = (phi->elts[3*1+2]-phi->elts[3*0+2])/(ds/2.0);
    for(unsigned long i=1 ; i<phi->N_element-1 ; i++)
    {
        DPHIX[i] = (phi->elts[3*(i+1)]-phi->elts[3*(i-1)])/ds;
        DPHIY[i] = (phi->elts[3*(i+1)+1]-phi->elts[3*(i-1)+1])/ds;
        DPHIZ[i] = (phi->elts[3*(i+1)+2]-phi->elts[3*(i-1)+2])/ds;
    }
    DPHIX[phi->N_element-1] = (phi->elts[3*(phi->N_element-1)]-phi->elts[3*(phi->N_element-2)])/(ds/2.0);
    DPHIY[phi->N_element-1] = (phi->elts[3*(phi->N_element-1)+1]-phi->elts[3*(phi->N_element-2)+1])/(ds/2.0);
    DPHIZ[phi->N_element-1] = (phi->elts[3*(phi->N_element-1)+2]-phi->elts[3*(phi->N_element-2)+2])/(ds/2.0);

    ///compute pointwise error
    double *EPS = new double[phi->N_element];
    if(EPS==NULL)
    {
        delete DPHIX;
        delete DPHIY;
        delete DPHIZ;
    }

    double dphi;
    double* T = new double[6];
    double Tx, Ty, Tz;
    for(unsigned long i=0 ; i<phi->N_element ; i++)
    {
        dphi = sqrt(SQR(DPHIX[i]) + SQR(DPHIY[i]) + SQR(DPHIZ[i]));
        syntheticTensor(phi->elts[3*i], phi->elts[3*i+1], phi->elts[3*i+2], T);
        Tx = T[0]*DPHIX[i] + T[1]*DPHIY[i] + T[2]*DPHIZ[i];
        Ty = T[1]*DPHIX[i] + T[3]*DPHIY[i] + T[4]*DPHIZ[i];
        Tz = T[2]*DPHIX[i] + T[4]*DPHIY[i] + T[5]*DPHIZ[i];
        EPS[i] = (1.0-(sqrt(SQR(Tx) + SQR(Ty) + SQR(Tz) )/(lambda1*dphi)))*dphi;
    }

    ///integrate
    C_toolbox_integrate* t = new C_toolbox_integrate();

    double eps = t->run(EPS, phi->N_element, 1.0, SIMPSON_RULE, NULL);

    ///destroy!!!!!!
    delete t;//->~C_integrate();
    delete DPHIX;
    delete DPHIY;
    delete DPHIZ;
    delete EPS;
    delete T;

    return eps/SQR(L);
}



double C_measure_synthetic::fiberCurvature(CtlStruct* phi)
{
    ///compute first derivative
    ///compute DPHI as it's done in fiberLength
    if(phi->N_element<=2)
    {
        return 0.0;
    }

    double *DPHIX = new double[phi->N_element];
    if(DPHIX==NULL)
    {
        return -1.0;
    }
    double *DPHIY = new double[phi->N_element];
    if(DPHIY==NULL)
    {
        delete DPHIX;
        return -1.0;
    }
    double *DPHIZ = new double[phi->N_element];
    if(DPHIZ==NULL)
    {
        delete DPHIY;
        delete DPHIX;
        return -1.0;
    }
    double ds = 2.0/((double) (phi->N_element-1));
    DPHIX[0] = (phi->elts[3*1]-phi->elts[3*0])/(ds/2.0);
    DPHIY[0] = (phi->elts[3*1+1]-phi->elts[3*0+1])/(ds/2.0);
    DPHIZ[0] = (phi->elts[3*1+2]-phi->elts[3*0+2])/(ds/2.0);
    for(unsigned long i=1 ; i<phi->N_element-1 ; i++)
    {
        DPHIX[i] = (phi->elts[3*(i+1)]-phi->elts[3*(i-1)])/ds;
        DPHIY[i] = (phi->elts[3*(i+1)+1]-phi->elts[3*(i-1)+1])/ds;
        DPHIZ[i] = (phi->elts[3*(i+1)+2]-phi->elts[3*(i-1)+2])/ds;
    }
    DPHIX[phi->N_element-1] = (phi->elts[3*(phi->N_element-1)]-phi->elts[3*(phi->N_element-2)])/(ds/2.0);
    DPHIY[phi->N_element-1] = (phi->elts[3*(phi->N_element-1)+1]-phi->elts[3*(phi->N_element-2)+1])/(ds/2.0);
    DPHIZ[phi->N_element-1] = (phi->elts[3*(phi->N_element-1)+2]-phi->elts[3*(phi->N_element-2)+2])/(ds/2.0);


    ///compute second derivative
    double *DDPHIX = new double[phi->N_element];
    if(DDPHIX==NULL)
    {
        delete DPHIX;
        delete DPHIY;
        delete DPHIZ;
        return -1.0;
    }
    double *DDPHIY = new double[phi->N_element];
    if(DDPHIY==NULL)
    {
        delete DDPHIX;
        delete DPHIX;
        delete DPHIY;
        delete DPHIZ;
        return -1.0;
    }
    double *DDPHIZ = new double[phi->N_element];
    if(DDPHIZ==NULL)
    {
        delete DDPHIY;
        delete DDPHIX;
        delete DPHIY;
        delete DPHIX;
        delete DPHIZ;
        return -1.0;
    }
    ds = 1.0/((double) (phi->N_element-1));
    DDPHIX[0] = (DPHIX[1]-DPHIX[0])/ds;
    DDPHIY[0] = (DPHIY[1]-DPHIY[0])/ds;
    DDPHIZ[0] = (DPHIZ[1]-DPHIZ[0])/ds;
    for(unsigned long i=1 ; i<phi->N_element-1 ; i++)
    {
        DDPHIX[i] = (phi->elts[3*(i+1)]-2.0*phi->elts[3*i]+phi->elts[3*(i-1)])/SQR(ds);
        DDPHIY[i] = (phi->elts[3*(i+1)+1]-2.0*phi->elts[3*i+1]+phi->elts[3*(i-1)+1])/SQR(ds);
        DDPHIZ[i] = (phi->elts[3*(i+1)+2]-2.0*phi->elts[3*i+2]+phi->elts[3*(i-1)+2])/SQR(ds);
    }
    DDPHIX[phi->N_element-1] = (DPHIX[phi->N_element-1]-DPHIX[phi->N_element-2])/ds;
    DDPHIY[phi->N_element-1] = (DPHIY[phi->N_element-1]-DPHIY[phi->N_element-2])/ds;
    DDPHIZ[phi->N_element-1] = (DPHIZ[phi->N_element-1]-DPHIZ[phi->N_element-2])/ds;

    ///compute curvature
    for(unsigned long i=1 ; i<phi->N_element-1 ; i++)
    {
        DDPHIX[i] = sqrt(SQR( DPHIY[i]*DDPHIZ[i] - DPHIZ[i]*DDPHIY[i] ) +\
                         SQR( DPHIZ[i]*DDPHIX[i] - DPHIX[i]*DDPHIZ[i] ) +\
                         SQR( DPHIX[i]*DDPHIY[i] - DPHIY[i]*DDPHIX[i] )) / (sqrt(SQR(DPHIX[i]) + SQR(DPHIY[i]) + SQR(DPHIZ[i]))* (SQR(DPHIX[i]) + SQR(DPHIY[i]) + SQR(DPHIZ[i])));
    }
    DDPHIX[0] = DDPHIX[1];
    DDPHIX[phi->N_element-1] = DDPHIX[phi->N_element-2];

    delete DDPHIY;
    delete DDPHIZ;
    delete DPHIY;
    delete DPHIX;
    delete DPHIZ;

    ///integrate curvature
    C_toolbox_integrate* t = new C_toolbox_integrate();
    //double K = t->run(DDPHIX, phi->N_element, 1.0, SIMPSON_RULE, NULL);
    double K = t->run(DDPHIX, phi->N_element, 1.0, SIMPSON_RULE, NULL);
    delete t;//->~C_integrate();
    ///return value
    return K;
}

double C_measure_synthetic::stdFiberCurvature(CtlStruct* phi, double muK)
{
    ///compute first derivative
    ///compute DPHI as it's done in fiberLength
    if(phi->N_element<=2)
    {
        return 0.0;
    }

    double *DPHIX = new double[phi->N_element];
    if(DPHIX==NULL)
    {
        return -1.0;
    }
    double *DPHIY = new double[phi->N_element];
    if(DPHIY==NULL)
    {
        delete DPHIX;
        return -1.0;
    }
    double *DPHIZ = new double[phi->N_element];
    if(DPHIZ==NULL)
    {
        delete DPHIY;
        delete DPHIX;
        return -1.0;
    }
    double ds = 2.0/((double) (phi->N_element-1));
    DPHIX[0] = (phi->elts[3*1]-phi->elts[3*0])/(ds/2.0);
    DPHIY[0] = (phi->elts[3*1+1]-phi->elts[3*0+1])/(ds/2.0);
    DPHIZ[0] = (phi->elts[3*1+2]-phi->elts[3*0+2])/(ds/2.0);
    for(unsigned long i=1 ; i<phi->N_element-1 ; i++)
    {
        DPHIX[i] = (phi->elts[3*(i+1)]-phi->elts[3*(i-1)])/ds;
        DPHIY[i] = (phi->elts[3*(i+1)+1]-phi->elts[3*(i-1)+1])/ds;
        DPHIZ[i] = (phi->elts[3*(i+1)+2]-phi->elts[3*(i-1)+2])/ds;
    }
    DPHIX[phi->N_element-1] = (phi->elts[3*(phi->N_element-1)]-phi->elts[3*(phi->N_element-2)])/(ds/2.0);
    DPHIY[phi->N_element-1] = (phi->elts[3*(phi->N_element-1)+1]-phi->elts[3*(phi->N_element-2)+1])/(ds/2.0);
    DPHIZ[phi->N_element-1] = (phi->elts[3*(phi->N_element-1)+2]-phi->elts[3*(phi->N_element-2)+2])/(ds/2.0);


    ///compute second derivative
    double *DDPHIX = new double[phi->N_element];
    if(DDPHIX==NULL)
    {
        delete DPHIX;
        delete DPHIY;
        delete DPHIZ;
        return -1.0;
    }
    double *DDPHIY = new double[phi->N_element];
    if(DDPHIY==NULL)
    {
        delete DDPHIX;
        delete DPHIX;
        delete DPHIY;
        delete DPHIZ;
        return -1.0;
    }
    double *DDPHIZ = new double[phi->N_element];
    if(DDPHIZ==NULL)
    {
        delete DDPHIY;
        delete DDPHIX;
        delete DPHIY;
        delete DPHIX;
        delete DPHIZ;
        return -1.0;
    }
    ds = 1.0/((double) (phi->N_element-1));
    DDPHIX[0] = (DPHIX[1]-DPHIX[0])/ds;
    DDPHIY[0] = (DPHIY[1]-DPHIY[0])/ds;
    DDPHIZ[0] = (DPHIZ[1]-DPHIZ[0])/ds;
    for(unsigned long i=1 ; i<phi->N_element-1 ; i++)
    {
        DDPHIX[i] = (phi->elts[3*(i+1)]-2.0*phi->elts[3*i]+phi->elts[3*(i-1)])/SQR(ds);
        DDPHIY[i] = (phi->elts[3*(i+1)+1]-2.0*phi->elts[3*i+1]+phi->elts[3*(i-1)+1])/SQR(ds);
        DDPHIZ[i] = (phi->elts[3*(i+1)+2]-2.0*phi->elts[3*i+2]+phi->elts[3*(i-1)+2])/SQR(ds);
    }
    DDPHIX[phi->N_element-1] = (DPHIX[phi->N_element-1]-DPHIX[phi->N_element-2])/ds;
    DDPHIY[phi->N_element-1] = (DPHIY[phi->N_element-1]-DPHIY[phi->N_element-2])/ds;
    DDPHIZ[phi->N_element-1] = (DPHIZ[phi->N_element-1]-DPHIZ[phi->N_element-2])/ds;

    ///compute curvature
    for(unsigned long i=1 ; i<phi->N_element-1 ; i++)
    {
        DDPHIX[i] = sqrt(SQR( DPHIY[i]*DDPHIZ[i] - DPHIZ[i]*DDPHIY[i] ) +\
                         SQR( DPHIZ[i]*DDPHIX[i] - DPHIX[i]*DDPHIZ[i] ) +\
                         SQR( DPHIX[i]*DDPHIY[i] - DPHIY[i]*DDPHIX[i] )) / sqrt(SQR(DPHIX[i]) + SQR(DPHIY[i]) + SQR(DPHIZ[i]));
        DDPHIX[i] = SQR(DDPHIX[i] - muK);
    }

    delete DDPHIY;
    delete DDPHIZ;
    delete DPHIY;
    delete DPHIX;
    delete DPHIZ;

    ///integrate curvature
    C_toolbox_integrate* t = new C_toolbox_integrate();
    double K = t->run(DDPHIX, phi->N_element, 1.0, SIMPSON_RULE, NULL);
    delete t;//->~C_integrate();
    ///return value
    return sqrt(K);
}

double C_measure_synthetic::fiberCurvatureCompare(CtlStruct* phi)
{
    if(phi->N_element<=2)
    {
        return 0.0;
    }

    double *DPHIX = new double[phi->N_element];
    if(DPHIX==NULL)
    {
        return -1.0;
    }
    double *DPHIY = new double[phi->N_element];
    if(DPHIY==NULL)
    {
        delete DPHIX;
        return -1.0;
    }
    double *DPHIZ = new double[phi->N_element];
    if(DPHIZ==NULL)
    {
        delete DPHIY;
        delete DPHIX;
        return -1.0;
    }
    double ds = 2.0/((double) (phi->N_element-1));
    DPHIX[0] = (phi->elts[3*1]-phi->elts[3*0])/(ds/2.0);
    DPHIY[0] = (phi->elts[3*1+1]-phi->elts[3*0+1])/(ds/2.0);
    DPHIZ[0] = (phi->elts[3*1+2]-phi->elts[3*0+2])/(ds/2.0);
    for(unsigned long i=1 ; i<phi->N_element-1 ; i++)
    {
        DPHIX[i] = (phi->elts[3*(i+1)]-phi->elts[3*(i-1)])/ds;
        DPHIY[i] = (phi->elts[3*(i+1)+1]-phi->elts[3*(i-1)+1])/ds;
        DPHIZ[i] = (phi->elts[3*(i+1)+2]-phi->elts[3*(i-1)+2])/ds;
    }
    DPHIX[phi->N_element-1] = (phi->elts[3*(phi->N_element-1)]-phi->elts[3*(phi->N_element-2)])/(ds/2.0);
    DPHIY[phi->N_element-1] = (phi->elts[3*(phi->N_element-1)+1]-phi->elts[3*(phi->N_element-2)+1])/(ds/2.0);
    DPHIZ[phi->N_element-1] = (phi->elts[3*(phi->N_element-1)+2]-phi->elts[3*(phi->N_element-2)+2])/(ds/2.0);


    ///compute second derivative
    double *DDPHIX = new double[phi->N_element];
    if(DDPHIX==NULL)
    {
        delete DPHIX;
        delete DPHIY;
        delete DPHIZ;
        return -1.0;
    }
    double *DDPHIY = new double[phi->N_element];
    if(DDPHIY==NULL)
    {
        delete DDPHIX;
        delete DPHIX;
        delete DPHIY;
        delete DPHIZ;
        return -1.0;
    }
    double *DDPHIZ = new double[phi->N_element];
    if(DDPHIZ==NULL)
    {
        delete DDPHIY;
        delete DDPHIX;
        delete DPHIY;
        delete DPHIX;
        delete DPHIZ;
        return -1.0;
    }
    ds = 1.0/((double) (phi->N_element-1));
    DDPHIX[0] = (DPHIX[1]-DPHIX[0])/ds;
    DDPHIY[0] = (DPHIY[1]-DPHIY[0])/ds;
    DDPHIZ[0] = (DPHIZ[1]-DPHIZ[0])/ds;
    for(unsigned long i=1 ; i<phi->N_element-1 ; i++)
    {
        DDPHIX[i] = (phi->elts[3*(i+1)]-2.0*phi->elts[3*i]+phi->elts[3*(i-1)])/SQR(ds);
        DDPHIY[i] = (phi->elts[3*(i+1)+1]-2.0*phi->elts[3*i+1]+phi->elts[3*(i-1)+1])/SQR(ds);
        DDPHIZ[i] = (phi->elts[3*(i+1)+2]-2.0*phi->elts[3*i+2]+phi->elts[3*(i-1)+2])/SQR(ds);
    }
    DDPHIX[phi->N_element-1] = (DPHIX[phi->N_element-1]-DPHIX[phi->N_element-2])/ds;
    DDPHIY[phi->N_element-1] = (DPHIY[phi->N_element-1]-DPHIY[phi->N_element-2])/ds;
    DDPHIZ[phi->N_element-1] = (DPHIZ[phi->N_element-1]-DPHIZ[phi->N_element-2])/ds;

    ///compute curvature
    for(unsigned long i=1 ; i<phi->N_element-1 ; i++)
    {
        DDPHIX[i] = sqrt(SQR( DPHIY[i]*DDPHIZ[i] - DPHIZ[i]*DDPHIY[i] ) +\
                         SQR( DPHIZ[i]*DDPHIX[i] - DPHIX[i]*DDPHIZ[i] ) +\
                         SQR( DPHIX[i]*DDPHIY[i] - DPHIY[i]*DDPHIX[i] )) / (sqrt(SQR(DPHIX[i]) + SQR(DPHIY[i]) + SQR(DPHIZ[i]))* (SQR(DPHIX[i]) + SQR(DPHIY[i]) + SQR(DPHIZ[i])));
        DDPHIX[i] = ABS(DDPHIX[i] - syntheticCurvature(phi->elts[3*i], phi->elts[3*i+1], phi->elts[3*i+2]));
    }

    DDPHIX[0] = DDPHIX[1];
    DDPHIX[0] = ABS(DDPHIX[0] - syntheticCurvature(phi->elts[0], phi->elts[1], phi->elts[2]));
    DDPHIX[phi->N_element-1] = DDPHIX[phi->N_element-2];
    DDPHIX[phi->N_element-1] = ABS(DDPHIX[phi->N_element-1] - syntheticCurvature(phi->elts[3*(phi->N_element-1)], phi->elts[3*(phi->N_element-1)+1], phi->elts[3*(phi->N_element-1)+2]));

    delete DDPHIY;
    delete DDPHIZ;
    delete DPHIY;
    delete DPHIX;
    delete DPHIZ;

    ///integrate curvature
    C_toolbox_integrate* t = new C_toolbox_integrate();
    //double K = t->run(DDPHIX, phi->N_element, 1.0, SIMPSON_RULE, NULL);
    double K = t->run(DDPHIX, phi->N_element, 1.0, SIMPSON_RULE, NULL);
    delete t;//->~C_integrate();
    ///return value
    return K;
}












double C_measure_synthetic::fiberLengthS(CtlStruct* dphi)
{
    double* DPHI = new double[dphi->N_element];
    if(DPHI==NULL)
    {
        return -1.0;
    }

    ///calculate speed norm
    for(unsigned long i=0 ; i<dphi->N_element ; i++)
    {
        DPHI[i] = sqrt(SQR(dphi->elts[3*i]) + SQR(dphi->elts[3*i+1]) + SQR(dphi->elts[3*i+2]));
    }

    ///integrate curvature
    C_toolbox_integrate* t = new C_toolbox_integrate();
    double L = t->run(DPHI, dphi->N_element, 1.0, SIMPSON_RULE, NULL);
    delete t;//->~C_integrate();

    delete DPHI;
    return L;
}
double C_measure_synthetic::fiberCurvatureS(CtlStruct* phi, CtlStruct* dphi, CtlStruct* ddphi)
{
    double* DPHI = new double[phi->N_element];
    if(DPHI==NULL)
    {
        return -1.0;
    }

    ///calculate speed norm
    for(unsigned long i=0 ; i<dphi->N_element ; i++)
    {
        DPHI[i] = sqrt(SQR(dphi->elts[3*i+1]*(ddphi->elts[3*i+2]) - dphi->elts[3*i+2]*(ddphi->elts[3*i+1]) ) +\
             SQR(dphi->elts[3*i+2]*(ddphi->elts[3*i]) - dphi->elts[3*i]*(ddphi->elts[3*i+2])) +\
             SQR(dphi->elts[3*i]*(ddphi->elts[3*i+1]) - dphi->elts[3*i+1]*(ddphi->elts[3*i])))/((SQR(dphi->elts[3*i]) + SQR(dphi->elts[3*i+1]) + SQR(dphi->elts[3*i+2]))*sqrt(SQR(dphi->elts[3*i]) + SQR(dphi->elts[3*i+1]) + SQR(dphi->elts[3*i+2])));
             //cout << DPHI[i] << endl;
    }
    //getchar();

    ///integrate curvature
    C_toolbox_integrate* t = new C_toolbox_integrate();
    double K = t->run(DPHI, dphi->N_element, 1.0, SIMPSON_RULE, NULL);
    delete t;//->~C_integrate();
    delete DPHI;
    return K;
}
double C_measure_synthetic::fiberCurvatureCompareS(CtlStruct* phi, CtlStruct* dphi, CtlStruct* ddphi)
{
    double* DPHI = new double[dphi->N_element];
    if(DPHI==NULL)
    {
        return -1.0;
    }

    ///calculate speed norm
    for(unsigned long i=0 ; i<dphi->N_element ; i++)
    {
        DPHI[i] = sqrt(SQR(dphi->elts[3*i+1]*(ddphi->elts[3*i+2]) - dphi->elts[3*i+2]*(ddphi->elts[3*i+1]) ) +\
             SQR(dphi->elts[3*i+2]*(ddphi->elts[3*i]) - dphi->elts[3*i]*(ddphi->elts[3*i+2])) +\
             SQR(dphi->elts[3*i]*(ddphi->elts[3*i+1]) - dphi->elts[3*i+1]*(ddphi->elts[3*i])))/\
             ((SQR(dphi->elts[3*i]) + SQR(dphi->elts[3*i+1]) + SQR(dphi->elts[3*i+2]))*sqrt(SQR(dphi->elts[3*i]) + SQR(dphi->elts[3*i+1]) + SQR(dphi->elts[3*i+2])));
        DPHI[i] = ABS(syntheticCurvature(phi->elts[3*i], phi->elts[3*i+1], phi->elts[3*i+2])-DPHI[i]);
    }

    ///integrate curvature
    C_toolbox_integrate* t = new C_toolbox_integrate();
    double K = t->run(DPHI, dphi->N_element, 1.0, SIMPSON_RULE, NULL);
    delete t;//->~C_integrate();
    delete DPHI;
    return K;
}



