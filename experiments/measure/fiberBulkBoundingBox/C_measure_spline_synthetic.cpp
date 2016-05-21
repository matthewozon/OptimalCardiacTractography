#include <C_measure_spline_synthetic.h>

C_measure_spline_synthetic::C_measure_spline_synthetic()
{
    //ctor
}

C_measure_spline_synthetic::~C_measure_spline_synthetic()
{
    //dtor
}

void C_measure_spline_synthetic::getLenAndErrPop(CtlCtlStruct* FIBERS, bool red)
{
    CtlCtlStruct* splineFIBERS;
    if(red)
    {
        unsigned long nb_sub_section = 5;
        C_spline* S = new C_spline();
        if(!S->computeFiberSpline(FIBERS, &splineFIBERS, nb_sub_section, CURVE))
        {
            delete S;//->~C_spline();
            return; ///?
        }
        delete S;//->~C_spline();
    }
    else
    {
        splineFIBERS = FIBERS;
    }

    FILE* f = fopen("fibPopulationErrAndLen.txt","w");
    if(f==NULL)
    {
        if(red)
        {
            delete splineFIBERS;//->~CtlCtlStruct();
        }
        return;
    }

    C_parametric_comparison* tool = new C_parametric_comparison();
    tool->alpha = helix_angle;
    for(unsigned long i=0 ; i<splineFIBERS->N ; i++)
    {
        fprintf(f,"%f, %f, %f, %f, %f\n", fiberLength(&(splineFIBERS->fibers[i])), fiberError(&(splineFIBERS->fibers[i])), tool->getDist(&(splineFIBERS->fibers[i])), tool->getDist3(&(splineFIBERS->fibers[i])), tool->getDist4(&(splineFIBERS->fibers[i])));
    }

    if(red)
    {
        delete splineFIBERS;//->~CtlCtlStruct();
    }
    delete tool;//->~C_parametric_comparison();
    fclose(f);

    return;
}

double C_measure_spline_synthetic::fiberLength(CtlStruct* phi)
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
    //cout << phi->elts[3*0] << " " << phi->elts[3*0+1] << " " << phi->elts[3*0+2] << endl;
    for(unsigned long i=1 ; i<phi->N_element-1 ; i++)
    {
        //cout << phi->elts[3*(i+1)] << " " << phi->elts[3*(i+1)+1] << " " << phi->elts[3*(i+1)+2] << endl;
        //cout << phi->elts[3*(i-1)] << " " << phi->elts[3*(i-1)+1] << " " << phi->elts[3*(i-1)+2] << endl;
        DPHI[i] = sqrt(SQR(phi->elts[3*(i+1)]-phi->elts[3*(i-1)]) + SQR(phi->elts[3*(i+1)+1]-phi->elts[3*(i-1)+1]) + SQR(phi->elts[3*(i+1)+2]-phi->elts[3*(i-1)+2]))/(ds);
        //cout << DPHI[i] << endl;
    }
    DPHI[phi->N_element-1] = sqrt(SQR(phi->elts[3*(phi->N_element-1)]-phi->elts[3*(phi->N_element-2)]) + SQR(phi->elts[3*(phi->N_element-1)+1]-phi->elts[3*(phi->N_element-2)+1]) + SQR(phi->elts[3*(phi->N_element-1)+2]-phi->elts[3*(phi->N_element-2)+2]) )/(ds/2.0);
    //cout << phi->elts[3*(phi->N_element-1)] << " " << phi->elts[3*(phi->N_element-1)+1] << " " << phi->elts[3*(phi->N_element-1)+2] << endl;
    ///integrate the norm of derivative parametric curve over the whole arc
    C_toolbox_integrate* t = new C_toolbox_integrate();

    double L = t->run(DPHI, phi->N_element, 1.0, SIMPSON_RULE, NULL);

    delete t;//->~C_integrate();
    delete DPHI;
    return L;
}

double C_measure_spline_synthetic::fiberError(CtlStruct* phi)
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

    return eps/( ((double) phi->N_element)*L );//SQR(L);
    ///return eps/SQR(L);
}

void C_measure_spline_synthetic::syntheticTensor(double x, double y, double z, double* T)
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
