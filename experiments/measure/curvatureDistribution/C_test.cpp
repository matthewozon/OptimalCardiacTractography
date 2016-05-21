#include <C_test.h>

C_test::C_test()
{
    //ctor
}

C_test::~C_test()
{
    //dtor
}

void C_test::testSynthetic(void)
{
    unsigned long Nx = 10;
    unsigned long Ny = 10;
    double dx = 1.0, dy=1.5;
    double x0 = -5.0, y0 = -7.5;
    double* T = new double[6];
    C_measure_synthetic* MS = new C_measure_synthetic();
    MS->helix_angle = Pi/2;
    MS->lambda1 = 1.0;
    MS->lambda2 = 0.2;
    MS->lambda3 = 0.2;
    FILE* f = fopen("angle90.data","w");
    fprintf(f, "%%Point, tensor component, first eigen vector, first eigen value, curvature\n");
    double** a = new double*[3];
    a[0] = new double[3];
    a[1] = new double[3];
    a[2] = new double[3];
    for(unsigned long i=0 ; i<Nx ; i++)
    {
        for(unsigned long j=0 ; j<Ny ; j++)
        {
            ///compute tensor field
            if(  !(((x0 + ((double) i)*dx)==0.0) && ((y0 + ((double) j)*dy)==0.0)))
            {
                MS->syntheticTensor(x0 + ((double) i)*dx, y0 + ((double) j)*dy, 0, T);
                a[0][0] = T[0]; a[0][1] = T[1]; a[0][2] = T[2];
                a[1][0] = T[1]; a[1][1] = T[3]; a[1][2] = T[4];
                a[2][0] = T[2]; a[2][1] = T[4]; a[2][2] = T[5];
                C_toolbox_eigen_sym* S = new C_toolbox_eigen_sym(a);
                fprintf(f, "%f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f\n", x0 + ((double) i)*dx, y0 + ((double) j)*dy, 0.0, T[0], T[1], T[2], T[3], T[4], T[5], S->z[0][0], S->z[1][0], S->z[2][0], S->d[0], MS->syntheticCurvature(x0 + ((double) i)*dx, y0 + ((double) j)*dy, 0.0));
                delete S;
            }

        }
    }

    fclose(f);
    delete [] a[0];
    delete [] a[1];
    delete [] a[2];
    delete [] a;
    delete [] T;
    delete MS;
    return;
}

void C_test::testSpline(void)
{
    unsigned long nbFib = 1;
    unsigned long N=1000;
    CtlCtlStruct* FIBER = new CtlCtlStruct(nbFib);

    double t;
    for(unsigned long j=0 ; j<nbFib ; j++)
    {
        FIBER->fibers[j].allocateArray(N);
        for(unsigned long i=0 ; i<N ; i++)
        {
            t = ((double) i)/((double) (N-1)); //0.25*
            FIBER->fibers[j].elts[3*i] = (1.0/(2*Pi))*cos(8*Pi*t);
            FIBER->fibers[j].elts[3*i+1] = (1.0/(2*Pi))*sin(8*Pi*t);
            FIBER->fibers[j].elts[3*i+2] = 2*t;
        }
    }
    FILE* f = fopen("fiber.data","w");
    C_spline* S = new C_spline();
    CtlCtlStruct* curveSpline;
    CtlCtlStruct* dphiSpline;
    CtlCtlStruct* ddphiSpline;
    S->computeFiberSpline(FIBER, &curveSpline, 1, CURVE);
    S->computeFiberSpline(FIBER, &dphiSpline, 1, SPEED);
    S->computeFiberSpline(FIBER, &ddphiSpline, 1, ACCELERATION);

    for(unsigned long j=0 ; j<nbFib ; j++)
    {
        FIBER->fibers[j].allocateArray(N);
        for(unsigned long i=0 ; i<curveSpline->fibers[j].N_element ; i++)
        {
            //cout << curveSpline->fibers[j].elts[3*i] << " " << curveSpline->fibers[j].elts[3*i+1] << " " << curveSpline->fibers[j].elts[3*i+2] << endl;
            fprintf(f,"%Lf, %Lf, %Lf, %Lf, %Lf, %Lf, %Lf, %Lf, %Lf\n", curveSpline->fibers[j].elts[3*i], curveSpline->fibers[j].elts[3*i+1], curveSpline->fibers[j].elts[3*i+2], dphiSpline->fibers[j].elts[3*i], dphiSpline->fibers[j].elts[3*i+1], dphiSpline->fibers[j].elts[3*i+2], ddphiSpline->fibers[j].elts[3*i], ddphiSpline->fibers[j].elts[3*i+1], ddphiSpline->fibers[j].elts[3*i+2]);
        }
        //fprintf(f, "%f, %f, %f, %f, %f, %f, %f, %f, %f\n", 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
    }



    delete curveSpline;
    delete dphiSpline;
    delete ddphiSpline;
    delete FIBER;
    fclose(f);
    return;
}

void C_test::testFiberLength(void)
{
    unsigned long nbFib = 1;
    unsigned long N=15;
    CtlCtlStruct* FIBER = new CtlCtlStruct(nbFib);

    double t;
    for(unsigned long j=0 ; j<nbFib ; j++)
    {
        FIBER->fibers[j].allocateArray(N);
        for(unsigned long i=0 ; i<N ; i++)
        {
            t = ((double) i)/((double) (N-1));
            FIBER->fibers[j].elts[3*i] = t;//(1.0/(2*pi))*cos(2*pi*t);//(1.0/(2*pi))* cos(8*pi*t);//t; //
            FIBER->fibers[j].elts[3*i+1] = 0.0;//(1.0/(2*pi))*sin(2*pi*t);//(1.0/(2*pi))*sin(8*pi*t); //0.0
            FIBER->fibers[j].elts[3*i+2] = 0.0;//4.0*t;
        }
    }
    C_spline* S = new C_spline();
    CtlCtlStruct* curveSpline;
    CtlCtlStruct* dphiSpline;
    CtlCtlStruct* ddphiSpline;
    S->computeFiberSpline(FIBER, &curveSpline, 1, CURVE);
    S->computeFiberSpline(FIBER, &dphiSpline, 1, SPEED);
    S->computeFiberSpline(FIBER, &ddphiSpline, 1, ACCELERATION);

    C_measure_synthetic* MS = new C_measure_synthetic();
    MS->helix_angle = 0.0*Pi/4.0;
    MS->lambda1 = 1.0;
    MS->lambda2 = 0.2;
    MS->lambda3 = 0.2;
    cout << "fiber length " << MS->meanFiberLength(FIBER) << endl;
    cout << "fiber length " << MS->meanFiberLength(curveSpline) << endl;
    cout << "fiber length " << MS->meanFiberLengthS(dphiSpline) << endl;

    ///destroy
    delete MS;
    delete curveSpline;
    delete dphiSpline;
    delete ddphiSpline;
    delete FIBER;
    return;
}


void C_test::testFiberCurvature(void)
{
    unsigned long nbFib = 1;
    unsigned long N=100;
    CtlCtlStruct* FIBER = new CtlCtlStruct(nbFib);

    double t;
    double alpha = Pi/4.0;
    for(unsigned long j=0 ; j<nbFib ; j++)
    {
        FIBER->fibers[j].allocateArray(N);
        for(unsigned long i=0 ; i<N ; i++)
        {
            t = ((long double) i)/((long double) (N-1));
            FIBER->fibers[j].elts[3*i] = cos(cos(alpha)*t);//t*t;// 10.0*   (1.0/(2.0*pi))*
            FIBER->fibers[j].elts[3*i+1] = -sin(cos(alpha)*t);//sin(2*pi*t); //t;// 10.0*    (1.0/(2.0*pi))*
            FIBER->fibers[j].elts[3*i+2] = sin(alpha)*t;
        }
    }
    C_spline* S = new C_spline();
    CtlCtlStruct* curveSpline;
    CtlCtlStruct* dphiSpline;
    CtlCtlStruct* ddphiSpline;
    S->computeFiberSpline(FIBER, &curveSpline, 1, CURVE);
    S->computeFiberSpline(FIBER, &dphiSpline, 1, SPEED);
    S->computeFiberSpline(FIBER, &ddphiSpline, 1, ACCELERATION);

    C_measure_synthetic* MS = new C_measure_synthetic();
    MS->helix_angle = Pi/4.0;
    MS->lambda1 = 1.0;
    MS->lambda2 = 0.2;
    MS->lambda3 = 0.2;
    cout << "fiber curvature " << MS->meanFiberCurvature(FIBER) << endl;
    cout << "fiber curvature " << MS->meanFiberCurvature(curveSpline) << endl;
    cout << "fiber curvature " << MS->meanFiberCurvatureS(curveSpline, dphiSpline, ddphiSpline) << endl;

    ///destroy
    delete MS;
    delete curveSpline;
    delete dphiSpline;
    delete ddphiSpline;
    delete FIBER;
    return;
}

void C_test::testFiberCurvatureCompare(void)
{
    unsigned long nbFib = 1;
    unsigned long N=15;
    CtlCtlStruct* FIBER = new CtlCtlStruct(nbFib);

    double t;
    double alpha=Pi/4.0;
    for(unsigned long j=0 ; j<nbFib ; j++)
    {
        FIBER->fibers[j].allocateArray(N);
        for(unsigned long i=0 ; i<N ; i++)
        {
            t = ((long double) i)/((long double) (N-1));
            FIBER->fibers[j].elts[3*i] = cos(cos(alpha)*t);//t*t;// 10.0* //t+1;//
            FIBER->fibers[j].elts[3*i+1] = sin(cos(alpha)*t); //t;// 10.0*0;//
            FIBER->fibers[j].elts[3*i+2] = sin(alpha)*t;
        }
    }
    C_spline* S = new C_spline();
    CtlCtlStruct* curveSpline;
    CtlCtlStruct* dphiSpline;
    CtlCtlStruct* ddphiSpline;
    S->computeFiberSpline(FIBER, &curveSpline, 1, CURVE);
    S->computeFiberSpline(FIBER, &dphiSpline, 1, SPEED);
    S->computeFiberSpline(FIBER, &ddphiSpline, 1, ACCELERATION);

    C_measure_synthetic* MS = new C_measure_synthetic();
    MS->helix_angle = Pi/4.0;
    MS->lambda1 = 1.0;
    MS->lambda2 = 0.2;
    MS->lambda3 = 0.2;
    cout << "fiber curvature compared " << MS->meanFiberCurvatureCompare(FIBER) << endl;
    cout << "fiber curvature compared " << MS->meanFiberCurvatureCompare(curveSpline) << endl;
    cout << "fiber curvature compared " << MS->meanFiberCurvatureCompareS(curveSpline, dphiSpline, ddphiSpline) << endl;

    ///destroy
    delete MS;
    delete curveSpline;
    delete dphiSpline;
    delete ddphiSpline;
    delete FIBER;
    return;
}


void C_test::testComparisonStreamPoint(void)
{
    C_parametric_comparison* tool = new C_parametric_comparison();
    tool->alpha = Pi/4.0;
    FILE* f = fopen("comparisonStream.data","w");

    ///create points from a streamline with r0 = 1, z0 = 1 at t = [0:0.01:1]*Dz*sqrt(1+cos(alpha)^2)/sin(alpha)
    double Dz = 3.0, z0 = 1, r0 = 1;
    double x, y, z;
    double dt = Dz*sqrt(1+SQR(cos(tool->alpha)))/sin(tool->alpha);
    for(unsigned long i=0 ; i<101 ; i++)
    {
        tool->getStreamPoint( ((double) i)*dt/100.0, r0, z0, &x, &y, &z);
        fprintf(f, "%f, %f, %f\n", x, y, z);
    }


    fclose(f);
    delete tool;
}

void C_test::testComparisonRadius(void)
{
    C_parametric_comparison* tool = new C_parametric_comparison();
    tool->alpha = Pi/4.0;
    FILE* f = fopen("comparisonStream.data","w");

    ///create points from a streamline with r0 = 1, z0 = 1 at t = [0:0.01:1]*Dz*sqrt(1+cos(alpha)^2)/sin(alpha)
    unsigned long N = 100;
    double Dz = 6.0, z0 = 100, r0 = 1;
    double x, y, z;
    double dt = Dz*sqrt(1+SQR(cos(tool->alpha)))/sin(tool->alpha);
    double* t = new double[N];
    CtlStruct* FIBER = new CtlStruct();
    FIBER->allocateArray(N);

    //.set up
    for(unsigned long i=0 ; i<N ; i++)
    {
        t[i] = ((double) i)*dt/((double) N-1);
        tool->getStreamPoint( t[i], r0, z0, &x, &y, &z);
        x += 0.1*(0.5-((double) rand())/((double) (RAND_MAX-1)));
        y += 0.1*(0.5-((double) rand())/((double) (RAND_MAX-1)));
        z += 0.1*(0.5-((double) rand())/((double) (RAND_MAX-1)));
        fprintf(f, "%f, %f, %f, %f\n", x, y, z, t[i]);
        FIBER->elts[3*i] = x;
        FIBER->elts[3*i+1] = y;
        FIBER->elts[3*i+2] = z;
    }

    //compute radius estimation
    r0 = 0.0;
    for(unsigned long i=0 ; i<FIBER->N_element ; i++)
    {
        r0 += sqrt(SQR(FIBER->elts[3*i]) + SQR(FIBER->elts[3*i+1]));
    }
    r0 = r0/((double) FIBER->N_element);
    r0 = tool->radiusEstimation(FIBER);//, t, r0);
    printf("radius: %f\n", r0);

    //compute z0 estimation
    z0 = 0.0;
    printf("altitude: %f\n", tool->altitudeEstimation(FIBER, t));

    delete FIBER;
    delete t;
    fclose(f);
    delete tool;
}


void C_test::testComparisonTime(void)
{
    C_parametric_comparison* tool = new C_parametric_comparison();
    tool->alpha = Pi/4.0;
    FILE* f = fopen("comparisonStream.data","w");

    ///create points from a streamline with r0 = 1, z0 = 1 at t = [0:0.01:1]*Dz*sqrt(1+cos(alpha)^2)/sin(alpha)
    unsigned long N = 100;
    double z0 = 0.0, r0 = 1;
    double x, y, z;
    double dt = 4.0*Pi/cos(tool->alpha) ;
    double* t = new double[N];
    CtlStruct* FIBER = new CtlStruct();
    FIBER->allocateArray(N);
    double s = 0.1;

    //.set up
    for(unsigned long i=0 ; i<N ; i++)
    {
        t[i] = ((double) i)*dt/((double) N-1);
        tool->getStreamPoint( t[i], r0, z0, &x, &y, &z);
        x += s*(0.5-((double) rand())/((double) (RAND_MAX-1)));
        y += s*(0.5-((double) rand())/((double) (RAND_MAX-1)));
        z += s*(0.5-((double) rand())/((double) (RAND_MAX-1)));
        FIBER->elts[3*i] = (long double) x;
        FIBER->elts[3*i+1] = (long double) y;
        FIBER->elts[3*i+2] = (long double) z;
        fprintf(f, "%f, %f, %f, %f\n", (double) (FIBER->elts[3*i]), (double) (FIBER->elts[3*i+1]), (double) (FIBER->elts[3*i+2]), t[i]);
    }

    //compute time estimation
    FILE* f2 = fopen("comparisonTime.data","w");
    tool->timeEstimation(FIBER, r0, z0, t);
    for(unsigned long i=0 ; i<N ; i++)
    {
        fprintf(f2, "%f\n", t[i]);
    }
    fclose(f2);

    delete FIBER;
    delete t;
    fclose(f);
    delete tool;
}



void C_test::testComparisonTimeAndParam(void)
{
    C_parametric_comparison* tool = new C_parametric_comparison();
    tool->alpha = Pi/4.0;
    FILE* f = fopen("comparisonStream.data","w");

    ///create points from a streamline with r0 = 1, z0 = 1 at t = [0:0.01:1]*Dz*sqrt(1+cos(alpha)^2)/sin(alpha)
    unsigned long N = 100;
    double z0 = 0.0, r0 = 1;
    double x, y, z;
    double dt = 4.0*Pi/cos(tool->alpha) ;
    double* t = new double[N];
    CtlStruct* FIBER = new CtlStruct();
    FIBER->allocateArray(N);
    double s = 0.1;

    //.set up
    for(unsigned long i=0 ; i<N ; i++)
    {
        t[i] = ((double) i)*dt/((double) N-1);
        tool->getStreamPoint( t[i], r0, z0, &x, &y, &z);
        fprintf(f, "%f, %f, %f, %f\n", x, y, z, t[i]);
        x += s*2.0*(0.5-((double) rand())/((double) (RAND_MAX-1)));
        y += s*2.0*(0.5-((double) rand())/((double) (RAND_MAX-1)));
        z += s*2.0*(0.5-((double) rand())/((double) (RAND_MAX-1)));
        FIBER->elts[3*i] = (long double) x;
        FIBER->elts[3*i+1] = (long double) y;
        FIBER->elts[3*i+2] = (long double) z;
        fprintf(f, "%f, %f, %f, %f\n", (double) (FIBER->elts[3*i]), (double) (FIBER->elts[3*i+1]), (double) (FIBER->elts[3*i+2]), t[i]);
    }

    //compute time estimation
    FILE* f2 = fopen("comparisonTime.data","w");
    tool->getTimeAndParam(FIBER, &r0, &z0, t);
    printf("%f, %f", r0, z0);
    for(unsigned long i=0 ; i<N ; i++)
    {
        tool->getStreamPoint( t[i], r0, z0, &x, &y, &z);
        fprintf(f2, "%f, %f, %f, %f\n", x, y, z, t[i]);
    }
    fclose(f2);

    delete FIBER;
    delete t;
    fclose(f);
    delete tool;
}

void C_test::testComparisonDist(void)
{
    C_parametric_comparison* tool = new C_parametric_comparison();
    tool->alpha = Pi/4.0;
    FILE* f = fopen("comparisonStream.data","w");

    ///create points from a streamline with r0 = 1, z0 = 1 at t = [0:0.01:1]*Dz*sqrt(1+cos(alpha)^2)/sin(alpha)
    unsigned long N = 100;
    double z0 = 0.0, r0 = 1;
    double x, y, z;
    double dt = 4.0*Pi/cos(tool->alpha) ;
    double* t = new double[N];
    CtlStruct* FIBER = new CtlStruct();
    FIBER->allocateArray(N);
    double s = 0.0;//0.005;

    //.set up
    for(unsigned long i=0 ; i<N ; i++)
    {
        t[i] = ((double) i)*dt/((double) N-1);
        tool->getStreamPoint( t[i], r0, z0, &x, &y, &z);
        fprintf(f, "%f, %f, %f, %f\n", x, y, z, t[i]);
        x += s*2.0*(0.5-((double) rand())/((double) (RAND_MAX-1)));
        y += s*2.0*(0.5-((double) rand())/((double) (RAND_MAX-1)));
        z += s*2.0*(0.5-((double) rand())/((double) (RAND_MAX-1)));
        FIBER->elts[3*i] = (long double) x;
        FIBER->elts[3*i+1] = (long double) y;
        FIBER->elts[3*i+2] = (long double) z;
        fprintf(f, "%f, %f, %f, %f\n", (double) (FIBER->elts[3*i]), (double) (FIBER->elts[3*i+1]), (double) (FIBER->elts[3*i+2]), t[i]);
    }

    //compute time estimation
    FILE* f2 = fopen("comparisonTime.data","w");
    tool->getTimeAndParam(FIBER, &r0, &z0, t);
    printf("%f, %f, %f", r0, z0, tool->getDist(FIBER));
    for(unsigned long i=0 ; i<N ; i++)
    {
        tool->getStreamPoint( t[i], r0, z0, &x, &y, &z);
        fprintf(f2, "%f, %f, %f, %f\n", x, y, z, t[i]);
    }
    fclose(f2);

    delete FIBER;
    delete t;
    fclose(f);
    delete tool;
}



void C_test::testComparisonTimeAndParamReal(void)
{
    C_parametric_comparison* tool = new C_parametric_comparison();
    tool->alpha = Pi/8.0;
    FILE* f = fopen("comparisonStream.data","w");

    ///create points from a streamline with r0 = 1, z0 = 1 at t = [0:0.01:1]*Dz*sqrt(1+cos(alpha)^2)/sin(alpha)
    C_load_fiber* LOAD = new C_load_fiber();
    CtlCtlStruct* FIBERS = LOAD->readFiber("SA_MOVE_1_U_15_alpha_0.2_beta_0.5_sigma_0.375_t0_1.5708_phi0_0.5_phi1_0.5_iter_96000_Tinit_0.9_Tend_0.04_M_100.fibSRC");
    C_spline* spline = new C_spline();
    CtlCtlStruct* m_ctlCtlStructSpline;
    if(!spline->computeFiberSpline(FIBERS, &m_ctlCtlStructSpline, 5, CURVE))
    {
        fclose(f);
        delete tool;
        delete LOAD;
        delete FIBERS;
        delete spline;
        cout << "false" << endl;
        return;
    }
    CtlStruct* FIBER = &(m_ctlCtlStructSpline->fibers[25]); //2 5 13 15 16 25 26!


    cout << FIBER->N_element << endl;
    unsigned long N = FIBER->N_element;
    double z0 = 0.0, r0 = 0;
    double x, y, z;
    double* t = new double[N];

    //

    //compute time estimation
    FILE* f2 = fopen("comparisonTime.data","w");
    tool->getTimeAndParam(FIBER, &r0, &z0, t);
    printf("%f, %f", r0, z0);
    for(unsigned long i=0 ; i<N ; i++)
    {
        fprintf(f, "%f, %f, %f, %f\n", (double) (FIBER->elts[3*i]), (double) (FIBER->elts[3*i+1]), (double) (FIBER->elts[3*i+2]), t[i]);
        tool->getStreamPoint( t[i], r0, z0, &x, &y, &z);
        fprintf(f2, "%f, %f, %f, %f\n", x, y, z, t[i]);
    }
    fclose(f2);

    delete spline;
    delete LOAD;
    delete FIBER;
    delete t;
    fclose(f);
    delete tool;
}

void C_test::testComparisonDistReal(void)
{
    C_parametric_comparison* tool = new C_parametric_comparison();
    tool->alpha = Pi/8.0;
    FILE* f = fopen("comparisonStream2.data","w");

    ///create points from a streamline with r0 = 1, z0 = 1 at t = [0:0.01:1]*Dz*sqrt(1+cos(alpha)^2)/sin(alpha)
    C_load_fiber* LOAD = new C_load_fiber();
    //CtlCtlStruct* FIBERS = LOAD->readFiber("SA_MOVE_1_U_15_alpha_0.2_beta_0.5_sigma_0.375_t0_1.5708_phi0_0.5_phi1_0.5_iter_96000_Tinit_0.9_Tend_0.04_M_100.fibSRC");
    CtlCtlStruct* FIBERS = LOAD->readFiber("/home/ozon/Documents/U15/synthetic/streamline/streamline_snr0phi22e/streamline_fiber_step_0.1_nbseed_836_minspeed_0.25_maxlength_50.fibSRC");

    CtlStruct* FIBER = &(FIBERS->fibers[7]);

    cout << FIBER->N_element << endl;
    unsigned long N = FIBER->N_element;
    double z0 = 0.0, r0 = 0;
    double x, y, z;
    double* t = new double[N];

    //

    //compute time estimation
    FILE* f2 = fopen("comparisonTime2.data","w");
    tool->getTimeAndParam(FIBER, &r0, &z0, t);
    printf("%f, %f, %f", r0, z0, tool->getDist(FIBER));
    for(unsigned long i=0 ; i<N ; i++)
    {
        fprintf(f, "%f, %f, %f, %f\n", (double) (FIBER->elts[3*i]), (double) (FIBER->elts[3*i+1]), (double) (FIBER->elts[3*i+2]), t[i]);
        tool->getStreamPoint( t[i], r0, z0, &x, &y, &z);
        fprintf(f2, "%f, %f, %f, %f\n", x, y, z, t[i]);
    }
    fclose(f2);

    delete LOAD;
    delete FIBER;
    delete t;
    fclose(f);
    delete tool;
}


void C_test::testHisto(void)
{
    C_load_fiber* LOAD = new C_load_fiber();
    bool red;

    red = true;
    //CtlCtlStruct* FIBERS = LOAD->readFiber("/home/ozon/Documents/U15/robustU15/SA_MOVE_1_U_15_alpha_0.2_beta_0.5_sigma_0.375_t0_1.5708_phi0_0.5_phi1_0.5_iter_96000_Tinit_0.9_Tend_0.04_M_100.fibSRC");
    //CtlCtlStruct* FIBERS = LOAD->readFiber("/home/ozon/Documents/U15/robustU15/SA_MOVE_1_U_15_alpha_0.2_beta_0.5_sigma_0.375_t0_1.5708_phi0_0.5_phi1_0.5_iter_96000_Tinit_0.9_Tend_0.045_M_20.fibSRC");
    //CtlCtlStruct* FIBERS = LOAD->readFiber("/home/ozon/Documents/U15/robustU15/SA_MOVE_1_U_15_alpha_0.2_beta_0.5_sigma_0.375_t0_1.5708_phi0_0.5_phi1_0.5_iter_96000_Tinit_1.1_Tend_0.05_M_10.fibSRC");
    //CtlCtlStruct* FIBERS = LOAD->readFiber("/home/ozon/Documents/U15/robustU15/SA_MOVE_1_U_15_alpha_0.2_beta_0.5_sigma_0.375_t0_1.5708_phi0_0.5_phi1_0.5_iter_96000_Tinit_1.2_Tend_0.025_M_0.fibSRC");

    //CtlCtlStruct* FIBERS = LOAD->readFiber("/home/ozon/Documents/U15/robustU15largeSNR100/SA_MOVE_1_U_15_alpha_0.2_beta_0.5_sigma_0.375_t0_1.5708_phi0_0.5_phi1_0.5_iter_96000_Tinit_0.576429_Tend_0.0133984_M__crop.fibSRC");
    //CtlCtlStruct* FIBERS = LOAD->readFiber("/home/ozon/Documents/U15/robustU15largeSNR20/SA_MOVE_1_U_15_alpha_0.2_beta_0.5_sigma_0.375_t0_1.5708_phi0_0.5_phi1_0.5_iter_96000_Tinit_0.54243_Tend_0.0158935_M__crop.fibSRC");
    //CtlCtlStruct* FIBERS = LOAD->readFiber("/home/ozon/Documents/U15/robustU15largeSNR10/SA_MOVE_1_U_15_alpha_0.2_beta_0.5_sigma_0.375_t0_1.5708_phi0_0.5_phi1_0.5_iter_96000_Tinit_0.653712_Tend_0.0165459_M__crop.fibSRC");
    //CtlCtlStruct* FIBERS = LOAD->readFiber("/home/ozon/Documents/U15/robustU15largeSNR0/SA_MOVE_1_U_15_alpha_0.2_beta_0.5_sigma_0.375_t0_1.5708_phi0_0.5_phi1_0.5_iter_96000_Tinit_0.915315_Tend_0.0135332_M__crop.fibSRC");

    //CtlCtlStruct* FIBERS = LOAD->readFiber("/home/ozon/Documents/U15/robustU1/SA_MOVE_1_U_1_alpha_1_phi0_0.5_phi1_1_iter_96000_Tinit_11.2_Tend_0.05_M_100.fibSRC");
    //CtlCtlStruct* FIBERS = LOAD->readFiber("/home/ozon/Documents/U15/robustU1/SA_MOVE_1_U_1_alpha_1_phi0_0.5_phi1_1_iter_96000_Tinit_11.2_Tend_0.06_M_20.fibSRC");
    //CtlCtlStruct* FIBERS = LOAD->readFiber("/home/ozon/Documents/U15/robustU1/SA_MOVE_1_U_1_alpha_1_phi0_0.5_phi1_1_iter_96000_Tinit_10.8_Tend_0.06_M_10.fibSRC");
    //CtlCtlStruct* FIBERS = LOAD->readFiber("/home/ozon/Documents/U15/robustU1/SA_MOVE_1_U_1_alpha_1_phi0_0.5_phi1_1_iter_96000_Tinit_10.5_Tend_0.1_M_0.fibSRC");

    //CtlCtlStruct* FIBERS = LOAD->readFiber("/home/ozon/Documents/U15/robustU1large/SA_MOVE_1_U_1_alpha_1_phi0_0.5_phi1_1_iter_96000_Tinit_8.77489_Tend_0.0939965_M_100_crop.fibSRC");
    //CtlCtlStruct* FIBERS = LOAD->readFiber("/home/ozon/Documents/U15/robustU1large/SA_MOVE_1_U_1_alpha_1_phi0_0.5_phi1_1_iter_96000_Tinit_8.77739_Tend_0.0828589_M_20_crop.fibSRC");
    //CtlCtlStruct* FIBERS = LOAD->readFiber("/home/ozon/Documents/U15/robustU1large/SA_MOVE_1_U_1_alpha_1_phi0_0.5_phi1_1_iter_96000_Tinit_8.48159_Tend_0.045248_M_10_crop.fibSRC");
    //CtlCtlStruct* FIBERS = LOAD->readFiber("/home/ozon/Documents/U15/robustU1large/SA_MOVE_1_U_1_alpha_1_phi0_0.5_phi1_1_iter_96000_Tinit_8.16043_Tend_0.04816_M_0_crop.fibSRC");

    //CtlCtlStruct* FIBERS = LOAD->readFiber("/home/ozon/Documents/U15/robustU1large2/SNR100/SA_MOVE_1_U_1_alpha_1_phi0_0.5_phi1_1_iter_96000_Tinit_8.82513_Tend_0.0693951_M__crop.fibSRC");
    //CtlCtlStruct* FIBERS = LOAD->readFiber("/home/ozon/Documents/U15/robustU1large2/SNR20/SA_MOVE_1_U_1_alpha_1_phi0_0.5_phi1_1_iter_96000_Tinit_8.82905_Tend_0.0466853_M__crop.fibSRC");
    //CtlCtlStruct* FIBERS = LOAD->readFiber("/home/ozon/Documents/U15/robustU1large2/SNR10/SA_MOVE_1_U_1_alpha_1_phi0_0.5_phi1_1_iter_96000_Tinit_8.54592_Tend_0.176652_M__crop.fibSRC");
    CtlCtlStruct* FIBERS = LOAD->readFiber("/home/ozon/Documents/U15/robustU1large2/SNR0/SA_MOVE_1_U_1_alpha_1_phi0_0.5_phi1_1_iter_96000_Tinit_8.10986_Tend_0.0483661_M__crop.fibSRC");

    //red = false;
    //CtlCtlStruct* FIBERS = LOAD->readFiber("/home/ozon/Documents/U15/synthetic/streamline/streamline_snr100phi22e/streamline_fiber_step_0.1_nbseed_836_minspeed_0.25_maxlength_50.fibSRC");
    //CtlCtlStruct* FIBERS = LOAD->readFiber("/home/ozon/Documents/U15/synthetic/streamline/streamline_snr20phi22e/streamline_fiber_step_0.1_nbseed_836_minspeed_0.25_maxlength_50.fibSRC");
    //CtlCtlStruct* FIBERS = LOAD->readFiber("/home/ozon/Documents/U15/synthetic/streamline/streamline_snr10phi22e/streamline_fiber_step_0.1_nbseed_836_minspeed_0.25_maxlength_50.fibSRC");
    //CtlCtlStruct* FIBERS = LOAD->readFiber("/home/ozon/Documents/U15/synthetic/streamline/streamline_snr0phi22e/streamline_fiber_step_0.1_nbseed_836_minspeed_0.25_maxlength_50.fibSRC");
    ///
    C_measure_spline_synthetic* tool = new C_measure_spline_synthetic();
    tool->helix_angle = Pi/8.0;
    tool->lambda1 = 1;
    tool->lambda2 = 0.2;
    tool->lambda3 = 0.2;

    cout << FIBERS->N << endl;

    tool->getLenAndErrPop(FIBERS, red);

    delete tool;
    delete LOAD;
    return;
}


#ifdef CROP

void C_test::testCrop(void)
{
    floating Rmin=8.5;
    floating Rmax=14.5;
    floating Zmin=-10.0;
    floating Zmax=10.0;
    C_crop* tool = new C_crop(Rmin, Rmax, Zmin, Zmax);
    //string savedEdge = "/home/ozon/Documents/U15/robustU15largeSNR100/SA_MOVE_1_U_15_alpha_0.2_beta_0.5_sigma_0.375_t0_1.5708_phi0_0.5_phi1_0.5_iter_96000_Tinit_0.576429_Tend_0.0133984_M_.cformat";
    //string savedEdge = "/home/ozon/Documents/U15/robustU15largeSNR30/SA_MOVE_1_U_15_alpha_0.2_beta_0.5_sigma_0.375_t0_1.5708_phi0_0.5_phi1_0.5_iter_96000_Tinit_0.580642_Tend_0.01331_M_.cformat";
    //string savedEdge = "/home/ozon/Documents/U15/robustU15largeSNR25/SA_MOVE_1_U_15_alpha_0.2_beta_0.5_sigma_0.375_t0_1.5708_phi0_0.5_phi1_0.5_iter_96000_Tinit_0.558773_Tend_0.0123249_M_.cformat";
    //string savedEdge = "/home/ozon/Documents/U15/robustU15largeSNR20/SA_MOVE_1_U_15_alpha_0.2_beta_0.5_sigma_0.375_t0_1.5708_phi0_0.5_phi1_0.5_iter_96000_Tinit_0.54243_Tend_0.0158935_M_.cformat";
    //string savedEdge = "/home/ozon/Documents/U15/robustU15largeSNR15/SA_MOVE_1_U_15_alpha_0.2_beta_0.5_sigma_0.375_t0_1.5708_phi0_0.5_phi1_0.5_iter_96000_Tinit_0.5525_Tend_0.0147515_M_.cformat";
    //string savedEdge = "/home/ozon/Documents/U15/robustU15largeSNR10/SA_MOVE_1_U_15_alpha_0.2_beta_0.5_sigma_0.375_t0_1.5708_phi0_0.5_phi1_0.5_iter_96000_Tinit_0.653712_Tend_0.0165459_M_.cformat";
    //string savedEdge = "/home/ozon/Documents/U15/robustU15largeSNR5/SA_MOVE_1_U_15_alpha_0.2_beta_0.5_sigma_0.375_t0_1.5708_phi0_0.5_phi1_0.5_iter_96000_Tinit_0.584744_Tend_0.017574_M_.cformat";
    //string savedEdge = "/home/ozon/Documents/U15/robustU15largeSNR0/SA_MOVE_1_U_15_alpha_0.2_beta_0.5_sigma_0.375_t0_1.5708_phi0_0.5_phi1_0.5_iter_96000_Tinit_0.915315_Tend_0.0135332_M_.cformat";

    //string savedEdge = "/home/ozon/Documents/U15/robustU1large/SA_MOVE_1_U_1_alpha_1_phi0_0.5_phi1_1_iter_96000_Tinit_8.74702_Tend_0.0822822_M_.cformat";
    //string savedEdge = "/home/ozon/Documents/U15/robustU1large/SA_MOVE_1_U_1_alpha_1_phi0_0.5_phi1_1_iter_96000_Tinit_8.77739_Tend_0.0828589_M_.cformat";
    //string savedEdge = "/home/ozon/Documents/U15/robustU1large/SA_MOVE_1_U_1_alpha_1_phi0_0.5_phi1_1_iter_96000_Tinit_8.77784_Tend_0.0850726_M_.cformat";
    //string savedEdge = "/home/ozon/Documents/U15/robustU1large/SA_MOVE_1_U_1_alpha_1_phi0_0.5_phi1_1_iter_96000_Tinit_8.77024_Tend_0.0795298_M_.cformat";
    //string savedEdge = "/home/ozon/Documents/U15/robustU1large/SA_MOVE_1_U_1_alpha_1_phi0_0.5_phi1_1_iter_96000_Tinit_8.77489_Tend_0.0939965_M_.cformat";
    //string savedEdge = "/home/ozon/Documents/U15/robustU1large/SA_MOVE_1_U_1_alpha_1_phi0_0.5_phi1_1_iter_96000_Tinit_8.48159_Tend_0.045248_M_.cformat";
    //string savedEdge = "/home/ozon/Documents/U15/robustU1large/SA_MOVE_1_U_1_alpha_1_phi0_0.5_phi1_1_iter_96000_Tinit_8.16043_Tend_0.04816_M_.cformat";
    //string savedEdge = "/home/ozon/Documents/U15/robustU1large/SA_MOVE_1_U_1_alpha_1_phi0_0.5_phi1_1_iter_96000_Tinit_8.22369_Tend_0.0739264_M_.cformat";


    //string savedEdge = "/home/ozon/Documents/U15/robustU1large2/SNR100/SA_MOVE_1_U_1_alpha_1_phi0_0.5_phi1_1_iter_96000_Tinit_8.82513_Tend_0.0693951_M_.cformat";
    //string savedEdge = "/home/ozon/Documents/U15/robustU1large2/SNR30/SA_MOVE_1_U_1_alpha_1_phi0_0.5_phi1_1_iter_96000_Tinit_8.83059_Tend_0.0692792_M_.cformat";
    //string savedEdge = "/home/ozon/Documents/U15/robustU1large2/SNR25/SA_MOVE_1_U_1_alpha_1_phi0_0.5_phi1_1_iter_96000_Tinit_8.83339_Tend_0.237809_M_.cformat";
    //string savedEdge = "/home/ozon/Documents/U15/robustU1large2/SNR20/SA_MOVE_1_U_1_alpha_1_phi0_0.5_phi1_1_iter_96000_Tinit_8.82905_Tend_0.0466853_M_.cformat";
    //string savedEdge = "/home/ozon/Documents/U15/robustU1large2/SNR15/SA_MOVE_1_U_1_alpha_1_phi0_0.5_phi1_1_iter_96000_Tinit_8.80117_Tend_0.197309_M_.cformat";
    //string savedEdge = "/home/ozon/Documents/U15/robustU1large2/SNR10/SA_MOVE_1_U_1_alpha_1_phi0_0.5_phi1_1_iter_96000_Tinit_8.54592_Tend_0.176652_M_.cformat";
    //string savedEdge = "/home/ozon/Documents/U15/robustU1large2/SNR5/SA_MOVE_1_U_1_alpha_1_phi0_0.5_phi1_1_iter_96000_Tinit_8.217_Tend_0.0487047_M_.cformat";
    string savedEdge = "/home/ozon/Documents/U15/robustU1large2/SNR0/SA_MOVE_1_U_1_alpha_1_phi0_0.5_phi1_1_iter_96000_Tinit_8.10986_Tend_0.0483661_M_.cformat";

    tool->cropHelixFiber(savedEdge);
    delete tool;
    return;
}

#endif
