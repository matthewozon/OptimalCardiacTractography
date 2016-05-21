#include <C_measure.h>

C_measure::C_measure()
{
    //ctor
    tool = new C_toolbox_integrate();
}

C_measure::~C_measure()
{
    //dtor
    delete tool;
}


//specific tools
void C_measure::fillArray(double x, double y, double z, double** a)
{
    if(m_rawTensor==NULL) return;
    signed long i = (signed long) (x/m_rawTensor->pixDimX);
    signed long j = (signed long) (y/m_rawTensor->pixDimY);
    signed long k = (signed long) (z/m_rawTensor->pixDimZ);

    if(i<0 || i>=(signed long)m_rawTensor->DimX || j<0 || j>=(signed long)m_rawTensor->DimY || k<0 || k>=(signed long)m_rawTensor->DimZ)
    {
        cout << "tensor out of bound " << i << " " << j << " " << k << endl;
        a[0][0] = 1.0;
        a[0][1] = 0.0;
        a[0][2] = 0.0;

        a[1][0] = 0.0;
        a[1][1] = 1.0;
        a[1][2] = 0.0;

        a[2][0] = 0.0;
        a[2][1] = 0.0;
        a[2][2] = 1.0;
    }
    else if(m_rawMask->raw3D[i][j][k]<0.5)
    {
        cout << "tensor out of mask" << i << " " << j << " " << k << endl;
        a[0][0] = 1.0;
        a[0][1] = 0.0;
        a[0][2] = 0.0;

        a[1][0] = 0.0;
        a[1][1] = 1.0;
        a[1][2] = 0.0;

        a[2][0] = 0.0;
        a[2][1] = 0.0;
        a[2][2] = 1.0;
    }
    else
    {
        a[0][0] = m_rawTensor->raw4D[i][j][k][0];
        a[0][1] = m_rawTensor->raw4D[i][j][k][1];
        a[0][2] = m_rawTensor->raw4D[i][j][k][2];

        a[1][0] = m_rawTensor->raw4D[i][j][k][1];
        a[1][1] = m_rawTensor->raw4D[i][j][k][3];
        a[1][2] = m_rawTensor->raw4D[i][j][k][4];

        a[2][0] = m_rawTensor->raw4D[i][j][k][2];
        a[2][1] = m_rawTensor->raw4D[i][j][k][4];
        a[2][2] = m_rawTensor->raw4D[i][j][k][5];
    }
    return;
}

double C_measure::TphiPrimeNorm2(double x, double y, double z, double dx, double dy, double dz)
{
    if(m_rawTensor==NULL) return -1.0;
    signed long i = (signed long) (x/m_rawTensor->pixDimX);
    signed long j = (signed long) (y/m_rawTensor->pixDimY);
    signed long k = (signed long) (z/m_rawTensor->pixDimZ);

    if(i<0 || i>=(signed long)m_rawTensor->DimX || j<0 || j>=(signed long)m_rawTensor->DimY || k<0 || k>=(signed long)m_rawTensor->DimZ) return -1.0;
    if(m_rawMask->raw3D[i][j][k]<0.5) return -1.0;

//    if(isnan(sqrt( SQR( (m_rawTensor->raw4D[i][j][k][0])*dx + (m_rawTensor->raw4D[i][j][k][1])*dy + (m_rawTensor->raw4D[i][j][k][2])*dz ) +\
//                   SQR( (m_rawTensor->raw4D[i][j][k][1])*dx + (m_rawTensor->raw4D[i][j][k][3])*dy + (m_rawTensor->raw4D[i][j][k][4])*dz ) +\
//                   SQR( (m_rawTensor->raw4D[i][j][k][2])*dx + (m_rawTensor->raw4D[i][j][k][4])*dy + (m_rawTensor->raw4D[i][j][k][5])*dz ) )))
//    {
//        cout << "sqrt isnan" << endl;
//        cout << SQR( (m_rawTensor->raw4D[i][j][k][0])*dx + (m_rawTensor->raw4D[i][j][k][1])*dy + (m_rawTensor->raw4D[i][j][k][2])*dz ) +\
//                SQR( (m_rawTensor->raw4D[i][j][k][1])*dx + (m_rawTensor->raw4D[i][j][k][3])*dy + (m_rawTensor->raw4D[i][j][k][4])*dz ) +\
//                SQR( (m_rawTensor->raw4D[i][j][k][2])*dx + (m_rawTensor->raw4D[i][j][k][4])*dy + (m_rawTensor->raw4D[i][j][k][5])*dz ) << endl;
//        cout << (m_rawTensor->raw4D[i][j][k][0])*dx + (m_rawTensor->raw4D[i][j][k][1])*dy + (m_rawTensor->raw4D[i][j][k][2])*dz << endl;
//        cout << (m_rawTensor->raw4D[i][j][k][1])*dx + (m_rawTensor->raw4D[i][j][k][3])*dy + (m_rawTensor->raw4D[i][j][k][4])*dz << endl;
//        cout << (m_rawTensor->raw4D[i][j][k][2])*dx + (m_rawTensor->raw4D[i][j][k][4])*dy + (m_rawTensor->raw4D[i][j][k][5])*dz << endl;
//        cout << "(" << dx << "," << dy << "," << dz << ")" << endl;
//    }

    return sqrt( SQR( (m_rawTensor->raw4D[i][j][k][0])*dx + (m_rawTensor->raw4D[i][j][k][1])*dy + (m_rawTensor->raw4D[i][j][k][2])*dz ) +\
          SQR( (m_rawTensor->raw4D[i][j][k][1])*dx + (m_rawTensor->raw4D[i][j][k][3])*dy + (m_rawTensor->raw4D[i][j][k][4])*dz ) +\
          SQR( (m_rawTensor->raw4D[i][j][k][2])*dx + (m_rawTensor->raw4D[i][j][k][4])*dy + (m_rawTensor->raw4D[i][j][k][5])*dz ) );
}


double C_measure::FA(double* eig)
{
    return sqrt( 0.5 * (SQR(eig[0] - eig[1]) + SQR(eig[1] - eig[2]) + SQR(eig[2] - eig[0])) / (SQR(eig[0]) + SQR(eig[1]) + SQR(eig[2])) );
}


double C_measure::LSplin(/**C_graph* G,*/ CtlStruct* phi)
{
    double dl = 1/((double) phi->N_element - 1);
    double* phiPrimeNorm = new double[phi->N_element];
    phiPrimeNorm[0] = sqrt( SQR((phi->elts[3*(1+1)] - phi->elts[3*1])/dl) +\
                           SQR((phi->elts[3*(1+1)+1] - phi->elts[3*1+1])/dl) +\
                           SQR((phi->elts[3*(1+1)+2] - phi->elts[3*1+2])/dl) );
    for(unsigned long n=1 ; n<phi->N_element-1 ; n++)
    {
        phiPrimeNorm[n] = sqrt( SQR((phi->elts[3*(n+1)] - phi->elts[3*(n-1)])/(2*dl)) +\
                               SQR((phi->elts[3*(n+1)+1] - phi->elts[3*(n-1)+1])/(2*dl)) +\
                               SQR((phi->elts[3*(n+1)+2] - phi->elts[3*(n-1)+2])/(2*dl)));
    }
    phiPrimeNorm[phi->N_element-1] = sqrt( SQR((phi->elts[3*(phi->N_element-1)] - phi->elts[3*(phi->N_element-2)])/dl) +\
                                               SQR((phi->elts[3*(phi->N_element-1)+1] - phi->elts[3*(phi->N_element-2)+1])/dl) +\
                                               SQR((phi->elts[3*(phi->N_element-1)+2] - phi->elts[3*(phi->N_element-2)+2])/dl) );

    //update mean
    double mean = tool->run(phiPrimeNorm, phi->N_element, 1, TRAPEZOIDAL_RULE, NULL);
    delete phiPrimeNorm;
    return mean;
}
double C_measure::meanLengthSplin(CtlCtlStruct* spline)
{
    double mean = 0;
    for(unsigned long i=0 ; i<spline->N ; i++)
    {
        mean += LSplin(&(spline->fibers[i]));
    }
    return mean/((double) spline->N);
}
double C_measure::stdLengthSplin(CtlCtlStruct* spline, double mean)
{
    double std = 0;
    for(unsigned long i=0 ; i<spline->N ; i++)
    {
        double delta = (LSplin(&(spline->fibers[i])) - mean);
        std += SQR(delta);
    }
    return sqrt(std/((double) spline->N));
}

C_point_array_data* C_measure::fillInterpolatePoints(CtlStruct* phi)
{
    C_point_array_data* P = new C_point_array_data(phi->N_element, 3, 3, false);

    //fill coordinates
    for(unsigned long i=0 ; i<phi->N_element ; i++)
    {
        P->pointCoordinate[i][0] = phi->elts[3*i];
        P->pointCoordinate[i][1] = phi->elts[3*i+1];
        P->pointCoordinate[i][2] = phi->elts[3*i+2];
    }

    //fill vector and FA
    C_toolbox_interpolate* t = new C_toolbox_interpolate(m_rawTensor, m_rawMask, 1.5*max(m_rawTensor->pixDimX, max(m_rawTensor->pixDimY, m_rawTensor->pixDimZ))/**sqrt(3)*(G->getLattice())*/);

    if(t->run(P,GAUSSIAN_PHI, TENSOR, false)!=_SUCCESS)
    {
        delete P;
        P = NULL;
    }

    delete t;
    return P;
}

double C_measure::meanLengthErrorSplinE2(CtlCtlStruct* spline, double* E2)
{
    double mean = 0.0, temp, temp2;
    unsigned long b = 0;
    *E2 = 0.0;
    for(unsigned long i=0 ; i<spline->N ; i++)
    {
        temp = lengthError(&(spline->fibers[i]), &temp2);
        if(temp>=0.0)
        {
            mean += temp;
            *E2 += temp2;
        }
        else
        {
            b++;
        }

    }

    if(spline->N > b)
    {
        *E2 = *E2/(((double) spline->N) - ((double) b));
        return mean/(((double) spline->N) - ((double) b));
    }
    else
    {
        *E2 = -1.0;
        return -1.0;
    }

}

double C_measure::lengthError(CtlStruct* phi, double* E2)
{

    double len = LSplin(phi);
    double dl = 1/((double) phi->N_element - 1);
    C_point_array_data* P = fillInterpolatePoints(phi); ///this allocate a C_point_array_data object: must be destroyed
    if(P==NULL)
    {
        *E2 = -1.0;
        return -1.0;
    }
    double* dX = new double[phi->N_element];
    double* dY = new double[phi->N_element];
    double* dZ = new double[phi->N_element];

    //compute derivatives
    dX[0] = (P->pointCoordinate[1][0] - P->pointCoordinate[0][0])/dl;
    dY[0] = (P->pointCoordinate[1][1] - P->pointCoordinate[0][1])/dl;
    dZ[0] = (P->pointCoordinate[1][2] - P->pointCoordinate[0][2])/dl;
    for(unsigned long i=1 ; i<phi->N_element-1 ; i++)
    {
        dX[i] = (P->pointCoordinate[i+1][0] - P->pointCoordinate[i-1][0])/(2*dl);
        dY[i] = (P->pointCoordinate[i+1][1] - P->pointCoordinate[i-1][1])/(2*dl);
        dZ[i] = (P->pointCoordinate[i+1][2] - P->pointCoordinate[i-1][2])/(2*dl);
    }
    dX[phi->N_element-1] = (P->pointCoordinate[phi->N_element-1][0] - P->pointCoordinate[phi->N_element-2][0])/dl;
    dY[phi->N_element-1] = (P->pointCoordinate[phi->N_element-1][1] - P->pointCoordinate[phi->N_element-2][1])/dl;
    dZ[phi->N_element-1] = (P->pointCoordinate[phi->N_element-1][2] - P->pointCoordinate[phi->N_element-2][2])/dl;

    //compute control points
    double* ctrlPoints = new double[phi->N_element];
    for(unsigned long i=0 ; i<phi->N_element ; i++)
    {
        ctrlPoints[i] = P->FA[i] *\
        sqrt(\
             SQR(P->vectorInterpolated[i][1]*dZ[i] - P->vectorInterpolated[i][2]*dY[i]) +\
             SQR(P->vectorInterpolated[i][2]*dX[i] - P->vectorInterpolated[i][0]*dZ[i]) +\
             SQR(P->vectorInterpolated[i][0]*dY[i] - P->vectorInterpolated[i][1]*dX[i])\
             );
    }

    //compute integrale value
    double r = tool->run(ctrlPoints, phi->N_element, 1.0, SIMPSON_RULE, NULL);

    //delete
    delete ctrlPoints;
    delete dZ;
    delete dY;
    delete dX;
    delete P;

    *E2 = r/len;
    return r/(SQR(len));
}

double C_measure::error_fiber(string fileNameFib)
{
    C_load_fiber* loadObj = new C_load_fiber();
    CtlCtlStruct* FIBERS = loadObj->readFiber(fileNameFib, true);
    if(FIBERS==NULL)
    {
        return -1;
    }
    double E = error_fiber(FIBERS);
    delete FIBERS;
    delete loadObj;
    return E;
}

double C_measure::error_fiber_tilde(string fileNameFib)
{
    C_load_fiber* loadObj = new C_load_fiber();
    CtlCtlStruct* FIBERS = loadObj->readFiber(fileNameFib, true);
    if(FIBERS==NULL)
    {
        return -1;
    }
    double E = error_fiber_tilde(FIBERS);
    delete FIBERS;
    delete loadObj;
    return E;
}

double  C_measure::error_fiber(CtlCtlStruct* FIBERS)
{
    double mean=0, tmp=0.0;
    unsigned long b=0;
    for(unsigned long i=0 ; i<FIBERS->N ; i++)
    {
        if(FIBERS->fibers[i].N_element>2)
        {
            tmp = error_fiber(&(FIBERS->fibers[i]));
            if(!isnan(tmp) && (tmp>=0.0)) mean +=tmp;
            else {b++; cout << "bug" << endl;}
        }
        else
        {
            b++;
        }
    }
    if(FIBERS->N>b)
    {
        #ifdef SHOULD_NOT_BE_USED
        return mean/(((double) FIBERS->N) - ((double) b));
        #else
        return 1.0-(mean/(((double) FIBERS->N) - ((double) b)));
        #endif
    }
    else
    {
        return -1.0;//NAN;//-1.0;
    }
}

double C_measure::error_fiber_tilde(CtlCtlStruct* FIBERS)
{
    double mean=0, tmp=0.0;
    unsigned long b=0;
    for(unsigned long i=0 ; i<FIBERS->N ; i++)
    {
        if(FIBERS->fibers[i].N_element>2)
        {
            tmp = error_fiber_tilde(&(FIBERS->fibers[i]));
            if(!isnan(tmp) && (tmp>=0.0)) mean +=tmp;
            else {b++; cout << "bug" << endl;}
        }
        else
        {
            b++;
        }
    }
    if(FIBERS->N>b)
    {
        return mean/(((double) FIBERS->N) - ((double) b));
    }
    else
    {
        return -1.0;//NAN;//-1.0;
    }
}

double C_measure::error_fiber(CtlStruct* phi)
{
#ifdef SHOULD_NOT_BE_USED
    unsigned long N = phi->N_element;
    //Point* p = G->getPoints();
    double TphiPrime;
    double** a = new double*[3];
    a[0] = new double[3];
    a[1] = new double[3];
    a[2] = new double[3];

    //
    C_toolbox_eigen_sym* v;

    double phiPrimeX, phiPrimeY, phiPrimeZ, Z;
    double E = 0;
    double dLength=0.0;
    unsigned long dN=0;

    //iter
    for(unsigned long i=1 ; i<N-1 ; i++)
    {
        fillArray(phi->elts[3*i], phi->elts[3*i+1], phi->elts[3*i+2], a);
        phiPrimeX = phi->elts[3*(i+1)] - phi->elts[3*(i-1)];
        phiPrimeY = phi->elts[3*(i+1)+1] - phi->elts[3*(i-1)+1];
        phiPrimeZ = phi->elts[3*(i+1)+2] - phi->elts[3*(i-1)+2];
        Z = sqrt(SQR(phiPrimeX) + SQR(phiPrimeY) + SQR(phiPrimeZ));
        phiPrimeX /= Z;
        phiPrimeY /= Z;
        phiPrimeZ /= Z;
        if(isnan(phiPrimeX) || isnan(phiPrimeY) || isnan(phiPrimeZ) || isinf(phiPrimeX) || isinf(phiPrimeY) || isinf(phiPrimeZ))
        {
            E = -1.0;
            i = N;
            //cout << "phiPrime" << endl;
        }
        else
        {
            TphiPrime = TphiPrimeNorm2(phi->elts[3*i], phi->elts[3*i+1], phi->elts[3*i+2], phiPrimeX, phiPrimeY, phiPrimeZ);
            if(TphiPrime<0.0)
            {
                dLength += 0.5*Z;//sqrt(SQR(phi->elts[3*(i+1)] - phi->elts[3*(i-1)]) + SQR(phi->elts[3*(i+1)+1] - phi->elts[3*(i-1)+1]) + SQR(phi->elts[3*(i+1)+2] - phi->elts[3*(i-1)+2]));
                //dLength += sqrt(SQR(phi->elts[3*(i)] - phi->elts[3*(i-1)]) + SQR(phi->elts[3*(i)+1] - phi->elts[3*(i-1)+1]) + SQR(phi->elts[3*(i)+2] - phi->elts[3*(i-1)+2]));
                dN++;
                //cout << "TphiPrime" << endl;
            }
            else
            {
                v = new C_toolbox_eigen_sym((/**const*/ double**) a, 3, false);
                if(isnan(1-(TphiPrime/(v->d[0]))) || isinf(1-(TphiPrime/(v->d[0]))))
                {
                    dLength += 0.5*Z;//sqrt(SQR(phi->elts[3*(i+1)] - phi->elts[3*(i-1)]) + SQR(phi->elts[3*(i+1)+1] - phi->elts[3*(i-1)+1]) + SQR(phi->elts[3*(i+1)+2] - phi->elts[3*(i-1)+2]));
                    //dLength += sqrt(SQR(phi->elts[3*(i)] - phi->elts[3*(i-1)]) + SQR(phi->elts[3*(i)+1] - phi->elts[3*(i-1)+1]) + SQR(phi->elts[3*(i)+2] - phi->elts[3*(i-1)+2]));
                    dN++;
                    //cout << "some nan" << endl;
                }
                else
                {
                    if(ABS(v->d[0])>SMALL_NUM)
                        E += (1-(TphiPrime/(v->d[0]))); /// because ||phi'||=1
                    else
                    {
                        dLength += 0.5*Z;//sqrt(SQR(phi->elts[3*(i+1)] - phi->elts[3*(i-1)]) + SQR(phi->elts[3*(i+1)+1] - phi->elts[3*(i-1)+1]) + SQR(phi->elts[3*(i+1)+2] - phi->elts[3*(i-1)+2]));
                        //dLength += sqrt(SQR(phi->elts[3*(i)] - phi->elts[3*(i-1)]) + SQR(phi->elts[3*(i)+1] - phi->elts[3*(i-1)+1]) + SQR(phi->elts[3*(i)+2] - phi->elts[3*(i-1)+2]));
                        dN++;
                    }

                    //if((1-(TphiPrime/(v->d[0])))<0.0) cout << (1-(TphiPrime/(v->d[0]))) << endl;
                }
                delete v;
            }
        }

    }


    //deletion
    delete a[2];
    delete a[1];
    delete a[0];
    delete a;

    if(dN==0)
    {
        //return E/(((double) N));
        return E/(((double) N) * len_fib(phi));
    }
    else
    {
        //cout << N << " " << dN << endl;
        if(dN>=N)
        {
            //cout << "error N" << N << endl;
            return -1.0;
        }
        else
        {
            if(dLength<len_fib(phi))
            {
                //return E/( (((double) N)-((double) dN)) );
                return E/( (((double) N)-((double) dN)) * (len_fib(phi) - dLength));
            }
            else
            {
                //cout << "returning -1.0: error length " << len_fib(phi) << " " << dLength << endl;
                //cout << len_fib(phi) << " " << dLength << endl;
                return -1.0;
            }
        }
    }
#else

    //compute the length of the fiber
    //compute Reimann integrale (approximation)
    unsigned long N=0;
    double S = 0.0, phiX, phiY, phiZ, sumDelta=0.0, phiXdev, phiYdev, phiZdev, tmpS, tmpL;
    double** D = new double*[3];
    D[0] = new double[3];
    D[1] = new double[3];
    D[2] = new double[3];

    for(unsigned long i=0 ; i<(phi->N_element-1) ; i++)
    {
        //calculate speed vector
        phiX = phi->elts[3*(i+1)]-phi->elts[3*i];
        phiY = phi->elts[3*(i+1)+1]-phi->elts[3*i+1];
        phiZ = phi->elts[3*(i+1)+2]-phi->elts[3*i+2];

        tmpL = sqrt(SQR(phiX) + SQR(phiY) + SQR(phiZ));

        ///check norm of speed vector///
        if(!(isnan(tmpL) || isinf(tmpL) || (tmpL<SMALL_NUM)))
        {

            //get tensor at point (phi->elts[3*i], phi->elts[3*i+2], phi->elts[3*i+2])
            fillArray(phi->elts[3*i], phi->elts[3*i+1], phi->elts[3*i+2], D);//must check what is in D

            //calculates tensor product
            phiXdev = D[0][0]*phiX + D[0][1]*phiY + D[0][2]*phiZ;
            phiYdev = D[1][0]*phiX + D[1][1]*phiY + D[1][2]*phiZ;
            phiZdev = D[2][0]*phiX + D[2][1]*phiY + D[2][2]*phiZ;

            //calculates eigen values
            C_toolbox_eigen_sym v((double**) D, 3, false);
            tmpS = sqrt(SQR(phiXdev) + SQR(phiYdev) + SQR(phiZdev))/v.d[0];
            if(!isnan(tmpS) && !isinf(tmpS))
            {
                //integrate
                S += tmpS;

                //sum the length
                sumDelta +=tmpL;
            }
            else
            {
                cout << "tmpS nan or inf " << tmpS << endl;
                N++;
            }
        }
        else
        {
            cout << "tmpL nan or inf " << tmpL << endl;
            N++;
        }
    }

    //deletion
    delete D[2];
    delete D[1];
    delete D[0];
    delete D;

    if(N<(phi->N_element-1)) return S/sumDelta;
    else return -1.0;

#endif

}

double C_measure::error_fiber_tilde(CtlStruct* phi)
{
    //compute the length of the fiber
    //compute Reimann integrale (approximation)
    unsigned long N=0;
    double S = 0.0, phiX, phiY, phiZ, sumDelta=0.0, phiXdev, phiYdev, phiZdev, tmpS, tmpL;
    double** D = new double*[3];
    D[0] = new double[3];
    D[1] = new double[3];
    D[2] = new double[3];

    for(unsigned long i=0 ; i<(phi->N_element-1) ; i++)
    {
        //calculate speed vector
        phiX = phi->elts[3*(i+1)]-phi->elts[3*i];
        phiY = phi->elts[3*(i+1)+1]-phi->elts[3*i+1];
        phiZ = phi->elts[3*(i+1)+2]-phi->elts[3*i+2];

        tmpL = sqrt(SQR(phiX) + SQR(phiY) + SQR(phiZ));

        ///check norm of speed vector///
        if(!(isnan(tmpL) || isinf(tmpL) || (tmpL<SMALL_NUM)))
        {

            //get tensor at point (phi->elts[3*i], phi->elts[3*i+2], phi->elts[3*i+2])
            fillArray(phi->elts[3*i], phi->elts[3*i+1], phi->elts[3*i+2], D);//must check what is in D

            //calculates tensor product
            phiXdev = D[0][0]*phiX + D[0][1]*phiY + D[0][2]*phiZ;
            phiYdev = D[1][0]*phiX + D[1][1]*phiY + D[1][2]*phiZ;
            phiZdev = D[2][0]*phiX + D[2][1]*phiY + D[2][2]*phiZ;

            //calculates eigen values
            C_toolbox_eigen_sym v((double**) D, 3, false);
            tmpS = sqrt(SQR(phiXdev) + SQR(phiYdev) + SQR(phiZdev))/v.d[0];
            if(!isnan(tmpS) && !isinf(tmpS))
            {
                //integrate
                S += tmpS;

                //sum the length
                sumDelta +=tmpL;
            }
            else
            {
                cout << "tmpS nan or inf " << tmpS << endl;
                N++;
            }
        }
        else
        {
            cout << "tmpL nan or inf " << tmpL << endl;
            N++;
        }
    }

    //deletion
    delete D[2];
    delete D[1];
    delete D[0];
    delete D;

    if(N<(phi->N_element-1)) return (sumDelta-S)/SQR(sumDelta);
    else return -1.0;

}

double C_measure::len_fib(CtlStruct* phi)
{
    double len = 0.0;
    for(unsigned long i=1 ; i<phi->N_element ; i++)
    {
        len += sqrt( SQR(phi->elts[3*i] - phi->elts[3*(i-1)]) + SQR(phi->elts[3*i+1] - phi->elts[3*(i-1)+1]) + SQR(phi->elts[3*i+2] - phi->elts[3*(i-1)+2]));
    }
    return len;
}

void C_measure::splinMeasure(string fileNameFib, double* mean, double* sigma, double* Etild2, double* E2)
{
    C_load_fiber* loadObj = new C_load_fiber();
    CtlCtlStruct* FIBERS = loadObj->readFiber(fileNameFib, false); ///interpolated fibers
    if(FIBERS==NULL)
    {
        return;
    }
    splinMeasure(FIBERS, mean, sigma, Etild2, E2);
    delete FIBERS;
    delete loadObj;
    return;
}
void C_measure::splinMeasure(CtlCtlStruct* FIBERS, double* mean, double* sigma, double* Etild2, double* E2)
{
    if(mean==NULL && sigma==NULL && Etild2==NULL && E2==NULL)
    {
        return;
    }

    if(mean!=NULL && sigma!=NULL)
    {
        *mean = meanLengthSplin(FIBERS);
        *sigma = stdLengthSplin(FIBERS, *mean);
    }
    else
    {
        if(mean!=NULL)
        {
            *mean = meanLengthSplin(FIBERS);
        }
        if(sigma!=NULL)
        {
            *sigma = stdLengthSplin(FIBERS, meanLengthSplin(FIBERS));
        }
    }

    if(Etild2!=NULL && E2!=NULL)
    {
        *Etild2 = meanLengthErrorSplinE2(FIBERS, E2);
    }
    else
    {
        if(Etild2!=NULL)
        {
            double dummy;
            *Etild2 = meanLengthErrorSplinE2(FIBERS, &dummy);
        }
        if(E2!=NULL)
        {
            meanLengthErrorSplinE2(FIBERS, E2);
        }
    }

    return;
}



long C_measure::getFibersStat(string fileNameFib, double* L, double* S)
{
    C_load_fiber* loadObj = new C_load_fiber();
    CtlCtlStruct* FIB = loadObj->readFiber(fileNameFib, true);
    if(FIB==NULL)
    {
        *L = 0;
        *S = 0;
        return 0;
    }
    *L = meanLengthSplin(FIB);
    *S = stdLengthSplin(FIB, *L);
    long N = FIB->N;
    delete FIB;
    delete loadObj;
    return N;
}
void C_measure::getLengthAndNormalizedTotalPerimeter(CtlCtlStruct* FIB, vector<double>* L, vector<double>* totPermim)
{
    double mL, mtotPermim;
    for(unsigned long i=0 ; i<FIB->N ; i++) //3->4
    {
        if(FIB->fibers[i].N_element>2)
        {
            getLengthAndNormalizedTotalPerimeter(&(FIB->fibers[i]), &mL, &mtotPermim);
            L->push_back(mL);
            totPermim->push_back(mtotPermim);
        }
        else if(FIB->fibers[i].N_element==2)
        {
            L->push_back(sqrt(SQR(FIB->fibers[i].elts[3]-FIB->fibers[i].elts[0]) + SQR(FIB->fibers[i].elts[4]-FIB->fibers[i].elts[1]) + SQR(FIB->fibers[i].elts[5]-FIB->fibers[i].elts[2])));
            totPermim->push_back(0.0);
        }
        else if(FIB->fibers[i].N_element==1)
        {
            L->push_back(0.0);
            totPermim->push_back(0.0);
        }
    }
    return;
}

void C_measure::getLengthAndNormalizedTotalPerimeter(CtlStruct* phi, double* L, double* totPerim)
{
    *L = LSplin(phi);///this is actualy not spline interpolation, but if phi is already spline interpolation


    double mL, mP, x1, x2, x3, y1, y2, y3, z1, z2, z3, x12, x13, y12, y13, z12, z13, nx, ny, nz, xp, yp, zp;
    int N;
//    for(int n=0 ; n<FIBERS->N ; n++)
//    {
        //cout << "fiber: " << n << " over " << FIBERS->N << endl;
        //generation des facettes de la n-eme fibre
        N = phi->N_element;
        //cout << N << " points in fiber #" << n << endl;
        vector<facette> F;
        for(int i=0 ; i<(N-2) ; i++)
        {
            for(int j=i+1 ; j<(N-1) ; j++)
            {
                for(int k=j+1 ; k<N ; k++)
                {
                    //verifie si la facette n'hexiste pas deja
//                    bool exist = false;
//                    for(int q=0 ; q< F.size() ; q++)
//                    {
//                        if(F.at(q).p1==i && F.at(q).p2==j && F.at(q).p3==k) exist=true;
//                        if(F.at(q).p2==i && F.at(q).p1==j && F.at(q).p3==k) exist=true;
//                        if(F.at(q).p1==i && F.at(q).p3==j && F.at(q).p2==k) exist=true;
//                        if(F.at(q).p2==i && F.at(q).p3==j && F.at(q).p1==k) exist=true;
//                        if(F.at(q).p3==i && F.at(q).p1==j && F.at(q).p2==k) exist=true;
//                        if(F.at(q).p3==i && F.at(q).p2==j && F.at(q).p1==k) exist=true;
//                        if(exist) q=F.size();
//                    }
//                    if(!exist)
//                    {
                        facette f(i,j,k);
                        F.push_back(f);
//                    }
                }
            }
        }

        //cout << "nombre de facttes : " << F.size() << endl;

        //extraction des facettes visible de l'exterieur
        vector<facette> Hull;
        bool currentSign, previousSign, changeSign, first;
        for(int i=0 ; i<F.size() ; i++)
        {
            //point de la facette
            x1 = phi->elts[3*F.at(i).p1];
            x2 = phi->elts[3*F.at(i).p2];
            x3 = phi->elts[3*F.at(i).p3];

            y1 = phi->elts[3*F.at(i).p1+1];
            y2 = phi->elts[3*F.at(i).p2+1];
            y3 = phi->elts[3*F.at(i).p3+1];

            z1 = phi->elts[3*F.at(i).p1+2];
            z2 = phi->elts[3*F.at(i).p2+2];
            z3 = phi->elts[3*F.at(i).p3+2];

            //vecteur generateur de la facette
            x12 = x2-x1;
            x13 = x3-x1;

            y12 = y2-y1;
            y13 = y3-y1;

            z12 = z2-z1;
            z13 = z3-z1;

            //vecteur normale au plan
            nx = y12*z13 - z12*y13;
            ny = z12*x13 - x12*z13;
            nz = x12*y13 - y12*x13;

            if( (SQR(nx)+SQR(ny)+SQR(nz))>0.0000001 )
            {
                //pour chaque point de la fibre
                first=true;
                changeSign = false;
                for(int j=0 ; j<phi->N_element ; j++)
                {
                    //si le point courant n'est pas sur la facette
                    if( (j!=F.at(i).p1) && (j!=F.at(i).p2) && (j!=F.at(i).p3))
                    {
                        //changement de repere : translation a l'origine de la facette (au point 1)
                        xp = phi->elts[3*j]-x1;
                        yp = phi->elts[3*j+1]-y1;
                        zp = phi->elts[3*j+2]-z1;

                        //de quel cote de la facette F.at(i) est le point j de FIBERS->fibers[n]
                        if(first)
                        {
                            //initialisation : pas de comparaison a faire
                            first=false;
                            currentSign = (xp*nx + yp*ny + zp*nz)>0.0;
                            //cout << "first" << endl;
                        }
                        else
                        {
                            //determine de quel cote se trouve le point et compare le cote avec le point precedent
                            previousSign = currentSign;
                            currentSign = (xp*nx + yp*ny + zp*nz)>0.0;
                            if(previousSign!=currentSign)
                            {
                                //si le point est de l'autre cote de la facette, la facette n'est pas visible de l'exterieur
                                changeSign = true;
                            }
                        }
                        if(changeSign)
                        {
                            //on peut donc arreter la boucle car la facette ne sera pas retenue pour la formation de l'enveloppe convex
                            j=phi->N_element;
                        }
                    }
                }

                //si la facette est visible de l'exterieure, on garde la facette dans l'enveloppe convexe.
                if(!changeSign) Hull.push_back(F.at(i));
            }


        }

        //calcul du perimetre total
        mP=0.0;
        for(int i=0 ; i<Hull.size() ; i++)
        {
            //point de la facette
            x1 = phi->elts[3*Hull.at(i).p1];
            x2 = phi->elts[3*Hull.at(i).p2];
            x3 = phi->elts[3*Hull.at(i).p3];

            y1 = phi->elts[3*Hull.at(i).p1+1];
            y2 = phi->elts[3*Hull.at(i).p2+1];
            y3 = phi->elts[3*Hull.at(i).p3+1];

            z1 = phi->elts[3*Hull.at(i).p1+2];
            z2 = phi->elts[3*Hull.at(i).p2+2];
            z3 = phi->elts[3*Hull.at(i).p3+2];

//            cout << "---------------------------------" << endl;
//            cout << "facette " << i << " p1 : (" << x1 << "," << y1 << "," << z1 << ")" << endl;
//            cout << "facette " << i << " p1 : (" << x2 << "," << y2 << "," << z2 << ")" << endl;
//            cout << "facette " << i << " p1 : (" << x3 << "," << y3 << "," << z3 << ")" << endl;

            //calcul du perimetre de la facette
            mP += sqrt(SQR(x1-x2) + SQR(y1-y2) + SQR(z1-z2));
            mP += sqrt(SQR(x2-x3) + SQR(y2-y3) + SQR(z2-z3));
            mP += sqrt(SQR(x1-x3) + SQR(y1-y3) + SQR(z1-z3));
        }
//        cout << "---------------------------------" << endl;


        //calcul la longueur de la fibre
        mL=0.0;
        for(int i=1 ; i<N ; i++)
        {
            mL += sqrt(SQR(phi->elts[3*i]-phi->elts[3*(i-1)]) +\
                       SQR(phi->elts[3*i+1]-phi->elts[3*(i-1)+1]) +\
                       SQR(phi->elts[3*i+2]-phi->elts[3*(i-1)+2]) );
        }
        *L = mL;

        //(*facetteNum).push_back(Hull.size());



        mP *=0.5/((double) Hull.size()); //chaque aretes des facettes a ete comptee deux fois
        *totPerim = mP;
//    }

    return;
}

void C_measure::getLengthAndBoxDimension(CtlCtlStruct* FIB, vector<double>* L, vector<double>* dU, vector<double>* dV, vector<double>* dW)
{
    double mL, mdU, mdV, mdW;
    for(unsigned long i=0 ; i<FIB->N ; i++) //3->4
    {
        if(FIB->fibers[i].N_element>2)
        {
            getLengthAndBoxDimension(&(FIB->fibers[i]), &mL, &mdU, &mdV, &mdW);
            L->push_back(mL);
            dU->push_back(mdU);
            dV->push_back(mdV);
            dW->push_back(mdW);
        }
        else if(FIB->fibers[i].N_element==2)
        {
            L->push_back(sqrt(SQR(FIB->fibers[i].elts[3]-FIB->fibers[i].elts[0]) + SQR(FIB->fibers[i].elts[4]-FIB->fibers[i].elts[1]) + SQR(FIB->fibers[i].elts[5]-FIB->fibers[i].elts[2])));
            dU->push_back(sqrt(SQR(FIB->fibers[i].elts[3]-FIB->fibers[i].elts[0]) + SQR(FIB->fibers[i].elts[4]-FIB->fibers[i].elts[1]) + SQR(FIB->fibers[i].elts[5]-FIB->fibers[i].elts[2])));
            dV->push_back(0.0);
            dW->push_back(0.0);
        }
        else if(FIB->fibers[i].N_element==1)
        {
            L->push_back(0.0);
            dU->push_back(0.0);
            dV->push_back(0.0);
            dW->push_back(0.0);
        }
    }
    return;
}
void C_measure::getLengthAndBoxDimension(string fileNameFib, vector<double>* L, vector<double>* dU, vector<double>* dV, vector<double>* dW, bool rawCond)
{
    C_load_fiber* loadObj = new C_load_fiber();
    CtlCtlStruct* FIB = loadObj->readFiber(fileNameFib, rawCond);
    double mL, mdU, mdV, mdW;
    for(unsigned long i=0 ; i<FIB->N ; i++) //3->4
    {
        if(FIB->fibers[i].N_element>2)
        {
            getLengthAndBoxDimension(&(FIB->fibers[i]), &mL, &mdU, &mdV, &mdW);
            L->push_back(mL);
            dU->push_back(mdU);
            dV->push_back(mdV);
            dW->push_back(mdW);
        }
        else if(FIB->fibers[i].N_element==2)
        {
            L->push_back(sqrt(SQR(FIB->fibers[i].elts[3]-FIB->fibers[i].elts[0]) + SQR(FIB->fibers[i].elts[4]-FIB->fibers[i].elts[1]) + SQR(FIB->fibers[i].elts[5]-FIB->fibers[i].elts[2])));
            dU->push_back(sqrt(SQR(FIB->fibers[i].elts[3]-FIB->fibers[i].elts[0]) + SQR(FIB->fibers[i].elts[4]-FIB->fibers[i].elts[1]) + SQR(FIB->fibers[i].elts[5]-FIB->fibers[i].elts[2])));
            dV->push_back(0.0);
            dW->push_back(0.0);
        }
        else if(FIB->fibers[i].N_element==1)
        {
            L->push_back(0.0);
            dU->push_back(0.0);
            dV->push_back(0.0);
            dW->push_back(0.0);
        }
    }
    delete FIB;
    delete loadObj;
    return;
}
void C_measure::getLengthAndBoxDimension(CtlStruct* phi, double* L, double* dU, double* dV, double* dW)
{
    *L = LSplin(phi);///this is actualy not spline interpolation, but if phi is already spline interpolation
    ///center the point distribution
    double** X = new double*[phi->N_element];
    for(unsigned long i=0 ; i<phi->N_element ; i++)
    {
        X[i] = new double[3];
        X[i][0] = (double) phi->elts[3*i+0];
        X[i][1] = (double) phi->elts[3*i+1];
        X[i][2] = (double) phi->elts[3*i+2];
    }
    //center data
    double* Xmu = new double[3];
    Xmu[0] = 0.0;
    Xmu[1] = 0.0;
    Xmu[2] = 0.0;
    for(unsigned long i=0 ; i<phi->N_element ; i++)
    {
        Xmu[0] += X[i][0];
        Xmu[1] += X[i][1];
        Xmu[2] += X[i][2];
    }
    Xmu[0] /= (double) phi->N_element;
    Xmu[1] /= (double) phi->N_element;
    Xmu[2] /= (double) phi->N_element;
    for(unsigned long i=0 ; i<phi->N_element ; i++)
    {
        X[i][0] -= Xmu[0];
        X[i][1] -= Xmu[1];
        X[i][2] -= Xmu[2];
    }
    delete Xmu;

    ///compute PCA
    C_toolbox_PCA_3D* toolPCA = new C_toolbox_PCA_3D();
    double** W = toolPCA->getPrincipalComponents3D(X, phi->N_element);
    delete toolPCA;

    ///rotate the data in frame W=(w1,w2,w3)
    double* Xnew = new double[3];
    for(unsigned long i=0 ; i<phi->N_element ; i++)
    {
        Xnew[0] = W[0][0]*X[i][0] + W[1][0]*X[i][1] + W[2][0]*X[i][2];
        Xnew[1] = W[0][1]*X[i][0] + W[1][1]*X[i][1] + W[2][1]*X[i][2];
        Xnew[2] = W[0][2]*X[i][0] + W[1][2]*X[i][1] + W[2][2]*X[i][2];
        X[i][0] = Xnew[0];
        X[i][1] = Xnew[1];
        X[i][2] = Xnew[2];
    }
    delete Xnew;
    for(unsigned long i=0 ; i<3 ; i++)
    {
        delete W[i];
    }
    delete W;

    ///find min and max of all components in new frame
    double xmin=X[0][0], xmax=X[0][0];
    double ymin=X[0][1], ymax=X[0][1];
    double zmin=X[0][2], zmax=X[0][2];
    for(unsigned long i=1 ; i<phi->N_element ; i++)
    {
        if(xmin>X[i][0]) xmin = X[i][0];
        if(xmax<X[i][0]) xmax = X[i][0];

        if(ymin>X[i][1]) ymin = X[i][1];
        if(ymax<X[i][1]) ymax = X[i][1];

        if(zmin>X[i][2]) zmin = X[i][2];
        if(zmax<X[i][2]) zmax = X[i][2];
    }

    *dU = xmax-xmin;
    *dV = ymax-ymin;
    *dW = zmax-zmin;

    ///free
    for(unsigned long i=0 ; i<phi->N_element ; i++)
    {
        delete X[i];
    }
    delete X;
    return;
}

void C_measure::getHelixAndTransversAngle(CtlCtlStruct* FIBERS, double xc, double yc, double zc, double dx, double dy, double dz,  vector<double>* helixAngle, vector<double>* transversAngle, vector<double>* r, vector<double>* theta, vector<double>* z)
{
    double X, Y, Z, U, V, W;
    //create rotation matrix
    double** R = new double*[3];
    R[0] = new double[3]; R[1] = new double[3]; R[2] = new double[3];

    R[2][0] = dx/sqrt(SQR(dx)+SQR(dy)+SQR(dz));        R[2][1] = dy/sqrt(SQR(dx)+SQR(dy)+SQR(dz));         R[2][2] = dz/sqrt(SQR(dx)+SQR(dy)+SQR(dz));
    R[0][0] = R[2][1]/sqrt(SQR(R[2][0])+SQR(R[2][1])); R[0][1] = -R[2][0]/sqrt(SQR(R[2][0])+SQR(R[2][1])); R[0][2] = 0.0;
    R[1][0] = R[2][1]*R[0][2] - R[2][2]*R[0][1];       R[1][1] = R[2][2]*R[0][0] - R[2][0]*R[0][2];        R[1][2] = R[2][0]*R[0][1] - R[2][1]*R[0][0];


    for(unsigned long i=0 ; i<FIBERS->N ; i++)
    {
        for(unsigned long j=1 ; j<FIBERS->fibers[i].N_element-1 ; j++)
        {
            X = FIBERS->fibers[i].elts[3*j+0]-xc;
            Y = FIBERS->fibers[i].elts[3*j+1]-yc;
            Z = FIBERS->fibers[i].elts[3*j+2]-zc;

            U = R[0][0]*X + R[0][1]*Y + R[0][2]*Z;
            V = R[1][0]*X + R[1][1]*Y + R[1][2]*Z;
            W = R[2][0]*X + R[2][1]*Y + R[2][2]*Z;

            r->push_back(sqrt(SQR(U)+SQR(V)));
            theta->push_back(atan2(V,U));
            z->push_back(W);

            helixAngle->push_back(getHelixAngle(dx, dy, dz, FIBERS->fibers[i].elts[3*(j+1)+0]-FIBERS->fibers[i].elts[3*(j-1)+0], FIBERS->fibers[i].elts[3*(j+1)+1]-FIBERS->fibers[i].elts[3*(j-1)+1], FIBERS->fibers[i].elts[3*(j+1)+2]-FIBERS->fibers[i].elts[3*(j-1)+2]));
            transversAngle->push_back(getTransversAngle(dx, dy, dz, X, Y, Z, FIBERS->fibers[i].elts[3*(j+1)+0]-FIBERS->fibers[i].elts[3*(j-1)+0], FIBERS->fibers[i].elts[3*(j+1)+1]-FIBERS->fibers[i].elts[3*(j-1)+1], FIBERS->fibers[i].elts[3*(j+1)+2]-FIBERS->fibers[i].elts[3*(j-1)+2]));
        }
    }
    delete R[0];
    delete R[1];
    delete R[2];
    delete R;
    return;
}

//
double C_measure::getHelixAngle(double dx, double dy, double dz, double dfx, double dfy, double dfz)
{
    return asin( (dx*dfx+dy*dfy+dz*dfz)/sqrt( (SQR(dx)+SQR(dy)+SQR(dz))*(SQR(dfx)+SQR(dfy)+SQR(dfz)) ) );
}
double C_measure::getTransversAngle(double dx, double dy, double dz, double dux, double duy, double duz, double dfx, double dfy, double dfz)
{
    double fxOrth = dfx - (dx*dfx+dy*dfy+dz*dfz)*dx/(SQR(dx)+SQR(dy)+SQR(dz));
    double fyOrth = dfy - (dx*dfx+dy*dfy+dz*dfz)*dy/(SQR(dx)+SQR(dy)+SQR(dz));
    double fzOrth = dfz - (dx*dfx+dy*dfy+dz*dfz)*dz/(SQR(dx)+SQR(dy)+SQR(dz));
    return asin((fxOrth*dux + fyOrth*duy + fzOrth*duz)/sqrt( (SQR(fxOrth)+SQR(fyOrth)+SQR(fzOrth))*(SQR(dux)+SQR(duy)+SQR(duz)) ));
}
