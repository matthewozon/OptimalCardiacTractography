#include <C_synthetic_DWI_patho.h>

C_synthetic_DWI_patho::C_synthetic_DWI_patho(double X, double Y, double Z, double xc, double yc, double zc, double r):C_synthetic_DWI()
{
    //ctor
    m_X = X;
    m_Y = Y;
    m_Z = Z;
    m_xc = xc;
    m_yc = yc;
    m_zc = zc;
    m_r = r;
    dimX = (unsigned long) X+1;
    dimY = (unsigned long) Y+1;
    dimZ = (unsigned long) Z+1;
}

C_synthetic_DWI_patho::~C_synthetic_DWI_patho()
{
    //dtor
}


rawData<double>* C_synthetic_DWI_patho::createB0Data(diffInfo* g) //, double xmin, double xmax, double ymin, double ymax, double zmin, double zmax, double rmin, double rmax
{

    rawData<double>* rawDataB0 = new rawData<double>(DICOM_DOUBLE, 3, dimX, dimY, dimZ);
    rawDataB0->pixDimX = 1.0;
    rawDataB0->pixDimY = 1.0;
    rawDataB0->pixDimZ = 1.0;
    rawDataB0->pixDimT = 1.0;
    (*rawDataB0) = 1.0;
    return rawDataB0;
}

rawData<double>* C_synthetic_DWI_patho::createMask(diffInfo* g) //, double xmin, double xmax, double ymin, double ymax, double zmin, double zmax, double rmin, double rmax
{
    rawData<double>* mask = createB0Data(g);
    *mask =  (*mask)>0.5;
    return mask;
}

vector< vector<rawData<double>*> > C_synthetic_DWI_patho::createDWIdata( vector< vector<diffInfo*> > dfi, vector<rawData<double>*> S0)
{
    vector<vector<rawData<double>*> > allRawData( dfi.size(), vector<rawData<double>*>(dfi.at(0).size(),NULL) );
    for(unsigned int i=0 ; i<dfi.size() ; i++)
    {
        for(unsigned int j=0 ; j<dfi.at(i).size() ; j++)
        {
            allRawData.at(i).at(j) = createOneDWI(dfi.at(i).at(j), S0.at(j));
        }
    }
    return allRawData;
}


rawData<double>* C_synthetic_DWI_patho::createOneDWI(diffInfo* g, rawData<double>* S0)
{
    rawData<double>* rawDataDWI = new rawData<double>(DICOM_DOUBLE, 3, dimX, dimY, dimZ);
    rawDataDWI->pixDimX = 1.0;
    rawDataDWI->pixDimY = 1.0;
    rawDataDWI->pixDimZ = 1.0;
    rawDataDWI->pixDimT = 1.0;
    (*rawDataDWI) = 0.0;

    double V1[3];
    double V2[3];
    double V3[3];
    double* D = new double[6];
    double vp1, vp2, vp3;

    ///for each pixel in rawDataDWI, set the 3 eigen vectors and compute DWI data
    double distSQRX, distSQRY, distSQRZ, SQRR = SQR(m_r);
    for(unsigned long i=0 ; i<dimX ; i++)
    {
        distSQRX = SQR(((double) i*g->voxelSize[0])-m_xc);
        for(unsigned long j=0 ; j<dimY ; j++)
        {
            distSQRY = SQR(((double) j*g->voxelSize[1])-m_yc);
            for(unsigned long k=0 ; k<dimZ ; k++)
            {
                distSQRZ = SQR(((double) k*g->voxelSize[2])-m_zc);
                if( (distSQRX+distSQRY+distSQRZ)<=SQRR)
                {
                    //draw a random v1
                    V1[0] = toolRand->doub()-0.5;
                    V1[1] = toolRand->doub()-0.5;
                    V1[2] = toolRand->doub()-0.5;
                    while(!(V1[0]!=0.0 && V1[1]!=0.0))
                    {
                        V1[0] = toolRand->doub()-0.5;
                        V1[1] = toolRand->doub()-0.5;
                        V1[2] = toolRand->doub()-0.5;
                    }

                    //calculate v2
                    V2[0] = -V1[1];
                    V2[1] = V1[0];
                    V2[2] = 0.0;

                    //calculate v3
                    V3[0] = V1[1]*V2[2] - V1[2]*V2[1];
                    V3[1] = V1[2]*V2[0] - V1[0]*V2[2];
                    V3[2] = V1[0]*V2[1] - V1[1]*V2[0];

                    //draw the eigen values
                    vp1 = 0.05+0.5*toolRand->doub();
                    vp2 = 0.05+0.5*toolRand->doub();
                    vp3 = 0.05+0.5*toolRand->doub();
                }
                else
                {
                    V1[0] = 1.0;
                    V1[1] = 0.0;
                    V1[2] = 0.0;

                    V2[0] = 0.0;
                    V2[1] = 1.0;
                    V2[2] = 0.0;

                    //calculate v3
                    V3[0] = 0.0;
                    V3[1] = 0.0;
                    V3[2] = 1.0;

                    vp1 = 1.0;
                    vp2 = 0.2;
                    vp3 = 0.2;
                }

                    ///create tensor from eigen vectors and eigen values
                    D[0] = vp1*SQR(V1[0])  + vp2*SQR(V2[0])  + vp3*SQR(V3[0]);
                    D[1] = vp1*V1[0]*V1[1] + vp2*V2[0]*V2[1] + vp3*V3[0]*V3[1];
                    D[2] = vp1*V1[0]*V1[2] + vp2*V2[0]*V2[2] + vp3*V3[0]*V3[2];
                    D[3] = vp1*SQR(V1[1])  + vp2*SQR(V2[1])  + vp3*SQR(V3[1]);
                    D[4] = vp1*V1[1]*V1[2] + vp2*V2[1]*V2[2] + vp3*V3[1]*V3[2];
                    D[5] = vp1*SQR(V1[2])  + vp2*SQR(V2[2])  + vp3*SQR(V3[2]);

                    //calculate D*g and store it in V1
                    V1[0] = D[0]*g->g[0] + D[1]*g->g[1] + D[2]*g->g[2];
                    V1[1] = D[1]*g->g[0] + D[3]*g->g[1] + D[4]*g->g[2];
                    V1[2] = D[2]*g->g[0] + D[4]*g->g[1] + D[5]*g->g[2];

                    //calculate g'*D*g
                    rawDataDWI->raw3D[i][j][k] =  g->g[0]*V1[0] + g->g[1]*V1[1] + g->g[2]*V1[2];

                    //calculate diffusion signal
                    rawDataDWI->raw3D[i][j][k] = S0->raw3D[i][j][k]*exp(-g->b*rawDataDWI->raw3D[i][j][k]);
            }
        }
    }
    delete D;
    return rawDataDWI;
}
