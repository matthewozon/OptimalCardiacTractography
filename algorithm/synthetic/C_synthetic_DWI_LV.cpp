#include <C_synthetic_DWI_LV.h>

C_synthetic_DWI_LV::C_synthetic_DWI_LV(double alphaMin, double alphaMax, double H, double a, double b, double A, double B):C_synthetic_DWI()
{
    //ctor
    m_alphamin = alphaMin;
    m_alphamax = alphaMax;
    m_H = H;
    m_a = a;
    m_b = b;
    m_A = A;
    m_B = B;
    m_xmin = 0.0; m_xmax = 2.0*a+4.0;
    m_ymin = 0.0; m_ymax = 2.0*a+4.0;
    m_zmin = -1.0; m_zmax = m_H+1.0;
}

C_synthetic_DWI_LV::~C_synthetic_DWI_LV()
{
    //dtor
}


rawData<double>* C_synthetic_DWI_LV::createB0Data(diffInfo* g)
{
    double Xc = 0.5*(m_xmin+m_xmax);
    double Yc = 0.5*(m_ymin+m_ymax);
    double r2;
    unsigned long dimX = 0, dimY = 0, dimZ = 0;
    for(double xx=m_xmin ; xx<=m_xmax ; xx+=g->voxelSize[0])
    {
        dimX++;
    }
    for(double yy=m_ymin ; yy<=m_ymax ; yy+=g->voxelSize[1])
    {
        dimY++;
    }
    for(double zz=m_zmin ; zz<=m_zmax ; zz+=g->voxelSize[2])
    {
        dimZ++;
    }
    rawData<double>* rawDataB0 = new rawData<double>(DICOM_DOUBLE, 3, dimX, dimY, dimZ);
    rawDataB0->pixDimX = 1.0;
    rawDataB0->pixDimY = 1.0;
    rawDataB0->pixDimZ = 1.0;
    rawDataB0->pixDimT = 1.0;
    (*rawDataB0) = 0.0;
    unsigned long i=0, j=0, k=0;

    ///for each pixel in rawDataDWI, set the 3 eigen vectors and compute DWI data
    i=0;
    for(double x=m_xmin ; x<=m_xmax ; x+=g->voxelSize[0], i++)
    {
        j=0;
        for(double y=m_ymin ; y<=m_ymax ; y+=g->voxelSize[1], j++)
        {
            r2 = SQR(x-Xc) + SQR(y-Yc);
            k=0;
            for(double z=m_zmin ; z<=m_zmax ; z+=g->voxelSize[2], k++)
            {
                if( ((r2/SQR(m_a))+SQR(z/m_b))<=1.0 &&  ((r2/SQR(m_A))+SQR(z/m_B))>=1.0 )
                {
                    rawDataB0->raw3D[i][j][k] = 1.0;
                }
            }
        }
    }
    return rawDataB0;
}

rawData<double>* C_synthetic_DWI_LV::createMask(diffInfo* g) //, double xmin, double xmax, double ymin, double ymax, double zmin, double zmax, double rmin, double rmax
{
    rawData<double>* mask = createB0Data(g);
    *mask =  (*mask)>0.5;
    return mask;
}

vector< vector<rawData<double>*> > C_synthetic_DWI_LV::createDWIdata( vector< vector<diffInfo*> > dfi, vector<rawData<double>*> S0)
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


rawData<double>* C_synthetic_DWI_LV::createOneDWI(diffInfo* g, rawData<double>* S0)
{
    double Xc = 0.5*(m_xmin+m_xmax);
    double Yc = 0.5*(m_ymin+m_ymax);
    double r2;
    unsigned long dimX = 0, dimY = 0, dimZ = 0;
    for(double xx=m_xmin ; xx<=m_xmax ; xx+=g->voxelSize[0])
    {
        dimX++;
    }
    for(double yy=m_ymin ; yy<=m_ymax ; yy+=g->voxelSize[1])
    {
        dimY++;
    }
    for(double zz=m_zmin ; zz<=m_zmax ; zz+=g->voxelSize[2])
    {
        dimZ++;
    }
    rawData<double>* rawDataDWI = new rawData<double>(DICOM_DOUBLE, 3, dimX, dimY, dimZ);
    rawDataDWI->pixDimX = 1.0;
    rawDataDWI->pixDimY = 1.0;
    rawDataDWI->pixDimZ = 1.0;
    rawDataDWI->pixDimT = 1.0;
    (*rawDataDWI) = 0.0;
    double dAlpha = m_alphamax - m_alphamin;
    double alpha;
    double dr;
    double rmin, rmax;
    double Z;

    double V1[3];// = {1.0, 0.0, 0.0};
    double V2[3];// = {0.0, 1.0, 0.0};
    double V3[3];// = {0.0, 0.0, 1.0};
    double* D = new double[6];
    double vp1 = 1.0, vp2 = 0.2, vp3 = 0.2;
    unsigned long i=0, j=0, k=0;

    ///for each pixel in rawDataDWI, set the 3 eigen vectors and compute DWI data
    i=0;
    for(double x=m_xmin ; x<=m_xmax ; x+=g->voxelSize[0], i++)
    {
        j=0;
        for(double y=m_ymin ; y<=m_ymax ; y+=g->voxelSize[1], j++)
        {
            r2 = SQR(x-Xc) + SQR(y-Yc);
            k=0;
            for(double z=m_zmin ; z<=m_zmax ; z+=g->voxelSize[2], k++)
            {
                //if we are in the ROI
                if( ((r2/SQR(m_a))+SQR(z/m_b))<=1.0 &&  ((r2/SQR(m_A))+SQR(z/m_B))>=1.0 )
                {
                    if(m_B>z) rmin = m_A*sqrt(1.0-SQR(z/m_B));
                    else rmin = 0.0;
                    rmax = m_a*sqrt(1.0-SQR(z/m_b));
                    dr = rmax-rmin;
                    alpha = dAlpha*(sqrt(r2)-rmin)/dr + m_alphamin;

                    //if r2==0 ==> isotropic tensor
                    if(r2<0.0001)
                    {
                        ///create isotropic
                        D[0] = vp3;
                        D[1] = 0.0;
                        D[2] = 0.0;
                        D[3] = vp3;
                        D[4] = 0.0;
                        D[5] = vp3;

                        //calculate D*g and store it in V1
                        V1[0] = D[0]*g->g[0] + D[1]*g->g[1] + D[2]*g->g[2];
                        V1[1] = D[1]*g->g[0] + D[3]*g->g[1] + D[4]*g->g[2];
                        V1[2] = D[2]*g->g[0] + D[4]*g->g[1] + D[5]*g->g[2];

                        //calculate g'*D*g
                        rawDataDWI->raw3D[i][j][k] =  g->g[0]*V1[0] + g->g[1]*V1[1] + g->g[2]*V1[2];

                        //calculate diffusion signal
                        rawDataDWI->raw3D[i][j][k] = S0->raw3D[i][j][k]*exp(-g->b*rawDataDWI->raw3D[i][j][k]);
                    }
                    else
                    {
                        //if angle to close to pi/2 %pi ==> special case
                        if(ABS(alpha-(Pi/2.0))<0.000001 || ABS(alpha+(Pi/2.0))<0.000001)
                        {
                            V1[0] = 0.0;
                            V1[1] = 0.0;
                            V1[2] = 1.0;

                            double r = sqrt(r2);
                            Z = sqrt( SQR((y-Yc)/r) + SQR(-(x-Xc)/r) );
                            V2[0] = (-(x-Xc)/r)/Z;
                            V2[1] = (-(y-Yc)/r)/Z;
                            V2[2] = 0.0;

                            V3[0] = V1[1]*V2[2] - V1[2]*V2[1];
                            V3[1] = V1[2]*V2[0] - V1[0]*V2[2];
                            V3[2] = V1[0]*V2[1] - V1[1]*V2[0];

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
                        else//usual case
                        {
                            double r = sqrt(r2);
                            Z = sqrt( SQR((y-Yc)/r) + SQR(-(x-Xc)/r) + SQR(alpha) );
                            V1[0] = ((y-Yc)/r)/Z;
                            V1[1] = (-(x-Xc)/r)/Z;
                            V1[2] = alpha/Z;

                            Z = sqrt( SQR((y-Yc)/r) + SQR(-(x-Xc)/r) );
                            V2[0] = (-(x-Xc)/r)/Z;
                            V2[1] = (-(y-Yc)/r)/Z;
                            V2[2] = 0.0;

                            V3[0] = V1[1]*V2[2] - V1[2]*V2[1];
                            V3[1] = V1[2]*V2[0] - V1[0]*V2[2];
                            V3[2] = V1[0]*V2[1] - V1[1]*V2[0];

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
            }
        }
    }
    delete D;
    return rawDataDWI;
}
