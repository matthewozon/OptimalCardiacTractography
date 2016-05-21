#include <C_synthetic_DWI_cylinder.h>

C_synthetic_DWI_cylinder::C_synthetic_DWI_cylinder(double xmin_, double xmax_, double ymin_, double ymax_, double zmin_, double zmax_, double rmin_, double rmax_, double alphamin_, double alphamax_):C_synthetic_DWI()
{
    //ctor
    xmin = xmin_;
    xmax = xmax_;
    ymin = ymin_;
    ymax = ymax_;
    zmin = zmin_;
    zmax = zmax_;
    rmin = rmin_;
    rmax = rmax_;
    alphamin = alphamin_;
    alphamax = alphamax_;
}

C_synthetic_DWI_cylinder::~C_synthetic_DWI_cylinder()
{
    //dtor
}


rawData<double>* C_synthetic_DWI_cylinder::createB0Data(diffInfo* g) //, double xmin, double xmax, double ymin, double ymax, double zmin, double zmax, double rmin, double rmax
{
    double Xc = 0.5*(xmin+xmax);
    double Yc = 0.5*(ymin+ymax);
    double r;
    unsigned long dimX = 0, dimY = 0, dimZ = 0;
    for(double xx=xmin ; xx<=xmax ; xx+=g->voxelSize[0])
    {
        dimX++;
    }
    for(double yy=ymin ; yy<=ymax ; yy+=g->voxelSize[1])
    {
        dimY++;
    }
    for(double zz=zmin ; zz<=zmax ; zz+=g->voxelSize[2])
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
    for(double x=xmin ; x<=xmax ; x+=g->voxelSize[0], i++)
    {
        j=0;
        for(double y=ymin ; y<=ymax ; y+=g->voxelSize[1], j++)
        {
            r = sqrt(SQR(x-Xc) + SQR(y-Yc));
            if(r<rmax && r>rmin)
            {
                k=0;
                for(double z=zmin ; z<=zmax ; z+=g->voxelSize[2], k++)
                {
                    //calculate diffusion signal
                    rawDataB0->raw3D[i][j][k] = 1.0;
                }
            }
        }
    }
    return rawDataB0;
}

rawData<double>* C_synthetic_DWI_cylinder::createMask(diffInfo* g) //, double xmin, double xmax, double ymin, double ymax, double zmin, double zmax, double rmin, double rmax
{
    rawData<double>* mask = createB0Data(g);
    *mask =  (*mask)>0.5;
//    mask->pixDimX = 1.0;
//    mask->pixDimY = 1.0;
//    mask->pixDimZ = 1.0;
//    mask->pixDimT = 1.0;
    return mask;
}

vector< vector<rawData<double>*> > C_synthetic_DWI_cylinder::createDWIdata( vector< vector<diffInfo*> > dfi, vector<rawData<double>*> S0)
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


rawData<double>* C_synthetic_DWI_cylinder::createOneDWI(diffInfo* g, rawData<double>* S0)
{
    double Xc = 0.5*(xmin+xmax);
    double Yc = 0.5*(ymin+ymax);
    double r;
    unsigned long dimX = 0, dimY = 0, dimZ = 0;
    for(double xx=xmin ; xx<=xmax ; xx+=g->voxelSize[0])
    {
        dimX++;
    }
    for(double yy=ymin ; yy<=ymax ; yy+=g->voxelSize[1])
    {
        dimY++;
    }
    for(double zz=zmin ; zz<=zmax ; zz+=g->voxelSize[2])
    {
        dimZ++;
    }
    rawData<double>* rawDataDWI = new rawData<double>(DICOM_DOUBLE, 3, dimX, dimY, dimZ);
    rawDataDWI->pixDimX = 1.0;
    rawDataDWI->pixDimY = 1.0;
    rawDataDWI->pixDimZ = 1.0;
    rawDataDWI->pixDimT = 1.0;
    (*rawDataDWI) = 0.0;
    double a = (alphamax - alphamin)/( rmax - rmin);
    double b = alphamin - a*rmin;
    double alpha;

    double Z;

    double V1[3];// = {1.0, 0.0, 0.0};
    double V2[3];// = {0.0, 1.0, 0.0};
    double V3[3];// = {0.0, 0.0, 1.0};
    double* D = new double[6];
    double vp1 = 1.0, vp2 = 0.2, vp3 = 0.2;
    unsigned long i=0, j=0, k=0;

    ///for each pixel in rawDataDWI, set the 3 eigen vectors and compute DWI data
    i=0;
    for(double x=xmin ; x<=xmax ; x+=g->voxelSize[0], i++)
    {
        j=0;
        for(double y=ymin ; y<=ymax ; y+=g->voxelSize[1], j++)
        {
            r = sqrt(SQR(x-Xc) + SQR(y-Yc));
            if(r<rmax && r>rmin)
            {
                if(ABS(a*r+b-(Pi/2))<0.000001 || ABS(a*r+b+(Pi/2))<0.000001)
                {
                    k=0;
                    for(double z=zmin ; z<=zmax ; z+=g->voxelSize[2], k++)
                    {
                        V1[0] = 0.0;
                        V1[1] = 0.0;
                        V1[2] = 1.0;

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
                else
                {
                    alpha = tan(a*r+b);
                    k=0;
                    for(double z=zmin ; z<=zmax ; z+=g->voxelSize[2], k++)
                    {
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
    delete D;
    return rawDataDWI;
}
