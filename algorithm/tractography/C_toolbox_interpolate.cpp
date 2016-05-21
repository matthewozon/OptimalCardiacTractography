#include <C_toolbox_interpolate.h>




//C_toolbox_interpolate::C_toolbox_interpolate(C_graph* G, double r):C_toolbox()
//{
//    //ctor
//    //paramOK = true;
//    R = r;
//    sigma = R/3.0;
//    m_graph = G;
//}
C_toolbox_interpolate::C_toolbox_interpolate(rawData<double>* dwiOrDtiData, rawData<double>* maskRawData, double r):C_toolbox()
{
    //ctor
    //paramOK = true;
    R = r;
    sigma = R/3.0;
    m_dwiOrDtiData = dwiOrDtiData;
    m_maskRawData = maskRawData;
}

C_toolbox_interpolate::~C_toolbox_interpolate()
{
    //dtor
}


///identically sample coordinate points
void C_toolbox_interpolate::regularInterpolatePoint(C_point_array_data* p, double* vStart, double* vEnd)
{
    //here it's assumed that vStart and vEnd are of dimension = dim
    if(p->nb_point>1)
    {
        double t = 0;
        for(unsigned long i=0 ; i<p->nb_point ; i++)
        {
            t = ((double) i)/((double) p->nb_point - 1.0);
            for(unsigned long j=0 ; j<p->dim ; j++)
                p->pointCoordinate[i][j] = vStart[j] + t * (vEnd[j] - vStart[j]);
        }
    }
    else
    {
        for(unsigned long j=0 ; j<p->dim ; j++)
            p->pointCoordinate[0][j] = 0.5 * ( vStart[j] + vEnd[j] );
    }


    return;
}
void C_toolbox_interpolate::regularInterpolatePoint(C_point_array_data* p, Point* vStart, Point* vEnd)
{
    //here it's assumed that vStart and vEnd are of dimension = dim
//    cout << "in interpolate" << endl;
    if(p->nb_point>1)
    {
        double t = 0;
        for(unsigned long i=0 ; i<p->nb_point ; i++)
        {
            t = ((double) i)/((double) p->nb_point - 1.0);
            p->pointCoordinate[i][0] = vStart->x + t * (vEnd->x - vStart->x);
            p->pointCoordinate[i][1] = vStart->y + t * (vEnd->y - vStart->y);
            p->pointCoordinate[i][2] = vStart->z + t * (vEnd->z - vStart->z);
        }
    }
    else
    {
//        cout << "interpolating" << endl;
//        cout << "p->pointCoordinate " << p->pointCoordinate << endl;
//        cout << "p->pointCoordinate[0] " << p->pointCoordinate[0] << endl;
        p->pointCoordinate[0][0] = 0.5 * ( vStart->x + vEnd->x );
        p->pointCoordinate[0][1] = 0.5 * ( vStart->y + vEnd->y );
        p->pointCoordinate[0][2] = 0.5 * ( vStart->z + vEnd->z );
//        cout << "done" << endl;
    }


    return;
}



//calling function
unsigned long C_toolbox_interpolate::run(C_point_array_data* p, unsigned long METHOD, unsigned long OBJECT_TYPE, bool yes, double** Hpsi, double b)
{
    if(p!=NULL)
    {
        if(p->pointCoordinate!=NULL && p->dim==3 && p->vectorInterpolated!=NULL && p->dimV==3 && p->nb_point>0)  //actually, dim and dimV must be 3!!!!
        {
            if( (METHOD==GAUSSIAN_PHI) | (METHOD==INV_SQRT) | (METHOD==CAUCHY_PHI) || (METHOD==NEAREST_NEIGHBOR) || (METHOD==LINEAR_INTERPOLATION))
            {
                if(OBJECT_TYPE==TENSOR)
                {
                    if(!interpolateFromRawDataDTI(p, METHOD, yes)) return ENDED_NOT_SUCCESS;
                }
                else if(OBJECT_TYPE==DWIS)
                {
                    if(!interpolateFromRawDataDWI(p, Hpsi, m_dwiOrDtiData->DimT-1 /**NgradDir*/, b /**b-value*/, METHOD, yes)) return ENDED_NOT_SUCCESS;
                }
                else
                {
                    return BAD_PARAM;
                }
            }
            else
            {
                //there is nothing else till now
                return BAD_PARAM;
            }
        }
        else
        {
            return BAD_PARAM;
        }
    }
    else
    {
        return NULL_POINTER;
    }

    return _SUCCESS;
}




bool C_toolbox_interpolate::interpolateFromRawDataDTI(C_point_array_data* p, unsigned long METHOD, bool yes) ///assume 2D or 3D +1D for tensor storage
{
    for(unsigned long i=0 ; i<p->nb_point ; i++)
    {
        ///run interpolation on one point
        interpolateFromRawDataDTI(p->pointCoordinate[i], METHOD, p->vectorInterpolated[i], &(p->FA[i]));
    }
    return true;
}
bool C_toolbox_interpolate::interpolateFromRawDataDWI(C_point_array_data* p, double** Hpsi, unsigned long NgradDir, double b, unsigned long METHOD, bool yes) ///assume 2D or 3D + 1D for diffusion storage (averaging alredy computed)
{
    for(unsigned long i=0 ; i<p->nb_point ; i++)
    {
        ///run interpolation on one point
        interpolateFromRawDataDWI(p->pointCoordinate[i], Hpsi, NgradDir, b, METHOD, p->vectorInterpolated[i], &(p->FA[i]));
    }
    return true;
}

bool C_toolbox_interpolate::interpolateFromRawDataDTI(double* p, unsigned long METHOD, double* unitVect, double* fa)
{
    if(m_dwiOrDtiData->DimT!=6) return false;
    ///interpolate on each slice (each component of the tensor)
    double** DTI = new double*[3];
    DTI[0] = new double[3];
    DTI[1] = new double[3];
    DTI[2] = new double[3];

    ///get current voxel indices ("nearest" neighbor)
    signed long i = (signed long) floor((p[0]/(m_dwiOrDtiData->pixDimX))); //(1.0/2.0) +
    signed long j = (signed long) floor((p[1]/(m_dwiOrDtiData->pixDimY))); //(1.0/2.0) +
    signed long k = (signed long) floor((p[2]/(m_dwiOrDtiData->pixDimZ))); //(1.0/2.0) +


    ///check if current point is in mask or close enough (it can happen that on a border the current nearest voxel falls right next to the mask, but it makes no sens to discard it)
    bool closeEnoughToMask = false;
    //cout << "coucou current indices are: " << i << " " << j << " " << k << " point coordinates: " << p[0] << " " << p[1] << " " << p[2] << endl;
    //cout << "mask dim " << m_maskRawData->DimX << " " << m_maskRawData->DimY << " " << m_maskRawData->DimZ << endl;
    if(m_maskRawData->raw3D[i][j][k]<=0.00001)
    {
        ///look for neighboring voxels
        for(signed long di=-1 ; di<=1 ; di++)
        {
            for(signed long dj=-1 ; dj<=1 ; dj++)
            {
                for(signed long dk=-1 ; dk<=1 ; dk++)
                {
                    if( ((i+di)>=0) && ((j+dj)>=0) && ((k+dk)>=0) && ((i+di)<(signed long) m_dwiOrDtiData->DimX) && ((j+dj)<(signed long) m_dwiOrDtiData->DimY) && ((k+dk)<(signed long) m_dwiOrDtiData->DimZ))
                    {
                        if(m_maskRawData->raw3D[i+di][j+dj][k+dk]>0.00001)
                        {
                            ///either set the new value of indices (or find the other voxel close in the mask and select the colsest one)
                            //i = i+di;
                            //j = j+dj;
                            //k = k+dk;
                            closeEnoughToMask = true;
                            di = 2;
                            dj = 2;
                            dk = 2;
                        }
                    }
                }
            }
        }
    }
    else
    {
        closeEnoughToMask = true;
    }
    //cout << "after checking mask" << endl;

    ///actual interpolation
    if(METHOD==NEAREST_NEIGHBOR && closeEnoughToMask)
    {
        ///compute "interpolation" and fill in the container
        DTI[0][0] = m_dwiOrDtiData->raw4D[i][j][k][0]; DTI[0][1] = m_dwiOrDtiData->raw4D[i][j][k][1]; DTI[0][2] = m_dwiOrDtiData->raw4D[i][j][k][2];
        DTI[1][0] = m_dwiOrDtiData->raw4D[i][j][k][1]; DTI[1][1] = m_dwiOrDtiData->raw4D[i][j][k][3]; DTI[1][2] = m_dwiOrDtiData->raw4D[i][j][k][4];
        DTI[2][0] = m_dwiOrDtiData->raw4D[i][j][k][2]; DTI[2][1] = m_dwiOrDtiData->raw4D[i][j][k][4]; DTI[2][2] = m_dwiOrDtiData->raw4D[i][j][k][5];
    }
    else if(METHOD==LINEAR_INTERPOLATION && closeEnoughToMask)
    {
        ///get the 7 other points (actually, just check if they are in the mask)
        unsigned long i0 = i, i1 = i+1;
        unsigned long j0 = j, j1 = j+1;
        unsigned long k0 = k, k1 = k+1;

        ///on upper boundaries, use nearest neighbor: should not happen in real case because heart does not fit the whole acquisition matrix
        if(i1>=m_dwiOrDtiData->DimX || j1>=m_dwiOrDtiData->DimY || k1>=m_dwiOrDtiData->DimZ)
        {
            DTI[0][0] = m_dwiOrDtiData->raw4D[i][j][k][0]; DTI[0][1] = m_dwiOrDtiData->raw4D[i][j][k][1]; DTI[0][2] = m_dwiOrDtiData->raw4D[i][j][k][2];
            DTI[1][0] = m_dwiOrDtiData->raw4D[i][j][k][1]; DTI[1][1] = m_dwiOrDtiData->raw4D[i][j][k][3]; DTI[1][2] = m_dwiOrDtiData->raw4D[i][j][k][4];
            DTI[2][0] = m_dwiOrDtiData->raw4D[i][j][k][2]; DTI[2][1] = m_dwiOrDtiData->raw4D[i][j][k][4]; DTI[2][2] = m_dwiOrDtiData->raw4D[i][j][k][5];
        }
        else //you're free to compute trilinear interpolation on each components
        {
            ///http://en.wikipedia.org/wiki/Trilinear_interpolation
            double x = p[0], y = p[1] , z = p[2];
            double x0 = (((double) i0))*(m_dwiOrDtiData->pixDimX); // + 0.5
            double y0 = (((double) j0))*(m_dwiOrDtiData->pixDimY); // + 0.5
            double z0 = (((double) k0))*(m_dwiOrDtiData->pixDimZ); // + 0.5
            double x1 = (((double) i1))*(m_dwiOrDtiData->pixDimX); // + 0.5
            double y1 = (((double) j1))*(m_dwiOrDtiData->pixDimY); // + 0.5
            double z1 = (((double) k1))*(m_dwiOrDtiData->pixDimZ); // + 0.5
            double xd = (x-x0)/(x1-x0);
            double yd = (y-y0)/(y1-y0);
            double zd = (z-z0)/(z1-z0);

            ///for all elements of the last dimension
            ///DTI[0][0] = m_dwiOrDtiData->raw4D[i][j][k][0]
            double* c = new double[6];
            for(unsigned long l=0 ; l<6 ; l++)
            {
                double c00 = m_dwiOrDtiData->raw4D[i0][j0][k0][l]*(1.0-xd) + m_dwiOrDtiData->raw4D[i1][j0][k0][l]*xd;
                double c01 = m_dwiOrDtiData->raw4D[i0][j0][k1][l]*(1.0-xd) + m_dwiOrDtiData->raw4D[i1][j0][k1][l]*xd;
                double c10 = m_dwiOrDtiData->raw4D[i0][j1][k0][l]*(1.0-xd) + m_dwiOrDtiData->raw4D[i1][j1][k0][l]*xd;
                double c11 = m_dwiOrDtiData->raw4D[i0][j1][k1][l]*(1.0-xd) + m_dwiOrDtiData->raw4D[i1][j1][k1][l]*xd;

                double c0 = c00*(1.0-yd) + c10*yd;
                double c1 = c01*(1.0-yd) + c11*yd;

                c[l] = c0*(1.0-zd) + c1*zd;
            }
            DTI[0][0] = c[0]; DTI[0][1] = c[1]; DTI[0][2] = c[2];
            DTI[1][0] = c[1]; DTI[1][1] = c[3]; DTI[1][2] = c[4];
            DTI[2][0] = c[2]; DTI[2][1] = c[4]; DTI[2][2] = c[5];
            delete c;
        }
    }
    else if( (METHOD==GAUSSIAN_PHI || METHOD==CAUCHY_PHI || METHOD==INV_SQRT) && closeEnoughToMask)
    {
        ///count number of point in the intersection of mask and ball centered on p of radius R
        //look for point that are actually reachable +/-0.5*R arround p
        signed long imin = (signed long) floor(((p[0]-R)/(m_dwiOrDtiData->pixDimX)));
        signed long jmin = (signed long) floor(((p[1]-R)/(m_dwiOrDtiData->pixDimY)));
        signed long kmin = (signed long) floor(((p[2]-R)/(m_dwiOrDtiData->pixDimZ)));
        signed long imax = (signed long) floor(((p[0]+R)/(m_dwiOrDtiData->pixDimX)));
        signed long jmax = (signed long) floor(((p[1]+R)/(m_dwiOrDtiData->pixDimY)));
        signed long kmax = (signed long) floor(((p[2]+R)/(m_dwiOrDtiData->pixDimZ)));
        unsigned long nbPointInIntersection = 0;
        for(signed long q=max(imin,(signed long) 0) ; q<=min(imax,(signed long) m_dwiOrDtiData->DimX-1) ; q++)
        {
            double x = (((double) q+ 0.5))*(m_dwiOrDtiData->pixDimX); //center of pixel
            x = SQR(p[0]-x);
            for(signed long r=max(jmin,(signed long) 0) ; r<=min(jmax,(signed long) m_dwiOrDtiData->DimY-1) ; r++)
            {
                double y = (((double) r+ 0.5))*(m_dwiOrDtiData->pixDimY); //center of pixel
                y = SQR(p[1]-y);
                for(signed long s=max(kmin,(signed long) 0) ; s<=min(kmax,(signed long) m_dwiOrDtiData->DimZ-1) ; s++)
                {
                    double z = (((double) s + 0.5))*(m_dwiOrDtiData->pixDimZ); //center of pixel
                    z = SQR(p[2]-z);
                    if( (x+y+z) < SQR(R) && (m_maskRawData->raw3D[q][r][s]>0.5)) ///it's faster to compute SQR than sqrt
                    {
                        ///indent number of points in intersection
                        nbPointInIntersection++;
                    }
                }
            }
        }

        if(nbPointInIntersection==0)
        {
            DTI[0][0] = m_dwiOrDtiData->raw4D[i][j][k][0]; DTI[0][1] = m_dwiOrDtiData->raw4D[i][j][k][1]; DTI[0][2] = m_dwiOrDtiData->raw4D[i][j][k][2];
            DTI[1][0] = m_dwiOrDtiData->raw4D[i][j][k][1]; DTI[1][1] = m_dwiOrDtiData->raw4D[i][j][k][3]; DTI[1][2] = m_dwiOrDtiData->raw4D[i][j][k][4];
            DTI[2][0] = m_dwiOrDtiData->raw4D[i][j][k][2]; DTI[2][1] = m_dwiOrDtiData->raw4D[i][j][k][4]; DTI[2][2] = m_dwiOrDtiData->raw4D[i][j][k][5];
        }
        else
        {
            ///allocate a weighting array
            double* Warray = new double[nbPointInIntersection];
            double Wg = 0.0;
            unsigned long u=0;
            ///initialize container
            DTI[0][0] = 0.0; DTI[0][1] = 0.0; DTI[0][2] = 0.0;
            DTI[1][0] = 0.0; DTI[1][1] = 0.0; DTI[1][2] = 0.0;
            DTI[2][0] = 0.0; DTI[2][1] = 0.0; DTI[2][2] = 0.0;
            ///compute interpolation using all points in intersection
            for(signed long q=max(imin,(signed long) 0) ; q<=min(imax,(signed long) m_dwiOrDtiData->DimX-1) ; q++)
            {
                double x = (((double) q+ 0.5))*(m_dwiOrDtiData->pixDimX); //center of pixel
                x = SQR(p[0]-x);
                for(signed long r=max(jmin,(signed long) 0) ; r<=min(jmax,(signed long) m_dwiOrDtiData->DimY-1) ; r++)
                {
                    double y = (((double) r+ 0.5))*(m_dwiOrDtiData->pixDimY); //center of pixel
                    y = SQR(p[1]-y);
                    for(signed long s=max(kmin,(signed long) 0) ; s<=min(kmax,(signed long) m_dwiOrDtiData->DimZ-1) ; s++)
                    {
                        double z = (((double) s + 0.5))*(m_dwiOrDtiData->pixDimZ); //center of pixel
                        z = SQR(p[2]-z);
                        if( (x+y+z) < SQR(R) && (m_maskRawData->raw3D[q][r][s]>0.5)) ///it's faster to compute SQR than sqrt
                        {
                            ///compute weight
                            if(METHOD==GAUSSIAN_PHI)
                            {
                                Warray[u] = gaussianPhi(sqrt(x+y+z));
                            }
                            else if(METHOD==CAUCHY_PHI)
                            {
                                Warray[u] = cauchyPhi(sqrt(x+y+z));
                            }
                            else //INV_SQRT
                            {
                                Warray[u] = invAbsPhi(sqrt(x+y+z));
                            }
                            ///add the weighted value to the global sum (+=)
                            DTI[0][0] += Warray[u]*m_dwiOrDtiData->raw4D[q][r][s][0]; DTI[0][1] += Warray[u]*m_dwiOrDtiData->raw4D[q][r][s][1]; DTI[0][2] += Warray[u]*m_dwiOrDtiData->raw4D[q][r][s][2];
                            DTI[1][0] += Warray[u]*m_dwiOrDtiData->raw4D[q][r][s][1]; DTI[1][1] += Warray[u]*m_dwiOrDtiData->raw4D[q][r][s][3]; DTI[1][2] += Warray[u]*m_dwiOrDtiData->raw4D[q][r][s][4];
                            DTI[2][0] += Warray[u]*m_dwiOrDtiData->raw4D[q][r][s][2]; DTI[2][1] += Warray[u]*m_dwiOrDtiData->raw4D[q][r][s][4]; DTI[2][2] += Warray[u]*m_dwiOrDtiData->raw4D[q][r][s][5];
                            ///add the weight to global weight (normalization factor)
                            Wg += Warray[u];
                            ///indent iterator
                            u++;
                            if(u==nbPointInIntersection)
                            {
                                q = imax+1;
                                r = jmax+1;
                                s = kmax+1;
                            }
                        }
                    }
                }
            }

            ///normalize he value
            DTI[0][0] /= Wg; DTI[0][1] /= Wg; DTI[0][2] /= Wg;
            DTI[1][0] /= Wg; DTI[1][1] /= Wg; DTI[1][2] /= Wg;
            DTI[2][0] /= Wg; DTI[2][1] /= Wg; DTI[2][2] /= Wg;

            ///delete
            delete Warray;
        }
    }
    else
    {
        delete DTI[0];
        delete DTI[1];
        delete DTI[2];
        delete DTI;
        return false;
    }

//    if(isnan(DTI[0][0]) || isnan(DTI[0][1]) || isnan(DTI[0][2]) || isnan(DTI[1][0]) || isnan(DTI[1][1]) || isnan(DTI[1][2]) || isnan(DTI[2][0]) || isnan(DTI[2][1]) || isnan(DTI[2][2]) )
//    {
//        cout << "NaN appears in DTI interpolation" << endl;
//    }

    ///extract eigen vectors and eigen values
    C_toolbox_eigen_sym* toolEigen2 = new C_toolbox_eigen_sym(DTI);
    ///store what must be stored
    *fa = FA(toolEigen2->d[0], toolEigen2->d[1], toolEigen2->d[2]);
    if(isnan(*fa))
    {
//        cout << DTI[0][0] << " " << DTI[0][1] << " " << DTI[0][2] << endl;
//        cout << DTI[1][0] << " " << DTI[1][1] << " " << DTI[1][2] << endl;
//        cout << DTI[2][0] << " " << DTI[2][1] << " " << DTI[2][2] << endl;
//        cout << "fa is NaN" << endl;
        *fa = 0.0;
        unitVect[0] = 1.0;
        unitVect[1] = 0.0;
        unitVect[2] = 0.0;
        delete toolEigen2;//->~C_toolbox_eigen_sym();
        delete DTI[0];
        delete DTI[1];
        delete DTI[2];
        delete[] DTI;
//        if(m_maskRawData->raw3D[i][j][k]<=0.00001) cout << "probably because the point was not in the mask" << endl;
        return false; ///check whether the point is inside the mask
    }
    double N = sqrt( SQR( toolEigen2->z[0][0] ) + SQR( toolEigen2->z[1][0] ) + SQR( toolEigen2->z[2][0] ) );
    unitVect[0] = toolEigen2->z[0][0]/N;
    unitVect[1] = toolEigen2->z[1][0]/N;
    unitVect[2] = toolEigen2->z[2][0]/N;
//    if(isnan(unitVect[0]) || isnan(unitVect[1]) || isnan(unitVect[2])) cout << "eigen vector is NaN" << endl;

    delete toolEigen2;//->~C_toolbox_eigen_sym();
    delete DTI[0];
    delete DTI[1];
    delete DTI[2];
    delete[] DTI;

    return true;
}
bool C_toolbox_interpolate::interpolateFromRawDataDWI(double* p, double** Hpsi, unsigned long NgradDir, double b, unsigned long METHOD, double* unitVect, double* fa)
{
    if(m_dwiOrDtiData->DimT<7) return false;
    ///interpolate on each slice (each diffusion direction)
    double* DWI = new double[m_dwiOrDtiData->DimT];

    ///get current voxel indices ("nearest" neighbor)
    signed long i = (signed long) floor((p[0]/(m_dwiOrDtiData->pixDimX))); //(1.0/2.0) +
    signed long j = (signed long) floor((p[1]/(m_dwiOrDtiData->pixDimY))); //(1.0/2.0) +
    signed long k = (signed long) floor((p[2]/(m_dwiOrDtiData->pixDimZ))); //(1.0/2.0) +

    ///check if current point is in mask or close enough (it can happen that on a border the current nearest voxel falls right next to the mask, but it makes no sens to discard it)
    bool closeEnoughToMask = false;
    if(m_maskRawData->raw3D[i][j][k]<=0.00001)
    {
        ///look for neighboring voxels
        for(signed long di=-1 ; di<=1 ; di++)
        {
            for(signed long dj=-1 ; dj<=1 ; dj++)
            {
                for(signed long dk=-1 ; dk<=1 ; dk++)
                {
                    if(i+di>=0 && j+dj>=0 && k+dk>=0 && i+di<(signed long) m_dwiOrDtiData->DimX && j+dj<(signed long) m_dwiOrDtiData->DimY && k+dk<(signed long) m_dwiOrDtiData->DimZ)
                    {
                        if(m_maskRawData->raw3D[i+di][j+dj][k+dk]>0.00001)
                        {
                            ///either set the new value of indices (or find the other voxel close in the mask and select the colsest one)
                            //i = i+di;
                            //j = j+dj;
                            //k = k+dk;
                            closeEnoughToMask = true;
                            di = 2;
                            dj = 2;
                            dk = 2;
                        }
                    }
                }
            }
        }
    }
    else
    {
        closeEnoughToMask = true;
    }

    if(METHOD==NEAREST_NEIGHBOR && closeEnoughToMask)
    {
        ///compute "interpolation" and fill in the container
        for(unsigned long l=0 ; l<m_dwiOrDtiData->DimT ; l++)
            DWI[l] = m_dwiOrDtiData->raw4D[i][j][k][l];
    }
    else if(METHOD==LINEAR_INTERPOLATION && closeEnoughToMask)
    {
        ///get the 7 other points (actually, just check if they are in the mask)
        unsigned long i0 = i, i1 = i+1;
        unsigned long j0 = j, j1 = j+1;
        unsigned long k0 = k, k1 = k+1;

        ///on upper boundaries, use nearest neighbor: should not happen in real case because heart does not fit the whole acquisition matrix
        if(i1>=m_dwiOrDtiData->DimX || j1>=m_dwiOrDtiData->DimY || k1>=m_dwiOrDtiData->DimZ)
        {
            for(unsigned long l=0 ; l<m_dwiOrDtiData->DimT ; l++)
                DWI[l] = m_dwiOrDtiData->raw4D[i][j][k][l];
        }
        else //you're free to compute trilinear interpolation on each components
        {
            ///http://en.wikipedia.org/wiki/Trilinear_interpolation
            double x = p[0], y = p[1] , z = p[2];
            double x0 = (((double) i0))*(m_dwiOrDtiData->pixDimX); // + 0.5
            double y0 = (((double) j0))*(m_dwiOrDtiData->pixDimY); // + 0.5
            double z0 = (((double) k0))*(m_dwiOrDtiData->pixDimZ); // + 0.5
            double x1 = (((double) i1))*(m_dwiOrDtiData->pixDimX); // + 0.5
            double y1 = (((double) j1))*(m_dwiOrDtiData->pixDimY); // + 0.5
            double z1 = (((double) k1))*(m_dwiOrDtiData->pixDimZ); // + 0.5
            double xd = (x-x0)/(x1-x0);
            double yd = (y-y0)/(y1-y0);
            double zd = (z-z0)/(z1-z0);

            ///for all elements of the last dimension
            for(unsigned long l=0 ; l<m_dwiOrDtiData->DimT ; l++)
            {
                double c00 = m_dwiOrDtiData->raw4D[i0][j0][k0][l]*(1.0-xd) + m_dwiOrDtiData->raw4D[i1][j0][k0][l]*xd;
                double c01 = m_dwiOrDtiData->raw4D[i0][j0][k1][l]*(1.0-xd) + m_dwiOrDtiData->raw4D[i1][j0][k1][l]*xd;
                double c10 = m_dwiOrDtiData->raw4D[i0][j1][k0][l]*(1.0-xd) + m_dwiOrDtiData->raw4D[i1][j1][k0][l]*xd;
                double c11 = m_dwiOrDtiData->raw4D[i0][j1][k1][l]*(1.0-xd) + m_dwiOrDtiData->raw4D[i1][j1][k1][l]*xd;

                double c0 = c00*(1.0-yd) + c10*yd;
                double c1 = c01*(1.0-yd) + c11*yd;

                DWI[l] = c0*(1.0-zd) + c1*zd;
            }
        }
    }
    else if((METHOD==GAUSSIAN_PHI || METHOD==CAUCHY_PHI || METHOD==INV_SQRT) && closeEnoughToMask)
    {
        ///count number of point in the intersection of mask and ball centered on p of radius R
        //look for point that are actually reachable +/-0.5*R arround p
        signed long imin = (signed long) floor(((p[0]-R)/(m_dwiOrDtiData->pixDimX)));
        signed long jmin = (signed long) floor(((p[1]-R)/(m_dwiOrDtiData->pixDimY)));
        signed long kmin = (signed long) floor(((p[2]-R)/(m_dwiOrDtiData->pixDimZ)));
        signed long imax = (signed long) floor(((p[0]+R)/(m_dwiOrDtiData->pixDimX)));
        signed long jmax = (signed long) floor(((p[1]+R)/(m_dwiOrDtiData->pixDimY)));
        signed long kmax = (signed long) floor(((p[2]+R)/(m_dwiOrDtiData->pixDimZ)));
        unsigned long nbPointInIntersection = 0;
        for(signed long q=max(imin,(signed long) 0) ; q<=min(imax,(signed long) m_dwiOrDtiData->DimX-1) ; q++)
        {
            double x = (((double) q) + 0.5)*(m_dwiOrDtiData->pixDimX);
            x = SQR(p[0]-x);
            for(signed long r=max(jmin,(signed long) 0) ; r<=min(jmax,(signed long) m_dwiOrDtiData->DimY-1) ; r++)
            {
                double y = (((double) r) + 0.5)*(m_dwiOrDtiData->pixDimY);
                y = SQR(p[1]-y);
                for(signed long s=max(kmin,(signed long) 0) ; s<=min(kmax,(signed long) m_dwiOrDtiData->DimZ-1) ; s++)
                {
                    double z = (((double) s) + 0.5)*(m_dwiOrDtiData->pixDimZ);
                    z = SQR(p[2]-z);
                    if( (x+y+z) < SQR(R) && (m_maskRawData->raw3D[q][r][s]>0.5)) ///it's faster to compute SQR than sqrt
                    {
                        ///indent number of points in intersection
                        nbPointInIntersection++;
                    }
                }
            }
        }

        if(nbPointInIntersection==0)
        {
            for(unsigned long l=0 ; l<m_dwiOrDtiData->DimT ; l++)
                DWI[l] = m_dwiOrDtiData->raw4D[i][j][k][l];
        }
        else
        {
            ///allocate a weighting array
            double* W = new double[nbPointInIntersection];
            double Wg = 0.0;
            unsigned long u=0;
            ///initialize container
            for(unsigned long l=0 ; l<m_dwiOrDtiData->DimT ; l++)
                DWI[l] = 0.0;

            ///compute interpolation using all points in intersection
            for(signed long q=max(imin,(signed long) 0) ; q<=min(imax,(signed long) m_dwiOrDtiData->DimX-1) ; q++)
            {
                double x = (((double) q) + 0.5)*(m_dwiOrDtiData->pixDimX);
                x = SQR(p[0]-x);
                for(signed long r=max(jmin,(signed long) 0) ; r<=min(jmax,(signed long) m_dwiOrDtiData->DimY-1) ; r++)
                {
                    double y = (((double) r) + 0.5)*(m_dwiOrDtiData->pixDimY);
                    y = SQR(p[1]-y);
                    for(signed long s=max(kmin,(signed long) 0) ; s<=min(kmax,(signed long) m_dwiOrDtiData->DimZ-1) ; s++)
                    {
                        double z = (((double) s) + 0.5)*(m_dwiOrDtiData->pixDimZ);
                        z = SQR(p[2]-z);
                        if( (x+y+z) < SQR(R) && (m_maskRawData->raw3D[q][r][s]>0.5)) ///it's faster to compute SQR than sqrt
                        {
                            ///compute weight
                            if(METHOD==GAUSSIAN_PHI)
                            {
                                W[u] = gaussianPhi(sqrt(x+y+z));
                            }
                            else if(METHOD==CAUCHY_PHI)
                            {
                                W[u] = cauchyPhi(sqrt(x+y+z));
                            }
                            else //INV_SQRT
                            {
                                W[u] = invAbsPhi(sqrt(x+y+z));
                            }
                            ///add the weighted value to the global sum (+=)
                            for(unsigned long l=0 ; l<m_dwiOrDtiData->DimT ; l++)
                                DWI[l] += W[u]*m_dwiOrDtiData->raw4D[q][r][s][l];

                            ///add the weight to global weight (normalization factor)
                            Wg += W[u];
                            ///indent iterator
                            u++;
                            if(u==nbPointInIntersection)
                            {
                                q = imax+1;
                                r = jmax+1;
                                s = kmax+1;
                            }
                        }
                    }
                }
            }
            delete W;
            ///normalize he value
            for(unsigned long l=0 ; l<m_dwiOrDtiData->DimT ; l++)
                DWI[l] /= Wg;
        }
    }
    else
    {
        delete DWI;
        return false;
    }

    ///compute DT from DW
    C_tensorMaker* toolSVD = new C_tensorMaker();
    double** DTI = toolSVD->makeTensor2(&(DWI[1]), Hpsi, b, DWI[0], m_dwiOrDtiData->DimT-1);
    delete DWI;
    delete toolSVD;

    ///extract eigen vectors and eigen values
    C_toolbox_eigen_sym* tool = new C_toolbox_eigen_sym(DTI);

    ///store what must be stored
    *fa = FA(tool->d[0], tool->d[1], tool->d[2]);
    double N = sqrt( SQR( tool->z[0][0] ) + SQR( tool->z[1][0] ) + SQR( tool->z[2][0] ) );
    unitVect[0] = tool->z[0][0]/N;
    unitVect[1] = tool->z[1][0]/N;
    unitVect[2] = tool->z[2][0]/N;

    delete tool;
    delete DTI[0];
    delete DTI[1];
    delete DTI[2];
    delete DTI;
    return true;
}



//phi function
double C_toolbox_interpolate::gaussianPhi(double r)
{
    return exp(-0.5*SQR(r/sigma));
}
double C_toolbox_interpolate::gaussianPhiNorm(double r)
{
    return gaussianPhi(r)/(sigma*sqrt(2*Pi));
}
double C_toolbox_interpolate::invAbsPhi(double r)
{
    return sqrt(cauchyPhi(r));
}
double C_toolbox_interpolate::cauchyPhi(double r)
{
    return 1.0/(1.0 + SQR(r/sigma));
}
double C_toolbox_interpolate::cauchyPhiNorm(double r)
{
    return cauchyPhi(r)/(sigma*Pi);
}
