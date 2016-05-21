#ifndef C_SAMPLING_H
#define C_SAMPLING_H

#include <struct.h>
#include <rawData.h>

template<class T> class C_sampling;

template<class T> class C_thread_NN : public C_thread
{
    friend class C_sampling<T>;
    public:
        /** Default constructor */
        C_thread_NN(rawData<T>* X, samplingFactor* sf, rawData<T>* Y, long idxStart, long idxEnd)
        {
            if(idxStart<=idxEnd)
            {
                m_idxStart = idxStart;
                m_idxEnd = idxEnd;//,  A->numel());
            }
            else
            {
                m_idxStart = idxEnd;//,  A->numel());
                m_idxEnd = idxStart;
            }
            m_X = X;
            m_Y = Y;
            m_sf = sf;
        };

        /** Default destructor */
        virtual ~C_thread_NN()
        {
            //
        };
        virtual void execute()
        {
            if(m_X->numDim==1)
            {
                long iOld;
                for(long idx=m_idxStart ; idx<m_idxEnd ; idx++)
                {
                    iOld = (long) ((double) idx)/(m_sf->Xfactor);
                    m_Y->raw1D[idx] = m_X->raw1D[iOld];
                }
            }
            if(m_X->numDim==2)
            {
                long i, j;
                long iOld, jOld;
                ///idx = i*DimY + j;
                for(long idx=m_idxStart ; idx<m_idxEnd ; idx++)
                {
                    i = floor(idx/m_Y->DimY);
                    j = idx - (floor(idx/m_Y->DimY)*m_Y->DimY);
                    iOld = (long) (((double) i)/(m_sf->Xfactor));
                    jOld = (long) (((double) j)/(m_sf->Yfactor));
                    m_Y->raw2D[i][j] = m_X->raw2D[iOld][jOld];
                }
            }
            if(m_X->numDim==3)
            {
                long i, j, k;
                ///DimY*i + j + DimX*DimY*k
                long iOld, jOld, kOld;
                for(long idx=m_idxStart ; idx<m_idxEnd ; idx++)
                {
                    k = floor(idx/(m_Y->DimX*m_Y->DimY));
                    i = floor((idx - (k*(m_Y->DimX*m_Y->DimY)))/m_Y->DimY);
                    j = idx - (k*(m_Y->DimX*m_Y->DimY)) - i*m_Y->DimY;
                    iOld = (long) (((double) i)/(m_sf->Xfactor));
                    jOld = (long) (((double) j)/(m_sf->Yfactor));
                    kOld = (long) ((double) k)/(m_sf->Zfactor);
//                    if(iOld<0 || iOld>=m_X->DimX) cout << "iOld pb " << iOld << " dimX = " << m_X->DimX << endl;
//                    if(jOld<0 || jOld>=m_X->DimY) cout << "jOld pb " << jOld << " dimY = " << m_X->DimY << endl;
//                    if(kOld<0 || kOld>=m_X->DimZ) cout << "kOld pb " << kOld << " dimZ = " << m_X->DimZ << endl;

                    m_Y->raw3D[i][j][k] = m_X->raw3D[iOld][jOld][kOld];
                }
            }
            if(m_X->numDim==4)
            {
                long i, j, k, l;
                ///DimY*i + j + DimX*DimY*k + DimX*DimY*DimZ*l
                long iOld, jOld, kOld, lOld;
                for(long idx=m_idxStart ; idx<m_idxEnd ; idx++)
                {
                    l = floor(idx/(m_Y->DimX*m_Y->DimY*m_Y->DimZ));
                    k = floor((idx - l*(m_Y->DimX*m_Y->DimY*m_Y->DimZ))/(m_Y->DimX*m_Y->DimY));
                    i = floor((idx - l*(m_Y->DimX*m_Y->DimY*m_Y->DimZ) - k*(m_Y->DimX*m_Y->DimY))/m_Y->DimY);
                    j = idx - l*(m_Y->DimX*m_Y->DimY*m_Y->DimZ) - k*(m_Y->DimX*m_Y->DimY) - i*m_Y->DimY;
                    iOld = (long) (((double) i)/(m_sf->Xfactor));
                    jOld = (long) (((double) j)/(m_sf->Yfactor));
                    kOld = (long) ((double) k)/(m_sf->Zfactor);
                    lOld = (long) ((double) l)/(m_sf->Tfactor);
                    m_Y->raw4D[i][j][k][l] = m_X->raw4D[iOld][jOld][kOld][lOld];
                }
            }
        };
    protected:
    private:
        long m_idxStart;
        long m_idxEnd;
        rawData<T>* m_X;
        samplingFactor* m_sf;
        rawData<T>* m_Y;
};







template<class T> class C_thread_linear : public C_thread
{
    friend class C_sampling<T>;
    public:
        /** Default constructor */
        C_thread_linear(rawData<T>* X, samplingFactor* sf, rawData<T>* Y, long idxStart, long idxEnd)
        {
            if(idxStart<=idxEnd)
            {
                m_idxStart = idxStart;
                m_idxEnd = idxEnd;//,  A->numel());
            }
            else
            {
                m_idxStart = idxEnd;//,  A->numel());
                m_idxEnd = idxStart;
            }
            m_X = X;
            m_Y = Y;
            m_sf = sf;
        };

        /** Default destructor */
        virtual ~C_thread_linear()
        {
            //
        };
        virtual void execute()
        {
            if(m_X->numDim==1)
            {
                unsigned long xUp, xDown;
                double xCenter;
                for(long idx=m_idxStart ; idx<m_idxEnd ; idx++)
                {

                    //get the actual position of the pixel in former signal
                    xCenter = ((double) idx)/(m_sf->Xfactor);

                    //get the surrounding pixels
                    xDown = (unsigned long) xCenter;
                    xUp = xDown + 1;

                    //compute linear interpolation
                    if(xDown==m_X->DimX-1)
                    {
                        xUp = xDown;
                        m_Y->raw1D[idx] = m_X->raw1D[xDown];
                    }
                    else
                    {
                        m_Y->raw1D[idx] = (((double) xUp) - xCenter)*m_X->raw1D[xDown] + (xCenter - ((double) xDown))*m_X->raw1D[xUp];
                    }
                }
            }
            if(m_X->numDim==2)
            {
                long i, j;
                long iOld, jOld;
                ///idx = i*DimY + j;
                unsigned long xUp, xDown, yUp, yDown;
                double xCenter, yCenter,w00, w10, w11, w01, dx, dy;
                for(long idx=m_idxStart ; idx<m_idxEnd ; idx++)
                {
                    i = floor(idx/m_Y->DimY);
                    j = idx - (floor(idx/m_Y->DimY)*m_Y->DimY);
                    //get the actual x-position of the pixel in former signal
                    xCenter = ((double) i)/(m_sf->Xfactor);

                    //get the surrounding pixels on x-dimension
                    xDown = (unsigned long) xCenter;
                    xUp = xDown + 1;
                    dx = xUp - xDown;

                    if(xUp==m_X->DimX)
                    {
                        xUp = xDown;
                        dx = 1.0;
                        //get the actual x-position of the pixel in former signal
                        yCenter = ((double) j)/(m_sf->Yfactor);

                        //get the surrounding pixels on x-dimension
                        yDown = (unsigned long) yCenter;
                        yUp = yDown + 1;
                        if(yUp==m_X->DimY)
                        {
                            yUp = yDown;
                            dy = 1.0;
                            m_Y->raw2D[i][j] = m_X->raw2D[xDown][yDown];
                        }
                        else
                        {
                            m_Y->raw2D[i][j] = (yUp-yCenter)*m_X->raw2D[xDown][yDown] + (yCenter-yDown)*m_X->raw2D[xDown][yUp];
                        }
                    }
                    else
                    {
                        //get the actual x-position of the pixel in former signal
                        yCenter = ((double) j)/(m_sf->Yfactor);

                        //get the surrounding pixels on x-dimension
                        yDown = (unsigned long) yCenter;
                        yUp = yDown + 1;
                        dy = yUp - yDown;
                        if(yUp==m_X->DimY)
                        {
                            yUp = yDown;
                            dy = 1.0;
                            m_Y->raw2D[i][j] = (xUp-xCenter)*m_X->raw2D[xDown][yDown] + (xCenter-xDown)*m_X->raw2D[xUp][yDown];
                        }
                        else
                        {
                            w00 = (xUp - xCenter)*(yUp - yCenter); //down dwon
                            w10 = (xCenter - xDown)*(yUp - yCenter); //up down
                            w01 = (xUp - xCenter)*(yCenter - yDown); //down up
                            w11 = (xCenter - xDown)*(yCenter - yDown); //up up
                            m_Y->raw2D[i][j] = (1.0/(dx*dy))*(w00*m_X->raw2D[xDown][yDown] + w10*m_X->raw2D[xUp][yDown] + w01*m_X->raw2D[xDown][yUp] + w11*m_X->raw2D[xUp][yUp]);
                        }
                    }

                }
            }
            if(m_X->numDim==3)
            {
                long i, j, k;
                ///DimY*i + j + DimX*DimY*k
                unsigned long xUp, xDown, yUp, yDown, zUp, zDown;
                double xCenter, yCenter, zCenter, w000, w100, w110, w010, w001, w101, w111, w011, dx, dy, dz;
                for(long idx=m_idxStart ; idx<m_idxEnd ; idx++)
                {
                    k = floor(idx/(m_Y->DimX*m_Y->DimY));
                    i = floor((idx - (k*(m_Y->DimX*m_Y->DimY)))/m_Y->DimY);
                    j = idx - (k*(m_Y->DimX*m_Y->DimY)) - i*m_Y->DimY;

                    xCenter = ((double) i)/(m_sf->Xfactor);

                    //get the surrounding pixels on x-dimension
                    xDown = (unsigned long) xCenter;
                    xUp = xDown + 1;
                    dx = xUp - xDown;
                    //compute linear interpolation
                    if(xUp==m_X->DimX)
                    {
                        //should handle boundary conditions on last YZ plane
                    }
                    else
                    {
                        //get the actual x-position of the pixel in former signal
                        yCenter = ((double) j)/(m_sf->Yfactor);

                        //get the surrounding pixels on x-dimension
                        yDown = (unsigned long) yCenter;
                        yUp = yDown + 1;
                        dy = yUp - yDown;
                        if(yUp==m_X->DimY)
                        {
                            //should handle boundary conditions on last XZ plane
                        }
                        else
                        {
                            //get the actual x-position of the pixel in former signal
                            zCenter = ((double) k)/(m_sf->Zfactor);

                            //get the surrounding pixels on x-dimension
                            zDown = (unsigned long) zCenter;
                            zUp = zDown + 1;
                            dz = zUp - zDown;

                            if(zUp==m_X->DimZ)
                            {
                                //should handle boundary conditions on last XY plane
                            }
                            else
                            {
                                ///x y z
                                w000 = (xUp - xCenter)*(yUp - yCenter)*(zUp-zCenter); //down dwon down
                                w100 = (xCenter - xDown)*(yUp - yCenter)*(zUp-zCenter); //up down down
                                w010 = (xUp - xCenter)*(yCenter - yDown)*(zUp-zCenter); //down up down
                                w110 = (xCenter - xDown)*(yCenter - yDown)*(zUp-zCenter); //up up down
                                w001 = (xUp - xCenter)*(yUp - yCenter)*(zCenter-zDown); //down dwon up
                                w101 = (xCenter - xDown)*(yUp - yCenter)*(zCenter-zDown); //up down up
                                w011 = (xUp - xCenter)*(yCenter - yDown)*(zCenter-zDown); //down up up
                                w111 = (xCenter - xDown)*(yCenter - yDown)*(zCenter-zDown); //up up up
                                m_Y->raw3D[i][j][k] = (1.0/(dx*dy*dz))*(w000*m_X->raw3D[xDown][yDown][zDown] +\
                                                                    w100*m_X->raw3D[xUp][yDown][zDown] +\
                                                                     w010*m_X->raw3D[xDown][yUp][zDown] +\
                                                                      w110*m_X->raw3D[xUp][yUp][zDown] +\
                                                                       w001*m_X->raw3D[xDown][yDown][zUp] +\
                                                                        w101*m_X->raw3D[xUp][yDown][zUp] +\
                                                                         w011*m_X->raw3D[xDown][yUp][zUp] +\
                                                                          w111*m_X->raw3D[xUp][yUp][zUp]);
                            }
                        }
                    }
                }
            }
            if(m_X->numDim==4)
            {
                long i, j, k, l;
                ///DimY*i + j + DimX*DimY*k + DimX*DimY*DimZ*l
                unsigned long xUp, xDown, yUp, yDown, zUp, zDown, tUp, tDown;
                double xCenter, yCenter, zCenter, tCenter, w0000, w1000, w1100, w0100, w0010, w1010, w1110, w0110, w0001, w1001, w1101, w0101, w0011, w1011, w1111, w0111, dx, dy, dz, dt;
                for(long idx=m_idxStart ; idx<m_idxEnd ; idx++)
                {
                    l = floor(idx/(m_Y->DimX*m_Y->DimY*m_Y->DimZ));
                    k = floor((idx - l*(m_Y->DimX*m_Y->DimY*m_Y->DimZ))/(m_Y->DimX*m_Y->DimY));
                    i = floor((idx - l*(m_Y->DimX*m_Y->DimY*m_Y->DimZ) - k*(m_Y->DimX*m_Y->DimY))/m_Y->DimY);
                    j = idx - l*(m_Y->DimX*m_Y->DimY*m_Y->DimZ) - k*(m_Y->DimX*m_Y->DimY) - i*m_Y->DimY;

                    //get the actual x-position of the pixel in former signal
                    xCenter = ((double) i)/(m_sf->Xfactor);

                    //get the surrounding pixels on x-dimension
                    xDown = (unsigned long) xCenter;
                    xUp = xDown + 1;
                    dx = xUp - xDown;
                    //compute linear interpolation
                    if(xUp==m_X->DimX)
                    {
                        //should handle boundary conditions on last YZT hyperplane
                    }
                    else
                    {
                        //get the actual x-position of the pixel in former signal
                        yCenter = ((double) j)/(m_sf->Yfactor);

                        //get the surrounding pixels on x-dimension
                        yDown = (unsigned long) yCenter;
                        yUp = yDown + 1;
                        dy = yUp - yDown;
                        if(yUp==m_X->DimY)
                        {
                            //should handle boundary conditions on last XZT hyperplane
                        }
                        else
                        {
                            //get the actual x-position of the pixel in former signal
                            zCenter = ((double) k)/(m_sf->Zfactor);

                            //get the surrounding pixels on x-dimension
                            zDown = (unsigned long) zCenter;
                            zUp = zDown + 1;
                            dz = zUp - zDown;

                            if(zUp==m_X->DimZ)
                            {
                                //should handle boundary conditions on last XYT hyperplane
                            }
                            else
                            {
                                //get the actual x-position of the pixel in former signal
                                tCenter = ((double) l)/(m_sf->Tfactor);

                                //get the surrounding pixels on x-dimension
                                tDown = (unsigned long) tCenter;
                                tUp = tDown + 1;
                                dt = tUp - tDown;

                                if(tUp==m_X->DimT)
                                {
                                    //should handle boundary conditions on last XYZ hyperplane
                                }
                                else
                                {
                                    ///x y z t
                                    w0000 = (xUp - xCenter)*(yUp - yCenter)*(zUp-zCenter)*(tUp-tCenter); //down dwon down down
                                    w1000 = (xCenter - xDown)*(yUp - yCenter)*(zUp-zCenter)*(tUp-tCenter); //up down down down
                                    w0100 = (xUp - xCenter)*(yCenter - yDown)*(zUp-zCenter)*(tUp-tCenter); //down up down down
                                    w1100 = (xCenter - xDown)*(yCenter - yDown)*(zUp-zCenter)*(tUp-tCenter); //up up down down
                                    w0010 = (xUp - xCenter)*(yUp - yCenter)*(zCenter-zDown)*(tUp-tCenter); //down dwon up down
                                    w1010 = (xCenter - xDown)*(yUp - yCenter)*(zCenter-zDown)*(tUp-tCenter); //up down up down
                                    w0110 = (xUp - xCenter)*(yCenter - yDown)*(zCenter-zDown)*(tUp-tCenter); //down up up down
                                    w1110 = (xCenter - xDown)*(yCenter - yDown)*(zCenter-zDown)*(tUp-tCenter); //up up up down
                                    w0001 = (xUp - xCenter)*(yUp - yCenter)*(zUp-zCenter)*(tCenter-tDown); //down dwon down up
                                    w1001 = (xCenter - xDown)*(yUp - yCenter)*(zUp-zCenter)*(tCenter-tDown); //up down down up
                                    w0101 = (xUp - xCenter)*(yCenter - yDown)*(zUp-zCenter)*(tCenter-tDown); //down up down up
                                    w1101 = (xCenter - xDown)*(yCenter - yDown)*(zUp-zCenter)*(tCenter-tDown); //up up down up
                                    w0011 = (xUp - xCenter)*(yUp - yCenter)*(zCenter-zDown)*(tCenter-tDown); //down dwon up up
                                    w1011 = (xCenter - xDown)*(yUp - yCenter)*(zCenter-zDown)*(tCenter-tDown); //up down up up
                                    w0111 = (xUp - xCenter)*(yCenter - yDown)*(zCenter-zDown)*(tCenter-tDown); //down up up up
                                    w1111 = (xCenter - xDown)*(yCenter - yDown)*(zCenter-zDown)*(tCenter-tDown); //up up up up
                                    m_Y->raw4D[i][j][k][l] = (1.0/(dx*dy*dz*dt))*(w0000*m_X->raw4D[xDown][yDown][zDown][tDown] +\
                                                                        w1000*m_X->raw4D[xUp][yDown][zDown][tDown] +\
                                                                         w0100*m_X->raw4D[xDown][yUp][zDown][tDown] +\
                                                                          w1100*m_X->raw4D[xUp][yUp][zDown][tDown] +\
                                                                           w0010*m_X->raw4D[xDown][yDown][zUp][tDown] +\
                                                                            w1010*m_X->raw4D[xUp][yDown][zUp][tDown] +\
                                                                             w0110*m_X->raw4D[xDown][yUp][zUp][tDown] +\
                                                                              w1110*m_X->raw4D[xUp][yUp][zUp][tDown] +\
                                                                               w0001*m_X->raw4D[xDown][yDown][zDown][tUp] +\
                                                                                w1001*m_X->raw4D[xUp][yDown][zDown][tUp] +\
                                                                                 w0101*m_X->raw4D[xDown][yUp][zDown][tUp] +\
                                                                                  w1101*m_X->raw4D[xUp][yUp][zDown][tUp] +\
                                                                                   w0011*m_X->raw4D[xDown][yDown][zUp][tUp] +\
                                                                                    w1011*m_X->raw4D[xUp][yDown][zUp][tUp] +\
                                                                                     w0111*m_X->raw4D[xDown][yUp][zUp][tUp] +\
                                                                                      w1111*m_X->raw4D[xUp][yUp][zUp][tUp]);
                                }
                            }
                        }
                    }
                }
            }
        };
    protected:
    private:
        long m_idxStart;
        long m_idxEnd;
        rawData<T>* m_X;
        samplingFactor* m_sf;
        rawData<T>* m_Y;
};







template<class T> class C_thread_gauss_trick : public C_thread
{
    friend class C_sampling<T>;
    public:
        /** Default constructor */
        C_thread_gauss_trick(rawData<T>* X, samplingFactor* sf, rawData<T>* Y, long idxStart, long idxEnd)
        {
            if(idxStart<=idxEnd)
            {
                m_idxStart = idxStart;
                m_idxEnd = idxEnd;//,  A->numel());
            }
            else
            {
                m_idxStart = idxEnd;//,  A->numel());
                m_idxEnd = idxStart;
            }
            m_X = X;
            m_Y = Y;
            m_sf = sf;
        };

        /** Default destructor */
        virtual ~C_thread_gauss_trick()
        {
            //
        };
        virtual void execute()
        {
            if(m_X->numDim==1)
            {
                unsigned long xUp, xDown;
                double xCenter;
                for(long idx=m_idxStart ; idx<m_idxEnd ; idx++)
                {

                    xCenter = ((double) idx)/(m_sf->Xfactor);
                    if(idx==0)
                    {
                        xDown=xCenter=0;
                        xUp = (m_sf->Rx) + 1;
                    }
                    else
                    {
                        xDown = (unsigned long) max( 0.0, xCenter - (double) (m_sf->Rx) );
                        xUp = (unsigned long) min((double) m_X->DimX, xCenter + (double) (m_sf->Rx) + 1.0 );
                        xCenter = 0.5*(double) (xDown+xUp);
                    }
                    double Z = 0.0, temp = 0.0, dX = 0.0;
                    for(unsigned long x=xDown ; x<xUp ; x++)
                    {
                        dX = (((double) x)-xCenter)*m_X->pixDimX;
                        temp += m_X->raw1D[x]*exp(-SQR(dX)/SQR(m_sf->sigx));
                        Z += exp(-SQR(dX)/SQR(m_sf->sigx));
                    }
                    m_Y->raw1D[idx] = (T) (temp/Z);
                }
            }
            if(m_X->numDim==2)
            {
                long i, j;
                ///idx = i*DimY + j;
                unsigned long xUp, xDown, yUp, yDown;
                double xCenter, yCenter;//, dx, dy;
                for(long idx=m_idxStart ; idx<m_idxEnd ; idx++)
                {
                    i = floor(idx/m_Y->DimY);
                    j = idx - (floor(idx/m_Y->DimY)*m_Y->DimY);

                    if(i==0)
                    {
                        xDown=xCenter=0;
                        xUp = (m_sf->Rx) +1;
                    }
                    else
                    {
                        xCenter = ((double) i)/(m_sf->Xfactor);//0.5*(double) (xDown+xUp);
                        xDown = (unsigned long) max(0.0, xCenter - (double) (m_sf->Rx) );
                        xUp = (unsigned long) min((double) m_X->DimX,  xCenter + (double) (m_sf->Rx) + 1.0);

                    }
                    if(j==0)
                    {
                        yDown=yCenter=0;
                        yUp = (unsigned long) (m_sf->Ry) +1;
                    }
                    else
                    {
                        yCenter = (((double) j)/(m_sf->Yfactor));//0.5*(double) (yDown+yUp);
                        yDown = (unsigned long) max( 0.0, yCenter - (double) (m_sf->Ry) );
                        yUp = (unsigned long) min((double) m_X->DimY,  yCenter + (double) (m_sf->Ry) + 1.0);
                    }
                    double Z = 0.0, temp = 0.0, tempZ, dX, dY;
                    for(unsigned long x=xDown ; x<xUp ; x++)
                    {
                        dX = (((double) x)-xCenter)*m_X->pixDimX;
                        for(unsigned long y=yDown ; y<yUp ; y++)
                        {
                            dY = (((double) y)-yCenter)*m_X->pixDimY;
                            tempZ = exp(-SQR(dX)/SQR(m_sf->sigx))*exp(-SQR(dY)/SQR(m_sf->sigy));
                            temp += m_X->raw2D[x][y]*tempZ;
                            Z += tempZ;
                        }
                    }
                    m_Y->raw2D[i][j] = (T) (temp/Z);
                }
            }
            if(m_X->numDim==3)
            {
                long i, j, k;
                ///DimY*i + j + DimX*DimY*k
                unsigned long xUp, xDown, yUp, yDown, zUp, zDown;
                double xCenter, yCenter, zCenter;
                for(long idx=m_idxStart ; idx<m_idxEnd ; idx++)
                {
                    k = floor(idx/(m_Y->DimX*m_Y->DimY));
                    i = floor((idx - (k*(m_Y->DimX*m_Y->DimY)))/m_Y->DimY);
                    j = idx - (k*(m_Y->DimX*m_Y->DimY)) - i*m_Y->DimY;

                    if(i==0)
                    {
                        xDown=xCenter=0;
                        xUp = (m_sf->Rx) +1;
                    }
                    else
                    {
                        xCenter = ((double) i)/(m_sf->Xfactor);//0.5*(double) (xDown+xUp);
                        xDown = (unsigned long) max(0.0, xCenter - (double) (m_sf->Rx) );
                        xUp = (unsigned long) min((double) m_X->DimX,  xCenter + (double) (m_sf->Rx) + 1.0);

                    }
                    if(j==0)
                    {
                        yDown=yCenter=0;
                        yUp = (unsigned long) (m_sf->Ry) +1;
                    }
                    else
                    {
                        yCenter = (((double) j)/(m_sf->Yfactor));//0.5*(double) (yDown+yUp);
                        yDown = (unsigned long) max( 0.0, yCenter - (double) (m_sf->Ry) );
                        yUp = (unsigned long) min((double) m_X->DimY,  yCenter + (double) (m_sf->Ry) + 1.0);
                    }
                    if(k==0)
                    {
                        zDown=zCenter=0;
                        zUp = (unsigned long) (m_sf->Rz) +1;
                    }
                    else
                    {
                        zCenter = (((double) k)/(m_sf->Zfactor));//0.5*(double) (yDown+yUp);
                        zDown = (unsigned long) max( 0.0, zCenter - (double) (m_sf->Rz) );
                        zUp = (unsigned long) min((double) m_X->DimZ,  zCenter + (double) (m_sf->Rz) + 1.0);
                    }
                    double Z = 0.0, temp = 0.0, tempZ, dX, dY, dZ;
                    for(unsigned long x=xDown ; x<xUp ; x++)
                    {
                        dX = (((double) x)-xCenter)*m_X->pixDimX;
                        for(unsigned long y=yDown ; y<yUp ; y++)
                        {
                            dY = (((double) y)-yCenter)*m_X->pixDimY;
                            for(unsigned long z=zDown ; z<zUp ; z++)
                            {
                                dZ = (((double) z)-zCenter)*m_X->pixDimZ;
                                tempZ = exp(-SQR(dX)/SQR(m_sf->sigx))*exp(-SQR(dY)/SQR(m_sf->sigy))*exp(-SQR(dZ)/SQR(m_sf->sigz));
                                temp += m_X->raw3D[x][y][z]*tempZ;
                                Z += tempZ;
                            }
                        }
                    }
                    m_Y->raw3D[i][j][k] = (T) (temp/Z);
                }
            }
            if(m_X->numDim==4)
            {
                long i, j, k, l;
                ///DimY*i + j + DimX*DimY*k + DimX*DimY*DimZ*l
                unsigned long xUp, xDown, yUp, yDown, zUp, zDown, tUp, tDown;
                double xCenter, yCenter, zCenter, tCenter;
                for(long idx=m_idxStart ; idx<m_idxEnd ; idx++)
                {
                    l = floor(idx/(m_Y->DimX*m_Y->DimY*m_Y->DimZ));
                    k = floor((idx - l*(m_Y->DimX*m_Y->DimY*m_Y->DimZ))/(m_Y->DimX*m_Y->DimY));
                    i = floor((idx - l*(m_Y->DimX*m_Y->DimY*m_Y->DimZ) - k*(m_Y->DimX*m_Y->DimY))/m_Y->DimY);
                    j = idx - l*(m_Y->DimX*m_Y->DimY*m_Y->DimZ) - k*(m_Y->DimX*m_Y->DimY) - i*m_Y->DimY;

                    if(i==0)
                    {
                        xDown=xCenter=0;
                        xUp = (m_sf->Rx) +1;
                    }
                    else
                    {
                        xCenter = ((double) i)/(m_sf->Xfactor);//0.5*(double) (xDown+xUp);
                        xDown = (unsigned long) max(0.0, xCenter - (double) (m_sf->Rx) );
                        xUp = (unsigned long) min((double) m_X->DimX,  xCenter + (double) (m_sf->Rx) + 1.0);

                    }
                    if(j==0)
                    {
                        yDown=yCenter=0;
                        yUp = (unsigned long) (m_sf->Ry) +1;
                    }
                    else
                    {
                        yCenter = (((double) j)/(m_sf->Yfactor));//0.5*(double) (yDown+yUp);
                        yDown = (unsigned long) max( 0.0, yCenter - (double) (m_sf->Ry) );
                        yUp = (unsigned long) min((double) m_X->DimY,  yCenter + (double) (m_sf->Ry) + 1.0);
                    }
                    if(k==0)
                    {
                        zDown=zCenter=0;
                        zUp = (unsigned long) (m_sf->Rz) +1;
                    }
                    else
                    {
                        zCenter = (((double) k)/(m_sf->Zfactor));//0.5*(double) (yDown+yUp);
                        zDown = (unsigned long) max( 0.0, zCenter - (double) (m_sf->Rz) );
                        zUp = (unsigned long) min((double) m_X->DimZ,  zCenter + (double) (m_sf->Rz) + 1.0);
                    }
                    if(l==0)
                    {
                        tDown=tCenter=0;
                        tUp = (unsigned long) (m_sf->Rt) +1;
                    }
                    else
                    {
                        tCenter = (((double) l)/(m_sf->Tfactor));//0.5*(double) (yDown+yUp);
                        tDown = (unsigned long) max( 0.0, tCenter - (double) (m_sf->Rt) );
                        tUp = (unsigned long) min((double) m_X->DimT,  tCenter + (double) (m_sf->Rt) + 1.0);
                    }
                    double Z = 0.0, temp = 0.0, tempZ, dX, dY, dZ, dT;
                    for(unsigned long x=xDown ; x<xUp ; x++)
                    {
                        dX = (((double) x)-xCenter)*m_X->pixDimX;
                        for(unsigned long y=yDown ; y<yUp ; y++)
                        {
                            dY = (((double) y)-yCenter)*m_X->pixDimY;
                            for(unsigned long z=zDown ; z<zUp ; z++)
                            {
                                dZ = (((double) z)-zCenter)*m_X->pixDimZ;
                                for(unsigned long t=tDown ; t<tUp ; t++)
                                {
                                    dT = (((double) t)-tCenter)*m_X->pixDimT;
                                    tempZ = exp(-SQR(dX)/SQR(m_sf->sigx))*exp(-SQR(dY)/SQR(m_sf->sigy))*exp(-SQR(dZ)/SQR(m_sf->sigz))*exp(-SQR(dT)/SQR(m_sf->sigt));
                                    temp += m_X->raw4D[x][y][z][t]*tempZ;
                                    Z += tempZ;
                                }
                            }
                        }
                    }
                    m_Y->raw4D[i][j][k][l] = (T) (temp/Z);
                }
            }
        };
    protected:
    private:
        long m_idxStart;
        long m_idxEnd;
        rawData<T>* m_X;
        samplingFactor* m_sf;
        rawData<T>* m_Y;
};


template<class T> class C_thread_gauss_conv3 : public C_thread
{
    friend class C_sampling<T>;
    public:
        /** Default constructor */
        C_thread_gauss_conv3(rawData<T>* X, rawData<T>* Y, long idxStart, long idxEnd)
        {
            if(idxStart<=idxEnd)
            {
                m_idxStart = idxStart;
                m_idxEnd = idxEnd;//,  A->numel());
            }
            else
            {
                m_idxStart = idxEnd;//,  A->numel());
                m_idxEnd = idxStart;
            }
            m_X = X;
            m_Y = Y;
            h = new double**[3];
            for(short i=0 ; i<3 ; i++)
            {
                h[i] = new double*[3];
                for(short j=0 ; j<3 ; j++)
                {
                    h[i][j] = new double[3];
                }
            }

            h[0][0][0] = 1.0/64.0; h[0][1][0] = 2.0/64.0; h[0][2][0] = 1.0/64.0;
            h[1][0][0] = 2.0/64.0; h[1][1][0] = 4.0/64.0; h[1][2][0] = 2.0/64.0;
            h[2][0][0] = 1.0/64.0; h[2][1][0] = 2.0/64.0; h[2][2][0] = 1.0/64.0;

            h[0][0][1] = 2.0/64.0; h[0][1][1] = 4.0/64.0; h[0][2][1] = 2.0/64.0;
            h[1][0][1] = 4.0/64.0; h[1][1][1] = 8.0/64.0; h[1][2][1] = 4.0/64.0;
            h[2][0][1] = 2.0/64.0; h[2][1][1] = 4.0/64.0; h[2][2][1] = 2.0/64.0;

            h[0][0][2] = 1.0/64.0; h[0][1][2] = 2.0/64.0; h[0][2][2] = 1.0/64.0;
            h[1][0][2] = 2.0/64.0; h[1][1][2] = 4.0/64.0; h[1][2][2] = 2.0/64.0;
            h[2][0][2] = 1.0/64.0; h[2][1][2] = 2.0/64.0; h[2][2][2] = 1.0/64.0;

        };

        /** Default destructor */
        virtual ~C_thread_gauss_conv3()
        {
            for(short i=0 ; i<3 ; i++)
            {
                for(short j=0 ; j<3 ; j++)
                {
                    delete h[i][j];
                }
                delete h[i];
            }
            delete h;
        };
        virtual void execute()
        {
            if(m_X->numDim==3)
            {
                long i, j, k;
                for(long idx=m_idxStart ; idx<m_idxEnd ; idx++)
                {
                    k = floor(idx/(m_Y->DimX*m_Y->DimY));
                    i = floor((idx - (k*(m_Y->DimX*m_Y->DimY)))/m_Y->DimY);
                    j = idx - (k*(m_Y->DimX*m_Y->DimY)) - i*m_Y->DimY;

                    if(i==0 || i==((long) m_Y->DimX-1) )
                    {
                        if(j==0 || j==((long) m_Y->DimY-1) )
                        {
                            if(k==0 || k==((long) m_Y->DimZ-1))
                            {
                                //
                            }
                            else
                            {
                                //
                            }
                        }
                        else
                        {
                            if(k==0 || k==((long) m_Y->DimZ-1) )
                            {
                                //
                            }
                            else
                            {
                                //
                            }
                        }
                    }
                    else if(j==0 || j==((long) m_Y->DimY-1))
                    {
                        if(k==0 || k==((long) m_Y->DimZ-1))
                        {
                            //
                        }
                        else
                        {
                            //
                        }
                    }
                    else if(k==0 || k==((long) m_Y->DimZ-1))
                    {
                        //
                    }
                    else
                    {
                        ///usual case
                        for(long di=-1 ; di<=1 ; di++)
                        {
                            for(long dj=-1 ; dj<=1 ; dj++)
                            {
                                for(long dk=-1 ; dk<=1 ; dk++)
                                {
//                                    if(i<0) cout << "i<0" << endl;
//                                    if(j<0) cout << "j<0" << endl;
//                                    if(k<0) cout << "k<0" << endl;
//                                    if(i>=m_Y->DimX) cout << "i<DimX" << endl;
//                                    if(j>=m_Y->DimY) cout << "j<DimY" << endl;
//                                    if(k>=m_Y->DimZ) cout << "k<DimZ" << endl;
//
//                                    if(i+di<0) cout << "i+di<0" << endl;
//                                    if(j+dj<0) cout << "j+dj<0" << endl;
//                                    if(k+dk<0) cout << "k+dk<0" << endl;
//                                    if(i+di>=m_Y->DimX) cout << "i+di>=DimX" << endl;
//                                    if(j+dj>=m_Y->DimY) cout << "j+dj>=DimY" << endl;
//                                    if(k+dk>=m_Y->DimZ) cout << "k+dk>=DimZ" << endl;
                                    m_Y->raw3D[i][j][k] += m_X->raw3D[i+di][j+dj][k+dk]*h[1-di][1-dj][1-dk];
                                }
                            }
                        }

                    }
                }
            }
            if(m_X->numDim==4)
            {
                long i, j, k, l;
                for(long idx=m_idxStart ; idx<m_idxEnd ; idx++)
                {
                    l = floor(idx/(m_Y->DimX*m_Y->DimY*m_Y->DimZ));
                    k = floor((idx - l*(m_Y->DimX*m_Y->DimY*m_Y->DimZ))/(m_Y->DimX*m_Y->DimY));
                    i = floor((idx - l*(m_Y->DimX*m_Y->DimY*m_Y->DimZ) - k*(m_Y->DimX*m_Y->DimY))/m_Y->DimY);
                    j = idx - l*(m_Y->DimX*m_Y->DimY*m_Y->DimZ) - k*(m_Y->DimX*m_Y->DimY) - i*m_Y->DimY;

                    if(i==0 || i==((long) m_Y->DimX-1) )
                    {
                        if(j==0 || j==((long) m_Y->DimY-1) )
                        {
                            if(k==0 || k==((long) m_Y->DimZ-1) )
                            {
                                //
                            }
                            else
                            {
                                //
                            }
                        }
                        else
                        {
                            if(k==0 || k==((long) m_Y->DimZ-1) )
                            {
                                //
                            }
                            else
                            {
                                //
                            }
                        }
                    }
                    else if(j==0 || j==((long) m_Y->DimY-1) )
                    {
                        if(k==0 || k==((long) m_Y->DimZ-1) )
                        {
                            //
                        }
                        else
                        {
                            //
                        }
                    }
                    else if(k==0 || k==((long) m_Y->DimZ-1) )
                    {
                        //
                    }
                    else
                    {
                        ///usual case
                        for(long di=-1 ; di<=1 ; di++)
                        {
                            for(long dj=-1 ; dj<=1 ; dj++)
                            {
                                for(long dk=-1 ; dk<=1 ; dk++)
                                {
//                                    if(i<0) cout << "i<0" << endl;
//                                    if(j<0) cout << "j<0" << endl;
//                                    if(k<0) cout << "k<0" << endl;
//                                    if(i>=m_Y->DimX) cout << "i<DimX" << endl;
//                                    if(j>=m_Y->DimY) cout << "j<DimY" << endl;
//                                    if(k>=m_Y->DimZ) cout << "k<DimZ" << endl;
//
//                                    if(i+di<0) cout << "i+di<0" << endl;
//                                    if(j+dj<0) cout << "j+dj<0" << endl;
//                                    if(k+dk<0) cout << "k+dk<0" << endl;
//                                    if(i+di>=m_Y->DimX) cout << "i+di>=DimX" << endl;
//                                    if(j+dj>=m_Y->DimY) cout << "j+dj>=DimY" << endl;
//                                    if(k+dk>=m_Y->DimZ) cout << "k+dk>=DimZ" << endl;

                                    m_Y->raw4D[i][j][k][l] += m_X->raw4D[i+di][j+dj][k+dk][l]*h[1-di][1-dj][1-dk];
                                }
                            }
                        }

                    }
                }
            }
        };
    private:
        long m_idxStart;
        long m_idxEnd;
        rawData<T>* m_X;
        rawData<T>* m_Y;
        double*** h;
};




template<class T> class C_sampling  //should be multithreaded
{
    public:
        /** Default constructor */
        C_sampling();
        /** Default destructor */
        virtual ~C_sampling();
        rawData<T>* samplingNearestNeighbor(rawData<T>* X, samplingFactor* sf);
        rawData<T>* samplingLinearInterpolation(rawData<T>* X, samplingFactor* sf); ///uses a current-neighborhood to calculate current value "low pass effect"
        //rawData<T>* samplingCubicInterpolation(rawData<T>* X, samplingFactor* sf); ///uses a current-neighborhood to calculate current value "sharpening effect"
        rawData<T>* samplingGaussTrick(rawData<T>* X, samplingFactor* sf); ///gaussian function depends on samplingFactor (kernel trick-like)

        ///gaussian pyramidal downsampling
        rawData<T>* gaussPyramidStep3D(rawData<T>* X);
};

template<class T> C_sampling<T>::C_sampling()
{
    //ctor
}

template<class T> C_sampling<T>::~C_sampling()
{
    //dtor
}


template<class T> rawData<T>* C_sampling<T>::samplingNearestNeighbor(rawData<T>* X, samplingFactor* sf)
{
    unsigned long newDimX, newDimY, newDimZ, newDimT;
    newDimX = ((unsigned long) ((double) X->DimX)*(sf->Xfactor));
    newDimY = ((unsigned long) ((double) X->DimY)*(sf->Yfactor));
    newDimZ = ((unsigned long) ((double) X->DimZ)*(sf->Zfactor));
    newDimT = ((unsigned long) ((double) X->DimT)*(sf->Tfactor));
    rawData<T>* Y = new rawData<T>(X->DICOM_TYPE, X->numDim, newDimX, newDimY, newDimZ, newDimT);
    *Y = 0;
    Y->pixDimX = X->pixDimX/(sf->Xfactor);
    Y->pixDimY = X->pixDimY/(sf->Yfactor);
    Y->pixDimZ = X->pixDimZ/(sf->Zfactor);
    Y->pixDimT = X->pixDimT/(sf->Tfactor);
    unsigned long NB_THREAD = 8;
    C_thread_NN<T>** THREAD_EQ = new C_thread_NN<T>*[NB_THREAD];
    for(unsigned long i=0 ; i<NB_THREAD ; i++)
    {
        unsigned long idxStart = (i*(Y->numel()))/NB_THREAD;
        unsigned long idxEnd = (((i+1)*(Y->numel()))/NB_THREAD);//-1;
        THREAD_EQ[i] = new C_thread_NN<T>(X, sf, Y, idxStart, idxEnd);
    }
    for(unsigned long i=0 ; i<NB_THREAD ; i++)
    {
        //cout << "start thread " << i << endl;
        THREAD_EQ[i]->start();
    }
    for(unsigned long i=0 ; i<NB_THREAD ; i++)
    {
        //cout << "wait thread " << i << endl;
        THREAD_EQ[i]->wait_for_exit();
    }
    for(unsigned long i=0 ; i<NB_THREAD ; i++)
    {
        delete THREAD_EQ[i];
    }
    delete THREAD_EQ;
    return Y;
}
template<class T> rawData<T>* C_sampling<T>::samplingLinearInterpolation(rawData<T>* X, samplingFactor* sf)
{
    unsigned long newDimX, newDimY, newDimZ, newDimT;
    newDimX = ((unsigned long) ((double) X->DimX)*(sf->Xfactor));
    newDimY = ((unsigned long) ((double) X->DimY)*(sf->Yfactor));
    newDimZ = ((unsigned long) ((double) X->DimZ)*(sf->Zfactor));
    newDimT = ((unsigned long) ((double) X->DimT)*(sf->Tfactor));
    rawData<T>* Y = new rawData<T>(X->DICOM_TYPE, X->numDim, newDimX, newDimY, newDimZ, newDimT);
    *Y = 0;
    Y->pixDimX = X->pixDimX/(sf->Xfactor);
    Y->pixDimY = X->pixDimY/(sf->Yfactor);
    Y->pixDimZ = X->pixDimZ/(sf->Zfactor);
    Y->pixDimT = X->pixDimT/(sf->Tfactor);
    unsigned long NB_THREAD = 8;
    C_thread_NN<T>** THREAD_EQ = new C_thread_NN<T>*[NB_THREAD];
    for(unsigned long i=0 ; i<NB_THREAD ; i++)
    {
        unsigned long idxStart = (i*(Y->numel()))/NB_THREAD;
        unsigned long idxEnd = (((i+1)*(Y->numel()))/NB_THREAD);//-1;
        THREAD_EQ[i] = new C_thread_NN<T>(X, sf, Y, idxStart, idxEnd);
    }
    for(unsigned long i=0 ; i<NB_THREAD ; i++)
    {
        THREAD_EQ[i]->start();
    }
    for(unsigned long i=0 ; i<NB_THREAD ; i++)
    {
        THREAD_EQ[i]->wait_for_exit();
    }
    for(unsigned long i=0 ; i<NB_THREAD ; i++)
    {
        delete THREAD_EQ[i];
    }
    delete THREAD_EQ;
    return Y;
}
template<class T> rawData<T>* C_sampling<T>::samplingGaussTrick(rawData<T>* X, samplingFactor* sf)
{
    unsigned long newDimX, newDimY, newDimZ, newDimT;
    newDimX = ((unsigned long) ((double) X->DimX)*(sf->Xfactor));
    newDimY = ((unsigned long) ((double) X->DimY)*(sf->Yfactor));
    newDimZ = ((unsigned long) ((double) X->DimZ)*(sf->Zfactor));
    newDimT = ((unsigned long) ((double) X->DimT)*(sf->Tfactor));
    rawData<T>* Y = new rawData<T>(X->DICOM_TYPE, X->numDim, newDimX, newDimY, newDimZ, newDimT);
    *Y = 0;
    Y->pixDimX = X->pixDimX/(sf->Xfactor);
    Y->pixDimY = X->pixDimY/(sf->Yfactor);
    Y->pixDimZ = X->pixDimZ/(sf->Zfactor);
    Y->pixDimT = X->pixDimT/(sf->Tfactor);
    unsigned long NB_THREAD = 8;
    C_thread_gauss_trick<T>** THREAD_EQ = new C_thread_gauss_trick<T>*[NB_THREAD];
    for(unsigned long i=0 ; i<NB_THREAD ; i++)
    {
        unsigned long idxStart = (i*(Y->numel()))/NB_THREAD;
        unsigned long idxEnd = (((i+1)*(Y->numel()))/NB_THREAD);//-1;
        THREAD_EQ[i] = new C_thread_gauss_trick<T>(X, sf, Y, idxStart, idxEnd);
    }
    for(unsigned long i=0 ; i<NB_THREAD ; i++)
    {
        THREAD_EQ[i]->start();
    }
    for(unsigned long i=0 ; i<NB_THREAD ; i++)
    {
        THREAD_EQ[i]->wait_for_exit();
    }
    for(unsigned long i=0 ; i<NB_THREAD ; i++)
    {
        delete THREAD_EQ[i];
    }
    delete THREAD_EQ;
    return Y;
}

template<class T> rawData<T>* C_sampling<T>::gaussPyramidStep3D(rawData<T>* X)
{
    if(X->numDim==1 || X->numDim==2) return NULL;
    rawData<T>* Y = new rawData<T>(X->DICOM_TYPE, X->numDim, X->DimX, X->DimY, X->DimZ, X->DimT);
    *Y = *X;
    Y->pixDimX = X->pixDimX;
    Y->pixDimY = X->pixDimY;
    Y->pixDimZ = X->pixDimZ;
    Y->pixDimT = X->pixDimT;
    unsigned long NB_THREAD = 8;
    C_thread_gauss_conv3<T>** THREAD_EQ = new C_thread_gauss_conv3<T>*[NB_THREAD];
    for(unsigned long i=0 ; i<NB_THREAD ; i++)
    {
        unsigned long idxStart = (i*(Y->numel()))/NB_THREAD;
        unsigned long idxEnd = (((i+1)*(Y->numel()))/NB_THREAD);//-1;
        THREAD_EQ[i] = new C_thread_gauss_conv3<T>(X, Y, idxStart, idxEnd);
    }
    for(unsigned long i=0 ; i<NB_THREAD ; i++)
    {
        THREAD_EQ[i]->start();
    }
    for(unsigned long i=0 ; i<NB_THREAD ; i++)
    {
        THREAD_EQ[i]->wait_for_exit();
    }
    for(unsigned long i=0 ; i<NB_THREAD ; i++)
    {
        delete THREAD_EQ[i];
    }
    delete THREAD_EQ;

//    return Y;

    //downsample
//    cout << "I am here" << endl;
    rawData<T>* Z = new rawData<T>(Y->DICOM_TYPE, Y->numDim, Y->DimX/2, Y->DimY/2, Y->DimZ/2, Y->DimT);
//    cout << "I am here?" << endl;
    *Z = 0;
    Z->pixDimX = 2*Y->pixDimX;
    Z->pixDimY = 2*Y->pixDimY;
    Z->pixDimZ = 2*Y->pixDimZ;
    Z->pixDimT = Y->pixDimT;
    if(Z->numDim==3)
    {
        for(long i=0 ; i<(long) Z->DimX ; i++)
        {
            for(long j=0 ; j<(long) Z->DimY ; j++)
            {
                for(long k=0 ; k<(long) Z->DimZ ; k++)
                {
                    Z->raw3D[i][j][k] = Y->raw3D[2*i][2*j][2*k];
                }
            }
        }
    }
    else if(Z->numDim==4)
    {
        for(long i=0 ; i<(long) Z->DimX ; i++)
        {
            for(long j=0 ; j<(long) Z->DimY ; j++)
            {
                for(long k=0 ; k<(long) Z->DimZ ; k++)
                {
                    for(long l=0 ; l<(long) Z->DimT ; l++)
                    {
                        Z->raw4D[i][j][k][l] = Y->raw4D[2*i][2*j][2*k][l];
                    }
                }
            }
        }
    }
    return Z;
}

#endif // C_SAMPLING_H
