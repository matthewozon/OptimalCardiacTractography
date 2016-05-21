#ifndef RAWDATA_H
#define RAWDATA_H

#include <C_thread.h>
#include <stdlib.h>
#include <math.h>
#include <struct.h>


template<class pixelType> class rawData;

template<class pixelType> class C_thread_equal : public C_thread
{
    friend class rawData<pixelType>;
    public:
        /** Default constructor */
        C_thread_equal(rawData<pixelType>* A, const pixelType b, long idxStart, long idxEnd)
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
            m_A = A;
            m_B = NULL;
            m_b = b;
            is_lvalue = true;
        };
        C_thread_equal(rawData<pixelType>* A, const rawData<pixelType>* B, long idxStart, long idxEnd)
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
            m_A = A;
            m_B = (rawData<pixelType>*) B;
            m_b = 0.0;
            is_lvalue = false;
        };

        /** Default destructor */
        virtual ~C_thread_equal()
        {
            //
        };
        virtual void execute()
        {
            if(is_lvalue)
            {
                if(m_A->numDim==1)
                {
                    for(long idx=m_idxStart ; idx<m_idxEnd ; idx++)
                    {
                        m_A->raw1D[idx] = m_b;
                    }
                }
                if(m_A->numDim==2)
                {
                    unsigned long i, j;
                    ///idx = i*DimY + j;
                    for(long idx=m_idxStart ; idx<m_idxEnd ; idx++)
                    {
                        i = m_A->getIndexI2D(idx);
                        j = m_A->getIndexJ2D(idx);
                        m_A->raw2D[i][j] = m_b;
                    }
                }
                if(m_A->numDim==3)
                {
                    unsigned long i, j, k;
                    ///DimY*i + j + DimX*DimY*k
                    for(long idx=m_idxStart ; idx<m_idxEnd ; idx++)
                    {
                        k = m_A->getIndexK3D(idx);
                        i = m_A->getIndexI3D(idx, k);
                        j = m_A->getIndexJ3D(idx, i, k);
                        m_A->raw3D[i][j][k] = m_b;
                    }
                }
                if(m_A->numDim==4)
                {
                    unsigned long i, j, k, l;
                    ///DimY*i + j + DimX*DimY*k + DimX*DimY*DimZ*l
                    for(long idx=m_idxStart ; idx<m_idxEnd ; idx++)
                    {
                        l = m_A->getIndexL4D(idx);
                        k = m_A->getIndexK4D(idx, l);
                        i = m_A->getIndexI4D(idx, k, l);
                        j = m_A->getIndexJ4D(idx, i, k, l);
                        m_A->raw4D[i][j][k][l] = m_b;
                    }
                }
            }
            else
            {
                if(m_A->numDim==1)
                {
                    for(long idx=m_idxStart ; idx<m_idxEnd ; idx++)
                    {
                        m_A->raw1D[idx] = m_B->raw1D[idx];
                    }
                }
                if(m_A->numDim==2)
                {
                    unsigned long i, j;
                    for(long idx=m_idxStart ; idx<m_idxEnd ; idx++)
                    {
                        i = m_A->getIndexI2D(idx);
                        j = m_A->getIndexJ2D(idx);
                        m_A->raw2D[i][j] = m_B->raw2D[i][j];
                    }
                }
                if(m_A->numDim==3)
                {
                    unsigned long i, j, k;
                    for(long idx=m_idxStart ; idx<m_idxEnd ; idx++)
                    {
                        k = m_A->getIndexK3D(idx);
                        i = m_A->getIndexI3D(idx, k);
                        j = m_A->getIndexJ3D(idx, i, k);
                        m_A->raw3D[i][j][k] = m_B->raw3D[i][j][k];
                    }
                }
                if(m_A->numDim==4)
                {
                    unsigned long i, j, k, l;
                    for(long idx=m_idxStart ; idx<m_idxEnd ; idx++)
                    {
                        l = m_A->getIndexL4D(idx);
                        k = m_A->getIndexK4D(idx, l);
                        i = m_A->getIndexI4D(idx, k, l);
                        j = m_A->getIndexJ4D(idx, i, k, l);
                        m_A->raw4D[i][j][k][l] = m_B->raw4D[i][j][k][l];
                    }
                }
                m_A->pixDimX = m_B->pixDimX;
                m_A->pixDimY = m_B->pixDimY;
                m_A->pixDimZ = m_B->pixDimZ;
                m_A->pixDimT = m_B->pixDimT;
            }
        };
    protected:
    private:
        long m_idxStart;
        long m_idxEnd;
        rawData<pixelType>* m_A;
        rawData<pixelType>* m_B;
        pixelType m_b;
        bool is_lvalue;
};

template<class pixelType> class C_thread_add : public C_thread
{
    friend class rawData<pixelType>;
    public:
        /** Default constructor */
        C_thread_add(rawData<pixelType>* A, const pixelType b, rawData<pixelType>* res, long idxStart, long idxEnd)
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
            m_A = A;
            m_B = NULL;
            m_b = b;
            m_R = res;
            is_lvalue = true;
        };
        C_thread_add(rawData<pixelType>* A, const rawData<pixelType>* B, rawData<pixelType>* res, long idxStart, long idxEnd)
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
            m_A = A;
            m_B = (rawData<pixelType>*) B;
            m_b = 0.0;
            m_R = res;
            is_lvalue = false;
        };

        /** Default destructor */
        virtual ~C_thread_add()
        {
            //
        };
        virtual void execute()
        {
            m_R->pixDimX = m_A->pixDimX;
            m_R->pixDimY = m_A->pixDimY;
            m_R->pixDimZ = m_A->pixDimZ;
            m_R->pixDimT = m_A->pixDimT;
            if(is_lvalue)
            {
                if(m_A->numDim==1)
                {
                    for(long idx=m_idxStart ; idx<m_idxEnd ; idx++)
                    {
                        m_R->raw1D[idx] = m_A->raw1D[idx] + m_b;
                    }
                }
                if(m_A->numDim==2)
                {
                    unsigned long i, j;
                    ///idx = i*DimY + j;
                    for(long idx=m_idxStart ; idx<m_idxEnd ; idx++)
                    {
                        i = m_A->getIndexI2D(idx);
                        j = m_A->getIndexJ2D(idx);
                        m_R->raw2D[i][j] = m_A->raw2D[i][j] + m_b;
                    }
                }
                if(m_A->numDim==3)
                {
                    unsigned long i, j, k;
                    ///DimY*i + j + DimX*DimY*k
                    for(long idx=m_idxStart ; idx<m_idxEnd ; idx++)
                    {
                        k = m_A->getIndexK3D(idx);
                        i = m_A->getIndexI3D(idx, k);
                        j = m_A->getIndexJ3D(idx, i, k);
                        m_R->raw3D[i][j][k] = m_A->raw3D[i][j][k] + m_b;
                    }
                }
                if(m_A->numDim==4)
                {
                    unsigned long i, j, k, l;
                    ///DimY*i + j + DimX*DimY*k + DimX*DimY*DimZ*l
                    for(long idx=m_idxStart ; idx<m_idxEnd ; idx++)
                    {
                        l = m_A->getIndexL4D(idx);
                        k = m_A->getIndexK4D(idx, l);
                        i = m_A->getIndexI4D(idx, k, l);
                        j = m_A->getIndexJ4D(idx, i, k, l);
                        m_R->raw4D[i][j][k][l] = m_A->raw4D[i][j][k][l] + m_b;
                    }
                }
            }
            else
            {
                if(m_A->numDim==1)
                {
                    for(long idx=m_idxStart ; idx<m_idxEnd ; idx++)
                    {
                        m_R->raw1D[idx] = m_A->raw1D[idx] + m_B->raw1D[idx];
                    }
                }
                if(m_A->numDim==2)
                {
                    unsigned long i, j;
                    for(long idx=m_idxStart ; idx<m_idxEnd ; idx++)
                    {
                        i = m_A->getIndexI2D(idx);
                        j = m_A->getIndexJ2D(idx);
                        m_R->raw2D[i][j] = m_A->raw2D[i][j] + m_B->raw2D[i][j];
                    }
                }
                if(m_A->numDim==3)
                {
                    unsigned long i, j, k;
                    for(long idx=m_idxStart ; idx<m_idxEnd ; idx++)
                    {
                        k = m_A->getIndexK3D(idx);
                        i = m_A->getIndexI3D(idx, k);
                        j = m_A->getIndexJ3D(idx, i, k);
                        m_R->raw3D[i][j][k] = m_A->raw3D[i][j][k] + m_B->raw3D[i][j][k];
                    }
                }
                if(m_A->numDim==4)
                {
                    unsigned long i, j, k, l;
                    for(long idx=m_idxStart ; idx<m_idxEnd ; idx++)
                    {
                        l = m_A->getIndexL4D(idx);
                        k = m_A->getIndexK4D(idx, l);
                        i = m_A->getIndexI4D(idx, k, l);
                        j = m_A->getIndexJ4D(idx, i, k, l);
                        m_R->raw4D[i][j][k][l] = m_A->raw4D[i][j][k][l] + m_B->raw4D[i][j][k][l];
                    }
                }
            }
        };
    protected:
    private:
        long m_idxStart;
        long m_idxEnd;
        rawData<pixelType>* m_A;
        rawData<pixelType>* m_B;
        rawData<pixelType>* m_R;
        pixelType m_b;
        bool is_lvalue;
};

template<class pixelType> class C_thread_sub : public C_thread
{
    friend class rawData<pixelType>;
    public:
        /** Default constructor */
        C_thread_sub(rawData<pixelType>* A, const pixelType b, rawData<pixelType>* res, long idxStart, long idxEnd)
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
            m_A = A;
            m_B = NULL;
            m_b = b;
            m_R = res;
            is_lvalue = true;
        };
        C_thread_sub(rawData<pixelType>* A, const rawData<pixelType>* B, rawData<pixelType>* res, long idxStart, long idxEnd)
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
            m_A = A;
            m_B = (rawData<pixelType>*) B;
            m_b = 0.0;
            m_R = res;
            is_lvalue = false;
        };

        /** Default destructor */
        virtual ~C_thread_sub()
        {
            //
        };
        virtual void execute()
        {
            m_R->pixDimX = m_A->pixDimX;
            m_R->pixDimY = m_A->pixDimY;
            m_R->pixDimZ = m_A->pixDimZ;
            m_R->pixDimT = m_A->pixDimT;
            if(is_lvalue)
            {
                if(m_A->numDim==1)
                {
                    for(unsigned long idx=m_idxStart ; idx<m_idxEnd ; idx++)
                    {
                        m_R->raw1D[idx] = m_A->raw1D[idx] - m_b;
                    }
                }
                if(m_A->numDim==2)
                {
                    unsigned long i, j;
                    ///idx = i*DimY + j;
                    for(unsigned long idx=m_idxStart ; idx<m_idxEnd ; idx++)
                    {
                        i = m_A->getIndexI2D(idx);
                        j = m_A->getIndexJ2D(idx);
                        m_R->raw2D[i][j] = m_A->raw2D[i][j] - m_b;
                    }
                }
                if(m_A->numDim==3)
                {
                    unsigned long i, j, k;
                    ///DimY*i + j + DimX*DimY*k
                    for(unsigned long idx=m_idxStart ; idx<m_idxEnd ; idx++)
                    {
                        k = m_A->getIndexK3D(idx);
                        i = m_A->getIndexI3D(idx, k);
                        j = m_A->getIndexJ3D(idx, i, k);
                        m_R->raw3D[i][j][k] = m_A->raw3D[i][j][k] - m_b;
                    }
                }
                if(m_A->numDim==4)
                {
                    unsigned long i, j, k, l;
                    ///DimY*i + j + DimX*DimY*k + DimX*DimY*DimZ*l
                    for(unsigned long idx=m_idxStart ; idx<m_idxEnd ; idx++)
                    {
                        l = m_A->getIndexL4D(idx);
                        k = m_A->getIndexK4D(idx, l);
                        i = m_A->getIndexI4D(idx, k, l);
                        j = m_A->getIndexJ4D(idx, i, k, l);
                        m_R->raw4D[i][j][k][l] = m_A->raw4D[i][j][k][l] - m_b;
                    }
                }
            }
            else
            {
                if(m_A->numDim==1)
                {
                    for(unsigned long idx=m_idxStart ; idx<m_idxEnd ; idx++)
                    {
                        m_R->raw1D[idx] = m_A->raw1D[idx] - m_B->raw1D[idx];
                    }
                }
                if(m_A->numDim==2)
                {
                    unsigned long i, j;
                    for(unsigned long idx=m_idxStart ; idx<m_idxEnd ; idx++)
                    {
                        i = m_A->getIndexI2D(idx);
                        j = m_A->getIndexJ2D(idx);
                        m_R->raw2D[i][j] = m_A->raw2D[i][j] - m_B->raw2D[i][j];
                    }
                }
                if(m_A->numDim==3)
                {
                    unsigned long i, j, k;
                    for(unsigned long idx=m_idxStart ; idx<m_idxEnd ; idx++)
                    {
                        k = m_A->getIndexK3D(idx);
                        i = m_A->getIndexI3D(idx, k);
                        j = m_A->getIndexJ3D(idx, i, k);
                        m_R->raw3D[i][j][k] = m_A->raw3D[i][j][k] - m_B->raw3D[i][j][k];
                    }
                }
                if(m_A->numDim==4)
                {
                    unsigned long i, j, k, l;
                    for(unsigned long idx=m_idxStart ; idx<m_idxEnd ; idx++)
                    {
                        l = m_A->getIndexL4D(idx);
                        k = m_A->getIndexK4D(idx, l);
                        i = m_A->getIndexI4D(idx, k, l);
                        j = m_A->getIndexJ4D(idx, i, k, l);
                        m_R->raw4D[i][j][k][l] = m_A->raw4D[i][j][k][l] - m_B->raw4D[i][j][k][l];
                    }
                }
            }
        };
    protected:
    private:
        long m_idxStart;
        long m_idxEnd;
        rawData<pixelType>* m_A;
        rawData<pixelType>* m_B;
        rawData<pixelType>* m_R;
        pixelType m_b;
        bool is_lvalue;
};

template<class pixelType> class C_thread_mult : public C_thread
{
    friend class rawData<pixelType>;
    public:
        /** Default constructor */
        C_thread_mult(rawData<pixelType>* A, const pixelType b, rawData<pixelType>* res, long idxStart, long idxEnd)
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
            m_A = A;
            m_B = NULL;
            m_b = b;
            m_R = res;
            is_lvalue = true;
        };
        C_thread_mult(rawData<pixelType>* A, const rawData<pixelType>* B, rawData<pixelType>* res, long idxStart, long idxEnd)
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
            m_A = A;
            m_B = (rawData<pixelType>*) B;
            m_b = 0.0;
            m_R = res;
            is_lvalue = false;
        };

        /** Default destructor */
        virtual ~C_thread_mult()
        {
            //
        };
        virtual void execute()
        {
            m_R->pixDimX = m_A->pixDimX;
            m_R->pixDimY = m_A->pixDimY;
            m_R->pixDimZ = m_A->pixDimZ;
            m_R->pixDimT = m_A->pixDimT;
            if(is_lvalue)
            {
                if(m_A->numDim==1)
                {
                    for(long idx=m_idxStart ; idx<m_idxEnd ; idx++)
                    {
                        m_R->raw1D[idx] = m_A->raw1D[idx] * m_b;
                    }
                }
                if(m_A->numDim==2)
                {
                    unsigned long i, j;
                    ///idx = i*DimY + j;
                    for(long idx=m_idxStart ; idx<m_idxEnd ; idx++)
                    {
                        i = m_A->getIndexI2D(idx);
                        j = m_A->getIndexJ2D(idx);
                        m_R->raw2D[i][j] = m_A->raw2D[i][j] * m_b;
                    }
                }
                if(m_A->numDim==3)
                {
                    unsigned long i, j, k;
                    ///DimY*i + j + DimX*DimY*k
                    for(long idx=m_idxStart ; idx<m_idxEnd ; idx++)
                    {
                        k = m_A->getIndexK3D(idx);
                        i = m_A->getIndexI3D(idx, k);
                        j = m_A->getIndexJ3D(idx, i, k);
                        m_R->raw3D[i][j][k] = m_A->raw3D[i][j][k] * m_b;
                    }
                }
                if(m_A->numDim==4)
                {
                    unsigned long i, j, k, l;
                    ///DimY*i + j + DimX*DimY*k + DimX*DimY*DimZ*l
                    for(long idx=m_idxStart ; idx<m_idxEnd ; idx++)
                    {
                        l = m_A->getIndexL4D(idx);
                        k = m_A->getIndexK4D(idx, l);
                        i = m_A->getIndexI4D(idx, k, l);
                        j = m_A->getIndexJ4D(idx, i, k, l);
                        m_R->raw4D[i][j][k][l] = m_A->raw4D[i][j][k][l] * m_b;
                    }
                }
            }
            else
            {
                if(m_A->numDim==1)
                {
                    for(long idx=m_idxStart ; idx<m_idxEnd ; idx++)
                    {
                        m_R->raw1D[idx] = m_A->raw1D[idx] * m_B->raw1D[idx];
                    }
                }
                if(m_A->numDim==2)
                {
                    unsigned long i, j;
                    for(long idx=m_idxStart ; idx<m_idxEnd ; idx++)
                    {
                        i = m_A->getIndexI2D(idx);
                        j = m_A->getIndexJ2D(idx);
                        m_R->raw2D[i][j] = m_A->raw2D[i][j] * m_B->raw2D[i][j];
                    }
                }
                if(m_A->numDim==3)
                {
                    unsigned long i, j, k;
                    for(long idx=m_idxStart ; idx<m_idxEnd ; idx++)
                    {
                        k = m_A->getIndexK3D(idx);
                        i = m_A->getIndexI3D(idx, k);
                        j = m_A->getIndexJ3D(idx, i, k);
                        m_R->raw3D[i][j][k] = m_A->raw3D[i][j][k] * m_B->raw3D[i][j][k];
                    }
                }
                if(m_A->numDim==4)
                {
                    unsigned long i, j, k, l;
                    for(long idx=m_idxStart ; idx<m_idxEnd ; idx++)
                    {
                        l = m_A->getIndexL4D(idx);
                        k = m_A->getIndexK4D(idx, l);
                        i = m_A->getIndexI4D(idx, k, l);
                        j = m_A->getIndexJ4D(idx, i, k, l);
                        m_R->raw4D[i][j][k][l] = m_A->raw4D[i][j][k][l] * m_B->raw4D[i][j][k][l];
                    }
                }
            }
        };
    protected:
    private:
        long m_idxStart;
        long m_idxEnd;
        rawData<pixelType>* m_A;
        rawData<pixelType>* m_B;
        rawData<pixelType>* m_R;
        pixelType m_b;
        bool is_lvalue;
};

template<class pixelType> class C_thread_div : public C_thread
{
    friend class rawData<pixelType>;
    public:
        /** Default constructor */
        C_thread_div(rawData<pixelType>* A, const pixelType b, rawData<pixelType>* res, long idxStart, long idxEnd)
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
            m_A = A;
            m_B = NULL;
            m_b = b;
            m_R = res;
            is_lvalue = true;
        };
        C_thread_div(rawData<pixelType>* A, const rawData<pixelType>* B, rawData<pixelType>* res, long idxStart, long idxEnd)
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
            m_A = A;
            m_B = (rawData<pixelType>*) B;
            m_b = 0.0;
            m_R = res;
            is_lvalue = false;
        };

        /** Default destructor */
        virtual ~C_thread_div()
        {
            //
        };
        virtual void execute()
        {
            m_R->pixDimX = m_A->pixDimX;
            m_R->pixDimY = m_A->pixDimY;
            m_R->pixDimZ = m_A->pixDimZ;
            m_R->pixDimT = m_A->pixDimT;
            if(is_lvalue)
            {
                if(m_A->numDim==1)
                {
                    for(long idx=m_idxStart ; idx<m_idxEnd ; idx++)
                    {
                        m_R->raw1D[idx] = m_A->raw1D[idx] / m_b;
                    }
                }
                if(m_A->numDim==2)
                {
                    unsigned long i, j;
                    ///idx = i*DimY + j;
                    for(long idx=m_idxStart ; idx<m_idxEnd ; idx++)
                    {
                        i = m_A->getIndexI2D(idx);
                        j = m_A->getIndexJ2D(idx);
                        m_R->raw2D[i][j] = m_A->raw2D[i][j] / m_b;
                    }
                }
                if(m_A->numDim==3)
                {
                    unsigned long i, j, k;
                    ///DimY*i + j + DimX*DimY*k
                    for(long idx=m_idxStart ; idx<m_idxEnd ; idx++)
                    {
                        k = m_A->getIndexK3D(idx);
                        i = m_A->getIndexI3D(idx, k);
                        j = m_A->getIndexJ3D(idx, i, k);
                        m_R->raw3D[i][j][k] = m_A->raw3D[i][j][k] / m_b;
                    }
                }
                if(m_A->numDim==4)
                {
                    unsigned long i, j, k, l;
                    ///DimY*i + j + DimX*DimY*k + DimX*DimY*DimZ*l
                    for(long idx=m_idxStart ; idx<m_idxEnd ; idx++)
                    {
                        l = m_A->getIndexL4D(idx);
                        k = m_A->getIndexK4D(idx, l);
                        i = m_A->getIndexI4D(idx, k, l);
                        j = m_A->getIndexJ4D(idx, i, k, l);
                        m_R->raw4D[i][j][k][l] = m_A->raw4D[i][j][k][l] / m_b;
                    }
                }
            }
            else
            {
                if(m_A->numDim==1)
                {
                    for(long idx=m_idxStart ; idx<m_idxEnd ; idx++)
                    {
                        m_R->raw1D[idx] = m_A->raw1D[idx] / m_B->raw1D[idx];
                    }
                }
                if(m_A->numDim==2)
                {
                    unsigned long i, j;
                    for(long idx=m_idxStart ; idx<m_idxEnd ; idx++)
                    {
                        i = m_A->getIndexI2D(idx);
                        j = m_A->getIndexJ2D(idx);
                        m_R->raw2D[i][j] = m_A->raw2D[i][j] / m_B->raw2D[i][j];
                    }
                }
                if(m_A->numDim==3)
                {
                    unsigned long i, j, k;
                    for(long idx=m_idxStart ; idx<m_idxEnd ; idx++)
                    {
                        k = m_A->getIndexK3D(idx);
                        i = m_A->getIndexI3D(idx, k);
                        j = m_A->getIndexJ3D(idx, i, k);
                        m_R->raw3D[i][j][k] = m_A->raw3D[i][j][k] / m_B->raw3D[i][j][k];
                    }
                }
                if(m_A->numDim==4)
                {
                    unsigned long i, j, k, l;
                    for(long idx=m_idxStart ; idx<m_idxEnd ; idx++)
                    {
                        l = m_A->getIndexL4D(idx);
                        k = m_A->getIndexK4D(idx, l);
                        i = m_A->getIndexI4D(idx, k, l);
                        j = m_A->getIndexJ4D(idx, i, k, l);
                        m_R->raw4D[i][j][k][l] = m_A->raw4D[i][j][k][l] / m_B->raw4D[i][j][k][l];
                    }
                }
            }
        };
    protected:
    private:
        long m_idxStart;
        long m_idxEnd;
        rawData<pixelType>* m_A;
        rawData<pixelType>* m_B;
        rawData<pixelType>* m_R;
        pixelType m_b;
        bool is_lvalue;
};

template<class pixelType> class C_thread_larger : public C_thread
{
    friend class rawData<pixelType>;
    public:
        /** Default constructor */
        C_thread_larger(rawData<pixelType>* A, const pixelType b, rawData<pixelType>* res, long idxStart, long idxEnd)
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
            m_A = A;
            m_B = NULL;
            m_b = b;
            m_R = res;
            is_lvalue = true;
        };
        C_thread_larger(rawData<pixelType>* A, const rawData<pixelType>* B, rawData<pixelType>* res, long idxStart, long idxEnd)
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
            m_A = A;
            m_B = (rawData<pixelType>*) B;
            m_b = 0.0;
            m_R = res;
            is_lvalue = false;
        };

        /** Default destructor */
        virtual ~C_thread_larger()
        {
            //
        };
        virtual void execute()
        {
            m_R->pixDimX = m_A->pixDimX;
            m_R->pixDimY = m_A->pixDimY;
            m_R->pixDimZ = m_A->pixDimZ;
            m_R->pixDimT = m_A->pixDimT;
            if(is_lvalue)
            {
                if(m_A->numDim==1)
                {
                    for(long idx=m_idxStart ; idx<m_idxEnd ; idx++)
                    {
                        m_R->raw1D[idx] = m_A->raw1D[idx] > m_b;
                    }
                }
                if(m_A->numDim==2)
                {
                    unsigned long i, j;
                    ///idx = i*DimY + j;
                    for(long idx=m_idxStart ; idx<m_idxEnd ; idx++)
                    {
                        i = m_A->getIndexI2D(idx);
                        j = m_A->getIndexJ2D(idx);
                        m_R->raw2D[i][j] = m_A->raw2D[i][j] > m_b;
                    }
                }
                if(m_A->numDim==3)
                {
                    unsigned long i, j, k;
                    ///DimY*i + j + DimX*DimY*k
                    for(long idx=m_idxStart ; idx<m_idxEnd ; idx++)
                    {
                        k = m_A->getIndexK3D(idx);
                        i = m_A->getIndexI3D(idx, k);
                        j = m_A->getIndexJ3D(idx, i, k);
                        m_R->raw3D[i][j][k] = m_A->raw3D[i][j][k] > m_b;
                    }
                }
                if(m_A->numDim==4)
                {
                    unsigned long i, j, k, l;
                    ///DimY*i + j + DimX*DimY*k + DimX*DimY*DimZ*l
                    for(long idx=m_idxStart ; idx<m_idxEnd ; idx++)
                    {
                        l = m_A->getIndexL4D(idx);
                        k = m_A->getIndexK4D(idx, l);
                        i = m_A->getIndexI4D(idx, k, l);
                        j = m_A->getIndexJ4D(idx, i, k, l);
                        m_R->raw4D[i][j][k][l] = m_A->raw4D[i][j][k][l] > m_b;
                    }
                }
            }
            else
            {
                if(m_A->numDim==1)
                {
                    for(long idx=m_idxStart ; idx<m_idxEnd ; idx++)
                    {
                        m_R->raw1D[idx] = m_A->raw1D[idx] > m_B->raw1D[idx];
                    }
                }
                if(m_A->numDim==2)
                {
                    unsigned long i, j;
                    for(long idx=m_idxStart ; idx<m_idxEnd ; idx++)
                    {
                        i = m_A->getIndexI2D(idx);
                        j = m_A->getIndexJ2D(idx);
                        m_R->raw2D[i][j] = m_A->raw2D[i][j] > m_B->raw2D[i][j];
                    }
                }
                if(m_A->numDim==3)
                {
                    unsigned long i, j, k;
                    for(long idx=m_idxStart ; idx<m_idxEnd ; idx++)
                    {
                        k = m_A->getIndexK3D(idx);
                        i = m_A->getIndexI3D(idx, k);
                        j = m_A->getIndexJ3D(idx, i, k);
                        m_R->raw3D[i][j][k] = m_A->raw3D[i][j][k] > m_B->raw3D[i][j][k];
                    }
                }
                if(m_A->numDim==4)
                {
                    unsigned long i, j, k, l;
                    for(long idx=m_idxStart ; idx<m_idxEnd ; idx++)
                    {
                        l = m_A->getIndexL4D(idx);
                        k = m_A->getIndexK4D(idx, l);
                        i = m_A->getIndexI4D(idx, k, l);
                        j = m_A->getIndexJ4D(idx, i, k, l);
                        m_R->raw4D[i][j][k][l] = m_A->raw4D[i][j][k][l] > m_B->raw4D[i][j][k][l];
                    }
                }
            }
        };
    protected:
    private:
        long m_idxStart;
        long m_idxEnd;
        rawData<pixelType>* m_A;
        rawData<pixelType>* m_B;
        rawData<pixelType>* m_R;
        pixelType m_b;
        bool is_lvalue;
};

template<class pixelType> class C_thread_larger_or_equal : public C_thread
{
    friend class rawData<pixelType>;
    public:
        /** Default constructor */
        C_thread_larger_or_equal(rawData<pixelType>* A, const pixelType b, rawData<pixelType>* res, long idxStart, long idxEnd)
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
            m_A = A;
            m_B = NULL;
            m_b = b;
            m_R = res;
            is_lvalue = true;
        };
        C_thread_larger_or_equal(rawData<pixelType>* A, const rawData<pixelType>* B, rawData<pixelType>* res, long idxStart, long idxEnd)
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
            m_A = A;
            m_B = (rawData<pixelType>*) B;
            m_b = 0.0;
            m_R = res;
            is_lvalue = false;
        };

        /** Default destructor */
        virtual ~C_thread_larger_or_equal()
        {
            //
        };
        virtual void execute()
        {
            m_R->pixDimX = m_A->pixDimX;
            m_R->pixDimY = m_A->pixDimY;
            m_R->pixDimZ = m_A->pixDimZ;
            m_R->pixDimT = m_A->pixDimT;
            if(is_lvalue)
            {
                if(m_A->numDim==1)
                {
                    for(unsigned long idx=m_idxStart ; idx<m_idxEnd ; idx++)
                    {
                        m_R->raw1D[idx] = m_A->raw1D[idx] >= m_b;
                    }
                }
                if(m_A->numDim==2)
                {
                    unsigned long i, j;
                    ///idx = i*DimY + j;
                    for(unsigned long idx=m_idxStart ; idx<m_idxEnd ; idx++)
                    {
                        i = m_A->getIndexI2D(idx);
                        j = m_A->getIndexJ2D(idx);
                        m_R->raw2D[i][j] = m_A->raw2D[i][j] >= m_b;
                    }
                }
                if(m_A->numDim==3)
                {
                    unsigned long i, j, k;
                    ///DimY*i + j + DimX*DimY*k
                    for(unsigned long idx=m_idxStart ; idx<m_idxEnd ; idx++)
                    {
                        k = m_A->getIndexK3D(idx);
                        i = m_A->getIndexI3D(idx, k);
                        j = m_A->getIndexJ3D(idx, i, k);
                        m_R->raw3D[i][j][k] = m_A->raw3D[i][j][k] >= m_b;
                    }
                }
                if(m_A->numDim==4)
                {
                    unsigned long i, j, k, l;
                    ///DimY*i + j + DimX*DimY*k + DimX*DimY*DimZ*l
                    for(unsigned long idx=m_idxStart ; idx<m_idxEnd ; idx++)
                    {
                        l = m_A->getIndexL4D(idx);
                        k = m_A->getIndexK4D(idx, l);
                        i = m_A->getIndexI4D(idx, k, l);
                        j = m_A->getIndexJ4D(idx, i, k, l);
                        m_R->raw4D[i][j][k][l] = m_A->raw4D[i][j][k][l] >= m_b;
                    }
                }
            }
            else
            {
                if(m_A->numDim==1)
                {
                    for(unsigned long idx=m_idxStart ; idx<m_idxEnd ; idx++)
                    {
                        m_R->raw1D[idx] = m_A->raw1D[idx] >= m_B->raw1D[idx];
                    }
                }
                if(m_A->numDim==2)
                {
                    unsigned long i, j;
                    for(unsigned long idx=m_idxStart ; idx<m_idxEnd ; idx++)
                    {
                        i = m_A->getIndexI2D(idx);
                        j = m_A->getIndexJ2D(idx);
                        m_R->raw2D[i][j] = m_A->raw2D[i][j] >= m_B->raw2D[i][j];
                    }
                }
                if(m_A->numDim==3)
                {
                    unsigned long i, j, k;
                    for(unsigned long idx=m_idxStart ; idx<m_idxEnd ; idx++)
                    {
                        k = m_A->getIndexK3D(idx);
                        i = m_A->getIndexI3D(idx, k);
                        j = m_A->getIndexJ3D(idx, i, k);
                        m_R->raw3D[i][j][k] = m_A->raw3D[i][j][k] >= m_B->raw3D[i][j][k];
                    }
                }
                if(m_A->numDim==4)
                {
                    unsigned long i, j, k, l;
                    for(unsigned long idx=m_idxStart ; idx<m_idxEnd ; idx++)
                    {
                        l = m_A->getIndexL4D(idx);
                        k = m_A->getIndexK4D(idx, l);
                        i = m_A->getIndexI4D(idx, k, l);
                        j = m_A->getIndexJ4D(idx, i, k, l);
                        m_R->raw4D[i][j][k][l] = m_A->raw4D[i][j][k][l] >= m_B->raw4D[i][j][k][l];
                    }
                }
            }
        };
    protected:
    private:
        long m_idxStart;
        long m_idxEnd;
        rawData<pixelType>* m_A;
        rawData<pixelType>* m_B;
        rawData<pixelType>* m_R;
        pixelType m_b;
        bool is_lvalue;
};

template<class pixelType> class C_thread_smaller : public C_thread
{
    friend class rawData<pixelType>;
    public:
        /** Default constructor */
        C_thread_smaller(rawData<pixelType>* A, const pixelType b, rawData<pixelType>* res, long idxStart, long idxEnd)
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
            m_A = A;
            m_B = NULL;
            m_b = b;
            m_R = res;
            is_lvalue = true;
        };
        C_thread_smaller(rawData<pixelType>* A, const rawData<pixelType>* B, rawData<pixelType>* res, long idxStart, long idxEnd)
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
            m_A = A;
            m_B = (rawData<pixelType>*) B;
            m_b = 0.0;
            m_R = res;
            is_lvalue = false;
        };

        /** Default destructor */
        virtual ~C_thread_smaller()
        {
            //
        };
        virtual void execute()
        {
            m_R->pixDimX = m_A->pixDimX;
            m_R->pixDimY = m_A->pixDimY;
            m_R->pixDimZ = m_A->pixDimZ;
            m_R->pixDimT = m_A->pixDimT;
            if(is_lvalue)
            {
                if(m_A->numDim==1)
                {
                    for(long idx=m_idxStart ; idx<m_idxEnd ; idx++)
                    {
                        m_R->raw1D[idx] = m_A->raw1D[idx] < m_b;
                    }
                }
                if(m_A->numDim==2)
                {
                    unsigned long i, j;
                    ///idx = i*DimY + j;
                    for(long idx=m_idxStart ; idx<m_idxEnd ; idx++)
                    {
                        i = m_A->getIndexI2D(idx);
                        j = m_A->getIndexJ2D(idx);
                        m_R->raw2D[i][j] = m_A->raw2D[i][j] < m_b;
                    }
                }
                if(m_A->numDim==3)
                {
                    unsigned long i, j, k;
                    ///DimY*i + j + DimX*DimY*k
                    for(long idx=m_idxStart ; idx<m_idxEnd ; idx++)
                    {
                        k = m_A->getIndexK3D(idx);
                        i = m_A->getIndexI3D(idx, k);
                        j = m_A->getIndexJ3D(idx, i, k);
                        m_R->raw3D[i][j][k] = m_A->raw3D[i][j][k] < m_b;
                    }
                }
                if(m_A->numDim==4)
                {
                    unsigned long i, j, k, l;
                    ///DimY*i + j + DimX*DimY*k + DimX*DimY*DimZ*l
                    for(long idx=m_idxStart ; idx<m_idxEnd ; idx++)
                    {
                        l = m_A->getIndexL4D(idx);
                        k = m_A->getIndexK4D(idx, l);
                        i = m_A->getIndexI4D(idx, k, l);
                        j = m_A->getIndexJ4D(idx, i, k, l);
                        m_R->raw4D[i][j][k][l] = m_A->raw4D[i][j][k][l] < m_b;
                    }
                }
            }
            else
            {
                if(m_A->numDim==1)
                {
                    for(long idx=m_idxStart ; idx<m_idxEnd ; idx++)
                    {
                        m_R->raw1D[idx] = m_A->raw1D[idx] < m_B->raw1D[idx];
                    }
                }
                if(m_A->numDim==2)
                {
                    unsigned long i, j;
                    for(long idx=m_idxStart ; idx<m_idxEnd ; idx++)
                    {
                        i = m_A->getIndexI2D(idx);
                        j = m_A->getIndexJ2D(idx);
                        m_R->raw2D[i][j] = m_A->raw2D[i][j] < m_B->raw2D[i][j];
                    }
                }
                if(m_A->numDim==3)
                {
                    unsigned long i, j, k;
                    for(long idx=m_idxStart ; idx<m_idxEnd ; idx++)
                    {
                        k = m_A->getIndexK3D(idx);
                        i = m_A->getIndexI3D(idx, k);
                        j = m_A->getIndexJ3D(idx, i, k);
                        m_R->raw3D[i][j][k] = m_A->raw3D[i][j][k] < m_B->raw3D[i][j][k];
                    }
                }
                if(m_A->numDim==4)
                {
                    unsigned long i, j, k, l;
                    for(long idx=m_idxStart ; idx<m_idxEnd ; idx++)
                    {
                        l = m_A->getIndexL4D(idx);
                        k = m_A->getIndexK4D(idx, l);
                        i = m_A->getIndexI4D(idx, k, l);
                        j = m_A->getIndexJ4D(idx, i, k, l);
                        m_R->raw4D[i][j][k][l] = m_A->raw4D[i][j][k][l] < m_B->raw4D[i][j][k][l];
                    }
                }
            }
        };
    protected:
    private:
        long m_idxStart;
        long m_idxEnd;
        rawData<pixelType>* m_A;
        rawData<pixelType>* m_B;
        rawData<pixelType>* m_R;
        pixelType m_b;
        bool is_lvalue;
};

template<class pixelType> class C_thread_smaller_or_equal : public C_thread
{
    friend class rawData<pixelType>;
    public:
        /** Default constructor */
        C_thread_smaller_or_equal(rawData<pixelType>* A, const pixelType b, rawData<pixelType>* res, long idxStart, long idxEnd)
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
            m_A = A;
            m_B = NULL;
            m_b = b;
            m_R = res;
            is_lvalue = true;
        };
        C_thread_smaller_or_equal(rawData<pixelType>* A, const rawData<pixelType>* B, rawData<pixelType>* res, long idxStart, long idxEnd)
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
            m_A = A;
            m_B = (rawData<pixelType>*) B;
            m_b = 0.0;
            m_R = res;
            is_lvalue = false;
        };

        /** Default destructor */
        virtual ~C_thread_smaller_or_equal()
        {
            //
        };
        virtual void execute()
        {
            m_R->pixDimX = m_A->pixDimX;
            m_R->pixDimY = m_A->pixDimY;
            m_R->pixDimZ = m_A->pixDimZ;
            m_R->pixDimT = m_A->pixDimT;
            if(is_lvalue)
            {
                if(m_A->numDim==1)
                {
                    for(unsigned long idx=m_idxStart ; idx<m_idxEnd ; idx++)
                    {
                        m_R->raw1D[idx] = m_A->raw1D[idx] <= m_b;
                    }
                }
                if(m_A->numDim==2)
                {
                    unsigned long i, j;
                    ///idx = i*DimY + j;
                    for(unsigned long idx=m_idxStart ; idx<m_idxEnd ; idx++)
                    {
                        i = m_A->getIndexI2D(idx);
                        j = m_A->getIndexJ2D(idx);
                        m_R->raw2D[i][j] = m_A->raw2D[i][j] <= m_b;
                    }
                }
                if(m_A->numDim==3)
                {
                    unsigned long i, j, k;
                    ///DimY*i + j + DimX*DimY*k
                    for(unsigned long idx=m_idxStart ; idx<m_idxEnd ; idx++)
                    {
                        k = m_A->getIndexK3D(idx);
                        i = m_A->getIndexI3D(idx, k);
                        j = m_A->getIndexJ3D(idx, i, k);
                        m_R->raw3D[i][j][k] = m_A->raw3D[i][j][k] <= m_b;
                    }
                }
                if(m_A->numDim==4)
                {
                    unsigned long i, j, k, l;
                    ///DimY*i + j + DimX*DimY*k + DimX*DimY*DimZ*l
                    for(unsigned long idx=m_idxStart ; idx<m_idxEnd ; idx++)
                    {
                        l = m_A->getIndexL4D(idx);
                        k = m_A->getIndexK4D(idx, l);
                        i = m_A->getIndexI4D(idx, k, l);
                        j = m_A->getIndexJ4D(idx, i, k, l);
                        m_R->raw4D[i][j][k][l] = m_A->raw4D[i][j][k][l] <= m_b;
                    }
                }
            }
            else
            {
                if(m_A->numDim==1)
                {
                    for(unsigned long idx=m_idxStart ; idx<m_idxEnd ; idx++)
                    {
                        m_R->raw1D[idx] = m_A->raw1D[idx] <= m_B->raw1D[idx];
                    }
                }
                if(m_A->numDim==2)
                {
                    unsigned long i, j;
                    for(unsigned long idx=m_idxStart ; idx<m_idxEnd ; idx++)
                    {
                        i = m_A->getIndexI2D(idx);
                        j = m_A->getIndexJ2D(idx);
                        m_R->raw2D[i][j] = m_A->raw2D[i][j] <= m_B->raw2D[i][j];
                    }
                }
                if(m_A->numDim==3)
                {
                    unsigned long i, j, k;
                    for(unsigned long idx=m_idxStart ; idx<m_idxEnd ; idx++)
                    {
                        k = m_A->getIndexK3D(idx);
                        i = m_A->getIndexI3D(idx, k);
                        j = m_A->getIndexJ3D(idx, i, k);
                        m_R->raw3D[i][j][k] = m_A->raw3D[i][j][k] <= m_B->raw3D[i][j][k];
                    }
                }
                if(m_A->numDim==4)
                {
                    unsigned long i, j, k, l;
                    for(unsigned long idx=m_idxStart ; idx<m_idxEnd ; idx++)
                    {
                        l = m_A->getIndexL4D(idx);
                        k = m_A->getIndexK4D(idx, l);
                        i = m_A->getIndexI4D(idx, k, l);
                        j = m_A->getIndexJ4D(idx, i, k, l);
                        m_R->raw4D[i][j][k][l] = m_A->raw4D[i][j][k][l] <= m_B->raw4D[i][j][k][l];
                    }
                }
            }
        };
    protected:
    private:
        long m_idxStart;
        long m_idxEnd;
        rawData<pixelType>* m_A;
        rawData<pixelType>* m_B;
        rawData<pixelType>* m_R;
        pixelType m_b;
        bool is_lvalue;
};


template<class pixelType> class C_thread_sum : public C_thread
{
    friend class rawData<pixelType>;
    public:
        pixelType partialSum;
        /** Default constructor */
        C_thread_sum(rawData<pixelType>* A, long idxStart, long idxEnd)
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
            m_A = A;
        };

        /** Default destructor */
        virtual ~C_thread_sum()
        {
            //
        };
        virtual void execute()
        {
            partialSum = (pixelType) 0.0;
            if(m_A->numDim==1)
            {
                for(long idx=m_idxStart ; idx<m_idxEnd ; idx++)
                {
                    partialSum += m_A->raw1D[idx];
                }
            }
            if(m_A->numDim==2)
            {
                unsigned long i, j;
                ///idx = i*DimY + j;
                for(long idx=m_idxStart ; idx<m_idxEnd ; idx++)
                {
                    i = m_A->getIndexI2D(idx);
                    j = m_A->getIndexJ2D(idx);
                    partialSum += m_A->raw2D[i][j];
                }
            }
            if(m_A->numDim==3)
            {
                unsigned long i, j, k;
                ///DimY*i + j + DimX*DimY*k
                for(long idx=m_idxStart ; idx<m_idxEnd ; idx++)
                {
                    k = m_A->getIndexK3D(idx);
                    i = m_A->getIndexI3D(idx, k);
                    j = m_A->getIndexJ3D(idx, i, k);
                    partialSum += m_A->raw3D[i][j][k];
                }
            }
            if(m_A->numDim==4)
            {
                unsigned long i, j, k, l;
                ///DimY*i + j + DimX*DimY*k + DimX*DimY*DimZ*l
                for(long idx=m_idxStart ; idx<m_idxEnd ; idx++)
                {
                    l = m_A->getIndexL4D(idx);
                    k = m_A->getIndexK4D(idx, l);
                    i = m_A->getIndexI4D(idx, k, l);
                    j = m_A->getIndexJ4D(idx, i, k, l);
                    partialSum += m_A->raw4D[i][j][k][l];
                }
            }
        };
    protected:
    private:
        long m_idxStart;
        long m_idxEnd;
        rawData<pixelType>* m_A;
};


template<class pixelType> class C_thread_remove_isolated : public C_thread
{
    friend class rawData<pixelType>;
    public:
        /** Default constructor */
        C_thread_remove_isolated(rawData<pixelType>* A, long idxStart, long idxEnd)
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
            m_A = A;
        };

        /** Default destructor */
        virtual ~C_thread_remove_isolated()
        {
            //
        };
        virtual void execute()
        {
            if(m_A->numDim==1)
            {
                if(m_idxStart==0)
                {
                    if(m_A->raw1D[0]!=0 && m_A->raw1D[1]==0) m_A->raw1D[0] = 0;
                    m_idxStart++;
                }
                if(m_idxEnd==m_A->numel())
                {
                    if(m_A->raw1D[m_idxEnd-1]!=0 && m_A->raw1D[m_idxEnd-2]==0) m_A->raw1D[m_idxEnd-1] = 0;
                    m_idxEnd--;
                }
                for(long idx=m_idxStart ; idx<m_idxEnd ; idx++)
                {
                    if(m_A->raw1D[idx]!=0)
                    {
                        if(m_A->raw1D[idx-1]==0 && m_A->raw1D[idx+1]==0) m_A->raw1D[idx]=0;
                    }
                }
            }
            if(m_A->numDim==2)
            {
                long i, j, diUp, djUp, diDown, djDown;
                for(long idx=m_idxStart ; idx<m_idxEnd ; idx++)
                {
                    i = m_A->getIndexI2D(idx);
                    j = m_A->getIndexJ2D(idx);
                    if(m_A->raw2D[i][j]!=0)
                    {
                        bool cond = true;
                        if(i==0) diDown=0;
                        else diDown=-1;
                        if(j==0) djDown=0;
                        else djDown=-1;

                        if(i==m_A->DimX-1) diUp=0;
                        else diUp=1;
                        if(j==m_A->DimY-1) djUp=0;
                        else djUp=1;

                        for(long ii=i+diDown ; ii<=i+diUp ; ii++)
                        {
                            for(long jj=j+djDown ; jj<=j+djUp ; jj++)
                            {
                                if(!(ii==i && jj==j))
                                {
                                    if(m_A->raw2D[ii][jj]!=0)
                                    {
                                        cond = false;
                                        ii=i+diUp+1;
                                        jj=j+djUp+1;
                                    }
                                }
                            }
                        }
                        if(cond) m_A->raw2D[i][j]=0;
                    }
                }
            }
            if(m_A->numDim==3)
            {
                long i, j, k, diUp, djUp, dkUp, diDown, djDown, dkDown;
                for(long idx=m_idxStart ; idx<m_idxEnd ; idx++)
                {
                    k = m_A->getIndexK3D(idx);
                    i = m_A->getIndexI3D(idx, k);
                    j = m_A->getIndexJ3D(idx, i, k);
                    if(m_A->raw3D[i][j][k]!=0)
                    {
                        bool cond = true;
                        if(i==0) diDown=0;
                        else diDown=-1;
                        if(j==0) djDown=0;
                        else djDown=-1;
                        if(k==0) dkDown=0;
                        else dkDown=-1;

                        if(i==m_A->DimX-1) diUp=0;
                        else diUp=1;
                        if(j==m_A->DimY-1) djUp=0;
                        else djUp=1;
                        if(k==m_A->DimZ-1) dkUp=0;
                        else dkUp=1;

                        for(long ii=i+diDown ; ii<=i+diUp ; ii++)
                        {
                            for(long jj=j+djDown ; jj<=j+djUp ; jj++)
                            {
                                for(long kk=k+dkDown ; kk<=k+dkUp ; kk++)
                                {
                                    if(!(ii==i && jj==j && kk==k))
                                    {
                                        if(m_A->raw3D[ii][jj][kk]!=0.0)
                                        {
                                            cond = false;
                                            ii=i+diUp+1;
                                            jj=j+djUp+1;
                                            kk=k+dkUp+1;
                                        }
                                    }
                                }
                            }
                        }
                        if(cond) m_A->raw3D[i][j][k]=0;
                    }
                }
            }
            if(m_A->numDim==4)
            {
                long i, j, k, l, diUp, djUp, dkUp, diDown, djDown, dkDown, dlUp, dlDown;
                for(long idx=m_idxStart ; idx<m_idxEnd ; idx++)
                {
                    l = m_A->getIndexL4D(idx);
                    k = m_A->getIndexK4D(idx, l);
                    i = m_A->getIndexI4D(idx, k, l);
                    j = m_A->getIndexJ4D(idx, i, k, l);
                    if(m_A->raw4D[i][j][k][l]!=0)
                    {
                        bool cond = true;
                        if(i==0) diDown=0;
                        else diDown=-1;
                        if(j==0) djDown=0;
                        else djDown=-1;
                        if(k==0) dkDown=0;
                        else dkDown=-1;
                        if(l==0) dlDown=0;
                        else dlDown=-1;

                        if(i==m_A->DimX-1) diUp=0;
                        else diUp=1;
                        if(j==m_A->DimY-1) djUp=0;
                        else djUp=1;
                        if(k==m_A->DimZ-1) dkUp=0;
                        else dkUp=1;
                        if(l==m_A->DimT-1) dlUp=0;
                        else dlUp=1;

                        for(long ii=i+diDown ; ii<=i+diUp ; ii++)
                        {
                            for(long jj=j+djDown ; jj<=j+djUp ; jj++)
                            {
                                for(long kk=k+dkDown ; kk<=k+dkUp ; kk++)
                                {
                                    for(long ll=l+dlDown ; ll<=l+dlUp ; ll++)
                                    {
                                        if(!(ii==i && jj==j && kk==k && ll==l))
                                        {
                                            if(m_A->raw4D[ii][jj][kk][ll]!=0.0)
                                            {
                                                cond = false;
                                                ii=i+diUp+1;
                                                jj=j+djUp+1;
                                                kk=k+dkUp+1;
                                                ll=l+dlUp+1;
                                            }
                                        }
                                    }
                                }
                            }
                        }
                        if(cond) m_A->raw4D[i][j][k][l]=0;
                    }
                }
            }
        };
    private:
        long m_idxStart;
        long m_idxEnd;
        rawData<pixelType>* m_A;
};






template<class pixelType> class rawData
{
    public:
        /** Default constructor */
        rawData(unsigned short myType /*type of data*/, unsigned short dim, unsigned int* dims);
        rawData(unsigned short myType /*type of data*/, unsigned short dim, unsigned long l=1, unsigned long m=1, unsigned long n=1, unsigned long o=1);
        void rawData_(unsigned long myType /*type of data*/, unsigned long dim, unsigned long l=1, unsigned long m=1, unsigned long n=1, unsigned long o=1);
        /** Default destructor */
        virtual ~rawData();

        ///operand
        rawData& operator= (const rawData& c); //OK
        rawData& operator= (const pixelType& x); //OK
        rawData& operator+ (const rawData& c); //OK
        rawData& operator+ (const pixelType& c); //OK
        rawData& operator- (const rawData& c); //OK
        rawData& operator- (const pixelType& c); //OK
        rawData& operator* (const rawData& c); //OK
        rawData& operator* (const pixelType& x);  //OK
        rawData& operator/ (const rawData& c); //OK
        rawData& operator/ (const pixelType& x);  //OK
        rawData& operator+= (const rawData& c); //OK
        rawData& operator+= (const pixelType& c);  //OK
        rawData& operator-= (const rawData& c); //OK
        rawData& operator-= (const pixelType& c);  //OK
        rawData& operator*= (const rawData& c); //OK
        rawData& operator*= (const pixelType& c);  //OK
        rawData& operator/= (const rawData& c); //OK
        rawData& operator/= (const pixelType& c); //OK

        rawData& operator> (const rawData& c); //OK
        rawData& operator> (const pixelType& c); //OK
        rawData& operator>= (const rawData& c); //OK
        rawData& operator>= (const pixelType& c); //OK
        rawData& operator< (const rawData& c); //OK
        rawData& operator< (const pixelType& c); //OK
        rawData& operator<= (const rawData& c); //OK
        rawData& operator<= (const pixelType& c); //OK

        ///attributs
        unsigned long numDim, DimX, DimY, DimZ, DimT, DICOM_TYPE;
        pixelType* raw1D;
        pixelType** raw2D;
        pixelType*** raw3D;
        pixelType**** raw4D;

        ///methods
        unsigned long getIndex1D(unsigned long i);
        unsigned long getIndex2D(unsigned long i, unsigned long j);
        unsigned long getIndex2DUnmosaic(unsigned long i, unsigned long j, unsigned long nb_slice);
        unsigned long getIndex3DUnmosaic(unsigned long i, unsigned long j, unsigned long k);
        unsigned long getIndex3D(unsigned long i, unsigned long j, unsigned long k);
        unsigned long getIndex4D(unsigned long i, unsigned long j, unsigned long k, unsigned long l);

        unsigned long getIndexI1D(unsigned long idx);
        unsigned long getIndexI2D(unsigned long idx);
        unsigned long getIndexI3D(unsigned long idx, unsigned long k);
        unsigned long getIndexI4D(unsigned long idx, unsigned long k, unsigned long l);

        unsigned long getIndexJ2D(unsigned long idx);
        unsigned long getIndexJ3D(unsigned long idx, unsigned long i, unsigned long k);
        unsigned long getIndexJ4D(unsigned long idx, unsigned long i, unsigned long k, unsigned long l);

        unsigned long getIndexK3D(unsigned long idx);
        unsigned long getIndexK4D(unsigned long idx, unsigned long l);

        unsigned long getIndexL4D(unsigned long idx);

        pixelType getMax();
        pixelType getMin();
        unsigned long numel();
        void removeSingle();

        double pixDimX;
        double pixDimY;
        double pixDimZ;
        double pixDimT;

        //may implement save to file
        bool write1Dtxt(char* name){return false;}
        bool write2Dtxt(char* name){return false;}
        bool write3Dtxt(char* name){return false;}
        bool write4Dtxt(char* name){return false;}

        bool write1Ddcm(char* name){return false;}
        bool write2Ddcm(char* name){return false;}
        bool write3Ddcm(char* name){return false;}
        bool write4Ddcm(char* name){return false;}
        //implement get one slice ((n-1)D)
        //should implement crop

        //new
        pixelType getValue(unsigned long idxPix);
        bool setValue(unsigned long idxPix, pixelType newValue);
        pixelType getSum(void);
        pixelType getMean(void);
        pixelType getVar(void);
        void crop1D(unsigned long xmin, unsigned long xmax);
        void crop2D(unsigned long xmin, unsigned long xmax, unsigned long ymin, unsigned long ymax);
        void crop3D(unsigned long xmin, unsigned long xmax, unsigned long ymin, unsigned long ymax, unsigned long zmin, unsigned long zmax);
        void crop4D(unsigned long xmin, unsigned long xmax, unsigned long ymin, unsigned long ymax, unsigned long zmin, unsigned long zmax, unsigned long tmin, unsigned long tmax);
        void swapXY2D(void);
        void swapXY3D(void);
        void swapXY4D(void);
};



template<class pixelType> rawData<pixelType>::rawData(unsigned short myType /*type of data*/, unsigned short dim, unsigned int* dims) //hypp dims = [rows, columns, dpth...]
{
    switch(dim)
    {
        case 1:
        {
            rawData_(myType, 1, dims[0]);
            break;
        }
        case 2:
        {
            rawData_(myType, 2, dims[1], dims[0]);
            break;
        }
        case 3:
        {
            rawData_(myType, 3, dims[1], dims[0], dims[2]);
            break;
        }
        case 4:
        {
            rawData_(myType, 4, dims[1], dims[0], dims[2], dims[3]);
            break;
        }
        default:
        {
            break;
        }
    }
}
template<class pixelType> rawData<pixelType>::rawData(unsigned short myType /*type of data*/, unsigned short dim, unsigned long l, unsigned long m, unsigned long n, unsigned long o)
{
    rawData_(myType, dim, l, m, n, o);
    pixDimX = -1.0;
    pixDimY = -1.0;
    pixDimZ = -1.0;
    pixDimT = -1.0;
};

template<class pixelType> rawData<pixelType>::~rawData()
{
    if(raw1D!=NULL)
    {
        delete raw1D;
    }
    if(raw2D!=NULL)
    {
        for(unsigned short i=0 ; i<DimX ; i++)
        {
            if(raw2D[i]!=NULL) delete raw2D[i];
        }
        delete raw2D;
    }
    if(raw3D!=NULL)
    {
        for(unsigned short i=0 ; i<DimX ; i++)
        {
            if(raw3D[i]!=NULL)
            {
                for(unsigned short j=0 ; j<DimY ; j++)
                {
                    if(raw3D[i][j]!=NULL) delete raw3D[i][j];
                }
                delete raw3D[i];
            }

        }
        delete raw3D;
    }
    if(raw4D!=NULL)
    {
        for(unsigned short i=0 ; i<DimX ; i++)
        {
            if(raw4D[i]!=NULL)
            {
                for(unsigned short j=0 ; j<DimY ; j++)
                {
                    if(raw4D[i][j]!=NULL)
                    {
                        for(unsigned short k=0 ; k<DimZ ; k++)
                        {
                            if(raw4D[i][j][k]!=NULL) delete raw4D[i][j][k];
                        }
                        delete raw4D[i][j];
                    }
                }
                delete raw4D[i];
            }
        }
        delete raw4D;

    }
}

template<class pixelType>void rawData<pixelType>::rawData_(unsigned long myType /*type of data*/, unsigned long dim, unsigned long l, unsigned long m, unsigned long n, unsigned long o)
{
    DICOM_TYPE = myType;
    if(DICOM_TYPE!=DICOM_ERROR)
    {
        numDim = dim;
        switch(numDim)
        {
            case 1:
            {
               raw1D = new pixelType[l];
               raw2D = NULL;
               raw3D = NULL;
               raw4D = NULL;
               DimX = l;
               DimY = 1;
               DimZ = 1;
               DimT = 1;
               break;
            }
            case 2:
            {
               raw2D = new pixelType*[l];
               for(unsigned short i=0 ; i<l ; i++)
               {
                   raw2D[i] = new pixelType[m];
               }
               raw1D = NULL;
               raw3D = NULL;
               raw4D = NULL;
               DimX = l;
               DimY = m;
               DimZ = 1;
               DimT = 1;
               break;
            }
            case 3:
            {
               raw3D = new pixelType**[l];
               for(unsigned short i=0 ; i<l ; i++)
               {
                   raw3D[i] = new pixelType*[m];
                   for(unsigned short j=0 ; j<m ; j++)
                   {
                       raw3D[i][j] = new pixelType[n];
                   }
               }
               raw1D = NULL;
               raw2D = NULL;
               raw4D = NULL;
               DimX = l;
               DimY = m;
               DimZ = n;
               DimT = 1;
               break;
            }
            case 4:
            {
               raw4D = new pixelType***[l];
               for(unsigned short i=0 ; i<l ; i++)
               {
                   raw4D[i] = new pixelType**[m];
                   for(unsigned short j=0 ; j<m ; j++)
                   {
                       raw4D[i][j] = new pixelType*[n];
                       for(unsigned short k=0 ; k<n ; k++)
                       {
                           raw4D[i][j][k] = new pixelType[o];
                       }
                   }
               }
               raw1D = NULL;
               raw2D = NULL;
               raw3D = NULL;
               DimX = l;
               DimY = m;
               DimZ = n;
               DimT = o;
               break;
            }
            default:
            {
                raw1D = NULL;
                raw2D = NULL;
                raw3D = NULL;
                raw4D = NULL;
                DimX = 0;
                DimY = 0;
                DimZ = 0;
                DimT = 0;
                break;
            }
        }
    }
    else
    {
        raw1D = NULL;
        raw2D = NULL;
        raw3D = NULL;
        raw4D = NULL;
        DimX = 0;
        DimY = 0;
        DimZ = 0;
        DimT = 0;
        numDim = 0;
    }

}


template<class pixelType> unsigned long rawData<pixelType>::numel()
{
    if(numDim==1) return DimX;
    if(numDim==2) return DimX*DimY;
    if(numDim==3) return DimX*DimY*DimZ;
    if(numDim==4) return DimX*DimY*DimZ*DimT;
    return 0;
}

template<class pixelType> rawData<pixelType>& rawData<pixelType>::operator= (const pixelType& x) //rawData::
{
    unsigned long NB_THREAD = 8;
    C_thread_equal<pixelType>** THREAD_EQ = new C_thread_equal<pixelType>*[NB_THREAD];
    for(unsigned long i=0 ; i<NB_THREAD ; i++)
    {
        unsigned long idxStart = (i*(this->numel()))/NB_THREAD;
        unsigned long idxEnd = (((i+1)*(this->numel()))/NB_THREAD);//-1;
        THREAD_EQ[i] = new C_thread_equal<pixelType>(this, x, idxStart, idxEnd);
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
    return *this;

};
template<class pixelType> rawData<pixelType>& rawData<pixelType>::operator= (const rawData& c) //rawData::
{
    if(this->numDim==c.numDim)
    {
        if(this->numDim==0)
        {
            //do nothing
            return *this;
        }
        else if(this->numDim>4)
        {
            //do nothing
            return *this;
        }
        else if(this->DimX==c.DimX && this->DimY==c.DimY && this->DimZ==c.DimZ && this->DimT==c.DimT)
        {
            unsigned long NB_THREAD = 8;
            C_thread_equal<pixelType>** THREAD_EQ = new C_thread_equal<pixelType>*[NB_THREAD];
            for(unsigned long i=0 ; i<NB_THREAD ; i++)
            {
                unsigned long idxStart = (i*(this->numel()))/NB_THREAD;
                unsigned long idxEnd = (((i+1)*(this->numel()))/NB_THREAD);//-1;
                THREAD_EQ[i] = new C_thread_equal<pixelType>(this, &c, idxStart, idxEnd);
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
            return *this;
        }
        else
        {
            //do nothing
            return *this;
        }
    }
    else
    {
        //do nothing
        return *this;
    }
}


template<class pixelType> rawData<pixelType>& rawData<pixelType>::operator+ (const rawData& c)
{
    rawData* result;
    if(this->numDim==c.numDim)
    {
        if(this->numDim==0)
        {
            result = new rawData(this->DICOM_TYPE, 0);
            return *result;
        }
        else if(this->numDim>4)
        {
            result = new rawData(this->DICOM_TYPE, 0);
            return *result;
        }
        else if(this->DimX==c.DimX && this->DimY==c.DimY && this->DimZ==c.DimZ && this->DimT==c.DimT)
        {
            result = new rawData(this->DICOM_TYPE, this->numDim, this->DimX, this->DimY, this->DimZ, this->DimT);
            unsigned long NB_THREAD = 8;
            C_thread_add<pixelType>** THREAD_EQ = new C_thread_add<pixelType>*[NB_THREAD];
            for(unsigned long i=0 ; i<NB_THREAD ; i++)
            {
                unsigned long idxStart = (i*(this->numel()))/NB_THREAD;
                unsigned long idxEnd = (((i+1)*(this->numel()))/NB_THREAD);//-1;
                THREAD_EQ[i] = new C_thread_add<pixelType>(this, &c, result, idxStart, idxEnd);
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
            return *result;
        }
        else
        {
            //return empty
            result = new rawData(this->DICOM_TYPE, 0);
            return *result;
        }
    }
    else
    {
        //return empty data
        result = new rawData(this->DICOM_TYPE, 0);
        return *result;
    }
}


template<class pixelType> rawData<pixelType>& rawData<pixelType>::operator+ (const pixelType& c)
{
    rawData* result = new rawData(this->DICOM_TYPE, this->numDim, this->DimX, this->DimY, this->DimZ, this->DimT);
    unsigned long NB_THREAD = 8;
    C_thread_add<pixelType>** THREAD_EQ = new C_thread_add<pixelType>*[NB_THREAD];
    for(unsigned long i=0 ; i<NB_THREAD ; i++)
    {
        unsigned long idxStart = (i*(this->numel()))/NB_THREAD;
        unsigned long idxEnd = (((i+1)*(this->numel()))/NB_THREAD);//-1;
        THREAD_EQ[i] = new C_thread_add<pixelType>(this, c, result, idxStart, idxEnd);
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
    return *result;
}

template<class pixelType> rawData<pixelType>& rawData<pixelType>::operator- (const rawData& c)
{
    rawData* result;
    if(this->numDim==c.numDim)
    {
        if(this->numDim==0)
        {
            result = new rawData(this->DICOM_TYPE, 0);
            return *result;
        }
        else if(this->numDim>4)
        {
            result = new rawData(this->DICOM_TYPE, 0);
            return *result;
        }
        else if(this->DimX==c.DimX && this->DimY==c.DimY && this->DimZ==c.DimZ && this->DimT==c.DimT)
        {
            result = new rawData(this->DICOM_TYPE, this->numDim, this->DimX, this->DimY, this->DimZ, this->DimT);
            unsigned long NB_THREAD = 8;
            C_thread_sub<pixelType>** THREAD_EQ = new C_thread_sub<pixelType>*[NB_THREAD];
            for(unsigned long i=0 ; i<NB_THREAD ; i++)
            {
                unsigned long idxStart = (i*(this->numel()))/NB_THREAD;
                unsigned long idxEnd = (((i+1)*(this->numel()))/NB_THREAD);//-1;
                THREAD_EQ[i] = new C_thread_sub<pixelType>(this, &c, result, idxStart, idxEnd);
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
            return *result;
        }
        else
        {
            //return empty
            result = new rawData(this->DICOM_TYPE, 0);
            return *result;
        }
    }
    else
    {
        //return empty data
        result = new rawData(this->DICOM_TYPE, 0);
        return *result;
    }
}


template<class pixelType> rawData<pixelType>& rawData<pixelType>::operator- (const pixelType& c)
{
    rawData* result = new rawData(this->DICOM_TYPE, this->numDim, this->DimX, this->DimY, this->DimZ, this->DimT);
    unsigned long NB_THREAD = 8;
    C_thread_sub<pixelType>** THREAD_EQ = new C_thread_sub<pixelType>*[NB_THREAD];
    for(unsigned long i=0 ; i<NB_THREAD ; i++)
    {
        unsigned long idxStart = (i*(this->numel()))/NB_THREAD;
        unsigned long idxEnd = (((i+1)*(this->numel()))/NB_THREAD);//-1;
        THREAD_EQ[i] = new C_thread_sub<pixelType>(this, c, result, idxStart, idxEnd);
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
    return *result;
}





template<class pixelType> rawData<pixelType>& rawData<pixelType>::operator+= (const rawData& c)
{
    *this = *this + c;
    return *this;
}
template<class pixelType> rawData<pixelType>& rawData<pixelType>::operator+= (const pixelType& c)
{
    *this = *this - c;
    return *this;
}
template<class pixelType> rawData<pixelType>& rawData<pixelType>::operator-= (const rawData& c)
{
    *this = *this - c;
    return *this;
}
template<class pixelType> rawData<pixelType>& rawData<pixelType>::operator-= (const pixelType& c)
{
    *this = *this - c;
    return *this;
}
template<class pixelType> rawData<pixelType>& rawData<pixelType>::operator*= (const rawData& c)
{
    *this = *this * c;
    return *this;
}
template<class pixelType> rawData<pixelType>& rawData<pixelType>::operator*= (const pixelType& c)
{
    *this = *this * c;
    return *this;
}
template<class pixelType> rawData<pixelType>& rawData<pixelType>::operator/= (const rawData& c)
{
    *this = *this / c;
    return *this;
}
template<class pixelType> rawData<pixelType>& rawData<pixelType>::operator/= (const pixelType& c)
{
    *this = *this / c;
    return *this;
}

template<class pixelType> rawData<pixelType>& rawData<pixelType>::operator* (const rawData& c)
{
    if(this->numDim==c.numDim)
    {
        if(this->numDim==0)
        {
            rawData* result = new rawData(this->DICOM_TYPE, 0);
            return *result;
        }
        else if(this->numDim>4)
        {
            rawData* result = new rawData(this->DICOM_TYPE, 0);
            return *result;
        }
        else if(this->DimX==c.DimX && this->DimY==c.DimY && this->DimZ==c.DimZ && this->DimT==c.DimT)
        {
            rawData* result = new rawData(this->DICOM_TYPE, this->numDim, this->DimX, this->DimY, this->DimZ, this->DimT);
            unsigned long NB_THREAD = 8;
            C_thread_mult<pixelType>** THREAD_EQ = new C_thread_mult<pixelType>*[NB_THREAD];
            for(unsigned long i=0 ; i<NB_THREAD ; i++)
            {
                unsigned long idxStart = (i*(this->numel()))/NB_THREAD;
                unsigned long idxEnd = (((i+1)*(this->numel()))/NB_THREAD);//-1;
                THREAD_EQ[i] = new C_thread_mult<pixelType>(this, &c, result, idxStart, idxEnd);
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
            return *result;
        }
        else
        {
            //return empty
            rawData* result = new rawData(this->DICOM_TYPE, 0);
            return *result;
        }
    }
    else
    {
        //return empty data
        rawData* result = new rawData(this->DICOM_TYPE, 0);
        return *result;
    }
}

template<class pixelType> rawData<pixelType>& rawData<pixelType>::operator* (const pixelType& x)
{
    rawData* result = new rawData(this->DICOM_TYPE, this->numDim, this->DimX, this->DimY, this->DimZ, this->DimT);
    unsigned long NB_THREAD = 8;
    C_thread_mult<pixelType>** THREAD_EQ = new C_thread_mult<pixelType>*[NB_THREAD];
    for(unsigned long i=0 ; i<NB_THREAD ; i++)
    {
        unsigned long idxStart = (i*(this->numel()))/NB_THREAD;
        unsigned long idxEnd = (((i+1)*(this->numel()))/NB_THREAD);//-1;
        THREAD_EQ[i] = new C_thread_mult<pixelType>(this, x, result, idxStart, idxEnd);
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
    return *result;

}


template<class pixelType> rawData<pixelType>& rawData<pixelType>::operator/ (const rawData& c)
{
    if(this->numDim==c.numDim)
    {
        if(this->numDim==0)
        {
            rawData* result = new rawData(this->DICOM_TYPE, 0);
            return *result;
        }
        else if(this->numDim>4)
        {
            rawData* result = new rawData(this->DICOM_TYPE, 0);
            return *result;
        }
        else if(this->DimX==c.DimX && this->DimY==c.DimY && this->DimZ==c.DimZ && this->DimT==c.DimT)
        {
            rawData* result = new rawData(this->DICOM_TYPE, this->numDim, this->DimX, this->DimY, this->DimZ, this->DimT);
            unsigned long NB_THREAD = 8;
            C_thread_div<pixelType>** THREAD_EQ = new C_thread_div<pixelType>*[NB_THREAD];
            for(unsigned long i=0 ; i<NB_THREAD ; i++)
            {
                unsigned long idxStart = (i*(this->numel()))/NB_THREAD;
                unsigned long idxEnd = (((i+1)*(this->numel()))/NB_THREAD);//-1;
                THREAD_EQ[i] = new C_thread_div<pixelType>(this, &c, result, idxStart, idxEnd);
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
            return *result;
        }
        else
        {
            //return empty
            rawData* result = new rawData(this->DICOM_TYPE, 0);
            return *result;
        }
    }
    else
    {
        //return empty data
        rawData* result = new rawData(this->DICOM_TYPE, 0);
        return *result;
    }
};

template<class pixelType> rawData<pixelType>& rawData<pixelType>::operator/ (const pixelType& x)
{
    rawData* result = new rawData(this->DICOM_TYPE, this->numDim, this->DimX, this->DimY, this->DimZ, this->DimT);
    unsigned long NB_THREAD = 8;
    C_thread_div<pixelType>** THREAD_EQ = new C_thread_div<pixelType>*[NB_THREAD];
    for(unsigned long i=0 ; i<NB_THREAD ; i++)
    {
        unsigned long idxStart = (i*(this->numel()))/NB_THREAD;
        unsigned long idxEnd = (((i+1)*(this->numel()))/NB_THREAD);//-1;
        THREAD_EQ[i] = new C_thread_div<pixelType>(this, x, result, idxStart, idxEnd);
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
    return *result;

};


template<class pixelType> rawData<pixelType>& rawData<pixelType>::operator> (const rawData& c)
{
    if(this->numDim==c.numDim)
    {
        if(this->numDim==0)
        {
            rawData* result = new rawData(this->DICOM_TYPE, 0);
            return *result;
        }
        else if(this->numDim>4)
        {
            rawData* result = new rawData(this->DICOM_TYPE, 0);
            return *result;
        }
        else if(this->DimX==c.DimX && this->DimY==c.DimY && this->DimZ==c.DimZ && this->DimT==c.DimT)
        {
            rawData* result = new rawData(this->DICOM_TYPE, this->numDim, this->DimX, this->DimY, this->DimZ, this->DimT);
            unsigned long NB_THREAD = 8;
            C_thread_larger<pixelType>** THREAD_EQ = new C_thread_larger<pixelType>*[NB_THREAD];
            for(unsigned long i=0 ; i<NB_THREAD ; i++)
            {
                unsigned long idxStart = (i*(this->numel()))/NB_THREAD;
                unsigned long idxEnd = (((i+1)*(this->numel()))/NB_THREAD);//-1;
                THREAD_EQ[i] = new C_thread_larger<pixelType>(this, &c, result, idxStart, idxEnd);
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
            return *result;
        }
        else
        {
            //return empty
            rawData* result = new rawData(this->DICOM_TYPE, 0);
            return *result;
        }
    }
    else
    {
        //return empty data
        rawData* result = new rawData(this->DICOM_TYPE, 0);
        return *result;
    }
};
template<class pixelType> rawData<pixelType>& rawData<pixelType>::operator> (const pixelType& c)
{
    if(this->numDim==0)
    {
        rawData* result = new rawData(this->DICOM_TYPE, 0);
        return *result;
    }
    if(this->numDim>4)
    {
        rawData* result = new rawData(this->DICOM_TYPE, 0);
        return *result;
    }

    rawData* result = new rawData(this->DICOM_TYPE, this->numDim, this->DimX, this->DimY, this->DimZ, this->DimT);
    unsigned long NB_THREAD = 8;
    C_thread_larger<pixelType>** THREAD_EQ = new C_thread_larger<pixelType>*[NB_THREAD];
    for(unsigned long i=0 ; i<NB_THREAD ; i++)
    {
        unsigned long idxStart = (i*(this->numel()))/NB_THREAD;
        unsigned long idxEnd = (((i+1)*(this->numel()))/NB_THREAD);//-1;
        THREAD_EQ[i] = new C_thread_larger<pixelType>(this, c, result, idxStart, idxEnd);
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
    return *result;
};
template<class pixelType> rawData<pixelType>& rawData<pixelType>::operator>= (const rawData& c)
{
    if(this->numDim==c.numDim)
    {
        if(this->numDim==0)
        {
            rawData* result = new rawData(this->DICOM_TYPE, 0);
            return *result;
        }
        else if(this->numDim>4)
        {
            rawData* result = new rawData(this->DICOM_TYPE, 0);
            return *result;
        }
        else if(this->DimX==c.DimX && this->DimY==c.DimY && this->DimZ==c.DimZ && this->DimT==c.DimT)
        {
            rawData* result = new rawData(this->DICOM_TYPE, this->numDim, this->DimX, this->DimY, this->DimZ, this->DimT);
            unsigned long NB_THREAD = 8;
            C_thread_larger_or_equal<pixelType>** THREAD_EQ = new C_thread_larger_or_equal<pixelType>*[NB_THREAD];
            for(unsigned long i=0 ; i<NB_THREAD ; i++)
            {
                unsigned long idxStart = (i*(this->numel()))/NB_THREAD;
                unsigned long idxEnd = (((i+1)*(this->numel()))/NB_THREAD);//-1;
                THREAD_EQ[i] = new C_thread_larger_or_equal<pixelType>(this, &c, result, idxStart, idxEnd);
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
            return *result;
        }
        else
        {
            //return empty
            rawData* result = new rawData(this->DICOM_TYPE, 0);
            return *result;
        }
    }
    else
    {
        //return empty data
        rawData* result = new rawData(this->DICOM_TYPE, 0);
        return *result;
    }
};
template<class pixelType> rawData<pixelType>& rawData<pixelType>::operator>= (const pixelType& c)
{
    if(this->numDim==0)
    {
        rawData* result = new rawData(this->DICOM_TYPE, 0);
        return *result;
    }
    if(this->numDim>4)
    {
        rawData* result = new rawData(this->DICOM_TYPE, 0);
        return *result;
    }

    rawData* result = new rawData(this->DICOM_TYPE, this->numDim, this->DimX, this->DimY, this->DimZ, this->DimT);
    unsigned long NB_THREAD = 8;
    C_thread_larger_or_equal<pixelType>** THREAD_EQ = new C_thread_larger_or_equal<pixelType>*[NB_THREAD];
    for(unsigned long i=0 ; i<NB_THREAD ; i++)
    {
        unsigned long idxStart = (i*(this->numel()))/NB_THREAD;
        unsigned long idxEnd = (((i+1)*(this->numel()))/NB_THREAD);//-1;
        THREAD_EQ[i] = new C_thread_larger_or_equal<pixelType>(this, c, result, idxStart, idxEnd);
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
    return *result;
};

template<class pixelType> rawData<pixelType>& rawData<pixelType>::operator< (const rawData& c)
{
    if(this->numDim==c.numDim)
    {
        if(this->numDim==0)
        {
            rawData* result = new rawData(this->DICOM_TYPE, 0);
            return *result;
        }
        else if(this->numDim>4)
        {
            rawData* result = new rawData(this->DICOM_TYPE, 0);
            return *result;
        }
        else if(this->DimX==c.DimX && this->DimY==c.DimY && this->DimZ==c.DimZ && this->DimT==c.DimT)
        {
            rawData* result = new rawData(this->DICOM_TYPE, this->numDim, this->DimX, this->DimY, this->DimZ, this->DimT);
            unsigned long NB_THREAD = 8;
            C_thread_smaller<pixelType>** THREAD_EQ = new C_thread_smaller<pixelType>*[NB_THREAD];
            for(unsigned long i=0 ; i<NB_THREAD ; i++)
            {
                unsigned long idxStart = (i*(this->numel()))/NB_THREAD;
                unsigned long idxEnd = (((i+1)*(this->numel()))/NB_THREAD);//-1;
                THREAD_EQ[i] = new C_thread_smaller<pixelType>(this, &c, result, idxStart, idxEnd);
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
            return *result;
        }
        else
        {
            //return empty
            rawData* result = new rawData(this->DICOM_TYPE, 0);
            return *result;
        }
    }
    else
    {
        //return empty data
        rawData* result = new rawData(this->DICOM_TYPE, 0);
        return *result;
    }
};

template<class pixelType> rawData<pixelType>& rawData<pixelType>::operator< (const pixelType& c)
{
    if(this->numDim==0)
    {
        rawData* result = new rawData(this->DICOM_TYPE, 0);
        return *result;
    }
    if(this->numDim>4)
    {
        rawData* result = new rawData(this->DICOM_TYPE, 0);
        return *result;
    }

    rawData* result = new rawData(this->DICOM_TYPE, this->numDim, this->DimX, this->DimY, this->DimZ, this->DimT);
    unsigned long NB_THREAD = 8;
    C_thread_smaller<pixelType>** THREAD_EQ = new C_thread_smaller<pixelType>*[NB_THREAD];
    for(unsigned long i=0 ; i<NB_THREAD ; i++)
    {
        unsigned long idxStart = (i*(this->numel()))/NB_THREAD;
        unsigned long idxEnd = (((i+1)*(this->numel()))/NB_THREAD);//-1;
        THREAD_EQ[i] = new C_thread_smaller<pixelType>(this, c, result, idxStart, idxEnd);
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
    return *result;
};

template<class pixelType> rawData<pixelType>& rawData<pixelType>::operator<= (const rawData& c)
{
    if(this->numDim==c.numDim)
    {
        if(this->numDim==0)
        {
            rawData* result = new rawData(this->DICOM_TYPE, 0);
            return *result;
        }
        else if(this->numDim>4)
        {
            rawData* result = new rawData(this->DICOM_TYPE, 0);
            return *result;
        }
        else if(this->DimX==c.DimX && this->DimY==c.DimY && this->DimZ==c.DimZ && this->DimT==c.DimT)
        {
            rawData* result = new rawData(this->DICOM_TYPE, this->numDim, this->DimX, this->DimY, this->DimZ, this->DimT);
            unsigned long NB_THREAD = 8;
            C_thread_smaller_or_equal<pixelType>** THREAD_EQ = new C_thread_smaller_or_equal<pixelType>*[NB_THREAD];
            for(unsigned long i=0 ; i<NB_THREAD ; i++)
            {
                unsigned long idxStart = (i*(this->numel()))/NB_THREAD;
                unsigned long idxEnd = (((i+1)*(this->numel()))/NB_THREAD);//-1;
                THREAD_EQ[i] = new C_thread_smaller_or_equal<pixelType>(this, &c, result, idxStart, idxEnd);
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
            return *result;
        }
        else
        {
            //return empty
            rawData* result = new rawData(this->DICOM_TYPE, 0);
            return *result;
        }
    }
    else
    {
        //return empty data
        rawData* result = new rawData(this->DICOM_TYPE, 0);
        return *result;
    }
};

template<class pixelType> rawData<pixelType>& rawData<pixelType>::operator<= (const pixelType& c)
{
    if(this->numDim==0)
    {
        rawData* result = new rawData(this->DICOM_TYPE, 0);
        return *result;
    }
    if(this->numDim>4)
    {
        rawData* result = new rawData(this->DICOM_TYPE, 0);
        return *result;
    }

    rawData* result = new rawData(this->DICOM_TYPE, this->numDim, this->DimX, this->DimY, this->DimZ, this->DimT);
    unsigned long NB_THREAD = 8;
    C_thread_smaller_or_equal<pixelType>** THREAD_EQ = new C_thread_smaller_or_equal<pixelType>*[NB_THREAD];
    for(unsigned long i=0 ; i<NB_THREAD ; i++)
    {
        unsigned long idxStart = (i*(this->numel()))/NB_THREAD;
        unsigned long idxEnd = (((i+1)*(this->numel()))/NB_THREAD);//-1;
        THREAD_EQ[i] = new C_thread_smaller_or_equal<pixelType>(this, c, result, idxStart, idxEnd);
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
    return *result;
};




template<class pixelType> pixelType rawData<pixelType>::getSum(void)
{
    unsigned long NB_THREAD = 8;
    C_thread_sum<pixelType>** THREAD_EQ = new C_thread_sum<pixelType>*[NB_THREAD];
    for(unsigned long i=0 ; i<NB_THREAD ; i++)
    {
        unsigned long idxStart = (i*(this->numel()))/NB_THREAD;
        unsigned long idxEnd = (((i+1)*(this->numel()))/NB_THREAD);//-1;
        THREAD_EQ[i] = new C_thread_sum<pixelType>(this, idxStart, idxEnd);
    }
    for(unsigned long i=0 ; i<NB_THREAD ; i++)
    {
        THREAD_EQ[i]->start();
    }
    for(unsigned long i=0 ; i<NB_THREAD ; i++)
    {
        THREAD_EQ[i]->wait_for_exit();
    }
    pixelType sum = 0.0;
    for(unsigned long i=0 ; i<NB_THREAD ; i++)
    {
        sum += THREAD_EQ[i]->partialSum;
        delete THREAD_EQ[i];
    }
    delete THREAD_EQ;
    return sum;
};

template<class pixelType> pixelType rawData<pixelType>::getMean(void)
{
    pixelType sum = getSum();
    return sum/((pixelType) this->numel());
}

template<class pixelType> pixelType rawData<pixelType>::getVar(void)
{
    pixelType mean_ = getMean();
    rawData<pixelType>* temp = new rawData<pixelType>(DICOM_TYPE, numDim, DimX, DimY, DimZ, DimT);
    *temp = *this;
    *temp -= mean_;
    *temp *= *temp;
    mean_ = temp->getSum();
    delete temp;
    return mean_/((pixelType) (this->numel()-1));
}


template<class pixelType> pixelType rawData<pixelType>::getValue(unsigned long idxPix)
{
    if(numDim==1)
    {
        return raw1D[idxPix];
    }
    else if(numDim==2)
    {
        return raw2D[getIndexI2D(idxPix)][getIndexJ2D(idxPix)];
    }
    else if(numDim==3)
    {
        unsigned long k = getIndexK3D(idxPix);
        unsigned long i = getIndexI3D(idxPix, k);
        unsigned long j = getIndexJ3D(idxPix, i, k);
        return raw3D[i][j][k];
    }
    else if(numDim==4)
    {
        unsigned long l = getIndexL4D(idxPix);
        unsigned long k = getIndexK4D(idxPix, l);
        unsigned long i = getIndexI4D(idxPix, k, l);
        unsigned long j = getIndexJ4D(idxPix, i, k, l);
        return raw4D[i][j][k][l];
    }
    else
    {
        return 0.0;
    }
}

template<class pixelType> bool rawData<pixelType>::setValue(unsigned long idxPix, pixelType newValue)
{
    if(numDim==1)
    {
        raw1D[idxPix] = newValue;
    }
    else if(numDim==2)
    {
        raw2D[getIndexI2D(idxPix)][getIndexJ2D(idxPix)] = newValue;
    }
    else if(numDim==3)
    {
        unsigned long k = getIndexK3D(idxPix);
        unsigned long i = getIndexI3D(idxPix, k);
        unsigned long j = getIndexJ3D(idxPix, i, k);
        raw3D[i][j][k] = newValue;
    }
    else if(numDim==4)
    {
        unsigned long l = getIndexL4D(idxPix);
        unsigned long k = getIndexK4D(idxPix, l);
        unsigned long i = getIndexI4D(idxPix, k, l);
        unsigned long j = getIndexJ4D(idxPix, i, k, l);
        raw4D[i][j][k][l] = newValue;
    }
    else
    {
        return false;
    }
    return true;
}


template<class pixelType> void rawData<pixelType>::removeSingle()
{
    unsigned long NB_THREAD = 8;
    C_thread_remove_isolated<pixelType>** THREAD_EQ = new C_thread_remove_isolated<pixelType>*[NB_THREAD];
    for(unsigned long i=0 ; i<NB_THREAD ; i++)
    {
        unsigned long idxStart = (i*(this->numel()))/NB_THREAD;
        unsigned long idxEnd = (((i+1)*(this->numel()))/NB_THREAD);//-1;
        THREAD_EQ[i] = new C_thread_remove_isolated<pixelType>(this, idxStart, idxEnd);
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
    return;
}




template<class pixelType>unsigned long rawData<pixelType>::getIndex1D(unsigned long i)
{
    return i;
}//OK
template<class pixelType>unsigned long rawData<pixelType>::getIndex2D(unsigned long i, unsigned long j)
{
    return i*DimY + j;
}//OK
template<class pixelType>unsigned long rawData<pixelType>::getIndex2DUnmosaic(unsigned long i, unsigned long j, unsigned long nb_slice)
{
    //get dimension of embeding square
    unsigned long subDimX, subDimY, squareDim;
    if(sqrt(nb_slice)==floor(sqrt(nb_slice)))
    {
        squareDim = (unsigned long) floor(sqrt(nb_slice));
    }
    else
    {
        squareDim = (unsigned long) floor(sqrt(nb_slice)) + 1;
    }
    subDimX = DimX/squareDim;
    subDimY = DimY/squareDim;

    //determine slice
    unsigned long Z, I, J;
    I = i/subDimX;
    J = j/subDimY;
    Z = I + squareDim*J;
    if(Z<nb_slice)
    {
        unsigned long i2 = i - I*subDimX;
        unsigned long j2 = j - J*subDimY;
        return Z*subDimX*subDimY + j2*subDimX + i2;
    }
    else
    {
        //return a number out of range
        return DimX*DimY;
    }
}
template<class pixelType>unsigned long rawData<pixelType>::getIndex3DUnmosaic(unsigned long i, unsigned long j, unsigned long k) //OK
{
    //get dimension of embeding square
    unsigned long /*supDimX,*/ supDimY, squareDim;
    if(sqrt(DimZ)==floor(sqrt(DimZ)))
    {
        squareDim = (unsigned long) floor(sqrt(DimZ));
    }
    else
    {
        squareDim = (unsigned long) floor(sqrt(DimZ)) + 1;
    }

    //dimension of the 2D image
    supDimY = DimY*squareDim;

    //convert i to I in 2D image
    unsigned long I = i + (floor(k/squareDim))*DimX;

    //convert j to J in 2D image
    unsigned long J = j + (k%squareDim)*DimY;

    return supDimY*I + J;
}
template<class pixelType> unsigned long rawData<pixelType>::getIndex3D(unsigned long i, unsigned long j, unsigned long k)
{
    return DimY*i + j + DimX*DimY*k;
}
template<class pixelType> unsigned long rawData<pixelType>::getIndex4D(unsigned long i, unsigned long j, unsigned long k, unsigned long l)
{
    return DimY*i + j + DimX*DimY*k + DimX*DimY*DimZ*l;
}




template<class pixelType> unsigned long rawData<pixelType>::getIndexI1D(unsigned long idx)
{
    return idx;
}
template<class pixelType> unsigned long rawData<pixelType>::getIndexI2D(unsigned long idx)
{
    return floor(idx/DimY);
}
template<class pixelType> unsigned long rawData<pixelType>::getIndexI3D(unsigned long idx, unsigned long k)
{
    return floor((idx - (k*(DimX*DimY)))/DimY);
}
template<class pixelType> unsigned long rawData<pixelType>::getIndexI4D(unsigned long idx, unsigned long k, unsigned long l)
{
    return floor((idx - l*(DimX*DimY*DimZ) - k*(DimX*DimY))/DimY);
}
template<class pixelType> unsigned long rawData<pixelType>::getIndexJ2D(unsigned long idx)
{
    return idx - (floor(idx/DimY)*DimY);
}
template<class pixelType> unsigned long rawData<pixelType>::getIndexJ3D(unsigned long idx, unsigned long i, unsigned long k)
{
    return idx - (k*(DimX*DimY)) - i*DimY;
}
template<class pixelType> unsigned long rawData<pixelType>::getIndexJ4D(unsigned long idx, unsigned long i, unsigned long k, unsigned long l)
{
    return idx - l*(DimX*DimY*DimZ) - k*(DimX*DimY) - i*DimY;
}
template<class pixelType> unsigned long rawData<pixelType>::getIndexK3D(unsigned long idx)
{
    return floor(idx/(DimX*DimY));
}
template<class pixelType> unsigned long rawData<pixelType>::getIndexK4D(unsigned long idx, unsigned long l)
{
    return floor((idx - l*(DimX*DimY*DimZ))/(DimX*DimY));
}
template<class pixelType> unsigned long rawData<pixelType>::getIndexL4D(unsigned long idx)
{
    return floor(idx/(DimX*DimY*DimZ));
}








template<class pixelType> pixelType rawData<pixelType>::getMax()
{
    pixelType m = (pixelType) 0;
    if(numDim==1 && raw1D!=NULL)
    {
        m = raw1D[0];
        for(unsigned long i=1 ; i<DimX ; i++)
        {
            if(m<raw1D[i]) m = raw1D[i];
        }
    }
    else if(numDim==2 && raw2D!=NULL)
    {
        m = raw2D[0][0];
        for(unsigned long i=0 ; i<DimX ; i++)
        {
            for(unsigned long j=0 ; j<DimY ; j++)
            {
                if(m<raw2D[i][j]) m = raw2D[i][j];
            }
        }
    }
    else if(numDim==3 && raw3D!=NULL)
    {
        m = raw3D[0][0][0];
        for(unsigned long i=0 ; i<DimX ; i++)
        {
            for(unsigned long j=0 ; j<DimY ; j++)
            {
                for(unsigned long k=0 ; k<DimZ ; k++)
                {
                    if(m<raw3D[i][j][k]) m = raw3D[i][j][k];
                }
            }
        }
    }
    else if(numDim==4 && raw4D!=NULL)
    {
        m = raw4D[0][0][0][0];
        for(unsigned long i=0 ; i<DimX ; i++)
        {
            for(unsigned long j=0 ; j<DimY ; j++)
            {
                for(unsigned long k=0 ; k<DimZ ; k++)
                {
                    for(unsigned long l=0 ; l<DimT ; l++)
                    {
                        if(m<raw4D[i][j][k][l]) m = raw4D[i][j][k][l];
                    }
                }
            }
        }
    }
    return m;
}
template<class pixelType> pixelType rawData<pixelType>::getMin()
{
    pixelType m = (pixelType) 0;
    if(numDim==1 && raw1D!=NULL)
    {
        m = raw1D[0];
        for(unsigned long i=1 ; i<DimX ; i++)
        {
            if(m>raw1D[i]) m = raw1D[i];
        }
    }
    else if(numDim==2 && raw2D!=NULL)
    {
        m = raw2D[0][0];
        for(unsigned long i=0 ; i<DimX ; i++)
        {
            for(unsigned long j=0 ; j<DimY ; j++)
            {
                if(m>raw2D[i][j]) m = raw2D[i][j];
            }
        }
    }
    else if(numDim==3 && raw3D!=NULL)
    {
        m = raw3D[0][0][0];
        for(unsigned long i=0 ; i<DimX ; i++)
        {
            for(unsigned long j=0 ; j<DimY ; j++)
            {
                for(unsigned long k=0 ; k<DimZ ; k++)
                {
                    if(m>raw3D[i][j][k]) m = raw3D[i][j][k];
                }
            }
        }
    }
    else if(numDim==4 && raw4D!=NULL)
    {
        m = raw4D[0][0][0][0];
        for(unsigned long i=0 ; i<DimX ; i++)
        {
            for(unsigned long j=0 ; j<DimY ; j++)
            {
                for(unsigned long k=0 ; k<DimZ ; k++)
                {
                    for(unsigned long l=0 ; l<DimT ; l++)
                    {
                        if(m>raw4D[i][j][k][l]) m = raw4D[i][j][k][l];
                    }
                }
            }
        }
    }
    return m;
}

template<class pixelType> void rawData<pixelType>::crop1D(unsigned long xmin, unsigned long xmax)
{
    unsigned long newDimX = xmax-xmin+1;
    pixelType* tempData = new pixelType[newDimX];
    if(tempData==NULL) return;
    DimX = newDimX;
    for(unsigned long i=0 ; i<newDimX ; i++)
    {
        tempData[i] = raw1D[xmin+i];
    }
    delete raw1D;
    raw1D = tempData;
    return;
}
template<class pixelType> void rawData<pixelType>::crop2D(unsigned long xmin, unsigned long xmax, unsigned long ymin, unsigned long ymax)
{
    unsigned long newDimX = xmax-xmin+1;
    unsigned long newDimY = ymax-ymin+1;

    //alloc new data array
    pixelType** tempData = new pixelType*[newDimX];


    for(unsigned long i=0 ; i<newDimX ; i++)
    {
        tempData[i] = new pixelType[newDimY];
        for(unsigned long j=0 ; j<newDimY ; j++)
        {
            tempData[i][j] = raw2D[xmin+i][ymin+j];
        }
    }
    for(unsigned long i=0 ; i<DimX ; i++)
    {
        delete raw2D[i];
    }
    delete raw2D;
    DimX = newDimX;
    DimY = newDimY;
    raw2D = tempData;
    return;
//    unsigned long newDimX = xmax-xmin+1;
//    unsigned long newDimY = ymax-ymin+1;
//
//    //alloc new data array
//    pixelType** tempData = new pixelType*[newDimX];
//    for(unsigned long i=0 ; i<xmin ; i++)
//    {
//        delete raw2D[i];
//    }
//    for(unsigned long i=xmax+1 ; i<DimY ; i++)
//    {
//        delete raw2D[i];
//    }
//
//
//    for(unsigned long i=0 ; i<newDimX ; i++)
//    {
//        tempData[i] = new pixelType[newDimY];
//        for(unsigned long j=0 ; j<newDimY ; j++)
//        {
//            tempData[i][j] = raw2D[xmin+i][ymin+j];
//        }
//        delete raw2D[xmin+i];
//    }
//    delete raw2D;
//    DimX = newDimX;
//    DimY = newDimY;
//    raw2D = tempData;
//    return;
}
template<class pixelType> void rawData<pixelType>::crop3D(unsigned long xmin, unsigned long xmax, unsigned long ymin, unsigned long ymax, unsigned long zmin, unsigned long zmax)
{
    unsigned long newDimX = xmax-xmin+1;
    unsigned long newDimY = ymax-ymin+1;
    unsigned long newDimZ = zmax-zmin+1;

    //alloc new data array
    pixelType*** tempData = new pixelType**[newDimX];


    for(unsigned long i=0 ; i<newDimX ; i++)
    {
        tempData[i] = new pixelType*[newDimY];
        for(unsigned long j=0 ; j<newDimY ; j++)
        {
            tempData[i][j] = new pixelType[newDimZ];
            for(unsigned long k=0 ; k<newDimZ ; k++)
            {
                tempData[i][j][k] = raw3D[xmin+i][ymin+j][zmin+k];
            }
        }
    }
    for(unsigned long i=0 ; i<DimX ; i++)
    {
        for(unsigned long j=0 ; j<DimY ; j++)
        {
            delete raw3D[i][j];
        }
        delete raw3D[i];
    }
    delete raw3D;
    DimX = newDimX;
    DimY = newDimY;
    DimZ = newDimZ;
    raw3D = tempData;
    return;
}
template<class pixelType> void rawData<pixelType>::crop4D(unsigned long xmin, unsigned long xmax, unsigned long ymin, unsigned long ymax, unsigned long zmin, unsigned long zmax, unsigned long tmin, unsigned long tmax)
{
    unsigned long newDimX = xmax-xmin+1;
    unsigned long newDimY = ymax-ymin+1;
    unsigned long newDimZ = zmax-zmin+1;
    unsigned long newDimT = tmax-tmin+1;

    //alloc new data array
    pixelType**** tempData = new pixelType***[newDimX];
    for(unsigned long i=0 ; i<newDimX ; i++)
    {
        tempData[i] = new pixelType**[newDimY];
        for(unsigned long j=0 ; j<newDimY ; j++)
        {
            tempData[i][j] = new pixelType*[newDimZ];
            for(unsigned long k=0 ; k<newDimZ ; k++)
            {
                tempData[i][j][k] = new pixelType[newDimT];
                for(unsigned long l=0 ; l<newDimT ; l++)
                {
                    tempData[i][j][k][l] = raw4D[xmin+i][ymin+j][zmin+k][tmin+l];
                }
            }
        }
    }
    for(unsigned long i=0 ; i<DimX ; i++)
    {
        for(unsigned long j=0 ; j<DimY ; j++)
        {
            for(unsigned long k=0 ; k<DimZ ; k++)
            {
                delete raw4D[i][j][k];
            }
            delete raw4D[i][j];
        }
        delete raw4D[i];
    }
    delete raw4D;
    DimX = newDimX;
    DimY = newDimY;
    DimZ = newDimZ;
    DimT = newDimT;
    raw4D = tempData;
    return;
}


template<class pixelType> void rawData<pixelType>::swapXY2D(void)
{
    pixelType** tempData = new pixelType*[DimY];
    for(unsigned long i=0 ; i<DimY ; i++)
        tempData[i] = new pixelType[DimX];

    for(unsigned long i=0 ; i<DimX ; i++ )
    {
        for(unsigned long j=0 ; j<DimY ; j++ )
        {
            tempData[j][i] = raw2D[i][j];
        }
        delete raw2D[i];
    }
    delete raw2D;
    raw2D = tempData;
    pixelType tempD;
    tempD = pixDimX;
    pixDimX = pixDimY;
    pixDimX = tempD;
    unsigned long tempUL;
    tempUL = DimX;
    DimX = DimY;
    DimY = tempUL;
    return;
}
template<class pixelType> void rawData<pixelType>::swapXY3D(void)
{
    pixelType*** tempData = new pixelType**[DimY];
    for(unsigned long i=0 ; i<DimY ; i++)
    {
        tempData[i] = new pixelType*[DimX];
        for(unsigned long j=0 ; j<DimX ; j++)
        {
            tempData[i][j] = new pixelType[DimZ];
        }
    }


    for(unsigned long i=0 ; i<DimX ; i++ )
    {
        for(unsigned long j=0 ; j<DimY ; j++ )
        {
            for(unsigned long k=0 ; k<DimZ ; k++ )
            {
                tempData[j][i][k] = raw3D[i][j][k];
            }
            delete raw3D[i][j];
        }
        delete raw3D[i];
    }
    delete raw3D;
    raw3D = tempData;
    pixelType tempD;
    tempD = pixDimX;
    pixDimX = pixDimY;
    pixDimX = tempD;
    unsigned long tempUL;
    tempUL = DimX;
    DimX = DimY;
    DimY = tempUL;
    return;
}
template<class pixelType> void rawData<pixelType>::swapXY4D(void)
{
    pixelType**** tempData = new pixelType***[DimY];
    for(unsigned long i=0 ; i<DimY ; i++)
    {
        tempData[i] = new pixelType**[DimX];
        for(unsigned long j=0 ; j<DimX ; j++)
        {
            tempData[i][j] = new pixelType*[DimZ];
            for(unsigned long k=0 ; k<DimZ ; k++)
            {
                tempData[i][j][k] = new pixelType[DimT];
            }
        }
    }


    for(unsigned long i=0 ; i<DimX ; i++ )
    {
        for(unsigned long j=0 ; j<DimY ; j++ )
        {
            for(unsigned long k=0 ; k<DimZ ; k++ )
            {
                for(unsigned long l=0 ; l<DimT ; l++ )
                {
                    tempData[j][i][k][l] = raw4D[i][j][k][l];
                }
                delete raw4D[i][j][k];
            }
            delete raw4D[i][j];
        }
        delete raw4D[i];
    }
    delete raw4D;
    raw4D = tempData;
    pixelType tempD;
    tempD = pixDimX;
    pixDimX = pixDimY;
    pixDimX = tempD;
    unsigned long tempUL;
    tempUL = DimX;
    DimX = DimY;
    DimY = tempUL;
    return;
}
#endif // RAWDATA_H
