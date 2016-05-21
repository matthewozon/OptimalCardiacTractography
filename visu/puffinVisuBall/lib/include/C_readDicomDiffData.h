#ifndef C_READDICOMDIFFDATA_H
#define C_READDICOMDIFFDATA_H

#include <C_ReadDicomDiffInfo.h>
#include <rawData.h>
#include <C_thread.h>


class C_readDicomDiffData;


template<class pixelType> class C_thread_load : public C_thread
{
    friend class rawData<pixelType>;
    public:
        /** Default constructor */
        C_thread_load(pixelType* buff, rawData<pixelType>* A, long idxStart, long idxEnd)
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
            m_buff = buff;
        };

        /** Default destructor */
        virtual ~C_thread_load()
        {
            //
        };
        virtual void execute()
        {
            if(m_A->numDim==1)
            {
                for(long idx=m_idxStart ; idx<m_idxEnd ; idx++)
                {
                    ///store data
                    m_A->raw1D[idx] = m_buff[idx];
                }
            }
            else if(m_A->numDim==2)
            {
                long i, j;
                for(long idx=m_idxStart ; idx<m_idxEnd ; idx++)
                {
                    i = floor(idx/m_A->DimY);
                    j = idx - (floor(idx/m_A->DimY)*m_A->DimY);
                    m_A->raw2D[i][j] = m_buff[m_A->getIndex2D(i,j)];
                }
            }
            else if(m_A->numDim==3)
            {
                long i, j, k;
                for(long idx=m_idxStart ; idx<m_idxEnd ; idx++)
                {
                    k = floor(idx/(m_A->DimX*m_A->DimY));
                    i = floor((idx - (k*(m_A->DimX*m_A->DimY)))/m_A->DimY);
                    j = idx - (k*(m_A->DimX*m_A->DimY)) - i*m_A->DimY;
                    m_A->raw3D[i][j][k] = m_buff[m_A->getIndex3D(i,j,k)];
                }
            }
            else if(m_A->numDim==4)
            {
                long i, j, k, l;
                for(long idx=m_idxStart ; idx<m_idxEnd ; idx++)
                {
                    l = floor(idx/(m_A->DimX*m_A->DimY*m_A->DimZ));
                    k = floor((idx - l*(m_A->DimX*m_A->DimY*m_A->DimZ))/(m_A->DimX*m_A->DimY));
                    i = floor((idx - l*(m_A->DimX*m_A->DimY*m_A->DimZ) - k*(m_A->DimX*m_A->DimY))/m_A->DimY);
                    j = idx - l*(m_A->DimX*m_A->DimY*m_A->DimZ) - k*(m_A->DimX*m_A->DimY) - i*m_A->DimY;
                    m_A->raw4D[i][j][k][l] = m_buff[m_A->getIndex4D(i,j,k,l)];
                }
            }
        };
    protected:
    private:
        long m_idxStart;
        long m_idxEnd;
        rawData<pixelType>* m_A;
        pixelType* m_buff;
};


///should create a template for this class?
class C_readDicomDiffData : public C_ReadDicomDiffInfo
{
    public:
        /** Default constructor */
        C_readDicomDiffData(string fileName);
        /** Default destructor */
        virtual ~C_readDicomDiffData();

        rawData<double>* getRawData();
        rawData<double>* getUnmosaicData(unsigned short nb_slice_per_mosa); //require 2D data and return 3D rawData

    protected:
        template<typename X> rawData<X>* getRawData_();
        template<typename X> rawData<X>* getUnmosaicData_(unsigned short nb_slice_per_mosa); //require 2D data and return 3D rawData
        template<typename X> rawData<double>* conver2double(rawData<X>* tempRaw);
};

#endif // C_READDICOMDIFFDATA_H
