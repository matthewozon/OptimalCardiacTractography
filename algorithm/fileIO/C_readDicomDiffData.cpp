#include <C_readDicomDiffData.h>

C_readDicomDiffData::C_readDicomDiffData(string fileName):C_ReadDicomDiffInfo(fileName)
{
    //ctor
}

C_readDicomDiffData::~C_readDicomDiffData()
{
    //dtor
}


rawData<double>* C_readDicomDiffData::getRawData()
{
    //get type of data
    rawData<double>* mRawData=NULL;

    switch(getPixelRepresentation())
    {
        case DICOM_BOOL:
        {
            rawData<bool>* tempRaw = getRawData_<bool>();
            if(tempRaw!=NULL)
            {
                mRawData = conver2double(tempRaw);
                delete tempRaw;
            }
            break;
        }
        case DICOM_CHAR:
        {
            rawData<char>* tempRaw = getRawData_<char>();
            if(tempRaw!=NULL)
            {
                mRawData = conver2double(tempRaw);
                delete tempRaw;
            }
            break;
        }
        case DICOM_SHORT:
        {
            rawData<short>* tempRaw = getRawData_<short>();
            if(tempRaw!=NULL)
            {
                mRawData = conver2double(tempRaw);
                delete tempRaw;
            }
            break;
        }
        case DICOM_LONG:
        {
            rawData<long>* tempRaw = getRawData_<long>();
            if(tempRaw!=NULL)
            {
                mRawData = conver2double(tempRaw);
                delete tempRaw;
            }
            break;
        }
        case DICOM_UCHAR:
        {
            rawData<unsigned char>* tempRaw = getRawData_<unsigned char>();
            if(tempRaw!=NULL)
            {
                mRawData = conver2double(tempRaw);
                delete tempRaw;
            }
            break;
        }
        case DICOM_USHORT:
        {
            rawData<unsigned short>* tempRaw = getRawData_<unsigned short>();
            if(tempRaw!=NULL)
            {
                mRawData = conver2double(tempRaw);
                delete tempRaw;
            }
            break;
        }
        case DICOM_ULONG:
        {
            rawData<unsigned long>* tempRaw = getRawData_<unsigned long>();
            if(tempRaw!=NULL)
            {
                mRawData = conver2double(tempRaw);
                delete tempRaw;
            }
            break;
        }
        case DICOM_FLOAT:
        {
            rawData<float>* tempRaw = getRawData_<float>();
            if(tempRaw!=NULL)
            {
                mRawData = conver2double(tempRaw);
                delete tempRaw;
            }
            break;
        }
        case DICOM_DOUBLE:
        {
            rawData<double>* tempRaw = getRawData_<double>();
            if(tempRaw!=NULL)
            {
                mRawData = conver2double(tempRaw);
                delete tempRaw;
            }
            break;
        }
        default:
        {
            break;
        }
    }

    return mRawData;
}
rawData<double>* C_readDicomDiffData::getUnmosaicData(unsigned short nb_slice_per_mosa)
{
    //get type of data
    rawData<double>* mRawData=NULL;

    switch(getPixelRepresentation())
    {
        case DICOM_BOOL:
        {
            rawData<bool>* tempRaw = getUnmosaicData_<bool>(nb_slice_per_mosa);
            if(tempRaw!=NULL)
            {
                mRawData = conver2double(tempRaw);
                delete tempRaw;
            }
            break;
        }
        case DICOM_CHAR:
        {
            rawData<char>* tempRaw = getUnmosaicData_<char>(nb_slice_per_mosa);
            if(tempRaw!=NULL)
            {
                mRawData = conver2double(tempRaw);
                delete tempRaw;
            }
            break;
        }
        case DICOM_SHORT:
        {
            rawData<short>* tempRaw = getUnmosaicData_<short>(nb_slice_per_mosa);
            if(tempRaw!=NULL)
            {
                mRawData = conver2double(tempRaw);
                delete tempRaw;
            }
            break;
        }
        case DICOM_LONG:
        {
            rawData<long>* tempRaw = getUnmosaicData_<long>(nb_slice_per_mosa);
            if(tempRaw!=NULL)
            {
                mRawData = conver2double(tempRaw);
                delete tempRaw;
            }
            break;
        }
        case DICOM_UCHAR:
        {
            rawData<unsigned char>* tempRaw = getUnmosaicData_<unsigned char>(nb_slice_per_mosa);
            if(tempRaw!=NULL)
            {
                mRawData = conver2double(tempRaw);
                delete tempRaw;
            }
            break;
        }
        case DICOM_USHORT:
        {
            rawData<unsigned short>* tempRaw = getUnmosaicData_<unsigned short>(nb_slice_per_mosa);
            if(tempRaw!=NULL)
            {
                mRawData = conver2double(tempRaw);
                delete tempRaw;
            }
            break;
        }
        case DICOM_ULONG:
        {
            rawData<unsigned long>* tempRaw = getUnmosaicData_<unsigned long>(nb_slice_per_mosa);
            if(tempRaw!=NULL)
            {
                mRawData = conver2double(tempRaw);
                delete tempRaw;
            }
            break;
        }
        case DICOM_FLOAT:
        {
            rawData<float>* tempRaw = getUnmosaicData_<float>(nb_slice_per_mosa);
            if(tempRaw!=NULL)
            {
                mRawData = conver2double(tempRaw);
                delete tempRaw;
            }
            break;
        }
        case DICOM_DOUBLE:
        {
            rawData<double>* tempRaw = getUnmosaicData_<double>(nb_slice_per_mosa);
            if(tempRaw!=NULL)
            {
                mRawData = conver2double(tempRaw);
                delete tempRaw;
            }
            break;
        }
        default:
        {
            break;
        }
    }

    return mRawData;
}





template<typename X> rawData<X>* C_readDicomDiffData::getRawData_()
{
    //returned pointer
    rawData<X>* rawDataX = NULL;

    //get number of byte in image buffer
    unsigned long buffLen = mImage.GetBufferLength();

    //allocate a temporary buffer of buffLen bytes
    char* buff = new char[buffLen];
    if(buff==NULL) return NULL;

    //get image buffer into temporary buffer
    if(mImage.GetBuffer(buff))
    {
        ///allocate storage
        rawDataX = new rawData<X>(getPixelRepresentation(), mImage.GetNumberOfDimensions(), (unsigned int*) mImage.GetDimensions());
        if(rawDataX==NULL)
        {
            delete buff;
            return NULL;
        }
        if(rawDataX->raw1D==NULL && rawDataX->raw2D==NULL && rawDataX->raw3D==NULL && rawDataX->raw4D==NULL)
        {
            delete buff;
            delete rawDataX;
            return NULL;
        }

        rawDataX->pixDimX = getPixelDimX();
        rawDataX->pixDimY = getPixelDimY();
        rawDataX->pixDimZ = getPixelDimZ();
        rawDataX->pixDimT = 1.0;

        //reorganize buffer into rawData structure
        X* bufferX = (X*) buff;

        unsigned long NB_THREAD = 8;
        C_thread_load<X>** THREAD_EQ = new C_thread_load<X>*[NB_THREAD];
        for(unsigned long i=0 ; i<NB_THREAD ; i++)
        {
            unsigned long idxStart = (i*(rawDataX->numel()))/NB_THREAD;
            unsigned long idxEnd = (((i+1)*(rawDataX->numel()))/NB_THREAD);//-1;
            THREAD_EQ[i] = new C_thread_load<X>(bufferX, rawDataX, idxStart, idxEnd);
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

//        if(rawDataX->numDim==1)
//        {
//            for(unsigned long i=0 ; i<rawDataX->DimX ; i++)
//            {
//                ///store data
//                rawDataX->raw1D[i] = bufferX[rawDataX->getIndex1D(i)];
//            }
//        }
//        else if(rawDataX->numDim==2)
//        {
//            for(unsigned long i=0 ; i<rawDataX->DimX ; i++)
//            {
//                for(unsigned long j=0 ; j<rawDataX->DimY ; j++)
//                {
//                    rawDataX->raw2D[i][j] = bufferX[rawDataX->getIndex2D(i,j)];
//
//                }
//            }
//        }
//        else if(rawDataX->numDim==3)
//        {
//            for(unsigned long i=0 ; i<rawDataX->DimX ; i++)
//            {
//                for(unsigned long j=0 ; j<rawDataX->DimY ; j++)
//                {
//                    for(unsigned long k=0 ; k<rawDataX->DimZ ; k++)
//                    {
//                        ///store data
//                        rawDataX->raw3D[i][j][k] = bufferX[rawDataX->getIndex3D(i,j,k)];
//                    }
//                }
//            }
//        }
//        else if(rawDataX->numDim==4)
//        {
//            for(unsigned long i=0 ; i<rawDataX->DimX ; i++)
//            {
//                for(unsigned long j=0 ; j<rawDataX->DimY ; j++)
//                {
//                    for(unsigned long k=0 ; k<rawDataX->DimZ ; k++)
//                    {
//                        for(unsigned long l=0 ; l<rawDataX->DimT ; l++)
//                        {
//                            ///store data
//                            rawDataX->raw4D[i][j][k][l] = bufferX[rawDataX->getIndex4D(i,j,k,l)];
//                        }
//                    }
//                }
//            }
//        }


    }
    delete buff;
    return rawDataX;
}
template<typename X> rawData<X>* C_readDicomDiffData::getUnmosaicData_(unsigned short nb_slice_per_mosa)
{
    //returned pointer
    rawData<X>* rawDataX = NULL;

    if(mImage.GetNumberOfDimensions()!=2) return NULL;

    //get number of byte in image buffer
    unsigned long buffLen = mImage.GetBufferLength();

    //allocate a temporary buffer of buffLen bytes
    char* buff = new char[buffLen];
    if(buff==NULL) return NULL;

    //get image buffer into temporary buffer
    if(mImage.GetBuffer(buff))
    {
        ///get subimage dimension
        unsigned long subDimX, subDimY, squareDim;
        if(sqrt(nb_slice_per_mosa)==floor(sqrt(nb_slice_per_mosa)))
        {
            squareDim = (unsigned short) floor(sqrt(nb_slice_per_mosa));
        }
        else
        {
            squareDim = (unsigned short) floor(sqrt(nb_slice_per_mosa)) + 1;
        }
        //cout << "before getNumberOfColumns" << endl;
        subDimX = getNumberOfColumns()/squareDim;
        //cout << "after getNumberOfColumns" << endl;
        subDimY = getNumberOfRows()/squareDim;
        //cout << "after getNumberOfRows" << endl;

        ///allocate storage
        rawDataX = new rawData<X>(getPixelRepresentation(), 2+1, subDimY, subDimX, nb_slice_per_mosa);
        if(rawDataX==NULL)
        {
            delete buff;
            return NULL;
        }
        if(rawDataX->raw1D==NULL && rawDataX->raw2D==NULL && rawDataX->raw3D==NULL && rawDataX->raw4D==NULL)
        {
            delete buff;
            delete rawDataX;
            return NULL;
        }

        //reorganize buffer into rawData structure
        X* bufferX = (X*) buff;
        for(unsigned long i=0 ; i<rawDataX->DimX ; i++)
        {
            for(unsigned long j=0 ; j<rawDataX->DimY ; j++)
            {
                for(unsigned long k=0 ; k<rawDataX->DimZ ; k++)
                {
                    ///store data
                    rawDataX->raw3D[i][j][k] = bufferX[rawDataX->getIndex3DUnmosaic(i,j,k)];
                }
            }
        }



    }

    //cout << "before getPixelDimX" << endl;
    rawDataX->pixDimX = getPixelDimX();
    //cout << "after getPixelDimX" << endl;
    rawDataX->pixDimY = getPixelDimY();
    rawDataX->pixDimZ = getPixelDimZ(); ///may not be right
    rawDataX->pixDimT = 1.0;
    delete buff;
    return rawDataX;
}



template<typename X> rawData<double>* C_readDicomDiffData::conver2double(rawData<X>* tempRaw)
{
    //copy data into new container

    //allocate rawData
    rawData<double>* mRawData = new rawData<double>(DICOM_DOUBLE, tempRaw->numDim, tempRaw->DimX, tempRaw->DimY, tempRaw->DimZ, tempRaw->DimT);
    if(mRawData==NULL)
    {
        return NULL;
    }
    if(mRawData->raw1D==NULL && mRawData->raw2D==NULL && mRawData->raw3D==NULL && mRawData->raw4D==NULL)
    {
        delete mRawData;
        return NULL;
    }

    if(tempRaw->numDim==1)//(mImage.GetNumberOfDimensions()==1)
    {

        //copy mRawData->raw1D
        for(unsigned int i=0 ; i<tempRaw->DimX ; i++)
        {
            mRawData->raw1D[i] = (double) tempRaw->raw1D[i];
        }
    }
    else if(tempRaw->numDim==2)//(mImage.GetNumberOfDimensions()==2)
    {
        //copy mRawData->raw2D
        for(unsigned int i=0 ; i<tempRaw->DimX ; i++)
        {
            for(unsigned int j=0 ; j<tempRaw->DimY ; j++)
            {
                mRawData->raw2D[i][j] = (double) tempRaw->raw2D[i][j];
            }

        }
    }
    else if(tempRaw->numDim==3)//if(mImage.GetNumberOfDimensions()==3)
    {
        //copy mRawData->raw3D
        for(unsigned int i=0 ; i<tempRaw->DimX ; i++)
        {
            for(unsigned int j=0 ; j<tempRaw->DimY ; j++)
            {
                for(unsigned int k=0 ; k<tempRaw->DimZ ; k++)
                {
                    mRawData->raw3D[i][j][k] = (double) tempRaw->raw3D[i][j][k];
                }
            }

        }
    }
    else if(tempRaw->numDim==4)//if(mImage.GetNumberOfDimensions()==4)
    {
        //copy mRawData->raw4D
        for(unsigned int i=0 ; i<tempRaw->DimX ; i++)
        {
            for(unsigned int j=0 ; j<tempRaw->DimY ; j++)
            {
                for(unsigned int k=0 ; tempRaw->DimZ ; k++)
                {
                    for(unsigned int l=0 ; l<tempRaw->DimT ; l++)
                    {
                        mRawData->raw4D[i][j][k][l] = (double) tempRaw->raw4D[i][j][k][l];
                    }
                }
            }

        }
    }
    else
    {
        delete mRawData;
        mRawData = NULL;
    }
    if(mRawData!=NULL)
    {
        mRawData->pixDimX = tempRaw->pixDimX;
        mRawData->pixDimY = tempRaw->pixDimY;
        mRawData->pixDimZ = tempRaw->pixDimZ;
        mRawData->pixDimT = tempRaw->pixDimT;
    }
    return mRawData;
}
