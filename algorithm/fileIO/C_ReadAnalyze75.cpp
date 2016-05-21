#include <C_ReadAnalyze75.h>


C_ReadAnalyze75::C_ReadAnalyze75()
{
    //ctor
}

C_ReadAnalyze75::~C_ReadAnalyze75()
{
    //dtor
}
unsigned short C_ReadAnalyze75::getDataRepresentation(dsr* dsrPtr)
{
    unsigned short datatype = dsrPtr->dime.datatype;
    //cout << datatype << endl;
    //cout << dsrPtr->dime.bitpix << endl;
    unsigned short DICOM_TYPE=DICOM_ERROR;
    if(datatype==DT_BINARY)
    {
        //readBinary data? must check bitdepth bitpix. Should be 8 (char) or 32 (int)
        if(dsrPtr->dime.bitpix==8) DICOM_TYPE = DICOM_CHAR;
        else if(dsrPtr->dime.bitpix==32) DICOM_TYPE = DICOM_LONG;
        else DICOM_TYPE = DICOM_ERROR;
    }
    else if(datatype==DT_UNSIGNED_CHAR)
    {
        // retrun an unsigned char array
        if(dsrPtr->dime.bitpix==8) DICOM_TYPE = DICOM_UCHAR;
        else DICOM_TYPE = DICOM_ERROR;
    }
    else if(datatype==DT_SIGNED_SHORT)
    {
        //check bitpix==2 short int or error
        if(dsrPtr->dime.bitpix==16) DICOM_TYPE = DICOM_SHORT;
        else DICOM_TYPE = DICOM_ERROR;
    }
    else if(datatype==DT_SIGNED_INT)
    {
        //check bitpix => 2 means short int, 4 means long int
        if(dsrPtr->dime.bitpix==16) DICOM_TYPE = DICOM_SHORT;
        else if(dsrPtr->dime.bitpix==32) DICOM_TYPE = DICOM_LONG;
        else DICOM_TYPE = DICOM_ERROR;
    }
    else if(datatype==DT_FLOAT) DICOM_TYPE = DICOM_FLOAT;
    else if(datatype==DT_COMPLEX) DICOM_TYPE = DICOM_FLOAT;
    else if(datatype==DT_DOUBLE) DICOM_TYPE = DICOM_DOUBLE;

    return DICOM_TYPE;
}

rawData<double>* C_ReadAnalyze75::getRawData(const char* filenameData)
{
    ///read analyze header
    dsr* dsrPtr = readAnalyzeInfo(getHdrName(filenameData).data());

//    cout << filenameData << " numDim = " << dsrPtr->dime.dim[0] << " " << dsrPtr->dime.dim[1] << " " << dsrPtr->dime.dim[2] << " " << dsrPtr->dime.dim[3] << " " << dsrPtr->dime.dim[4] << endl;
//    getchar();

    ///get dicom type (pixel representation for rawData)
    unsigned short DICOM_TYPE=getDataRepresentation(dsrPtr);
//    cout << "data representation: " << DICOM_TYPE << endl;

    if(DICOM_TYPE==DICOM_ERROR)
    {
        delete dsrPtr;
        return NULL;
    }
    rawData<double>* mRawData=NULL;
    switch(DICOM_TYPE)
    {
        case DICOM_BOOL:
        {
            rawData<bool>* tempRaw = getRawData_<bool>(filenameData, dsrPtr);
            if(tempRaw!=NULL)
            {
                mRawData = conver2double(tempRaw);
                delete tempRaw;
            }
            break;
        }
        case DICOM_CHAR:
        {
            rawData<char>* tempRaw = getRawData_<char>(filenameData, dsrPtr);
            if(tempRaw!=NULL)
            {
                mRawData = conver2double(tempRaw);
                delete tempRaw;
            }
            break;
        }
        case DICOM_SHORT:
        {
            rawData<short>* tempRaw = getRawData_<short>(filenameData, dsrPtr);
            if(tempRaw!=NULL)
            {
                mRawData = conver2double(tempRaw);
                delete tempRaw;
            }
            break;
        }
        case DICOM_LONG:
        {
            rawData<long>* tempRaw = getRawData_<long>(filenameData, dsrPtr);
            if(tempRaw!=NULL)
            {
                mRawData = conver2double(tempRaw);
                delete tempRaw;
            }
            break;
        }
        case DICOM_UCHAR:
        {
            rawData<unsigned char>* tempRaw = getRawData_<unsigned char>(filenameData, dsrPtr);
            if(tempRaw!=NULL)
            {
                mRawData = conver2double(tempRaw);
                delete tempRaw;
            }
            break;
        }
        case DICOM_USHORT:
        {
            rawData<unsigned short>* tempRaw = getRawData_<unsigned short>(filenameData, dsrPtr);
            if(tempRaw!=NULL)
            {
                mRawData = conver2double(tempRaw);
                delete tempRaw;
            }
            break;
        }
        case DICOM_ULONG:
        {
            rawData<unsigned long>* tempRaw = getRawData_<unsigned long>(filenameData, dsrPtr);
            if(tempRaw!=NULL)
            {
                mRawData = conver2double(tempRaw);
                delete tempRaw;
            }
            break;
        }
        case DICOM_FLOAT:
        {
            rawData<float>* tempRaw = getRawData_<float>(filenameData, dsrPtr);
            //cout << "tempRaw " << tempRaw << endl;
            if(tempRaw!=NULL)
            {
                mRawData = conver2double(tempRaw);
                delete tempRaw;
            }
            break;
        }
        case DICOM_DOUBLE:
        {
            rawData<double>* tempRaw = getRawData_<double>(filenameData, dsrPtr);
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

    delete dsrPtr;

    return mRawData;

}

rawData<double>* C_ReadAnalyze75::getUnmosaicData(const char* filenameData, unsigned short nb_slice_per_mosa)
{
    //cout << "1" << endl;
    ///read analyze header
    dsr* dsrPtr = readAnalyzeInfo(getHdrName(filenameData).data());

    ///get dicom type (pixel representation for rawData)
    //cout << "2" << endl;
    unsigned short DICOM_TYPE=getDataRepresentation(dsrPtr);
    //cout << "data representation: " << DICOM_TYPE << endl;

    if(DICOM_TYPE==DICOM_ERROR)
    {
        //cout << "3" << endl;
        delete dsrPtr;
        return NULL;
    }
    rawData<double>* mRawData=NULL;
    switch(DICOM_TYPE)
    {
        case DICOM_BOOL:
        {
            rawData<bool>* tempRaw = getUnmosaicData_<bool>(filenameData, dsrPtr, nb_slice_per_mosa);
            if(tempRaw!=NULL)
            {
                mRawData = conver2double(tempRaw);
                delete tempRaw;
            }
            break;
        }
        case DICOM_CHAR:
        {
            rawData<char>* tempRaw = getUnmosaicData_<char>(filenameData, dsrPtr, nb_slice_per_mosa);
            if(tempRaw!=NULL)
            {
                mRawData = conver2double(tempRaw);
                delete tempRaw;
            }
            break;
        }
        case DICOM_SHORT:
        {
            rawData<short>* tempRaw = getUnmosaicData_<short>(filenameData, dsrPtr, nb_slice_per_mosa);
            if(tempRaw!=NULL)
            {
                mRawData = conver2double(tempRaw);
                delete tempRaw;
            }
            break;
        }
        case DICOM_LONG:
        {
            rawData<long>* tempRaw = getUnmosaicData_<long>(filenameData, dsrPtr, nb_slice_per_mosa);
            if(tempRaw!=NULL)
            {
                mRawData = conver2double(tempRaw);
                delete tempRaw;
            }
            break;
        }
        case DICOM_UCHAR:
        {
            rawData<unsigned char>* tempRaw = getUnmosaicData_<unsigned char>(filenameData, dsrPtr, nb_slice_per_mosa);
            if(tempRaw!=NULL)
            {
                mRawData = conver2double(tempRaw);
                delete tempRaw;
            }
            break;
        }
        case DICOM_USHORT:
        {
            rawData<unsigned short>* tempRaw = getUnmosaicData_<unsigned short>(filenameData, dsrPtr, nb_slice_per_mosa);
            if(tempRaw!=NULL)
            {
                mRawData = conver2double(tempRaw);
                delete tempRaw;
            }
            break;
        }
        case DICOM_ULONG:
        {
            rawData<unsigned long>* tempRaw = getUnmosaicData_<unsigned long>(filenameData, dsrPtr, nb_slice_per_mosa);
            if(tempRaw!=NULL)
            {
                mRawData = conver2double(tempRaw);
                delete tempRaw;
            }
            break;
        }
        case DICOM_FLOAT:
        {
            rawData<float>* tempRaw = getUnmosaicData_<float>(filenameData, dsrPtr, nb_slice_per_mosa);
            //cout << "tempRaw " << tempRaw << endl;
            if(tempRaw!=NULL)
            {
                mRawData = conver2double(tempRaw);
                delete tempRaw;
            }
            break;
        }
        case DICOM_DOUBLE:
        {
            rawData<double>* tempRaw = getUnmosaicData_<double>(filenameData, dsrPtr, nb_slice_per_mosa);
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

    delete dsrPtr;

    return mRawData;

}

dsr* C_ReadAnalyze75::readAnalyzeInfo(const char* filenameHDR)
{
    FILE *f=fopen(filenameHDR,"rb");
    if(f==NULL)
    {
        return (dsr*) NULL;
    }

    //fread();
	fseek(f,0,SEEK_END); //set current positon to the end of the file
	unsigned long n = ftell(f)/sizeof(dsr); //tell the size of the file divide by length of structure = nb of structure
	if(n!=1)
	{
	    cout << "hdr pb" << endl;
	    fclose(f);
	    return (dsr*) NULL;
	}
	fseek(f,0,SEEK_SET); //set current possition to the begining of the file
	dsr* dsrPtr = new dsr[n];//*sizeof(PointWithNeighbors)];
	fread(dsrPtr,sizeof(dsr),n,f); //read into buffer
    fclose(f);
    ///check for a special case
    if(dsrPtr->dime.dim[0]==4)
    {
        if(dsrPtr->dime.dim[4]<2) dsrPtr->dime.dim[0]=3;
    }
    return dsrPtr;
}

char* C_ReadAnalyze75::readAnalyzeData(const char* filenameData, dsr* dsrPtr, long int offsetElement, long int nbElement)
{
    FILE* f = fopen(filenameData, "rb");
    if(f==NULL)
    {
        return NULL;
    }

    ///get size of elements to read
    unsigned short datatype = dsrPtr->dime.datatype;
    unsigned short SIZE_OF_ELTS = 0;
    if(datatype==DT_BINARY)
    {
        //readBinary data? must check bitdepth bitpix. Should be 8 (char) or 32 (int)
        if(dsrPtr->dime.bitpix==8) SIZE_OF_ELTS = 1;
        else if(dsrPtr->dime.bitpix==32) SIZE_OF_ELTS = 4;
        else SIZE_OF_ELTS = 0;
    }
    else if(datatype==DT_UNSIGNED_CHAR)
    {
        // retrun an unsigned char array
        if(dsrPtr->dime.bitpix==8) SIZE_OF_ELTS = 1;
        else SIZE_OF_ELTS = 0;
    }
    else if(datatype==DT_SIGNED_SHORT)
    {
        //check bitpix==2 short int or error
        if(dsrPtr->dime.bitpix==16) SIZE_OF_ELTS = 2;
        else SIZE_OF_ELTS = 0;
    }
    else if(datatype==DT_SIGNED_INT)
    {
        //check bitpix => 2 means short int, 4 means long int
        if(dsrPtr->dime.bitpix==16) SIZE_OF_ELTS = 2;
        else if(dsrPtr->dime.bitpix==32) SIZE_OF_ELTS = 4;
        else SIZE_OF_ELTS = 0;
    }
    else if(datatype==DT_FLOAT) SIZE_OF_ELTS = 4;
    else if(datatype==DT_COMPLEX) SIZE_OF_ELTS = 4;
    else if(datatype==DT_DOUBLE) SIZE_OF_ELTS = 8;

    //check dim
    fseek(f,0,SEEK_END);
    unsigned long n = ftell(f)/SIZE_OF_ELTS;
    unsigned long dim = getNumel(dsrPtr, dsrPtr->dime.dim[0]);
    if(n!= dim)
    {
        cout << "dimension pb numel should be " << dim << " but numel in file is " << n << " file: " << filenameData << " is corrupted!" << endl;
        fclose(f);
        return NULL;
    }
    if(nbElement!=-1)
    {
        if( ((unsigned long) offsetElement)+((unsigned long) nbElement)>n)
        {
            cout << "element to read are out of range" << endl;
            fclose(f);
            return NULL;
        }
        n = nbElement;
    }
    fseek(f,offsetElement*SIZE_OF_ELTS,SEEK_SET);//fseek(f,sizeof(unsigned char)*offsetElement,SEEK_SET);

    //read data
    char* data = new char[n*SIZE_OF_ELTS];
	fread(data,SIZE_OF_ELTS,n,f); //read into buffer

	//close file
    fclose(f);
    return data;
}



template<typename T> rawData<T>* C_ReadAnalyze75::getRawData_(const char* filenameData, dsr* dsrPtr)
{
    ///read analyze header
    //dsr* dsrPtr = readAnalyzeInfo(getHdrName(filenameData).data());

    ///get dicom type (pixel representation for rawData)
    unsigned short DICOM_TYPE=getDataRepresentation(dsrPtr);
    //cout << "re-DICOM_TYPE" << DICOM_TYPE << endl;

    if(DICOM_TYPE==DICOM_ERROR)
    {
        delete dsrPtr;
        return NULL;
    }

    rawData<T>* rawDataT = NULL;
    if(dsrPtr->dime.dim[0]==2)
    {
        ///load data into char buffer
        char* buff = readAnalyzeData(getImgName(filenameData).data(), dsrPtr, 0, (dsrPtr->dime.dim[1])*(dsrPtr->dime.dim[2]));
        if(buff==NULL) return NULL;

        ///cast to T-buffer
        T* bufferT = (T*) buff;

        ///allocate rawData
        rawDataT = new rawData<T>(DICOM_TYPE, 2, dsrPtr->dime.dim[1], dsrPtr->dime.dim[2]);
        rawDataT->pixDimX = dsrPtr->dime.pixdim[1];
        rawDataT->pixDimY = dsrPtr->dime.pixdim[2];
        rawDataT->pixDimZ = dsrPtr->dime.pixdim[3];
        rawDataT->pixDimT = dsrPtr->dime.pixdim[4];

        ///fill in rawData
        for(unsigned long i=0 ; i<(unsigned long) dsrPtr->dime.dim[1] ; i++) //X dim
        {
            for(unsigned long j=0 ; j<(unsigned long) dsrPtr->dime.dim[2] ; j++) //Y dim
            {
                rawDataT->raw2D[i][j] = bufferT[j*((unsigned long) dsrPtr->dime.dim[1]) + i];
            }
        }

        ///destroy
        delete buff;
    }
    else if(dsrPtr->dime.dim[0]==3)
    {
        ///load data into char buffer
        char* buff = readAnalyzeData(getImgName(filenameData).data(), dsrPtr, 0, (dsrPtr->dime.dim[1])*(dsrPtr->dime.dim[2])*(dsrPtr->dime.dim[3]));
        if(buff==NULL) return NULL;

        ///cast to T-buffer
        T* bufferT = (T*) buff;

        ///allocate rawData
        rawDataT = new rawData<T>(DICOM_TYPE, 3, dsrPtr->dime.dim[1], dsrPtr->dime.dim[2], dsrPtr->dime.dim[3]);
        rawDataT->pixDimX = dsrPtr->dime.pixdim[1];
        rawDataT->pixDimY = dsrPtr->dime.pixdim[2];
        rawDataT->pixDimZ = dsrPtr->dime.pixdim[3];
        rawDataT->pixDimT = dsrPtr->dime.pixdim[4];

        ///fill in rawData
        for(unsigned long i=0 ; i<(unsigned long) dsrPtr->dime.dim[1] ; i++) //X dim
        {
            for(unsigned long j=0 ; j<(unsigned long) dsrPtr->dime.dim[2] ; j++) //Y dim
            {
                for(unsigned long k=0 ; k<(unsigned long) dsrPtr->dime.dim[3] ; k++) //Z dim
                {
                    rawDataT->raw3D[i][j][k] = bufferT[k*((unsigned long) dsrPtr->dime.dim[2])*((unsigned long) dsrPtr->dime.dim[1]) + j*((unsigned long) dsrPtr->dime.dim[1]) + i];
                }
            }
        }

        ///destroy
        delete buff;
    }
    else if(dsrPtr->dime.dim[0]==4)
    {
        ///allocate rawData
        rawDataT = new rawData<T>(DICOM_TYPE, 4, dsrPtr->dime.dim[1], dsrPtr->dime.dim[2], dsrPtr->dime.dim[3], dsrPtr->dime.dim[4]);
        rawDataT->pixDimX = dsrPtr->dime.pixdim[1];
        rawDataT->pixDimY = dsrPtr->dime.pixdim[2];
        rawDataT->pixDimZ = dsrPtr->dime.pixdim[3];
        rawDataT->pixDimT = dsrPtr->dime.pixdim[4];

        ///get size of elements to read
//        unsigned short datatype = dsrPtr->dime.datatype;
//        unsigned short SIZE_OF_ELTS = 0;
//        if(datatype==DT_BINARY)
//        {
//            //readBinary data? must check bitdepth bitpix. Should be 8 (char) or 32 (int)
//            if(dsrPtr->dime.bitpix==8) SIZE_OF_ELTS = 1;
//            else if(dsrPtr->dime.bitpix==32) SIZE_OF_ELTS = 4;
//            else SIZE_OF_ELTS = 0;
//        }
//        else if(datatype==DT_UNSIGNED_CHAR)
//        {
//            // retrun an unsigned char array
//            if(dsrPtr->dime.bitpix==8) SIZE_OF_ELTS = 1;
//            else SIZE_OF_ELTS = 0;
//        }
//        else if(datatype==DT_SIGNED_SHORT)
//        {
//            //check bitpix==2 short int or error
//            if(dsrPtr->dime.bitpix==16) SIZE_OF_ELTS = 2;
//            else SIZE_OF_ELTS = 0;
//        }
//        else if(datatype==DT_SIGNED_INT)
//        {
//            //check bitpix => 2 means short int, 4 means long int
//            if(dsrPtr->dime.bitpix==16) SIZE_OF_ELTS = 2;
//            else if(dsrPtr->dime.bitpix==32) SIZE_OF_ELTS = 4;
//            else SIZE_OF_ELTS = 0;
//        }
//        else if(datatype==DT_FLOAT) SIZE_OF_ELTS = 4;
//        else if(datatype==DT_COMPLEX) SIZE_OF_ELTS = 4;
//        else if(datatype==DT_DOUBLE) SIZE_OF_ELTS = 8;

//        cout << "SIZE_OF_ELTS: " << SIZE_OF_ELTS << endl;

        ///fill in rawData
        for(unsigned long t=0 ; t<(unsigned long) dsrPtr->dime.dim[4] ; t++) //grad dim
        {
            ///load data into char buffer
            unsigned long volumeData = (dsrPtr->dime.dim[1])*(dsrPtr->dime.dim[2])*(dsrPtr->dime.dim[3]);
            char* buff = readAnalyzeData(getImgName(filenameData).data(), dsrPtr, t*volumeData, volumeData);
            if(buff==NULL) return NULL;

            ///cast to T-buffer
            T* bufferT = (T*) buff;
            for(unsigned long i=0 ; i<(unsigned long) dsrPtr->dime.dim[1] ; i++) //X dim
            {
                for(unsigned long j=0 ; j<(unsigned long) dsrPtr->dime.dim[2] ; j++) //Y dim
                {
                    for(unsigned long k=0 ; k<(unsigned long) dsrPtr->dime.dim[3] ; k++) //Z dim
                    {
                        rawDataT->raw4D[i][j][k][t] = bufferT[k*((unsigned long) dsrPtr->dime.dim[2])*((unsigned long) dsrPtr->dime.dim[1]) + j*((unsigned long) dsrPtr->dime.dim[1]) + i];
                    }
                }
            }
            delete buff;
        }
    }
    return rawDataT;
}


template<typename T> rawData<T>* C_ReadAnalyze75::getUnmosaicData_(const char* filenameData, dsr* dsrPtr, unsigned short nb_slice_per_mosa)
{
    ///get dicom type (pixel representation for rawData)
    unsigned short DICOM_TYPE=getDataRepresentation(dsrPtr);
//    cout << "DICOM_TYPE" << DICOM_TYPE << endl;

    if(DICOM_TYPE==DICOM_ERROR)
    {
        delete dsrPtr;
        return NULL;
    }

    rawData<T>* rawDataT = NULL;
    if(dsrPtr->dime.dim[0]==2)
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
        rawDataT = getRawData_<T>(filenameData, dsrPtr);
        rawDataT->swapXY2D();
        T tempT;
        for(unsigned long i=0 ; i<rawDataT->DimX ; i++)
        {
            for(unsigned long j=0 ; j<rawDataT->DimY/2 ; j++)
            {
                tempT = rawDataT->raw2D[i][rawDataT->DimY-1-j];
                rawDataT->raw2D[i][rawDataT->DimY-1-j] = rawDataT->raw2D[i][j];
                rawDataT->raw2D[i][j] = tempT;
            }
        }
        for(unsigned long i=0 ; i<rawDataT->DimX/2 ; i++)
        {
            for(unsigned long j=0 ; j<rawDataT->DimY ; j++)
            {
                tempT = rawDataT->raw2D[rawDataT->DimX-1-i][j];
                rawDataT->raw2D[rawDataT->DimX-1-i][j] = rawDataT->raw2D[i][j];
                rawDataT->raw2D[i][j] = tempT;
            }
        }
        for(unsigned long i=0 ; i<rawDataT->DimX ; i++)
        {
            for(unsigned long j=0 ; j<rawDataT->DimY/2 ; j++)
            {
                tempT = rawDataT->raw2D[i][rawDataT->DimY-1-j];
                rawDataT->raw2D[i][rawDataT->DimY-1-j] = rawDataT->raw2D[i][j];
                rawDataT->raw2D[i][j] = tempT;
            }
        }
//        FILE* f = fopen("blu.txt","w");
//        cout << filenameData << endl;
//        for(unsigned long i=0 ; i<rawDataT->DimX ; i++)
//        {
//            for(unsigned long j=0 ; j<rawDataT->DimY ; j++)
//            {
//                fprintf(f,"%i ", rawDataT->raw2D[i][j]);
//            }
//            fprintf(f,"\n");
//        }
//        fclose(f);
//        getchar();

        //create one copy for each slice
        rawData<T>** copySliceRawData = new rawData<T>*[nb_slice_per_mosa];
        for(unsigned short i=0 ; i<nb_slice_per_mosa ; i++)
        {
            copySliceRawData[i] = new rawData<T>(DICOM_TYPE, 2, rawDataT->DimX, rawDataT->DimY, rawDataT->DimZ, rawDataT->DimT);
            *(copySliceRawData[i]) = *rawDataT;
        }

        subDimX = rawDataT->DimX/squareDim;
        subDimY = rawDataT->DimY/squareDim;

        delete rawDataT;
        rawDataT = new rawData<T>(DICOM_TYPE, 3, subDimX, subDimY, nb_slice_per_mosa);

        unsigned short zIndex = 0;
        for(unsigned long i=0 ; i<squareDim ; i++)
        {
            for(unsigned long j=0 ; j<squareDim ; j++)
            {
                copySliceRawData[zIndex]->crop2D(i*subDimX, (i+1)*subDimX-1, j*subDimY, (j+1)*subDimY-1);
                zIndex++;
                if(zIndex==nb_slice_per_mosa)
                {
                    i=squareDim;
                    j=squareDim;
                }
            }
        }
        for(unsigned short k=0 ; k<nb_slice_per_mosa ; k++)
        {
            for(unsigned long i=0 ; i<subDimX ; i++)
            {
                for(unsigned long j=0 ; j<subDimY ; j++)
                {
                    rawDataT->raw3D[i][j][k] = copySliceRawData[k]->raw2D[i][j];
                }
            }
            delete copySliceRawData[k];
        }
        delete copySliceRawData;

//        rawDataT = rawDataT;


//        ///load data into char buffer
//        char* buff = readAnalyzeData(getImgName(filenameData).data(), dsrPtr, 0, (dsrPtr->dime.dim[1])*(dsrPtr->dime.dim[2]));
//        if(buff==NULL) return NULL;
//
//        ///cast to T-buffer
//        T* bufferT = (T*) buff;
//
//        ///allocate rawData
//        rawDataT = new rawData<T>(DICOM_TYPE, 2+1, subDimX, subDimY, nb_slice_per_mosa);
//        if(rawDataT==NULL)
//        {
//            delete buff;
//            return NULL;
//        }
//        if(rawDataT->raw1D==NULL && rawDataT->raw2D==NULL && rawDataT->raw3D==NULL && rawDataT->raw4D==NULL)
//        {
//            delete buff;
//            delete rawDataT;
//            return NULL;
//        }
//        rawDataT->pixDimX = dsrPtr->dime.pixdim[1];
//        rawDataT->pixDimY = dsrPtr->dime.pixdim[2];
//        rawDataT->pixDimZ = dsrPtr->dime.pixdim[3];
//        rawDataT->pixDimT = dsrPtr->dime.pixdim[4];
//
//        ///fill in rawData
//        for(unsigned long i=0 ; i<rawDataT->DimX ; i++) //X dim
//        {
//            for(unsigned long j=0 ; j<rawDataT->DimY ; j++) //Y dim
//            {
//                for(unsigned long k=0 ; k<rawDataT->DimZ ; k++)
//                {
//                    rawDataT->raw3D[i][j][k] = bufferT[(dsrPtr->dime.dim[1])*(dsrPtr->dime.dim[2])-1-rawDataT->getIndex3DUnmosaic(i,j,k)];
//                }
//            }
//        }
//
//        ///destroy
//        delete buff;
    }
    return rawDataT;
}



template<typename T> rawData<double>* C_ReadAnalyze75::conver2double(rawData<T>* tempRaw)
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
                for(unsigned int k=0 ; k<tempRaw->DimZ ; k++)
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

vector< vector<double> > C_ReadAnalyze75::getGradientDirectionFormTextFile(string gradFileName)
{
    ///assume that gradient directions are stored on lines (1 line = 3 double). The components are separated by one space
    FILE* file = fopen(gradFileName.data(), "r");
    char line[256];

    vector<double> dGradDir(3,0.0);
    vector< vector<double> > gradDir;
    while (fgets(line, 256, file) != NULL)
    {
        string tempBuf = line;
        unsigned int found1 = tempBuf.find_first_of(" ");
        unsigned int found2 = tempBuf.find_first_of(" ", found1+1);
        if(found1!=string::npos && found2!=string::npos)
        {
            dGradDir.at(0) = atof(tempBuf.substr(0, found1).data());
            dGradDir.at(1) = atof(tempBuf.substr(found1+1, found2 -found1).data());
            dGradDir.at(2) = atof(tempBuf.substr(found2).data());
            gradDir.push_back(dGradDir);
        }
    }
    fclose(file);
    return gradDir;
}


bool C_ReadAnalyze75::saveAna(const char* filenameAna, dsr* ptrDSR, rawData<double>* data)
{
    //check or fill dimensions
    bool r = writeAnalyzeInfo(getHdrName(filenameAna).data(), ptrDSR);
    return (r && writeAnalyzeData(getImgName(filenameAna).data(), data));
}
bool C_ReadAnalyze75::writeAnalyzeInfo(const char* filenameHDR, dsr* ptrDSR)
{
    FILE *f=fopen(filenameHDR,"wb");
    if(f==NULL)
    {
        return false;
    }
	size_t ret = fwrite((char*) ptrDSR, sizeof(char), sizeof(dsr), f);
	fclose(f);
	if(ret!=sizeof(dsr))
	{
	    return false;
	}
	else
	{
	    return true;
	}
}
bool C_ReadAnalyze75::writeAnalyzeData(const char* filenameData, rawData<double>* data)
{
    FILE *f=fopen(filenameData,"wb");
    if(f==NULL)
    {
        return false;
    }
    size_t ret = 0.0;
    if(data->numDim==1)
    {
        ret = fwrite((char*) data->raw1D, sizeof(char), sizeof(double)*data->DimX, f);
    }
    else if(data->numDim==2)
    {
        double* tempBuff = new double[data->DimX];
        ret = 0;
        for(unsigned long j=0 ; j<data->DimY ; j++)
        {
            for(unsigned long i=0 ; i<data->DimX ; i++)
            {
                tempBuff[i] = data->raw2D[i][j];
            }
            ret += fwrite((char*) tempBuff, sizeof(char), sizeof(double)*data->DimX, f);
        }
        delete tempBuff;
    }
    else if(data->numDim==3)
    {
        double* tempBuff = new double[data->DimX];
        ret = 0;
        for(unsigned long k=0 ; k<data->DimZ ; k++)
        {
            for(unsigned long j=0 ; j<data->DimY ; j++)
            {
                for(unsigned long i=0 ; i<data->DimX ; i++)
                {
                    tempBuff[i] = data->raw3D[i][j][k];
                }
                ret += fwrite((char*) tempBuff, sizeof(char), sizeof(double)*data->DimX, f);
            }
        }
        delete tempBuff;
    }
    else if(data->numDim==4)
    {
        double* tempBuff = new double[data->DimX];
        ret = 0;
        for(unsigned long l=0 ; l<data->DimT ; l++)
        {
            for(unsigned long k=0 ; k<data->DimZ ; k++)
            {
                for(unsigned long j=0 ; j<data->DimY ; j++)
                {
                    for(unsigned long i=0 ; i<data->DimX ; i++)
                    {
                        tempBuff[i] = data->raw4D[i][j][k][l];
                    }
                    ret += fwrite((char*) tempBuff, sizeof(char), sizeof(double)*data->DimX, f);
                }
            }
        }
        delete tempBuff;
    }
	fclose(f);
	if(ret!=data->numel())
	{
	    return false;
	}
	else
	{
	    return true;
	}
    return true;
}
bool C_ReadAnalyze75::getGradientDirectionFormTextFile(string gradFileName, vector< vector<double> > grad)
{
    FILE* file = fopen(getDatName(gradFileName).data(), "w");
    if(file==NULL) return false;
    for(unsigned int i=0 ; i<grad.size() ; i++)
    {
        for(unsigned int j=0 ; j<grad.at(i).size() ; j++)
        {
            fprintf(file, "%f ", grad.at(i).at(j));
        }
        fprintf(file,"\n");
    }
    fclose(file);
    return true;
}

short int C_ReadAnalyze75::getDim(dsr* dsrPtr)
{
    if(dsrPtr==NULL)
    {
        return 0;
    }
    return dsrPtr->dime.dim[0];
}

long int C_ReadAnalyze75::getNumel(dsr* dsrPtr, short int numDim)
{
    long int Nel = 1;
    //short int N = getDim(dsrPtr);
    for(short int i=1 ; i<=numDim ; i++)
    {
        if(dsrPtr->dime.dim[i]>0)
        {
            Nel *= dsrPtr->dime.dim[i];
        }
    }
    return Nel;
}


string C_ReadAnalyze75::getImgName(string name)
{
    if(getSuffix(name).compare(".img")==0)
    {
        string newName = name;
        return newName;
    }
    string base =  getBaseName(name);
    base.append(".img");
    return base;
}
string C_ReadAnalyze75::getHdrName(string name)
{
    if(getSuffix(name).compare(".hdr")==0)
    {
        string newName = name;
        return newName;
    }
    string base =  getBaseName(name);
    base.append(".hdr");
    return base;
}
string C_ReadAnalyze75::getDatName(string name)
{
    if(getSuffix(name).compare(".dat")==0)
    {
        string newName = name;
        return newName;
    }
    string base =  getBaseName(name);
    base.append(".dat");
    return base;
}
string C_ReadAnalyze75::getBaseName(string name)
{
    string newName = name;
    size_t n = newName.find_last_of(".");
    if(n==string::npos)
    {
        return newName;
    }
    return newName.substr(0,n);
}

string C_ReadAnalyze75::getSuffix(string name)
{
    if(true)
    {
        size_t n = name.find_last_of(".");
        if(n==string::npos)
        {
            return "";
        }
        return name.substr(n);
    }
    else
    {
        return name.substr(name.length()-5);
    }
}
