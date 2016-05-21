#include "C_maskMaker.h"

C_maskMaker::C_maskMaker()
{
    //ctor
}

C_maskMaker::~C_maskMaker()
{
    //dtor
}


rawData<double>* C_maskMaker::createMaskFromDWIs(vector<vector<rawData<double>*> > allRawData, samplingFactor* SF, unsigned long SAMPLING_TYPE)
{
    ///check dimension
    unsigned long numDim = allRawData.at(0).at(0)->numDim;
    if(numDim!=2 && numDim!=3) return NULL;

    ///instanciate sampling tool
    C_sampling<double>* samplingTool = new C_sampling<double>();

    ///creating mask from several gradient directions
    unsigned long DX;
    unsigned long DY;
    unsigned long DZ;
    if(SF==NULL || SAMPLING_TYPE==NO_SAMPLING)
    {
        DX = allRawData.at(0).at(0)->DimX;
        DY = allRawData.at(0).at(0)->DimY;
        DZ = allRawData.at(0).at(0)->DimZ;
    }
    else
    {
        DX = ((unsigned long) ((double) allRawData.at(0).at(0)->DimX)*(SF->Xfactor));
        DY = ((unsigned long) ((double) allRawData.at(0).at(0)->DimY)*(SF->Yfactor));
        DZ = ((unsigned long) ((double) allRawData.at(0).at(0)->DimZ)*(SF->Zfactor));
    }
    rawData<double>* segRawData;// = new rawData<double>(DICOM_DOUBLE,3,DX,DY,DZ);
    if(numDim==2)
    {
        segRawData = new rawData<double>(DICOM_DOUBLE,2,DX,DY);
    }
    else
    {
        segRawData = new rawData<double>(DICOM_DOUBLE,3,DX,DY,DZ);
    }
    if(SF==NULL || SAMPLING_TYPE==NO_SAMPLING)
    {
        segRawData->pixDimX = allRawData.at(0).at(0)->pixDimX;
        segRawData->pixDimY = allRawData.at(0).at(0)->pixDimY;
        segRawData->pixDimZ = allRawData.at(0).at(0)->pixDimZ;
        segRawData->pixDimT = allRawData.at(0).at(0)->pixDimT;
    }
    else
    {
        segRawData->pixDimX = allRawData.at(0).at(0)->pixDimX*(SF->Xfactor);
        segRawData->pixDimY = allRawData.at(0).at(0)->pixDimY*(SF->Yfactor);
        segRawData->pixDimZ = allRawData.at(0).at(0)->pixDimZ*(SF->Zfactor);
        segRawData->pixDimT = allRawData.at(0).at(0)->pixDimT*(SF->Tfactor);
    }

    unsigned long NgradientDirection = allRawData.size();
    unsigned long Nrepetition = allRawData.at(0).size();
    *(segRawData) = 0.0;
    for(unsigned long i=0 ; i<NgradientDirection ; i++)
    {
        for(unsigned long j=0 ; j<Nrepetition ; j++)
        {
            if(SAMPLING_TYPE==NO_SAMPLING)
            {
                (*segRawData) += *(allRawData.at(i).at(j));
            }
            else
            {
                if(SAMPLING_TYPE==NEAREST_NEIGHBOR)
                {
                    rawData<double>* temp = samplingTool->samplingNearestNeighbor(allRawData.at(i).at(j), SF);
                    (*segRawData) += *(temp);
                    delete temp;//->~rawData();
                }
                else if(SAMPLING_TYPE==LINEAR)
                {
                    rawData<double>* temp = samplingTool->samplingLinearInterpolation(allRawData.at(i).at(j), SF);
                    (*segRawData) += *(temp);
                    delete temp;//->~rawData();
                }
                else if(SAMPLING_TYPE==KERNEL_TRICK)
                {
                    rawData<double>* temp = samplingTool->samplingGaussTrick(allRawData.at(i).at(j), SF);
                    (*segRawData) += *(temp);
                    delete temp;//->~rawData();
                }
            }
        }
    }
    (*segRawData) = (*segRawData)/((double) (NgradientDirection*Nrepetition));
    delete samplingTool;

    ///Compute kmeans
    double th = 0.0;
    if(numDim==2)
    {
        C_kmean2D* KMEAN = new C_kmean2D(segRawData->raw2D, segRawData->DimX, segRawData->DimY);
        KMEAN->run(3);
        th = KMEAN->mu.at(2) + KMEAN->sig.at(2);
        delete KMEAN;
    }
    else
    {
        C_kmean3D* KMEAN = new C_kmean3D(segRawData->raw3D, segRawData->DimX, segRawData->DimY, segRawData->DimZ);
        KMEAN->run(3);
        th = KMEAN->mu.at(2) + 0.9*KMEAN->sig.at(2);
        delete KMEAN;
    }
    *segRawData = (*segRawData)>th;

    ///remove isolated pixels (bondary condition=0)
    if(numDim==2)
    {
        //first and last row
        for(unsigned long j=0 ; j<segRawData->DimY ; j++)
        {
            segRawData->raw2D[0][j] = 0;
            segRawData->raw2D[segRawData->DimX-1][j] = 0;
        }
        //first and last column
        for(unsigned long i=1 ; i<segRawData->DimX-1 ; i++)
        {
            segRawData->raw2D[i][0] = 0;
            segRawData->raw2D[i][segRawData->DimY-1] = 0;
        }

        //inside
        unsigned long neighborCounter;
        for(unsigned long i=1 ; i<segRawData->DimX-1 ; i++)
        {
            for(unsigned long j=1 ; j<segRawData->DimY-1 ; j++)
            {
                ///if the current pixel is part of the mask
                if(segRawData->raw2D[i][j]>0.5)
                {
                    neighborCounter = 0;
                    ///count number of neighbors
                    if(segRawData->raw2D[i+1][j]>0.5) neighborCounter++;
                    if(segRawData->raw2D[i-1][j]>0.5) neighborCounter++;
                    if(segRawData->raw2D[i][j+1]>0.5) neighborCounter++;
                    if(segRawData->raw2D[i][j-1]>0.5) neighborCounter++;

                    if(segRawData->raw2D[i+1][j+1]>0.5) neighborCounter++;
                    if(segRawData->raw2D[i+1][j-1]>0.5) neighborCounter++;
                    if(segRawData->raw2D[i-1][j-1]>0.5) neighborCounter++;
                    if(segRawData->raw2D[i-1][j+1]>0.5) neighborCounter++;

                    ///if not enough neighbors, delete point from mask
                    if(neighborCounter<2) segRawData->raw2D[i][j] = 0;
                }
            }
        }
    }
    else
    {

        //plans YZ (first and last)
        for(unsigned long j=0 ; j<segRawData->DimY ; j++)
        {
            for(unsigned long k=0 ; k<segRawData->DimZ ; k++)
            {
                segRawData->raw3D[0][j][k] = 0;
                segRawData->raw3D[segRawData->DimX-1][j][k] = 0;
            }
        }
        //plans XZ (first and last)
        for(unsigned long i=0 ; i<segRawData->DimX ; i++)
        {
            for(unsigned long k=1 ; k<segRawData->DimZ-1 ; k++)
            {
                segRawData->raw3D[i][0][k] = 0;
                segRawData->raw3D[i][segRawData->DimY-1][k] = 0;
            }
        }
        //plans XY (first and last)
        for(unsigned long i=1 ; i<segRawData->DimX-1 ; i++)
        {
            for(unsigned long j=1 ; j<segRawData->DimY-1 ; j++)
            {
                segRawData->raw3D[i][j][0] = 0;
                segRawData->raw3D[i][j][segRawData->DimZ-1] = 0;
            }
        }
        //inside volume
        unsigned long neighborCounter;
        for(unsigned long i=1 ; i<segRawData->DimX-1 ; i++)
        {
            for(unsigned long j=1 ; j<segRawData->DimY-1 ; j++)
            {
                for(unsigned long k=1 ; k<segRawData->DimZ-1 ; k++)
                {
                    ///if the current pixel is part of the mask
                    if(segRawData->raw3D[i][j][k]>0.5)
                    {
                        neighborCounter = 0;
                        ///count number of neighbors
                        if(segRawData->raw3D[i+1][j][k+1]>0.5) neighborCounter++;
                        if(segRawData->raw3D[i-1][j][k+1]>0.5) neighborCounter++;
                        if(segRawData->raw3D[i][j+1][k+1]>0.5) neighborCounter++;
                        if(segRawData->raw3D[i][j-1][k+1]>0.5) neighborCounter++;
                        if(segRawData->raw3D[i][j][k+1]>0.5) neighborCounter++;
                        if(segRawData->raw3D[i+1][j+1][k+1]>0.5) neighborCounter++;
                        if(segRawData->raw3D[i+1][j-1][k+1]>0.5) neighborCounter++;
                        if(segRawData->raw3D[i-1][j-1][k+1]>0.5) neighborCounter++;
                        if(segRawData->raw3D[i-1][j+1][k+1]>0.5) neighborCounter++;

                        if(segRawData->raw3D[i+1][j][k]>0.5) neighborCounter++;
                        if(segRawData->raw3D[i-1][j][k]>0.5) neighborCounter++;
                        if(segRawData->raw3D[i][j+1][k]>0.5) neighborCounter++;
                        if(segRawData->raw3D[i][j-1][k]>0.5) neighborCounter++;
                        if(segRawData->raw3D[i+1][j+1][k]>0.5) neighborCounter++;
                        if(segRawData->raw3D[i+1][j-1][k]>0.5) neighborCounter++;
                        if(segRawData->raw3D[i-1][j-1][k]>0.5) neighborCounter++;
                        if(segRawData->raw3D[i-1][j+1][k]>0.5) neighborCounter++;

                        if(segRawData->raw3D[i+1][j][k-1]>0.5) neighborCounter++;
                        if(segRawData->raw3D[i-1][j][k-1]>0.5) neighborCounter++;
                        if(segRawData->raw3D[i][j+1][k-1]>0.5) neighborCounter++;
                        if(segRawData->raw3D[i][j-1][k-1]>0.5) neighborCounter++;
                        if(segRawData->raw3D[i][j][k-1]>0.5) neighborCounter++;
                        if(segRawData->raw3D[i+1][j+1][k-1]>0.5) neighborCounter++;
                        if(segRawData->raw3D[i+1][j-1][k-1]>0.5) neighborCounter++;
                        if(segRawData->raw3D[i-1][j-1][k-1]>0.5) neighborCounter++;
                        if(segRawData->raw3D[i-1][j+1][k-1]>0.5) neighborCounter++;

                        ///if not enough neighbors, delete point from mask
                        if(neighborCounter<2) segRawData->raw3D[i][j][k] = 0;
                    }
                }
            }
        }
    }

    return segRawData;
}


rawData<double>* C_maskMaker::createMaskFromB0(vector<rawData<double>*> allRawData, samplingFactor* SF, unsigned long SAMPLING_TYPE)
{
    ///check dimension
    unsigned long numDim = allRawData.at(0)->numDim;
    if(numDim!=2 && numDim!=3) return NULL;

    ///instanciate sampling tool
    C_sampling<double>* samplingTool = new C_sampling<double>();

    ///creating mask from several gradient directions
    unsigned long DX;
    unsigned long DY;
    unsigned long DZ;
    if(SF==NULL || SAMPLING_TYPE==NO_SAMPLING)
    {
        DX = allRawData.at(0)->DimX;
        DY = allRawData.at(0)->DimY;
        DZ = allRawData.at(0)->DimZ;
    }
    else
    {
        DX = ((unsigned long) ((double) allRawData.at(0)->DimX)*(SF->Xfactor));
        DY = ((unsigned long) ((double) allRawData.at(0)->DimY)*(SF->Yfactor));
        DZ = ((unsigned long) ((double) allRawData.at(0)->DimZ)*(SF->Zfactor));
    }
    rawData<double>* segRawData;
    if(numDim==2)
    {
        segRawData = new rawData<double>(DICOM_DOUBLE,2,DX,DY);
    }
    else
    {
        segRawData = new rawData<double>(DICOM_DOUBLE,3,DX,DY,DZ);
    }
    if(SF==NULL || SAMPLING_TYPE==NO_SAMPLING)
    {
        segRawData->pixDimX = allRawData.at(0)->pixDimX;
        segRawData->pixDimY = allRawData.at(0)->pixDimY;
        segRawData->pixDimZ = allRawData.at(0)->pixDimZ;
        segRawData->pixDimT = allRawData.at(0)->pixDimT;
    }
    else
    {
        segRawData->pixDimX = allRawData.at(0)->pixDimX*(SF->Xfactor);
        segRawData->pixDimY = allRawData.at(0)->pixDimY*(SF->Yfactor);
        segRawData->pixDimZ = allRawData.at(0)->pixDimZ*(SF->Zfactor);
        segRawData->pixDimT = allRawData.at(0)->pixDimT*(SF->Tfactor);
    }

    unsigned long Nrepetition = allRawData.size();
    *(segRawData) = 0.0;
    for(unsigned long j=0 ; j<Nrepetition ; j++)
    {
        if(SAMPLING_TYPE==NO_SAMPLING)
        {
            (*segRawData) += *(allRawData.at(j));
        }
        else
        {
            if(SAMPLING_TYPE==NEAREST_NEIGHBOR)
            {
                rawData<double>* temp = samplingTool->samplingNearestNeighbor(allRawData.at(j), SF);
                (*segRawData) += *(temp);
                delete temp;//->~rawData();
            }
            else if(SAMPLING_TYPE==LINEAR)
            {
                rawData<double>* temp = samplingTool->samplingLinearInterpolation(allRawData.at(j), SF);
                (*segRawData) += *(temp);
                delete temp;//->~rawData();
            }
            else if(SAMPLING_TYPE==KERNEL_TRICK)
            {
                rawData<double>* temp = samplingTool->samplingGaussTrick(allRawData.at(j), SF);
                (*segRawData) += *(temp);
                delete temp;//->~rawData();
            }
        }
    }
    (*segRawData) = (*segRawData)/((double) Nrepetition);
    delete samplingTool;

    ///Compute kmeans
    double thLow = segRawData->getMin(), thHigh = segRawData->getMax();
    if(numDim==2)
    {
        C_kmean2D* KMEAN = new C_kmean2D(segRawData->raw2D, segRawData->DimX, segRawData->DimY);
        KMEAN->run(3);
        thLow = KMEAN->mu.at(0);
        thHigh = KMEAN->mu.at(1);
        delete KMEAN;
    }
    else
    {
        C_kmean3D* KMEAN = new C_kmean3D(segRawData->raw3D, segRawData->DimX, segRawData->DimY, segRawData->DimZ);
        KMEAN->run(3);
//        cout << KMEAN->mu.at(0) << " " << KMEAN->mu.at(1) << " " << KMEAN->mu.at(2) << endl;
        thLow = KMEAN->mu.at(0);
        thHigh = KMEAN->mu.at(1);
        delete KMEAN;
    }
    *segRawData = (*segRawData)*((*segRawData)<thHigh);
    *segRawData = (*segRawData)>thLow;

    ///remove isolated pixels (bondary condition=0)
    if(numDim==2)
    {
        //first and last row
        for(unsigned long j=0 ; j<segRawData->DimY ; j++)
        {
            segRawData->raw2D[0][j] = 0;
            segRawData->raw2D[segRawData->DimX-1][j] = 0;
        }
        //first and last column
        for(unsigned long i=1 ; i<segRawData->DimX-1 ; i++)
        {
            segRawData->raw2D[i][0] = 0;
            segRawData->raw2D[i][segRawData->DimY-1] = 0;
        }

        //inside
        unsigned long neighborCounter;
        for(unsigned long i=1 ; i<segRawData->DimX-1 ; i++)
        {
            for(unsigned long j=1 ; j<segRawData->DimY-1 ; j++)
            {
                ///if the current pixel is part of the mask
                if(segRawData->raw2D[i][j]>0.5)
                {
                    neighborCounter = 0;
                    ///count number of neighbors
                    if(segRawData->raw2D[i+1][j]>0.5) neighborCounter++;
                    if(segRawData->raw2D[i-1][j]>0.5) neighborCounter++;
                    if(segRawData->raw2D[i][j+1]>0.5) neighborCounter++;
                    if(segRawData->raw2D[i][j-1]>0.5) neighborCounter++;

                    if(segRawData->raw2D[i+1][j+1]>0.5) neighborCounter++;
                    if(segRawData->raw2D[i+1][j-1]>0.5) neighborCounter++;
                    if(segRawData->raw2D[i-1][j-1]>0.5) neighborCounter++;
                    if(segRawData->raw2D[i-1][j+1]>0.5) neighborCounter++;

                    ///if not enough neighbors, delete point from mask
                    if(neighborCounter<2) segRawData->raw2D[i][j] = 0;
                }
            }
        }
    }
    else
    {

        //plans YZ (first and last)
        for(unsigned long j=0 ; j<segRawData->DimY ; j++)
        {
            for(unsigned long k=0 ; k<segRawData->DimZ ; k++)
            {
                segRawData->raw3D[0][j][k] = 0;
                segRawData->raw3D[segRawData->DimX-1][j][k] = 0;
            }
        }
        //plans XZ (first and last)
        for(unsigned long i=0 ; i<segRawData->DimX ; i++)
        {
            for(unsigned long k=1 ; k<segRawData->DimZ-1 ; k++)
            {
                segRawData->raw3D[i][0][k] = 0;
                segRawData->raw3D[i][segRawData->DimY-1][k] = 0;
            }
        }
        //plans XY (first and last)
        for(unsigned long i=1 ; i<segRawData->DimX-1 ; i++)
        {
            for(unsigned long j=1 ; j<segRawData->DimY-1 ; j++)
            {
                segRawData->raw3D[i][j][0] = 0;
                segRawData->raw3D[i][j][segRawData->DimZ-1] = 0;
            }
        }
        //inside volume
        unsigned long neighborCounter;
        for(unsigned long i=1 ; i<segRawData->DimX-1 ; i++)
        {
            for(unsigned long j=1 ; j<segRawData->DimY-1 ; j++)
            {
                for(unsigned long k=1 ; k<segRawData->DimZ-1 ; k++)
                {
                    ///if the current pixel is part of the mask
                    if(segRawData->raw3D[i][j][k]>0.5)
                    {
                        neighborCounter = 0;
                        ///count number of neighbors
                        if(segRawData->raw3D[i+1][j][k+1]>0.5) neighborCounter++;
                        if(segRawData->raw3D[i-1][j][k+1]>0.5) neighborCounter++;
                        if(segRawData->raw3D[i][j+1][k+1]>0.5) neighborCounter++;
                        if(segRawData->raw3D[i][j-1][k+1]>0.5) neighborCounter++;
                        if(segRawData->raw3D[i][j][k+1]>0.5) neighborCounter++;
                        if(segRawData->raw3D[i+1][j+1][k+1]>0.5) neighborCounter++;
                        if(segRawData->raw3D[i+1][j-1][k+1]>0.5) neighborCounter++;
                        if(segRawData->raw3D[i-1][j-1][k+1]>0.5) neighborCounter++;
                        if(segRawData->raw3D[i-1][j+1][k+1]>0.5) neighborCounter++;

                        if(segRawData->raw3D[i+1][j][k]>0.5) neighborCounter++;
                        if(segRawData->raw3D[i-1][j][k]>0.5) neighborCounter++;
                        if(segRawData->raw3D[i][j+1][k]>0.5) neighborCounter++;
                        if(segRawData->raw3D[i][j-1][k]>0.5) neighborCounter++;
                        if(segRawData->raw3D[i+1][j+1][k]>0.5) neighborCounter++;
                        if(segRawData->raw3D[i+1][j-1][k]>0.5) neighborCounter++;
                        if(segRawData->raw3D[i-1][j-1][k]>0.5) neighborCounter++;
                        if(segRawData->raw3D[i-1][j+1][k]>0.5) neighborCounter++;

                        if(segRawData->raw3D[i+1][j][k-1]>0.5) neighborCounter++;
                        if(segRawData->raw3D[i-1][j][k-1]>0.5) neighborCounter++;
                        if(segRawData->raw3D[i][j+1][k-1]>0.5) neighborCounter++;
                        if(segRawData->raw3D[i][j-1][k-1]>0.5) neighborCounter++;
                        if(segRawData->raw3D[i][j][k-1]>0.5) neighborCounter++;
                        if(segRawData->raw3D[i+1][j+1][k-1]>0.5) neighborCounter++;
                        if(segRawData->raw3D[i+1][j-1][k-1]>0.5) neighborCounter++;
                        if(segRawData->raw3D[i-1][j-1][k-1]>0.5) neighborCounter++;
                        if(segRawData->raw3D[i-1][j+1][k-1]>0.5) neighborCounter++;

                        ///if not enough neighbors, delete point from mask
                        if(neighborCounter<2) segRawData->raw3D[i][j][k] = 0;
                    }
                }
            }
        }
    }

    return segRawData;
}


rawData<double>* C_maskMaker::getMaskFromANAFile(const char* filenameData)
{
    C_ReadAnalyze75* toolRead = new C_ReadAnalyze75();
    dsr* dsrPtr = toolRead->readAnalyzeInfo(filenameData);
    rawData<double>* temp = toolRead->getRawData(filenameData);
    temp->pixDimX = dsrPtr->dime.pixdim[1];
    temp->pixDimY = dsrPtr->dime.pixdim[2];
    temp->pixDimZ = dsrPtr->dime.pixdim[3];
    temp->pixDimT = dsrPtr->dime.pixdim[4];
    delete toolRead;

    if(temp->numDim==3)
    {
        delete dsrPtr;
        return temp;
    }
    if(temp->numDim<3)
    {
        delete dsrPtr;
        delete temp;
        return NULL;
    }
    if(temp->numDim==4)
    {
        ///re-write the container as a 3D container with the first hyperplan of the 4D container
        rawData<double>* m = new rawData<double>(DICOM_DOUBLE, 3, temp->DimX, temp->DimY, temp->DimZ);
        for(unsigned long i=0 ; i<temp->DimX ; i++)
        {
            for(unsigned long j=0 ; j<temp->DimY ; j++)
            {
                for(unsigned long k=0 ; k<temp->DimZ ; k++)
                {
                    m->raw3D[i][j][k] = temp->raw4D[i][j][k][0];
                }
            }
        }
        delete temp;
        m->pixDimX = dsrPtr->dime.pixdim[1];
        m->pixDimY = dsrPtr->dime.pixdim[2];
        m->pixDimZ = dsrPtr->dime.pixdim[3];
        m->pixDimT = dsrPtr->dime.pixdim[4];
        delete dsrPtr;
        return m;
    }
    return NULL;
}

rawData<double>* C_maskMaker::createMaskFromTensors(rawData<double>* rawTensor, double FAmin, double FAmax)
{
    rawData<double>* segRawData = new rawData<double>(DICOM_DOUBLE,3, rawTensor->DimX, rawTensor->DimY, rawTensor->DimZ);
    double FA_;
    for(unsigned long i=0 ; i<rawTensor->DimX ; i++)
    {
        for(unsigned long j=0 ; j<rawTensor->DimY ; j++)
        {
            for(unsigned long k=0 ; k<rawTensor->DimZ ; k++)
            {
                C_toolbox_eigen_sym* toolEIG = new C_toolbox_eigen_sym(rawTensor->raw4D[i][j][k], false);
                FA_ = 0.5*(SQR(toolEIG->d[1]-toolEIG->d[0]) + SQR(toolEIG->d[2]-toolEIG->d[1]) + SQR(toolEIG->d[0]-toolEIG->d[2]))/(SQR(toolEIG->d[0])+SQR(toolEIG->d[1])+SQR(toolEIG->d[2]));
                if(FA_>=SQR(FAmin) && FA_<=SQR(FAmax))
                {
                    segRawData->raw3D[i][j][k] = 1.0;
                }
                else
                {
                    segRawData->raw3D[i][j][k] = 0.0;
                }
                delete toolEIG;
            }
        }
    }
    return segRawData;
}
