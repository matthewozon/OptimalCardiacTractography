#include <C_getDicomFromDirectory.h>

C_getDicomFromDirectory::C_getDicomFromDirectory(string directoryName)
{
    //ctor
    m_directoryName = directoryName;
    C_sort_files* sortFiles = new C_sort_files();
    dataFile = sortFiles->getSiFileNameDCM(directoryName);
    SoFile = sortFiles->getS0FileNameDCM(directoryName);
    roiFileName = sortFiles->getROIFileNameDCM(directoryName);
    NgradientDirection = dataFile.size();
    Nrepetition = SoFile.size(); //only for S0: it will be averaged
    C_readDicomDiffData* m = new C_readDicomDiffData(SoFile.at(0).data());
    m_nb_slice_in_mosa = m->getNumberOfSlice();
    delete m;
    delete sortFiles;
}

C_getDicomFromDirectory::~C_getDicomFromDirectory()
{
    //dtor
}

vector<diffInfo*> C_getDicomFromDirectory::getDiffInfoS0(void)
{
    ///load S0 info and data
    vector<diffInfo*> myB0InfoVect(Nrepetition,NULL);
    double** tempDoubleTableB0;
    double* tempDoubleArrayB0;
    for(unsigned int j=0 ; j<SoFile.size() ; j++)
    {
        C_readDicomDiffData* m = new C_readDicomDiffData(SoFile.at(j).data());
        myB0InfoVect.at(j) = new diffInfo;
        tempDoubleArrayB0 = m->getGradientDirectionSliceBasis();
        if(tempDoubleArrayB0==NULL)
        {
            myB0InfoVect.at(j)->g[0] = 0.0;    myB0InfoVect.at(j)->g[1] = 0.0;    myB0InfoVect.at(j)->g[2] = 1.0;
        }
        else
        {
            myB0InfoVect.at(j)->g[0] = tempDoubleArrayB0[0];    myB0InfoVect.at(j)->g[1] = tempDoubleArrayB0[1];    myB0InfoVect.at(j)->g[2] = tempDoubleArrayB0[2];
        }
        tempDoubleTableB0 = m->getChangeBasisMatrix();
        if(tempDoubleTableB0==NULL)
        {
            myB0InfoVect.at(j)->X_slice[0] = 1.0;   myB0InfoVect.at(j)->X_slice[1] = 0.0;   myB0InfoVect.at(j)->X_slice[2] = 0.0;
            myB0InfoVect.at(j)->Y_slice[0] = 0.0;   myB0InfoVect.at(j)->Y_slice[1] = 1.0;   myB0InfoVect.at(j)->Y_slice[2] = 0.0;
            myB0InfoVect.at(j)->Z_slice[0] = 0.0;   myB0InfoVect.at(j)->Z_slice[1] = 0.0;   myB0InfoVect.at(j)->Z_slice[2] = 1.0;
        }
        else
        {
            myB0InfoVect.at(j)->X_slice[0] = tempDoubleTableB0[0][0];   myB0InfoVect.at(j)->X_slice[1] = tempDoubleTableB0[0][1];   myB0InfoVect.at(j)->X_slice[2] = tempDoubleTableB0[0][2];
            myB0InfoVect.at(j)->Y_slice[0] = tempDoubleTableB0[1][0];   myB0InfoVect.at(j)->Y_slice[1] = tempDoubleTableB0[1][1];   myB0InfoVect.at(j)->Y_slice[2] = tempDoubleTableB0[1][2];
            myB0InfoVect.at(j)->Z_slice[0] = tempDoubleTableB0[2][0];   myB0InfoVect.at(j)->Z_slice[1] = tempDoubleTableB0[2][1];   myB0InfoVect.at(j)->Z_slice[2] = tempDoubleTableB0[2][2];
        }

        myB0InfoVect.at(j)->b = m->getBValue();
        myB0InfoVect.at(j)->echoTime = m->getEchoTime();
        myB0InfoVect.at(j)->repTime = m->getRepetitionTime();
        myB0InfoVect.at(j)->voxelSize[0] = m->getPixelDimX();
        myB0InfoVect.at(j)->voxelSize[1] = m->getPixelDimY();
        myB0InfoVect.at(j)->voxelSize[2] = m->getPixelDimZ();
        if(tempDoubleArrayB0!=NULL)
        {
            delete tempDoubleArrayB0;
        }
        if(tempDoubleTableB0!=NULL)
        {
            delete tempDoubleTableB0[0];
            delete tempDoubleTableB0[1];
            delete tempDoubleTableB0[2];
            delete tempDoubleTableB0;
        }
        delete m;
    }
    return myB0InfoVect;
}
vector<rawData<double>*> C_getDicomFromDirectory::getDataS0(void)
{
    ///load S0 info and data
    vector<rawData<double>*> allB0Data(Nrepetition,NULL);
    for(unsigned int j=0 ; j<SoFile.size() ; j++)
    {
        C_readDicomDiffData* m = new C_readDicomDiffData(SoFile.at(j).data());
        allB0Data.at(j) = m->getUnmosaicData(m_nb_slice_in_mosa/*m->getNumberOfSlice()*/);
        delete m;
    }
    return allB0Data;
}

rawData<double>* C_getDicomFromDirectory::getDataROI(void)
{
    rawData<double>* ROI_DATA = NULL;
    //cout << roiFileName << endl;
    C_readDicomDiffData* m = new C_readDicomDiffData(roiFileName);
    ROI_DATA = m->getUnmosaicData(m_nb_slice_in_mosa/*m->getNumberOfSlice()*/);
    //cout << "coucou" << endl;
    delete m;
    return ROI_DATA;
}

vector<vector<diffInfo*> > C_getDicomFromDirectory::getDiffInfoSi(void)
{
    ///load Si info
    vector< vector<diffInfo*> > myDiffInfoVect( NgradientDirection, vector<diffInfo*>(1/*Nrepetition*/,NULL) );
    double** tempDoubleTable;
    double* tempDoubleArray;
    for(unsigned int i=0 ; i<dataFile.size() ; i++)
    {
        for(unsigned int j=0 ; j<dataFile.at(i).size() ; j++)
        {
            C_readDicomDiffData* m = new C_readDicomDiffData(dataFile.at(i).at(j).data());
            myDiffInfoVect.at(i).at(j) = new diffInfo;
            tempDoubleArray = m->getGradientDirectionSliceBasis();
            myDiffInfoVect.at(i).at(j)->g[0] = tempDoubleArray[0];    myDiffInfoVect.at(i).at(j)->g[1] = tempDoubleArray[1];    myDiffInfoVect.at(i).at(j)->g[2] = tempDoubleArray[2];
            tempDoubleTable = m->getChangeBasisMatrix();
            myDiffInfoVect.at(i).at(j)->X_slice[0] = tempDoubleTable[0][0];   myDiffInfoVect.at(i).at(j)->X_slice[1] = tempDoubleTable[0][1];   myDiffInfoVect.at(i).at(j)->X_slice[2] = tempDoubleTable[0][2];
            myDiffInfoVect.at(i).at(j)->Y_slice[0] = tempDoubleTable[1][0];   myDiffInfoVect.at(i).at(j)->Y_slice[1] = tempDoubleTable[1][1];   myDiffInfoVect.at(i).at(j)->Y_slice[2] = tempDoubleTable[1][2];
            myDiffInfoVect.at(i).at(j)->Z_slice[0] = tempDoubleTable[2][0];   myDiffInfoVect.at(i).at(j)->Z_slice[1] = tempDoubleTable[2][1];   myDiffInfoVect.at(i).at(j)->Z_slice[2] = tempDoubleTable[2][2];
            myDiffInfoVect.at(i).at(j)->b = m->getBValue();
            myDiffInfoVect.at(i).at(j)->echoTime = m->getEchoTime();
            myDiffInfoVect.at(i).at(j)->repTime = m->getRepetitionTime();
            myDiffInfoVect.at(i).at(j)->voxelSize[0] = m->getPixelDimX();
            myDiffInfoVect.at(i).at(j)->voxelSize[1] = m->getPixelDimY();
            myDiffInfoVect.at(i).at(j)->voxelSize[2] = m->getPixelDimZ();
            delete tempDoubleArray;
            delete tempDoubleTable[0];
            delete tempDoubleTable[1];
            delete tempDoubleTable[2];
            delete tempDoubleTable;
            delete m;
        }
    }
    return myDiffInfoVect;
}
vector<vector<rawData<double>*> > C_getDicomFromDirectory::getDataSi(void)
{
    ///load Si data
    vector<vector<rawData<double>*> > allRawData( NgradientDirection, vector<rawData<double>*>(1/*Nrepetition*/,NULL) );
    for(unsigned int i=0 ; i<dataFile.size() ; i++)
    {
        for(unsigned int j=0 ; j<dataFile.at(i).size() ; j++)
        {
            C_readDicomDiffData* m = new C_readDicomDiffData(dataFile.at(i).at(j).data());
            allRawData.at(i).at(j) = m->getUnmosaicData(m_nb_slice_in_mosa/*m->getNumberOfSlice()*/);
            delete m;
        }
    }
    return allRawData;
}
