#include <C_thread_multiscale.h>

C_thread_multiscale::C_thread_multiscale(vector< vector<diffInfo*> > dfi_ /**info about dwi b!=0*/,\
                                          vector< vector<rawData<double>*> > dwiData_ /**dwi b!=0*/,\
                                           vector<diffInfo*> b0Info_ /**info about dwi b=0*/,\
                                            vector<rawData<double>*> b0Data_ /**dwi with b=0*/,\
                                             rawData<double>* dwiMask_ /**a mask for ROI*/):C_thread()
{
    //ctor
    dfi = vector< vector<diffInfo*> >( dfi_.size(), vector<diffInfo*>(dfi_.at(0).size(),NULL) );
    dwiData = vector< vector<rawData<double>*> >(dwiData_.size(), vector<rawData<double>*>(dwiData_.at(0).size(), NULL) );
    b0Info = vector<diffInfo*>(b0Info_.size(), NULL);
    b0Data = vector<rawData<double>*>(b0Data_.size(), NULL);

    dwiMask = new rawData<double>(dwiMask_->DICOM_TYPE, dwiMask_->numDim, dwiMask_->DimX, dwiMask_->DimY, dwiMask_->DimZ, dwiMask_->DimT);
    *dwiMask = *dwiMask_;

    for(unsigned long i=0 ; i<b0Data_.size() ; i++ )
    {
        b0Data.at(i) = new rawData<double>(b0Data_.at(i)->DICOM_TYPE, b0Data_.at(i)->numDim, b0Data_.at(i)->DimX, b0Data_.at(i)->DimY, b0Data_.at(i)->DimZ, b0Data_.at(i)->DimT);
        *(b0Data.at(i)) = *(b0Data_.at(i));
        b0Info.at(i) = new diffInfo;
        *(b0Info.at(i)) = *(b0Info_.at(i));

    }
    for(unsigned long i=0 ; i<dwiData_.size() ; i++)
    {
        for(unsigned long j=0 ; j<dwiData_.at(i).size() ; j++)
        {
            dwiData.at(i).at(j) = new rawData<double>(dwiData_.at(i).at(j)->DICOM_TYPE, dwiData_.at(i).at(j)->numDim, dwiData_.at(i).at(j)->DimX, dwiData_.at(i).at(j)->DimY, dwiData_.at(i).at(j)->DimZ, dwiData_.at(i).at(j)->DimT);
            *(dwiData.at(i).at(j)) = *(dwiData_.at(i).at(j));
            dfi.at(i).at(j) = new diffInfo;
            *(dfi.at(i).at(j)) = *(dfi_.at(i).at(j));
        }
    }
}

C_thread_multiscale::~C_thread_multiscale()
{
    //dtor
    for(unsigned int i=0 ; i<b0Data.size() ; i++)
    {
        delete b0Data.at(i);
        delete b0Info.at(i);
    }
    b0Data.clear();
    b0Info.clear();

    ///should delete dfi
    ///should delete dwiData
    for(unsigned int i=0 ; i<dwiData.size() ; i++)
    {
        for(unsigned int j=0 ; j<dwiData.at(i).size() ; j++)
        {
            delete dwiData.at(i).at(j);
            delete dfi.at(i).at(j);
        }
    }
    dwiData.clear();
    dfi.clear();
}


void C_thread_multiscale::execute()
{
    C_sampling<double>* toolSAMPLING = new C_sampling<double>();
    rawData<double>* maskTemp;
    rawData<double>* tensorData;
    //cout << "in execute: " << downSamplingType << " " << numPyramidLevel << endl;
    switch(downSamplingType)
    {
        case DWI_DOWN:
        {
            ///should check the number of diff dir and b value and...
            ///down sample DWI
                //b0
            rawData<double>* tempTensorData;
            for(unsigned long i=0 ; i<b0Data.size() ; i++ )
            {
                //down sampling for each b0
                for(unsigned long ii=0 ; ii<numPyramidLevel ; ii++)
                {
                    tempTensorData = toolSAMPLING->gaussPyramidStep3D(b0Data.at(i));
                    delete b0Data.at(i);
                    b0Data.at(i) = tempTensorData;
                }
            }
                //b!=0
            for(unsigned long i=0 ; i<dwiData.size() ; i++)
            {
                for(unsigned long j=0 ; j<dwiData.at(i).size() ; j++)
                {
                    //for each dwi signal
                    for(unsigned long ii=0 ; ii<numPyramidLevel ; ii++)
                    {
                        tempTensorData = toolSAMPLING->gaussPyramidStep3D(dwiData.at(i).at(j));
                        delete dwiData.at(i).at(j);
                        dwiData.at(i).at(j) = tempTensorData;
                    }
                }
            }
                //mask
            rawData<double>* tempDwiMask;
            //cout << dwiMask->numDim << " " << dwiMask->DimX << " " << dwiMask->DimY << " " << dwiMask->DimZ << " " << dwiMask->DimT << endl;
            for(unsigned long i=0 ; i<numPyramidLevel ; i++)
            {
                tempDwiMask = toolSAMPLING->gaussPyramidStep3D(dwiMask);
                delete dwiMask;
                dwiMask = tempDwiMask;
                (*dwiMask) = (*dwiMask)>0.5;
            }
            maskTemp = dwiMask;
            ///compute tensor //the method (makeTensors) computing tensor is implemented in library
            C_tensorMaker* toolTensor = new C_tensorMaker();
            tensorData = toolTensor->makeTensors(dfi /**info about dwi b!=0*/,\
                                          dwiData /**dwi b!=0*/,\
                                           b0Info /**info about dwi b=0*/,\
                                            b0Data /**dwi with b=0*/,\
                                             dwiMask /**a mask for ROI*/);
            delete toolTensor;
            break;
        }
        case DTI_DOWN:
        {
            ///check if 4 dimension
            ///compute tensor //the method (makeTensors) computing tensor is implemented in library
            C_tensorMaker* toolTensor = new C_tensorMaker();
            tensorData = toolTensor->makeTensors(dfi /**info about dwi b!=0*/,\
                                          dwiData /**dwi b!=0*/,\
                                           b0Info /**info about dwi b=0*/,\
                                            b0Data /**dwi with b=0*/,\
                                             dwiMask /**a mask for ROI*/);
            delete toolTensor;
            ///down sample DTI in usual space
            rawData<double>* tempTensorData;
            rawData<double>* tempDwiMask;
            for(unsigned long i=0 ; i<numPyramidLevel ; i++)
            {
                tempTensorData = toolSAMPLING->gaussPyramidStep3D(tensorData);
                delete tensorData;
                tensorData = tempTensorData;
                tempDwiMask = toolSAMPLING->gaussPyramidStep3D(dwiMask);
                delete dwiMask;
                dwiMask = tempDwiMask;
                (*dwiMask) = (*dwiMask)>0.5;
            }
            maskTemp = dwiMask;
            break;
        }
        case DTI_DOWN_LOG:
        {
            try
            {
                ///check if 4 dimension
                ///compute tensor //the method (makeTensors) computing tensor is implemented in library
                C_tensorMaker* toolTensor = new C_tensorMaker();
                tensorData = toolTensor->makeTensors(dfi /**info about dwi b!=0*/,\
                                              dwiData /**dwi b!=0*/,\
                                               b0Info /**info about dwi b=0*/,\
                                                b0Data /**dwi with b=0*/,\
                                                 dwiMask /**a mask for ROI*/);
                delete toolTensor;
                ///down sample DTI in logM space
                    //go to log Euclidean metrics
                 logM(tensorData, dwiMask);
                    //compute actual down sampling
                rawData<double>* tempTensorData;
                rawData<double>* tempDwiMask;
                for(unsigned long i=0 ; i<numPyramidLevel ; i++)
                {
                    tempTensorData = toolSAMPLING->gaussPyramidStep3D(tensorData);
                    delete tensorData;
                    tensorData = tempTensorData;
                    tempDwiMask = toolSAMPLING->gaussPyramidStep3D(dwiMask);
                    delete dwiMask;
                    dwiMask = tempDwiMask;
                    (*dwiMask) = (*dwiMask)>0.5;
                }
                    //go back to straight metrics
                expM(tensorData, dwiMask);
                maskTemp = dwiMask;
            }
            catch(char const* msg)
            {
                cout << msg << endl;
            }
            break;
        }
        default:
        {
            ///compute tensor //the method (makeTensors) computing tensor is implemented in library
            C_tensorMaker* toolTensor = new C_tensorMaker();
            tensorData = toolTensor->makeTensors(dfi /**info about dwi b!=0*/,\
                                          dwiData /**dwi b!=0*/,\
                                           b0Info /**info about dwi b=0*/,\
                                            b0Data /**dwi with b=0*/,\
                                             dwiMask /**a mask for ROI*/);
            delete toolTensor;
            ///dwon sample only graph: check dim also
            ///fill sampling factors so that you can have all graph resolution
            samplingFactor* m_SFGraph = new samplingFactor();
            m_SFGraph->Xfactor = 1.0/(pow((double) 2.0, (int) numPyramidLevel));
            m_SFGraph->Yfactor = 1.0/(pow((double) 2.0, (int) numPyramidLevel));
            m_SFGraph->Zfactor = 1.0/(pow((double) 2.0, (int) numPyramidLevel));
            m_SFGraph->Tfactor = 1.0;
            ///generate mask for the given resolution
            maskTemp = toolSAMPLING->samplingNearestNeighbor(dwiMask, m_SFGraph);
            (*maskTemp) = (*maskTemp)>0.5;
            delete m_SFGraph;
            break;
        }
    }
    delete toolSAMPLING;



    ///save data
    try
    {
        C_ReadAnalyze75* toolSave = new C_ReadAnalyze75();
        dsr* dsrPtr = new dsr;
        dsrPtr->hk.extents = 16384;
        //dsrPtr->hk.db_name = (char*) "private";
        dsrPtr->hk.regular = (char) 'r';
        dsrPtr->dime.datatype = DT_DOUBLE;
        dsrPtr->dime.bitpix = 8*sizeof(double);
        dsrPtr->dime.dim[0] = 4;
        dsrPtr->dime.dim[1] = tensorData->DimX;
        dsrPtr->dime.dim[2] = tensorData->DimY;
        dsrPtr->dime.dim[3] = tensorData->DimZ;
        dsrPtr->dime.dim[4] = tensorData->DimT;
        dsrPtr->dime.pixdim[0] = 4;
        dsrPtr->dime.pixdim[1] = tensorData->pixDimX;
        dsrPtr->dime.pixdim[2] = tensorData->pixDimY;
        dsrPtr->dime.pixdim[3] = tensorData->pixDimZ;
        dsrPtr->dime.pixdim[4] = tensorData->pixDimT;
        ostringstream tensorName;
        tensorName << "tensors_xx_xy_xz_yy_yz_zz";
        switch(downSamplingType)
        {
            case DWI_DOWN:
            {
                tensorName << "_DWI";
                break;
            }
            case DTI_DOWN:
            {
                tensorName << "_DTI";
                break;
            }
            case DTI_DOWN_LOG:
            {
                tensorName << "_DTILOG";
                break;
            }
            default:
            {
                tensorName << "_GRAPH";
                break;
            }
        }
        tensorName << "_pyramidalLevel_" << numPyramidLevel;
        toolSave->saveAna(tensorName.str().data(), dsrPtr, tensorData);


        dsrPtr->dime.dim[0] = 3;
        dsrPtr->dime.dim[1] = dwiMask->DimX;
        dsrPtr->dime.dim[2] = dwiMask->DimY;
        dsrPtr->dime.dim[3] = dwiMask->DimZ;
        dsrPtr->dime.dim[4] = 0;
        dsrPtr->dime.pixdim[0] = 3;
        dsrPtr->dime.pixdim[1] = dwiMask->pixDimX;
        dsrPtr->dime.pixdim[2] = dwiMask->pixDimY;
        dsrPtr->dime.pixdim[3] = dwiMask->pixDimZ;
        dsrPtr->dime.pixdim[4] = 0;
        ostringstream maskName;
        maskName << "tensors_mask";
        switch(downSamplingType)
        {
            case DWI_DOWN:
            {
                maskName << "_DWI";
                break;
            }
            case DTI_DOWN:
            {
                maskName << "_DTI";
                break;
            }
            case DTI_DOWN_LOG:
            {
                maskName << "_DTILOG";
                break;
            }
            default:
            {
                maskName << "_GRAPH";
                break;
            }
        }
        maskName << "_pyramidalLevel_" << numPyramidLevel;
        toolSave->saveAna(maskName.str().data(), dsrPtr, dwiMask);
        delete dsrPtr;
        delete toolSave;
    }
    catch(...)
    {
        cout << "tensors and/or mask could not be saved" << endl;
    }



    cout << "after switch in execute: " << downSamplingType << " " << numPyramidLevel << endl;
    cout << maskTemp->pixDimX << " " << maskTemp->pixDimY << " " << maskTemp->pixDimZ << " " << maskTemp->pixDimT << endl;
    cout << tensorData->pixDimX << " " << tensorData->pixDimY << " " << tensorData->pixDimZ << " " << tensorData->pixDimT << endl;
    cout << dwiMask->pixDimX << " " << dwiMask->pixDimY << " " << dwiMask->pixDimZ << " " << dwiMask->pixDimT << endl;

    ///create graph at given resolution
    time_t startTime = time(0);
    C_graph* G = new C_graph(maskTemp);
    time_t endTime = time(0);
    cout << "time for graph creation: " << difftime(endTime, startTime) * 1000.0 << " at pyramid level: " << numPyramidLevel  << endl;

    C_energy* U = new C_energy(G, tensorData, dwiMask);

    ///energy set up
    U->METHOD_INTERPOLATION = GAUSSIAN_PHI;
    U->OBJECT_INTERPOLATED = TENSOR;
    U->METHOD_INTEGRATION = SIMPSON_RULE;
    U->m_nbInterpolatingPoint = 1;
    U->m_alpha = 0.2;
    U->m_beta = 0.5;
    U->m_sigma = 0.375;
    U->m_phi0 = PHI0;
    U->m_phi1 = PHI1;
    U->m_t0 = Pi/2.0;
    U->m_R = 1.5*tensorData->pixDimX; //must be fix: regarding the resolution 1

    ///parameter
    unsigned long energyKind = U15;

    ///init energy
    unsigned long nb_thread = 8;
    //cout << "init energy in thread for resolution: " << m_resolution_factor << endl;
    startTime = time(0);
    U->initEnergy(energyKind, nb_thread);
    endTime = time(0);
    cout << "time for energy init: " << difftime(endTime, startTime) * 1000.0 << " at pyramid level: " << numPyramidLevel  << endl;
    //cout << "finish energy in thread for resolution: " << m_resolution_factor << endl;
    int endMem = getValue();
    cout << "vMem = " << endMem << " at pyramid level: " << numPyramidLevel << endl;
    //cout << "diff vMem = " << endMem -curMem << endl;
    //return;

    U->initPointData(energyKind);

    ///go for minimization
    unsigned long coolingKind = SA_EXP_BLOCK_COOLING;
    unsigned long kindOfMove = SA_MOVE_1;//SA_MOVE_1

    ///temperature assessment
    C_temperature* T = new C_temperature(G, U, energyKind);
    T->m_alphaMixture = 0.975;
    T->m_moveFromCenter = false;
    T->m_varSigma = 0.2*tensorData->pixDimX;
    T->m_voxelConstraint = true;
    double startingRate;
    double endRate;
    if(energyKind==U1)
    {
        startingRate = 0.3;
        endRate = 0.0030;
    }
    else
    {
        startingRate = 0.1;//2;//7;
        endRate = 0.0003;//0.001
    }
    T->runTemperature(kindOfMove, kindOfMove, startingRate, endRate, 1000);
    double Tinit = T->getTinit();
    double Tend = T->getTend();
    delete T;


    ///run minimization
    unsigned long nbIter = 64000;
    C_SA* minimization_by_SA = new C_SA(G, U, Tinit, Tend, nbIter, coolingKind, energyKind, kindOfMove);
    minimization_by_SA->alphaMixture = 0.975;
    minimization_by_SA->m_varSigma = 0.2*tensorData->pixDimX;
    minimization_by_SA->m_voxelConstraint = true;
    minimization_by_SA->m_moveFromCenter = false;
    ostringstream expName;
    expName << "SA_MOVE_1_SYNTHETIC_HEART_U15_Tinit_" << Tinit << "_Tend_" << Tend;
    expName << "_nbiter_" << nbIter << "_alpha_" << U->m_alpha << "_beta_" << U->m_beta << "_sigma_" << U->m_sigma;
    expName << "_Niterpol_" << U->m_nbInterpolatingPoint << "_radius_" << U->m_R << "_pyramidLevel_" << numPyramidLevel;
    switch(downSamplingType)
    {
        case DWI_DOWN:
        {
            expName << "_DWI";
            break;
        }
        case DTI_DOWN:
        {
            expName << "_DTI";
            break;
        }
        case DTI_DOWN_LOG:
        {
            expName << "_DTILOG";
            break;
        }
        default:
        {
            expName << "_GRAPH";
            break;
        }
    }
    minimization_by_SA->baseName = expName.str();
    minimization_by_SA->m_expTAG = "_";

    //init graph
    minimization_by_SA->randInit(0.2);

    //run
    minimization_by_SA->flipRate = true;
    minimization_by_SA->energyMeasure = true;
    startTime = time(0);
    minimization_by_SA->run();
    endTime = time(0);
    cout << "time per iteration: " << difftime(endTime, startTime) * 1000.0/((double) (nbIter+MAX_EDGES_PER_NODE)) << endl;

    G->saveEdges((char*) expName.str().data());

    ///save fibers
    C_bundle* s = new C_bundle(G, true); //extract fibers
    minimization_by_SA->saveFiber((char*) expName.str().data(), s->FIBERS); //save fiber in file

    ///free memory
    delete s;
    delete minimization_by_SA;

    ///delete mask
    if(downSamplingType==GRAPH_DOWN) delete maskTemp;

    delete tensorData;

    delete U;
    delete G;
    return;
}


int C_thread_multiscale::parseLine(char* line)
{
    int i = strlen(line);
    while (*line < '0' || *line > '9') line++;
    line[i-3] = '\0';
    i = atoi(line);
    return i;
}

int C_thread_multiscale::getValue() //Note: this value is in KB!
{
    FILE* file = fopen("/proc/self/status", "r");
    int result = -1;
    char line[128];


    while (fgets(line, 128, file) != NULL){
        if (strncmp(line, "VmSize:", 7) == 0){
            result = parseLine(line);
            break;
        }
    }
    fclose(file);
    return result;
}

void C_thread_multiscale::logM(rawData<double>* rawDTI, rawData<double>* rawMask)
{
    ///check dimensions
    if(rawDTI->numDim!=4) return;
    if(rawDTI->DimT!=6) return;
    if(rawMask->numDim!=3) return;
    //tool for matrix product
    C_toolbox_SVD* toolProdMat = new C_toolbox_SVD();
    ///for each tensor
    double** D;
    double** PtimesD;
    double** Pt;
    double** M;
    C_toolbox_eigen_sym* toolEIG;
    for(unsigned long i=0 ; i<rawDTI->DimX ; i++)
    {
        for(unsigned long j=0 ; j<rawDTI->DimY ; j++)
        {
            for(unsigned long k=0 ; k<rawDTI->DimZ ; k++)
            {
                if(rawMask->raw3D[i][j][k]>0.00001)
                {
                    ///tensor diagnalization
                    toolEIG = new C_toolbox_eigen_sym(rawDTI->raw4D[i][j][k]);
                    ///compute log of each eigen value
                    toolEIG->d[0] = log(toolEIG->d[0]); //this is actually ln
                    toolEIG->d[1] = log(toolEIG->d[1]);
                    toolEIG->d[2] = log(toolEIG->d[2]);
                    ///form back the tensor and store it into argument;
                    D = toolProdMat->diag(toolEIG->d,3);
                    PtimesD = toolProdMat->matProd(toolEIG->z, 3, 3, D, 3, 3);
                    Pt = toolProdMat->transpose(toolEIG->z, 3, 3);
                    M = toolProdMat->matProd(PtimesD, 3, 3, Pt, 3, 3);
                    rawDTI->raw4D[i][j][k][0] = M[0][0];
                    rawDTI->raw4D[i][j][k][1] = M[0][1];
                    rawDTI->raw4D[i][j][k][2] = M[0][2];
                    rawDTI->raw4D[i][j][k][3] = M[1][1];
                    rawDTI->raw4D[i][j][k][4] = M[1][2];
                    rawDTI->raw4D[i][j][k][5] = M[2][2];
                    ///delete what was allocated in loop
                    for(short u=0 ; u<3 ; u++)
                    {
                        delete M[u];
                        delete Pt[u];
                        delete PtimesD[u];
                        delete D[u];
                    }
                    delete M;
                    delete Pt;
                    delete PtimesD;
                    delete D;
                    delete toolEIG;
                }
            }
        }
    }
    delete toolProdMat;
}

void C_thread_multiscale::expM(rawData<double>* rawDTI, rawData<double>* rawMask)
{
    ///check dimensions
    if(rawDTI->numDim!=4) return;
    if(rawDTI->DimT!=6) return;
    if(rawMask->numDim!=3) return;
    //tool for matrix product
    C_toolbox_SVD* toolProdMat = new C_toolbox_SVD();
    ///for each tensor
    double** D;
    double** PtimesD;
    double** Pt;
    double** M;
    C_toolbox_eigen_sym* toolEIG;
    for(unsigned long i=0 ; i<rawDTI->DimX ; i++)
    {
        for(unsigned long j=0 ; j<rawDTI->DimY ; j++)
        {
            for(unsigned long k=0 ; k<rawDTI->DimZ ; k++)
            {
                if(rawMask->raw3D[i][j][k]>0.00001)
                {
                    ///tensor diagnalization
                    toolEIG = new C_toolbox_eigen_sym(rawDTI->raw4D[i][j][k]);
                    ///compute log of each eigen value
                    toolEIG->d[0] = exp(toolEIG->d[0]);
                    toolEIG->d[1] = exp(toolEIG->d[1]);
                    toolEIG->d[2] = exp(toolEIG->d[2]);
                    ///form back the tensor and store it into argument;
                    D = toolProdMat->diag(toolEIG->d,3);
                    PtimesD = toolProdMat->matProd(toolEIG->z, 3, 3, D, 3, 3);
                    Pt = toolProdMat->transpose(toolEIG->z, 3, 3);
                    M = toolProdMat->matProd(PtimesD, 3, 3, Pt, 3, 3);
                    rawDTI->raw4D[i][j][k][0] = M[0][0];
                    rawDTI->raw4D[i][j][k][1] = M[0][1];
                    rawDTI->raw4D[i][j][k][2] = M[0][2];
                    rawDTI->raw4D[i][j][k][3] = M[1][1];
                    rawDTI->raw4D[i][j][k][4] = M[1][2];
                    rawDTI->raw4D[i][j][k][5] = M[2][2];
                    ///delete what was allocated in loop
                    for(short u=0 ; u<3 ; u++)
                    {
                        delete M[u];
                        delete Pt[u];
                        delete PtimesD[u];
                        delete D[u];
                    }
                    delete M;
                    delete Pt;
                    delete PtimesD;
                    delete D;
                    delete toolEIG;
                }

            }
        }
    }
    delete toolProdMat;
}
