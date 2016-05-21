#include <C_tensorMaker.h>

C_tensorMaker::C_tensorMaker()
{
    //ctor
    toolSVD = new C_toolbox_SVD();
}

C_tensorMaker::~C_tensorMaker()
{
    //dtor
    delete toolSVD;
}


rawData<double>* C_tensorMaker::makeTensors(vector<vector<diffInfo*> > dwiInfo, vector<vector<rawData<double>*> > dwiData, vector<diffInfo*>  b0Info, vector<rawData<double>*> b0Data, rawData<double>* dwiMask)
{
    //get number of direction
    unsigned long NgradientDirection = dwiData.size();
    //get number of repetition
    unsigned long Nrepetition = dwiData.at(0).size();
    //get number of dimension in rawData
    unsigned long Ndim = dwiData.at(0).at(0)->numDim;
    //get b value of DWI
    double b = dwiInfo.at(0).at(0)->b;
    //get b value of b0 (may be different from 0)
    double b0 = b0Info.at(0)->b;
    double gx0 = b0Info.at(0)->g[0];
    double gy0 = b0Info.at(0)->g[1];
    double gz0 = b0Info.at(0)->g[2];

    //allocate a rawData of dimension + 1: rawDataTensor (last dim has 6 components)
    rawData<double>* rawDataTensor = NULL;
    if(Ndim==1)
    {
        rawDataTensor = new rawData<double>(DICOM_DOUBLE, Ndim+1, dwiData.at(0).at(0)->DimX, 6);
    }
    else if(Ndim==2)
    {
        rawDataTensor = new rawData<double>(DICOM_DOUBLE, Ndim+1, dwiData.at(0).at(0)->DimX, dwiData.at(0).at(0)->DimY, 6);
    }
    else if(Ndim==3)
    {
        rawDataTensor = new rawData<double>(DICOM_DOUBLE, Ndim+1, dwiData.at(0).at(0)->DimX, dwiData.at(0).at(0)->DimY, dwiData.at(0).at(0)->DimZ, 6);
    }
    else
    {
        cout << "Will not be able to produce tensors" << endl;
        return rawDataTensor;
    }
    rawDataTensor->pixDimX = dwiData.at(0).at(0)->pixDimX;
    rawDataTensor->pixDimY = dwiData.at(0).at(0)->pixDimY;
    rawDataTensor->pixDimZ = dwiData.at(0).at(0)->pixDimZ;
    rawDataTensor->pixDimT = dwiData.at(0).at(0)->pixDimT;

    //if several repetition, average all and store the results in myAverageDwiDataVector, otherwise copy pointer in myAverageDwiDataVector
    //the same for b0 data
    vector<rawData<double>*> myAverageDwiDataVector = averageRawData(dwiData);
    rawData<double>* myAverageB0DatVector = new rawData<double>(DICOM_DOUBLE, b0Data.at(0)->numDim, b0Data.at(0)->DimX, b0Data.at(0)->DimY, b0Data.at(0)->DimZ, b0Data.at(0)->DimT);
    (*myAverageB0DatVector) = (*b0Data.at(0));
    if(Nrepetition!=1)
    {
        for(unsigned long j=1 ; j<Nrepetition ; j++)
        {
            (*myAverageB0DatVector) += (*b0Data.at(j));
        }
        (*myAverageB0DatVector) /= ((double) Nrepetition);
    }

    //create H matrix
    double** H = new double*[NgradientDirection]; //should check allocation result
    double gx, gy, gz;
    //cout << "computing H matrix" << endl;
    for(unsigned long i=0 ; i<NgradientDirection ; i++)
    {
        H[i] = new double[6];
        gx = dwiInfo.at(i).at(0)->g[0];
        gy = dwiInfo.at(i).at(0)->g[1];
        gz = dwiInfo.at(i).at(0)->g[2];

        //create a row of the raw matrix (before modification if needed)
        H[i][0] = SQR(gx); H[i][1] = SQR(gy); H[i][2] = SQR(gz); H[i][3] = 2*gx*gy; H[i][4] = 2*gx*gz; H[i][5] = 2*gy*gz;

        if(b0>0.0001)
        {
            //compute correction
            H[i][0] -= (b0/b)*SQR(gx0); H[i][1] -= (b0/b)*SQR(gy0); H[i][2] -= (b0/b)*SQR(gz0); H[i][3] -= 2*(b0/b)*gx0*gy0; H[i][4] -= 2*(b0/b)*gx0*gz0; H[i][5] -= 2*(b0/b)*gy0*gz0;
        }
    }



    //double** DTI = T->matProd(HpseudoInv, 6, NgradientDirection, Y, NgradientDirection, 1);
    //for each voxel in mask, fill Y, calculate Hpsi*Y, store the result in rawDataTensor

    *rawDataTensor = 0.0;

    computeTensors(myAverageDwiDataVector, myAverageB0DatVector, dwiMask, b, H, rawDataTensor);

    //destroy
    for(unsigned int i=0 ; i<myAverageDwiDataVector.size() ; i++)
    {
        delete myAverageDwiDataVector.at(i);//->~rawData();
    }
    myAverageDwiDataVector.clear();
    delete myAverageB0DatVector;//->~rawData();
    for(unsigned long i=0 ; i<NgradientDirection ; i++)
    {
        delete H[i];
    }
    delete H;

    return rawDataTensor;
}

bool C_tensorMaker::computeTensors(vector<rawData<double>*> myAverageDwiDataVector, rawData<double>* myAverageB0DatVector, rawData<double>* dwiMask, double b, double**H, rawData<double>* tensorOutput)
{
    //create Hpsi
    unsigned long NgradientDirection = myAverageDwiDataVector.size();
    double** Hpsi = toolSVD->pseudoInv(H,NgradientDirection, 6);
//    cout << "H" << endl;
//    for(unsigned long i=0 ; i<NgradientDirection ; i++)
//    {
//        cout << H[i][0] << ", " << H[i][1] << ", " << H[i][2] << ", " << H[i][3] << ", " << H[i][4] << ", " << H[i][5] << endl;
//    }
//    cout << "Hpsi" << endl;
//    for(unsigned long i=0 ; i<NgradientDirection ; i++)
//    {
//        cout << Hpsi[0][i] << ", " << Hpsi[1][i] << ", " << Hpsi[2][i] << ", " << Hpsi[3][i] << ", " << Hpsi[4][i] << ", " << Hpsi[5][i] << endl;
//    }
//    getchar();

    //allocate Y
    double** Y = new double*[NgradientDirection];
    for(unsigned long i=0 ; i<NgradientDirection ; i++)
    {
        Y[i] = new double[1];
    }

    //for each element
    if(myAverageB0DatVector->numDim==1)
    {
        for(unsigned long i=0 ; i<myAverageB0DatVector->DimX ; i++ )
        {
            if(dwiMask->raw1D[i]>0.5) //0=not in mask, 1=in mask
            {
                //fill Y
                for(unsigned long l=0 ; l<NgradientDirection ; l++)
                {
                    Y[l][0] = (1.0/b)*log(myAverageB0DatVector->raw1D[i]/myAverageDwiDataVector.at(l)->raw1D[i]);
                    if(Y[l][0]!=Y[l][0] || isinf(Y[l][0])) //in case the operation produces a NaN, try something else
                    {
                        Y[l][0] = (-1.0/b)*log(myAverageDwiDataVector.at(l)->raw1D[i]/myAverageB0DatVector->raw1D[i]);
                        if(Y[l][0]!=Y[l][0] || isinf(Y[l][0]))//remove the contribution of this element
                        {
                            cout << "may produce biased tensor" << endl;
                            Y[l][0] = 0.0;
                        }
                    }
                }
                //compute tensor
                double** DTI = toolSVD->matProd(Hpsi, 6, NgradientDirection, Y, NgradientDirection, 1); //produce a matrix 6 rows by 1 column [Dxx Dyy Dzz Dxy Dxz Dyz]
                //store it as [Dxx Dxy Dxz Dyy Dyz Dzz]
                tensorOutput->raw2D[i][0] = DTI[0][0];
                tensorOutput->raw2D[i][1] = DTI[3][0];
                tensorOutput->raw2D[i][2] = DTI[4][0];
                tensorOutput->raw2D[i][3] = DTI[1][0];
                tensorOutput->raw2D[i][4] = DTI[5][0];
                tensorOutput->raw2D[i][5] = DTI[2][0];

                ///normalization by main eigen value: usefull only for U1, but there is no problem for the other energy function
                C_toolbox_eigen_sym* toolEIGEN = new C_toolbox_eigen_sym(DTI);
                if(toolEIGEN->d[0]>0.0)
                {
                    tensorOutput->raw2D[i][0] /= toolEIGEN->d[0];
                    tensorOutput->raw2D[i][1] /= toolEIGEN->d[0];
                    tensorOutput->raw2D[i][2] /= toolEIGEN->d[0];
                    tensorOutput->raw2D[i][3] /= toolEIGEN->d[0];
                    tensorOutput->raw2D[i][4] /= toolEIGEN->d[0];
                    tensorOutput->raw2D[i][5] /= toolEIGEN->d[0];
                }
                delete toolEIGEN;
                for(unsigned long m=0 ; m<6 ; m++)
                {
                    delete DTI[m];
                }
                delete DTI;
            }

        }
        return true;
    }
    else if(myAverageB0DatVector->numDim==2)
    {
        for(unsigned long i=0 ; i<myAverageB0DatVector->DimX ; i++ )
        {
            for(unsigned long j=0 ; j<myAverageB0DatVector->DimY ; j++ )
            {
                if(dwiMask->raw2D[i][j]>0.5) //0=not in mask, 1=in mask
                {
                    //fill Y
                    for(unsigned long l=0 ; l<NgradientDirection ; l++)
                    {
                        Y[l][0] = (1.0/b)*log(myAverageB0DatVector->raw2D[i][j]/myAverageDwiDataVector.at(l)->raw2D[i][j]);
                        if(Y[l][0]!=Y[l][0] || isinf(Y[l][0])) //in case the operation produces a NaN, try something else
                        {
                            Y[l][0] = (-1.0/b)*log(myAverageDwiDataVector.at(l)->raw2D[i][j]/myAverageB0DatVector->raw2D[i][j]);
                            if(Y[l][0]!=Y[l][0] || isinf(Y[l][0]))//remove the contribution of this element
                            {
                                cout << "may produce biased tensor" << endl;
                                Y[l][0] = 0.0;
                            }
                        }
                    }
                    //compute tensor
                    double** DTI = toolSVD->matProd(Hpsi, 6, NgradientDirection, Y, NgradientDirection, 1); //produce a matrix 6 rows by 1 column [Dxx Dyy Dzz Dxy Dxz Dyz]
                    //store it as [Dxx Dxy Dxz Dyy Dyz Dzz]
                    tensorOutput->raw3D[i][j][0] = DTI[0][0];
                    tensorOutput->raw3D[i][j][1] = DTI[3][0];
                    tensorOutput->raw3D[i][j][2] = DTI[4][0];
                    tensorOutput->raw3D[i][j][3] = DTI[1][0];
                    tensorOutput->raw3D[i][j][4] = DTI[5][0];
                    tensorOutput->raw3D[i][j][5] = DTI[2][0];

                    ///normalization by main eigen value: usefull only for U1, but there is no problem for the other energy function
                    C_toolbox_eigen_sym* toolEIGEN = new C_toolbox_eigen_sym(DTI);//should try to handle errors
                    if(toolEIGEN->d[0]>0.0)
                    {
                        tensorOutput->raw3D[i][j][0] /= toolEIGEN->d[0];
                        tensorOutput->raw3D[i][j][1] /= toolEIGEN->d[0];
                        tensorOutput->raw3D[i][j][2] /= toolEIGEN->d[0];
                        tensorOutput->raw3D[i][j][3] /= toolEIGEN->d[0];
                        tensorOutput->raw3D[i][j][4] /= toolEIGEN->d[0];
                        tensorOutput->raw3D[i][j][5] /= toolEIGEN->d[0];
                    }
                    delete toolEIGEN;
                    for(unsigned long m=0 ; m<6 ; m++)
                    {
                        delete DTI[m];
                    }
                    delete DTI;
                }
            }
        }
        return true;
    }
    else
    {
        for(unsigned long i=0 ; i<myAverageB0DatVector->DimX ; i++ )
        {
            for(unsigned long j=0 ; j<myAverageB0DatVector->DimY ; j++ )
            {
                for(unsigned long k=0 ; k<myAverageB0DatVector->DimZ ; k++ )
                {
                    if(dwiMask->raw3D[i][j][k]>0.5) //0=not in mask, 1=in mask
                    {
                        //fill Y
                        for(unsigned long l=0 ; l<NgradientDirection ; l++)
                        {
                            Y[l][0] = (1.0/b)*log(myAverageB0DatVector->raw3D[i][j][k]/myAverageDwiDataVector.at(l)->raw3D[i][j][k]);
                            if(Y[l][0]!=Y[l][0] || isinf(Y[l][0])) //in case the operation produces a NaN, try something else
                            {
                                Y[l][0] = (-1.0/b)*log(myAverageDwiDataVector.at(l)->raw3D[i][j][k]/myAverageB0DatVector->raw3D[i][j][k]);
                                if(Y[l][0]!=Y[l][0] || isinf(Y[l][0]))//remove the contribution of this element
                                {
                                    cout << "may produce biased tensor" << endl;
                                    Y[l][0] = 0.0;
                                }
                            }
//                            if(isinf(Y[l][0]))
//                            {
//                                if(Y[l][0]!=Y[l][0]) cout << "Y[l][0]!=Y[l][0] even for Inf values" << endl;
//                                cout << Y[l][0] << endl;
//                                cout << "b value: " << b << " B0 map: " << myAverageB0DatVector->raw3D[i][j][k] << " DWI: " << myAverageDwiDataVector.at(l)->raw3D[i][j][k] << endl;
//                            }
                        }
                        //compute tensor
                        double** DTI = toolSVD->matProd(Hpsi, 6, NgradientDirection, Y, NgradientDirection, 1); //produce a matrix 6 rows by 1 column [Dxx Dyy Dzz Dxy Dxz Dyz]
                        //store it as [Dxx Dxy Dxz Dyy Dyz Dzz]
                        tensorOutput->raw4D[i][j][k][0] = DTI[0][0];
                        tensorOutput->raw4D[i][j][k][1] = DTI[3][0];
                        tensorOutput->raw4D[i][j][k][2] = DTI[4][0];
                        tensorOutput->raw4D[i][j][k][3] = DTI[1][0];
                        tensorOutput->raw4D[i][j][k][4] = DTI[5][0];
                        tensorOutput->raw4D[i][j][k][5] = DTI[2][0];

                        ///normalization by main eigen value: usefull only for U1, but there is no problem for the other energy function
//                        cout << endl;
//                        cout << tensorOutput->raw4D[i][j][k][0] << " " << tensorOutput->raw4D[i][j][k][1] << " " << tensorOutput->raw4D[i][j][k][2] << endl;
//                        cout << tensorOutput->raw4D[i][j][k][1] << " " << tensorOutput->raw4D[i][j][k][3] << " " << tensorOutput->raw4D[i][j][k][4] << endl;
//                        cout << tensorOutput->raw4D[i][j][k][2] << " " << tensorOutput->raw4D[i][j][k][4] << " " << tensorOutput->raw4D[i][j][k][5] << endl;
//                        cout << endl;
                        C_toolbox_eigen_sym* toolEIGEN = new C_toolbox_eigen_sym(DTI);
                        if(toolEIGEN->d[0]>0.0)
                        {
                            tensorOutput->raw4D[i][j][k][0] /= toolEIGEN->d[0];
                            tensorOutput->raw4D[i][j][k][1] /= toolEIGEN->d[0];
                            tensorOutput->raw4D[i][j][k][2] /= toolEIGEN->d[0];
                            tensorOutput->raw4D[i][j][k][3] /= toolEIGEN->d[0];
                            tensorOutput->raw4D[i][j][k][4] /= toolEIGEN->d[0];
                            tensorOutput->raw4D[i][j][k][5] /= toolEIGEN->d[0];
                        }
                        delete toolEIGEN;
//                        double tr = DTI[0][0] + DTI[1][0] + DTI[2][0];
//                        cout << "tensor in mask" << endl;
//                        cout << "\t" << tensorOutput->raw4D[i][j][k][0]/tr << "\t" << tensorOutput->raw4D[i][j][k][1]/tr << "\t" << tensorOutput->raw4D[i][j][k][2]/tr << endl;
//                        cout << "\t" << tensorOutput->raw4D[i][j][k][1]/tr << "\t" << tensorOutput->raw4D[i][j][k][3]/tr << "\t" << tensorOutput->raw4D[i][j][k][4]/tr << endl;
//                        cout << "\t" << tensorOutput->raw4D[i][j][k][2]/tr << "\t" << tensorOutput->raw4D[i][j][k][4]/tr << "\t" << tensorOutput->raw4D[i][j][k][5]/tr << endl;
//                        getchar();
                        for(unsigned long m=0 ; m<6 ; m++)
                        {
                            delete DTI[m];
                        }
                        delete DTI;
                    }
                }
            }
        }
        return true;
    }


    //destroy garbage
    for(unsigned long i=0 ; i<NgradientDirection ; i++)
    {
        delete Y[i];
    }
    delete Y;
    for(unsigned long i=0 ; i<6 ; i++)
    {
        delete Hpsi[i];
    }
    delete Hpsi;
}

rawData<double>* C_tensorMaker::makeTensors(vector<diffInfo*> dwiInfo, rawData<double>* dwiDataS, rawData<double>* dwiMask)
{
    if(dwiDataS->numDim!=4) return NULL;

    //get number of direction
    unsigned long NgradientDirection = dwiDataS->DimT-1;

    //get number of dimension in rawData
    //unsigned long Ndim = dwiDataS->numDim;
    //get b value of DWI
    double b = dwiInfo.at(1)->b;
    //get b value of b0 (may be different from 0)
    double b0 = dwiInfo.at(0)->b;
    double gx0 = dwiInfo.at(0)->g[0];
    double gy0 = dwiInfo.at(0)->g[1];
    double gz0 = dwiInfo.at(0)->g[2];

    //allocate a rawData of dimension + 1: rawDataTensor (last dim has 6 components)
    rawData<double>* rawDataTensor = new rawData<double>(DICOM_DOUBLE, 4, dwiDataS->DimX, dwiDataS->DimY, dwiDataS->DimZ, 6);
    rawDataTensor->pixDimX = dwiDataS->pixDimX;
    rawDataTensor->pixDimY = dwiDataS->pixDimY;
    rawDataTensor->pixDimZ = dwiDataS->pixDimZ;
    rawDataTensor->pixDimT = dwiDataS->pixDimT;

    //create H matrix
    double** H = new double*[NgradientDirection]; //should check allocation result
    double gx, gy, gz;
    //cout << "computing H matrix" << endl;
    for(unsigned long i=0 ; i<NgradientDirection ; i++)
    {
        H[i] = new double[6];
        gx = dwiInfo.at(i+1)->g[0];
        gy = dwiInfo.at(i+1)->g[1];
        gz = dwiInfo.at(i+1)->g[2];

        //create a row of the raw matrix (before modification if needed)
        H[i][0] = SQR(gx); H[i][1] = SQR(gy); H[i][2] = SQR(gz); H[i][3] = 2*gx*gy; H[i][4] = 2*gx*gz; H[i][5] = 2*gy*gz;

        if(b0>0.0001)
        {
            //compute correction
            H[i][0] -= (b0/b)*SQR(gx0); H[i][1] -= (b0/b)*SQR(gy0); H[i][2] -= (b0/b)*SQR(gz0); H[i][3] -= 2*(b0/b)*gx0*gy0; H[i][4] -= 2*(b0/b)*gx0*gz0; H[i][5] -= 2*(b0/b)*gy0*gz0;
        }
        //cout << H[i][0] << " " << H[i][1] << " " << H[i][2] << " " << H[i][3] << " " << H[i][4] << " " << H[i][5] << endl;
    }
    //getchar();

    double** Hpsi = toolSVD->pseudoInv(H,NgradientDirection, 6);


    //double** DTI = T->matProd(HpseudoInv, 6, NgradientDirection, Y, NgradientDirection, 1);
    //for each voxel in mask, fill Y, calculate Hpsi*Y, store the result in rawDataTensor
    double** Y = new double*[NgradientDirection];
    for(unsigned long i=0 ; i<NgradientDirection ; i++)
    {
        Y[i] = new double[1];
    }
    for(unsigned long i=0 ; i<rawDataTensor->DimX ; i++)
    {
        for(unsigned long j=0 ; j<rawDataTensor->DimY ; j++)
        {
            for(unsigned long k=0 ; k<rawDataTensor->DimZ ; k++)
            {
                if(dwiMask->raw3D[i][j][k]>0.5)
                {
                    ///fill Y vector
                    for(unsigned long t=0 ; t<NgradientDirection ; t++)
                    {
//                        if(dwiDataS->raw4D[i][j][k][0]==0 || dwiDataS->raw4D[i][j][k][t+1]==0)
//                        {
//                            cout << "will generate a bad number" << endl;
//                        }
                        Y[t][0] = (1.0/b)*log(dwiDataS->raw4D[i][j][k][0]/dwiDataS->raw4D[i][j][k][t+1]);
                        if(Y[t][0]!=Y[t][0] || isinf(Y[t][0])) //in case the operation produces a NaN, try something else
                        {
                            Y[t][0] = (-1.0/b)*log(dwiDataS->raw4D[i][j][k][t+1]/dwiDataS->raw4D[i][j][k][0]);
                            if(Y[t][0]!=Y[t][0] || isinf(Y[t][0]))//remove the contribution of this element
                            {
                                cout << "may produce biased tensor" << endl;
                                Y[t][0] = 0.0;
                            }
                        }
                        //cout << Y[t][0] << " ";
                    }
//                    cout << endl;
//                    getchar();

                    ///compute DTI
                    double** DTI = toolSVD->matProd(Hpsi, 6, NgradientDirection, Y, NgradientDirection, 1);
                    ///store
                    rawDataTensor->raw4D[i][j][k][0] = DTI[0][0];
                    rawDataTensor->raw4D[i][j][k][1] = DTI[3][0];
                    rawDataTensor->raw4D[i][j][k][2] = DTI[4][0];
                    rawDataTensor->raw4D[i][j][k][3] = DTI[1][0];
                    rawDataTensor->raw4D[i][j][k][4] = DTI[5][0];
                    rawDataTensor->raw4D[i][j][k][5] = DTI[2][0];

                    ///delete
                    for(unsigned long m=0 ; m<6 ; m++)
                    {
                        delete DTI[m];
                    }
                    delete DTI;
                }
            }
        }
    }


    //destroy
    for(unsigned long i=0 ; i<NgradientDirection ; i++)
    {
        delete Y[i];
    }
    delete Y;
    for(unsigned long i=0 ; i<6 ; i++)
    {
        delete Hpsi[i];
    }
    delete Hpsi;
    for(unsigned long i=0 ; i<NgradientDirection ; i++)
    {
        delete H[i];
    }
    delete H;

    return rawDataTensor;
}
double** C_tensorMaker::makeHpsi(vector<vector<diffInfo*> > dwiInfo, vector<diffInfo*> b0Info)
{
    //get number of direction
    unsigned long N = dwiInfo.size();
    //get b value of DWI
    double b = dwiInfo.at(0).at(0)->b;
    //get b value of b0 (may be different from 0)
    double b0 = b0Info.at(0)->b;
    double gx0 = b0Info.at(0)->g[0];
    double gy0 = b0Info.at(0)->g[1];
    double gz0 = b0Info.at(0)->g[2];

    //create H matrix
    double** H = new double*[N]; //should check allocation result
    double gx, gy, gz;

    for(unsigned long i=0 ; i<N ; i++)
    {
        H[i] = new double[6];
        gx = dwiInfo.at(i).at(0)->g[0];
        gy = dwiInfo.at(i).at(0)->g[1];
        gz = dwiInfo.at(i).at(0)->g[2];

        //create a row of the raw matrix (before modification if needed)
        H[i][0] = SQR(gx); H[i][1] = SQR(gy); H[i][2] = SQR(gz); H[i][3] = 2*gx*gy; H[i][4] = 2*gx*gz; H[i][5] = 2*gy*gz;

        if(b0>0.0001)
        {
            //compute correction
            H[i][0] -= (b0/b)*SQR(gx0); H[i][1] -= (b0/b)*SQR(gy0); H[i][2] -= (b0/b)*SQR(gz0); H[i][3] -= 2*(b0/b)*gx0*gy0; H[i][4] -= 2*(b0/b)*gx0*gz0; H[i][5] -= 2*(b0/b)*gy0*gz0;
        }
    }

    return toolSVD->pseudoInv(H,N, 6);
}
double* C_tensorMaker::makeTensor(double* dwiArray/***N rows*/, double** Hpsi /**6 rows by N columns*/, double b /**b-value for dwiArray*/, double S0, unsigned long N)
{
    double** Y = new double*[N];
    for(unsigned long i=0 ; i<N ; i++)
    {
        Y[i] = new double[1];
    }

    for(unsigned long t=0 ; t<N ; t++)
    {
        Y[t][0] = (1.0/b)*log(S0/dwiArray[t]);
        if(Y[t][0]!=Y[t][0] || isinf(Y[t][0])) //in case the operation produces a NaN, try something else
        {
            Y[t][0] = (-1.0/b)*log(dwiArray[t]/S0);
            if(Y[t][0]!=Y[t][0] || isinf(Y[t][0]))//remove the contribution of this element
            {
                cout << "may produce biased tensor" << endl;
                Y[t][0] = 0.0;
            }
        }
    }

    ///compute DTI
    double** tempDTI = toolSVD->matProd(Hpsi, 6, N, Y, N, 1);
    ///store
    double* DTI = new double[6];
    DTI[0] = tempDTI[0][0];
    DTI[1] = tempDTI[3][0];
    DTI[2] = tempDTI[4][0];
    DTI[3] = tempDTI[1][0];
    DTI[4] = tempDTI[5][0];
    DTI[5] = tempDTI[2][0];

    ///delete
    for(unsigned long m=0 ; m<6 ; m++)
    {
        delete tempDTI[m];
    }
    delete tempDTI;
    for(unsigned long i=0 ; i<N ; i++)
    {
        delete Y[i];
    }
    delete Y;
    return DTI;
}

double** C_tensorMaker::makeTensor2(double* dwiArray/***N rows*/, double** Hpsi /**6 rows by N columns*/, double b /**b-value for dwiArray*/, double S0, unsigned long N)
{
    double* tempDTI = makeTensor(dwiArray, Hpsi, b, S0, N);
    ///store
    double** DTI = new double*[3];
    DTI[0] = new double[3];
    DTI[1] = new double[3];
    DTI[2] = new double[3];

    DTI[0][0] = tempDTI[0]; DTI[0][1] = tempDTI[1]; DTI[0][2] = tempDTI[2];
    DTI[1][0] = tempDTI[1]; DTI[1][1] = tempDTI[3]; DTI[1][2] = tempDTI[4];
    DTI[2][0] = tempDTI[2]; DTI[2][1] = tempDTI[4]; DTI[2][2] = tempDTI[5];

    delete tempDTI;
    return DTI;
}

rawData<double>* C_tensorMaker::getTensors(rawData<double>* dwiDataS, rawData<double>* dwiMask)
{
//    cout << "dwiDataS dimX: " << dwiDataS->DimX << " dwiMask dimX " << dwiMask->DimX << endl;
//    cout << "dwiDataS dimY: " << dwiDataS->DimY << " dwiMask dimY " << dwiMask->DimY << endl;
//    cout << "dwiDataS dimZ: " << dwiDataS->DimZ << " dwiMask dimZ " << dwiMask->DimZ << endl;

    if(dwiDataS->numDim!=4) return NULL;

    //allocate a rawData of dimension + 1: rawDataTensor (last dim has 6 components)
    rawData<double>* rawDataTensor = new rawData<double>(DICOM_DOUBLE, 4, dwiDataS->DimX, dwiDataS->DimY, dwiDataS->DimZ, 6);
    rawDataTensor->pixDimX = dwiDataS->pixDimX;
    rawDataTensor->pixDimY = dwiDataS->pixDimY;
    rawDataTensor->pixDimZ = dwiDataS->pixDimZ;
    rawDataTensor->pixDimT = dwiDataS->pixDimT;
    for(unsigned long i=0 ; i<rawDataTensor->DimX ; i++)
    {
        for(unsigned long j=0 ; j<rawDataTensor->DimY ; j++)
        {
            for(unsigned long k=0 ; k<rawDataTensor->DimZ ; k++)
            {
                if(dwiMask->raw3D[i][j][k]>0.5)
                {
                    ///fill Y vector

                    ///store
                    rawDataTensor->raw4D[i][j][k][0] = dwiDataS->raw4D[i][j][k][0]; //xx
                    rawDataTensor->raw4D[i][j][k][1] = dwiDataS->raw4D[i][j][k][1]; //xy
                    rawDataTensor->raw4D[i][j][k][2] = dwiDataS->raw4D[i][j][k][2]; //xz
                    rawDataTensor->raw4D[i][j][k][3] = dwiDataS->raw4D[i][j][k][3]; //yy
                    rawDataTensor->raw4D[i][j][k][4] = dwiDataS->raw4D[i][j][k][4]; //yz
                    rawDataTensor->raw4D[i][j][k][5] = dwiDataS->raw4D[i][j][k][5]; //zz
                }
                else
                {
                    rawDataTensor->raw4D[i][j][k][0] = 0.0;
                    rawDataTensor->raw4D[i][j][k][1] = 0.0;
                    rawDataTensor->raw4D[i][j][k][2] = 0.0;
                    rawDataTensor->raw4D[i][j][k][3] = 0.0;
                    rawDataTensor->raw4D[i][j][k][4] = 0.0;
                    rawDataTensor->raw4D[i][j][k][5] = 0.0;
                }
            }
        }
    }


    return rawDataTensor;
}


vector<rawData<double>*> C_tensorMaker::averageRawData(vector<vector<rawData<double>*> > myData)
{
    unsigned long N = myData.size();
    unsigned long Nrep = myData.at(0).size();
    unsigned long Ndim = myData.at(0).at(0)->numDim;
    unsigned long DimX = myData.at(0).at(0)->DimX;
    unsigned long DimY = myData.at(0).at(0)->DimY;
    unsigned long DimZ = myData.at(0).at(0)->DimZ;
    unsigned long DimT = myData.at(0).at(0)->DimT;
    vector<rawData<double>*> averageData(N,NULL);
    if(N==1)
    {
        for(unsigned long i=0 ; i<N ; i++)
        {
            averageData.at(i) = myData.at(i).at(0);
        }
    }
    else
    {
        for(unsigned long i=0 ; i<N ; i++)
        {
            averageData.at(i) = new rawData<double>(DICOM_DOUBLE, Ndim, DimX, DimY, DimZ, DimT);
            (*averageData.at(i)) = (*myData.at(i).at(0));
            for(unsigned long j=1 ; j<Nrep ; j++)
            {
                (*averageData.at(i)) += (*myData.at(i).at(j));
            }
            (*averageData.at(i)) /= ((double) Nrep);
        }
    }
    return averageData;
}
