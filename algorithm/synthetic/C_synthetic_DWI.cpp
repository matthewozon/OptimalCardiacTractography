#include "C_synthetic_DWI.h"

C_synthetic_DWI::C_synthetic_DWI()
{
    //ctor
    unsigned long long int idum = ((unsigned long long int) time(NULL))*((unsigned long long int) time(NULL));
    while(idum==THE_NUMBER_WHICH_MUST_NOT_BE_NAMED)
    {
        idum = ((unsigned long long int) time(NULL))*((unsigned long long int) time(NULL));
    }
    toolRand = new C_toolbox_rand2(idum);
}

C_synthetic_DWI::~C_synthetic_DWI()
{
    //dtor
    delete toolRand;
}


diffInfo* C_synthetic_DWI::createB0Info()
{
    //
    diffInfo* dfi = new diffInfo;
    dfi->b = 0.0;
    dfi->echoTime = 80;
    dfi->repTime = 860;
    dfi->voxelSize[0] = 1.0; dfi->voxelSize[1] = 1.0; dfi->voxelSize[2] = 1.0;
    dfi->g[0] = 1.0; dfi->g[1] = 0.0; dfi->g[2] = 0.0;
    dfi->X_slice[0] = 1.0; dfi->X_slice[1] = 0.0; dfi->X_slice[2] = 0.0;
    dfi->Y_slice[0] = 0.0; dfi->Y_slice[1] = 1.0; dfi->Y_slice[2] = 0.0;
    dfi->Z_slice[0] = 0.0; dfi->Z_slice[1] = 0.0; dfi->Z_slice[2] = 1.0;
    return dfi;
}


vector< vector<diffInfo*> > C_synthetic_DWI::createDirectionsWithRepetion()
{
    ///1 1 0
    ///1 -1 0
    ///0 1 1
    ///0 1 -1
    ///1 0 1
    ///-1 0 1
    double staticDir[6][3] = {{1.0, 1.0, 0.0}, {1.0, -1.0, 0.0}, {0.0, 1.0, 1.0}, {0.0, 1.0, -1.0}, {1.0, 0.0, 1.0}, {-1.0, 0.0, 1.0}};
    short NgradientDirection = 6;
    short Nrepetition = 1;
    vector< vector<diffInfo*> > myDiffInfoVect( NgradientDirection, vector<diffInfo*>(Nrepetition,NULL) );
    for(short i=0 ; i<NgradientDirection ; i++)
    {
        for(short j=0 ; j<Nrepetition ; j++)
        {
            myDiffInfoVect.at(i).at(j) = new diffInfo;
            myDiffInfoVect.at(i).at(j)->b = 1.0;
            myDiffInfoVect.at(i).at(j)->echoTime = 80;
            myDiffInfoVect.at(i).at(j)->repTime = 860;
            myDiffInfoVect.at(i).at(j)->voxelSize[0] = 1.0; myDiffInfoVect.at(i).at(j)->voxelSize[1] = 1.0; myDiffInfoVect.at(i).at(j)->voxelSize[2] = 1.0;
            myDiffInfoVect.at(i).at(j)->g[0] = staticDir[i][0]/sqrt(SQR(staticDir[i][0]) + SQR(staticDir[i][1]) + SQR(staticDir[i][2]));
            myDiffInfoVect.at(i).at(j)->g[1] = staticDir[i][1]/sqrt(SQR(staticDir[i][0]) + SQR(staticDir[i][1]) + SQR(staticDir[i][2]));
            myDiffInfoVect.at(i).at(j)->g[2] = staticDir[i][2]/sqrt(SQR(staticDir[i][0]) + SQR(staticDir[i][1]) + SQR(staticDir[i][2]));
            myDiffInfoVect.at(i).at(j)->X_slice[0] = 1.0; myDiffInfoVect.at(i).at(j)->X_slice[1] = 0.0; myDiffInfoVect.at(i).at(j)->X_slice[2] = 0.0;
            myDiffInfoVect.at(i).at(j)->Y_slice[0] = 0.0; myDiffInfoVect.at(i).at(j)->Y_slice[1] = 1.0; myDiffInfoVect.at(i).at(j)->Y_slice[2] = 0.0;
            myDiffInfoVect.at(i).at(j)->Z_slice[0] = 0.0; myDiffInfoVect.at(i).at(j)->Z_slice[1] = 0.0; myDiffInfoVect.at(i).at(j)->Z_slice[2] = 1.0;
            //cout << myDiffInfoVect.at(i).at(j)->g[0] << " " << myDiffInfoVect.at(i).at(j)->g[1] << " " << myDiffInfoVect.at(i).at(j)->g[2] << endl;
        }
    }
    return myDiffInfoVect;
}








double C_synthetic_DWI::gaussDist(double mu, double sigma)
{
    double u = toolRand->doub();
    return ((u) > SMALL_NUM ? (mu + ( sigma * sqrt(-2*log(u)) ) * cos(2*Pi*toolRand->doub())) : (0));
}

double C_synthetic_DWI::riceDist(double mu1, double mu2, double sigma)
{
    double x = gaussDist(mu1, sigma);
    double y = gaussDist(mu2, sigma);
    double r = sqrt(SQR(x) + SQR(y));
    return r;
}

double C_synthetic_DWI::getSignalEnergy(rawData<double>* S, rawData<double>* mask)
{
    double E = 0.0;//, mu=0.0;
    unsigned long n = 0;
//    for(unsigned long i=0 ; i<S->DimX ; i++)
//    {
//        for(unsigned long j=0 ; j<S->DimY ; j++)
//        {
//            for(unsigned long k=0 ; k<S->DimZ ; k++)
//            {
//                if(mask->raw3D[i][j][k]>0.5)
//                {
//                    mu += S->raw3D[i][j][k];
//                    n++;
//                }
//            }
//        }
//    }
//    mu /= (double) n;
    for(unsigned long i=0 ; i<S->DimX ; i++)
    {
        for(unsigned long j=0 ; j<S->DimY ; j++)
        {
            for(unsigned long k=0 ; k<S->DimZ ; k++)
            {
                if(mask->raw3D[i][j][k]>0.5)
                {
                    //E += SQR(S->raw3D[i][j][k]-mu);
                    E += S->raw3D[i][j][k];//SQR(S->raw3D[i][j][k]);
                    n++;
                }
            }
        }
    }
    E /= (double) n;//(n-1);
    return E;//sqrt(E);
}


double C_synthetic_DWI::getSTDEVfromDWI(rawData<double>* S, rawData<double>* mask, double SNR)
{
    return sqrt(2/(4-Pi))*getSignalEnergy(S, mask)/sqrt(SNR);
    //return sqrt(getSignalEnergy(S, mask)/SNR)*sqrt(2/(4-Pi));
}


void C_synthetic_DWI::addNoise(rawData<double>* S, rawData<double>* mask, double SNR)
{
    double sigma = getSTDEVfromDWI(S, mask, SNR);
    for(unsigned long i=0 ; i<S->DimX ; i++)
    {
        for(unsigned long j=0 ; j<S->DimY ; j++)
        {
            for(unsigned long k=0 ; k<S->DimZ ; k++)
            {
                if(mask->raw3D[i][j][k]>0.5)
                {
                    //cout << S->raw3D[i][j][k] << ", " << riceDist(0.0, 0.0, sigma) << endl;
                    S->raw3D[i][j][k] += riceDist(0.0, 0.0, sigma);
                }
            }
        }
    }
    return;
}
