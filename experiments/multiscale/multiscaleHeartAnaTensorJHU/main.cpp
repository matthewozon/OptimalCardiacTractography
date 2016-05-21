#include <iostream>
#include <dirent.h> //get file in directory
#include <algorithm> //sorting algorithm


////for memory check
#include "stdlib.h"
#include "stdio.h"
#include "string.h"

int parseLine(char* line)
{
    int i = strlen(line);
    while (*line < '0' || *line > '9') line++;
    line[i-3] = '\0';
    i = atoi(line);
    return i;
}

 int getValue() //Note: this value is in KB!
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


#include <C_ReadAnalyze75.h>
#include <C_tensorMaker.h>
#include <C_maskMaker.h>

#include <C_toolbox_interpolate.h>

#include <C_energy.h>
#include <C_temperature.h>
#include <C_SA.h>
#include <C_bundle.h>

#include <C_sampling.h>

#include <C_thread_multiscale.h>

#define ISOTROPIC_SAMPLING
//#define TEST_GRAPH_ISOTROPIC
#ifdef TEST_GRAPH_ISOTROPIC
    #define NTHREAD 1
#else
    #define NTHREAD 6
#endif
//#define ON_MY_COMPUTER

using namespace std;

int main(int argc, char *argv[])
{
    if(argc==1) ///ANA test
    {
        int curMem = getValue();
        cout << endl;
        cout << "vMem = " << curMem << endl;

        ///get gradient direction
        C_ReadAnalyze75* ANAtool = new C_ReadAnalyze75();
        rawData<double>* myTensos = ANAtool->getRawData("DT080803");


        ///creating mask from several gradient directions
        C_maskMaker* maskTool = new C_maskMaker();
        rawData<double>* segRawData;
        segRawData = ANAtool->getRawData("ROI");
        (*segRawData) = (*segRawData)>0.5;
        segRawData->removeSingle();
        delete maskTool;

        #ifdef ISOTROPIC_SAMPLING
        C_sampling<double>* toolISOTROPE = new C_sampling<double>();
        samplingFactor* m_SF = new samplingFactor();
        ///fill sampling factors so that you can have all graph resolution
        m_SF->Xfactor = 1.0;//128.0/((double) segRawData->DimX);
        m_SF->Yfactor = 1.0;//128.0/((double) segRawData->DimY);
        m_SF->Zfactor = segRawData->pixDimZ*(m_SF->Xfactor)/segRawData->pixDimX; //allB0Data.at(0)->pixDimZ/allB0Data.at(0)->pixDimX;//generate a bug in graph generation: maybe a comparason < instead of <=
        m_SF->Tfactor = 1.0;
        cout << "Xfactor: " << m_SF->Xfactor << " Yfacto: " << m_SF->Yfactor << " Zfactor: " << m_SF->Zfactor << " Tfactor: " << m_SF->Tfactor << endl;

        rawData<double>* tempRawData = toolISOTROPE->samplingNearestNeighbor(segRawData, m_SF);
        delete segRawData;
        segRawData = new rawData<double>(DICOM_DOUBLE, tempRawData->numDim, tempRawData->DimX, tempRawData->DimY, tempRawData->DimZ, tempRawData->DimT);
        (*segRawData) = (*tempRawData)>0.5;
        segRawData->removeSingle();

        delete tempRawData;
        tempRawData = toolISOTROPE->samplingNearestNeighbor(myTensos, m_SF);
        delete myTensos;
        myTensos = new rawData<double>(DICOM_DOUBLE, tempRawData->numDim, tempRawData->DimX, tempRawData->DimY, tempRawData->DimZ, tempRawData->DimT);
        *myTensos = *tempRawData;

        dsr* ptrDSR = new dsr;
        ptrDSR->hk.extents = 16384;
        //ptrDSR->hk.db_name = (char*) "private";
        ptrDSR->hk.regular = (char) 'r';
        ptrDSR->dime.datatype = DT_DOUBLE;
        ptrDSR->dime.bitpix = 8*sizeof(double);
        ptrDSR->dime.dim[0] = 3;
        ptrDSR->dime.dim[1] = segRawData->DimX; ptrDSR->dime.dim[2] = segRawData->DimY; ptrDSR->dime.dim[3] = segRawData->DimZ;
        ptrDSR->dime.pixdim[0] = 3;
        ptrDSR->dime.pixdim[1] = segRawData->pixDimX; ptrDSR->dime.pixdim[2] = segRawData->pixDimY; ptrDSR->dime.pixdim[3] = segRawData->pixDimZ;
        ANAtool->saveAna( (const char*) "test_isotrope_data", ptrDSR, segRawData);

        //return 3;
        delete tempRawData;
        delete m_SF;
        delete toolISOTROPE;
        #endif
        delete ANAtool;

        ///fill sampling factors so that you can have all graph resolution
        #ifdef TEST_GRAPH_ISOTROPIC
        //unsigned short NTHREAD = 1;
        double resolutionFactor[NTHREAD] = {0.1875};
        #else
        //unsigned short NTHREAD = 7;
            #ifdef ON_MY_COMPUTER
            double resolutionFactor[NTHREAD] = {0.1, 0.11574, 0.125, 0.15, 0.175, 0.2};//{1.0/4.0, 1.0/3.5, 1.0/3.0, 1.0/2.5, 1.0/2.0, 1.0/1.5, 1.0/1.0};
            #else
            double resolutionFactor[NTHREAD] = {0.25, 0.3, 0.35, 0.4, 0.45, 0.5};
            #endif
        #endif

        C_thread_multiscale** THREAD_RUN = new C_thread_multiscale*[NTHREAD];
        for(unsigned short i=0 ; i<NTHREAD ; i++)
        {
            ///change resolution at will
            THREAD_RUN[i] = new C_thread_multiscale(resolutionFactor[i], segRawData, myTensos, 0, 0);
            THREAD_RUN[i]->curMem = curMem;
        }

        //launch threads
        for(unsigned short i=0 ; i<NTHREAD ; i++)
        {
            //cout << "start thread " << i << endl;
            THREAD_RUN[i]->start();
        }
        //wait for all threads
        for(unsigned short i=0 ; i<NTHREAD ; i++)
        {
            //cout << "wait thread " << i << endl;
            THREAD_RUN[i]->wait_for_exit();
        }
        //delete threads
        for(unsigned short i=0 ; i<NTHREAD ; i++)
        {
            delete THREAD_RUN[i];
        }
        delete THREAD_RUN;
        delete segRawData;
        delete myTensos;
    }
    return 60;
}
