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
#include <C_synthetic_DWI_LV.h>

#include <C_toolbox_interpolate.h>

#include <C_energy.h>
#include <C_temperature.h>
#include <C_SA.h>
#include <C_bundle.h>

#include <C_sampling.h>

#include <C_thread_multiscale.h>


struct tractData
{
    rawData<double>* dwiMask;
    rawData<double>* tensorData;
};

tractData createDataAndComputeTensor(char* SNRc);


#define NTHREAD 4

using namespace std;

int main(int argc, char *argv[])
{
    if(argc!=2)
    {
        cout << "command: " << argv[0] << " should be run as: " << argv[0] << " #pyramid level<unsigned int>" << endl;
        return -1;
    }

    ///allocate tool for synthetic data creation
    C_synthetic_DWI_LV toolDWI(-Pi/4.0, Pi/4.0, 4.0+2.0*64.0/2.0, 2.0*16.0/*/2.0*/, 2.0*64.0/2.0, 2.0*9.0/*/2.0*/, 2.0*48.0/2.0);
    //C_synthetic_DWI_cylinder toolDWI(XMIN, XMAX, YMIN, YMAX, ZMIN, ZMAX, RMIN, RMAX, ALPHAMIN, ALPHAMAX);

    ///create structures to store DWI information (gradient direction, b-value... cf struct definition)
    vector<diffInfo*> b0Info;
    b0Info.push_back(toolDWI.createB0Info());
    vector< vector<diffInfo*> > dfi = toolDWI.createDirectionsWithRepetion(); //create six direction
    ///actually create data
    vector<rawData<double>*> b0Data;
    b0Data.push_back(toolDWI.createB0Data(b0Info.at(0))); // create only one b0, because there is no need for more in synthetic case
    rawData<double>* dwiMask = new rawData<double>(DICOM_DOUBLE, 3, b0Data.at(0)->DimX, b0Data.at(0)->DimY, b0Data.at(0)->DimZ);// allocate mask
    (*dwiMask) = (*(b0Data.at(0))>0.0); ///for real data the creation of the mask is different. Here, is S0 is strictly positive, it is in mask
    dwiMask->pixDimX = 1.0;
    dwiMask->pixDimY = 1.0;
    dwiMask->pixDimZ = 1.0;
    dwiMask->pixDimT = 1.0;
    vector< vector<rawData<double>*> > dwiData = toolDWI.createDWIdata(dfi, b0Data);


    ///fill sampling factors so that you can have all graph resolution
    //short int numPyramidLevel[NTHREAD] = {0, 0, 0, 0, 1, 1, 1, 1};//, 2, 2, 2, 2};
    short int downSamplingType[NTHREAD] = {DWI_DOWN, DTI_DOWN, DTI_DOWN_LOG, GRAPH_DOWN};//, DWI_DOWN, DTI_DOWN, DTI_DOWN_LOG, GRAPH_DOWN};

    //short int numPyramidLevel[NTHREAD] = {2};//, 1, 1, 1, 2, 2, 2, 2};
    //short int downSamplingType[NTHREAD] = {DTI_DOWN_LOG};//DWI_DOWN, DTI_DOWN, DTI_DOWN_LOG, GRAPH_DOWN, DWI_DOWN, DTI_DOWN, DTI_DOWN_LOG, GRAPH_DOWN};

    C_thread_multiscale** THREAD_RUN = new C_thread_multiscale*[NTHREAD];
    for(unsigned short i=0 ; i<NTHREAD ; i++)
    {
        ///change resolution at will
        THREAD_RUN[i] = new C_thread_multiscale(dfi, dwiData, b0Info, b0Data, dwiMask);
        THREAD_RUN[i]->downSamplingType = downSamplingType[i];
        THREAD_RUN[i]->numPyramidLevel = atoi(argv[1]);//numPyramidLevel[i];
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
    delete [] THREAD_RUN;



    ///should delete b0Data;
    ///should delete b0Info
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

    return 60;
}
