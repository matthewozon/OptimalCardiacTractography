#ifndef C_THREAD_MULTISCALE_H
#define C_THREAD_MULTISCALE_H

#include <C_thread.h>

#include <C_toolbox_eigen_sym.h>

#include <C_energy.h>
#include <C_temperature.h>
#include <C_SA.h>
#include <C_bundle.h>

#include <C_sampling.h>

#include <rawData.h>

#include <C_ReadAnalyze75.h>

#include <time.h>

#include <cmath> //for pow function

//for memory check
#include "stdlib.h"
#include "stdio.h"
#include "string.h"

#define DWI_DOWN 0
#define DTI_DOWN 1
#define DTI_DOWN_LOG 2
#define GRAPH_DOWN 3


class C_thread_multiscale : public C_thread
{
    public:
        /** Default constructor */
        C_thread_multiscale(vector< vector<diffInfo*> > dfi_ /**info about dwi b!=0*/,\
                                          vector< vector<rawData<double>*> > dwiData_ /**dwi b!=0*/,\
                                           vector<diffInfo*> b0Info_ /**info about dwi b=0*/,\
                                            vector<rawData<double>*> b0Data_ /**dwi with b=0*/,\
                                             rawData<double>* dwiMask_ /**a mask for ROI*/);
        /** Default destructor */
        virtual ~C_thread_multiscale();
        virtual void execute();
        unsigned short downSamplingType;
        unsigned short numPyramidLevel;


    protected:
        //double m_resolution_factor; //!< Member variable "resolution_factor"
        int getValue(void);
        int parseLine(char* line);
        void logM(rawData<double>* rawDTI, rawData<double>* rawMask);
        void expM(rawData<double>* rawDTI, rawData<double>* rawMask);

        vector<diffInfo*> b0Info;
        vector< vector<diffInfo*> > dfi;
        vector<rawData<double>*> b0Data;
        rawData<double>* dwiMask;
        vector< vector<rawData<double>*> > dwiData;
};

#endif // C_THREAD_MULTISCALE_H
