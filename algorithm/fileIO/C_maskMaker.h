#ifndef C_MASKMAKER_H
#define C_MASKMAKER_H

#include <C_kmean2D.h>
#include <C_kmean3D.h>
#include <C_sampling.h>
#include <C_ReadAnalyze75.h>
#include <C_toolbox_eigen_sym.h>
#include <settings_def.h>

#include <iostream>

using namespace std;



class C_maskMaker
{
    public:
        /** Default constructor */
        C_maskMaker();
        /** Default destructor */
        virtual ~C_maskMaker();
        rawData<double>* createMaskFromDWIs(vector<vector<rawData<double>*> > allRawData, samplingFactor* SF, unsigned long SAMPLING_TYPE);
        rawData<double>* createMaskFromB0(vector<rawData<double>*> allRawData, samplingFactor* SF, unsigned long SAMPLING_TYPE);
        rawData<double>* createMaskFromTensors(rawData<double>* rawTensor, double FAmin, double FAmax);

        rawData<double>* getMaskFromANAFile(const char* filenameData);
};

#endif // C_MASKMAKER_H
