#ifndef C_SYNTHETIC_DWI_H
#define C_SYNTHETIC_DWI_H

#include <rawData.h>
#include <vector>
#include <struct.h>
#include <MACRO_DEF.h>
#include <C_toolbox_rand2.h>

using namespace std;

class C_synthetic_DWI
{
    public:
        /** Default constructor */
        C_synthetic_DWI();
        /** Default destructor */
        virtual ~C_synthetic_DWI();

        //gradient info
        vector< vector<diffInfo*> > createDirectionsWithRepetion(/**actually, there is only on repetition*/);

        //mask creation
        virtual rawData<double>* createMask(diffInfo* g)=0; //, double xmin, double xmax, double ymin, double ymax, double zmin, double zmax, double rmin, double rmax

        //B0
        virtual rawData<double>* createB0Data(diffInfo* g)=0; //, double xmin, double xmax, double ymin, double ymax, double zmin, double zmax, double rmin, double rmax
        diffInfo* createB0Info();

        //data
        virtual vector< vector<rawData<double>*> > createDWIdata(vector< vector<diffInfo*> > dfi, vector<rawData<double>*> S0)=0; //, double xmin, double xmax, double ymin, double ymax, double zmin, double zmax, double rmin, double rmax, double alphamin, double alphamax
        virtual rawData<double>* createOneDWI(diffInfo* g, rawData<double>* S0)=0; //, double xmin, double xmax, double ymin, double ymax, double zmin, double zmax, double rmin, double rmax, double alphamin, double alphamax

        //noise addition
        void addNoise(rawData<double>* S, rawData<double>* mask, double SNR);
        //noise creation
    protected:
        double getSignalEnergy(rawData<double>* S, rawData<double>* mask);
        double getSTDEVfromDWI(rawData<double>* S, rawData<double>* mask, double SNR);

        C_toolbox_rand2* toolRand;
        double gaussDist(double mu, double sigma);
        double riceDist(double mu1, double mu2, double sigma);
};

#endif // C_SYNTHETIC_DWI_H
