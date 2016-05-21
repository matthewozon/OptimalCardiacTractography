#ifndef C_SYNTHETIC_DWI_PATHO_H
#define C_SYNTHETIC_DWI_PATHO_H

#include <C_synthetic_DWI.h>


class C_synthetic_DWI_patho : public C_synthetic_DWI
{
    public:
        /** Default constructor */
        C_synthetic_DWI_patho(double X, double Y, double Z, double xc, double yc, double zc, double r);
        /** Default destructor */
        virtual ~C_synthetic_DWI_patho();

        //mask
        virtual rawData<double>* createMask(diffInfo* g); //

        //B0
        virtual rawData<double>* createB0Data(diffInfo* g); //

        //data
        virtual vector< vector<rawData<double>*> > createDWIdata(vector< vector<diffInfo*> > dfi, vector<rawData<double>*> S0); //
    protected:
        virtual rawData<double>* createOneDWI(diffInfo* g, rawData<double>* S0);

        double m_X, m_Y, m_Z;
        double m_xc, m_yc, m_zc;
        double m_r;
        unsigned long dimX, dimY, dimZ;
};

#endif // C_SYNTHETIC_DWI_PATHO_H
