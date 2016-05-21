#ifndef C_SYNTHETIC_DWI_LV_H
#define C_SYNTHETIC_DWI_LV_H

#include <C_synthetic_DWI.h>


class C_synthetic_DWI_LV : public C_synthetic_DWI
{
    public:
        /** Default constructor */
        C_synthetic_DWI_LV(double alphaMin/*=-Pi/4.0*/, double alphaMax/*=Pi/4.0*/, double H/*=64.0*/, double a/*=16.0*/, double b/*=64.0*/, double A/*=9.0*/, double B/*=48.0*/);
        /** Default destructor */
        virtual ~C_synthetic_DWI_LV();

        virtual rawData<double>* createMask(diffInfo* g); //

        //B0
        virtual rawData<double>* createB0Data(diffInfo* g); //

        //data
        virtual vector< vector<rawData<double>*> > createDWIdata(vector< vector<diffInfo*> > dfi, vector<rawData<double>*> S0); //
        virtual rawData<double>* createOneDWI(diffInfo* g, rawData<double>* S0);
    private:
        double m_alphamin;
        double m_alphamax;
        double m_H;//height
        double m_a;//outsider ellipsoid
        double m_b;
        double m_A;//insider ellipsoid
        double m_B;
        double m_xmin, m_xmax;
        double m_ymin, m_ymax;
        double m_zmin, m_zmax;
};

#endif // C_SYNTHETIC_DWI_LV_H
