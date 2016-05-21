#ifndef C_SYNTHETIC_DWI_CYLINDER_H
#define C_SYNTHETIC_DWI_CYLINDER_H

#include <C_synthetic_DWI.h>


class C_synthetic_DWI_cylinder : public C_synthetic_DWI
{
    public:
        /** Default constructor */
        C_synthetic_DWI_cylinder(double xmin_, double xmax_, double ymin_, double ymax_, double zmin_, double zmax_, double rmin_, double rmax_, double alphamin_, double alphamax_);
        /** Default destructor */
        virtual ~C_synthetic_DWI_cylinder();
        virtual rawData<double>* createMask(diffInfo* g); //

        //B0
        virtual rawData<double>* createB0Data(diffInfo* g); //

        //data
        virtual vector< vector<rawData<double>*> > createDWIdata(vector< vector<diffInfo*> > dfi, vector<rawData<double>*> S0); //
        virtual rawData<double>* createOneDWI(diffInfo* g, rawData<double>* S0);
    protected:
        double xmin;
        double xmax;
        double ymin;
        double ymax;
        double zmin;
        double zmax;
        double rmin;
        double rmax;
        double alphamin;
        double alphamax;
};

#endif // C_SYNTHETIC_DWI_CYLINDER_H
