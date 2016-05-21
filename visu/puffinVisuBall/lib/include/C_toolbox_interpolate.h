/**
In this class should occur all interpolations, either on DWI or DTI and the result should be return in a interpolatePoints structure with a result in vector array

Case that must be implemented:
    - interpolation from tensor (no matter if the tensor map is down/up-sampled or not, or if the DWI data are down/up-sampled)
    - NOT: interpolation from DWI must not be implemented in this class because it must include tensor calculation which is implemented in another class
*/


#ifndef C_INTERPOLATE_H
#define C_INTERPOLATE_H

#include <C_toolbox.h>
#include <C_graph.h>
#include <rawData.h>
#include <C_point_array_data.h>
#include <C_toolbox_eigen_sym.h>
#include <C_tensorMaker.h>




class C_toolbox_interpolate : public C_toolbox  //ok for allocations
{
    public:
        /** Default constructor */
        //C_toolbox_interpolate(C_graph* G, double r=1.0);
        C_toolbox_interpolate(rawData<double>* dwiOrDtiData, rawData<double>* maskRawData, double r=1.0);
        /** Default destructor */
        virtual ~C_toolbox_interpolate();

        //the called function that run a given interpolation on a given set of points
        void regularInterpolatePoint(C_point_array_data* p, double* vStart, double* vEnd); ///init the point on which to compute the interpolation (liner interpolation between two points)
        void regularInterpolatePoint(C_point_array_data* p, Point* vStart, Point* vEnd);

        unsigned long run(C_point_array_data* p, unsigned long METHOD, unsigned long OBJECT_TYPE, bool yes, double** Hpsi=NULL, double b=1000.0); ///if yes==true => C_point_array_data must initialize vector3Interpolated, eig1, eig2 and eig3

        //set settings
        void setSigma(double newSigma){sigma = newSigma;}
    protected:
        double R;  //radius of embeding/surounding ball
        double sigma;  //parameter for phi functions (careful, this is not obiously the same meaning for all functions)
        //double x0, y0, z0; //patient origine (often taken to be (0,0,0))
        //double dx, dy, dz; //voxel size

    private:
        C_graph* m_graph;
        rawData<double>* m_dwiOrDtiData;
        rawData<double>* m_maskRawData;
        bool interpolateFromRawDataDTI(C_point_array_data* p, unsigned long METHOD, bool yes); ///assume 3D +1D for tensor storage
        bool interpolateFromRawDataDWI(C_point_array_data* p, double** Hpsi, unsigned long NgradDir, double b, unsigned long METHOD, bool yes); ///assume 3D + 1D for diffusion storage (averaging alredy computed)
        bool interpolateFromRawDataDTI(double* p, unsigned long METHOD, double* unitVect, double* FA);
        bool interpolateFromRawDataDWI(double* p, double** Hpsi, unsigned long NgradDir, double b, unsigned long METHOD, double* unitVect, double* FA);


        //phi functions
        double gaussianPhi(double r);
        double gaussianPhiNorm(double r);
        double invAbsPhi(double r);
        double cauchyPhi(double r);
        double cauchyPhiNorm(double r);
};

#endif // C_INTERPOLATE_H
