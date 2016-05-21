#ifndef C_STREAMLINE_H
#define C_STREAMLINE_H

#include <rawData.h>
#include <C_toolbox_eigen_sym.h>
#include <C_bundle.h>
#include <C_toolbox_interpolate.h>
#include <C_point_array_data.h>
#include <sstream>

#define MINIMAL   0
#define EVOLUTION 1


class C_streamline
{
    public:
        /** Default constructor */
        C_streamline(rawData<double>* rawTensors, rawData<double>* rawMask, double wpunct=0.2);
        /** Default destructor */
        virtual ~C_streamline();
        bool run2(double LEVEL_SAMPLING=1, double dt=0.1, double MINSPEED=0.25, double MAXLENGTH=20);
        bool run2(unsigned short NB_SEED=1, double dt=0.1, double MINSPEED=0.25, double MAXLENGTH=20); //OK without measure
        rawData<double>* m_rawTensorsCOmpare;
    protected:
        int id;
        double m_wpunct;

        rawData<double>* m_rawTensors;
        rawData<double>* m_rawMask;
        unsigned long nbPoint;

        bool saveFiber(std::string fileNameFib, vector< vector<BundlePoint> > fibers); //OK


        ///real streamline
        bool isPartOfVolume(double x, double y, double z); //re-implemented
        double* F(double x, double y, double z, double t); ///dr(x,y,z,t)/dt = F(x,y,z,t) //re-implemented
        double* EULERstep(double x, double y, double z, double t, double dt, bool reverse); //OK
        double* RK4step(double x, double y, double z, double t, double dt, bool reverse); //OK
        double* RK4K1(double x, double y, double z, double t, double dt); //OK
        double* RK4K2(double x, double y, double z, double k1x, double k1y, double k1z, double t, double dt); //OK
        double* RK4K3(double x, double y, double z, double k2x, double k2y, double k2z, double t, double dt); //OK
        double* RK4K4(double x, double y, double z, double k3x, double k3y, double k3z, double t, double dt); //OK
        vector<BundlePoint> streamLine(double x0, double y0, double z0, double dt, bool reverse, double MINSPEED, double MAXLENGTH); //OK
        vector<BundlePoint> streamLine(double x0, double y0, double z0, double dt, double MINSPEED, double MAXLENGTH); //OK
        vector< vector<BundlePoint> > streamLines(unsigned short NB_SEED, double dt, double MINSPEED, double MAXLENGTH); //to be modified and multithread
        vector< vector<BundlePoint> > streamLines(double LEVEL_SAMPLING, double dt, double MINSPEED, double MAXLENGTH);

        //measures must be re-implemented without graph
        void measure(vector< vector<BundlePoint> > V, double* L/**mean length*/, double* S/**mean standard deviavtion of length*/, double* E/**mean fiber error*/);
        void measure(vector<BundlePoint>  v, double* L/**length*/, double* E/**fiber error*/);
};

#endif // C_STREAMLINE_H
