#ifndef C_TOOLBOX_SPLINE_H
#define C_TOOLBOX_SPLINE_H

#include <C_toolbox.h>
#include <struct.h>

#include <C_thread_spline.h>

class C_thread_spline;

class C_toolbox_spline : public C_toolbox
{
    friend class C_thread_spline; //
    public:
        /** Default constructor */
        C_toolbox_spline();
        /** Default destructor */
        virtual ~C_toolbox_spline();

        bool computeFiberSpline(CtlCtlStruct* m_ctlCtlStruct, CtlCtlStruct** m_ctlCtlStructSpline, unsigned long nb_sub_section=5);
    protected:
    private:
        //B-spline
        bool computeFiberSpline(CtlStruct* srcBuff, CtlStruct* destBuff, unsigned long nb_sub_section);
        bool computeBernsteinCoef(float u, unsigned long M, float* BernsteinCoef);
};

#endif // C_TOOLBOX_SPLINE_H
