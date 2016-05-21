#ifndef C_THREAD_SPLINE_H
#define C_THREAD_SPLINE_H

#include <C_thread.h>
#include <struct.h>
#include <C_toolbox_spline.h>

class C_toolbox_spline;

class C_thread_spline : public C_thread
{
    public:
        /** Default constructor */
        C_thread_spline(C_toolbox_spline* _Spline);
        /** Default destructor */
        virtual ~C_thread_spline();
        unsigned long idxStart;
        unsigned long idxEnd;

        CtlCtlStruct* m_ctlCtlStruct;
        CtlCtlStruct** m_ctlCtlStructSpline;
        unsigned long nb_sub_section;

        virtual void execute();

    protected:
        C_toolbox_spline* m_Spline;
};

#endif // C_THREAD_SPLINE_H
