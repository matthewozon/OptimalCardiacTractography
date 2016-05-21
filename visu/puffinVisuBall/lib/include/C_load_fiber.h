#ifndef C_LOAD_FIBER_H
#define C_LOAD_FIBER_H

#include <struct.h>
#include <C_toolbox_spline.h>
#include <fstream>
#include <settings_def.h>
//#include <structs.h>
#include <struct.h>


class C_load_fiber
{
    public:
        /** Default constructor */
        C_load_fiber();
        /** Default destructor */
        virtual ~C_load_fiber();

        CtlCtlStruct* readFiber(string fileNameFib);///, bool rawData = true);
        CtlCtlStruct* readFiber(string fileNameFib, bool rawData);
};

#endif // C_LOAD_FIBER_H
