#ifndef C_THREAD_MULTISCALE_H
#define C_THREAD_MULTISCALE_H

#include <C_thread.h>

#include <C_energy.h>
#include <C_temperature.h>
#include <C_SA.h>
#include <C_bundle.h>

#include <C_sampling.h>

#include <rawData.h>

#include <time.h>

//for memory check
#include "stdlib.h"
#include "stdio.h"
#include "string.h"



class C_thread_multiscale : public C_thread
{
    public:
        /** Default constructor */
        C_thread_multiscale(double resolution_factor, rawData<double>* maskDTI, rawData<double>* dataDTI, unsigned long NgradientDirection, unsigned long Nrepetition);
        /** Default destructor */
        virtual ~C_thread_multiscale();
        int curMem;
        virtual void execute();
    protected:
        double m_resolution_factor; //!< Member variable "resolution_factor"
        rawData<double>* m_maskDTI;
        rawData<double>* m_dataDTI;
        unsigned long m_NgradientDirection;
        unsigned long m_Nrepetition;
        int getValue(void);
        int parseLine(char* line);
};

#endif // C_THREAD_MULTISCALE_H
