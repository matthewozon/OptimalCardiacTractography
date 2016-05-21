#ifndef C_COOLING_SCHEDULE_H
#define C_COOLING_SCHEDULE_H
#include <settings_def.h>
#include <iostream>
#include <math.h>

using namespace std;


class C_cooling_schedule
{
    public:
        /** Default constructor */
        C_cooling_schedule(double Tinit, double Tend, unsigned long nbIter, unsigned long stepLength=1);
        /** Default destructor */
        ~C_cooling_schedule();

        /** methods to modify attributes (use only if necessary!!!) */
        bool setTinit(double newTinit);
        bool setTend(double newTend);
        bool setnbIter(unsigned long newnbIter);
        bool setCurrentIter(unsigned long newIter);
        bool resetCurrentIter(void);
        bool setStepLength(unsigned long newStepLen);
        double getTinit(void);
        double getTend(void);
        unsigned long getnbIter(void);
        unsigned long getCurrentIter(void);
        unsigned long getStepLen(void);

        /**  calling function*/
        double schedule(unsigned long iter, unsigned long kind);
        double schedule(unsigned long kind);

    protected:
    private:
        /** cooling function */
        double expSchedule(unsigned long iter);
        double expSchedule(void);
        double expBlockSchedule(unsigned long iter);
        double expBlockSchedule(void);
        //double linRandSchedule(unsigned long iter); //maybe useless
        //double linRandSchedule(void); //maybe useless

        //attribute concening the cooling schedule
        double m_Tinit; //initial temperature
        double m_Tend; //final temperature
        unsigned long m_nbIter; //number of iteration/step to reach the final temperature (Tend) starting from Tinit
        unsigned long m_current_step;
        unsigned long m_stepLength;
        double step_exp;
        double step_exp_block;
};

#endif // C_COOLING_SCHEDULE_H
