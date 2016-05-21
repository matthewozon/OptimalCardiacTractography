#include <C_cooling_schedule.h>

C_cooling_schedule::C_cooling_schedule(double Tinit, double Tend, unsigned long nbIter, unsigned long stepLength)
{
    //ctor
    m_Tinit = Tinit;
    m_Tend = Tend;
    m_nbIter = nbIter;
    m_current_step = 0;
    m_stepLength = stepLength;
    //step_exp = Tinit*(Tend/Tinit)^(i/n)
    step_exp = pow(m_Tend/m_Tinit,1/((double)(m_nbIter-1)));
    step_exp_block = pow(m_Tend/m_Tinit,1/((double) floor( ((double) m_nbIter-1)/((double) m_stepLength)) ));
}

C_cooling_schedule::~C_cooling_schedule()
{
    //dtor
}


/** utilities */
bool C_cooling_schedule::setTinit(double newTinit) //should compute new exp_step and exp_step_block and...
{
    m_Tinit = newTinit;
    step_exp = pow(m_Tend/m_Tinit,1/((double)(m_nbIter-1)));
    step_exp_block = pow(m_Tend/m_Tinit,1/((double) floor( ((double) m_nbIter-1)/((double) m_stepLength)) ));
    return true;
}
bool C_cooling_schedule::setTend(double newTend) //should compute new exp_step and exp_step_block and...
{
    m_Tend = newTend;
    step_exp = pow(m_Tend/m_Tinit,1/((double)(m_nbIter-1)));
    step_exp_block = pow(m_Tend/m_Tinit,1/((double) floor( ((double) m_nbIter-1)/((double) m_stepLength)) ));
    return true;
}
bool C_cooling_schedule::setnbIter(unsigned long newnbIter) //should compute new exp_step and exp_step_block and...
{
    m_nbIter = newnbIter;
    step_exp = pow(m_Tend/m_Tinit,1/((double)(m_nbIter-1)));
    step_exp_block = pow(m_Tend/m_Tinit,1/((double) floor( ((double) m_nbIter-1)/((double) m_stepLength)) ));
    return true;
}
bool C_cooling_schedule::setCurrentIter(unsigned long newIter)
{
    m_current_step = newIter;
    return true;
}
bool C_cooling_schedule::resetCurrentIter(void)
{
    m_current_step = 0;
    return true;
}
bool C_cooling_schedule::setStepLength(unsigned long newStepLen) //should compute new exp_step_block
{
    m_stepLength = newStepLen;
    step_exp_block = pow(m_Tend/m_Tinit,1/((double) floor( ((double) m_nbIter-1)/((double) m_stepLength)) ));
    return true;
}
double C_cooling_schedule::getTinit(void)
{
    return m_Tinit;
}
double C_cooling_schedule::getTend(void)
{
    return m_Tend;
}
unsigned long C_cooling_schedule::getnbIter(void)
{
    return m_nbIter;
}
unsigned long C_cooling_schedule::getCurrentIter(void)
{
    return m_current_step;
}
unsigned long C_cooling_schedule::getStepLen(void)
{
    return m_stepLength;
}

/** global call*/
double C_cooling_schedule::schedule(unsigned long iter, unsigned long kind)
{
    double T;
    if(kind==COOLING_EXP){
        T = expSchedule(iter);
    }else{//if(kind==COOLING_EXP_BLOCK)
        T = expBlockSchedule(iter);
    }
    return T;
}
double C_cooling_schedule::schedule(unsigned long kind)
{
    double T;
    if(kind==COOLING_EXP){
        T = expSchedule();
    }else{//if(kind==COOLING_EXP_BLOCK)
        T = expBlockSchedule();
    }
    return T;
}

/** cooling function */
double C_cooling_schedule::expSchedule(unsigned long iter)
{
    return m_Tinit*pow(step_exp,(double) iter);
}
double C_cooling_schedule::expSchedule(void)
{
    //update first current step (use current_step-1 as current step right after)
    m_current_step++;
    return m_Tinit*pow(step_exp,(double) m_current_step-1);
}
double C_cooling_schedule::expBlockSchedule(unsigned long iter)
{
    double tmp = floor( ((double) iter)/((double) m_stepLength));
    return m_Tinit*pow(step_exp_block,tmp);
}
double C_cooling_schedule::expBlockSchedule(void)
{
    m_current_step++;
    double tmp = floor( ((double) (m_current_step-1))/((double) m_stepLength));
    return m_Tinit*pow(step_exp_block,tmp);
}
