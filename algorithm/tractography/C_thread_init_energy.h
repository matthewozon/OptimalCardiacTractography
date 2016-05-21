#ifndef C_THREAD_INIT_ENERGY_H
#define C_THREAD_INIT_ENERGY_H

#include <C_thread.h>
#include <C_energy.h>

class C_energy;

class C_thread_init_energy : public C_thread
{
    //friend class C_energy;
    public:
        /** Default constructor */
        C_thread_init_energy(C_energy* U, unsigned long idxStart, unsigned long idxEnd, unsigned long kind);
        /** Default destructor */
        virtual ~C_thread_init_energy();
        virtual void execute();
        double subTotU;
    protected:
    private:
        unsigned long m_idxStart;
        unsigned long m_idxEnd;
        unsigned long m_kind;
        C_energy* m_U;
};

#endif // C_THREAD_INIT_ENERGY_H
