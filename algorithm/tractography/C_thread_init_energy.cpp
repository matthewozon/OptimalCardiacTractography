#include "C_thread_init_energy.h"

C_thread_init_energy::C_thread_init_energy(C_energy* U, unsigned long idxStart, unsigned long idxEnd, unsigned long kind):C_thread()
{
    m_U = U;
    m_idxStart = idxStart;
    m_idxEnd = idxEnd;
    m_kind = kind;
    subTotU = 0.0;
}

C_thread_init_energy::~C_thread_init_energy()
{
    //dtor
}


void C_thread_init_energy::execute()
{
    subTotU = m_U->initEnergy(m_kind, m_idxStart, m_idxEnd);
    return;
}
