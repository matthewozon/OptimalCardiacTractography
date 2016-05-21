#ifndef C_ICM_H
#define C_ICM_H

#include <C_minimization.h>


class C_ICM : public C_minimization
{
    public:
        C_ICM(C_graph* graph/**a pointer on an existing graph*/,\
              C_energy* U /**pointer to the energy related to the C_graph graph: assuming it's already set up*/,\
                 unsigned short energy_func=ICM_U1 /**kind of energy function: see settings_def.h*/,\
                  unsigned short move_kind=ICM_MOVE_1 /**kind of moves: see settings_def.h*/,\
                     bool rand_sweep=true);
        /** Default destructor */
        virtual ~C_ICM();

        /** working methods */
        virtual int run(void);
        virtual int iteration(void);

        //settings
        unsigned short m_energy_function;
        unsigned short m_move;
        bool m_rand_sweep;

        std::string m_fileName;
        std::string m_TAG;
    private:


        //energy related
        C_energy* m_U;
        double m_relative_energy;
        double globU;

        //ICM implementation
        int icm_iteration(void);
        int icm_iteration_v2(void);
};

#endif // C_ICM_H

