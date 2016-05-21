#ifndef C_SA_H
#define C_SA_H

#include <C_minimization.h>
#include <string.h>

#include <fstream>
#include <sstream>

class C_SA : public C_minimization
{
    public:
        /** Default constructor */
        C_SA(C_graph* graph /**a pointer on an existing graph*/,\
             C_energy* U /**pointer to the energy related to the C_graph graph: assuming it's already set up*/,\
              double Tinit, double Tend  /**initial and final temperature*/,\
                unsigned long nbIter /**number of iterations, which is multiply by the number of edges in the graph=> actual nb of iter = nbIter*nb_edges*/,\
                    unsigned short cooling_sch=SA_EXP_BLOCK_COOLING /**kind of cooling schedule: see settings_def.h*/,\
                     unsigned short energy_func=SA_U1 /**kind of energy function: see settings_def.h*/,\
                      unsigned short move_kind=SA_MOVE_1 /**kind of moves: see settings_def.h*/);
        /** Default destructor */
        virtual ~C_SA();

        /** utilities */
        bool setStepLen(unsigned long stepLength);
        double getTinit(void);
        double getTend(void);
        unsigned long getnbIter(void);
        //unsigned short getEnergyFunction(void);

        std::string baseName;
        /** working methods */
        virtual int run(void);
        virtual int iteration(void);

        bool flipRate;
        bool energyMeasure;
        double alphaMixture;
        double m_varSigma;
        bool m_voxelConstraint;
        bool m_moveFromCenter;
        string m_expTAG;


        unsigned short m_cooling_schedule;
        unsigned short m_energy_function;
        unsigned short m_move;

        long int nbFlip;

    protected:
        //random parameter
        unsigned long long int idum;
        std::string m_fileName;

    private:
        double m_relative_energy;
        double globU;

        //temperature of the current iteration
        double T;
        double m_Tinit;
        double m_Tend;
        unsigned long m_nbIter;

        //tools
        C_toolbox_rand2 *m_rand;
        C_move *m_moves;
        C_energy *m_U;
        C_cooling_schedule *Tcooling;

        //SA implementation
        int iterationM1(void);
        int iterationM2(void);
        int iterationM1M2(void);
        int iterationM4(void);
        int iterationP1(void);
        int iterationP1M1(void);
        int iterationP1M2(void);

        //int iterationSC(void);
        //int icm_iteration(void);

        //tool
        string baseSaveName(void);
        void writeMfile(double startU, double endU);
        void writeEnergiesAndFlipRate(fstream* filestr, double countflips, double modulo);


};

#endif // C_SA_H
