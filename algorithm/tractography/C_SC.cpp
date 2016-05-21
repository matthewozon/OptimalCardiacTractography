#include <C_SC.h>

C_SC::C_SC(C_graph* graph /**a pointer on an existing graph*/,\
             C_energy* U /**pointer to the energy related to the C_graph graph: assuming it's already set up*/,\
              double Tinit, double Tend  /**initial and final temperature*/,\
                unsigned long nbIter /**number of iterations, which is multiply by the number of edges in the graph=> actual nb of iter = nbIter*nb_edges*/,\
                 unsigned short cooling_sch /**kind of cooling schedule: see settings_def.h*/,\
                  unsigned short energy_func/**kind of energy function: see settings_def.h*/,\
                   unsigned short move_kind/**kind of moves: see settings_def.h*/):C_minimization(graph)
{
    //ctor
    m_Tinit = Tinit;
    m_Tend = Tend;
    m_nbIter = nbIter;
    m_expTAG = "";

    m_cooling_schedule = cooling_sch;
    m_energy_function = energy_func;
    m_move = move_kind;
    //m_rand_sweep = rand_sweep;
    m_moves = new C_move(m_graph);
    m_U = U;


    idum = (unsigned long long int) time(NULL);
    while(idum==THE_NUMBER_WHICH_MUST_NOT_BE_NAMED)
    {
        idum = (unsigned long long int) time(NULL);
    }
    m_rand = new C_toolbox_rand2(idum);
    m_relative_energy = 0;

    Tcooling = new C_cooling_schedule(m_Tinit,m_Tend, m_nbIter*m_graph->getLenEdges(), m_graph->getLenEdges());
    m_fileName = "new_file.cformat";
    alphaMixture = 0.5;
    flipRate = true;
    energyMeasure = true;
    nbFlip = 0;
}

C_SC::~C_SC()
{
    //dtor
    delete m_moves;
    delete m_rand;
    delete Tcooling;
}



/** utilities */
bool C_SC::setStepLen(unsigned long stepLength)
{
    return Tcooling->setStepLength(stepLength);
}
double C_SC::getTinit(void)
{
    return m_Tinit;
}
double C_SC::getTend(void)
{
    return m_Tend;
}
unsigned long C_SC::getnbIter(void)
{
    return m_nbIter;
}
unsigned short C_SC::getEnergyFunction(void)
{
    return m_energy_function;
}

/** working methods */
int C_SC::run(void)
{
    ostringstream oss;
    oss << baseSaveName() << ".data";


    fstream filestr;
    if(flipRate || energyMeasure)
    {
        filestr.open ((char*) oss.str().data(), fstream::in | fstream::out | fstream::app);
    }

    //initialize the points
    m_U->initPointData(m_energy_function);

    //compute total energy at first
    double Uglobal=0.0, Udata=0.0, Utopo=0.0, startUglobal=0.0;
    startUglobal = m_U->getGlobalEnergy(m_energy_function, &Udata, &Utopo);

	int flips=-1, tp, proposedUpHillMove, acceptedUpHillMove;
	m_relative_energy=0;
	unsigned long tries=m_graph->getLenEdges();


	//initialize m_move
	m_moves->initMove(m_move);
    m_moves->Tinit = m_Tinit;
    m_moves->Tend = m_Tend;


    if(m_move==SC_MOVE_V1 || m_move==SC_MOVE_V1V || m_move==SC_MOVE_V1PT)
	{
	    m_moves->m_voxelConstraint = m_voxelConstraint;
	    m_moves->m_moveFromCenter = m_moveFromCenter;
	    m_moves->m_mix->alpha1 = alphaMixture;
	    m_moves->m_mix->alpha2 = 1.0-alphaMixture;
	    m_U->initThreads();
	}

	for (unsigned long i=0 ; i<m_nbIter+MAX_EDGE_MOVED_AT_ONCE ; i++) // for the cooling schedule
	{
	    ///save current fibers if we are close to the end
	    if(i>=m_nbIter)
        {
            if(m_move==SC_MOVE_V1V)
            {
                m_moves->m_varSigma = 10*m_moves->m_varSigma;
                //save current graph configuration;
                C_bundle* extractFibers = new C_bundle(m_graph,false);
                ostringstream oss;
                oss << "fibers_iteration_" << i << "_SC_MOVE_V1V_" << m_expTAG;
                saveFiber(oss.str(), extractFibers->FIBERS);
                delete extractFibers;
            }
        }
	    proposedUpHillMove = 0;
	    acceptedUpHillMove = 0;
        for (unsigned long count=0; count < tries ; count++) // for all edges
	    {
	        ///update temperature
	        if(i < m_nbIter)
	        {
	            ///as long as we are running SA
	            T = Tcooling->schedule(m_cooling_schedule);
	        }
	        else
	        {
	            ///once SA is run, compute 10 iteration more with very low temperature: almost ICM
	            T = min(0.0000001,m_Tend/1000);
	            if(m_move==SC_MOVE_V1V)
	            {
	                m_moves->m_mix->alpha1=0.0;
	            }
	        }

            //for the cases of the energies U17 and U18, the reltive energy should initialized at each iteration so that it corresponds to the new cost function
                if(m_energy_function==SC_U18)
	        {
	            m_U->computeBetaDecrease(T, m_Tinit, m_Tend);
	        }

            m_moves->T = T;

            ///compute an iteration at the given temperature
            tp = iteration();
            flips += tp;
            if(nbFlip>0)
            {
                proposedUpHillMove++;
                if(tp>0)
                {
                    acceptedUpHillMove++;
                }
            }
	    }

	    //write measures
	    writeEnergiesAndFlipRate(&filestr, acceptedUpHillMove, proposedUpHillMove);
	}

	///really important: it frees memory allocated for moves
    m_moves->endMove();
    if(m_move==SC_MOVE_V1 || m_move==SC_MOVE_V1V || m_move==SC_MOVE_V1PT)
	{
	    m_U->deleteThreads();
	}

    ///the end!!!
    if(flipRate || energyMeasure)
    {
        filestr.close();
    }

    ///write M-file
    if(flipRate || energyMeasure)
    {
        ///measure
        Uglobal = m_U->getGlobalEnergy(m_energy_function, &Udata, &Utopo);
        writeMfile(startUglobal, Uglobal);
    }

    return (int) flips;
}
void C_SC::writeEnergiesAndFlipRate(fstream* filestr, double countflips, double modulo)
{
    //other measures: energies.
    if(energyMeasure || flipRate)
    {
        if(energyMeasure && flipRate)
        {
            double Uglobal, Udata, Utopo;
            Uglobal = m_U->getGlobalEnergy(m_energy_function, &Udata, &Utopo);
            *filestr<< T << " " << countflips/(double) modulo << " " << m_relative_energy/**-delta*/ << " " << Uglobal << " " << Udata << " " << Utopo<<endl;
        }
        else if(energyMeasure)
        {
            double Uglobal, Udata, Utopo;
            Uglobal = m_U->getGlobalEnergy(m_energy_function, &Udata, &Utopo);
            *filestr<< T << " " << m_relative_energy << " " << Uglobal << " " << Udata << " " << Utopo<<endl;
        }
        else
        {
            *filestr<< T << " " << countflips/(double) modulo << " " << m_relative_energy <<endl;
        }
    }
}
void C_SC::writeMfile(double startU, double endU)
{
    unsigned long sample = (unsigned long) (m_nbIter/5000);
    //write M-file
    if(flipRate && energyMeasure)
    {
        //
        ostringstream oss, ossMfile;
        oss << baseSaveName() << ".data";
        ossMfile << baseSaveName() << ".m";
        std::string Mfile = ossMfile.str();
        size_t found = Mfile.find(".");
        while(found!=string::npos)
        {
            //replace "." by "_"
            Mfile.replace(found, 1, "_");
            //find instance of "."
            found = Mfile.find(".");
        }
        found = Mfile.rfind("_");
        if(found!=string::npos)
        {
            Mfile.replace(found, 1, ".");
        }
        fstream filestrMfile;
        filestrMfile.open ((char*) Mfile.data(), fstream::out | fstream::app);
        filestrMfile << "clc" << endl;
        filestrMfile << "clear all" << endl;
        filestrMfile << "close all" << endl;
        filestrMfile << endl;
        filestrMfile << "%load data" << endl;
        filestrMfile << "SA = load('" << oss.str() << "');" << endl;
        if(sample>1)
        {
            filestrMfile << "SA = SA(1:" << sample << ":end,:);" << endl;
        }
        filestrMfile << endl;
        filestrMfile << "%set variables" << endl;
        filestrMfile << "Ustart = " << startU << ";" << endl;
        filestrMfile << "Uend = " << endU << ";" << endl;
        filestrMfile << "alpha = " << m_U->m_alpha << ";" << endl;
        filestrMfile << endl;
        filestrMfile << "%compute difference: abs(Uestimated-Ucomputed)"<< endl;
        filestrMfile << "deltaU = abs((Ustart+SA(:,3)) - SA(:,4));" << endl;
        filestrMfile << "%show estimated energy (blue), computed energy (red), difference (green), data energy (magenta), topological energy (yellow), starting and ending energies"<< endl;
        filestrMfile << "figure"<< endl;
        filestrMfile << "plot(1:size(SA,1), repmat(Ustart, 1, size(SA,1)))" << endl;
        filestrMfile << "hold on" << endl;
        filestrMfile << "plot(1:size(SA,1), repmat(Uend, 1, size(SA,1)))" << endl;
        filestrMfile << "plot(1:size(SA,1), Ustart+SA(:,3),'ob')" << endl;
        filestrMfile << "plot(1:size(SA,1), SA(:,4), '.r')" << endl;
        filestrMfile << "plot(1:size(SA,1), deltaU, 'g')" << endl;
        filestrMfile << "plot(1:size(SA,1), SA(:,5), 'm')" << endl;
        filestrMfile << "plot(1:size(SA,1), SA(:,6), 'y')" << endl;
        filestrMfile << "hold off" << endl;
        filestrMfile << "title('" << oss.str() << "')" << endl;
        filestrMfile << "xlabel('iteration')" << endl;
        filestrMfile << "ylabel('energy')" << endl;
        filestrMfile << endl;
        filestrMfile << "%show acceptance rate"<< endl;
        filestrMfile << "figure"<< endl;
        filestrMfile << "plot(1:size(SA,1), SA(:,2), '.b')" << endl;
        filestrMfile << "title('" << oss.str() << "')" << endl;
        filestrMfile << "xlabel('iteration')" << endl;
        filestrMfile << "ylabel('acceptance rate')" << endl;
        filestrMfile.close();
    }
    else
    {
        if(flipRate)
        {
            ostringstream oss, ossMfile;
            oss << baseSaveName() << ".data";
            ossMfile << baseSaveName() << ".m";
            std::string Mfile = ossMfile.str();
            size_t found = Mfile.find(".");
            while(found!=string::npos)
            {
                //replace "." by "_"
                Mfile.replace(found, 1, "_");
                //find instance of "."
                found = Mfile.find(".");
            }
            found = Mfile.rfind("_");
            if(found!=string::npos)
            {
                Mfile.replace(found, 1, ".");
            }
            fstream filestrMfile;
            filestrMfile.open ((char*) Mfile.data(), fstream::out | fstream::app);
            filestrMfile << "clc" << endl;
            filestrMfile << "clear all" << endl;
            filestrMfile << "close all" << endl;
            filestrMfile << endl;
            filestrMfile << "%load data" << endl;
            filestrMfile << "SA = load('" << oss.str() << "');" << endl;
            if(sample>1)
            {
                filestrMfile << "SA = SA(1:" << sample << ":end,:);" << endl;
            }
            filestrMfile << endl;
            filestrMfile << "%set variables" << endl;
            filestrMfile << "Ustart = " << startU << ";" << endl;
            filestrMfile << "Uend = " << endU << ";" << endl;
            filestrMfile << "alpha = " << m_U->m_alpha << ";" << endl;
            filestrMfile << endl;
            filestrMfile << "%show estimated (red), starting and ending energies"<< endl;
            filestrMfile << "figure"<< endl;
            filestrMfile << "plot(1:size(SA,1), repmat(Ustart, 1, size(SA,1)))" << endl;
            filestrMfile << "hold on" << endl;
            filestrMfile << "plot(1:size(SA,1), repmat(Uend, 1, size(SA,1)))" << endl;
            filestrMfile << "plot(1:size(SA,1), Ustart+SA(:,3), '.r')" << endl;
            filestrMfile << "hold off" << endl;
            filestrMfile << "title('" << oss.str() << "')" << endl;
            filestrMfile << "xlabel('iteration')" << endl;
            filestrMfile << "ylabel('energy')" << endl;
            filestrMfile << endl;
            filestrMfile << "%show acceptance rate"<< endl;
            filestrMfile << "figure"<< endl;
            filestrMfile << "plot(1:size(SA,1), SA(:,2), '.b')" << endl;
            filestrMfile << "title('" << oss.str() << "')" << endl;
            filestrMfile << "xlabel('iteration')" << endl;
            filestrMfile << "ylabel('acceptance rate')" << endl;
            filestrMfile.close();
        }
        if(energyMeasure)
        {
            ///filestr<< 0.0 << " " << m_relative_energy << " " << Uglobal << " " << Udata << " " << Utopo <<endl;
            ostringstream oss, ossMfile;
            oss << baseSaveName() << ".data";
            ossMfile << baseSaveName() << ".m";
            std::string Mfile = ossMfile.str();
            size_t found = Mfile.find(".");
            while(found!=string::npos)
            {
                //replace "." by "_"
                Mfile.replace(found, 1, "_");
                //find instance of "."
                found = Mfile.find(".");
            }
            found = Mfile.rfind("_");
            if(found!=string::npos)
            {
                Mfile.replace(found, 1, ".");
            }
            fstream filestrMfile;
            filestrMfile.open ((char*) Mfile.data(), fstream::out | fstream::app);
            filestrMfile << "clc" << endl;
            filestrMfile << "clear all" << endl;
            filestrMfile << "close all" << endl;
            filestrMfile << endl;
            filestrMfile << "%load data" << endl;
            filestrMfile << "SA = load('" << oss.str() << "');" << endl;
            if(sample>1)
            {
                filestrMfile << "SA = SA(1:" << sample << ":end,:);" << endl;
            }
            filestrMfile << endl;
            filestrMfile << "%set variables" << endl;
            filestrMfile << "Ustart = " << startU << ";" << endl;
            filestrMfile << "Uend = " << endU << ";" << endl;
            filestrMfile << "alpha = " << m_U->m_alpha << ";" << endl;
            filestrMfile << endl;
            filestrMfile << "%compute difference: abs(Uestimated-Ucomputed)"<< endl;
            filestrMfile << "deltaU = abs((Ustart+SA(:,2)) - SA(:,3));" << endl;
            filestrMfile << "%show estimated energy (blue), computed energy (red), difference (green), data energy (magenta), topological energy (yellow), starting and ending energies"<< endl;
            filestrMfile << "figure"<< endl;
            filestrMfile << "plot(1:size(SA,1), repmat(Ustart, 1, size(SA,1)))" << endl;
            filestrMfile << "hold on" << endl;
            filestrMfile << "plot(1:size(SA,1), repmat(Uend, 1, size(SA,1)))" << endl;
            filestrMfile << "plot(1:size(SA,1), Ustart+SA(:,2),'ob')" << endl;
            filestrMfile << "plot(1:size(SA,1), SA(:,3), '.r')" << endl;
            filestrMfile << "plot(1:size(SA,1), deltaU, 'g')" << endl;
            filestrMfile << "plot(1:size(SA,1), SA(:,4), 'm')" << endl;
            filestrMfile << "plot(1:size(SA,1), SA(:,5), 'y')" << endl;
            filestrMfile << "hold off" << endl;
            filestrMfile << "title('" << oss.str() << "')" << endl;
            filestrMfile << "xlabel('iteration')" << endl;
            filestrMfile << "ylabel('energy')" << endl;
            filestrMfile << endl;
            filestrMfile.close();
        }
    }

}
int C_SC::iteration(void)
{
    int flip;
    if(m_move==SC_MOVE_1_2)
    {
        flip = iteration12();
    }
    else if(m_move==SC_MOVE_1_2_4)
    {
        flip = iteration124();
    }
    else if(m_move==SC_MOVE_V1 || m_move==SC_MOVE_V1V || m_move==SC_MOVE_V1PT)
	{
	    flip = iterationP1M1();
	}
    else
    {
        //Da ist ein Problem :-(
        flip = -2;
    }

    return flip;
}

int C_SC::iteration12(void)
{
    //check/be carefull with alpha_topo
    int flip=0;

    //draw two edges at random
    if(!m_moves->move(SC_MOVE_1_2))
    {
        return -1;
    } //draw a configuration in pre-defined solution landscape
    Edge * edge1=(m_moves->aukMove->e)[0];


    if(m_moves->aukMove->nb_edge==1)
    {
        //change 10
        double partial_e1p1, partial_e1p2;
        double delta_energy = m_U->deltaEdgeU(edge1, &partial_e1p1, &partial_e1p2, m_energy_function);
        if(delta_energy>=0.0)
        {
            nbFlip = 1;
        }
        else
        {
            nbFlip = 0;
        }
        if( delta_energy < -SMALL_NUM || m_rand->doub() < exp(-delta_energy/T))
        {
            flip++;
            flip_edge(edge1,partial_e1p1,partial_e1p2);
            m_relative_energy+=delta_energy;

        }
        return flip;
    }


    if(m_moves->aukMove->nb_edge==2)
    {
        double partial_e2p1,partial_e2p2;
        Edge * edge2=(m_moves->aukMove->e)[1];
        //change 11
        double partial_e1p1, partial_e1p2;
        double delta_energy = m_U->deltaEdgeU(edge1, &partial_e1p1, &partial_e1p2, m_energy_function);

        flip_edge(edge1,partial_e1p1,partial_e1p2);
        delta_energy += m_U->deltaEdgeU(edge2, &partial_e2p1, &partial_e2p2, m_energy_function);
        if(delta_energy>=0.0)
        {
            nbFlip = 2;
        }
        else
        {
            nbFlip = 0;
        }

        if( delta_energy < -SMALL_NUM || m_rand->doub() < exp(-delta_energy/T))
        {
            flip+=2;
            flip_edge(edge2,partial_e2p1,partial_e2p2);
            m_relative_energy+=delta_energy;
        }
        else
        {
            flip_edge(edge1,partial_e1p1,partial_e1p2);
        }
        return flip;
    }
    return flip;
}

int C_SC::iteration124(void)
{
    //check/be carefull with alpha_topo
    int flip=0;

    //draw two edges at random
    if(!m_moves->move(SC_MOVE_1_2_4))
    {
        return -1;
    } //draw a configuration in pre-defined solution landscape
    Edge * edge1=(m_moves->aukMove->e)[0];


    if(m_moves->aukMove->nb_edge==1)
    {
        //change 10
        double partial_e1p1, partial_e1p2;
        double delta_energy = m_U->deltaEdgeU(edge1, &partial_e1p1, &partial_e1p2, m_energy_function);
        if(delta_energy>=0.0)
        {
            nbFlip = 1;
        }
        else
        {
            nbFlip = 0;
        }
        if( delta_energy < -SMALL_NUM || m_rand->doub() < exp(-delta_energy/T))
        {
            flip++;
            flip_edge(edge1,partial_e1p1,partial_e1p2);
            m_relative_energy+=delta_energy;

        }
        return flip;
    }


    if(m_moves->aukMove->nb_edge==2)
    {
        double partial_e2p1,partial_e2p2;
        Edge * edge2=(m_moves->aukMove->e)[1];
        //change 11
        double partial_e1p1, partial_e1p2;
        double delta_energy = m_U->deltaEdgeU(edge1, &partial_e1p1, &partial_e1p2, m_energy_function);

        flip_edge(edge1,partial_e1p1,partial_e1p2);
        delta_energy += m_U->deltaEdgeU(edge2, &partial_e2p1, &partial_e2p2, m_energy_function);
        if(delta_energy>=0.0)
        {
            nbFlip = 2;
        }
        else
        {
            nbFlip = 0;
        }

        if( delta_energy < -SMALL_NUM || m_rand->doub() < exp(-delta_energy/T))
        {
            flip+=2;
            flip_edge(edge2,partial_e2p1,partial_e2p2);
            m_relative_energy+=delta_energy;
        }
        else
        {
            flip_edge(edge1,partial_e1p1,partial_e1p2);
        }
        return flip;
    }

    if(m_moves->aukMove->nb_edge==4)
    {
        Edge * edge2=(m_moves->aukMove->e)[1];
        Edge * edge3=(m_moves->aukMove->e)[2];
        Edge * edge4=(m_moves->aukMove->e)[3];
        //change 11
        double partial_e1p1, partial_e1p2;
        double delta_energy = m_U->deltaEdgeU(edge1, &partial_e1p1, &partial_e1p2, m_energy_function);
        flip_edge(edge1,partial_e1p1,partial_e1p2);

        double partial_e2p1,partial_e2p2;
        delta_energy += m_U->deltaEdgeU(edge2, &partial_e2p1, &partial_e2p2, m_energy_function);
        flip_edge(edge2,partial_e2p1,partial_e2p2);

        double partial_e3p1,partial_e3p2;
        delta_energy += m_U->deltaEdgeU(edge3, &partial_e3p1, &partial_e3p2, m_energy_function);
        flip_edge(edge3,partial_e3p1,partial_e3p2);

        double partial_e4p1,partial_e4p2;
        delta_energy += m_U->deltaEdgeU(edge4, &partial_e4p1, &partial_e4p2, m_energy_function);

        if(delta_energy>=0.0)
        {
            nbFlip = 4;
        }
        else
        {
            nbFlip = 0;
        }

        if( delta_energy < -SMALL_NUM || m_rand->doub() < exp(-delta_energy/T))
        {
            flip+=4;
            flip_edge(edge4,partial_e4p1,partial_e4p2);
            m_relative_energy+=delta_energy;
        }
        else
        {
            flip_edge(edge3,partial_e3p1,partial_e3p2);
            flip_edge(edge2,partial_e2p1,partial_e2p2);
            flip_edge(edge1,partial_e1p1,partial_e1p2);
        }
        return flip;
    }
    return flip;
}


int C_SC::iterationP1M1(void)
{
    //check/be carefull with alpha_topo
    int flip=0;

    //draw two edges at random
    m_moves->move(m_move/*SC_MOVE_V1*/); //draw a configuration in pre-defined solution landscape


    if(m_moves->aukMove->nb_edge==1)
    {
        //change 10
        Edge * edge1=(m_moves->aukMove->e)[0];
        double partial_e1p1, partial_e1p2;
        double delta_energy = m_U->deltaEdgeU(edge1, &partial_e1p1, &partial_e1p2, m_energy_function);

        double dU = delta_energy;

        if(delta_energy>=0.0)
        {
            nbFlip = 1;
        }
        else
        {
            nbFlip = 0;
        }

        if( dU < -SMALL_NUM || m_rand->doub() < exp(-dU/T))
        {
            flip++;
            flip_edge(edge1,partial_e1p1,partial_e1p2);
            m_relative_energy+=delta_energy;
        }
        return flip;
    }


    if(m_moves->aukMove->nb_edge==0)
    {
        //int flip=0;
        double delta_energy;//, partial_sum1,partial_sum2;

        //draw a random site in graph
        Point* pp = m_graph->getPoints();
        Point* V = &(pp[m_moves->aukMove->idxPoint]);
        Point* newV = m_moves->aukMove->p;
        vector<long> idxPointVector;
        vector< vector<long> > idxPairVector;
        vector< vector<double> > dataPairVector;
        //compute energy differnce
        delta_energy = m_U->deltaVertexMoveU(V,newV, m_energy_function, &idxPointVector, &idxPairVector, &dataPairVector);


        if(delta_energy>0.0)
        {
            nbFlip = 1;
        }
        else
        {
            nbFlip = 0;
        }

        //check edge has to be flipped
        if( delta_energy < -SMALL_NUM || m_rand->doub() < exp(-delta_energy/T))
        {
            //flip count flip and update local variables
            flip++;
            m_relative_energy+=delta_energy;
            long idxP = -1;

            //modify point coordinates
            V->x = newV->x;
            V->y = newV->y;
            V->z = newV->z;
            //modify vectors
            for(int k=0 ; k<V->n_edges ; k++)
            {
                Point*  q = V->edges[k].other_point;
                //way forward
                V->edges[k].dx = q->x - V->x;
                V->edges[k].dy = q->y - V->y;
                V->edges[k].dz = q->z - V->z;

                //way back
                for(int kq=0 ; kq<q->n_edges ; kq++)
                {
                    if(q->edges[kq].other_point==V)
                    {
                        q->edges[kq].dx = V->x - q->x;
                        q->edges[kq].dy = V->y - q->y;
                        q->edges[kq].dz = V->z - q->z;
                    }
                }
            }
            //update graph: error_sum (=-sum(data)), and data
            for(long i=0 ; i<(long) idxPointVector.size() ; i++)
            {
                //get point index
                idxP = idxPointVector.at(i);
                //get the edge binding V and pp[idxP] and modify the vector coordinates
                Edge* e2 = NULL;
                //long idxE=-1;
                for(int k=0 ; k<pp[idxP].n_edges ; k++)
                {
                    if(pp[idxP].edges[k].other_point==V)
                    {
                        e2 = pp[idxP].edges[k].edge;
                        k = pp[idxP].n_edges;
                    }
                }
                //remove contribution of current data
                if(e2->connected)
                {
                    for(int k=0 ; k<pp[idxP].n_edges ; k++)
                    {
                        if(pp[idxP].edges[k].connected)
                        {
                            if(m_energy_function==U15 || m_energy_function==U18)
                            {
                                pp[idxP].error_sum -= m_U->U15data(e2, pp[idxP].edges[k].edge);
                            }
                            else
                            {
                                pp[idxP].error_sum -= m_U->U16data(e2, pp[idxP].edges[k].edge);
                            }
                        }
                    }
                }

                //update energy info
                for(long j=0 ; j<(long) idxPairVector.at(i).size() ; j++)
                {
                    //update data
                    pp[idxP].data[idxPairVector.at(i).at(j)] = dataPairVector.at(i).at(j);
                }

                //add contribution of current data
                if(e2->connected)
                {
                    for(int k=0 ; k<pp[idxP].n_edges ; k++)
                    {
                        if(pp[idxP].edges[k].connected)
                        {
                            if(m_energy_function==U15 || m_energy_function==U18)
                            {
                                pp[idxP].error_sum += m_U->U15data(e2, pp[idxP].edges[k].edge);
                            }
                            else
                            {
                                pp[idxP].error_sum += m_U->U16data(e2, pp[idxP].edges[k].edge);
                            }
                        }
                    }
                }
            }
        }
        return flip;
    }
    return flip;
}


string C_SC::baseSaveName(void)
{
    ostringstream oss;

    oss << baseName << "_" << m_expTAG ;
    return oss.str();
}