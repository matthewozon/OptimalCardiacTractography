#include <C_temperature.h>

C_temperature::C_temperature(C_graph* graph, C_energy* U, unsigned short energy_kind)
{
    //ctor
    m_energy_kind = energy_kind;
    m_U = U;
    m_graph = graph;

    idum = (unsigned long long int) time(NULL);
    while(idum==THE_NUMBER_WHICH_MUST_NOT_BE_NAMED)
    {
        idum = (unsigned long long int) time(NULL);
    }
    m_rand = new C_toolbox_rand2(idum);

    Tinit = 999999999999999999999.0;
    Tend = 0.00000000000000000001;
    m_varSigma = 0.2;
    m_alphaMixture = 0.975;
    m_voxelConstraint = true;
    m_moveFromCenter = true;
}
C_temperature::~C_temperature()
{
    //dtor
    delete m_rand;
}

double C_temperature::rootFinder(double* deltaU /*an array of positive value*/, unsigned long M /*number of elements in deltaU*/, double tau /*acceptance rate*/, double epsilon /*precision on beta*/)
{
    double beta = 0.0;
    double f = 1.0;
    long int i=0;
    while(f>=epsilon)
    {
        beta -= (fxtau(beta, deltaU, (double) M, tau)/fxtau_deriv(beta, deltaU, (double) M));
        f = fxtau(beta, deltaU, (double) M, tau);
        i++;
        if(i>100000)
        {
            throw ENDED_NOT_SUCCESS;
        }
    }
    return 1.0/beta;
}

double C_temperature::fxtau(double beta, double* deltaU, double M, double tau)
{
    double res = -M*tau;
    for(unsigned long k=0; k<(unsigned long) M; k++)
        res += exp(-beta*deltaU[k]);
    return res;
}
double C_temperature::fxtau_deriv(double beta, double* deltaU, double M)
{
    double res = 0.0;
    for(unsigned long k=0; k<(unsigned long) M; k++)
        res -= deltaU[k]*exp(-beta*deltaU[k]);

    return res;
}


bool C_temperature::runTemperature(unsigned long kind_move_init, unsigned long kind_move_end, double tauInit, double tauEnd, unsigned long M)
{
    unsigned long N = M;
    m_graph->randInit(0.2);
    startingMove = true;
    if(m_energy_kind==U18)
    {
        m_U->m_beta0 = m_U->m_beta;
        m_U->computeBetaDecrease(1.0,1.0,0.5);
    }
    double* deltaU = uphill_energy(kind_move_init, true, &N); //ok
    if(deltaU!=NULL)
    {
        if(N==0)
        {
            delete deltaU;
            return false;
        }
        //select delta energy by generating a markov chain according to kernel communication
        Tinit = rootFinder(deltaU, N, tauInit,0.00000000001);
        delete deltaU;

        //select delta energy by generating a markov chain according to kernel communication
        N = M;
        m_graph->randInit(0.15);
        startingMove = false;
        if(m_energy_kind==U18)
        {
            m_U->m_beta0 = m_U->m_beta;
            m_U->computeBetaDecrease(0.5,1.0,0.5);
        }
        deltaU = uphill_energy(kind_move_end, true/**false*/, &N);
        if(m_energy_kind==U18)
        {
            m_U->m_beta = m_U->m_beta0;
        }
        if(deltaU!=NULL)
        {
            if(N==0)
            {
                delete deltaU;
                return false;
            }
            Tend = rootFinder(deltaU, N, tauEnd,0.00000000001);

            delete deltaU;
        }else
        {
//            cout << "did not manage to find second deltaU" << endl;
            return false;
        }

    }else
    {
//        cout << "did not manage to find first deltaU" << endl;
        return false;
    }
    return true;
}




double* C_temperature::uphill_energy(unsigned long kind_move, bool init, unsigned long* N)
{
    //usual init
    unsigned long M = *N;
    double* deltaU = new double[M];
    if(deltaU==NULL) return deltaU;

    //init a move
    C_move* temperature_move = new C_move(m_graph);
    temperature_move->initMove(kind_move);

    //generate Markov chain according to communication mechanism till having M uphill moves
    unsigned long k = 0;
    double delta;

    if((kind_move==MOVE_T1+MOVE_T2) || (kind_move==MOVE_V) || (kind_move==MOVE_V+MOVE_T1) || (kind_move==MOVE_V+MOVE_T2) || (kind_move==SC_MOVE_V1) || (kind_move==SC_MOVE_V1PT) || (kind_move==SC_MOVE_V1V) )
    {
//        cout << "HERE? " << temperature_move->m_mix << endl;
        temperature_move->m_mix->alpha1 = m_alphaMixture;
        temperature_move->m_mix->alpha2 = 1.0-m_alphaMixture;
        temperature_move->m_varSigma = m_varSigma;
        temperature_move->m_voxelConstraint = m_voxelConstraint;
        temperature_move->m_moveFromCenter = m_moveFromCenter;
        ///force temperatures of move so that it knows which move to generate
        if(startingMove)
        {
            temperature_move->Tinit = 2.0;
            temperature_move->Tend = 1.0;
            temperature_move->T = 2.0;
        }
        else
        {
            temperature_move->Tinit = 2.0;
            temperature_move->Tend = 1.0;
            temperature_move->T = 1.0;
        }
    }

    //if we are playing with multithreading
    if(kind_move==MOVE_V || (kind_move==MOVE_V+MOVE_T1) || (kind_move==MOVE_V+MOVE_T2) || (kind_move==SC_MOVE_V1) || (kind_move==SC_MOVE_V1PT) || (kind_move==SC_MOVE_V1V) )
    {
        m_U->initThreads();
    }

    unsigned long i = 0;
    while(k<M && i<500*M)
    {
        i++;
        //draw a new configuration

        //compute deltaUTmp and update graph
        //this part could be move to iteration of SA by setting
        //T as large as possible, and getting m_relative_energy before and after each move,
        //delta = m_relative_energy_after-m_relative_energy_before, but I have to load a new graph...
        //let's see if it's really usefull to have an ICM object there, otherwise, delete ICM and use SA
        if(kind_move==MOVE_T1)
        {
            //cout << "before move generation" << endl;
            temperature_move->move(kind_move);
            //cout << "before delta generation" << endl;
            delta = uphill_energyM1(temperature_move->aukMove->e, init);
            //cout << "after delta generation delta = " << delta << endl;
        }
        else if(kind_move==MOVE_T2)
        {
            temperature_move->move(kind_move);
            delta = uphill_energyM2(temperature_move->aukMove->e, init);
        }
        else if(kind_move==MOVE_T4)
        {
            temperature_move->move(kind_move);
            delta = uphill_energyM4(temperature_move->aukMove->e, init);
        }
        else if(kind_move==MOVE_T1+MOVE_T2)
        {
            temperature_move->move(MOVE_T1+MOVE_T2);
            delta = uphill_energyM1M2(temperature_move->aukMove, init);
        }
        else if(kind_move==MOVE_V)
        {
            temperature_move->move(MOVE_V);
            delta = uphill_energyP1(temperature_move->aukMove->idxPoint, temperature_move->aukMove->p);
        }
        else if(kind_move==MOVE_V+MOVE_T1)
        {
            temperature_move->move(MOVE_V+MOVE_T1);
            if(temperature_move->aukMove->nb_edge==0)
            {
                delta = uphill_energyP1(temperature_move->aukMove->idxPoint, temperature_move->aukMove->p);
            }
            else
            {
                delta = uphill_energyM1(temperature_move->aukMove->e, init);
            }
        }
        else if( (kind_move==SC_MOVE_V1) || (kind_move==SC_MOVE_V1V) || (kind_move==SC_MOVE_V1PT) )
        {
            temperature_move->move(kind_move);
            if(temperature_move->aukMove->nb_edge==0)
            {
                delta = uphill_energyP1(temperature_move->aukMove->idxPoint, temperature_move->aukMove->p);
            }
            else
            {
                delta = uphill_energyM1(temperature_move->aukMove->e, init);
            }
        }
        else //if(kind_move==MOVE_V+MOVE_T2)
        {
            temperature_move->move(MOVE_V+MOVE_T2);
            if(temperature_move->aukMove->nb_edge==0)
            {
                delta = uphill_energyP1(temperature_move->aukMove->idxPoint, temperature_move->aukMove->p);
            }
            else
            {
                delta = uphill_energyM2(temperature_move->aukMove->e, init);
            }
        }

        if(delta>0)
        {
            deltaU[k] = delta;
            k++;
        }
    }
    *N = k;
    temperature_move->endMove();
    delete temperature_move;

    //if we are playing with multithreading
    if(kind_move==MOVE_V || (kind_move==MOVE_V+MOVE_T1) || (kind_move==MOVE_V+MOVE_T2) || (kind_move==SC_MOVE_V1) || (kind_move==SC_MOVE_V1V) || (kind_move==SC_MOVE_V1PT))
    {
        m_U->deleteThreads();
    }

    return deltaU;
}


void C_temperature::flip_edge(Edge * edge,double psum_p1,double psum_p2)
{

	int was_connected=edge->connected;
	edge->connected = !was_connected;
	edge->this_in_p1->connected = !was_connected;
	edge->this_in_p2->connected = !was_connected;

    if (was_connected)
    {
        edge->p1->degree--;
        edge->p2->degree--;
        edge->p1->error_sum-=psum_p1;
        edge->p2->error_sum-=psum_p2;
    }
    else
    {
        edge->p1->degree++;
        edge->p2->degree++;
        edge->p1->error_sum+=psum_p1;
        edge->p2->error_sum+=psum_p2;
    }
}

double C_temperature::uphill_energyM1(Edge** edges, bool init)
{
    double partial_sum1 = 0, partial_sum2 = 0;
    double delta = m_U->deltaEdgeU(edges[0], &partial_sum1, &partial_sum2, m_energy_kind);
    if(init){
        flip_edge(edges[0], partial_sum1, partial_sum2);
    }
    return delta;
}
double C_temperature::uphill_energyP1(long idxPoint, Point* newV)
{
    vector<long> idxPointVector;
    vector< vector<long> > idxPairVector;
    vector< vector<double> > dataPairVector;
    if(newV->x<0.0) newV->x=0.0;
    if(newV->y<0.0) newV->y=0.0;
    if(newV->z<0.0) newV->z=0.0;
    double dU = m_U->deltaVertexMoveU(&(m_graph->getPoints()[idxPoint]),newV, m_energy_kind, &idxPointVector, &idxPairVector, &dataPairVector);
    idxPointVector.clear();
    idxPairVector.clear();
    dataPairVector.clear();
    return dU;
}
double C_temperature::uphill_energyM2(Edge** edges, bool init)
{
    double delta = 0;
    double partial_sum1, partial_sum2;
    Edge * edge1=edges[0];
    Edge * edge2=edges[1];

    delta = m_U->deltaEdgeU(edge1, &partial_sum1, &partial_sum2, m_energy_kind);
    flip_edge(edge1,partial_sum1,partial_sum2);
    double tmp_partial_1=partial_sum1, tmp_partial_2=partial_sum2;
    delta += m_U->deltaEdgeU(edge2, &partial_sum1, &partial_sum2, m_energy_kind);
    if(init){
        //update graph
        flip_edge(edge2,partial_sum1,partial_sum2);
    }else{
        //go back to older configuration
        flip_edge(edge1,tmp_partial_1,tmp_partial_2);
    }

    return delta;

}

double C_temperature::uphill_energyM4(Edge** edges, bool init)
{
    Edge* edge1 = edges[0];
    Edge* edge2 = edges[1];
    Edge* edge3 = edges[2];
    Edge* edge4 = edges[3];

    //change 11
    double partial_e1p1, partial_e1p2;
    double delta = m_U->deltaEdgeU(edge1, &partial_e1p1, &partial_e1p2, m_energy_kind);
    flip_edge(edge1,partial_e1p1,partial_e1p2);

    double partial_e2p1,partial_e2p2;
    delta += m_U->deltaEdgeU(edge2, &partial_e2p1, &partial_e2p2, m_energy_kind);
    flip_edge(edge2,partial_e2p1,partial_e2p2);

    double partial_e3p1,partial_e3p2;
    delta += m_U->deltaEdgeU(edge3, &partial_e3p1, &partial_e3p2, m_energy_kind);
    flip_edge(edge3,partial_e3p1,partial_e3p2);

    double partial_e4p1,partial_e4p2;
    delta += m_U->deltaEdgeU(edge4, &partial_e4p1, &partial_e4p2, m_energy_kind);

    if(init)
    {
        flip_edge(edge4,partial_e4p1,partial_e4p2);
    }
    else
    {
        flip_edge(edge3,partial_e3p1,partial_e3p2);
        flip_edge(edge2,partial_e2p1,partial_e2p2);
        flip_edge(edge1,partial_e1p1,partial_e1p2);
    }
    return delta;
}


double C_temperature::uphill_energyM1M2(MovingEdge* aukMove, bool init)
{
    double delta = 0;
    double partial_sum1, partial_sum2;

    Edge * edge1=aukMove->e[0];


    if(aukMove->nb_edge==2)
    {
        Edge * edge2=aukMove->e[1];
        delta = m_U->deltaEdgeU(edge1, &partial_sum1, &partial_sum2, m_energy_kind);
        flip_edge(edge1,partial_sum1,partial_sum2);
        double tmp_partial_1=partial_sum1, tmp_partial_2=partial_sum2;
        delta += m_U->deltaEdgeU(edge2, &partial_sum1, &partial_sum2, m_energy_kind);
        if(init){
            //update graph
            flip_edge(edge2,partial_sum1,partial_sum2);
        }else{
            //go back to older configuration
            flip_edge(edge1,tmp_partial_1,tmp_partial_2);
        }
    }
    else
    {
        delta = m_U->deltaEdgeU(edge1, &partial_sum1, &partial_sum2, m_energy_kind);
        if(init)
        {
            //update graph
            flip_edge(edge1,partial_sum1,partial_sum2);
        }
    }
    return delta;
}

