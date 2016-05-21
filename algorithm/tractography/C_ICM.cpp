#include <C_ICM.h>


C_ICM::C_ICM(C_graph* graph,\
             C_energy* U /**pointer to the energy related to the C_graph graph: assuming it's already set up*/,\
                 unsigned short energy_func /**kind of energy function: see settings_def.h*/,\
                  unsigned short move_kind /**kind of moves: see settings_def.h*/,\
                     bool rand_sweep):C_minimization(graph)
{
    //ctor
    m_energy_function = energy_func;
    m_move = move_kind;
    m_rand_sweep = rand_sweep;
    m_U = U;
    m_fileName = "new_file.cformat";
}

C_ICM::~C_ICM()
{
}


/** working methods */
int C_ICM::run(void)
{
    fstream filestr;
    string dataName = (m_fileName+m_TAG)+".data";
    filestr.open ((char*) dataName.data(), fstream::out);// | fstream::app);

    int flips;

    m_U->initPointData(m_energy_function);

    //compute total energy at first
    double Edata=0.0, Etopo=0.0;
    double totU = m_U->getGlobalEnergy(m_energy_function, &Edata, &Etopo);
    globU = totU;
    filestr<< flips << ", " << globU << ", " << totU << ", " << Edata << ", " << Etopo << endl;

    //compute iterations
    unsigned long regFlip[3] = {10, 10, 10};
    unsigned long iter = 0;
    do{
        flips = iteration();
        cout << "iteration " << iter << " #flips " << flips << endl;
        iter++;
        totU = m_U->getGlobalEnergy(m_energy_function, &Edata, &Etopo);
        filestr<< flips << ", " << globU << ", " << totU << ", " << Edata << ", " << Etopo << endl;
        //decalage
        regFlip[2] = regFlip[1];
        regFlip[1] = regFlip[0];
        regFlip[0] = (unsigned int) flips;
        flips = (int) ((regFlip[0] + regFlip[1] + regFlip[2])/3);
        cout << flips << endl;
    }while(flips>1);

    totU = m_U->getGlobalEnergy(m_energy_function, &Edata, &Etopo);
    filestr<< flips << ", " << globU << ", " << totU << ", " << Edata << ", " << Etopo << endl;
    filestr.close();

    return (int) 1;
}

int C_ICM::iteration(void)
{
    //according to m_move and m_marc_option
    if(m_move==ICM_MOVE_1)
    {
        return icm_iteration();
    } else //if(m_move==ICM_MOVE_2)
    {
        return icm_iteration_v2();
    }
}




/**Carole's methods*/
int C_ICM::icm_iteration(void)
{
    Edge * edges = m_graph->getEdges();
    int n_edges = m_graph->getLenEdges();
    int count_flips=0;
    double partial_sum1, partial_sum2, delta_energy;
    unsigned long * index_array=0;
    if (m_rand_sweep) index_array= sweep_idx_edges();

    for(int i=0;i<n_edges;i++)
    {
        struct Edge* current_edge;
        if (m_rand_sweep)
                current_edge=&edges[index_array[i]];
        else
                current_edge=&edges[i];

        delta_energy = m_U->deltaEdgeU(current_edge,&partial_sum1, &partial_sum2, m_energy_function);
        if( delta_energy < -SMALL_NUM)
        {
            count_flips++;
            flip_edge(current_edge,partial_sum1,partial_sum2);
            m_relative_energy+=delta_energy;
            globU += delta_energy;
        }
    }

    if(m_rand_sweep) delete index_array;

    return count_flips;
}

int C_ICM::icm_iteration_v2(void)
{
    Point * points = m_graph->getPoints();
    int count_flips=0;
    double partial_e1p1,partial_e1p2;
    double partial_e2p1,partial_e2p2;
    double delta;
    Edge* edge1;
    Edge* edge2;
    Point* current_point;
    unsigned long * index_array=0;
    if (m_rand_sweep) index_array=sweep_idx_points();

    int n_points = m_graph->getLenPoints();
    for(int i=0;i<n_points;i++)
    {
        if (m_rand_sweep)
                current_point=&points[index_array[i]];
        else
                current_point=&points[i];
        for (int e1=0;e1<current_point->n_edges;e1++) //For all combinations of adjacent edges
        {
            for(int e2=0;e2<e1;e2++)
            {
                edge1 = current_point->edges[e1].edge;
                edge2 = current_point->edges[e2].edge;
                delta = m_U->deltaEdgeU(edge1,&partial_e1p1, &partial_e1p2, m_energy_function);
                flip_edge(edge1,partial_e1p1,partial_e1p2);
                delta += m_U->deltaEdgeU(edge2,&partial_e2p1, &partial_e2p2, m_energy_function);
                if( delta< -SMALL_NUM)
                {
                        count_flips++;
                        flip_edge(edge2,partial_e2p1,partial_e2p2);
                        m_relative_energy+=delta;
                        globU += delta;
                }
                else
                        flip_edge(edge1,partial_e1p1,partial_e1p2);
            }
        }
    }

    if (m_rand_sweep) delete index_array;

    return count_flips;
}
