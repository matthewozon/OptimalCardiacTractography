#include "C_thread_deltaU_move_vertex.h"

C_thread_deltaU_move_vertex::C_thread_deltaU_move_vertex(C_energy* U):C_thread()
{
    //ctor
    m_U = U;
}

C_thread_deltaU_move_vertex::~C_thread_deltaU_move_vertex()
{
    //dtor
}

void C_thread_deltaU_move_vertex::resetAttributes()
{
    u = NULL;
    v = NULL;
    v_new = NULL;
    idxPair.clear();
    dataPair.clear();
    dU_data = 0.0;
//    dU_topo = 0.0;
}


void C_thread_deltaU_move_vertex::execute()
{
    ///assuming that all attributes have been set

    ///find the edge, and its index, binding u and the current v
    Edge* e2 = NULL;
    long idxEdge2 = -1;
    for(long i=0 ; i<u->n_edges ; i++) ///go through all edges of u
    {
        if(u->edges[i].other_point==v)
        {
            e2 = u->edges[i].edge;
            idxEdge2 = i;
            i = u->n_edges;
        }
    }

    ///calculate new energy for pairs including the point u and v, with the new position of v(v_new)
    dU_data = 0.0;
//    dU_topo = 0.0;
    vector<long> idxEdgeE1;
    for(long i=0 ; i<u->n_edges ; i++) ///go through all edges of u
    {
        if(i!=idxEdge2)
        {
            ///store pair idx
            idxPair.push_back(m_U->m_graph->getIndexOfPair(u, i, idxEdge2));
            idxEdgeE1.push_back(i);

            ///calculate and store data fidelity with U15dataComputation(v_new, other_point)
            if(m_kind==U15 || m_kind==U18)
            {
                dataPair.push_back(m_U->U15dataComputation(u->edges[i].other_point,v_new)); //should be the time consuming part
            }
            else
            {
                dataPair.push_back(m_U->U16dataComputation(u->edges[i].other_point,v_new)); //should be the time consuming part
            }
        }
    }

    ///if at least two edges are connected, calculate energy variation
    if(u->degree>1)
    {
        ///sumup with normalization
        if(e2->connected)
        {
            for(unsigned int i=0 ; i<dataPair.size() ; i++)
            {
                //if the other edge is connected, add it up
                if(u->edges[idxEdgeE1.at(i)].connected)
                {
                    dU_data -= (dataPair.at(i)-u->data[idxPair.at(i)]);
                }
            }
            dU_data /= (u->degree*(u->degree-1))/2;
        }
    }
    return;
}
