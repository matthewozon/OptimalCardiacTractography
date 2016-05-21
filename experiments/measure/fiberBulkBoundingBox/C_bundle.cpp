#include "C_bundle.h"

C_bundle::C_bundle(C_graph* graph, bool cutAcuteEdges)
{
    //ctor
    m_graph = new C_graph(*graph);
    m_point = m_graph->getPoints();
    nb_points = m_graph->getLenPoints();
    m_edge = m_graph->getEdges();
    nb_edges = m_graph->getLenEdges();
    FIBERS.clear();
    try{
        if(!extractFibers(cutAcuteEdges)) throw 4;
    }
    catch(...)
    {
        cout << "An error had been thrown while extracting fibers: no fiber had been extracted" << endl;
        //throw 1;
    }

    if(false) ///wait before setting it to true: wait fora good functional that does not allow folding fibers
    {
        try
        {
            fiberProc();// throw 4;
        }
        catch(...)
        {
            cout << "An error had been thrown while extracting fibers" << endl;
        }
    }


}

C_bundle::~C_bundle()
{
    //dtor
    FIBERS.clear();
    if(m_graph!=NULL)
    {
        delete m_graph;//->~C_graph();
    }
}





bool C_bundle::extractFibers(bool cutAcuteEdges)
{

    ///graph corrections
    ///get rid of points where d(v)>2
    //clean graph: for each point which is connected to more than 2 edge, set all edges to 0
    try{
        if(!cleanGraph()) throw 2; //ok
    }catch(bad_alloc& ba)
    {
        cout << "In cleanGraph: " << ba.what() << endl;
        throw 2;
    }
    catch(...)
    {
        cout << "cleanGraph generated an exeption: fibers won't be extracted" << endl;
        throw 2;
    }
    //check if there are still something in graph
    if(m_graph->getNbConnectedEdge()==0)
    {
        cout << "no edge left after cleanGraph" << endl;
        FIBERS.clear();
        return false;
    }

    ///cute acute edges: delete points where d(v)=2 & angle<pi/2
    //it could be interesting to shrink the space :-): means one can store all connected edges and point in arrays to find it faster after
    vector< vector<unsigned long> > fibers;
    vector<unsigned long> tempIdx;

    if(cutAcuteEdges)
    {
        //cut forbiden patern
        try{
            if(!cutAcuteEdge()) throw 1;
        }catch(bad_alloc& ba)
        {
            cout << "In cutAcuteEdge: " << ba.what() << endl;
            throw 1;
        }
        catch(...)
        {
            cout << "cutAcuteEdge generated an exeption: fibers won't be extracted" << endl;
            throw 1;
        }

        //check if there are still something in graph
        if(m_graph->getNbConnectedEdge()==0)
        {
            cout << "no edge left after cutAcuteEdge" << endl;
            FIBERS.clear();
            return false;
        }
    }

    ///clean the mess: single edges
    //delete single edge
    try
    {
        if(!deleteSingle()) throw 7; ///ok
    }
    catch(bad_alloc& ba)
    {
        cout << "In deleteSingle: " << ba.what() << endl;
        throw 7;
    }
    catch(...)
    {
        cout << "deleteSingle generated an exeption: fibers won't be extracted" << endl;
        throw 7;
    }
    //check if there are still something in graph
    if(m_graph->getNbConnectedEdge()==0)
    {
        cout << "no edge left after deleteSingle" << endl;
        FIBERS.clear();
        return false;
    }


    ///fiber extraction


    ///fiber corection: get rid of folding fibers

    ///IT SHALL RETURN










    //extract not-loop-shaped fibers
    try{
        if(!extractLinearFibers())
        {
            throw 3;
        }
    }catch(bad_alloc& ba)
    {
        cout << "In extractLinearFibers: " << ba.what() << endl;
        throw 3;
    }
    catch(...)
    {
        cout << "extractLinearFibers generated an exeption: fibers won't be extracted" << endl;
        throw 3;
    }

    //check if there are still something in graph
    if(m_graph->getNbConnectedEdge()==0)
    {
        return true;
    }

    //get number of connected edges
    unsigned long nb_connected_edge = 0; //must be the same as nb_vertices :-) because there are only loops
    for(unsigned long k=0 ; k<nb_edges ; k++)
    {
        //check if edge is connected
        if(m_edge[k].connected)
        {
            nb_connected_edge++;
        }
    }




    try{
        if(!extractLoopFibers()) throw 5;
        nb_connected_edge = 0; //must be the same as nb_vertices :-) because there are only loops
        for(unsigned long k=0 ; k<nb_edges ; k++)
        {
            //check if edge is connected
            if(m_edge[k].connected)
            {
                nb_connected_edge++;
            }
        }
    }catch(bad_alloc& ba)
    {
        cout << "In large loop fucking bad_alloc: " << ba.what() << endl;
        throw 4;
    }

    return true;
}


bool C_bundle::extractLinearFibers()
{
    //get all starting point
    vector<unsigned long> idx_start;

    for(unsigned long i=0 ; i<nb_points ; i++)
    {
        if(m_point[i].degree==1)
        {
            idx_start.push_back(i);
        }
    }

    //store it in usual array
    unsigned long nb_idx_start = idx_start.size();
    unsigned long** idx_S = new unsigned long*[nb_idx_start];
    for(unsigned long i=0 ; i<idx_start.size() ; i++)
    {
        idx_S[i] = new unsigned long[2];
        idx_S[i][0] = idx_start.at(i);
        idx_S[i][1] = 1;
    }
    idx_start.erase(idx_start.begin(), idx_start.end());

    //while there is at least one starting point
    unsigned long start_idx , end_idx, enable_idx;
    enable_idx = nb_idx_start;
    start_idx = nb_points;
    while(enable_idx>1)
    {
        //draw the first starting point and disable it in idx_start array
        for(unsigned long i=0 ; i<nb_idx_start ; i++)
        {
            if(idx_S[i][1]==1)
            {
                idx_S[i][1] = 0;
                start_idx = idx_S[i][0];
                i = nb_idx_start;
            }
        }
        enable_idx--;

        //track the fiber
        extractLinearFiber(start_idx);


        //hope it's OK
        unsigned long nb_fibers = FIBERS.size();
        unsigned long nb_elt_last_fiber = FIBERS.at(nb_fibers-1).size();
        end_idx = (FIBERS.at(nb_fibers-1)).at(nb_elt_last_fiber-1);

        //disable ending index in idx_start array
        for(unsigned long i=0 ; i<nb_idx_start ; i++)
        {
            if(idx_S[i][0]==end_idx)
            {
                if(idx_S[i][1]==1)
                {
                    idx_S[i][1] = 0;
                    i = nb_idx_start;
                }
                else
                {
                    enable_idx = 0;
                }

            }
        }


        //check how many starting point are enabled
        enable_idx--;
    }

    for(unsigned long i=0 ; i<idx_start.size() ; i++)
    {
        delete idx_S[i];
    }
    delete idx_S;


    return true;
}


bool C_bundle::extractLinearFiber(unsigned long idx_start)
{
    vector<unsigned long> tempIdx;
    tempIdx.push_back(idx_start);
    unsigned long new_idx;

    //init
    Edge* e = NULL; //really important to set it to NULL for the first point
    Point* p = &(m_point[idx_start]);

    for(unsigned i=0 ; i<(unsigned long) p->n_edges ; i++)
    {
        if((p->edges[i]).edge->connected)
        {
            e = (p->edges[i]).edge;
            i = p->n_edges;
        }
    }

    disconnectEdge(e);

    if(p!=e->p2 && p!=e->p1) //means that p != p1 and p != p2 => does not belong to the given edge!!!!! should never happen
    {
        p = NULL;
    }
    else
    {
        if(p==e->p1)
        {
            p = e->p2;
        }
        else
        {
            p = e->p1;
        }
    }

    new_idx = p - m_point;
    tempIdx.push_back(new_idx);

    //iter
    bool cond = true;
    unsigned long k=0;
    unsigned long NN = m_graph->getLenPoints();
    while(cond && k<NN+1){
        k++;
        p = &(m_point[new_idx]);

        switch(p->degree)
        {
            case 0:
            {
                //means we reach the end
                cond = false;
                break;
            }
            case 1:
            {
                //usual case in fiber
                for(unsigned i=0 ; i<(unsigned long) p->n_edges ; i++)
                {
                    if((p->edges[i]).edge->connected)
                    {
                        //didItCrossIfWay = true;
                        e = (p->edges[i]).edge;
                        i = p->n_edges;
                    }
                }

                disconnectEdge(e);

                if(p!=e->p2 && p!=e->p1) //means that p != p1 and p != p2 => does not belong to the given edge!!!!! should never happen
                {
                    p = NULL;
                }
                else
                {
                    if(p==e->p1)
                    {
                        p = e->p2;
                    }
                    else
                    {
                        p = e->p1;
                    }
                }


                new_idx = p - m_point;

                tempIdx.push_back(new_idx);

                break;
            }
            default:
            {
                //means there is a pb!!!! or this is the end of a bundle
                cond = false;
                break;
            }
        }

    }

    FIBERS.push_back(tempIdx);
//    if(k==NN)
//    {
//        return false;
//    }
    return true;
}


bool C_bundle::extractLoopFibers(void)
{
    //count how many point with two connected edge are there
    unsigned long count = 0;
    for(unsigned long i=0 ; i<nb_points ; i++)
    {
        if(m_point[i].degree==2)
        {
            count++;
        }
    }

    unsigned long idx;
    while(count>2)
    {
        //find first point with two connected edges
        idx = 0;
        for(unsigned long i=0 ; i<nb_points ; i++)
        {
            if(m_point[i].degree==2)
            {
                idx = i;
                i = nb_points;
            }
            if(m_point[i].degree!=2 && i==nb_points-1 && idx!=0)
            {
                idx=nb_points;
            }
        }

        if(idx<nb_points)
        {
            //find the first edge
            for(unsigned long j=0 ; j<(unsigned long) m_point[idx].n_edges ; j++)
            {
                if(m_point[idx].edges[j].edge->connected)
                {
                    //disconnect it
                    disconnectEdge(m_point[idx].edges[j].edge);
                    j = m_point[idx].n_edges;
                }
            }

            //run loop extraction
            extractLoopFiber(idx);

            //update count
            count = 0;
            for(unsigned long i=0 ; i<nb_points ; i++)
            {
                if(m_point[i].degree==2)
                {
                    count++;
                }
            }
        }
        else
        {
            count = 0;
        }

    }
    return true;
}

bool C_bundle::extractLoopFiber(unsigned long idx)
{
    return extractLinearFiber(idx);
}



bool C_bundle::cutAcuteEdge()
{
    //Assumption: for all points , d(v)<=2
    double ex1, ex2, ey1, ey2, ez1, ez2, cosine;
    for(unsigned long i=0 ; i<m_graph->getLenPoints() ; i++)
    {
        if(m_point[i].degree==2)
        {
            ///find the two connected edges
            short int idxe1 = -1;
            short int idxe2 = -1;
            for(long j=0 ; j<m_point[i].n_edges ; j++)
            {
                if(m_point[i].edges[j].connected)
                {
                    if(idxe1==-1)
                    {
                        idxe1 = j;
                    }
                    else
                    {
                        idxe2 = j;
                        j = m_point[i].n_edges;
                    }
                }
            }

            ///compute cosine
            ex1 = m_point[i].edges[idxe1].other_point->x - m_point[i].x;
            ey1 = m_point[i].edges[idxe1].other_point->y - m_point[i].y;
            ez1 = m_point[i].edges[idxe1].other_point->z - m_point[i].z;

            ex2 = m_point[i].edges[idxe2].other_point->x - m_point[i].x;
            ey2 = m_point[i].edges[idxe2].other_point->y - m_point[i].y;
            ez2 = m_point[i].edges[idxe2].other_point->z - m_point[i].z;

            cosine = (ex1*ex2 + ey1*ey2 + ez1*ez2)/(sqrt(SQR(ex1) + SQR(ey1) + SQR(ez1)) * sqrt(SQR(ex2) + SQR(ey2) + SQR(ez2)));

            if(cosine>0)
            {
                disconnectEdge(m_point[i].edges[idxe1].edge);
            }

        }
    }

    //unsigned long N_edge = m_graph->getLenEdges();
//    for(unsigned long i=0 ; i<nb_edges ; i++)
//    {
//        if(m_edge[i].connected)
//        {
//            ex1 = m_edge[i].p2->x - m_edge[i].p1->x;
//            ey1 = m_edge[i].p2->y - m_edge[i].p1->y;
//            ez1 = m_edge[i].p2->z - m_edge[i].p1->z;
//            //check all edges around and delete if one form acute angle with the current one (then go to following edge) (check only edges with larger indices)
//            for(unsigned long j=0 ; j< (unsigned long) m_edge[i].p1->n_edges ; j++)
//            {
//                //check if edge is connected and if it's different from the current one
//                if(m_edge[i].p1->edges[j].edge != &(m_edge[i]) && m_edge[i].p1->edges[j].edge->connected &&  (m_edge[i].p1->edges[j].edge)>m_edge+i)
//                {
//                    if(m_edge[i].p1!=m_edge[i].p1->edges[j].edge->p1)
//                    {
//                        ex2 = m_edge[i].p1->edges[j].edge->p1->x - m_edge[i].p1->x;
//                        ey2 = m_edge[i].p1->edges[j].edge->p1->y - m_edge[i].p1->y;
//                        ez2 = m_edge[i].p1->edges[j].edge->p1->z - m_edge[i].p1->z;
//                    }
//                    else
//                    {
//                        ex2 = m_edge[i].p1->edges[j].edge->p2->x - m_edge[i].p1->x;
//                        ey2 = m_edge[i].p1->edges[j].edge->p2->y - m_edge[i].p1->y;
//                        ez2 = m_edge[i].p1->edges[j].edge->p2->z - m_edge[i].p1->z;
//                    }
//                    //check angle
//                    if(ex1*ex2 + ey1*ey2 + ez1*ez2>=0)
//                    {
//                        m_edge[i].connected = !m_edge[i].connected;
//                        m_edge[i].this_in_p1->connected = m_edge[i].connected;
//                        m_edge[i].this_in_p2->connected = m_edge[i].connected;
//                        m_edge[i].p1->degree--;
//                        m_edge[i].p2->degree--;
//                    }
//                }
//            }
//
//            if(m_edge[i].connected)
//            {
//                ex1 = m_edge[i].p1->x - m_edge[i].p2->x;
//                ey1 = m_edge[i].p1->y - m_edge[i].p2->y;
//                ez1 = m_edge[i].p1->z - m_edge[i].p2->z;
//                //check all edges around and delete if one form acute angle with the current one (then go to following edge) (check only edges with larger indices)
//                for(unsigned long j=0 ; j< (unsigned long) m_edge[i].p2->n_edges ; j++)
//                {
//                    //check if edge is connected and if it's different from the current one
//                    if(m_edge[i].p2->edges[j].edge != &(m_edge[i]) && m_edge[i].p2->edges[j].edge->connected &&  (m_edge[i].p2->edges[j].edge)>m_edge+i)
//                    {
//                        if(m_edge[i].p2!=m_edge[i].p2->edges[j].edge->p1)
//                        {
//                            ex2 = m_edge[i].p2->edges[j].edge->p1->x - m_edge[i].p2->x;
//                            ey2 = m_edge[i].p2->edges[j].edge->p1->y - m_edge[i].p2->y;
//                            ez2 = m_edge[i].p2->edges[j].edge->p1->z - m_edge[i].p2->z;
//                        }
//                        else
//                        {
//                            ex2 = m_edge[i].p2->edges[j].edge->p2->x - m_edge[i].p2->x;
//                            ey2 = m_edge[i].p2->edges[j].edge->p2->y - m_edge[i].p2->y;
//                            ez2 = m_edge[i].p2->edges[j].edge->p2->z - m_edge[i].p2->z;
//                        }
//
//                        //check angle
//                        if(ex1*ex2 + ey1*ey2 + ez1*ez2>=0)
//                        {
//                            m_edge[i].connected = !m_edge[i].connected;
//                            m_edge[i].this_in_p1->connected = m_edge[i].connected;
//                            m_edge[i].this_in_p2->connected = m_edge[i].connected;
//                            m_edge[i].p1->degree--;
//                            m_edge[i].p2->degree--;
//                        }
//                    }
//                }
//            }
//        }
//    }
    //bool cond=false;
    int bbb;
    for(unsigned long i=0 ; i<nb_points ; i++)
    {
        bbb = 0;
        //m_point[i].degree = 0;
        for(unsigned long j=0 ; j<(unsigned long) m_point[i].n_edges ; j++)
        {
            if(m_point[i].edges[j].edge->connected)
            {
                bbb++;
            }
        }
        if(bbb!=m_point[i].degree)
        {
            //cond = true;
            m_point[i].degree = bbb;
        }
    }

    return true;
}



bool C_bundle::cleanGraph()
{
    Edge* eTemp;

    //clean graph: for each point which is connected to more than 2 edge, set all edges to 0
    for(unsigned long i=0 ; i<nb_points ; i++)
    {
        if(m_point[i].degree>2)
        {
            //should list all pair
            //compute angle for all pair
            //keep only the pair with smaller cosin
            for(unsigned long j=0 ; j<(unsigned long) m_point[i].n_edges ; j++)
            {
                eTemp = m_point[i].edges[j].edge;
                if(eTemp->connected)
                {
                    disconnectEdge(eTemp);
                }
            }
        }
    }
    return true;

}


vector<unsigned long> C_bundle::getStartingIdx()
{
    //Point* my_points = my_graph->getPoints();
    //unsigned long nb_points = my_graph->getLenPoints();
    unsigned long nb_connected_edges;
    vector<unsigned long> idx_start;
    for(unsigned long i=0 ; i<nb_points ; i++)
    {
        nb_connected_edges = 0;
        for(unsigned long j=0 ; j<(unsigned long) m_point[i].n_edges ; j++)
        {
            if(m_point[i].edges[j].connected)
            {
                nb_connected_edges++;
                if(nb_connected_edges>1)
                {
                    j=m_point[i].n_edges;
                }
            }
        }
        if(nb_connected_edges==1)
        {
            idx_start.push_back(i);
        }
    }
    return idx_start;

}

void C_bundle::disconnectEdge(Edge* e)
{
    if(e->connected)
    {
        (e->p2->degree)--;
        (e->p1->degree)--;
        e->connected = 0;
        e->this_in_p1->connected = 0;
        e->this_in_p2->connected = 0;
    }

}


unsigned long C_bundle::findFollowingPoint(Point* p, Edge* e)
{
    unsigned long r = 0;
    if(e!=NULL)
    {
        for(unsigned long k=0 ; k<(unsigned long) p->n_edges ; k++)
        {
            if(p->edges[k].connected)
            {
                if(e!=p->edges[k].edge)
                {
                    e = p->edges[k].edge;
                    //update i
                    if(e->p1 == p)
                    {
                        //choose p2 as new
                        r = e->p2 - m_point;
                    }else{
                        //choose p1 as next point
                        r = e->p1 - m_point;
                    }

                    disconnectEdge(e);

                    //breaking condition
                    k = p->n_edges;
                }
            }
        }
    }
    else
    {
        for(unsigned long k=0 ; k<(unsigned long) p->n_edges ; k++)
        {
            if(p->edges[k].connected)
            {
                e = p->edges[k].edge;
                //update i
                if(e->p1 == p)
                {
                    //choose p2 as new
                    r = e->p2 - m_point;
                }else{
                    //choose p1 as next point
                    r = e->p1 - m_point;
                }

                disconnectEdge(e);

                //breaking condition
                k = p->n_edges;
            }
        }
    }
    return r;

}

bool C_bundle::deleteSingle()
{
    bool cond1, cond2;
    for(unsigned long i=0 ; i<nb_edges ; i++)
    {
        if(m_edge[i].connected)
        {
            cond1 = (m_edge[i].p1->degree==1);
            cond2 = (m_edge[i].p2->degree==1);
            if(cond1 && cond2)
            {
                disconnectEdge(&(m_edge[i]));
            }
        }
    }
    return true;
}










bool C_bundle::fiberProc(void)
{
    vector< vector<unsigned long> > tempFIBERS;
    unsigned long N = FIBERS.size();
    if(N==0)
    {
        cout << "no fiber found" << endl;
        return false;
    }
    else
    {
        cout << N << " fibers have been found" << endl;
    }

    for(unsigned long i=0 ; i<N ; i++)
    {
        if(FIBERS.at(i).size()>4)
        {
            vector< vector<unsigned long> > tt = fiberProc(FIBERS.at(i));
            //tt.push_back(FIBERS.at(i));
            for(unsigned long n=0 ; n<tt.size() ; n++)
            {
                tempFIBERS.push_back(tt.at(n));
            }
        }
    }

    FIBERS = tempFIBERS;

    return true;
}

vector< vector<unsigned long> > C_bundle::fiberProc(vector<unsigned long> fib)
{
    vector< unsigned long > idxToGetRidOf; /// indices of points that are to be deleted: index in fiber
    idxToGetRidOf.clear();
    vector< vector<unsigned long> > returnedFiber;
    returnedFiber.clear();
    double dphi = 0;
    double dphi1 = 0;
    double dphi2 = 0;
    for(unsigned long j=0 ; j<fib.size()-4 ; j++)
    {
        ///compute curvature for voxel j+2
        ///compute distance between j and j+4
        dphi = sqrt(SQR(m_point[fib.at(j)].x - m_point[fib.at(j+4)].x) + SQR(m_point[fib.at(j)].y - m_point[fib.at(j+4)].y) + SQR(m_point[fib.at(j)].z - m_point[fib.at(j+4)].z));
        dphi1 = sqrt(SQR(m_point[fib.at(j+2)].x - m_point[fib.at(j+4)].x) + SQR(m_point[fib.at(j+2)].y - m_point[fib.at(j+4)].y) + SQR(m_point[fib.at(j+2)].z - m_point[fib.at(j+4)].z));
        dphi2 = sqrt(SQR(m_point[fib.at(j)].x - m_point[fib.at(j+2)].x) + SQR(m_point[fib.at(j)].y - m_point[fib.at(j+2)].y) + SQR(m_point[fib.at(j)].z - m_point[fib.at(j+2)].z));
        if(SQR(dphi)<=SQR(dphi1) + SQR(dphi2))
        {
            idxToGetRidOf.push_back(j+2);
        }
    }
    //returnedFiber.push_back(fib);
    //return returnedFiber;
    if(!idxToGetRidOf.empty())
    {
        //returnedFiber.push_back(fib);
        //return returnedFiber;
        ///compute spitting and deletion if necesary

        ///compute length of each subpart
        unsigned long n_sub = idxToGetRidOf.size();//+1;
        vector< unsigned long > tt;
        unsigned long sizeSubPart;
        for(unsigned long i=0 ; i<n_sub ; i++)
        {
            if(i==0)
            {
                sizeSubPart = idxToGetRidOf.at(i);
            }
            else if(i==n_sub-1)
            {
                sizeSubPart = (fib.size()-1) - idxToGetRidOf.at(i);
            }
            else
            {
                sizeSubPart = idxToGetRidOf.at(i+1) - idxToGetRidOf.at(i) - 1;
            }
            if(sizeSubPart>4)
            {
                if(i==0)
                {
                    for(unsigned long j=0 ; j<idxToGetRidOf.at(i) ; j++)
                    {
                        tt.push_back(fib.at(j));
                    }
                    returnedFiber.push_back(tt);
                    tt.clear();
                }
                else if(i==n_sub-1)
                {
                    for(unsigned long j=idxToGetRidOf.at(i)+1 ; j<fib.size() ; j++)
                    {
                        tt.push_back(fib.at(j));
                    }
                    returnedFiber.push_back(tt);
                    tt.clear();
                }
                else
                {
                    for(unsigned long j=idxToGetRidOf.at(i)+1 ; j<idxToGetRidOf.at(i+1) ; j++)
                    {
                        tt.push_back(fib.at(j));
                    }
                    returnedFiber.push_back(tt);
                    tt.clear();
                }
            } ///else do nothing it's going to get rid of it after by not including it
        }
    }
    else
    {
        ///return the original fib
        returnedFiber.push_back(fib);
    }
    return returnedFiber;
}

