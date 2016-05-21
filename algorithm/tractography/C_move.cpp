#include <C_move.h>

C_move::C_move(C_graph *graph)
{
    //ctor
    m_graph = graph;
    aukMove = new MovingEdge;
    aukMove->e = NULL;
    aukMove->p = NULL;
    aukMove->idxPoint = 0;
    aukMove->nb_edge = 0;
    myLandScape.e = NULL;
    myLandScape.dim = 0;
    myLandScape.sizeOfLandScape = 0;
    unsigned long long int idum = ((unsigned long long int) time(NULL))*((unsigned long long int) time(NULL));
    while(idum==THE_NUMBER_WHICH_MUST_NOT_BE_NAMED)
    {
        idum = ((unsigned long long int) time(NULL))*((unsigned long long int) time(NULL));
    }
    t = new C_toolbox_rand2(idum);
    m_mix = NULL; //new Mixture;
    T = 0;
    Tinit = 0;
    Tend = 0;
    n = 0;
    m_varSigma = 0.2;
    m_voxelConstraint = true;
    m_moveFromCenter = true;
}

C_move::~C_move()
{
    //dtor
    deleteLandScape();
    if(t!=NULL) delete t;

    if(m_mix!=NULL)
    {
        delete m_mix;
    }
    if(aukMove!=(MovingEdge*) NULL)
        delete aukMove;
}

bool C_move::initMove(long kind) //should init a large array of with the whole solution landscape
{
    if(m_mix!=NULL)
    {
        delete m_mix;
    }
    if(aukMove->e!=NULL)
    {
        delete aukMove->e;
        aukMove->nb_edge = 0;
    }
    aukMove->e = new Edge*[MAX_EDGE_MOVED_AT_ONCE];
    if(aukMove->e==NULL) return false;

    if(aukMove->p!=NULL)
    {
        delete aukMove->p;
        aukMove->idxPoint = 0;
    }
    aukMove->p = new Point[1];
    if(aukMove->p==NULL) return false;

    if(kind==MOVE_T1)
    {
        aukMove->nb_edge = 1;
    }
    else if(kind==MOVE_T2 || kind==MOVE_T1+MOVE_T2 || kind==SC_MOVE_1_2 || kind==SC_MOVE_1_2_4)
    {
        if(kind==MOVE_T1+MOVE_T2) //add it will it's no longer a test
        {
            m_mix = new Mixture;
            m_mix->alpha1 = 0.5;
            m_mix->alpha2 = 0.5;
        }
        aukMove->nb_edge = 2;
        getLandScape2();
    }
    else if(kind==MOVE_T4)
    {
        aukMove->nb_edge = 4;
    }
    else if(kind==MOVE_V) //means
    {
        //aukMove->p = new Point;
    }
    else if( (kind==(MOVE_V+MOVE_T1)) || kind==SC_MOVE_V1 || kind==SC_MOVE_V1V || kind==SC_MOVE_V1PT)
    {
        m_mix = new Mixture;
        m_mix->alpha1 = 0.5;
        m_mix->alpha2 = 0.5;
        aukMove->nb_edge = 1;
    }
    else// if(kind==MOVE_V+MOVE_T2)
    {
        m_mix = new Mixture;
        m_mix->alpha1 = 0.5;
        m_mix->alpha2 = 0.5;
        aukMove->nb_edge = 2;
        getLandScape2();
    }

    return true;
}
bool C_move::initMoveICM(long kind) //should init a large array of with the whole solution landscape
{
    if(aukMove->e!=NULL)
    {
        delete aukMove->e;
        aukMove->nb_edge = 0;
    }
    aukMove->e = new Edge*[MAX_EDGE_MOVED_AT_ONCE];
    if(aukMove->e==NULL) return false;

    if(kind==MOVE_T1)
    {
        aukMove->nb_edge = 1;
    }
    else if(kind==MOVE_T2)
    {
        aukMove->nb_edge = 2;
        getLandScape2();
    }
    else //means MOVE_T4
    {
        aukMove->nb_edge = 0;
        return false;//aukMove->nb_edge = 4;
    }
    return true;
}

bool C_move::endMove(void)
{
    if(aukMove->p!=NULL)
    {
        delete aukMove->p;
        aukMove->p = NULL;
    }
    if(aukMove->e!=NULL)
    {
        delete aukMove->e;
        aukMove->e = NULL;
    }
    deleteLandScape();
    return true;
}

bool C_move::endMoveICM(void)
{
    if(aukMove->e!=NULL)
    {
        delete aukMove->e;
        aukMove->e = NULL;
    }
    deleteLandScape();
    return true;
}

bool C_move::move(long kind) //, Mixture* mix
{
    if(kind==MOVE_T1)
    {
        return move1();
    }
    else if(kind==MOVE_T2)
    {
        return move2(); //set it back after test
    }
    else if(kind==MOVE_T4)
    {
        return move4();
    }
    else if(kind==MOVE_T1+MOVE_T2)
    {
        return move12(true);
    }
    else if(kind==SC_MOVE_1_2)
    {
        return move12(false);
    }
    else if(kind==SC_MOVE_1_2_4)
    {
        return move124();
    }
    else if(kind==MOVE_V)
    {
        return moveVertex();
    }
    else if(kind==MOVE_V+MOVE_T1)
    {
        return movev1(true, false);
    }
    else if(kind==SC_MOVE_V1)
    {
        return movev1(false, false);
    }
    else if(kind==SC_MOVE_V1V)
    {
        return movev1(false, true);
    }
    else if(kind==SC_MOVE_V1PT)
    {
        return moveRatioVertexEdgeT();
    }
    else if(kind==MOVE_V+MOVE_T2)
    {
        return movev2(true);
    }
    else
    {
        return false;
    }
}
bool C_move::moveICM(long kind) //, Mixture* mix
{
    if(kind==MOVE_T1)
    {
        return move1ICM();
    }
    else if(kind==MOVE_T2)
    {
        return move2ICM(); //set it back after test
    }
    else
    {
        return false;
    }
}
bool C_move::move1(void)
{
    //chose at random an edge
    aukMove->nb_edge = 1;
    unsigned long idx = (unsigned long) ((double) m_graph->getLenEdges()-1) * (t->doub()) ;
    aukMove->e[0] = &((m_graph->getEdges())[idx]);
    if(aukMove->e[0])
    {
        return false;
    }
    return true;
}
bool C_move::move1ICM(void)
{
    //chose at random an edge
    aukMove->nb_edge = 1;
    aukMove->e[0] = &((m_graph->getEdges())[n]);
    n++;
    if(n>=m_graph->getLenEdges())
    {
        n=0;
    }
    if(aukMove->e[0])
    {
        return false;
    }
    return true;
}
bool C_move::move2(void)
{
    if(myLandScape.sizeOfLandScape!=0)
    {
        //draw at random a pair of edge in lanscape
        aukMove->nb_edge = 2;
        unsigned long i = (unsigned long) ((t->doub())*((double) (myLandScape.sizeOfLandScape-1) ));
        aukMove->e[0] = myLandScape.e[i][0];
        aukMove->e[1] = myLandScape.e[i][1];
        if(aukMove->e[0]==NULL || aukMove->e[1]==NULL)
        {
            return false;
        }
    }
    else
    {
        return false;
    }

    return true;
}
bool C_move::move2ICM(void)
{
    if(myLandScape.sizeOfLandScape!=0)
    {
        //draw at random a pair of edge in lanscape
        aukMove->nb_edge = 2;
        ///unsigned long i = (unsigned long) ((t->doub())*((double) (myLandScape.sizeOfLandScape-1) ));
        aukMove->e[0] = myLandScape.e[n][0];
        aukMove->e[1] = myLandScape.e[n][1];
        n++;
        if(n>=myLandScape.sizeOfLandScape)
        {
            n=0;
        }
        if(aukMove->e[0]==NULL || aukMove->e[1]==NULL)
        {
            return false;
        }
    }
    else
    {
        return false;
    }

    return true;
}
bool C_move::move4(void)
{
    aukMove->nb_edge = 4;
    aukMove->e[0] = &((m_graph->getEdges())[ ( (unsigned long) ( ((double) m_graph->getLenEdges()-1) * (t->doub()) ) ) ]);

    Point* p2;
    if(t->doub()>=0.5)
    {
        ///go in one way
        p2 = aukMove->e[0]->p1;
    }
    else
    {
        ///go the other way
        p2 = aukMove->e[0]->p2;
    }

    ///draw second edge
    unsigned long idx = ( (unsigned long) ( ((double) p2->n_edges-1) * (t->doub()) ) );
    aukMove->e[1] = p2->edges[idx].edge;
    unsigned long k = 0;
    while(aukMove->e[0]==aukMove->e[1] && k<(p2->n_edges*10))
    {
        k++;
        idx = ( (unsigned long) ( ((double) p2->n_edges-1) * (t->doub()) ) );
        aukMove->e[1] = p2->edges[idx].edge;
    }
    if(k==(p2->n_edges*10))
    {
        return false;
    }

    ///draw third edge
    Point* p3 = aukMove->e[1]->p1;
    if(p3==p2)
    {
        p3 = aukMove->e[1]->p2;
    }
    idx = ( (unsigned long) ( ((double) p3->n_edges-1) * (t->doub()) ) );
    aukMove->e[2] = p3->edges[idx].edge;
    k = 0;
    while(aukMove->e[1]==aukMove->e[2] && k<(p3->n_edges*10))
    {
        k++;
        idx = ( (unsigned long) ( ((double) p3->n_edges-1) * (t->doub()) ) );
        aukMove->e[2] = p3->edges[idx].edge;
    }
    if(k==(p3->n_edges*10))
    {
        return false;
    }

    ///draw fourth edge
    Point* p4 = aukMove->e[2]->p1;
    if(p4==p3)
    {
        p4 = aukMove->e[2]->p2;
    }
    idx = ( (unsigned long) ( ((double) p4->n_edges-1) * (t->doub()) ) );
    aukMove->e[3] = p4->edges[idx].edge;
    k = 0;
    while( (aukMove->e[2]==aukMove->e[3] || aukMove->e[0]==aukMove->e[3]) && k<(p4->n_edges*10))
    {
        k++;
        idx = ( (unsigned long) ( ((double) p4->n_edges-1) * (t->doub()) ) );
        aukMove->e[3] = p4->edges[idx].edge;
    }
    if(k==(p4->n_edges*10))
    {
        return false;
    }
    return true;

}


bool C_move::move12(bool sa)
{
    if(myLandScape.sizeOfLandScape!=0)
    {
        //draw at random a pair of edge in lanscape

        if(sa)
        {
            if(t->doub()<m_mix->alpha1)
            {
                move1();
            }
            else
            {
                move2();
            }
        }
        else
        {
            if(Tend!=Tinit)
            {
                //double threshold = (T-Tend)/(Tinit-Tend);
                double threshold = 1-log(T/Tinit)/log(Tend/Tinit);/// (T-Tend)/(Tinit-Tend);
                if(t->doub()<=threshold)
                {
                    move2();
                }
                else
                {
                    move1();
                }
            }
            else
            {
                return false;
            }
        }

    }
    else
    {
        return false;
    }
    return true;
}



bool C_move::move124(void)
{
    if(myLandScape.sizeOfLandScape!=0)
    {
        if(Tend!=Tinit)
        {
            double x = log(T/Tinit)/log(Tend/Tinit); //progression bar
            //if(T<((Tinit+Tend)/2.0))
            if(x>0.5)
            {
                ///1 or 2 edges
                //double threshold = ((Tinit+Tend)-2.0*T)/(Tinit-Tend);
                double threshold = 2.0 * (x - 0.5);
                if(t->doub()>threshold)
                {
                    ///draw 2 edges
                    move2();
                }
                else
                {
                    ///draw 1 edge
                    move1();
                }
            }
            else
            {
                ///2 or 4 edges
                //double threshold = (2.0*T-(Tinit+Tend))/(Tinit-Tend);
                double threshold = 2.0 * (0.5 - x);
                if(t->doub()>threshold)
                {
                    ///draw 2 edges
                    move2();
                }
                else
                {
                    ///draw 4 edges
                    move4();
                }
            }

        }
        else
        {
            return false;
        }


    }
    else
    {
        return false;
    }
    return true;
}

bool C_move::moveVertex(void)
{
    ///notify that you are moving a vertex, not an edge
    aukMove->nb_edge=0;
    ///draw a random radius
    double mu = 0;
    double sigma = m_varSigma;//0.2;
    double u = t->doub();
    double dx = ((u) > SMALL_NUM ? (mu + ( sigma * sqrt(-2*log(u)) ) * cos(2*Pi*t->doub())) : (0));
    u = t->doub();
    double dy = ((u) > SMALL_NUM ? (mu + ( sigma * sqrt(-2*log(u)) ) * cos(2*Pi*t->doub())) : (0));
    u = t->doub();
    double dz = ((u) > SMALL_NUM ? (mu + ( sigma * sqrt(-2*log(u)) ) * cos(2*Pi*t->doub())) : (0));
    aukMove->idxPoint = (unsigned long) ((double) m_graph->getLenPoints()-1) * (t->doub()) ;
    if(m_voxelConstraint)
    {
        if(m_moveFromCenter)
        {
            aukMove->p->x = m_graph->getPoints()[aukMove->idxPoint].x0+dx;
            aukMove->p->y = m_graph->getPoints()[aukMove->idxPoint].y0+dy;
            aukMove->p->z = m_graph->getPoints()[aukMove->idxPoint].z0+dz;
        }
        else
        {
            aukMove->p->x = m_graph->getPoints()[aukMove->idxPoint].x+dx;
            aukMove->p->y = m_graph->getPoints()[aukMove->idxPoint].y+dy;
            aukMove->p->z = m_graph->getPoints()[aukMove->idxPoint].z+dz;
        }
        while(ABS(aukMove->p->x-m_graph->getPoints()[aukMove->idxPoint].x0)>=0.4999*m_graph->getPixDX())
        {
            u = t->doub();
            dx = ((u) > SMALL_NUM ? (mu + ( sigma * sqrt(-2*log(u)) ) * cos(2*Pi*t->doub())) : (0));
            if(m_moveFromCenter)
            {
                aukMove->p->x = m_graph->getPoints()[aukMove->idxPoint].x0+dx;
            }
            else
            {
                aukMove->p->x = m_graph->getPoints()[aukMove->idxPoint].x+dx;
            }
        }
        while(ABS(aukMove->p->y-m_graph->getPoints()[aukMove->idxPoint].y0)>=0.4999*m_graph->getPixDY())
        {
            u = t->doub();
            dy = ((u) > SMALL_NUM ? (mu + ( sigma * sqrt(-2*log(u)) ) * cos(2*Pi*t->doub())) : (0));
            if(m_moveFromCenter)
            {
                aukMove->p->y = m_graph->getPoints()[aukMove->idxPoint].y0+dy;
            }
            else
            {
                aukMove->p->y = m_graph->getPoints()[aukMove->idxPoint].y+dy;
            }
        }
        while(ABS(aukMove->p->z-m_graph->getPoints()[aukMove->idxPoint].z0)>=0.4999*m_graph->getPixDZ())
        {
            u = t->doub();
            dz = ((u) > SMALL_NUM ? (mu + ( sigma * sqrt(-2*log(u)) ) * cos(2*Pi*t->doub())) : (0));
            if(m_moveFromCenter)
            {
                aukMove->p->z = m_graph->getPoints()[aukMove->idxPoint].z0+dz;
            }
            else
            {
                aukMove->p->z = m_graph->getPoints()[aukMove->idxPoint].z+dz;
            }
        }
    }
    else
    {
        if(m_moveFromCenter)
        {
            aukMove->p->x = m_graph->getPoints()[aukMove->idxPoint].x0+dx;
            aukMove->p->y = m_graph->getPoints()[aukMove->idxPoint].y0+dy;
            aukMove->p->z = m_graph->getPoints()[aukMove->idxPoint].z0+dz;
        }
        else
        {
            aukMove->p->x = m_graph->getPoints()[aukMove->idxPoint].x+dx;
            aukMove->p->y = m_graph->getPoints()[aukMove->idxPoint].y+dy;
            aukMove->p->z = m_graph->getPoints()[aukMove->idxPoint].z+dz;
        }

        while(!belongToROI(aukMove->p->x, aukMove->p->y, aukMove->p->z))
        {
            u = t->doub();
            dx = ((u) > SMALL_NUM ? (mu + ( sigma * sqrt(-2*log(u)) ) * cos(2*Pi*t->doub())) : (0));
            u = t->doub();
            dy = ((u) > SMALL_NUM ? (mu + ( sigma * sqrt(-2*log(u)) ) * cos(2*Pi*t->doub())) : (0));
            u = t->doub();
            dz = ((u) > SMALL_NUM ? (mu + ( sigma * sqrt(-2*log(u)) ) * cos(2*Pi*t->doub())) : (0));
            if(m_moveFromCenter)
            {
                aukMove->p->x = m_graph->getPoints()[aukMove->idxPoint].x0+dx;
                aukMove->p->y = m_graph->getPoints()[aukMove->idxPoint].y0+dy;
                aukMove->p->z = m_graph->getPoints()[aukMove->idxPoint].z0+dz;
            }
            else
            {
                aukMove->p->x = m_graph->getPoints()[aukMove->idxPoint].x+dx;
                aukMove->p->y = m_graph->getPoints()[aukMove->idxPoint].y+dy;
                aukMove->p->z = m_graph->getPoints()[aukMove->idxPoint].z+dz;
            }
        }
    }
    return true;
}

bool C_move::belongToROI(double x, double y, double z)
{
    signed long i= (signed long) (x/m_graph->getMask()->pixDimX);
    signed long j= (signed long) (y/m_graph->getMask()->pixDimY);
    signed long k= (signed long) (z/m_graph->getMask()->pixDimZ);
    if(i<0) return false;
    if(i>=(signed long)m_graph->getMask()->DimX) return false;
    if(j<0) return false;
    if(j>=(signed long)m_graph->getMask()->DimY) return false;
    if(k<0) return false;
    if(k>=(signed long)m_graph->getMask()->DimZ) return false;

    if(m_graph->getMask()->raw3D[i][j][k]>0.5) return true;

    return false;
}

bool C_move::movev1(bool sa, bool addLargeMove)
{
    //draw at random a pair of edge in lanscape

    if(sa)
    {
        if(t->doub()<m_mix->alpha1)
        {
            move1();
        }
        else
        {
            moveVertex();
        }
    }
    else
    {
        if(Tend!=Tinit)
        {
            if(t->doub()<m_mix->alpha1)
            {
                move1();
            }
            else
            {
                ///update m_varSigma for the current iteration
                double varSigmaOpt = 0.3*m_graph->getPixDX();
                if(addLargeMove)
                {
                    if(T>Tinit) m_varSigma = 1.5*varSigmaOpt;
                    else if(T<Tend) m_varSigma = 0.1*varSigmaOpt;
                    else m_varSigma = 1.5*varSigmaOpt - 1.4*(log(T/Tinit)/log(Tend/Tinit))*varSigmaOpt;
                }
                else
                {
                    if(T>Tinit) m_varSigma = 1.5*varSigmaOpt;
                    else if(T<Tend) m_varSigma = varSigmaOpt;
                    else m_varSigma = 1.5*varSigmaOpt - 1.4*(log(T/Tinit)/log(Tend/Tinit))*varSigmaOpt;
                }
                moveVertex();
            }
        }
        else
        {
            return false;
        }
    }
    return true;
}

bool C_move::moveRatioVertexEdgeT(void)
{
    if(Tend!=Tinit)
    {
        if(false)
        {
            if(T>Tinit) m_mix->alpha1 = 1.0;
            else if(T<Tend) m_mix->alpha1 = 0.975;
            else m_mix->alpha1 = 1.0 - 0.025*(log(T/Tinit)/log(Tend/Tinit));
        }
        else
        {
            if(T>Tinit) m_mix->alpha1 = 1.0;
            else if(T<Tend) m_mix->alpha1 = 0.95;
            else
            {
                double x = log(T/Tinit)/log(Tend/Tinit);
                if(x<0.5)
                {
                    m_mix->alpha1 = 1.0;
                }
                else
                {
                    m_mix->alpha1 = 1.0 - 0.05*(2.0*x-1.0);
                }

            }
        }


        m_mix->alpha2 = 1.0 - m_mix->alpha1;
        if(t->doub()<m_mix->alpha1)
        {
            move1();
        }
        else
        {
            ///update m_varSigma for the current iteration
            double varSigmaOpt = 0.3*m_graph->getPixDX();
            //m_varSigma = 1.5*varSigmaOpt - 1.4*(log(T/Tinit)/log(Tend/Tinit))*varSigmaOpt;
            if(T>Tinit) m_varSigma = 1.5*varSigmaOpt;
            else if(T<Tend) m_varSigma = 0.1*varSigmaOpt;
            else m_varSigma = 1.5*varSigmaOpt - 1.4*(log(T/Tinit)/log(Tend/Tinit))*varSigmaOpt;
            moveVertex();
        }
    }
    else
    {
        return false;
    }
}


bool C_move::movev2(bool sa)
{
    if(myLandScape.sizeOfLandScape!=0)
    {
        //draw at random a pair of edge in lanscape

        if(sa)
        {
            if(t->doub()<m_mix->alpha1)
            {
                move2();
            }
            else
            {
                moveVertex();
            }
        }
        else
        {
            ///to be implemented later
            return false;
        }

    }
    else
    {
        return false;
    }
    return true;
}







bool C_move::getLandScape2(void)
{
    //get all posssible edge pairs
    vector<unsigned long> idx1;
    vector<unsigned long> idx2;
    unsigned long nbEdge = m_graph->getLenEdges();
    unsigned long idxTemp;
    Edge* tr = m_graph->getEdges();
    Point* p;
    for(unsigned long i=0 ; i<nbEdge ; i++)
    {
        //check all edges neibouring current edge
        //on first side
        p = tr[i].p1;
        for(unsigned long j=0 ; j<(unsigned long) p->n_edges ; j++)
        {
            //compute current edge idx
            idxTemp = (p->edges[j]).edge - tr;
            if(idxTemp>i)
            {
                //store edge
                idx1.push_back(i);
                idx2.push_back(idxTemp);
            }
        }
        //on second side
        p = tr[i].p2;
        for(unsigned long j=0 ; j<(unsigned long) p->n_edges ; j++)
        {
            //compute current edge idx
            idxTemp = (p->edges[j]).edge - tr;
            if(idxTemp>i)
            {
                //store edge
                idx1.push_back(i);
                idx2.push_back(idxTemp);
            }
        }

    }

    //allocate a 2D array of edge addresses
    if(myLandScape.e!=NULL)
    {
        for(unsigned long i=0 ; i<myLandScape.sizeOfLandScape ; i++)
        {
            delete myLandScape.e[i];
        }
        delete myLandScape.e;
        myLandScape.e = NULL;
    }
    myLandScape.dim = 2;
    myLandScape.sizeOfLandScape = idx1.size();
    myLandScape.e = new Edge**[myLandScape.sizeOfLandScape];

    for(unsigned long i=0 ; i<myLandScape.sizeOfLandScape ; i++)
    {
        myLandScape.e[i] = new Edge*[2];
        myLandScape.e[i][0] = &(tr[idx1.at(i)]);
        myLandScape.e[i][1] = &(tr[idx2.at(i)]);
    }
    idx1.erase(idx1.begin(), idx1.end());
    idx2.erase(idx2.begin(), idx2.end());
    return true;
}


bool C_move::deleteLandScape(void)
{
    if(myLandScape.e!=NULL)
    {
        for(unsigned long i=0 ; i<myLandScape.sizeOfLandScape ; i++)
        {
            delete myLandScape.e[i];
        }
        delete myLandScape.e;
        myLandScape.e = NULL;
    }
    myLandScape.dim = 0;
    myLandScape.sizeOfLandScape = 0;
    return true;
}
