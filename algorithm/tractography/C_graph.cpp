#include <C_graph.h>

C_graph::C_graph(rawData<double>* mask)
{
    //ctor
    m_points = NULL;
    m_edges = NULL;
    nb_points=0;
    nb_edges=0;
    m_mask = mask;
    pix_dx = mask->pixDimX;
    pix_dy = mask->pixDimY;
    pix_dz = mask->pixDimZ;
    DX = pix_dx*mask->DimX;
    DY = pix_dy*mask->DimY;
    DZ = pix_dz*mask->DimZ;
    //create the graph from a mask: rawData
    //nb_points = countNbPointInMask(mask);
    m_points = createPointArrayFromMask();
    nb_edges = countNbEdgesInMask();
    m_edges = createEdgeArray();
    bindPointsAndEdges();//, mask->pixDimX, mask->pixDimY, mask->pixDimZ);

}
C_graph::C_graph(C_graph &G)
{
    this->nb_points = G.getLenPoints();
    this->nb_edges = G.getLenEdges();
    copyGraph(&(this->m_points), &(this->m_edges), G.getPoints(), G.getEdges(), this->nb_points, this->nb_edges);
}
bool C_graph::copyGraph(Point** points, Edge** edges, Point* m_points_source, Edge* m_edges_source, unsigned long Np, unsigned long Ne)
{
    unsigned long idx_f = 0;
    unsigned long idx_e = 0;
    *points = new Point[Np];
    *edges = new Edge[Ne];

    //copy points
    for (unsigned long i=0;i<Np;i++) //Read points, count connections
	{
		//copy coordinate of point into new point set
		(*points)[i].x=m_points_source[i].x;
		(*points)[i].y=m_points_source[i].y;
		(*points)[i].z=m_points_source[i].z;
		#ifdef HANDLE_VERTEX_MOVE
		(*points)[i].x0=m_points_source[i].x0;
		(*points)[i].y0=m_points_source[i].y0;
		(*points)[i].z0=m_points_source[i].z0;
                #endif
		//determine the number of possible edges and the number of connected edge at point i
		(*points)[i].n_edges=m_points_source[i].n_edges;
		(*points)[i].degree=m_points_source[i].degree;
		((*points)[i]).edges= new Vector[(*points)[i].n_edges];
	}

	//copy edges and update points
	for (unsigned long i=0;i<Np;i++) //Store edges and connections
    {
        for(int j=0;j<(*points)[i].n_edges;j++)
        {
            //calculate coordinates of Vector edge of point i
            (*points)[i].edges[j].dx=m_points_source[i].edges[j].dx;
            (*points)[i].edges[j].dy=m_points_source[i].edges[j].dy;
            (*points)[i].edges[j].dz=m_points_source[i].edges[j].dz;

            //store the address of the jth edge of point i
            idx_f = ((m_points_source[i].edges[j].edge) - m_edges_source);
            (*points)[i].edges[j].edge=&((*edges)[idx_f]);


            //store the address of the point connected by this edge
            idx_e = ((m_points_source[i].edges[j].other_point) - m_points_source);
            (*points)[i].edges[j].other_point=&((*points)[idx_e]);


            (*points)[i].edges[j].connected=m_points_source[i].edges[j].connected;
            //fill edge structure
            if (idx_e > (long)i) //If not, we have already processed this edge
            {
                (*edges)[idx_f].this_in_p1=&((*points)[i].edges[j]);
                (*edges)[idx_f].p1=&((*points)[i]);
                (*edges)[idx_f].p2=&((*points)[idx_e]);
                (*edges)[idx_f].connected=m_edges_source[idx_f].connected;
                //(*edges)[idx_f].data_energy=m_edges_source[idx_f].data_energy;
            }
            else
                (*edges)[idx_f].this_in_p2=&((*points)[i].edges[j]);

        }
    }
    return true;
}

C_graph::~C_graph()
{
    //dtor
    if(m_points!=NULL)
    {
        delete m_points;
    }
    if(m_edges!=NULL)
    {
        delete m_edges;
    }
}



bool C_graph::saveEdges(char * filename)
{
    SaveEdges *save_edges = new SaveEdges[nb_edges];
	for (unsigned long i=0;i< nb_edges;i++)
	{
		save_edges[i].x1=m_edges[i].p1->x;
		save_edges[i].y1=m_edges[i].p1->y;
		save_edges[i].z1=m_edges[i].p1->z;
		save_edges[i].x2=m_edges[i].p2->x;
		save_edges[i].y2=m_edges[i].p2->y;
		save_edges[i].z2=m_edges[i].p2->z;
		save_edges[i].connected=m_edges[i].connected;
	}
	FILE *f=fopen(filename,"wb");
	fwrite(save_edges,sizeof(SaveEdges),nb_edges,f);
	fclose(f);
	delete save_edges;
    return true;
}


Edge* C_graph::findEdge(Point* p1, Point* p2)
{
    for(int i=0 ; i<p1->n_edges ; i++)
    {
        if(p1->edges[i].other_point==p2)
        {
            return p1->edges[i].edge;
        }
    }
    return NULL;
}
Point* C_graph::findOtherPoint(Edge* e, Point* p)
{
    if(p!=e->p1)
    {
        return e->p1;
    }
    if(p!=e->p2)
    {
        return e->p2;
    }
    return NULL;
}
Point* C_graph::findCommunPoint(Edge* e1, Edge* e2)
{
    if(e1->p1==e2->p1) return e1->p1;
    if(e1->p1==e2->p2) return e1->p1;
    if(e1->p2==e2->p1) return e1->p2;
    if(e1->p2==e2->p2) return e1->p2;

    return NULL;
}

unsigned long C_graph::getNbConnectedEdge(void)
{
    unsigned long l = 0;
    for(unsigned long i=0 ; i<nb_edges ; i++)
    {
        if(m_edges[i].connected)
        {
            l++;
        }
    }

    return l;
}
void C_graph::randInit(double alphaConnected)
{
    C_toolbox_rand2 *p = new C_toolbox_rand2(rand());
    for(unsigned long i=0 ; i<nb_edges ; i++)
    {
        if(p->doub()<alphaConnected)
        {
            m_edges[i].connected = true;
            m_edges[i].this_in_p1->connected = true;
            m_edges[i].this_in_p2->connected = true;
        }
        else
        {
            m_edges[i].connected = false;
            m_edges[i].this_in_p1->connected = false;
            m_edges[i].this_in_p2->connected = false;
        }
    }
    //set the actual degree for each point
    for(unsigned long i=0 ; i<nb_points ; i++)
    {
        m_points[i].degree = 0;
        for(unsigned long n=0 ; n<(unsigned long) m_points[i].n_edges ; n++)
        {
            if(m_points[i].edges[n].connected)
            {
                m_points[i].degree++;
            }
        }
    }

    delete p;
    return;
}



unsigned long C_graph::countNbPointInMask(/**rawData<double>* maskor a segmented rawData*/)
{
    ///check dimension
    if(m_mask->numDim==0) return 0;
    if(m_mask->numDim>4) return 0;

    unsigned long N=0;

    if(m_mask->numDim==2)
    {
        ///mosa case?
        for(unsigned long i=0 ; i<m_mask->DimX ; i++)
        {
            for(unsigned long j=0 ; j<m_mask->DimY ; j++)
            {
                if(m_mask->raw2D[i][j]>0.5)
                {
                    N++;
                }
            }
        }
    }

    if(m_mask->numDim==3)
    {
        ///usual case (unmosa)
        for(unsigned long i=0 ; i<m_mask->DimX ; i++)
        {
            for(unsigned long j=0 ; j<m_mask->DimY ; j++)
            {
                for(unsigned long k=0 ; k<m_mask->DimZ ; k++)
                {
                    if(m_mask->raw3D[i][j][k]>0.5)
                    {
                        N++;
                    }
                }
            }
        }
    }
    return N;
}

Point* C_graph::createPointArrayFromMask(/**rawData<double>* maskor a segmented rawData*/)
{
    nb_points = countNbPointInMask();
    if(nb_points==0) return NULL;

    if(m_mask->numDim!=3) return NULL; //maybe later handle other cases  && mask->numDim!=2

    Point* p = new Point[nb_points];

    ///for each element of the container, check whether it is in mask and store it in point array if true
    unsigned long n=0;
    for(signed long k=0 ; k<(signed long) m_mask->DimZ ; k++)
    {
        for(signed long j=0 ; j<(signed long) m_mask->DimY ; j++)
        {
            for(signed long i=0 ; i<(signed long) m_mask->DimX ; i++)
            {
                if(m_mask->raw3D[i][j][k]>0)
                {
                    //point location
                    p[n].x = (((double) i)+0.5)*m_mask->pixDimX;
                    p[n].y = (((double) j)+0.5)*m_mask->pixDimY;
                    p[n].z = (((double) k)+0.5)*m_mask->pixDimZ;
                    #ifdef HANDLE_VERTEX_MOVE
                    p[n].x0 = p[n].x;
                    p[n].y0 = p[n].y;
                    p[n].z0 = p[n].z;
                    #endif
                    n++;
                    if(n==nb_points)
                    {
                        i = m_mask->DimX;
                        j = m_mask->DimY;
                        k = m_mask->DimZ;
                    }
                }
            }
        }
    }

    ///count number of neighbors
    for(unsigned long i=0 ; i<nb_points ; i++)
    {
        //anotherCounter += p[i].n_edges;
        ///check all later point
        signed long nbVector = 0;
        for(unsigned long ii=0 ; ii<nb_points ; ii++)
        {
            if(ii!=i)
            {
                if( (ABS(p[i].x-p[ii].x)<=1.0001*pix_dx)/**1.001*/ && (ABS(p[i].y-p[ii].y)<=1.0001*pix_dy)/**1.001*/ && (ABS(p[i].z-p[ii].z)<=1.0001*pix_dz)/**1.001*/)
                {
                    //indenet nbVector
                    nbVector++;
                    if(nbVector==MAX_EDGES_PER_NODE)
                    {
                        ii=nb_points;
                    }
                }
            }
        }

        p[i].n_edges = nbVector;
        p[i].degree = nbVector; //all edges are connected
        //allocate connector to edge
        p[i].edges = new Vector[nbVector];
        for(int d=0 ; d<p[i].n_edges ; d++)
        {
            p[i].edges[d].other_point=NULL;
            p[i].edges[d].edge = NULL;
        }

        ///define data vector: for new graph
        if(nbVector>1)
        {
            p[i].data = new double[(nbVector*(nbVector-1))/2];
            ///fill the array of data as follow: for the edge i and the edge j, with i<j, the data is stored at index (i+1)*nn - (i+1)*(i+2)/2 + 1 - (nn - (j+1)) -2
            unsigned long idx;
            for(unsigned long ii=0 ; ii<(unsigned long) p[i].n_edges ; ii++)
            {
                for(unsigned long jj=ii+1 ; jj<(unsigned long) p[i].n_edges ; jj++)
                {
                    idx = (ii+1)*p[i].n_edges - (ii+1)*(ii+2)/2 + 1 - (p[i].n_edges - (jj+1)) - 2;
                    p[i].data[idx] = idx;
                }
            }
        }
        else
        {
            cout << "WILL GENERATE A NULL ARRAY AT POINT: " << n << endl;
            p[i].data = NULL;
        }
    }

    return p;
}


unsigned long C_graph::countNbEdgesInMask(void)
{
    unsigned long Ne = 0;
    for(unsigned long i=0 ; i<nb_points ; i++)
        Ne += m_points[i].n_edges; //each edge is counted exactly twice (once for each points that share the edge)

    return (Ne/2);
}

Edge* C_graph::createEdgeArray()
{
    return new Edge[nb_edges];
}


bool C_graph::bindPointsAndEdges(void) /**Point* p, unsigned long Np, Edge* e, unsigned long Ne, double pixSizeX, double pixSizeY, double pixSizeZ*/
{
    unsigned long edgeCounter = 0;//, anotherCounter = 0;

    ///for each point, connect vectors to other points
//    cout << "runing on points" << endl;
    for(unsigned long i=0 ; i<nb_points ; i++)
    {
        //anotherCounter += p[i].n_edges;
        ///check all later point
        signed long nbVector = 0;
        for(unsigned long ii=0 ; ii<nb_points ; ii++)
        {
            if(ii!=i)
            {
                if( (ABS(m_points[i].x-m_points[ii].x)<=1.0001*pix_dx)/**1.001*/ && (ABS(m_points[i].y-m_points[ii].y)<=1.0001*pix_dy)/**1.001*/ && (ABS(m_points[i].z-m_points[ii].z)<=1.0001*pix_dz)/**1.001*/)
                {
                    m_points[i].edges[nbVector].other_point = &(m_points[ii]);
                    m_points[i].edges[nbVector].connected = 1;
                    m_points[i].edges[nbVector].dx = m_points[i].edges[nbVector].other_point->x - m_points[i].x;
                    m_points[i].edges[nbVector].dy = m_points[i].edges[nbVector].other_point->y - m_points[i].y;
                    m_points[i].edges[nbVector].dz = m_points[i].edges[nbVector].other_point->z - m_points[i].z;
                    m_points[i].edges[nbVector].edge = NULL;

                    //indenet nbVector
                    nbVector++;
                    if(nbVector==m_points[i].n_edges)
                    {
                        ii=nb_points;
                    }
                }
            }
        }
    }

    ///for each point in graph
    edgeCounter = 0;
    for(unsigned long i=0 ; i<nb_points ; i++)
    {
        ///for each vector
        for(int j=0 ; j<m_points[i].n_edges ; j++ )
        {
            Point* p1 = &(m_points[i]);
            Point* p2 = m_points[i].edges[j].other_point;

            ///find the vector of p2 that bind p1
            int k = -1;
            for(int jj=0 ; jj<p2->n_edges ; jj++)
            {
                if(p2->edges[jj].other_point==p1) //((abs(p1->x-p2->x)<=0.001*pixSizeX) && (abs(p1->y-p2->y)<=0.001*pixSizeY) && (abs(p1->z-p2->z)<=0.001*pixSizeZ))
                {
                    k = jj;
                    jj = p2->n_edges;
                }
            }
            if(k==-1)
            {
                cout << "throwing an error" << endl;
                throw "graph corrupted";
            }

            if( p1->edges[j].edge==NULL && p2->edges[k].edge==NULL)
            {
                //set a new common edge
                p1->edges[j].edge = &(m_edges[edgeCounter]);
                p2->edges[k].edge = &(m_edges[edgeCounter]);
                p1->edges[j].edge->connected = 1;
                p1->edges[j].edge->p1 = p1;
                p1->edges[j].edge->p2 = p2;
                p1->edges[j].edge->this_in_p1 = &(p1->edges[j]);
                p1->edges[j].edge->this_in_p2 = &(p2->edges[k]);
                edgeCounter++;
            }
            else
            {
                if( (p1->edges[j].edge!=NULL) && (p2->edges[k].edge!=NULL))
                {
                    //cout << "pb occured with edge dedication" << endl;
                }
                else if(p1->edges[j].edge!=NULL)
                {
                    p2->edges[k].edge = p1->edges[j].edge;
                }
                else if(p2->edges[k].edge!=NULL)
                {
                    p1->edges[j].edge = p2->edges[k].edge;
                }
            }
        }
    }
    return true;
}



long C_graph::getIndexOfPair(Point* p, long idx1, long idx2)
{
    long i = idx1;
    long j = idx2;
    if(i>j) swap(i,j);
    return (i+1)*p->n_edges - (i+1)*(i+2)/2 + 1 - (p->n_edges - (j+1)) - 2;
}

