#include <C_energy.h>

C_energy::C_energy(C_graph *graph, rawData<double>* _rawData, rawData<double>* _rawDataMask)
{
    //ctor
    m_nbInterpolatingPoint = 100;
    m_graph = graph;
    m_alpha = 0.2;
    m_beta = 0.5;
    m_beta0 = m_beta;
    m_sigma = 0.375;
    m_t0 = Pi/2.0;
    m_phi0 = PHI0;
    m_phi1 = PHI1;
    METHOD_INTERPOLATION = GAUSSIAN_PHI;
    OBJECT_INTERPOLATED = TENSOR;
    METHOD_INTEGRATION = SIMPSON_RULE;
    m_rawData = _rawData;
    m_rawDataMask = _rawDataMask;
    Hpsi=NULL;
    b = 0.0;

    THREAD_DELTA = NULL;
    m_R = 1.5;
}

C_energy::~C_energy()
{
    //dtor
}
#ifdef MY_TEST_MULTITHREAD
void C_energy::initEnergy(long kind, unsigned long nb_thread)
#else
void C_energy::initEnergy(long kind)
#endif
{
    #ifdef MY_TEST_MULTITHREAD
    initPairEdgeData(kind, nb_thread);
    #else
    initEnergy(kind, 0, m_graph->getLenPoints());
    #endif
    return;
}


double C_energy::initEnergy(unsigned long kind, unsigned long idxStart, unsigned long idxEnd)
{
    if(kind==U15 || kind==U1 || kind==U18)
    {
        //init data container
        return initPairEdgeDataMultiThread(kind, idxStart, idxEnd);
    }
    else
    {
        cout << "initialization failed: energy function does not match any defined function" << endl;
    }
    return 0.0;
}

bool C_energy::initPointData(long kind)
{
    if(m_graph!=NULL)
    {
        unsigned long nb_points = m_graph->getLenPoints();
        if(nb_points==0)
        {
            cout << "no point in graph" << endl;
            return false;
        }
        Point* points = m_graph->getPoints();
        if(points==NULL)
        {
            cout << "NULL pointer for points in graph" << endl;
            return false;
        }
        if(kind==U1)
        {
            for (unsigned long i=0;i<nb_points;i++)
            {
                points[i].error_sum = 0;
                for (unsigned long j=0;j< (unsigned long)points[i].n_edges;j++)
                {
                    if(points[i].edges[j].connected)
                    {
                        for(unsigned long k=0;k<j;k++)
                        {
                            if (points[i].edges[k].connected)
                            {
                                points[i].error_sum+=U1data(&points[i],&points[i].edges[j],&points[i].edges[k]);///why aren't we adding contribution of 0/1 edge?
                            }
                        }
                    }
                }
            }
        }
        else if(kind==U15 || kind==U18)
        {
            for (unsigned long i=0; i<nb_points ;i++)
            {
                points[i].error_sum = 0;
                if(points[i].degree>1)
                {
                    for(unsigned long j=1 ; j<(unsigned long) points[i].n_edges ; j++) //starts at 1 instead of 0 because there is nothing to compute at 0
                    {
                        if(points[i].edges[j].connected)
                        {
                            for(unsigned long k=0 ; k<j ; k++)
                            {
                                if(points[i].edges[k].connected)
                                {
//                                    cout << "point: " << i << " edges: " << j << " " << k << endl;
//                                    cout << "point: " << &points[i] << " edges: " << points[i].edges[j].edge << " " << points[i].edges[k].edge << endl;
//                                    cout << "point sum: " << points[i].error_sum << endl;
                                    if(points[i].edges[j].edge==NULL || points[i].edges[k].edge==NULL) cout << " an edge is NULL, there are " << points[i].n_edges << " edges on point " << i << endl;
                                    if(kind==U15 || kind==U18)
                                    {
                                        points[i].error_sum += U15data(points[i].edges[j].edge, points[i].edges[k].edge);
                                    }
                                    else
                                    {
                                        points[i].error_sum += U16data(points[i].edges[j].edge, points[i].edges[k].edge);
                                    }
//                                    cout << "point sum: " << points[i].error_sum << endl;
                                }
                            }
                        }
                    }
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

double C_energy::getGlobalEnergy(long kind, double* Edata, double* Etopo)
{
    *Edata = 0;
    *Etopo = 0;
    unsigned long nb_points = m_graph->getLenPoints();
    Point *points = m_graph->getPoints();

    if(kind==U1)
    {
        double Utemp;
        for (unsigned long i=0;i<nb_points;i++)
        {
            if(points[i].degree>=2){
                if(points[i].degree>2)
                {
                    *Etopo += (points[i].degree-2);
                }
                Utemp = 0;
                for (unsigned long j=0;j< (unsigned long)points[i].n_edges;j++)
                {
                    if(points[i].edges[j].connected)
                    {
                        for(unsigned long k=0;k<j;k++)
                        {
                            if(points[i].edges[k].connected)
                            {
                                Utemp += U1data(&points[i],&points[i].edges[j],&points[i].edges[k]);
                            }
                        }
                    }
                }
                Utemp = 2.0*Utemp/((double) points[i].degree * (double) (points[i].degree-1));
            }
            else if(points[i].degree==1)
            {
                Utemp = F1;
            }
            else
            {
                Utemp = F0;
            }
            *Edata += Utemp;
        }
        return (*Edata) + m_alpha*(*Etopo);
    }
    else if(kind==U15 || kind==U18)
    {
        //
        double dataEnergy;
        for (unsigned long i=0; i<nb_points ;i++)
        {
            dataEnergy = 0;
            *Etopo += U15topo(&points[i]);

            if(points[i].degree>1)
            {
                for(unsigned long j=1 ; j<(unsigned long) points[i].n_edges ; j++)
                {
                    if(points[i].edges[j].connected)
                    {
                        for(unsigned long k=0 ; k<j ; k++)
                        {
                            if(points[i].edges[k].connected)
                            {
                                dataEnergy += U15data(points[i].edges[j].edge, points[i].edges[k].edge);
                            }
                        }
                    }
                }
                dataEnergy = 2.0*dataEnergy/((double) (points[i].degree*(points[i].degree-1)));
            }
            *Edata += dataEnergy;
        }
        //double r = (*Edata) + (*Etopo);
//        if(r!=r) cout << "NaN is coming" << endl;
        return (*Edata) + (*Etopo);
    }
    else
    {
        return 0;
    }

}


//calling method
double C_energy::deltaEdgeU(Edge* e, double* partial_sum_p1, double* partial_sum_p2, long kind)
{
    if(kind==U1){
        return deltaEdgeU1(e, partial_sum_p1, partial_sum_p2);
    }
    else if(kind==U15 || kind==U18)
    {
        double s = deltaEdgeU15(e, partial_sum_p1, partial_sum_p2);
        return s;
    }
    else
    {
        throw "energy is not defined";
    }
}

//private methods: compute calculation

double C_energy::U1data(Point * point, Vector * edge1, Vector * edge2)
{
    ///look for the commun point
    Point* p_e1_e2 = m_graph->findCommunPoint(edge1->edge, edge2->edge);
    ///get edge indices (local indices)
    long i=-1, j=-1, N=(signed long) p_e1_e2->n_edges;
    for(long n=0 ; n<N ; n++)
    {
        if(p_e1_e2->edges[n].edge==edge1->edge) i=n;
        if(p_e1_e2->edges[n].edge==edge2->edge) j=n;
    }
    if(i==-1 || j==-1 || i==j) return 0.0;
    ///calculate index for data array
    if(j<i) swap(i,j); //we must have i<j
    return p_e1_e2->data[(i+1)*N - (i+1)*(i+2)/2 + 1 - (N - (j+1)) - 2];
}

double C_energy::deltaEdgeU1(Edge* e, double* partial_sum_p1, double* partial_sum_p2)
{
    return deltaU1(e->p1, e->this_in_p1, partial_sum_p1) + deltaU1(e->p2, e->this_in_p2, partial_sum_p2);
}
double C_energy::deltaU1(Point *point, Vector * edge, double * partial_sum)
{
    return deltaU1data(point, edge, partial_sum) + deltaU1topo(point, edge);
}
double C_energy::deltaU1data(Point *point, Vector * edge, double * partial_sum)
{
    //
    double delta;
	*partial_sum=0.;
	//for all possible edges of the conexity around the point *point (but those which are connected if any)
	for (int i=0,j=0;j<point->degree-edge->connected;i++){
	    //if the current edge is connected and if it is not the edge that is being considered, then computed the data fidefity term and add it to the partial_sum
		if (point->edges[i].connected && &point->edges[i] != edge )
		{
			*partial_sum+=U1data(point,edge,&point->edges[i]);
			j++;
		}
	}

	//depending on the conection of the considered edge, compute the defference of energy: delta. But I can't get how is it computed (caculus)
	if(edge->connected)
	{
		int degree=point->degree; //degree in connected state //number of edges coming out from the vertex point
		if (degree>2) //if more than 2 edges at vertex point
			delta=2.0*(( ((double) point->error_sum)-*partial_sum)/((((double) degree)-1.0)*(((double) degree)-2.0)) - ((double) point->error_sum)/(((double) degree)*(((double) degree)-1.0)));// -m_alpha; // I need to know what is error_sum
		else if (degree==2) //if exactly 2 edges //delta=1. - *partial_sum;
			delta = F1 - *partial_sum;
		else delta = -F0; //delta = -0.5; //if no edge or the end of a bundle
	} else
	{
		int degree=point->degree+1; //degree in connected state
		if (degree>2)
			delta=-2.0*( ((double) point->error_sum)/((((double) degree)-1.0)*(((double) degree)-2.0)) - (((double) point->error_sum)+*partial_sum)/(((double) degree)*(((double) degree)-1.0)));// +m_alpha;
		else if (degree==2)
			delta=(*partial_sum) - F1; //delta=(*partial_sum) - 1.;
		else delta = F0; //delta = 0.5;
	}
	return delta;
}
double C_energy::deltaU1topo(Point *point, Vector * edge)
{
    if(edge->connected)
    {
        int degree=point->degree; //degree in connected state //number of edges coming out from the vertex point
        if (degree>2){
            return -m_alpha;
        }else{
            return 0;
        }
    }else
    {
        int degree=point->degree+1; //degree in connected state //number of edges coming out from the vertex point
        if (degree>2){
            return m_alpha;
        }else{
            return 0;
        }
    }
}







double C_energy::U15data(Edge* e1, Edge* e2)
{
    if(e1==NULL || e2==NULL) throw -1.0;
    ///look for the commun point
    Point* p_e1_e2 = m_graph->findCommunPoint(e1, e2);
    if(p_e1_e2==NULL)
    {
        cout << "error in U15data: commun point does not exist!" << endl;
        return 0.0;
    }
    ///get edge indices (local indices)
    long i=-1, j=-1, N=(signed long) p_e1_e2->n_edges;
    for(long n=0 ; n<N ; n++)
    {
        if(p_e1_e2->edges[n].edge==e1) i=n;
        if(p_e1_e2->edges[n].edge==e2) j=n;
    }
    if(i==-1 || j==-1 || i==j) return 0.0;
    ///calculate index for data array
    if(j<i) //swap(i,j); //we must have i<j
    {
        long temp = i;
        i=j;
        j=temp;
    }
    return -p_e1_e2->data[(i+1)*N - (i+1)*(i+2)/2 + 1 - (N - (j+1)) - 2];
}
double C_energy::U15topo(Point* p)
{
    double r;// = phi(p->degree);;
    if(p->degree==2)
    {
        //find both connected edges
        Edge** e = new Edge*[2];
        unsigned long j=0;
        for(unsigned long i=0 ; i<(unsigned long) p->n_edges ; i++)
        {
            if(p->edges[i].edge->connected)
            {
                e[j] = p->edges[i].edge;
                j++;
                if(j==2)
                {
                    i = p->n_edges;
                }
            }
        }
        //return a function of cos angle
        r = m_beta*f3(angleFunc2(e[0], e[1]), m_t0, m_sigma);
        delete e;
    }
    else
    {
        r = m_alpha*phi_2(p->degree);
    }
    //if(isinf(r) || r!=r) cout << "pb in topo" << endl;
    return r;
}
double C_energy::deltaEdgeU15(Edge* e, double* partial_sum_p1, double* partial_sum_p2)
{
    return deltaU15(e->p1, e->this_in_p1, partial_sum_p1) + deltaU15(e->p2, e->this_in_p2, partial_sum_p2);
}
double C_energy::deltaU15(Point *point, Vector * edge, double * partial_sum)
{
    return deltaU15data(point, edge, partial_sum) + deltaU15topo(point, edge);
}
double C_energy::deltaU15data(Point *point, Vector * edge, double * partial_sum)
{
    if(point->degree==0)
    {
        *partial_sum = 0;
        return 0.0;
    }
    if(point->degree==1 && edge->connected)
    {
        *partial_sum = 0;
        return 0.0;
    }

    if(point->degree==1 && !edge->connected)
    {
        //do something special while going from 1 -> 2
        *partial_sum = 0;
        for(unsigned long n=0 ; n<(unsigned long) point->n_edges ; n++) //sum over n
        {
            if(point->edges[n].edge!=edge->edge) //n different from k
            {
                if(point->edges[n].edge->connected)
                {
                    *partial_sum += U15data(point->edges[n].edge, edge->edge); /// positive //UrawData(point->edges[n].edge, edge->edge, U15);
                    n = point->n_edges;
                }
            }
        }
        return *partial_sum;
    }

    if(point->degree==2 && edge->connected)
    {
        //do something special while going from 2 -> 1
        *partial_sum = 0;
        ///n=1 n+1=0
        double temp = 0;
        for(unsigned long n=0 ; n<(unsigned long) point->n_edges ; n++) //sum over n
        {
            if(point->edges[n].edge!=edge->edge) //n different from k
            {
                if(point->edges[n].edge->connected)
                {
                    temp -= U15data(point->edges[n].edge, edge->edge); /// negative //UrawData(point->edges[n].edge, edge->edge, U15);
                }
            }
        }
        *partial_sum = -temp;
        return temp;
    }


    ///compute delta data and store it in *partial_sum
    *partial_sum = 0;

    //compute current state
    if(edge->edge->connected)
    {
        ///n=1 n+1=0
        for(unsigned long n=0 ; n<(unsigned long) point->n_edges ; n++) //sum over n
        {
            if(point->edges[n].edge!=edge->edge) //n different from k
            {
                if(point->edges[n].edge->connected)
                {
                    *partial_sum -= U15data(point->edges[n].edge, edge->edge); /// negative //UrawData(point->edges[n].edge, edge->edge, U15);
                }
            }
        }
    }
    else
    {
        ///n=0 n+1=1
        for(unsigned long n=0 ; n<(unsigned long) point->n_edges ; n++) //sum over n
        {
            if(point->edges[n].edge!=edge->edge) //n different from k
            {
                if(point->edges[n].edge->connected)
                {
                    *partial_sum += U15data(point->edges[n].edge, edge->edge); /// positive //UrawData(point->edges[n].edge, edge->edge, U15);
                }
            }
        }
    }


    ///delta energy
    double delta=0;
    double Zcurrent, Znext;
    Zcurrent = point->degree * (point->degree-1) / 2;
    if(edge->connected)
    {
        Znext = (point->degree-1) * (point->degree-2) / 2;
    }
    else
    {
        Znext = (point->degree+1) * (point->degree-1+1) / 2;
    }

    delta = ((1.0/Znext)-(1.0/Zcurrent)) * point->error_sum + (*partial_sum/Znext);


    if(edge->connected)
    {
        *partial_sum = -*partial_sum;
    }
    return delta;
}
double C_energy::deltaU15topo(Point *point, Vector * edge)
{
    if(edge->edge->connected)
    {
        int degree=point->degree; //degree in connected state //number of edges coming out from the vertex point
        double res = 0, resDegree = 0;
        if(degree==3)
        {
            //find both remaining edges
            Edge* e1 = NULL;
            Edge* e2 = NULL;
            for(unsigned long i=0 ; i<(unsigned long)point->n_edges ; i++)
            {
                if(point->edges[i].connected)
                {
                    if(point->edges[i].edge != edge->edge)
                    {
                        if(e1==NULL)
                        {
                            e1 = point->edges[i].edge;
                        }
                        else
                        {
                            e2 = point->edges[i].edge;
                            i = point->n_edges;
                        }
                    }
                }
            }
            if(e1!=NULL && e2!=NULL)
            {
                res = f3(angleFunc2(e1, e2), m_t0, m_sigma);
            }
            else
            {
                //throw an error
                cout << "throw an error because degree didn't match the actual number of edges: degree=3" << endl;
                throw -1;
            }
        }

        if(degree==2)
        {
            //find both remaining edges
            Edge* e1 = NULL;
            for(unsigned long i=0 ; i<(unsigned long) point->n_edges ; i++)
            {
                if(point->edges[i].connected)
                {
                    if(point->edges[i].edge != edge->edge)
                    {
                        e1 = point->edges[i].edge;
                        i = point->n_edges;
                    }
                }
            }
            if(e1!=NULL)
            {
                res = -f3(angleFunc2(e1, edge->edge), m_t0, m_sigma);
            }
            else
            {
                //throw an error
                cout << "throw an error because degree didn't match the actual number of edges: degree=2" << endl;
                throw -1;
            }
        }

        resDegree = phi_2(degree-1) - phi_2(degree);

        ///return separate or not
        return m_alpha*resDegree + m_beta*res;
    }
    else
    {
        int degree=point->degree; //degree in connected state //number of edges coming out from the vertex point
        double res = 0, resDegree = 0;
        if(degree==1)
        {
            Edge* e1 = NULL;
            for(unsigned long i=0 ; i<(unsigned long) point->n_edges ; i++)
            {
                if(point->edges[i].connected)
                {
                    e1 = point->edges[i].edge;
                    i = point->n_edges;
                }
            }
            if(e1!=NULL)
            {
                res = f3(angleFunc2(e1, edge->edge), m_t0, m_sigma);
            }
            else
            {
                //throw an error
                cout << "else throw an error because degree didn't match the actual number of edges: degree=1" << endl;
                throw -1;
            }
        }
        if(degree==2)
        {
            Edge* e1 = NULL;
            Edge* e2 = NULL;
            for(unsigned long i=0 ; i<(unsigned long) point->n_edges ; i++)
            {
                if(point->edges[i].connected)
                {
                    if(e1==NULL)
                    {
                        e1 = point->edges[i].edge;
                    }
                    else
                    {
                        e2 = point->edges[i].edge;
                        i = point->n_edges;
                    }
                }
            }
            if(e1!=NULL && e2!=NULL)
            {
                res = -f3(angleFunc2(e1, e2), m_t0, m_sigma);
            }
            else
            {
                //throw an error
                cout << "else throw an error because degree didn't match the actual number of edges: degree=2" << endl;
                throw -1;
            }
        }
        resDegree = phi_2(degree+1) - phi_2(degree);

        ///return separate or not
        return m_alpha*resDegree + m_beta*res;
    }
}













double C_energy::U16data(Edge* e1, Edge* e2)
{
    return U15data(e1, e2);
}
double C_energy::U16topo(Point* p)
{
    double qwer;
    if(p->n_edges<EDGES_PER_NODE)
    {
        //modif m_phi1 and m_beta
        double m_phi1_temp = m_phi1;
        m_phi1 = 0.0;
        m_beta *= 2.0;
        //compute delta
        qwer = U15topo(p);
        //reset m_phi1 and m_beta
        m_phi1 = m_phi1_temp;
        m_beta /= 2.0;
    }
    else
    {
        //compute delta
        qwer = U15topo(p);
    }
    return qwer;
}
double C_energy::deltaEdgeU16(Edge* e, double* partial_sum_p1, double* partial_sum_p2)
{
    return deltaU16(e->p1, e->this_in_p1, partial_sum_p1) + deltaU16(e->p2, e->this_in_p2, partial_sum_p2);
}
double C_energy::deltaU16(Point *point, Vector * edge, double * partial_sum)
{
    return deltaU16data(point, edge, partial_sum) + deltaU16topo(point, edge);
}
double C_energy::deltaU16data(Point *point, Vector * edge, double * partial_sum)
{
    return deltaU15data(point, edge, partial_sum);
}
double C_energy::deltaU16topo(Point *point, Vector * edge)
{
    double qwer;
    if(point->n_edges<EDGES_PER_NODE)
    {
        //modif m_phi1 and m_beta
        double m_phi1_temp = m_phi1;
        //m_beta0 = m_beta;
        m_beta *= 2.0;
        m_phi1 = 0.0;
        //compute delta
        qwer = deltaU15topo(point, edge);
        //reset m_phi1 and m_beta
        m_phi1 = m_phi1_temp;
        m_beta /= 2.0;
    }
    else
    {
        //compute delta
        qwer = deltaU15topo(point, edge);
    }
    return qwer;
}












///tools
double C_energy::angleFunc(Edge* e1, Edge* e2)
{
    Point *p1, *p2, *p12;

    if(e1->p1==e2->p1)
    {
        //
        p1 = e1->p2;
        p12 = e1->p1;
        p2 = e2->p2;
    }
    else if(e1->p1==e2->p2)
    {
        //
        p1 = e1->p2;
        p12 = e1->p1;
        p2 = e2->p1;
    }
    else if(e1->p2==e2->p1)
    {
        //
        p1 = e1->p1;
        p12 = e1->p2;
        p2 = e2->p2;
    }
    else if(e1->p2==e2->p2)
    {
        //
        p1 = e1->p1;
        p12 = e1->p2;
        p2 = e2->p1;
    }
    else
    {
        cout << "will throw an exception -1.0" << endl;
        throw -1.0;
    }

    return angleFunc(p1, p2, p12);
}
double C_energy::angleFunc(Point* p1, Point* p2, Point* p12)
{
    double du_x, du_y, du_z, dv_x, dv_y, dv_z;
    du_x = p1->x - p12->x;
    du_y = p1->y - p12->y;
    du_z = p1->z - p12->z;
    dv_x = p2->x - p12->x;
    dv_y = p2->y - p12->y;
    dv_z = p2->z - p12->z;
    return (du_x*dv_x + du_y*dv_y + du_z*dv_z) / (  sqrt( SQR(du_x) + SQR(du_y) + SQR(du_z) )  *   sqrt( SQR(dv_x) + SQR(dv_y) + SQR(dv_z) )  );
}
double C_energy::angleFunc2(Edge* e1, Edge* e2)
{
    double COS = angleFunc(e1, e2);

    if(COS>=1)
    {
        return 0.0;
    }

    if(COS<=-1)
    {
        return Pi;
    }

    return acos(COS);
}


double C_energy::U1dataComputation(Point * point, Vector * edge1, Vector * edge2)
{
    if (edge1->dx*edge2->dx+edge1->dy*edge2->dy+edge1->dz*edge2->dz <-1e-12) //do not tolerate right angles
    {
        double ex,ey,ez;
        ex= (double) edge1->dx - (double) edge2->dx;
        ey= (double) edge1->dy - (double) edge2->dy;
        ez= (double) edge1->dz - (double) edge2->dz;
        double rvalue = 0.0;
        ///if DTI data
        if(m_rawData->DimT==6 && m_rawData->numDim==4)
        {
            ///find current pixel
            signed long i = (signed long) floor((point->x/(m_rawData->pixDimX)));
            signed long j = (signed long) floor((point->y/(m_rawData->pixDimY)));
            signed long k = (signed long) floor((point->z/(m_rawData->pixDimZ)));
            double* DTI = new double[6];
            DTI[0] = m_rawData->raw4D[i][j][k][0];
            DTI[1] = m_rawData->raw4D[i][j][k][1];
            DTI[2] = m_rawData->raw4D[i][j][k][2];
            DTI[3] = m_rawData->raw4D[i][j][k][3];
            DTI[4] = m_rawData->raw4D[i][j][k][4];
            DTI[5] = m_rawData->raw4D[i][j][k][5];
            rvalue = 1. - sqrt((SQR( ((double) DTI[0]) *ex+ ((double) DTI[1])*ey+ ((double) DTI[2])*ez)+\
                SQR(((double) DTI[1])*ex+ ((double) DTI[3])*ey+ ((double) DTI[4])*ez)+\
                SQR(((double) DTI[2])*ex+ ((double) DTI[4])*ey+ ((double) DTI[5])*ez))/
                (SQR(ex)+SQR(ey)+SQR(ez)));
            delete DTI;
        }
        else if(m_rawData->DimT>6 && Hpsi!=NULL && b!=0.0 && m_rawData->numDim==4) ///if DWI data
        {
            ///find current pixel
            signed long i = (signed long) floor((point->x/(m_rawData->pixDimX)));
            signed long j = (signed long) floor((point->y/(m_rawData->pixDimY)));
            signed long k = (signed long) floor((point->z/(m_rawData->pixDimZ)));
            C_tensorMaker* toolTensor = new C_tensorMaker();
            double* DTI = toolTensor->makeTensor(&(m_rawData->raw4D[i][j][k][1]), Hpsi, b, m_rawData->raw4D[i][j][k][0], m_rawData->DimT-1);
            rvalue = 1. - sqrt((SQR( ((double) DTI[0]) *ex+ ((double) DTI[1])*ey+ ((double) DTI[2])*ez)+\
                SQR(((double) DTI[1])*ex+ ((double) DTI[3])*ey+ ((double) DTI[4])*ez)+\
                SQR(((double) DTI[2])*ex+ ((double) DTI[4])*ey+ ((double) DTI[5])*ez))/
                (SQR(ex)+SQR(ey)+SQR(ez)));
            delete DTI;
            delete toolTensor;
        }
        else
        {
            throw NAN;
        }
        if (rvalue > SMALL_NUM)
                return cbrtf(rvalue);
        else return 0.;
    }
    else return 1.;
}

double C_energy::U15dataComputation(Point* p1, Point* p2)//, unsigned long METHOD_INTERPOLATION, unsigned long OBJECT_INTERPOLATED, unsigned long METHOD_INTEGRATION
{
    C_toolbox_interpolate* interpolateTool = new C_toolbox_interpolate(m_rawData, m_rawDataMask, m_R/**physical distance*//**1.5*//**double r, should depend on lattice parameters?*/);

    //in case of thread, there's a need for several tools
    C_point_array_data* ctrlPointThread = new C_point_array_data(m_nbInterpolatingPoint, 3, 3, false);
    interpolateTool->regularInterpolatePoint(ctrlPointThread, p1, p2);

    //init container: interpolated data
    double dX = p2->x - p1->x;
    double dY = p2->y - p1->y;
    double dZ = p2->z - p1->z;
    double NORM = sqrt( SQR(dX) + SQR(dY) + SQR(dZ) ); //integrale over s\in[0,1] of ||phiprime(s)||


    //fill interpolated data
    interpolateTool->run(ctrlPointThread, METHOD_INTERPOLATION, OBJECT_INTERPOLATED, false, Hpsi, b);
//    cout << "after run interpolateTool" << endl;
    double r;
    if(m_nbInterpolatingPoint==1)
    {
        r = ctrlPointThread->FA[0]*ABS(ctrlPointThread->vectorInterpolated[0][0]*dX+\
            ctrlPointThread->vectorInterpolated[0][1]*dY+\
            ctrlPointThread->vectorInterpolated[0][2]*dZ)/NORM;
    }
    else
    {
        //integrate over the 3D
        double* ctrlPoints = new double[m_nbInterpolatingPoint];
        for(unsigned long i=0 ; i<ctrlPointThread->nb_point ; i++)
        {
            ctrlPoints[i] = ctrlPointThread->FA[i]*ABS(ctrlPointThread->vectorInterpolated[i][0]*dX+\
            ctrlPointThread->vectorInterpolated[i][1]*dY+\
            ctrlPointThread->vectorInterpolated[i][2]*dZ)/NORM;
        }
        C_toolbox_integrate* integrateToolThread = new C_toolbox_integrate();
        r = integrateToolThread->run(ctrlPoints, ctrlPointThread->nb_point, 1.0, METHOD_INTEGRATION);

        delete integrateToolThread;
        delete ctrlPoints;
    }
//    cout << "at the end" << endl;
    delete ctrlPointThread;
    delete interpolateTool;
//    if(isnan(r)) cout << "NaN in energy init" << endl;
//    if(isinf(r)) cout << "Infty in energy init" << endl;
    return r;
}

double C_energy::U16dataComputation(Point* p1, Point* p2)
{
    return U15dataComputation(p1, p2);
}

double C_energy::phi(int n)
{
    switch(n)
    {
        case 0:
            return m_phi0;
        break;
        case 1:
            return m_phi1;
        break;
        case 2:
            return 0.0;
        break;
        default:
            return (double) n -1.0;
        break;
    }
}

double C_energy::phi_2(int n)
{
    switch(n)
    {
        case 0:
            return m_phi0;
        break;
        case 1:
            return m_phi1;
        break;
        case 2:
            return 0.0;
        break;
        default:
            return (double) n;
        break;
    }
}

double C_energy::f(double t, double gamma)
{
    return pow(0.5 * (1 + t), gamma);
}

double C_energy::f2(double t, double t0, double sigma)
{
    return 1.0/(1.0+exp(-(t-t0)/sigma));
}
double C_energy::f3(double t, double t0, double sigma)
{
    return 1.0/(1.0+exp((t-t0)/sigma));
}

void C_energy::computeBetaIncrease(double T, double Tinit, double Tend) //beta increase with iteration
{
    if(T>Tinit) m_beta = 0.0;
    else if(T<Tend) m_beta = 0.1;//0.5;
    else m_beta = 0.1*log(T/Tinit)/log(Tend/Tinit); //0.5*log(T/Tinit)/log(Tend/Tinit);
    return;
}
void C_energy::computeBetaDecrease(double T, double Tinit, double Tend) //beta increase with iteration
{
    if(T>Tinit) m_beta = 0.5;
    else if(T<Tend) m_beta = 0.1;//0.0;
    else m_beta = 0.4*(1.0-(log(T/Tinit)/log(Tend/Tinit))) + 0.1;//0.5*(1.0-(log(T/Tinit)/log(Tend/Tinit)));
    return;
}

#ifdef MY_TEST_MULTITHREAD
double C_energy::initPairEdgeData(unsigned long kind, unsigned long nb_thread)
#else
double C_energy::initPairEdgeData(unsigned long kind)
#endif
{
    #ifdef MY_TEST_MULTITHREAD
    //unsigned long nb_thread=8;
    C_thread_init_energy** THREADS_U = new C_thread_init_energy*[nb_thread];
    for(unsigned long i=0 ; i<nb_thread ; i++)
    {
        unsigned long idxStart = (i*(m_graph->getLenPoints()))/nb_thread;
        unsigned long idxEnd = (((i+1)*(m_graph->getLenPoints()))/nb_thread);//-1;
        THREADS_U[i] = new C_thread_init_energy(this,idxStart, idxEnd, kind);
//        cout << "start thread: " << i << endl;
        THREADS_U[i]->start();
    }

    ///wait for all threads to finish
    for(unsigned long i=0 ; i<nb_thread ; i++)
    {
//        cout << "wait for thread: " << i << endl;
        THREADS_U[i]->wait_for_exit();
    }

    ///sum up energy
    double totU = 0.0;
    for(unsigned long i=0 ; i<nb_thread ; i++)
    {
        totU += THREADS_U[i]->subTotU;
//        cout << i << " " << totU  << endl;
    }

    ///delete all threads
    for(unsigned long i=0 ; i<nb_thread ; i++)
    {
        delete THREADS_U[i];
    }
    delete THREADS_U;
    return totU;
    #else
    return initPairEdgeDataMultiThread(kind, 0, m_graph->getLenPoints());
    #endif
}


double C_energy::initPairEdgeDataMultiThread(unsigned long kind, unsigned long idxStart, unsigned long idxEnd)
{
    if(kind==U15 || kind==U1 || kind==U18)
    {
        if(idxEnd>m_graph->getLenPoints()) idxEnd = m_graph->getLenPoints();
        if(idxEnd<idxStart || idxStart>=m_graph->getLenPoints() )
        {
            return 0.0;
        }
        Point* p = m_graph->getPoints();
        Point* p1;
        Point* p2;
        double totU = 0;
        for(unsigned long n=idxStart/**0*/ ; n<idxEnd/**m_graph->getLenPoints()*/ ; n++)
        {
            if(p[n].n_edges>1)
            {
                unsigned long idx;
                for(unsigned long ii=0 ; ii<(unsigned long) p[n].n_edges ; ii++)
                {
                    p1 = p[n].edges[ii].other_point;
                    for(unsigned long jj=ii+1 ; jj<(unsigned long) p[n].n_edges ; jj++)
                    {
                        p2 = p[n].edges[jj].other_point;
                        idx = (ii+1)*p[n].n_edges - (ii+1)*(ii+2)/2 + 1 - (p[n].n_edges - (jj+1)) - 2;
                        if(kind==U15 || kind==U18)
                        {
                            if(p1!=NULL && p2!=NULL) p[n].data[idx] = U15dataComputation(p1, p2);
                        }
                        else //U1
                        {
                            if(p1!=NULL && p2!=NULL) p[n].data[idx] = U1dataComputation(&(p[n]), &(p[n].edges[ii]), &(p[n].edges[jj]));
                        }
                        if(p[n].edges[jj].connected && p[n].edges[ii].connected) totU -= p[n].data[idx];
                    }
                }
            }
        }
        return totU;
    }
    else
    {
        return 0.0;
    }
}


bool C_energy::initThreads()
{
    deleteThreads();
    THREAD_DELTA = new C_thread_deltaU_move_vertex*[EDGES_PER_NODE]; //C_thread_deltaU_move_vertex**
    if(THREAD_DELTA==NULL) return false;
    for(int i=0 ; i<EDGES_PER_NODE ; i++)
    {
        THREAD_DELTA[i] = new C_thread_deltaU_move_vertex(this);
        if(THREAD_DELTA[i]==NULL)
        {
            deleteThreads();
            return false;
        }
    }
    return true;
}
void C_energy::deleteThreads()
{
    if(THREAD_DELTA!=NULL)
    {
        for(int i=0 ; i<EDGES_PER_NODE ; i++)
        {
            if(THREAD_DELTA[i]!=NULL)
            {
                delete THREAD_DELTA[i];
                THREAD_DELTA[i] = NULL;
            }
        }
        delete THREAD_DELTA;
        THREAD_DELTA = NULL;
    }
}

double C_energy::deltaVertexMoveU(Point* V /**point of the graph*/, Point* newV /**not in graph*/, long kind, vector<long>* idxPoint, vector<vector<long> >* idxPair, vector<vector<double> >* dataPair)
{
    idxPoint->clear();
    idxPair->clear();
    dataPair->clear();
    Point* V0 = m_graph->getPoints();

    double dU = 0.0;

    if(m_nbInterpolatingPoint<=3)
    {
        ///
        *idxPoint = vector<long>(V->n_edges,0);
        *idxPair = vector< vector<long> >(V->n_edges, vector<long>(0,0));
        *dataPair = vector< vector<double> >(V->n_edges, vector<double>(0.0,0));
        for(int n=0 ; n<V->n_edges ; n++)
        {
            idxPoint->at(n) = V->edges[n].other_point-V0;
            Point* u = V->edges[n].other_point;
            Point* v = V;
            Point* v_new = newV;
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
            double dU_data = 0.0;
        //    dU_topo = 0.0;
            vector<long> idxEdgeE1;
            for(long i=0 ; i<u->n_edges ; i++) ///go through all edges of u
            {
                if(i!=idxEdge2)
                {
                    ///store pair idx
                    idxPair->at(n).push_back(m_graph->getIndexOfPair(u, i, idxEdge2));
                    idxEdgeE1.push_back(i);

                    ///calculate and store data fidelity with U15dataComputation(v_new, other_point)
                    if(kind==U15 || kind==U18)
                    {
                        dataPair->at(n).push_back(U15dataComputation(u->edges[i].other_point,v_new)); //should be the time consuming part
                    }
                    else
                    {
                        dataPair->at(n).push_back(U16dataComputation(u->edges[i].other_point,v_new)); //should be the time consuming part
                    }
                }
            }

            ///if at least two edges are connected, calculate energy variation
            if(u->degree>1)
            {
                ///sumup with normalization
                if(e2->connected)
                {
                    for(unsigned int i=0 ; i<dataPair->at(n).size() ; i++)
                    {
                        //if the other edge is connected, add it up
                        if(u->edges[idxEdgeE1.at(i)].connected)
                        {
                            dU_data -= (dataPair->at(n).at(i)-u->data[idxPair->at(n).at(i)]);
                        }
                    }
                    dU_data /= (u->degree*(u->degree-1))/2;
                }
            }
            dU += dU_data;
        }
    }
    else
    {
        for(int i=0 ; i<V->n_edges ; i++)
        {
            THREAD_DELTA[i]->resetAttributes();
            ///set attributes of thread
            THREAD_DELTA[i]->u = V->edges[i].other_point;
            THREAD_DELTA[i]->v = V;
            THREAD_DELTA[i]->v_new = newV;
            THREAD_DELTA[i]->m_kind = kind;
            ///start thread
            THREAD_DELTA[i]->start();
        }

        ///wait for all thread to finish job
        for(int i=0 ; i<V->n_edges ; i++)
        {
            THREAD_DELTA[i]->wait_for_exit();
        }

        ///retreive data
        *idxPoint = vector<long>(V->n_edges,0);
        *idxPair = vector< vector<long> >(V->n_edges, vector<long>(0,0));
        *dataPair = vector< vector<double> >(V->n_edges, vector<double>(0.0,0));
        for(int i=0 ; i<V->n_edges ; i++) ///for all point in Nv
        {
            ///calculate energy variation from the neighborhood of V
            dU += THREAD_DELTA[i]->dU_data;
            ///fill vector that maps point that may be modified
            idxPoint->at(i) = V->edges[i].other_point-V0;
            ///for the ith point, get information about pairs of edges htat may be modified
            idxPair->at(i) = THREAD_DELTA[i]->idxPair;
            dataPair->at(i) = THREAD_DELTA[i]->dataPair;
        }
    }


    ///compute topo
    double tempX = V->x;
    double tempY = V->y;
    double tempZ = V->z;
    for(int i=0 ; i<V->n_edges ; i++) ///for all point in Nv
    {
        ///if at least two edges are connected, calculate energy variation
        if(V->edges[i].other_point->degree==2)
        {
            ///sumup with normalization
            if(V->edges[i].edge->connected)
            {
                if(kind==U15 || kind==U18)
                {
                    dU -= U15topo(V->edges[i].other_point);
                }
                else
                {
                    dU -= U16topo(V->edges[i].other_point);
                }
                V->x = newV->x;
                V->y = newV->y;
                V->z = newV->z;
                if(kind==U15 || kind==U18)
                {
                    dU += U15topo(V->edges[i].other_point);
                }
                else
                {
                    dU += U16topo(V->edges[i].other_point);
                }
                V->x = tempX;
                V->y = tempY;
                V->z = tempZ;
            }
        }
    }

    //delta topo of the considered vertex
    if(V->degree==2) //topological term will be modified
    {
        //compute the difference between the possible and current one
        if(kind==U15 || kind==U18)
        {
            dU -= U15topo(V);
        }
        else
        {
            dU -= U16topo(V);
        }
        V->x = newV->x;
        V->y = newV->y;
        V->z = newV->z;
        if(kind==U15 || kind==U18)
        {
            dU += U15topo(V);
        }
        else
        {
            dU += U16topo(V);
        }
        V->x = tempX;
        V->y = tempY;
        V->z = tempZ;
    }
    return dU;
}
