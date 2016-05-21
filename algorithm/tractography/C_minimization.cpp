#include "C_minimization.h"

extern "C"
{
C_minimization::C_minimization(C_graph* graph)
{
    //m_graph_filename = "";
    m_graph = graph;
}
C_minimization::~C_minimization()
{
    //dtor
}


void C_minimization::flip_edge(Edge * edge,double psum_p1,double psum_p2)
{
	int was_connected=edge->connected;
	edge->connected = !was_connected;
	edge->this_in_p1->connected = !was_connected;
	edge->this_in_p2->connected = !was_connected;

	if(true)
	{
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
	else
	{
	    if(was_connected)
	    {
	        edge->p1->degree--;
            edge->p2->degree--;
	    }
	    else
	    {
	        edge->p1->degree++;
            edge->p2->degree++;
	    }
		edge->p1->error_sum+=psum_p1;
		edge->p2->error_sum+=psum_p2;
	}
}
unsigned long* C_minimization::sweep_idx_edges(void)
{
	return sweep_idx(m_graph->getLenEdges());
}
unsigned long* C_minimization::sweep_idx_points(void)
{
	return sweep_idx(m_graph->getLenPoints());
}
unsigned long* C_minimization::sweep_idx(unsigned long N)
{
    unsigned long * index_array= new unsigned long[N];
	for (unsigned long i=0;i<N;i++)
		index_array[i]=i;
	shuffle(index_array,N);
	return index_array;
}
void C_minimization::shuffle(unsigned long *array, size_t n)
{
    if (n > 1) {
        size_t i;
		for (i = 0; i < n - 1; i++) {
		  size_t j = i + rand() / (RAND_MAX / (n - i) + 1);
		  unsigned long t = array[j];
		  array[j] = array[i];
		  array[i] = t;
		}
    }
}
bool C_minimization::saveFiber(std::string fileNameFib, vector< vector<unsigned long> > fibers)
{
    //file names
    string headerFile = fileNameFib + ".fibHDR";
    string sourceFile = fileNameFib + ".fibSRC";

    //filestreams
    fstream filestrSRC, filestrHDR;
    filestrHDR.open (headerFile.data(), ofstream::binary | ofstream::out); //should check if ok
    filestrSRC.open (sourceFile.data(), ofstream::binary | ofstream::out); //should check if ok

    //get number of bundles
    unsigned long NB = fibers.size();
    filestrHDR.write((char*) &NB, sizeof(unsigned long)); //write it in header file

    //declare a buffer to store and write
    BundlePoint* b_fibers;

    //get the graph point addresses
    Point* my_points = m_graph->getPoints();

    unsigned long maxN;
    //for each bundle
    for(unsigned long u=0 ; u<NB ; u++)
    {
        //get the size of the current bundle and save it in header file
        maxN = fibers.at(u).size();
        filestrHDR.write((char*) &maxN, sizeof(unsigned long));

        //allocate bundle
        b_fibers = new BundlePoint[maxN]; //should check if allocation doesn't fail

        //fill in the allocated bundle
        for(unsigned long r=0 ; r<maxN ; r++)
        {
            //fill fibers
            b_fibers[r].x = my_points[fibers.at(u).at(r)].x;
            b_fibers[r].y = my_points[fibers.at(u).at(r)].y;
            b_fibers[r].z = my_points[fibers.at(u).at(r)].z;
            b_fibers[r].idx = fibers.at(u).at(r);
        }

        //write data in source file
        filestrSRC.write((char*) b_fibers, maxN*sizeof(BundlePoint));
        delete b_fibers;
    }

    //close files
    filestrSRC.close();
    filestrHDR.close();


    return true;
}
void C_minimization::disconnectEdge(Edge* e)
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

void C_minimization::randInit(double alphaConnected)
{
    C_toolbox_rand2 *p = new C_toolbox_rand2(rand());

    Edge* e = m_graph->getEdges();
    unsigned long Ne = m_graph->getLenEdges();
//    cout << "start initialization of " << Ne << " edges" << endl;
    for(unsigned long i=0 ; i<Ne ; i++)
    {
        if(p->doub()<alphaConnected)
        {
            e[i].connected = true;
            e[i].this_in_p1->connected = true;
            e[i].this_in_p2->connected = true;
        }
        else
        {
            e[i].connected = false;
            e[i].this_in_p1->connected = false;
            e[i].this_in_p2->connected = false;
        }
    }
    //set the actual degree for each point
    unsigned long Np = m_graph->getLenPoints();
    Point* p3 = m_graph->getPoints();
//    cout << "start initialization of " << Np << " points" << endl;
    for(unsigned long i=0 ; i<Np ; i++)
    {
        p3[i].degree = 0;
        for(unsigned long n=0 ; n<(unsigned long) p3[i].n_edges ; n++)
        {
            if(p3[i].edges[n].connected)
            {
                p3[i].degree++;
            }
        }
    }
//    cout << "rand init done" << endl;

    delete p;//->~C_toolbox_rand2();
    return;
}

}
