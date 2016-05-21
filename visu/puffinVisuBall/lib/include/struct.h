#ifndef STRUCT_H_INCLUDED
#define STRUCT_H_INCLUDED

//#include <ludcmp.h>
#include <math.h>
#include <stdlib.h>
#include <vector>
#include <settings_def.h>


typedef float floating;

struct Vector { //This is used to represent neighboring edges to a point
	struct Edge* edge; //Link to the edge in the edge array //should be protected //address of one edge
    double dx; //Vector coordinates of the edge for metric evaluation
	double dy;
	double dz;
	char connected; //no bool type in C
	struct Point * other_point; //This is there for faster neighborhood operations, useful for icm_v4 //should be protected
};


struct Point { //This is to store the points of the graph in the points array
	double x; //Coordinates of the point
	double y;
	double z;
	int degree; //Cached for faster access // d(v) = sum(we)over e belonging to the set of adjacent edges, v is a given vertex (point) and we weight of edge e (0 or 1) (unsigned long?)
	int n_edges; //(unsigned long?)
	double error_sum; //tangent error
	struct Vector * edges; //Neighboring edges to the point //should be protected
	double* data;
	#ifdef _OLD_GRAPH
	double matrix[6]; //DTI Tensor //maybe we can add eigenvalues of the matrix
	#endif
	#ifdef HANDLE_VERTEX_MOVE
	double x0;
	double y0;
	double z0;
	#endif
	//floating matrix[6]; //DTI Tensor //maybe we can add eigenvalues of the matrix
	/* symmetric matrix format:
	 * ( 0 1 2 )
	 * (   3 4 )
	 * (     5 ) */
};

struct Edge {
	struct Point * p1; //should be protected
	struct Point * p2; //should be protected
	struct Vector * this_in_p1; //Pointer to the edge in the point's neighbor list //should be protected
	struct Vector * this_in_p2; //should be protected
	char connected; //no bool type in C
};



typedef struct _PointWithNeighbors { //This structure is used to read the graph
	floating x;
	floating y;
	floating z;
	int edges[MAX_EDGES_PER_NODE]; //global index in edge array //Array of ids connected neighbors, negative if disconnected,0 for no edge //an array containing the id of each neighboor point connected to this point
	int edge_ids[MAX_EDGES_PER_NODE]; //local index in vector //Array of edge ids, for easier interpretation of the graph // an array containing the id of each edge coming out from this point
	floating matrix[6];
	/* symmetric matrix format:
	 * ( 0 1 2 )
	 * (   3 4 )
	 * (     5 ) */
} PointWithNeighbors;


typedef struct _SaveEdges{
	floating x1;
	floating y1;
	floating z1;
	floating x2;
	floating y2;
	floating z2;
	int connected;
} SaveEdges;


struct diffInfo
{
    double g[3]; //unit gradient direction expressed in patient basis
    double echoTime; //echo time
    double repTime;  //repetition time
    double b; //b-value
    double X_slice[3]; //first vector of patient basis expressed in machine basis
    double Y_slice[3]; //second vector...
    double Z_slice[3]; //third  vector...
    double voxelSize[3]; //size of voxels (X dim, Ydim, Z dim)
};

struct samplingFactor ///down or up sampling factors: set to 1 if no sampling to do
{
    ///sampling factors
    double Xfactor;
    double Yfactor;
    double Zfactor;
    double Tfactor;

    ///dimension of sampling ball (the "ball" that defines the elements to take into account for interpolation)
    unsigned long Rx;
    unsigned long Ry;
    unsigned long Rz;
    unsigned long Rt;

    ///variance of gaussian
    double sigx;
    double sigy;
    double sigz;
    double sigt;

    samplingFactor()
    {
        Xfactor = 1.0;
        Yfactor = 1.0;
        Zfactor = 1.0;
        Tfactor = 1.0;
        Rx = 0.0;
        Ry = 0.0;
        Rz = 0.0;
        Rt = 0.0;
        sigx = 0.3;
        sigy = 0.3;
        sigz = 0.3;
        sigt = 0.3;
    };
    ~samplingFactor()
    {
        //
    };
};

/* Acceptable values for datatype defined in settings_def.h*/

typedef struct _BundlePoint{
    double x;
	double y;
	double z;
	unsigned long idx; //corresponding index in graph (equivalent to the one in .cformat file)
} BundlePoint;

struct Mixture{
    //unsigned long nb_mix; //number of mixture
    //double* alpha; //proportion of each mixture. sum alpha = 1
    double alpha1;
    double alpha2;
};

struct MovingEdge{
    unsigned long nb_edge; //number of edge that can be moved
    Edge** e;//[MAX_EDGE_MOVED_AT_ONCE]; //array containing the addresses of the concerned edges
    Point* p;//not in graph
    unsigned long idxPoint;
};

struct LandScape{
    Edge*** e; //array of sizeOfLandScape line and dim raw
    unsigned long sizeOfLandScape; //tell how many solution are in array e
    unsigned long dim; //tell how many edge are used to do a new solution
    double data1;
    double data2;
    //IDX* idx; //idx[i].idx is an array containing all address where edge i can be found
};


struct CtlStruct{ //this struct is to send data back to python
    unsigned long N_element; //means there are actualy 3*N_element GLfloat (because of x, y, z)
	double* elts;
	unsigned long* idx;
	CtlStruct()
	{
	    N_element = 0;
	    elts = NULL;
	    idx = NULL;
	}
	~CtlStruct()
	{
	    N_element = 0;
	    deleteArray();
	}
	void allocateArray(unsigned long N)
	{
	    deleteArray();
	    N_element = N;
	    elts = new double[3*N_element];
	    idx = new unsigned long[N_element];
	}

	private:
	void deleteArray(void)
	{
	    if(elts!=NULL)
	    {
	        delete elts;
	        elts = NULL;
	    }
	    if(idx!=NULL)
	    {
	        delete idx;
	        idx = NULL;
	    }
	}
};

struct CtlCtlStruct{ //this struct is to send data back to python
    unsigned long N;
	CtlStruct* fibers;
	CtlCtlStruct(unsigned long n)
	{
	    N = n;
	    fibers = NULL;
	    allocateCtrlStruct(N);
	}
	~CtlCtlStruct()
	{
	    N = 0;
	    deleteCtrlStruct();
	}

	private:
	void allocateCtrlStruct(unsigned long n)
	{
	    if(fibers!=NULL)
	    {
	        deleteCtrlStruct();
	    }
	    fibers = new CtlStruct[n];
	}
	void deleteCtrlStruct(void)
	{
	    for(unsigned long i=0 ; i<N ; i++)
	    {
	        fibers[i].~CtlStruct();// is it really necessary?
	    }
	    delete[] fibers; //may be not?
	    fibers = NULL;

	}
};

//typedef struct _CtlStruct{ //this struct is to send data back to python
//    unsigned long N_element; //means there are actualy 3*N_element GLfloat (because of x, y, z)
//	float* elts;
//	unsigned long* idx;
//} CtlStruct;
//
//typedef struct _CtlCtlStruct{ //this struct is to send data back to python
//    unsigned long N;
//	CtlStruct* fibers;
//} CtlCtlStruct;

#endif // STRUCT_H_INCLUDED
