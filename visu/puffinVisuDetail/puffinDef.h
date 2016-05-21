#ifndef PUFFIN_DEFS
#define PUFFIN_DEFS

#include <MACRO_DEF.h>
#include <struct.h>

#define M_TENSOR 1
#define M_EDGE   2
#define M_FIBERS 4
#define M_RAW_FIBERS 8

//#define LOWZ 0//20//15//15//7//16//30//20//8//58//1.38
//#define HIGHZ -10//30//18//9//20//17.5//0//31.4//21.4//0//21.4//23//21.4//-8//64//46.5
#define LOWR 0
#define HIGHR 200000

typedef float floating;
#define MAX_EDGES_PER_NODE 128
#define MIN_PHYSICAL_LENGTH 2 //60//0//0


struct TensorData { //This is to store the points of the graph in the points array
	double x; //Coordinates of the point
	double y;
	double z;
	//double matrix[6]; //DTI Tensor //maybe we can add eigenvalues of the matrix matrix can be computed from eigenvalues and eigen vectors
	/* symmetric matrix format:
	 * ( 0 1 2 )
	 * (   3 4 )
	 * (     5 ) */
	 double eigenValue[3];
	 double eigenVector1[3];
	 double eigenVector2[3];
	 double eigenVector3[3];
	 //rotation angle to draw the ellipse corresponding to the tensor (apply the rotation theta around z0, phi around y1 and psi around z2, where 0 1 and 2 correspond to the first coordinate system, the second and the third rep.)
	 double theta;
	 double phi;
	 double psi;

};




struct plotData { //This is to store the points of the graph in the points array
	TensorData* T;
    unsigned long nbTensor;
    CtlCtlStruct* Fib;
    int rawFibers;
    plotData()
    {
        T = (TensorData*) 0;
        nbTensor = 0;
        Fib = (CtlCtlStruct*) 0;
        rawFibers = 0;
    }
    ~plotData()
    {
        if(T!=(TensorData*) 0)
        {
            delete T;
        }
        if(Fib!=(CtlCtlStruct*) 0)
        {
            delete Fib;
        }
    }
};

#endif
