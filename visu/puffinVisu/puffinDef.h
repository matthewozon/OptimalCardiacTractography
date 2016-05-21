#ifndef PUFFIN_DEFS
#define PUFFIN_DEFS

#include <struct.h>
#include <vtkPolyData.h>
#include <vtkSmartPointer.h>

#define M_TENSOR 1
#define M_EDGE   2
#define M_FIBERS 4
#define M_LARGER_FIBER_AND_CONTROL_POINTS 8
#define M_LARGER_FIBER 16
#define M_RAW_FIBERS 32
#define M_BOTH_LARGEST_FIBER 64
#define M_SAVE 128
#define M_VTK_POLYDATA 256
#define M_HDR_READER 512

#define LOWZ 58//8//58
#define HIGHZ 64//-8//64
#define LOWR 0
#define HIGHR 200000

typedef float floating;
#define MAX_EDGES_PER_NODE 128

#define SQR(x) ((x)*(x))

//#define LATER

typedef struct _CtlStruct{ //this struct is to send data back to python
    unsigned long N_element; //means there are actualy 3*N_element GLfloat (because of x, y, z)
	float* elts;
	unsigned long* idx;
} CtlStruct_;

typedef struct __CtlCtlStruct{ //this struct is to send data back to python
    unsigned long N;
        CtlStruct_* fibers;
} CtlCtlStruct_;


//typedef struct _BundlePoint{ //this struct is to send data back to python
//	double x;
//	double y;
//	double z;
//	unsigned long idx;
//} BundlePoint;

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

//typedef struct _PointWithNeighbors { //This structure is used to read the graph
//	floating x;
//	floating y;
//	floating z;
//	int edges[MAX_EDGES_PER_NODE]; //Array of ids connected neighbors, negative if disconnected,0 for no edge //an array containing the id of each neighboor point connected to this point
//	int edge_ids[MAX_EDGES_PER_NODE]; //Array of edge ids, for easier interpretation of the graph // an array containing the id of each edge coming out from this point
//	floating matrix[6];
//	/* symmetric matrix format:
//	 * ( 0 1 2 )
//	 * (   3 4 )
//	 * (     5 ) */
//} PointWithNeighbors;



#ifdef LATER
typedef struct _SaveEdges{ //this struct is to send data back to python
	floating x1;
	floating y1;
	floating z1;
	floating x2;
	floating y2;
	floating z2;
	int connected;
} SaveEdges;
#endif

struct plotData { //This is to store the points of the graph in the points array
	TensorData* T;
    unsigned long nbTensor;
    CtlCtlStruct* Fib;
    int rawFibers;
    vtkSmartPointer<vtkPolyData> polyDataV;
    #ifdef LATER
        SaveEdges* E;
        unsigned long nbEdges;
    #endif
    plotData()
    {
        T = (TensorData*) 0;
        nbTensor = 0;
        Fib = (CtlCtlStruct*) 0;
        rawFibers = 0;
        polyDataV = 0;
        #ifdef LATER
            E = (SaveEdges*) 0;
            nbEdges = 0;
        #endif
    }
    ~plotData()
    {
        if(T!=(TensorData*) 0)
        {
            delete T;
        }
        if(Fib!=(CtlCtlStruct*) 0)
        {
            if(Fib->fibers!=(CtlStruct*) 0)
            {
                for(unsigned long i=0 ; i<Fib->N ; i++)
                {
                    if(Fib->fibers[i].elts!=(double*) 0)
                    {
                        delete Fib->fibers[i].elts;
                    }
                    if(Fib->fibers[i].idx!=(unsigned long*) 0)
                    {
                        delete Fib->fibers[i].idx;
                    }
                }
                delete Fib->fibers;
            }
            delete Fib;
        }
        #ifdef LATER
            if(E!=(SaveEdges*) 0)
            {
                delete E;
            }
        #endif
    }
};



//    ANALYZE
//TM
// Header File Format
//*
//*  (c) Copyright, 1986-1995
//*  Biomedical Imaging Resource
//*  Mayo Foundation
//*
//*  dbh.h
//*
//*  databse sub-definitions
//*/
//struct header_key   /* header key    */
//{                                 /* off + size */
//      int sizeof_hdr;   /* 0 + 4            */ //byte size of header
//      char data_type[10];  /* 4 + 10           */
//      char db_name[18];  /* 14 + 18          */
//      int extents;   /* 32 + 4           */ //should be 16384
//      short int session_error;  /* 36 + 2           */
//      char regular;    /* 38 + 1           */  //must be r to indicated that all images/volumes are of the same size
//      char hkey_un0;                   /* 39 + 1   */
//};                               /* total=40 bytes  */

//struct image_dimension
//{                                  /* off + size       */
//      short int dim[8];                 /* 0 + 16           */
//                        /** dim[0] = number of dim in db (usually 4)
//                            dim[1] =  Image X dimension; number of pixels in an image row
//                            dim[2]     Image Y dimension; number of pixel rows in slice
//                            dim[3]     Volume Z dimension; number of slices in a volume
//                            dim[4]     Time points, number of volumes in database
//                        */
//      short int unused8;                /* 16 + 2           */
//      short int unused9;                /* 18 + 2           */
//      short int unused10;               /* 20 + 2           */
//      short int unused11;               /* 22 + 2           */
//      short int unused12;               /* 24 + 2           */
//      short int unused13;               /* 26 + 2           */
//      short int unused14;               /* 28 + 2           */
//      short int datatype;               /* 30 + 2           */ //datatype for this image set (ex: DT_FLOAT)
//      short int bitpix;                 /* 32 + 2           */ //number of bits per pixel; 1, 8, 16, 32, or 64.
//      short int dim_un0;                /* 34 + 2           */ //unused
//      float pixdim[8];                  /* 36 + 32          */
//                      /**
//                           pixdim[] specifies the voxel dimensitons: Parallel array to dim[], giving real world measurements in mm. and ms.
//                           pixdim[1] - voxel width in mm.
//                           pixdim[2] - voxel height in mm.
//                           pixdim[3] - interslice distance / slice thickness in mm.
//                               ...etc
//                      */
//      float vox_offset;                 /* 68 + 4            */ //byte offset in the .img file at which voxels start. This value can negative to specify that the absolute value is applied for every image in the file.
//      float funused1;                   /* 72 + 4            */
//      float funused2;                   /* 76 + 4            */
//      float funused3;                   /* 80 + 4            */
//      float cal_max;                    /* 84 + 4            */ //specify the range of calibration values
//      float cal_min;                    /* 88 + 4            */ //specify the range of calibration values
//      float compressed;                 /* 92 + 4            */
//      float verified;                   /* 96 + 4            */
//      int glmax,glmin;                  /* 100 + 8           */ //The maximum and minimum pixel values for the entire database.
//};                                 /* total=108 bytes  */

//struct data_history
//{                                  /* off + size       */
//      char descrip[80];                 /* 0 + 80           */
//      char aux_file[24];                /* 80 + 24         */
//      char orient;                      /* 104 + 1          */
//      char originator[10];              /* 105 + 10         */
//      char generated[10];              /* 115 + 10         */
//      char scannum[10];   /* 125 + 10        */
//      char patient_id[10];             /* 135 + 10         */
//      char exp_date[10];   /* 145 + 10        */
//      char exp_time[10];   /* 155 + 10        */
//      char hist_un0[3];    /* 165 + 3          */
//      int views;                         /* 168 + 4          */
//      int vols_added;                   /* 172 + 4          */
//      int start_field;                  /* 176 + 4          */
//      int field_skip;                   /* 180 + 4         */
//      int omax, omin;                   /* 184 + 8         */
//      int smax, smin;                   /* 192 + 8         */
//};

//struct dsr
//{
//      struct header_key hk;              /* 0 + 40            */
//      struct image_dimension dime;      /* 40 + 108          */
//      struct data_history hist;          /* 148 + 200         */
//};                                  /* total= 348 bytes
/* Acceptable values for datatype */
#define DT_NONE                    0
#define DT_UNKNOWN                0
#define DT_BINARY                  1
#define DT_UNSIGNED_CHAR          2
#define DT_SIGNED_SHORT           4
#define DT_SIGNED_INT              8
#define DT_FLOAT                   16
#define DT_COMPLEX                 32
#define DT_DOUBLE                  64
#define DT_RGB                     128
#define DT_ALL                     255


//#define DATA_READ
#ifdef DATA_READ
struct rawPoint
{
    double x;
    double y;
    double z;
    double M[6];
};

struct rawPointStruct
{
    rawPoint* rawP;
    unsigned long nbRawP;
};
#endif

#endif
