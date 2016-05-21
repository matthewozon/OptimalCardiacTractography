#ifndef C_VISU_IO_H
#define C_VISU_IO_H

#include <puffinDef.h>
#ifdef LATER
    #include <vector>
#endif
#include <C_toolbox_eigen_sym.h>
//#include <C_analyze.h>
#include <C_ReadAnalyze75.h>


#include <fstream>

#include <iostream>

#include <vtkSmartPointer.h>
#include <vtkPolyData.h>
#include <vtkPoints.h>
#include <vtkCellArray.h>
#include <vtkPolyDataWriter.h>
#include <vtkPolyDataReader.h>
using namespace std;

struct sliceIdx
{
    unsigned long fibIdx;
    unsigned long nbSubFib;
    vector<unsigned long> idxStart;
    vector<unsigned long> idxStop;
};

#define EXTRACT_DWI
#undef EXTRACT_DWI

class C_visu_IO
{
    public:
        /** Default constructor */
        C_visu_IO();
        /** Default destructor */
        virtual ~C_visu_IO();

        //calling func
        bool load(string fileName, long kind);

        void save(string filename);

        //attributes
        plotData* data;
        double Zmax;
        double Zmin;
//        TensorData* m_t;
//        unsigned long m_nb_tensor;
//        CtlCtlStruct* m_ctlCtlStruct;
//        #ifdef LATER
//            SaveEdges* m_e;
//            unsigned long m_nb_edges;
//        #endif

    private:

        //reading func
        bool readFiber(string fileNameFib);
        bool readLargestFiber(string fileNameFib);
        bool readVTKfibers(string fileNameFib);
        bool loadTensorFromPointData(string fileNamePoints);
        #ifdef LATER
            bool loadEdge(string fileNameEdge);
        #endif
        void analyzeInfo(const char* filenameHDR);

        //tools
        bool getTensorsFromRawPoints(PointWithNeighbors** p, unsigned long nb_point, TensorData** t); //compute the tensor data from raw points
        bool getTensorFromRawPoint(PointWithNeighbors* p, TensorData* t); //compute a tensor data from a raw point address
        string getBaseName(string fieName);

        #ifdef EXTRACT_DWI
        double getDWIfromDTI(PointWithNeighbors* p, double gx, double gy, double gz, double b);
        #endif
};

#endif // C_VISU_IO_H
