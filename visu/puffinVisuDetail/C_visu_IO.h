#ifndef C_VISU_IO_H
#define C_VISU_IO_H

#include <puffinDef.h>
#include <C_toolbox_eigen_sym.h>

#include <fstream>

#include <C_measure.h>

#include <iostream>

//#include <vtkSmartPointer.h>
//#include <vtkPolyData.h>
//#include <vtkPoints.h>
//#include <vtkCellArray.h>
//#include <vtkPolyDataWriter.h>
//#include <vtkPolyDataReader.h>
using namespace std;

struct sliceIdx
{
    unsigned long fibIdx;
    unsigned long nbSubFib;
    vector<unsigned long> idxStart;
    vector<unsigned long> idxStop;
};

class C_visu_IO
{
    public:
        /** Default constructor */
        C_visu_IO();
        /** Default destructor */
        virtual ~C_visu_IO();

        //calling func
        bool load(string fileName, long kind);

        bool selectFiberDetail(double detailMin, double detailMax);

        bool saveFiber(std::string fileNameFib);

        //attributes
        plotData* data;
        double Zmax;
        double Zmin;

    private:

        //reading func
        bool readFiber(string fileNameFib);
        bool loadTensorFromPointData(string fileNamePoints);

        //tools
        bool getTensorsFromRawPoints(PointWithNeighbors** p, unsigned long nb_point, TensorData** t); //compute the tensor data from raw points
        bool getTensorFromRawPoint(PointWithNeighbors* p, TensorData* t); //compute a tensor data from a raw point address
        string getBaseName(string fieName);
};

#endif // C_VISU_IO_H
