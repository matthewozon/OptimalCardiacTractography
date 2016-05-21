#ifndef C_VISU_DRAW_H
#define C_VISU_DRAW_H

#include <vtkPolyData.h>
#include <vtkPoints.h>
#include <vtkCellArray.h>
#include <vtkFloatArray.h>
#include <vtkDoubleArray.h>
#include <vtkPointData.h>

#include <vtkSmartPointer.h>
#include <vtkTubeFilter.h>

#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkCamera.h>
#include <vtkInteractorStyleTrackballCamera.h>

#include <vtkPolyDataAlgorithm.h>
#include <vtkVertexGlyphFilter.h>


#include <vtkKochanekSpline.h>
#include <vtkCardinalSpline.h>
#include <vtkParametricFunctionSource.h>
#include <vtkParametricSpline.h>
#include <vtkParametricEllipsoid.h>

#include <vtkTransform.h>
#include <vtkTransformPolyDataFilter.h>
#include <vtkProperty.h>

#include <puffinDef.h>
#include <C_toolbox_spline.h>




using namespace std;


class C_visu_draw
{
    public:
        /** Default constructor */
        C_visu_draw();
        /** Default destructor */
        virtual ~C_visu_draw();

        void draw(plotData* data, unsigned long nbSubDiv=1, double radius=0.2, float R=0, float G=0, float B=0);
        //void draw(vector< vtkSmartPointer<vtkPoints> > points, long nbSubDiv, float R, float G, float B, double radius);

    private:
        void draw(CtlCtlStruct* m_ctlCtlStruct, long nbSubDiv, float R, float G, float B, double radius);
        void draw(vector< vtkSmartPointer<vtkPolyDataMapper> > vectorMap, float R, float G, float B);

        void computeSpline(vtkSmartPointer<vtkPoints> points, long nbSubDiv, double radius, vtkSmartPointer<vtkTubeFilter>);

        vtkSmartPointer<vtkPolyDataMapper> computeMapper(vtkSmartPointer<vtkTubeFilter> tube);

        vector< vtkSmartPointer<vtkPolyDataMapper> > run(vector< vtkSmartPointer<vtkPoints> > P, long nbSubDiv, double radius);
        vector< vtkSmartPointer<vtkPolyDataMapper> > run(CtlCtlStruct* m_ctlCtlStruct, long nbSubDiv, double radius);

        void computeColor(vtkSmartPointer<vtkPolyData> polydata);

        //draw
        vector< vtkSmartPointer< vtkActor > > drawTensors(TensorData* T, unsigned long nbTensor);
        vtkSmartPointer< vtkActor > drawTensor(TensorData* T);

        vector< vtkSmartPointer< vtkActor > > drawFibers(CtlCtlStruct* Fib, unsigned long nbSubDiv, double radius);
        vector< vtkSmartPointer< vtkActor > > drawRawFibers(CtlCtlStruct* Fib, double radius);
        vtkSmartPointer< vtkActor > drawFiber(CtlStruct* Fib, double radius);


        double fibLength(CtlStruct* fibers);

};

#endif // C_VISU_DRAW_H
