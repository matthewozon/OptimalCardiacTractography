#include <iostream>
#include <C_visu_IO.h>
#include <C_visu_draw.h>

//#include <stdio.h>
//#include <string.h>
//#include <math.h>
//#include <stdlib.h>
//#include <time.h>
//#include <unistd.h>
//#include <sys/swap.h>
//#include <analyze.h>
//#include <substitutions.h>
//libtpcimgio
///usr/include/libtpcimgio

using namespace std;

///**for test*/
//#include <vtkParametricEllipsoid.h>
//#include <vtkSmartPointer.h>
//#include <vtkTransformPolyDataFilter.h>
//#include <vtkSphereSource.h>
//#include <vtkPolyData.h>
//#include <vtkTransform.h>
//#include <vtkPolyDataMapper.h>
//#include <vtkActor.h>
//#include <vtkRenderWindow.h>
//#include <vtkRenderer.h>
//#include <vtkRenderWindowInteractor.h>
//#include <vtkProperty.h>
////#include <math.h>
#include <vtkPolyDataWriter.h>

vtkSmartPointer<vtkActor> ellipsoid(double Xc, double Yc, double Zc, double l1, double l2, double l3, double theta, double phi, double psi, double vp1[3]);
void showTensor(TensorData* T);

int main(int argc, char **argv)
{

    if(argc>2)
    {
        C_visu_IO* ioPuffin = new C_visu_IO();
        bool canDraw = false;
        if(strcmp(argv[1],"TENSOR")==0 || atoi(argv[1])==1)
        {
            if(ioPuffin->load(argv[2], M_TENSOR))
            {
                canDraw = true;
            }
            else
            {
                canDraw = false;
            }
        }
        else if(strcmp(argv[1],"EDGE")==0 || atoi(argv[1])==2)
        {
            #ifdef LATER
                if(ioPuffin->load(argv[2], M_EDGE))
                {
                    canDraw = true;
                }
                else
                {
                    canDraw = false;
                }
            #else
                canDraw = false;
                cout << "Not yet available" << endl;
            #endif
        }
        else if(strcmp(argv[1],"FIBER")==0 || atoi(argv[1])==4)
        {
            ioPuffin->Zmax = 10;
            ioPuffin->Zmin = 1;
            if(ioPuffin->load(argv[2], M_FIBERS))
            {
                canDraw = true;
            }
            else
            {
                canDraw = false;
            }
        }
        else if(strcmp(argv[1],"LARGER_FIBER_AND_CONTROL_POINTS")==0 || atoi(argv[1])==8)
        {
            if(ioPuffin->load(argv[2], M_LARGER_FIBER_AND_CONTROL_POINTS))
            {
                canDraw = true;
            }
            else
            {
                canDraw = false;
            }
        }
        else if(strcmp(argv[1],"M_LARGER_FIBER")==0 || atoi(argv[1])==16)
        {
            if(ioPuffin->load(argv[2], M_LARGER_FIBER))
            {
                canDraw = true;
            }
            else
            {
                canDraw = false;
            }
        }
        else if(strcmp(argv[1],"M_RAW_FIBERS")==0 || atoi(argv[1])==32)
        {
            ioPuffin->Zmax = 10;
            ioPuffin->Zmin = 1;
            if(ioPuffin->load(argv[2], M_RAW_FIBERS))
            {
                canDraw = true;
            }
            else
            {
                canDraw = false;
            }
        }
        else if(strcmp(argv[1],"M_BOTH_LARGEST_FIBER")==0 || atoi(argv[1])==64)
        {
            if(ioPuffin->load(argv[2], M_BOTH_LARGEST_FIBER))
            {
                canDraw = true;
            }
            else
            {
                canDraw = false;
            }
        }
        else if(strcmp(argv[1],"M_SAVE")==0 || atoi(argv[1])==128)
        {
            canDraw = false;
            ioPuffin->save(argv[2]);
        }
        else if(strcmp(argv[1],"M_VTK_POLYDATA")==0 || atoi(argv[1])==256)
        {
            canDraw = true;
            ioPuffin->load(argv[2], M_VTK_POLYDATA);
        }
        else if(strcmp(argv[1],"M_HDR_READER")==0 || atoi(argv[1])==M_HDR_READER)
        {
            canDraw = false;
            ioPuffin->load(argv[2], M_HDR_READER);
        }
        else
        {
            canDraw = false;
            cout << "unhandle mode" << endl;
        }



        if(canDraw)
        {
            //draw it
            C_visu_draw* drawPuffin = new C_visu_draw();
            drawPuffin->draw(ioPuffin->data, 1, 0.5*3*0.25*4*0.2, 1.0, 1.0, 1.0);
            delete drawPuffin;
        }
        else
        {
            cout << "nothing can be drawn... together :-)" << endl;
        }

        //delete ioPuffin; //FIXME: big big error in the deletion of memory... garbage collector help me...
    }
    else if(argc==1)
    {
//        vector< vtkSmartPointer<vtkPoints> > points;// = vtkPoints::New(VTK_FLOAT);
//    // Create a line
//        unsigned int i=0;
//        for(float dx=0 ; dx<2.0 ; dx+=0.05)
//        {
//            points.push_back(vtkPoints::New(VTK_FLOAT));
//            points.at(i)->InsertNextPoint(0.0+dx,0.0,0.0);
//            points.at(i)->InsertNextPoint(0.0+dx,0.0,1.0);
//            points.at(i)->InsertNextPoint(1.0+dx,1.0,1.0);
//            points.at(i)->InsertNextPoint(1.0+dx,0.0,0.0);
//            i++;
//
//        }

//        vector< vtkSmartPointer<vtkPoints> > points;// = vtkPoints::New(VTK_FLOAT);
    // Create a line
//        unsigned int i=0;
//        points.push_back(vtkPoints::New(VTK_FLOAT));
//        points.at(i)->InsertNextPoint(0.0,0.0,0.0);
//        points.at(i)->InsertNextPoint(1.0,1.0,0.0);
//        i++;
//        points.push_back(vtkPoints::New(VTK_FLOAT));
//        points.at(i)->InsertNextPoint(0.0,1.0,0.0);
//        points.at(i)->InsertNextPoint(1.0,0.0,0.0);
//        i++;
//        points.push_back(vtkPoints::New(VTK_FLOAT));
//        points.at(i)->InsertNextPoint(0.0,0.0,0.0);
//        points.at(i)->InsertNextPoint(1.0,0.0,0.0);
//        points.at(i)->InsertNextPoint(2.0,0.0,0.0);
//        points.at(i)->InsertNextPoint(3.0,0.1,0.0);
//        points.at(i)->InsertNextPoint(2.0,0.2,0.0);
//        points.at(i)->InsertNextPoint(1.0,0.2,0.0);
//        points.at(i)->InsertNextPoint(0.0,0.2,0.0);
//        i++;
//
        C_visu_draw* drawPuffin = new C_visu_draw();
//        drawPuffin->draw(points, 9, 0.0+1, 0.0+1, 0.0+1, 0.03);



        CtlCtlStruct* Fib = new CtlCtlStruct(29);
        //Fib->N = 29;
        Fib->fibers = new CtlStruct[Fib->N];
        for(int i=0 ; i<Fib->N ; i++)
        {
            Fib->fibers[i].N_element = 2;
            Fib->fibers[i].elts = new double[3*Fib->fibers[i].N_element];
        }

        Fib->fibers[0].elts[0] = 0.0; Fib->fibers[0].elts[1] = 0.0; Fib->fibers[0].elts[2] = 0.0;
        Fib->fibers[0].elts[3] = 0.0; Fib->fibers[0].elts[4] = 1.0; Fib->fibers[0].elts[5] = 0.0;

        Fib->fibers[1].elts[0] = 0.0; Fib->fibers[1].elts[1] = 0.0; Fib->fibers[1].elts[2] = 0.0;
        Fib->fibers[1].elts[3] = 1.0; Fib->fibers[1].elts[4] = 0.0; Fib->fibers[1].elts[5] = 0.0;

        Fib->fibers[2].elts[0] = 0.0; Fib->fibers[2].elts[1] = 0.0; Fib->fibers[2].elts[2] = 0.0;
        Fib->fibers[2].elts[3] = 1.0; Fib->fibers[2].elts[4] = 1.0; Fib->fibers[2].elts[5] = 0.0;

        Fib->fibers[3].elts[0] = 0.0; Fib->fibers[3].elts[1] = 1.0; Fib->fibers[3].elts[2] = 0.0;
        Fib->fibers[3].elts[3] = 1.0; Fib->fibers[3].elts[4] = 0.0; Fib->fibers[3].elts[5] = 0.0;

        Fib->fibers[4].elts[0] = 0.0; Fib->fibers[4].elts[1] = 1.0; Fib->fibers[4].elts[2] = 0.0;
        Fib->fibers[4].elts[3] = 1.0; Fib->fibers[4].elts[4] = 1.0; Fib->fibers[4].elts[5] = 0.0;

        Fib->fibers[5].elts[0] = 0.0; Fib->fibers[5].elts[1] = 1.0; Fib->fibers[5].elts[2] = 0.0;
        Fib->fibers[5].elts[3] = 1.0; Fib->fibers[5].elts[4] = 2.0; Fib->fibers[5].elts[5] = 0.0;

        Fib->fibers[6].elts[0] = 0.0; Fib->fibers[6].elts[1] = 1.0; Fib->fibers[6].elts[2] = 0.0;
        Fib->fibers[6].elts[3] = 0.0; Fib->fibers[6].elts[4] = 2.0; Fib->fibers[6].elts[5] = 0.0;

        Fib->fibers[7].elts[0] = 0.0; Fib->fibers[7].elts[1] = 2.0; Fib->fibers[7].elts[2] = 0.0;
        Fib->fibers[7].elts[3] = 1.0; Fib->fibers[7].elts[4] = 1.0; Fib->fibers[7].elts[5] = 0.0;

        Fib->fibers[8].elts[0] = 0.0; Fib->fibers[8].elts[1] = 2.0; Fib->fibers[8].elts[2] = 0.0;
        Fib->fibers[8].elts[3] = 1.0; Fib->fibers[8].elts[4] = 2.0; Fib->fibers[8].elts[5] = 0.0;

        Fib->fibers[9].elts[0] = 1.0; Fib->fibers[9].elts[1] = 0.0; Fib->fibers[9].elts[2] = 0.0;
        Fib->fibers[9].elts[3] = 2.0; Fib->fibers[9].elts[4] = 0.0; Fib->fibers[9].elts[5] = 0.0;

        Fib->fibers[10].elts[0] = 1.0; Fib->fibers[10].elts[1] = 0.0; Fib->fibers[10].elts[2] = 0.0;
        Fib->fibers[10].elts[3] = 2.0; Fib->fibers[10].elts[4] = 1.0; Fib->fibers[10].elts[5] = 0.0;

        Fib->fibers[11].elts[0] = 1.0; Fib->fibers[11].elts[1] = 0.0; Fib->fibers[11].elts[2] = 0.0;
        Fib->fibers[11].elts[3] = 1.0; Fib->fibers[11].elts[4] = 1.0; Fib->fibers[11].elts[5] = 0.0;

        Fib->fibers[12].elts[0] = 1.0; Fib->fibers[12].elts[1] = 1.0; Fib->fibers[12].elts[2] = 0.0;
        Fib->fibers[12].elts[3] = 2.0; Fib->fibers[12].elts[4] = 0.0; Fib->fibers[12].elts[5] = 0.0;

        Fib->fibers[13].elts[0] = 1.0; Fib->fibers[13].elts[1] = 1.0; Fib->fibers[13].elts[2] = 0.0;
        Fib->fibers[13].elts[3] = 2.0; Fib->fibers[13].elts[4] = 1.0; Fib->fibers[13].elts[5] = 0.0;

        Fib->fibers[14].elts[0] = 1.0; Fib->fibers[14].elts[1] = 1.0; Fib->fibers[14].elts[2] = 0.0;
        Fib->fibers[14].elts[3] = 2.0; Fib->fibers[14].elts[4] = 2.0; Fib->fibers[14].elts[5] = 0.0;

        Fib->fibers[15].elts[0] = 1.0; Fib->fibers[15].elts[1] = 1.0; Fib->fibers[15].elts[2] = 0.0;
        Fib->fibers[15].elts[3] = 1.0; Fib->fibers[15].elts[4] = 2.0; Fib->fibers[15].elts[5] = 0.0;

        Fib->fibers[16].elts[0] = 1.0; Fib->fibers[16].elts[1] = 2.0; Fib->fibers[16].elts[2] = 0.0;
        Fib->fibers[16].elts[3] = 2.0; Fib->fibers[16].elts[4] = 1.0; Fib->fibers[16].elts[5] = 0.0;

        Fib->fibers[17].elts[0] = 1.0; Fib->fibers[17].elts[1] = 2.0; Fib->fibers[17].elts[2] = 0.0;
        Fib->fibers[17].elts[3] = 2.0; Fib->fibers[17].elts[4] = 2.0; Fib->fibers[17].elts[5] = 0.0;

        Fib->fibers[18].elts[0] = 2.0; Fib->fibers[18].elts[1] = 0.0; Fib->fibers[18].elts[2] = 0.0;
        Fib->fibers[18].elts[3] = 3.0; Fib->fibers[18].elts[4] = 0.0; Fib->fibers[18].elts[5] = 0.0;

        Fib->fibers[19].elts[0] = 2.0; Fib->fibers[19].elts[1] = 0.0; Fib->fibers[19].elts[2] = 0.0;
        Fib->fibers[19].elts[3] = 3.0; Fib->fibers[19].elts[4] = 1.0; Fib->fibers[19].elts[5] = 0.0;

        Fib->fibers[20].elts[0] = 2.0; Fib->fibers[20].elts[1] = 0.0; Fib->fibers[20].elts[2] = 0.0;
        Fib->fibers[20].elts[3] = 2.0; Fib->fibers[20].elts[4] = 1.0; Fib->fibers[20].elts[5] = 0.0;

        Fib->fibers[21].elts[0] = 2.0; Fib->fibers[21].elts[1] = 1.0; Fib->fibers[21].elts[2] = 0.0;
        Fib->fibers[21].elts[3] = 3.0; Fib->fibers[21].elts[4] = 0.0; Fib->fibers[21].elts[5] = 0.0;

        Fib->fibers[22].elts[0] = 2.0; Fib->fibers[22].elts[1] = 1.0; Fib->fibers[22].elts[2] = 0.0;
        Fib->fibers[22].elts[3] = 3.0; Fib->fibers[22].elts[4] = 1.0; Fib->fibers[22].elts[5] = 0.0;

        Fib->fibers[23].elts[0] = 2.0; Fib->fibers[23].elts[1] = 1.0; Fib->fibers[23].elts[2] = 0.0;
        Fib->fibers[23].elts[3] = 3.0; Fib->fibers[23].elts[4] = 2.0; Fib->fibers[23].elts[5] = 0.0;

        Fib->fibers[24].elts[0] = 2.0; Fib->fibers[24].elts[1] = 1.0; Fib->fibers[24].elts[2] = 0.0;
        Fib->fibers[24].elts[3] = 2.0; Fib->fibers[24].elts[4] = 2.0; Fib->fibers[24].elts[5] = 0.0;

        Fib->fibers[25].elts[0] = 2.0; Fib->fibers[25].elts[1] = 2.0; Fib->fibers[25].elts[2] = 0.0;
        Fib->fibers[25].elts[3] = 3.0; Fib->fibers[25].elts[4] = 1.0; Fib->fibers[25].elts[5] = 0.0;

        Fib->fibers[26].elts[0] = 2.0; Fib->fibers[26].elts[1] = 2.0; Fib->fibers[26].elts[2] = 0.0;
        Fib->fibers[26].elts[3] = 3.0; Fib->fibers[26].elts[4] = 2.0; Fib->fibers[26].elts[5] = 0.0;

        Fib->fibers[27].elts[0] = 3.0; Fib->fibers[27].elts[1] = 0.0; Fib->fibers[27].elts[2] = 0.0;
        Fib->fibers[27].elts[3] = 3.0; Fib->fibers[27].elts[4] = 1.0; Fib->fibers[27].elts[5] = 0.0;

        Fib->fibers[28].elts[0] = 3.0; Fib->fibers[28].elts[1] = 1.0; Fib->fibers[28].elts[2] = 0.0;
        Fib->fibers[28].elts[3] = 3.0; Fib->fibers[28].elts[4] = 2.0; Fib->fibers[28].elts[5] = 0.0;

        //drawPuffin->draw(Fib, 1, 1.0, 1.0, 1.0, 0.05);
        vector< vtkSmartPointer< vtkActor > > fibActor =  drawPuffin->drawFibers(Fib, 1, 0.03);

        TensorData* tens = new TensorData[12];
        for(int i=0 ; i<4 ; i++)
        {
            double mX = (double) 2*(i+2);
            for(int j=1 ; j<4 ; j++)
            {
                double mY = (double) 2*j+2;
                tens[i*3+j-1].x = (double) i;
                tens[i*3+j-1].y = (double) j-1;
                tens[i*3+j-1].z = 0.0;

                tens[i*3+j-1].eigenValue[0] = 0.5;
                tens[i*3+j-1].eigenValue[1] = 0.1;
                tens[i*3+j-1].eigenValue[2] = 0.1;

                tens[i*3+j-1].eigenVector1[0] = -sin(atan2(mY,mX));
                tens[i*3+j-1].eigenVector1[1] = cos(atan2(mY,mX));
                tens[i*3+j-1].eigenVector1[2] = 0.0;

                tens[i*3+j-1].eigenVector2[0] = cos(atan2(mY,mX));
                tens[i*3+j-1].eigenVector2[1] = sin(atan2(mY,mX));
                tens[i*3+j-1].eigenVector2[2] = 0.0;

                tens[i*3+j-1].eigenVector3[0] = 0.0;
                tens[i*3+j-1].eigenVector3[1] = 0.0;
                tens[i*3+j-1].eigenVector3[2] = 1.0;

                tens[i*3+j-1].phi = 0.0;
                tens[i*3+j-1].psi = 0.0;
                tens[i*3+j-1].theta = atan2(mY,mX);
            }
        }
        vtkSmartPointer<vtkRenderer> renderer = vtkSmartPointer<vtkRenderer>::New();
        vector< vtkSmartPointer< vtkActor > > actorTensor = drawPuffin->drawTensors(tens, 12);
        for(unsigned int j=0 ; j<actorTensor.size() ; j++ )
        {
            renderer->AddActor(actorTensor.at(j));
        }

        for(unsigned int j=0 ; j<fibActor.size() ; j++ )
        {
            renderer->AddActor(fibActor.at(j));
        }
        renderer->SetBackground(1.0, 1.0, 1.0);
        renderer->ResetCamera();

        vtkSmartPointer<vtkRenderWindow> renWin = vtkSmartPointer<vtkRenderWindow>::New();
        vtkSmartPointer<vtkRenderWindowInteractor> iren = vtkSmartPointer<vtkRenderWindowInteractor>::New();
        iren->SetRenderWindow(renWin);
        renWin->AddRenderer(renderer);
        renWin->SetSize(500, 500);
        renWin->Render();

        vtkSmartPointer<vtkInteractorStyleTrackballCamera> style = vtkSmartPointer<vtkInteractorStyleTrackballCamera>::New();
        iren->SetInteractorStyle(style);

        iren->Start();


        plotData* data = new plotData();
        data->Fib = Fib;
        data->T = tens;
        data->nbTensor = 12;
        data->rawFibers = 2;
        //drawPuffin->draw(data, 1, 0.02, 1.0, 1.0, 1.0);

        delete drawPuffin;
        delete data;
    }
    else
    {
        cout << "no file or mode specified." << endl;
        return -1;
    }
    return 0;
}

void showTensor(TensorData* T)
 {
     cout << "(x,y,z) = (" << T->x << "," << T->y << "," << T->z << ")" << endl;
     cout << "eigen value (vp1,vp2,vp3) = (" << T->eigenValue[0] << "," << T->eigenValue[1] << "," << T->eigenValue[2] << ")" << endl;
     cout << "main direction (Vx,Vy,Vz) = (" << T->eigenVector1[0] << "," << T->eigenVector1[1] << "," << T->eigenVector1[2] << ")" << endl;
     cout << "(tetha,phi,psi) = (" << T->theta << "," << T->phi << "," << T->psi << ")" << endl;
 }

vtkSmartPointer<vtkActor> ellipsoid(double Xc, double Yc, double Zc, double l1/** > */, double l2/** > */, double l3, double theta, double phi, double psi, double vp1[3])
{
    vtkSmartPointer<vtkParametricEllipsoid> parametricObject = vtkSmartPointer<vtkParametricEllipsoid>::New();

    vtkSmartPointer<vtkParametricFunctionSource> ellipsoidSource = vtkSmartPointer<vtkParametricFunctionSource>::New();
    parametricObject->SetXRadius(l1);
    parametricObject->SetYRadius(l2);
    parametricObject->SetZRadius(l3);
    ellipsoidSource->SetParametricFunction(parametricObject);
    ellipsoidSource->Update();

    vtkSmartPointer<vtkTransform> translation = vtkSmartPointer<vtkTransform>::New();

    translation->Translate(Xc, Yc, Zc);
    translation->RotateZ(180*theta/3.14161);
    translation->RotateY(180*phi/3.14161);
    translation->RotateX(180*psi/3.14161);
    double FA = sqrt( (SQR(l1-l2) + SQR(l2-l3) + SQR(l3-l1) ) / ( SQR(l1) +SQR(l2) + SQR(l3) ) );
    translation->Scale(FA+0.1, FA+0.1, FA+0.1);

    vtkSmartPointer<vtkTransformPolyDataFilter> transformFilter = vtkSmartPointer<vtkTransformPolyDataFilter>::New();
    transformFilter->SetInputConnection(ellipsoidSource->GetOutputPort());
    transformFilter->SetTransform(translation);
    transformFilter->Update();

    // Create a mapper and actor
    vtkSmartPointer<vtkPolyDataMapper> transformedMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    transformedMapper->SetInputConnection(transformFilter->GetOutputPort());

    vtkSmartPointer<vtkActor> transformedActor = vtkSmartPointer<vtkActor>::New();
    transformedActor->SetMapper(transformedMapper);
    double Z =  sqrt( SQR(vp1[0]) +SQR(vp1[1]) + SQR(vp1[2]));
    transformedActor->GetProperty()->SetColor(ABS(vp1[0])/Z,ABS(vp1[1])/Z,ABS(vp1[2])/Z);
    return transformedActor;
}
