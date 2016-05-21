#include <C_visu_draw.h>

C_visu_draw::C_visu_draw()
{
    //ctor
//    Zmax = HIGHZ;
//    Zmin = LOWZ;
    m_FAth = 0.227;
}

C_visu_draw::~C_visu_draw()
{
    //dtor
}
void C_visu_draw::draw(plotData* data, unsigned long nbSubDiv, double radius, float R, float G, float B)
{
    vtkSmartPointer<vtkRenderer> renderer = vtkSmartPointer<vtkRenderer>::New();

    ///the tensor showing part must be modified! may help: http://www.cmake.org/Wiki/VTK/Examples/Cxx/Filtering/vtkVertexGlyphFilter
    ///this is something similar to what I did http://www.cmake.org/Wiki/VTK/Examples/Visualization/VectorOfActors
    vector< vtkSmartPointer< vtkActor > > actorTensor;
    if(data->nbTensor!=0)
    {
        actorTensor = drawTensors(data->T, data->nbTensor);
        for(unsigned int j=0 ; j<actorTensor.size() ; j++ )
        {
            renderer->AddActor(actorTensor.at(j));
        }
    }

    vector< vtkSmartPointer< vtkActor > > actorFib;
    if(data->Fib!=NULL)
    {
        if(data->Fib->N!=0)
        {
            if(data->rawFibers==1)
            {
                actorFib = drawRawFibers(data->Fib, radius);
            }
            else if(data->rawFibers==2)
            {
                actorFib = drawFibers(data->Fib, nbSubDiv, radius);
                actorFib.push_back(drawRawFibers(data->Fib, radius).at(0));
            }
            else
            {
                actorFib = drawFibers(data->Fib, nbSubDiv, radius);
            }
            for(unsigned int j=0 ; j<actorFib.size() ; j++ )
            {
                renderer->AddActor(actorFib.at(j));
            }
        }
    }


    renderer->SetBackground(R, G, B);

    // Make an oblique view
    //renderer->GetActiveCamera()->Azimuth(30);
    //renderer->GetActiveCamera()->Azimuth(0);
    //renderer->GetActiveCamera()->Pitch(0);
    //renderer->GetActiveCamera()->Yaw(0);
    //renderer->GetActiveCamera()->Roll(0);
    renderer->GetActiveCamera()->Elevation(30);
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
    return;
}



void C_visu_draw::draw(CtlCtlStruct* m_ctlCtlStruct, long nbSubDiv, float R, float G, float B, double radius)
{
    draw(run(m_ctlCtlStruct, nbSubDiv, radius), R, G, B);
    return;
}

//void C_visu_draw::draw(vector< vtkSmartPointer<vtkPoints> > points, long nbSubDiv, float R, float G, float B, double radius)
//{
//    draw(run(points, nbSubDiv, radius), R, G, B);
//    return;
//}

void C_visu_draw::draw(vector< vtkSmartPointer<vtkPolyDataMapper> > vectorMap, float R, float G, float B)
{
    vtkSmartPointer<vtkRenderer> renderer = vtkSmartPointer<vtkRenderer>::New();
    vector< vtkSmartPointer<vtkActor> > actor;// = vtkSmartPointer<vtkActor>::New();

    for(unsigned int j=0 ; j<vectorMap.size() ; j++ )
    {
        actor.push_back(vtkSmartPointer<vtkActor>::New());
        actor.at(j)->SetMapper(vectorMap.at(j));
        renderer->AddActor(actor.at(j));
    }


    //renderer->SetBackground(.2, .3, .4);
    renderer->SetBackground(R, G, B);
    // Make an oblique view
    renderer->GetActiveCamera()->Azimuth(0);//30
    renderer->GetActiveCamera()->Elevation(30);
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
    return;
}
void C_visu_draw::computeSpline(vtkSmartPointer<vtkPoints> points, long nbSubDiv, double radius, vtkSmartPointer<vtkTubeFilter> tube)
{
    int nb = points->GetNumberOfPoints();
    vtkSmartPointer<vtkParametricSpline> parametricObject = vtkSmartPointer<vtkParametricSpline>::New();
    parametricObject->SetPoints(points);
    vtkSmartPointer<vtkParametricFunctionSource> lineSource = vtkSmartPointer<vtkParametricFunctionSource>::New();
    lineSource->SetUResolution(nb*nbSubDiv);
    lineSource->SetParametricFunction(parametricObject);
    lineSource->Update();

    //get polyData from spline
    vtkSmartPointer<vtkPolyData> polyData = lineSource->GetOutput();
    computeColor(polyData);

    tube->SetInput(polyData);
    tube->SetNumberOfSides(16);
    tube->SetRadius(radius);

    return;
}
vtkSmartPointer<vtkPolyDataMapper> C_visu_draw::computeMapper(vtkSmartPointer<vtkTubeFilter> tube)
{
    vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    mapper->SetInputConnection(tube->GetOutputPort());
    mapper->ScalarVisibilityOn();
    mapper->SetScalarModeToUsePointFieldData();
    mapper->SelectColorArray("Colors");
    return mapper;
}
vector< vtkSmartPointer<vtkPolyDataMapper> > C_visu_draw::run(vector< vtkSmartPointer<vtkPoints> > P, long nbSubDiv, double radius)
{
    vector< vtkSmartPointer<vtkPolyDataMapper> > vectMap;
    for(unsigned int i=0 ; i<P.size() ; i++)
    {
        vtkSmartPointer<vtkTubeFilter> tube = vtkSmartPointer<vtkTubeFilter>::New();
        computeSpline(P.at(i), nbSubDiv, radius, tube);
        vectMap.push_back(computeMapper(tube));
    }
    return vectMap;
}
vector< vtkSmartPointer<vtkPolyDataMapper> > C_visu_draw::run(CtlCtlStruct* m_ctlCtlStruct, long nbSubDiv, double radius)
{
    vector< vtkSmartPointer<vtkPoints> > P;
    for(unsigned long i=0 ; i<m_ctlCtlStruct->N ; i++)
    {
        P.push_back(vtkPoints::New(VTK_FLOAT));
        //insert points
        for(unsigned long j=0 ; j<m_ctlCtlStruct->fibers[i].N_element ; j++)
        {
            P.at(i)->InsertNextPoint(m_ctlCtlStruct->fibers[i].elts[3*j],m_ctlCtlStruct->fibers[i].elts[3*j+1],m_ctlCtlStruct->fibers[i].elts[3*j+2]);
        }
    }
    return run(P, nbSubDiv, radius);
}

void C_visu_draw::computeColor(vtkSmartPointer<vtkPolyData> polyData)
{
    //get number of point in polyData
    int N = polyData->GetNumberOfPoints();

    //get pointer to point set
    vtkPoints* P = polyData->GetPoints();


    //compute color mapping
    vtkSmartPointer<vtkDoubleArray> colors = vtkSmartPointer<vtkDoubleArray>::New();
    colors->SetName("Colors");
    colors->SetNumberOfComponents(3);
    colors->SetNumberOfTuples(N);

    double r1[3];
    double r2[3];
    double R=1.0, G=1.0, B=0.0;
    P->GetPoint(0,r1);
    P->GetPoint(1,r2);
    R = /*0.0*/ABS(r2[0] - r1[0]);
    G = /*0.0*/ABS(r2[1] - r1[1]);
    B = /*1+0.0*/ABS(r2[2] - r1[2]);

    B = /*0.0*/ABS(r2[0] - r1[0]);
    R = /*0.0*/ABS(r2[1] - r1[1]);
    G = /*1+0.0*/ABS(r2[2] - r1[2]);
    double N2 = sqrt(R*R + B*B + G*G);
    R /= N2;
    G /= N2;
    B /= N2;
    colors->InsertTuple3(0,R,G,B);
    for(int t=1 ; t<N-1 ; t++)
    {
        P->GetPoint(t-1,r1);
        P->GetPoint(t+1,r2);
        R = /*0.0*/ABS(r2[0] - r1[0]);
        G = /*0.0*/ABS(r2[1] - r1[1]);
        B = /*1+0.0*/ABS(r2[2] - r1[2]);
        N2 = sqrt(R*R + B*B + G*G);
        R /= N2;
        G /= N2;
        B /= N2;
        colors->InsertTuple3(t,R,G,B);
    }
    P->GetPoint(N-2,r1);
    P->GetPoint(N-1,r2);
    R = /*0.0*/ABS(r2[0] - r1[0]);
    G = /*0.0*/ABS(r2[1] - r1[1]);
    B = /*1+0.0*/ABS(r2[2] - r1[2]);

    B = /*0.0*/ABS(r2[0] - r1[0]);
    R = /*0.0*/ABS(r2[1] - r1[1]);
    G = /*1+0.0*/ABS(r2[2] - r1[2]);

    N2 = sqrt(R*R + B*B + G*G);
    R /= N2;
    G /= N2;
    B /= N2;
    colors->InsertTuple3(N-1,R,G,B);
    polyData->GetPointData()->AddArray(colors);
    return;
}











vector< vtkSmartPointer< vtkActor > > C_visu_draw::drawTensors(TensorData* T, unsigned long nbTensor)
{
    vector< vtkSmartPointer< vtkActor > > v;
    for(unsigned long i=0 ; i<nbTensor ; i++)
    {
        if( sqrt( 0.5*(SQR(T[i].eigenValue[0] - T[i].eigenValue[1]) + SQR(T[i].eigenValue[1] - T[i].eigenValue[2]) + SQR(T[i].eigenValue[2] - T[i].eigenValue[0]))/(SQR(T[i].eigenValue[0]) + SQR(T[i].eigenValue[1]) + SQR(T[i].eigenValue[2])) ) >m_FAth )
        {
           v.push_back(drawTensor(&T[i]));
        }
    }
    return v;
}
vtkSmartPointer< vtkActor > C_visu_draw::drawTensor(TensorData* T)
{
    vtkSmartPointer<vtkParametricEllipsoid> parametricObject = vtkSmartPointer<vtkParametricEllipsoid>::New();

    vtkSmartPointer<vtkParametricFunctionSource> ellipsoidSource = vtkSmartPointer<vtkParametricFunctionSource>::New();
    double Zeig = 1.2*sqrt(SQR(T->eigenValue[0]) + SQR(T->eigenValue[1]) + SQR(T->eigenValue[2]));
    parametricObject->SetXRadius(T->eigenValue[0]/Zeig);
    parametricObject->SetYRadius(T->eigenValue[1]/Zeig);
    parametricObject->SetZRadius(T->eigenValue[2]/Zeig);
    ellipsoidSource->SetParametricFunction(parametricObject);
    ellipsoidSource->Update();

    vtkSmartPointer<vtkTransform> translation = vtkSmartPointer<vtkTransform>::New();

    translation->Translate(T->x, T->y, T->z);
    translation->RotateZ(180*(T->theta)/3.14161);
    translation->RotateY(180*(T->phi)/3.14161);
    translation->RotateX(180*(T->psi)/3.14161);
    //double FA = sqrt( 0.5* (SQR(T->eigenValue[0]-T->eigenValue[1]) + SQR(T->eigenValue[1]-T->eigenValue[2]) + SQR(T->eigenValue[2]-T->eigenValue[0]) ) / ( SQR(T->eigenValue[0]) +SQR(T->eigenValue[1]) + SQR(T->eigenValue[2]) ) );
    //FA = FA/2;
    //translation->Scale(FA+0.1, FA+0.1, FA+0.1);

    vtkSmartPointer<vtkTransformPolyDataFilter> transformFilter = vtkSmartPointer<vtkTransformPolyDataFilter>::New();
    transformFilter->SetInputConnection(ellipsoidSource->GetOutputPort());
    transformFilter->SetTransform(translation);
    transformFilter->Update();

    // Create a mapper and actor
    vtkSmartPointer<vtkPolyDataMapper> transformedMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    transformedMapper->SetInputConnection(transformFilter->GetOutputPort());

    vtkSmartPointer<vtkActor> transformedActor = vtkSmartPointer<vtkActor>::New();
    transformedActor->SetMapper(transformedMapper);
    double Z =  sqrt( SQR(T->eigenVector1[0]) +SQR(T->eigenVector1[1]) + SQR(T->eigenVector1[2]));
    transformedActor->GetProperty()->SetColor(ABS(T->eigenVector1[0])/Z,ABS(T->eigenVector1[1])/Z,ABS(T->eigenVector1[2])/Z);

    return transformedActor;
}







double C_visu_draw::fibLength(CtlStruct* fibers)
{
    if(fibers->N_element<2)
    {
        return 0;
    }
    double L = 0;
    for(unsigned long i=0 ; i<fibers->N_element-1 ; i++)
    {
        L += sqrt( SQR(fibers->elts[3*i] - fibers->elts[3*(i+1)]) + SQR(fibers->elts[3*i+1] - fibers->elts[3*(i+1)+1]) + SQR(fibers->elts[3*i+2] - fibers->elts[3*(i+1)+2]) );
    }
    return L;
}

vector< vtkSmartPointer< vtkActor > > C_visu_draw::drawFibers(CtlCtlStruct* Fib, unsigned long nbSubDiv, double radius)
{
    //compute spline approximation
    C_toolbox_spline* spline = new C_toolbox_spline();
    CtlCtlStruct** splineStruct = new CtlCtlStruct*;
    spline->computeFiberSpline(Fib, splineStruct, nbSubDiv);
    vector< vtkSmartPointer< vtkActor > > v;
    for(unsigned long i=0 ; i<Fib->N ; i++)
    {
        if(fibLength(&((*splineStruct)->fibers[i]))>MIN_PHYSICAL_LENGTH)
        {
            v.push_back(drawFiber( &((*splineStruct)->fibers[i]), radius ));
        }
    }
    delete spline;//->~C_spline();
    return v;
}
vector< vtkSmartPointer< vtkActor > > C_visu_draw::drawRawFibers(CtlCtlStruct* Fib, double radius)
{
    //compute spline approximation
    vector< vtkSmartPointer< vtkActor > > v;
    for(unsigned long i=0 ; i<Fib->N ; i++)
    {
        if(fibLength(&(Fib->fibers[i]))>MIN_PHYSICAL_LENGTH)
        {
            v.push_back(drawFiber( &(Fib->fibers[i]), radius ));
        }
    }
    return v;
}
vtkSmartPointer< vtkActor > C_visu_draw::drawFiber(CtlStruct* Fib, double radius)
{
    //move fiber point in vtkpoints
    vtkSmartPointer<vtkPoints> P = vtkPoints::New(VTK_FLOAT);
    vtkSmartPointer<vtkCellArray> lineVect = vtkSmartPointer<vtkCellArray>::New();
    vtkSmartPointer<vtkPolyData> polyDataVect = vtkPolyData::New();

    lineVect->InsertNextCell( Fib->N_element-1 );  //should I return to the version that remove the last point? yes, but it's due to the file format
    for(unsigned long j=0 ; j<Fib->N_element-1 ; j++)
    {
        //insert point
        P->InsertNextPoint(Fib->elts[3*j],Fib->elts[3*j+1],Fib->elts[3*j+2]);
        lineVect->InsertCellPoint(j);
    }
    lineVect->UpdateCellCount( Fib->N_element-1 );

    //fill polydata
    polyDataVect->SetPoints(P);
    polyDataVect->SetLines(lineVect);

    //compute a color for each point of the line
    computeColor(polyDataVect);

    //compute spline estimation and store result in tube
    vtkSmartPointer<vtkTubeFilter> tube = vtkSmartPointer<vtkTubeFilter>::New();

    //set tube
    tube->SetInput(polyDataVect);
    tube->SetNumberOfSides(25);//16
    tube->SetRadius(radius);
    tube->Update();

    //computeSpline(P, nbSubDiv, radius, tube);

    //set a mapper
    vtkSmartPointer<vtkPolyDataMapper> vectMap;
    vectMap = computeMapper(tube);

    //set an actor from mapper
    vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
    actor->SetMapper(vectMap);
    return actor;
}

