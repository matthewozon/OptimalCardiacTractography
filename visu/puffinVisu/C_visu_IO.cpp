#include <C_visu_IO.h>


C_visu_IO::C_visu_IO()
{
    //ctor
    data = new plotData;
    data->rawFibers = 0;
    Zmax = -1;
    Zmin = 1;
    //m_ctlCtlStruct = NULL;
}

C_visu_IO::~C_visu_IO()
{
    //dtor
    delete data;
//    if(m_ctlCtlStruct!=NULL)
//    {
//        delete m_ctlCtlStruct;
//    }
}

bool C_visu_IO::load(string fileName, long kind)
{
    if(kind==M_TENSOR)
    {
        return loadTensorFromPointData(fileName);
    }
    else if(kind==M_EDGE)
    {
        #ifdef LATER
            return loadEdge(fileName);
        #else
            return false;
        #endif
    }
    else if(kind==M_FIBERS || kind==M_RAW_FIBERS)
    {
        if(kind==M_RAW_FIBERS)
        {
            data->rawFibers = 1;
        }
        return readFiber(fileName);
    }
    else if(kind==M_LARGER_FIBER_AND_CONTROL_POINTS || kind==M_LARGER_FIBER || kind==M_BOTH_LARGEST_FIBER)
    {
        if(kind==M_LARGER_FIBER_AND_CONTROL_POINTS)
        {
            data->rawFibers = 1;
        }
        if(kind==M_BOTH_LARGEST_FIBER)
        {
            data->rawFibers = 2;
        }
        return readLargestFiber(fileName);
    }
    else if(kind==M_VTK_POLYDATA)
    {
        vtkSmartPointer< vtkPolyDataReader > reader = vtkSmartPointer< vtkPolyDataReader >::New();
        reader->SetFileName( "/home/ozon/Documents/codeblocks_workspace/test_algo_workspace/visu/puffinVisu/dummy.vtk" ); //(fileName.append(".vtk")).data()
        if(reader->OpenVTKFile()==0)
        {
            cout << "pb occured while reading data" << endl;
        }
        cout << reader->GetOutput()->GetNumberOfPoints() << endl;;
        if(reader->GetOutput()==NULL)
        {
            cout << "no way!" << endl;
        }
        reader->GetOutput()->Print(cout);
        //data->polyDataV->Print(cout);
        return true;
    }
    else if(kind==M_HDR_READER)
    {
        ///show full hdr file
        analyzeInfo(fileName.data());
        return false;
    }
    return false;
}


void C_visu_IO::analyzeInfo(const char* filenameHDR)
{
    C_ReadAnalyze75* ANALYSE = new C_ReadAnalyze75();
    //load informations in INFO
    dsr* INFO = ANALYSE->readAnalyzeInfo(filenameHDR);
    delete ANALYSE;

    ///show everythings in INFO, but what is empty
    cout << "information about: " << filenameHDR << endl;
    cout << "header_key" << endl;
    cout << endl;
    cout << "\t" << "size of header: " << INFO->hk.sizeof_hdr << endl;
    cout << "\t" << "data type: " << INFO->hk.data_type << endl;
    cout << "\t" << "data base name: " << INFO->hk.db_name << endl;
    cout << "\t" << "extents (usually 16384): " << INFO->hk.extents << endl;
    cout << "\t" << "session error: " << INFO->hk.session_error << endl;
    cout << "\t" << "regular (r indicated that all images are of same size): " << INFO->hk.regular << endl;
    cout << "\t" << "unused char: " << INFO->hk.hkey_un0 << endl;
    cout << endl;
    cout << endl;
    cout << "image_dimension" << endl;
    cout << endl;
    cout << "\t" << "number of dimension in DB (usually 4): " << INFO->dime.dim[0] << endl;
    cout << "\t" << "Image X dimension; number of pixels in an image row: " << INFO->dime.dim[1] << endl;
    cout << "\t" << "Image Y dimension; number of pixel rows in slice: " << INFO->dime.dim[2] << endl;
    cout << "\t" << "Volume Z dimension; number of slices in a volume: " << INFO->dime.dim[3] << endl;
    cout << "\t" << "Time points, number of volumes in database: " << INFO->dime.dim[4] << endl;
    cout << "\t" << "other dimension: " << INFO->dime.dim[5] << endl;
    cout << "\t" << "other dimension: " << INFO->dime.dim[6] << endl;
    cout << "\t" << "other dimension: " << INFO->dime.dim[7] << endl;
    cout << "\t" << "unused 8: " << INFO->dime.unused8 << endl;
    cout << "\t" << "unused 9: " << INFO->dime.unused9 << endl;
    cout << "\t" << "unused 10: " << INFO->dime.unused10 << endl;
    cout << "\t" << "unused 11: " << INFO->dime.unused11 << endl;
    cout << "\t" << "unused 12: " << INFO->dime.unused12 << endl;
    cout << "\t" << "unused 13: " << INFO->dime.unused13 << endl;
    cout << "\t" << "unused 14: " << INFO->dime.unused14 << endl;
//    if(==)
    switch(INFO->dime.datatype)
    {
        case DT_NONE:
        {
            cout << "\t" << "no data or unknown type" << endl;
            break;
        }
        case DT_BINARY:
        {
            cout << "\t" << "binary data" << endl;
            break;
        }
        case DT_UNSIGNED_CHAR:
        {
            cout << "\t" << "unsigned char data" << endl;
            break;
        }
        case DT_SIGNED_SHORT:
        {
            cout << "\t" << "signed short data" << endl;
            break;
        }
        case DT_SIGNED_INT:
        {
            cout << "\t" << "signed int data" << endl;
            break;
        }
        case DT_FLOAT:
        {
            cout << "\t" << "float data" << endl;
            break;
        }
        case DT_COMPLEX:
        {
            cout << "\t" << "complex data" << endl;
            break;
        }
        case DT_DOUBLE:
        {
            cout << "\t" << "double data" << endl;
            break;
        }
        case DT_RGB:
        {
            cout << "\t" << "RGB data" << endl;
            break;
        }
        case DT_ALL:
        {
            cout << "\t" << "ALL data" << endl;
            break;
        }
        default:
        {
            cout << "\t" << "no data" << endl;
            break;
        }
    }
    cout << "\t" << "bit per voxel: " << INFO->dime.bitpix << endl;
    cout << "\t" << "voxel dimensitons: giving real world measurements in mm. and ms.: " << INFO->dime.pixdim[0] << endl;
    cout << "\t" << "voxel width in mm: " << INFO->dime.pixdim[1] << endl;
    cout << "\t" << "voxel height in mm: " << INFO->dime.pixdim[2] << endl;
    cout << "\t" << "interslice in mm: " << INFO->dime.pixdim[3] << endl;
    cout << "\t" << "interframe delay in ms: " << INFO->dime.pixdim[4] << endl;
    cout << "\t" << "other measurement: " << INFO->dime.pixdim[5] << endl;
    cout << "\t" << "other measurement: " << INFO->dime.pixdim[6] << endl;
    cout << "\t" << "other measurement: " << INFO->dime.pixdim[7] << endl;
    cout << "\t" << "voxel offset: " << INFO->dime.vox_offset << endl;
    cout << "\t" << "funused1: " << INFO->dime.funused1 << endl;
    cout << "\t" << "funused2: " << INFO->dime.funused2 << endl;
    cout << "\t" << "funused3: " << INFO->dime.funused3 << endl;
    cout << "\t" << "max calibration value: " << INFO->dime.cal_max << endl;
    cout << "\t" << "min calibration value: " << INFO->dime.cal_min << endl;
    cout << "\t" << "compression rate: " << INFO->dime.compressed << endl;
    cout << "\t" << "verified: " << INFO->dime.verified << endl;
    cout << "\t" << "max value in db: " << INFO->dime.glmax << endl;
    cout << "\t" << "min value in db: " << INFO->dime.glmin << endl;

    cout << endl;
    cout << endl;
    cout << "data_history" << endl;
    cout << endl;
    cout << "\t" << "description: " << INFO->hist.descrip << endl;
    cout << "\t" << "auxiliar files: " << INFO->hist.aux_file << endl;
    cout << "\t" << "orient: " << INFO->hist.orient << endl;
    cout << "\t" << "originator: " << INFO->hist.originator << endl;
    cout << "\t" << "generated: " << INFO->hist.generated << endl;
    cout << "\t" << "scannum: " << INFO->hist.scannum << endl;
    cout << "\t" << "patient_id: " << INFO->hist.patient_id << endl;
    cout << "\t" << "date: " << INFO->hist.exp_date << endl;
    cout << "\t" << "time: " << INFO->hist.exp_time << endl;
    cout << "\t" << "number of views: " << INFO->hist.views << endl;
    cout << "\t" << "vols_added: " << INFO->hist.vols_added << endl;
    cout << "\t" << "start_field: " << INFO->hist.start_field << endl;
    cout << "\t" << "field_skip: " << INFO->hist.field_skip << endl;
    cout << "\t" << "omax: " << INFO->hist.omax << " omin: " << INFO->hist.omin << endl;
    cout << "\t" << "smax: " << INFO->hist.smax << " smin: " << INFO->hist.smin << endl;
    return;
}


string C_visu_IO::getBaseName(string fileName)
{
    size_t found = fileName.find_last_of(".");
    string str;
    if(found==string::npos)
    {
        str = fileName;
    }
    else
    {
        size_t S = fileName.size();
        if(found<S-7)
        {
            str = fileName;
        }
        else
        {
            str = fileName.substr(0,found);
        }
    }
    return str;

}
bool C_visu_IO::readFiber(string fileNameFib)//, CtlCtlStruct* fibersCtl)
{
    //file names
    string headerFile = getBaseName(fileNameFib) + ".fibHDR";
    string sourceFile = getBaseName(fileNameFib) + ".fibSRC";

    //
    fstream fhdr;
    fhdr.open (headerFile.data(), ifstream::in | ifstream::binary ); //
    if(fhdr.is_open())
    {
        //
        unsigned long nbBundlePoints = 0;

        //get number of data = nb_bundle+1
        fhdr.seekg (0, ios::end);
        unsigned long length = fhdr.tellg();
        fhdr.seekg (0, ios::beg);
        unsigned long N = (unsigned long) (length/sizeof(unsigned long));

        unsigned long* hdrData = new unsigned long[N];
        fhdr.read ((char*) hdrData,length);
        if(fhdr.bad())
        {
            fhdr.close();
            delete hdrData;
            return false;
        }else{

            data->Fib = new CtlCtlStruct(hdrData[0]);
            //data->Fib->N = hdrData[0];
            data->Fib->fibers = new CtlStruct[data->Fib->N];
            for(unsigned long o=1 ; o< hdrData[0]+1 ; o++)
            {
                (data->Fib->fibers)[o-1].N_element = hdrData[o];
                (data->Fib->fibers)[o-1].elts = new double[3*(hdrData[o])];
                (data->Fib->fibers)[o-1].idx = new unsigned long[hdrData[o]];
                nbBundlePoints += hdrData[o];
            }
        }
        fhdr.close();
        if(nbBundlePoints==0)
        {
            return false;
        }


        //data file
        fstream fsrc;
        fsrc.open (sourceFile.data(), ifstream::in | ifstream::binary ); //
        if(fsrc.is_open())
        {
            //get length of file
            unsigned long lengthB = nbBundlePoints*sizeof(BundlePoint);
            BundlePoint* srcData = new BundlePoint[nbBundlePoints];
            fsrc.read ((char*) srcData,lengthB);
            if(fsrc.bad())
            {
                fsrc.close();
                delete srcData;
                return false;
            }else{
                //for each bundle
                unsigned long offset = 0;
                for(unsigned long o=1 ; o< hdrData[0]+1 ; o++)
                {
                    for(unsigned long g=0 ; g<hdrData[o] ; g++)
                    {
                        (data->Fib->fibers)[o-1].elts[3*g] = srcData[g+offset].x;
                        (data->Fib->fibers)[o-1].elts[3*g+1] = srcData[g+offset].y;
                        (data->Fib->fibers)[o-1].elts[3*g+2] = -srcData[g+offset].z;  ///add a minus on this line to show the human heart in the right way
                        (data->Fib->fibers)[o-1].idx[g] = srcData[g+offset].idx;
                    }
                    offset += hdrData[o];
                }
            }
            delete srcData;
            delete hdrData;
            fsrc.close();





            ///if only a slice is desired
//            Zmin = 15; //28  base
//            Zmax = 20; //38

//            Zmin = 55; //28  middle
//            Zmax = 60; //38

//            Zmin = 85; //28  apex
//            Zmax = -90; //38

//            Zmin = 70; //28
//            Zmax = 75; //38
            Zmin = LOWZ;
            Zmax = -HIGHZ;
//            Zmin = -HIGHZ;
//            Zmax = -LOWZ;
//            double Xmax = 140.0, Xmin=120.0;
            //double Ymax = 118.0, Ymin=86.0;
            double Xmax = 140.0, Xmin=120.0;
            double Ymax = 110.0, Ymin=86.0;
//            double Xmin = -140.0, Xmax = -120.0;
//            double Ymin = -110.0, Ymax = -86.0;
            if(Zmax>Zmin+0.000000001)
            {
                vector< sliceIdx > A;
                sliceIdx B;
                bool wasInSlice;
                double R;
                for(unsigned long i=0 ; i<data->Fib->N ; i++ )
                {
                    B.fibIdx = i;
                    B.nbSubFib = 0;
                    B.idxStart.clear();
                    B.idxStop.clear();
                    wasInSlice = false;
                    for(unsigned long j=0 ; j<(data->Fib->fibers)[i].N_element ; j++ )
                    {
                        R = sqrt( SQR( (data->Fib->fibers)[i].elts[3*j]) + SQR( (data->Fib->fibers)[i].elts[3*j+1]));
                        bool condX = true || (((data->Fib->fibers)[i].elts[3*j]>Xmin) && ((data->Fib->fibers)[i].elts[3*j]<Xmax));
                        bool condY = true || (((data->Fib->fibers)[i].elts[3*j+1]>Ymin) && ((data->Fib->fibers)[i].elts[3*j+1]<Ymax));
                        if( (data->Fib->fibers)[i].elts[3*j+2]<=Zmax &&  (data->Fib->fibers)[i].elts[3*j+2]>=Zmin && R>=LOWR && R<=HIGHR && condX && condY)
                        {
                            if(wasInSlice==false)
                            {
                                B.nbSubFib++;
                                B.idxStart.push_back(j);
                                wasInSlice = true;
                            }
                        }
                        else
                        {
                            if(wasInSlice==true)
                            {
                                B.idxStop.push_back(j);
                                wasInSlice = false;
                            }
                        }
                    }
                    if(wasInSlice==true)
                    {
                        B.idxStop.push_back((data->Fib->fibers)[i].N_element-1);
                    }
                    if(B.nbSubFib!=0)
                    {
                        ///check if there are at least 2 points in each sub-fiber
                        A.push_back(B);
                    }
                }

                ///allocate new fiber container if there are enough fibers, otherwise return false
                unsigned long NFIB = 0;
                for(unsigned int i=0 ; i<A.size() ; i++)
                {
                    NFIB += A.at(i).nbSubFib;
                }

                if(NFIB==0)
                {
                    return false;
                }

                //cout << NFIB << " detected" << endl;

                CtlCtlStruct* tempFib = new CtlCtlStruct(NFIB);
                //tempFib->N = NFIB;
                tempFib->fibers = new CtlStruct[NFIB];
                unsigned long tempIdx = 0;
                //cout << "There are " << A.size() << " fibers involved" << endl;
                for(unsigned int i=0 ; i<A.size() ; i++)
                {
                    //cout << "Fiber: " << i << " there are " << A.at(i).nbSubFib << " sub fibers"<< endl;
                    for(unsigned int j=0 ; j<A.at(i).nbSubFib ; j++)
                    {
                        //cout << "Subfiber: " << j << endl;
                        ///allocate fiber
                        (tempFib->fibers)[tempIdx].N_element = A.at(i).idxStop.at(j) - A.at(i).idxStart.at(j) + 1;
                        (tempFib->fibers)[tempIdx].elts = new double[3*((tempFib->fibers)[tempIdx].N_element)];
                        (tempFib->fibers)[tempIdx].idx = new unsigned long[(tempFib->fibers)[tempIdx].N_element];

                        ///fill in fiber
                        for(unsigned long k=0 ; k<(tempFib->fibers)[tempIdx].N_element ; k++)
                        {
                            (tempFib->fibers)[tempIdx].elts[3*k] = (data->Fib->fibers)[A.at(i).fibIdx].elts[3*(k+A.at(i).idxStart.at(j))];
                            (tempFib->fibers)[tempIdx].elts[3*k+1] = (data->Fib->fibers)[A.at(i).fibIdx].elts[3*(k+A.at(i).idxStart.at(j))+1];
                            (tempFib->fibers)[tempIdx].elts[3*k+2] = (data->Fib->fibers)[A.at(i).fibIdx].elts[3*(k+A.at(i).idxStart.at(j))+2];
                        }

                        ///update idx
                        tempIdx++;
                    }
                }
                A.clear();
                //cout << "Still in :-p" << endl;
                if(false) ///try to measure helix angle
                {
                    FILE* f = fopen("sheetAngle.data","a");
                    for(unsigned long i=0 ; i<tempFib->N ; i++)
                    {
                        double Xmean = 0.5*((tempFib->fibers)[i].elts[0] + (tempFib->fibers)[i].elts[3*((tempFib->fibers)[i].N_element-1)]) ;
                        double Ymean = 0.5*((tempFib->fibers)[i].elts[0+1] + (tempFib->fibers)[i].elts[3*((tempFib->fibers)[i].N_element-1)+1]) ;
                        //double Angle = (3.14159265/2.0) - atan2((tempFib->fibers)[i].elts[0+2] - (tempFib->fibers)[i].elts[3*((tempFib->fibers)[i].N_element-1)+2], (tempFib->fibers)[i].elts[0+1] - (tempFib->fibers)[i].elts[3*((tempFib->fibers)[i].N_element-1)+1]);
                        double dZ = (tempFib->fibers)[i].elts[0+2] - (tempFib->fibers)[i].elts[3*((tempFib->fibers)[i].N_element-1)+2];
                        double dX = (tempFib->fibers)[i].elts[0] - (tempFib->fibers)[i].elts[3*((tempFib->fibers)[i].N_element-1)];
//                        double Angle = atan(dX/dZ);
//                        if(Angle>0)
//                        {
//                            Angle -= 0.5*3.14159265;
//                        }
//                        else
//                        {
//                            Angle += 0.5*3.14159265;
//                        }
                        double Angle = atan(dZ/dX);
                        if(Angle==Angle)
                        {
                            fprintf(f,"%f, %f, %f\n", Xmean, Ymean, Angle);
                        }
                        else
                        {
                            fprintf(f,"%f, %f, %f\n", Xmean, Ymean, 0.5*3.14159265);
                        }
                    }
                    fclose(f);
                }


                ///fill in fibers
                delete data->Fib;
                data->Fib = tempFib;
            }


        }else{
            return false;
        }


    }
    else
    {
        return false;
    }
    return true;
}


bool C_visu_IO::readLargestFiber(string fileNameFib)//, CtlCtlStruct* fibersCtl)
{
    //file names
    string headerFile = fileNameFib + ".fibHDR";
    string sourceFile = fileNameFib + ".fibSRC";

    //
    fstream fhdr;
    fhdr.open (headerFile.data(), ifstream::in | ifstream::binary ); //
    if(fhdr.is_open())
    {
        //
        unsigned long nbBundlePoints = 0;

        //get number of data = nb_bundle+1
        fhdr.seekg (0, ios::end);
        unsigned long length = fhdr.tellg();
        fhdr.seekg (0, ios::beg);
        unsigned long N = (unsigned long) (length/sizeof(unsigned long));

        unsigned long* hdrData = new unsigned long[N];
        fhdr.read ((char*) hdrData,length);
        unsigned long idxLargest = 1;
        if(fhdr.bad())
        {
            fhdr.close();
            delete hdrData;
            return false;
        }else{

            data->Fib = new CtlCtlStruct(1);
            //data->Fib->N = 1;
            data->Fib->fibers = new CtlStruct[data->Fib->N];
            unsigned long NB = 0;
            for(unsigned long o=1 ; o< hdrData[0]+1 ; o++)
            {
                nbBundlePoints += hdrData[o];
                if(hdrData[o]>NB)
                {
                    NB = hdrData[o];
                    idxLargest = o;
                }
            }
            if(nbBundlePoints!=0)
            {
                (data->Fib->fibers)[0].N_element = hdrData[idxLargest];
                (data->Fib->fibers)[0].elts = new double[3*(hdrData[idxLargest])];
                (data->Fib->fibers)[0].idx = new unsigned long[hdrData[idxLargest]];
            }

        }
        fhdr.close();
        if(nbBundlePoints==0)
        {
            return false;
        }


        //data file
        fstream fsrc;
        fsrc.open (sourceFile.data(), ifstream::in | ifstream::binary ); //
        if(fsrc.is_open())
        {
            //get length of file
            unsigned long lengthB = nbBundlePoints*sizeof(BundlePoint);
            BundlePoint* srcData = new BundlePoint[nbBundlePoints];
            fsrc.read ((char*) srcData,lengthB);
            if(fsrc.bad())
            {
                fsrc.close();
                delete srcData;
                return false;
            }else{
                //for each bundle
                unsigned long offset = 0;
                for(unsigned long o=1 ; o< hdrData[0]+1 ; o++)
                {
                    if(o==idxLargest)
                    {
                        for(unsigned long g=0 ; g<hdrData[o] ; g++)
                        {
                            (data->Fib->fibers)[0].elts[3*g] = srcData[g+offset].x;
                            (data->Fib->fibers)[0].elts[3*g+1] = srcData[g+offset].y;
                            (data->Fib->fibers)[0].elts[3*g+2] = srcData[g+offset].z;
                            (data->Fib->fibers)[0].idx[g] = srcData[g+offset].idx;
                        }
                    }
                    offset += hdrData[o];
                }
            }
            delete srcData;
            delete hdrData;
            fsrc.close();


        }else{
            return false;
        }


    }
    else
    {
        return false;
    }
    return true;
}

bool C_visu_IO::loadTensorFromPointData(string fileNamePoints) //this might be useless
{
    fstream f;
    f.open (fileNamePoints.data(), ifstream::in);
    if(f.is_open())
    {
        //get length of file
        f.seekg (0, ios::end);
        unsigned long length = f.tellg();
        f.seekg (0, ios::beg);
        unsigned long nb_point = (unsigned long) (length/sizeof(PointWithNeighbors));

        PointWithNeighbors* m_p = new PointWithNeighbors[nb_point];
        if(m_p!=NULL)
        {
            f.read ((char*) m_p,length);
            if(f.bad())
            {
                f.close();
                delete m_p;
                m_p = NULL;
                return false;
            }
            f.close();

            //convert data into tensor data
            data->T = new TensorData[nb_point];

            getTensorFromRawPoint(m_p, data->T);
            #ifdef EXTRACT_DWI
            fstream fdwi;
            fdwi.open ("dwi_6dir.data", ifstream::out);
            #endif

            for(unsigned long i=1 ; i<nb_point ; i++)
            {
                (m_p)++;
                (data->T)++;
                getTensorFromRawPoint(m_p, data->T);
                #ifdef EXTRACT_DWI
                if(m_p->z<=HIGHZ && m_p->z>=LOWZ)
                {
                    fdwi << m_p->x << ", " << m_p->y << ", " << m_p->z << ", " << getDWIfromDTI(m_p, 1.0/sqrt(2.0), 1.0/sqrt(2.0), 0.0, 1.0) << ", " << getDWIfromDTI(m_p, 1.0/sqrt(2.0), -1.0/sqrt(2.0), 0.0, 1.0) << ", " << getDWIfromDTI(m_p, 0.0, 1.0/sqrt(2.0), 1.0/sqrt(2.0), 1.0) << ", " << getDWIfromDTI(m_p, 0.0, 1.0/sqrt(2.0), -1.0/sqrt(2.0), 1.0) << ", " << getDWIfromDTI(m_p, 1.0/sqrt(2.0), 0.0, 1.0/sqrt(2.0), 1.0) << ", " << getDWIfromDTI(m_p, -1.0/sqrt(2.0), 0.0, 1.0/sqrt(2.0), 1.0) << endl;
                }
                #endif
            }
            (m_p) -= (nb_point-1);
            (data->T) -= (nb_point-1);
            data->nbTensor = nb_point;
            delete m_p;
            #ifdef EXTRACT_DWI
            fdwi.close();
            #endif

        }else
        {
            f.close();
            return false;
        }

    }
    else
    {
        return false;
    }
    return true;
}

#ifdef EXTRACT_DWI
double C_visu_IO::getDWIfromDTI(PointWithNeighbors* p, double gx, double gy, double gz, double b)
{
    double gtDg = gx*(((double) (p->matrix[0]))*gx + ((double) (p->matrix[1]))*gy + ((double) (p->matrix[2]))*gz) +\
                gy*(((double) (p->matrix[1]))*gx + ((double) (p->matrix[3]))*gy + ((double) (p->matrix[4]))*gz) +\
                gz*(((double) (p->matrix[2]))*gx + ((double) (p->matrix[4]))*gy + ((double) (p->matrix[5]))*gz);
    return exp(-b*gtDg);
}
#endif

bool C_visu_IO::getTensorsFromRawPoints(PointWithNeighbors** p, unsigned long nb_point, TensorData** t) //compute the tensor data from raw points
{
    if(t!=NULL && *p!=NULL)
    {
        *t = new TensorData[nb_point];

        getTensorFromRawPoint(*p, *t);
        for(unsigned long i=1 ; i<nb_point ; i++)
        {
            (*p)++;
            (*t)++;
            getTensorFromRawPoint(*p, *t);
        }
        (*p) -= (nb_point-1);
        (*t) -= (nb_point-1);
    }
    else
    {
        *t = NULL;
        return false;
    }

    return true;
}

bool C_visu_IO::getTensorFromRawPoint(PointWithNeighbors* p, TensorData* t) //compute a tensor data from a raw point address
{
    double** matrix = new double*[3];
    for(int i=0 ; i<3 ; i++)
    {
        matrix[i] = new double[3];
    }
    matrix[0][0] = (double) (p->matrix[0]);
    matrix[1][0] = (double) (p->matrix[1]);
    matrix[2][0] = (double) (p->matrix[2]);
    matrix[0][1] = (double) (p->matrix[1]);
    matrix[1][1] = (double) (p->matrix[3]);
    matrix[2][1] = (double) (p->matrix[4]);
    matrix[0][2] = (double) (p->matrix[2]);
    matrix[1][2] = (double) (p->matrix[4]);
    matrix[2][2] = (double) (p->matrix[5]);

    C_toolbox_eigen_sym* eig = new C_toolbox_eigen_sym(matrix, 3, true);


    t->eigenValue[0] = eig->d[0];
    t->eigenValue[1] = eig->d[1];
    t->eigenValue[2] = eig->d[2];

    t->eigenVector1[0] = eig->z[0][0];
    t->eigenVector1[1] = eig->z[1][0];
    t->eigenVector1[2] = eig->z[2][0];

    t->eigenVector2[0] = eig->z[0][1];
    t->eigenVector2[1] = eig->z[1][1];
    t->eigenVector2[2] = eig->z[2][1];

    t->eigenVector3[0] = eig->z[0][2];
    t->eigenVector3[1] = eig->z[1][2];
    t->eigenVector3[2] = eig->z[2][2];

    delete eig;

    double theta = atan2(t->eigenVector1[1], t->eigenVector1[0]);
    double phi = atan2(t->eigenVector1[2], t->eigenVector1[0]*cos(theta)+t->eigenVector1[1]*sin(theta));
    double psi = atan2(t->eigenVector2[0]*cos(theta)*sin(phi) + t->eigenVector2[1]*sin(theta)*sin(phi) + t->eigenVector2[2]*cos(phi), t->eigenVector2[0]*(sin(theta)) + t->eigenVector2[1]*cos(theta));
    t->theta = theta;
    t->phi = phi;
    t->psi = psi;

    t->x = (double) (p->x);
    t->y = (double) (p->y);
    t->z = (double) (p->z);

    for(int i=0 ; i<3 ; i++)
    {
        if(matrix[i]!=NULL)
            delete [] matrix[i];
    }
    if(matrix!=NULL)
        delete [] matrix;

    return true;
}


#ifdef LATER
bool C_visu_IO::loadEdge(std::string fileNameEdge)
{
    fstream f;
    f.open (fileNameEdge.data(), ifstream::in);
    if(f.is_open())
    {
        //get length of file
        f.seekg (0, ios::end);
        unsigned long length = f.tellg();
        f.seekg (0, ios::beg);
        data->nbEdges = (unsigned long) (length/sizeof(SaveEdges));
//        if(e!=NULL)
        //{
            SaveEdges* eTemp = new SaveEdges[data->nbEdges];
            f.read ((char*) eTemp,length);

            if(f.bad())
            {
                f.close();
                delete eTemp;
                data->E = NULL;
                return false;
            }
            f.close();

            //select only connected edges
            vector<unsigned long> idx;
            for(unsigned long i=0 ; i<data->nbEdges ; i++)
            {
                if(eTemp[i].connected)
                {
                    //push_back
                    idx.push_back(i);
                }
            }

            //update number of connected edges
            data->nbEdges = idx.size();

            //allocate array of Saved edges
            data->E = new SaveEdges[data->nbEdges];

            //store connected edges
            for(unsigned long i=0 ; i<data->nbEdges ; i++)
            {
                data->E[i].x1 = eTemp[idx.at(i)].x1;
                data->E[i].y1 = eTemp[idx.at(i)].y1;
                data->E[i].z1 = eTemp[idx.at(i)].z1;
                data->E[i].x2 = eTemp[idx.at(i)].x2;
                data->E[i].y2 = eTemp[idx.at(i)].y2;
                data->E[i].z2 = eTemp[idx.at(i)].z2;
                data->E[i].connected = eTemp[idx.at(i)].connected;
            }
            idx.clear();
            delete eTemp;
        //}

    }
    else
    {
        data->E = NULL;
        return false;
    }
    return true;
}
#endif







///write
void C_visu_IO::save(string filename)
{
    load(filename, M_RAW_FIBERS);


    vtkSmartPointer<vtkPoints> P = vtkPoints::New(VTK_FLOAT);
    vtkSmartPointer<vtkCellArray> lineVect = vtkSmartPointer<vtkCellArray>::New();
    vtkSmartPointer<vtkPolyData> polyDataVect = vtkPolyData::New();
    unsigned long n=0;
    for(unsigned long i=0 ; i<data->Fib->N ; i++ )
    {
        lineVect->InsertNextCell( data->Fib->fibers[i].N_element );
        for(unsigned long j=0 ; j<data->Fib->fibers[i].N_element ; j++)
        {
            P->InsertNextPoint(data->Fib->fibers[i].elts[3*j], data->Fib->fibers[i].elts[3*j+1], data->Fib->fibers[i].elts[3*j+2]);
            lineVect->InsertCellPoint(n);
            n++;
        }
        lineVect->UpdateCellCount( data->Fib->fibers[i].N_element );
    }
    polyDataVect->SetPoints(P);
    polyDataVect->SetLines(lineVect);
    vtkSmartPointer<vtkPolyDataWriter> writer4 = vtkSmartPointer<vtkPolyDataWriter>::New();
    writer4->SetFileName( (filename.append(".vtk")).data() );
    writer4->SetInput( polyDataVect );
    writer4->Write();
    return;
}
