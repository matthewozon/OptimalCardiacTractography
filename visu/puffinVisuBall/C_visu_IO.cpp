#include <C_visu_IO.h>


C_visu_IO::C_visu_IO()
{
    //ctor
    data = new plotData;
    data->rawFibers = 0;
    Zmax = HIGHZ;
    Zmin = LOWZ;
}

C_visu_IO::~C_visu_IO()
{
    //dtor
    delete data;
}

bool C_visu_IO::load(string fileName, long kind)
{
    if(kind==M_TENSOR)
    {
        return loadTensorFromPointData(fileName);
    }
    else if(kind==M_FIBERS || kind==M_RAW_FIBERS)
    {
        if(kind==M_RAW_FIBERS)
        {
            data->rawFibers = 1;
        }
        return readFiber(fileName);
    }
    return false;
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
bool C_visu_IO::readFiber(string fileNameFib)
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
            //data->Fib->fibers = new CtlStruct[data->Fib->N];// done in constructor
            for(unsigned long o=1 ; o< hdrData[0]+1 ; o++)
            {
                (data->Fib->fibers)[o-1].allocateArray(hdrData[o]);
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
                        (data->Fib->fibers)[o-1].elts[3*g+2] = srcData[g+offset].z;
                        (data->Fib->fibers)[o-1].idx[g] = srcData[g+offset].idx;
                    }
                    offset += hdrData[o];
                }
            }
            delete srcData;
            delete hdrData;
            fsrc.close();





            ///if only a slice is desired
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
                        if( (data->Fib->fibers)[i].elts[3*j+2]<=Zmax &&  (data->Fib->fibers)[i].elts[3*j+2]>=Zmin && R>=LOWR && R<=HIGHR)
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


                CtlCtlStruct* tempFib = new CtlCtlStruct(NFIB);
                //tempFib->fibers = new CtlStruct[NFIB];
                unsigned long tempIdx = 0;
                for(unsigned int i=0 ; i<A.size() ; i++)
                {
                    for(unsigned int j=0 ; j<A.at(i).nbSubFib ; j++)
                    {
                        ///allocate fiber
                        (tempFib->fibers)[tempIdx].allocateArray(A.at(i).idxStop.at(j) - A.at(i).idxStart.at(j) + 1);

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

            ///if only some part is de-rised
            if(Zmax>Zmin+0.000000001)
            {
                unsigned long NN=0;
                for(unsigned long i=0 ; i<nb_point ; i++)
                {
                    if(m_p[i].z>Zmin &&  m_p[i].z<Zmax)
                    {
                        NN++;
                    }
                }
                cout << NN << endl;
                PointWithNeighbors* new_m_p = new PointWithNeighbors[NN];
                if(new_m_p==NULL)
                {
                    delete m_p;
                    m_p = NULL;
                    return false;
                }
                unsigned long j=0;
                for(unsigned long i=0 ; i<nb_point ; i++)
                {
                    if(m_p[i].z>Zmin &&  m_p[i].z<Zmax)
                    {
                        new_m_p[j].x = m_p[i].x;
                        new_m_p[j].y = m_p[i].y;
                        new_m_p[j].z = m_p[i].z;
                        new_m_p[j].matrix[0] = m_p[i].matrix[0];
                        new_m_p[j].matrix[1] = m_p[i].matrix[1];
                        new_m_p[j].matrix[2] = m_p[i].matrix[2];
                        new_m_p[j].matrix[3] = m_p[i].matrix[3];
                        new_m_p[j].matrix[4] = m_p[i].matrix[4];
                        new_m_p[j].matrix[5] = m_p[i].matrix[5];
                        j++;
                    }
                    if(j==NN)
                    {
                        i=nb_point;
                    }
                }
                nb_point = NN;
                delete m_p;
                m_p = new_m_p;
            }

            //convert data into tensor data
            data->T = new TensorData[nb_point];

            getTensorFromRawPoint(m_p, data->T);

            for(unsigned long i=1 ; i<nb_point ; i++)
            {
                (m_p)++;
                (data->T)++;
                getTensorFromRawPoint(m_p, data->T);
            }
            (m_p) -= (nb_point-1);
            (data->T) -= (nb_point-1);
            data->nbTensor = nb_point;
            delete m_p;


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

    C_toolbox_eigen_sym* eig = new C_toolbox_eigen_sym((double**) matrix, 3, true);


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

    delete eig;//->~C_eigen_sym();

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
            delete matrix[i];
    }
    if(matrix!=NULL)
        delete matrix;

    return true;
}


bool C_visu_IO::selectFiberDetail(double x, double y, double z, double R)
{
    ///count fibers in range and store their indices
    vector<unsigned long> idx;
    bool dummyCondition;
    for(unsigned long i=0 ; i<data->Fib->N ; i++)
    {
        //if at least one point of the fiber crosses the ball
        dummyCondition = false;
        for(unsigned long j=0 ; j<data->Fib->fibers[i].N_element ; j++)
        {
            if( SQR(data->Fib->fibers[i].elts[3*j+0]-x) + SQR(data->Fib->fibers[i].elts[3*j+1]-y) + SQR(data->Fib->fibers[i].elts[3*j+2]-z) < SQR(R) )
            {
                dummyCondition = true;
            }
        }
        if( dummyCondition )
        {
            idx.push_back(i);
        }
    }
    cout << "there are: " << idx.size() << " fibers remaining" << endl;
    if(idx.size()==0) return false;
    ///create a new fiber structure
    CtlCtlStruct* tempFibers = new CtlCtlStruct(idx.size());
    ///store only fibers in the range
    for(unsigned int i=0 ; i<idx.size() ; i++)
    {
        tempFibers->fibers[i].allocateArray(data->Fib->fibers[idx.at(i)].N_element);
        //copy elements
        for(unsigned long j=0 ; j<tempFibers->fibers[i].N_element ; j++ )
        {
            tempFibers->fibers[i].elts[3*j+0] = data->Fib->fibers[idx.at(i)].elts[3*j+0];
            tempFibers->fibers[i].elts[3*j+1] = data->Fib->fibers[idx.at(i)].elts[3*j+1];
            tempFibers->fibers[i].elts[3*j+2] = data->Fib->fibers[idx.at(i)].elts[3*j+2];
        }
    }
    ///delete old fibers
    delete data->Fib;
    ///replace with new fibers
    data->Fib = tempFibers;
    return true;
}
