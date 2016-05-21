#include <C_streamline.h>

C_streamline::C_streamline(rawData<double>* rawTensors, rawData<double>* rawMask, double wpunct)
{
    //ctor
    m_wpunct = wpunct;
    m_rawTensorsCOmpare = NULL;
    m_rawTensors = rawTensors;
    m_rawMask = rawMask;
}

C_streamline::~C_streamline()
{
    //dtor
}

bool C_streamline::run2(double LEVEL_SAMPLING, double dt, double MINSPEED, double MAXLENGTH)
{
    //cout << "HERE0?" << endl;
    vector< vector<BundlePoint> > V = streamLines(LEVEL_SAMPLING, dt, MINSPEED, MAXLENGTH);
    //cout << "HERE?" << endl;

    ///compute measures
    try
    {
        //cout << "after stream" << endl;
        if(m_rawTensorsCOmpare==NULL) m_rawTensorsCOmpare = m_rawTensors;
        double L=0.0, S=0.0, E=0.0;
        measure(V, &L/**mean length*/, &S/**mean standard deviavtion of length*/, &E/**mean fiber error*/);
        //cout << "after measure" << endl;

        ///write to file
        ostringstream ossdata;
        ossdata << "streamline_result_step_" << dt << "_nbseed_" << LEVEL_SAMPLING << "_minspeed_" << MINSPEED << "_maxlength_" << MAXLENGTH << ".data";
        FILE* f = fopen(ossdata.str().data() ,"a");
        ///step size, nb seed, min speed, maxlength, mean length, std length, mean fiber error
        if(f!=NULL)
        {
            fprintf(f, "%f, %f, %f, %f, %f, %f, %f\n", dt, LEVEL_SAMPLING, MINSPEED, MAXLENGTH, L, S, E);
            fclose(f);
        }
    }
    catch(...)
    {
        cout << "an error occured while computing measures: integration step " << dt << " sampling level: " << LEVEL_SAMPLING << " min speed " << MINSPEED << " max length " << MAXLENGTH << endl;
    }

    ostringstream ossfiber;
    ossfiber << "streamline_fiber_step_" << dt << "_nbseed_" << LEVEL_SAMPLING << "_minspeed_" << MINSPEED << "_maxlength_" << MAXLENGTH;
    return saveFiber(ossfiber.str().data(), V);
}
bool C_streamline::run2(unsigned short NB_SEED, double dt, double MINSPEED, double MAXLENGTH)
{
    vector< vector<BundlePoint> > V = streamLines(NB_SEED, dt, MINSPEED, MAXLENGTH);

    ///compute measures
    try
    {
        //cout << "after stream" << endl;
        if(m_rawTensorsCOmpare==NULL) m_rawTensorsCOmpare = m_rawTensors;
        double L=0.0, S=0.0, E=0.0;
        measure(V, &L/**mean length*/, &S/**mean standard deviavtion of length*/, &E/**mean fiber error*/);
        //cout << "after measure" << endl;

        ///write to file
        ostringstream ossdata;
        ossdata << "streamline_result_step_" << dt << "_nbseed_" << NB_SEED << "_minspeed_" << MINSPEED << "_maxlength_" << MAXLENGTH << ".data";
        FILE* f = fopen(ossdata.str().data() ,"a");
        ///step size, nb seed, min speed, maxlength, mean length, std length, mean fiber error
        if(f!=NULL)
        {
            fprintf(f, "%f, %i, %f, %f, %f, %f, %f\n", dt, NB_SEED, MINSPEED, MAXLENGTH, L, S, E);
            fclose(f);
        }
    }
    catch(...)
    {
        cout << "an error occured while computing measures: integration step " << dt << " nb seed: " << NB_SEED << " min speed " << MINSPEED << " max length " << MAXLENGTH << endl;
    }

    ostringstream ossfiber;
    ossfiber << "streamline_fiber_step_" << dt << "_nbseed_" << NB_SEED << "_minspeed_" << MINSPEED << "_maxlength_" << MAXLENGTH;
    return saveFiber(ossfiber.str().data(), V);
}

///do NOT delete measures
void C_streamline::measure(vector< vector<BundlePoint> > V, double* L/**mean length*/, double* S/**mean standard deviavtion of length*/, double* E/**mean fiber error*/)
{
    double* LL = new double[V.size()];
    double* EE = new double[V.size()];
    for(unsigned int i=0 ; i<V.size() ; i++)
    {
        measure(V.at(i), &(LL[i])/**length*/, &(EE[i])/**fiber error*/);
    }

    *L = 0;
    for(unsigned int i=0 ; i<V.size() ; i++)
    {
        (*L) += LL[i];
    }
    *L = (*L)/((double) V.size());

    *S = 0;
    for(unsigned int i=0 ; i<V.size() ; i++)
    {
        (*S) += SQR(LL[i]-(*L));
    }
    *S = sqrt((*S)/((double) V.size()));

    *E = 0;
    for(unsigned int i=0 ; i<V.size() ; i++)
    {
        (*E) += EE[i];
    }
    *E = (*E)/((double) V.size());
    delete EE;
    delete LL;
    return;
}
void C_streamline::measure(vector<BundlePoint>  v, double* L/**length*/, double* E/**fiber error*/)
{
    ///create interpolate point with all array
    if(v.size()>10)
    {
        ///compute dr (array of length v.size()-1)
        double* dr = new double[3];
        *L = 0;
        *E = 0;
        double Tdr = 0;
        C_toolbox_eigen_sym* toolEIG;
        //check if all points are in volume
        for(unsigned int i=0 ; i<v.size() ; i++)
        {
            if(!isPartOfVolume(v.at(i).x, v.at(i).y, v.at(i).z))
            {
                //remove point from fiber
//                cout << "point #" << i+1 << " is not in volume. There are " << v.size() << " points in fiber" << endl;
                v.erase(v.begin()+i);
                i--;
//                cout << "point #" << i+1 << " is not in volume. There are " << v.size() << " points in fiber" << endl;
            }
        }
        for(unsigned int i=0 ; i<v.size()-1 ; i++)
        {
            //get coordinates of the point in raw data
            unsigned long idxX = (unsigned long) (v.at(i).x/m_rawTensors->pixDimX);
            unsigned long idxY = (unsigned long) (v.at(i).y/m_rawTensors->pixDimY);
            unsigned long idxZ = (unsigned long) (v.at(i).z/m_rawTensors->pixDimZ);

            //get position variation
            dr[0] = v.at(i+1).x - v.at(i).x;
            dr[1] = v.at(i+1).y - v.at(i).y;
            dr[2] = v.at(i+1).z - v.at(i).z;
            *L += sqrt(SQR(dr[0]) + SQR(dr[1]) + SQR(dr[2]));

            //compute tensor
            Tdr = sqrt( SQR(m_rawTensors->raw4D[idxX][idxY][idxZ][0]*dr[0] + m_rawTensors->raw4D[idxX][idxY][idxZ][1]*dr[1] + m_rawTensors->raw4D[idxX][idxY][idxZ][2]*dr[2]) +\
                       SQR(m_rawTensors->raw4D[idxX][idxY][idxZ][1]*dr[0] + m_rawTensors->raw4D[idxX][idxY][idxZ][3]*dr[1] + m_rawTensors->raw4D[idxX][idxY][idxZ][4]*dr[2]) +\
                       SQR(m_rawTensors->raw4D[idxX][idxY][idxZ][2]*dr[0] + m_rawTensors->raw4D[idxX][idxY][idxZ][4]*dr[1] + m_rawTensors->raw4D[idxX][idxY][idxZ][5]*dr[2]) );
            toolEIG = new C_toolbox_eigen_sym(m_rawTensors->raw4D[idxX][idxY][idxZ], false);
            *E += (1 - (Tdr/(toolEIG->d[0]*sqrt(SQR(dr[0]) + SQR(dr[1]) + SQR(dr[2])))) );
            delete toolEIG;
        }
        *E = (*E)/(((double) v.size()-1)*(*L));

        delete dr;
    }
    else
    {
        *L = 0;
        *E = 0;
    }
    return;
}

bool C_streamline::saveFiber(std::string fileNameFib, vector< vector<BundlePoint> > fibers)
{
    //file names
    string headerFile = fileNameFib + ".fibHDR";
    string sourceFile = fileNameFib + ".fibSRC";

    //filestreams
    fstream filestrSRC, filestrHDR;
    filestrHDR.open (headerFile.data(), ofstream::binary | ofstream::out); //should check if ok
    filestrSRC.open (sourceFile.data(), ofstream::binary | ofstream::out); //should check if ok

    //get number of bundles
    unsigned long NB = fibers.size();
    filestrHDR.write((char*) &NB, sizeof(unsigned long)); //write it in header file

    //declare a buffer to store and write
    BundlePoint* b_fibers;

    unsigned long maxN;
    //for each bundle
    for(unsigned long u=0 ; u<NB ; u++)
    {
        //get the size of the current bundle and save it in header file
        maxN = fibers.at(u).size();
        filestrHDR.write((char*) &maxN, sizeof(unsigned long));

        //allocate bundle
        b_fibers = new BundlePoint[maxN]; //should check if allocation doesn't fail

        //fill in the allocated bundle
        for(unsigned long r=0 ; r<maxN ; r++)
        {
            //fill fibers
            b_fibers[r].x = fibers.at(u).at(r).x;
            b_fibers[r].y = fibers.at(u).at(r).y;
            b_fibers[r].z = fibers.at(u).at(r).z;
            b_fibers[r].idx = fibers.at(u).at(r).idx;
        }

        //write data in source file
        filestrSRC.write((char*) b_fibers, maxN*sizeof(BundlePoint));
        delete b_fibers;
    }

    //close files
    filestrSRC.close();
    filestrHDR.close();


    return true;
}



bool C_streamline::isPartOfVolume(double x, double y, double z)
{
    if(x<0.0 || y<0.0 || z<0.0) return false;
    if(x>=m_rawMask->pixDimX*m_rawMask->DimX || y>=m_rawMask->pixDimY*m_rawMask->DimY || z>=m_rawMask->pixDimZ*m_rawMask->DimZ) return false;
    long i = (long) (x/m_rawMask->pixDimX);
    if(i<0 || i>=(long) m_rawMask->DimX) return false;
    long j = (long) (y/m_rawMask->pixDimY);
    if(j<0 || j>=(long) m_rawMask->DimY) return false;
    long k = (long) (z/m_rawMask->pixDimZ);
    if(k<0 || k>=(long) m_rawMask->DimZ) return false;
    return (m_rawMask->raw3D[i][j][k]>0.5);
}
double* C_streamline::F(double x, double y, double z, double t) ///dr(x,y,z,t)/dt = F(x,y,z,t)
{
    double* v = new double[3];
    if(isPartOfVolume(x,y,z))
    {
        C_toolbox_interpolate* INTERPOLATE = new C_toolbox_interpolate(m_rawTensors, m_rawMask); //, double r=1.0
        C_point_array_data* p = new C_point_array_data(1 /**one point*/, 3, 3, false/**do not init all array*/);
        p->pointCoordinate[0][0] = x;
        p->pointCoordinate[0][1] = y;
        p->pointCoordinate[0][2] = z;

        ///compute tensor, main eigen vector and FA at point (x,y,z)
        //cout << "before interpolation at point " << x << " " << y << " " << z << "idx: " << (unsigned long) (x/m_rawMask->pixDimX) << " " << (unsigned long) (y/m_rawMask->pixDimY) << " " << (unsigned long) (z/m_rawMask->pixDimZ) << endl;
        INTERPOLATE->run(p, NEAREST_NEIGHBOR /*LINEAR_INTERPOLATION GAUSSIAN_PHI*/, TENSOR, false);
        //cout << "after interpolation" << endl;

        ///compute the return value
        v[0] = (p->FA[0])*(p->vectorInterpolated[0][0]);
        v[1] = (p->FA[0])*(p->vectorInterpolated[0][1]);
        v[2] = (p->FA[0])*(p->vectorInterpolated[0][2]);
        //cout << v[0] << " " << v[1] << " " << v[2] << " vector in field"<< endl;
        delete p;
        delete INTERPOLATE;
    }
    else
    {
        v[0] = 0.0;
        v[1] = 0.0;
        v[2] = 0.0;
        //cout << "not in volume: SPEED=0" << endl;
    }

    return v;
}
double* C_streamline::EULERstep(double x, double y, double z, double t, double dt, bool reverse)
{
    double* dr = F(x,y,z,t);
    if(reverse)
    {
        dr[0] = -dt*dr[0];
        dr[1] = -dt*dr[1];
        dr[2] = -dt*dr[2];
    }
    else
    {
        dr[0] = dt*dr[0];
        dr[1] = dt*dr[1];
        dr[2] = dt*dr[2];
    }
    return dr;
}
double* C_streamline::RK4step(double x, double y, double z, double t, double dt, bool reverse)
{
    double* K1 = RK4K1(x, y, z, t, dt);
    if(reverse) ///reverse the global direction
    {
        K1[0] = -K1[0];
        K1[1] = -K1[1];
        K1[2] = -K1[2];
    }
    double* K2 = RK4K2(x, y, z, K1[0], K1[1], K1[2], t, dt);
    if(K2[0]*K1[0] + K2[1]*K1[1] + K2[2]*K1[2]<0) ///assume that we should go in almost the same direction
    {
        K2[0] = -K2[0];
        K2[1] = -K2[1];
        K2[2] = -K2[2];
    }
    double* K3 = RK4K3(x, y, z, K2[0], K2[1], K2[2], t, dt);
    if(K3[0]*K2[0] + K3[1]*K2[1] + K3[2]*K2[2]<0) ///assume that we should go in almost the same direction
    {
        K3[0] = -K3[0];
        K3[1] = -K3[1];
        K3[2] = -K3[2];
    }
    double* K4 = RK4K4(x, y, z, K3[0], K3[1], K3[2], t, dt);
    if(K4[0]*K3[0] + K4[1]*K3[1] + K4[2]*K3[2]<0) ///assume that we should go in almost the same direction
    {
        K4[0] = -K4[0];
        K4[1] = -K4[1];
        K4[2] = -K4[2];
    }
    double* dr = new double[3];
    dr[0] = (1.0/6.0)*(K1[0] + 2.0*K2[0] + 2*K3[0] + K4[0]);
    dr[1] = (1.0/6.0)*(K1[1] + 2.0*K2[1] + 2*K3[1] + K4[1]);
    dr[2] = (1.0/6.0)*(K1[2] + 2.0*K2[2] + 2*K3[2] + K4[2]);
    delete K1;
    delete K2;
    delete K3;
    delete K4;
    return dr;
}
double* C_streamline::RK4K1(double x, double y, double z, double t, double dt)
{
    double* K1 = F(x,y,z,t);
    K1[0] *= dt;
    K1[1] *= dt;
    K1[2] *= dt;
    return K1;
}
double* C_streamline::RK4K2(double x, double y, double z, double k1x, double k1y, double k1z, double t, double dt)
{
    double* K2 = F(x+0.5*k1x,y+0.5*k1y,z+0.5*k1z,t+0.5*dt);
    K2[0] *= dt;
    K2[1] *= dt;
    K2[2] *= dt;
    return K2;
}
double* C_streamline::RK4K3(double x, double y, double z, double k2x, double k2y, double k2z, double t, double dt)
{
    double* K3 = F(x+0.5*k2x,y+0.5*k2y,z+0.5*k2z,t+0.5*dt);
    K3[0] *= dt;
    K3[1] *= dt;
    K3[2] *= dt;
    return K3;
}
double* C_streamline::RK4K4(double x, double y, double z, double k3x, double k3y, double k3z, double t, double dt)
{
    double* K4 = F(x+k3x,y+k3y,z+k3z,t+dt);
    K4[0] *= dt;
    K4[1] *= dt;
    K4[2] *= dt;
    return K4;
}

vector<BundlePoint> C_streamline::streamLine(double x0, double y0, double z0, double dt, bool reverse, double MINSPEED, double MAXLENGTH)
{
    vector< BundlePoint > V;
    BundlePoint v;
    v.x = x0;
    v.y = y0;
    v.z = z0;
    v.idx = 0;
    V.push_back(v);
    double t = 0;
    double length = 0;
    double x=x0, y=y0, z=z0;
    unsigned long idx = 0;

    ///init
    double speed;
    //cout << "coucou" << endl;
    double* dr = RK4step(x0, y0, z0, t, dt, reverse);
    //double* dr = EULERstep(x, y, z, t, dt, reverse);
    //cout << "bye" << endl;
    double* dr2;
    t += dt;
    length += sqrt(SQR(dr[0]) + SQR(dr[1]) + SQR(dr[2]));
    speed = length/dt;
    x += dr[0];
    y += dr[1];
    z += dr[2];
    idx++;
    v.x = x;
    v.y = y;
    v.z = z;
    v.idx = idx;
    V.push_back(v);

    ///iter
    while(length<MAXLENGTH)
    {
        //cout << length << endl;
        dr2 = RK4step(x, y, z, t, dt, false);
        //dr2 = EULERstep(x, y, z, t, dt, false);
        //cout << "should change sens? " << dr[0]*dr2[0]+dr[1]*dr2[1]+dr[2]*dr2[2] << endl;
        if(dr[0]*dr2[0]+dr[1]*dr2[1]+dr[2]*dr2[2]<0) /// if we're going in the same direction
        {
            //cout << "change sens " << x << " " << y << " " << z << endl;
            //if(!isPartOfVolume(x,y,z)) cout << "How is that even possible?" << endl;
            delete dr2;
            dr2 = RK4step(x, y, z, t, dt, true);
            //dr2 = EULERstep(x, y, z, t, dt, true);
        }
        //cout << dr2 << endl;
        //cout << dr2[0] << " " << dr2[1] << " " << dr2[2] << endl;

        ///stopping criteria
        speed = sqrt(SQR(dr2[0]) + SQR(dr2[1]) + SQR(dr2[2]))/dt;
        //cout << "speed: " << speed << endl;
        if(speed<MINSPEED)
        {
            //cout << "min speed" << endl;
            length = 2*MAXLENGTH;
        }
        else if(acos((dr[0]*dr2[0]+dr[1]*dr2[1]+dr[2]*dr2[2])/(sqrt(SQR(dr[0]) + SQR(dr[1]) + SQR(dr[2]))*sqrt(SQR(dr2[0]) + SQR(dr2[1]) + SQR(dr2[2]))))>(Pi/2.5) )
        {
            //cout << "angle" << endl;
            length = 2*MAXLENGTH;
        }
        else if(!isPartOfVolume(x+dr2[0], y+dr2[1], z+dr2[2]))
        {
            //cout << "not in mask" << endl;
            length = 2*MAXLENGTH;
        }
        else
        {
            //cout << "going on" << endl;
            t += dt;
            dr[0] = dr2[0];
            dr[1] = dr2[1];
            dr[2] = dr2[2];
            x += dr[0];
            y += dr[1];
            z += dr[2];
            idx++;
            v.x = x;
            v.y = y;
            v.z = z;
            v.idx = idx;
            V.push_back(v);

        }
        length += dt*speed;
        delete dr2;
    }
    //cout << "end fiber" << endl;
    delete dr;
    return V;
}

vector<BundlePoint> C_streamline::streamLine(double x0, double y0, double z0, double dt, double MINSPEED, double MAXLENGTH)
{
    //cout << "coucou 1" << endl;
    vector<BundlePoint> V1 = streamLine(x0, y0, z0, dt, false, MINSPEED, MAXLENGTH);
    vector<BundlePoint> V2 = streamLine(x0, y0, z0, dt, true, MINSPEED, MAXLENGTH);

    ///merge both streamlines

    ///reverse V1
    vector<BundlePoint> V;
    for(unsigned int i=0 ; i<V1.size() ; i++)
    {
        V.push_back(V1.at(V1.size()-i-1));
    }
    ///cat V1 wiht V2\{first element}
    for(unsigned int i=1 ; i<V2.size() ; i++)
    {
        V.push_back(V2.at(i));
    }
    return V;
}
vector< vector<BundlePoint> > C_streamline::streamLines(double LEVEL_SAMPLING, double dt, double MINSPEED, double MAXLENGTH)
{
    //list all points in mask
    unsigned long idxX, idxY, idxZ;
    vector< vector<BundlePoint> > V;
    //cout << m_rawTensors << endl;
    //cout << m_rawTensors->DimX << " " <<m_rawTensors->DimY << " " << m_rawTensors->DimZ << " " << m_rawTensors->DimT << endl;
    //cout << m_rawTensors->pixDimX << " " <<m_rawTensors->pixDimY << " " << m_rawTensors->pixDimZ << " " << m_rawTensors->pixDimT << endl;
    for(double x=0.5*m_rawTensors->pixDimX ; x<m_rawTensors->DimX*m_rawTensors->pixDimX ; x+= LEVEL_SAMPLING*m_rawTensors->pixDimX)
    {
        idxX = (unsigned long) floor(x/m_rawTensors->pixDimX);
        for(double y=0.5*m_rawTensors->pixDimY ; y<m_rawTensors->DimY*m_rawTensors->pixDimY ; y+= LEVEL_SAMPLING*m_rawTensors->pixDimY)
        {
            idxY = (unsigned long) floor(y/m_rawTensors->pixDimY);
            for(double z=0.5*m_rawTensors->pixDimZ ; z<m_rawTensors->DimZ*m_rawTensors->pixDimZ ; z+= LEVEL_SAMPLING*m_rawTensors->pixDimZ)
            {
                idxZ = (unsigned long) floor(z/m_rawTensors->pixDimZ);
                if(idxX<m_rawTensors->DimX && idxY<m_rawTensors->DimY && idxZ<m_rawTensors->DimZ)
                {
                    //cout << idxX << " " << idxY << " " << idxZ << endl;
                    try
                    {
                        if(m_rawMask->raw3D[idxX][idxY][idxZ]>0.5)
                        {
                            vector<BundlePoint> v = streamLine(x, y, z, dt, MINSPEED, MAXLENGTH);
                            if(v.size()>2.10)
                            {
                                V.push_back(v);
                            }
                            else
                            {
                                v.clear();
                            }
                        }
                    }
                    catch(char const* A)
                    {
                        cout << A << " in streamLines(double LEVEL_SAMPLING, double dt, double MINSPEED, double MAXLENGTH), but will continue streamlining." << endl;
                    }

                }

            }
        }
    }
    return V;
}

vector< vector<BundlePoint> > C_streamline::streamLines(unsigned short NB_SEED, double dt, double MINSPEED, double MAXLENGTH) //to be modified and multithread
{
//    cout << "I am here and pointers are " << m_rawMask << " and " << m_rawTensors << endl;
//    cout << "numdim mask " << m_rawMask->numDim << endl;
//    cout << "pointer raw3D mask " << m_rawMask->raw3D << endl;
//    cout << "dim mask " << m_rawMask->DimX << " " << m_rawMask->DimY << " " << m_rawMask->DimZ << endl;
//    cout << "pix dim mask " << m_rawMask->pixDimX << " " << m_rawMask->pixDimY << " " << m_rawMask->pixDimZ << endl;
//    cout << "numdim tensor " << m_rawTensors->numDim << endl;
//    cout << "pointer raw4D tensor " << m_rawTensors->raw4D << endl;
//    cout << "dim tensor " << m_rawTensors->DimX << " " << m_rawTensors->DimY << " " << m_rawTensors->DimZ << endl;
//    cout << "pix dim tensor " << m_rawTensors->pixDimX << " " << m_rawTensors->pixDimY << " " << m_rawTensors->pixDimZ << endl;
    //srand ( time(NULL) );
    //list all points in mask
    vector<unsigned long> idxX;
    vector<unsigned long> idxY;
    vector<unsigned long> idxZ;
    for(unsigned long i=0 ; i<m_rawMask->DimX ; i++)
    {
        for(unsigned long j=0 ; j<m_rawMask->DimY ; j++)
        {
            for(unsigned long k=0 ; k<m_rawMask->DimZ ; k++)
            {
                if(m_rawMask->raw3D[i][j][k]>0.5)
                {
                    idxX.push_back(i);
                    idxY.push_back(j);
                    idxZ.push_back(k);
                }
            }
        }
    }
    double x, y, z;
    unsigned long idx;
    nbPoint = idxX.size();
    vector< vector<BundlePoint> > V;
    C_toolbox_rand2* toolRAND = new C_toolbox_rand2(time(NULL));
    for(unsigned short i=0 ; i<NB_SEED ; i++)
    {
        idx = ((unsigned long) (nbPoint*toolRAND->doub()));
        x = ((double) idxX.at(idx))*m_rawTensors->pixDimX + m_rawTensors->pixDimX*( 0.5 - 0.0*toolRAND->doub() );
        if(x<0.0) x = 0.0;
        y = ((double) idxY.at(idx))*m_rawTensors->pixDimY + m_rawTensors->pixDimY*( 0.5 - 0.0*toolRAND->doub() );
        if(y<0.0) y = 0.0;
        z = ((double) idxZ.at(idx))*m_rawTensors->pixDimZ + m_rawTensors->pixDimZ*( 0.5 - 0.0*toolRAND->doub() );
        if(z<0.0) z = 0.0;
//        x = ((double) idxX.at(idx))*m_rawTensors->pixDimX + m_rawTensors->pixDimX*( 0.5 - (((double) rand())/((double) RAND_MAX+1)) );
//        y = ((double) idxY.at(idx))*m_rawTensors->pixDimY + m_rawTensors->pixDimY*( 0.5 - (((double) rand())/((double) RAND_MAX+1)) );
//        z = ((double) idxZ.at(idx))*m_rawTensors->pixDimZ + m_rawTensors->pixDimZ*( 0.5 - (((double) rand())/((double) RAND_MAX+1)) );
//        cout << "seed: " << i << " point: " << idx << " " << x << " " << y << " " << z << " mask at point " << idxX.at(idx) << " " << idxY.at(idx) << " " << idxZ.at(idx) << " "<< m_rawMask->raw3D[idxX.at(idx)][idxY.at(idx)][idxZ.at(idx)] << endl;
        vector<BundlePoint> v = streamLine(x, y, z, dt, MINSPEED, MAXLENGTH);
        if(v.size()>10)
        {
            V.push_back(v);
        }
        else
        {
            v.clear();
        }
    }
    delete toolRAND;
    return V;
}

