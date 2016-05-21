#include <C_load_fiber.h>

C_load_fiber::C_load_fiber()
{
    //ctor
}

C_load_fiber::~C_load_fiber()
{
    //dtor
}

CtlCtlStruct* C_load_fiber::readFiber(string fileNameFib, bool rawData)//, CtlCtlStruct* fibersCtl)
{
    //file names
    string headerFile = fileNameFib;/// + ".fibHDR";
    string sourceFile = fileNameFib;/// + ".fibSRC";
    if(fileNameFib.find(".fibHDR")==fileNameFib.npos && fileNameFib.find(".fibSRC")==fileNameFib.npos && fileNameFib.find(".fib")==fileNameFib.npos)
    {
        ///means there is no ext
        headerFile += ".fibHDR";
        sourceFile += ".fibSRC";
    }
    else
    {
        if(fileNameFib.find(".fibHDR")!=fileNameFib.npos)
        {
            sourceFile.replace(fileNameFib.find("HDR"), 3, "SRC");
        }
        else if(fileNameFib.find(".fibSRC")!=fileNameFib.npos)
        {
            headerFile.replace(fileNameFib.find("SRC"), 3, "HDR");
        }
        else if(fileNameFib.find(".fib")!=fileNameFib.npos)
        {
            headerFile.append("HDR");
            sourceFile.append("SRC");
        }
        else
        {
            return NULL;
        }
    }

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

        if(N==0)
        {
            fhdr.close();
            return NULL;
        }

        unsigned long* hdrData = new unsigned long[N];
        fhdr.read ((char*) hdrData,length);
        CtlCtlStruct* Fib;
        if(fhdr.bad())
        {
            fhdr.close();
            delete hdrData;
            return NULL;
        }else{
            if(hdrData[0]==0)
            {
                fhdr.close();
                delete hdrData;
                return NULL;
            }
            Fib = new CtlCtlStruct(hdrData[0]);
            for(unsigned long o=1 ; o< hdrData[0]+1 ; o++)
            {
                Fib->fibers[o-1].allocateArray(hdrData[o]);
                nbBundlePoints += hdrData[o];
            }
        }
        fhdr.close();
        if(nbBundlePoints==0)
        {
            delete Fib;//->~CtlCtlStruct();
            delete hdrData;
            return NULL;
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
                delete Fib;//->~CtlCtlStruct();
                delete hdrData;
                return NULL;
            }else{
                //for each bundle
                unsigned long offset = 0;
                for(unsigned long o=1 ; o< hdrData[0]+1 ; o++)
                {
                    for(unsigned long g=0 ; g<hdrData[o] ; g++)
                    {
                        (Fib->fibers)[o-1].elts[3*g] = srcData[g+offset].x;
                        (Fib->fibers)[o-1].elts[3*g+1] = srcData[g+offset].y;
                        (Fib->fibers)[o-1].elts[3*g+2] = srcData[g+offset].z;
                        (Fib->fibers)[o-1].idx[g] = srcData[g+offset].idx;
                    }
                    offset += hdrData[o];
                }
            }
            delete srcData;
            delete hdrData;
            fsrc.close();

            if(rawData)
            {
                return Fib;
            }
            else
            {
                //compute B-splines
                CtlCtlStruct* temp;// = new CtlCtlStruct*;
                C_toolbox_spline* s = new C_toolbox_spline();
                s->computeFiberSpline(Fib, &temp, 11);
                delete s;//->~C_spline();
                delete Fib;//->~CtlCtlStruct();
                return temp;
            }




        }else{
            delete Fib;
            delete hdrData;
            return NULL;
        }


    }
    else
    {
        return NULL;
    }
    return NULL;
}

CtlCtlStruct* C_load_fiber::readFiber(string fileNameFib)///, bool rawData)//, CtlCtlStruct* fibersCtl)
{
    //file names
    string headerFile = fileNameFib;/// + ".fibHDR";
    string sourceFile = fileNameFib;/// + ".fibSRC";
    if(fileNameFib.find(".fibHDR")==fileNameFib.npos && fileNameFib.find(".fibSRC")==fileNameFib.npos && fileNameFib.find(".fib")==fileNameFib.npos)
    {
        ///means there is no ext
        headerFile += ".fibHDR";
        sourceFile += ".fibSRC";
    }
    else
    {
        if(fileNameFib.find(".fibHDR")!=fileNameFib.npos)
        {
            sourceFile.replace(fileNameFib.find("HDR"), 3, "SRC");
        }
        else if(fileNameFib.find(".fibSRC")!=fileNameFib.npos)
        {
            headerFile.replace(fileNameFib.find("SRC"), 3, "HDR");
        }
        else if(fileNameFib.find(".fib")!=fileNameFib.npos)
        {
            headerFile.append("HDR");
            sourceFile.append("SRC");
        }
        else
        {
            return NULL;
        }
    }

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
        CtlCtlStruct* Fib;
        if(fhdr.bad())
        {
            fhdr.close();
            delete hdrData;
            return NULL;
        }else{

            Fib = new CtlCtlStruct(hdrData[0]);
            for(unsigned long o=1 ; o< hdrData[0]+1 ; o++)
            {
                Fib->fibers[o-1].allocateArray(hdrData[o]);
                nbBundlePoints += hdrData[o];
            }
        }
        fhdr.close();
        if(nbBundlePoints==0)
        {
            delete Fib;
            delete hdrData;
            return NULL;
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
                delete Fib;
                delete hdrData;
                return NULL;
            }else{
                //for each bundle
                unsigned long offset = 0;
                for(unsigned long o=1 ; o< hdrData[0]+1 ; o++)
                {
                    for(unsigned long g=0 ; g<hdrData[o] ; g++)
                    {
                        (Fib->fibers)[o-1].elts[3*g] = (long double) srcData[g+offset].x;
                        (Fib->fibers)[o-1].elts[3*g+1] = (long double) srcData[g+offset].y;
                        (Fib->fibers)[o-1].elts[3*g+2] = (long double) srcData[g+offset].z;
                        (Fib->fibers)[o-1].idx[g] = srcData[g+offset].idx;
                    }
                    offset += hdrData[o];
                }
            }
            delete srcData;
            delete hdrData;
            fsrc.close();

//            if(rawData)
//            {
                return Fib;
//            }
//            else
//            {
//                //compute B-splines
//                CtlCtlStruct* temp;// = new CtlCtlStruct*;
//                C_spline* s = new C_spline();
//                s->computeFiberSpline(Fib, &temp, 11);
//                s->~C_spline();
//                Fib->~CtlCtlStruct();
//                return temp;
//            }




        }else{
            delete Fib;//->~CtlCtlStruct();
            delete hdrData;
            return NULL;
        }


    }
    else
    {
        return NULL;
    }
    return NULL;
}



