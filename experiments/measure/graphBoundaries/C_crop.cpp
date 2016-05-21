#include <C_crop.h>

#ifdef CROP
C_crop::C_crop(floating Rmin, floating Rmax, floating Zmin, floating Zmax):C_measure()
{
    //ctor
    m_Rmin = Rmin;
    m_Rmax = Rmax;
    m_Zmin = Zmin;
    m_Zmax = Zmax;
}

C_crop::~C_crop()
{
    //dtor
}

void C_crop::cropHelixFiber(string fileNameSaveEdges)
{
    string fileNameFib = fileNameSaveEdges.substr(0, fileNameSaveEdges.length()-8);
    fileNameFib.append("_crop");
    unsigned long nb_edge;
    SaveEdges* e = readSavedEdges(fileNameSaveEdges, &nb_edge);
    C_graph* G = getGraphFromSaveEdges(e, nb_edge);
    delete [] e;

    //extract fibers from graph
    #ifdef NO_ACUTE
        C_bundle* s = new C_bundle(G, true);
    #else
        C_bundle* s = new C_bundle(G, false);
    #endif
    //save fibers
    saveFiber(fileNameFib, s->FIBERS, G);
    delete G;
    delete s;
    return;
}


SaveEdges* C_crop::readSavedEdges(std::string fileNameEdge, unsigned long* nb_edge)
{
    fstream f;
    f.open (fileNameEdge.data(), ifstream::in);
    if(f.is_open())
    {
        //get length of file
        f.seekg (0, ios::end);
        unsigned long length = f.tellg();
        f.seekg (0, ios::beg);
        *nb_edge = (unsigned long) (length/sizeof(SaveEdges));

        SaveEdges* eTemp = new SaveEdges[*nb_edge];
        f.read ((char*) eTemp,length);

        if(f.bad())
        {
            f.close();
            delete [] eTemp;
            return NULL;
        }
        f.close();

        //select only connected edges
        vector<unsigned long> idx;
        floating R1, R2;
        for(unsigned long i=0 ; i<*nb_edge ; i++)
        {
            if(eTemp[i].connected)
            {
                ///if the whole edge is inside the crop area
                if(eTemp[i].z1<=m_Zmax && eTemp[i].z2<=m_Zmax && eTemp[i].z1>=m_Zmin && eTemp[i].z2>=m_Zmin)
                {
                    R1 = sqrt(SQR(eTemp[i].x1) + SQR(eTemp[i].y1));
                    R2 = sqrt(SQR(eTemp[i].x2) + SQR(eTemp[i].y2));
                    if(R1<=m_Rmax && R2<=m_Rmax && R1>=m_Rmin && R2>=m_Rmin)
                    {
                        //push_back
                        idx.push_back(i);
                    }

                }
            }
        }

        //update number of connected edges
        *nb_edge = idx.size();

        //allocate array of Saved edges
        SaveEdges* e = new SaveEdges[*nb_edge];

        //store connected edges
        for(unsigned long i=0 ; i<*nb_edge ; i++)
        {
            e[i].x1 = eTemp[idx.at(i)].x1;
            e[i].y1 = eTemp[idx.at(i)].y1;
            e[i].z1 = eTemp[idx.at(i)].z1;

            e[i].x2 = eTemp[idx.at(i)].x2;
            e[i].y2 = eTemp[idx.at(i)].y2;
            e[i].z2 = eTemp[idx.at(i)].z2;

            e[i].connected = 1;
        }
        idx.clear();
        delete [] eTemp;
        return e;

    }
    else
    {
        return NULL;
    }
}

#endif
