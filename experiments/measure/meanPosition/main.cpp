#include <iostream>
#include <sstream>
#include <vector>
#include <stdio.h>
#include <C_graph.h>
#include <C_measure.h>
#include <C_parametric_comparison.h>
#include <C_ReadAnalyze75.h>


using namespace std;

string getAlpha(string filename);
string getBeta(string filename);
string getSigma(string filename);
string getAmplitude(string filename);
string getSNR(string filename);
string getProbaMove(string filename);
string getStdMove(string filename);
string getPhi0(string filename);
string getPhi1(string filename);

//#define FIBER_SIMILARITY
//#define MEAN_FIBER_ERROR
#define MEAN_POSITION
//#define CURVATURE
//#define GRAPH

int main(int argc, char **argv)
{
#ifdef GRAPH
    cout << "Command must be run as: " << argv[0] << " mask.hdr<string>" << endl;
    if(argc==2)
    {
        C_ReadAnalyze75 toolANA;
        rawData<double>* roiData = toolANA.getRawData(argv[1]);
        C_graph G(roiData);
        cout << argv[1] << " has " << G.getLenPoints() << " points and " << G.getLenEdges() << " edges" << endl;
        //get min and max position in ROI
        Point* P = G.getPoints();
        double minX = P[0].x, maxX = P[0].x;
        double minY = P[0].y, maxY = P[0].y;
        double minZ = P[0].z, maxZ = P[0].z;
        double meanX = P[0].x, meanY = P[0].y, meanZ = P[0].z;
        for(unsigned long i=1 ; i<G.getLenPoints() ; i++)
        {
            meanX += P[i].x;
            meanY += P[i].y;
            meanZ += P[i].z;

            if(minX>P[i].x) minX = P[i].x;
            if(maxX<P[i].x) maxX = P[i].x;

            if(minY>P[i].y) minY = P[i].y;
            if(maxY<P[i].y) maxY = P[i].y;

            if(minZ>P[i].y) minZ = P[i].z;
            if(maxZ<P[i].y) maxZ = P[i].z;
        }
        cout << "min (x,y,z) (" << minX << ", " << minY << ", " << minZ << ")" << endl;
        cout << "max (x,y,z) (" << maxX << ", " << maxY << ", " << maxZ << ")" << endl;
        cout << "mean (x,y,z) (" << meanX << ", " << meanY << ", " << meanZ << ")" << endl;
    }

#endif
#ifdef MEAN_POSITION
    cout << argv[1] << endl;
    C_load_fiber myLoadObject;
    CtlCtlStruct* FIBERS = myLoadObject.readFiber(argv[1],false);
    double xc = 0.0;
    double yc = 0.0;
    double zc = 0.0;
    double minx=FIBERS->fibers[0].elts[0], miny=FIBERS->fibers[0].elts[1], minz=FIBERS->fibers[0].elts[2];
    double maxx=FIBERS->fibers[0].elts[0], maxy=FIBERS->fibers[0].elts[1], maxz=FIBERS->fibers[0].elts[2];
    unsigned long nbPoints = 0;
    for(unsigned long n=0 ; n<FIBERS->N ; n++ )
    {
        for(unsigned long i=0 ; i<FIBERS->fibers[n].N_element ; i++ )
        {
            nbPoints++;
            xc += FIBERS->fibers[n].elts[3*i];
            yc += FIBERS->fibers[n].elts[3*i+1];
            zc += FIBERS->fibers[n].elts[3*i+2];
            if(minx>FIBERS->fibers[n].elts[3*i]) minx = FIBERS->fibers[n].elts[3*i];
            if(maxx<FIBERS->fibers[n].elts[3*i]) maxx = FIBERS->fibers[n].elts[3*i];

            if(miny>FIBERS->fibers[n].elts[3*i+1]) miny = FIBERS->fibers[n].elts[3*i+1];
            if(maxy<FIBERS->fibers[n].elts[3*i+1]) maxy = FIBERS->fibers[n].elts[3*i+1];

            if(minz>FIBERS->fibers[n].elts[3*i+2]) minz = FIBERS->fibers[n].elts[3*i+2];
            if(maxz<FIBERS->fibers[n].elts[3*i+2]) maxz = FIBERS->fibers[n].elts[3*i+2];
        }
    }
    xc = xc/((double) nbPoints);
    yc = yc/((double) nbPoints);
    zc = zc/((double) nbPoints);
    double sxc = 0.0;
    double syc = 0.0;
    double szc = 0.0;
    for(unsigned long n=0 ; n<FIBERS->N ; n++ )
    {
        for(unsigned long i=0 ; i<FIBERS->fibers[n].N_element ; i++ )
        {
            sxc += SQR(FIBERS->fibers[n].elts[3*i]-xc);
            syc += SQR(FIBERS->fibers[n].elts[3*i+1]-yc);
            szc += SQR(FIBERS->fibers[n].elts[3*i+2]-zc);
        }
    }
    sxc = sxc/((double) (nbPoints-1));
    syc = syc/((double) (nbPoints-1));
    szc = szc/((double) (nbPoints-1));
    cout << "center of fibers is: (" << xc << "," << yc << "," << zc << ")" << endl;
    cout << "std of fibers is: (" << sqrt(sxc) << "," << sqrt(syc) << "," << sqrt(szc) << ")" << endl;
    cout << "min (x,y,z) (" << minx << ", " << miny << ", " << minz << ")" << endl;
    cout << "max (x,y,z) (" << maxx << ", " << maxy << ", " << maxz << ")" << endl;
    return 9;
#endif


#ifdef MEAN_FIBER_ERROR
    cout << "Command must be run as: " << argv[0] << " tensor.hdr<string> mask.hdr<string> tx<double> ty<double> tz<double> file1.fibSRC<string>... fileN.fibSRC" << endl;
    cout << "only for fibers centered on the center defined by the mask" << endl;
    if(argc>6)
    {
        //C_ReadAnalyze75.h
        //ostringstream oss;
        FILE* f = fopen("measureErrorFiber.data","a");
        if(f==NULL)
        {
            cout << "Couldn't open data file." << endl;
            return -1;
        }

        C_load_fiber loadObj;
        C_measure M;

        ///don't forget to add the tensor data and the mask
        M.m_rawTensor = NULL;///two possibilities: load it from an analyze file or create it
        M.m_rawMask = NULL;///id as the previous line
        C_ReadAnalyze75 toolANA;
        M.m_rawTensor = toolANA.getRawData(argv[1]);
        M.m_rawMask = toolANA.getRawData(argv[2]);

        if (M.m_rawTensor==NULL)
        {
            cout << "tensor field is not loaded" << endl;
            return -1;
        }
        if (M.m_rawMask==NULL)
        {
            cout << "mask is not loaded" << endl;
            return -1;
        }

        C_parametric_comparison tool;
        tool.alpha = Pi/8.0;
        double L, S, epsilon = 0;


        string alpha = "";
        string beta = "";
        string sigma = "";
        string SNR = "";
        string probaMove = "";
        string stdMove = "";
        string PHI0 = "";
        string PHI1 = "";

        //THIS BOOLEAN DEFINES WHETHER OR NOT THE DATA ARE LOADED RAW (TRUE) OR APPROXIMATED (FALSE)
        bool rawDataCond = true;//false;

        for(int i=6 ; i<argc ; i++)
        {
            try
            {
                cout << argv[i] << endl;
                if(rawDataCond)
                {
                    cout << "raw fibers" << endl;
                }
                else
                {
                    cout << "approximated fibers" << endl;
                }

                CtlCtlStruct* FIBERS = loadObj.readFiber(argv[i],rawDataCond); ///load interpolated fibers

                //center the fibers
                //double xc = 0.0;
                //double yc = 0.0;
                //unsigned long nbPoints = 0;
                for(unsigned long n=0 ; n<FIBERS->N ; n++ )
                {
                    for(unsigned long i=0 ; i<FIBERS->fibers[n].N_element ; i++ )
                    {
                        FIBERS->fibers[n].elts[3*i] += atof(argv[3]);
                        FIBERS->fibers[n].elts[(3*i)+1] += atof(argv[4]);
                        FIBERS->fibers[n].elts[(3*i)+2] += atof(argv[5]);
                    }
                }
                //cout << "center of fibers is: (" << xc/((double) nbPoints) << "," << yc/((double) nbPoints) << ")" << endl;
                if(FIBERS==NULL)
                {
                    cout << "no fiber found" << endl;
                    epsilon = -1;
                    L = 0;
                    S = 0;
                }
                else
                {
                    try
                    {
                        //must be raw fiber
                        epsilon = M.error_fiber(FIBERS);//(argv[i]);
                        M.getFibersStat(argv[i], &L, &S);
                        cout << "err: " << epsilon << endl;
                        cout << "mean length " << L << " sdt length " << S << endl;
                    }
                    catch(...)
                    {
                        epsilon = -1;
                        L = 0;
                        S = 0;
                    }
                }

                ///if possible, get alpha, beta and sigma, otherwise set it to 0
                alpha = getAlpha(argv[i]);
                beta = getBeta(argv[i]);
                sigma = getSigma(argv[i]);
                SNR = getSNR(argv[i]);
                probaMove = getProbaMove(argv[i]);
                stdMove = getStdMove(argv[i]);
                PHI0 = getPhi0(argv[i]);
                PHI1 = getPhi1(argv[i]);
                if(SNR.empty())
                {
                    SNR = "100";
                }

                ///write parameters and measure: SNR, alpha, beta, sigma, mean length, std length, carole error, mean distance, mean angle, my measure
                //fprintf(f, "%s, %s, %s, %s, %s, %s, %s, %s, %f, %f, %f, %f, %f, %f\n", SNR.data(), alpha.data(), beta.data(), sigma.data(), PHI0.data(), PHI1.data(), probaMove.data(), stdMove.data(), L, S, epsilon, m, A, m2);
                fprintf(f, "%s, %s, %s, %s, %s, %s, %s, %s, %f, %f, %f\n", SNR.data(), alpha.data(), beta.data(), sigma.data(), PHI0.data(), PHI1.data(), probaMove.data(), stdMove.data(), L, S, epsilon);

                if(FIBERS!=NULL)
                {
                    delete FIBERS;
                }

            }
            catch(...)
            {
                cout << argv[i] << " failed computing errors" << endl;
            }

        }
        fclose(f);
    }
#endif
#ifdef FIBER_SIMILARITY
    cout << "Command must be run as: " << argv[0] << " tx<double> ty<double> tz<double> file1.fibSRC<string>... fileN.fibSRC" << endl;
    cout << "only for fibers centered on (0,0,z) with helix angle pi/8" << endl;
    if(argc>4)
    {
        //C_ReadAnalyze75.h
        ostringstream oss;
        FILE* f = fopen("measureSimilarityRaw.data","a");
        if(f==NULL)
        {
            cout << "Couldn't open data file." << endl;
            return -1;
        }

        C_load_fiber* loadObj = new C_load_fiber();
        C_measure* M = new C_measure();

        C_parametric_comparison* tool = new C_parametric_comparison();
        tool->alpha = Pi/8.0;
        double m, m2, A;


        string alpha = "";
        string beta = "";
        string sigma = "";
        string SNR = "";
        string probaMove = "";
        string stdMove = "";
        string PHI0 = "";
        string PHI1 = "";

        //THIS BOOLEAN DEFINES WHETHER OR NOT THE DATA ARE LOADED RAW (TRUE) OR APPROXIMATED (FALSE)
        bool rawDataCond = true;
        for(int i=4 ; i<argc ; i++)
        {
            try
            {
                cout << argv[i] << endl;
                if(rawDataCond)
                {
                    cout << "raw fibers" << endl;
                }
                else
                {
                    cout << "approximated fibers" << endl;
                }

                CtlCtlStruct* FIBERS = loadObj->readFiber(argv[i],rawDataCond);
                for(unsigned long n=0 ; n<FIBERS->N ; n++ )
                {
                    for(unsigned long i=0 ; i<FIBERS->fibers[n].N_element ; i++ )
                    {
                        FIBERS->fibers[n].elts[3*i] += atof(argv[1]);
                        FIBERS->fibers[n].elts[3*i+1] += atof(argv[2]);
                        FIBERS->fibers[n].elts[3*i+2] += atof(argv[3]);
                    }
                }
                //cout << "center of fibers is: (" << xc/((double) nbPoints) << "," << yc/((double) nbPoints) << ")" << endl;
                if(FIBERS==NULL)
                {
                    cout << "no fiber found" << endl;
                    m = -1;
                    m2 = -1;
                    A = -1.0;
                }
                else
                {

                    try
                    {
                        m = tool->getMeanDist(FIBERS);    ///mean distance
                        A = tool->getMeanDist4(FIBERS);   ///mean angle
                        m2 = tool->getMeanDist3(FIBERS);  ///actual measure
                        cout << "error m = " << m << " error m2 = " << m2 << endl;
                    }
                    catch(...)
                    {
                        m = -1;
                        m2 = -1;
                        A = -1.0;
                        cout << "an error has been caught while calculating mean fiber error m or m2" << endl;
                    }
                }

                ///if possible, get alpha, beta and sigma, otherwise set it to 0
                alpha = getAlpha(argv[i]);
                beta = getBeta(argv[i]);
                sigma = getSigma(argv[i]);
                SNR = getSNR(argv[i]);
                probaMove = getProbaMove(argv[i]);
                stdMove = getStdMove(argv[i]);
                PHI0 = getPhi0(argv[i]);
                PHI1 = getPhi1(argv[i]);
                if(SNR.empty())
                {
                    SNR = "100";
                }

                ///write parameters and measure: SNR, alpha, beta, sigma, mean length, std length, carole error, mean distance, mean angle, my measure
                //fprintf(f, "%s, %s, %s, %s, %s, %s, %s, %s, %f, %f, %f, %f, %f, %f\n", SNR.data(), alpha.data(), beta.data(), sigma.data(), PHI0.data(), PHI1.data(), probaMove.data(), stdMove.data(), L, S, epsilon, m, A, m2);
                fprintf(f, "%s, %s, %s, %s, %s, %s, %s, %s, %f, %f, %f\n", SNR.data(), alpha.data(), beta.data(), sigma.data(), PHI0.data(), PHI1.data(), probaMove.data(), stdMove.data(), m, A, m2);

                if(FIBERS!=NULL)
                {
                    delete FIBERS;
                }

            }
            catch(...)
            {
                cout << argv[i] << " failed computing errors" << endl;
            }

        }

        delete tool;
        delete M;
        delete loadObj;
        fclose(f);
    }
    else
    {
        cout << "graph file, fiber files" << endl;
    }
#endif

#ifdef CURVATURE

    if(argc>1)
    {
        C_load_fiber loadObj;

        //THIS BOOLEAN DEFINES WHETHER OR NOT THE DATA ARE LOADED RAW (TRUE) OR APPROXIMATED (FALSE)
        bool rawDataCond = true;//may try with false
        for(int i=1 ; i<argc ; i++)
        {
            try
            {
                //open file
                ostringstream oss;
                oss << "measureCurvatureDistributionMeanLength_" << i << ".data";
                //oss << "measureU2DistributionMeanLength_" << i << ".data";
                FILE* f = fopen(oss.str().data(),"w"); //"a"

                //load fibers
                cout << argv[i] << endl;
                if(rawDataCond)
                {
                    cout << "raw fibers" << endl;
                }
                else
                {
                    cout << "approximated fibers" << endl;
                }
                CtlCtlStruct* FIBERS = loadObj.readFiber(argv[i],rawDataCond);

                //calculate mean curvature for each fiber
                double cur;
                double dFx, dFy, dFz;
                double ddFx, ddFy, ddFz;
                double tmpX, tmpY, tmpZ;
                double L;
                for(unsigned long n=0 ; n<FIBERS->N ; n++ )
                {

                    if(FIBERS->fibers[n].N_element>3)
                    {
                        cur = 0.0;
                        if(false)
                        {
                            for(unsigned long i=1 ; i<FIBERS->fibers[n].N_element-1 ; i++ )
                            {
                                dFx = FIBERS->fibers[n].elts[3*(i+1)]-FIBERS->fibers[n].elts[3*(i-1)];
                                dFy = FIBERS->fibers[n].elts[3*(i+1)+1]-FIBERS->fibers[n].elts[3*(i-1)+1];
                                dFz = FIBERS->fibers[n].elts[3*(i+1)+2]-FIBERS->fibers[n].elts[3*(i-1)+2];
                                ddFx = FIBERS->fibers[n].elts[3*(i+1)]-2*FIBERS->fibers[n].elts[3*i]+FIBERS->fibers[n].elts[3*(i-1)];
                                ddFy = FIBERS->fibers[n].elts[3*(i+1)+1]-2*FIBERS->fibers[n].elts[3*i+1]+FIBERS->fibers[n].elts[3*(i-1)+1];
                                ddFz = FIBERS->fibers[n].elts[3*(i+1)+2]-2*FIBERS->fibers[n].elts[3*i+2]+FIBERS->fibers[n].elts[3*(i-1)+2];
                                tmpX = dFy*ddFz - dFz*ddFy;
                                tmpY = dFz*ddFx - dFx*ddFz;
                                tmpZ = dFx*ddFy - dFy*ddFx;
                                cur += sqrt(tmpX*tmpX + tmpY*tmpY + tmpZ*tmpZ)/pow(dFx*dFx + dFy*dFy + dFz*dFz,1.5);
                            }
                            cur /= ((double) FIBERS->fibers[n].N_element);
                        }
                        else if(false)
                        {
                            for(unsigned long i=1 ; i<FIBERS->fibers[n].N_element-1 ; i++ )
                            {
                                //calculate vector's coordinates
                                dFx = FIBERS->fibers[n].elts[3*(i+1)]-FIBERS->fibers[n].elts[3*i];
                                dFy = FIBERS->fibers[n].elts[3*(i+1)+1]-FIBERS->fibers[n].elts[3*i+1];
                                dFz = FIBERS->fibers[n].elts[3*(i+1)+2]-FIBERS->fibers[n].elts[3*i+2];

                                ddFx = FIBERS->fibers[n].elts[3*(i-1)] - FIBERS->fibers[n].elts[3*i];
                                ddFy = FIBERS->fibers[n].elts[3*(i-1)+1] - FIBERS->fibers[n].elts[3*i+1];
                                ddFz = FIBERS->fibers[n].elts[3*(i-1)+2] - FIBERS->fibers[n].elts[3*i+2];

                                //calculate cos of angle
                                tmpX = (dFx*ddFx+dFy*ddFy+dFz*ddFz)/(sqrt(dFx*dFx+dFy*dFy+dFz*dFz)*sqrt(ddFx*ddFx+ddFy*ddFy+ddFz*ddFz));


                                //get angle
                                if (tmpX < -1.0) tmpX = -1.0 ;
                                else if (tmpX > 1.0) tmpX = 1.0 ;
                                tmpY = acos(tmpX);

                                //compute curvature terme
                                if(tmpY==tmpY) cur +=  1.0/(1.0+exp( (tmpY-(Pi/2.0))/0.375 ));//checking for nan
                            }
                        }
                        else if(false)
                        {
                            tmpZ = 0.0;
                            for(unsigned long i=1 ; i<FIBERS->fibers[n].N_element-1 ; i++ )
                            {
                                //calculate vector's coordinates
                                dFx = FIBERS->fibers[n].elts[3*(i+1)]-FIBERS->fibers[n].elts[3*i];
                                dFy = FIBERS->fibers[n].elts[3*(i+1)+1]-FIBERS->fibers[n].elts[3*i+1];
                                dFz = FIBERS->fibers[n].elts[3*(i+1)+2]-FIBERS->fibers[n].elts[3*i+2];

                                ddFx = FIBERS->fibers[n].elts[3*(i-1)] - FIBERS->fibers[n].elts[3*i];
                                ddFy = FIBERS->fibers[n].elts[3*(i-1)+1] - FIBERS->fibers[n].elts[3*i+1];
                                ddFz = FIBERS->fibers[n].elts[3*(i-1)+2] - FIBERS->fibers[n].elts[3*i+2];

                                //calculate cos of angle
                                tmpX = (dFx*ddFx+dFy*ddFy+dFz*ddFz)/(sqrt(dFx*dFx+dFy*dFy+dFz*dFz)*sqrt(ddFx*ddFx+ddFy*ddFy+ddFz*ddFz));


                                //get angle
                                if (tmpX < -1.0) tmpX = -1.0 ;
                                else if (tmpX > 1.0) tmpX = 1.0 ;
                                tmpY = acos(tmpX);

                                //compute curvature terme
                                if(tmpY==tmpY) //checking for nan
                                {
                                    cur +=  (1.0/(1.0+exp( (tmpY-(Pi/2.0))/0.375 )))*sqrt(dFx*dFx+dFy*dFy+dFz*dFz);
                                    tmpZ += sqrt(dFx*dFx+dFy*dFy+dFz*dFz);
                                }
                            }
                        }
                        else
                        {
                            //this version makes sens: curve integral of the curvature normalised by the length.
                            L = 0.0;
                            for(unsigned long i=1 ; i<FIBERS->fibers[n].N_element-1 ; i++ )
                            {
                                dFx = FIBERS->fibers[n].elts[3*(i+1)]-FIBERS->fibers[n].elts[3*(i-1)];
                                dFy = FIBERS->fibers[n].elts[3*(i+1)+1]-FIBERS->fibers[n].elts[3*(i-1)+1];
                                dFz = FIBERS->fibers[n].elts[3*(i+1)+2]-FIBERS->fibers[n].elts[3*(i-1)+2];
                                ddFx = FIBERS->fibers[n].elts[3*(i+1)]-2*FIBERS->fibers[n].elts[3*i]+FIBERS->fibers[n].elts[3*(i-1)];
                                ddFy = FIBERS->fibers[n].elts[3*(i+1)+1]-2*FIBERS->fibers[n].elts[3*i+1]+FIBERS->fibers[n].elts[3*(i-1)+1];
                                ddFz = FIBERS->fibers[n].elts[3*(i+1)+2]-2*FIBERS->fibers[n].elts[3*i+2]+FIBERS->fibers[n].elts[3*(i-1)+2];
                                tmpX = dFy*ddFz - dFz*ddFy;
                                tmpY = dFz*ddFx - dFx*ddFz;
                                tmpZ = dFx*ddFy - dFy*ddFx;
                                //compute curvature terme
                                if(tmpX==tmpX && tmpY==tmpY && tmpZ==tmpZ) //checking for nan
                                {
                                    cur += sqrt(tmpX*tmpX + tmpY*tmpY + tmpZ*tmpZ)/pow(dFx*dFx + dFy*dFy + dFz*dFz,1.5);
                                    L += sqrt(dFx*dFx+dFy*dFy+dFz*dFz);
                                }

                            }
                        }

                        if(cur<0.0) cout << "negative curvature?" << endl;
                        if(L<0.0) cout << "negative length?" << endl;

                        if(!isnan(cur/L)) fprintf(f, "%f\n", cur/L);
                    }
                }

                fclose(f);

                if(FIBERS!=NULL)
                {
                    delete FIBERS;
                }

            }
            catch(...)
            {
                cout << argv[i] << " failed computing errors" << endl;
            }

        }
    }
    else
    {
        cout << "Command must be run as: " << argv[0] << " file1.fibSRC<string>... fileN.fibSRC" << endl;
    }
#endif

    return argc;
}


string getAlpha(string filename)
{
    size_t n = filename.find("alpha");
    if(n==filename.npos)
    {
        return (char*) "0";
    }
    n = filename.find("_", n);
    size_t m = filename.find("_", n+1);
    if(m==n || m==filename.npos)
    {
        return (char*) "0";
    }
    return filename.substr(n+1, m-(n+1));
}
string getBeta(string filename)
{
    size_t n = filename.find("beta");
    if(n==filename.npos)
    {
        return (char*) "0";
    }
    n = filename.find("_", n);
    size_t m = filename.find("_", n+1);
    if(m==n || m==filename.npos)
    {
        return (char*) "0";
    }
    return filename.substr(n+1, m-(n+1));
}
string getSigma(string filename)
{
    size_t n = filename.find("sigma");
    if(n==filename.npos)
    {
        return (char*) "0";
    }
    n = filename.find("_", n);
    size_t m = filename.find("_", n+1);
    if(m==n || m==filename.npos)
    {
        return (char*) "0";
    }
    return filename.substr(n+1, m-(n+1));
}

string getAmplitude(string filename)
{
    size_t n = filename.find("_A_");
    if(n==filename.npos)
    {
        return (char*) "0";
    }
    n = filename.find("_", n+1);
    size_t m = filename.find("_", n+1);
    if(m==n || m==filename.npos)
    {
        return (char*) "0";
    }
    return filename.substr(n+1, m-(n+1));
}

string getSNR(string filename)
{
    size_t n = filename.find("SNR");
    if(n==filename.npos)
    {
        return (char*) "INF";
    }
    n = filename.find("_", n);
    size_t m = filename.find("_", n+1);
    if(m==n || m==filename.npos)
    {
        return (char*) "INF";
    }
    return filename.substr(n+1, m-(n+1));
//    size_t n = filename.find("_M_");
//    if(n==filename.npos)
//    {
//        return (char*) "0";
//    }
//    n = filename.find("_", n+1);
//    size_t m = filename.find(".", n+1);
//    if(m==n || m==filename.npos)
//    {
//        return (char*) "0";
//    }
//    return filename.substr(n+1, m-(n+1));
}

string getProbaMove(string filename)
{
    size_t n = filename.find("p_");
    if(n==filename.npos)
    {
        return (char*) "0";
    }
    n = filename.find("_", n);
    size_t m = filename.find("_", n+1);
    if(m==n || m==filename.npos)
    {
        return (char*) "0";
    }
    return filename.substr(n+1, m-(n+1));
}

string getStdMove(string filename)
{
    size_t n = filename.find("stdmove");
    if(n==filename.npos)
    {
        return (char*) "0";
    }
    n = filename.find("_", n);
    size_t m = filename.find("_", n+1);
    if(m==n || m==filename.npos)
    {
        return (char*) "0";
    }
    return filename.substr(n+1, m-(n+1));
}

string getPhi0(string filename)
{
    size_t n = filename.find("phi0");
    if(n==filename.npos)
    {
        return (char*) "0.5";
    }
    n = filename.find("_", n);
    size_t m = filename.find("_", n+1);
    if(m==n || m==filename.npos)
    {
        return (char*) "0.5";
    }
    return filename.substr(n+1, m-(n+1));
}

string getPhi1(string filename)
{
    size_t n = filename.find("phi1");
    if(n==filename.npos)
    {
        return (char*) "0.5";
    }
    n = filename.find("_", n);
    size_t m = filename.find("_", n+1);
    if(m==n || m==filename.npos)
    {
        return (char*) "0.5";
    }
    return filename.substr(n+1, m-(n+1));
}
