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


#define CONVEX_HULL
//#define BOUNDING_BOX

int main(int argc, char **argv)
{

    if(argc!=2)
    {
        cout << "Command must be run as: " << argv[0] << " file.fibSRC <string>" << endl;
    }
    else
    {
#ifdef CONVEX_HULL
        FILE* f = fopen("measureFiberBulkConvex.data","a");
#else
        FILE* f = fopen("measureFiberBulkBox.data","a");
#endif
        if(f==NULL)
        {
            cout << "Couldn't open data file." << endl;
            return -1;
        }

        bool rawDataCond = true; //only for raw data
        try
        {
            cout << argv[1] << endl;
            C_load_fiber loadObj;
            C_measure M;
#ifdef CONVEX_HULL
            vector<double> L, totPerimeter;
#else
            vector<double> L, dU,dV,dW;
#endif

            CtlCtlStruct* FIBERS = loadObj.readFiber(argv[1],rawDataCond);

            if(FIBERS!=NULL)
            {
#ifdef CONVEX_HULL
                //try a real convex hull algorithm
                //delete FIBERS;
//                    CtlCtlStruct* FIBERS = new CtlCtlStruct(1);
//                    FIBERS->fibers[0].N_element = 14;
//                    FIBERS->fibers[0].elts = new double[FIBERS->fibers[0].N_element*3];
//                    FIBERS->fibers[0].elts[3*0] = 0.0; FIBERS->fibers[0].elts[3*0+1] = 0.0; FIBERS->fibers[0].elts[3*0+2] = 0.0;
//                    FIBERS->fibers[0].elts[3*1] = 1.0; FIBERS->fibers[0].elts[3*1+1] = 0.0; FIBERS->fibers[0].elts[3*1+2] = 0.0;
//                    FIBERS->fibers[0].elts[3*2] = 0.0; FIBERS->fibers[0].elts[3*2+1] = 1.0; FIBERS->fibers[0].elts[3*2+2] = 0.0;
//                    FIBERS->fibers[0].elts[3*3] = 1.0; FIBERS->fibers[0].elts[3*3+1] = 1.0; FIBERS->fibers[0].elts[3*3+2] = 0.0;

//                    FIBERS->fibers[0].elts[3*4] = 0.0; FIBERS->fibers[0].elts[3*4+1] = 0.0; FIBERS->fibers[0].elts[3*4+2] = 1.0;
//                    FIBERS->fibers[0].elts[3*5] = 1.0; FIBERS->fibers[0].elts[3*5+1] = 0.0; FIBERS->fibers[0].elts[3*5+2] = 1.0;
//                    FIBERS->fibers[0].elts[3*6] = 0.0; FIBERS->fibers[0].elts[3*6+1] = 1.0; FIBERS->fibers[0].elts[3*6+2] = 1.0;
//                    FIBERS->fibers[0].elts[3*7] = 1.0; FIBERS->fibers[0].elts[3*7+1] = 1.0; FIBERS->fibers[0].elts[3*7+2] = 1.0;


//                    FIBERS->fibers[0].elts[3*8] = 0.5; FIBERS->fibers[0].elts[3*8+1] = 0.5; FIBERS->fibers[0].elts[3*8+2] = 0.5;
//                    FIBERS->fibers[0].elts[3*9] = 0.5; FIBERS->fibers[0].elts[3*9+1] = 0.5; FIBERS->fibers[0].elts[3*9+2] = 0.51;
//                    FIBERS->fibers[0].elts[3*10] = 0.6; FIBERS->fibers[0].elts[3*10+1] = 0.4; FIBERS->fibers[0].elts[3*10+2] = 0.2;
//                    FIBERS->fibers[0].elts[3*11] = 0.7; FIBERS->fibers[0].elts[3*11+1] = 0.5; FIBERS->fibers[0].elts[3*11+2] = 0.11;
//                    FIBERS->fibers[0].elts[3*12] = 0.6; FIBERS->fibers[0].elts[3*12+1] = 0.6; FIBERS->fibers[0].elts[3*12+2] = 0.01;
//                    FIBERS->fibers[0].elts[3*13] = 0.8; FIBERS->fibers[0].elts[3*13+1] = 0.7; FIBERS->fibers[0].elts[3*13+2] = 0.51;
                M.getLengthAndNormalizedTotalPerimeter(FIBERS,&L, &totPerimeter);
                for(unsigned int i=0 ; i<FIBERS->N ; i++)
                {
                    fprintf(f, "%f, %f\n", L.at(i), totPerimeter.at(i));
                }
#else
                M.getLengthAndBoxDimension(FIBERS,&L, &dU, &dV, &dW);
                for(unsigned int i=0 ; i<FIBERS->N ; i++)
                {
                    fprintf(f, "%f, %f, %f, %f\n", L.at(i), dU.at(i), dV.at(i), dW.at(i));
                }
#endif
                delete FIBERS;
            }

        }
        catch(...)
        {
            cout << argv[1] << " failed computing errors" << endl;
        }

        fclose(f);
    }
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
