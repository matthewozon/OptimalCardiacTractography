#ifndef C_SPLINE_H
#define C_SPLINE_H


/**this class need to be checked
    generate a set of control point circle shaped with 100 points.
    export in file .dat the curve, speed and accelaration
    DO NOT DELETE THIS CLASS: IT COMPUTES SPLINE CURVE, SPEED AND ACCELERATION
*/

//#include <GL/gl.h>
//#include <structs.h>
//#include <puffinDef.h>
#include <vector>
//#include <math.h>
#include <struct.h>

#define SQR(x) ((x)*(x))
#define CUB(x) ((x)*(x)*(x))

#include <iostream>
#include <stdio.h>
//struct CtlStruct{ //this struct is to send data back to python
//    unsigned long N_element; //means there are actualy 3*N_element GLfloat (because of x, y, z)
//	long double* elts;
//	unsigned long* idx;
//	CtlStruct()
//	{
//	    N_element = 0;
//	    elts = NULL;
//	    idx = NULL;
//	}
//	~CtlStruct()
//	{
//	    N_element = 0;
//	    deleteArray();
//	}
//	void allocateArray(unsigned long N)
//	{
//	    deleteArray();
//	    N_element = N;
//	    elts = new long double[3*N_element];
//	    idx = new unsigned long[N_element];
//	}
//
//	private:
//	void deleteArray(void)
//	{
//	    if(elts!=NULL)
//	    {
//	        delete elts;
//	        elts = NULL;
//	    }
//	    if(idx!=NULL)
//	    {
//	        delete idx;
//	        idx = NULL;
//	    }
//	}
//};

//struct CtlCtlStruct{ //this struct is to send data back to python
//    unsigned long N;
//	CtlStruct* fibers;
//	CtlCtlStruct(unsigned long n)
//	{
//	    N = n;
//	    fibers = NULL;
//	    allocateCtrlStruct(N);
//	}
//	~CtlCtlStruct()
//	{
//	    N = 0;
//	    deleteCtrlStruct();
//	}
//
//	private:
//	void allocateCtrlStruct(unsigned long n)
//	{
//	    if(fibers!=NULL)
//	    {
//	        deleteCtrlStruct();
//	    }
//	    fibers = new CtlStruct[n];
//	}
//	void deleteCtrlStruct(void)
//	{
//	    for(unsigned long i=0 ; i<N ; i++)
//	    {
//	        fibers[i].~CtlStruct();
//	    }
//	    delete[] fibers; //may be not?
//	    fibers = NULL;
//
//	}
//};

#define CURVE 0
#define SPEED 1
#define ACCELERATION 2

using namespace std;
class C_spline
{
    public:
        /** Default constructor */
        C_spline();
        /** Default destructor */
        virtual ~C_spline();

        bool computeFiberSpline(CtlCtlStruct* m_ctlCtlStruct, CtlCtlStruct** m_ctlCtlStructSpline, unsigned long nb_sub_section, unsigned short MODE=CURVE);
    protected:
    private:
        //B-spline
        bool computeFiberSpline(CtlStruct* srcBuff, CtlStruct* destBuff, unsigned short MODE);
        bool computeBernsteinCoef(long double u, unsigned long M, long double* BernsteinCoef);

};

#endif // C_SPLINE_H
