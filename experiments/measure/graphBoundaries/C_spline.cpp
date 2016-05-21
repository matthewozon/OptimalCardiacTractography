#include <C_spline.h>

C_spline::C_spline()
{
    //ctor
}

C_spline::~C_spline()
{
    //dtor
}



//Bernstein
bool C_spline::computeFiberSpline(CtlCtlStruct* m_ctlCtlStruct, CtlCtlStruct** m_ctlCtlStructSpline, unsigned long nb_sub_section, unsigned short MODE)
{
    *m_ctlCtlStructSpline = new CtlCtlStruct(m_ctlCtlStruct->N);

    for(unsigned long i=0 ; i<(*m_ctlCtlStructSpline)->N ; i++)
    {
        ((*m_ctlCtlStructSpline)->fibers[i]).allocateArray(((m_ctlCtlStruct->fibers[i]).N_element-1)*nb_sub_section+1);

        for(unsigned long j=0 ; j<((*m_ctlCtlStructSpline)->fibers[i]).N_element ; j++ )
        {
            ((*m_ctlCtlStructSpline)->fibers[i]).idx[j] = (unsigned long) j/nb_sub_section;
        }

        //compute B-spline with the given control points and store it in dest buff
        computeFiberSpline(&(m_ctlCtlStruct->fibers[i]), &((*m_ctlCtlStructSpline)->fibers[i]), MODE);

    }
    return true;

}

bool C_spline::computeFiberSpline(CtlStruct* srcBuff, CtlStruct* destBuff, unsigned short MODE)
{
    long double du = 1.0/((double) (destBuff->N_element-1));

    long double* BernsteinCoef;
    BernsteinCoef = new long double[srcBuff->N_element-MODE];
    //unsigned long i=0;

    //for half
    //for(float u=0 ; u<=1 ; u+=du)
    long double u=0.0;
    for(unsigned long i=0 ; i<destBuff->N_element ; i++)
    {
        //compute bernstein coef for a given u and return them in BernsteinCoef
        computeBernsteinCoef(u, srcBuff->N_element-MODE, BernsteinCoef);

        //compute vector value for u and 1-u
        long double rx = 0;
        long double ry = 0;
        long double rz = 0;
        for(unsigned long j=0 ; j<srcBuff->N_element-MODE ; j++)
        {
            if(MODE==CURVE)
            {
                rx += BernsteinCoef[j]*srcBuff->elts[3*j];
                ry += BernsteinCoef[j]*srcBuff->elts[3*j+1];
                rz += BernsteinCoef[j]*srcBuff->elts[3*j+2];
            }
            else if(MODE==SPEED)
            {
                rx += (srcBuff->N_element-1)*BernsteinCoef[j]*(srcBuff->elts[3*(j+1)]-srcBuff->elts[3*j]);
                ry += (srcBuff->N_element-1)*BernsteinCoef[j]*(srcBuff->elts[3*(j+1)+1]-srcBuff->elts[3*j+1]);
                rz += (srcBuff->N_element-1)*BernsteinCoef[j]*(srcBuff->elts[3*(j+1)+2]-srcBuff->elts[3*j+2]);
            }
            else  ///ACCELERATION
            {
                rx += (srcBuff->N_element-1)*(srcBuff->N_element-2)*BernsteinCoef[j]*(srcBuff->elts[3*(j+2)]-2*(srcBuff->elts[3*(j+1)])+srcBuff->elts[3*j]);
                ry += (srcBuff->N_element-1)*(srcBuff->N_element-2)*BernsteinCoef[j]*(srcBuff->elts[3*(j+2)+1]-2*(srcBuff->elts[3*(j+1)+1])+srcBuff->elts[3*j+1]);
                rz += (srcBuff->N_element-1)*(srcBuff->N_element-2)*BernsteinCoef[j]*(srcBuff->elts[3*(j+2)+2]-2*(srcBuff->elts[3*(j+1)+2])+srcBuff->elts[3*j+2]);
            }
        }

        //fill the destination at two
        destBuff->elts[3*i] = rx;
        destBuff->elts[3*i+1] = ry;
        destBuff->elts[3*i+2] = rz;
        destBuff->idx[i] = (unsigned long) (srcBuff->N_element-1)*u;

        //update i
        //i++;

        //cout << "u=" << u << " rx=" << rx << " ry=" << ry << " rz=" << rz << endl;
        u+=du;
    }

    delete [] BernsteinCoef;
    return true;
}

bool C_spline::computeBernsteinCoef(long double u, unsigned long M, long double* BernsteinCoef)
{
    //init
    long double* tempB = new long double[M];
    long double Bprevious, Bcurrent;

    for(unsigned long i=0 ; i<M ; i++)
    {
        if(i==0)
        {
            //
            tempB[M-1-i] = 1.0;
        }
        else
        {
            tempB[M-1-i] = tempB[M-i]*(1.0-u);
        }
    }
    BernsteinCoef[0] = tempB[0];

    //loop
    unsigned long idxCurrent;
    for(unsigned long i=1 ; i<M ; i++)
    {
        Bprevious = 0;
        for(unsigned long k=i ; k<M ; k++)
        {
            //compute the current value
            idxCurrent = M-1-k;
            if(i==k)
            {
                //compute the first value
                Bcurrent = u*tempB[idxCurrent+1];
            }
            else
            {
                //compute a normal value
                Bcurrent = (1-u)*tempB[idxCurrent+1] + u*Bprevious;
            }
            //temporary store the previous value
            Bprevious = tempB[idxCurrent] ;

            //store the current value
            tempB[idxCurrent]  = Bcurrent;

        }

        //store it
        BernsteinCoef[i] = tempB[0];
    }

    return true;
}
