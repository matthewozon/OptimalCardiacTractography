#include <C_toolbox_spline.h>

C_toolbox_spline::C_toolbox_spline():C_toolbox()
{
    //ctor
}

C_toolbox_spline::~C_toolbox_spline()
{
    //dtor
}


//Bernstein
bool C_toolbox_spline::computeFiberSpline(CtlCtlStruct* m_ctlCtlStruct, CtlCtlStruct** m_ctlCtlStructSpline, unsigned long nb_sub_section)
{
    *m_ctlCtlStructSpline = new CtlCtlStruct(m_ctlCtlStruct->N);

    unsigned long NB_THREAD = 8;
    C_thread_spline** THREAD_SPLINE = new C_thread_spline*[NB_THREAD];
    for(unsigned long i=0 ; i<NB_THREAD ; i++)
    {
        THREAD_SPLINE[i] = new C_thread_spline(this);
        ///set attributes of thread
        THREAD_SPLINE[i]->idxStart = (i*((*m_ctlCtlStructSpline)->N))/NB_THREAD;
        THREAD_SPLINE[i]->idxEnd = (((i+1)*((*m_ctlCtlStructSpline)->N))/NB_THREAD);
        THREAD_SPLINE[i]->m_ctlCtlStruct = m_ctlCtlStruct;
        THREAD_SPLINE[i]->m_ctlCtlStructSpline = m_ctlCtlStructSpline;
        THREAD_SPLINE[i]->nb_sub_section = nb_sub_section;
        ///start thread
        THREAD_SPLINE[i]->start();
    }

    ///wait for all thread to finish job
    for(unsigned long i=0 ; i<NB_THREAD ; i++)
    {
        THREAD_SPLINE[i]->wait_for_exit();
    }

    //delete all threads
    for(unsigned long i=0 ; i<NB_THREAD ; i++)
    {
        if(THREAD_SPLINE[i]!=NULL)
        {
            delete THREAD_SPLINE[i];
        }
    }
    delete THREAD_SPLINE;
    return true;

}

bool C_toolbox_spline::computeFiberSpline(CtlStruct* srcBuff, CtlStruct* destBuff, unsigned long nb_sub_section)
{
    unsigned long m = srcBuff->N_element-1;
    unsigned long M = m+1;
    unsigned long N = nb_sub_section;
    float du = 1.0/((float) (N*(M-1)));


    //GLfloat* BernsteinCoef;
    float* BernsteinCoef = new float[M];
    unsigned long i=0;

    //for half
    for(float u=0 ; u<=1 ; u+=du)
    {
        //compute bernstein coef for a given u and return them in BernsteinCoef
        computeBernsteinCoef(u, M, BernsteinCoef);

        //compute vector value for u and 1-u
        float rx = 0;
        float ry = 0;
        float rz = 0;
        for(unsigned long j=0 ; j<M ; j++)
        {
            rx += BernsteinCoef[j]*srcBuff->elts[3*j];
            ry += BernsteinCoef[j]*srcBuff->elts[3*j+1];
            rz += BernsteinCoef[j]*srcBuff->elts[3*j+2];
        }

        //fill the destination at two
        destBuff->elts[3*i] = rx;
        destBuff->elts[3*i+1] = ry;
        destBuff->elts[3*i+2] = rz;
        destBuff->idx[i] = (unsigned long) (M-1)*du;

        //update i
        i++;


    }

    delete BernsteinCoef;
    return true;
}

//bool C_spline::computeBernsteinCoef(float u, unsigned long M, GLfloat* BernsteinCoef)
bool C_toolbox_spline::computeBernsteinCoef(float u, unsigned long M, float* BernsteinCoef)
{
    //init
    float* tempB = new float[M];
    float Bprevious, Bcurrent;

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
