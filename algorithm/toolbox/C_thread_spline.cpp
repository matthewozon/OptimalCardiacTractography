#include <C_thread_spline.h>

C_thread_spline::C_thread_spline(C_toolbox_spline* _Spline)
{
    //ctor
    m_Spline = _Spline;
}

C_thread_spline::~C_thread_spline()
{
    //dtor
}

void C_thread_spline::execute()
{

    for(unsigned long i=idxStart ; i<idxEnd ; i++)
    {
        ((*m_ctlCtlStructSpline)->fibers[i]).allocateArray(((m_ctlCtlStruct->fibers[i]).N_element-1)*nb_sub_section+1);

        for(unsigned long j=0 ; j<((*m_ctlCtlStructSpline)->fibers[i]).N_element ; j++ )
        {
            ((*m_ctlCtlStructSpline)->fibers[i]).idx[j] = (unsigned long) j/nb_sub_section;
        }

        //compute B-spline with the given control points and store it in dest buff
        m_Spline->computeFiberSpline(&(m_ctlCtlStruct->fibers[i]), &((*m_ctlCtlStructSpline)->fibers[i]), nb_sub_section);

    }
//    for(unsigned long i=idxStart ; i<idxEnd/**(*m_ctlCtlStructSpline)->N*/ ; i++)
//    {
//        //allocate a new fiber
//        ((*m_ctlCtlStructSpline)->fibers[i]).N_element = ((m_ctlCtlStruct->fibers[i]).N_element-1)*nb_sub_section+1;
//        //((*m_ctlCtlStructSpline)->fibers[i]).elts = new GLfloat[3*((*m_ctlCtlStructSpline)->fibers[i]).N_element];
//        ((*m_ctlCtlStructSpline)->fibers[i]).elts = new float[3*((*m_ctlCtlStructSpline)->fibers[i]).N_element];
//        ((*m_ctlCtlStructSpline)->fibers[i]).idx = new unsigned long[((*m_ctlCtlStructSpline)->fibers[i]).N_element];
//
//        for(unsigned long j=0 ; j<((*m_ctlCtlStructSpline)->fibers[i]).N_element ; j++ )
//        {
//            ((*m_ctlCtlStructSpline)->fibers[i]).idx[j] = (unsigned long) j/nb_sub_section;
//        }
//
//        //compute B-spline with the given control points and store it in dest buff
//        m_Spline->computeFiberSpline(&(m_ctlCtlStruct->fibers[i]), &((*m_ctlCtlStructSpline)->fibers[i]), nb_sub_section);
//
//    }
    return;
}
