#ifndef C_CROP_H
#define C_CROP_H

#ifdef CROP
#include <C_measure.h>


class C_crop : public C_measure
{
    public:
        /** Default constructor */
        C_crop(floating Rmin, floating Rmax, floating Zmin, floating Zmax);
        /** Default destructor */
        virtual ~C_crop();
        void cropHelixFiber(string fileNameSavedEdges);
    protected:
        SaveEdges* readSavedEdges(std::string fileNameEdge, unsigned long* nb_edge);
    private:
        floating m_Rmin;
        floating m_Rmax;
        floating m_Zmin;
        floating m_Zmax;
};

#endif

#endif // C_CROP_H
