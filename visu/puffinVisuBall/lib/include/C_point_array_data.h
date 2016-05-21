#ifndef C_POINT_ARRAY_DATA_H
#define C_POINT_ARRAY_DATA_H

#include <stdlib.h>


class C_point_array_data
{
    public:
        /** Default constructor */
        C_point_array_data(unsigned long N, unsigned long spaceDim/**=3*/, unsigned long vectorDim/**=3*/, bool yes);
        /** Default destructor */
        virtual ~C_point_array_data();

        ///allocator
        void allocateArrays(unsigned long N, unsigned long spaceDim=3, unsigned long vectorDim=3, bool yes=false);

        ///delete
        void deleteArrays(void);

        ///data
        unsigned long nb_point; //number of points (lenth of interpolation)
        unsigned long dim; //dimension of space always <=3 is anyway = 3
        unsigned long dimV; //vector dimension <=3
        double** pointCoordinate; // an array of size = nb_point * dim
        double** vectorInterpolated; // an array of size = nb_point * dimV
        double** vector3Interpolated; // an array of size = nb_point * dimV
        double* FA;
        double* eig1;
        double* eig2;
        double* eig3;
};

#endif // C_POINT_ARRAY_DATA_H
