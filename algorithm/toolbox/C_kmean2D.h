#ifndef C_KMEAN2D_H
#define C_KMEAN2D_H

#include <stdlib.h>
#include <stdio.h>
#include <cmath>
#include <vector>
#include <iostream>
#include <MACRO_DEF.h>

using namespace std;

class C_kmean2D ///should create the same in 3D
{
    public:
        /** Default constructor */
        C_kmean2D(double** mIM, unsigned long mL, unsigned long mC, bool seg=false); //do not allocate memory for segmented image by default
        /** Default destructor */
        virtual ~C_kmean2D();

        bool run(unsigned long N /*number of centroids to create*/);

        vector<double> mu; //array of mean value
        vector<double> sig; //array of std
        unsigned long nbCentroids;
    protected:
        vector<unsigned long> /*vector of pixel index closest to X[x]*/ pixelCloseTo(vector<double> X /*array of current centroids*/, unsigned long x /*index of the centroid doncidered*/);
        unsigned long /*label of pixelValue*/ pixelCloseTo(vector<double> X /*array of current centroids*/, double pixelValue);

        double** IM; //do not destroy this pointer (copy of image pointer)
        bool isIMSEG;
        unsigned long** IMSEG; //labels of segmented image refer to centroid indices
        unsigned long L; //nb rows
        unsigned long C; //nb columns
        double getMeanOFValuesAtIndex(vector<unsigned long> IDX); //mu = sum(IM(IDX))/(length(IDX))
        double getStandardDeviationOFValuesAtIndex(vector<unsigned long> IDX); //sig = sqrt( sum(SQR(IM(IDX)-mu))/(length(IDX)-1) )
        unsigned long getRowFromIdx(unsigned long idx);
        unsigned long getColumnFromIdx(unsigned long idx);
        unsigned long getIdxFromRowAndColumn(unsigned long l, unsigned long c);
    public:
        double getMaxValue(void);
        double getMinValue(void);
};

#endif // C_KMEAN_H
