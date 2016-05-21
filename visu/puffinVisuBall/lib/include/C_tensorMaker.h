#ifndef C_TENSORMAKER_H
#define C_TENSORMAKER_H

#include <C_toolbox_SVD.h>
//#include <C_readDicomDiffData.h>
#include <rawData.h>
#include <struct.h>
#include <MACRO_DEF.h>

class C_tensorMaker
{
    public:
        /** Default constructor */
        C_tensorMaker();
        /** Default destructor */
        virtual ~C_tensorMaker();
        /**OK*/rawData<double>* makeTensors(vector<vector<diffInfo*> > dwiInfo, vector<vector<rawData<double>*> > dwiData, vector<diffInfo*> b0Info, vector<rawData<double>*> b0Data);
        /**OK*/rawData<double>* makeTensors(vector<vector<diffInfo*> > dwiInfo, vector<vector<rawData<double>*> > dwiData, vector<diffInfo*> b0Info, vector<rawData<double>*> b0Data, rawData<double>* dwiMask);

        ///dwiData must be 4D with at least 6 elements in last dimension and dwiMask must be 3D. We assume that there is no average in DWI so that the la dimension has Ndir+1 elements
        /**OK*/rawData<double>* makeTensors(vector<diffInfo*> dwiInfo, rawData<double>* dwiDataS, rawData<double>* dwiMask);

        ///retrieve tensors that were already computed and stored in ANA file (last dimension must be 6 for the 6 entries of tensor xx (0), xy (1), xz (2), yy(3), yz(4) and zz(5))
        /**some cases are not handled*/rawData<double>* getTensors(rawData<double>* dwiDataS, rawData<double>* dwiMask);

        ///compute a single tensor
        /**to be checked*/double* makeTensor(double* dwiArray/***N rows*/, double** Hpsi /**6 rows by N columns*/, double b /**b-value for dwiArray*/, double S0, unsigned long N);
        /**to be checked*/double** makeTensor2(double* dwiArray/***N rows*/, double** Hpsi /**6 rows by N columns*/, double b /**b-value for dwiArray*/, double S0, unsigned long N);
        /**to be checked*/double** makeHpsi(vector<vector<diffInfo*> > dwiInfo, vector<diffInfo*> b0Info);
    protected:
        C_toolbox_SVD* toolSVD;
        //double** m_Hpsi; //pseudo inverse storage 6 rows by NgradDirection columns
        //double** Y;    //something like (1/b)ln(S0/Si): NgradDirection rows by 1 column
        /**OK*/vector<rawData<double>*> averageRawData(vector<vector<rawData<double>*> > myData);
        /**OK*/bool computeTensors(vector<rawData<double>*> myAverageDwiDataVector, rawData<double>* myAverageB0DatVector, rawData<double>* dwiMask, double b, double**H, rawData<double>* tensorOutput);
};

#endif // C_TENSORMAKER_H
