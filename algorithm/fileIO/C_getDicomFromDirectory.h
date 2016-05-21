#ifndef C_GETDICOMFROMDIRECTORY_H
#define C_GETDICOMFROMDIRECTORY_H

#include <C_sort_files.h>
#include <C_ReadDicomDiffInfo.h>
#include <C_readDicomDiffData.h>

#include <stdio.h>

using namespace std;

class C_getDicomFromDirectory
{
    public:
        /** Default constructor */
        C_getDicomFromDirectory(string directoryName);
        /** Default destructor */
        virtual ~C_getDicomFromDirectory();


        vector<diffInfo*> getDiffInfoS0(void);
        vector<rawData<double>*> getDataS0(void);
        vector<vector<diffInfo*> > getDiffInfoSi(void);
        vector<vector<rawData<double>*> > getDataSi(void);
        rawData<double>* getDataROI(void);
        unsigned long getNumDirection(void){return NgradientDirection;};
        unsigned long getNumRepetition(void){return Nrepetition;};

    protected:
        string m_directoryName;
        vector< /**gradient direction*/ vector< /**averaging*/string > > dataFile;
        vector< string > SoFile;
        string roiFileName;
        unsigned long NgradientDirection;
        unsigned long Nrepetition;
        unsigned long m_nb_slice_in_mosa;
};

#endif // C_GETDICOMFROMDIRECTORY_H
