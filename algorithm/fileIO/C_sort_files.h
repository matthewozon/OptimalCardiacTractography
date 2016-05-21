#ifndef C_SORT_FILES_H
#define C_SORT_FILES_H

#include <vector>
#include <string>
#include <math.h>
#include <dirent.h> //get file in directory
#include <algorithm> //sorting algorithm

#include <C_ReadDicomDiffInfo.h>

#include <iostream>
#include <stdio.h>

using namespace std;


class C_sort_files
{
    public:
        /** Default constructor */
        C_sort_files();
        /** Default destructor */
        virtual ~C_sort_files();

        ///DCM
        vector<string> getDicomFileName(string dirName);
        vector<string> getS0FileNameDCM(vector<string> fileNames, unsigned long NgradientDirection, unsigned long Nrepetition);
        vector<vector<string> > getSiFileNameDCM(vector<string> fileNames, unsigned long NgradientDirection, unsigned long Nrepetition);
        vector<string> getS0FileNameDCM(string directoryName);
        vector<vector<string> > getSiFileNameDCM(string directoryName);
        string getROIFileNameDCM(string directoryName);

        ///ANA
        vector<string> getANAFileName(string dirName);
        vector<string> getS0FileNameANA(vector<string> fileNames, unsigned long NgradientDirection, unsigned long Nrepetition, string patternInFileNames);
        vector<vector<string> > getSiFileNameANA(vector<string> fileNames, unsigned long NgradientDirection, unsigned long Nrepetition, string patternInFileNames);

};


#endif // C_SORT_FILES_H
