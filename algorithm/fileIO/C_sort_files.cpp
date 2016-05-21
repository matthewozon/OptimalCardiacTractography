#include <C_sort_files.h>

C_sort_files::C_sort_files()
{
    //ctor
}

C_sort_files::~C_sort_files()
{
    //dtor
}

//"callback function
bool numericalSorting(string str1,string str2){
    unsigned found1 = str1.find_last_of(".");
    unsigned found2 = str2.find_last_of(".");
    unsigned short num1 = atoi(str1.substr(found1-4,4).data());
    unsigned short num2 = atoi(str2.substr(found2-4,4).data());
    return (num1<num2);
}
bool numericalSortingNoSuffix(string str1,string str2){
    unsigned found1 = str1.length();
    unsigned found2 = str2.length();
    unsigned short num1 = atoi(str1.substr(found1-4,4).data());
    unsigned short num2 = atoi(str2.substr(found2-4,4).data());
    return (num1<num2);
}
bool vectorEquality(vector<double> A, vector<double> B)
{
    if(A.size()!=B.size()) return false;
    double x = 0.0;
    for(unsigned long i=0 ; i<A.size() ; i++)
    {
        x += A.at(i)*B.at(i);
    }
    return (x*x>0.999*0.999);
}

vector<string> C_sort_files::getDicomFileName(string dirName)
{
    vector<string> fileNames;
    DIR *d;
    dirent *dir;
    d = opendir(dirName.data());
    if (d)
    {
        while ((dir = readdir(d)) != NULL)
        {
            string str = dir->d_name;
            if(str.length()>4)
            {
                unsigned found = str.find_last_of(".");
                string suffix = str.substr(found+1);
                string prefix = str.substr(0,found);
                if(prefix.length()>0 && suffix.length()>0)
                {
                    if(suffix.compare("dcm")==0)
                    {
                        if(prefix.substr(0,2).compare("IM")==0) //special selection
                        {
                            //unsigned short num = atoi(str.substr(found-4,4).data());
                            //if(prefix.compare("ROI")!=0 && prefix.compare("ROI2")!=0)
                            //{
                                fileNames.push_back(dir->d_name);
                            //}
                        }

                    }
                }
            }
        }
        closedir(d);
    }

    ///sort filenames
    sort(fileNames.begin(), fileNames.end(), numericalSorting);
    return fileNames;
}
vector<string> C_sort_files::getS0FileNameDCM(vector<string> fileNames, unsigned long NgradientDirection, unsigned long Nrepetition)
{
    vector< string > SoFile(Nrepetition,"");
    unsigned long j=0;
    for(unsigned int i=0 ; i<fileNames.size() ; i++)
    {
        if((i%(NgradientDirection+1))==0)
        {
            SoFile.at(j) = fileNames.at(i);
            j++;
        }
        if(j==Nrepetition)
        {
            i = fileNames.size();
        }
    }
    return SoFile;
}
vector<vector<string> > C_sort_files::getSiFileNameDCM(vector<string> fileNames, unsigned long NgradientDirection, unsigned long Nrepetition)
{
    vector< /**gradient direction*/ vector< /**averaging*/string > > dataFile(NgradientDirection, vector< /**averaging*/string >(Nrepetition,"") );
    unsigned long j=0;
    for(unsigned int i=0 ; i<fileNames.size() ; i++)
    {
        if((i%(NgradientDirection+1))!=0)
        {
            unsigned int nr = floor(((double)i)/((double)(NgradientDirection+1)));
            unsigned int ng = i-nr*(NgradientDirection+1)-1;
            dataFile.at(ng).at(nr) = fileNames.at(i);
            j++;
        }
        if(j==(NgradientDirection*Nrepetition))
        {
            i = fileNames.size();
        }
    }
    return dataFile;
}
vector<string> C_sort_files::getS0FileNameDCM(string directoryName)
{
    vector<string> fileNames;
    DIR *d;
    dirent *dir;
    d = opendir(directoryName.data());
    if (d)
    {
        while ((dir = readdir(d)) != NULL)
        {
            string str = dir->d_name;
            if(str.length()>4)
            {
                unsigned found = str.find_last_of(".");
                string suffix = str.substr(found+1);
                string prefix = str.substr(0,found);
                if(prefix.length()>0 && suffix.length()>0)
                {
                    if(suffix.compare("dcm")==0)
                    {
                        fileNames.push_back(dir->d_name);
                    }
                }
            }
        }
        closedir(d);
    }
    ///list and store
    vector< string > SoFile;
    C_ReadDicomDiffInfo* toolREAD_INFO_DCM;
    double b;
    for(unsigned int i=0 ; i<fileNames.size() ; i++)
    {
        toolREAD_INFO_DCM = new C_ReadDicomDiffInfo(fileNames.at(i));
        b = toolREAD_INFO_DCM->getBValue();
        if(b<20 && b>=0.0)
        {
            SoFile.push_back(fileNames.at(i));
        }
        delete toolREAD_INFO_DCM;
    }
    fileNames.clear();
    return SoFile;
}
vector<vector<string> > C_sort_files::getSiFileNameDCM(string directoryName)
{
    vector<string> fileNames;
    DIR *d;
    dirent *dir;
    d = opendir(directoryName.data());
    if (d)
    {
        while ((dir = readdir(d)) != NULL)
        {
            string str = dir->d_name;
            if(str.length()>4)
            {
                unsigned found = str.find_last_of(".");
                string suffix = str.substr(found+1);
                string prefix = str.substr(0,found);
                if(prefix.length()>0 && suffix.length()>0)
                {
                    if(suffix.compare("dcm")==0)
                    {
                        fileNames.push_back(dir->d_name);

                    }
                }
            }
        }
        closedir(d);
    }
    ///count and store file b!=0 (also store all directions)
    vector< string > SiFileArray;
//    vector< double > gradX;
//    vector< double > gradY;
//    vector< double > gradZ;
//    vector< unsigned long > idxGrad;
//    double* gradVect;
    C_ReadDicomDiffInfo* toolREAD_INFO_DCM;
    for(unsigned int i=0 ; i<fileNames.size() ; i++)
    {
        toolREAD_INFO_DCM = new C_ReadDicomDiffInfo(fileNames.at(i));
        if(toolREAD_INFO_DCM->getBValue()>=20)
        {
            SiFileArray.push_back(fileNames.at(i));
//            cout << fileNames.at(i) << endl;
//            gradVect = toolREAD_INFO_DCM->getGradientDirection();
//            gradX.push_back(gradVect[0]);
//            gradY.push_back(gradVect[1]);
//            gradZ.push_back(gradVect[2]);
//            delete gradVect;
        }
        delete toolREAD_INFO_DCM;
    }
    fileNames.clear();

    ///do not compute unique: 1) there are no pairs of vectors which are strictly the same (mechanical constraints), 2) svd will do the job
//    ///make a uniq on directions
//    vector< double > gradXUniq;
//    vector< double > gradYUniq;
//    vector< double > gradZUniq;
//    gradXUniq.push_back(gradX.at(0));
//    gradYUniq.push_back(gradY.at(0));
//    gradZUniq.push_back(gradZ.at(0));
//    idxGrad.push_back(0);
//    bool belongToUniqDirection;
//    for(unsigned int i=1 ; i<gradX.size() ; i++)
//    {
//        belongToUniqDirection = false;
//        for(unsigned int j=0 ; j<gradXUniq.size() ; j++)
//        {
//            if(gradX.at(i)==gradXUniq.at(j) && gradY.at(i)==gradYUniq.at(j) && gradZ.at(i)==gradZUniq.at(j))
//            //if( (gradX.at(i)*gradXUniq.at(j) + gradY.at(i)*gradYUniq.at(j) + gradZ.at(i)*gradZUniq.at(j))*(gradX.at(i)*gradXUniq.at(j) + gradY.at(i)*gradYUniq.at(j) + gradZ.at(i)*gradZUniq.at(j))>0.999*0.999 ) //assuming that gradient direction are normalized
//            {
//                //turn condition true
//                belongToUniqDirection = true;
//                //set the gradient idex
//                idxGrad.push_back(j);
//                //set j to final value
//                j = gradXUniq.size();
//            }
//        }
//        if(!belongToUniqDirection)
//        {
//            //update gradient list
//            gradXUniq.push_back(gradX.at(i));
//            gradYUniq.push_back(gradY.at(i));
//            gradZUniq.push_back(gradZ.at(i));
//            //set a new index (because it was not existing in the previous loop, no entry was created)
//            idxGrad.push_back(gradXUniq.size()-1);
//        }
//    }
    ///determine the number of repetition
    unsigned long Ndirection = SiFileArray.size();
    unsigned long Nrepetition = 1;//gradX.size()/Ndirection;
    vector< /**gradient direction*/ vector< /**averaging*/string > > SiFile(Ndirection, vector< /**averaging*/string >(0/*Nrepetition*/,"") );
    for(unsigned int i=0 ; i<Ndirection/*SiFileArray.size()*/ ; i++)
    {
        SiFile.at(i).push_back(SiFileArray.at(i));
    }
    SiFileArray.clear();
//    idxGrad.clear();
//    gradX.clear();
//    gradY.clear();
//    gradZ.clear();
//    gradXUniq.clear();
//    gradYUniq.clear();
//    gradZUniq.clear();
    return SiFile;
}

string C_sort_files::getROIFileNameDCM(string directoryName)
{
    string ROI_DCM = "";
    vector<string> fileNames;
    DIR *d;
    dirent *dir;
    d = opendir(directoryName.data());
    if (d)
    {
        while ((dir = readdir(d)) != NULL)
        {
            string str = dir->d_name;
            if(str.length()>4)
            {
                unsigned found = str.find_last_of(".");
                string suffix = str.substr(found+1);
                string prefix = str.substr(0,found);
                if(prefix.length()>0 && suffix.length()>0)
                {
                    if(suffix.compare("dcm")==0)
                    {
                        fileNames.push_back(dir->d_name);
                    }
                }
            }
        }
        closedir(d);
    }

    ///find first file with negative b (if no B_value in dicom header, the method getgetBValue returns -1)
    C_ReadDicomDiffInfo* toolREAD_INFO_DCM;
    double b;
    for(unsigned int i=0 ; i<fileNames.size() ; i++)
    {
        toolREAD_INFO_DCM = new C_ReadDicomDiffInfo(fileNames.at(i));
        b = toolREAD_INFO_DCM->getBValue();
        if(b<0.1)
        {
            ROI_DCM = fileNames.at(i);
        }
        delete toolREAD_INFO_DCM;
    }
    fileNames.clear();

    return ROI_DCM;
}






///ANA
vector<string> C_sort_files::getANAFileName(string dirName)
{
    vector<string> fileNames;
    DIR *d;
    dirent *dir;
    d = opendir(dirName.data());
    if (d)
    {
        while ((dir = readdir(d)) != NULL)
        {
            string str = dir->d_name;
            if(str.length()>4)
            {
                unsigned found = str.find_last_of(".");
                string suffix = str.substr(found+1);
                string prefix = str.substr(0,found);
                if(prefix.length()>0 && suffix.length()>0)
                {
                    if(suffix.compare("img")==0)
                    {
                        //if(prefix.substr(0,2).compare("IM")==0) ///special selection: try to find a special patern
                        //{
                            ///special patern 389_09_20090630_001_012_ep2d_diff_stejskal_mddw_12_p2_b0_1000_ep_b if
                            //unsigned short num = atoi(str.substr(found-4,4).data());
                            //if(prefix.compare("ROI")!=0 && prefix.compare("ROI2")!=0)
                            //{
                                fileNames.push_back(prefix); /**store only base name*/  //(dir->d_name);
                            //}
                        //}

                    }
                }
            }
        }
        closedir(d);
    }

    ///sort filenames those files do not obviously contain numbers
    //sort(fileNames.begin(), fileNames.end(), numericalSorting);
    return fileNames;
}
vector<string> C_sort_files::getS0FileNameANA(vector<string> fileNames, unsigned long NgradientDirection, unsigned long Nrepetition, string patternInFileNames)
{
    ///we assume that the final part of the names are of form: XXX_YYYY where XXX is a 3 digit number (number of repetition) and YYYY a 4 digit number (YYYYth acquiisition)
    ///and the first part of the name is patternInFileNames

    ///filter names
    vector<string> fileNamesWithPattern;
    for(unsigned int i=0 ; i<fileNames.size() ; i++)
    {
        if( fileNames.at(i).compare(0, patternInFileNames.length(), patternInFileNames)==0 ) ///if the first part of the filename is patternInFileNames
        {
            fileNamesWithPattern.push_back(fileNames.at(i));
        }
    }

    ///sort filenames
    sort(fileNamesWithPattern.begin(), fileNamesWithPattern.end(), numericalSortingNoSuffix);

    ///create storage
    vector< string > SoFile(Nrepetition,"");

    unsigned long j=0;
    for(unsigned int i=0 ; i<fileNamesWithPattern.size() ; i++)
    {
        if((i%(NgradientDirection+1))==0)
        {
            SoFile.at(j) = fileNamesWithPattern.at(i);
            j++;
        }
        if(j==Nrepetition)
        {
            i = fileNamesWithPattern.size();
        }
    }
    return SoFile;
}
vector<vector<string> > C_sort_files::getSiFileNameANA(vector<string> fileNames, unsigned long NgradientDirection, unsigned long Nrepetition, string patternInFileNames)
{
    ///we assume that the final part of the names are of form: XXX_YYYY where XXX is a 3 digit number (number of repetition) and YYYY a 4 digit number (YYYYth acquiisition)
    ///and the first part of the name is patternInFileNames

    ///filter names
    vector<string> fileNamesWithPattern;
    for(unsigned int i=0 ; i<fileNames.size() ; i++)
    {
        if( fileNames.at(i).compare(0, patternInFileNames.length(), patternInFileNames)==0 ) ///if the first part of the filename is patternInFileNames
        {
            fileNamesWithPattern.push_back(fileNames.at(i));
        }
    }

    ///sort filenames
    sort(fileNamesWithPattern.begin(), fileNamesWithPattern.end(), numericalSortingNoSuffix);

    ///create storage
    vector< /**gradient direction*/ vector< /**averaging*/string > > dataFile(NgradientDirection, vector< /**averaging*/string >(Nrepetition,"") );
    unsigned long j=0;
    for(unsigned int i=0 ; i<fileNamesWithPattern.size() ; i++)
    {
        if((i%(NgradientDirection+1))!=0)
        {
            unsigned int nr = floor(((double)i)/((double)(NgradientDirection+1)));
            unsigned int ng = i-nr*(NgradientDirection+1)-1;
            dataFile.at(ng).at(nr) = fileNamesWithPattern.at(i);
            j++;
        }
        if(j==(NgradientDirection*Nrepetition))
        {
            i = fileNamesWithPattern.size();
        }
    }
    return dataFile;
}


