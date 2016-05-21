#ifndef C_WRITEDICOMDIFFINFO_H
#define C_WRITEDICOMDIFFINFO_H

#include <gdcmImage.h>

#include <struct.h>
#include <settings_def.h>
#include <iostream>
//#include <gdcmReader.h>
//#include <gdcmImageReader.h>
#include <gdcmWriter.h>
#include <gdcmDataSet.h>
#include <gdcmAttribute.h>
#include <gdcmCSAHeader.h>
#include <gdcmDictPrinter.h>



using namespace std;


///example from
///http://gdcm.sourceforge.net/html/GenFakeImage_8cxx-example.html

//http://gdcm.sourceforge.net/html/rle2img_8cxx-example.html
//http://gdcm.sourceforge.net/html/classgdcm_1_1StreamImageWriter.html#a9b6dd175ff5cfd0d875eec0b0931ed1d
//https://github.com/malaterre/GDCM/blob/master/Examples/Cxx/StreamImageReaderTest.cxx


class C_WriteDicomDiffInfo
{
    public:
        /** Default constructor */
        C_WriteDicomDiffInfo(/*string fileName*/);
        /** Default destructor */
        virtual ~C_WriteDicomDiffInfo();

        //dimension of the acquisition matrix: external
        void setNumberOfDimensions(unsigned short N);
        void setDimX(unsigned long dimx);
        void setDimY(unsigned long dimy);
        void setDimZ(unsigned long dimz);
        void setDimT(unsigned long dimt);

        //pixel encoding
        void setBitsAllocated(unsigned short n); ///ok
        void setPixelRepresentation(unsigned short DICOM_REP); ///ok
        //void signedPixelRepresentation(bool signedRepresentation); ///implemented

        //MONOCHROME OR RGB OR ...
        void setPhotometricInterpretation(string photo); ///ok


        //diffusion related data
        void setDiffInfo(diffInfo* dfi);


        protected:
        //acquisition set up
        void setPatientBasisVector(double* X_slice, double* Y_slice); ///pb with tag
        void setGradientDirectionSliceBasis(double* G); ///pb with tag

        //double** getBMatrix();
        void setBValue(double b); ///pb with tag
        void setEchoTime(double Te); ///ok
        void setRepetitionTime(double Tr); ///ok


        //distance between pixels in x, y and z direction
        void setPixelDimX(double dx); ///pb
        void setPixelDimY(double dy); ///pb
        void setPixelDimZ(double dz); ///ok

        gdcm::DataSet* ds; ///the first object to create and fill
        gdcm::File* mFile; ///SetDataSet(&ds)
        gdcm::Image* mImage; ///SetFile(mFile)
        gdcm::FileMetaInformation *mMetaInfo;

};

#endif // C_WRITEDICOMDIFFINFO_H
