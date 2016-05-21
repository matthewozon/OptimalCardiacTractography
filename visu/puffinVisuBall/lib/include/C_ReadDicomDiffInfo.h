#ifndef C_READDICOMDIFFINFO_H
#define C_READDICOMDIFFINFO_H

//#include <struct.h>
#include <settings_def.h>
#include <iostream>
#include <gdcmReader.h>
#include <gdcmImageReader.h>
//#include <gdcmWriter.h>
#include <gdcmDataSet.h>
#include <gdcmAttribute.h>
#include <gdcmCSAHeader.h>
#include <gdcmDictPrinter.h>

//// relevant GE private tags
//const gdcm::DictEntry GEDictBValue( 0x0043, 0x1039, "IS", "1", "B Value of diffusion weighting" );
//const gdcm::DictEntry GEDictXGradient( 0x0019, 0x10bb, "DS", "1", "X component of gradient direction" );
//const gdcm::DictEntry GEDictYGradient( 0x0019, 0x10bc, "DS", "1", "Y component of gradient direction" );
//const gdcm::DictEntry GEDictZGradient( 0x0019, 0x10bd, "DS", "1", "Z component of gradient direction" );
//
//// relevant Siemens private tags
//const gdcm::DictEntry SiemensMosiacParameters( 0x0051, 0x100b, "IS", "1", "Mosiac Matrix Size" );
//const gdcm::DictEntry SiemensDictNMosiac( 0x0019, 0x100a, "US", "1", "Number of Images In Mosaic" );
//const gdcm::DictEntry SiemensDictBValue( 0x0019, 0x100c, "IS", "1", "B Value of diffusion weighting" );
//const gdcm::DictEntry SiemensDictDiffusionDirection( 0x0019, 0x100e, "FD", "3", "Diffusion Gradient Direction" );
//const gdcm::DictEntry SiemensDictDiffusionMatrix( 0x0019, 0x1027, "FD", "6", "Diffusion Matrix" );



using namespace std;

class C_ReadDicomDiffInfo
{
    public:
        /** Default constructor */
        C_ReadDicomDiffInfo(string fileName); //OK
        /** Default destructor */
        virtual ~C_ReadDicomDiffInfo(); //OK

        //acquisition set up
        double** getChangeBasisMatrix(); //OK
        double* getGradientDirection(); //OK in machine frame
        double* getGradientDirectionSliceBasis(); //OK

        double getBValue(); //OK
        double getEchoTime(); //OK
        double getRepetitionTime(); //OK

        //pixel encoding
        unsigned short getBitsAllocated(); //OK
        unsigned short getPixelRepresentation(); //OK
        string getSequenceVariant(); //OK but useless?
        bool signedPixelRepresentation();//OK

        //MONOCHROME OR RGB OR ...
        string getPhotometricInterpretation(); //OK

        //dimension of the images in file
        unsigned long getNumberOfRows(); //OK
        unsigned long getNumberOfColumns(); //OK

        //distance between pixels in x, y and z direction
        double* getPixelDim(); //OK
        double getPixelDimX(); //OK
        double getPixelDimY(); //OK
        double getPixelDimZ(); //OK

        //dimension of the acquisition matrix
        unsigned long getAcquisitionMatrixNumberOfRows(); //OK
        unsigned long getAcquisitionMatrixNumberOfColumns(); //OK
        unsigned short getNumberOfSlice(); //OK
        //void printMosaicMatrixSize();


    protected:
        ///spcific attributes for reading information elements.

        //readers
        gdcm::Reader reader;
        gdcm::ImageReader mImageReader;

        //containers
        gdcm::DataSet ds;
        gdcm::File mFile;
        gdcm::Image mImage;

};

//(0018,0021)	CS	SequenceVariant	1-n
//(0028,0103)	US	PixelRepresentation 1
//(0028,0106)	US/SS	SmallestImagePixelValue	1
//(0028,0107)	US/SS	LargestImagePixelValue	1
//(7FE0,0000)	UL	PixelDataGroupLength	1
//(7FE0,0010)	OW/OB	PixelData	1

#endif // C_READDICOMDIFFINFO_H
