#ifndef C_WRITEDICOMDIFFDATA_H
#define C_WRITEDICOMDIFFDATA_H

#include <C_WriteDicomDiffInfo.h>
#include <gdcmImageWriter.h>
#include <rawData.h>

//#define CSA_WRITE
//http://gdcm.sourceforge.net/html/HelloVizWorld_8cxx-example.html#_a10
//http://gdcm.sourceforge.net/html/classgdcm_1_1ImageWriter.html

class C_writeDicomDiffData : public C_WriteDicomDiffInfo
{
    public:
        /** Default constructor */
        C_writeDicomDiffData(string fileName);
        /** Default destructor */
        virtual ~C_writeDicomDiffData();

        //fill image with with rawData
        void fillImage(rawData<unsigned short>* imRawData);

        //write data
        bool writeData();//mImage
    protected:
        gdcm::ImageWriter* mImageWriter;
};

#endif // C_WRITEDICOMDIFFDATA_H
