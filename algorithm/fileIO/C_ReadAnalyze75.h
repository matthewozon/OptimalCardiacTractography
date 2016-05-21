#ifndef C_READANALYZE75_H
#define C_READANALYZE75_H

#include <rawData.h>
#include <string>
#include <stdio.h>
#include <iostream>

/*    ANALYZE
TM
 Header File Format
*
*  (c) Copyright, 1986-1995
*  Biomedical Imaging Resource
*  Mayo Foundation
*
*  dbh.h
*
*  databse sub-definitions
*/
struct header_key   /* header key    */
{                                 /* off + size */
      int sizeof_hdr;   /* 0 + 4            */ //byte size of header
      char data_type[10];  /* 4 + 10           */
      char db_name[18];  /* 14 + 18          */
      int extents;   /* 32 + 4           */ //should be 16384
      short int session_error;  /* 36 + 2           */
      char regular;    /* 38 + 1           */  //must be r to indicated that all images/volumes are of the same size
      char hkey_un0;                   /* 39 + 1   */
};                               /* total=40 bytes  */

struct image_dimension
{                                  /* off + size       */
      short int dim[8];                 /* 0 + 16           */
                        /** dim[0] = number of dim in db (usually 4)
                            dim[1] =  Image X dimension; number of pixels in an image row
                            dim[2]     Image Y dimension; number of pixel rows in slice
                            dim[3]     Volume Z dimension; number of slices in a volume
                            dim[4]     Time points, number of volumes in database
                        */
      short int unused8;                /* 16 + 2           */
      short int unused9;                /* 18 + 2           */
      short int unused10;               /* 20 + 2           */
      short int unused11;               /* 22 + 2           */
      short int unused12;               /* 24 + 2           */
      short int unused13;               /* 26 + 2           */
      short int unused14;               /* 28 + 2           */
      short int datatype;               /* 30 + 2           */ //datatype for this image set (ex: DT_FLOAT)
      short int bitpix;                 /* 32 + 2           */ //number of bits per pixel; 1, 8, 16, 32, or 64.
      short int dim_un0;                /* 34 + 2           */ //unused
      float pixdim[8];                  /* 36 + 32          */
                      /**
                           pixdim[] specifies the voxel dimensitons: Parallel array to dim[], giving real world measurements in mm. and ms.
                           pixdim[1] - voxel width in mm.
                           pixdim[2] - voxel height in mm.
                           pixdim[3] - interslice distance / slice thickness in mm.
                               ...etc
                      */
      float vox_offset;                 /* 68 + 4            */ //byte offset in the .img file at which voxels start. This value can negative to specify that the absolute value is applied for every image in the file.
      float funused1;                   /* 72 + 4            */
      float funused2;                   /* 76 + 4            */
      float funused3;                   /* 80 + 4            */
      float cal_max;                    /* 84 + 4            */ //specify the range of calibration values
      float cal_min;                    /* 88 + 4            */ //specify the range of calibration values
      float compressed;                 /* 92 + 4            */
      float verified;                   /* 96 + 4            */
      int glmax,glmin;                  /* 100 + 8           */ //The maximum and minimum pixel values for the entire database.
};                                 /* total=108 bytes  */

struct data_history
{                                  /* off + size       */
      char descrip[80];                 /* 0 + 80           */
      char aux_file[24];                /* 80 + 24         */
      char orient;                      /* 104 + 1          */
      char originator[10];              /* 105 + 10         */
      char generated[10];              /* 115 + 10         */
      char scannum[10];   /* 125 + 10        */
      char patient_id[10];             /* 135 + 10         */
      char exp_date[10];   /* 145 + 10        */
      char exp_time[10];   /* 155 + 10        */
      char hist_un0[3];    /* 165 + 3          */
      int views;                         /* 168 + 4          */
      int vols_added;                   /* 172 + 4          */
      int start_field;                  /* 176 + 4          */
      int field_skip;                   /* 180 + 4         */
      int omax, omin;                   /* 184 + 8         */ //TE/TR
      int smax, smin;                   /* 192 + 8         */ //b value, magnetic field magniyude (mT)
};

struct dsr
{
      struct header_key hk;              /* 0 + 40            */
      struct image_dimension dime;      /* 40 + 108          */
      struct data_history hist;          /* 148 + 200         */
};                                  /* total= 348 bytes */

using namespace std;

class C_ReadAnalyze75
{
    public:
        /** Default constructor */
        C_ReadAnalyze75();
        /** Default destructor */
        virtual ~C_ReadAnalyze75();

        rawData<double>* getRawData(const char* filenameData);
        rawData<double>* getUnmosaicData(const char* filenameData, unsigned short nb_slice_per_mosa);

        ///get the gradient directions stored in a file. Each row is one direction, elements of the row are separated by a single space
        vector< vector<double> > getGradientDirectionFormTextFile(string gradFileName);

        dsr* readAnalyzeInfo(const char* filenameHDR);

        bool saveAna(const char* filenameAna, dsr* ptrDSR, rawData<double>* data);

    protected:
        char* readAnalyzeData(const char* filenameData, dsr* dsrPtr, long int offsetElement, long int nbElement);
        template<typename T> rawData<T>* getRawData_(const char* filenameData, dsr* dsrPtr);
        template<typename T> rawData<T>* getUnmosaicData_(const char* filenameData, dsr* dsrPtr, unsigned short nb_slice_per_mosa);
        unsigned short getDataRepresentation(dsr* dsrPtr);
        template<typename T> rawData<double>* conver2double(rawData<T>* tempRaw);

        bool writeAnalyzeInfo(const char* filenameHDR, dsr* ptrDSR);
        bool writeAnalyzeData(const char* filenameData, rawData<double>* data);
        bool getGradientDirectionFormTextFile(string gradFileName, vector< vector<double> > grad);

        //tools
        short int getDim(dsr* dsrPtr);
        long int getNumel(dsr* dsrPtr, short int numDim);
        string getImgName(string name);
        string getHdrName(string name);
        string getDatName(string name);
        string getBaseName(string name);
        string getSuffix(string name);
};




#endif // C_READANALYZE75_H
