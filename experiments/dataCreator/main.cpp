#define _NOTDEF
#ifdef _NOTDEF
#include <string.h>
#include <sstream>
#include <iostream>
#include <C_tensorMaker.h>
#include <C_synthetic_DWI_cylinder.h>
#include <settings_volume.h>
#include <C_ReadAnalyze75.h>

//#define SAVE_DWI

using namespace std;


int main(int argc, char *argv[])
{
    if(argc!=3)
    {
        cout << "the command should be run like: " << argv[0] << " SNR<double> expTag<string>" << endl;
        return 25;
    }
    ///allocate tool for synthetic data creation
    C_synthetic_DWI_cylinder* toolDWI = new C_synthetic_DWI_cylinder(XMIN, XMAX, YMIN, YMAX, ZMIN, ZMAX, RMIN, RMAX, ALPHAMIN, ALPHAMAX);

    ///create structures to store DWI information (gradient direction, b-value... cf struct definition)
    vector<diffInfo*> b0Info;
    b0Info.push_back(toolDWI->createB0Info());
    vector< vector<diffInfo*> > dfi = toolDWI->createDirectionsWithRepetion(); //create six direction
    ///actually create data
    vector<rawData<double>*> b0Data;
    b0Data.push_back(toolDWI->createB0Data(b0Info.at(0))); // create only one b0, because there is no need for more in synthetic case
    rawData<double>* segRawData = new rawData<double>(DICOM_DOUBLE, 3, b0Data.at(0)->DimX, b0Data.at(0)->DimY, b0Data.at(0)->DimZ);// allocate mask
    (*segRawData) = (*(b0Data.at(0))>0.0); ///for real data the creation of the mask is different. Here, is S0 is strictly positive, it is in mask
    segRawData->pixDimX = 1.0;
    segRawData->pixDimY = 1.0;
    segRawData->pixDimZ = 1.0;
    segRawData->pixDimT = 1.0;


    vector< vector<rawData<double>*> > dwiData = toolDWI->createDWIdata(dfi, b0Data);

    ///add, or not, rician noise
    if(strcmp(argv[1],"INF")!=0)
    {
        double SNR = atof(argv[1]);
        for(unsigned int i=0 ; i<dwiData.size() ; i++)
        {
            for(unsigned int j=0 ; j<dwiData.at(i).size() ; j++)
            {
                toolDWI->addNoise(dwiData.at(i).at(j) /**slice to be corrupted*/, segRawData/**mask ROI*/, SNR);
            }
        }
        for(unsigned int i=0 ; i<b0Data.size() ; i++)
        {
            toolDWI->addNoise(b0Data.at(i) /**slice to be corrupted*/, segRawData/**mask ROI*/, SNR);
        }
    }

#ifdef SAVE_DWI
    C_ReadAnalyze75* toolSaveDWI = new C_ReadAnalyze75();
    dsr* dsrPtrDWI = new dsr;
    dsrPtrDWI->hk.extents = 16384;
    //dsrPtr->hk.db_name = (char*) "private";
    dsrPtrDWI->hk.regular = (char) 'r';
    dsrPtrDWI->dime.datatype = DT_DOUBLE;
    dsrPtrDWI->dime.bitpix = 8*sizeof(double);
    dsrPtrDWI->dime.dim[0] = 3;
    dsrPtrDWI->dime.dim[1] = dwiData.at(0).at(0)->DimX;
    dsrPtrDWI->dime.dim[2] = dwiData.at(0).at(0)->DimY;
    dsrPtrDWI->dime.dim[3] = dwiData.at(0).at(0)->DimZ;
    dsrPtrDWI->dime.pixdim[0] = 3;
    dsrPtrDWI->dime.pixdim[1] = dwiData.at(0).at(0)->pixDimX;
    dsrPtrDWI->dime.pixdim[2] = dwiData.at(0).at(0)->pixDimY;
    dsrPtrDWI->dime.pixdim[3] = dwiData.at(0).at(0)->pixDimZ;
    if(strcmp(argv[1],"INF")==0)
    {
        toolSaveDWI->saveAna("DWI_SNR_INF", dsrPtrDWI, dwiData.at(0).at(0));
    }
    else
    {
        ostringstream NAME_DWI;
        NAME_DWI << "DWI_SNR_" << argv[1] << "_alpha_22_5_rmin_" << (int) RMIN << "_rmax_" << (int) RMAX << "_N_" << argv[2];
        toolSaveDWI->saveAna(NAME_DWI.str().data(), dsrPtrDWI, dwiData.at(0).at(0));
    }
#endif



    ///compute tensor //the method (makeTensors) computing tensor is implemented in library
    C_tensorMaker toolTensor;// = C_tensorMaker();
    rawData<double>* myTensos = toolTensor.makeTensors(dfi /**info about dwi b!=0*/,\
                                  dwiData /**dwi b!=0*/,\
                                   b0Info /**info about dwi b=0*/,\
                                    b0Data /**dwi with b=0*/,\
                                     segRawData /**a mask for ROI*/);

    ///should delete b0Data;
    ///should delete b0Info
    for(unsigned int i=0 ; i<b0Data.size() ; i++)
    {
        delete b0Data.at(i);
        delete b0Info.at(i);
    }
    b0Data.clear();
    b0Info.clear();


    ///should delete dfi
    ///should delete dwiData
    for(unsigned int i=0 ; i<dwiData.size() ; i++)
    {
        for(unsigned int j=0 ; j<dwiData.at(i).size() ; j++)
        {
            delete dwiData.at(i).at(j);
            delete dfi.at(i).at(j);
        }
    }
    dwiData.clear();
    dfi.clear();

    delete toolDWI;

    ///save data
    try
    {
        C_ReadAnalyze75* toolSave = new C_ReadAnalyze75();
        dsr* dsrPtr = new dsr;
        dsrPtr->hk.extents = 16384;
        //dsrPtr->hk.db_name = (char*) "private";
        dsrPtr->hk.regular = (char) 'r';
        dsrPtr->dime.datatype = DT_DOUBLE;
        dsrPtr->dime.bitpix = 8*sizeof(double);
        dsrPtr->dime.dim[0] = 4;
        dsrPtr->dime.dim[1] = myTensos->DimX;
        dsrPtr->dime.dim[2] = myTensos->DimY;
        dsrPtr->dime.dim[3] = myTensos->DimZ;
        dsrPtr->dime.dim[4] = myTensos->DimT;
        dsrPtr->dime.pixdim[0] = 4;
        dsrPtr->dime.pixdim[1] = myTensos->pixDimX;
        dsrPtr->dime.pixdim[2] = myTensos->pixDimY;
        dsrPtr->dime.pixdim[3] = myTensos->pixDimZ;
        dsrPtr->dime.pixdim[4] = myTensos->pixDimT;
        if(strcmp(argv[1],"INF")==0)
        {
            toolSave->saveAna("tensors_xx_xy_xz_yy_yz_zz_SNR_INF", dsrPtr, myTensos);
        }
        else
        {
            ostringstream NAME;
            NAME << "DTI_cylinder_6_grad_1_rep_SNR_" << argv[1] << "_alpha_22_5_rmin_" << (int) RMIN << "_rmax_" << (int) RMAX << "_N_" << argv[2];
            toolSave->saveAna(NAME.str().data(), dsrPtr, myTensos);
        }



        dsrPtr->dime.dim[0] = 3;
        dsrPtr->dime.dim[1] = segRawData->DimX;
        dsrPtr->dime.dim[2] = segRawData->DimY;
        dsrPtr->dime.dim[3] = segRawData->DimZ;
        dsrPtr->dime.dim[4] = 0;
        dsrPtr->dime.pixdim[0] = 3;
        dsrPtr->dime.pixdim[1] = segRawData->pixDimX;
        dsrPtr->dime.pixdim[2] = segRawData->pixDimY;
        dsrPtr->dime.pixdim[3] = segRawData->pixDimZ;
        dsrPtr->dime.pixdim[4] = 0;
        if(strcmp(argv[1],"INF")==0)
        {
            toolSave->saveAna("tensors_mask_SNR_INF", dsrPtr, segRawData);
        }
        else
        {
            ostringstream NAME;
     	    NAME << "mask_cylinder_6_grad_1_rep_SNR_" << argv[1] << "_alpha_22_5_rmin_" << (int) RMIN << "_rmax_" << (int) RMAX << "_N_" << argv[2];
            toolSave->saveAna(NAME.str().data(), dsrPtr, segRawData);
        }
        delete dsrPtr;
        delete toolSave;
    }
    catch(...)
    {
        cout << "tensors and/or mask could not be saved" << endl;
    }


    ///delete
    delete myTensos;
    delete segRawData;

    return 25;
}
#endif
