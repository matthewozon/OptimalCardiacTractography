//for memory check
#include "stdlib.h"
#include "stdio.h"
#include "string.h"
#include <iostream>
#include <dirent.h> //get file in directory
#include <algorithm> //sorting algorithm

//for the tractography
#include <C_readDicomDiffData.h> //#lib fileIO
#include <C_tensorMaker.h> //#lib syntheticData
#include <C_sort_files.h>  //#lib fileIO
#include <C_maskMaker.h>   //#lib syntheticData
#include <C_toolbox_interpolate.h> //#lib tractPuff
#include <C_energy.h> //#lib tractPuff
#include <C_temperature.h> //#lib tractPuff
#include <C_SA.h> //#lib tractPuff
#include <C_bundle.h> //#lib tractPuff


//for more clarity in the code:
//loading and computation of the data
struct tractData
{
    rawData<double>* segRawData;
    rawData<double>* myTensos;
};
tractData loadDataAndComputeTensorOfNewBornFromDicom(std::string);
//for the virtual memory allocated for the process
int parseLine(char*);
int getValue(void);


//it is shamefull, but it happens that I use that:
using namespace std;


//the code aims to load data from a directory (it is supposed that there are 192 diffusion directions and only one repetition per direction) and compute the tensors and mask for the ROI
//then, it initializes the energy and minimization process, and finaly compute the minimzation and saves the fibers from the optimal configuration of the graph
int main(int argc, char *argv[])
{
    if(argc==1) ///DICOM
    {
        //load tensors and ROI
	tractData myTractData = loadDataAndComputeTensorOfNewBornFromDicom(".");//, segRawData, myTensos);
	rawData<double>* segRawData = myTractData.segRawData;
	rawData<double>* myTensos = myTractData.myTensos;
	cout << "segRawData " << segRawData << endl;
	cout << "myTensos " << myTensos << endl;

        if(myTensos!=NULL && segRawData!=NULL)
        {
	    cout << "HELLO" << endl;
            ///get memory before allocating and initiating graph+energy
            int curMem = getValue();
            cout << endl;
            cout << "vMem = " << curMem << endl;
            C_graph* G = new C_graph(segRawData);
            C_energy* U = new C_energy(G, myTensos, segRawData);

            ///energy set up
            U->METHOD_INTERPOLATION = GAUSSIAN_PHI;
            U->OBJECT_INTERPOLATED = TENSOR;
            U->METHOD_INTEGRATION = SIMPSON_RULE;
            U->m_nbInterpolatingPoint = 100;
            U->m_alpha = 0.2;
            U->m_beta = 0.5;
            U->m_sigma = 0.375;
            U->m_phi0 = PHI0;
            U->m_phi1 = PHI1;
            U->m_t0 = Pi/2.0;
            U->m_R = 3.0;

            ///parameter
            unsigned long energyKind = U15; //if U1 is used, the tensors must be normalized by the largest eigenvalue

            ///init energy
            unsigned long nb_thread = 8; //defines the number of thread that will calculate the data fidelity term for each pair of the graph (use as much as possible), but it is not used in the SA process
            U->initEnergy(energyKind, nb_thread);//compute data fidelity term for each pair of edges

            U->initPointData(energyKind);//set the actual energy for each vertex (using what was calculated for the pair of edges)
            ///get memory before allocating and initiating graph+energy
            int endMem = getValue();
            cout << "vMem = " << endMem << endl;
            cout << "diff vMem = " << endMem -curMem << endl; //show the total memory that is used by the process


            ///go for minimization
            unsigned long coolingKind = SA_EXP_BLOCK_COOLING; //this is the default value for the cooling schedule: exponential decrease with constant steps
            unsigned long kindOfMove = SA_MOVE_1; //what type of communication kernel is used

            ///temperature assessment
            C_temperature* T = new C_temperature(G, U, energyKind);
	    //the acceptance rates are slightly different regarding the energy that is used
            double startingRate;
            double endRate;
            if(energyKind==U1)
            {
                startingRate = 0.3;
                endRate = 0.0030;
            }
            else
            {
                startingRate = 0.7;
                endRate = 0.001;
            }
	    //compute the initial and final temperatures for the given communication kernels and rates
            T->runTemperature(kindOfMove, kindOfMove, startingRate, endRate, 1000/*number of steps for the uphill generation*/);
            double Tinit = T->getTinit();
            double Tend = T->getTend();
            delete T;

	    //if the temperature are suitable for the process (should check for NaN (isnan) and Inf (isinf))
	    if(Tinit>0.0 && Tend>0.0)
            {
		///run minimization
                unsigned long nbIter = 64;//000;//which is actually the number of constant steps (number of iteration=nbStep*nbEdges)
                C_SA minimization_by_SA(G, U, Tinit, Tend, nbIter, coolingKind, energyKind, kindOfMove);// = new C_SA(G, U, Tinit, Tend, nbIter, coolingKind, energyKind, kindOfMove);

		//some tag that will appear at the end of the filenames
                minimization_by_SA.m_expTAG = "TEST_1";

                //init graph (draw at random a configuration with a given percentage of connected edges)
                minimization_by_SA.randInit(0.2);

		//those two variables indicates if the acceptance rate and the energy are recorded during the process. It should always be set to true because it is the only experimental proof of convergence of the algorithm. If the acceptance rate does not reach the expected value or reaches it too early, the minimization is probably not optimal. For the enregy, it is interesting to record the values calculated by the formula and by variation because if both are the same, it ensures that both methods calculates the same energy (mostly used for debugging, but always recorded just in case...)
                minimization_by_SA.flipRate = true;
                minimization_by_SA.energyMeasure = true;
                //actually comoute the minimization
                minimization_by_SA.run();

		//save the graph configuration in case we want to use it for a warm restart. But I never used it in 4 years, though it might be usefull for simulation on the cluster.
                ostringstream oss;
                oss << "the_experiment_name" << ".cformat";
                G->saveEdges((char*) oss.str().data());

                ///save fibers
                C_bundle s(G, true);// = new C_bundle(G, true); //extract fibers from the graph configuration and arrange in fiber structures (they should be modified because of the complexity regarding their use)
                minimization_by_SA.saveFiber("fibers_", s.FIBERS); //save fiber in file (on file for the header .fibHDR and one for the sources .fibSRC)

                ///free memory
                //delete s;
                //delete minimization_by_SA;
                delete segRawData;
	        delete myTensos;
                delete U;
                delete G;
            }
        }
    }
    return 60;
}




tractData loadDataAndComputeTensorOfNewBornFromDicom(std::string dirName) //, rawData<double>* segRawData, rawData<double>* myTensos)
{
	//get all files needed for DTI computation
        C_sort_files sortFiles;// = new C_sort_files();//just a tool that helps for sorting files in directories

        ///store filenames
        unsigned long NgradientDirection = 192, Nrepetition = 1;
        //string dirName = ".";
        vector< /**gradient direction*/ vector< /**averaging*/string > > dataFile = sortFiles.getSiFileNameDCM(sortFiles.getDicomFileName(dirName), NgradientDirection, Nrepetition);
        vector< string > SoFile = sortFiles.getS0FileNameDCM(sortFiles.getDicomFileName(dirName), NgradientDirection, Nrepetition);
        //delete sortFiles;


        ///load S0 info and data
        vector<rawData<double>*> allB0Data(Nrepetition,NULL);// = new rawData<double>*[Nrepetition];
        vector<diffInfo*> myB0InfoVect(Nrepetition,NULL);
        double** tempDoubleTableB0;//temporary array for the rotation matrix (from the MRI machine reference to the patient/image reference)
        double* tempDoubleArrayB0;//temporary array for the gradient directions
        for(unsigned int j=0 ; j<SoFile.size() ; j++)
        {
            C_readDicomDiffData m(SoFile.at(j).data());// = new C_readDicomDiffData(SoFile.at(j).data());
            allB0Data.at(j) = m.getUnmosaicData(m.getNumberOfSlice());
            myB0InfoVect.at(j) = new diffInfo;
            tempDoubleArrayB0 = m.getGradientDirectionSliceBasis();
            myB0InfoVect.at(j)->g[0] = tempDoubleArrayB0[0];    myB0InfoVect.at(j)->g[1] = tempDoubleArrayB0[1];    myB0InfoVect.at(j)->g[2] = tempDoubleArrayB0[2];
            tempDoubleTableB0 = m.getChangeBasisMatrix();
            myB0InfoVect.at(j)->X_slice[0] = tempDoubleTableB0[0][0];   myB0InfoVect.at(j)->X_slice[1] = tempDoubleTableB0[0][1];   myB0InfoVect.at(j)->X_slice[2] = tempDoubleTableB0[0][2];
            myB0InfoVect.at(j)->Y_slice[0] = tempDoubleTableB0[1][0];   myB0InfoVect.at(j)->Y_slice[1] = tempDoubleTableB0[1][1];   myB0InfoVect.at(j)->Y_slice[2] = tempDoubleTableB0[1][2];
            myB0InfoVect.at(j)->Z_slice[0] = tempDoubleTableB0[2][0];   myB0InfoVect.at(j)->Z_slice[1] = tempDoubleTableB0[2][1];   myB0InfoVect.at(j)->Z_slice[2] = tempDoubleTableB0[2][2];
            myB0InfoVect.at(j)->b = m.getBValue();
            myB0InfoVect.at(j)->echoTime = m.getEchoTime();
            myB0InfoVect.at(j)->repTime = m.getRepetitionTime();
            myB0InfoVect.at(j)->voxelSize[0] = m.getPixelDimX();
            myB0InfoVect.at(j)->voxelSize[1] = m.getPixelDimY();
            myB0InfoVect.at(j)->voxelSize[2] = m.getPixelDimZ();
            delete tempDoubleArrayB0;
            delete tempDoubleTableB0[0];
            delete tempDoubleTableB0[1];
            delete tempDoubleTableB0[2];
            delete tempDoubleTableB0;
            //delete m;
        }




        ///load Si info and data
        vector<vector<rawData<double>*> > allRawData( NgradientDirection, vector<rawData<double>*>(Nrepetition,NULL) );// = new rawData<double>*[NgradientDirection*Nrepetition];
        vector< vector<diffInfo*> > myDiffInfoVect( NgradientDirection, vector<diffInfo*>(Nrepetition,NULL) );
        double** tempDoubleTable;
        double* tempDoubleArray;
        for(unsigned int i=0 ; i<dataFile.size() ; i++)
        {
            for(unsigned int j=0 ; j<dataFile.at(i).size() ; j++)
            {
                C_readDicomDiffData m(dataFile.at(i).at(j).data());// = new C_readDicomDiffData(dataFile.at(i).at(j).data());
                allRawData.at(i).at(j) = m.getUnmosaicData(m.getNumberOfSlice());
                myDiffInfoVect.at(i).at(j) = new diffInfo;
                tempDoubleArray = m.getGradientDirectionSliceBasis();
                myDiffInfoVect.at(i).at(j)->g[0] = tempDoubleArray[0];    myDiffInfoVect.at(i).at(j)->g[1] = tempDoubleArray[1];    myDiffInfoVect.at(i).at(j)->g[2] = tempDoubleArray[2];
                tempDoubleTable = m.getChangeBasisMatrix();
                myDiffInfoVect.at(i).at(j)->X_slice[0] = tempDoubleTable[0][0];   myDiffInfoVect.at(i).at(j)->X_slice[1] = tempDoubleTable[0][1];   myDiffInfoVect.at(i).at(j)->X_slice[2] = tempDoubleTable[0][2];
                myDiffInfoVect.at(i).at(j)->Y_slice[0] = tempDoubleTable[1][0];   myDiffInfoVect.at(i).at(j)->Y_slice[1] = tempDoubleTable[1][1];   myDiffInfoVect.at(i).at(j)->Y_slice[2] = tempDoubleTable[1][2];
                myDiffInfoVect.at(i).at(j)->Z_slice[0] = tempDoubleTable[2][0];   myDiffInfoVect.at(i).at(j)->Z_slice[1] = tempDoubleTable[2][1];   myDiffInfoVect.at(i).at(j)->Z_slice[2] = tempDoubleTable[2][2];
                myDiffInfoVect.at(i).at(j)->b = m.getBValue();
                myDiffInfoVect.at(i).at(j)->echoTime = m.getEchoTime();
                myDiffInfoVect.at(i).at(j)->repTime = m.getRepetitionTime();
                myDiffInfoVect.at(i).at(j)->voxelSize[0] = m.getPixelDimX();
                myDiffInfoVect.at(i).at(j)->voxelSize[1] = m.getPixelDimY();
                myDiffInfoVect.at(i).at(j)->voxelSize[2] = m.getPixelDimZ();
                delete tempDoubleArray;
                delete tempDoubleTable[0];
                delete tempDoubleTable[1];
                delete tempDoubleTable[2];
                delete tempDoubleTable;
                //delete m;
            }
        }

	//create container for the pair mask and tensors
	tractData myTractData;
        ///creating mask from several gradient directions
        C_maskMaker maskTool;// = new C_maskMaker();
        myTractData.segRawData = maskTool.createMaskFromDWIs(allRawData, NULL, NO_SAMPLING);
        //delete maskTool;




        ///calculate tensors from the arranged data: directions, b-value, S0 and diffusion weighted measurements
        C_tensorMaker myTensorMaker;// = new C_tensorMaker();
        myTractData.myTensos = myTensorMaker.makeTensors(myDiffInfoVect, allRawData, myB0InfoVect, allB0Data, myTractData.segRawData);
        if(myTractData.myTensos->raw4D!=NULL)
        {
            cout << "4D data were created with dimx=" << myTractData.myTensos->DimX <<  " dimy=" << myTractData.myTensos->DimY <<  " dimz=" << myTractData.myTensos->DimZ <<  " dimt=" << myTractData.myTensos->DimT << endl;
        }
        //delete myTensorMaker;

        ///can delete myDiffInfoVect, allRawData, myB0InfoVect, allB0Data
        //cout << "deleting" << endl;
        for(unsigned int i=0 ; i<dataFile.size() ; i++)
        {
            for(unsigned int j=0 ; j<dataFile.at(i).size() ; j++)
            {
                delete myDiffInfoVect.at(i).at(j);
                delete allRawData.at(i).at(j);
            }
        }
        //cout << "still deleting" << endl;
        for(unsigned int j=0 ; j<allB0Data.size() ; j++)
        {
            delete allB0Data.at(j);
            delete myB0InfoVect.at(j);
        }
        //cout << "deleted" << endl;

        for(unsigned long i=0 ; i<myTractData.segRawData->DimX ; i++ )
        {
            for(unsigned long j=0 ; j<myTractData.segRawData->DimY ; j++ )
            {
                for(unsigned long k=0 ; k<myTractData.segRawData->DimZ ; k++ )
                {
                    swap(myTractData.segRawData->raw3D[i][j][k],myTractData.segRawData->raw3D[j][i][k]);
                    swap(myTractData.myTensos->raw4D[i][j][k][0],myTractData.myTensos->raw4D[j][i][k][0]);
                    swap(myTractData.myTensos->raw4D[i][j][k][1],myTractData.myTensos->raw4D[j][i][k][1]);
                    swap(myTractData.myTensos->raw4D[i][j][k][2],myTractData.myTensos->raw4D[j][i][k][2]);
                    swap(myTractData.myTensos->raw4D[i][j][k][3],myTractData.myTensos->raw4D[j][i][k][3]);
                    swap(myTractData.myTensos->raw4D[i][j][k][4],myTractData.myTensos->raw4D[j][i][k][4]);
                    swap(myTractData.myTensos->raw4D[i][j][k][5],myTractData.myTensos->raw4D[j][i][k][5]);
                }
            }
        }
	cout << "segRawData " << myTractData.segRawData << endl;
	cout << "myTensos " << myTractData.myTensos << endl;
    return myTractData;
}


//create functions that can read the allocated memory of a process in the process
int parseLine(char* line)
{
    int i = strlen(line);
    while (*line < '0' || *line > '9') line++;
    line[i-3] = '\0';
    i = atoi(line);
    return i;
}

int getValue(void) //Note: this value is in KB!
{
    FILE* file = fopen("/proc/self/status", "r");
    int result = -1;
    char line[128];


    while (fgets(line, 128, file) != NULL){
        if (strncmp(line, "VmSize:", 7) == 0){
            result = parseLine(line);
            break;
        }
    }
    fclose(file);
    return result;
}
