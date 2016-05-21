#include "C_thread_multiscale.h"

C_thread_multiscale::C_thread_multiscale(double resolution_factor, rawData<double>* maskDTI, rawData<double>* dataDTI, unsigned long NgradientDirection, unsigned long Nrepetition):C_thread()
{
    //ctor
    m_resolution_factor = resolution_factor;
    m_maskDTI = maskDTI;
    m_dataDTI = dataDTI;
    m_NgradientDirection = NgradientDirection;
    m_Nrepetition = Nrepetition;
}

C_thread_multiscale::~C_thread_multiscale()
{
    //dtor
}


void C_thread_multiscale::execute()
{
    C_sampling<double>* toolSAMPLING = new C_sampling<double>();
    samplingFactor* m_SFGraph = new samplingFactor();
    ///fill sampling factors so that you can have all graph resolution
    m_SFGraph->Xfactor = m_resolution_factor;
    m_SFGraph->Yfactor = m_resolution_factor;
    m_SFGraph->Zfactor = m_resolution_factor;
    m_SFGraph->Tfactor = 1.0;
    ///generate mask for the given resolution
    rawData<double>* maskTemp = toolSAMPLING->samplingNearestNeighbor(m_maskDTI, m_SFGraph);
    (*maskTemp) = (*maskTemp)>0.5;
    delete m_SFGraph;
    ///create graph at given resolution
//    cout << "DTI: " << m_dataDTI->DimX << " " << m_dataDTI->DimY << " " << m_dataDTI->DimZ << " " << m_dataDTI->DimT << endl;
//    cout << "pixDim: " << m_dataDTI->pixDimX << " " << m_dataDTI->pixDimY << " " << m_dataDTI->pixDimZ << " " << m_dataDTI->pixDimT << endl;
//    cout << "makskDTI: " << m_maskDTI->DimX << " " << m_maskDTI->DimY << " " << m_maskDTI->DimZ << " " << m_maskDTI->DimT << endl;
//    cout << "pixDim: " << m_maskDTI->pixDimX << " " << m_maskDTI->pixDimY << " " << m_maskDTI->pixDimZ << " " << m_maskDTI->pixDimT << endl;
//    cout << "maksk: " << maskTemp->DimX << " " << maskTemp->DimY << " " << maskTemp->DimZ << " " << maskTemp->DimT << endl;
//    cout << "pixDim: " << maskTemp->pixDimX << " " << maskTemp->pixDimY << " " << maskTemp->pixDimZ << " " << maskTemp->pixDimT << endl;
    //return;
    time_t startTime = time(0);
    C_graph* G = new C_graph(maskTemp);
    time_t endTime = time(0);
    cout << "time for graph creation: " << difftime(endTime, startTime) * 1000.0 << " at resolution: " << m_resolution_factor  << endl;

    C_energy* U = new C_energy(G, m_dataDTI, m_maskDTI);

    ///energy set up
    U->METHOD_INTERPOLATION = GAUSSIAN_PHI;
    U->OBJECT_INTERPOLATED = TENSOR;
    U->METHOD_INTEGRATION = SIMPSON_RULE;
    U->m_nbInterpolatingPoint = 7;
    U->m_alpha = 0.2;
    U->m_beta = 0.5;
    U->m_sigma = 0.375;
    U->m_phi0 = PHI0;
    U->m_phi1 = PHI1;
    U->m_t0 = Pi/2.0;
    U->m_R = 1.5*m_dataDTI->pixDimX; //must be fix: regarding the resolution 1

    ///parameter
    unsigned long energyKind = U15;

    ///init energy
    unsigned long nb_thread = 8;
    //cout << "init energy in thread for resolution: " << m_resolution_factor << endl;
    startTime = time(0);
    U->initEnergy(energyKind, nb_thread);
    endTime = time(0);
    cout << "time for energy init: " << difftime(endTime, startTime) * 1000.0 << " at resolution: " << m_resolution_factor  << endl;
    //cout << "finish energy in thread for resolution: " << m_resolution_factor << endl;
    int endMem = getValue();
    cout << "vMem = " << endMem << " at resolution: " << m_resolution_factor << endl;
    //cout << "diff vMem = " << endMem -curMem << endl;
    //return;

    U->initPointData(energyKind);

    ///go for minimization
    unsigned long coolingKind = SA_EXP_BLOCK_COOLING;
    unsigned long kindOfMove = SA_MOVE_1;

    ///temperature assessment
    C_temperature* T = new C_temperature(G, U, energyKind);
    double startingRate;
    double endRate;
    if(energyKind==U1)
    {
        startingRate = 0.3;
        endRate = 0.0030;
    }
    else
    {
        startingRate = 0.4;//7;
        endRate = 0.001;
    }
    T->runTemperature(kindOfMove, kindOfMove, startingRate, endRate, 1000);
    double Tinit = T->getTinit();
    double Tend = T->getTend();
    delete T;


    ///run minimization
    unsigned long nbIter = 96000;
    C_SA* minimization_by_SA = new C_SA(G, U, Tinit, Tend, nbIter, coolingKind, energyKind, kindOfMove);
    ostringstream expName;
    expName << "SA_MOVE_1_FULL_HEART_ANA_JHU_Ng_" << m_NgradientDirection << "_Nr_" << m_Nrepetition << "_U15_Tinit_" << Tinit << "_Tend_" << Tend;
    expName << "_nbiter_" << nbIter << "_alpha_" << U->m_alpha << "_beta_" << U->m_beta << "_sigma_" << U->m_sigma;
    expName << "_Niterpol_" << U->m_nbInterpolatingPoint << "_radius_" << U->m_R << "_res_" << m_resolution_factor;
    minimization_by_SA->baseName = expName.str();
    minimization_by_SA->m_expTAG = "_";

    //init graph
    minimization_by_SA->randInit(0.2);

    //run
    minimization_by_SA->flipRate = true;
    minimization_by_SA->energyMeasure = true;
    startTime = time(0);
    minimization_by_SA->run();
    endTime = time(0);
    cout << "time per iteration: " << difftime(endTime, startTime) * 1000.0/((double) (nbIter+MAX_EDGES_PER_NODE)) << endl;

    G->saveEdges((char*) expName.str().data());

    ///save fibers
    C_bundle* s = new C_bundle(G, true); //extract fibers
    minimization_by_SA->saveFiber((char*) expName.str().data(), s->FIBERS); //save fiber in file

    ///free memory
    delete s;
    delete minimization_by_SA;

    ///delete mask
    delete maskTemp;
    delete toolSAMPLING;

    delete U;
    delete G;
    return;
}


int C_thread_multiscale::parseLine(char* line)
{
    int i = strlen(line);
    while (*line < '0' || *line > '9') line++;
    line[i-3] = '\0';
    i = atoi(line);
    return i;
}

int C_thread_multiscale::getValue() //Note: this value is in KB!
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
