#include <C_toolbox_PCA_3D.h>

C_toolbox_PCA_3D::C_toolbox_PCA_3D():C_toolbox()
{
    //ctor
}

C_toolbox_PCA_3D::~C_toolbox_PCA_3D()
{
    //dtor
}


double** C_toolbox_PCA_3D::getPrincipalComponents3D(double** X, unsigned long N)
{
    //compute mean
    double* Xmu = new double[3];
    Xmu[0] = 0.0;
    Xmu[1] = 0.0;
    Xmu[2] = 0.0;
    for(unsigned long i=0 ; i<N ; i++)
    {
        Xmu[0] += X[i][0];
        Xmu[1] += X[i][1];
        Xmu[2] += X[i][2];
    }
    Xmu[0] /= (double) N;
    Xmu[1] /= (double) N;
    Xmu[2] /= (double) N;
    //compute cov of X => 3*3 matrix
    double* covMAT = new double[6];
    covMAT[0] = 0.0; covMAT[1] = 0.0; covMAT[2] = 0.0; covMAT[3] = 0.0; covMAT[4] = 0.0; covMAT[5] = 0.0;
    for(unsigned long i=0 ; i<N ; i++)
    {
        covMAT[0] += SQR(X[i][0]-Xmu[0]); covMAT[1] += (X[i][0]-Xmu[0])*(X[i][1]-Xmu[1]); covMAT[2] += (X[i][0]-Xmu[0])*(X[i][2]-Xmu[2]);
                                          covMAT[3] += SQR(X[i][1]-Xmu[1]);               covMAT[4] += (X[i][1]-Xmu[1])*(X[i][2]-Xmu[2]);
                                                                                          covMAT[5] += SQR(X[i][2]-Xmu[2]);
    }
    delete Xmu;
    covMAT[0] /= (double) N;
    covMAT[1] /= (double) N;
    covMAT[2] /= (double) N;
    covMAT[3] /= (double) N;
    covMAT[4] /= (double) N;
    covMAT[5] /= (double) N;

//    FILE* f = fopen("covMAT.data", "w");
//    fprintf(f, "%f, %f, %f\n", covMAT[0], covMAT[1], covMAT[2]);
//    fprintf(f, "%f, %f, %f\n", covMAT[1], covMAT[3], covMAT[4]);
//    fprintf(f, "%f, %f, %f\n", covMAT[2], covMAT[4], covMAT[5]);
//    fclose(f);
    //get rotation matrix P
    C_toolbox_eigen_sym* toolDIAG = new C_toolbox_eigen_sym(covMAT);
    double** P = new double*[3];
    P[0] = new double[3]; P[1] = new double[3]; P[2] = new double[3];
    P[0][0] = toolDIAG->z[0][0]; P[0][1] = toolDIAG->z[0][1]; P[0][2] = toolDIAG->z[0][2];
    P[1][0] = toolDIAG->z[1][0]; P[1][1] = toolDIAG->z[1][1]; P[1][2] = toolDIAG->z[1][2];
    P[2][0] = toolDIAG->z[2][0]; P[2][1] = toolDIAG->z[2][1]; P[2][2] = toolDIAG->z[2][2];
    delete toolDIAG;
    delete covMAT;
    return P;
}
