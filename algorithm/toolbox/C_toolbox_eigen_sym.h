#ifndef C_TOOLBOX_EIGEN_SYM_H
#define C_TOOLBOX_EIGEN_SYM_H

#include <C_toolbox.h>

//#include <stdlib.h>
#include <cmath>
#include <float.h>
//#include <settings_def.h>



class C_toolbox_eigen_sym : public C_toolbox
{
    public:
        int n;
        double** z; //each column correspond to an eigenvector
        double* d, *e; //contain the eigenvalues sorted descendingly
        bool yesvecs;
        /** Default constructor */
        C_toolbox_eigen_sym(/**const*/ double** a, const int nb_rows=3, bool yesvec=true);
        C_toolbox_eigen_sym(/**const*/ double* a, bool yesvec=true); //special case for 3*3 symetric matrix
        C_toolbox_eigen_sym(const double* dd, int dim_dd, const double* ee/*, int dim_ee, bool yesvec=true*/);
        /** Default destructor */
        virtual ~C_toolbox_eigen_sym();

        void sortEig(void);
        void tred2(void);
        void tqli(void);
        double pythag(const double a, const double b);
        void eigsrt(double* d, int n, double** v=NULL);
};

#endif // C_TOOLBOX_EIGEN_SYM_H
