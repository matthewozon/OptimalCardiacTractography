#ifndef C_TOOLBOX_PCA_3D_H
#define C_TOOLBOX_PCA_3D_H

#include <C_toolbox.h>
#include <C_toolbox_eigen_sym.h>


class C_toolbox_PCA_3D : public C_toolbox
{
    public:
        /** Default constructor */
        C_toolbox_PCA_3D();
        /** Default destructor */
        virtual ~C_toolbox_PCA_3D();

        double** getPrincipalComponents3D(double** X, unsigned long N); //gaussian assumption

    protected:
    private:
};

#endif // C_TOOLBOX_PCA_3D_H
