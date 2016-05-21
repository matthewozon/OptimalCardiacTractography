#include <C_point_array_data.h>

C_point_array_data::C_point_array_data(unsigned long N, unsigned long spaceDim/**=3*/, unsigned long vectorDim/**=3*/, bool yes)
{
    //ctor
    FA = NULL;
    pointCoordinate = NULL;
    vectorInterpolated = NULL;
    vector3Interpolated = NULL;
    eig1 = NULL;
    eig2 = NULL;
    eig3 = NULL;
    nb_point = 0;
    dim = 0;
    dimV = 0;
    allocateArrays(N, spaceDim, vectorDim, yes);
}

C_point_array_data::~C_point_array_data()
{
    //dtor
    deleteArrays();
}

void C_point_array_data::allocateArrays(unsigned long N, unsigned long spaceDim, unsigned long vectorDim, bool yes)
{
    deleteArrays();
    nb_point = N;
    dim = spaceDim;
    dimV = vectorDim;

    pointCoordinate = new double*[nb_point];
    if(pointCoordinate!=NULL)
    {
        if(dim>0)
        {
            for(unsigned long i=0 ; i<nb_point ; i++)
            {
                pointCoordinate[i] = new double[dim];
            }
        }
        else
        {
            for(unsigned long i=0 ; i<nb_point ; i++)
            {
                pointCoordinate[i] = NULL;
            }

        }
    }
    else
    {
        dim = 0;
    }

    //vecto init
    vectorInterpolated = new double*[nb_point];
    if(vectorInterpolated!=NULL)
    {
        if(dimV>0)
        {
            for(unsigned long i=0 ; i<nb_point ; i++)
            {
                vectorInterpolated[i] = new double[dimV];
            }
        }
        else
        {
            for(unsigned long i=0 ; i<nb_point ; i++)
            {
                vectorInterpolated[i] = NULL;
            }

        }
    }
    else
    {
        dimV = 0;
    }

    if(yes)
    {
        vector3Interpolated = new double*[nb_point];
        if(vector3Interpolated!=NULL)
        {
            if(dimV>0)
            {
                for(unsigned long i=0 ; i<nb_point ; i++)
                {
                    vector3Interpolated[i] = new double[dimV];
                }
            }
            else
            {
                for(unsigned long i=0 ; i<nb_point ; i++)
                {
                    vector3Interpolated[i] = NULL;
                }

            }
        }
    }

    //FA init
    FA = new double[nb_point];
    if(FA!=NULL)
    {
        for(unsigned long i=0 ; i<nb_point ; i++)
        {
            FA[i] = 0.0;  //init as water diffusion FA=0
        }
    }

    if(yes)
    {
        eig1 = new double[nb_point];
        if(eig1!=NULL)
        {
            for(unsigned long i=0 ; i<nb_point ; i++)
            {
                eig1[i] = 0.0;  //init as water diffusion FA=0
            }
        }

        eig2 = new double[nb_point];
        if(eig2!=NULL)
        {
            for(unsigned long i=0 ; i<nb_point ; i++)
            {
                eig2[i] = 0.0;  //init as water diffusion FA=0
            }
        }

        eig3 = new double[nb_point];
        if(eig3!=NULL)
        {
            for(unsigned long i=0 ; i<nb_point ; i++)
            {
                eig3[i] = 0.0;  //init as water diffusion FA=0
            }
        }
    }
}

void C_point_array_data::deleteArrays(void)
{
    if(pointCoordinate!=NULL)
    {
        for(unsigned long i=0 ; i<nb_point ; i++)
        {
            if(pointCoordinate[i]!=NULL)
            {
                delete pointCoordinate[i];
            }
        }
        delete pointCoordinate;
    }

    if(vectorInterpolated!=NULL)
    {
        for(unsigned long i=0 ; i<nb_point ; i++)
        {
            if(vectorInterpolated[i]!=NULL)
            {
                delete vectorInterpolated[i];
            }
        }
        delete vectorInterpolated;
    }

    if(vector3Interpolated!=NULL)
    {
        for(unsigned long i=0 ; i<nb_point ; i++)
        {
            if(vector3Interpolated[i]!=NULL)
            {
                delete vector3Interpolated[i];
            }
        }
        delete vector3Interpolated;
    }


    if(FA!=NULL)
    {
        delete FA;
    }

    if(eig1!=NULL)
    {
        delete eig1;
    }

    if(eig2!=NULL)
    {
        delete eig2;
    }

    if(eig3!=NULL)
    {
        delete eig3;
    }
}
