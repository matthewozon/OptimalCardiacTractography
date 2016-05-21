#include <C_kmean2D.h>


C_kmean2D::C_kmean2D(double** mIM, unsigned long mL, unsigned long mC, bool seg)
{
    //ctor
    L = mL;
    C = mC;
    IM = mIM;
    IMSEG = NULL;
    isIMSEG = seg;
    if(isIMSEG)
    {
        IMSEG = new unsigned long*[L];
        if(IMSEG!=NULL)
        {
            for(unsigned long i=0 ; i<L ; i++)
            {
                IMSEG[i] = new unsigned long[C]; //should check if all allocations went well
            }
        }
    }

    mu = vector<double>(0,0);
    sig = vector<double>(0,0);
}

C_kmean2D::~C_kmean2D()
{
    //dtor
    if(isIMSEG)
    {
        if(IMSEG!=NULL)
        {
            for(unsigned long i=0 ; i<L ; i++)
            {
                if(IMSEG[i]!=NULL) delete IMSEG[i];
            }
            delete IMSEG;
        }
    }
    mu.clear();
    sig.clear();
}

//tools
unsigned long C_kmean2D::getRowFromIdx(unsigned long idx)
{
    //idx = c + l*C
    //return idx %L;
    return floor(idx/C);
}
unsigned long C_kmean2D::getColumnFromIdx(unsigned long idx)
{
    //idx = c + l*C
    return idx - floor(idx/C)*C;
    //return idx %L;
}
unsigned long C_kmean2D::getIdxFromRowAndColumn(unsigned long l, unsigned long c)
{
    return c + l*C;
}

double C_kmean2D::getMeanOFValuesAtIndex(vector<unsigned long> IDX) //mu = sum(IM(IDX))/(length(IDX))
{
    double m = 0.0;
    if(IDX.size()!=0)
    {
        for(unsigned long i=0 ; i<IDX.size() ; i++)
        {
//            if(getRowFromIdx(IDX.at(i))>L)
//            {
//                cout << "rows out of bound " << getRowFromIdx(IDX.at(i)) << endl;
//                getchar();
//            }
//            if(getColumnFromIdx(IDX.at(i))>C)
//            {
//                cout << "column out of bound " << getColumnFromIdx(IDX.at(i)) << endl;
//                getchar();
//            }
            m += IM[getRowFromIdx(IDX.at(i))][getColumnFromIdx(IDX.at(i))];
        }
        m = m/((double) IDX.size());

    }
    else
    {
        for(unsigned long l=0 ; l<L ; l++)
        {
            for(unsigned long c=0 ; c<C ; c++)
            {
                m += IM[l][c];
            }
        }
        m = m/((double) L*C);
    }
    if(m!=m) m = 0.0;
    if(isinf(m)) m = 0.0;
    return m;
}
double C_kmean2D::getMaxValue(void)
{
    double m = IM[0][0];
    for(unsigned long l=0 ; l<L ; l++)
    {
        for(unsigned long c=0 ; c<C ; c++)
        {
            if(m<IM[l][c])
            {
                m = IM[l][c];
            }
        }
    }
    return m;
}
double C_kmean2D::getMinValue(void)
{
    double m = IM[0][0];
    for(unsigned long l=0 ; l<L ; l++)
    {
        for(unsigned long c=0 ; c<C ; c++)
        {
            if(m>IM[l][c])
            {
                m = IM[l][c];
            }
        }
    }
    return m;
}
double C_kmean2D::getStandardDeviationOFValuesAtIndex(vector<unsigned long> IDX) //sig = sqrt( sum(SQR(IM(IDX)-mu))/(length(IDX)-1) )
{
    double m = getMeanOFValuesAtIndex(IDX);

    //getchar();
    double s = 0.0;
    if(IDX.size()!=0)
    {
        if(IDX.size()==1) return 0.0;
        for(unsigned long i=0 ; i<IDX.size() ; i++)
        {
            s += SQR(IM[getRowFromIdx(IDX.at(i))][getColumnFromIdx(IDX.at(i))] -m);
        }
        if(isinf(s) || s!=s) s = 0.0;
        s = sqrt(s/((double) IDX.size()-1));
    }
    else
    {
        for(unsigned long l=0 ; l<L ; l++)
        {
            for(unsigned long c=0 ; c<C ; c++)
            {
                s += SQR(IM[l][c]-m);
            }
        }
        s = sqrt(s/((double) L*C -1));
        if(isinf(s) || s!=s) s = 0.0;
    }

    return s;
}

//working methods
unsigned long /*label of pixelValue*/ C_kmean2D::pixelCloseTo(vector<double> X /*vector of current centroids*/, double pixelValue)
{
    unsigned long label = 0;
    double m = abs(X.at(0)-pixelValue);
    for(unsigned long i=1 ; i<X.size() ; i++)
    {
        if(abs(X.at(i)-pixelValue)<m)
        {
            m = abs(X.at(i)-pixelValue);
            label = i;
        }
    }
    return label;
}
vector<unsigned long> /*vector of pixel index closest to X[x]*/ C_kmean2D::pixelCloseTo(vector<double> X /*vector of current centroids*/, unsigned long x /*index of the centroid doncidered*/)
{
    vector<unsigned long> IDX;
    for(unsigned long l=0 ; l<L ; l++)
    {
        for(unsigned long c=0 ; c<C ; c++)
        {
            if(pixelCloseTo(X, IM[l][c])==x)
            {
                IDX.push_back(getIdxFromRowAndColumn(l,c));
            }
        }
    }
    return IDX;
}

bool C_kmean2D::run(unsigned long N /*number of centroids to create*/)
{
    if(N==0) return false;
    if(mu.size()!=0) mu.clear();
    if(sig.size()!=0) sig.clear();
    mu = vector<double>(N,0.0);
    sig = vector<double>(N,0.0);
    if(N==1)
    {
        if(isIMSEG)
        {
            for(unsigned long l=0 ; l<L ; l++)
            {
                for(unsigned long c=0 ; c<C ; c++)
                {
                    IMSEG[l][c] = 0;
                }
            }
        }

        mu.at(0) = getMeanOFValuesAtIndex(vector<unsigned long>(0,0));
        sig.at(0) = getStandardDeviationOFValuesAtIndex(vector<unsigned long>(0,0));
    }
    else
    {
        ///Kmeans
        //init centroids and index vector
        vector< vector<unsigned long> > IDX(N, vector<unsigned long>(0,0) );
        double m = getMeanOFValuesAtIndex(vector<unsigned long>(0,0));
        double s = getStandardDeviationOFValuesAtIndex(vector<unsigned long>(0,0));
        //cout << m << " " << s << endl;
        //getchar();
        //if(m+s>getMaxValue() || m-s<getMinValue())
        //{
            double sInf = max(m-s, getMinValue());
            double sUp = min(m+s, getMaxValue());
            s = sUp-sInf;
        //}
        for(unsigned long i=0 ; i<N ; i++)
        {
            mu.at(i) = m + 2.0*((((double) i)/ ((double) N-1) )-0.5)*s; //should check if it's not out of bound
            //cout << mu.at(i) << endl;
        }
        //getchar();


        //iterations
        for(unsigned short iter=0 ; iter<50 ; iter++)
        {
            //for each class, get indices of the closest pixels
            for(unsigned long i=0 ; i<N ; i++)
            {
                IDX.at(i) = pixelCloseTo(mu, i);
            }

            //for each class, updatecentroids (should check if this is still sorted and classes should be created)
            for(unsigned long i=0 ; i<N ; i++)
            {
                mu.at(i) = getMeanOFValuesAtIndex(IDX.at(i));
                sig.at(i) = getStandardDeviationOFValuesAtIndex(IDX.at(i));
                //std::cout << mu.at(i) << " " << sig.at(i) << std::endl;
            }
        }

        if(isIMSEG)
        {
            ///store segmented image
            for(unsigned int i=0 ; i<IDX.size() ; i++)
            {
                for(unsigned int k=0 ; k<IDX.at(i).size() ; k++)
                {
                    IMSEG[getRowFromIdx(IDX.at(i).at(k))][getColumnFromIdx(IDX.at(i).at(k))] = i;
                }
            }
        }

    }
    return true;

}
