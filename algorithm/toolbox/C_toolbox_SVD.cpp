#include <C_toolbox_SVD.h>

C_toolbox_SVD::C_toolbox_SVD():C_toolbox()
{
    //ctor
}

C_toolbox_SVD::~C_toolbox_SVD()
{
    //dtor
}


double C_toolbox_SVD::PYTHAG(double a, double b)
{
    double at = fabs(a), bt = fabs(b), ct, result;

    if (at > bt)       { ct = bt / at; result = at * sqrt(1.0 + ct * ct); }
    else if (bt > 0.0) { ct = at / bt; result = bt * sqrt(1.0 + ct * ct); }
    else result = 0.0;
    return(result);
}


long int C_toolbox_SVD::dsvd(double **a, long int m, long int n, double *w, double **v)
{
    long int flag, i, its, j, jj, k, l, nm;
    double c, f, h, s, x, y, z;
    double anorm = 0.0, g = 0.0, scale = 0.0;
    double *rv1;

    if (m < n)
    {
        //fprintf(stderr, "#rows must be > #cols \n");
        cout << "#rows must be > #cols \n";
        return 0;
    }

    rv1 = new double[n];//(double *)malloc((unsigned int) n*sizeof(double));

/* Householder reduction to bidiagonal form */
    for (i = 0; i < n; i++)
    {
        /* left-hand reduction */
        l = i + 1;
        rv1[i] = scale * g;
        g = s = scale = 0.0;
        if (i < m)
        {
            for (k = i; k < m; k++)
                scale += fabs((double)a[k][i]);
            if (scale)
            {
                for (k = i; k < m; k++)
                {
                    a[k][i] = (double)((double)a[k][i]/scale);
                    s += ((double)a[k][i] * (double)a[k][i]);
                }
                f = (double)a[i][i];
                g = -SIGN(sqrt(s), f);
                h = f * g - s;
                a[i][i] = (double)(f - g);
                if (i != n - 1)
                {
                    for (j = l; j < n; j++)
                    {
                        for (s = 0.0, k = i; k < m; k++)
                            s += ((double)a[k][i] * (double)a[k][j]);
                        f = s / h;
                        for (k = i; k < m; k++)
                            a[k][j] += (double)(f * (double)a[k][i]);
                    }
                }
                for (k = i; k < m; k++)
                    a[k][i] = (double)((double)a[k][i]*scale);
            }
        }
        w[i] = (double)(scale * g);

        /* right-hand reduction */
        g = s = scale = 0.0;
        if (i < m && i != n - 1)
        {
            for (k = l; k < n; k++)
                scale += fabs((double)a[i][k]);
            if (scale)
            {
                for (k = l; k < n; k++)
                {
                    a[i][k] = (double)((double)a[i][k]/scale);
                    s += ((double)a[i][k] * (double)a[i][k]);
                }
                f = (double)a[i][l];
                g = -SIGN(sqrt(s), f);
                h = f * g - s;
                a[i][l] = (double)(f - g);
                for (k = l; k < n; k++)
                    rv1[k] = (double)a[i][k] / h;
                if (i != m - 1)
                {
                    for (j = l; j < m; j++)
                    {
                        for (s = 0.0, k = l; k < n; k++)
                            s += ((double)a[j][k] * (double)a[i][k]);
                        for (k = l; k < n; k++)
                            a[j][k] += (double)(s * rv1[k]);
                    }
                }
                for (k = l; k < n; k++)
                    a[i][k] = (double)((double)a[i][k]*scale);
            }
        }
        anorm = MAX(anorm, (fabs((double)w[i]) + fabs(rv1[i])));
    }

    /* accumulate the right-hand transformation */
    for (i = n - 1; i >= 0; i--)
    {
        if (i < n - 1)
        {
            if (g)
            {
                for (j = l; j < n; j++)
                    v[j][i] = (double)(((double)a[i][j] / (double)a[i][l]) / g);
                    /* double division to avoid underflow */
                for (j = l; j < n; j++)
                {
                    for (s = 0.0, k = l; k < n; k++)
                        s += ((double)a[i][k] * (double)v[k][j]);
                    for (k = l; k < n; k++)
                        v[k][j] += (double)(s * (double)v[k][i]);
                }
            }
            for (j = l; j < n; j++)
                v[i][j] = v[j][i] = 0.0;
        }
        v[i][i] = 1.0;
        g = rv1[i];
        l = i;
    }

    /* accumulate the left-hand transformation */
    for (i = n - 1; i >= 0; i--)
    {
        l = i + 1;
        g = (double)w[i];
        if (i < n - 1)
            for (j = l; j < n; j++)
                a[i][j] = 0.0;
        if (g)
        {
            g = 1.0 / g;
            if (i != n - 1)
            {
                for (j = l; j < n; j++)
                {
                    for (s = 0.0, k = l; k < m; k++)
                        s += ((double)a[k][i] * (double)a[k][j]);
                    f = (s / (double)a[i][i]) * g;
                    for (k = i; k < m; k++)
                        a[k][j] += (double)(f * (double)a[k][i]);
                }
            }
            for (j = i; j < m; j++)
                a[j][i] = (double)((double)a[j][i]*g);
        }
        else
        {
            for (j = i; j < m; j++)
                a[j][i] = 0.0;
        }
        ++a[i][i];
    }

    /* diagonalize the bidiagonal form */
    for (k = n - 1; k >= 0; k--)
    {                             /* loop over singular values */
        for (its = 0; its < 30; its++)
        {                         /* loop over allowed iterations */
            flag = 1;
            for (l = k; l >= 0; l--)
            {                     /* test for splitting */
                nm = l - 1;
                if (fabs(rv1[l]) + anorm == anorm)
                {
                    flag = 0;
                    break;
                }
                if (fabs((double)w[nm]) + anorm == anorm)
                    break;
            }
            if (flag)
            {
                c = 0.0;
                s = 1.0;
                for (i = l; i <= k; i++)
                {
                    f = s * rv1[i];
                    if (fabs(f) + anorm != anorm)
                    {
                        g = (double)w[i];
                        h = PYTHAG(f, g);
                        w[i] = (double)h;
                        h = 1.0 / h;
                        c = g * h;
                        s = (- f * h);
                        for (j = 0; j < m; j++)
                        {
                            y = (double)a[j][nm];
                            z = (double)a[j][i];
                            a[j][nm] = (double)(y * c + z * s);
                            a[j][i] = (double)(z * c - y * s);
                        }
                    }
                }
            }
            z = (double)w[k];
            if (l == k)
            {                  /* convergence */
                if (z < 0.0)
                {              /* make singular value nonnegative */
                    w[k] = (double)(-z);
                    for (j = 0; j < n; j++)
                        v[j][k] = (-v[j][k]);
                }
                break;
            }
            if (its >= 30) {
                //free((void*) rv1);
                delete rv1;
                ///fprintf(stderr, "No convergence after 30,000! iterations \n");
                cout << "No convergence after 30,000! iterations \n";
                return(0);
            }

            /* shift from bottom 2 x 2 minor */
            x = (double)w[l];
            nm = k - 1;
            y = (double)w[nm];
            g = rv1[nm];
            h = rv1[k];
            f = ((y - z) * (y + z) + (g - h) * (g + h)) / (2.0 * h * y);
            g = PYTHAG(f, 1.0);
            f = ((x - z) * (x + z) + h * ((y / (f + SIGN(g, f))) - h)) / x;

            /* next QR transformation */
            c = s = 1.0;
            for (j = l; j <= nm; j++)
            {
                i = j + 1;
                g = rv1[i];
                y = (double)w[i];
                h = s * g;
                g = c * g;
                z = PYTHAG(f, h);
                rv1[j] = z;
                c = f / z;
                s = h / z;
                f = x * c + g * s;
                g = g * c - x * s;
                h = y * s;
                y = y * c;
                for (jj = 0; jj < n; jj++)
                {
                    x = (double)v[jj][j];
                    z = (double)v[jj][i];
                    v[jj][j] = (double)(x * c + z * s);
                    v[jj][i] = (double)(z * c - x * s);
                }
                z = PYTHAG(f, h);
                w[j] = (double)z;
                if (z)
                {
                    z = 1.0 / z;
                    c = f * z;
                    s = h * z;
                }
                f = (c * g) + (s * y);
                x = (c * y) - (s * g);
                for (jj = 0; jj < m; jj++)
                {
                    y = (double)a[jj][j];
                    z = (double)a[jj][i];
                    a[jj][j] = (double)(y * c + z * s);
                    a[jj][i] = (double)(z * c - y * s);
                }
            }
            rv1[l] = 0.0;
            rv1[k] = f;
            w[k] = (double)x;
        }
    }
    //free((void*) rv1);
    delete rv1;
    return 1;
}




double** C_toolbox_SVD::transpose(double** a, long int m, long int n)
{
    double** at = new double*[n];
    for(long int i=0 ; i<n ; i++)
    {
        at[i] = new double[m];
        for(long int j=0; j<m ; j++)
        {
            at[i][j] = a[j][i];
        }
    }
    return at;
}


double** C_toolbox_SVD::matProd(double** A, long int ma, long int na, double** B, long int mb, long int nb)
{
    double** C = NULL;
    if(na==mb)
    {
        long int K = na;
        ///ma*nb
        C = new double*[ma];
        for(long int i=0 ; i< ma ; i++)
        {
            C[i] = new double[nb];
            for(long int j=0 ; j<nb ; j++)
            {
                C[i][j] = 0;
                for(long int k=0 ; k<K ; k++)
                {
                    C[i][j] += A[i][k]*B[k][j];
                }
            }
        }

    }
    return C;
}

double** C_toolbox_SVD::diag(double*A, long int n)
{
    double** v = new double*[n];
    for(long int i=0 ; i<n ; i++)
    {
        v[i] = new double[n];
        for(long int j=0 ; j<n ; j++)
        {
            v[i][j]=0;
        }
        v[i][i]=A[i];
    }
    return v;
}

double** C_toolbox_SVD::pseudoInv(double** A, long int m, long int n)
{
    double *w = new double[n];
    double** v = new double*[n];
    for(long int i=0 ; i<n ; i++)
    {
        v[i] = new double[n];
    }
    long int r = dsvd(A, m, n, w, v);
    if(r!=1)
    {
        delete w; //double *w = new double[n];
        for(long int i=0 ; i<n ; i++)
        {
            delete v[i];// = new double[n];
        }
        delete v;//double** v = new double*[n];
        return NULL;
    }

    ///Vmatrix*diag(1./Wdiag)*(Umatrix')
    for(long int i=0 ; i<n ; i++)
    {
        w[i] = 1.0/w[i];
    }
    double** Winv = diag(w, n);
    double** Ut = transpose(A, m, n);
//cout << " go fuck yourself!!!!!!!!!!!!" << endl;
    double** TEMP = matProd(v, n, n, Winv, n, n);
    double** Ainv = matProd(TEMP, n, n, Ut, n, m);

    ///
    delete TEMP;
    delete w; //double *w = new double[n];
    for(long int i=0 ; i<n ; i++)
    {
        delete v[i];// = new double[n];
    }
    delete v;//double** v = new double*[n];
    return Ainv;
}

long int* C_toolbox_SVD::sortVW(double*m_w, long int n)
{
    long int* idx = new long int[n];
    double* w = new double[n];
    for(long int i=0 ; i<n ; i++)
    {
        idx[i] = i;
        w[i] = m_w[i];
    }
    double tempD;
    long int tempI;
    long int nbSwap;
    for(long int i=0 ; i<n-1 ; i++)
    {
        nbSwap = 0;
        for(long int j=0 ; j<n-1 ; j++)
        {
            if(w[j]>w[j+1])
            {
                ///swap values
                tempD = w[j];
                w[j] = w[j+1];
                w[j+1] = tempD;

                ///swap indices
                tempI = idx[j];
                idx[j] = idx[j+1];
                idx[j+1] = tempI;

                ///increment nbSwap
                nbSwap++;
            }
        }
        if(nbSwap==0)
        {
            //cout << i << "/" << n << endl;
            i=n;
        }
    }
    delete w;
    return idx;
}
