#include <C_toolbox_eigen_sym.h>


C_toolbox_eigen_sym::C_toolbox_eigen_sym(/**const*/ double** a, const int nb_rows, bool yesvec):C_toolbox()
{
    n = nb_rows;
    z = new double*[n];
    for(int i=0; i<n ; i++)
    {
        z[i] = new double[n];
        for(int j=0 ; j<n ; j++)
        {
            z[i][j] = a[i][j];
        }
    }
    d = new double[n];
    e = new double[n];
    yesvecs = yesvec;
    tred2();
    tqli();
    sortEig();
}
C_toolbox_eigen_sym::C_toolbox_eigen_sym(/**const*/ double* a, bool yesvec):C_toolbox()
{
    n = 3;
    z = new double*[n];
    z[0] = new double[n];
    z[1] = new double[n];
    z[2] = new double[n];
    z[0][0] = a[0]; z[0][1] = a[1]; z[0][2] = a[2];
    z[1][0] = a[1]; z[1][1] = a[3]; z[1][2] = a[4];
    z[2][0] = a[2]; z[2][1] = a[4]; z[2][2] = a[5];
    d = new double[n];
    e = new double[n];
    yesvecs = yesvec;
    tred2();
    tqli();
    sortEig();
}

C_toolbox_eigen_sym::C_toolbox_eigen_sym(const double* dd, int dim_dd, const double* ee/*, int dim_ee, bool yesvec*/):C_toolbox()
{
    n = dim_dd;
    d = new double[n];
    for(int i=0 ; i<n ; i++) d[i] = dd[i];
    e = new double[n];
    for(int i=0 ; i<n ; i++) e[i] = ee[i];
    z = new double*[n];
    for(int i=0; i<n ; i++)
    {
        z[i] = new double[n];
        for(int j=0 ; j<n ; j++){
            if(i!=j)
                z[i][j] = 0;
            else
                z[i][i]=1.0;
        }

    }
    tqli();
    sortEig();
}
C_toolbox_eigen_sym::~C_toolbox_eigen_sym()
{
    //dtor
    if(d!=NULL)
    {
        delete d;
    }
    if(e!=NULL)
    {
        delete e;
    }
    if(z!=NULL)
    {
        for(int i=0 ; i<n ; i++)
        {
            if(z[i]!=NULL)
            {
                delete z[i];
            }
        }
        delete z;
    }
}


void C_toolbox_eigen_sym::sortEig(void) { //sort descending order
    if (yesvecs)
        eigsrt(d,n,z);
    else
        eigsrt(d,n);
}

void C_toolbox_eigen_sym::eigsrt(double* d, int n, double** v)
{
	int k;

	for (int i=0;i<n-1;i++) {
		double p=d[k=i];
		for (int j=i;j<n;j++)
			if (d[j] >= p) p=d[k=j];
		if (k != i) {
			d[k]=d[i];
			d[i]=p;
			if (v != NULL)
				for (int j=0;j<n;j++) {
					p=v[j][i];
					v[j][i]=v[j][k];
					v[j][k]=p;
				}
		}
	}
}

void C_toolbox_eigen_sym::tred2(void)
{
	int l,k,j,i;
	double scale,hh,h,g,f;
	for (i=n-1;i>0;i--) {
		l=i-1;
		h=scale=0.0;
		if (l > 0) {
			for (k=0;k<i;k++)
				scale += abs(z[i][k]);
			if (scale == 0.0)
				e[i]=z[i][l];
			else {
				for (k=0;k<i;k++) {
					z[i][k] /= scale;
					h += z[i][k]*z[i][k];
				}
				f=z[i][l];
				g=(f >= 0.0 ? -sqrt(h) : sqrt(h));
				e[i]=scale*g;
				h -= f*g;
				z[i][l]=f-g;
				f=0.0;
				for (j=0;j<i;j++) {
					if (yesvecs)
						z[j][i]=z[i][j]/h;
					g=0.0;
					for (k=0;k<j+1;k++)
						g += z[j][k]*z[i][k];
					for (k=j+1;k<i;k++)
						g += z[k][j]*z[i][k];
					e[j]=g/h;
					f += e[j]*z[i][j];
				}
				hh=f/(h+h);
				for (j=0;j<i;j++) {
					f=z[i][j];
					e[j]=g=e[j]-hh*f;
					for (k=0;k<j+1;k++)
						z[j][k] -= (f*e[k]+g*z[i][k]);
				}
			}
		} else
			e[i]=z[i][l];
		d[i]=h;
	}
	if (yesvecs) d[0]=0.0;
	e[0]=0.0;
	for (i=0;i<n;i++) {
		if (yesvecs) {
			if (d[i] != 0.0) {
				for (j=0;j<i;j++) {
					g=0.0;
					for (k=0;k<i;k++)
						g += z[i][k]*z[k][j];
					for (k=0;k<i;k++)
						z[k][j] -= g*z[k][i];
				}
			}
			d[i]=z[i][i];
			z[i][i]=1.0;
			for (j=0;j<i;j++) z[j][i]=z[i][j]=0.0;
		} else {
			d[i]=z[i][i];
		}
	}
}

void C_toolbox_eigen_sym::tqli()
{
	int m,l,iter,i,k;
	double s,r,p,g,f,dd,c,b;
	const double EPS=DBL_EPSILON;
	for (i=1;i<n;i++) e[i-1]=e[i];
	e[n-1]=0.0;
	for (l=0;l<n;l++) {
		iter=0;
		do {
			for (m=l;m<n-1;m++) {
				dd=abs(d[m])+abs(d[m+1]);
				if (abs(e[m]) <= EPS*dd) break;
			}
			if (m != l) {
				if (iter++ == 90) throw("Too many iterations in tqli"); //used to be 30
				g=(d[l+1]-d[l])/(2.0*e[l]);
				r=pythag(g,1.0);
				g=d[m]-d[l]+e[l]/(g+SIGN(r,g));
				s=c=1.0;
				p=0.0;
				for (i=m-1;i>=l;i--) {
					f=s*e[i];
					b=c*e[i];
					e[i+1]=(r=pythag(f,g));
					if (r == 0.0) {
						d[i+1] -= p;
						e[m]=0.0;
						break;
					}
					s=f/r;
					c=g/r;
					g=d[i+1]-p;
					r=(d[i]-g)*s+2.0*c*b;
					d[i+1]=g+(p=s*r);
					g=c*r-b;
					if (yesvecs) {
						for (k=0;k<n;k++) {
							f=z[k][i+1];
							z[k][i+1]=s*z[k][i]+c*f;
							z[k][i]=c*z[k][i]-s*f;
						}
					}
				}
				if (r == 0.0 && i >= l) continue;
				d[l] -= p;
				e[l]=g;
				e[m]=0.0;
			}
		} while (m != l);
	}
}




double C_toolbox_eigen_sym::pythag(const double a, const double b) {
	double absa=abs(a), absb=abs(b);
	return (absa > absb ? absa*sqrt(1.0+SQR(absb/absa)) :
		(absb == 0.0 ? 0.0 : absb*sqrt(1.0+SQR(absa/absb))));
}
