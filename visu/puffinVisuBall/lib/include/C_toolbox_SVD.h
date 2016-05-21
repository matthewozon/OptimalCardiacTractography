#ifndef C_TOOLBOX_SVD_H
#define C_TOOLBOX_SVD_H

#include <C_toolbox.h>
//#include <math.h>
#include <iostream>
//#include <settings_def.h>

/************************************************************
 *                                                          *
 *  Permission is hereby granted  to  any  individual   or  *
 *  institution   for  use,  copying, or redistribution of  *
 *  this code and associated documentation,  provided       *
 *  that   such  code  and documentation are not sold  for  *
 *  profit and the  following copyright notice is retained  *
 *  in the code and documentation:                          *
 *     Copyright (c) held by Dianne Cook                    *
 *  All Rights Reserved.                                    *
 *                                                          *
 *  Questions and comments are welcome, and I request       *
 *  that you share any modifications with me.               *
 *                                                          *
 *                Dianne Cook                               *
 *             dicook@iastate.edu                           *
 *                                                          *
 ************************************************************/
/**http://www.public.iastate.edu/~dicook/JSS/paper/code/defs_and_types.h*/

#define PRECISION1 32768
#define PRECISION2 16384
/*#define PI 3.1415926535897932*/
//#define MIN(x,y) ( (x) < (y) ? (x) : (y) )
//#define MAX(x,y) ((x)>(y)?(x):(y))
//#define SIGN(a, b) ((b) >= 0.0 ? fabs(a) : -fabs(a))
//#define ABS(x) ((x)>(0)?(x):(-x))
#define MAXINT 2147483647
#define ASCII_TEXT_BORDER_WIDTH 4
#define MAXHIST 100
#define STEP0 0.01
#define FORWARD 1
#define BACKWARD -1
#define PROJ_DIM 5
#define True 1
#define False 0


typedef struct {
	float x, y, z;
} fcoords;

typedef struct {
	long x, y, z;
} lcoords;

typedef struct {
	int x, y, z;
} icoords;

typedef struct {
	float min, max;
} lims;


/* grand tour history */
typedef struct hist_rec {
  struct hist_rec *prev, *next;
  float *basis[3];
  int pos;
} hist_rec;

using namespace std;

class C_toolbox_SVD : public C_toolbox
{
    public:
        /** Default constructor */
        C_toolbox_SVD();
        /** Default destructor */
        virtual ~C_toolbox_SVD();
        double PYTHAG(double a, double b);
        long int dsvd(double **a, long int m, long int n, double *w, double **v);
        double** transpose(double** a, long int m, long int n);
        double** matProd(double** A, long int ma, long int na, double** B, long int mb, long int nb);
        double** diag(double*A, long int n);
        double** pseudoInv(double** A, long int m, long int n);
        long int* sortVW(double*m_w, long int n);
};

#endif // C_TOOLBOX_SVD_H
