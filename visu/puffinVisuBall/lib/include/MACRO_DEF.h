#ifndef MACRO_DEF_H_INCLUDED
#define MACRO_DEF_H_INCLUDED

#define MAX(a,b) ((a) > (b) ? (a) : (b))
#define MIN(a,b) ((a) >= (b) ? (b) : (a))
#define MAXABS(a,b,c) MAX(ABS(a),MAX(ABS(b),ABS(c)))
#define ABS(a) ((a) > 0 ? (a) : -(a))
//#define ABS(a) ((a) > 0 ? (a) : -(a))
#define SQR(x) ((x)*(x))
#define CUB(x) ((x)*(x)*(x))
#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))

#endif // MACRO_DEF_H_INCLUDED
