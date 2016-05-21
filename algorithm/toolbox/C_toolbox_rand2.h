#ifndef C_TOOLBOX_RAND2_H
#define C_TOOLBOX_RAND2_H

#include <C_toolbox.h>

#define THE_NUMBER_WHICH_MUST_NOT_BE_NAMED 4101842887655102017LL


class C_toolbox_rand2 : public C_toolbox
{
    public:
        /** Default constructor */
        //C_rand2();
        /** Default destructor */
        virtual ~C_toolbox_rand2();


        unsigned long long int u,v,w;

        C_toolbox_rand2(unsigned long long int j);/* : v(4101842887655102017LL), w(1) {
        u = j ^ v; int64();
        v = u; int64();
        w = v; int64();
        }*/
        /*inline*/ unsigned long long int int64(); /*{
        u = u * 2862933555777941757LL + 7046029254386353087LL;
        v ^= v >> 17; v ^= v << 31; v ^= v >> 8;
        w = 4294957665U*(w & 0xffffffff) + (w >> 32);
        unsigned long long int x = u ^ (u << 21); x ^= x >> 35; x ^= x << 4;
        return (x + v) ^ w;
        }*/
        /*inline*/ double doub(); /*{ return 5.42101086242752217E-20 * int64(); }*/
        /*inline*/ unsigned int int32(); /*{ return (unsigned int)int64(); }*/
};

#endif // C_TOOLBOX_RAND2_H
