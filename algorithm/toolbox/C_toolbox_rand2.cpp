#include <C_toolbox_rand2.h>


C_toolbox_rand2::C_toolbox_rand2(unsigned long long int j) : C_toolbox(), v(THE_NUMBER_WHICH_MUST_NOT_BE_NAMED), w(1)
{
    u = j ^ v; int64();
    v = u; int64();
    w = v; int64();
}

C_toolbox_rand2::~C_toolbox_rand2()
{
    //dtor
}

/*inline*/ unsigned long long int C_toolbox_rand2::int64()
{
    u = u * 2862933555777941757LL + 7046029254386353087LL;
    v ^= v >> 17; v ^= v << 31; v ^= v >> 8;
    w = 4294957665U*(w & 0xffffffff) + (w >> 32);
    unsigned long long int x = u ^ (u << 21); x ^= x >> 35; x ^= x << 4;
    return (x + v) ^ w;
}

/*inline*/ double C_toolbox_rand2::doub()
{
    return 5.42101086242752217E-20 * int64();
}


/*inline*/ unsigned int C_toolbox_rand2::int32()
{
    return (unsigned int)int64();
}
