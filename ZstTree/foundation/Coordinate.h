
#ifndef _FOUNDATION_COORDINATE_H
#define _FOUNDATION_COORDINATE_H

#include "BiStates.h"
#include "BaseDefine.h"

template <class NN, bool DOCHECK = true>//todo 为何使用NN作为通用类型，delay和x,y不是同一个类型啊？
class F_Point
{
public:
    F_Point() {
        x = y = 0;
        z = 0;
        delay = 0;
        t = 0 ;
    }
    NN x, y;    /* coordinate */ 
    int z;    /* number of tier */
    NN delay;
    NN t;    /* used for sorting */
};

#endif
