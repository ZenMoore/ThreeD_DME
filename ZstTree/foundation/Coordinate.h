
#ifndef _FOUNDATION_COORDINATE_H
#define _FOUNDATION_COORDINATE_H

#include "BiStates.h"
#include "BaseDefine.h"

template <class NN, bool DOCHECK = true>
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
    NN t;    /* usded for sorting */
};

#endif
