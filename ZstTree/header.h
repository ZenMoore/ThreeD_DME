#ifndef _HEADER_H
#define _HEADER_H

#include <iostream>
#include <fstream>
#include <cstring>
#include <vector>
#include <assert.h>
#include <limits>
#include <map>
#include "foundation/BaseDefine.h"
#include "foundation/BiStates.h"
#include "foundation/Trr.h"
#include "foundation/DOUBLE.h"
using namespace std;

#define CHECK YES
#define SHOW_INFO NO

#define DBL_MAX __DBL_MAX__//todo 这是什么，下面的PUCAP_SCALE是什么
#define NIL -1

#define BUCKET_CAPACITY   240
#define MAX_N_SINKS       3101
#define MAX_N_NODES       (2*MAX_N_SINKS)
#define N_Neighbor 1

#define PUCAP_SCALE (1E+12)

typedef F_Point < DOUBLE, true > PointType ;
typedef F_Trr < DOUBLE, true > TrrType ;

class LineType
{
public:
    LineType(void) { }
    LineType(PointType p1, PointType p2) {
        assert(p1.z == p2.z);
        endpoint[0] = p1;
        endpoint[1] = p2;
    }
    bool Is_Point(void) {
        if ( (endpoint[0].x - endpoint[1].x)._ABS() + (endpoint[0].y - endpoint[1].y) < FUZZ ) 
            return YES;
        else
            return NO;
    }

public:
    PointType endpoint[2];
};

class MsType 
{
public:
    MsType() {
        npts = 0;
        dist = DOUBLE(DBL_MAX);
        side = uod = NIL;
        segment_L = segment_R = NULL;
    }
    // ~MsType() {
    //     if (segment_L) {
    //         segment_L;
    //     }
    //     if (segment_R) {
    //         delete segment_R;
    //     }
    // }

    int npts;//todo npts, uod, side 是什么？
    PointType vertex[2];
    DOUBLE L_EdgeLen, R_EdgeLen;
    DOUBLE dist, capac, subtree_cost;
    int side, uod;
    MsType *segment_L, *segment_R;
    TrrType L_Trr, R_Trr;
    PointType L_line[2], R_line[2];
};

class NodeType {
public:
    NodeType() {
        id = L = R = PAR = root_id = NIL;
    }
    PointType m_mergePt;
    MsType segment;
    TrrType trr;
    int id, L, R, PAR, root_id;//todo 这里为什么是int, par是什么
};

typedef struct buckettype {
    buckettype() {
        num = 0;
    }
    int num;    /* number of elements in the bucket */   
    int element[BUCKET_CAPACITY];
} BucketType;

class PairType {
public:
    int x, y, x_root_id, y_root_id;
    DOUBLE cost;
};

#include "utils.h"

#endif