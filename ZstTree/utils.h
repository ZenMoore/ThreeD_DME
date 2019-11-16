#ifndef _UTILS_H
#define _UTILS_H

#include "header.h"

/* types of a line segment */ 
#define VERTICAL          0
#define HORIZONTAL        1
#define MANHATTAN_ARC     2           /* include single points */
#define FLAT              3
#define TILT              4

// Trr
void build_NodeTRR(NodeType *node);
void build_trr(TrrType *ms, DOUBLE d, TrrType *trr);
bool trrContain(TrrType *t1, TrrType *t2);
int TRR2segment(TrrType *trr, PointType pt[]);
// Line
DOUBLE linedist(PointType lpt0, PointType lpt1, PointType lpt2, PointType lpt3, PointType ans[2]);
int calc_line_type(PointType pt1, PointType pt2);
void line2ms(TrrType *ms, PointType p1, PointType p2);
DOUBLE ms_distance(TrrType *ms1, TrrType *ms2);
void ms2line(TrrType *ms, PointType *p1, PointType *p2);
bool parallel_line(PointType p1,PointType p2,PointType p3,PointType p4);
bool Manhattan_arc(PointType p1, PointType p2);
bool a_segment_TRR(TrrType *trr);
// Point
DOUBLE Point_dist(PointType p1, PointType p2);
bool same_Point(PointType p1, PointType p2);
int pt_on_line_segment(PointType pt,PointType pt1,PointType pt2);
int PT_on_line_segment(PointType *pt,PointType pt1,PointType pt2);

DOUBLE QuadraticEquation(DOUBLE a, DOUBLE b, DOUBLE c);
DOUBLE pt2linedist(PointType p1, PointType p2, PointType p3, PointType *ans);

#endif