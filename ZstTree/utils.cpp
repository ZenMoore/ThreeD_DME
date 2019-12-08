#include "utils.h"

/* types of interseciton of two line segments */
#define P1             1//交于端点1
#define P2             2//交于端点2
#define P3             3//交于端点3
#define P4             4//交于端点4
#define XING           5//交叉
#define OVERLAPPING    6//重合
#define SAME_EDGE      7//完全重合


void TRR2pt(TrrType *trr, PointType *pt) {
    pt->x = (trr->ylow + trr->xhi )/2;
    pt->y = (trr->ylow - trr->xhi )/2;
}
bool empty_trr(TrrType *t) {
    t->SelfDOCHECK();
    if (t->xlow > t->xhi || t->ylow > t->yhi) {
        return(YES);
    } else {
        return(NO);
    }
}
void check_trr(TrrType *t) {
    empty_trr(t);
    assert(t->xlow <= t->xhi);
    assert(t->ylow <= t->yhi);
}
int TRR2segment(TrrType *trr, PointType pt[]) {
    TrrType ms;
    int n;

    DOUBLE a = trr->xhi - trr->xlow;
    DOUBLE b = trr->yhi - trr->ylow;

    check_trr(trr);

    if (a == 0 && b == 0 ) {
        TRR2pt(trr, &(pt[0]));
        n = 1;
    } else if (a == 0 || b == 0 ) {
        ms2line(trr, &(pt[0]), &(pt[1]));//todo 我的理解是：这是需要detour的情况，我们通过将线段延长后，将线连到延长段上完成detour
        n = 2;
    } else {
        assert(0);
    }
    return(n);
}
bool Manhattan_arc(PointType p1, PointType p2) {
    DOUBLE a = (p2.y-p1.y)._ABS();
    DOUBLE b = (p2.x-p1.x)._ABS();
    return (a == b) ? YES : NO;
}
bool same_Point(PointType p1, PointType p2) {
    DOUBLE dist;
    dist = (p1.x - p2.x)._ABS() + (p1.y - p2.y)._ABS();
    return dist == 0;
}
bool parallel_line(PointType p1,PointType p2,PointType p3,PointType p4) {
    DOUBLE dx1, dx2, dy1, dy2;

    if (same_Point(p1,p2) || same_Point(p3,p4)) return(NO);

    dx1 = p1.x-p2.x; 
    dy1 = p1.y-p2.y;
    dx2 = p3.x-p4.x; 
    dy2 = p3.y-p4.y;

    if ( dx1 == 0 && dx2 == 0) { /* parallel vertical lines */
        return(YES);
    } else if ( dx1 != 0 && dx2 != 0 && dy1/dx1 == dy2/dx2) {
        return(YES);  /* other or horizontal parallel lines */
    }
    return(NO);
}
void build_trr(TrrType *ms, DOUBLE d, TrrType *trr) {//todo TrrType能是斜矩形吗？
    trr->xlow = ms->xlow - d;
    trr->xhi  = ms->xhi + d;

    trr->ylow = ms->ylow - d;
    trr->yhi  = ms->yhi + d;

    trr->z = ms->z;

} /* build_trr */
void build_NodeTRR(NodeType *node) {//todo 解释一下这个trr生成算法
    int j;
    DOUBLE  x, y;
    unsigned z;
 
    node->trr.MakeDiamond(node->segment.vertex[0], 0);
    assert(node->segment.npts == 1 || node->segment.npts == 2);
    if (node->segment.npts == 2) {
        x = node->segment.vertex[1].x;
        y = node->segment.vertex[1].y;
        z = node->segment.vertex[1].z;
        node->trr.xlow = tMIN(node->trr.xlow, x-y);
        node->trr.xhi  = tMAX(node->trr.xhi , x-y);
        node->trr.ylow = tMIN(node->trr.ylow, x+y);
        node->trr.yhi  = tMAX(node->trr.yhi , x+y);
        assert(node->trr.z == z);
    }
}

DOUBLE QuadraticEquation(DOUBLE a, DOUBLE b, DOUBLE c) {
    DOUBLE derta = b*b - a*c*4;
    assert(derta >= 0);
    return (derta._SQRT() - b) / (a * 2);
}

DOUBLE Point_dist(PointType p1, PointType p2) {
    return((p1.x - p2.x)._ABS() + (p1.y - p2.y)._ABS());
}
bool trrContain(TrrType *t1, TrrType *t2) {

    if ( t1->xhi >= t2->xhi - FUZZ && t1->xlow <= t2->xlow + FUZZ  &&
        t1->yhi >= t2->yhi - FUZZ && t1->ylow <= t2->ylow + FUZZ) {
        return(YES);
    } else { 
        return(NO); 
    }
}
/********************************************************************/
/*  merging segment to line segment                                 */
/********************************************************************/
bool a_segment_TRR(TrrType *trr) {
    /* check if *trr is a valid merging segment */
    if (trr->xlow == trr->xhi || trr->ylow == trr->yhi) {
        return(YES);
    } else {
        return(NO);
    }
}
void check_ms(TrrType *ms) {
    /* *ms must be a valid trr */
    check_trr(ms);

    /* check if *ms is a valid merging segment */
    assert(a_segment_TRR(ms));

    /* remove the epsilon error */
    ms->SelfDOCHECK();
}
void ms_to_line(TrrType *ms, DOUBLE *x1, DOUBLE *y1, int *z1, DOUBLE *y2, DOUBLE *x2, int *z2) {//todo 这个怎么就变成了line了？还是和make diamond一样的问题
    /* *ms must be a valide merging segement */
    check_ms(ms);
    *x1 = (ms->ylow + ms->xhi ) / 2;
    *y1 = (ms->ylow - ms->xhi ) / 2;
    *z1 = ms->z;
    *x2 = (ms->yhi  + ms->xlow) / 2;
    *y2 = (ms->yhi  - ms->xlow) / 2;
    *z2 = ms->z;
    /* p2 must be higher than p1 */
    assert(*y1 <= *y2);
}
void ms2line(TrrType *ms, PointType *p1, PointType *p2) {
    ms_to_line(ms, &(p1->x), &(p1->y), &(p1->z), &(p2->x), &(p2->y), &(p2->z)); 
}
/****************************************************************************/
/*  calculate line type                                                     */
/****************************************************************************/
int calc_line_type(PointType pt1, PointType pt2) {
    DOUBLE a = (pt1.x-pt2.x)._ABS();
    DOUBLE b = (pt1.y-pt2.y)._ABS();

    if (a == b) {
        return(MANHATTAN_ARC);
    } else if (a == 0) {
        return(VERTICAL);
    } else if (b == 0) {
        return(HORIZONTAL);
    } else if (a > b) {
        return(FLAT);
    } else if (a < b) {
        return(TILT);
    }
    assert ( 0 ) ;
}

/****************************************************************************/
/*  check if two bounding boxes intesect.                                   */
/****************************************************************************/
int bbox_overlap(DOUBLE x1, DOUBLE y1, DOUBLE x2, DOUBLE y2,
                 DOUBLE x3, DOUBLE y3, DOUBLE x4, DOUBLE y4) {
    if (tMIN(x1,x2) >= tMAX(x3,x4) + FUZZ ) return(NO);
    if (tMIN(x3,x4) >= tMAX(x1,x2) + FUZZ ) return(NO);
    if (tMIN(y1,y2) >= tMAX(y3,y4) + FUZZ ) return(NO);
    if (tMIN(y3,y4) >= tMAX(y1,y2) + FUZZ ) return(NO);
    return(YES);
}

/****************************************************************************/
/* check if (x,y) is in line_segment p1,p2                                  */
/****************************************************************************/
int ON_line_segment(DOUBLE *x,DOUBLE *y,DOUBLE x1,DOUBLE y1,DOUBLE x2,DOUBLE y2)
{
    DOUBLE new_x, new_y;
    DOUBLE t = DBL_MAX;
    int ans = NO;

    DOUBLE d = (x1 - x2)._ABS() + (y1 - y2)._ABS();      
    DOUBLE d1 = (x1 - *x)._ABS() + (y1 - *y)._ABS();  
    DOUBLE d2 = (x2 - *x)._ABS() + (y2 - *y)._ABS();

    //cout << "不知道要不要变x和y" << endl;
    //cout << "不知道要不要删前两个if" << endl;
    if ( (d1 + d2 - d)._ABS() <= 2 * FUZZ ) {
        if ( d1 == 0 ) {
            *x = x1;     
            *y = y1;      
            ans = YES;
        } else if ( d2 == 0 ) {
            *x = x2;     
            *y = y2;      
            ans = YES;
        } else {
            new_x = *x; 
            new_y = (y2-y1)*(*x-x1)/(x2-x1) + y1;
            t = (new_y - *y)._ABS();  
            if (t < FUZZ) {
                ans = YES;
                *x = new_x;
                *y = new_y;
            }
        }
    }
    return(ans);
}
int on_line_segment(DOUBLE x,DOUBLE y,DOUBLE x1,DOUBLE y1,DOUBLE x2,DOUBLE y2)
{
    return(ON_line_segment(&x,&y,x1,y1,x2,y2));
}

/***********************************************************************/
/* check if (x,y) is in the box defined by max_x,min_x,max_y,min_y     */
/***********************************************************************/
int in_bbox(DOUBLE x, DOUBLE y, DOUBLE x1,DOUBLE y1, DOUBLE x2,DOUBLE y2) {
    if (x <= tMAX(x1,x2) + FUZZ && x >= tMIN(x1,x2) - FUZZ && y <= tMAX(y1,y2)
        + FUZZ && y >= tMIN(y1,y2) - FUZZ )
    { return(YES); }
    return(NO);
}

/****************************************************************************/
/*  check if two line segments intersect.                                   */
/****************************************************************************/
int L_intersect(DOUBLE *x, DOUBLE *y, DOUBLE x1, DOUBLE y1,DOUBLE x2, 
DOUBLE y2, DOUBLE x3, DOUBLE y3, DOUBLE x4, DOUBLE y4) {

    assert(x1 != x2 && x3 != x4);
    DOUBLE mm=(y2-y1)/(x2-x1);
    DOUBLE nn=(y4-y3)/(x4-x3);
    assert(mm == 1 || mm == -1);
    assert(nn == 1 || nn == -1);

    /* zero-length_edge */ 
    if ((x1-x2)._ABS() + (y1-y2)._ABS() < FUZZ || (x3-x4)._ABS() + (y3-y4)._ABS() < FUZZ) 
        assert(0);             
    
    if (!bbox_overlap(x1,y1,x2,y2,x3,y3,x4,y4)) return(NO);

    //cout << "确认x1和x3的对应关系" << endl;
    if ((x1-x3)._ABS() + (y1-y3)._ABS() + (x2-x4)._ABS() + (y2-y4)._ABS() < FUZZ)  return(SAME_EDGE);

    int count = 0;
    *x = DBL_MAX; *y = DBL_MAX;
    if (on_line_segment(x1,y1,x3,y3,x4,y4)) {
        *x = x1;
        *y = y1;
        count++;
    }
    if (on_line_segment(x2,y2,x3,y3,x4,y4)) {
        *x = x2;
        *y = y2;
        count++;
    }
    if (on_line_segment(x3,y3,x1,y1,x2,y2)) {
        *x = x3;
        *y = y3;
        count++;
    }
    if (on_line_segment(x4,y4,x1,y1,x2,y2)) {
        *x = x4;
        *y = y4;
        count++;
    }
    if ((x1-x3)._ABS() + (y1-y3)._ABS() < FUZZ) count--;        
    if ((x1-x4)._ABS() + (y1-y4)._ABS() < FUZZ) count--;        
    if ((x2-x3)._ABS() + (y2-y3)._ABS() < FUZZ) count--;        
    if ((x2-x4)._ABS() + (y2-y4)._ABS() < FUZZ) count--;        
    if (count >= 2) return(OVERLAPPING);


    if (mm == nn) {    
        return(NO);         /* disjoint parallel line segments */
    } else {               /* Non-parallel line segments */
        *x = ((y3-nn*x3) - (y1-mm*x1)) / (mm-nn);
        *y = ((y1 + mm*(*x - x1)) + (y3 + nn*(*x - x3))) / 2.0;
    }//求出两线交点

    if (in_bbox(*x,*y,x1,y1,x2,y2) && in_bbox(*x,*y,x3,y3,x4,y4) ) {

        if ( *x == x1 && *y == y1 ) {
            return(P1);
        } else if ( *x == x2 && *y == y2 ) {
            return(P2);
        } else if ( *x == x3 && *y == y3 ) {
            return(P3);
        } else if ( *x == x4 && *y == y4 ) {
            return(P4);
        } else {
            return(XING);
        }
    }
    return(NO);
}

int lineIntersect(PointType *p, PointType p1, PointType p2, PointType p3, PointType p4) {
    DOUBLE x, y; 
    int ans; 

    ans = L_intersect(&x, &y, p1.x, p1.y, p2.x, p2.y, p3.x, p3.y, p4.x, p4.y);
    p->x = x;
    p->y = y;
    return ans;
}

/****************************************************************************/
/* check if point pt is in line_segment p1, p2                           */
/****************************************************************************/
int pt_on_line_segment(PointType pt,PointType pt1,PointType pt2)
{
  return(ON_line_segment(&(pt.x),&(pt.y),pt1.x,pt1.y,pt2.x,pt2.y));
}
/****************************************************************************/
/* check if (x,y) is in line_segment p1,p2                                  */
/* if yes, eliminate the epsilon error. */
/****************************************************************************/
int PT_on_line_segment(PointType *pt,PointType pt1,PointType pt2) {
    return ON_line_segment(&(pt->x), &(pt->y), pt1.x, pt1.y, pt2.x, pt2.y);
}
/****************************************************************************/
/* return the shortest distance between point p1 and line p2-p3 */
/****************************************************************************/
/********************************************************************/
/* check if *ms is a valid merging segment */
/********************************************************************/
void line_to_ms(TrrType *ms, DOUBLE x1, DOUBLE y1, DOUBLE x2, DOUBLE y2, int z) {
    DOUBLE a,b;

    if (y1 > y2) { /* p2 must be higher than p1 */
        a = x1; 
        b = y1; 
        x1 = x2; 
        y1 = y2; 
        x2 = a; 
        y2 = b; 
    }
    ms->ylow = (x1+y1); 
    ms->yhi  = (x2+y2);
    ms->xlow = (x2-y2); 
    ms->xhi  = (x1-y1);
    ms->z = z;
    /* *ms must be a Manhattan arc. */
    check_ms(ms); 
}
void line2ms(TrrType *ms, PointType p1, PointType p2) {
    assert(p1.z == p2.z);
    line_to_ms(ms, p1.x, p1.y, p2.x, p2.y, p1.z);  
}
DOUBLE ms_distance(TrrType *ms1, TrrType *ms2) {
    DOUBLE x1a,x1b,y1a,y1b,x2a,x2b,y2a,y2b;
    DOUBLE d1,d2,d3,d4;
    DOUBLE t1,t2,t3;

    /* remember that ms1 and ms2 are represented in the Linfinity metric */
    /* first check that these are valid merging segments */
    if (CHECK == YES) check_ms(ms1);
    if (CHECK == YES) check_ms(ms2);

    x1a = ms1->xhi; y1a = ms1->yhi;
    x1b = ms1->xlow; y1b = ms1->ylow;

    x2a = ms2->xhi; y2a = ms2->yhi;
    x2b = ms2->xlow; y2b = ms2->ylow;

    /*there are four cases to consider:
        (1) there is no intersetion between x-coords of ms1 and x-coords of ms2
            and no intersetion tween y-coords of ms1 and y-coords of ms2
        (2) no intersection tween x-coords, but some intersection in y-coords
        (3) no intersection tween y-coords, but some intersection in x-coords
        (4) intersection in both coords (distance=0)

        (1) find min distance between endpoints
        (2) find min distance in x-coords
        (3) find min distance in y-coords
        (4) return distance = 0
        */

        /* if no intersection between x & y coordinates at all, take
        the min distance between endpoints */
    if ( ((x1a < x2b) || (x2a < x1b))  && ((y1a < y2b) || (y2a < y1b)) ) {

        /* find the distance between all 4 pairs of endpoints */
        d1 = tMAX((x1a-x2a)._ABS(), (y1a-y2a)._ABS());
        d2 = tMAX((x1b-x2a)._ABS(), (y1b-y2a)._ABS());
        d3 = tMAX((x1a-x2b)._ABS(), (y1a-y2b)._ABS());
        d4 = tMAX((x1b-x2b)._ABS(), (y1b-y2b)._ABS());

        t1 =  tMIN(d1,d2);
        t2 =  tMIN(d3,d4);
        t3 =  tMIN(t1,t2);
        return( t3 );
    }
    else if ((x1a < x2b) || (x2a < x1b)) {
        t1 = x2b - x1a;
        t2 = x1b - x2a;
        t3 =  tMAX(t1,t2);  /* take max here because one will be negative */
        return( t3 );
    }
    else if ((y1a < y2b) || (y2a < y1b)) {
        t1 = y2b - y1a;
        t2 = y1b - y2a;
        t3 =  tMAX(t1,t2); /* take max here because one will be negative */
        return( t3 );
    }
    else return( 0 );
  
} /* ms_distance() */
void make_intersect_sub( TrrType *trr1, TrrType *trr2, TrrType *t ) {

    t->xlow = tMAX(trr1->xlow,trr2->xlow);
    t->xhi  = tMIN(trr1->xhi, trr2->xhi );

    t->ylow = tMAX(trr1->ylow,trr2->ylow);
    t->yhi  = tMIN(trr1->yhi, trr2->yhi );

    if (trr1->z == trr2->z) {
        t->z = trr1->z;
    } else {
        t->z = 0;
    }
}
/********************************************************************/
/*                                                                  */
/********************************************************************/
void make_intersect( TrrType *trr1, TrrType *trr2, TrrType *t ) {

    make_intersect_sub(trr1, trr2, t);
    /* check if intersection is non-empty */
    check_trr(t); 

} /* make_intersect */
/********************************************************************/
/* find the mid-point of the core of trr */
/********************************************************************/
void core_mid_point(TrrType *trr, PointType *p) {
    DOUBLE tx,ty;

    tx = (trr->xlow + trr->xhi) / 2;
    ty = (trr->ylow + trr->yhi) / 2;

    p->x = (tx+ty)/2;
    p->y = (ty-tx)/2;
    p->z = trr->z;

} /* core_mid_point */
DOUBLE pt2linedist_case_ms(PointType p1, PointType p2, PointType p3, PointType *ans) {
    TrrType ms0, ms1, trr0, ms2;
    DOUBLE dist;

    ms0.MakeDiamond(p1, 0);
    line2ms(&ms1, p2, p3);
    dist = ms_distance(&ms0, &ms1);
    trr0.MakeDiamond(p1, dist);
    make_intersect(&trr0,&ms1,&ms2);
    core_mid_point(&ms2, ans);
    assert(p2.z == p3.z);
    ans->z = p2.z;

    return(dist);
}
DOUBLE pt2linedist(PointType p1, PointType p2, PointType p3, PointType *ans) {
    int i,n;
    PointType candidate_pt[4];
    DOUBLE d, min_d, a, b;

    a = (p3.x - p2.x)._ABS(); 
    b = (p3.y - p2.y)._ABS(); 
    if (same_Point(p2,p3)) {
        *ans = p2;
        min_d = Point_dist(p1,p2);
    } else if (pt_on_line_segment(p1,p2,p3)) {
        *ans = p1;
        min_d = 0;
    } else {  /* p1, p2 is Manhattan arc */
        assert(a == b);
        min_d = pt2linedist_case_ms(p1,p2,p3,ans);
    }
    assert(pt_on_line_segment(*ans,p2,p3));
    return(min_d);
}

DOUBLE linedist(PointType lpt0, PointType lpt1, PointType lpt2, PointType lpt3, PointType ans[2]) {
    DOUBLE dist, mindist = DBL_MAX;
    PointType intersect, pt[2][2];
    int k; 

    pt[0][0] = lpt0;
    pt[0][1] = lpt1;
    pt[1][0] = lpt2;
    pt[1][1] = lpt3;

    int n0 = same_Point(lpt0, lpt1) ? 1 : 2;
    int n1 = same_Point(lpt2, lpt3) ? 1 : 2;

    if (n0 == 2 && n1 == 2) {
        k = lineIntersect(&intersect, lpt0, lpt1, lpt2, lpt3);
        if (k > 0) {
            assert(lpt0.z == lpt1.z);
            assert(lpt2.z == lpt3.z);
            ans[0] = ans[1] = intersect; 
            ans[0].z = lpt0.z;
            ans[1].z = lpt2.z;
            return(0.0);
        }
    }

    for (int i = 0; i < n0; ++i) {
        k = (i+1)%2; 
        for (int j = 0; j < n1; ++j) {
            dist = pt2linedist(pt[i][j], pt[k][0], pt[k][1], &intersect);
            if (dist < mindist ) {
                mindist = dist;
                ans[i] = pt[i][j]; 
                ans[k] = intersect;
            }
        }
    }
    return(mindist);
}