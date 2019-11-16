#include "ZstTree.h"

void ZstTree::calc_merge_pt_delay(MsType *mergesegment, DOUBLE d0, DOUBLE d1, bool same_tier, int lor) {
    mergesegment->vertex[0].delay = calc_merge_pt_delay_sub(mergesegment, d0, d1, same_tier, lor, mergesegment->uod);
    mergesegment->vertex[1].delay = mergesegment->vertex[0].delay;
}
/*************************************************************************/
/*  calculate edge delay                                                 */
/*************************************************************************/
DOUBLE ZstTree::calc_delay_increase(DOUBLE cap, DOUBLE leng, DOUBLE pures, DOUBLE pucap) {
    DOUBLE t, t0, t1;

    /* Elmore delay model */
    t0 = pures * leng * (pucap * leng / 2 + cap);
    return(t0);
}
DOUBLE ZstTree::_pt_delay_increase(DOUBLE leng, DOUBLE cap, PointType *q0,PointType *q1, int lor, int uod, int i, bool same_tier) {
    DOUBLE t;

    DOUBLE h = (q0->x - q1->x)._ABS();
    DOUBLE v = (q0->y - q1->y)._ABS();
    assert((leng == h + v) || (leng >= h + v));
    if (lor == -1) {
        t = calc_delay_increase(cap, leng, m_pures, m_pucap);
        if (same_tier) {
            
        } else {
            assert(m_topoMode == MMMMODE);
            t += calc_delay_increase(cap + leng * m_pucap, 1, m_TsvR, m_TsvC);
        }
    } else if(lor == 0) {
        if (uod == -1) {
            assert(h + v == leng);
            if (i == 0) {
                t = calc_delay_increase(cap, leng, m_pures, m_pucap);
            } else {
                t = calc_delay_increase(cap, leng, m_pures, m_pucap);
                t = t + calc_delay_increase(cap + leng*m_pucap, 1, m_TsvR, m_TsvC);
            }
        } else if (uod / 2 == 0) {
            if (uod % 2 == 0) {
                if (i == 0) {
                    assert(leng == 0);
                    t = calc_delay_increase(cap, leng, m_pures, m_pucap);
                } else {
                    t = calc_delay_increase(cap, h + v, m_pures, m_pucap);
                    t = t + calc_delay_increase(cap + (h+v) * m_pucap, 1, m_TsvR, m_TsvC);
                    t = t + calc_delay_increase(cap + (h+v) * m_pucap + m_TsvC, leng-h-v, m_pures, m_pucap);
                }
            } else {
                if (i == 0) {
                    assert(leng == 0);
                    t = calc_delay_increase(cap, leng, m_pures, m_pucap);
                } else {
                    t = calc_delay_increase(cap, leng, m_pures, m_pucap);
                    t = t + calc_delay_increase(cap + leng * m_pucap, 1, m_TsvR, m_TsvC);
                }
            }
        } else {
            if (i == 0) {
                t = calc_delay_increase(cap, leng, m_pures, m_pucap);
            } else {
                assert(leng == 0);
                t = calc_delay_increase(cap, leng, m_pures, m_pucap);
                t = t + calc_delay_increase(cap + leng * m_pucap, 1, m_TsvR, m_TsvC);
            }
        }
    } else if (lor == 1) {
        if (uod == -1) {
            assert(h + v == leng);
            if (i == 0) {
                t = calc_delay_increase(cap, leng, m_pures, m_pucap);
                t = t + calc_delay_increase(cap + leng*m_pucap, 1, m_TsvR, m_TsvC);
            } else {
                t = calc_delay_increase(cap, leng, m_pures, m_pucap);
            }
        } else if (uod / 2 == 0) {
            if (i == 0) {
                assert(leng == 0);
                t = calc_delay_increase(cap, leng, m_pures, m_pucap);
                t = t + calc_delay_increase(cap + leng * m_pucap, 1, m_TsvR, m_TsvC);
            } else {
                t = calc_delay_increase(cap, leng, m_pures, m_pucap);
            }
        } else {
            if (uod % 2 == 0) {
                if (i == 0) {
                    t = calc_delay_increase(cap, leng, m_pures, m_pucap);
                    t = t + calc_delay_increase(cap + leng * m_pucap, 1, m_TsvR, m_TsvC);
                } else {
                    assert(leng == 0);
                    t = calc_delay_increase(cap, leng, m_pures, m_pucap);
                }
            } else {
                if (i == 0) {
                    t = calc_delay_increase(cap, h + v, m_pures, m_pucap);
                    t = t + calc_delay_increase(cap + (h+v) * m_pucap, 1, m_TsvR, m_TsvC);
                    t = t + calc_delay_increase(cap + (h+v) * m_pucap + m_TsvC, leng-h-v, m_pures, m_pucap);
                } else {
                    assert(leng == 0);
                    t = calc_delay_increase(cap, leng, m_pures, m_pucap);
                }
            }
        }
    } else { assert(0); }
    return(t);
}
DOUBLE ZstTree::calc_edge_delay(int i, int par, int lor, int uod) {
    
    int v;
    if (i == 0) {
        v = Node[par].L;
    } else {
        v = Node[par].R;
    }
    PointType &pt = Node[ par ].m_mergePt ;
    PointType &qt = Node[ v   ].m_mergePt ;
    DOUBLE cap = Node[v].segment.capac;
    bool same_tier = (Node[Node[par].L].segment.vertex[0].z == Node[Node[par].R].segment.vertex[0].z);
    
    DOUBLE t = _pt_delay_increase(EdgeLength[v], cap, &(qt), &(pt), lor, uod, i, same_tier);
    
    return(t);
}
/*********************************************************************/
/*  calculate ZST tree cost                                             */
/*********************************************************************/
void ZstTree::calc_TreeCost_sub(int par, int v, DOUBLE *Tcost, DOUBLE *Tdist, int *n_detour) {
    DOUBLE dist, length;

    if (v<0) return;
    PointType &pt = Node[par].m_mergePt ;
    PointType &qt = Node[v].m_mergePt ;
    dist = Point_dist(pt, qt );
    length = EdgeLength[v];
    assert(length == dist || dist < length);
    if (dist < length) {
        (*n_detour)++;
    }
    *Tcost += length;
    *Tdist += dist;
    calc_TreeCost_sub(v,Node[v].L, Tcost, Tdist, n_detour);
    calc_TreeCost_sub(v,Node[v].R, Tcost, Tdist, n_detour);
}
int ZstTree::calc_TreeCost(int v, DOUBLE *Tcost, DOUBLE *Tdist) {
    int n_detour;

    n_detour = 0; 
    *Tcost = *Tdist = 0;

    calc_TreeCost_sub(v, Node[v].L, Tcost, Tdist, &n_detour);
    calc_TreeCost_sub(v, Node[v].R, Tcost, Tdist, &n_detour);
    return(n_detour);
}
/*********************************************************************/
/*  calculate ZST delay                                              */
/*********************************************************************/
void ZstTree::calc_ZST_delay_sub(int v) {
 
    int L = Node[v].L;
    int R = Node[v].R;
    
    PointType &ptL = Node[L].m_mergePt;
    PointType &ptR = Node[R].m_mergePt;
    if (ptL.delay == NIL) {
        calc_ZST_delay_sub(L);
    }
    
    if (ptR.delay == NIL) {
        calc_ZST_delay_sub(R);
    }
    
    DOUBLE tL, tR;
    int side = Node[v].segment.side;
    int uod = Node[v].segment.uod;
    tL = calc_edge_delay(0, v, side, uod);
    tR = calc_edge_delay(1, v, side, uod);
    

    PointType &pt = Node[ v ].m_mergePt ;
    if ( (ptL.delay+tL) - (ptR.delay+tR) < 1.1 ) {
        // assert(ptL.delay+tL == ptR.delay+tR);
        pt.delay = ptL.delay + tL;
    } else {
        assert(0);
    }

    //print_ZST_delay_error(v, L, R, tL, tR);
}
void ZstTree::calc_ZST_delay(int v) {

    for (int i = m_nTerms; i < m_nPoints; ++i) {
        PointType &pt = Node[i].m_mergePt;
        pt.delay = NIL;
    }
    calc_ZST_delay_sub(v);
}
/*********************************************************************/
/*  calculate Elomore merge distance                                 */
/*********************************************************************/
DOUBLE ZstTree::calc_Elmore_merge_distance_sub(DOUBLE cap1,DOUBLE delay1, DOUBLE cap2, DOUBLE delay2, DOUBLE d, int lor) {

    DOUBLE r = m_pures, c = m_pucap;
    DOUBLE rt = m_TsvR, ct = m_TsvC;
    DOUBLE x, y, z;
    if (lor == -1) {  /* same tier*/
        y = (delay2 - delay1 + r*d*(cap2 + c*d/2));
        z = (r*(cap1 + cap2 + c*d));
    } else if (lor == 0) {
        y = (delay2 - delay1 + rt*(c*d + cap2 + ct/2) + r*d*(cap2 + c*d/2));
        z = ((r*(cap1 + cap2 + c*d)) + rt*c);
    } else if (lor == 1) {
        y = (delay2 - delay1 - rt*(ct/2 + cap1) + r*d*(cap2 + c*d/2));
        //y = (delay2 - delay1 + rt*(ct/2 - cap1) + r*d*(cap2 + c*d/2));
        z = ((r*(cap1 + cap2 + c*d)) + rt*c);
    } else { assert(0); }

    x = y / z;
    return x;
}
DOUBLE ZstTree::calc_Elmore_merge_distance_detour1(DOUBLE cap1,DOUBLE delay1, DOUBLE cap2, DOUBLE delay2, 
DOUBLE d, DOUBLE *d1, DOUBLE *d2, int lor, int uod) {

    DOUBLE r = m_pures, c = m_pucap;
    DOUBLE rt = m_TsvR, ct = m_TsvC;
    DOUBLE x, h;
    DOUBLE A_, B_, C_;
    if (lor == -1) {
        cap2 *= r;
        x = ((cap2*cap2 + r*c* (delay1-delay2)*2)._SQRT() - cap2) / (r*c);
    } else if (lor == 0) {
        // h = r*cap2 + rt*c;
        // x = ((h*h - r*c*(delay2 - delay1 + rt*ct/2 + rt*cap2)*2)._SQRT() - h) / (r*c);
        if (uod == 0) {
            A_ = r*c/2;
            B_ = r*cap2 + r*ct + r*c*d;
            C_ = delay2 - delay1 + r*c*d*d/2 + r*d*cap2 + rt*ct/2 + rt*c*d + rt*cap2;
            x = d + QuadraticEquation(A_, B_, C_);
        } else if (uod == 1) {
            A_ = r*c/2;
            B_ = r*c*d + r*cap2 + rt*c;
            C_ = r*c*d*d/2 + r*d*cap2 + rt*ct/2 + rt*c*d + rt*cap2 + delay2 - delay1;
            x = d + QuadraticEquation(A_, B_, C_);
        } else { assert(0); }
    } else if (lor == 1) {
        A_ = r*c/2;
        B_ = r*c*d + r*cap2;
        C_ = r*c*d*d/2 + r*d*cap2 - rt*ct/2 - rt*cap1 + delay2 - delay1;
        x = d + QuadraticEquation(A_, B_, C_);
        // h = r*cap2;
        // x = ((h*h - r*c*(delay2 - delay1 - rt*ct/2 - rt*cap1)*2)._SQRT() - h) / (r*c);
    } else { assert(0); }
    *d1 = 0;
    *d2 = x;
    return x;
}
DOUBLE ZstTree::calc_Elmore_merge_distance_detour2(DOUBLE cap1,DOUBLE delay1, DOUBLE cap2, DOUBLE delay2, 
DOUBLE d, DOUBLE *d1, DOUBLE *d2, int lor, int uod) {

    DOUBLE r = m_pures, c = m_pucap;
    DOUBLE rt = m_TsvR, ct = m_TsvC;
    DOUBLE x, h;
    DOUBLE A_, B_, C_;
    if (lor == -1) {
        cap1 *= r;
        x = ((cap1*cap1 + r*c* (delay2-delay1)*2)._SQRT() - cap1) / (r*c);
    } else if (lor == 0) {
        A_ = r*c/2;
        B_ = r*c*d + r*cap1;
        C_ = r*c*d*d/2 + r*d*cap1 - rt*ct/2 - rt*cap2 + delay1 - delay2;
        x = d + QuadraticEquation(A_, B_, C_);
        // h = r*cap1;
        // x = ((h*h - r*c*(delay1 - delay2 - rt*ct/2 - rt*cap2)*2)._SQRT() - h) / (r*c);
    } else if (lor == 1) {
        // h = r*cap1 + rt*c;
        // x = ((h*h - r*c*(delay1 - delay2 + (rt*ct/2 + rt)*cap1)*2)._SQRT() - h) / (r*c);
        if (uod == 0) {
            A_ = r*c/2;
            B_ = r*c*d + r*cap1 + rt*c;
            C_ = r*c*d*d/2 + r*d*cap1 + rt*ct/2 + rt*c*d + rt*cap1 + delay1 - delay2;
            x = d + QuadraticEquation(A_, B_, C_);
        } else if (uod == 1) {
            A_ = r*c/2;
            B_ = r*ct + r*c*d + r*cap1;
            C_ = r*c*d*d/2 + r*d*cap1 + rt*ct/2 + rt*c*d + rt*cap1 + delay1 - delay2;
            x = d + QuadraticEquation(A_, B_, C_);
        } else { assert(0); }
    } else { assert(0); }
    *d1 = x;
    *d2 = 0;
    return x;
}
DOUBLE ZstTree::calc_Elmore_merge_distance(DOUBLE cap1,DOUBLE delay1, DOUBLE cap2, DOUBLE delay2,
DOUBLE d,DOUBLE *d1,DOUBLE *d2, int lor, int *uod) {

    DOUBLE x = calc_Elmore_merge_distance_sub(cap1, delay1, cap2, delay2, d, lor);

    DOUBLE distance[2];
    DOUBLE ea[2], eb[2];
    if (x < 0) {
        for (int i = 0; i < 2; i++) {
            distance[i] = calc_Elmore_merge_distance_detour1(cap1, delay1, cap2, delay2, d, &ea[i], &eb[i], lor, i);
        }
        *uod = (distance[0] <= distance[1]) ? 0 : 1; 
        *d1 = ea[*uod];
        *d2 = eb[*uod];
        assert(distance[*uod] >= d);
        return distance[*uod];
    } else if (x > d) {
        for (int i = 0; i < 2; i++) {
            distance[i] = calc_Elmore_merge_distance_detour2(cap1, delay1, cap2, delay2, d, &ea[i], &eb[i], lor, i);
        }
        *uod = (distance[0] <= distance[1]) ? 0 : 1; 
        *uod += 2;
        *d1 = ea[*uod%2];
        *d2 = eb[*uod%2];
        assert(distance[*uod%2] >= d);
        return distance[*uod%2];
    } else {
        *d1 = x;
        *d2 = d-x;
        *uod = -1;
        return x;
    }
}  
/* calc_merge_distance */
int ZstTree::calc_merge_distance(DOUBLE cap1, DOUBLE delay1, DOUBLE cap2, DOUBLE delay2, 
DOUBLE d, DOUBLE *d1,DOUBLE *d2, bool same_tier, int *uod, int force_lor) {
                              
    if (same_tier) {
        calc_Elmore_merge_distance(cap1, delay1, cap2, delay2, d, d1, d2, -1, uod);
        return -1;
    } else {
        if (force_lor == 0 || force_lor == 1) {
            DOUBLE a = calc_Elmore_merge_distance(cap1, delay1, cap2, 
                                    delay2, d, d1, d2, force_lor, uod);
            return force_lor;                                            
        }

        int lor;
        DOUBLE da1, da2, db1, db2;
        int uod2[2];
        DOUBLE a = calc_Elmore_merge_distance(cap1, delay1, cap2, 
                                    delay2, d, &da1, &da2, 0, &uod2[0]);
        DOUBLE b = calc_Elmore_merge_distance(cap1, delay1, cap2, 
                                    delay2, d, &db1, &db2, 1, &uod2[1]); 
        assert(a >= 0 && b >= 0);                                
        if (a >= 0 && a <= d) {
            assert(a == da1 && d - a == da2);
            *d1 = da1;
            *d2 = da2;
            lor = 0;
        } else if (b >= 0 && b <= d) {
            assert(b == db1 && d - b == db2);
            *d1 = db1;
            *d2 = db2;
            lor = 1;
        } else if (a > b) {
            *d1 = db1;
            *d2 = db2;
            lor = 1;
        } else {
            *d1 = da1;
            *d2 = da2;
            lor = 0;
        }
        *uod = uod2[lor];
        return lor;
    }
}
DOUBLE ZstTree::calc_merge_pt_delay_sub(MsType *segment, DOUBLE d0, DOUBLE d1, bool same_tier, int lor, int uod) {
    DOUBLE t0, t1;
    
    DOUBLE r = m_pures, c = m_pucap;
    DOUBLE rt = m_TsvR, ct = m_TsvC;
    DOUBLE cap1 = segment->segment_L->capac;
    DOUBLE cap2 = segment->segment_R->capac;
    DOUBLE d = segment->dist, x = d0 + d1 - d;
    if (same_tier) {
        assert(lor == -1);
        t0 = r*d0*(c*d0/2+cap1) + segment->L_line[0].delay;
        t1 = r*d1*(c*d1/2+cap2) + segment->R_line[0].delay;
    } else {
        assert(lor == 0 || lor == 1);
        if (lor == 0) {
            if (uod == -1) {
                t0=r*d0*(c*d0/2+cap1)+ segment->L_line[0].delay;
                t1=r*d1*(c*d1/2+cap2) + m_TsvR*(m_TsvC/2+c*d1+cap2)+segment->R_line[0].delay;
            } else if(uod / 2 == 0) {
                if (uod % 2 == 0) {
                    t0 = segment->L_line[0].delay;
                    t1 = r*d*(c*d/2 + cap2) + rt*(ct/2 + c*d + cap2) + r*x*(c*x/2 + cap2 + ct + c*d) + segment->R_line[0].delay;
                } else {
                    t0 = segment->L_line[0].delay;
                    t1 = r * (d+x) * (c*(d+x)/2+cap2) + rt*(ct/2 + c*(d+x) + cap2) + segment->R_line[0].delay;
                }
            } else {
                t0 = r * (d+x) * (c*(d+x)/2+cap1) + segment->L_line[0].delay;
                t1 = rt*(ct/2 + cap2) + segment->R_line[0].delay;
            }
        } else if (lor == 1) {
            if (uod == -1) {
                t0=r*d0*(c*d0/2+cap1) + m_TsvR*(m_TsvC/2+c*d0+cap1)+segment->L_line[0].delay;
                t1=r*d1*(c*d1/2+cap2)+segment->R_line[0].delay;
            } else if (uod / 2 == 0) {
                t0 = rt*(ct/2 + cap1) + segment->L_line[0].delay;
                t1 = r * (d+x) * (c*(d+x)/2 + cap2) + segment->R_line[0].delay;
            } else {
                if (uod % 2 == 0) {
                    t0 = r * (d+x) * (c*(d+x)/2 + cap1) + rt*(ct/2 + c*(d+x) + cap1) + segment->L_line[0].delay;
                    t1 = segment->R_line[0].delay;
                } else {
                    t0 = r*d*(c*d/2 + cap1) + rt*(ct/2 + c*d + cap1) + r*x*(c*x/2 + ct + c*d + cap1) + segment->L_line[0].delay;
                    t1 = segment->R_line[0].delay;
                }
            }
        }
    }
    
    if (t0 != t1) {
        printf("lor: %d  uod: %d  same_tier: %d\n", lor, uod, same_tier);
    }
    assert(t0 == t1);
    return(t0);
}
int ZstTree::calc_segment_EdgeLen(MsType *segment, DOUBLE *d0, DOUBLE *d1, int force_lor) {

    DOUBLE cap0 = segment->segment_L->capac;
    DOUBLE cap1 = segment->segment_R->capac;
    int z0 = segment->segment_L->vertex[0].z;
    int z1 = segment->segment_R->vertex[0].z;
    bool same_tier = (z0 == z1);

    // assert(segment->L_line[0].delay == segment->L_line[1].delay);
    // assert(segment->R_line[0].delay == segment->R_line[1].delay);
    DOUBLE tL = segment->L_line[0].delay;
    DOUBLE tR = segment->R_line[0].delay;

    int uod;
    int lor = calc_merge_distance(cap0, tL, cap1, tR, segment->dist, d0, d1, same_tier, &uod, force_lor);
    segment->side = lor;
    segment->uod = uod;
    /* just for check */
    calc_merge_pt_delay_sub(segment, *d0, *d1, same_tier, lor, uod);
    
    segment->L_EdgeLen = *d0;
    segment->R_EdgeLen = *d1;

    return lor;
}
DOUBLE ZstTree::calc_merging_cost(NodeType *node_L, NodeType *node_R) {
    MsType *segment_L = &(node_L->segment);
    MsType *segment_R = &(node_R->segment);
    MsType segment;

    calc_MS(&segment, segment_L, segment_R, -1, -1);
    MergeSegment(&segment, segment_L, segment_R, -1); 
    return(segment.subtree_cost);
}