#include "ZstTree.h"


/*********************************************************************/
/*  check merge segment                                              */
/*********************************************************************/
void ZstTree::check_MergeSegment(MsType *segment) {
    MsType *segment_L = segment->segment_L;
    MsType *segment_R = segment->segment_R;
    assert(segment_L->npts == 1 || segment_L->npts == 2);
    assert(segment_R->npts == 1 || segment_R->npts == 2);
}
/*********************************************************************/
/*  check updated MS                                                 */
/*********************************************************************/
void ZstTree::check_JS_MS(MsType segment) {
    TrrType ms0, ms1;

    line2ms(&ms0, segment.L_line[0], segment.L_line[1]);
    line2ms(&ms1, segment.R_line[0], segment.R_line[1]);

    assert(trrContain(&segment.L_Trr, &ms0));
    assert(trrContain(&segment.R_Trr, &ms1));
}
void ZstTree::check_update_MS(MsType *segment, PointType pt[4], int line0type) {
    int parallel;
    DOUBLE d;
    PointType tmp_pt[2],q0,q1,q2,q3;

    parallel = parallel_line(pt[0],pt[1],pt[2],pt[3]);

    if (parallel) {
        if (line0type == FLAT) { /* not consider non-rectilinear  parallel lines */
        assert(0);                                     
        } else if (line0type == TILT) {
        assert(0);                                     
        }
    }

    q0 = segment->L_line[0];
    q1 = segment->L_line[1];
    q2 = segment->R_line[0];
    q3 = segment->R_line[1];

    d=linedist(q0,q1,q2,q3,tmp_pt);
    if (CHECK == YES) {
        assert(same_Point(q0,q1) || same_Point(q2,q3) || 
            parallel_line(q0,q1,q2,q3));
        assert(d == segment->dist) ;
        assert(pt_on_line_segment(q0, pt[0], pt[1]));
        assert(pt_on_line_segment(q1, pt[0], pt[1]));
        assert(pt_on_line_segment(q2, pt[2], pt[3]));
        assert(pt_on_line_segment(q3, pt[2], pt[3]));
    }

}
/*********************************************************************/
/*  check compare neighbors                                          */
/*********************************************************************/
void ZstTree::check_compare_neighbors(int i, int j) {
    assert(i >= 0 && j >= 0);
    int k1 = Node[i].root_id;
    int k2 = Node[j].root_id;
    
    assert(k1 >= 0 && k2 >= 0 && k1 < m_nPoints &&  k2 < m_nPoints);
    assert(k1 != k2 && !Marked[i] && !Marked[j]);
}
/*********************************************************************/
/*  check root id                                                    */
/*********************************************************************/
void ZstTree::count_tree_nodes(int root, int v, int *n) {
    if (v < 0) return;
    (*n)++;
    assert(Node[v].root_id == root); 
    count_tree_nodes(root, Node[v].L, n);
    count_tree_nodes(root, Node[v].R, n);
}
void ZstTree::check_root_id(int root) {
    int n = 0;
    count_tree_nodes(root, root, &n);
}