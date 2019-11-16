#include "ZstTree.h"




/*********************************************************************/
/*  build nearest neighbor graph                                     */
/*********************************************************************/
int ZstTree::xlow_index(NodeType node) {
    DOUBLE t = (node.trr.xlow - m_ctrr.xlow) * N_Index / m_ctrr.Width(X);
    int i = tMIN( (int) t.value, N_Index - 1); 
    assert(i >= 0 && i < N_Index);
    return i;
}
int ZstTree::xhi_index(NodeType node) {
    DOUBLE t = (node.trr.xhi - m_ctrr.xlow) * N_Index / m_ctrr.Width(X);
    int i = tMIN( (int) t.value, N_Index - 1); 
    assert (i >= 0 && i < N_Index);
    return i; 
}
int ZstTree::ylow_index(NodeType node) {
    DOUBLE t = (node.trr.ylow - m_ctrr.ylow) * N_Index / m_ctrr.Width(Y) ;         
    int i = tMIN( (int) t.value, N_Index - 1); 
    assert (i >= 0 && i < N_Index);
    return i; 
}
int ZstTree::yhi_index(NodeType node) {
    DOUBLE t = (node.trr.yhi - m_ctrr.ylow) * N_Index / m_ctrr.Width(Y) ;         
    int i = tMIN( (int) t.value, N_Index - 1); 
    assert (i >= 0 && i < N_Index);
    return i; 
}
void ZstTree::BuildMainTrr(int n) {
    m_ctrr = TrrType();
    int j = 0;
    for (int i = 0; i < n; i++) {
        if (!Marked[i]) {
            m_ctrr.Enclose(Node[i].trr);
            j++;
        }
    }
    double t = sqrt((double) j) ;
    N_Index = tMAX(1, (int) t);

    assert(N_Index <= MAX_N_Index);
}
void ZstTree::bucket_partitioning_sub(int v) {
    int j1 = xlow_index( Node[v] );
    int j2 = xhi_index( Node[v] );
    int k1 = ylow_index( Node[v] );
    int k2 = yhi_index( Node[v] );

    int n;
    for (int j = j1; j <= j2; j++) {
        for (int k = k1; k <= k2; k++) {
            n = Bucket[j][k].num;
            Bucket[j][k].element[n] = v;
            Bucket[j][k].num = n + 1;
            if (Bucket[j][k].num > BUCKET_CAPACITY) {
                printf("Bucket[%d][%d].num = %d > %d \n", 
                    j, k, Bucket[j][k].num, BUCKET_CAPACITY);
            }
            assert(Bucket[j][k].num <= BUCKET_CAPACITY);
        }
    }
}
void ZstTree::bucket_partitioning(int n) {
    BuildMainTrr(n);
    for (int i = 0; i < N_Index; i++) {
        for (int j = 0; j < N_Index; j++) { 
            Bucket[i][j].num = 0;
        } 
    }
    for (int i = 0; i < n; i++) {
        if (!Marked[i]) {
            bucket_partitioning_sub(i);
        }
    }
    if (SHOW_INFO) print_Bucket();
}
void ZstTree::update_neighbors(int x, int y, DOUBLE cost) {

    int n = N_neighbors[x];
    int j;
    for (j = 0; j < n; ++j) {
        assert(Node[x].root_id != Node[The_NEIghors[x][j]].root_id);
        if ( Neighbor_Cost[x][j] > cost ) {
            break;
        }
    }
    if (j>=N_Neighbor) return; 
    N_neighbors[x]= n = tMIN(n+1, N_Neighbor);
    
    for (int i = n-2; i >= j;--i) {
        assert(0);
        The_NEIghors[x][i+1] = The_NEIghors[x][i];
        Neighbor_Cost[x][i+1] = Neighbor_Cost[x][i];
    }
    The_NEIghors[x][j] = y;
    Neighbor_Cost[x][j] = cost;
    //cout << x << " " << y << " " << cost.value << endl;
}
void ZstTree::compare_neighbors(int i, int j) {

    if (CHECK == YES) check_compare_neighbors(i, j); 
    
    PointType &pt = Node[i].m_mergePt ;
    pt.t += 1.0;  /* number of comparisons for mr[i]; */ 

    DOUBLE cost = calc_merging_cost(&(Node[i]), &(Node[j]));
    int k1 = Node[i].root_id;
    int k2 = Node[j].root_id;
    assert(k1 == i && k2 == j);
    DOUBLE old_cost = Node[k1].segment.subtree_cost + Node[k2].segment.subtree_cost;

    //DOUBLE cost_inc = cost - old_cost;
    DOUBLE cost_inc = cost;
    update_neighbors(i, j, cost_inc);
}
void ZstTree::do_compare_neighbors(int i,int j) {
    int k1 = Node[i].root_id;
    int k2 = Node[j].root_id;
    if (k1 != k2)  {
        compare_neighbors(i, j);
    }
}
void ZstTree::compare_neighbors_in_bucket(int v, int x, int y) {
    assert(v != NIL);
    assert(x >= 0 && y >= 0 && x < N_Index && y < N_Index);
    for (int j = 0; j < Bucket[x][y].num; ++j) {
        do_compare_neighbors(v, Bucket[x][y].element[j]);
    }
}
void ZstTree::init_nearest_neighbor(int v) {
    PointType &pt = Node[v].m_mergePt ;
    pt.t = 0.0;   /* number of comparisons */ 
    N_neighbors[v] = 0;

    int x1 = xlow_index( Node[v]);
    int x2 = xhi_index( Node[v]);
    int y1 = ylow_index( Node[v]);
    int y2 = yhi_index( Node[v]);
    for (int i = x1; i <= x2; ++i) {
        for (int j = y1; j <= y2; ++j) {
            compare_neighbors_in_bucket(v, i, j);
        }
    }
}
void ZstTree::calc_nearest_neighbor_sub(int v,int out_xlow,int out_xhi,int out_ylow, int out_yhi) {

    int n1 = tMAX(0, out_ylow);
    int n2 = tMIN(N_Index-1, out_yhi);
    for (int j = n1; j <= n2; ++j) {
        if (out_xlow >= 0)     compare_neighbors_in_bucket(v,out_xlow,j);
        if (out_xhi < N_Index) compare_neighbors_in_bucket(v,out_xhi,j);
    }
    n1 = tMAX(0,out_xlow+1);
    n2 = tMIN(N_Index-1,out_xhi-1);
    for (int j=n1;j<=n2;++j) {
        if (out_ylow>=0)       compare_neighbors_in_bucket(v,j, out_ylow); 
        if (out_yhi < N_Index) compare_neighbors_in_bucket(v,j, out_yhi);
    }
}
int ZstTree::calc_nearest_neighbor(int v, int inc) {
    int  x1,x2,y1,y2;

    x1 = xlow_index( Node[v] ) - inc;
    x2 = xhi_index( Node[v] )  + inc;
    y1 = ylow_index( Node[v] ) - inc;
    y2 = yhi_index( Node[v] )  + inc;

    if (x1>=0 || y1>=0 || x2 <N_Index || y2 < N_Index) {
        calc_nearest_neighbor_sub(v,x1,x2,y1,y2);
        return(YES);
    } else {
        return(NO);
    }
    /* assert(x1>=0 || y1>=0 || x2 <N_Index || y2 < N_Index); */
}
int ZstTree::construct_NNG(int n) {
    int graphSize = 0; // control search range for nearest neighbors

    assert(n > 1);
    int inc = 0, k1, j;
    do {
        ++inc;
        k1 = 0;
        for (int i = 0; i < n; ++i) {
            j = UnMarkedNodes[i];
            calc_nearest_neighbor(j, inc);
            if (N_neighbors[j] > 0) { k1++; }
        }
        if (CHECK) assert(inc <= N_Index);
    } while (k1  == 0 || inc < graphSize);
    assert(k1 > 0);
    if (SHOW_INFO) print_NNG(n);
    return(k1);
}
void ZstTree::build_nearest_neighbor_graph(int n, int show_info) {

    bucket_partitioning(n);
    for (int i = 0; i < n; ++i) {
        if (!Marked[i]) {
            init_nearest_neighbor(i);
        }
    }

    int j = 0;
    for (int i = 0; i < n; ++i) {
        if (!Marked[i] ) { 
            UnMarkedNodes[j++] = i; 
        }
    }

    int k1 = construct_NNG(j);

    if (show_info) {
        show_NNG_infomation(j, k1);
    }
}
/*********************************************************************/
/*  initialize                                                       */
/*********************************************************************/
void ZstTree::init_calc_whole_ZST() {
    // init_all_Node
    Curr_Npoints = m_nTerms + 1;
    for (int i = 0; i < m_nPoints; i++) {
        N_neighbors[i] = 0;
        Marked[i] = YES;
        if (i < m_nTerms) {
            build_NodeTRR(&Node[i]);
        } else {
            Node[i].segment.vertex[0].delay = NIL;
        }
    }
}
int ZstTree::init_ExG_DME() {
    int n = 0;
    for (int i = 0; i < m_nTerms; i++) {
        Marked[i] = NO;
        if (CHECK == YES) check_root_id(i);
        n++;
    }
    assert(n > 0);
    return n;
}
/*********************************************************************/
/*  select the cloest pair of merging segments.                      */
/*********************************************************************/
int pair_compare_inc(const void *p1, const void *q1) {

    PairType *p = (PairType*) p1, *q = (PairType*) q1;

    assert(p->x == p->x_root_id);
    assert(p->y == p->y_root_id);
    assert(q->x == q->x_root_id);
    assert(q->y == q->y_root_id);

    return (p->cost > q->cost) ? YES: NO;
}
int ZstTree::count_merge_pairs(int n_nodes) { 
    
    int n = 0, k;
    for (int i = 0; i < n_nodes; i++) {
        k = N_neighbors[i];
        if (!Marked[i] && k >0 ) {
            for (int j = 0; j < k; ++j) {
                Best_Pair[n].x = i;
                Best_Pair[n].y = The_NEIghors[i][j];
                Best_Pair[n].x_root_id = Node[Best_Pair[n].x].root_id;
                Best_Pair[n].y_root_id = Node[Best_Pair[n].y].root_id;
                Best_Pair[n].cost = Neighbor_Cost[i][j];
                n++;
            }
        }
    }
    assert(n > 1);
    return(n);
}
int ZstTree::init_pairs_to_merge(int n_nodes) { 
    int i, n;

    n = count_merge_pairs(n_nodes);

    qsort(Best_Pair, n, sizeof(PairType), pair_compare_inc); 

    if (SHOW_INFO) print_pairs_to_merge(n);

    return(n);
}
int ZstTree::pairs_to_merge(int n_nodes) { 
    int besta1, besta2, i1, i2;
    char *Flag;

    int n = init_pairs_to_merge(n_nodes);
    Flag = (char *) calloc(m_nPoints, sizeof(char));
    assert(Flag != NULL);

    for (int i = 0; i < n_nodes; i++) { Flag[i] = NO; }

    int j = 0;
    for (int i = 0; i < n; i++) {
        besta1 = Best_Pair[i].x;
        besta2 = Best_Pair[i].y;
        assert(besta1 >=0 && besta2>=0);
        i1 = Node[besta1].root_id;
        i2 = Node[besta2].root_id;
        assert(i1 == besta1);
        assert(i2 == besta2);
        if (!Flag[i1] && !Flag[i2]) {
            Flag[i1] = Flag[i2] = YES;
            Best_Pair[j] = Best_Pair[i];
            j++;
        /* Best_Pair[j++]=Best_Pair[i];  */
        }
    }
    assert(j > 0); 
    assert(n_nodes + j <= (int) m_nPoints); 
    free(Flag);

    if (SHOW_INFO) print_pairs_to_merge(j);
    return(j);
}
/*********************************************************************/
/*  merging process                                                  */
/*********************************************************************/
void ZstTree::update_MS(MsType *segment, PointType pt[4], PointType pts[2]) {
    DOUBLE max_x,max_y,min_x,min_y, t; 
    TrrType ms0,ms1,trr0,trr1,t0,t1;

    int line0type = calc_line_type(pt[0],pt[1]);
    int line1type = calc_line_type(pt[2],pt[3]);
    assert(line0type == MANHATTAN_ARC);
    assert(line1type == MANHATTAN_ARC);

    line2ms(&ms0, pt[0], pt[1]);
    line2ms(&ms1, pt[2], pt[3]);
    
    t = ms_distance(&ms0, &ms1);
    assert(t == segment->dist);  //assert(equivalent(t, area->dist, 1E-6)); 
    segment->dist = t;
    segment->L_Trr = ms0;
    segment->R_Trr = ms1;
    build_trr(&ms0, segment->dist, &trr0);
    build_trr(&ms1, segment->dist, &trr1);
    t0.MakeIntersect(trr1, ms0);
    t0.z = ms0.z;
    t1.MakeIntersect(trr0, ms1);
    t1.z = ms1.z;
    ms2line(&t0, &(segment->L_line[0]), &(segment->L_line[1]) );
    ms2line(&t1, &(segment->R_line[0]), &(segment->R_line[1]) );

    if (CHECK == YES) check_update_MS(segment, pt, line0type);
}
void ZstTree::calc_MS_sub2(MsType *segment, PointType pt[4]) {
    DOUBLE ldist, old_dist, old_area, new_area;
    PointType pts[2],  old_JS[2][2]; 
    TrrType old_L_MS,  old_R_MS; 

    ldist = linedist(pt[0], pt[1], pt[2], pt[3], pts);
    if (ldist == segment->dist) {
        assert(0);
    } else if (ldist < segment->dist )  {
        segment->dist = ldist;
        update_MS(segment, pt, pts);
    }
    assert(calc_line_type(segment->L_line[0], segment->L_line[1]) == MANHATTAN_ARC );
    assert(calc_line_type(segment->R_line[0], segment->R_line[1]) == MANHATTAN_ARC );
    check_JS_MS(*segment);
}
void calc_BS_located(PointType *pt, MsType *segment, PointType *p1, PointType *p2) {
    int found;
    int n = segment->npts;
    for (int i = 0; i < n; ++i) {
        *p1 = segment->vertex[i];
        *p2 = segment->vertex[(i+1)%n];
        found = PT_on_line_segment(pt, *p1, *p2);
        if (found) return;
    }
    //print_Point(stdout, *pt);
    //print_area(area);
    assert(0);
}
void calc_pt_delays(MsType *segment, PointType *q1, PointType q0, PointType q2) {
    
    assert(pt_on_line_segment(*q1, q0, q2));
    DOUBLE d1 = Point_dist(q0,*q1);
    DOUBLE a = (q0.x - q2.x)._ABS();
    DOUBLE b = (q0.y - q2.y)._ABS();
    DOUBLE d2 = a + b;
    
    if (d1 == 0) {
        q1->delay = q0.delay;
    } else if (same_Point(*q1, q2)) {
        q1->delay = q2.delay;
    } else if (a == b) { /* p0p1 is a Manhattan arc, a point */
        assert(q0.delay == q2.delay);
        q1->delay = q0.delay = q2.delay;
        q1->delay = q0.delay = q2.delay;
    } else {
        assert(0);
    }
}
void ZstTree::calc_MS_delay(MsType *segment, MsType *segment_L, MsType *segment_R) {
    PointType p1,p2;
    for (int i = 0; i < 2; ++i) {
        calc_BS_located(&(segment->L_line[i]), segment_L, &p1, &p2);
        calc_pt_delays(segment_L, &(segment->L_line[i]), p1, p2);

        calc_BS_located(&(segment->R_line[i]), segment_R, &p1, &p2);
        calc_pt_delays(segment_R, &(segment->R_line[i]), p1, p2);
    }
}
void ZstTree::calc_MS_sub1(MsType *segment, MsType *segment_L, MsType *segment_R) {
    PointType pt[4], q[2];

    int mn1 = segment_L->npts;
    int mn2 = segment_R->npts;
    int k0 = 0;
    int k1 = 1 % mn1; 
    pt[0] = segment_L->vertex[k0];
    pt[1] = segment_L->vertex[k1];
    int k2 = 0;
    int k3 = 1 % mn2;
    pt[2] = segment_R->vertex[k2];
    pt[3] = segment_R->vertex[k3];
    calc_MS_sub2(segment, pt); 
    calc_MS_delay(segment, segment_L, segment_R);
}
void ZstTree::calc_MS(MsType *segment, MsType *segment_L, MsType *segment_R, int v1, int v2) {
    int n1,n2;
    int j0,j1,j2,j3;
    PointType pt[4];

    segment->dist = DBL_MAX;
    n1 = segment_L->npts;
    n2 = segment_R->npts;
    if (n1 == 1  && n2 == 1) {
        segment->L_line[0] = segment->L_line[1] = segment_L->vertex[0];  
        segment->R_line[0] = segment->R_line[1] = segment_R->vertex[0];  
        segment->L_Trr.MakeDiamond(segment_L->vertex[0], 0);
        segment->R_Trr.MakeDiamond(segment_R->vertex[0], 0);
        segment->dist = Point_dist(segment_L->vertex[0], segment_R->vertex[0]) ; 
    } else {
        if (v1 < 0 && v2 < 0) {
            calc_MS_sub1(segment, segment_L, segment_R); 
        } else {
            assert(0);
        }
    }
    assert(calc_line_type(segment->L_line[0], segment->L_line[1]) == MANHATTAN_ARC);
    assert(calc_line_type(segment->R_line[0], segment->R_line[1]) == MANHATTAN_ARC);
    check_JS_MS(*segment);
}
void ZstTree::align_MS(PointType pt[2]) {
    PointType p;
    if (pt[0].y == pt[1].y) {
        if (pt[0].x < pt[1].x ) {
            p = pt[0];
            pt[0] = pt[1];
            pt[1] = p;
        }
    } else if ( pt[0].y < pt[1].y ) {
        p = pt[0];
        pt[0] = pt[1];
        pt[1] = p;
    }
}
void ZstTree::MS_processing(MsType *segment) {
    int n_MS_L, n_MS_R;
    if (same_Point(segment->L_line[0], segment->L_line[1])) {
        n_MS_L = 1;
    } else {
        n_MS_L = 2;
    }
    if (same_Point(segment->R_line[0], segment->R_line[1])) {
        n_MS_R = 1;
    } else {
        n_MS_R = 2;
    }
    if (n_MS_L==2) align_MS(segment->L_line);
    if (n_MS_R==2) align_MS(segment->R_line);
}
void construct_TRR_mr(MsType *segment, TrrType *ms0, TrrType *ms1, DOUBLE d0,
                     DOUBLE d1, int lor) {
    int i;
    TrrType trr0, trr1, trr;

    build_trr(ms0, d0, &trr0);
    build_trr(ms1, d1, &trr1);
    trr.MakeIntersect(trr0, trr1);
    if (lor == 0) {
        trr.z = trr0.z;
    } else if (lor == 1) {
        trr.z = trr1.z;
    } else if (lor == -1) {
        assert(trr1.z == trr0.z);
        trr.z = trr1.z;
    } else {
        assert(0);
    }
    segment->npts = TRR2segment(&trr, segment->vertex);
    // vertex只赋值了第一个点
}
void ZstTree::MergeManhattanArc(MsType *segment, MsType *segment_L, MsType *segment_R, int force_lor) {
    DOUBLE d0, d1;

    assert(segment->dist == ms_distance(&segment->L_Trr, &segment->R_Trr));
    
    int z1 = segment_L->vertex[0].z;
    int z2 = segment_R->vertex[0].z;
    bool same_tier = (z1 == z2);
    int lor = calc_segment_EdgeLen(segment, &d0, &d1, force_lor);
    if (same_tier) {
        if (m_topoMode == MMMMODE) {
            
        } else {
            segment->vertex[0].z = z2;
        }
    } else {
        if (lor == 0) {
            segment->vertex[0].z = z1;
        } else if (lor == 1) {
            segment->vertex[0].z = z2;
        }
    }
    
    assert(d0 >= 0);
    assert(d1 >= 0);

    calc_merge_pt_delay(segment, d0, d1, same_tier, lor);
             
    construct_TRR_mr(segment, &segment->L_Trr, &segment->R_Trr, segment->L_EdgeLen, segment->R_EdgeLen, lor);
}
void ZstTree::MS_processing_sub2(MsType *segment) {
    if (a_segment_TRR(&segment->L_Trr)) {
        ms2line(&segment->L_Trr,  &(segment->L_line[0]), &(segment->L_line[1]) );
    }
    if (a_segment_TRR(&segment->R_Trr)) {
        ms2line(&segment->R_Trr,  &(segment->R_line[0]), &(segment->R_line[1]) );
    }
}
DOUBLE segment_merge_cost(MsType *segment) {

    DOUBLE t = segment->L_EdgeLen + segment->R_EdgeLen;
    if (t >= 0) {
        assert( t == segment->dist || t >= segment->dist);
        return(t);
    } else {
        return(segment->dist);
    }
}
void ZstTree::MergeSegment_sub(MsType *segment) {
    
    DOUBLE merge_cost = segment_merge_cost(segment);
    int z1 = segment->segment_L->vertex[0].z;
    int z2 = segment->segment_R->vertex[0].z;
    bool same_tier = (z1 == z2);

    segment->subtree_cost = merge_cost + segment->segment_L->subtree_cost + segment->segment_R->subtree_cost;
   
    if (same_tier) {
        segment->capac =  segment->segment_L->capac + segment->segment_R->capac + merge_cost * m_pucap;
        if (m_topoMode == MMMMODE && segment->vertex[0].z != segment->segment_L->vertex[0].z) {
            assert(segment->segment_L->vertex[0].z == segment->segment_R->vertex[0].z);
            segment->capac += m_TsvC;
        }
    } else {
        segment->capac =
            segment->segment_L->capac + segment->segment_R->capac + merge_cost * m_pucap+m_TsvC;
    }

    DOUBLE t = segment->L_EdgeLen + segment->R_EdgeLen;
    assert(t  >= segment->dist - FUZZ || t==-2.0);

}
void ZstTree::MergeSegment(MsType *segment, MsType *segment_L, MsType *segment_R, int force_lor) {

    segment->segment_L = segment_L;
    segment->segment_R = segment_R;
    segment->L_EdgeLen = segment->R_EdgeLen = NIL;
    MS_processing(segment); 

    check_MergeSegment(segment);

    MergeManhattanArc(segment, segment_L, segment_R, force_lor);

  	assert(Manhattan_arc(segment->L_line[0], segment->L_line[1]));
	assert(Manhattan_arc(segment->R_line[0], segment->R_line[1]));
    MS_processing_sub2(segment);

  	MergeSegment_sub(segment);
}
void ZstTree::Merge2Segmentcase1(MsType *segment, MsType *segment_L, MsType *segment_R, int position) {
    int a = segment_L->vertex->z;
    int b = segment_R->vertex->z;
    if (a == b){
        MergeSegment(segment, segment_L, segment_R, -1);
    }else{
        MergeSegment(segment, segment_L, segment_R, position);
    }
}
void ZstTree::Merge2Nodes(NodeType *node, NodeType *node_L,NodeType *node_R, int force_lor) 
{
    MsType *segment   = &(node->segment);
    MsType *segment_L = &(node_L->segment);
    MsType *segment_R = &(node_R->segment);
    if (!LOOK_AHEAD) {

        calc_MS(segment, segment_L, segment_R, -1, -1);
        MergeSegment(segment, segment_L, segment_R, force_lor);

    } else {
        assert(m_topoMode != MMMMODE);

        DOUBLE a,b,c,d;
        MsType *segment_L_L, *segment_L_R;
        MsType *segment_R_L, *segment_R_R;
        assert(node_L->PAR == node_R->PAR);
        int v = node_L->PAR;
        
        bool have_l = node_L->L != NIL && node_L->R != NIL;
        bool have_r = node_R->L != NIL && node_R->R != NIL;
        if (have_l)
        {
            a = segment_L->segment_L->vertex->z;
            b = segment_L->segment_R->vertex->z;
            segment_L_L = &(Node[node_L->L].segment);
            segment_L_R = &(Node[node_L->R].segment);
        }
        if (have_r)
        {
            c = segment_R->segment_L->vertex->z;
            d = segment_R->segment_R->vertex->z;
            segment_R_L=&(Node[node_R->L].segment);
            segment_R_R=&(Node[node_R->R].segment);
        }
        DOUBLE cost[4];
        for (int i = 0; i < 2; i++)
        {
            for (int j = 0; j < 2; j++)
            {
                if (have_l)
                {
                    calc_MS(segment_L, segment_L_L, segment_L_R,  -1,-1);
                    Merge2Segmentcase1(segment_L, segment_L_L, segment_L_R, i);
                }
                if (have_r)
                {
                    calc_MS(segment_R, segment_R_L, segment_R_R,  -1,-1);
                    Merge2Segmentcase1(segment_R, segment_R_L, segment_R_R, j);
                }
                calc_MS(segment, segment_L, segment_R,  -1,-1);
                MergeSegment(segment, segment_L, segment_R, -1);
                cost[i + (j << 1)] = skewchecking(v);
            }
        }
        int ct = 0;
        for (int i = 0; i < 4; i++)
        {
            if (cost[i] < cost[ct])
                ct = i;
        }
        
        if (have_l)
        {
            calc_MS(segment_L, segment_L_L, segment_L_R,  -1,-1);
            Merge2Segmentcase1(segment_L, segment_L_L, segment_L_R, ct % 2);
        }
        if (have_r)
        {
            calc_MS(segment_R, segment_R_L, segment_R_R,  -1,-1);
            Merge2Segmentcase1(segment_R, segment_R_L, segment_R_R, ct >> 1);
        }
        calc_MS(segment, segment_L, segment_R,  -1,-1);
        MergeSegment(segment, segment_L, segment_R, -1);
    }
}
void ZstTree::Merge2Trees(int v, int L, int R) {
    Node[L].PAR = v;
    Node[R].PAR = v;
    Node[v].L = L;
    Node[v].R = R;
    Merge2Nodes(&(Node[v]), &(Node[L]), &(Node[R]), -1); 
    build_NodeTRR(&(Node[v]));
}
/*********************************************************************/
/*  embedding                                                        */
/*********************************************************************/
void check_embedding(PointType *p0, PointType *p1, MsType *area) {
    PointType tmp_pt0, tmp_pt1;
    calc_BS_located(p0, area, &tmp_pt0, &tmp_pt1);
    calc_BS_located(p1, area, &tmp_pt0, &tmp_pt1);
}
void ZstTree::embedding_sub(int p, int v, PointType pt1, PointType pt2, DOUBLE edgelen,
                   MsType *segment) {
    DOUBLE dist, a, b, old_dist;
    int unblocked; 
    TrrType trr1, trr2,ms;

    assert(p >= 0);

    PointType &pt = Node[p].m_mergePt;
    PointType &qt = Node[v].m_mergePt;

    a = (pt1.x-pt2.x)._ABS();
    b = (pt1.y-pt2.y)._ABS();
    if (a == 0 && b == 0) {  /* pt1 and pt2 are smae points */
        qt = pt1;  
    } else if (a == 0) {  /* vertical line */
        assert(0);
        qt       .x = pt1.x;
        qt       .y = pt       .y;
    } else if (b == 0) {  /* horizontal line */
        assert(0);
        qt       .x = pt       .x;
        qt       .y = pt1.y;
    } else {  /* other types of line segment */
        pt2linedist(pt, pt1, pt2, &(qt));
    }
    assert(pt_on_line_segment(qt, pt1,pt2));

    old_dist = dist = Point_dist(pt, qt);

    assert(edgelen >= 0);

    /* detour wiring required */
    if (edgelen < dist - FUZZ) {
    printf("warning: p=%d:v= %d, edgelen = %.9f, dist = %.9f\n", 
            p,v, edgelen.value, dist.value);
    }
    assert(edgelen >= dist - 100*FUZZ);
    EdgeLength[v] = tMAX(edgelen, dist);
}
void ZstTree::embedding(int p, int child) {
    PointType p0, p1;
    MsType *segment;
    DOUBLE t, edgelen;
    int v;

    if (child==0) { v = Node[p].L; } else { v = Node[p].R; }
    segment = &(Node[p].segment);
    if (child==0) {
        edgelen = segment->L_EdgeLen;
        p0 = segment->L_line[0];
        p1 = segment->L_line[1];
    } else {
        edgelen = segment->R_EdgeLen;
        p0 = segment->R_line[0];
        p1 = segment->R_line[1];
    }

    segment = &(Node[v].segment);

    check_embedding(&p0, &p1, segment);

    embedding_sub(p, v, p0, p1, edgelen, segment);

    if (v >= m_nTerms) embedding(v,0);
    if (v >= m_nTerms) embedding(v,1);
}
void calc_segment_center(MsType *segment, PointType *pt) {
    DOUBLE x,y;

    x = y = 0.0;
    int n = segment->npts;
    for (int i = 0; i < n; ++i) {
        x = x + segment->vertex[i].x;
        y = y + segment->vertex[i].y;
    }
    pt->x = x / n;
    pt->y = y / n;
    pt->z = segment->vertex[0].z;
}
void ZstTree::TopDown_Embedding(int v) {

    PointType &pt = Node[v].m_mergePt ;
    NodeType *from  = SuperRootNode () ;
    int superRoot = SuperRootNodeIndex () ;

    PointType &qt = from->m_mergePt;
    MsType *segment = &(Node[v].segment);
    calc_segment_center(segment, &pt);
    qt.x = pt.x;
    qt.y = pt.y;
    qt.z = pt.z;
    EdgeLength[superRoot] = 0;

    if (v >= m_nTerms) embedding(v,0);
    if (v >= m_nTerms) embedding(v,1);
}
/*********************************************************************/
/*  main functions of building ZST                                   */
/*********************************************************************/
void ZstTree::updateRootId(int root_id,int v) {
    if (v < 0)  return;

    NodeType *node = &Node[v] ;
    node->root_id = root_id;
    Neighbor_Cost[v][0] = DBL_MAX;
    
    updateRootId(root_id, node->L);
    updateRootId(root_id, node->R);
}
void ZstTree::init_marked(int v) {
    if (v < 0) return;
    Marked[v] = YES;
    init_marked(Node[v].L);
    init_marked(Node[v].R);
}
int ZstTree::ExG_DME(int &n_trees, int v, int show_info) {
    int n, i, j;
    while (n_trees > 1) {
        build_nearest_neighbor_graph(v, show_info);
        n = pairs_to_merge(v);
        if (show_info) {
            printf("\n#_trees=%d, #_nodes=%d, #_mering_pairs=%d\n",n_trees, v, n);
        }
        n_trees -= n;
        for (int k = 0; k < n; k++,v++) {
            i = Best_Pair[k].x;
            j = Best_Pair[k].y;
            Merge2Trees(v, i, j);
            updateRootId(v, v);
            Marked[v] = NO;
            init_marked(Node[v].L);
            init_marked(Node[v].R);
        }
    }
    return v;
}
void ZstTree::ZSTs_at_level(int show_info) {
    int n_trees;
    n_trees = init_ExG_DME();
    if (show_info) print_run(n_trees);
    if (n_trees > 1) {
        Curr_Npoints = ExG_DME(n_trees, Curr_Npoints, show_info);
    }
    assert(n_trees == 1);
}
void ZstTree::findzMinMax(const vector<NodeType> nodes, int &min1, int &max1) {
    // find the min&max z of a set
    int n = nodes.size();

    if (!n) {
        cerr << "Try to find z min&max for an empty node set!" << endl;
        return;
    }
    // initialization
    min1 = numeric_limits<int>::max();
    max1 = 0;
    for (int i = 0; i <= n-1; i++) {
        int iz = nodes[i].m_mergePt.z;
        if (iz < min1)
            min1 = iz;
        if (iz > max1)
            max1 = iz;
    }
    assert(min1 <= max1);
}
NodeType ZstTree::centerofset(const vector<NodeType > nodes) {
    // return the among the nodes
    int n = nodes.size();
    DOUBLE x = 0;
    DOUBLE y = 0;
    DOUBLE z = 0;
    NodeType cnode;
    PointType coor;
    if (n == 1)
        return nodes[0];
    else if (n == 0)
        assert(0);  //return nodes[NIL] ;
    for (int i = 0; i <= n - 1; i++) {
        x = x + nodes[i].m_mergePt.x;
        y = y + nodes[i].m_mergePt.y;
    }
    coor.x=x/n;
    coor.y=y/n;
    coor.z = -1; // to be replaced
    cnode.m_mergePt=coor;

    return cnode;
}
NodeType ZstTree::Zcut(const vector<NodeType> nodes,vector<NodeType>& nodes1, vector<NodeType>& nodes2) {
// Z cut for 3D nodes, return sroot, nodes1, nodes2
    NodeType sroot;
    int n = nodes.size();
    int zmin = numeric_limits<int>::max(), zmax = 0;
    int zstd = 0;
    if (n <= 1) {
        assert(0);
        // nodes1.clear();
        // nodes2.clear();
        // return sroot; // cannot be divided
    }
    // generate a sroot and mark it as merge point if true
    sroot = centerofset(nodes);

    findzMinMax(nodes, zmin, zmax);

    // determine which z should be used to divide nodes and sroot should be
    if (zmin > m_SrcZ)
        zstd = zmin;
    else if (zmax < m_SrcZ)
        zstd = zmax;
    else
        zstd = m_SrcZ;

    sroot.m_mergePt.z = zstd;
    // divides nodes
    for (int i = 0; i <= n - 1; i++) {
        int iz = nodes[i].m_mergePt.z;
        if (iz == zstd)
            nodes1.push_back(nodes[i]);
        else
            nodes2.push_back(nodes[i]); // assuming there are already nodes in srcz tier
    }
    // there is no node at zstd, divide in another way
    if (nodes1.empty()) {
        nodes2.clear();
        for (int i = 0; i <= n - 1; i++) {
            int iz = nodes[i].m_mergePt.z;
            if (iz < zstd)
                nodes1.push_back(nodes[i]);
            else
                nodes2.push_back(nodes[i]); // assuming there are already nodes in srcz tier
        }
    }

    return sroot;
}
vector<NodeType> ZstTree::subSet(const vector<NodeType>& nodes,
                         NodeType center, char direction) {
    // divide a set into two sets
    vector<NodeType> snodes;
    DOUBLE x0 = center.m_mergePt.x;
    DOUBLE y0 = center.m_mergePt.y;
    int i;
    int n = nodes.size();
    if (n <= 1)
        assert(0);
        //return snodes;
    for (i = 0; i <= n-1; i++) {
        switch (direction) {
            case 'l':
                if (nodes[i].m_mergePt.x <= x0)
                    snodes.push_back(nodes[i]);
                break;
            case 'r':
                if (nodes[i].m_mergePt.x> x0)
                    snodes.push_back(nodes[i]);
                break;
            case 'b':
                if (nodes[i].m_mergePt.y <= y0)
                    snodes.push_back(nodes[i]);
                break;
            case 't':
                if (nodes[i].m_mergePt.y> y0)
                    snodes.push_back(nodes[i]);
                break;
            default:
                cerr << "Wrong MMM subset!" << endl;
                break;
        }
    }
    return snodes;
}
int ZstTree::minN(const vector<NodeType>& nodes) {
    // min nsink within the same tier in nodes
    int n = numeric_limits<int>::max(), nodesize = nodes.size(), z;
    map<int, int> ns; // nsinks in different tiers
    map<int, int>::iterator it;

    // determine n for each tier
    for (int i = 0; i <= nodesize-1; i++) {
        z = nodes[i].m_mergePt.z;
        if (ns.count(z)) {
            // similar sinks have been found
            ns[z] = ns[z] + 1;
        } else
            // insert new tier
            ns.insert(pair<int, int>(z, 1));
    }

    // find the min n
    for (it = ns.begin(); it != ns.end(); it++) {
        if ((*it).second < n)
            n = (*it).second;
    }

    return n;

}
NodeType ZstTree::XYcut(const vector<NodeType> nodes,vector<NodeType>& nodes1, vector<NodeType>& nodes2, DOUBLE &TSVratio, char direction) {
    NodeType sroot;
    int zmin = 100, zmax = 0;
    int nmin1 = 0, nmin2 = 0;;

    if (nodes.size() == 1) {
        assert(0);
        // only one node
        // sroot = nodes[0];
        // return sroot;
    } else if (nodes.size() == 0) {
        assert(0);
        // empty nodes
        // return  nodes[NIL];
    }
    // generate a sroot and mark it as merge point if true
    sroot = centerofset(nodes);

    //determine z of sroot
    findzMinMax(nodes, zmin, zmax);
    if (zmin > m_SrcZ)
        sroot.m_mergePt.z = zmin;
    else if (zmax < m_SrcZ)
        sroot.m_mergePt.z = zmax;
    else
        sroot.m_mergePt.z = m_SrcZ;
    if (direction == 'x') {
        // use horizontal division for current set
        nodes1 = subSet(nodes, sroot, 'l'); // left subset
        nodes2 = subSet(nodes, sroot, 'r'); // right subset
        if (nodes1.size() == 0 || nodes2.size() == 0) {
            // x of all nodes are the same, y-cut should be tried
            nodes1 = subSet(nodes, sroot, 't'); // top subset
            nodes2 = subSet(nodes, sroot, 'b'); // bottom subset
        }
    }else {
        // use vertical division for current set
        nodes1 = subSet(nodes, sroot, 't'); // top subset
        nodes2 = subSet(nodes, sroot, 'b'); // bottom subset
        if (nodes1.size() == 0 || nodes2.size() == 0) {
            // y of all nodes are the same, x-cut should be tried
            nodes1 = subSet(nodes, sroot, 'l'); // top subset
            nodes2 = subSet(nodes, sroot, 'r'); // bottom subset
        }
    }
    // determine TSV bound ratio
    nmin1 = minN(nodes1);
    nmin2 = minN(nodes2);
    if (nmin1 == 0 || nmin2 == 0) {
        cerr << "Min # of nodes is 0 for a set of divided nodes!" << endl;
        return nodes[NIL];
    }
    TSVratio = DOUBLE(nmin1) / DOUBLE(nmin2);

    return sroot;
}
NodeType ZstTree::MMM3D(vector<NodeType> nodes, int ntsv, char direction, int &ii) {

    NodeType sroot, sroot1, sroot2;
    vector<NodeType> nodes1, nodes2; // subsets of nodes
    int TSVbound1 = 1, TSVbound2 = 1,m;
    DOUBLE TSVratio, n;
    char direction2 = direction;
    int zmin, zmax;
    
    if (nodes.empty()) {
        throw(" Error: MMM3D for EMPTY nodes!");
    } else if (nodes.size() == 1) {
        return nodes[0];
    }

    ii -= 1;
    // check z of nodes
    findzMinMax(nodes, zmin, zmax);
    if (ntsv == 0)
        throw(" Error: TSV bound for MMM3D is 0!");
    else if (ntsv == 1 && zmin != zmax) {
        //only one tsv can be used and the nodes are in different tiers, Z cut for 3D nodes
        sroot = Zcut(nodes, nodes1, nodes2);
        sroot.id = ii;
        TSVbound1 = 1;
        TSVbound2 = 1;
    }else {
        // XY cut, 2D MMM
        sroot = XYcut(nodes, nodes1, nodes2, TSVratio, direction); // ratio = TSVbound1/TSVbound2
        sroot.id = ii;
        //from->m_stnPt=sroot->m_stnPt;
        // from->m_stnPt=sroot->m_stnPt;
        n = (TSVratio / (TSVratio + 1)) * DOUBLE(ntsv);
        m = n._ROUND();
        TSVbound1 = max(m,1); // must be larger or equal to 1
        TSVbound2 = ntsv - TSVbound1;
        if (TSVbound2 == 0) {
            TSVbound2 = 1;
            TSVbound1 = max(ntsv - TSVbound2, 1); // TSVbound2 has been reset, TSVbound1 needs to be re-calculated again
        }
        if (direction == 'x')
            direction2 = 'y';
        else
            direction2 = 'x';
    }
    sroot2 = MMM3D(nodes2, TSVbound2, direction2,ii);
    sroot1 = MMM3D(nodes1, TSVbound1, direction2,ii);

    // generate TSVs if needed between sroot and sroot1
    sroot.L = sroot1.id;
    Node[sroot.id].L = sroot1.id;
    sroot.R = sroot2.id;
    Node[sroot.id].R = sroot2.id;
    Node[sroot1.id].PAR=sroot.id;
    Node[sroot2.id].PAR=sroot.id;

    Node[sroot.id].m_mergePt.z = sroot.m_mergePt.z;
    Node[sroot.id].segment.vertex[0].z = sroot.m_mergePt.z;

    return sroot;

}
void ZstTree::DME3D(int v) {
    int L,R;

    assert(v>=0);
    assert( (unsigned) v>=m_nTerms);
    assert(Node[v].segment.npts <= 0);

    L = Node[v].L;
    R = Node[v].R;
    assert( L>=0 && R>=0 && L < (int) m_nPoints && R < (int) m_nPoints);
    if (Node[L].segment.npts<=0) {
        DME3D(L);
    }
    if (Node[R].segment.npts<=0) {
        DME3D(R);
    }
    int force_lor;
    if (Node[v].m_mergePt.z == Node[L].m_mergePt.z) {
        force_lor = 0;
    } else if (Node[v].m_mergePt.z == Node[R].m_mergePt.z) {
        force_lor = 1;
    } else { assert(0); }
    Merge2Nodes(&(Node[v]), &(Node[L]), &(Node[R]), force_lor);
}
int ZstTree::MMM_DME() {
    int ii = m_nPoints;
    NodeType root = MMM3D(Snode, m_TsvBound, 'x', ii);
    cout << "TsvNum: " << TotalNumTsv(RootNodeIndex()) << endl;
    //root.m_mergePt.z = m_SrcZ; // change to source z
    DME3D(RootNodeIndex());
}
DOUBLE ZstTree::skewchecking(int v) {
    DOUBLE t, cost, cost2;

    TopDown_Embedding(v);
    calc_ZST_delay(v);
    
    int root = RootNodeIndex () ;
    if (v==root) set_SuperRoot();

    calc_TreeCost(v, &cost, &t);
    cost2 = Node[v].segment.subtree_cost;
    if (cost != cost2) {
        assert(cost == cost2);  //assert(equivalent(cost, cost2, 1E-5));
    }

    PointType &pt = Node[v].m_mergePt;

    return(cost);
}
void ZstTree::calc_whole_ZST(TwoStates show_info) {
    init_calc_whole_ZST();
    if (m_topoMode == MMMMODE) {
        MMM_DME();
    } else {
        ZSTs_at_level(show_info);
    }
    skewchecking( RootNodeIndex() );
    if (show_info) 
        print_answer(SinksFileName().c_str(), RootNodeIndex());
}
/*********************************************************************/
/*  fundamental functions to build tree                              */
/*********************************************************************/
DOUBLE ZstTree::TotalPower() {
    return Node[RootNodeIndex()].segment.capac * 1.2 * 1.2 * 1E9;
}
int ZstTree::sibling(int p, int par) {
    if (p==Node[par].L)
        return(Node[par].R);
    if (p==Node[par].R)
        return(Node[par].L);
    assert(0);
}
int ZstTree::TotalNumTsv(int v) {
    if (v < m_nTerms) return 0;

    int totalNum = 0;
    if (Node[v].L != NIL) totalNum += TotalNumTsv(Node[v].L);
    if (Node[v].R != NIL) totalNum += TotalNumTsv(Node[v].R);
    bool b1 = (Node[v].m_mergePt.z != Node[Node[v].L].m_mergePt.z);
    bool b2 = (Node[v].m_mergePt.z != Node[Node[v].R].m_mergePt.z);
    return totalNum + (b1 | b2);
}
DOUBLE ZstTree::TotalLength() {
    DOUBLE Tcost = 0;
    DOUBLE Tdist=0 ;
  
    int root = RootNodeIndex () ;
    calc_TreeCost(root, &Tcost, &Tdist);
    return Tcost;
}
ZstTree::ZstTree(const string inputSinksFileName, TopologyMode topoMode, TwoStates LOOK_AHEAD)
{
    this->m_inputSinksFileName = inputSinksFileName;
    this->m_topoMode = topoMode;
    this->LOOK_AHEAD = LOOK_AHEAD;
}
ZstTree::~ZstTree() {/* need something here... */};
void ZstTree::ConstructTree(void)
{
    read_input_file();
    print_header();
    calc_whole_ZST(NO);
}