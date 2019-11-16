#include "ZstTree.h"

/*********************************************************************/
/*  print best pairs                                                 */
/*********************************************************************/
void ZstTree::print_pairs_to_merge(int n) {
    printf("\n");
    printf("Best Pairs (total: %d):\n", n);
    for (int i = 0; i < n; i++) {
        printf("x: %d  y: %d  cost: %.0lf\n", Best_Pair[i].x, Best_Pair[i].y, Best_Pair[i].cost.value);
    }
}
/*********************************************************************/
/*  print the nearest neighbors                                      */
/*********************************************************************/
void ZstTree::print_NNG(int n) {
    printf("\n");
    int j;
    for (int i = 0; i < n; i++) {
        j = UnMarkedNodes[i];
        assert(N_neighbors[j] == 1);
        printf("id: %d  The_NEIghor: %d\n", j, The_NEIghors[j][0]);
    }
}
/*********************************************************************/
/*  show NNG information                                             */
/*********************************************************************/
void ZstTree::show_NNG_infomation(int n, int k1) {
    int i, j,k2,k3;
    DOUBLE t;

    for (i=0,t=0,k2=k3=0;i<n;++i) {
        j = UnMarkedNodes[i];
        if (N_neighbors[j]>0) {
        k2 += N_neighbors[j];
        k3 = tMAX(k3, N_neighbors[j]);
        }
        PointType &pt = Node[j].m_mergePt ;
        t = t + pt.t;
    }
    assert(t > 0);    
    printf("\n");
    printf("comparisons/nodes = %d/%d = %.1f \n", (int) t.value, n,  (t/n).value);
    printf("nodes_with_neighbor/nodes = %d/%d = %.2f\n", k1, n, ((DOUBLE)k1/n).value);
    printf("neighbors/nodes = %d/%d = %.2f\n", k2,n, ((DOUBLE)k2/n).value);
    printf("max_neighbors = %d \n", k3);
    printf("\n\n");
    fflush(stdout);
}
/*********************************************************************/
/*  print bucket                                                     */
/*********************************************************************/
void ZstTree::print_Bucket() {
    int s1 = 0, k, max_size, min_size;
    max_size = min_size = Bucket[0][0].num;
    for (int i = 0; i < N_Index; ++i) {
        for (int j = 0; j < N_Index; ++j) {
            k = Bucket[i][j].num;
            s1 += k;
            max_size = tMAX(max_size, k); 
            min_size = tMIN(min_size, k); 
        }
    }
    int s2 = 0;
    for (unsigned i = 0; i < m_nPoints; ++i) {
        if (!Marked[i])
            s2++;
    }

    int n_bucket = N_Index * N_Index; 
    double t1 = (double) n_bucket;
    double t2 = (double) s1/ (double) s2; 
    printf("Nodes/Bucket:ave:%f(max:%d,min:%d) --> Redundancy:%f \n", 
        s1/t1, max_size, min_size, t2); 
    printf("n_Bucket:%d    n_unmarked_Root:%d \n", n_bucket, s2); 
    printf("xlow: %.8lf  xhi: %.8lf  ylow: %.8lf  yhi: %.8lf\n", 
        m_ctrr.xlow.value, m_ctrr.xhi.value, m_ctrr.ylow.value, m_ctrr.yhi.value);
    for (int i = 0; i < N_Index; i++) {
        for (int j = 0; j < N_Index; j++) {
            printf("Bucket[%d][%d].num = %d\n", i, j, Bucket[i][j].num);
        }
    }
    if (s1 < s2)  {
        printf("Buckets contains %d elements \n", s1);
        printf("There are %d elements \n", s2);
    }
    assert(s1 >= s2);
}
/*********************************************************************/
/*  print run                                                        */
/*********************************************************************/
void ZstTree::print_run(int n_trees) {
    printf("============================================================\n");
    printf("=  run ExG_DME \n");
    printf("=  n_trees = %d \n", n_trees);
    printf("============================================================\n");
}
/*********************************************************************/
/*  printf sink nodes                                                */
/*********************************************************************/
void ZstTree::print_sink_nodes(FILE *f) {
    printf("\n");
    for (int i = 0; i < m_nTerms; i++) {
        fprintf(f, "id: %d  L: %d  R: %d  PAR: %d  root_id: %d\n", 
            Node[i].id, Node[i].L, Node[i].R, Node[i].PAR, Node[i].root_id);
    }
    printf("\n");
}
/*********************************************************************/
/*  print header                                                     */
/*********************************************************************/
void ZstTree::print_header() {

    printf("\n--------------------------------------------------------\n");
    printf("3D_DME %s (nterms=%d, Npoints=%d) ", m_inputSinksFileName.c_str(), m_nTerms, m_nPoints);
    printf("\n");
    if ( m_topoMode == GREEDYMODE ) {
        printf("Greedy mode");
    } else {
        printf("MMM mode");
    }
    printf("\n");
    printf("NumPins : %d\n", m_nTerms);
    printf("PerUnitResistance : %lf\n", m_pures.value);
    printf("PerUnitCapacitance : %.18lf\n", m_pucap.value);
    printf("TsvBound : %d\n", m_TsvBound);
    printf("SrcZ : %d\n", m_SrcZ);
    printf("TsvR : %lf\n", m_TsvR.value);
    printf("TsvC : %.16lf\n", m_TsvC.value);
    printf("--------------------------------------------------------\n");
    printf("\n"); 
    fflush(stdout);
}

/*********************************************************************/
/*  print answer                                                     */
/*********************************************************************/
DOUBLE ZstTree::print_answer(const char fn[], int v) {
    PointType p;
    DOUBLE Tcost = 0, Tdist=0, t, cap; 
    int n_detour_nodes = calc_TreeCost(v, &Tcost, &Tdist); 
    DOUBLE sCost = Node[v].segment.subtree_cost ;

    assert(Tcost == sCost);  //assert( equivalent(Tcost, sCost, 1.0));
    printf("\n%s:WL=%10.0f B=", fn, Tcost.value);

    printf("ps ");
    cap = Node[v].segment.capac;
    PointType &pt = Node[v].m_mergePt ;
    printf("t:%.1f (%.1f)(%.3f)", t.value, pt.delay.value, cap.value);
    t = 1E15/PUCAP_SCALE;
    printf("\n\n\n");
    printf("Total wirelength: %f\n", Tcost.value);
    printf("R: %f Oh\n", m_pures.value);
    printf("C: %f FF\n", (m_pucap*t).value);
    printf("\n");
    
    if (Tdist != Tcost) {
        t = (Tcost- Tdist)*100.0/Tdist;
        printf("Treewire = %.1f (detour:%.1f%%) ***", Tdist.value, t.value); 
        t = (double) (n_detour_nodes)*100.0/m_nTerms;
        printf("n_detour_nodes = %d (%f%%) ***\n", n_detour_nodes, t.value); 
    }
    fflush(stdout);

    return(Tcost);

}