#ifndef _ZSTTREE_H
#define _ZSTTREE_H

#include "header.h"

class ZstTree
{
public:
    enum TopologyMode {
        GREEDYMODE,
        MMMMODE
    };
public:
    // fundamental
    ZstTree(const string inputSinksFileName, TopologyMode topoMode, TwoStates LOOK_AHEAD);
    ~ZstTree();
    void ConstructTree(void);

    // read file
    void set_K(void);
    void Ex_DME_memory_allocation(void);
    void ExG_DME_memory_allocation(void);
    bool tReadInputHeader(const string fn);

    void AddInputNode ( int i , DOUBLE x, DOUBLE y, int z, DOUBLE cap, DOUBLE delay);
    void AddSinkNode ( int i , DOUBLE x, DOUBLE y, int z, DOUBLE cap, DOUBLE delay);
    bool tReadInputNodes (const string fn);
    void read_input_file(void);

    // main steps
    void init_calc_whole_ZST();
    void calc_whole_ZST(TwoStates show_info);
    void ZSTs_at_level(int show_info);
    int init_ExG_DME(void);
    int ExG_DME(int &n_trees, int v, int show_info);

    // build nearest neighbor graph
    int xlow_index(NodeType node);
    int xhi_index(NodeType node);
    int ylow_index(NodeType node);
    int yhi_index(NodeType node);
    void BuildMainTrr(int n);
    void bucket_partitioning_sub(int v);
    void bucket_partitioning(int n);
    void build_nearest_neighbor_graph(int n, int show_info);
    void init_nearest_neighbor(int v);
    void compare_neighbors_in_bucket(int v,int x,int y);
    void do_compare_neighbors(int i, int j);  
    void compare_neighbors(int i, int j);
    void update_neighbors(int x,int y, DOUBLE cost);

    void calc_nearest_neighbor_sub(int v,int out_xlow,int out_xhi,int out_ylow, int out_yhi);
    int calc_nearest_neighbor(int v,int inc);
    int construct_NNG(int n);
    
    // calculate
    DOUBLE calc_merging_cost(NodeType *node_L, NodeType *node_R);
    int calc_segment_EdgeLen(MsType *mergesegment, DOUBLE *d0, DOUBLE *d1, int lor);  
    DOUBLE calc_Elmore_merge_distance(DOUBLE cap1,DOUBLE delay1, DOUBLE cap2, DOUBLE delay2,DOUBLE d,DOUBLE *d1,DOUBLE *d2, int lor, int *uod);  
    DOUBLE calc_Elmore_merge_distance_sub(DOUBLE cap1,DOUBLE delay1, DOUBLE cap2, DOUBLE delay2, DOUBLE d, int lor);
    DOUBLE calc_merge_pt_delay_sub(MsType *mergesegment, DOUBLE d0, DOUBLE d1, bool same_tier, int lor, int uod);
    void calc_merge_pt_delay(MsType *segment, DOUBLE d0, DOUBLE d1, bool same_tier, int lor);
    int calc_merge_distance(DOUBLE cap1, DOUBLE delay1, DOUBLE cap2, DOUBLE delay2, DOUBLE d, DOUBLE *d1,DOUBLE *d2, bool same_tier, int *uod, int force_lor);
    DOUBLE calc_Elmore_merge_distance_detour1(DOUBLE cap1,DOUBLE delay1, DOUBLE cap2, DOUBLE delay2, 
                                              DOUBLE d, DOUBLE *d1, DOUBLE *d2, int lor, int uod);
    DOUBLE calc_Elmore_merge_distance_detour2(DOUBLE cap1,DOUBLE delay1, DOUBLE cap2, DOUBLE delay2, 
                                              DOUBLE d, DOUBLE *d1, DOUBLE *d2, int lor, int uod);
    void calc_ZST_delay(int v);
    void calc_ZST_delay_sub(int v);
    int calc_TreeCost(int v, DOUBLE *Tcost, DOUBLE *Tdist);
    void calc_TreeCost_sub(int par, int v, DOUBLE *Tcost, DOUBLE *Tdist, int *n_detour);

    DOUBLE calc_edge_delay(int v, int par, int lor, int uod);
    DOUBLE _pt_delay_increase(DOUBLE leng, DOUBLE cap, PointType *q0,PointType *q1, int lor, int uod, int i, bool same_tier);
    DOUBLE calc_delay_increase(DOUBLE cap, DOUBLE leng, DOUBLE pures, DOUBLE pucap);
    
    // select pairs
    int pairs_to_merge(int n_nodes);
    int init_pairs_to_merge(int n_nodes);
    int count_merge_pairs(int n_nodes);

    // merge
    void Merge2Trees(int v, int L, int R);
    void Merge2Nodes(NodeType *node, NodeType *node_L, NodeType *node_R, int force_lor);
    void MergeSegment(MsType  *segment, MsType *segment_L, MsType *segment_R, int force_lor);
    void MergeSegment_sub(MsType *segment);
    void Merge2Segmentcase1(MsType *segment, MsType *segment_L, MsType *segment_R, int position);

    void calc_MS(MsType *segment, MsType *segment_L, MsType *segment_R, int v1, int v2);
    void calc_MS_sub2(MsType *segment, PointType pt[4]);
    void calc_MS_sub1(MsType *segment, MsType *segment_L, MsType *segment_R);
    void update_MS(MsType *segment, PointType pt[4], PointType pts[2]);
    void calc_MS_delay(MsType *segment, MsType *segment_L, MsType *segment_R);

    void MS_processing(MsType *segment);
    void align_MS(PointType pt[2]);
    void MergeManhattanArc(MsType *mergesegment, MsType *mergesegment_L, MsType *mergesegment_R, int force_lor);
    void MS_processing_sub2(MsType *mergesegment);

    void updateRootId(int root_id, int v);
    void init_marked(int v);
    
    DOUBLE skewchecking(int v);
    void TopDown_Embedding(int v);
    void embedding_sub(int p, int v, PointType pt1, PointType pt2, DOUBLE edgelen, MsType *mergesegment);
    void embedding(int p, int child);
    

    
    
    
    // void calc_MS_sub2(MsType *mergesegment, PointType pt[4]);
    // void calc_MS_sub1(MsType *mergesegment, MsType *mergesegment_L, MsType *mergesegment_R);
    

    
    
    
  
    // check
    void count_tree_nodes(int root, int v, int *n);
    void check_root_id(int root);
    void check_compare_neighbors(int i, int j);
    void check_update_MS(MsType *segment, PointType pt[4], int line0type);
    void check_JS_MS(MsType segment);
    void check_MergeSegment(MsType *segment);

    // print
    void print_header(void);
    void print_sink_nodes(FILE *f);
    void print_run(int n_trees);
    void print_Bucket();
    void print_NNG(int n);
    void show_NNG_infomation(int n, int k1);
    void print_pairs_to_merge(int j);
    DOUBLE print_answer(const char fn[], int v);

    DOUBLE TotalLength(void);
    DOUBLE TotalPower();
    int sibling(int p, int q);
    int TotalNumTsv(int v);

    // MMM_DME
    void DME3D(int v);
    NodeType MMM3D(vector <NodeType > nodes  , int ntsv, char direction,int &ii);
    void findzMinMax(const vector<NodeType> nodes, int& min1,int& max1);
    NodeType centerofset(const vector<NodeType > nodes);
    int minN(const vector<NodeType>& nodes);
    vector<NodeType> subSet(const vector<NodeType>& nodes, NodeType center, char direction);
    NodeType XYcut(const vector<NodeType> nodes,vector<NodeType>& nodes1, vector<NodeType>& nodes2, DOUBLE &TSVratio, char direction);
    NodeType Zcut(const vector<NodeType> nodes,vector<NodeType>& nodes1, vector<NodeType>& nodes2);
    int MMM_DME(void);
    
    // void print_ZST_delay_error(int v, int L, int R, DOUBLE tL, DOUBLE tR);

    string SinksFileName() const {
        return m_inputSinksFileName;
    }
    TopologyMode getTopologyMode() const {
        return m_topoMode;
    }
    void SetNterms(int n) {
        m_nTerms = n;
        m_nPoints = n << 1;
    }
    int Nterms() const {
        return m_nTerms;
    }
    int Npoints() const {
        return m_nPoints;
    }
    int RootNodeIndex () const {
        return 2 * m_nTerms - 1 ;
    }
    NodeType* SuperRootNode ( ) {
        return &Node[ m_nTerms ] ;
    }
    int SuperRootNodeIndex () const {
        return m_nTerms ;
    }
    void set_SuperRoot() {
        int superRoot = SuperRootNodeIndex () ;
        int root = RootNodeIndex () ;
        
        Node[superRoot] = Node[root];
        Node[superRoot].PAR = Node[superRoot].R = NIL;
        Node[superRoot].L = root;
        Node[root].PAR = superRoot;
        Node[superRoot].id = superRoot;
        Node[superRoot].root_id = superRoot;
    }     
    void SetPerUnitResistance(DOUBLE pures) {
        m_pures = pures;
    }
    void SetPerUnitCapacitance(DOUBLE pucap) {
        m_pucap = pucap;
    }
    void SetTsvBound(int TsvBound) {
        m_TsvBound = TsvBound;
    }
    void SetSrcZ(int SrcZ) {
        m_SrcZ = SrcZ;
    }
    void SetTsvR(DOUBLE TsvR) {
        m_TsvR = TsvR;
    }
    void SetTsvC(DOUBLE TsvC) {
        m_TsvC = TsvC;
    }

private:
    enum COEFFICIENT {
        INTERCONNECT = 0,
        TSV
    };

    TwoStates LOOK_AHEAD = NO;

    // fundamental variables
    string m_inputSinksFileName;
    TopologyMode m_topoMode;

    // input variables
    int m_nTerms, m_nPoints;
    DOUBLE m_pures, m_pucap;
    int m_TsvBound, m_SrcZ;
    DOUBLE m_TsvR, m_TsvC;
    DOUBLE m_K[2];

    // "global" variables
    vector <NodeType> Node;
    vector <NodeType> Snode;
    vector <DOUBLE> EdgeLength;
    int *N_neighbors;
    DOUBLE **Neighbor_Cost;
    int **The_NEIghors;

    PairType *Best_Pair;
    char *Marked;
    int *UnMarkedNodes;
    BucketType **Bucket;

    int Curr_Npoints;

    TrrType m_ctrr;
    int N_Index, MAX_N_Index;
};

#endif