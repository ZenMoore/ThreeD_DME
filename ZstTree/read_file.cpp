#include "ZstTree.h"

/*********************************************************************/
/*  utils                                    */
/*********************************************************************/
static void tSplit ( char *buf, vector< char * > tokv, const char *sep = 0 ) {
    if( sep == 0 ) {
        sep = "\t ";
    }
    tokv.clear();
    char *temp ;
    char *s  = buf ;
    
    while( ( temp = strtok( s, sep )) != 0 ) {
        tokv.push_back ( temp ) ;
        // cerr << temp  << " " ;
        s = 0 ;
    }
    // cerr << endl ;
}
/*********************************************************************/
/*  read header                                                      */
/*********************************************************************/
void ZstTree::set_K() {
    m_TsvC *= PUCAP_SCALE;
    m_pucap *= PUCAP_SCALE;
    m_K[INTERCONNECT] = m_pures * m_pucap * 0.5;
    m_K[TSV] = m_TsvR * m_TsvC * 0.5;
}
void ZstTree::Ex_DME_memory_allocation(void) {
    for (int i = 0; i < m_nPoints; i++) {
        Node.push_back( NodeType() );
        EdgeLength.push_back( DOUBLE(0) );
    }
    for (int i = 0; i < m_nTerms; i++) {
        Snode.push_back ( NodeType() ) ;
    }

    N_neighbors = (int *) calloc(m_nPoints, sizeof(int));
    assert(N_neighbors != NULL);

    Neighbor_Cost = (DOUBLE **) calloc(m_nPoints, sizeof(DOUBLE *));
    assert(Neighbor_Cost != NULL);

    for (int  i = 0; i < m_nPoints; i++) {
        Neighbor_Cost[i] = (DOUBLE *) calloc(N_Neighbor, sizeof(DOUBLE));
        assert(Neighbor_Cost[i] != NULL);
    }

    The_NEIghors = (int **) calloc(m_nPoints, sizeof(int *));
    assert(The_NEIghors != NULL);

    for (int i = 0; i < m_nPoints; i++) {
        The_NEIghors[i] = (int *) calloc(N_Neighbor, sizeof(int));
        assert(The_NEIghors[i] != NULL);
    }
}
void ZstTree::ExG_DME_memory_allocation() {
    Best_Pair = (PairType *) calloc(m_nPoints*N_Neighbor, sizeof(PairType));
    assert(Best_Pair != NULL);

    Marked = (char *) calloc(m_nPoints, sizeof(char));
    assert(Marked != NULL);

    UnMarkedNodes = (int *) calloc(m_nPoints, sizeof(int));
    assert(UnMarkedNodes != NULL);

    double t = sqrt((double) m_nPoints);
    int n = tMAX(1, (int) t);
    Bucket = (BucketType **) calloc(n, sizeof(BucketType *));
    assert(Bucket != NULL);

    for ( int i =0; i<n; i++) {
        Bucket[i] = (BucketType *) calloc(n, sizeof(BucketType));
        assert(Bucket[i] != NULL);
    }
    MAX_N_Index = n;
}
bool ZstTree::tReadInputHeader(const string fn) {
    const int MAX_STRING_LENGTH = 1024;
    char commentChar = '#';
    char line[MAX_STRING_LENGTH];

    ifstream inFile;
    inFile.open(fn.c_str());
    if (!inFile ) return false;

    double pures, pucap, TsvR, TsvC;
    int nterms, TsvBound, SrcZ;
    while ( inFile.getline(line, MAX_STRING_LENGTH) ) {
        char buf[ MAX_STRING_LENGTH ];
        strcpy ( buf , line ) ;
        vector <char*> tokv ;
        tSplit ( buf, tokv, " " ) ;
        char token [ MAX_STRING_LENGTH ] ;
        sscanf(line,"%s \n", token ) ;
        string key = token ;

        if (  token[0] == commentChar)  {
        } else if ( key=="NumPins" ) {
            sscanf( line,"%s : %d \n", token , &nterms ) ;
            SetNterms ( nterms ) ;
        } else if ( key=="PerUnitResistance" ) {
            sscanf( line,"%s : %lf \n", token , &pures ) ;
            SetPerUnitResistance ( pures ) ;
        } else if ( key=="PerUnitCapacitance" ) {
        sscanf( line,"%s : %lf \n", token , &pucap ) ;
        SetPerUnitCapacitance ( pucap ) ;
        } else if ( key=="TsvBound" ) {
            sscanf( line,"%s : %d \n", token , &TsvBound ) ;
            SetTsvBound ( TsvBound ) ;
        } else if ( key=="SrcZ" ) {
            sscanf( line,"%s : %d \n", token , &SrcZ ) ;
            SetSrcZ ( SrcZ ) ;
        } else if ( key=="TsvR" ) {
            sscanf( line,"%s : %lf \n", token , &TsvR ) ;
            SetTsvR ( TsvR ) ;
        } else if ( key=="TsvC" ) {
            sscanf( line,"%s : %lf \n", token , &TsvC ) ;
            SetTsvC ( TsvC ) ;
        }
    }
    set_K() ;

    Ex_DME_memory_allocation();
    ExG_DME_memory_allocation();

    inFile.close() ;
    return true ;
}
/*********************************************************************/
/*  read sink nodes                                                  */
/*********************************************************************/
void ZstTree::AddInputNode ( int i , DOUBLE x, DOUBLE y, int z, DOUBLE cap, DOUBLE delay) {
    NodeType &from = Node[i];
    from.id = from.root_id = i;

    MsType &segment = from.segment;
    segment.npts = 1;
    segment.capac = cap;

    PointType &pt = segment.vertex[0];
    pt.x = x;
    pt.y = y;
    pt.z = z;
    pt.delay = delay;
    from.m_mergePt = pt;
}
void ZstTree::AddSinkNode ( int i , DOUBLE x, DOUBLE y, int z, DOUBLE cap, DOUBLE delay) {
    NodeType &from = Snode[i];
    from.id = from.root_id = i;

    MsType &segment = from.segment;
    segment.npts = 1;
    segment.capac = cap;

    PointType &pt = segment.vertex[0];
    pt.x = x;
    pt.y = y;
    pt.z = z;
    pt.delay = delay;
    from.m_mergePt = pt;
}
bool ZstTree::tReadInputNodes ( const string fn ) {
    const int MAX_STRING_LENGTH = 1024;
    char line[ MAX_STRING_LENGTH ] ;

    ifstream inFile;
    inFile.open ( fn.c_str() ) ;
    if (!inFile ) return false ;

    int i = 0 ;
    int id = 0, z = 0;
    DOUBLE x = 0, y = 0, cap = 0, delay = 0;
    double a, b;

    while ( inFile.getline (line, MAX_STRING_LENGTH ) ) {

        if ( strlen( line ) == 0 ) continue ;
        char token [ MAX_STRING_LENGTH ] ;
        sscanf(line,"%s", token ) ;
        string key = token ;
        
        if (key=="Sink") {
            if ( i ) {
                AddInputNode ( id, x, y, z, cap, delay ) ;
                AddSinkNode ( id, x, y, z, cap, delay ) ;
                x = y = cap = delay = 0 ;
                id = z = 0;
            }
            sscanf(line,"%s :  %d \n", token , &id ) ;
            i++ ;
        } else if (  key == "Coordinate" ) {
            sscanf(line,"%s :  %lf %lf %d  \n", token , &a, &b, &z) ;
            x = a;
            y = b;
        } else if ( key=="Capacitive" ) {
            char token2 [ MAX_STRING_LENGTH ] ;
            sscanf(line,"%s %s : %lf \n", token, token2, &a ) ;
            cap = a;
            cap = cap*PUCAP_SCALE;
        } else if (  key == "Downstream_Delay" ) {
            sscanf(line,"%s :  %lf \n", token , &a ) ;
            delay = a;
            delay = delay*PUCAP_SCALE; // scale up 
        } else if (  key == "#" ) {
        } else if (  key == "NumPins" ) {
        } else if (  key == "PerUnitResistance" ) {
        } else if (  key == "PerUnitCapacitance" ) {
        } else if (  key == "TsvBound" ) {
        } else if (  key == "SrcZ" ) {
        } else if (  key == "TsvR" ) {
        } else if (  key == "TsvC" ) {
        } else {
            printf ( "Unknown keyword  %s \n", key.c_str() ) ;
            assert(0);
        }
    }
    if ( i ) {
        AddInputNode ( id, x, y, z, cap, delay ) ;
        AddSinkNode ( id, x, y, z, cap, delay ) ;
    }

    assert (i == m_nTerms);
    if (SHOW_INFO) print_sink_nodes(stdout);

    inFile.close() ;
    return true ;
}
/*********************************************************************/
/*  main funtion for reading file                                    */
/*********************************************************************/
void ZstTree::read_input_file(void)
{
    string fn = m_inputSinksFileName;
    bool ok = tReadInputHeader(fn) && tReadInputNodes(fn);
    if (!ok) {
        cerr << "cannot open file " << fn << endl;
        exit(0);
    }
}