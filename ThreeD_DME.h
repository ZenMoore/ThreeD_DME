#ifndef _THREED_DME_H
#define _THREED_DME_H

#include "ZstTree/ZstTree.h"
using namespace std;

class ThreeD_DME
{

public:
    ThreeD_DME(const string inputSinksFileName, ZstTree::TopologyMode topoMode, TwoStates LOOK_AHEAD) {
        m_tree = new ZstTree(inputSinksFileName, topoMode, LOOK_AHEAD);
    }
    ~ThreeD_DME() {
        delete m_tree; 
    }
    void ConstructTree(void) {
        m_tree->ConstructTree();
    }
    DOUBLE TotalLength(void) {
        return m_tree->TotalLength();
    }
    int TotalNumTsv(void) {
        return m_tree->TotalNumTsv(m_tree->RootNodeIndex());
    }
    DOUBLE TotalPower(void) {
        return m_tree->TotalPower();
    }

private:
    ZstTree *m_tree;
};

#endif