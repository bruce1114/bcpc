#ifndef BICLIQUETREE
#define BICLIQUETREE
#include<vector>
#include<stdint.h>

struct biCliqueTreeNode{
    std::vector<uint32_t> sons;
    int inValue,uv;
    int resMark;
    bool qpMark;
};

struct biCliqueTree{
    std::vector<biCliqueTreeNode> treeNodePool;
    uint32_t root;
    std::vector<int> firstlayerson;

    inline uint32_t getNewNode(){
        treeNodePool.emplace_back(biCliqueTreeNode());
        uint32_t size=treeNodePool.size();
        return size-1;
    }
};


#endif