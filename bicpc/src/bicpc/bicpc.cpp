#ifndef BICPC_H
#define BICPC_H

#include "../bigraph/bigraph.cpp"
#include "../tools/linearSetThree.hpp"
#include "../tools/unionfind.cpp"
#include "../tools/utils.cpp"
#include "bicliquetree.cpp"
#include<iostream>
#include<chrono>
#include<unordered_map>
// #define MAX_N 10000
#ifdef LARGELOCAL
#define MAX_N 10000000
#else
#define MAX_N 10000
#endif
#define MAX_K 20

uint64_t combination_cache[MAX_N+1][MAX_K+1] = {0};

class BICPC{
public:
    bigraph* g;
    uint32_t p,q;
    uint32_t threshold[2];
    linearSetThree S[2];
    std::vector<std::vector<uint32_t>> ws;
    std::vector<std::vector<uint32_t>> ws2[2];
    biCliqueTree runtree;
    std::vector<uint32_t> validLeafNodes;
    std::vector<uint32_t> deg;

    std::vector<std::vector<uint32_t> > cliqueuv[2];
    std::vector<uint32_t> RUV[2];
    std::vector<int> RUVlook[2];
    std::vector<int> xStackNode;
    std::vector<int> curStackNode;
    std::vector<int> connectResults;
    uint32_t firstlayert;
    unf fans;


    //debug
    int connectBiCliqueInXzoneCnt=0;
    int bbranchcnt=0;

    int anchor,opAnchor,anchorThreshold,opThreshold;
    uint64_t qpcliqueNum;
    std::vector<std::vector<int>> opNeighbors;
    std::vector<uint32_t> RP;
    std::vector<std::vector<uint32_t>> two_hop_pos_end_layer;
    std::vector<std::vector<uint32_t>> op_vertices;
    std::vector<std::vector<uint32_t>> anchorCandLayer;
    std::vector<uint32_t> seed_vertices;
    std::vector<uint32_t> temp_comb;
    std::vector<int> anchorlab;

    std::vector<std::vector<int>> v2mbc[2];
    std::vector<int> commonMbclique;
    std::vector<int> commonMbcliquePos;
    std::vector<bool> visGroup;
    std::vector<int> visGroups;

    int mbeTreeNodeCnt=0;
    int pqnodecnt=0;
    int vpqnodecnt=0;
    uint64_t pqnodedepth_total=0;
    uint64_t bicliquelistNodecnt=0;


    void baseline();
    void qpbcl();
    void qbcpcLeafRes();
    void qbcpcLeafResforabnode();// for virtual abnode, depth of real abnode
    void qbcpcLeafResbase();
    void qbcpcLeafRespoor();
    void qbcpcLeafRespoorCntpqnode();
    void conductQpbcl(int left);
    void listMBCliques();
    void listMBCliquesTreeFullPivot();
    void listMBCliquesTree();
    void bbranchpureFullTree(uint32_t deep,uint32_t side,uint32_t inValue,int fid);
    void bbranchpure(uint32_t deep,uint32_t side,uint32_t inValue);
    void outputBCPCs();
    void quasi2finalBCPC();
    void connectBiCliqueInXzone(int side);
    int connectBiCliqueInCurZone(int side);
    int connectAllBiClique(int side,int thisNode,bool first);
    int connectAllBiCliqueNew(int side,int thisNode);
    int connectBiCliqueNewNoX(int side,int thisNode);
    int connectAllBiCliqueLeafRes();
    int connectAllBiCliqueLeafResbase();
    int getGroupNum();
    int bbranch_simplep(uint32_t deep,uint32_t side,uint32_t inValue);
    int bbranchLeafRes_simplep(uint32_t deep,uint32_t side,uint32_t inValue);
    int bbranchLeafResforabnode_simplep(uint32_t deep,uint32_t side,uint32_t inValue);
    int bbranchLeafResbase_simplep(uint32_t deep,uint32_t side,uint32_t inValue);
    int bbranchLeafRespoor_simplep(uint32_t deep,uint32_t side,uint32_t inValue);
    int bbranchLeafRespoorcntpqnode_simplep(uint32_t deep,uint32_t side,uint32_t inValue);
    void justConnect(uint32_t side);
    void reportTreeNode(int thisNodeId,int fid,int pivot,int pivotside);
    // void qpcoreReduction();

    void prepareAndListqpCliques();
    void listQPCliques(int left);
    double estimate_cost(int side);
    void find_combinations(uint64_t &comb_count, uint32_t beg_offset, uint32_t curr_depth, uint32_t total_depth);
    uint64_t count_combinations(uint64_t _n,uint64_t _k);

    void listAndConnectMBC();
    void listAndConnectMBCbase();
    void prepareV2MBC();
    void conductListAndConnectMBC(int left,int validPosR);
    void conductListAndConnectMBCbase(int left,int validPosR);
    bool checkCommonGroup(int validPosR);
    void prune_find_combinations(uint64_t &comb_count, uint32_t beg_offset, uint32_t curr_depth, uint32_t total_depth,int validPosR);
    void prune_find_combinations_base(uint64_t &comb_count, uint32_t beg_offset, uint32_t curr_depth, uint32_t total_depth,int validPosR);

    BICPC(const std::string filename,uint32_t rp,uint32_t rq){
        g=new bigraph(filename);
        std::cerr<<"read done"<<std::endl;
        p=rp,q=rq;
        threshold[0]=p;
        threshold[1]=q;

        S[0].resize(g->n[0]);
        S[1].resize(g->n[1]);
        deg.resize(std::max(g->n[0],g->n[1]),0);

        ws.resize(g->maxDu + g->maxDv);
        ws2[0].resize(g->maxDu + g->maxDv);
        ws2[1].resize(g->maxDu + g->maxDv);
        RUVlook[0].resize(g->n[0],0);
        RUVlook[1].resize(g->n[1],0);

        // fans.init(1281832);

        // std::cerr<<g->n[0]<<" "<<g->n[1]<<std::endl;
        // std::cerr<<S[0].sz<<" "<<S[1].sz<<std::endl;

    }

    BICPC(const std::string filename,uint32_t scala,uint32_t rp,uint32_t rq){
        g=new bigraph(filename,scala);
        std::cerr<<"read done"<<std::endl;
        p=rp,q=rq;
        threshold[0]=p;
        threshold[1]=q;

        S[0].resize(g->n[0]);
        S[1].resize(g->n[1]);
        deg.resize(std::max(g->n[0],g->n[1]),0);

        ws.resize(g->maxDu + g->maxDv);
        ws2[0].resize(g->maxDu + g->maxDv);
        ws2[1].resize(g->maxDu + g->maxDv);
        RUVlook[0].resize(g->n[0],0);
        RUVlook[1].resize(g->n[1],0);

        // fans.init(1281832);

        // std::cerr<<g->n[0]<<" "<<g->n[1]<<std::endl;
        // std::cerr<<S[0].sz<<" "<<S[1].sz<<std::endl;

    }

    ~BICPC(){
        delete g;
    }



    bool cmp(const int & a, const int & b){
        for(int i=0;i<cliqueuv[0][a].size()&&i<cliqueuv[0][b].size();++i){
            if(cliqueuv[0][a][i]!=cliqueuv[0][b][i]){
                return cliqueuv[0][a][i]<cliqueuv[0][b][i];
            }
        }
        return cliqueuv[0][a].size()<cliqueuv[0][b].size();
    }

    void outputBicliques(std::vector<uint32_t> * theOldLabels){
        //new labels to old labels
        for(int i=0;i<cliqueuv[0].size();++i){
            for(int j=0;j<cliqueuv[0][i].size();++j){
                cliqueuv[0][i][j]=theOldLabels[0][cliqueuv[0][i][j]];
            }
        }
        for(int i=0;i<cliqueuv[1].size();++i){
            for(int j=0;j<cliqueuv[1][i].size();++j){
                cliqueuv[1][i][j]=theOldLabels[1][cliqueuv[1][i][j]];
            }
        }

        for(int i=0;i<cliqueuv[0].size();++i){
            sort(cliqueuv[0][i].begin(),cliqueuv[0][i].end());
        }
        for(int i=0;i<cliqueuv[1].size();++i){
            sort(cliqueuv[1][i].begin(),cliqueuv[1][i].end());
        }

        int bicliquenum=cliqueuv[0].size();

        std::vector<int> indexarr;
        indexarr.resize(bicliquenum);
        for(int i=0;i<bicliquenum;++i){
            indexarr[i]=i;
        }

        sort(indexarr.begin(),indexarr.end());

        //print final results
        for(int i=0;i<indexarr.size();++i){
            int theindex=indexarr[i];
            for(int j=0;j<cliqueuv[0][theindex].size();++j){
                std::cout<<cliqueuv[0][theindex][j]<<" ";
            }
            std::cout<<std::endl;

            for(int j=0;j<cliqueuv[1][theindex].size();++j){
                std::cout<<cliqueuv[1][theindex][j]<<" ";
            }
            std::cout<<std::endl;
        }
    }

    int bbranch(uint32_t deep,uint32_t side,uint32_t inValue);
    int bbranchLeafRes(uint32_t deep,uint32_t side,uint32_t inValue);
    int bbranchLeafResforabnode(uint32_t deep,uint32_t side,uint32_t inValue);
    int bbranchLeafResbase(uint32_t deep,uint32_t side,uint32_t inValue);
    int bbranchLeafRespoor(uint32_t deep,uint32_t side,uint32_t inValue);
    int bbranchLeafRespoorcntpqnode(uint32_t deep,uint32_t side,uint32_t inValue);

    void qbcpc(){
        auto mt1 = std::chrono::steady_clock::now();

        uint32_t t=0, z=1;
        if(g->maxDu<g->maxDv) {t=1,z=0;}
        firstlayert=t;

        // std::cerr<<"t: "<<t<<" "<<"z: "<<z<<std::endl;
        // std::cerr<<"missing clique:"<<std::endl;
        // std::cerr<<g->labelsL[1]<<" "<<g->labelsL[184683]<<std::endl;
        // std::cerr<<g->labelsR[133490]<<std::endl;

        runtree.root=runtree.getNewNode();
        runtree.treeNodePool[runtree.root].inValue=-1;
        runtree.treeNodePool[runtree.root].resMark=-1;
        runtree.treeNodePool[runtree.root].qpMark=false;
        runtree.treeNodePool[runtree.root].uv=-1;

        runtree.firstlayerson.resize(g->n[t],-1);
        xStackNode.reserve(g->n[0]+g->n[1]);

        //first layer no pivot
        for(uint32_t u=0;u<g->n[t];u++){
            S[t].c=S[t].r=g->n[t];
            S[z].c=S[z].r=g->n[z];

            

            //form C set for z
            for(uint32_t i=g->pos[t][u];i<g->pos[t][u+1];++i){
                uint32_t v=g->e[t][i];
                S[z].swapByPos(--S[z].c,S[z].pos(v));
            }

            S[z].x=S[z].c;

            //form C set for t
            for(uint32_t i = g->pos[t][u]; i < g->pos[t][u + 1]; i++) {
                uint32_t v = g->e[t][i];
                if(g->pos[z][v + 1] > 0){
                    for(uint32_t j = g->pos[z][v + 1] - 1; j >= g->pos[z][v]; j--) {
                        uint32_t w = g->e[z][j];
                        
                        if(w > u) {
                            uint32_t pw = S[t].pos(w);
                            if(pw < S[t].c) {
                                S[t].swapByPos(--S[t].c, pw);
                            }
                        }
                        else break;

                        if(j == 0) break;
                    }
                }
            }
            S[t].x = S[t].c;

            //形成t的X集合
            for(uint32_t i = g->pos[t][u]; i < g->pos[t][u + 1]; i++) {
                uint32_t v = g->e[t][i];
                if(g->pos[z][v + 1] > 0){
                    for(uint32_t j = g->pos[z][v]; j < g->pos[z][v + 1]; j++) {
                        uint32_t w = g->e[z][j];
                        
                        if(w < u) {
                            uint32_t pw = S[t].pos(w);
                            if(pw < S[t].x) {
                                S[t].swapByPos(--S[t].x, pw);
                            }
                        }
                        else break;

                        if(j == 0) break;
                    }
                }
            }

            // if(u==1859991){
            //     for(int i=S[t].c;i<S[t].r;++i){
            //         std::cerr<<g->old_lables[0][S[t][i]]<<" ";
            //     }
            //     std::cerr<<std::endl;
            //     for(int i=S[z].c;i<S[z].r;++i){
            //         std::cerr<<g->old_lables[1][S[z][i]]<<" ";
            //     }
            //     std::cerr<<std::endl;
            // }


            if(S[z].CIsEmpty()) continue;

            RUV[t].push_back(u);
            RUVlook[t][u]=1;
            int res=bbranch(0,t,u);
            RUV[t].pop_back();
            RUVlook[t][u]=0;

            if(res!=-1){
                runtree.firstlayerson[u]=res;
                runtree.treeNodePool[runtree.root].sons.push_back(res);
            }

        }

        std::cerr<<"maxBiCliqueCount: "<<cliqueuv[0].size()<<std::endl;
        auto mt2 = std::chrono::steady_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(mt2 - mt1);
        std::cerr << "qbcpc time: " << duration.count() << " ms" << std::endl;
    }

    void collectBiClique(int startNode,std::vector<int>& holdCliqueId);
    void collectBiCliqueLeafRes(int startNode,std::vector<int>& holdCliqueId);
    void collectBiCliqueLeafResbase(int startNode,std::vector<int>& holdCliqueId);
    
};

void BICPC::prepareV2MBC(){

}

void BICPC::conductQpbcl(int left){
    if(left==2){
        for(int i=0;i<anchorCandLayer[left].size();++i){
            int u=anchorCandLayer[left][i];

            RP.push_back(u);
            op_vertices[2].clear();
            if(anchorThreshold<=2){
                for(int j=g->pos[anchor][u];j<g->pos[anchor][u+1];++j){
                    int v=g->e[anchor][j];
                    op_vertices[2].push_back(v);
                }
            }else{
                int j=0,k=g->pos[anchor][u];
                while(j<op_vertices[3].size()&&k<g->pos[anchor][u+1]){
                    if(op_vertices[3][j]==g->e[anchor][k]){
                        op_vertices[2].push_back(op_vertices[3][j]);
                        j++;
                        k++;
                    }else if(op_vertices[3][j]<g->e[anchor][k]){
                        j++;
                    }else{
                        k++;
                    }
                }
            }

            if(op_vertices[2].size()<opThreshold){
                RP.pop_back();
                continue;
            }

            for(int j=g->two_hop_pos[u];j<two_hop_pos_end_layer[2][u];++j){
                int v=g->two_hop_e[j];

                RP.push_back(v);
                op_vertices[1].clear();
                int k=0,o=g->pos[anchor][v];
                while(k<op_vertices[2].size()&&o<g->pos[anchor][v+1]){
                    if(op_vertices[2][k]==g->e[anchor][o]){
                        op_vertices[1].push_back(op_vertices[2][k]);
                        k++;
                        o++;
                    }else if(op_vertices[2][k]<g->e[anchor][o]){
                        k++;
                    }else{
                        o++;
                    }
                }
                if(op_vertices[1].size()<opThreshold){
                    RP.pop_back();
                    continue;
                }

                #ifdef COUNTQPONLY
                qpcliqueNum+=count_combinations(op_vertices[1].size(),opThreshold);
                #else
                seed_vertices.clear();
                for(int k=0;k<op_vertices[1].size();++k){
                    seed_vertices.push_back(op_vertices[1][k]);
                }
                uint64_t comb_count = 0;
                // std::vector< std::vector<int> > combinations;
                uint32_t comb_size = opThreshold;
                // vector<int> temp_comb(comb_size);
                // find_combinations(comb_count, 0, comb_size, comb_size);

                //shrink v2mbc
                for(int k=0;k<seed_vertices.size();++k){
                    int w=seed_vertices[k];
                    int tempnr=0;
                    for(int o=0;o<v2mbc[opAnchor][w].size();++o){
                        int bc=v2mbc[opAnchor][w][o];
                        // if(commonMbcliquePos[bc]<=nr){
                        //     exchange(v2mbc[opAnchor][w],o,tempnr++);
                        // }
                    }
                }

                //debug
                // if(RP[0]==1826028&&RP[1]==1934355&&RP[2]==1952117){
                //     std::cerr<<"here"<<std::endl;
                // }

                // prune_find_combinations(comb_count, 0, comb_size, comb_size,nr);
                #endif

                RP.pop_back();
            }

            RP.pop_back();
        }
        return;
    }
}

void BICPC::qpbcl(){
    qpcliqueNum=0;
    anchor=0;
    if(estimate_cost(0)>estimate_cost(1)) anchor=1;
    anchorThreshold=(anchor==0?p:q);
    opThreshold=(anchor==0?q:p);
    std::cerr<<"anchor: "<<anchor<<std::endl;

    g->prepareTowHopNeighbor(anchor,opThreshold);
    two_hop_pos_end_layer.resize(anchorThreshold+1);
    anchorCandLayer.resize(anchorThreshold+1);
    op_vertices.resize(anchorThreshold+1);
    seed_vertices.reserve(g->n[anchor^1]);
    temp_comb.resize(opThreshold);
    anchorlab.resize(g->n[anchor],anchorThreshold);

    for(int i=0;i<=anchorThreshold;++i){
        two_hop_pos_end_layer[i].resize(g->n[anchor]);
        op_vertices[i].reserve((anchor==0?g->maxDu:g->maxDv));
        anchorCandLayer[i].reserve(g->n[anchor]);
    }
    for(int u=0;u<g->n[anchor];++u){
        two_hop_pos_end_layer[anchorThreshold][u]=g->two_hop_pos[u+1];
        anchorCandLayer[anchorThreshold].push_back(u);
    }

    // listQPCliques(anchorThreshold);
    conductQpbcl(anchorThreshold);
}

void BICPC::justConnect(uint32_t side){
    if(RUV[0].size()>=threshold[0]&&RUV[1].size()>=threshold[1]){
        int res=connectAllBiCliqueNew(side,-1);
        connectResults.clear();
    }
}

int BICPC::bbranchLeafResforabnode_simplep(uint32_t deep,uint32_t side,uint32_t inValue){
    // if(S[1].sz!=g->n[1]){
    //     std::cerr<<"error"<<std::endl;
    // }
    // if(S[z].CIsEmpty()) return;

    int thisNode=-1;//temporally no tree node is built here
    int oriSide=side;

    uint32_t t=side;
    uint32_t z=side ^ 1;

    if(S[t].CIsEmpty()&&S[z].CIsEmpty()){
        if(S[t].XIsEmpty()&&S[z].XIsEmpty()){
            if(RUV[t].size()>=threshold[t]&&RUV[z].size()>=threshold[z]&&RUV[t].size()+RUV[z].size()>threshold[t]+threshold[z]){
                thisNode=runtree.getNewNode();
                runtree.treeNodePool[thisNode].inValue=inValue;
                runtree.treeNodePool[thisNode].uv=side;
                runtree.treeNodePool[thisNode].resMark=-1;
                runtree.treeNodePool[thisNode].qpMark=false;

                cliqueuv[t].emplace_back(RUV[t]);
                cliqueuv[z].emplace_back(RUV[z]);

                runtree.treeNodePool[thisNode].resMark=cliqueuv[t].size()-1;

                //check pqnode
                if(RUV[t].size()==threshold[t]&&t==oriSide||RUV[z].size()==threshold[z]&&z==oriSide){
                    pqnodecnt++;
                    pqnodedepth_total+=RUV[t].size()+RUV[z].size();
                }
            }


            return thisNode;
        }
        //check pqnode
        if(RUV[t].size()>=threshold[t]&&RUV[z].size()>=threshold[z]){
            if(RUV[t].size()==threshold[t]&&t==oriSide||RUV[z].size()==threshold[z]&&z==oriSide){
                // connectAllBiCliqueLeafResbase();
                vpqnodecnt++;
            }
        }
        return -1;
    }

    if(RUV[t].size()+S[t].cSize()<threshold[t]||RUV[z].size()+S[z].cSize()<threshold[z]){
        return thisNode;
    }

    //find pivot

    //debug
    // if(S[1].sz!=g->n[1]){
    //     std::cerr<<"error"<<std::endl;
    // }

    //initialize copies of xcr
    uint32_t x[2];
    uint32_t c[2];
    uint32_t r[2];
    x[t] = S[t].x; x[z] = S[z].x;
    c[t] = S[t].c; c[z] = S[z].c;
    r[t] = S[t].r; r[z] = S[z].r;

    int oriXstackSize=xStackNode.size();

    //find pivot
    // uint32_t maxI[2] = {0}, maxE[2] = {0}, maxUV[2] = {0};
    // maxI[0] = S[0].c;
    // maxI[1] = S[1].c;
    // maxUV[0] = S[0][S[0].c];
    // maxUV[1] = S[1][S[1].c];
    // t=0,z=1;
    uint32_t wsSize;


    //find pivot in C of t
    // for(uint32_t i = S[t].c; i < S[t].r; i++) {//scan all vertices in C
    // // for(uint32_t i = S[t].c; i <= S[t].c; i++) {
    //     uint32_t u = S[t][i];
    //     uint32_t e = 0;

    //     if(g->deg(u, t) > S[z].r - S[z].c) {
    //         for(uint32_t j = S[z].c; j < S[z].r; j++) {
    //             uint32_t v = S[z][j];
    //             if(g->connect(u, v, t)) {
    //                 e++;
    //                 deg[v]++;
    //             }
    //         }
    //     }
    //     else {
    //         for(uint32_t j = g->pos[t][u]; j < g->pos[t][u + 1]; j++) {
    //             uint32_t v = g->e[t][j];
    //             uint32_t pv = S[z].pos(v);
    //             if(S[z].c <= pv && pv < S[z].r) {//in C
    //                 e++;
    //                 deg[v]++;
    //             }
    //         }
    //     }

    //     if(e > maxE[z]) {
    //         maxE[z] = e;
    //         maxI[t] = i;
    //         maxUV[t] = u;
    //     }
    // }
    // //find pivot in C of z
    // for(uint32_t i = S[1].c; i < S[1].r; i++) {//scan all vertices in C
    // // for(uint32_t i = S[t].c; i <= S[t].c; i++) {
    //     uint32_t u = S[1][i];
    //     uint32_t e = deg[u];
    //     deg[u] = 0;
    //     if(e > maxE[0]) {
    //         maxE[0] = e;
    //         maxI[1] = i;
    //         maxUV[1] = u;
    //     }
    // }

    // //check if pivot in X of both sides
    // bool isPivotInX[2] = {false, false};
    // for(t = 0; t < 2; t++) {
    //     z = t ^ 1;//t + z = 1

    //     for(uint32_t i = S[t].x; i < S[t].c; i++) {//scan all vertices in X
    //     // for(uint32_t i = S[t].c; i <= S[t].c; i++) {
    //         uint32_t u = S[t][i];
    //         uint32_t e = 0;

    //         if(g->deg(u, t) > S[z].r - S[z].c) {
    //             for(uint32_t j = S[z].c; j < S[z].r; j++) {
    //                 uint32_t v = S[z][j];
    //                 if(g->connect(u, v, t)) {
    //                     e++;
    //                 }
    //             }
    //         }
    //         else {
    //             for(uint32_t j = g->pos[t][u]; j < g->pos[t][u + 1]; j++) {
    //                 uint32_t v = g->e[t][j];
    //                 uint32_t pv = S[z].pos(v);
    //                 if(S[z].c <= pv && pv < S[z].r) {//in C
    //                     e++;
    //                 }
    //             }
    //         }

    //         if(e > maxE[z]) {
    //             isPivotInX[t] = true;
    //             maxE[z] = e;
    //             maxI[t] = i;
    //             maxUV[t] = u;
    //         }
    //     }
    // }

    // //pivot found. pivot must have at least one neighbor in this algorithm
    // t = 1;
    // if(S[0].cSize() - maxE[0] >= S[1].cSize() - maxE[1]) {
    //     t = 0;
    // }
    // z = t ^ 1;


    // //r is more like the right limit of C
    // //r[z] doesnt follow S[z].r
    // //pivot could be in X, and C of z is changed, but never mind, it will recover by r[z]
    // if(g->deg(maxUV[t], t) > S[z].r - S[z].c) {
    //     if(S[z].r>0){
    //         for(uint32_t i = S[z].r - 1; i >= S[z].c; i--) {
    //             if(!g->connect(maxUV[t], S[z][i], t)) {
    //                 S[z].swapByPos(--S[z].r, i);
    //             }
    
    //             if(i == 0) break;
    //         }
    //     }
    // }else{
    //     uint32_t u = maxUV[t];
    //     uint32_t tmpR = S[z].c;
    //     for(uint32_t j = g->pos[t][u]; j < g->pos[t][u + 1]; j++) {
    //         uint32_t v = g->e[t][j];
    //         uint32_t pv = S[z].pos(v);
    //         if(S[z].c <= pv && pv < S[z].r) {
    //             S[z].swapByPos(tmpR++, pv);
    //         }
    //     }
    //     S[z].r = tmpR;
    // }
    // wsSize = r[z] - S[z].r;
    // if(ws[deep].size() < wsSize) {
    //     ws[deep].resize(wsSize * 2);
    // }
    // memcpy(ws[deep].data(), S[z].begin() + S[z].r, sizeof(uint32_t) * wsSize);


    // //iterate on pivot
    // if(isPivotInX[t]==false){
        
    //     S[t].swapByPos(maxI[t], --r[t]);
    //     S[t].r = r[t];

    //     //x[z] doesnt follow S[z].x, used to recover S[z].x
    //     if(g->deg(maxUV[t], t) > S[z].c - S[z].x) {
    //         for(uint32_t i = S[z].x; i < S[z].c; i++) {
    //             if(!g->connect(maxUV[t], S[z][i], t)) {
    //                 S[z].swapByPos(S[z].x++, i);
    //             }
    //         }
    //     }else{
    //         uint32_t u = maxUV[t];
    //         uint32_t tmpX = S[z].c;
    //         for(uint32_t j = g->pos[t][u]; j < g->pos[t][u + 1]; j++) {
    //             uint32_t v = g->e[t][j];
    //             uint32_t pv = S[z].pos(v);
    //             if(S[z].x <= pv && pv < S[z].c) {
    //                 S[z].swapByPos(--tmpX, pv);
    //             }
    //         }
    //         S[z].x = tmpX;
    //     }



    //     RUV[t].push_back(maxUV[t]);
    //     RUVlook[t][maxUV[t]]=1;

    //     //debug
    //     // if(S[1].sz!=g->n[1]){
    //     //     std::cerr<<"error"<<std::endl;
    //     // }
    //     // if(thisNode==401812){
    //     //     std::cerr<<"here"<<std::endl;
    //     // }

    //     // #ifdef PARTIALPIVOT
    //     // int res;
    //     // if(S[z].cSize()+S[t].cSize()<=20){
    //     //     res=bbranch_simplep(deep+1,t,maxUV[t]);
    //     // }else{
    //     //     res=bbranch(deep+1,t,maxUV[t]);
    //     // }
    //     // #else
    //     //C of z must not be empty
    //     int res=bbranchLeafRes(deep+1,t,maxUV[t]);
    //     // #endif

    //     //debug
    //     // if(S[1].sz!=g->n[1]){
    //     //     std::cerr<<"error"<<std::endl;
    //     // }

    //     RUV[t].pop_back();
    //     RUVlook[t][maxUV[t]]=0;
    //     if(res!=-1){
    //         if(thisNode==-1){
    //             thisNode=runtree.getNewNode();
    //             runtree.treeNodePool[thisNode].inValue=inValue;
    //             runtree.treeNodePool[thisNode].uv=side;
    //             runtree.treeNodePool[thisNode].resMark=-1;
    //             runtree.treeNodePool[thisNode].qpMark=false;
    //         }
    //         runtree.treeNodePool[thisNode].sons.push_back(res);
    //         xStackNode.push_back(res);
    //     }

    //     //debug
    //     // if(S[1].sz!=g->n[1]){
    //     //     std::cerr<<"error"<<std::endl;
    //     // }

    //     //put pivot into X
    //     S[t].swapByPos(c[t]++, r[t]++);
    // }




    //iterate on both sides
    int tt=side;
    uint32_t wsize[2];
    for(int o=0;o<=1;++o){
        tt=tt^1;
        int zz=tt^1;
        wsSize=r[tt]-c[tt];
        wsize[tt]=wsSize;
        if(ws2[tt][deep].size() < wsSize) {
            ws2[tt][deep].resize(wsSize * 2);
        }
        memcpy(ws2[tt][deep].data(), S[tt].begin() + c[tt], sizeof(uint32_t) * wsSize);

        for(uint32_t j = 0; j < wsSize; j++){
            uint32_t w = ws2[tt][deep][j];
            S[tt].swapByPos(S[tt].pos(w), --r[tt]);

            //recover next level parameter
            S[0].x = x[0]; S[1].x = x[1];
            S[0].c = c[0]; S[1].c = c[1];
            S[0].r = r[0]; S[1].r = r[1];

            if(g->deg(w, tt) > S[zz].c - S[zz].x) {
                for(uint32_t i = S[zz].x; i < S[zz].c; i++) {
                    if(!g->connect(w, S[zz][i], tt)) {
                        S[zz].swapByPos(S[zz].x++, i);
                    }
                }
            }
            else {
                uint32_t tmpX = S[zz].c;
                for(uint32_t i = g->pos[tt][w]; i < g->pos[tt][w + 1]; i++) {
                    uint32_t u = g->e[tt][i];
                    uint32_t pu = S[zz].pos(u);
                    if(S[zz].x <= pu && pu < S[zz].c) {
                        S[zz].swapByPos(--tmpX, pu);
                    }
                }
                S[zz].x = tmpX;
            }

            if(g->deg(w, tt) > S[zz].r - S[zz].c) {
                if(S[zz].r > 0){
                    for(uint32_t i = S[zz].r - 1; i >= S[zz].c; i--) {
                        if(!g->connect(w, S[zz][i], tt)) {
                            S[zz].swapByPos(--S[zz].r, i);
                        }

                        if(i == 0) break;
                    }
                }
            }else{
                uint32_t tmpR = S[zz].c;
                for(uint32_t j = g->pos[tt][w]; j < g->pos[tt][w + 1]; j++) {
                    uint32_t u = g->e[tt][j];
                    uint32_t pu = S[zz].pos(u);
                    if(S[zz].c <= pu && pu < S[zz].r) {
                        S[zz].swapByPos(tmpR++, pu);
                    }
                }
                S[zz].r = tmpR;
            }

            // //debug
            // if(deep==0&&inValue==1859991&&w==1911392){
            //     // std::cerr<<<<std::endl;
            //     for(int ii=S[z].c;ii<S[z].r;++ii){
            //         std::cerr<<g->old_lables[z][S[z][ii]]<<" ";
            //     }
            //     std::cerr<<std::endl;
                
            // }

            // if(S[zz].cSize()>0){
            RUV[tt].push_back(w);
            RUVlook[tt][w]=1;
            //debug
            // if(S[1].sz!=g->n[1]){
            //     std::cerr<<"error"<<std::endl;
            // }

            int res=bbranchLeafResforabnode_simplep(deep+1,tt,w);
            RUV[tt].pop_back();
            RUVlook[tt][w]=0;

            if(res!=-1){
                if(thisNode==-1){
                    thisNode=runtree.getNewNode();
                    runtree.treeNodePool[thisNode].inValue=inValue;
                    runtree.treeNodePool[thisNode].uv=side;
                    runtree.treeNodePool[thisNode].resMark=-1;
                    runtree.treeNodePool[thisNode].qpMark=false;
                }
                runtree.treeNodePool[thisNode].sons.push_back(res);
                xStackNode.push_back(res);
            }
            // }

            S[tt].swapByPos(c[tt]++,r[tt]++);//move w to X
        }
    }

    for(int o=0;o<=1;++o){
        for(uint32_t j = 0; j < wsize[o]; j++) {
            uint32_t w = ws2[o][deep][j];
            S[o].swapByPos(S[o].pos(w), --c[o]);
        }
    }

    //recover x,c,r
    S[t].x = x[t]; S[z].x = x[z];
    S[t].c = c[t]; S[z].c = c[z];
    S[t].r = r[t]; S[z].r = r[z];

    //connect that should be connected
    //almost the same as before
    if(RUV[t].size()>=threshold[t]&&RUV[z].size()>=threshold[z]){
        int res=-1;
        // if(RUV[t].size()==threshold[t]&&t==oriSide||RUV[z].size()==threshold[z]&&z==oriSide){// only the first one
        //     res=connectAllBiCliqueLeafResbase();
        //     if(thisNode!=-1){
        //         runtree.treeNodePool[thisNode].resMark=res;
        //     }
        // }else{
        //     if(thisNode!=-1){
        //         fans.updateSize(cliqueuv[0].size());
        //         for(int i=0;i<runtree.treeNodePool[thisNode].sons.size();++i){
        //             int sonNode=runtree.treeNodePool[thisNode].sons[i];
        //             connectResults.push_back(runtree.treeNodePool[sonNode].resMark);
        //         }

        //         for(int i=1;i<connectResults.size();++i){
        //             fans.merge(connectResults[i],connectResults[0]);
        //         }
        //         runtree.treeNodePool[thisNode].resMark=fans.find(connectResults[0]);
        //         connectResults.clear();
        //     }

        // }

        if(thisNode!=-1){
            fans.updateSize(cliqueuv[0].size());
            for(int i=0;i<runtree.treeNodePool[thisNode].sons.size();++i){
                int sonNode=runtree.treeNodePool[thisNode].sons[i];
                connectResults.push_back(runtree.treeNodePool[sonNode].resMark);
            }

            for(int i=1;i<connectResults.size();++i){
                fans.merge(connectResults[i],connectResults[0]);
            }
            runtree.treeNodePool[thisNode].resMark=fans.find(connectResults[0]);
            connectResults.clear();

            //check pqnode
            if(RUV[t].size()==threshold[t]&&t==oriSide||RUV[z].size()==threshold[z]&&z==oriSide){
                pqnodecnt++;
                pqnodedepth_total+=RUV[t].size()+RUV[z].size();
            }
        }else{
            if(RUV[t].size()==threshold[t]&&t==oriSide||RUV[z].size()==threshold[z]&&z==oriSide){
                vpqnodecnt++;
            }
        }


    }

    xStackNode.resize(oriXstackSize);
    return thisNode;
}

int BICPC::bbranchLeafRespoorcntpqnode_simplep(uint32_t deep,uint32_t side,uint32_t inValue){
    // if(S[1].sz!=g->n[1]){
    //     std::cerr<<"error"<<std::endl;
    // }
    // if(S[z].CIsEmpty()) return;

    int thisNode=-1;//temporally no tree node is built here
    int oriSide=side;

    uint32_t t=side;
    uint32_t z=side ^ 1;

    if(S[t].CIsEmpty()&&S[z].CIsEmpty()){
        if(S[t].XIsEmpty()&&S[z].XIsEmpty()){
            if(RUV[t].size()>=threshold[t]&&RUV[z].size()>=threshold[z]&&RUV[t].size()+RUV[z].size()>threshold[t]+threshold[z]){
                thisNode=runtree.getNewNode();
                runtree.treeNodePool[thisNode].inValue=inValue;
                runtree.treeNodePool[thisNode].uv=side;
                runtree.treeNodePool[thisNode].resMark=-1;
                runtree.treeNodePool[thisNode].qpMark=false;

                cliqueuv[t].emplace_back(RUV[t]);
                cliqueuv[z].emplace_back(RUV[z]);

                runtree.treeNodePool[thisNode].resMark=cliqueuv[t].size()-1;

                //check pqnode
                if(RUV[t].size()==threshold[t]&&t==oriSide||RUV[z].size()==threshold[z]&&z==oriSide){
                    pqnodecnt++;
                }
            }


            return thisNode;
        }
        // if(RUV[t].size()>=threshold[t]&&RUV[z].size()>=threshold[z]){
        //     if(RUV[t].size()==threshold[t]&&t==oriSide||RUV[z].size()==threshold[z]&&z==oriSide){
        //         connectAllBiCliqueLeafResbase();
        //     }
        // }
        return -1;
    }

    if(RUV[t].size()+S[t].cSize()<threshold[t]||RUV[z].size()+S[z].cSize()<threshold[z]){
        return thisNode;
    }

    //find pivot

    //debug
    // if(S[1].sz!=g->n[1]){
    //     std::cerr<<"error"<<std::endl;
    // }

    //initialize copies of xcr
    uint32_t x[2];
    uint32_t c[2];
    uint32_t r[2];
    x[t] = S[t].x; x[z] = S[z].x;
    c[t] = S[t].c; c[z] = S[z].c;
    r[t] = S[t].r; r[z] = S[z].r;

    int oriXstackSize=xStackNode.size();

    //find pivot
    // uint32_t maxI[2] = {0}, maxE[2] = {0}, maxUV[2] = {0};
    // maxI[0] = S[0].c;
    // maxI[1] = S[1].c;
    // maxUV[0] = S[0][S[0].c];
    // maxUV[1] = S[1][S[1].c];
    // t=0,z=1;
    uint32_t wsSize;


    //find pivot in C of t
    // for(uint32_t i = S[t].c; i < S[t].r; i++) {//scan all vertices in C
    // // for(uint32_t i = S[t].c; i <= S[t].c; i++) {
    //     uint32_t u = S[t][i];
    //     uint32_t e = 0;

    //     if(g->deg(u, t) > S[z].r - S[z].c) {
    //         for(uint32_t j = S[z].c; j < S[z].r; j++) {
    //             uint32_t v = S[z][j];
    //             if(g->connect(u, v, t)) {
    //                 e++;
    //                 deg[v]++;
    //             }
    //         }
    //     }
    //     else {
    //         for(uint32_t j = g->pos[t][u]; j < g->pos[t][u + 1]; j++) {
    //             uint32_t v = g->e[t][j];
    //             uint32_t pv = S[z].pos(v);
    //             if(S[z].c <= pv && pv < S[z].r) {//in C
    //                 e++;
    //                 deg[v]++;
    //             }
    //         }
    //     }

    //     if(e > maxE[z]) {
    //         maxE[z] = e;
    //         maxI[t] = i;
    //         maxUV[t] = u;
    //     }
    // }
    // //find pivot in C of z
    // for(uint32_t i = S[1].c; i < S[1].r; i++) {//scan all vertices in C
    // // for(uint32_t i = S[t].c; i <= S[t].c; i++) {
    //     uint32_t u = S[1][i];
    //     uint32_t e = deg[u];
    //     deg[u] = 0;
    //     if(e > maxE[0]) {
    //         maxE[0] = e;
    //         maxI[1] = i;
    //         maxUV[1] = u;
    //     }
    // }

    // //check if pivot in X of both sides
    // bool isPivotInX[2] = {false, false};
    // for(t = 0; t < 2; t++) {
    //     z = t ^ 1;//t + z = 1

    //     for(uint32_t i = S[t].x; i < S[t].c; i++) {//scan all vertices in X
    //     // for(uint32_t i = S[t].c; i <= S[t].c; i++) {
    //         uint32_t u = S[t][i];
    //         uint32_t e = 0;

    //         if(g->deg(u, t) > S[z].r - S[z].c) {
    //             for(uint32_t j = S[z].c; j < S[z].r; j++) {
    //                 uint32_t v = S[z][j];
    //                 if(g->connect(u, v, t)) {
    //                     e++;
    //                 }
    //             }
    //         }
    //         else {
    //             for(uint32_t j = g->pos[t][u]; j < g->pos[t][u + 1]; j++) {
    //                 uint32_t v = g->e[t][j];
    //                 uint32_t pv = S[z].pos(v);
    //                 if(S[z].c <= pv && pv < S[z].r) {//in C
    //                     e++;
    //                 }
    //             }
    //         }

    //         if(e > maxE[z]) {
    //             isPivotInX[t] = true;
    //             maxE[z] = e;
    //             maxI[t] = i;
    //             maxUV[t] = u;
    //         }
    //     }
    // }

    // //pivot found. pivot must have at least one neighbor in this algorithm
    // t = 1;
    // if(S[0].cSize() - maxE[0] >= S[1].cSize() - maxE[1]) {
    //     t = 0;
    // }
    // z = t ^ 1;


    // //r is more like the right limit of C
    // //r[z] doesnt follow S[z].r
    // //pivot could be in X, and C of z is changed, but never mind, it will recover by r[z]
    // if(g->deg(maxUV[t], t) > S[z].r - S[z].c) {
    //     if(S[z].r>0){
    //         for(uint32_t i = S[z].r - 1; i >= S[z].c; i--) {
    //             if(!g->connect(maxUV[t], S[z][i], t)) {
    //                 S[z].swapByPos(--S[z].r, i);
    //             }
    
    //             if(i == 0) break;
    //         }
    //     }
    // }else{
    //     uint32_t u = maxUV[t];
    //     uint32_t tmpR = S[z].c;
    //     for(uint32_t j = g->pos[t][u]; j < g->pos[t][u + 1]; j++) {
    //         uint32_t v = g->e[t][j];
    //         uint32_t pv = S[z].pos(v);
    //         if(S[z].c <= pv && pv < S[z].r) {
    //             S[z].swapByPos(tmpR++, pv);
    //         }
    //     }
    //     S[z].r = tmpR;
    // }
    // wsSize = r[z] - S[z].r;
    // if(ws[deep].size() < wsSize) {
    //     ws[deep].resize(wsSize * 2);
    // }
    // memcpy(ws[deep].data(), S[z].begin() + S[z].r, sizeof(uint32_t) * wsSize);


    // //iterate on pivot
    // if(isPivotInX[t]==false){
        
    //     S[t].swapByPos(maxI[t], --r[t]);
    //     S[t].r = r[t];

    //     //x[z] doesnt follow S[z].x, used to recover S[z].x
    //     if(g->deg(maxUV[t], t) > S[z].c - S[z].x) {
    //         for(uint32_t i = S[z].x; i < S[z].c; i++) {
    //             if(!g->connect(maxUV[t], S[z][i], t)) {
    //                 S[z].swapByPos(S[z].x++, i);
    //             }
    //         }
    //     }else{
    //         uint32_t u = maxUV[t];
    //         uint32_t tmpX = S[z].c;
    //         for(uint32_t j = g->pos[t][u]; j < g->pos[t][u + 1]; j++) {
    //             uint32_t v = g->e[t][j];
    //             uint32_t pv = S[z].pos(v);
    //             if(S[z].x <= pv && pv < S[z].c) {
    //                 S[z].swapByPos(--tmpX, pv);
    //             }
    //         }
    //         S[z].x = tmpX;
    //     }



    //     RUV[t].push_back(maxUV[t]);
    //     RUVlook[t][maxUV[t]]=1;

    //     //debug
    //     // if(S[1].sz!=g->n[1]){
    //     //     std::cerr<<"error"<<std::endl;
    //     // }
    //     // if(thisNode==401812){
    //     //     std::cerr<<"here"<<std::endl;
    //     // }

    //     // #ifdef PARTIALPIVOT
    //     // int res;
    //     // if(S[z].cSize()+S[t].cSize()<=20){
    //     //     res=bbranch_simplep(deep+1,t,maxUV[t]);
    //     // }else{
    //     //     res=bbranch(deep+1,t,maxUV[t]);
    //     // }
    //     // #else
    //     //C of z must not be empty
    //     int res=bbranchLeafRes(deep+1,t,maxUV[t]);
    //     // #endif

    //     //debug
    //     // if(S[1].sz!=g->n[1]){
    //     //     std::cerr<<"error"<<std::endl;
    //     // }

    //     RUV[t].pop_back();
    //     RUVlook[t][maxUV[t]]=0;
    //     if(res!=-1){
    //         if(thisNode==-1){
    //             thisNode=runtree.getNewNode();
    //             runtree.treeNodePool[thisNode].inValue=inValue;
    //             runtree.treeNodePool[thisNode].uv=side;
    //             runtree.treeNodePool[thisNode].resMark=-1;
    //             runtree.treeNodePool[thisNode].qpMark=false;
    //         }
    //         runtree.treeNodePool[thisNode].sons.push_back(res);
    //         xStackNode.push_back(res);
    //     }

    //     //debug
    //     // if(S[1].sz!=g->n[1]){
    //     //     std::cerr<<"error"<<std::endl;
    //     // }

    //     //put pivot into X
    //     S[t].swapByPos(c[t]++, r[t]++);
    // }




    //iterate on both sides
    int tt=side;
    uint32_t wsize[2];
    for(int o=0;o<=1;++o){
        tt=tt^1;
        int zz=tt^1;
        wsSize=r[tt]-c[tt];
        wsize[tt]=wsSize;
        if(ws2[tt][deep].size() < wsSize) {
            ws2[tt][deep].resize(wsSize * 2);
        }
        memcpy(ws2[tt][deep].data(), S[tt].begin() + c[tt], sizeof(uint32_t) * wsSize);

        for(uint32_t j = 0; j < wsSize; j++){
            uint32_t w = ws2[tt][deep][j];
            S[tt].swapByPos(S[tt].pos(w), --r[tt]);

            //recover next level parameter
            S[0].x = x[0]; S[1].x = x[1];
            S[0].c = c[0]; S[1].c = c[1];
            S[0].r = r[0]; S[1].r = r[1];

            if(g->deg(w, tt) > S[zz].c - S[zz].x) {
                for(uint32_t i = S[zz].x; i < S[zz].c; i++) {
                    if(!g->connect(w, S[zz][i], tt)) {
                        S[zz].swapByPos(S[zz].x++, i);
                    }
                }
            }
            else {
                uint32_t tmpX = S[zz].c;
                for(uint32_t i = g->pos[tt][w]; i < g->pos[tt][w + 1]; i++) {
                    uint32_t u = g->e[tt][i];
                    uint32_t pu = S[zz].pos(u);
                    if(S[zz].x <= pu && pu < S[zz].c) {
                        S[zz].swapByPos(--tmpX, pu);
                    }
                }
                S[zz].x = tmpX;
            }

            if(g->deg(w, tt) > S[zz].r - S[zz].c) {
                if(S[zz].r > 0){
                    for(uint32_t i = S[zz].r - 1; i >= S[zz].c; i--) {
                        if(!g->connect(w, S[zz][i], tt)) {
                            S[zz].swapByPos(--S[zz].r, i);
                        }

                        if(i == 0) break;
                    }
                }
            }else{
                uint32_t tmpR = S[zz].c;
                for(uint32_t j = g->pos[tt][w]; j < g->pos[tt][w + 1]; j++) {
                    uint32_t u = g->e[tt][j];
                    uint32_t pu = S[zz].pos(u);
                    if(S[zz].c <= pu && pu < S[zz].r) {
                        S[zz].swapByPos(tmpR++, pu);
                    }
                }
                S[zz].r = tmpR;
            }

            // //debug
            // if(deep==0&&inValue==1859991&&w==1911392){
            //     // std::cerr<<<<std::endl;
            //     for(int ii=S[z].c;ii<S[z].r;++ii){
            //         std::cerr<<g->old_lables[z][S[z][ii]]<<" ";
            //     }
            //     std::cerr<<std::endl;
                
            // }

            // if(S[zz].cSize()>0){
            RUV[tt].push_back(w);
            RUVlook[tt][w]=1;
            //debug
            // if(S[1].sz!=g->n[1]){
            //     std::cerr<<"error"<<std::endl;
            // }

            int res=bbranchLeafRespoorcntpqnode_simplep(deep+1,tt,w);
            RUV[tt].pop_back();
            RUVlook[tt][w]=0;

            if(res!=-1){
                if(thisNode==-1){
                    thisNode=runtree.getNewNode();
                    runtree.treeNodePool[thisNode].inValue=inValue;
                    runtree.treeNodePool[thisNode].uv=side;
                    runtree.treeNodePool[thisNode].resMark=-1;
                    runtree.treeNodePool[thisNode].qpMark=false;
                }
                runtree.treeNodePool[thisNode].sons.push_back(res);
                xStackNode.push_back(res);
            }
            // }

            S[tt].swapByPos(c[tt]++,r[tt]++);//move w to X
        }
    }

    for(int o=0;o<=1;++o){
        for(uint32_t j = 0; j < wsize[o]; j++) {
            uint32_t w = ws2[o][deep][j];
            S[o].swapByPos(S[o].pos(w), --c[o]);
        }
    }

    //recover x,c,r
    S[t].x = x[t]; S[z].x = x[z];
    S[t].c = c[t]; S[z].c = c[z];
    S[t].r = r[t]; S[z].r = r[z];

    //connect that should be connected
    //almost the same as before
    if(RUV[t].size()>=threshold[t]&&RUV[z].size()>=threshold[z]){
        int res=-1;
        // if(RUV[t].size()==threshold[t]&&t==oriSide||RUV[z].size()==threshold[z]&&z==oriSide){// only the first one
        //     res=connectAllBiCliqueLeafResbase();
        //     if(thisNode!=-1){
        //         runtree.treeNodePool[thisNode].resMark=res;
        //     }
        // }else{
        //     if(thisNode!=-1){
        //         fans.updateSize(cliqueuv[0].size());
        //         for(int i=0;i<runtree.treeNodePool[thisNode].sons.size();++i){
        //             int sonNode=runtree.treeNodePool[thisNode].sons[i];
        //             connectResults.push_back(runtree.treeNodePool[sonNode].resMark);
        //         }

        //         for(int i=1;i<connectResults.size();++i){
        //             fans.merge(connectResults[i],connectResults[0]);
        //         }
        //         runtree.treeNodePool[thisNode].resMark=fans.find(connectResults[0]);
        //         connectResults.clear();
        //     }

        // }

        if(thisNode!=-1){
            fans.updateSize(cliqueuv[0].size());
            for(int i=0;i<runtree.treeNodePool[thisNode].sons.size();++i){
                int sonNode=runtree.treeNodePool[thisNode].sons[i];
                connectResults.push_back(runtree.treeNodePool[sonNode].resMark);
            }

            for(int i=1;i<connectResults.size();++i){
                fans.merge(connectResults[i],connectResults[0]);
            }
            runtree.treeNodePool[thisNode].resMark=fans.find(connectResults[0]);
            connectResults.clear();

            //check pqnode
            if(RUV[t].size()==threshold[t]&&t==oriSide||RUV[z].size()==threshold[z]&&z==oriSide){
                pqnodecnt++;
            }
        }


    }

    xStackNode.resize(oriXstackSize);
    return thisNode;
}

int BICPC::bbranchLeafRespoor_simplep(uint32_t deep,uint32_t side,uint32_t inValue){
    // if(S[1].sz!=g->n[1]){
    //     std::cerr<<"error"<<std::endl;
    // }
    // if(S[z].CIsEmpty()) return;

    int thisNode=-1;//temporally no tree node is built here
    int oriSide=side;

    uint32_t t=side;
    uint32_t z=side ^ 1;

    if(S[t].CIsEmpty()&&S[z].CIsEmpty()){
        if(S[t].XIsEmpty()&&S[z].XIsEmpty()){
            if(RUV[t].size()>=threshold[t]&&RUV[z].size()>=threshold[z]&&RUV[t].size()+RUV[z].size()>threshold[t]+threshold[z]){
                thisNode=runtree.getNewNode();
                runtree.treeNodePool[thisNode].inValue=inValue;
                runtree.treeNodePool[thisNode].uv=side;
                runtree.treeNodePool[thisNode].resMark=-1;
                runtree.treeNodePool[thisNode].qpMark=false;

                cliqueuv[t].emplace_back(RUV[t]);
                cliqueuv[z].emplace_back(RUV[z]);

                runtree.treeNodePool[thisNode].resMark=cliqueuv[t].size()-1;
            }
            return thisNode;
        }
        // if(RUV[t].size()>=threshold[t]&&RUV[z].size()>=threshold[z]){
        //     if(RUV[t].size()==threshold[t]&&t==oriSide||RUV[z].size()==threshold[z]&&z==oriSide){
        //         connectAllBiCliqueLeafResbase();
        //     }
        // }
        return -1;
    }

    if(RUV[t].size()+S[t].cSize()<threshold[t]||RUV[z].size()+S[z].cSize()<threshold[z]){
        return thisNode;
    }

    //find pivot

    //debug
    // if(S[1].sz!=g->n[1]){
    //     std::cerr<<"error"<<std::endl;
    // }

    //initialize copies of xcr
    uint32_t x[2];
    uint32_t c[2];
    uint32_t r[2];
    x[t] = S[t].x; x[z] = S[z].x;
    c[t] = S[t].c; c[z] = S[z].c;
    r[t] = S[t].r; r[z] = S[z].r;

    int oriXstackSize=xStackNode.size();

    //find pivot
    // uint32_t maxI[2] = {0}, maxE[2] = {0}, maxUV[2] = {0};
    // maxI[0] = S[0].c;
    // maxI[1] = S[1].c;
    // maxUV[0] = S[0][S[0].c];
    // maxUV[1] = S[1][S[1].c];
    // t=0,z=1;
    uint32_t wsSize;


    //find pivot in C of t
    // for(uint32_t i = S[t].c; i < S[t].r; i++) {//scan all vertices in C
    // // for(uint32_t i = S[t].c; i <= S[t].c; i++) {
    //     uint32_t u = S[t][i];
    //     uint32_t e = 0;

    //     if(g->deg(u, t) > S[z].r - S[z].c) {
    //         for(uint32_t j = S[z].c; j < S[z].r; j++) {
    //             uint32_t v = S[z][j];
    //             if(g->connect(u, v, t)) {
    //                 e++;
    //                 deg[v]++;
    //             }
    //         }
    //     }
    //     else {
    //         for(uint32_t j = g->pos[t][u]; j < g->pos[t][u + 1]; j++) {
    //             uint32_t v = g->e[t][j];
    //             uint32_t pv = S[z].pos(v);
    //             if(S[z].c <= pv && pv < S[z].r) {//in C
    //                 e++;
    //                 deg[v]++;
    //             }
    //         }
    //     }

    //     if(e > maxE[z]) {
    //         maxE[z] = e;
    //         maxI[t] = i;
    //         maxUV[t] = u;
    //     }
    // }
    // //find pivot in C of z
    // for(uint32_t i = S[1].c; i < S[1].r; i++) {//scan all vertices in C
    // // for(uint32_t i = S[t].c; i <= S[t].c; i++) {
    //     uint32_t u = S[1][i];
    //     uint32_t e = deg[u];
    //     deg[u] = 0;
    //     if(e > maxE[0]) {
    //         maxE[0] = e;
    //         maxI[1] = i;
    //         maxUV[1] = u;
    //     }
    // }

    // //check if pivot in X of both sides
    // bool isPivotInX[2] = {false, false};
    // for(t = 0; t < 2; t++) {
    //     z = t ^ 1;//t + z = 1

    //     for(uint32_t i = S[t].x; i < S[t].c; i++) {//scan all vertices in X
    //     // for(uint32_t i = S[t].c; i <= S[t].c; i++) {
    //         uint32_t u = S[t][i];
    //         uint32_t e = 0;

    //         if(g->deg(u, t) > S[z].r - S[z].c) {
    //             for(uint32_t j = S[z].c; j < S[z].r; j++) {
    //                 uint32_t v = S[z][j];
    //                 if(g->connect(u, v, t)) {
    //                     e++;
    //                 }
    //             }
    //         }
    //         else {
    //             for(uint32_t j = g->pos[t][u]; j < g->pos[t][u + 1]; j++) {
    //                 uint32_t v = g->e[t][j];
    //                 uint32_t pv = S[z].pos(v);
    //                 if(S[z].c <= pv && pv < S[z].r) {//in C
    //                     e++;
    //                 }
    //             }
    //         }

    //         if(e > maxE[z]) {
    //             isPivotInX[t] = true;
    //             maxE[z] = e;
    //             maxI[t] = i;
    //             maxUV[t] = u;
    //         }
    //     }
    // }

    // //pivot found. pivot must have at least one neighbor in this algorithm
    // t = 1;
    // if(S[0].cSize() - maxE[0] >= S[1].cSize() - maxE[1]) {
    //     t = 0;
    // }
    // z = t ^ 1;


    // //r is more like the right limit of C
    // //r[z] doesnt follow S[z].r
    // //pivot could be in X, and C of z is changed, but never mind, it will recover by r[z]
    // if(g->deg(maxUV[t], t) > S[z].r - S[z].c) {
    //     if(S[z].r>0){
    //         for(uint32_t i = S[z].r - 1; i >= S[z].c; i--) {
    //             if(!g->connect(maxUV[t], S[z][i], t)) {
    //                 S[z].swapByPos(--S[z].r, i);
    //             }
    
    //             if(i == 0) break;
    //         }
    //     }
    // }else{
    //     uint32_t u = maxUV[t];
    //     uint32_t tmpR = S[z].c;
    //     for(uint32_t j = g->pos[t][u]; j < g->pos[t][u + 1]; j++) {
    //         uint32_t v = g->e[t][j];
    //         uint32_t pv = S[z].pos(v);
    //         if(S[z].c <= pv && pv < S[z].r) {
    //             S[z].swapByPos(tmpR++, pv);
    //         }
    //     }
    //     S[z].r = tmpR;
    // }
    // wsSize = r[z] - S[z].r;
    // if(ws[deep].size() < wsSize) {
    //     ws[deep].resize(wsSize * 2);
    // }
    // memcpy(ws[deep].data(), S[z].begin() + S[z].r, sizeof(uint32_t) * wsSize);


    // //iterate on pivot
    // if(isPivotInX[t]==false){
        
    //     S[t].swapByPos(maxI[t], --r[t]);
    //     S[t].r = r[t];

    //     //x[z] doesnt follow S[z].x, used to recover S[z].x
    //     if(g->deg(maxUV[t], t) > S[z].c - S[z].x) {
    //         for(uint32_t i = S[z].x; i < S[z].c; i++) {
    //             if(!g->connect(maxUV[t], S[z][i], t)) {
    //                 S[z].swapByPos(S[z].x++, i);
    //             }
    //         }
    //     }else{
    //         uint32_t u = maxUV[t];
    //         uint32_t tmpX = S[z].c;
    //         for(uint32_t j = g->pos[t][u]; j < g->pos[t][u + 1]; j++) {
    //             uint32_t v = g->e[t][j];
    //             uint32_t pv = S[z].pos(v);
    //             if(S[z].x <= pv && pv < S[z].c) {
    //                 S[z].swapByPos(--tmpX, pv);
    //             }
    //         }
    //         S[z].x = tmpX;
    //     }



    //     RUV[t].push_back(maxUV[t]);
    //     RUVlook[t][maxUV[t]]=1;

    //     //debug
    //     // if(S[1].sz!=g->n[1]){
    //     //     std::cerr<<"error"<<std::endl;
    //     // }
    //     // if(thisNode==401812){
    //     //     std::cerr<<"here"<<std::endl;
    //     // }

    //     // #ifdef PARTIALPIVOT
    //     // int res;
    //     // if(S[z].cSize()+S[t].cSize()<=20){
    //     //     res=bbranch_simplep(deep+1,t,maxUV[t]);
    //     // }else{
    //     //     res=bbranch(deep+1,t,maxUV[t]);
    //     // }
    //     // #else
    //     //C of z must not be empty
    //     int res=bbranchLeafRes(deep+1,t,maxUV[t]);
    //     // #endif

    //     //debug
    //     // if(S[1].sz!=g->n[1]){
    //     //     std::cerr<<"error"<<std::endl;
    //     // }

    //     RUV[t].pop_back();
    //     RUVlook[t][maxUV[t]]=0;
    //     if(res!=-1){
    //         if(thisNode==-1){
    //             thisNode=runtree.getNewNode();
    //             runtree.treeNodePool[thisNode].inValue=inValue;
    //             runtree.treeNodePool[thisNode].uv=side;
    //             runtree.treeNodePool[thisNode].resMark=-1;
    //             runtree.treeNodePool[thisNode].qpMark=false;
    //         }
    //         runtree.treeNodePool[thisNode].sons.push_back(res);
    //         xStackNode.push_back(res);
    //     }

    //     //debug
    //     // if(S[1].sz!=g->n[1]){
    //     //     std::cerr<<"error"<<std::endl;
    //     // }

    //     //put pivot into X
    //     S[t].swapByPos(c[t]++, r[t]++);
    // }




    //iterate on both sides
    int tt=side;
    uint32_t wsize[2];
    for(int o=0;o<=1;++o){
        tt=tt^1;
        int zz=tt^1;
        wsSize=r[tt]-c[tt];
        wsize[tt]=wsSize;
        if(ws2[tt][deep].size() < wsSize) {
            ws2[tt][deep].resize(wsSize * 2);
        }
        memcpy(ws2[tt][deep].data(), S[tt].begin() + c[tt], sizeof(uint32_t) * wsSize);

        for(uint32_t j = 0; j < wsSize; j++){
            uint32_t w = ws2[tt][deep][j];
            S[tt].swapByPos(S[tt].pos(w), --r[tt]);

            //recover next level parameter
            S[0].x = x[0]; S[1].x = x[1];
            S[0].c = c[0]; S[1].c = c[1];
            S[0].r = r[0]; S[1].r = r[1];

            if(g->deg(w, tt) > S[zz].c - S[zz].x) {
                for(uint32_t i = S[zz].x; i < S[zz].c; i++) {
                    if(!g->connect(w, S[zz][i], tt)) {
                        S[zz].swapByPos(S[zz].x++, i);
                    }
                }
            }
            else {
                uint32_t tmpX = S[zz].c;
                for(uint32_t i = g->pos[tt][w]; i < g->pos[tt][w + 1]; i++) {
                    uint32_t u = g->e[tt][i];
                    uint32_t pu = S[zz].pos(u);
                    if(S[zz].x <= pu && pu < S[zz].c) {
                        S[zz].swapByPos(--tmpX, pu);
                    }
                }
                S[zz].x = tmpX;
            }

            if(g->deg(w, tt) > S[zz].r - S[zz].c) {
                if(S[zz].r > 0){
                    for(uint32_t i = S[zz].r - 1; i >= S[zz].c; i--) {
                        if(!g->connect(w, S[zz][i], tt)) {
                            S[zz].swapByPos(--S[zz].r, i);
                        }

                        if(i == 0) break;
                    }
                }
            }else{
                uint32_t tmpR = S[zz].c;
                for(uint32_t j = g->pos[tt][w]; j < g->pos[tt][w + 1]; j++) {
                    uint32_t u = g->e[tt][j];
                    uint32_t pu = S[zz].pos(u);
                    if(S[zz].c <= pu && pu < S[zz].r) {
                        S[zz].swapByPos(tmpR++, pu);
                    }
                }
                S[zz].r = tmpR;
            }

            // //debug
            // if(deep==0&&inValue==1859991&&w==1911392){
            //     // std::cerr<<<<std::endl;
            //     for(int ii=S[z].c;ii<S[z].r;++ii){
            //         std::cerr<<g->old_lables[z][S[z][ii]]<<" ";
            //     }
            //     std::cerr<<std::endl;
                
            // }

            // if(S[zz].cSize()>0){
            RUV[tt].push_back(w);
            RUVlook[tt][w]=1;
            //debug
            // if(S[1].sz!=g->n[1]){
            //     std::cerr<<"error"<<std::endl;
            // }

            int res=bbranchLeafRespoor_simplep(deep+1,tt,w);
            RUV[tt].pop_back();
            RUVlook[tt][w]=0;

            if(res!=-1){
                if(thisNode==-1){
                    thisNode=runtree.getNewNode();
                    runtree.treeNodePool[thisNode].inValue=inValue;
                    runtree.treeNodePool[thisNode].uv=side;
                    runtree.treeNodePool[thisNode].resMark=-1;
                    runtree.treeNodePool[thisNode].qpMark=false;
                }
                runtree.treeNodePool[thisNode].sons.push_back(res);
                xStackNode.push_back(res);
            }
            // }

            S[tt].swapByPos(c[tt]++,r[tt]++);//move w to X
        }
    }

    for(int o=0;o<=1;++o){
        for(uint32_t j = 0; j < wsize[o]; j++) {
            uint32_t w = ws2[o][deep][j];
            S[o].swapByPos(S[o].pos(w), --c[o]);
        }
    }

    //recover x,c,r
    S[t].x = x[t]; S[z].x = x[z];
    S[t].c = c[t]; S[z].c = c[z];
    S[t].r = r[t]; S[z].r = r[z];

    //connect that should be connected
    //almost the same as before
    if(RUV[t].size()>=threshold[t]&&RUV[z].size()>=threshold[z]){
        int res=-1;
        // if(RUV[t].size()==threshold[t]&&t==oriSide||RUV[z].size()==threshold[z]&&z==oriSide){// only the first one
        //     res=connectAllBiCliqueLeafResbase();
        //     if(thisNode!=-1){
        //         runtree.treeNodePool[thisNode].resMark=res;
        //     }
        // }else{
        //     if(thisNode!=-1){
        //         fans.updateSize(cliqueuv[0].size());
        //         for(int i=0;i<runtree.treeNodePool[thisNode].sons.size();++i){
        //             int sonNode=runtree.treeNodePool[thisNode].sons[i];
        //             connectResults.push_back(runtree.treeNodePool[sonNode].resMark);
        //         }

        //         for(int i=1;i<connectResults.size();++i){
        //             fans.merge(connectResults[i],connectResults[0]);
        //         }
        //         runtree.treeNodePool[thisNode].resMark=fans.find(connectResults[0]);
        //         connectResults.clear();
        //     }

        // }

        if(thisNode!=-1){
            fans.updateSize(cliqueuv[0].size());
            for(int i=0;i<runtree.treeNodePool[thisNode].sons.size();++i){
                int sonNode=runtree.treeNodePool[thisNode].sons[i];
                connectResults.push_back(runtree.treeNodePool[sonNode].resMark);
            }

            for(int i=1;i<connectResults.size();++i){
                fans.merge(connectResults[i],connectResults[0]);
            }
            runtree.treeNodePool[thisNode].resMark=fans.find(connectResults[0]);
            connectResults.clear();
        }


    }

    xStackNode.resize(oriXstackSize);
    return thisNode;
}

int BICPC::bbranchLeafResbase_simplep(uint32_t deep,uint32_t side,uint32_t inValue){
    // if(S[1].sz!=g->n[1]){
    //     std::cerr<<"error"<<std::endl;
    // }
    // if(S[z].CIsEmpty()) return;

    int thisNode=-1;//temporally no tree node is built here
    int oriSide=side;

    uint32_t t=side;
    uint32_t z=side ^ 1;

    if(S[t].CIsEmpty()&&S[z].CIsEmpty()){
        if(S[t].XIsEmpty()&&S[z].XIsEmpty()){
            if(RUV[t].size()>=threshold[t]&&RUV[z].size()>=threshold[z]&&RUV[t].size()+RUV[z].size()>threshold[t]+threshold[z]){
                thisNode=runtree.getNewNode();
                runtree.treeNodePool[thisNode].inValue=inValue;
                runtree.treeNodePool[thisNode].uv=side;
                runtree.treeNodePool[thisNode].resMark=-1;
                runtree.treeNodePool[thisNode].qpMark=false;

                cliqueuv[t].emplace_back(RUV[t]);
                cliqueuv[z].emplace_back(RUV[z]);

                runtree.treeNodePool[thisNode].resMark=cliqueuv[t].size()-1;
            }
            return thisNode;
        }
        if(RUV[t].size()>=threshold[t]&&RUV[z].size()>=threshold[z]){
            if(RUV[t].size()==threshold[t]&&t==oriSide||RUV[z].size()==threshold[z]&&z==oriSide){
                connectAllBiCliqueLeafResbase();
            }
        }
        return -1;
    }

    if(RUV[t].size()+S[t].cSize()<threshold[t]||RUV[z].size()+S[z].cSize()<threshold[z]){
        return thisNode;
    }

    //find pivot

    //debug
    // if(S[1].sz!=g->n[1]){
    //     std::cerr<<"error"<<std::endl;
    // }

    //initialize copies of xcr
    uint32_t x[2];
    uint32_t c[2];
    uint32_t r[2];
    x[t] = S[t].x; x[z] = S[z].x;
    c[t] = S[t].c; c[z] = S[z].c;
    r[t] = S[t].r; r[z] = S[z].r;

    int oriXstackSize=xStackNode.size();

    //find pivot
    // uint32_t maxI[2] = {0}, maxE[2] = {0}, maxUV[2] = {0};
    // maxI[0] = S[0].c;
    // maxI[1] = S[1].c;
    // maxUV[0] = S[0][S[0].c];
    // maxUV[1] = S[1][S[1].c];
    // t=0,z=1;
    uint32_t wsSize;


    //find pivot in C of t
    // for(uint32_t i = S[t].c; i < S[t].r; i++) {//scan all vertices in C
    // // for(uint32_t i = S[t].c; i <= S[t].c; i++) {
    //     uint32_t u = S[t][i];
    //     uint32_t e = 0;

    //     if(g->deg(u, t) > S[z].r - S[z].c) {
    //         for(uint32_t j = S[z].c; j < S[z].r; j++) {
    //             uint32_t v = S[z][j];
    //             if(g->connect(u, v, t)) {
    //                 e++;
    //                 deg[v]++;
    //             }
    //         }
    //     }
    //     else {
    //         for(uint32_t j = g->pos[t][u]; j < g->pos[t][u + 1]; j++) {
    //             uint32_t v = g->e[t][j];
    //             uint32_t pv = S[z].pos(v);
    //             if(S[z].c <= pv && pv < S[z].r) {//in C
    //                 e++;
    //                 deg[v]++;
    //             }
    //         }
    //     }

    //     if(e > maxE[z]) {
    //         maxE[z] = e;
    //         maxI[t] = i;
    //         maxUV[t] = u;
    //     }
    // }
    // //find pivot in C of z
    // for(uint32_t i = S[1].c; i < S[1].r; i++) {//scan all vertices in C
    // // for(uint32_t i = S[t].c; i <= S[t].c; i++) {
    //     uint32_t u = S[1][i];
    //     uint32_t e = deg[u];
    //     deg[u] = 0;
    //     if(e > maxE[0]) {
    //         maxE[0] = e;
    //         maxI[1] = i;
    //         maxUV[1] = u;
    //     }
    // }

    // //check if pivot in X of both sides
    // bool isPivotInX[2] = {false, false};
    // for(t = 0; t < 2; t++) {
    //     z = t ^ 1;//t + z = 1

    //     for(uint32_t i = S[t].x; i < S[t].c; i++) {//scan all vertices in X
    //     // for(uint32_t i = S[t].c; i <= S[t].c; i++) {
    //         uint32_t u = S[t][i];
    //         uint32_t e = 0;

    //         if(g->deg(u, t) > S[z].r - S[z].c) {
    //             for(uint32_t j = S[z].c; j < S[z].r; j++) {
    //                 uint32_t v = S[z][j];
    //                 if(g->connect(u, v, t)) {
    //                     e++;
    //                 }
    //             }
    //         }
    //         else {
    //             for(uint32_t j = g->pos[t][u]; j < g->pos[t][u + 1]; j++) {
    //                 uint32_t v = g->e[t][j];
    //                 uint32_t pv = S[z].pos(v);
    //                 if(S[z].c <= pv && pv < S[z].r) {//in C
    //                     e++;
    //                 }
    //             }
    //         }

    //         if(e > maxE[z]) {
    //             isPivotInX[t] = true;
    //             maxE[z] = e;
    //             maxI[t] = i;
    //             maxUV[t] = u;
    //         }
    //     }
    // }

    // //pivot found. pivot must have at least one neighbor in this algorithm
    // t = 1;
    // if(S[0].cSize() - maxE[0] >= S[1].cSize() - maxE[1]) {
    //     t = 0;
    // }
    // z = t ^ 1;


    // //r is more like the right limit of C
    // //r[z] doesnt follow S[z].r
    // //pivot could be in X, and C of z is changed, but never mind, it will recover by r[z]
    // if(g->deg(maxUV[t], t) > S[z].r - S[z].c) {
    //     if(S[z].r>0){
    //         for(uint32_t i = S[z].r - 1; i >= S[z].c; i--) {
    //             if(!g->connect(maxUV[t], S[z][i], t)) {
    //                 S[z].swapByPos(--S[z].r, i);
    //             }
    
    //             if(i == 0) break;
    //         }
    //     }
    // }else{
    //     uint32_t u = maxUV[t];
    //     uint32_t tmpR = S[z].c;
    //     for(uint32_t j = g->pos[t][u]; j < g->pos[t][u + 1]; j++) {
    //         uint32_t v = g->e[t][j];
    //         uint32_t pv = S[z].pos(v);
    //         if(S[z].c <= pv && pv < S[z].r) {
    //             S[z].swapByPos(tmpR++, pv);
    //         }
    //     }
    //     S[z].r = tmpR;
    // }
    // wsSize = r[z] - S[z].r;
    // if(ws[deep].size() < wsSize) {
    //     ws[deep].resize(wsSize * 2);
    // }
    // memcpy(ws[deep].data(), S[z].begin() + S[z].r, sizeof(uint32_t) * wsSize);


    // //iterate on pivot
    // if(isPivotInX[t]==false){
        
    //     S[t].swapByPos(maxI[t], --r[t]);
    //     S[t].r = r[t];

    //     //x[z] doesnt follow S[z].x, used to recover S[z].x
    //     if(g->deg(maxUV[t], t) > S[z].c - S[z].x) {
    //         for(uint32_t i = S[z].x; i < S[z].c; i++) {
    //             if(!g->connect(maxUV[t], S[z][i], t)) {
    //                 S[z].swapByPos(S[z].x++, i);
    //             }
    //         }
    //     }else{
    //         uint32_t u = maxUV[t];
    //         uint32_t tmpX = S[z].c;
    //         for(uint32_t j = g->pos[t][u]; j < g->pos[t][u + 1]; j++) {
    //             uint32_t v = g->e[t][j];
    //             uint32_t pv = S[z].pos(v);
    //             if(S[z].x <= pv && pv < S[z].c) {
    //                 S[z].swapByPos(--tmpX, pv);
    //             }
    //         }
    //         S[z].x = tmpX;
    //     }



    //     RUV[t].push_back(maxUV[t]);
    //     RUVlook[t][maxUV[t]]=1;

    //     //debug
    //     // if(S[1].sz!=g->n[1]){
    //     //     std::cerr<<"error"<<std::endl;
    //     // }
    //     // if(thisNode==401812){
    //     //     std::cerr<<"here"<<std::endl;
    //     // }

    //     // #ifdef PARTIALPIVOT
    //     // int res;
    //     // if(S[z].cSize()+S[t].cSize()<=20){
    //     //     res=bbranch_simplep(deep+1,t,maxUV[t]);
    //     // }else{
    //     //     res=bbranch(deep+1,t,maxUV[t]);
    //     // }
    //     // #else
    //     //C of z must not be empty
    //     int res=bbranchLeafRes(deep+1,t,maxUV[t]);
    //     // #endif

    //     //debug
    //     // if(S[1].sz!=g->n[1]){
    //     //     std::cerr<<"error"<<std::endl;
    //     // }

    //     RUV[t].pop_back();
    //     RUVlook[t][maxUV[t]]=0;
    //     if(res!=-1){
    //         if(thisNode==-1){
    //             thisNode=runtree.getNewNode();
    //             runtree.treeNodePool[thisNode].inValue=inValue;
    //             runtree.treeNodePool[thisNode].uv=side;
    //             runtree.treeNodePool[thisNode].resMark=-1;
    //             runtree.treeNodePool[thisNode].qpMark=false;
    //         }
    //         runtree.treeNodePool[thisNode].sons.push_back(res);
    //         xStackNode.push_back(res);
    //     }

    //     //debug
    //     // if(S[1].sz!=g->n[1]){
    //     //     std::cerr<<"error"<<std::endl;
    //     // }

    //     //put pivot into X
    //     S[t].swapByPos(c[t]++, r[t]++);
    // }




    //iterate on both sides
    int tt=side;
    uint32_t wsize[2];
    for(int o=0;o<=1;++o){
        tt=tt^1;
        int zz=tt^1;
        wsSize=r[tt]-c[tt];
        wsize[tt]=wsSize;
        if(ws2[tt][deep].size() < wsSize) {
            ws2[tt][deep].resize(wsSize * 2);
        }
        memcpy(ws2[tt][deep].data(), S[tt].begin() + c[tt], sizeof(uint32_t) * wsSize);

        for(uint32_t j = 0; j < wsSize; j++){
            uint32_t w = ws2[tt][deep][j];
            S[tt].swapByPos(S[tt].pos(w), --r[tt]);

            //recover next level parameter
            S[0].x = x[0]; S[1].x = x[1];
            S[0].c = c[0]; S[1].c = c[1];
            S[0].r = r[0]; S[1].r = r[1];

            if(g->deg(w, tt) > S[zz].c - S[zz].x) {
                for(uint32_t i = S[zz].x; i < S[zz].c; i++) {
                    if(!g->connect(w, S[zz][i], tt)) {
                        S[zz].swapByPos(S[zz].x++, i);
                    }
                }
            }
            else {
                uint32_t tmpX = S[zz].c;
                for(uint32_t i = g->pos[tt][w]; i < g->pos[tt][w + 1]; i++) {
                    uint32_t u = g->e[tt][i];
                    uint32_t pu = S[zz].pos(u);
                    if(S[zz].x <= pu && pu < S[zz].c) {
                        S[zz].swapByPos(--tmpX, pu);
                    }
                }
                S[zz].x = tmpX;
            }

            if(g->deg(w, tt) > S[zz].r - S[zz].c) {
                if(S[zz].r > 0){
                    for(uint32_t i = S[zz].r - 1; i >= S[zz].c; i--) {
                        if(!g->connect(w, S[zz][i], tt)) {
                            S[zz].swapByPos(--S[zz].r, i);
                        }

                        if(i == 0) break;
                    }
                }
            }else{
                uint32_t tmpR = S[zz].c;
                for(uint32_t j = g->pos[tt][w]; j < g->pos[tt][w + 1]; j++) {
                    uint32_t u = g->e[tt][j];
                    uint32_t pu = S[zz].pos(u);
                    if(S[zz].c <= pu && pu < S[zz].r) {
                        S[zz].swapByPos(tmpR++, pu);
                    }
                }
                S[zz].r = tmpR;
            }

            // //debug
            // if(deep==0&&inValue==1859991&&w==1911392){
            //     // std::cerr<<<<std::endl;
            //     for(int ii=S[z].c;ii<S[z].r;++ii){
            //         std::cerr<<g->old_lables[z][S[z][ii]]<<" ";
            //     }
            //     std::cerr<<std::endl;
                
            // }

            // if(S[zz].cSize()>0){
            RUV[tt].push_back(w);
            RUVlook[tt][w]=1;
            //debug
            // if(S[1].sz!=g->n[1]){
            //     std::cerr<<"error"<<std::endl;
            // }

            int res=bbranchLeafResbase_simplep(deep+1,tt,w);
            RUV[tt].pop_back();
            RUVlook[tt][w]=0;

            if(res!=-1){
                if(thisNode==-1){
                    thisNode=runtree.getNewNode();
                    runtree.treeNodePool[thisNode].inValue=inValue;
                    runtree.treeNodePool[thisNode].uv=side;
                    runtree.treeNodePool[thisNode].resMark=-1;
                    runtree.treeNodePool[thisNode].qpMark=false;
                }
                runtree.treeNodePool[thisNode].sons.push_back(res);
                xStackNode.push_back(res);
            }
            // }

            S[tt].swapByPos(c[tt]++,r[tt]++);//move w to X
        }
    }

    for(int o=0;o<=1;++o){
        for(uint32_t j = 0; j < wsize[o]; j++) {
            uint32_t w = ws2[o][deep][j];
            S[o].swapByPos(S[o].pos(w), --c[o]);
        }
    }

    //recover x,c,r
    S[t].x = x[t]; S[z].x = x[z];
    S[t].c = c[t]; S[z].c = c[z];
    S[t].r = r[t]; S[z].r = r[z];

    //connect that should be connected
    //almost the same as before
    if(RUV[t].size()>=threshold[t]&&RUV[z].size()>=threshold[z]){
        int res=-1;
        if(RUV[t].size()==threshold[t]&&t==oriSide||RUV[z].size()==threshold[z]&&z==oriSide){// only the first one
            res=connectAllBiCliqueLeafResbase();
            if(thisNode!=-1){
                runtree.treeNodePool[thisNode].resMark=res;
            }
        }else{
            if(thisNode!=-1){
                fans.updateSize(cliqueuv[0].size());
                for(int i=0;i<runtree.treeNodePool[thisNode].sons.size();++i){
                    int sonNode=runtree.treeNodePool[thisNode].sons[i];
                    connectResults.push_back(runtree.treeNodePool[sonNode].resMark);
                }

                for(int i=1;i<connectResults.size();++i){
                    fans.merge(connectResults[i],connectResults[0]);
                }
                runtree.treeNodePool[thisNode].resMark=fans.find(connectResults[0]);
                connectResults.clear();
            }

        }


    }

    xStackNode.resize(oriXstackSize);
    return thisNode;
}

int BICPC::bbranchLeafRes_simplep(uint32_t deep,uint32_t side,uint32_t inValue){
    // if(S[1].sz!=g->n[1]){
    //     std::cerr<<"error"<<std::endl;
    // }
    // if(S[z].CIsEmpty()) return;

    int thisNode=-1;//temporally no tree node is built here
    int oriSide=side;

    uint32_t t=side;
    uint32_t z=side ^ 1;

    if(S[t].CIsEmpty()&&S[z].CIsEmpty()){
        if(S[t].XIsEmpty()&&S[z].XIsEmpty()){
            if(RUV[t].size()>=threshold[t]&&RUV[z].size()>=threshold[z]&&RUV[t].size()+RUV[z].size()>threshold[t]+threshold[z]){
                thisNode=runtree.getNewNode();
                runtree.treeNodePool[thisNode].inValue=inValue;
                runtree.treeNodePool[thisNode].uv=side;
                runtree.treeNodePool[thisNode].resMark=-1;
                runtree.treeNodePool[thisNode].qpMark=false;

                cliqueuv[t].emplace_back(RUV[t]);
                cliqueuv[z].emplace_back(RUV[z]);

                runtree.treeNodePool[thisNode].resMark=cliqueuv[t].size()-1;
            }
            return thisNode;
        }
        if(RUV[t].size()>=threshold[t]&&RUV[z].size()>=threshold[z]){
            if(RUV[t].size()==threshold[t]&&t==oriSide||RUV[z].size()==threshold[z]&&z==oriSide){
                connectAllBiCliqueLeafRes();
            }
        }
        return -1;
    }

    if(RUV[t].size()+S[t].cSize()<threshold[t]||RUV[z].size()+S[z].cSize()<threshold[z]){
        return thisNode;
    }

    //find pivot

    //debug
    // if(S[1].sz!=g->n[1]){
    //     std::cerr<<"error"<<std::endl;
    // }

    //initialize copies of xcr
    uint32_t x[2];
    uint32_t c[2];
    uint32_t r[2];
    x[t] = S[t].x; x[z] = S[z].x;
    c[t] = S[t].c; c[z] = S[z].c;
    r[t] = S[t].r; r[z] = S[z].r;

    int oriXstackSize=xStackNode.size();

    //find pivot
    // uint32_t maxI[2] = {0}, maxE[2] = {0}, maxUV[2] = {0};
    // maxI[0] = S[0].c;
    // maxI[1] = S[1].c;
    // maxUV[0] = S[0][S[0].c];
    // maxUV[1] = S[1][S[1].c];
    // t=0,z=1;
    uint32_t wsSize;


    //find pivot in C of t
    // for(uint32_t i = S[t].c; i < S[t].r; i++) {//scan all vertices in C
    // // for(uint32_t i = S[t].c; i <= S[t].c; i++) {
    //     uint32_t u = S[t][i];
    //     uint32_t e = 0;

    //     if(g->deg(u, t) > S[z].r - S[z].c) {
    //         for(uint32_t j = S[z].c; j < S[z].r; j++) {
    //             uint32_t v = S[z][j];
    //             if(g->connect(u, v, t)) {
    //                 e++;
    //                 deg[v]++;
    //             }
    //         }
    //     }
    //     else {
    //         for(uint32_t j = g->pos[t][u]; j < g->pos[t][u + 1]; j++) {
    //             uint32_t v = g->e[t][j];
    //             uint32_t pv = S[z].pos(v);
    //             if(S[z].c <= pv && pv < S[z].r) {//in C
    //                 e++;
    //                 deg[v]++;
    //             }
    //         }
    //     }

    //     if(e > maxE[z]) {
    //         maxE[z] = e;
    //         maxI[t] = i;
    //         maxUV[t] = u;
    //     }
    // }
    // //find pivot in C of z
    // for(uint32_t i = S[1].c; i < S[1].r; i++) {//scan all vertices in C
    // // for(uint32_t i = S[t].c; i <= S[t].c; i++) {
    //     uint32_t u = S[1][i];
    //     uint32_t e = deg[u];
    //     deg[u] = 0;
    //     if(e > maxE[0]) {
    //         maxE[0] = e;
    //         maxI[1] = i;
    //         maxUV[1] = u;
    //     }
    // }

    // //check if pivot in X of both sides
    // bool isPivotInX[2] = {false, false};
    // for(t = 0; t < 2; t++) {
    //     z = t ^ 1;//t + z = 1

    //     for(uint32_t i = S[t].x; i < S[t].c; i++) {//scan all vertices in X
    //     // for(uint32_t i = S[t].c; i <= S[t].c; i++) {
    //         uint32_t u = S[t][i];
    //         uint32_t e = 0;

    //         if(g->deg(u, t) > S[z].r - S[z].c) {
    //             for(uint32_t j = S[z].c; j < S[z].r; j++) {
    //                 uint32_t v = S[z][j];
    //                 if(g->connect(u, v, t)) {
    //                     e++;
    //                 }
    //             }
    //         }
    //         else {
    //             for(uint32_t j = g->pos[t][u]; j < g->pos[t][u + 1]; j++) {
    //                 uint32_t v = g->e[t][j];
    //                 uint32_t pv = S[z].pos(v);
    //                 if(S[z].c <= pv && pv < S[z].r) {//in C
    //                     e++;
    //                 }
    //             }
    //         }

    //         if(e > maxE[z]) {
    //             isPivotInX[t] = true;
    //             maxE[z] = e;
    //             maxI[t] = i;
    //             maxUV[t] = u;
    //         }
    //     }
    // }

    // //pivot found. pivot must have at least one neighbor in this algorithm
    // t = 1;
    // if(S[0].cSize() - maxE[0] >= S[1].cSize() - maxE[1]) {
    //     t = 0;
    // }
    // z = t ^ 1;


    // //r is more like the right limit of C
    // //r[z] doesnt follow S[z].r
    // //pivot could be in X, and C of z is changed, but never mind, it will recover by r[z]
    // if(g->deg(maxUV[t], t) > S[z].r - S[z].c) {
    //     if(S[z].r>0){
    //         for(uint32_t i = S[z].r - 1; i >= S[z].c; i--) {
    //             if(!g->connect(maxUV[t], S[z][i], t)) {
    //                 S[z].swapByPos(--S[z].r, i);
    //             }
    
    //             if(i == 0) break;
    //         }
    //     }
    // }else{
    //     uint32_t u = maxUV[t];
    //     uint32_t tmpR = S[z].c;
    //     for(uint32_t j = g->pos[t][u]; j < g->pos[t][u + 1]; j++) {
    //         uint32_t v = g->e[t][j];
    //         uint32_t pv = S[z].pos(v);
    //         if(S[z].c <= pv && pv < S[z].r) {
    //             S[z].swapByPos(tmpR++, pv);
    //         }
    //     }
    //     S[z].r = tmpR;
    // }
    // wsSize = r[z] - S[z].r;
    // if(ws[deep].size() < wsSize) {
    //     ws[deep].resize(wsSize * 2);
    // }
    // memcpy(ws[deep].data(), S[z].begin() + S[z].r, sizeof(uint32_t) * wsSize);


    // //iterate on pivot
    // if(isPivotInX[t]==false){
        
    //     S[t].swapByPos(maxI[t], --r[t]);
    //     S[t].r = r[t];

    //     //x[z] doesnt follow S[z].x, used to recover S[z].x
    //     if(g->deg(maxUV[t], t) > S[z].c - S[z].x) {
    //         for(uint32_t i = S[z].x; i < S[z].c; i++) {
    //             if(!g->connect(maxUV[t], S[z][i], t)) {
    //                 S[z].swapByPos(S[z].x++, i);
    //             }
    //         }
    //     }else{
    //         uint32_t u = maxUV[t];
    //         uint32_t tmpX = S[z].c;
    //         for(uint32_t j = g->pos[t][u]; j < g->pos[t][u + 1]; j++) {
    //             uint32_t v = g->e[t][j];
    //             uint32_t pv = S[z].pos(v);
    //             if(S[z].x <= pv && pv < S[z].c) {
    //                 S[z].swapByPos(--tmpX, pv);
    //             }
    //         }
    //         S[z].x = tmpX;
    //     }



    //     RUV[t].push_back(maxUV[t]);
    //     RUVlook[t][maxUV[t]]=1;

    //     //debug
    //     // if(S[1].sz!=g->n[1]){
    //     //     std::cerr<<"error"<<std::endl;
    //     // }
    //     // if(thisNode==401812){
    //     //     std::cerr<<"here"<<std::endl;
    //     // }

    //     // #ifdef PARTIALPIVOT
    //     // int res;
    //     // if(S[z].cSize()+S[t].cSize()<=20){
    //     //     res=bbranch_simplep(deep+1,t,maxUV[t]);
    //     // }else{
    //     //     res=bbranch(deep+1,t,maxUV[t]);
    //     // }
    //     // #else
    //     //C of z must not be empty
    //     int res=bbranchLeafRes(deep+1,t,maxUV[t]);
    //     // #endif

    //     //debug
    //     // if(S[1].sz!=g->n[1]){
    //     //     std::cerr<<"error"<<std::endl;
    //     // }

    //     RUV[t].pop_back();
    //     RUVlook[t][maxUV[t]]=0;
    //     if(res!=-1){
    //         if(thisNode==-1){
    //             thisNode=runtree.getNewNode();
    //             runtree.treeNodePool[thisNode].inValue=inValue;
    //             runtree.treeNodePool[thisNode].uv=side;
    //             runtree.treeNodePool[thisNode].resMark=-1;
    //             runtree.treeNodePool[thisNode].qpMark=false;
    //         }
    //         runtree.treeNodePool[thisNode].sons.push_back(res);
    //         xStackNode.push_back(res);
    //     }

    //     //debug
    //     // if(S[1].sz!=g->n[1]){
    //     //     std::cerr<<"error"<<std::endl;
    //     // }

    //     //put pivot into X
    //     S[t].swapByPos(c[t]++, r[t]++);
    // }




    //iterate on both sides
    int tt=side;
    uint32_t wsize[2];
    for(int o=0;o<=1;++o){
        tt=tt^1;
        int zz=tt^1;
        wsSize=r[tt]-c[tt];
        wsize[tt]=wsSize;
        if(ws2[tt][deep].size() < wsSize) {
            ws2[tt][deep].resize(wsSize * 2);
        }
        memcpy(ws2[tt][deep].data(), S[tt].begin() + c[tt], sizeof(uint32_t) * wsSize);

        for(uint32_t j = 0; j < wsSize; j++){
            uint32_t w = ws2[tt][deep][j];
            S[tt].swapByPos(S[tt].pos(w), --r[tt]);

            //recover next level parameter
            S[0].x = x[0]; S[1].x = x[1];
            S[0].c = c[0]; S[1].c = c[1];
            S[0].r = r[0]; S[1].r = r[1];

            if(g->deg(w, tt) > S[zz].c - S[zz].x) {
                for(uint32_t i = S[zz].x; i < S[zz].c; i++) {
                    if(!g->connect(w, S[zz][i], tt)) {
                        S[zz].swapByPos(S[zz].x++, i);
                    }
                }
            }
            else {
                uint32_t tmpX = S[zz].c;
                for(uint32_t i = g->pos[tt][w]; i < g->pos[tt][w + 1]; i++) {
                    uint32_t u = g->e[tt][i];
                    uint32_t pu = S[zz].pos(u);
                    if(S[zz].x <= pu && pu < S[zz].c) {
                        S[zz].swapByPos(--tmpX, pu);
                    }
                }
                S[zz].x = tmpX;
            }

            if(g->deg(w, tt) > S[zz].r - S[zz].c) {
                if(S[zz].r > 0){
                    for(uint32_t i = S[zz].r - 1; i >= S[zz].c; i--) {
                        if(!g->connect(w, S[zz][i], tt)) {
                            S[zz].swapByPos(--S[zz].r, i);
                        }

                        if(i == 0) break;
                    }
                }
            }else{
                uint32_t tmpR = S[zz].c;
                for(uint32_t j = g->pos[tt][w]; j < g->pos[tt][w + 1]; j++) {
                    uint32_t u = g->e[tt][j];
                    uint32_t pu = S[zz].pos(u);
                    if(S[zz].c <= pu && pu < S[zz].r) {
                        S[zz].swapByPos(tmpR++, pu);
                    }
                }
                S[zz].r = tmpR;
            }

            // //debug
            // if(deep==0&&inValue==1859991&&w==1911392){
            //     // std::cerr<<<<std::endl;
            //     for(int ii=S[z].c;ii<S[z].r;++ii){
            //         std::cerr<<g->old_lables[z][S[z][ii]]<<" ";
            //     }
            //     std::cerr<<std::endl;
                
            // }

            // if(S[zz].cSize()>0){
            RUV[tt].push_back(w);
            RUVlook[tt][w]=1;
            //debug
            // if(S[1].sz!=g->n[1]){
            //     std::cerr<<"error"<<std::endl;
            // }

            int res=bbranchLeafRes_simplep(deep+1,tt,w);
            RUV[tt].pop_back();
            RUVlook[tt][w]=0;

            if(res!=-1){
                if(thisNode==-1){
                    thisNode=runtree.getNewNode();
                    runtree.treeNodePool[thisNode].inValue=inValue;
                    runtree.treeNodePool[thisNode].uv=side;
                    runtree.treeNodePool[thisNode].resMark=-1;
                    runtree.treeNodePool[thisNode].qpMark=false;
                }
                runtree.treeNodePool[thisNode].sons.push_back(res);
                xStackNode.push_back(res);
            }
            // }

            S[tt].swapByPos(c[tt]++,r[tt]++);//move w to X
        }
    }

    for(int o=0;o<=1;++o){
        for(uint32_t j = 0; j < wsize[o]; j++) {
            uint32_t w = ws2[o][deep][j];
            S[o].swapByPos(S[o].pos(w), --c[o]);
        }
    }

    //recover x,c,r
    S[t].x = x[t]; S[z].x = x[z];
    S[t].c = c[t]; S[z].c = c[z];
    S[t].r = r[t]; S[z].r = r[z];

    //connect that should be connected
    //almost the same as before
    if(RUV[t].size()>=threshold[t]&&RUV[z].size()>=threshold[z]){
        int res=-1;
        if(RUV[t].size()==threshold[t]&&t==oriSide||RUV[z].size()==threshold[z]&&z==oriSide){// only the first one
            res=connectAllBiCliqueLeafRes();
            if(thisNode!=-1){
                runtree.treeNodePool[thisNode].resMark=res;
            }
        }else{
            if(thisNode!=-1){
                fans.updateSize(cliqueuv[0].size());
                for(int i=0;i<runtree.treeNodePool[thisNode].sons.size();++i){
                    int sonNode=runtree.treeNodePool[thisNode].sons[i];
                    connectResults.push_back(runtree.treeNodePool[sonNode].resMark);
                }

                for(int i=1;i<connectResults.size();++i){
                    fans.merge(connectResults[i],connectResults[0]);
                }
                runtree.treeNodePool[thisNode].resMark=fans.find(connectResults[0]);
                connectResults.clear();
            }

        }


    }

    xStackNode.resize(oriXstackSize);
    return thisNode;
}

int BICPC::bbranch_simplep(uint32_t deep,uint32_t side,uint32_t inValue){
    // if(S[1].sz!=g->n[1]){
    //     std::cerr<<"error"<<std::endl;
    // }
    // if(S[z].CIsEmpty()) return;

    int thisNode=-1;//temporally no tree node is built here
    int oriSide=side;
    uint32_t t=side;
    uint32_t z=side ^ 1;

    if(S[z].cSize()==0){
        if(S[t].x==S[t].c&&S[t].c==S[t].r){
            if(S[z].x==S[z].c&&S[z].c==S[z].r){
                return thisNode;
            }
        }
        justConnect(oriSide);
        return thisNode;
    }

    //TASK1
    //Check if the current R forms a answer with the CAND on the opposite side of newly chosen vertex
    bool found=true;
    int curCliqueId=-1;
    if(S[z].XIsEmpty()&&RUV[t].size()>=threshold[t]&&(RUV[z].size()+S[z].r-S[z].c>=threshold[z])){
        if(RUV[t].size()+RUV[z].size()+S[z].r-S[z].c>threshold[t]+threshold[z]){
            for(uint32_t i = S[t].x; i < S[t].r; i++) {
                //zy
                // printf("here\n");

                bool connectAll = true;
                for(uint32_t j = S[z].c; j < S[z].r; j++) {
                    if(!g->connect(S[t][i], S[z][j], t)) {
                        connectAll = false;
                        break;
                    }
                }
                if(connectAll) {
                    found = false;
                    break;
                }
            }

            if(found){//found a maximal biclique: RUV and cand of z
                thisNode=runtree.getNewNode();
                runtree.treeNodePool[thisNode].inValue=inValue;
                runtree.treeNodePool[thisNode].uv=side;
                runtree.treeNodePool[thisNode].resMark=-1;
                runtree.treeNodePool[thisNode].qpMark=false;

                cliqueuv[t].emplace_back(RUV[t]);
                int oldsize=RUV[z].size();
                RUV[z].reserve(RUV[z].size()+S[z].r-S[z].c);
                for(int i=S[z].c;i<S[z].r;++i){
                    RUV[z].push_back(S[z][i]);
                }
                cliqueuv[z].emplace_back(RUV[z]);
                // runtree.treeNodePool[thisNode].resMark=cliqueuv[t].size()-1;

                //reset RUV
                RUV[z].resize(oldsize);

                runtree.treeNodePool[thisNode].resMark=cliqueuv[t].size()-1;
                curCliqueId=cliqueuv[t].size()-1;
            }
        }
    }


    //TASK2, make sure C of t,z have edge(s), by shrink C of z
    // uint32_t newR = S[z].r;
    // uint32_t oldR = S[z].r;
    // for(uint32_t i = S[z].c; i < newR; ) {
    //     bool degZero = true;
    //     for(uint32_t j = S[t].c; j < S[t].r; j++) {
    //         if(g->connect(S[z][i], S[t][j], z)) {
    //             degZero = false;
    //             break;
    //         }
    //     }
    //     if(degZero) {
    //         S[z].swapByPos(i, --newR);
    //     }
    //     else i++;
    // }
    // S[z].r = newR;
    bool emptyc=(S[t].cSize()==0||S[z].cSize()==0);

    if(emptyc){
        if(RUV[0].size()>=threshold[0]&&RUV[1].size()>=threshold[1]){
            // S[z].r=oldR;
            // if(thisNode!=-1){
            //     connectResults.push_back(runtree.treeNodePool[thisNode].resMark);
            // }
            // int res=connectBiCliqueInCurZone(oriSide);
            // if(res!=-1){
            //     connectResults.push_back(res);
            // }
            // if((S[t].x<S[t].c)||(S[z].x<S[z].c)){
            //     // fans.updateSize(cliqueuv[0].size());

            //     // if(curCliqueId!=-1) curStackNode.push_back(curCliqueId);
            //     connectBiCliqueInXzone(oriSide);

            //     // if(curCliqueId!=-1) curStackNode.pop_back();

            // }

            // // if(connectResults.size()==0){
            // //     std::cerr<<"final connect id error"<<std::endl;
            // //     exit(0);
            // // }

            // fans.updateSize(cliqueuv[0].size());
            // if(connectResults.size()>1){
            //     for(int i=1;i<connectResults.size();++i){
            //         fans.merge(connectResults[0],connectResults[i]);
            //     }
            // }

            int res=connectAllBiCliqueNew(oriSide,thisNode);

            if(thisNode!=-1){
                // if(res==-1){
                //     std::cerr<<"error definitely"<<std::endl;
                //     exit(0);
                // }
                // if(runtree.treeNodePool[thisNode].resMark==-1){
                //     std::cerr<<"resmark error"<<std::endl;
                //     exit(0);
                // }

                // if(res!=-1){
                //     fans.merge(runtree.treeNodePool[thisNode].resMark,res);
                //     runtree.treeNodePool[thisNode].resMark=fans.find(res);
                // }

                runtree.treeNodePool[thisNode].qpMark=true;
                runtree.treeNodePool[thisNode].resMark=res;
            }

            connectResults.clear();
            // S[z].r=newR;
        }
        return thisNode;
    }

    //TASK3, find pivot and start new recursions
    // if(S[t].cSize()==0||S[z].cSize()==0){
    //     //to connect
    //     if(RUV[0].size()>=p&&RUV[1].size()>=q){
    //         if(thisNode!=-1){
    //             connectResults.push_back(runtree.treeNodePool[thisNode].resMark);
    //         }
    //         int res=connectBiCliqueInCurZone(oriSide);
    //         if(res!=-1){
    //             connectResults.push_back(res);
    //         }
    //         if((S[0].x<S[0].c||S[1].x<S[1].c)){
    //             // fans.updateSize(cliqueuv[0].size());

    //             // if(curCliqueId!=-1) curStackNode.push_back(curCliqueId);
    //             connectBiCliqueInXzone(oriSide);

    //             // if(curCliqueId!=-1) curStackNode.pop_back();

    //         }

    //         if(connectResults.size()==0){
    //             std::cerr<<"final connect id error"<<std::endl;
    //             exit(0);
    //         }

    //         fans.updateSize(cliqueuv[0].size());
    //         if(connectResults.size()>1){
    //             for(int i=1;i<connectResults.size();++i){
    //                 fans.merge(connectResults[0],connectResults[i]);
    //             }
    //         }

    //         if(thisNode!=-1){
    //             // if(runtree.treeNodePool[thisNode].resMark==-1){
    //             //     std::cerr<<"resmark error"<<std::endl;
    //             //     exit(0);
    //             // }

    //             // if(res!=-1){
    //             //     fans.merge(runtree.treeNodePool[thisNode].resMark,res);
    //             //     runtree.treeNodePool[thisNode].resMark=fans.find(res);
    //             // }

    //             runtree.treeNodePool[thisNode].qpMark=true;
    //             runtree.treeNodePool[thisNode].resMark=fans.find(connectResults[0]);
    //         }

    //         connectResults.clear();

    //     }
    //     return thisNode;
    // }

    if(RUV[t].size()+S[t].cSize()<threshold[t]||RUV[z].size()+S[z].cSize()<threshold[z]){
        return thisNode;
    }

    //debug
    // if(S[1].sz!=g->n[1]){
    //     std::cerr<<"error"<<std::endl;
    // }

    // #ifdef NOSHRINK
    // S[z].r=oldR;
    // #endif


    //initialize copies of xcr
    uint32_t x[2];
    uint32_t c[2];
    uint32_t r[2];
    x[t] = S[t].x; x[z] = S[z].x;
    c[t] = S[t].c; c[z] = S[z].c;
    r[t] = S[t].r; r[z] = S[z].r;

    int oriXstackSize=xStackNode.size();

    //find pivot
    // uint32_t maxI[2] = {0}, maxE[2] = {0}, maxUV[2] = {0};
    // maxI[0] = S[0].c;
    // maxI[1] = S[1].c;
    // maxUV[0] = S[0][S[0].c];
    // maxUV[1] = S[1][S[1].c];
    // t=0,z=1;
    uint32_t wsSize;


    //find pivot in C of t
    // for(uint32_t i = S[t].c; i < S[t].r; i++) {//scan all vertices in C
    // // for(uint32_t i = S[t].c; i <= S[t].c; i++) {
    //     uint32_t u = S[t][i];
    //     uint32_t e = 0;

    //     if(g->deg(u, t) > S[z].r - S[z].c) {
    //         for(uint32_t j = S[z].c; j < S[z].r; j++) {
    //             uint32_t v = S[z][j];
    //             if(g->connect(u, v, t)) {
    //                 e++;
    //                 deg[v]++;
    //             }
    //         }
    //     }
    //     else {
    //         for(uint32_t j = g->pos[t][u]; j < g->pos[t][u + 1]; j++) {
    //             uint32_t v = g->e[t][j];
    //             uint32_t pv = S[z].pos(v);
    //             if(S[z].c <= pv && pv < S[z].r) {//in C
    //                 e++;
    //                 deg[v]++;
    //             }
    //         }
    //     }

    //     if(e > maxE[z]) {
    //         maxE[z] = e;
    //         maxI[t] = i;
    //         maxUV[t] = u;
    //     }
    // }
    // //find pivot in C of z
    // for(uint32_t i = S[1].c; i < S[1].r; i++) {//scan all vertices in C
    // // for(uint32_t i = S[t].c; i <= S[t].c; i++) {
    //     uint32_t u = S[1][i];
    //     uint32_t e = deg[u];
    //     deg[u] = 0;
    //     if(e > maxE[0]) {
    //         maxE[0] = e;
    //         maxI[1] = i;
    //         maxUV[1] = u;
    //     }
    // }

    // //check if pivot in X of both sides
    // bool isPivotInX[2] = {false, false};
    // for(t = 0; t < 2; t++) {
    //     z = t ^ 1;//t + z = 1

    //     for(uint32_t i = S[t].x; i < S[t].c; i++) {//scan all vertices in X
    //     // for(uint32_t i = S[t].c; i <= S[t].c; i++) {
    //         uint32_t u = S[t][i];
    //         uint32_t e = 0;

    //         if(g->deg(u, t) > S[z].r - S[z].c) {
    //             for(uint32_t j = S[z].c; j < S[z].r; j++) {
    //                 uint32_t v = S[z][j];
    //                 if(g->connect(u, v, t)) {
    //                     e++;
    //                 }
    //             }
    //         }
    //         else {
    //             for(uint32_t j = g->pos[t][u]; j < g->pos[t][u + 1]; j++) {
    //                 uint32_t v = g->e[t][j];
    //                 uint32_t pv = S[z].pos(v);
    //                 if(S[z].c <= pv && pv < S[z].r) {//in C
    //                     e++;
    //                 }
    //             }
    //         }

    //         if(e > maxE[z]) {
    //             isPivotInX[t] = true;
    //             maxE[z] = e;
    //             maxI[t] = i;
    //             maxUV[t] = u;
    //         }
    //     }
    // }

    // //pivot found. pivot must have at least one neighbor in this algorithm
    // t = 1;
    // if(S[0].cSize() - maxE[0] >= S[1].cSize() - maxE[1]) {
    //     t = 0;
    // }
    // z = t ^ 1;


    // //r is more like the right limit of C
    // //r[z] doesnt follow S[z].r
    // //pivot could be in X, and C of z is changed, but never mind, it will recover by r[z]
    // if(g->deg(maxUV[t], t) > S[z].r - S[z].c) {
    //     for(uint32_t i = S[z].r - 1; i >= S[z].c; i--) {
    //         if(!g->connect(maxUV[t], S[z][i], t)) {
    //             S[z].swapByPos(--S[z].r, i);
    //         }

    //         if(i == 0) break;
    //     }
    // }else{
    //     uint32_t u = maxUV[t];
    //     uint32_t tmpR = S[z].c;
    //     for(uint32_t j = g->pos[t][u]; j < g->pos[t][u + 1]; j++) {
    //         uint32_t v = g->e[t][j];
    //         uint32_t pv = S[z].pos(v);
    //         if(S[z].c <= pv && pv < S[z].r) {
    //             S[z].swapByPos(tmpR++, pv);
    //         }
    //     }
    //     S[z].r = tmpR;
    // }
    // wsSize = r[z] - S[z].r;
    // if(ws[deep].size() < wsSize) {
    //     ws[deep].resize(wsSize * 2);
    // }
    // memcpy(ws[deep].data(), S[z].begin() + S[z].r, sizeof(uint32_t) * wsSize);


    //iterate on pivot
    // if(isPivotInX[t]==false){
        
    //     S[t].swapByPos(maxI[t], --r[t]);
    //     S[t].r = r[t];

    //     //x[z] doesnt follow S[z].x, used to recover S[z].x
    //     if(g->deg(maxUV[t], t) > S[z].c - S[z].x) {
    //         for(uint32_t i = S[z].x; i < S[z].c; i++) {
    //             if(!g->connect(maxUV[t], S[z][i], t)) {
    //                 S[z].swapByPos(S[z].x++, i);
    //             }
    //         }
    //     }else{
    //         uint32_t u = maxUV[t];
    //         uint32_t tmpX = S[z].c;
    //         for(uint32_t j = g->pos[t][u]; j < g->pos[t][u + 1]; j++) {
    //             uint32_t v = g->e[t][j];
    //             uint32_t pv = S[z].pos(v);
    //             if(S[z].x <= pv && pv < S[z].c) {
    //                 S[z].swapByPos(--tmpX, pv);
    //             }
    //         }
    //         S[z].x = tmpX;
    //     }



    //     RUV[t].push_back(maxUV[t]);
    //     RUVlook[t][maxUV[t]]=1;
    //     if(curCliqueId!=-1) curStackNode.push_back(curCliqueId);

    //     //debug
    //     // if(S[1].sz!=g->n[1]){
    //     //     std::cerr<<"error"<<std::endl;
    //     // }
    //     // if(thisNode==401812){
    //     //     std::cerr<<"here"<<std::endl;
    //     // }

    //     //C of z must not be empty
    //     int res=bbranch(deep+1,t,maxUV[t]);

    //     //debug
    //     // if(S[1].sz!=g->n[1]){
    //     //     std::cerr<<"error"<<std::endl;
    //     // }

    //     RUV[t].pop_back();
    //     RUVlook[t][maxUV[t]]=0;
    //     if(curCliqueId!=-1) curStackNode.pop_back();
    //     if(res!=-1){
    //         if(thisNode==-1){
    //             thisNode=runtree.getNewNode();
    //             runtree.treeNodePool[thisNode].inValue=inValue;
    //             runtree.treeNodePool[thisNode].uv=side;
    //             runtree.treeNodePool[thisNode].resMark=-1;
    //             runtree.treeNodePool[thisNode].qpMark=false;
    //         }
    //         runtree.treeNodePool[thisNode].sons.push_back(res);
    //         xStackNode.push_back(res);
    //     }

    //     //debug
    //     // if(S[1].sz!=g->n[1]){
    //     //     std::cerr<<"error"<<std::endl;
    //     // }

    //     //put pivot into X
    //     S[t].swapByPos(c[t]++, r[t]++);
    // }


    //iterate on both sides
    int tt=side;
    uint32_t wsize[2];
    for(int o=0;o<=1;++o){
        tt=tt^1;
        int zz=tt^1;
        wsSize=r[tt]-c[tt];
        wsize[tt]=wsSize;
        if(ws2[tt][deep].size() < wsSize) {
            ws2[tt][deep].resize(wsSize * 2);
        }
        memcpy(ws2[tt][deep].data(), S[tt].begin() + c[tt], sizeof(uint32_t) * wsSize);

        for(uint32_t j = 0; j < wsSize; j++){
            uint32_t w = ws2[tt][deep][j];
            S[tt].swapByPos(S[tt].pos(w), --r[tt]);

            //recover next level parameter
            S[0].x = x[0]; S[1].x = x[1];
            S[0].c = c[0]; S[1].c = c[1];
            S[0].r = r[0]; S[1].r = r[1];

            if(g->deg(w, tt) > S[zz].c - S[zz].x) {
                for(uint32_t i = S[zz].x; i < S[zz].c; i++) {
                    if(!g->connect(w, S[zz][i], tt)) {
                        S[zz].swapByPos(S[zz].x++, i);
                    }
                }
            }
            else {
                uint32_t tmpX = S[zz].c;
                for(uint32_t i = g->pos[tt][w]; i < g->pos[tt][w + 1]; i++) {
                    uint32_t u = g->e[tt][i];
                    uint32_t pu = S[zz].pos(u);
                    if(S[zz].x <= pu && pu < S[zz].c) {
                        S[zz].swapByPos(--tmpX, pu);
                    }
                }
                S[zz].x = tmpX;
            }

            if(g->deg(w, tt) > S[zz].r - S[zz].c) {
                if(S[zz].r > 0){
                    for(uint32_t i = S[zz].r - 1; i >= S[zz].c; i--) {
                        if(!g->connect(w, S[zz][i], tt)) {
                            S[zz].swapByPos(--S[zz].r, i);
                        }

                        if(i == 0) break;
                    }
                }
            }else{
                uint32_t tmpR = S[zz].c;
                for(uint32_t j = g->pos[tt][w]; j < g->pos[tt][w + 1]; j++) {
                    uint32_t u = g->e[tt][j];
                    uint32_t pu = S[zz].pos(u);
                    if(S[zz].c <= pu && pu < S[zz].r) {
                        S[zz].swapByPos(tmpR++, pu);
                    }
                }
                S[zz].r = tmpR;
            }

            // //debug
            // if(deep==0&&inValue==1859991&&w==1911392){
            //     // std::cerr<<<<std::endl;
            //     for(int ii=S[z].c;ii<S[z].r;++ii){
            //         std::cerr<<g->old_lables[z][S[z][ii]]<<" ";
            //     }
            //     std::cerr<<std::endl;
                
            // }

            // if(S[zz].cSize()>0){
            RUV[tt].push_back(w);
            RUVlook[tt][w]=1;
            if(curCliqueId!=-1) curStackNode.push_back(curCliqueId);
            //debug
            // if(S[1].sz!=g->n[1]){
            //     std::cerr<<"error"<<std::endl;
            // }

            int res=bbranch_simplep(deep+1,tt,w);
            if(curCliqueId!=-1) curStackNode.pop_back();
            RUV[tt].pop_back();
            RUVlook[tt][w]=0;

            if(res!=-1){
                if(thisNode==-1){
                    thisNode=runtree.getNewNode();
                    runtree.treeNodePool[thisNode].inValue=inValue;
                    runtree.treeNodePool[thisNode].uv=side;
                    runtree.treeNodePool[thisNode].resMark=-1;
                    runtree.treeNodePool[thisNode].qpMark=false;
                }
                runtree.treeNodePool[thisNode].sons.push_back(res);
                xStackNode.push_back(res);
            }
            // }

            S[tt].swapByPos(c[tt]++,r[tt]++);//move w to X
        }
    }

    for(int o=0;o<=1;++o){
        for(uint32_t j = 0; j < wsize[o]; j++) {
            uint32_t w = ws2[o][deep][j];
            S[o].swapByPos(S[o].pos(w), --c[o]);
        }
    }

    // if(c[t]<r[t]){//S[t].c>0
    //     //iterate on vertices beyound pivot
    //     for(uint32_t j = 0; j < wsSize; j++){
    //         uint32_t w = ws[deep][j];
    //         S[z].swapByPos(S[z].pos(w), --r[z]);

    //         //recover next level parameter
    //         S[t].x = x[t]; S[z].x = x[z];
    //         S[t].c = c[t]; S[z].c = c[z];
    //         S[t].r = r[t]; S[z].r = r[z];

    //         if(g->deg(w, z) > S[t].c - S[t].x) {
    //             for(uint32_t i = S[t].x; i < S[t].c; i++) {
    //                 if(!g->connect(w, S[t][i], z)) {
    //                     S[t].swapByPos(S[t].x++, i);
    //                 }
    //             }
    //         }
    //         else {
    //             uint32_t tmpX = S[t].c;
    //             for(uint32_t j = g->pos[z][w]; j < g->pos[z][w + 1]; j++) {
    //                 uint32_t u = g->e[z][j];
    //                 uint32_t pu = S[t].pos(u);
    //                 if(S[t].x <= pu && pu < S[t].c) {
    //                     S[t].swapByPos(--tmpX, pu);
    //                 }
    //             }
    //             S[t].x = tmpX;
    //         }

    //         if(g->deg(w, z) > S[t].c - S[t].r) {
    //             if(S[t].r > 0){
    //                 for(uint32_t i = S[t].r - 1; i >= S[t].c; i--) {
    //                     if(!g->connect(w, S[t][i], z)) {
    //                         S[t].swapByPos(--S[t].r, i);
    //                     }

    //                     if(i == 0) break;
    //                 }
    //             }
    //         }else{
    //             uint32_t tmpR = S[t].c;
    //             for(uint32_t j = g->pos[z][w]; j < g->pos[z][w + 1]; j++) {
    //                 uint32_t u = g->e[z][j];
    //                 uint32_t pu = S[t].pos(u);
    //                 if(S[t].c <= pu && pu < S[t].r) {
    //                     S[t].swapByPos(tmpR++, pu);
    //                 }
    //             }
    //             S[t].r = tmpR;
    //         }

    //         // //debug
    //         // if(deep==0&&inValue==1859991&&w==1911392){
    //         //     // std::cerr<<<<std::endl;
    //         //     for(int ii=S[z].c;ii<S[z].r;++ii){
    //         //         std::cerr<<g->old_lables[z][S[z][ii]]<<" ";
    //         //     }
    //         //     std::cerr<<std::endl;
                
    //         // }

    //         if(S[t].cSize()>0){
    //             RUV[z].push_back(w);
    //             RUVlook[z][w]=1;
    //             if(curCliqueId!=-1) curStackNode.push_back(curCliqueId);
    //             //debug
    //             // if(S[1].sz!=g->n[1]){
    //             //     std::cerr<<"error"<<std::endl;
    //             // }

    //             int res=bbranch(deep+1,z,w);
    //             if(curCliqueId!=-1) curStackNode.pop_back();
    //             RUV[z].pop_back();
    //             RUVlook[z][w]=0;

    //             if(res!=-1){
    //                 if(thisNode==-1){
    //                     thisNode=runtree.getNewNode();
    //                     runtree.treeNodePool[thisNode].inValue=inValue;
    //                     runtree.treeNodePool[thisNode].uv=side;
    //                     runtree.treeNodePool[thisNode].resMark=-1;
    //                     runtree.treeNodePool[thisNode].qpMark=false;
    //                 }
    //                 runtree.treeNodePool[thisNode].sons.push_back(res);
    //                 xStackNode.push_back(res);
    //             }
    //         }

    //         S[z].swapByPos(c[z]++,r[z]++);//move w to X
    //     }

    //     for(uint32_t j = 0; j < wsSize; j++) {
    //         uint32_t w = ws[deep][j];
    //         S[z].swapByPos(S[z].pos(w), --c[z]);
    //     }
    // }

    //recover pivot
    // if(isPivotInX[t] == false) {
    //     S[t].swapByPos(S[t].pos(maxUV[t]), --c[t]);
    // }

    //recover x,c,r
    S[t].x = x[t]; S[z].x = x[z];
    S[t].c = c[t]; S[z].c = c[z];
    S[t].r = r[t]; S[z].r = r[z];

    //connect that should be connected
    //almost the same as before
    if(RUV[0].size()>=p&&RUV[1].size()>=q){
        // fans.updateSize(cliqueuv[0].size());
        // S[oriSide^1].r=oldR;
        // if(thisNode!=-1&&runtree.treeNodePool[thisNode].resMark!=-1){
        //     connectResults.push_back(runtree.treeNodePool[thisNode].resMark);
        // }
        // int res=connectBiCliqueInCurZone(oriSide);
        // if(res!=-1){
        //     connectResults.push_back(res);
        // }

        // if(S[oriSide^1].x<S[oriSide^1].c||)
        // connectBiCliqueInXzone(oriSide);

        // if(connectResults.size()==0){
        //     std::cerr<<"final connect id error"<<std::endl;
        //     exit(0);
        // }

        // fans.updateSize(cliqueuv[0].size());
        // if(connectResults.size()>1){
        //     for(int i=1;i<connectResults.size();++i){
        //         fans.merge(connectResults[0],connectResults[i]);
        //     }
        // }
        int res=connectAllBiCliqueNew(oriSide,thisNode);

        if(thisNode!=-1){
            // if(res==-1){
            //     std::cerr<<"error definitely"<<std::endl;
            //     exit(0);
            // }
            runtree.treeNodePool[thisNode].qpMark=true;
            runtree.treeNodePool[thisNode].resMark=res;
        }

        connectResults.clear();
        // S[oriSide^1].r=newR;
    }

    xStackNode.resize(oriXstackSize);
    return thisNode;
}

void BICPC::prune_find_combinations_base(uint64_t &comb_count, uint32_t beg_offset, uint32_t curr_depth, uint32_t total_depth,int validPosR){
    #ifdef COUNTBICLIQUETREENODE
    bicliquelistNodecnt++;
    #endif

    uint32_t N = seed_vertices.size();
    if (curr_depth == 0) {
        // cliqueuv[opAnchor].emplace_back(temp_comb);
        // cliqueuv[anchor].emplace_back(RP);
        // combs.emplace_back(comb);
        comb_count++;

        for(int i=1;i<=validPosR;++i){
            fans.merge(commonMbclique[i],commonMbclique[0]);
        }

        return;
    }

    // if(N-beg_offset<curr_depth) return;

    for (uint32_t i = beg_offset; i < N; i++){

        if(N-i<curr_depth) break;

        //check common mbc
        int v=seed_vertices[i];
        int nvalidPosR=0;
        for(int j=0;j<v2mbc[opAnchor][v].size();++j){
            int bc=v2mbc[opAnchor][v][j];
            if(commonMbcliquePos[bc]<=validPosR){
                exchange(commonMbclique,commonMbcliquePos[bc],nvalidPosR++,commonMbcliquePos);
            }else break;
        }
        nvalidPosR--;
        // if(checkCommonGroup(nvalidPosR)==false) continue;

        temp_comb[total_depth-curr_depth] = v;

        //shrink v2mbc
        if(curr_depth-1>0){
            for(int j=i+1;j<N;++j){
                int w=seed_vertices[j];
                int tempnr=0;
                for(int k=0;k<v2mbc[opAnchor][w].size();++k){
                    int bc=v2mbc[opAnchor][w][k];
                    if(commonMbcliquePos[bc]>validPosR) break;
                    if(commonMbcliquePos[bc]<=nvalidPosR){
                        exchange(v2mbc[opAnchor][w],k,tempnr++);
                    }
                }
            }
        }

        prune_find_combinations_base(comb_count, i + 1, curr_depth - 1, total_depth,nvalidPosR);
    }
}

//shrink v2mbc before use
void BICPC::prune_find_combinations(uint64_t &comb_count, uint32_t beg_offset, uint32_t curr_depth, uint32_t total_depth,int validPosR){

    #ifdef COUNTBICLIQUETREENODE
    bicliquelistNodecnt++;
    #endif

    uint32_t N = seed_vertices.size();
    if (curr_depth == 0) {
        // cliqueuv[opAnchor].emplace_back(temp_comb);
        // cliqueuv[anchor].emplace_back(RP);
        // combs.emplace_back(comb);
        comb_count++;

        for(int i=1;i<=validPosR;++i){
            fans.merge(commonMbclique[i],commonMbclique[0]);
        }

        return;
    }

    // if(N-beg_offset<curr_depth) return;

    for (uint32_t i = beg_offset; i < N; i++){

        if(N-i<curr_depth) break;

        //check common mbc
        int v=seed_vertices[i];
        int nvalidPosR=0;
        for(int j=0;j<v2mbc[opAnchor][v].size();++j){
            int bc=v2mbc[opAnchor][v][j];
            if(commonMbcliquePos[bc]<=validPosR){
                exchange(commonMbclique,commonMbcliquePos[bc],nvalidPosR++,commonMbcliquePos);
            }else break;
        }
        nvalidPosR--;
        if(checkCommonGroup(nvalidPosR)==false) continue;

        temp_comb[total_depth-curr_depth] = v;

        //shrink v2mbc
        if(curr_depth-1>0){
            for(int j=i+1;j<N;++j){
                int w=seed_vertices[j];
                int tempnr=0;
                for(int k=0;k<v2mbc[opAnchor][w].size();++k){
                    int bc=v2mbc[opAnchor][w][k];
                    if(commonMbcliquePos[bc]>validPosR) break;
                    if(commonMbcliquePos[bc]<=nvalidPosR){
                        exchange(v2mbc[opAnchor][w],k,tempnr++);
                    }
                }
            }
        }

        prune_find_combinations(comb_count, i + 1, curr_depth - 1, total_depth,nvalidPosR);
    }
}

bool BICPC::checkCommonGroup(int validPosR){
    for(int i=0;i<=validPosR;++i){
        int bc=commonMbclique[i];
        int father=fans.find(bc);
        if(visGroup[father]==false){
            visGroup[father]=true;
            visGroups.push_back(father);
        }
    }

    for(int i=0;i<visGroups.size();++i){
        visGroup[visGroups[i]]=false;
    }
    int num=visGroups.size();
    visGroups.clear();

    return num>1;
}

void BICPC::conductListAndConnectMBCbase(int left,int validPosR){

    #ifdef COUNTBICLIQUETREENODE
    bicliquelistNodecnt++;
    #endif

    if(left==2){
        for(int i=0;i<anchorCandLayer[left].size();++i){
            int u=anchorCandLayer[left][i];

            // //debug
            // if(u==1826028||u==1934355||u==1952117){
            //     std::cerr<<"here"<<std::endl;
            // }

            //check common mbc
            int nvalidPosR=0;
            for(int j=0;j<v2mbc[anchor][u].size();++j){
                int bc=v2mbc[anchor][u][j];
                if(commonMbcliquePos[bc]<=validPosR){
                    exchange(commonMbclique,commonMbcliquePos[bc],nvalidPosR++,commonMbcliquePos);
                }else break;
            }
            nvalidPosR--;
            // if(checkCommonGroup(nvalidPosR)==false) continue;


            RP.push_back(u);
            op_vertices[2].clear();
            if(anchorThreshold<=2){
                for(int j=g->pos[anchor][u];j<g->pos[anchor][u+1];++j){
                    int v=g->e[anchor][j];
                    op_vertices[2].push_back(v);
                }
            }else{
                int j=0,k=g->pos[anchor][u];
                while(j<op_vertices[3].size()&&k<g->pos[anchor][u+1]){
                    if(op_vertices[3][j]==g->e[anchor][k]){
                        op_vertices[2].push_back(op_vertices[3][j]);
                        j++;
                        k++;
                    }else if(op_vertices[3][j]<g->e[anchor][k]){
                        j++;
                    }else{
                        k++;
                    }
                }
            }

            if(op_vertices[2].size()<opThreshold){
                RP.pop_back();
                continue;
            }

            for(int j=g->two_hop_pos[u];j<two_hop_pos_end_layer[2][u];++j){
                int v=g->two_hop_e[j];

                // //debug
                // if(v==1826028||v==1934355||v==1952117){
                //     std::cerr<<"here"<<std::endl;
                // }

                //check common mbc
                int nr=0;
                for(int k=0;k<v2mbc[anchor][v].size();++k){
                    int bc=v2mbc[anchor][v][k];
                    if(commonMbcliquePos[bc]<=nvalidPosR){
                        exchange(commonMbclique,commonMbcliquePos[bc],nr++,commonMbcliquePos);
                    }else if(commonMbcliquePos[bc]>validPosR){// because no shrink
                        break;
                    }
                }
                nr--;
                // if(checkCommonGroup(nr)==false) continue;


                RP.push_back(v);
                op_vertices[1].clear();
                int k=0,o=g->pos[anchor][v];
                while(k<op_vertices[2].size()&&o<g->pos[anchor][v+1]){
                    if(op_vertices[2][k]==g->e[anchor][o]){
                        op_vertices[1].push_back(op_vertices[2][k]);
                        k++;
                        o++;
                    }else if(op_vertices[2][k]<g->e[anchor][o]){
                        k++;
                    }else{
                        o++;
                    }
                }
                if(op_vertices[1].size()<opThreshold){
                    RP.pop_back();
                    continue;
                }

                #ifdef COUNTQPONLY
                qpcliqueNum+=count_combinations(op_vertices[1].size(),opThreshold);
                #else
                seed_vertices.clear();
                for(int k=0;k<op_vertices[1].size();++k){
                    seed_vertices.push_back(op_vertices[1][k]);
                }
                uint64_t comb_count = 0;
                // std::vector< std::vector<int> > combinations;
                uint32_t comb_size = opThreshold;
                // vector<int> temp_comb(comb_size);
                // find_combinations(comb_count, 0, comb_size, comb_size);

                //shrink v2mbc
                for(int k=0;k<seed_vertices.size();++k){
                    int w=seed_vertices[k];
                    int tempnr=0;
                    for(int o=0;o<v2mbc[opAnchor][w].size();++o){
                        int bc=v2mbc[opAnchor][w][o];
                        if(commonMbcliquePos[bc]<=nr){
                            exchange(v2mbc[opAnchor][w],o,tempnr++);
                        }
                    }
                }

                //debug
                // if(RP[0]==1826028&&RP[1]==1934355&&RP[2]==1952117){
                //     std::cerr<<"here"<<std::endl;
                // }

                prune_find_combinations_base(comb_count, 0, comb_size, comb_size,nr);
                #endif

                RP.pop_back();
            }

            RP.pop_back();
        }
        return;
    }

    for(int i=0;i<anchorCandLayer[left].size();++i){
        int u=anchorCandLayer[left][i];

        // //debug
        // if(u==1826028||u==1934355||u==1952117){
        //     std::cerr<<"here"<<std::endl;
        // }

        //check common mbc
        int nvalidPosR=0;
        for(int j=0;j<v2mbc[anchor][u].size();++j){
            int bc=v2mbc[anchor][u][j];
            if(commonMbcliquePos[bc]<=validPosR){
                exchange(commonMbclique,commonMbcliquePos[bc],nvalidPosR++,commonMbcliquePos);
            }else break;
        }
        nvalidPosR--;
        // if(checkCommonGroup(nvalidPosR)==false) continue;



        RP.push_back(u);
        op_vertices[left].clear();
        if(left==anchorThreshold){
            for(int j=g->pos[anchor][u];j<g->pos[anchor][u+1];++j){
                int v=g->e[anchor][j];
                op_vertices[left].push_back(v);
            }
        }else{
            int j=0,k=g->pos[anchor][u];
            while(j<op_vertices[left+1].size()&&k<g->pos[anchor][u+1]){
                if(op_vertices[left+1][j]==g->e[anchor][k]){
                    op_vertices[left].push_back(op_vertices[left+1][j]);
                    j++;
                    k++;
                }else if(op_vertices[left+1][j]<g->e[anchor][k]){
                    j++;
                }else{
                    k++;
                }
            }
        }

        if(op_vertices[left].size()<opThreshold){
            RP.pop_back();
            continue;
        }

        anchorCandLayer[left-1].clear();
        for(int j=g->two_hop_pos[u];j<two_hop_pos_end_layer[left][u];++j){
            int v=g->two_hop_e[j];
            if(anchorlab[v]==left){
                anchorlab[v]=left-1;
                anchorCandLayer[left-1].push_back(v);
                two_hop_pos_end_layer[left-1][v]=g->two_hop_pos[v];
            }
        }
        for(int j=0;j<anchorCandLayer[left-1].size();++j){
            int v=anchorCandLayer[left-1][j];
            int end=two_hop_pos_end_layer[left][v];
            for(int k=g->two_hop_pos[v];k<end;++k){
                int w=g->two_hop_e[k];
                if(anchorlab[w]==left-1){
                    two_hop_pos_end_layer[left-1][v]++;
                }else{
                    g->two_hop_e[k--]=g->two_hop_e[--end];
                    g->two_hop_e[end]=w;
                }
            }


            //shrink v2mbc
            int tempnr=0;
            for(int k=0;k<v2mbc[anchor][v].size();++k){
                int bc=v2mbc[anchor][v][k];
                if(commonMbcliquePos[bc]>validPosR) break;
                if(commonMbcliquePos[bc]<=nvalidPosR){
                    exchange(v2mbc[anchor][v],k,tempnr++);
                }
            }
        }

        conductListAndConnectMBCbase(left-1,nvalidPosR);
        for(int j=0;j<anchorCandLayer[left-1].size();++j){
            int v=anchorCandLayer[left-1][j];
            anchorlab[v]=left;
        }

        RP.pop_back();
    }
}

void BICPC::conductListAndConnectMBC(int left,int validPosR){

    #ifdef COUNTBICLIQUETREENODE
    bicliquelistNodecnt++;
    #endif

    if(left==2){
        for(int i=0;i<anchorCandLayer[left].size();++i){
            int u=anchorCandLayer[left][i];

            // //debug
            // if(u==1826028||u==1934355||u==1952117){
            //     std::cerr<<"here"<<std::endl;
            // }

            //check common mbc
            int nvalidPosR=0;
            for(int j=0;j<v2mbc[anchor][u].size();++j){
                int bc=v2mbc[anchor][u][j];
                if(commonMbcliquePos[bc]<=validPosR){
                    exchange(commonMbclique,commonMbcliquePos[bc],nvalidPosR++,commonMbcliquePos);
                }else break;
            }
            nvalidPosR--;
            if(checkCommonGroup(nvalidPosR)==false) continue;


            RP.push_back(u);
            op_vertices[2].clear();
            if(anchorThreshold<=2){
                for(int j=g->pos[anchor][u];j<g->pos[anchor][u+1];++j){
                    int v=g->e[anchor][j];
                    op_vertices[2].push_back(v);
                }
            }else{
                int j=0,k=g->pos[anchor][u];
                while(j<op_vertices[3].size()&&k<g->pos[anchor][u+1]){
                    if(op_vertices[3][j]==g->e[anchor][k]){
                        op_vertices[2].push_back(op_vertices[3][j]);
                        j++;
                        k++;
                    }else if(op_vertices[3][j]<g->e[anchor][k]){
                        j++;
                    }else{
                        k++;
                    }
                }
            }

            if(op_vertices[2].size()<opThreshold){
                RP.pop_back();
                continue;
            }

            for(int j=g->two_hop_pos[u];j<two_hop_pos_end_layer[2][u];++j){
                int v=g->two_hop_e[j];

                // //debug
                // if(v==1826028||v==1934355||v==1952117){
                //     std::cerr<<"here"<<std::endl;
                // }

                //check common mbc
                int nr=0;
                for(int k=0;k<v2mbc[anchor][v].size();++k){
                    int bc=v2mbc[anchor][v][k];
                    if(commonMbcliquePos[bc]<=nvalidPosR){
                        exchange(commonMbclique,commonMbcliquePos[bc],nr++,commonMbcliquePos);
                    }else if(commonMbcliquePos[bc]>validPosR){// because no shrink
                        break;
                    }
                }
                nr--;
                if(checkCommonGroup(nr)==false) continue;


                RP.push_back(v);
                op_vertices[1].clear();
                int k=0,o=g->pos[anchor][v];
                while(k<op_vertices[2].size()&&o<g->pos[anchor][v+1]){
                    if(op_vertices[2][k]==g->e[anchor][o]){
                        op_vertices[1].push_back(op_vertices[2][k]);
                        k++;
                        o++;
                    }else if(op_vertices[2][k]<g->e[anchor][o]){
                        k++;
                    }else{
                        o++;
                    }
                }
                if(op_vertices[1].size()<opThreshold){
                    RP.pop_back();
                    continue;
                }

                #ifdef COUNTQPONLY
                qpcliqueNum+=count_combinations(op_vertices[1].size(),opThreshold);
                #else
                seed_vertices.clear();
                for(int k=0;k<op_vertices[1].size();++k){
                    seed_vertices.push_back(op_vertices[1][k]);
                }
                uint64_t comb_count = 0;
                // std::vector< std::vector<int> > combinations;
                uint32_t comb_size = opThreshold;
                // vector<int> temp_comb(comb_size);
                // find_combinations(comb_count, 0, comb_size, comb_size);

                //shrink v2mbc
                for(int k=0;k<seed_vertices.size();++k){
                    int w=seed_vertices[k];
                    int tempnr=0;
                    for(int o=0;o<v2mbc[opAnchor][w].size();++o){
                        int bc=v2mbc[opAnchor][w][o];
                        if(commonMbcliquePos[bc]<=nr){
                            exchange(v2mbc[opAnchor][w],o,tempnr++);
                        }
                    }
                }

                //debug
                // if(RP[0]==1826028&&RP[1]==1934355&&RP[2]==1952117){
                //     std::cerr<<"here"<<std::endl;
                // }

                prune_find_combinations(comb_count, 0, comb_size, comb_size,nr);
                #endif

                RP.pop_back();
            }

            RP.pop_back();
        }
        return;
    }

    for(int i=0;i<anchorCandLayer[left].size();++i){
        int u=anchorCandLayer[left][i];

        // //debug
        // if(u==1826028||u==1934355||u==1952117){
        //     std::cerr<<"here"<<std::endl;
        // }

        //check common mbc
        int nvalidPosR=0;
        for(int j=0;j<v2mbc[anchor][u].size();++j){
            int bc=v2mbc[anchor][u][j];
            if(commonMbcliquePos[bc]<=validPosR){
                exchange(commonMbclique,commonMbcliquePos[bc],nvalidPosR++,commonMbcliquePos);
            }else break;
        }
        nvalidPosR--;
        if(checkCommonGroup(nvalidPosR)==false) continue;



        RP.push_back(u);
        op_vertices[left].clear();
        if(left==anchorThreshold){
            for(int j=g->pos[anchor][u];j<g->pos[anchor][u+1];++j){
                int v=g->e[anchor][j];
                op_vertices[left].push_back(v);
            }
        }else{
            int j=0,k=g->pos[anchor][u];
            while(j<op_vertices[left+1].size()&&k<g->pos[anchor][u+1]){
                if(op_vertices[left+1][j]==g->e[anchor][k]){
                    op_vertices[left].push_back(op_vertices[left+1][j]);
                    j++;
                    k++;
                }else if(op_vertices[left+1][j]<g->e[anchor][k]){
                    j++;
                }else{
                    k++;
                }
            }
        }

        if(op_vertices[left].size()<opThreshold){
            RP.pop_back();
            continue;
        }

        anchorCandLayer[left-1].clear();
        for(int j=g->two_hop_pos[u];j<two_hop_pos_end_layer[left][u];++j){
            int v=g->two_hop_e[j];
            if(anchorlab[v]==left){
                anchorlab[v]=left-1;
                anchorCandLayer[left-1].push_back(v);
                two_hop_pos_end_layer[left-1][v]=g->two_hop_pos[v];
            }
        }
        for(int j=0;j<anchorCandLayer[left-1].size();++j){
            int v=anchorCandLayer[left-1][j];
            int end=two_hop_pos_end_layer[left][v];
            for(int k=g->two_hop_pos[v];k<end;++k){
                int w=g->two_hop_e[k];
                if(anchorlab[w]==left-1){
                    two_hop_pos_end_layer[left-1][v]++;
                }else{
                    g->two_hop_e[k--]=g->two_hop_e[--end];
                    g->two_hop_e[end]=w;
                }
            }


            //shrink v2mbc
            int tempnr=0;
            for(int k=0;k<v2mbc[anchor][v].size();++k){
                int bc=v2mbc[anchor][v][k];
                if(commonMbcliquePos[bc]>validPosR) break;
                if(commonMbcliquePos[bc]<=nvalidPosR){
                    exchange(v2mbc[anchor][v],k,tempnr++);
                }
            }
        }

        conductListAndConnectMBC(left-1,nvalidPosR);
        for(int j=0;j<anchorCandLayer[left-1].size();++j){
            int v=anchorCandLayer[left-1][j];
            anchorlab[v]=left;
        }

        RP.pop_back();
    }
}

void BICPC::listAndConnectMBCbase(){
    auto mt1 = std::chrono::steady_clock::now();

    //prepare v2mbc
    v2mbc[0].resize(g->n[0]);
    v2mbc[1].resize(g->n[1]);

    int mbcNum=cliqueuv[0].size();
    for(int i=0;i<mbcNum;++i){
        for(int j=0;j<cliqueuv[0][i].size();++j){
            int u=cliqueuv[0][i][j];
            v2mbc[0][u].push_back(i);
        }
        for(int j=0;j<cliqueuv[1][i].size();++j){
            int v=cliqueuv[1][i][j];
            v2mbc[1][v].push_back(i);
        }
    }

    //prepare common mbclique
    commonMbclique.resize(mbcNum);
    commonMbcliquePos.resize(mbcNum);
    for(int i=0;i<mbcNum;++i){
        commonMbclique[i]=commonMbcliquePos[i]=i;
    }

    //prepare groups
    fans.updateSize(mbcNum);
    visGroup.resize(mbcNum,false);
    visGroups.reserve(mbcNum);

    //prepare original qplist
    qpcliqueNum=0;
    anchor=0;
    if(estimate_cost(0)>estimate_cost(1)) anchor=1;
    anchorThreshold=(anchor==0?p:q);
    opThreshold=(anchor==0?q:p);
    opAnchor=anchor^1;
    std::cerr<<"anchor: "<<anchor<<std::endl;

    g->prepareTowHopNeighbor(anchor,opThreshold);
    two_hop_pos_end_layer.resize(anchorThreshold+1);
    anchorCandLayer.resize(anchorThreshold+1);
    op_vertices.resize(anchorThreshold+1);
    seed_vertices.reserve(g->n[anchor^1]);
    temp_comb.resize(opThreshold);
    anchorlab.resize(g->n[anchor],anchorThreshold);

    for(int i=0;i<=anchorThreshold;++i){
        two_hop_pos_end_layer[i].resize(g->n[anchor]);
        op_vertices[i].reserve((anchor==0?g->maxDu:g->maxDv));
        anchorCandLayer[i].reserve(g->n[anchor]);
    }
    for(int u=0;u<g->n[anchor];++u){
        two_hop_pos_end_layer[anchorThreshold][u]=g->two_hop_pos[u+1];
        anchorCandLayer[anchorThreshold].push_back(u);
    }

    conductListAndConnectMBCbase(anchorThreshold,mbcNum-1);





    auto mt2 = std::chrono::steady_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(mt2 - mt1);
    std::cerr << "list time: " << duration.count() << " ms" << std::endl;

    #ifdef COUNTBICLIQUETREENODE
    std::cerr << "list node num: " << bicliquelistNodecnt << std::endl;
    #endif
}

void BICPC::listAndConnectMBC(){
    auto mt1 = std::chrono::steady_clock::now();

    //prepare v2mbc
    v2mbc[0].resize(g->n[0]);
    v2mbc[1].resize(g->n[1]);

    int mbcNum=cliqueuv[0].size();
    for(int i=0;i<mbcNum;++i){
        for(int j=0;j<cliqueuv[0][i].size();++j){
            int u=cliqueuv[0][i][j];
            v2mbc[0][u].push_back(i);
        }
        for(int j=0;j<cliqueuv[1][i].size();++j){
            int v=cliqueuv[1][i][j];
            v2mbc[1][v].push_back(i);
        }
    }

    //prepare common mbclique
    commonMbclique.resize(mbcNum);
    commonMbcliquePos.resize(mbcNum);
    for(int i=0;i<mbcNum;++i){
        commonMbclique[i]=commonMbcliquePos[i]=i;
    }

    //prepare groups
    fans.updateSize(mbcNum);
    visGroup.resize(mbcNum,false);
    visGroups.reserve(mbcNum);

    //prepare original qplist
    qpcliqueNum=0;
    anchor=0;
    if(estimate_cost(0)>estimate_cost(1)) anchor=1;
    anchorThreshold=(anchor==0?p:q);
    opThreshold=(anchor==0?q:p);
    opAnchor=anchor^1;
    std::cerr<<"anchor: "<<anchor<<std::endl;

    g->prepareTowHopNeighbor(anchor,opThreshold);
    two_hop_pos_end_layer.resize(anchorThreshold+1);
    anchorCandLayer.resize(anchorThreshold+1);
    op_vertices.resize(anchorThreshold+1);
    seed_vertices.reserve(g->n[anchor^1]);
    temp_comb.resize(opThreshold);
    anchorlab.resize(g->n[anchor],anchorThreshold);

    for(int i=0;i<=anchorThreshold;++i){
        two_hop_pos_end_layer[i].resize(g->n[anchor]);
        op_vertices[i].reserve((anchor==0?g->maxDu:g->maxDv));
        anchorCandLayer[i].reserve(g->n[anchor]);
    }
    for(int u=0;u<g->n[anchor];++u){
        two_hop_pos_end_layer[anchorThreshold][u]=g->two_hop_pos[u+1];
        anchorCandLayer[anchorThreshold].push_back(u);
    }

    conductListAndConnectMBC(anchorThreshold,mbcNum-1);

    auto mt2 = std::chrono::steady_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(mt2 - mt1);
    std::cerr << "list time: " << duration.count() << " ms" << std::endl;
    #ifdef COUNTBICLIQUETREENODE
    std::cerr << "list node num: " << bicliquelistNodecnt << std::endl;
    #endif

}

uint64_t BICPC::count_combinations(uint64_t _n,uint64_t _k){
    if(combination_cache[_n][_k] > 0){
        return combination_cache[_n][_k];
    }
    else{
        if (_k > _n) {
            return 0;
        }
        uint64_t backup_n = _n;
        uint64_t backup_k = _k;
        uint64_t r = 1;
        for (uint64_t d = 1; d <= _k; ++d) {
            r *= _n--;
            r /= d;
        }
        if(_n <= MAX_N && _k <= MAX_K){
            combination_cache[backup_n][backup_k] = r;
        }
        return r;
    }
}

void BICPC::find_combinations(uint64_t &comb_count, uint32_t beg_offset, uint32_t curr_depth, uint32_t total_depth){
    uint32_t N = seed_vertices.size();
    if (curr_depth == 0) {
        cliqueuv[anchor^1].emplace_back(temp_comb);
        cliqueuv[anchor].emplace_back(RP);
        // combs.emplace_back(comb);
        comb_count++;
        return;
    }
    for (uint32_t i = beg_offset; i < N; i++){
        temp_comb[total_depth-curr_depth] = seed_vertices[i];
        find_combinations(comb_count, i + 1, curr_depth - 1, total_depth);
    }
}

double BICPC::estimate_cost(int side){
    srand(time(NULL));
    int n1=g->n[0];
    int n2=g->n[1];
    int opside=side^1;

    int num_vertices=n1+n2;
    int vertices_cap=std::max(n1,n2);
    int num_rounds=ceil(num_vertices*0.01);

    uint64_t total_two_hop_deg=0;
    uint64_t max_two_hop_deg=0;

    std::vector<uint32_t> common_neig_map, aux_array_two_neig;
    common_neig_map.resize(vertices_cap,0);
    aux_array_two_neig.resize(vertices_cap);
    // uint32_t offset=(side==0?0:n1);
    uint32_t common_neig_threshold=(side==0?q:p);

    for(int r=0;r<num_rounds;++r){
        uint64_t estimated_two_hop_deg=0;
        uint32_t i=rand()%g->n[side];
        uint32_t idx=0;
        for(int j=g->pos[side][i];j<g->pos[side][i+1];++j){
            int v=g->e[side][j];
            for(int k=g->pos[opside][v];k<g->pos[opside][v+1];++k){
                int w=g->e[opside][k];
                if(w<i){
                    common_neig_map[w]++;
                    if(common_neig_map[w]==1){
                        aux_array_two_neig[idx++]=w;
                    }
                }else{
                    break;
                }
            }
        }

        for(int j=0;j<idx;++j){
            int w=aux_array_two_neig[j];
            if(common_neig_map[w]>=common_neig_threshold){
                estimated_two_hop_deg++;
            }
            common_neig_map[w]=0;
        }

        max_two_hop_deg=std::max(max_two_hop_deg,estimated_two_hop_deg);
        total_two_hop_deg+=estimated_two_hop_deg*g->n[side];
    }

    total_two_hop_deg/=num_rounds;
    uint64_t avg_two_hop_deg=std::max((uint64_t)2,total_two_hop_deg/g->n[side]);
    uint64_t pq_value=(side==0?p:q);
    double totalCost=total_two_hop_deg*std::pow(avg_two_hop_deg,pq_value-2);
    return totalCost;
}

void BICPC::listQPCliques(int left){
    if(left==2){
        for(int i=0;i<anchorCandLayer[left].size();++i){
            int u=anchorCandLayer[left][i];
            RP.push_back(u);
            op_vertices[2].clear();
            if(anchorThreshold<=2){
                for(int j=g->pos[anchor][u];j<g->pos[anchor][u+1];++j){
                    int v=g->e[anchor][j];
                    op_vertices[2].push_back(v);
                }
            }else{
                int j=0,k=g->pos[anchor][u];
                while(j<op_vertices[3].size()&&k<g->pos[anchor][u+1]){
                    if(op_vertices[3][j]==g->e[anchor][k]){
                        op_vertices[2].push_back(op_vertices[3][j]);
                        j++;
                        k++;
                    }else if(op_vertices[3][j]<g->e[anchor][k]){
                        j++;
                    }else{
                        k++;
                    }
                }
            }

            if(op_vertices[2].size()<opThreshold){
                RP.pop_back();
                continue;
            }

            for(int j=g->two_hop_pos[u];j<two_hop_pos_end_layer[2][u];++j){
                int v=g->two_hop_e[j];
                RP.push_back(v);
                op_vertices[1].clear();
                int k=0,o=g->pos[anchor][v];
                while(k<op_vertices[2].size()&&o<g->pos[anchor][v+1]){
                    if(op_vertices[2][k]==g->e[anchor][o]){
                        op_vertices[1].push_back(op_vertices[2][k]);
                        k++;
                        o++;
                    }else if(op_vertices[2][k]<g->e[anchor][o]){
                        k++;
                    }else{
                        o++;
                    }
                }
                if(op_vertices[1].size()<opThreshold){
                    RP.pop_back();
                    continue;
                }

                #ifdef COUNTQPONLY
                qpcliqueNum+=count_combinations(op_vertices[1].size(),opThreshold);
                #else
                seed_vertices.clear();
                for(int k=0;k<op_vertices[1].size();++k){
                    seed_vertices.push_back(op_vertices[1][k]);
                }
                uint64_t comb_count = 0;
                // std::vector< std::vector<int> > combinations;
                uint32_t comb_size = opThreshold;
                // vector<int> temp_comb(comb_size);
                find_combinations(comb_count, 0, comb_size, comb_size);
                #endif

                RP.pop_back();
            }

            RP.pop_back();
        }
        return;
    }

    for(int i=0;i<anchorCandLayer[left].size();++i){
        int u=anchorCandLayer[left][i];
        RP.push_back(u);
        op_vertices[left].clear();
        if(left==anchorThreshold){
            for(int j=g->pos[anchor][u];j<g->pos[anchor][u+1];++j){
                int v=g->e[anchor][j];
                op_vertices[left].push_back(v);
            }
        }else{
            int j=0,k=g->pos[anchor][u];
            while(j<op_vertices[left+1].size()&&k<g->pos[anchor][u+1]){
                if(op_vertices[left+1][j]==g->e[anchor][k]){
                    op_vertices[left].push_back(op_vertices[left+1][j]);
                    j++;
                    k++;
                }else if(op_vertices[left+1][j]<g->e[anchor][k]){
                    j++;
                }else{
                    k++;
                }
            }
        }

        if(op_vertices[left].size()<opThreshold){
            RP.pop_back();
            continue;
        }

        anchorCandLayer[left-1].clear();
        for(int j=g->two_hop_pos[u];j<two_hop_pos_end_layer[left][u];++j){
            int v=g->two_hop_e[j];
            if(anchorlab[v]==left){
                anchorlab[v]=left-1;
                anchorCandLayer[left-1].push_back(v);
                two_hop_pos_end_layer[left-1][v]=g->two_hop_pos[v];
            }
        }
        for(int j=0;j<anchorCandLayer[left-1].size();++j){
            int v=anchorCandLayer[left-1][j];
            int end=two_hop_pos_end_layer[left][v];
            for(int k=g->two_hop_pos[v];k<end;++k){
                int w=g->two_hop_e[k];
                if(anchorlab[w]==left-1){
                    two_hop_pos_end_layer[left-1][v]++;
                }else{
                    g->two_hop_e[k--]=g->two_hop_e[--end];
                    g->two_hop_e[end]=w;
                }
            }
        }

        listQPCliques(left-1);
        for(int j=0;j<anchorCandLayer[left-1].size();++j){
            int v=anchorCandLayer[left-1][j];
            anchorlab[v]=left;
        }

        RP.pop_back();
    }
}

void BICPC::prepareAndListqpCliques(){
    qpcliqueNum=0;
    anchor=0;
    if(estimate_cost(0)>estimate_cost(1)) anchor=1;
    anchorThreshold=(anchor==0?p:q);
    opThreshold=(anchor==0?q:p);
    std::cerr<<"anchor: "<<anchor<<std::endl;

    g->prepareTowHopNeighbor(anchor,opThreshold);
    two_hop_pos_end_layer.resize(anchorThreshold+1);
    anchorCandLayer.resize(anchorThreshold+1);
    op_vertices.resize(anchorThreshold+1);
    seed_vertices.reserve(g->n[anchor^1]);
    temp_comb.resize(opThreshold);
    anchorlab.resize(g->n[anchor],anchorThreshold);

    for(int i=0;i<=anchorThreshold;++i){
        two_hop_pos_end_layer[i].resize(g->n[anchor]);
        op_vertices[i].reserve((anchor==0?g->maxDu:g->maxDv));
        anchorCandLayer[i].reserve(g->n[anchor]);
    }
    for(int u=0;u<g->n[anchor];++u){
        two_hop_pos_end_layer[anchorThreshold][u]=g->two_hop_pos[u+1];
        anchorCandLayer[anchorThreshold].push_back(u);
    }

    listQPCliques(anchorThreshold);
}

int BICPC::getGroupNum(){
    int bicliqueNum=cliqueuv[0].size();
    int groupnum=0;
    fans.updateSize(bicliqueNum);
    for(int i=0;i<bicliqueNum;++i){
        int father=fans.find(i);
        if(father==i) groupnum++;
    }

    return groupnum;

}

int BICPC::connectAllBiCliqueLeafResbase(){

    for(int i=S[firstlayert].x;i<S[firstlayert].c;++i){
        if(runtree.firstlayerson[S[firstlayert][i]]!=-1){
            int startnode=runtree.firstlayerson[S[firstlayert][i]];
            collectBiCliqueLeafResbase(startnode,connectResults);
            if(RUVlook[firstlayert][runtree.treeNodePool[startnode].inValue]==1){
                std::cerr<<"possible error"<<std::endl;
                exit(0);
            }
        }
    }

    for(int i=0;i<xStackNode.size();++i){
        int inValue=runtree.treeNodePool[xStackNode[i]].inValue;
        int t=runtree.treeNodePool[xStackNode[i]].uv;

        int pvalue=S[t].pos(inValue);
        if((pvalue>=S[t].x&&pvalue<S[t].r)||(RUVlook[t][inValue]==1)){
            collectBiCliqueLeafResbase(xStackNode[i],connectResults);
            if(RUVlook[t][inValue]==1){
                std::cerr<<"possible error"<<std::endl;
                exit(0);
            }
        }
    }

    if(connectResults.size()==0){
        std::cerr<<"connect result empty"<<std::endl;
        exit(0);
    }

    fans.updateSize(cliqueuv[0].size());
    if(connectResults.size()>1){
        for(int i=1;i<connectResults.size();++i){
            fans.merge(connectResults[0],connectResults[i]);
        }
    }

    int res=fans.find(connectResults[0]);
    connectResults.clear();
    return res;
}

int BICPC::connectAllBiCliqueLeafRes(){

    for(int i=S[firstlayert].x;i<S[firstlayert].c;++i){
        if(runtree.firstlayerson[S[firstlayert][i]]!=-1){
            int startnode=runtree.firstlayerson[S[firstlayert][i]];
            collectBiCliqueLeafRes(startnode,connectResults);
            if(RUVlook[firstlayert][runtree.treeNodePool[startnode].inValue]==1){
                std::cerr<<"possible error"<<std::endl;
                exit(0);
            }
        }
    }

    for(int i=0;i<xStackNode.size();++i){
        int inValue=runtree.treeNodePool[xStackNode[i]].inValue;
        int t=runtree.treeNodePool[xStackNode[i]].uv;

        int pvalue=S[t].pos(inValue);
        if((pvalue>=S[t].x&&pvalue<S[t].r)||(RUVlook[t][inValue]==1)){
            collectBiCliqueLeafRes(xStackNode[i],connectResults);
            if(RUVlook[t][inValue]==1){
                std::cerr<<"possible error"<<std::endl;
                exit(0);
            }
        }
    }

    if(connectResults.size()==0){
        std::cerr<<"connect result empty"<<std::endl;
        exit(0);
    }

    fans.updateSize(cliqueuv[0].size());
    if(connectResults.size()>1){
        for(int i=1;i<connectResults.size();++i){
            fans.merge(connectResults[0],connectResults[i]);
        }
    }

    int res=fans.find(connectResults[0]);
    connectResults.clear();
    return res;
}

int BICPC::connectBiCliqueNewNoX(int side,int thisNode){
    int t=side;
    int z=side^1;

    if(thisNode==-1) return -1;

    if(runtree.treeNodePool[thisNode].resMark!=-1){
        connectResults.push_back(runtree.treeNodePool[thisNode].resMark);
    }

    int res=connectBiCliqueInCurZone(side);
    if(res!=-1){
        connectResults.push_back(res);
    }

    for(int i=0;i<runtree.treeNodePool[thisNode].sons.size();++i){
        int sonNode=runtree.treeNodePool[thisNode].sons[i];
        connectResults.push_back(runtree.treeNodePool[sonNode].resMark);
    }

    if(connectResults.size()==0){
        std::cerr<<"no results error"<<std::endl;
        exit(0);
    }

    fans.updateSize(cliqueuv[0].size());
    if(connectResults.size()>1){
        for(int i=1;i<connectResults.size();++i){
            fans.merge(connectResults[0],connectResults[i]);
        }
    }

    res=fans.find(connectResults[0]);
    connectResults.clear();
    return res;
}

int BICPC::connectAllBiCliqueNew(int side,int thisNode){
    int t=side;
    int z=side^1;
    bool needX=false;

    if(thisNode!=-1&&runtree.treeNodePool[thisNode].resMark!=-1){
        connectResults.push_back(runtree.treeNodePool[thisNode].resMark);
    }

    int res=connectBiCliqueInCurZone(side);
    if(res!=-1){
        connectResults.push_back(res);
    }

    #ifdef XZONEDIRECT
    needX=true;
    #endif

    if(needX==false){
        for(uint32_t i = S[t].x; i < S[t].r; i++) {
            for(uint32_t j = S[z].x; j < S[z].r; j++) {
                if(g->connect(S[t][i], S[z][j], t)) {
                    needX=true;
                    break;
                }
            }
            if(needX){
                break;
            }
        }
    }
    
    if(needX) connectBiCliqueInXzone(side);

    if(connectResults.size()==0){
        std::cerr<<"final connect id error"<<std::endl;
        exit(0);
        // return -1;
    }

    fans.updateSize(cliqueuv[0].size());
    if(connectResults.size()>1){
        for(int i=1;i<connectResults.size();++i){
            fans.merge(connectResults[0],connectResults[i]);
        }
    }

    return fans.find(connectResults[0]);
}

int BICPC::connectAllBiClique(int side,int thisNode,bool first){
    int t=side;
    int z=side^1;
    bool needX=false;

    if(thisNode!=-1&&runtree.treeNodePool[thisNode].resMark!=-1){
        connectResults.push_back(runtree.treeNodePool[thisNode].resMark);
    }else{
        needX=true;
    }

    int res=connectBiCliqueInCurZone(side);
    if(res!=-1){
        connectResults.push_back(res);
    }else{
        needX=true;
    }

    #ifdef XZONEDIRECT
    needX=true;
    #endif

    if(needX==false){
        if(first==false){
            needX=true;
        }else{
            for(uint32_t i = S[t].x; i < S[t].r; i++) {
                for(uint32_t j = S[z].x; j < S[z].r; j++) {
                    if(g->connect(S[t][i], S[z][j], t)) {
                        needX=true;
                        break;
                    }
                }
                if(needX){
                    break;
                }
            }
        }
    }
    
    if(needX) connectBiCliqueInXzone(side);

    if(connectResults.size()==0){
        std::cerr<<"final connect id error"<<std::endl;
        exit(0);
        // return -1;
    }

    fans.updateSize(cliqueuv[0].size());
    if(connectResults.size()>1){
        for(int i=1;i<connectResults.size();++i){
            fans.merge(connectResults[0],connectResults[i]);
        }
    }

    return fans.find(connectResults[0]);

}

int BICPC::connectBiCliqueInCurZone(int side){
    for(int i=curStackNode.size()-1;i>=0;--i){
        int cliqueId=curStackNode[i];
        if(cliqueuv[side^1][cliqueId].size()<RUV[side^1].size()){
            return -1;
        }
        if(cliqueuv[side^1][cliqueId].size()==RUV[side^1].size()){
            return cliqueId;
        }
    }
    return -1;
}

void BICPC::quasi2finalBCPC(){
    auto mt1 = std::chrono::steady_clock::now();

    int bicliqueNum=cliqueuv[0].size();
    std::vector<std::vector<int>> fa2list;
    fa2list.resize(bicliqueNum);
    std::vector<int> falist;

    fans.updateSize(bicliqueNum);
    for(int c=0;c<bicliqueNum;++c){
        int father=fans.find(c);
        fa2list[father].push_back(c);
        if(father==c) falist.push_back(c);
    }

    std::cerr<<"group num: "<<falist.size()<<std::endl;

    std::vector<std::vector<int>> node2bicliques[2];
    node2bicliques[0].resize(g->n[0]+1);
    node2bicliques[1].resize(g->n[1]+1);
    std::vector<int> clique2group;
    clique2group.resize(bicliqueNum);

    //prepare node to biclique
    for(int i=0;i<falist.size();++i){
        int father=falist[i];
        for(int j=0;j<fa2list[father].size();++j){
            int c=fa2list[father][j];
            clique2group[c]=i;
            for(int k=0;k<cliqueuv[0][c].size();++k){
                int u=cliqueuv[0][c][k];
                node2bicliques[0][u].push_back(c);
            }
            for(int k=0;k<cliqueuv[1][c].size();++k){
                int v=cliqueuv[1][c][k];
                node2bicliques[1][v].push_back(c);
            }
        }
    }

    //prepare jump array
    std::vector<std::vector<int>> node2bicliqueJump[2];
    node2bicliqueJump[0].resize(g->n[0]+1);
    node2bicliqueJump[1].resize(g->n[1]+1);
    for(int u=0;u<g->n[0];++u){
        int next=node2bicliques[0][u].size();
        if(next==0) continue;

        node2bicliqueJump[0][u].resize(next);
        int curGroup=clique2group[node2bicliques[0][u][next-1]];
        for(int j=node2bicliques[0][u].size()-1;j>=0;--j){
            int c=node2bicliques[0][u][j];
            if(clique2group[c]==curGroup) node2bicliqueJump[0][u][j]=next;
            else{
                curGroup=clique2group[c];
                next=j+1;
                node2bicliqueJump[0][u][j]=next;
            }
        }
    }
    for(int v=0;v<g->n[1];++v){
        int next=node2bicliques[1][v].size();
        if(next==0) continue;

        node2bicliqueJump[1][v].resize(next);
        int curGroup=clique2group[node2bicliques[1][v][next-1]];
        for(int j=node2bicliques[1][v].size()-1;j>=0;--j){
            int c=node2bicliques[1][v][j];
            if(clique2group[c]==curGroup) node2bicliqueJump[1][v][j]=next;
            else{
                curGroup=clique2group[c];
                next=j+1;
                node2bicliqueJump[1][v][j]=next;
            }
        }
    }

    std::vector<int> overlap[2];
    std::vector<int> computed;

    overlap[0].resize(bicliqueNum,0);
    overlap[1].resize(bicliqueNum,0);
    computed.reserve(bicliqueNum);

    for(int c=0;c<bicliqueNum;++c){
        int curGroup=clique2group[c];
        for(int j=0;j<cliqueuv[0][c].size();++j){
            int u=cliqueuv[0][c][j];
            for(int k=0;k<node2bicliques[0][u].size();){
                int nc=node2bicliques[0][u][k];
                if(clique2group[nc]==curGroup){
                    k=node2bicliqueJump[0][u][k];
                }else if(fans.find(nc)==fans.find(c)){
                    k=node2bicliqueJump[0][u][k];
                }else{
                    if(overlap[0][nc]++ == 0){
                        computed.push_back(nc);
                    }
                    k++;
                }
            }
        }

        for(int j=0;j<cliqueuv[1][c].size();++j){
            int v=cliqueuv[1][c][j];
            for(int k=0;k<node2bicliques[1][v].size();){
                int nc=node2bicliques[1][v][k];
                if(clique2group[nc]==curGroup){
                    k=node2bicliqueJump[1][v][k];
                }else if(fans.find(nc)==fans.find(c)){
                    k=node2bicliqueJump[1][v][k];
                }else if(overlap[0][nc]<p){
                    k++;
                }else{
                    overlap[1][nc]++;
                    if(overlap[1][nc]>=q){
                        fans.merge(nc,c);
                        k=node2bicliqueJump[1][v][k];
                    }else{
                        k++;
                    }
                }
            }
        }

        for(int j=0;j<computed.size();++j){
            int c=computed[j];
            overlap[0][c]=overlap[1][c]=0;
        }
        computed.clear();
    }

    auto mt2 = std::chrono::steady_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(mt2 - mt1);
    std::cerr << "traverse time: " << duration.count() << " ms" << std::endl;
}

void BICPC::outputBCPCs(){
    //new labels to old labels
    std::vector<uint32_t> * theOldLabels=g->old_lables;
    for(int i=0;i<cliqueuv[0].size();++i){
        for(int j=0;j<cliqueuv[0][i].size();++j){
            cliqueuv[0][i][j]=theOldLabels[0][cliqueuv[0][i][j]];
        }
    }
    for(int i=0;i<cliqueuv[1].size();++i){
        for(int j=0;j<cliqueuv[1][i].size();++j){
            cliqueuv[1][i][j]=theOldLabels[1][cliqueuv[1][i][j]];
        }
    }

    for(int i=0;i<cliqueuv[0].size();++i){
        sort(cliqueuv[0][i].begin(),cliqueuv[0][i].end());
    }
    for(int i=0;i<cliqueuv[1].size();++i){
        sort(cliqueuv[1][i].begin(),cliqueuv[1][i].end());
    }

    // std::vector<std::vector<int>> bicpcResults;
    int bicliqueNum=cliqueuv[0].size();

    fans.updateSize(bicliqueNum);
    std::unordered_map<int,std::vector<int> > bicpcResults;
    for(int i=0;i<bicliqueNum;++i){
        if(cliqueuv[0][i].size()==p&&cliqueuv[1][i].size()==q) continue;
        int father=fans.find(i);
        bicpcResults[father].push_back(i);
    }

    

    std::vector<std::vector<int>> results;
    std::unordered_map<int,std::vector<int> >::iterator it=bicpcResults.begin();
    while(it!=bicpcResults.end()){
        results.push_back(it->second);
        it++;
    }

    for(int i=0;i<results.size();++i){
        sort(results[i].begin(),results[i].end());
    }
    sort(results.begin(),results.end());


    //final print
    for(int i=0;i<results.size();++i){
        std::cout<<"qpcpc:"<<std::endl;
        for(int j=0;j<results[i].size();++j){
            int c=results[i][j];

            for(int k=0;k<cliqueuv[0][c].size();++k){
                std::cout<<cliqueuv[0][c][k]<<" ";
            }
            std::cout<<std::endl;
            for(int k=0;k<cliqueuv[1][c].size();++k){
                std::cout<<cliqueuv[1][c][k]<<" ";
            }
            std::cout<<std::endl;
        }

    }


}

void BICPC::qbcpcLeafRespoorCntpqnode(){
    auto mt1 = std::chrono::steady_clock::now();

    uint32_t t=0, z=1;
    if(g->maxDu<g->maxDv) {t=1,z=0;}
    firstlayert=t;

    // std::cerr<<"t: "<<t<<" "<<"z: "<<z<<std::endl;
    // std::cerr<<"missing clique:"<<std::endl;
    // std::cerr<<g->labelsL[1]<<" "<<g->labelsL[184683]<<std::endl;
    // std::cerr<<g->labelsR[133490]<<std::endl;

    runtree.root=runtree.getNewNode();
    runtree.treeNodePool[runtree.root].inValue=-1;
    runtree.treeNodePool[runtree.root].resMark=-1;
    runtree.treeNodePool[runtree.root].qpMark=false;
    runtree.treeNodePool[runtree.root].uv=-1;

    runtree.firstlayerson.resize(g->n[t],-1);
    xStackNode.reserve(g->n[0]+g->n[1]);

    //first layer no pivot
    for(uint32_t u=0;u<g->n[t];u++){
        S[t].c=S[t].r=g->n[t];
        S[z].c=S[z].r=g->n[z];

        

        //form C set for z
        for(uint32_t i=g->pos[t][u];i<g->pos[t][u+1];++i){
            uint32_t v=g->e[t][i];
            S[z].swapByPos(--S[z].c,S[z].pos(v));
        }

        S[z].x=S[z].c;

        //form C set for t
        for(uint32_t i = g->pos[t][u]; i < g->pos[t][u + 1]; i++) {
            uint32_t v = g->e[t][i];
            if(g->pos[z][v + 1] > 0){
                for(uint32_t j = g->pos[z][v + 1] - 1; j >= g->pos[z][v]; j--) {
                    uint32_t w = g->e[z][j];
                    
                    if(w > u) {
                        uint32_t pw = S[t].pos(w);
                        if(pw < S[t].c) {
                            S[t].swapByPos(--S[t].c, pw);
                        }
                    }
                    else break;

                    if(j == 0) break;
                }
            }
        }
        S[t].x = S[t].c;

        //形成t的X集合
        for(uint32_t i = g->pos[t][u]; i < g->pos[t][u + 1]; i++) {
            uint32_t v = g->e[t][i];
            if(g->pos[z][v + 1] > 0){
                for(uint32_t j = g->pos[z][v]; j < g->pos[z][v + 1]; j++) {
                    uint32_t w = g->e[z][j];
                    
                    if(w < u) {
                        uint32_t pw = S[t].pos(w);
                        if(pw < S[t].x) {
                            S[t].swapByPos(--S[t].x, pw);
                        }
                    }
                    else break;

                    if(j == 0) break;
                }
            }
        }


        if(S[z].CIsEmpty()) continue;

        RUV[t].push_back(u);
        RUVlook[t][u]=1;
        // int res=bbranch(0,t,u);
        int res=bbranchLeafRespoorcntpqnode(0,t,u);
        RUV[t].pop_back();
        RUVlook[t][u]=0;

        if(res!=-1){
            runtree.firstlayerson[u]=res;
            runtree.treeNodePool[runtree.root].sons.push_back(res);
        }

    }

    std::cerr<<"maxBiCliqueCount: "<<cliqueuv[0].size()<<std::endl;
    std::cerr<<"alpha,beta node cnt: "<<pqnodecnt<<std::endl;
    auto mt2 = std::chrono::steady_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(mt2 - mt1);
    std::cerr << "qbcpc time: " << duration.count() << " ms" << std::endl;
}

void BICPC::qbcpcLeafRespoor(){
    auto mt1 = std::chrono::steady_clock::now();

    uint32_t t=0, z=1;
    if(g->maxDu<g->maxDv) {t=1,z=0;}
    firstlayert=t;

    // std::cerr<<"t: "<<t<<" "<<"z: "<<z<<std::endl;
    // std::cerr<<"missing clique:"<<std::endl;
    // std::cerr<<g->labelsL[1]<<" "<<g->labelsL[184683]<<std::endl;
    // std::cerr<<g->labelsR[133490]<<std::endl;

    runtree.root=runtree.getNewNode();
    runtree.treeNodePool[runtree.root].inValue=-1;
    runtree.treeNodePool[runtree.root].resMark=-1;
    runtree.treeNodePool[runtree.root].qpMark=false;
    runtree.treeNodePool[runtree.root].uv=-1;

    runtree.firstlayerson.resize(g->n[t],-1);
    xStackNode.reserve(g->n[0]+g->n[1]);

    //first layer no pivot
    for(uint32_t u=0;u<g->n[t];u++){
        S[t].c=S[t].r=g->n[t];
        S[z].c=S[z].r=g->n[z];

        

        //form C set for z
        for(uint32_t i=g->pos[t][u];i<g->pos[t][u+1];++i){
            uint32_t v=g->e[t][i];
            S[z].swapByPos(--S[z].c,S[z].pos(v));
        }

        S[z].x=S[z].c;

        //form C set for t
        for(uint32_t i = g->pos[t][u]; i < g->pos[t][u + 1]; i++) {
            uint32_t v = g->e[t][i];
            if(g->pos[z][v + 1] > 0){
                for(uint32_t j = g->pos[z][v + 1] - 1; j >= g->pos[z][v]; j--) {
                    uint32_t w = g->e[z][j];
                    
                    if(w > u) {
                        uint32_t pw = S[t].pos(w);
                        if(pw < S[t].c) {
                            S[t].swapByPos(--S[t].c, pw);
                        }
                    }
                    else break;

                    if(j == 0) break;
                }
            }
        }
        S[t].x = S[t].c;

        //形成t的X集合
        for(uint32_t i = g->pos[t][u]; i < g->pos[t][u + 1]; i++) {
            uint32_t v = g->e[t][i];
            if(g->pos[z][v + 1] > 0){
                for(uint32_t j = g->pos[z][v]; j < g->pos[z][v + 1]; j++) {
                    uint32_t w = g->e[z][j];
                    
                    if(w < u) {
                        uint32_t pw = S[t].pos(w);
                        if(pw < S[t].x) {
                            S[t].swapByPos(--S[t].x, pw);
                        }
                    }
                    else break;

                    if(j == 0) break;
                }
            }
        }


        if(S[z].CIsEmpty()) continue;

        RUV[t].push_back(u);
        RUVlook[t][u]=1;
        // int res=bbranch(0,t,u);
        int res=bbranchLeafRespoor(0,t,u);
        RUV[t].pop_back();
        RUVlook[t][u]=0;

        if(res!=-1){
            runtree.firstlayerson[u]=res;
            runtree.treeNodePool[runtree.root].sons.push_back(res);
        }

    }

    std::cerr<<"maxBiCliqueCount: "<<cliqueuv[0].size()<<std::endl;
    auto mt2 = std::chrono::steady_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(mt2 - mt1);
    std::cerr << "qbcpc time: " << duration.count() << " ms" << std::endl;
}

void BICPC::qbcpcLeafResbase(){
    auto mt1 = std::chrono::steady_clock::now();

    uint32_t t=0, z=1;
    if(g->maxDu<g->maxDv) {t=1,z=0;}
    firstlayert=t;

    // std::cerr<<"t: "<<t<<" "<<"z: "<<z<<std::endl;
    // std::cerr<<"missing clique:"<<std::endl;
    // std::cerr<<g->labelsL[1]<<" "<<g->labelsL[184683]<<std::endl;
    // std::cerr<<g->labelsR[133490]<<std::endl;

    runtree.root=runtree.getNewNode();
    runtree.treeNodePool[runtree.root].inValue=-1;
    runtree.treeNodePool[runtree.root].resMark=-1;
    runtree.treeNodePool[runtree.root].qpMark=false;
    runtree.treeNodePool[runtree.root].uv=-1;

    runtree.firstlayerson.resize(g->n[t],-1);
    xStackNode.reserve(g->n[0]+g->n[1]);

    //first layer no pivot
    for(uint32_t u=0;u<g->n[t];u++){
        S[t].c=S[t].r=g->n[t];
        S[z].c=S[z].r=g->n[z];

        

        //form C set for z
        for(uint32_t i=g->pos[t][u];i<g->pos[t][u+1];++i){
            uint32_t v=g->e[t][i];
            S[z].swapByPos(--S[z].c,S[z].pos(v));
        }

        S[z].x=S[z].c;

        //form C set for t
        for(uint32_t i = g->pos[t][u]; i < g->pos[t][u + 1]; i++) {
            uint32_t v = g->e[t][i];
            if(g->pos[z][v + 1] > 0){
                for(uint32_t j = g->pos[z][v + 1] - 1; j >= g->pos[z][v]; j--) {
                    uint32_t w = g->e[z][j];
                    
                    if(w > u) {
                        uint32_t pw = S[t].pos(w);
                        if(pw < S[t].c) {
                            S[t].swapByPos(--S[t].c, pw);
                        }
                    }
                    else break;

                    if(j == 0) break;
                }
            }
        }
        S[t].x = S[t].c;

        //形成t的X集合
        for(uint32_t i = g->pos[t][u]; i < g->pos[t][u + 1]; i++) {
            uint32_t v = g->e[t][i];
            if(g->pos[z][v + 1] > 0){
                for(uint32_t j = g->pos[z][v]; j < g->pos[z][v + 1]; j++) {
                    uint32_t w = g->e[z][j];
                    
                    if(w < u) {
                        uint32_t pw = S[t].pos(w);
                        if(pw < S[t].x) {
                            S[t].swapByPos(--S[t].x, pw);
                        }
                    }
                    else break;

                    if(j == 0) break;
                }
            }
        }


        if(S[z].CIsEmpty()) continue;

        RUV[t].push_back(u);
        RUVlook[t][u]=1;
        // int res=bbranch(0,t,u);
        int res=bbranchLeafResbase(0,t,u);
        RUV[t].pop_back();
        RUVlook[t][u]=0;

        if(res!=-1){
            runtree.firstlayerson[u]=res;
            runtree.treeNodePool[runtree.root].sons.push_back(res);
        }

    }

    std::cerr<<"maxBiCliqueCount: "<<cliqueuv[0].size()<<std::endl;
    auto mt2 = std::chrono::steady_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(mt2 - mt1);
    std::cerr << "qbcpc time: " << duration.count() << " ms" << std::endl;
}

void BICPC::qbcpcLeafResforabnode(){
    auto mt1 = std::chrono::steady_clock::now();

    uint32_t t=0, z=1;
    if(g->maxDu<g->maxDv) {t=1,z=0;}
    firstlayert=t;

    // std::cerr<<"t: "<<t<<" "<<"z: "<<z<<std::endl;
    // std::cerr<<"missing clique:"<<std::endl;
    // std::cerr<<g->labelsL[1]<<" "<<g->labelsL[184683]<<std::endl;
    // std::cerr<<g->labelsR[133490]<<std::endl;

    runtree.root=runtree.getNewNode();
    runtree.treeNodePool[runtree.root].inValue=-1;
    runtree.treeNodePool[runtree.root].resMark=-1;
    runtree.treeNodePool[runtree.root].qpMark=false;
    runtree.treeNodePool[runtree.root].uv=-1;

    runtree.firstlayerson.resize(g->n[t],-1);
    xStackNode.reserve(g->n[0]+g->n[1]);

    //first layer no pivot
    for(uint32_t u=0;u<g->n[t];u++){
        S[t].c=S[t].r=g->n[t];
        S[z].c=S[z].r=g->n[z];

        

        //form C set for z
        for(uint32_t i=g->pos[t][u];i<g->pos[t][u+1];++i){
            uint32_t v=g->e[t][i];
            S[z].swapByPos(--S[z].c,S[z].pos(v));
        }

        S[z].x=S[z].c;

        //form C set for t
        for(uint32_t i = g->pos[t][u]; i < g->pos[t][u + 1]; i++) {
            uint32_t v = g->e[t][i];
            if(g->pos[z][v + 1] > 0){
                for(uint32_t j = g->pos[z][v + 1] - 1; j >= g->pos[z][v]; j--) {
                    uint32_t w = g->e[z][j];
                    
                    if(w > u) {
                        uint32_t pw = S[t].pos(w);
                        if(pw < S[t].c) {
                            S[t].swapByPos(--S[t].c, pw);
                        }
                    }
                    else break;

                    if(j == 0) break;
                }
            }
        }
        S[t].x = S[t].c;

        //形成t的X集合
        for(uint32_t i = g->pos[t][u]; i < g->pos[t][u + 1]; i++) {
            uint32_t v = g->e[t][i];
            if(g->pos[z][v + 1] > 0){
                for(uint32_t j = g->pos[z][v]; j < g->pos[z][v + 1]; j++) {
                    uint32_t w = g->e[z][j];
                    
                    if(w < u) {
                        uint32_t pw = S[t].pos(w);
                        if(pw < S[t].x) {
                            S[t].swapByPos(--S[t].x, pw);
                        }
                    }
                    else break;

                    if(j == 0) break;
                }
            }
        }


        if(S[z].CIsEmpty()) continue;

        RUV[t].push_back(u);
        RUVlook[t][u]=1;
        // int res=bbranch(0,t,u);
        int res=bbranchLeafResforabnode(0,t,u);
        RUV[t].pop_back();
        RUVlook[t][u]=0;

        if(res!=-1){
            runtree.firstlayerson[u]=res;
            runtree.treeNodePool[runtree.root].sons.push_back(res);
        }

    }

    std::cerr<<"maxBiCliqueCount: "<<cliqueuv[0].size()<<std::endl;
    auto mt2 = std::chrono::steady_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(mt2 - mt1);
    std::cerr << "qbcpc time: " << duration.count() << " ms" << std::endl;
    std::cerr<<"real abnode: "<<pqnodecnt<<std::endl;
    std::cerr<<"virtual abnode: "<<vpqnodecnt<<std::endl;
    std::cerr<<"total abnode: "<<pqnodecnt+vpqnodecnt<<std::endl;
    std::cerr<<"ave abnode: "<<pqnodedepth_total/pqnodecnt<<std::endl;
}

void BICPC::qbcpcLeafRes(){
    auto mt1 = std::chrono::steady_clock::now();

    uint32_t t=0, z=1;
    if(g->maxDu<g->maxDv) {t=1,z=0;}
    firstlayert=t;

    // std::cerr<<"t: "<<t<<" "<<"z: "<<z<<std::endl;
    // std::cerr<<"missing clique:"<<std::endl;
    // std::cerr<<g->labelsL[1]<<" "<<g->labelsL[184683]<<std::endl;
    // std::cerr<<g->labelsR[133490]<<std::endl;

    runtree.root=runtree.getNewNode();
    runtree.treeNodePool[runtree.root].inValue=-1;
    runtree.treeNodePool[runtree.root].resMark=-1;
    runtree.treeNodePool[runtree.root].qpMark=false;
    runtree.treeNodePool[runtree.root].uv=-1;

    runtree.firstlayerson.resize(g->n[t],-1);
    xStackNode.reserve(g->n[0]+g->n[1]);

    //first layer no pivot
    for(uint32_t u=0;u<g->n[t];u++){
        S[t].c=S[t].r=g->n[t];
        S[z].c=S[z].r=g->n[z];

        

        //form C set for z
        for(uint32_t i=g->pos[t][u];i<g->pos[t][u+1];++i){
            uint32_t v=g->e[t][i];
            S[z].swapByPos(--S[z].c,S[z].pos(v));
        }

        S[z].x=S[z].c;

        //form C set for t
        for(uint32_t i = g->pos[t][u]; i < g->pos[t][u + 1]; i++) {
            uint32_t v = g->e[t][i];
            if(g->pos[z][v + 1] > 0){
                for(uint32_t j = g->pos[z][v + 1] - 1; j >= g->pos[z][v]; j--) {
                    uint32_t w = g->e[z][j];
                    
                    if(w > u) {
                        uint32_t pw = S[t].pos(w);
                        if(pw < S[t].c) {
                            S[t].swapByPos(--S[t].c, pw);
                        }
                    }
                    else break;

                    if(j == 0) break;
                }
            }
        }
        S[t].x = S[t].c;

        //形成t的X集合
        for(uint32_t i = g->pos[t][u]; i < g->pos[t][u + 1]; i++) {
            uint32_t v = g->e[t][i];
            if(g->pos[z][v + 1] > 0){
                for(uint32_t j = g->pos[z][v]; j < g->pos[z][v + 1]; j++) {
                    uint32_t w = g->e[z][j];
                    
                    if(w < u) {
                        uint32_t pw = S[t].pos(w);
                        if(pw < S[t].x) {
                            S[t].swapByPos(--S[t].x, pw);
                        }
                    }
                    else break;

                    if(j == 0) break;
                }
            }
        }


        if(S[z].CIsEmpty()) continue;

        RUV[t].push_back(u);
        RUVlook[t][u]=1;
        // int res=bbranch(0,t,u);
        int res=bbranchLeafRes(0,t,u);
        RUV[t].pop_back();
        RUVlook[t][u]=0;

        if(res!=-1){
            runtree.firstlayerson[u]=res;
            runtree.treeNodePool[runtree.root].sons.push_back(res);
        }

    }

    std::cerr<<"maxBiCliqueCount: "<<cliqueuv[0].size()<<std::endl;
    auto mt2 = std::chrono::steady_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(mt2 - mt1);
    std::cerr << "qbcpc time: " << duration.count() << " ms" << std::endl;
}

void BICPC::baseline(){
    listMBCliques();

    auto mt1 = std::chrono::steady_clock::now();

    std::vector<std::vector<int>> node2bicliques[2];
    node2bicliques[0].resize(g->n[0]+1);
    node2bicliques[1].resize(g->n[1]+1);

    //prepare node to biclique
    for(int i=0;i<cliqueuv[0].size();++i){
        for(int j=0;j<cliqueuv[0][i].size();++j){
            int u=cliqueuv[0][i][j];
            node2bicliques[0][u].push_back(i);
        }
        for(int j=0;j<cliqueuv[1][i].size();++j){
            int v=cliqueuv[1][i][j];
            node2bicliques[1][v].push_back(i);
        }
    }

    std::vector<int> overlap[2];
    std::vector<int> computed;
    int bicliqueNum=cliqueuv[0].size();
    fans.updateSize(bicliqueNum);

    uint64_t cliqueGraphEdgeNum=0;

    overlap[0].resize(bicliqueNum,0);
    overlap[1].resize(bicliqueNum,0);
    computed.reserve(bicliqueNum);
    for(int c=0;c<bicliqueNum;++c){
        for(int j=0;j<cliqueuv[0][c].size();++j){
            int u=cliqueuv[0][c][j];
            for(int k=0;k<node2bicliques[0][u].size();++k){
                int nc=node2bicliques[0][u][k];
                if(fans.find(nc)==fans.find(c)) continue;

                if(overlap[0][nc]++ == 0){
                    computed.push_back(nc);
                }
            }
        }

        for(int j=0;j<cliqueuv[1][c].size();++j){
            int v=cliqueuv[1][c][j];
            for(int k=0;k<node2bicliques[1][v].size();++k){
                int nc=node2bicliques[1][v][k];

                if(overlap[0][nc]<p) continue;
                if(fans.find(nc)==fans.find(c)) continue;

                overlap[1][nc]++;
                if(overlap[1][nc]>=q){
                    fans.merge(nc,c);
                }
            }
        }

        for(int j=0;j<computed.size();++j){
            int nc=computed[j];
            overlap[0][nc]=overlap[1][nc]=0;
        }
        cliqueGraphEdgeNum+=computed.size();
        computed.clear();
    }
    
    std::cerr<<"biclique overlap number: "<<cliqueGraphEdgeNum<<std::endl;
    auto mt2 = std::chrono::steady_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(mt2 - mt1);
    std::cerr << "traverse time: " << duration.count() << " ms" << std::endl;

}

void BICPC::collectBiCliqueLeafResbase(int startNode,std::vector<int>& holdCliqueId){
    // if(runtree.treeNodePool[startNode].resMark!=-1){
    //     holdCliqueId.push_back(runtree.treeNodePool[startNode].resMark);
    //     return;
    // }
    if(runtree.treeNodePool[startNode].sons.size()==0){
        holdCliqueId.push_back(runtree.treeNodePool[startNode].resMark);
        return;
    }

    std::vector<uint32_t>& sons=runtree.treeNodePool[startNode].sons;
    for(int i=0;i<sons.size();++i){
        int inValue=runtree.treeNodePool[sons[i]].inValue;
        int t=runtree.treeNodePool[sons[i]].uv;
        int z=t ^ 1;

        int pvalue=S[t].pos(inValue);
        if((pvalue>=S[t].x&&pvalue<S[t].r)||(RUVlook[t][inValue]==1)){
            collectBiCliqueLeafResbase(sons[i],holdCliqueId);
            // if(RUVlook[t][inValue]==1) break;
        }
    }
}

void BICPC::collectBiCliqueLeafRes(int startNode,std::vector<int>& holdCliqueId){
    if(runtree.treeNodePool[startNode].resMark!=-1){
        holdCliqueId.push_back(runtree.treeNodePool[startNode].resMark);
        return;
    }

    std::vector<uint32_t>& sons=runtree.treeNodePool[startNode].sons;
    for(int i=0;i<sons.size();++i){
        int inValue=runtree.treeNodePool[sons[i]].inValue;
        int t=runtree.treeNodePool[sons[i]].uv;
        int z=t ^ 1;

        int pvalue=S[t].pos(inValue);
        if((pvalue>=S[t].x&&pvalue<S[t].r)||(RUVlook[t][inValue]==1)){
            collectBiCliqueLeafRes(sons[i],holdCliqueId);
            if(RUVlook[t][inValue]==1) break;
        }
    }
}

void BICPC::collectBiClique(int startNode,std::vector<int>& holdCliqueId){
    if(runtree.treeNodePool[startNode].resMark!=-1){
        holdCliqueId.push_back(runtree.treeNodePool[startNode].resMark);
    }

    if(runtree.treeNodePool[startNode].qpMark==true) return;

    std::vector<uint32_t>& sons=runtree.treeNodePool[startNode].sons;
    for(int i=0;i<sons.size();++i){
        int inValue=runtree.treeNodePool[sons[i]].inValue;
        int t=runtree.treeNodePool[sons[i]].uv;
        int z=t ^ 1;

        int pvalue=S[t].pos(inValue);
        if((pvalue>=S[t].x&&pvalue<S[t].r)||(RUVlook[t][inValue]==1)){
            collectBiClique(sons[i],holdCliqueId);
            if(RUVlook[t][inValue]==1) break;
        }
    }
}

void BICPC::connectBiCliqueInXzone(int side){
    // std::cerr<<"here"<<std::endl;
    //debug
    // connectBiCliqueInXzoneCnt++;
    // if(connectBiCliqueInXzoneCnt==682){
    //     connectBiCliqueInXzoneCnt=0;
    // }
    int oriConnectResSize=connectResults.size();
    for(int i=S[firstlayert].x;i<S[firstlayert].c;++i){
        if(runtree.firstlayerson[S[firstlayert][i]]!=-1){
            int startnode=runtree.firstlayerson[S[firstlayert][i]];
            collectBiClique(startnode,connectResults);
            if(RUVlook[firstlayert][runtree.treeNodePool[startnode].inValue]==1){
                std::cerr<<"possible error"<<std::endl;
                exit(0);
            }
        }
    }

    for(int i=0;i<xStackNode.size();++i){
        int inValue=runtree.treeNodePool[xStackNode[i]].inValue;
        int t=runtree.treeNodePool[xStackNode[i]].uv;

        int pvalue=S[t].pos(inValue);
        if((pvalue>=S[t].x&&pvalue<S[t].r)||(RUVlook[t][inValue]==1)){
            collectBiClique(xStackNode[i],connectResults);
            if(RUVlook[t][inValue]==1){
                std::cerr<<"possible error"<<std::endl;
                exit(0);
            }
        }
    }

    return;
}

int BICPC::bbranchLeafRespoorcntpqnode(uint32_t deep,uint32_t side,uint32_t inValue){
    // if(S[1].sz!=g->n[1]){
    //     std::cerr<<"error"<<std::endl;
    // }
    // if(S[z].CIsEmpty()) return;

    int thisNode=-1;//temporally no tree node is built here
    int oriSide=side;

    uint32_t t=side;
    uint32_t z=side ^ 1;

    if(S[t].CIsEmpty()&&S[z].CIsEmpty()){
        if(S[t].XIsEmpty()&&S[z].XIsEmpty()){
            if(RUV[t].size()>=threshold[t]&&RUV[z].size()>=threshold[z]&&RUV[t].size()+RUV[z].size()>threshold[t]+threshold[z]){
                thisNode=runtree.getNewNode();
                runtree.treeNodePool[thisNode].inValue=inValue;
                runtree.treeNodePool[thisNode].uv=side;
                runtree.treeNodePool[thisNode].resMark=-1;
                runtree.treeNodePool[thisNode].qpMark=false;

                cliqueuv[t].emplace_back(RUV[t]);
                cliqueuv[z].emplace_back(RUV[z]);

                runtree.treeNodePool[thisNode].resMark=cliqueuv[t].size()-1;

                //check pqnode
                if(RUV[t].size()==threshold[t]&&t==oriSide||RUV[z].size()==threshold[z]&&z==oriSide){
                    pqnodecnt++;
                }
            }
            

            return thisNode;
        }
        // if(RUV[t].size()>=threshold[t]&&RUV[z].size()>=threshold[z]){
        //     if(RUV[t].size()==threshold[t]&&t==oriSide||RUV[z].size()==threshold[z]&&z==oriSide){
        //         connectAllBiCliqueLeafResbase();
        //     }
        // }
        return -1;
    }

    if(RUV[t].size()+S[t].cSize()<threshold[t]||RUV[z].size()+S[z].cSize()<threshold[z]){
        return thisNode;
    }

    //find pivot

    //debug
    // if(S[1].sz!=g->n[1]){
    //     std::cerr<<"error"<<std::endl;
    // }

    //initialize copies of xcr
    uint32_t x[2];
    uint32_t c[2];
    uint32_t r[2];
    x[t] = S[t].x; x[z] = S[z].x;
    c[t] = S[t].c; c[z] = S[z].c;
    r[t] = S[t].r; r[z] = S[z].r;

    int oriXstackSize=xStackNode.size();

    //find pivot
    uint32_t maxI[2] = {0}, maxE[2] = {0}, maxUV[2] = {0};
    maxI[0] = S[0].c;
    maxI[1] = S[1].c;
    maxUV[0] = S[0][S[0].c];
    maxUV[1] = S[1][S[1].c];
    t=0,z=1;
    uint32_t wsSize;


    //find pivot in C of t
    for(uint32_t i = S[t].c; i < S[t].r; i++) {//scan all vertices in C
    // for(uint32_t i = S[t].c; i <= S[t].c; i++) {
        uint32_t u = S[t][i];
        uint32_t e = 0;

        if(g->deg(u, t) > S[z].r - S[z].c) {
            for(uint32_t j = S[z].c; j < S[z].r; j++) {
                uint32_t v = S[z][j];
                if(g->connect(u, v, t)) {
                    e++;
                    deg[v]++;
                }
            }
        }
        else {
            for(uint32_t j = g->pos[t][u]; j < g->pos[t][u + 1]; j++) {
                uint32_t v = g->e[t][j];
                uint32_t pv = S[z].pos(v);
                if(S[z].c <= pv && pv < S[z].r) {//in C
                    e++;
                    deg[v]++;
                }
            }
        }

        if(e > maxE[z]) {
            maxE[z] = e;
            maxI[t] = i;
            maxUV[t] = u;
        }
    }
    //find pivot in C of z
    for(uint32_t i = S[1].c; i < S[1].r; i++) {//scan all vertices in C
    // for(uint32_t i = S[t].c; i <= S[t].c; i++) {
        uint32_t u = S[1][i];
        uint32_t e = deg[u];
        deg[u] = 0;
        if(e > maxE[0]) {
            maxE[0] = e;
            maxI[1] = i;
            maxUV[1] = u;
        }
    }

    //check if pivot in X of both sides
    bool isPivotInX[2] = {false, false};
    for(t = 0; t < 2; t++) {
        z = t ^ 1;//t + z = 1

        for(uint32_t i = S[t].x; i < S[t].c; i++) {//scan all vertices in X
        // for(uint32_t i = S[t].c; i <= S[t].c; i++) {
            uint32_t u = S[t][i];
            uint32_t e = 0;

            if(g->deg(u, t) > S[z].r - S[z].c) {
                for(uint32_t j = S[z].c; j < S[z].r; j++) {
                    uint32_t v = S[z][j];
                    if(g->connect(u, v, t)) {
                        e++;
                    }
                }
            }
            else {
                for(uint32_t j = g->pos[t][u]; j < g->pos[t][u + 1]; j++) {
                    uint32_t v = g->e[t][j];
                    uint32_t pv = S[z].pos(v);
                    if(S[z].c <= pv && pv < S[z].r) {//in C
                        e++;
                    }
                }
            }

            if(e > maxE[z]) {
                isPivotInX[t] = true;
                maxE[z] = e;
                maxI[t] = i;
                maxUV[t] = u;
            }
        }
    }

    //pivot found. pivot must have at least one neighbor in this algorithm
    t = 1;
    if(S[0].cSize() - maxE[0] >= S[1].cSize() - maxE[1]) {
        t = 0;
    }
    z = t ^ 1;


    //r is more like the right limit of C
    //r[z] doesnt follow S[z].r
    //pivot could be in X, and C of z is changed, but never mind, it will recover by r[z]
    if(g->deg(maxUV[t], t) > S[z].r - S[z].c) {
        if(S[z].r>0){
            for(uint32_t i = S[z].r - 1; i >= S[z].c; i--) {
                if(!g->connect(maxUV[t], S[z][i], t)) {
                    S[z].swapByPos(--S[z].r, i);
                }
    
                if(i == 0) break;
            }
        }
    }else{
        uint32_t u = maxUV[t];
        uint32_t tmpR = S[z].c;
        for(uint32_t j = g->pos[t][u]; j < g->pos[t][u + 1]; j++) {
            uint32_t v = g->e[t][j];
            uint32_t pv = S[z].pos(v);
            if(S[z].c <= pv && pv < S[z].r) {
                S[z].swapByPos(tmpR++, pv);
            }
        }
        S[z].r = tmpR;
    }
    wsSize = r[z] - S[z].r;
    if(ws[deep].size() < wsSize) {
        ws[deep].resize(wsSize * 2);
    }
    memcpy(ws[deep].data(), S[z].begin() + S[z].r, sizeof(uint32_t) * wsSize);


    //iterate on pivot
    if(isPivotInX[t]==false){
        
        S[t].swapByPos(maxI[t], --r[t]);
        S[t].r = r[t];

        //x[z] doesnt follow S[z].x, used to recover S[z].x
        if(g->deg(maxUV[t], t) > S[z].c - S[z].x) {
            for(uint32_t i = S[z].x; i < S[z].c; i++) {
                if(!g->connect(maxUV[t], S[z][i], t)) {
                    S[z].swapByPos(S[z].x++, i);
                }
            }
        }else{
            uint32_t u = maxUV[t];
            uint32_t tmpX = S[z].c;
            for(uint32_t j = g->pos[t][u]; j < g->pos[t][u + 1]; j++) {
                uint32_t v = g->e[t][j];
                uint32_t pv = S[z].pos(v);
                if(S[z].x <= pv && pv < S[z].c) {
                    S[z].swapByPos(--tmpX, pv);
                }
            }
            S[z].x = tmpX;
        }



        RUV[t].push_back(maxUV[t]);
        RUVlook[t][maxUV[t]]=1;

        //debug
        // if(S[1].sz!=g->n[1]){
        //     std::cerr<<"error"<<std::endl;
        // }
        // if(thisNode==401812){
        //     std::cerr<<"here"<<std::endl;
        // }

        #ifdef PARTIALPIVOT
        int res;
        if(S[z].cSize()+S[t].cSize()<=PARTIALPIVOTVAL){
            res=bbranchLeafRespoorcntpqnode_simplep(deep+1,t,maxUV[t]);
        }else{
            res=bbranchLeafRespoorcntpqnode(deep+1,t,maxUV[t]);
        }
        #else
        //C of z must not be empty
        int res=bbranchLeafRespoorcntpqnode(deep+1,t,maxUV[t]);
        #endif

        //debug
        // if(S[1].sz!=g->n[1]){
        //     std::cerr<<"error"<<std::endl;
        // }

        RUV[t].pop_back();
        RUVlook[t][maxUV[t]]=0;
        if(res!=-1){
            if(thisNode==-1){
                thisNode=runtree.getNewNode();
                runtree.treeNodePool[thisNode].inValue=inValue;
                runtree.treeNodePool[thisNode].uv=side;
                runtree.treeNodePool[thisNode].resMark=-1;
                runtree.treeNodePool[thisNode].qpMark=false;
            }
            runtree.treeNodePool[thisNode].sons.push_back(res);
            xStackNode.push_back(res);
        }

        //debug
        // if(S[1].sz!=g->n[1]){
        //     std::cerr<<"error"<<std::endl;
        // }

        //put pivot into X
        S[t].swapByPos(c[t]++, r[t]++);
    }




    // if(c[t]<r[t]){//S[t].c>0
        //iterate on vertices beyound pivot
        for(uint32_t j = 0; j < wsSize; j++){
            uint32_t w = ws[deep][j];
            S[z].swapByPos(S[z].pos(w), --r[z]);

            //recover next level parameter
            S[t].x = x[t]; S[z].x = x[z];
            S[t].c = c[t]; S[z].c = c[z];
            S[t].r = r[t]; S[z].r = r[z];

            if(g->deg(w, z) > S[t].c - S[t].x) {
                for(uint32_t i = S[t].x; i < S[t].c; i++) {
                    if(!g->connect(w, S[t][i], z)) {
                        S[t].swapByPos(S[t].x++, i);
                    }
                }
            }
            else {
                uint32_t tmpX = S[t].c;
                for(uint32_t j = g->pos[z][w]; j < g->pos[z][w + 1]; j++) {
                    uint32_t u = g->e[z][j];
                    uint32_t pu = S[t].pos(u);
                    if(S[t].x <= pu && pu < S[t].c) {
                        S[t].swapByPos(--tmpX, pu);
                    }
                }
                S[t].x = tmpX;
            }

            if(g->deg(w, z) > S[t].r - S[t].c) {
                if(S[t].r > 0){
                    for(uint32_t i = S[t].r - 1; i >= S[t].c; i--) {
                        if(!g->connect(w, S[t][i], z)) {
                            S[t].swapByPos(--S[t].r, i);
                        }

                        if(i == 0) break;
                    }
                }
            }else{
                uint32_t tmpR = S[t].c;
                for(uint32_t j = g->pos[z][w]; j < g->pos[z][w + 1]; j++) {
                    uint32_t u = g->e[z][j];
                    uint32_t pu = S[t].pos(u);
                    if(S[t].c <= pu && pu < S[t].r) {
                        S[t].swapByPos(tmpR++, pu);
                    }
                }
                S[t].r = tmpR;
            }

            // //debug
            // if(deep==0&&inValue==1859991&&w==1911392){
            //     // std::cerr<<<<std::endl;
            //     for(int ii=S[z].c;ii<S[z].r;++ii){
            //         std::cerr<<g->old_lables[z][S[z][ii]]<<" ";
            //     }
            //     std::cerr<<std::endl;
                
            // }

            // if(S[t].cSize()>0){
                RUV[z].push_back(w);
                RUVlook[z][w]=1;
                //debug
                // if(S[1].sz!=g->n[1]){
                //     std::cerr<<"error"<<std::endl;
                // }

                #ifdef PARTIALPIVOT
                int res;
                if(S[z].cSize()+S[t].cSize()<=PARTIALPIVOTVAL){
                    res=bbranchLeafRespoorcntpqnode_simplep(deep+1,z,w);
                }else{
                    res=bbranchLeafRespoorcntpqnode(deep+1,z,w);
                }
                #else
                //C of z must not be empty
                int res=bbranchLeafRespoorcntpqnode(deep+1,z,w);
                #endif

                RUV[z].pop_back();
                RUVlook[z][w]=0;

                if(res!=-1){
                    if(thisNode==-1){
                        thisNode=runtree.getNewNode();
                        runtree.treeNodePool[thisNode].inValue=inValue;
                        runtree.treeNodePool[thisNode].uv=side;
                        runtree.treeNodePool[thisNode].resMark=-1;
                        runtree.treeNodePool[thisNode].qpMark=false;
                    }
                    runtree.treeNodePool[thisNode].sons.push_back(res);
                    xStackNode.push_back(res);
                }
            // }

            S[z].swapByPos(c[z]++,r[z]++);//move w to X
        }

        for(uint32_t j = 0; j < wsSize; j++) {
            uint32_t w = ws[deep][j];
            S[z].swapByPos(S[z].pos(w), --c[z]);
        }
    // }

    //recover pivot
    if(isPivotInX[t] == false) {
        S[t].swapByPos(S[t].pos(maxUV[t]), --c[t]);
    }

    //recover x,c,r
    S[t].x = x[t]; S[z].x = x[z];
    S[t].c = c[t]; S[z].c = c[z];
    S[t].r = r[t]; S[z].r = r[z];

    //connect that should be connected
    //almost the same as before
    if(RUV[t].size()>=threshold[t]&&RUV[z].size()>=threshold[z]){
        int res=-1;
        // if(RUV[t].size()==threshold[t]&&t==oriSide||RUV[z].size()==threshold[z]&&z==oriSide){// only the first one
        //     res=connectAllBiCliqueLeafResbase();
        //     if(thisNode!=-1){
        //         runtree.treeNodePool[thisNode].resMark=res;
        //     }
        // }else{
        //     if(thisNode!=-1){
        //         fans.updateSize(cliqueuv[0].size());
        //         for(int i=0;i<runtree.treeNodePool[thisNode].sons.size();++i){
        //             int sonNode=runtree.treeNodePool[thisNode].sons[i];
        //             connectResults.push_back(runtree.treeNodePool[sonNode].resMark);
        //         }

        //         for(int i=1;i<connectResults.size();++i){
        //             fans.merge(connectResults[i],connectResults[0]);
        //         }
        //         runtree.treeNodePool[thisNode].resMark=fans.find(connectResults[0]);
        //         connectResults.clear();
        //     }

        // }

        if(thisNode!=-1){
            fans.updateSize(cliqueuv[0].size());
            for(int i=0;i<runtree.treeNodePool[thisNode].sons.size();++i){
                int sonNode=runtree.treeNodePool[thisNode].sons[i];
                connectResults.push_back(runtree.treeNodePool[sonNode].resMark);
            }

            for(int i=1;i<connectResults.size();++i){
                fans.merge(connectResults[i],connectResults[0]);
            }
            runtree.treeNodePool[thisNode].resMark=fans.find(connectResults[0]);
            connectResults.clear();

            //check pqnode
            if(RUV[t].size()==threshold[t]&&t==oriSide||RUV[z].size()==threshold[z]&&z==oriSide){
                pqnodecnt++;
            }
        }


    }

    xStackNode.resize(oriXstackSize);
    return thisNode;
}

int BICPC::bbranchLeafRespoor(uint32_t deep,uint32_t side,uint32_t inValue){
    // if(S[1].sz!=g->n[1]){
    //     std::cerr<<"error"<<std::endl;
    // }
    // if(S[z].CIsEmpty()) return;

    int thisNode=-1;//temporally no tree node is built here
    int oriSide=side;

    uint32_t t=side;
    uint32_t z=side ^ 1;

    if(S[t].CIsEmpty()&&S[z].CIsEmpty()){
        if(S[t].XIsEmpty()&&S[z].XIsEmpty()){
            if(RUV[t].size()>=threshold[t]&&RUV[z].size()>=threshold[z]&&RUV[t].size()+RUV[z].size()>threshold[t]+threshold[z]){
                thisNode=runtree.getNewNode();
                runtree.treeNodePool[thisNode].inValue=inValue;
                runtree.treeNodePool[thisNode].uv=side;
                runtree.treeNodePool[thisNode].resMark=-1;
                runtree.treeNodePool[thisNode].qpMark=false;

                cliqueuv[t].emplace_back(RUV[t]);
                cliqueuv[z].emplace_back(RUV[z]);

                runtree.treeNodePool[thisNode].resMark=cliqueuv[t].size()-1;
            }
            return thisNode;
        }
        // if(RUV[t].size()>=threshold[t]&&RUV[z].size()>=threshold[z]){
        //     if(RUV[t].size()==threshold[t]&&t==oriSide||RUV[z].size()==threshold[z]&&z==oriSide){
        //         connectAllBiCliqueLeafResbase();
        //     }
        // }
        return -1;
    }

    if(RUV[t].size()+S[t].cSize()<threshold[t]||RUV[z].size()+S[z].cSize()<threshold[z]){
        return thisNode;
    }

    //find pivot

    //debug
    // if(S[1].sz!=g->n[1]){
    //     std::cerr<<"error"<<std::endl;
    // }

    //initialize copies of xcr
    uint32_t x[2];
    uint32_t c[2];
    uint32_t r[2];
    x[t] = S[t].x; x[z] = S[z].x;
    c[t] = S[t].c; c[z] = S[z].c;
    r[t] = S[t].r; r[z] = S[z].r;

    int oriXstackSize=xStackNode.size();

    //find pivot
    uint32_t maxI[2] = {0}, maxE[2] = {0}, maxUV[2] = {0};
    maxI[0] = S[0].c;
    maxI[1] = S[1].c;
    maxUV[0] = S[0][S[0].c];
    maxUV[1] = S[1][S[1].c];
    t=0,z=1;
    uint32_t wsSize;


    //find pivot in C of t
    for(uint32_t i = S[t].c; i < S[t].r; i++) {//scan all vertices in C
    // for(uint32_t i = S[t].c; i <= S[t].c; i++) {
        uint32_t u = S[t][i];
        uint32_t e = 0;

        if(g->deg(u, t) > S[z].r - S[z].c) {
            for(uint32_t j = S[z].c; j < S[z].r; j++) {
                uint32_t v = S[z][j];
                if(g->connect(u, v, t)) {
                    e++;
                    deg[v]++;
                }
            }
        }
        else {
            for(uint32_t j = g->pos[t][u]; j < g->pos[t][u + 1]; j++) {
                uint32_t v = g->e[t][j];
                uint32_t pv = S[z].pos(v);
                if(S[z].c <= pv && pv < S[z].r) {//in C
                    e++;
                    deg[v]++;
                }
            }
        }

        if(e > maxE[z]) {
            maxE[z] = e;
            maxI[t] = i;
            maxUV[t] = u;
        }
    }
    //find pivot in C of z
    for(uint32_t i = S[1].c; i < S[1].r; i++) {//scan all vertices in C
    // for(uint32_t i = S[t].c; i <= S[t].c; i++) {
        uint32_t u = S[1][i];
        uint32_t e = deg[u];
        deg[u] = 0;
        if(e > maxE[0]) {
            maxE[0] = e;
            maxI[1] = i;
            maxUV[1] = u;
        }
    }

    //check if pivot in X of both sides
    bool isPivotInX[2] = {false, false};
    for(t = 0; t < 2; t++) {
        z = t ^ 1;//t + z = 1

        for(uint32_t i = S[t].x; i < S[t].c; i++) {//scan all vertices in X
        // for(uint32_t i = S[t].c; i <= S[t].c; i++) {
            uint32_t u = S[t][i];
            uint32_t e = 0;

            if(g->deg(u, t) > S[z].r - S[z].c) {
                for(uint32_t j = S[z].c; j < S[z].r; j++) {
                    uint32_t v = S[z][j];
                    if(g->connect(u, v, t)) {
                        e++;
                    }
                }
            }
            else {
                for(uint32_t j = g->pos[t][u]; j < g->pos[t][u + 1]; j++) {
                    uint32_t v = g->e[t][j];
                    uint32_t pv = S[z].pos(v);
                    if(S[z].c <= pv && pv < S[z].r) {//in C
                        e++;
                    }
                }
            }

            if(e > maxE[z]) {
                isPivotInX[t] = true;
                maxE[z] = e;
                maxI[t] = i;
                maxUV[t] = u;
            }
        }
    }

    //pivot found. pivot must have at least one neighbor in this algorithm
    t = 1;
    if(S[0].cSize() - maxE[0] >= S[1].cSize() - maxE[1]) {
        t = 0;
    }
    z = t ^ 1;


    //r is more like the right limit of C
    //r[z] doesnt follow S[z].r
    //pivot could be in X, and C of z is changed, but never mind, it will recover by r[z]
    if(g->deg(maxUV[t], t) > S[z].r - S[z].c) {
        if(S[z].r>0){
            for(uint32_t i = S[z].r - 1; i >= S[z].c; i--) {
                if(!g->connect(maxUV[t], S[z][i], t)) {
                    S[z].swapByPos(--S[z].r, i);
                }
    
                if(i == 0) break;
            }
        }
    }else{
        uint32_t u = maxUV[t];
        uint32_t tmpR = S[z].c;
        for(uint32_t j = g->pos[t][u]; j < g->pos[t][u + 1]; j++) {
            uint32_t v = g->e[t][j];
            uint32_t pv = S[z].pos(v);
            if(S[z].c <= pv && pv < S[z].r) {
                S[z].swapByPos(tmpR++, pv);
            }
        }
        S[z].r = tmpR;
    }
    wsSize = r[z] - S[z].r;
    if(ws[deep].size() < wsSize) {
        ws[deep].resize(wsSize * 2);
    }
    memcpy(ws[deep].data(), S[z].begin() + S[z].r, sizeof(uint32_t) * wsSize);


    //iterate on pivot
    if(isPivotInX[t]==false){
        
        S[t].swapByPos(maxI[t], --r[t]);
        S[t].r = r[t];

        //x[z] doesnt follow S[z].x, used to recover S[z].x
        if(g->deg(maxUV[t], t) > S[z].c - S[z].x) {
            for(uint32_t i = S[z].x; i < S[z].c; i++) {
                if(!g->connect(maxUV[t], S[z][i], t)) {
                    S[z].swapByPos(S[z].x++, i);
                }
            }
        }else{
            uint32_t u = maxUV[t];
            uint32_t tmpX = S[z].c;
            for(uint32_t j = g->pos[t][u]; j < g->pos[t][u + 1]; j++) {
                uint32_t v = g->e[t][j];
                uint32_t pv = S[z].pos(v);
                if(S[z].x <= pv && pv < S[z].c) {
                    S[z].swapByPos(--tmpX, pv);
                }
            }
            S[z].x = tmpX;
        }



        RUV[t].push_back(maxUV[t]);
        RUVlook[t][maxUV[t]]=1;

        //debug
        // if(S[1].sz!=g->n[1]){
        //     std::cerr<<"error"<<std::endl;
        // }
        // if(thisNode==401812){
        //     std::cerr<<"here"<<std::endl;
        // }

        #ifdef PARTIALPIVOT
        int res;
        if(S[z].cSize()+S[t].cSize()<=PARTIALPIVOTVAL){
            res=bbranchLeafRespoor_simplep(deep+1,t,maxUV[t]);
        }else{
            res=bbranchLeafRespoor(deep+1,t,maxUV[t]);
        }
        #else
        //C of z must not be empty
        int res=bbranchLeafRespoor(deep+1,t,maxUV[t]);
        #endif

        //debug
        // if(S[1].sz!=g->n[1]){
        //     std::cerr<<"error"<<std::endl;
        // }

        RUV[t].pop_back();
        RUVlook[t][maxUV[t]]=0;
        if(res!=-1){
            if(thisNode==-1){
                thisNode=runtree.getNewNode();
                runtree.treeNodePool[thisNode].inValue=inValue;
                runtree.treeNodePool[thisNode].uv=side;
                runtree.treeNodePool[thisNode].resMark=-1;
                runtree.treeNodePool[thisNode].qpMark=false;
            }
            runtree.treeNodePool[thisNode].sons.push_back(res);
            xStackNode.push_back(res);
        }

        //debug
        // if(S[1].sz!=g->n[1]){
        //     std::cerr<<"error"<<std::endl;
        // }

        //put pivot into X
        S[t].swapByPos(c[t]++, r[t]++);
    }




    // if(c[t]<r[t]){//S[t].c>0
        //iterate on vertices beyound pivot
        for(uint32_t j = 0; j < wsSize; j++){
            uint32_t w = ws[deep][j];
            S[z].swapByPos(S[z].pos(w), --r[z]);

            //recover next level parameter
            S[t].x = x[t]; S[z].x = x[z];
            S[t].c = c[t]; S[z].c = c[z];
            S[t].r = r[t]; S[z].r = r[z];

            if(g->deg(w, z) > S[t].c - S[t].x) {
                for(uint32_t i = S[t].x; i < S[t].c; i++) {
                    if(!g->connect(w, S[t][i], z)) {
                        S[t].swapByPos(S[t].x++, i);
                    }
                }
            }
            else {
                uint32_t tmpX = S[t].c;
                for(uint32_t j = g->pos[z][w]; j < g->pos[z][w + 1]; j++) {
                    uint32_t u = g->e[z][j];
                    uint32_t pu = S[t].pos(u);
                    if(S[t].x <= pu && pu < S[t].c) {
                        S[t].swapByPos(--tmpX, pu);
                    }
                }
                S[t].x = tmpX;
            }

            if(g->deg(w, z) > S[t].r - S[t].c) {
                if(S[t].r > 0){
                    for(uint32_t i = S[t].r - 1; i >= S[t].c; i--) {
                        if(!g->connect(w, S[t][i], z)) {
                            S[t].swapByPos(--S[t].r, i);
                        }

                        if(i == 0) break;
                    }
                }
            }else{
                uint32_t tmpR = S[t].c;
                for(uint32_t j = g->pos[z][w]; j < g->pos[z][w + 1]; j++) {
                    uint32_t u = g->e[z][j];
                    uint32_t pu = S[t].pos(u);
                    if(S[t].c <= pu && pu < S[t].r) {
                        S[t].swapByPos(tmpR++, pu);
                    }
                }
                S[t].r = tmpR;
            }

            // //debug
            // if(deep==0&&inValue==1859991&&w==1911392){
            //     // std::cerr<<<<std::endl;
            //     for(int ii=S[z].c;ii<S[z].r;++ii){
            //         std::cerr<<g->old_lables[z][S[z][ii]]<<" ";
            //     }
            //     std::cerr<<std::endl;
                
            // }

            // if(S[t].cSize()>0){
                RUV[z].push_back(w);
                RUVlook[z][w]=1;
                //debug
                // if(S[1].sz!=g->n[1]){
                //     std::cerr<<"error"<<std::endl;
                // }

                #ifdef PARTIALPIVOT
                int res;
                if(S[z].cSize()+S[t].cSize()<=PARTIALPIVOTVAL){
                    res=bbranchLeafRespoor_simplep(deep+1,z,w);
                }else{
                    res=bbranchLeafRespoor(deep+1,z,w);
                }
                #else
                //C of z must not be empty
                int res=bbranchLeafRespoor(deep+1,z,w);
                #endif

                RUV[z].pop_back();
                RUVlook[z][w]=0;

                if(res!=-1){
                    if(thisNode==-1){
                        thisNode=runtree.getNewNode();
                        runtree.treeNodePool[thisNode].inValue=inValue;
                        runtree.treeNodePool[thisNode].uv=side;
                        runtree.treeNodePool[thisNode].resMark=-1;
                        runtree.treeNodePool[thisNode].qpMark=false;
                    }
                    runtree.treeNodePool[thisNode].sons.push_back(res);
                    xStackNode.push_back(res);
                }
            // }

            S[z].swapByPos(c[z]++,r[z]++);//move w to X
        }

        for(uint32_t j = 0; j < wsSize; j++) {
            uint32_t w = ws[deep][j];
            S[z].swapByPos(S[z].pos(w), --c[z]);
        }
    // }

    //recover pivot
    if(isPivotInX[t] == false) {
        S[t].swapByPos(S[t].pos(maxUV[t]), --c[t]);
    }

    //recover x,c,r
    S[t].x = x[t]; S[z].x = x[z];
    S[t].c = c[t]; S[z].c = c[z];
    S[t].r = r[t]; S[z].r = r[z];

    //connect that should be connected
    //almost the same as before
    if(RUV[t].size()>=threshold[t]&&RUV[z].size()>=threshold[z]){
        int res=-1;
        // if(RUV[t].size()==threshold[t]&&t==oriSide||RUV[z].size()==threshold[z]&&z==oriSide){// only the first one
        //     res=connectAllBiCliqueLeafResbase();
        //     if(thisNode!=-1){
        //         runtree.treeNodePool[thisNode].resMark=res;
        //     }
        // }else{
        //     if(thisNode!=-1){
        //         fans.updateSize(cliqueuv[0].size());
        //         for(int i=0;i<runtree.treeNodePool[thisNode].sons.size();++i){
        //             int sonNode=runtree.treeNodePool[thisNode].sons[i];
        //             connectResults.push_back(runtree.treeNodePool[sonNode].resMark);
        //         }

        //         for(int i=1;i<connectResults.size();++i){
        //             fans.merge(connectResults[i],connectResults[0]);
        //         }
        //         runtree.treeNodePool[thisNode].resMark=fans.find(connectResults[0]);
        //         connectResults.clear();
        //     }

        // }

        if(thisNode!=-1){
            fans.updateSize(cliqueuv[0].size());
            for(int i=0;i<runtree.treeNodePool[thisNode].sons.size();++i){
                int sonNode=runtree.treeNodePool[thisNode].sons[i];
                connectResults.push_back(runtree.treeNodePool[sonNode].resMark);
            }

            for(int i=1;i<connectResults.size();++i){
                fans.merge(connectResults[i],connectResults[0]);
            }
            runtree.treeNodePool[thisNode].resMark=fans.find(connectResults[0]);
            connectResults.clear();
        }


    }

    xStackNode.resize(oriXstackSize);
    return thisNode;
}

int BICPC::bbranchLeafResbase(uint32_t deep,uint32_t side,uint32_t inValue){
    // if(S[1].sz!=g->n[1]){
    //     std::cerr<<"error"<<std::endl;
    // }
    // if(S[z].CIsEmpty()) return;

    int thisNode=-1;//temporally no tree node is built here
    int oriSide=side;

    uint32_t t=side;
    uint32_t z=side ^ 1;

    if(S[t].CIsEmpty()&&S[z].CIsEmpty()){
        if(S[t].XIsEmpty()&&S[z].XIsEmpty()){
            if(RUV[t].size()>=threshold[t]&&RUV[z].size()>=threshold[z]&&RUV[t].size()+RUV[z].size()>threshold[t]+threshold[z]){
                thisNode=runtree.getNewNode();
                runtree.treeNodePool[thisNode].inValue=inValue;
                runtree.treeNodePool[thisNode].uv=side;
                runtree.treeNodePool[thisNode].resMark=-1;
                runtree.treeNodePool[thisNode].qpMark=false;

                cliqueuv[t].emplace_back(RUV[t]);
                cliqueuv[z].emplace_back(RUV[z]);

                runtree.treeNodePool[thisNode].resMark=cliqueuv[t].size()-1;
            }
            return thisNode;
        }
        if(RUV[t].size()>=threshold[t]&&RUV[z].size()>=threshold[z]){
            if(RUV[t].size()==threshold[t]&&t==oriSide||RUV[z].size()==threshold[z]&&z==oriSide){
                connectAllBiCliqueLeafResbase();
            }
        }
        return -1;
    }

    if(RUV[t].size()+S[t].cSize()<threshold[t]||RUV[z].size()+S[z].cSize()<threshold[z]){
        return thisNode;
    }

    //find pivot

    //debug
    // if(S[1].sz!=g->n[1]){
    //     std::cerr<<"error"<<std::endl;
    // }

    //initialize copies of xcr
    uint32_t x[2];
    uint32_t c[2];
    uint32_t r[2];
    x[t] = S[t].x; x[z] = S[z].x;
    c[t] = S[t].c; c[z] = S[z].c;
    r[t] = S[t].r; r[z] = S[z].r;

    int oriXstackSize=xStackNode.size();

    //find pivot
    uint32_t maxI[2] = {0}, maxE[2] = {0}, maxUV[2] = {0};
    maxI[0] = S[0].c;
    maxI[1] = S[1].c;
    maxUV[0] = S[0][S[0].c];
    maxUV[1] = S[1][S[1].c];
    t=0,z=1;
    uint32_t wsSize;


    //find pivot in C of t
    for(uint32_t i = S[t].c; i < S[t].r; i++) {//scan all vertices in C
    // for(uint32_t i = S[t].c; i <= S[t].c; i++) {
        uint32_t u = S[t][i];
        uint32_t e = 0;

        if(g->deg(u, t) > S[z].r - S[z].c) {
            for(uint32_t j = S[z].c; j < S[z].r; j++) {
                uint32_t v = S[z][j];
                if(g->connect(u, v, t)) {
                    e++;
                    deg[v]++;
                }
            }
        }
        else {
            for(uint32_t j = g->pos[t][u]; j < g->pos[t][u + 1]; j++) {
                uint32_t v = g->e[t][j];
                uint32_t pv = S[z].pos(v);
                if(S[z].c <= pv && pv < S[z].r) {//in C
                    e++;
                    deg[v]++;
                }
            }
        }

        if(e > maxE[z]) {
            maxE[z] = e;
            maxI[t] = i;
            maxUV[t] = u;
        }
    }
    //find pivot in C of z
    for(uint32_t i = S[1].c; i < S[1].r; i++) {//scan all vertices in C
    // for(uint32_t i = S[t].c; i <= S[t].c; i++) {
        uint32_t u = S[1][i];
        uint32_t e = deg[u];
        deg[u] = 0;
        if(e > maxE[0]) {
            maxE[0] = e;
            maxI[1] = i;
            maxUV[1] = u;
        }
    }

    //check if pivot in X of both sides
    bool isPivotInX[2] = {false, false};
    for(t = 0; t < 2; t++) {
        z = t ^ 1;//t + z = 1

        for(uint32_t i = S[t].x; i < S[t].c; i++) {//scan all vertices in X
        // for(uint32_t i = S[t].c; i <= S[t].c; i++) {
            uint32_t u = S[t][i];
            uint32_t e = 0;

            if(g->deg(u, t) > S[z].r - S[z].c) {
                for(uint32_t j = S[z].c; j < S[z].r; j++) {
                    uint32_t v = S[z][j];
                    if(g->connect(u, v, t)) {
                        e++;
                    }
                }
            }
            else {
                for(uint32_t j = g->pos[t][u]; j < g->pos[t][u + 1]; j++) {
                    uint32_t v = g->e[t][j];
                    uint32_t pv = S[z].pos(v);
                    if(S[z].c <= pv && pv < S[z].r) {//in C
                        e++;
                    }
                }
            }

            if(e > maxE[z]) {
                isPivotInX[t] = true;
                maxE[z] = e;
                maxI[t] = i;
                maxUV[t] = u;
            }
        }
    }

    //pivot found. pivot must have at least one neighbor in this algorithm
    t = 1;
    if(S[0].cSize() - maxE[0] >= S[1].cSize() - maxE[1]) {
        t = 0;
    }
    z = t ^ 1;


    //r is more like the right limit of C
    //r[z] doesnt follow S[z].r
    //pivot could be in X, and C of z is changed, but never mind, it will recover by r[z]
    if(g->deg(maxUV[t], t) > S[z].r - S[z].c) {
        if(S[z].r>0){
            for(uint32_t i = S[z].r - 1; i >= S[z].c; i--) {
                if(!g->connect(maxUV[t], S[z][i], t)) {
                    S[z].swapByPos(--S[z].r, i);
                }
    
                if(i == 0) break;
            }
        }
    }else{
        uint32_t u = maxUV[t];
        uint32_t tmpR = S[z].c;
        for(uint32_t j = g->pos[t][u]; j < g->pos[t][u + 1]; j++) {
            uint32_t v = g->e[t][j];
            uint32_t pv = S[z].pos(v);
            if(S[z].c <= pv && pv < S[z].r) {
                S[z].swapByPos(tmpR++, pv);
            }
        }
        S[z].r = tmpR;
    }
    wsSize = r[z] - S[z].r;
    if(ws[deep].size() < wsSize) {
        ws[deep].resize(wsSize * 2);
    }
    memcpy(ws[deep].data(), S[z].begin() + S[z].r, sizeof(uint32_t) * wsSize);


    //iterate on pivot
    if(isPivotInX[t]==false){
        
        S[t].swapByPos(maxI[t], --r[t]);
        S[t].r = r[t];

        //x[z] doesnt follow S[z].x, used to recover S[z].x
        if(g->deg(maxUV[t], t) > S[z].c - S[z].x) {
            for(uint32_t i = S[z].x; i < S[z].c; i++) {
                if(!g->connect(maxUV[t], S[z][i], t)) {
                    S[z].swapByPos(S[z].x++, i);
                }
            }
        }else{
            uint32_t u = maxUV[t];
            uint32_t tmpX = S[z].c;
            for(uint32_t j = g->pos[t][u]; j < g->pos[t][u + 1]; j++) {
                uint32_t v = g->e[t][j];
                uint32_t pv = S[z].pos(v);
                if(S[z].x <= pv && pv < S[z].c) {
                    S[z].swapByPos(--tmpX, pv);
                }
            }
            S[z].x = tmpX;
        }



        RUV[t].push_back(maxUV[t]);
        RUVlook[t][maxUV[t]]=1;

        //debug
        // if(S[1].sz!=g->n[1]){
        //     std::cerr<<"error"<<std::endl;
        // }
        // if(thisNode==401812){
        //     std::cerr<<"here"<<std::endl;
        // }

        #ifdef PARTIALPIVOT
        int res;
        if(S[z].cSize()+S[t].cSize()<=PARTIALPIVOTVAL){
            res=bbranchLeafResbase_simplep(deep+1,t,maxUV[t]);
        }else{
            res=bbranchLeafResbase(deep+1,t,maxUV[t]);
        }
        #else
        //C of z must not be empty
        int res=bbranchLeafResbase(deep+1,t,maxUV[t]);
        #endif

        //debug
        // if(S[1].sz!=g->n[1]){
        //     std::cerr<<"error"<<std::endl;
        // }

        RUV[t].pop_back();
        RUVlook[t][maxUV[t]]=0;
        if(res!=-1){
            if(thisNode==-1){
                thisNode=runtree.getNewNode();
                runtree.treeNodePool[thisNode].inValue=inValue;
                runtree.treeNodePool[thisNode].uv=side;
                runtree.treeNodePool[thisNode].resMark=-1;
                runtree.treeNodePool[thisNode].qpMark=false;
            }
            runtree.treeNodePool[thisNode].sons.push_back(res);
            xStackNode.push_back(res);
        }

        //debug
        // if(S[1].sz!=g->n[1]){
        //     std::cerr<<"error"<<std::endl;
        // }

        //put pivot into X
        S[t].swapByPos(c[t]++, r[t]++);
    }




    // if(c[t]<r[t]){//S[t].c>0
        //iterate on vertices beyound pivot
        for(uint32_t j = 0; j < wsSize; j++){
            uint32_t w = ws[deep][j];
            S[z].swapByPos(S[z].pos(w), --r[z]);

            //recover next level parameter
            S[t].x = x[t]; S[z].x = x[z];
            S[t].c = c[t]; S[z].c = c[z];
            S[t].r = r[t]; S[z].r = r[z];

            if(g->deg(w, z) > S[t].c - S[t].x) {
                for(uint32_t i = S[t].x; i < S[t].c; i++) {
                    if(!g->connect(w, S[t][i], z)) {
                        S[t].swapByPos(S[t].x++, i);
                    }
                }
            }
            else {
                uint32_t tmpX = S[t].c;
                for(uint32_t j = g->pos[z][w]; j < g->pos[z][w + 1]; j++) {
                    uint32_t u = g->e[z][j];
                    uint32_t pu = S[t].pos(u);
                    if(S[t].x <= pu && pu < S[t].c) {
                        S[t].swapByPos(--tmpX, pu);
                    }
                }
                S[t].x = tmpX;
            }

            if(g->deg(w, z) > S[t].r - S[t].c) {
                if(S[t].r > 0){
                    for(uint32_t i = S[t].r - 1; i >= S[t].c; i--) {
                        if(!g->connect(w, S[t][i], z)) {
                            S[t].swapByPos(--S[t].r, i);
                        }

                        if(i == 0) break;
                    }
                }
            }else{
                uint32_t tmpR = S[t].c;
                for(uint32_t j = g->pos[z][w]; j < g->pos[z][w + 1]; j++) {
                    uint32_t u = g->e[z][j];
                    uint32_t pu = S[t].pos(u);
                    if(S[t].c <= pu && pu < S[t].r) {
                        S[t].swapByPos(tmpR++, pu);
                    }
                }
                S[t].r = tmpR;
            }

            // //debug
            // if(deep==0&&inValue==1859991&&w==1911392){
            //     // std::cerr<<<<std::endl;
            //     for(int ii=S[z].c;ii<S[z].r;++ii){
            //         std::cerr<<g->old_lables[z][S[z][ii]]<<" ";
            //     }
            //     std::cerr<<std::endl;
                
            // }

            // if(S[t].cSize()>0){
                RUV[z].push_back(w);
                RUVlook[z][w]=1;
                //debug
                // if(S[1].sz!=g->n[1]){
                //     std::cerr<<"error"<<std::endl;
                // }

                #ifdef PARTIALPIVOT
                int res;
                if(S[z].cSize()+S[t].cSize()<=PARTIALPIVOTVAL){
                    res=bbranchLeafResbase_simplep(deep+1,z,w);
                }else{
                    res=bbranchLeafResbase(deep+1,z,w);
                }
                #else
                //C of z must not be empty
                int res=bbranchLeafResbase(deep+1,z,w);
                #endif

                RUV[z].pop_back();
                RUVlook[z][w]=0;

                if(res!=-1){
                    if(thisNode==-1){
                        thisNode=runtree.getNewNode();
                        runtree.treeNodePool[thisNode].inValue=inValue;
                        runtree.treeNodePool[thisNode].uv=side;
                        runtree.treeNodePool[thisNode].resMark=-1;
                        runtree.treeNodePool[thisNode].qpMark=false;
                    }
                    runtree.treeNodePool[thisNode].sons.push_back(res);
                    xStackNode.push_back(res);
                }
            // }

            S[z].swapByPos(c[z]++,r[z]++);//move w to X
        }

        for(uint32_t j = 0; j < wsSize; j++) {
            uint32_t w = ws[deep][j];
            S[z].swapByPos(S[z].pos(w), --c[z]);
        }
    // }

    //recover pivot
    if(isPivotInX[t] == false) {
        S[t].swapByPos(S[t].pos(maxUV[t]), --c[t]);
    }

    //recover x,c,r
    S[t].x = x[t]; S[z].x = x[z];
    S[t].c = c[t]; S[z].c = c[z];
    S[t].r = r[t]; S[z].r = r[z];

    //connect that should be connected
    //almost the same as before
    if(RUV[t].size()>=threshold[t]&&RUV[z].size()>=threshold[z]){
        int res=-1;
        if(RUV[t].size()==threshold[t]&&t==oriSide||RUV[z].size()==threshold[z]&&z==oriSide){// only the first one
            res=connectAllBiCliqueLeafResbase();
            if(thisNode!=-1){
                runtree.treeNodePool[thisNode].resMark=res;
            }
        }else{
            if(thisNode!=-1){
                fans.updateSize(cliqueuv[0].size());
                for(int i=0;i<runtree.treeNodePool[thisNode].sons.size();++i){
                    int sonNode=runtree.treeNodePool[thisNode].sons[i];
                    connectResults.push_back(runtree.treeNodePool[sonNode].resMark);
                }

                for(int i=1;i<connectResults.size();++i){
                    fans.merge(connectResults[i],connectResults[0]);
                }
                runtree.treeNodePool[thisNode].resMark=fans.find(connectResults[0]);
                connectResults.clear();
            }

        }


    }

    xStackNode.resize(oriXstackSize);
    return thisNode;
}

int BICPC::bbranchLeafResforabnode(uint32_t deep,uint32_t side,uint32_t inValue){
    // if(S[1].sz!=g->n[1]){
    //     std::cerr<<"error"<<std::endl;
    // }
    // if(S[z].CIsEmpty()) return;

    int thisNode=-1;//temporally no tree node is built here
    int oriSide=side;

    uint32_t t=side;
    uint32_t z=side ^ 1;

    if(S[t].CIsEmpty()&&S[z].CIsEmpty()){
        if(S[t].XIsEmpty()&&S[z].XIsEmpty()){
            if(RUV[t].size()>=threshold[t]&&RUV[z].size()>=threshold[z]&&RUV[t].size()+RUV[z].size()>threshold[t]+threshold[z]){
                thisNode=runtree.getNewNode();
                runtree.treeNodePool[thisNode].inValue=inValue;
                runtree.treeNodePool[thisNode].uv=side;
                runtree.treeNodePool[thisNode].resMark=-1;
                runtree.treeNodePool[thisNode].qpMark=false;

                cliqueuv[t].emplace_back(RUV[t]);
                cliqueuv[z].emplace_back(RUV[z]);

                runtree.treeNodePool[thisNode].resMark=cliqueuv[t].size()-1;

                //check pqnode
                if(RUV[t].size()==threshold[t]&&t==oriSide||RUV[z].size()==threshold[z]&&z==oriSide){
                    pqnodecnt++;
                    pqnodedepth_total+=RUV[t].size()+RUV[z].size();
                }
            }
            

            return thisNode;
        }
        //check pqnode
        if(RUV[t].size()>=threshold[t]&&RUV[z].size()>=threshold[z]){
            if(RUV[t].size()==threshold[t]&&t==oriSide||RUV[z].size()==threshold[z]&&z==oriSide){
                // connectAllBiCliqueLeafResbase();
                vpqnodecnt++;
            }
        }
        return -1;
    }

    if(RUV[t].size()+S[t].cSize()<threshold[t]||RUV[z].size()+S[z].cSize()<threshold[z]){
        return thisNode;
    }

    //find pivot

    //debug
    // if(S[1].sz!=g->n[1]){
    //     std::cerr<<"error"<<std::endl;
    // }

    //initialize copies of xcr
    uint32_t x[2];
    uint32_t c[2];
    uint32_t r[2];
    x[t] = S[t].x; x[z] = S[z].x;
    c[t] = S[t].c; c[z] = S[z].c;
    r[t] = S[t].r; r[z] = S[z].r;

    int oriXstackSize=xStackNode.size();

    //find pivot
    uint32_t maxI[2] = {0}, maxE[2] = {0}, maxUV[2] = {0};
    maxI[0] = S[0].c;
    maxI[1] = S[1].c;
    maxUV[0] = S[0][S[0].c];
    maxUV[1] = S[1][S[1].c];
    t=0,z=1;
    uint32_t wsSize;


    //find pivot in C of t
    for(uint32_t i = S[t].c; i < S[t].r; i++) {//scan all vertices in C
    // for(uint32_t i = S[t].c; i <= S[t].c; i++) {
        uint32_t u = S[t][i];
        uint32_t e = 0;

        if(g->deg(u, t) > S[z].r - S[z].c) {
            for(uint32_t j = S[z].c; j < S[z].r; j++) {
                uint32_t v = S[z][j];
                if(g->connect(u, v, t)) {
                    e++;
                    deg[v]++;
                }
            }
        }
        else {
            for(uint32_t j = g->pos[t][u]; j < g->pos[t][u + 1]; j++) {
                uint32_t v = g->e[t][j];
                uint32_t pv = S[z].pos(v);
                if(S[z].c <= pv && pv < S[z].r) {//in C
                    e++;
                    deg[v]++;
                }
            }
        }

        if(e > maxE[z]) {
            maxE[z] = e;
            maxI[t] = i;
            maxUV[t] = u;
        }
    }
    //find pivot in C of z
    for(uint32_t i = S[1].c; i < S[1].r; i++) {//scan all vertices in C
    // for(uint32_t i = S[t].c; i <= S[t].c; i++) {
        uint32_t u = S[1][i];
        uint32_t e = deg[u];
        deg[u] = 0;
        if(e > maxE[0]) {
            maxE[0] = e;
            maxI[1] = i;
            maxUV[1] = u;
        }
    }

    //check if pivot in X of both sides
    bool isPivotInX[2] = {false, false};
    for(t = 0; t < 2; t++) {
        z = t ^ 1;//t + z = 1

        for(uint32_t i = S[t].x; i < S[t].c; i++) {//scan all vertices in X
        // for(uint32_t i = S[t].c; i <= S[t].c; i++) {
            uint32_t u = S[t][i];
            uint32_t e = 0;

            if(g->deg(u, t) > S[z].r - S[z].c) {
                for(uint32_t j = S[z].c; j < S[z].r; j++) {
                    uint32_t v = S[z][j];
                    if(g->connect(u, v, t)) {
                        e++;
                    }
                }
            }
            else {
                for(uint32_t j = g->pos[t][u]; j < g->pos[t][u + 1]; j++) {
                    uint32_t v = g->e[t][j];
                    uint32_t pv = S[z].pos(v);
                    if(S[z].c <= pv && pv < S[z].r) {//in C
                        e++;
                    }
                }
            }

            if(e > maxE[z]) {
                isPivotInX[t] = true;
                maxE[z] = e;
                maxI[t] = i;
                maxUV[t] = u;
            }
        }
    }

    //pivot found. pivot must have at least one neighbor in this algorithm
    t = 1;
    if(S[0].cSize() - maxE[0] >= S[1].cSize() - maxE[1]) {
        t = 0;
    }
    z = t ^ 1;


    //r is more like the right limit of C
    //r[z] doesnt follow S[z].r
    //pivot could be in X, and C of z is changed, but never mind, it will recover by r[z]
    if(g->deg(maxUV[t], t) > S[z].r - S[z].c) {
        if(S[z].r>0){
            for(uint32_t i = S[z].r - 1; i >= S[z].c; i--) {
                if(!g->connect(maxUV[t], S[z][i], t)) {
                    S[z].swapByPos(--S[z].r, i);
                }
    
                if(i == 0) break;
            }
        }
    }else{
        uint32_t u = maxUV[t];
        uint32_t tmpR = S[z].c;
        for(uint32_t j = g->pos[t][u]; j < g->pos[t][u + 1]; j++) {
            uint32_t v = g->e[t][j];
            uint32_t pv = S[z].pos(v);
            if(S[z].c <= pv && pv < S[z].r) {
                S[z].swapByPos(tmpR++, pv);
            }
        }
        S[z].r = tmpR;
    }
    wsSize = r[z] - S[z].r;
    if(ws[deep].size() < wsSize) {
        ws[deep].resize(wsSize * 2);
    }
    memcpy(ws[deep].data(), S[z].begin() + S[z].r, sizeof(uint32_t) * wsSize);


    //iterate on pivot
    if(isPivotInX[t]==false){
        
        S[t].swapByPos(maxI[t], --r[t]);
        S[t].r = r[t];

        //x[z] doesnt follow S[z].x, used to recover S[z].x
        if(g->deg(maxUV[t], t) > S[z].c - S[z].x) {
            for(uint32_t i = S[z].x; i < S[z].c; i++) {
                if(!g->connect(maxUV[t], S[z][i], t)) {
                    S[z].swapByPos(S[z].x++, i);
                }
            }
        }else{
            uint32_t u = maxUV[t];
            uint32_t tmpX = S[z].c;
            for(uint32_t j = g->pos[t][u]; j < g->pos[t][u + 1]; j++) {
                uint32_t v = g->e[t][j];
                uint32_t pv = S[z].pos(v);
                if(S[z].x <= pv && pv < S[z].c) {
                    S[z].swapByPos(--tmpX, pv);
                }
            }
            S[z].x = tmpX;
        }



        RUV[t].push_back(maxUV[t]);
        RUVlook[t][maxUV[t]]=1;

        //debug
        // if(S[1].sz!=g->n[1]){
        //     std::cerr<<"error"<<std::endl;
        // }
        // if(thisNode==401812){
        //     std::cerr<<"here"<<std::endl;
        // }

        #ifdef PARTIALPIVOT
        int res;
        if(S[z].cSize()+S[t].cSize()<=PARTIALPIVOTVAL){
            res=bbranchLeafResforabnode_simplep(deep+1,t,maxUV[t]);
        }else{
            res=bbranchLeafResforabnode(deep+1,t,maxUV[t]);
        }
        #else
        //C of z must not be empty
        int res=bbranchLeafResforabnode(deep+1,t,maxUV[t]);
        #endif

        //debug
        // if(S[1].sz!=g->n[1]){
        //     std::cerr<<"error"<<std::endl;
        // }

        RUV[t].pop_back();
        RUVlook[t][maxUV[t]]=0;
        if(res!=-1){
            if(thisNode==-1){
                thisNode=runtree.getNewNode();
                runtree.treeNodePool[thisNode].inValue=inValue;
                runtree.treeNodePool[thisNode].uv=side;
                runtree.treeNodePool[thisNode].resMark=-1;
                runtree.treeNodePool[thisNode].qpMark=false;
            }
            runtree.treeNodePool[thisNode].sons.push_back(res);
            xStackNode.push_back(res);
        }

        //debug
        // if(S[1].sz!=g->n[1]){
        //     std::cerr<<"error"<<std::endl;
        // }

        //put pivot into X
        S[t].swapByPos(c[t]++, r[t]++);
    }




    // if(c[t]<r[t]){//S[t].c>0
        //iterate on vertices beyound pivot
        for(uint32_t j = 0; j < wsSize; j++){
            uint32_t w = ws[deep][j];
            S[z].swapByPos(S[z].pos(w), --r[z]);

            //recover next level parameter
            S[t].x = x[t]; S[z].x = x[z];
            S[t].c = c[t]; S[z].c = c[z];
            S[t].r = r[t]; S[z].r = r[z];

            if(g->deg(w, z) > S[t].c - S[t].x) {
                for(uint32_t i = S[t].x; i < S[t].c; i++) {
                    if(!g->connect(w, S[t][i], z)) {
                        S[t].swapByPos(S[t].x++, i);
                    }
                }
            }
            else {
                uint32_t tmpX = S[t].c;
                for(uint32_t j = g->pos[z][w]; j < g->pos[z][w + 1]; j++) {
                    uint32_t u = g->e[z][j];
                    uint32_t pu = S[t].pos(u);
                    if(S[t].x <= pu && pu < S[t].c) {
                        S[t].swapByPos(--tmpX, pu);
                    }
                }
                S[t].x = tmpX;
            }

            if(g->deg(w, z) > S[t].r - S[t].c) {
                if(S[t].r > 0){
                    for(uint32_t i = S[t].r - 1; i >= S[t].c; i--) {
                        if(!g->connect(w, S[t][i], z)) {
                            S[t].swapByPos(--S[t].r, i);
                        }

                        if(i == 0) break;
                    }
                }
            }else{
                uint32_t tmpR = S[t].c;
                for(uint32_t j = g->pos[z][w]; j < g->pos[z][w + 1]; j++) {
                    uint32_t u = g->e[z][j];
                    uint32_t pu = S[t].pos(u);
                    if(S[t].c <= pu && pu < S[t].r) {
                        S[t].swapByPos(tmpR++, pu);
                    }
                }
                S[t].r = tmpR;
            }

            // //debug
            // if(deep==0&&inValue==1859991&&w==1911392){
            //     // std::cerr<<<<std::endl;
            //     for(int ii=S[z].c;ii<S[z].r;++ii){
            //         std::cerr<<g->old_lables[z][S[z][ii]]<<" ";
            //     }
            //     std::cerr<<std::endl;
                
            // }

            // if(S[t].cSize()>0){
                RUV[z].push_back(w);
                RUVlook[z][w]=1;
                //debug
                // if(S[1].sz!=g->n[1]){
                //     std::cerr<<"error"<<std::endl;
                // }

                #ifdef PARTIALPIVOT
                int res;
                if(S[z].cSize()+S[t].cSize()<=PARTIALPIVOTVAL){
                    res=bbranchLeafResforabnode_simplep(deep+1,z,w);
                }else{
                    res=bbranchLeafResforabnode(deep+1,z,w);
                }
                #else
                //C of z must not be empty
                int res=bbranchLeafResforabnode(deep+1,z,w);
                #endif

                RUV[z].pop_back();
                RUVlook[z][w]=0;

                if(res!=-1){
                    if(thisNode==-1){
                        thisNode=runtree.getNewNode();
                        runtree.treeNodePool[thisNode].inValue=inValue;
                        runtree.treeNodePool[thisNode].uv=side;
                        runtree.treeNodePool[thisNode].resMark=-1;
                        runtree.treeNodePool[thisNode].qpMark=false;
                    }
                    runtree.treeNodePool[thisNode].sons.push_back(res);
                    xStackNode.push_back(res);
                }
            // }

            S[z].swapByPos(c[z]++,r[z]++);//move w to X
        }

        for(uint32_t j = 0; j < wsSize; j++) {
            uint32_t w = ws[deep][j];
            S[z].swapByPos(S[z].pos(w), --c[z]);
        }
    // }

    //recover pivot
    if(isPivotInX[t] == false) {
        S[t].swapByPos(S[t].pos(maxUV[t]), --c[t]);
    }

    //recover x,c,r
    S[t].x = x[t]; S[z].x = x[z];
    S[t].c = c[t]; S[z].c = c[z];
    S[t].r = r[t]; S[z].r = r[z];

    //connect that should be connected
    //almost the same as before
    if(RUV[t].size()>=threshold[t]&&RUV[z].size()>=threshold[z]){
        int res=-1;
        // if(RUV[t].size()==threshold[t]&&t==oriSide||RUV[z].size()==threshold[z]&&z==oriSide){// only the first one
        //     res=connectAllBiCliqueLeafResbase();
        //     if(thisNode!=-1){
        //         runtree.treeNodePool[thisNode].resMark=res;
        //     }
        // }else{
        //     if(thisNode!=-1){
        //         fans.updateSize(cliqueuv[0].size());
        //         for(int i=0;i<runtree.treeNodePool[thisNode].sons.size();++i){
        //             int sonNode=runtree.treeNodePool[thisNode].sons[i];
        //             connectResults.push_back(runtree.treeNodePool[sonNode].resMark);
        //         }

        //         for(int i=1;i<connectResults.size();++i){
        //             fans.merge(connectResults[i],connectResults[0]);
        //         }
        //         runtree.treeNodePool[thisNode].resMark=fans.find(connectResults[0]);
        //         connectResults.clear();
        //     }

        // }

        if(thisNode!=-1){
            fans.updateSize(cliqueuv[0].size());
            for(int i=0;i<runtree.treeNodePool[thisNode].sons.size();++i){
                int sonNode=runtree.treeNodePool[thisNode].sons[i];
                connectResults.push_back(runtree.treeNodePool[sonNode].resMark);
            }

            for(int i=1;i<connectResults.size();++i){
                fans.merge(connectResults[i],connectResults[0]);
            }
            runtree.treeNodePool[thisNode].resMark=fans.find(connectResults[0]);
            connectResults.clear();

            //check pqnode
            if(RUV[t].size()==threshold[t]&&t==oriSide||RUV[z].size()==threshold[z]&&z==oriSide){
                pqnodecnt++;
                pqnodedepth_total+=RUV[t].size()+RUV[z].size();
            }
        }else{
            if(RUV[t].size()==threshold[t]&&t==oriSide||RUV[z].size()==threshold[z]&&z==oriSide){
                vpqnodecnt++;
            }
        }


    }

    xStackNode.resize(oriXstackSize);
    return thisNode;
}

int BICPC::bbranchLeafRes(uint32_t deep,uint32_t side,uint32_t inValue){
    // if(S[1].sz!=g->n[1]){
    //     std::cerr<<"error"<<std::endl;
    // }
    // if(S[z].CIsEmpty()) return;

    int thisNode=-1;//temporally no tree node is built here
    int oriSide=side;

    uint32_t t=side;
    uint32_t z=side ^ 1;

    if(S[t].CIsEmpty()&&S[z].CIsEmpty()){
        if(S[t].XIsEmpty()&&S[z].XIsEmpty()){
            if(RUV[t].size()>=threshold[t]&&RUV[z].size()>=threshold[z]&&RUV[t].size()+RUV[z].size()>threshold[t]+threshold[z]){
                thisNode=runtree.getNewNode();
                runtree.treeNodePool[thisNode].inValue=inValue;
                runtree.treeNodePool[thisNode].uv=side;
                runtree.treeNodePool[thisNode].resMark=-1;
                runtree.treeNodePool[thisNode].qpMark=false;

                cliqueuv[t].emplace_back(RUV[t]);
                cliqueuv[z].emplace_back(RUV[z]);

                runtree.treeNodePool[thisNode].resMark=cliqueuv[t].size()-1;
            }
            return thisNode;
        }
        if(RUV[t].size()>=threshold[t]&&RUV[z].size()>=threshold[z]){
            if(RUV[t].size()==threshold[t]&&t==oriSide||RUV[z].size()==threshold[z]&&z==oriSide){
                connectAllBiCliqueLeafRes();
            }
        }
        return -1;
    }

    if(RUV[t].size()+S[t].cSize()<threshold[t]||RUV[z].size()+S[z].cSize()<threshold[z]){
        return thisNode;
    }

    //find pivot

    //debug
    // if(S[1].sz!=g->n[1]){
    //     std::cerr<<"error"<<std::endl;
    // }

    //initialize copies of xcr
    uint32_t x[2];
    uint32_t c[2];
    uint32_t r[2];
    x[t] = S[t].x; x[z] = S[z].x;
    c[t] = S[t].c; c[z] = S[z].c;
    r[t] = S[t].r; r[z] = S[z].r;

    int oriXstackSize=xStackNode.size();

    //find pivot
    uint32_t maxI[2] = {0}, maxE[2] = {0}, maxUV[2] = {0};
    maxI[0] = S[0].c;
    maxI[1] = S[1].c;
    maxUV[0] = S[0][S[0].c];
    maxUV[1] = S[1][S[1].c];
    t=0,z=1;
    uint32_t wsSize;


    //find pivot in C of t
    for(uint32_t i = S[t].c; i < S[t].r; i++) {//scan all vertices in C
    // for(uint32_t i = S[t].c; i <= S[t].c; i++) {
        uint32_t u = S[t][i];
        uint32_t e = 0;

        if(g->deg(u, t) > S[z].r - S[z].c) {
            for(uint32_t j = S[z].c; j < S[z].r; j++) {
                uint32_t v = S[z][j];
                if(g->connect(u, v, t)) {
                    e++;
                    deg[v]++;
                }
            }
        }
        else {
            for(uint32_t j = g->pos[t][u]; j < g->pos[t][u + 1]; j++) {
                uint32_t v = g->e[t][j];
                uint32_t pv = S[z].pos(v);
                if(S[z].c <= pv && pv < S[z].r) {//in C
                    e++;
                    deg[v]++;
                }
            }
        }

        if(e > maxE[z]) {
            maxE[z] = e;
            maxI[t] = i;
            maxUV[t] = u;
        }
    }
    //find pivot in C of z
    for(uint32_t i = S[1].c; i < S[1].r; i++) {//scan all vertices in C
    // for(uint32_t i = S[t].c; i <= S[t].c; i++) {
        uint32_t u = S[1][i];
        uint32_t e = deg[u];
        deg[u] = 0;
        if(e > maxE[0]) {
            maxE[0] = e;
            maxI[1] = i;
            maxUV[1] = u;
        }
    }

    //check if pivot in X of both sides
    bool isPivotInX[2] = {false, false};
    for(t = 0; t < 2; t++) {
        z = t ^ 1;//t + z = 1

        for(uint32_t i = S[t].x; i < S[t].c; i++) {//scan all vertices in X
        // for(uint32_t i = S[t].c; i <= S[t].c; i++) {
            uint32_t u = S[t][i];
            uint32_t e = 0;

            if(g->deg(u, t) > S[z].r - S[z].c) {
                for(uint32_t j = S[z].c; j < S[z].r; j++) {
                    uint32_t v = S[z][j];
                    if(g->connect(u, v, t)) {
                        e++;
                    }
                }
            }
            else {
                for(uint32_t j = g->pos[t][u]; j < g->pos[t][u + 1]; j++) {
                    uint32_t v = g->e[t][j];
                    uint32_t pv = S[z].pos(v);
                    if(S[z].c <= pv && pv < S[z].r) {//in C
                        e++;
                    }
                }
            }

            if(e > maxE[z]) {
                isPivotInX[t] = true;
                maxE[z] = e;
                maxI[t] = i;
                maxUV[t] = u;
            }
        }
    }

    //pivot found. pivot must have at least one neighbor in this algorithm
    t = 1;
    if(S[0].cSize() - maxE[0] >= S[1].cSize() - maxE[1]) {
        t = 0;
    }
    z = t ^ 1;


    //r is more like the right limit of C
    //r[z] doesnt follow S[z].r
    //pivot could be in X, and C of z is changed, but never mind, it will recover by r[z]
    if(g->deg(maxUV[t], t) > S[z].r - S[z].c) {
        if(S[z].r>0){
            for(uint32_t i = S[z].r - 1; i >= S[z].c; i--) {
                if(!g->connect(maxUV[t], S[z][i], t)) {
                    S[z].swapByPos(--S[z].r, i);
                }
    
                if(i == 0) break;
            }
        }
    }else{
        uint32_t u = maxUV[t];
        uint32_t tmpR = S[z].c;
        for(uint32_t j = g->pos[t][u]; j < g->pos[t][u + 1]; j++) {
            uint32_t v = g->e[t][j];
            uint32_t pv = S[z].pos(v);
            if(S[z].c <= pv && pv < S[z].r) {
                S[z].swapByPos(tmpR++, pv);
            }
        }
        S[z].r = tmpR;
    }
    wsSize = r[z] - S[z].r;
    if(ws[deep].size() < wsSize) {
        ws[deep].resize(wsSize * 2);
    }
    memcpy(ws[deep].data(), S[z].begin() + S[z].r, sizeof(uint32_t) * wsSize);


    //iterate on pivot
    if(isPivotInX[t]==false){
        
        S[t].swapByPos(maxI[t], --r[t]);
        S[t].r = r[t];

        //x[z] doesnt follow S[z].x, used to recover S[z].x
        if(g->deg(maxUV[t], t) > S[z].c - S[z].x) {
            for(uint32_t i = S[z].x; i < S[z].c; i++) {
                if(!g->connect(maxUV[t], S[z][i], t)) {
                    S[z].swapByPos(S[z].x++, i);
                }
            }
        }else{
            uint32_t u = maxUV[t];
            uint32_t tmpX = S[z].c;
            for(uint32_t j = g->pos[t][u]; j < g->pos[t][u + 1]; j++) {
                uint32_t v = g->e[t][j];
                uint32_t pv = S[z].pos(v);
                if(S[z].x <= pv && pv < S[z].c) {
                    S[z].swapByPos(--tmpX, pv);
                }
            }
            S[z].x = tmpX;
        }



        RUV[t].push_back(maxUV[t]);
        RUVlook[t][maxUV[t]]=1;

        //debug
        // if(S[1].sz!=g->n[1]){
        //     std::cerr<<"error"<<std::endl;
        // }
        // if(thisNode==401812){
        //     std::cerr<<"here"<<std::endl;
        // }

        #ifdef PARTIALPIVOT
        int res;
        if(S[z].cSize()+S[t].cSize()<=PARTIALPIVOTVAL){
            res=bbranchLeafRes_simplep(deep+1,t,maxUV[t]);
        }else{
            res=bbranchLeafRes(deep+1,t,maxUV[t]);
        }
        #else
        //C of z must not be empty
        int res=bbranchLeafRes(deep+1,t,maxUV[t]);
        #endif

        //debug
        // if(S[1].sz!=g->n[1]){
        //     std::cerr<<"error"<<std::endl;
        // }

        RUV[t].pop_back();
        RUVlook[t][maxUV[t]]=0;
        if(res!=-1){
            if(thisNode==-1){
                thisNode=runtree.getNewNode();
                runtree.treeNodePool[thisNode].inValue=inValue;
                runtree.treeNodePool[thisNode].uv=side;
                runtree.treeNodePool[thisNode].resMark=-1;
                runtree.treeNodePool[thisNode].qpMark=false;
            }
            runtree.treeNodePool[thisNode].sons.push_back(res);
            xStackNode.push_back(res);
        }

        //debug
        // if(S[1].sz!=g->n[1]){
        //     std::cerr<<"error"<<std::endl;
        // }

        //put pivot into X
        S[t].swapByPos(c[t]++, r[t]++);
    }




    // if(c[t]<r[t]){//S[t].c>0
        //iterate on vertices beyound pivot
        for(uint32_t j = 0; j < wsSize; j++){
            uint32_t w = ws[deep][j];
            S[z].swapByPos(S[z].pos(w), --r[z]);

            //recover next level parameter
            S[t].x = x[t]; S[z].x = x[z];
            S[t].c = c[t]; S[z].c = c[z];
            S[t].r = r[t]; S[z].r = r[z];

            if(g->deg(w, z) > S[t].c - S[t].x) {
                for(uint32_t i = S[t].x; i < S[t].c; i++) {
                    if(!g->connect(w, S[t][i], z)) {
                        S[t].swapByPos(S[t].x++, i);
                    }
                }
            }
            else {
                uint32_t tmpX = S[t].c;
                for(uint32_t j = g->pos[z][w]; j < g->pos[z][w + 1]; j++) {
                    uint32_t u = g->e[z][j];
                    uint32_t pu = S[t].pos(u);
                    if(S[t].x <= pu && pu < S[t].c) {
                        S[t].swapByPos(--tmpX, pu);
                    }
                }
                S[t].x = tmpX;
            }

            if(g->deg(w, z) > S[t].r - S[t].c) {
                if(S[t].r > 0){
                    for(uint32_t i = S[t].r - 1; i >= S[t].c; i--) {
                        if(!g->connect(w, S[t][i], z)) {
                            S[t].swapByPos(--S[t].r, i);
                        }

                        if(i == 0) break;
                    }
                }
            }else{
                uint32_t tmpR = S[t].c;
                for(uint32_t j = g->pos[z][w]; j < g->pos[z][w + 1]; j++) {
                    uint32_t u = g->e[z][j];
                    uint32_t pu = S[t].pos(u);
                    if(S[t].c <= pu && pu < S[t].r) {
                        S[t].swapByPos(tmpR++, pu);
                    }
                }
                S[t].r = tmpR;
            }

            // //debug
            // if(deep==0&&inValue==1859991&&w==1911392){
            //     // std::cerr<<<<std::endl;
            //     for(int ii=S[z].c;ii<S[z].r;++ii){
            //         std::cerr<<g->old_lables[z][S[z][ii]]<<" ";
            //     }
            //     std::cerr<<std::endl;
                
            // }

            // if(S[t].cSize()>0){
                RUV[z].push_back(w);
                RUVlook[z][w]=1;
                //debug
                // if(S[1].sz!=g->n[1]){
                //     std::cerr<<"error"<<std::endl;
                // }

                #ifdef PARTIALPIVOT
                int res;
                if(S[z].cSize()+S[t].cSize()<=PARTIALPIVOTVAL){
                    res=bbranchLeafRes_simplep(deep+1,z,w);
                }else{
                    res=bbranchLeafRes(deep+1,z,w);
                }
                #else
                //C of z must not be empty
                int res=bbranchLeafRes(deep+1,z,w);
                #endif

                RUV[z].pop_back();
                RUVlook[z][w]=0;

                if(res!=-1){
                    if(thisNode==-1){
                        thisNode=runtree.getNewNode();
                        runtree.treeNodePool[thisNode].inValue=inValue;
                        runtree.treeNodePool[thisNode].uv=side;
                        runtree.treeNodePool[thisNode].resMark=-1;
                        runtree.treeNodePool[thisNode].qpMark=false;
                    }
                    runtree.treeNodePool[thisNode].sons.push_back(res);
                    xStackNode.push_back(res);
                }
            // }

            S[z].swapByPos(c[z]++,r[z]++);//move w to X
        }

        for(uint32_t j = 0; j < wsSize; j++) {
            uint32_t w = ws[deep][j];
            S[z].swapByPos(S[z].pos(w), --c[z]);
        }
    // }

    //recover pivot
    if(isPivotInX[t] == false) {
        S[t].swapByPos(S[t].pos(maxUV[t]), --c[t]);
    }

    //recover x,c,r
    S[t].x = x[t]; S[z].x = x[z];
    S[t].c = c[t]; S[z].c = c[z];
    S[t].r = r[t]; S[z].r = r[z];

    //connect that should be connected
    //almost the same as before
    if(RUV[t].size()>=threshold[t]&&RUV[z].size()>=threshold[z]){
        int res=-1;
        if(RUV[t].size()==threshold[t]&&t==oriSide||RUV[z].size()==threshold[z]&&z==oriSide){// only the first one
            res=connectAllBiCliqueLeafRes();
            if(thisNode!=-1){
                runtree.treeNodePool[thisNode].resMark=res;
            }
        }else{
            if(thisNode!=-1){
                fans.updateSize(cliqueuv[0].size());
                for(int i=0;i<runtree.treeNodePool[thisNode].sons.size();++i){
                    int sonNode=runtree.treeNodePool[thisNode].sons[i];
                    connectResults.push_back(runtree.treeNodePool[sonNode].resMark);
                }

                for(int i=1;i<connectResults.size();++i){
                    fans.merge(connectResults[i],connectResults[0]);
                }
                runtree.treeNodePool[thisNode].resMark=fans.find(connectResults[0]);
                connectResults.clear();
            }

        }


    }

    xStackNode.resize(oriXstackSize);
    return thisNode;
}

int BICPC::bbranch(uint32_t deep,uint32_t side,uint32_t inValue){
    // if(S[1].sz!=g->n[1]){
    //     std::cerr<<"error"<<std::endl;
    // }
    // if(S[z].CIsEmpty()) return;

    int thisNode=-1;//temporally no tree node is built here
    int oriSide=side;

    //TASK1
    //Check if the current R forms a answer with the CAND on the opposite side of newly chosen vertex
    bool found=true;
    int curCliqueId=-1;
    uint32_t t=side;
    uint32_t z=side ^ 1;
    if(S[z].XIsEmpty()&&RUV[t].size()>=threshold[t]&&(RUV[z].size()+S[z].r-S[z].c>=threshold[z])){
        if(RUV[t].size()+RUV[z].size()+S[z].r-S[z].c>threshold[t]+threshold[z]){
            for(uint32_t i = S[t].x; i < S[t].r; i++) {
                //zy
                // printf("here\n");

                bool connectAll = true;
                for(uint32_t j = S[z].c; j < S[z].r; j++) {
                    if(!g->connect(S[t][i], S[z][j], t)) {
                        connectAll = false;
                        break;
                    }
                }
                if(connectAll) {
                    found = false;
                    break;
                }
            }

            if(found){//found a maximal biclique: RUV and cand of z
                thisNode=runtree.getNewNode();
                runtree.treeNodePool[thisNode].inValue=inValue;
                runtree.treeNodePool[thisNode].uv=side;
                runtree.treeNodePool[thisNode].resMark=-1;
                runtree.treeNodePool[thisNode].qpMark=false;

                cliqueuv[t].emplace_back(RUV[t]);
                int oldsize=RUV[z].size();
                RUV[z].reserve(RUV[z].size()+S[z].r-S[z].c);
                for(int i=S[z].c;i<S[z].r;++i){
                    RUV[z].push_back(S[z][i]);
                }
                cliqueuv[z].emplace_back(RUV[z]);
                // runtree.treeNodePool[thisNode].resMark=cliqueuv[t].size()-1;

                //reset RUV
                RUV[z].resize(oldsize);

                runtree.treeNodePool[thisNode].resMark=cliqueuv[t].size()-1;
                curCliqueId=cliqueuv[t].size()-1;
            }
        }
    }


    //TASK2, make sure C of t,z have edge(s), by shrink C of z
    uint32_t newR = S[z].r;
    uint32_t oldR = S[z].r;
    for(uint32_t i = S[z].c; i < newR; ) {
        bool degZero = true;
        for(uint32_t j = S[t].c; j < S[t].r; j++) {
            if(g->connect(S[z][i], S[t][j], z)) {
                degZero = false;
                break;
            }
        }
        if(degZero) {
            S[z].swapByPos(i, --newR);
        }
        else i++;
    }
    S[z].r = newR;
    bool emptyc=(S[t].cSize()==0||S[z].cSize()==0);

    if(emptyc){
        // #ifdef LESSCONNECT
        // if(RUV[t].size()>=threshold[t]&&RUV[z].size()>=threshold[z]){
        //     if(RUV[t].size()==threshold[t]&&t==oriSide||RUV[z].size()==threshold[z]&&z==oriSide){
        //         int res=connectAllBiCliqueNew(oriSide,thisNode);
        //         if(thisNode!=-1){
        //             runtree.treeNodePool[thisNode].qpMark=true;
        //             runtree.treeNodePool[thisNode].resMark=res;
        //         }
        //         connectResults.clear();
        //     }
        // }

        // if(thisNode!=-1){
        //     if(RUV[t].size()>=threshold[t]&&RUV[z].size()>=threshold[z]){
        //         if(RUV[t].size()==threshold[t]&&t==oriSide||RUV[z].size()==threshold[z]&&z==oriSide){
        //             int res=connectAllBiCliqueNew(oriSide,thisNode);
        //         }
        //     }
        // }
        // #else
        if(RUV[0].size()>=threshold[0]&&RUV[1].size()>=threshold[1]){
            S[z].r=oldR;

            int res=connectAllBiCliqueNew(oriSide,thisNode);

            if(thisNode!=-1){
                runtree.treeNodePool[thisNode].qpMark=true;
                runtree.treeNodePool[thisNode].resMark=res;
            }

            connectResults.clear();
            S[z].r=newR;
        }
        // #endif
        return thisNode;
    }

    //TASK3, find pivot and start new recursions
    // if(S[t].cSize()==0||S[z].cSize()==0){
    //     //to connect
    //     if(RUV[0].size()>=p&&RUV[1].size()>=q){
    //         if(thisNode!=-1){
    //             connectResults.push_back(runtree.treeNodePool[thisNode].resMark);
    //         }
    //         int res=connectBiCliqueInCurZone(oriSide);
    //         if(res!=-1){
    //             connectResults.push_back(res);
    //         }
    //         if((S[0].x<S[0].c||S[1].x<S[1].c)){
    //             // fans.updateSize(cliqueuv[0].size());

    //             // if(curCliqueId!=-1) curStackNode.push_back(curCliqueId);
    //             connectBiCliqueInXzone(oriSide);

    //             // if(curCliqueId!=-1) curStackNode.pop_back();

    //         }

    //         if(connectResults.size()==0){
    //             std::cerr<<"final connect id error"<<std::endl;
    //             exit(0);
    //         }

    //         fans.updateSize(cliqueuv[0].size());
    //         if(connectResults.size()>1){
    //             for(int i=1;i<connectResults.size();++i){
    //                 fans.merge(connectResults[0],connectResults[i]);
    //             }
    //         }

    //         if(thisNode!=-1){
    //             // if(runtree.treeNodePool[thisNode].resMark==-1){
    //             //     std::cerr<<"resmark error"<<std::endl;
    //             //     exit(0);
    //             // }

    //             // if(res!=-1){
    //             //     fans.merge(runtree.treeNodePool[thisNode].resMark,res);
    //             //     runtree.treeNodePool[thisNode].resMark=fans.find(res);
    //             // }

    //             runtree.treeNodePool[thisNode].qpMark=true;
    //             runtree.treeNodePool[thisNode].resMark=fans.find(connectResults[0]);
    //         }

    //         connectResults.clear();

    //     }
    //     return thisNode;
    // }

    if(RUV[t].size()+S[t].cSize()<threshold[t]||RUV[z].size()+S[z].cSize()<threshold[z]){
        return thisNode;
    }

    //debug
    // if(S[1].sz!=g->n[1]){
    //     std::cerr<<"error"<<std::endl;
    // }

    #ifdef NOSHRINK
    S[z].r=oldR;
    #endif


    //initialize copies of xcr
    uint32_t x[2];
    uint32_t c[2];
    uint32_t r[2];
    x[t] = S[t].x; x[z] = S[z].x;
    c[t] = S[t].c; c[z] = S[z].c;
    r[t] = S[t].r; r[z] = S[z].r;

    int oriXstackSize=xStackNode.size();

    //find pivot
    uint32_t maxI[2] = {0}, maxE[2] = {0}, maxUV[2] = {0};
    maxI[0] = S[0].c;
    maxI[1] = S[1].c;
    maxUV[0] = S[0][S[0].c];
    maxUV[1] = S[1][S[1].c];
    t=0,z=1;
    uint32_t wsSize;


    //find pivot in C of t
    for(uint32_t i = S[t].c; i < S[t].r; i++) {//scan all vertices in C
    // for(uint32_t i = S[t].c; i <= S[t].c; i++) {
        uint32_t u = S[t][i];
        uint32_t e = 0;

        if(g->deg(u, t) > S[z].r - S[z].c) {
            for(uint32_t j = S[z].c; j < S[z].r; j++) {
                uint32_t v = S[z][j];
                if(g->connect(u, v, t)) {
                    e++;
                    deg[v]++;
                }
            }
        }
        else {
            for(uint32_t j = g->pos[t][u]; j < g->pos[t][u + 1]; j++) {
                uint32_t v = g->e[t][j];
                uint32_t pv = S[z].pos(v);
                if(S[z].c <= pv && pv < S[z].r) {//in C
                    e++;
                    deg[v]++;
                }
            }
        }

        if(e > maxE[z]) {
            maxE[z] = e;
            maxI[t] = i;
            maxUV[t] = u;
        }
    }
    //find pivot in C of z
    for(uint32_t i = S[1].c; i < S[1].r; i++) {//scan all vertices in C
    // for(uint32_t i = S[t].c; i <= S[t].c; i++) {
        uint32_t u = S[1][i];
        uint32_t e = deg[u];
        deg[u] = 0;
        if(e > maxE[0]) {
            maxE[0] = e;
            maxI[1] = i;
            maxUV[1] = u;
        }
    }

    //check if pivot in X of both sides
    bool isPivotInX[2] = {false, false};
    for(t = 0; t < 2; t++) {
        z = t ^ 1;//t + z = 1

        for(uint32_t i = S[t].x; i < S[t].c; i++) {//scan all vertices in X
        // for(uint32_t i = S[t].c; i <= S[t].c; i++) {
            uint32_t u = S[t][i];
            uint32_t e = 0;

            if(g->deg(u, t) > S[z].r - S[z].c) {
                for(uint32_t j = S[z].c; j < S[z].r; j++) {
                    uint32_t v = S[z][j];
                    if(g->connect(u, v, t)) {
                        e++;
                    }
                }
            }
            else {
                for(uint32_t j = g->pos[t][u]; j < g->pos[t][u + 1]; j++) {
                    uint32_t v = g->e[t][j];
                    uint32_t pv = S[z].pos(v);
                    if(S[z].c <= pv && pv < S[z].r) {//in C
                        e++;
                    }
                }
            }

            if(e > maxE[z]) {
                isPivotInX[t] = true;
                maxE[z] = e;
                maxI[t] = i;
                maxUV[t] = u;
            }
        }
    }

    //pivot found. pivot must have at least one neighbor in this algorithm
    t = 1;
    if(S[0].cSize() - maxE[0] >= S[1].cSize() - maxE[1]) {
        t = 0;
    }
    z = t ^ 1;


    //r is more like the right limit of C
    //r[z] doesnt follow S[z].r
    //pivot could be in X, and C of z is changed, but never mind, it will recover by r[z]
    if(g->deg(maxUV[t], t) > S[z].r - S[z].c) {
        for(uint32_t i = S[z].r - 1; i >= S[z].c; i--) {
            if(!g->connect(maxUV[t], S[z][i], t)) {
                S[z].swapByPos(--S[z].r, i);
            }

            if(i == 0) break;
        }
    }else{
        uint32_t u = maxUV[t];
        uint32_t tmpR = S[z].c;
        for(uint32_t j = g->pos[t][u]; j < g->pos[t][u + 1]; j++) {
            uint32_t v = g->e[t][j];
            uint32_t pv = S[z].pos(v);
            if(S[z].c <= pv && pv < S[z].r) {
                S[z].swapByPos(tmpR++, pv);
            }
        }
        S[z].r = tmpR;
    }
    wsSize = r[z] - S[z].r;
    if(ws[deep].size() < wsSize) {
        ws[deep].resize(wsSize * 2);
    }
    memcpy(ws[deep].data(), S[z].begin() + S[z].r, sizeof(uint32_t) * wsSize);


    //iterate on pivot
    if(isPivotInX[t]==false){
        
        S[t].swapByPos(maxI[t], --r[t]);
        S[t].r = r[t];

        //x[z] doesnt follow S[z].x, used to recover S[z].x
        if(g->deg(maxUV[t], t) > S[z].c - S[z].x) {
            for(uint32_t i = S[z].x; i < S[z].c; i++) {
                if(!g->connect(maxUV[t], S[z][i], t)) {
                    S[z].swapByPos(S[z].x++, i);
                }
            }
        }else{
            uint32_t u = maxUV[t];
            uint32_t tmpX = S[z].c;
            for(uint32_t j = g->pos[t][u]; j < g->pos[t][u + 1]; j++) {
                uint32_t v = g->e[t][j];
                uint32_t pv = S[z].pos(v);
                if(S[z].x <= pv && pv < S[z].c) {
                    S[z].swapByPos(--tmpX, pv);
                }
            }
            S[z].x = tmpX;
        }



        RUV[t].push_back(maxUV[t]);
        RUVlook[t][maxUV[t]]=1;
        if(curCliqueId!=-1) curStackNode.push_back(curCliqueId);

        //debug
        // if(S[1].sz!=g->n[1]){
        //     std::cerr<<"error"<<std::endl;
        // }
        // if(thisNode==401812){
        //     std::cerr<<"here"<<std::endl;
        // }

        #ifdef PARTIALPIVOT
        int res;
        if(S[z].cSize()+S[t].cSize()<=PARTIALPIVOTVAL){
            res=bbranch_simplep(deep+1,t,maxUV[t]);
        }else{
            res=bbranch(deep+1,t,maxUV[t]);
        }
        #else
        //C of z must not be empty
        int res=bbranch(deep+1,t,maxUV[t]);
        #endif

        //debug
        // if(S[1].sz!=g->n[1]){
        //     std::cerr<<"error"<<std::endl;
        // }

        RUV[t].pop_back();
        RUVlook[t][maxUV[t]]=0;
        if(curCliqueId!=-1) curStackNode.pop_back();
        if(res!=-1){
            if(thisNode==-1){
                thisNode=runtree.getNewNode();
                runtree.treeNodePool[thisNode].inValue=inValue;
                runtree.treeNodePool[thisNode].uv=side;
                runtree.treeNodePool[thisNode].resMark=-1;
                runtree.treeNodePool[thisNode].qpMark=false;
            }
            runtree.treeNodePool[thisNode].sons.push_back(res);
            xStackNode.push_back(res);
        }

        //debug
        // if(S[1].sz!=g->n[1]){
        //     std::cerr<<"error"<<std::endl;
        // }

        //put pivot into X
        S[t].swapByPos(c[t]++, r[t]++);
    }




    if(c[t]<r[t]){//S[t].c>0
        //iterate on vertices beyound pivot
        for(uint32_t j = 0; j < wsSize; j++){
            uint32_t w = ws[deep][j];
            S[z].swapByPos(S[z].pos(w), --r[z]);

            //recover next level parameter
            S[t].x = x[t]; S[z].x = x[z];
            S[t].c = c[t]; S[z].c = c[z];
            S[t].r = r[t]; S[z].r = r[z];

            if(g->deg(w, z) > S[t].c - S[t].x) {
                for(uint32_t i = S[t].x; i < S[t].c; i++) {
                    if(!g->connect(w, S[t][i], z)) {
                        S[t].swapByPos(S[t].x++, i);
                    }
                }
            }
            else {
                uint32_t tmpX = S[t].c;
                for(uint32_t j = g->pos[z][w]; j < g->pos[z][w + 1]; j++) {
                    uint32_t u = g->e[z][j];
                    uint32_t pu = S[t].pos(u);
                    if(S[t].x <= pu && pu < S[t].c) {
                        S[t].swapByPos(--tmpX, pu);
                    }
                }
                S[t].x = tmpX;
            }

            if(g->deg(w, z) > S[t].r - S[t].c) {
                if(S[t].r > 0){
                    for(uint32_t i = S[t].r - 1; i >= S[t].c; i--) {
                        if(!g->connect(w, S[t][i], z)) {
                            S[t].swapByPos(--S[t].r, i);
                        }

                        if(i == 0) break;
                    }
                }
            }else{
                uint32_t tmpR = S[t].c;
                for(uint32_t j = g->pos[z][w]; j < g->pos[z][w + 1]; j++) {
                    uint32_t u = g->e[z][j];
                    uint32_t pu = S[t].pos(u);
                    if(S[t].c <= pu && pu < S[t].r) {
                        S[t].swapByPos(tmpR++, pu);
                    }
                }
                S[t].r = tmpR;
            }

            // //debug
            // if(deep==0&&inValue==1859991&&w==1911392){
            //     // std::cerr<<<<std::endl;
            //     for(int ii=S[z].c;ii<S[z].r;++ii){
            //         std::cerr<<g->old_lables[z][S[z][ii]]<<" ";
            //     }
            //     std::cerr<<std::endl;
                
            // }

            if(S[t].cSize()>0){
                RUV[z].push_back(w);
                RUVlook[z][w]=1;
                if(curCliqueId!=-1) curStackNode.push_back(curCliqueId);
                //debug
                // if(S[1].sz!=g->n[1]){
                //     std::cerr<<"error"<<std::endl;
                // }

                #ifdef PARTIALPIVOT
                int res;
                if(S[z].cSize()+S[t].cSize()<=PARTIALPIVOTVAL){
                    res=bbranch_simplep(deep+1,z,w);
                }else{
                    res=bbranch(deep+1,z,w);
                }
                #else
                //C of z must not be empty
                int res=bbranch(deep+1,z,w);
                #endif

                if(curCliqueId!=-1) curStackNode.pop_back();
                RUV[z].pop_back();
                RUVlook[z][w]=0;

                if(res!=-1){
                    if(thisNode==-1){
                        thisNode=runtree.getNewNode();
                        runtree.treeNodePool[thisNode].inValue=inValue;
                        runtree.treeNodePool[thisNode].uv=side;
                        runtree.treeNodePool[thisNode].resMark=-1;
                        runtree.treeNodePool[thisNode].qpMark=false;
                    }
                    runtree.treeNodePool[thisNode].sons.push_back(res);
                    xStackNode.push_back(res);
                }
            }

            S[z].swapByPos(c[z]++,r[z]++);//move w to X
        }

        for(uint32_t j = 0; j < wsSize; j++) {
            uint32_t w = ws[deep][j];
            S[z].swapByPos(S[z].pos(w), --c[z]);
        }
    }

    //recover pivot
    if(isPivotInX[t] == false) {
        S[t].swapByPos(S[t].pos(maxUV[t]), --c[t]);
    }

    //recover x,c,r
    S[t].x = x[t]; S[z].x = x[z];
    S[t].c = c[t]; S[z].c = c[z];
    S[t].r = r[t]; S[z].r = r[z];

    //connect that should be connected
    //almost the same as before
    if(RUV[0].size()>=p&&RUV[1].size()>=q){
        // fans.updateSize(cliqueuv[0].size());
        S[oriSide^1].r=oldR;
        // if(thisNode!=-1&&runtree.treeNodePool[thisNode].resMark!=-1){
        //     connectResults.push_back(runtree.treeNodePool[thisNode].resMark);
        // }
        // int res=connectBiCliqueInCurZone(oriSide);
        // if(res!=-1){
        //     connectResults.push_back(res);
        // }

        // if(S[oriSide^1].x<S[oriSide^1].c||)
        // connectBiCliqueInXzone(oriSide);

        // if(connectResults.size()==0){
        //     std::cerr<<"final connect id error"<<std::endl;
        //     exit(0);
        // }

        // fans.updateSize(cliqueuv[0].size());
        // if(connectResults.size()>1){
        //     for(int i=1;i<connectResults.size();++i){
        //         fans.merge(connectResults[0],connectResults[i]);
        //     }
        // }
        int res=connectAllBiCliqueNew(oriSide,thisNode);

        if(thisNode!=-1){
            // if(res==-1){
            //     std::cerr<<"error definitely"<<std::endl;
            //     exit(0);
            // }
            runtree.treeNodePool[thisNode].qpMark=true;
            runtree.treeNodePool[thisNode].resMark=res;
        }

        connectResults.clear();
        S[oriSide^1].r=newR;
    }

    xStackNode.resize(oriXstackSize);
    return thisNode;

    
}

void BICPC::bbranchpureFullTree(uint32_t deep,uint32_t side,uint32_t inValue,int fid){
    // if(S[z].CIsEmpty()) return;

    // int thisNode=-1;//temporally no tree node is built here
    uint32_t t=side;
    uint32_t z=side ^ 1;

    int thisTreeId=mbeTreeNodeCnt++;


    if(S[t].CIsEmpty()&&S[z].CIsEmpty()){
        if(S[t].XIsEmpty()&&S[z].XIsEmpty()&&RUV[t].size()>=threshold[t]&&RUV[z].size()>=threshold[z]){
            if(RUV[t].size()+RUV[z].size()>threshold[t]+threshold[z]){
                cliqueuv[t].emplace_back(RUV[t]);
                cliqueuv[z].emplace_back(RUV[z]);
            }
        }
        
        reportTreeNode(thisTreeId,fid,-1,-1);
        return;
    }

    //TASK3, find pivot and start new recursions

    if(RUV[t].size()+S[t].cSize()<threshold[t]||RUV[z].size()+S[z].cSize()<threshold[z]){
        reportTreeNode(thisTreeId,fid,-1,-1);
        return;
    }


    // #ifdef NOSHRINK
    // S[z].r=oldR;
    // #endif


    //find pivot
    uint32_t maxI[2] = {0}, maxE[2] = {0}, maxUV[2] = {0};
    maxI[0] = S[0].c;
    maxI[1] = S[1].c;
    maxUV[0] = S[0][S[0].c];
    maxUV[1] = S[1][S[1].c];
    t=0,z=1;
    //find pivot in C of t
    for(uint32_t i = S[t].c; i < S[t].r; i++) {//scan all vertices in C
    // for(uint32_t i = S[t].c; i <= S[t].c; i++) {
        uint32_t u = S[t][i];
        uint32_t e = 0;

        if(g->deg(u, t) > S[z].r - S[z].c) {
            for(uint32_t j = S[z].c; j < S[z].r; j++) {
                uint32_t v = S[z][j];
                if(g->connect(u, v, t)) {
                    e++;
                    deg[v]++;
                }
            }
        }
        else {
            for(uint32_t j = g->pos[t][u]; j < g->pos[t][u + 1]; j++) {
                uint32_t v = g->e[t][j];
                uint32_t pv = S[z].pos(v);
                if(S[z].c <= pv && pv < S[z].r) {//in C
                    e++;
                    deg[v]++;
                }
            }
        }

        if(e > maxE[z]) {
            maxE[z] = e;
            maxI[t] = i;
            maxUV[t] = u;
        }
    }
    //find pivot in C of z
    for(uint32_t i = S[1].c; i < S[1].r; i++) {//scan all vertices in C
    // for(uint32_t i = S[t].c; i <= S[t].c; i++) {
        uint32_t u = S[1][i];
        uint32_t e = deg[u];
        deg[u] = 0;
        if(e > maxE[0]) {
            maxE[0] = e;
            maxI[1] = i;
            maxUV[1] = u;
        }
    }

    //check if pivot in X of both sides
    bool isPivotInX[2] = {false, false};
    for(t = 0; t < 2; t++) {
        z = t ^ 1;//t + z = 1

        for(uint32_t i = S[t].x; i < S[t].c; i++) {//scan all vertices in X
        // for(uint32_t i = S[t].c; i <= S[t].c; i++) {
            uint32_t u = S[t][i];
            uint32_t e = 0;

            if(g->deg(u, t) > S[z].r - S[z].c) {
                for(uint32_t j = S[z].c; j < S[z].r; j++) {
                    uint32_t v = S[z][j];
                    if(g->connect(u, v, t)) {
                        e++;
                    }
                }
            }
            else {
                for(uint32_t j = g->pos[t][u]; j < g->pos[t][u + 1]; j++) {
                    uint32_t v = g->e[t][j];
                    uint32_t pv = S[z].pos(v);
                    if(S[z].c <= pv && pv < S[z].r) {//in C
                        e++;
                    }
                }
            }

            if(e > maxE[z]) {
                isPivotInX[t] = true;
                maxE[z] = e;
                maxI[t] = i;
                maxUV[t] = u;
            }
        }
    }

    //pivot found.
    t = 1;
    if(S[0].cSize() - maxE[0] >= S[1].cSize() - maxE[1]) {
        t = 0;
    }
    z = t ^ 1;


    //report tree node
    reportTreeNode(thisTreeId,fid,maxUV[t],t);


    //initialize copies of xcr
    uint32_t x[2];
    uint32_t c[2];
    uint32_t r[2];
    x[t] = S[t].x; x[z] = S[z].x;
    c[t] = S[t].c; c[z] = S[z].c;
    r[t] = S[t].r; r[z] = S[z].r;

    //r is more like the right limit of C
    //r[z] doesnt follow S[z].r
    //pivot could be in X, and C of z is changed, but never mind, it will recover by r[z]
    if(g->deg(maxUV[t], t) > S[z].r - S[z].c) {
        if(S[z].r>0){
            for(uint32_t i = S[z].r - 1; i >= S[z].c; i--) {
                if(!g->connect(maxUV[t], S[z][i], t)) {
                    S[z].swapByPos(--S[z].r, i);
                }
    
                if(i == 0) break;
            }
        }
    }else{
        uint32_t u = maxUV[t];
        uint32_t tmpR = S[z].c;
        for(uint32_t j = g->pos[t][u]; j < g->pos[t][u + 1]; j++) {
            uint32_t v = g->e[t][j];
            uint32_t pv = S[z].pos(v);
            if(S[z].c <= pv && pv < S[z].r) {
                S[z].swapByPos(tmpR++, pv);
            }
        }
        S[z].r = tmpR;
    }

    uint32_t wsSize = r[z] - S[z].r;
    if(ws[deep].size() < wsSize) {
        ws[deep].resize(wsSize * 2);
    }
    memcpy(ws[deep].data(), S[z].begin() + S[z].r, sizeof(uint32_t) * wsSize);

    // //debug
    // if(deep==0&&inValue==1859991){
    //     std::cerr<<"t: "<<t<<" pivot: "<<g->old_lables[t][maxUV[t]]<<std::endl;
    //     std::cerr<<wsSize<<std::endl;
    //     for(int i=0;i<ws[deep].size();++i){
    //         std::cerr<<g->old_lables[z][ws[deep][i]]<<" ";
    //     }
    //     std::cerr<<std::endl;
    // }
    // //debug
    // if(deep==0&&inValue==1859991){
    //     for(int ii=S[t].c;ii<S[t].r;++ii){
    //         std::cerr<<g->old_lables[t][S[t][ii]]<<" ";
    //     }
    //     std::cerr<<std::endl;
    //     for(int ii=S[z].c;ii<S[z].r;++ii){
    //         std::cerr<<g->old_lables[z][S[z][ii]]<<" ";
    //     }
    //     std::cerr<<std::endl;
    // }

    // int oriXstackSize=xStackNode.size();
    //iterate on pivot
    if(isPivotInX[t]==false){
        
        S[t].swapByPos(maxI[t], --r[t]);
        S[t].r = r[t];

        //x[z] doesnt follow S[z].x, used to recover S[z].x
        if(g->deg(maxUV[t], t) > S[z].c - S[z].x) {
            for(uint32_t i = S[z].x; i < S[z].c; i++) {
                if(!g->connect(maxUV[t], S[z][i], t)) {
                    S[z].swapByPos(S[z].x++, i);
                }
            }
        }else{
            uint32_t u = maxUV[t];
            uint32_t tmpX = S[z].c;
            for(uint32_t j = g->pos[t][u]; j < g->pos[t][u + 1]; j++) {
                uint32_t v = g->e[t][j];
                uint32_t pv = S[z].pos(v);
                if(S[z].x <= pv && pv < S[z].c) {
                    S[z].swapByPos(--tmpX, pv);
                }
            }
            S[z].x = tmpX;
        }

        RUV[t].push_back(maxUV[t]);
        // RUVlook[t][maxUV[t]]=1;
        
        bbranchpureFullTree(deep+1,t,maxUV[t],thisTreeId);
        RUV[t].pop_back();
        // RUVlook[t][maxUV[t]]=0;
        // if(res!=-1){
        //     if(thisNode==-1){
        //         thisNode=runtree.getNewNode();
        //         runtree.treeNodePool[thisNode].inValue=inValue;
        //         runtree.treeNodePool[thisNode].uv=side;
        //         runtree.treeNodePool[thisNode].resMark=-1;
        //     }
        //     runtree.treeNodePool[thisNode].sons.push_back(res);
        //     xStackNode.push_back(res);
        // }

        //put pivot into X
        S[t].swapByPos(c[t]++, r[t]++);
        

    }

    // if(deep==0&&inValue==1859991){
    //     std::cerr<<S[t].cSize()<<std::endl;
    // }

    // if(c[t]<r[t]){//S[t].c>0
        //iterate on vertices beyound pivot
        for(uint32_t j = 0; j < wsSize; j++){
            uint32_t w = ws[deep][j];
            S[z].swapByPos(S[z].pos(w), --r[z]);

            //recover next level parameter
            S[t].x = x[t]; S[z].x = x[z];
            S[t].c = c[t]; S[z].c = c[z];
            S[t].r = r[t]; S[z].r = r[z];

            if(g->deg(w, z) > S[t].c - S[t].x) {
                for(uint32_t i = S[t].x; i < S[t].c; i++) {
                    if(!g->connect(w, S[t][i], z)) {
                        S[t].swapByPos(S[t].x++, i);
                    }
                }
            }
            else {
                uint32_t tmpX = S[t].c;
                for(uint32_t j = g->pos[z][w]; j < g->pos[z][w + 1]; j++) {
                    uint32_t u = g->e[z][j];
                    uint32_t pu = S[t].pos(u);
                    if(S[t].x <= pu && pu < S[t].c) {
                        S[t].swapByPos(--tmpX, pu);
                    }
                }
                S[t].x = tmpX;
            }

            if(g->deg(w, z) > S[t].c - S[t].r) {
                if(S[t].r > 0){
                    for(uint32_t i = S[t].r - 1; i >= S[t].c; i--) {
                        if(!g->connect(w, S[t][i], z)) {
                            S[t].swapByPos(--S[t].r, i);
                        }

                        if(i == 0) break;
                    }
                }
            }else{
                uint32_t tmpR = S[t].c;
                for(uint32_t j = g->pos[z][w]; j < g->pos[z][w + 1]; j++) {
                    uint32_t u = g->e[z][j];
                    uint32_t pu = S[t].pos(u);
                    if(S[t].c <= pu && pu < S[t].r) {
                        S[t].swapByPos(tmpR++, pu);
                    }
                }
                S[t].r = tmpR;
            }

            // //debug
            // if(deep==0&&inValue==1859991&&w==1911392){
            //     // std::cerr<<<<std::endl;
            //     for(int ii=S[z].c;ii<S[z].r;++ii){
            //         std::cerr<<g->old_lables[z][S[z][ii]]<<" ";
            //     }
            //     std::cerr<<std::endl;
                
            // }

            // if(S[t].cSize()>0){
                RUV[z].push_back(w);
                // RUVlook[z][w]=1;
                bbranchpureFullTree(deep+1,z,w,thisTreeId);
                RUV[z].pop_back();
                // RUVlook[z][w]=0;

                // if(res!=-1){
                //     if(thisNode==-1){
                //         thisNode=runtree.getNewNode();
                //         runtree.treeNodePool[thisNode].inValue=inValue;
                //         runtree.treeNodePool[thisNode].uv=side;
                //         runtree.treeNodePool[thisNode].resMark=-1;
                //     }
                //     runtree.treeNodePool[thisNode].sons.push_back(res);
                //     xStackNode.push_back(res);
                // }
            // }

            S[z].swapByPos(c[z]++,r[z]++);//move w to X
        }

        for(uint32_t j = 0; j < wsSize; j++) {
            uint32_t w = ws[deep][j];
            S[z].swapByPos(S[z].pos(w), --c[z]);
        }
    // }

    //recover pivot
    if(isPivotInX[t] == false) {
        S[t].swapByPos(S[t].pos(maxUV[t]), --c[t]);
    }

    //recover x,c,r
    S[t].x = x[t]; S[z].x = x[z];
    S[t].c = c[t]; S[z].c = c[z];
    S[t].r = r[t]; S[z].r = r[z];

    //connect that should be connected
    // if(RUV[0].size()>=p&&RUV[1].size()>=q){
    //     fans.updateSize(cliqueuv[0].size());

    //     int res=connectBiCliques();
    //     if(thisNode!=-1){
    //         if(runtree.treeNodePool[thisNode].resMark==-1){
    //             runtree.treeNodePool[thisNode].resMark=res;
    //         }else{
    //             fans.merge(runtree.treeNodePool[thisNode].resMark,res);
    //             runtree.treeNodePool[thisNode].resMark=fans.find(res);
    //         }
    //     }
    // }

    // xStackNode.resize(oriXstackSize);
    return;
}

void BICPC::bbranchpure(uint32_t deep,uint32_t side,uint32_t inValue){
    // if(S[z].CIsEmpty()) return;

    // int thisNode=-1;//temporally no tree node is built here

    //TASK1
    //Check if the current R forms a answer with the CAND on the opposite side of newly chosen vertex
    bool found=true;
    uint32_t t=side;
    uint32_t z=side ^ 1;
    if(S[z].XIsEmpty()&&RUV[t].size()>=threshold[t]&&(RUV[z].size()+S[z].r-S[z].c>=threshold[z])){
        if(RUV[t].size()+RUV[z].size()+S[z].r-S[z].c>threshold[t]+threshold[z]){
            for(uint32_t i = S[t].x; i < S[t].r; i++) {
                //zy
                // printf("here\n");

                bool connectAll = true;
                for(uint32_t j = S[z].c; j < S[z].r; j++) {
                    if(!g->connect(S[t][i], S[z][j], t)) {
                        connectAll = false;
                        break;
                    }
                }
                if(connectAll) {
                    found = false;
                    break;
                }
            }

            if(found){//found a maximal biclique: RUV and cand of z
                // thisNode=runtree.getNewNode();
                // runtree.treeNodePool[thisNode].inValue=inValue;
                // runtree.treeNodePool[thisNode].uv=side;
                // runtree.treeNodePool[thisNode].resMark=-1;

                cliqueuv[t].emplace_back(RUV[t]);
                int oldsize=RUV[z].size();
                RUV[z].reserve(RUV[z].size()+S[z].r-S[z].c);
                for(int i=S[z].c;i<S[z].r;++i){
                    RUV[z].push_back(S[z][i]);
                }
                cliqueuv[z].emplace_back(RUV[z]);
                // runtree.treeNodePool[thisNode].resMark=cliqueuv[t].size()-1;

                //reset RUV
                RUV[z].resize(oldsize);

                // if(RUV[0].size()>=p&&RUV[1].size()>=q){
                //     runtree.treeNodePool[thisNode].resMark=thisNode;
                // }
            }
        }
    }


    //TASK2, make sure C of t,z have edge(s), by shrink C of z
    uint32_t newR = S[z].r;
    uint32_t oldR=S[z].r;
    for(uint32_t i = S[z].c; i < newR; ) {
        bool degZero = true;
        for(uint32_t j = S[t].c; j < S[t].r; j++) {
            if(g->connect(S[z][i], S[t][j], z)) {
                degZero = false;
                break;
            }
        }
        if(degZero) {
            S[z].swapByPos(i, --newR);
        }
        else i++;
    }
    S[z].r = newR;

    //TASK3, find pivot and start new recursions
    if(S[t].cSize()==0||S[z].cSize()==0){
        //to connect
        // if(RUV[0].size()>=p&&RUV[1].size()>=q&&(S[0].x<S[0].c||S[1].x<S[1].c)){
        //     fans.updateSize(cliqueuv[0].size());

        //     int res=connectBiCliques();
        //     if(thisNode!=-1){
        //         if(runtree.treeNodePool[thisNode].resMark==-1){
        //             std::cerr<<"resmark error"<<std::endl;
        //             exit(0);
        //         }

        //         fans.merge(runtree.treeNodePool[thisNode].resMark,res);
        //         runtree.treeNodePool[thisNode].resMark=fans.find(res);
        //     }
        // }
        return;
    }

    if(RUV[t].size()+S[t].cSize()<threshold[t]||RUV[z].size()+S[z].cSize()<threshold[z]){
        return;
    }


    // #ifdef NOSHRINK
    // S[z].r=oldR;
    // #endif


    //find pivot
    uint32_t maxI[2] = {0}, maxE[2] = {0}, maxUV[2] = {0};
    maxI[0] = S[0].c;
    maxI[1] = S[1].c;
    maxUV[0] = S[0][S[0].c];
    maxUV[1] = S[1][S[1].c];
    t=0,z=1;
    //find pivot in C of t
    for(uint32_t i = S[t].c; i < S[t].r; i++) {//scan all vertices in C
    // for(uint32_t i = S[t].c; i <= S[t].c; i++) {
        uint32_t u = S[t][i];
        uint32_t e = 0;

        if(g->deg(u, t) > S[z].r - S[z].c) {
            for(uint32_t j = S[z].c; j < S[z].r; j++) {
                uint32_t v = S[z][j];
                if(g->connect(u, v, t)) {
                    e++;
                    deg[v]++;
                }
            }
        }
        else {
            for(uint32_t j = g->pos[t][u]; j < g->pos[t][u + 1]; j++) {
                uint32_t v = g->e[t][j];
                uint32_t pv = S[z].pos(v);
                if(S[z].c <= pv && pv < S[z].r) {//in C
                    e++;
                    deg[v]++;
                }
            }
        }

        if(e > maxE[z]) {
            maxE[z] = e;
            maxI[t] = i;
            maxUV[t] = u;
        }
    }
    //find pivot in C of z
    for(uint32_t i = S[1].c; i < S[1].r; i++) {//scan all vertices in C
    // for(uint32_t i = S[t].c; i <= S[t].c; i++) {
        uint32_t u = S[1][i];
        uint32_t e = deg[u];
        deg[u] = 0;
        if(e > maxE[0]) {
            maxE[0] = e;
            maxI[1] = i;
            maxUV[1] = u;
        }
    }

    //check if pivot in X of both sides
    bool isPivotInX[2] = {false, false};
    for(t = 0; t < 2; t++) {
        z = t ^ 1;//t + z = 1

        for(uint32_t i = S[t].x; i < S[t].c; i++) {//scan all vertices in X
        // for(uint32_t i = S[t].c; i <= S[t].c; i++) {
            uint32_t u = S[t][i];
            uint32_t e = 0;

            if(g->deg(u, t) > S[z].r - S[z].c) {
                for(uint32_t j = S[z].c; j < S[z].r; j++) {
                    uint32_t v = S[z][j];
                    if(g->connect(u, v, t)) {
                        e++;
                    }
                }
            }
            else {
                for(uint32_t j = g->pos[t][u]; j < g->pos[t][u + 1]; j++) {
                    uint32_t v = g->e[t][j];
                    uint32_t pv = S[z].pos(v);
                    if(S[z].c <= pv && pv < S[z].r) {//in C
                        e++;
                    }
                }
            }

            if(e > maxE[z]) {
                isPivotInX[t] = true;
                maxE[z] = e;
                maxI[t] = i;
                maxUV[t] = u;
            }
        }
    }

    //pivot found. pivot must have at least one neighbor in this algorithm
    t = 1;
    if(S[0].cSize() - maxE[0] >= S[1].cSize() - maxE[1]) {
        t = 0;
    }
    z = t ^ 1;

    //initialize copies of xcr
    uint32_t x[2];
    uint32_t c[2];
    uint32_t r[2];
    x[t] = S[t].x; x[z] = S[z].x;
    c[t] = S[t].c; c[z] = S[z].c;
    r[t] = S[t].r; r[z] = S[z].r;

    //r is more like the right limit of C
    //r[z] doesnt follow S[z].r
    //pivot could be in X, and C of z is changed, but never mind, it will recover by r[z]
    if(g->deg(maxUV[t], t) > S[z].r - S[z].c) {
        for(uint32_t i = S[z].r - 1; i >= S[z].c; i--) {
            if(!g->connect(maxUV[t], S[z][i], t)) {
                S[z].swapByPos(--S[z].r, i);
            }

            if(i == 0) break;
        }
    }else{
        uint32_t u = maxUV[t];
        uint32_t tmpR = S[z].c;
        for(uint32_t j = g->pos[t][u]; j < g->pos[t][u + 1]; j++) {
            uint32_t v = g->e[t][j];
            uint32_t pv = S[z].pos(v);
            if(S[z].c <= pv && pv < S[z].r) {
                S[z].swapByPos(tmpR++, pv);
            }
        }
        S[z].r = tmpR;
    }

    uint32_t wsSize = r[z] - S[z].r;
    if(ws[deep].size() < wsSize) {
        ws[deep].resize(wsSize * 2);
    }
    memcpy(ws[deep].data(), S[z].begin() + S[z].r, sizeof(uint32_t) * wsSize);

    // //debug
    // if(deep==0&&inValue==1859991){
    //     std::cerr<<"t: "<<t<<" pivot: "<<g->old_lables[t][maxUV[t]]<<std::endl;
    //     std::cerr<<wsSize<<std::endl;
    //     for(int i=0;i<ws[deep].size();++i){
    //         std::cerr<<g->old_lables[z][ws[deep][i]]<<" ";
    //     }
    //     std::cerr<<std::endl;
    // }
    // //debug
    // if(deep==0&&inValue==1859991){
    //     for(int ii=S[t].c;ii<S[t].r;++ii){
    //         std::cerr<<g->old_lables[t][S[t][ii]]<<" ";
    //     }
    //     std::cerr<<std::endl;
    //     for(int ii=S[z].c;ii<S[z].r;++ii){
    //         std::cerr<<g->old_lables[z][S[z][ii]]<<" ";
    //     }
    //     std::cerr<<std::endl;
    // }

    // int oriXstackSize=xStackNode.size();
    //iterate on pivot
    if(isPivotInX[t]==false){
        
        S[t].swapByPos(maxI[t], --r[t]);
        S[t].r = r[t];

        //x[z] doesnt follow S[z].x, used to recover S[z].x
        if(g->deg(maxUV[t], t) > S[z].c - S[z].x) {
            for(uint32_t i = S[z].x; i < S[z].c; i++) {
                if(!g->connect(maxUV[t], S[z][i], t)) {
                    S[z].swapByPos(S[z].x++, i);
                }
            }
        }else{
            uint32_t u = maxUV[t];
            uint32_t tmpX = S[z].c;
            for(uint32_t j = g->pos[t][u]; j < g->pos[t][u + 1]; j++) {
                uint32_t v = g->e[t][j];
                uint32_t pv = S[z].pos(v);
                if(S[z].x <= pv && pv < S[z].c) {
                    S[z].swapByPos(--tmpX, pv);
                }
            }
            S[z].x = tmpX;
        }

        RUV[t].push_back(maxUV[t]);
        // RUVlook[t][maxUV[t]]=1;
        
        //C of z must not be empty
        bbranchpure(deep+1,t,maxUV[t]);
        RUV[t].pop_back();
        // RUVlook[t][maxUV[t]]=0;
        // if(res!=-1){
        //     if(thisNode==-1){
        //         thisNode=runtree.getNewNode();
        //         runtree.treeNodePool[thisNode].inValue=inValue;
        //         runtree.treeNodePool[thisNode].uv=side;
        //         runtree.treeNodePool[thisNode].resMark=-1;
        //     }
        //     runtree.treeNodePool[thisNode].sons.push_back(res);
        //     xStackNode.push_back(res);
        // }

        //put pivot into X
        S[t].swapByPos(c[t]++, r[t]++);
        

    }

    // if(deep==0&&inValue==1859991){
    //     std::cerr<<S[t].cSize()<<std::endl;
    // }

    if(c[t]<r[t]){//S[t].c>0
        //iterate on vertices beyound pivot
        for(uint32_t j = 0; j < wsSize; j++){
            uint32_t w = ws[deep][j];
            S[z].swapByPos(S[z].pos(w), --r[z]);

            //recover next level parameter
            S[t].x = x[t]; S[z].x = x[z];
            S[t].c = c[t]; S[z].c = c[z];
            S[t].r = r[t]; S[z].r = r[z];

            if(g->deg(w, z) > S[t].c - S[t].x) {
                for(uint32_t i = S[t].x; i < S[t].c; i++) {
                    if(!g->connect(w, S[t][i], z)) {
                        S[t].swapByPos(S[t].x++, i);
                    }
                }
            }
            else {
                uint32_t tmpX = S[t].c;
                for(uint32_t j = g->pos[z][w]; j < g->pos[z][w + 1]; j++) {
                    uint32_t u = g->e[z][j];
                    uint32_t pu = S[t].pos(u);
                    if(S[t].x <= pu && pu < S[t].c) {
                        S[t].swapByPos(--tmpX, pu);
                    }
                }
                S[t].x = tmpX;
            }

            if(g->deg(w, z) > S[t].c - S[t].r) {
                if(S[t].r > 0){
                    for(uint32_t i = S[t].r - 1; i >= S[t].c; i--) {
                        if(!g->connect(w, S[t][i], z)) {
                            S[t].swapByPos(--S[t].r, i);
                        }

                        if(i == 0) break;
                    }
                }
            }else{
                uint32_t tmpR = S[t].c;
                for(uint32_t j = g->pos[z][w]; j < g->pos[z][w + 1]; j++) {
                    uint32_t u = g->e[z][j];
                    uint32_t pu = S[t].pos(u);
                    if(S[t].c <= pu && pu < S[t].r) {
                        S[t].swapByPos(tmpR++, pu);
                    }
                }
                S[t].r = tmpR;
            }

            // //debug
            // if(deep==0&&inValue==1859991&&w==1911392){
            //     // std::cerr<<<<std::endl;
            //     for(int ii=S[z].c;ii<S[z].r;++ii){
            //         std::cerr<<g->old_lables[z][S[z][ii]]<<" ";
            //     }
            //     std::cerr<<std::endl;
                
            // }

            if(S[t].cSize()>0){
                RUV[z].push_back(w);
                // RUVlook[z][w]=1;
                bbranchpure(deep+1,z,w);
                RUV[z].pop_back();
                // RUVlook[z][w]=0;

                // if(res!=-1){
                //     if(thisNode==-1){
                //         thisNode=runtree.getNewNode();
                //         runtree.treeNodePool[thisNode].inValue=inValue;
                //         runtree.treeNodePool[thisNode].uv=side;
                //         runtree.treeNodePool[thisNode].resMark=-1;
                //     }
                //     runtree.treeNodePool[thisNode].sons.push_back(res);
                //     xStackNode.push_back(res);
                // }
            }

            S[z].swapByPos(c[z]++,r[z]++);//move w to X
        }

        for(uint32_t j = 0; j < wsSize; j++) {
            uint32_t w = ws[deep][j];
            S[z].swapByPos(S[z].pos(w), --c[z]);
        }
    }

    //recover pivot
    if(isPivotInX[t] == false) {
        S[t].swapByPos(S[t].pos(maxUV[t]), --c[t]);
    }

    //recover x,c,r
    S[t].x = x[t]; S[z].x = x[z];
    S[t].c = c[t]; S[z].c = c[z];
    S[t].r = r[t]; S[z].r = r[z];

    //connect that should be connected
    // if(RUV[0].size()>=p&&RUV[1].size()>=q){
    //     fans.updateSize(cliqueuv[0].size());

    //     int res=connectBiCliques();
    //     if(thisNode!=-1){
    //         if(runtree.treeNodePool[thisNode].resMark==-1){
    //             runtree.treeNodePool[thisNode].resMark=res;
    //         }else{
    //             fans.merge(runtree.treeNodePool[thisNode].resMark,res);
    //             runtree.treeNodePool[thisNode].resMark=fans.find(res);
    //         }
    //     }
    // }

    // xStackNode.resize(oriXstackSize);
    return;

    
}

void BICPC::reportTreeNode(int thisNodeId,int fid,int pivot,int pivotside){
    std::cout<<"idx "<<thisNodeId<<" fidx "<<fid<<std::endl;
    for(int i=0;i<RUV[0].size();++i){
        std::cout<<g->old_lables[0][RUV[0][i]]<<" ";
    }
    std::cout<<"; ";
    for(int i=S[0].c;i<S[0].r;++i){
        std::cout<<g->old_lables[0][S[0][i]]<<" ";
    }
    std::cout<<"; ";
    for(int i=S[0].x;i<S[0].c;++i){
        std::cout<<g->old_lables[0][S[0][i]]<<" ";
    }
    std::cout<<",";
    for(int i=0;i<RUV[1].size();++i){
        std::cout<<g->old_lables[1][RUV[1][i]]<<" ";
    }
    std::cout<<"; ";
    for(int i=S[1].c;i<S[1].r;++i){
        std::cout<<g->old_lables[1][S[1][i]]<<" ";
    }
    std::cout<<"; ";
    for(int i=S[1].x;i<S[1].c;++i){
        std::cout<<g->old_lables[1][S[1][i]]<<" ";
    }

    if(pivot>=0){
        std::cout<<",";
        char sidemark='U';
        if(pivotside==1) sidemark='V';
        std::cout<<sidemark<<" "<<g->old_lables[pivotside][pivot];
    }

    std::cout<<std::endl;
}

void BICPC::listMBCliquesTree(){
    auto mt1 = std::chrono::steady_clock::now();

    mbeTreeNodeCnt=0;
    int thisNodeIdx=mbeTreeNodeCnt++;
    uint32_t t=0, z=1;
    if(g->maxDu<g->maxDv) {t=1,z=0;}
    // firstlayert=t;

    // std::cerr<<"t: "<<t<<" "<<"z: "<<z<<std::endl;
    // std::cerr<<"missing clique:"<<std::endl;
    // std::cerr<<g->labelsL[1]<<" "<<g->labelsL[184683]<<std::endl;
    // std::cerr<<g->labelsR[133490]<<std::endl;

    // runtree.root=runtree.getNewNode();
    // runtree.treeNodePool[runtree.root].inValue=-1;
    // runtree.treeNodePool[runtree.root].resMark=-1;
    // runtree.treeNodePool[runtree.root].uv=-1;

    // runtree.firstlayerson.resize(g->n[t],-1);
    // xStackNode.reserve(g->n[0]+g->n[1]);

    //first layer no pivot
    for(uint32_t u=0;u<g->n[t];u++){
        S[t].c=S[t].r=g->n[t];
        S[z].c=S[z].r=g->n[z];


        // //debug
        // if(u==25341){
        //     std::cerr<<"here"<<std::endl;
        // }
        

        //form C set for z
        for(uint32_t i=g->pos[t][u];i<g->pos[t][u+1];++i){
            uint32_t v=g->e[t][i];
            S[z].swapByPos(--S[z].c,S[z].pos(v));
        }

        S[z].x=S[z].c;

        //form C set for t
        for(uint32_t i = g->pos[t][u]; i < g->pos[t][u + 1]; i++) {
            uint32_t v = g->e[t][i];
            if(g->pos[z][v + 1] > 0){
                for(uint32_t j = g->pos[z][v + 1] - 1; j >= g->pos[z][v]; j--) {
                    uint32_t w = g->e[z][j];
                    
                    if(w > u) {
                        uint32_t pw = S[t].pos(w);
                        if(pw < S[t].c) {
                            S[t].swapByPos(--S[t].c, pw);
                        }
                    }
                    else break;

                    if(j == 0) break;
                }
            }
        }
        S[t].x = S[t].c;

        //形成t的X集合
        for(uint32_t i = g->pos[t][u]; i < g->pos[t][u + 1]; i++) {
            uint32_t v = g->e[t][i];
            if(g->pos[z][v + 1] > 0){
                for(uint32_t j = g->pos[z][v]; j < g->pos[z][v + 1]; j++) {
                    uint32_t w = g->e[z][j];
                    
                    if(w < u) {
                        uint32_t pw = S[t].pos(w);
                        if(pw < S[t].x) {
                            S[t].swapByPos(--S[t].x, pw);
                        }
                    }
                    else break;

                    if(j == 0) break;
                }
            }
        }

        // if(u==1859991){
        //     for(int i=S[t].c;i<S[t].r;++i){
        //         std::cerr<<g->old_lables[0][S[t][i]]<<" ";
        //     }
        //     std::cerr<<std::endl;
        //     for(int i=S[z].c;i<S[z].r;++i){
        //         std::cerr<<g->old_lables[1][S[z][i]]<<" ";
        //     }
        //     std::cerr<<std::endl;
        // }


        if(S[z].CIsEmpty()) continue;

        RUV[t].push_back(u);
        // RUVlook[t][u]=1;
        bbranchpureFullTree(0,t,u,thisNodeIdx);
        RUV[t].pop_back();
        // RUVlook[t][u]=0;

        // if(res!=-1){
        //     runtree.firstlayerson[u]=res;
        //     runtree.treeNodePool[runtree.root].sons.push_back(res);
        // }

    }

    std::cerr<<"maxBiCliqueCount: "<<cliqueuv[0].size()<<std::endl;
    auto mt2 = std::chrono::steady_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(mt2 - mt1);
    std::cerr << "list mbclique time: " << duration.count() << " ms" << std::endl;
}

void BICPC::listMBCliquesTreeFullPivot(){
    auto mt1 = std::chrono::steady_clock::now();

    mbeTreeNodeCnt=0;
    int thisNodeIdx=mbeTreeNodeCnt++;
    uint32_t t=0, z=1;
    if(g->maxDu<g->maxDv) {t=1,z=0;}
    // firstlayert=t;

    // std::cerr<<"t: "<<t<<" "<<"z: "<<z<<std::endl;
    // std::cerr<<"missing clique:"<<std::endl;
    // std::cerr<<g->labelsL[1]<<" "<<g->labelsL[184683]<<std::endl;
    // std::cerr<<g->labelsR[133490]<<std::endl;

    // runtree.root=runtree.getNewNode();
    // runtree.treeNodePool[runtree.root].inValue=-1;
    // runtree.treeNodePool[runtree.root].resMark=-1;
    // runtree.treeNodePool[runtree.root].uv=-1;

    // runtree.firstlayerson.resize(g->n[t],-1);
    // xStackNode.reserve(g->n[0]+g->n[1]);

    S[0].r=g->n[0];S[0].c=S[0].x=0;
    S[1].r=g->n[1];S[1].c=S[1].x=0;

    bbranchpureFullTree(0,0,0,thisNodeIdx);

    // //first layer no pivot
    // for(uint32_t u=0;u<g->n[t];u++){
    //     S[t].c=S[t].r=g->n[t];
    //     S[z].c=S[z].r=g->n[z];


    //     // //debug
    //     // if(u==25341){
    //     //     std::cerr<<"here"<<std::endl;
    //     // }
        

    //     //form C set for z
    //     for(uint32_t i=g->pos[t][u];i<g->pos[t][u+1];++i){
    //         uint32_t v=g->e[t][i];
    //         S[z].swapByPos(--S[z].c,S[z].pos(v));
    //     }

    //     S[z].x=S[z].c;

    //     //form C set for t
    //     for(uint32_t i = g->pos[t][u]; i < g->pos[t][u + 1]; i++) {
    //         uint32_t v = g->e[t][i];
    //         if(g->pos[z][v + 1] > 0){
    //             for(uint32_t j = g->pos[z][v + 1] - 1; j >= g->pos[z][v]; j--) {
    //                 uint32_t w = g->e[z][j];
                    
    //                 if(w > u) {
    //                     uint32_t pw = S[t].pos(w);
    //                     if(pw < S[t].c) {
    //                         S[t].swapByPos(--S[t].c, pw);
    //                     }
    //                 }
    //                 else break;

    //                 if(j == 0) break;
    //             }
    //         }
    //     }
    //     S[t].x = S[t].c;

    //     //形成t的X集合
    //     for(uint32_t i = g->pos[t][u]; i < g->pos[t][u + 1]; i++) {
    //         uint32_t v = g->e[t][i];
    //         if(g->pos[z][v + 1] > 0){
    //             for(uint32_t j = g->pos[z][v]; j < g->pos[z][v + 1]; j++) {
    //                 uint32_t w = g->e[z][j];
                    
    //                 if(w < u) {
    //                     uint32_t pw = S[t].pos(w);
    //                     if(pw < S[t].x) {
    //                         S[t].swapByPos(--S[t].x, pw);
    //                     }
    //                 }
    //                 else break;

    //                 if(j == 0) break;
    //             }
    //         }
    //     }

    //     // if(u==1859991){
    //     //     for(int i=S[t].c;i<S[t].r;++i){
    //     //         std::cerr<<g->old_lables[0][S[t][i]]<<" ";
    //     //     }
    //     //     std::cerr<<std::endl;
    //     //     for(int i=S[z].c;i<S[z].r;++i){
    //     //         std::cerr<<g->old_lables[1][S[z][i]]<<" ";
    //     //     }
    //     //     std::cerr<<std::endl;
    //     // }


    //     if(S[z].CIsEmpty()) continue;

    //     RUV[t].push_back(u);
    //     // RUVlook[t][u]=1;
    //     bbranchpureFullTree(0,t,u,thisNodeIdx);
    //     RUV[t].pop_back();
    //     // RUVlook[t][u]=0;

    //     // if(res!=-1){
    //     //     runtree.firstlayerson[u]=res;
    //     //     runtree.treeNodePool[runtree.root].sons.push_back(res);
    //     // }

    // }

    std::cerr<<"maxBiCliqueCount: "<<cliqueuv[0].size()<<std::endl;
    auto mt2 = std::chrono::steady_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(mt2 - mt1);
    std::cerr << "list mbclique time: " << duration.count() << " ms" << std::endl;
}

void BICPC::listMBCliques(){
    auto mt1 = std::chrono::steady_clock::now();

    uint32_t t=0, z=1;
    if(g->maxDu<g->maxDv) {t=1,z=0;}
    // firstlayert=t;

    // std::cerr<<"t: "<<t<<" "<<"z: "<<z<<std::endl;
    // std::cerr<<"missing clique:"<<std::endl;
    // std::cerr<<g->labelsL[1]<<" "<<g->labelsL[184683]<<std::endl;
    // std::cerr<<g->labelsR[133490]<<std::endl;

    // runtree.root=runtree.getNewNode();
    // runtree.treeNodePool[runtree.root].inValue=-1;
    // runtree.treeNodePool[runtree.root].resMark=-1;
    // runtree.treeNodePool[runtree.root].uv=-1;

    // runtree.firstlayerson.resize(g->n[t],-1);
    // xStackNode.reserve(g->n[0]+g->n[1]);

    //first layer no pivot
    for(uint32_t u=0;u<g->n[t];u++){
        S[t].c=S[t].r=g->n[t];
        S[z].c=S[z].r=g->n[z];


        // //debug
        // if(u==25341){
        //     std::cerr<<"here"<<std::endl;
        // }
        

        //form C set for z
        for(uint32_t i=g->pos[t][u];i<g->pos[t][u+1];++i){
            uint32_t v=g->e[t][i];
            S[z].swapByPos(--S[z].c,S[z].pos(v));
        }

        S[z].x=S[z].c;

        //form C set for t
        for(uint32_t i = g->pos[t][u]; i < g->pos[t][u + 1]; i++) {
            uint32_t v = g->e[t][i];
            if(g->pos[z][v + 1] > 0){
                for(uint32_t j = g->pos[z][v + 1] - 1; j >= g->pos[z][v]; j--) {
                    uint32_t w = g->e[z][j];
                    
                    if(w > u) {
                        uint32_t pw = S[t].pos(w);
                        if(pw < S[t].c) {
                            S[t].swapByPos(--S[t].c, pw);
                        }
                    }
                    else break;

                    if(j == 0) break;
                }
            }
        }
        S[t].x = S[t].c;

        //形成t的X集合
        for(uint32_t i = g->pos[t][u]; i < g->pos[t][u + 1]; i++) {
            uint32_t v = g->e[t][i];
            if(g->pos[z][v + 1] > 0){
                for(uint32_t j = g->pos[z][v]; j < g->pos[z][v + 1]; j++) {
                    uint32_t w = g->e[z][j];
                    
                    if(w < u) {
                        uint32_t pw = S[t].pos(w);
                        if(pw < S[t].x) {
                            S[t].swapByPos(--S[t].x, pw);
                        }
                    }
                    else break;

                    if(j == 0) break;
                }
            }
        }

        // if(u==1859991){
        //     for(int i=S[t].c;i<S[t].r;++i){
        //         std::cerr<<g->old_lables[0][S[t][i]]<<" ";
        //     }
        //     std::cerr<<std::endl;
        //     for(int i=S[z].c;i<S[z].r;++i){
        //         std::cerr<<g->old_lables[1][S[z][i]]<<" ";
        //     }
        //     std::cerr<<std::endl;
        // }


        if(S[z].CIsEmpty()) continue;

        RUV[t].push_back(u);
        // RUVlook[t][u]=1;
        bbranchpure(0,t,u);
        RUV[t].pop_back();
        // RUVlook[t][u]=0;

        // if(res!=-1){
        //     runtree.firstlayerson[u]=res;
        //     runtree.treeNodePool[runtree.root].sons.push_back(res);
        // }

    }

    std::cerr<<"maxBiCliqueCount: "<<cliqueuv[0].size()<<std::endl;
    auto mt2 = std::chrono::steady_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(mt2 - mt1);
    std::cerr << "list mbclique time: " << duration.count() << " ms" << std::endl;
}

#endif