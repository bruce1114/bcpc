#ifndef BIGRAPH_CPP
#define BIGRAPH_CPP

#include <vector>
#include <algorithm>
#include <iostream>
#include<fstream>

#include "../tools/fastIO.hpp"
#include "../tools/listLinearHeap.hpp"
#include "../tools/hopstotchHash.hpp"


struct bigraph{
    uint32_t n1,n2,m,maxDu,maxDv;
    uint32_t core[2];
    uint32_t n[2];

    struct Edge{
        uint32_t u,v;
    };
    std::vector<Edge> edges;

    std::vector<uint32_t> posU,e1,posV,e2;
    std::vector<uint32_t> pos[2];
    std::vector<uint32_t> e[2];
    std::vector<uint32_t> cores[2];

    std::vector<CuckooHash> cuhash[2];

    std::vector<uint32_t> old_lables[2];
    std::vector<uint32_t> labelsL;
    std::vector<uint32_t> labelsR;

    // std::vector<std::vector<uint32_t>> two_hop_neighbor;
    std::vector<uint32_t> two_hop_pos,two_hop_e;

    void prepareTowHopNeighbor(int anchor,int qpthreshold);

    bigraph(const std::string& filename){
        read(filename);

        // //debug
        // for(int i=0;i<old_lables[0].size();++i){
        //     std::cout<<old_lables[0][i]<<std::endl;
        // }
        // for(int i=0;i<old_lables[1].size();++i){
        //     std::cout<<old_lables[1][i]<<std::endl;
        // }
        // exit(0);

        // cuhash[0].resize(n[0])
        // cuhash[1].resize(n[1]);
        // for(uint32_t t = 0; t <= 1; t++) {
        //     for(uint32_t i = 0; i < n[t]; i++) {
        //         uint32_t d = pos[t][i + 1] - pos[t][i];
        //         cuhash[t][i].reserve(d+1);
        //         for (uint32_t j = pos[t][i]; j < pos[t][i + 1];++j)
        //             cuhash[t][i].insert(e[t][j]);
        //     }
        // }
    }

    bigraph(const std::string& filename,uint32_t scala){
        read(filename,scala);

        // //debug
        // for(int i=0;i<old_lables[0].size();++i){
        //     std::cout<<old_lables[0][i]<<std::endl;
        // }
        // for(int i=0;i<old_lables[1].size();++i){
        //     std::cout<<old_lables[1][i]<<std::endl;
        // }
        // exit(0);

        // cuhash[0].resize(n[0])
        // cuhash[1].resize(n[1]);
        // for(uint32_t t = 0; t <= 1; t++) {
        //     for(uint32_t i = 0; i < n[t]; i++) {
        //         uint32_t d = pos[t][i + 1] - pos[t][i];
        //         cuhash[t][i].reserve(d+1);
        //         for (uint32_t j = pos[t][i]; j < pos[t][i + 1];++j)
        //             cuhash[t][i].insert(e[t][j]);
        //     }
        // }
    }

    void read(const std::string& filename,uint32_t scala){
        fastIO in(filename,"r");

        n1 = in.getUInt();
        n2 = in.getUInt();
        m = in.getUInt();

        m=m*scala/100;

        edges.resize(m);
        e1.resize(m);
        e2.resize(m);
        posU.resize(n1 + 5);
        posV.resize(n2 + 5);

        for(uint32_t i = 0; i < m; i++) {
            edges[i].u = in.getUInt();
            edges[i].v = in.getUInt();
        }

        cores[0].resize(n1);
        cores[1].resize(n2);
        old_lables[0].resize(n1);
        old_lables[1].resize(n2);
        changeToTwoHopCoreOrder();

        n[0] = n1;
        n[1] = n2;
        pos[0] = std::move(posU);
        pos[1] = std::move(posV);
        e[0] = std::move(e1);
        e[1] = std::move(e2);
    }

    void read(const std::string& filename){
        fastIO in(filename,"r");

        n1 = in.getUInt();
        n2 = in.getUInt();
        m = in.getUInt();

        edges.resize(m);
        e1.resize(m);
        e2.resize(m);
        posU.resize(n1 + 5);
        posV.resize(n2 + 5);

        for(uint32_t i = 0; i < m; i++) {
            edges[i].u = in.getUInt();
            edges[i].v = in.getUInt();
        }

        cores[0].resize(n1);
        cores[1].resize(n2);
        old_lables[0].resize(n1);
        old_lables[1].resize(n2);
        changeToTwoHopCoreOrder();

        n[0] = n1;
        n[1] = n2;
        pos[0] = std::move(posU);
        pos[1] = std::move(posV);
        e[0] = std::move(e1);
        e[1] = std::move(e2);
    }

    void changeToTwoHopCoreOrder() {
        // printf("here two\n");fflush(stdout);
        std::vector<uint32_t> d1, d2;
        
        d1.resize(n1,0);
        d2.resize(n2,0);

        for(uint32_t i = 0; i < m; i++) {
            ++d1[edges[i].u];
            ++d2[edges[i].v];
        }

        maxDu = 0;
        for(uint32_t i = 0; i < n1; i++) {
            maxDu = std::max(maxDu, d1[i]);
        }
        maxDv = 0;
        for(uint32_t i = 0; i < n2; i++) {
            maxDv = std::max(maxDv, d2[i]);
        }

        posU[0] = 0;
        for(uint32_t u = 0; u < n1; u++) {
            posU[u + 1] = d1[u] + posU[u];
        }
        for(uint32_t i = 0; i < m; i++) {
            e1[posU[edges[i].u]++] = edges[i].v; 
        }
        posU[0] = 0;
        for(uint32_t u = 0; u < n1; u++) {
            posU[u + 1] = d1[u] + posU[u];
        } 
        
        posV[0] = 0;
        for(uint32_t v = 0; v < n2; v++) {
            posV[v + 1] = d2[v] + posV[v];
        }
        for(uint32_t i = 0; i < m; i++) {
            e2[posV[edges[i].v]++] = edges[i].u; 
        }
        posV[0] = 0;
        for(uint32_t v = 0; v < n2; v++) {
            posV[v + 1] = d2[v] + posV[v];
        }

        uint32_t N = n1+n2;
        std::vector<uint32_t> dequeue;
        std::vector<uint32_t> cdeg;
        dequeue.resize(N);
        cdeg.resize(N, 0);
        uint32_t tsize = 0;
        
        n[0] = n1; n[1] = n2;

        //将度为0的节点添加到dequeue中
        for(uint32_t tt = 0; tt <= 1; tt++) {
            for(uint32_t i = 0; i < n[tt]; i++) {
                uint32_t d = 0;
                if (tt == 0) {
                    d = d1[i];
                    cdeg[i] = d;
                    if (d == 0) dequeue[tsize++] = i;
                }
                else  {
                    d = d2[i];
                    cdeg[i+n1] = d;
                    if (d == 0) dequeue[tsize++] = i+n1;
                }
            }
        }
     
        //连坐
        uint32_t rsize = 0;
        while (tsize > rsize) {
            uint32_t s = rsize;
            rsize = tsize;
            for (uint32_t i = s; i < rsize; ++i) {
                int v = dequeue[i];
                if (v >= n1) {
                    v -= n1;
                    for (uint32_t j = posV[v]; j < posV[v+1]; ++j) {
                        uint32_t u = e2[j];
                        if (cdeg[u] > 0) {
                            cdeg[u] --;
                            if (cdeg[u] == 0) dequeue[tsize++] = u;
                        }
                    }
                }
                else {
                    for (uint32_t j = posU[v]; j < posU[v+1]; ++j) {
                        uint32_t u = e1[j];
                        u += n1;
                        if (cdeg[u] > 0) {
                            cdeg[u] --;
                            if (cdeg[u] == 0) dequeue[tsize++] = u;
                        }
                    }
                }
            }
        //     break;
        }

        uint32_t maxTwoHopDegreeU = 0, maxTwoHopDegreeV = 0;
        // d1.clear();
        d1.resize(d1.size(),0);//zy
        d2.resize(d2.size(),0);//zy
        std::vector<uint32_t> stk;
        uint32_t n = std::max(n1, n2);
        std::vector<uint32_t> ids(n);
        std::vector<uint32_t> keys(n);
        labelsL.resize(n1);
        labelsR.resize(n2);
        uint32_t l = 0;
        
        for(uint32_t u = 0; u < n1; u++) ids[u] = u;

        //统计U中节点的1、2跳邻居总数
        for(uint32_t u = 0; u < n1; u++) {
            if (cdeg[u] == 0) {
                keys[u] = 0;
                continue;
            }
            stk.clear();
            uint32_t tmp = 0;

            for(uint32_t i = posU[u]; i < posU[u + 1]; i++) {
                uint32_t v = e1[i];
                if (cdeg[v+n1] == 0) continue;
                tmp++;
                for(uint32_t j = posV[v]; j < posV[v + 1]; j++) {
                    uint32_t w = e2[j];
                    if(d1[w] == 0 && cdeg[w] > 0) {
                        d1[w] = 1;
                        stk.push_back(w);
                    }
                }
            }
            tmp += stk.size();

            maxTwoHopDegreeU = std::max(maxTwoHopDegreeU, tmp);
            keys[u] = tmp;

            for(auto w : stk) d1[w] = 0;
        }

        ListLinearHeap heap(n1, maxTwoHopDegreeU + 1);
        heap.init(n1, maxTwoHopDegreeU + 1, ids.data(), keys.data());
        core[0] = 0;
        core[1] = 0;

        
        for(uint32_t i = 0; i < n1; i++) {
            uint32_t u, degU;

            if(!heap.pop_min(u, degU)) printf("errorLheap\n");
            old_lables[0][l] = u;
            labelsL[u] = l++;//record the place where u is (in core number sequece)
            core[0] = std::max(core[0], degU);
            cores[0][l-1] = core[0];//record the core number order
            if (core[0] == 0) continue;
            stk.clear();

            //对u的所有二跳邻居减1
            for(uint32_t i = posU[u]; i < posU[u + 1]; i++) {
                uint32_t v = e1[i];
                if (cdeg[v+n1] == 0) continue;
                for(uint32_t j = posV[v]; j < posV[v + 1]; j++) {
                    uint32_t w = e2[j];
                    if(d1[w] == 0 && cdeg[w] > 0) {
                        d1[w] = 1;
                        stk.push_back(w);
                        heap.decrement(w, 1);
                    }
                }
            }

            for(auto w : stk) d1[w] = 0;
        }

        d2.clear();
        l = 0;
        
        //统计V中节点的1、2跳邻居总数
        for(uint32_t u = 0; u < n2; u++) ids[u] = u;
        for(uint32_t u = 0; u < n2; u++) {
            if (cdeg[u+n1] == 0) {
                keys[u] = 0;
                continue;
            }
            stk.clear();
            uint32_t tmp = 0;

            for(uint32_t i = posV[u]; i < posV[u + 1]; i++) {
                uint32_t v = e2[i];
                if (cdeg[v] == 0) continue;
                tmp++;
                for(uint32_t j = posU[v]; j < posU[v + 1]; j++) {
                    uint32_t w = e1[j];
                    if(d2[w] == 0 && cdeg[w+n1] > 0) {
                        d2[w] = 1;
                        stk.push_back(w);
                    }
                }
            }
            tmp += stk.size();

            maxTwoHopDegreeV = std::max(maxTwoHopDegreeV, tmp);
            keys[u] = tmp;

            for(auto w : stk) d2[w] = 0;
        }

        // printf("core %u %u\n", maxTwoHopDegreeU, maxTwoHopDegreeV);

        ListLinearHeap rheap(n2, maxTwoHopDegreeV + 1);
        rheap.init(n2, maxTwoHopDegreeV + 1, ids.data(), keys.data());

        for(uint32_t i = 0; i < n2; i++) {
            uint32_t u, degU;

            if(!rheap.pop_min(u, degU)) printf("errorLheap\n");
            old_lables[1][l] = u;
            labelsR[u] = l++;
            core[1] = std::max(core[1], degU);
            cores[1][l-1] = core[1];
            if (core[1] == 0) continue;
            stk.clear();
            for(uint32_t i = posV[u]; i < posV[u + 1]; i++) {
                uint32_t v = e2[i];
                if (cdeg[v] == 0) continue;
                for(uint32_t j = posU[v]; j < posU[v + 1]; j++) {
                    uint32_t w = e1[j];
                    if(d2[w] == 0 && cdeg[w+n1] > 0) {
                        d2[w] = 1;
                        stk.push_back(w);
                        heap.decrement(w, 1);
                    }
                }
            }

            for(auto w : stk) d2[w] = 0;
        }
        
        //??
        for (uint32_t i = 0; i < n1; ++i) {
            int v = labelsL[i];
            cores[0][v] = cdeg[i];
        }
        for (uint32_t i = 0; i < n2; ++i) {
            int v = labelsR[i];
            cores[1][v] = cdeg[i+n1];
        }

        for(uint32_t i = 0; i < m; i++) {
            edges[i].u = labelsL[edges[i].u];
            edges[i].v = labelsR[edges[i].v];
        }

        std::fill(d1.begin(), d1.begin() + n1, 0);
        std::fill(d2.begin(), d2.begin() + n2, 0);
        std::fill(posU.begin(), posU.begin() + n1 + 1, 0);
        std::fill(posV.begin(), posV.begin() + n2 + 1, 0);

        for(uint32_t i = 0; i < m; i++) {
            ++d1[edges[i].u];
            ++d2[edges[i].v];
        }

        for(uint32_t i = 0; i < n1; i++) {
            posU[i + 1] = posU[i] + d1[i];
        }
        for(uint32_t i = 0; i < n2; i++) {
            posV[i + 1] = posV[i] + d2[i];
        }

        for(uint32_t i = 0; i < m; i++) {
            e1[ posU[edges[i].u]++ ] = edges[i].v;
        }
        for(uint32_t i = 0; i < m; i++) {
            e2[ posV[edges[i].v]++ ] = edges[i].u;
        }

        posU[0] = posV[0] = 0;
        for(uint32_t i = 0; i < n1; i++) {
            posU[i + 1] = posU[i] + d1[i];
        }
        for(uint32_t i = 0; i < n2; i++) {
            posV[i + 1] = posV[i] + d2[i];
        }

        for(uint32_t i = 0; i < n1; i++) {
            std::sort(e1.begin() + posU[i], e1.begin() + posU[i + 1]);
        }
        for(uint32_t i = 0; i < n2; i++) {
            std::sort(e2.begin() + posV[i], e2.begin() + posV[i + 1]);
        }
    }

    bool connect(uint32_t u,uint32_t v,uint32_t t){
        return cuhash[t][u].find(v);
    }

    uint32_t deg(uint32_t u, uint32_t t) {
        return pos[t][u + 1] - pos[t][u];
    }

    void qpcoreReduction(int p,int q);

    void getcoreSize(int& usize,int& vsize);

    void initializeHash();

    void writeOldLabels();
};

void bigraph::writeOldLabels(){
    std::ofstream fout("old_label.txt",std::ios::out);
    int maxn=std::max(n[0],n[1]);

    for(int i=0;i<maxn;++i){
        
        int leftres,rightres;

        if(i<n[0]) leftres=old_lables[0][i];
        else leftres=-1;
        if(i<n[1]) rightres=old_lables[1][i];
        else rightres=-1;

        fout<<i<<" "<<leftres<<" "<<rightres<<std::endl;
    }
}

void bigraph::initializeHash(){
    cuhash[0].resize(n[0]);
    cuhash[1].resize(n[1]);
    for(uint32_t t = 0; t <= 1; t++) {
        for(uint32_t i = 0; i < n[t]; i++) {
            uint32_t d = pos[t][i + 1] - pos[t][i];
            cuhash[t][i].reserve(d+1);
            for (uint32_t j = pos[t][i]; j < pos[t][i + 1];++j)
                cuhash[t][i].insert(e[t][j]);
        }
    }
}

void bigraph::getcoreSize(int& usize,int& vsize){
    usize=0;
    vsize=0;
    for(int u=0;u<n[0];++u){
        if(pos[0][u]<pos[0][u+1]) usize++;
    }
    for(int v=0;v<n[1];++v){
        if(pos[1][v]<pos[1][v+1]) vsize++;
    }
}

void bigraph::qpcoreReduction(int p,int q){
    std::vector<uint32_t> d1, d2;
    std::vector<uint32_t> ids1,ids2;
    
    d1.resize(n1,0);
    d2.resize(n2,0);

    ids1.resize(n1);
    ids2.resize(n2);

    for(int i=0;i<n1;++i) ids1[i]=i;
    for(int i=0;i<n2;++i) ids2[i]=i;

    for(uint32_t i = 0; i < m; i++) {
        ++d1[edges[i].u];
        ++d2[edges[i].v];
    }

    uint32_t maxd1=0;
    uint32_t maxd2=0;
    for(int i=0;i<n1;++i){
        maxd1=std::max(maxd1,d1[i]);
    }
    for(int i=0;i<n2;++i){
        maxd2=std::max(maxd2,d2[i]);
    }

    ListLinearHeap heap1(n1,maxd1+1);
    ListLinearHeap heap2(n2,maxd2+1);
    heap1.init(n1,maxd1+1,ids1.data(),d1.data());
    heap2.init(n2,maxd2+1,ids2.data(),d2.data());

    while(true){
        int upop_num=0;
        int vpop_num=0;

        //pop U
        while(true){
            uint32_t u,du;
            if(heap1.get_min(u,du)==false){
                break;
            }else if(du>=q){
                break;
            }

            heap1.pop_min(u,du);
            upop_num++;
            d1[u]=0;

            for(int i=pos[0][u];i<pos[0][u+1];++i){
                int v=e[0][i];
                if(d2[v]==0) continue;
                heap2.decrement(v,1);
            }
        }

        //pop V
        while(true){
            uint32_t v,dv;
            if(heap2.get_min(v,dv)==false){
                break;
            }else if(dv>=p){
                break;
            }

            heap2.pop_min(v,dv);
            vpop_num++;
            d2[v]=0;

            for(int i=pos[1][v];i<pos[1][v+1];++i){
                int u=e[1][i];
                if(d1[u]==0) continue;
                heap1.decrement(u,1);
            }
        }

        if(upop_num==0||vpop_num) break;
    }

    std::vector<uint32_t> d1c=d1;
    std::vector<uint32_t> d2c=d2;

    posU=std::vector<uint32_t>(n1+5);
    posV=std::vector<uint32_t>(n2+5);
    std::fill(d1.begin(), d1.begin() + n1, 0);
    std::fill(d2.begin(), d2.begin() + n2, 0);
    std::fill(posU.begin(), posU.begin() + n1 + 1, 0);
    std::fill(posV.begin(), posV.begin() + n2 + 1, 0);


    std::vector<Edge> nedges;
    nedges.reserve(edges.size());
    for(int i=0;i<m;++i){
        uint32_t u=edges[i].u;
        uint32_t v=edges[i].v;
        if(d1c[u]==0||d2c[v]==0) continue;

        nedges.push_back(edges[i]);
        ++d1[edges[i].u];
        ++d2[edges[i].v];
    }
    edges=nedges;
    m=edges.size();

    e1=std::vector<uint32_t>(m);
    e2=std::vector<uint32_t>(m);


    for(uint32_t i = 0; i < n1; i++) {
        posU[i + 1] = posU[i] + d1[i];
    }
    for(uint32_t i = 0; i < n2; i++) {
        posV[i + 1] = posV[i] + d2[i];
    }

    for(uint32_t i = 0; i < m; i++) {
        e1[ posU[edges[i].u]++ ] = edges[i].v;
    }
    for(uint32_t i = 0; i < m; i++) {
        e2[ posV[edges[i].v]++ ] = edges[i].u;
    }

    posU[0] = posV[0] = 0;
    for(uint32_t i = 0; i < n1; i++) {
        posU[i + 1] = posU[i] + d1[i];
    }
    for(uint32_t i = 0; i < n2; i++) {
        posV[i + 1] = posV[i] + d2[i];
    }

    for(uint32_t i = 0; i < n1; i++) {
        std::sort(e1.begin() + posU[i], e1.begin() + posU[i + 1]);
    }
    for(uint32_t i = 0; i < n2; i++) {
        std::sort(e2.begin() + posV[i], e2.begin() + posV[i + 1]);
    }

    pos[0] = std::move(posU);
    pos[1] = std::move(posV);
    e[0] = std::move(e1);
    e[1] = std::move(e2);

}

void bigraph::prepareTowHopNeighbor(int anchor,int qpthreshold){
    int nanchor=n[anchor];
    int t=anchor;
    int z=anchor^1;

    std::vector<uint32_t> common_neig_map,aux_array_two_neig;
    common_neig_map.resize(nanchor,0);
    aux_array_two_neig.resize(nanchor);

    two_hop_pos.resize(nanchor+1,0);
    std::vector<std::pair<int,int>> edges;
    std::vector<int> two_hop_deg;
    two_hop_deg.resize(nanchor+1,0);

    for(int u=0;u<nanchor;++u){
        int idx=0;
        for(int j=pos[t][u];j<pos[t][u+1];++j){
            int v=e[t][j];
            for(int k=pos[z][v];k<pos[z][v+1];++k){
                int w=e[z][k];
                if(w<u){
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
            if(common_neig_map[w]>=qpthreshold){
                // two_hop_neighbor[w].push_back(u);//direct edge
                edges.push_back(std::pair<int,int>(w,u));
            }
            common_neig_map[w]=0;
        }
    }

    for(int i=0;i<edges.size();++i){
        int u=edges[i].first;
        int v=edges[i].second;

        two_hop_pos[u]++;
    }

    for(int u=1;u<nanchor;++u){
        two_hop_pos[u]+=two_hop_pos[u-1];
    }

    two_hop_e.resize(two_hop_pos[nanchor-1]);
    two_hop_pos[nanchor]=two_hop_pos[nanchor-1];

    for(int i=0;i<edges.size();++i){
        int u=edges[i].first;
        int v=edges[i].second;

        two_hop_e[--two_hop_pos[u]]=v;
    }

    //no need to sort?
    // for(int u=0;u<nanchor;++u){
    //     sort(two_hop_e.begin()+two_hop_pos[u],two_hop_e.begin()+two_hop_pos[u+1]);
    // }


    
}


#endif