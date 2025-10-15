#ifndef UNIONFINDCPP
#define UNIONFINDCPP
#include<vector>
using std::vector;

struct unf{
    vector<int> fa;

    inline void init(int n){
        fa.resize(n);
        for(int i=0;i<fa.size();++i) fa[i]=i;
    }

    inline void reinit(int l,int r){
        for(int i=l;i<=r;++i){
            fa[i]=i;
        }
    }

    inline void updateSize(int n){
        int oriSize=fa.size();
        if(oriSize>=n) return;
        fa.resize(n);
        for(int i=oriSize;i<n;++i) fa[i]=i;
    }

    inline int find(int x){
        if(x<0||x>=fa.size()){
            std::cerr<<"mem error"<<std::endl;
        }
        if(x==fa[x]) return x;
        else return (fa[x]=find(fa[x]));
    }
    
    inline void merge(int x,int y){
        fa[find(x)]=find(y);
    }
};

#endif