#ifndef UTILSCPP
#define UTILSCPP


#include<vector>

inline void exchange(std::vector<int>& arr,int aIndex,int bIndex,std::vector<int>& arrIndex){
    int va=arr[aIndex];
    int vb=arr[bIndex];
    arr[aIndex]=vb;
    arr[bIndex]=va;
    arrIndex[va]=bIndex;
    arrIndex[vb]=aIndex;
}

inline void exchange(vector<int>& arr,int aIndex,int bIndex){
    int va=arr[aIndex];
    int vb=arr[bIndex];
    arr[aIndex]=vb;
    arr[bIndex]=va;

    // if(arr[aIndex]<0||arr[bIndex]<0){
    //     std::cerr<<"error"<<std::endl;
    // }
}

#endif