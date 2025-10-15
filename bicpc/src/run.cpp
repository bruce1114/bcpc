#include"tools/getArgs.hpp"
#include<iostream>
#include<chrono>
#include"bicpc/bicpc.cpp"

std::vector<std::vector<uint32_t> > ordered_cliqueuv[2];//zy

bool comp(int a,int b){//zy

    for(int i=0;i<ordered_cliqueuv[0][a].size()&&i<ordered_cliqueuv[0][b].size();++i){
        if(ordered_cliqueuv[0][a][i]!=ordered_cliqueuv[0][b][i]){
            return ordered_cliqueuv[0][a][i]<ordered_cliqueuv[0][b][i];
        }
    }
    return ordered_cliqueuv[0][a].size()<ordered_cliqueuv[0][b].size();
}

void outputbiclique(std::vector<std::vector<uint32_t> > * pcliqueuv, std::vector<uint32_t> * pOldLabels){//zy

    for(int i=0;i<pcliqueuv[0].size();++i){
        for(int j=0;j<pcliqueuv[0][i].size();++j){
            #ifdef RECOVEROLDLABEL
            pcliqueuv[0][i][j]=pOldLabels[0][pcliqueuv[0][i][j]];
            #else
            pcliqueuv[0][i][j]=pcliqueuv[0][i][j];
            #endif
        }
    }
    for(int i=0;i<pcliqueuv[1].size();++i){
        for(int j=0;j<pcliqueuv[1][i].size();++j){
            #ifdef RECOVEROLDLABEL
            pcliqueuv[1][i][j]=pOldLabels[1][pcliqueuv[1][i][j]];
            #else
            pcliqueuv[1][i][j]=pcliqueuv[1][i][j];
            #endif
        }
    }

    #ifdef SORTOUTPUT
    for(int i=0;i<pcliqueuv[0].size();++i){
        sort(pcliqueuv[0][i].begin(),pcliqueuv[0][i].end());
    }
    for(int i=0;i<pcliqueuv[1].size();++i){
        sort(pcliqueuv[1][i].begin(),pcliqueuv[1][i].end());
    }
    #endif

    int bicliquenum=pcliqueuv[0].size();
    std::vector<int> indexvec;
    indexvec.resize(bicliquenum);
    for(int i=0;i<bicliquenum;++i){
        indexvec[i]=i;
    }

    ordered_cliqueuv[0]=pcliqueuv[0];
    ordered_cliqueuv[1]=pcliqueuv[1];
    
    #ifdef SORTOUTPUT
    sort(indexvec.begin(),indexvec.end(),comp);
    #endif

    //print final results
    for(int i=0;i<indexvec.size();++i){
        int theindex=indexvec[i];
        for(int j=0;j<ordered_cliqueuv[0][theindex].size();++j){
            std::cout<<ordered_cliqueuv[0][theindex][j]<<" ";
        }
        std::cout<<std::endl;

        for(int j=0;j<ordered_cliqueuv[1][theindex].size();++j){
            std::cout<<ordered_cliqueuv[1][theindex][j]<<" ";
        }
        std::cout<<std::endl;
    }
}

void outputBCPCs(std::vector<std::vector<uint32_t> > * pcliqueuv,std::vector<uint32_t> * pOldLabels,unf& fans){
    for(int i=0;i<pcliqueuv[0].size();++i){
        for(int j=0;j<pcliqueuv[0][i].size();++j){
            #ifdef RECOVEROLDLABEL
            pcliqueuv[0][i][j]=pOldLabels[0][pcliqueuv[0][i][j]];
            #else
            pcliqueuv[0][i][j]=pcliqueuv[0][i][j];
            #endif
        }
    }
    for(int i=0;i<pcliqueuv[1].size();++i){
        for(int j=0;j<pcliqueuv[1][i].size();++j){
            #ifdef RECOVEROLDLABEL
            pcliqueuv[1][i][j]=pOldLabels[1][pcliqueuv[1][i][j]];
            #else
            pcliqueuv[1][i][j]=pcliqueuv[1][i][j];
            #endif
        }
    }

    for(int i=0;i<pcliqueuv[0].size();++i){
        sort(pcliqueuv[0][i].begin(),pcliqueuv[0][i].end());
    }
    for(int i=0;i<pcliqueuv[1].size();++i){
        sort(pcliqueuv[1][i].begin(),pcliqueuv[1][i].end());
    }

    int bicliquenum=pcliqueuv[0].size();
    std::vector<int> indexvec;
    indexvec.resize(bicliquenum);
    for(int i=0;i<bicliquenum;++i){
        indexvec[i]=i;
    }

    ordered_cliqueuv[0]=pcliqueuv[0];
    ordered_cliqueuv[1]=pcliqueuv[1];
    
    sort(indexvec.begin(),indexvec.end(),comp);

    std::vector<int> reindexvec;
    reindexvec.resize(bicliquenum);
    for(int i=0;i<bicliquenum;++i){
        reindexvec[indexvec[i]]=i;
    }

    std::unordered_map<int,std::vector<int> > bicpcResults;
    for(int i=0;i<bicliquenum;++i){
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
        for(int j=0;j<results[i].size();++j){
            results[i][j]=reindexvec[results[i][j]];
        }
    }

    for(int i=0;i<results.size();++i){
        sort(results[i].begin(),results[i].end());
    }
    sort(results.begin(),results.end());

    for(int i=0;i<results.size();++i){
        for(int j=0;j<results[i].size();++j){
            results[i][j]=indexvec[results[i][j]];
        }
    }

    //final print
    for(int i=0;i<results.size();++i){
        std::cout<<"qpcpc:"<<std::endl;
        for(int j=0;j<results[i].size();++j){
            int c=results[i][j];

            for(int k=0;k<pcliqueuv[0][c].size();++k){
                std::cout<<pcliqueuv[0][c][k]<<" ";
            }
            std::cout<<std::endl;
            for(int k=0;k<pcliqueuv[1][c].size();++k){
                std::cout<<pcliqueuv[1][c][k]<<" ";
            }
            std::cout<<std::endl;
        }

    }
}


int main(int argc, char* argv[]) {

    argsController ac(argc,argv);

    if(!ac.exist("-f")||!ac.exist("-p")||!ac.exist("-q")||!ac.exist("-a")||!ac.exist("-o")) {
        std::cerr<<"parameter error"<<std::endl;
        return 0;
    }

    std::string filename=ac["-f"];
    uint32_t p=std::stoi(ac["-p"]);
    uint32_t q=std::stoi(ac["-q"]);
    std::string alg=ac["-a"];
    std::cerr<<argv[0]<<" "<<filename<<" p "<<p<<" q "<<q<<" a "<<alg<<std::endl;
    int outputmark=std::stoi(ac["-o"]);

    BICPC bicpc(filename,p,q);

    // //debug
    // bicpc.g->writeOldLabels();
    // return 0;

    auto t1 = std::chrono::steady_clock::now();

    #ifdef COREREDUCTION
    bicpc.g->qpcoreReduction(p,q);
    // int usize,vsize;
    // bicpc.g->getcoreSize(usize,vsize);
    // std::cerr<<usize<<" "<<vsize<<std::endl;
    #endif
    bicpc.g->initializeHash();

    if(alg=="baseline"){
        bicpc.baseline();
    }else if(alg=="qbicpc"){
        bicpc.qbcpc();
        bicpc.quasi2finalBCPC();
    }else if(alg=="qbicpcL"){
        bicpc.qbcpcLeafRes();
        bicpc.quasi2finalBCPC();
    }else if(alg=="qbicpcLforabnode"){
        bicpc.qbcpcLeafResforabnode();
    }else if(alg=="qbicpcLbase"){
        bicpc.qbcpcLeafResbase();
        bicpc.quasi2finalBCPC();
    }else if(alg=="qbicpcLpoor"){
        bicpc.qbcpcLeafRespoor();
        bicpc.quasi2finalBCPC();
    }else if(alg=="qbicpcLpoors1"){
        bicpc.qbcpcLeafRespoor();
    }else if(alg=="pqnodecnt"){
        bicpc.qbcpcLeafRespoorCntpqnode();
    }else if(alg=="qplist"){
        bicpc.prepareAndListqpCliques();
    }else if(alg=="mbc"){
        bicpc.listMBCliques();
    }else if(alg=="mbctreefullp"){
        bicpc.listMBCliquesTreeFullPivot();
    }else if(alg=="mbctree"){
        bicpc.listMBCliquesTree();
    }else if(alg=="qpbcl"){
        bicpc.listMBCliques();
        bicpc.listAndConnectMBCbase();
    }else if(alg=="qpbcl_mbc"){
        bicpc.listMBCliques();
        bicpc.listAndConnectMBC();
    }else if(alg=="qpbcl_qcpc"){
        bicpc.qbcpc();
        std::cerr<<"group num: "<<bicpc.getGroupNum()<<std::endl;
        bicpc.listAndConnectMBC();
    }else if(alg=="qpbcl_qcpcL"){
        bicpc.qbcpcLeafRes();
        std::cerr<<"group num: "<<bicpc.getGroupNum()<<std::endl;
        bicpc.listAndConnectMBC();
    }else{
        std::cerr<<"alg error"<<std::endl;
        return 0;
    }
    auto t2 = std::chrono::steady_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
    std::cerr << "time: " << duration.count() << " ms" << std::endl;

    if(alg=="baseline"||alg=="qbicpc"||alg=="qpbcl_qcpc"||alg=="qpbcl_qcpcL"||alg=="qpbcl_mbc"||alg=="qbicpcL"||alg=="qpbcl"||alg=="qbicpcLbase"||alg=="qbicpcLpoor"){
        std::cerr<<"BCPC num: "<<bicpc.getGroupNum()<<std::endl;
    }else if(alg=="qplist"){
        #ifndef COUNTQPONLY
        bicpc.qpcliqueNum=bicpc.cliqueuv[0].size();
        #endif
        std::cerr<<"qpclique num: "<<bicpc.qpcliqueNum<<std::endl;
    }

    if(outputmark==0) return 0;

    if(alg=="mbc"||alg=="qplist"||alg=="mbctree"||alg=="mbctreefullp"){
        outputbiclique(bicpc.cliqueuv,bicpc.g->old_lables);
    }else if(alg!="pqnodecnt"&&alg!="qbicpcLpoors1"&&alg!="qbicpcLforabnode"){
        // bicpc.outputBCPCs();
        outputBCPCs(bicpc.cliqueuv,bicpc.g->old_lables,bicpc.fans);
    }
    
    return 0;
    

}