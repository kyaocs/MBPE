//
//  main.cpp
//  maximal_balanced_kplex_enum
//
//  Created by kai on 2023/3/21.
//

#include "Utility.h"
#include "LinearHeap.h"

void CTprune();
void comp_process_order();

void load_graph(string input_graph)
{
    string buffer;
    ifstream input_file(input_graph, ios::in);
    if (!input_file.is_open()){cout << "cannot open file : "<<input_graph<<endl;exit(1);}
    else{
        input_file >> n >> m;
        map<ui, int> * s_G = new map<ui, int>[n];
        ui tu, tv;
        int flag;
        while (input_file >> tu >> tv >> flag){
            if(tu == tv) continue;
            assert(tu >= 0 && tu < n); assert(tv >= 0 && tv < n);
            assert(flag == 1 || flag == -1);
            s_G[tu].insert(make_pair(tv, flag));
            s_G[tv].insert(make_pair(tu, flag));
        }
        m = 0; pm = 0; nm = 0;
        for(ui i = 0; i < n; i++){
            const map<ui, int> & nei = s_G[i];
            for(auto e : nei){
                if(e.second == 1) ++ pm;
                else ++ nm;
            }
            m += nei.size();
        }
        assert(m%2 == 0); assert(pm%2 == 0); assert(nm%2 == 0);
        m /= 2; pm /= 2; nm /= 2;
#ifdef _onlyShowG_
        cout<<"n:"<<n<<", m:"<<m<<" ("<<pm<<"+"<<nm<<")"<<endl;
        double r = (double)nm/m;
        cout<<"nm/m = "<<r<<endl;
#endif
        input_file.close();
        
        p_pstart = new ui[n+1];
        p_edges = new ui[2*pm];
        p_pend = new ui[n+1];
        
        n_pstart = new ui[n+1];
        n_edges = new ui[2*nm];
        n_pend = new ui[n+1];
        
        degree = new ui[n];
        p_degree = new ui[n];
        n_degree = new ui[n];
        
        //construct positive edges
        p_pstart[0] = 0;
        for(ui i = 0; i < n; i++){
            const map<ui, int> & nei = s_G[i];
            ui start_idx = p_pstart[i];
            int t_d = 0;
            for(auto e : nei){
                if(e.second == 1){
                    p_edges[start_idx ++] = e.first;
                    ++ t_d;
                }
            }
            p_pstart[i+1] = start_idx;
            p_degree[i] = t_d;
        }
        assert(p_pstart[n] == 2*pm);
        
        //construct negative edges
        n_pstart[0] = 0;
        for(ui i = 0; i < n; i++){
            const map<ui, int> & nei = s_G[i];
            ui start_idx = n_pstart[i];
            int t_d = 0;
            for(auto e : nei){
                if(e.second == -1){
                    n_edges[start_idx ++] = e.first;
                    ++ t_d;
                }
            }
            n_pstart[i+1] = start_idx;
            n_degree[i] = t_d;
        }
        assert(n_pstart[n] == 2*nm);
        max_deg = 0;
        for(ui i = 0; i < n; i++) {
            p_pend[i] = p_pstart[i+1];
            n_pend[i] = n_pstart[i+1];
            degree[i] = p_degree[i] + n_degree[i];
            if(degree[i] > max_deg) max_deg = degree[i];
        }
        true_id = new ui[n];
        original_id = new ui[n];
        for(ui i = 0; i < n; i++) original_id[i] = i;
        delete [] s_G;
    }
    
#ifdef _CaseStudy_
    input_graph.erase(input_graph.end() - 4, input_graph.end());
    input_graph.append("_map.txt"); //the mapping file should be named in this way.
    input_file.open(input_graph);
    if(!input_file.is_open()){ cout<<"cannot open map file !"<<endl; exit(1); }
    ui vid;
    string content;
    while (input_file >> content >> vid) id2str[vid] = content;
    input_file.close();
#endif
#ifdef _onlyShowG_
    process_order = new ui[n];
    for(ui i = 0; i < n; i++) process_order[i] = i;
    if(n==0) return;
    comp_process_order();
    cout<<"max_deg = "<<max_deg<<", max_core = "<<max_core<<endl;
    cout<<"only show the statistic of the graph! now will exit!"<<endl;
    exit(1);
#endif
}

void init_hash()
{
    sup_pp = new dense_hash_map<ui, int>[n];
    sup_nn = new dense_hash_map<ui, int>[n];
    sup_np = new dense_hash_map<ui, int>[n];
    e_sign = new dense_hash_map<ui, int>[n];
    e_del = new dense_hash_map<ui, bool>[n];
    
    
    for(ui i = 0; i < n; i++){
        sup_pp[i].resize(degree[i]);
        sup_nn[i].resize(degree[i]);
        sup_np[i].resize(degree[i]);
        e_sign[i].resize(degree[i]);
        e_del[i].resize(degree[i]);
    }
    for(ui i = 0; i < n; i++){
        sup_pp[i].set_empty_key(INF);
        sup_nn[i].set_empty_key(INF);
        sup_np[i].set_empty_key(INF);
        e_sign[i].set_empty_key(INF);
        e_del[i].set_empty_key(INF);
    }
    for(ui i = 0; i < n; i++){
        for(ui j = p_pstart[i]; j < p_pstart[i+1]; j++){
            sup_pp[i][p_edges[j]] = 0; sup_pp[p_edges[j]][i] = 0;
            sup_nn[i][p_edges[j]] = 0; sup_nn[p_edges[j]][i] = 0;
            sup_np[i][p_edges[j]] = 0; sup_np[p_edges[j]][i] = 0;
            e_sign[i][p_edges[j]] = 1; e_sign[p_edges[j]][i] = 1;
            e_del[i][p_edges[j]] = 0; e_del[p_edges[j]][i] = 0;
        }
        for(ui j = n_pstart[i]; j < n_pstart[i+1]; j++){
            sup_pp[i][n_edges[j]] = 0; sup_pp[n_edges[j]][i] = 0;
            sup_nn[i][n_edges[j]] = 0; sup_nn[n_edges[j]][i] = 0;
            sup_np[i][n_edges[j]] = 0; sup_np[n_edges[j]][i] = 0;
            e_sign[i][n_edges[j]] = -1; e_sign[n_edges[j]][i] = -1;
            e_del[i][n_edges[j]] = 0; e_del[n_edges[j]][i] = 0;
        }
    }
}

inline void reord_deg(ui & v1, ui & v2)
{
    if(degree[v1]>degree[v2]){
        ui t = v2;
        v2 = v1;
        v1 = t;
    }
}

void edge_reduction(ui & dele_count)
{
    int pp_thre = max(tau - 2*k, 0);
    int nn_thre = max(tau - 2*k + 2, 0);
    int np_thre = max(tau - 2*k + 1, 0);
    
    cout<<"ppthre = "<<pp_thre<<endl;
    cout<<"nn_thre = "<<nn_thre<<endl;
    cout<<"np_thre = "<<np_thre<<endl;
    
//    Timer t;
    //counting sort by degree
    ui * vs = new ui[n];
    ui * C = new ui[max_deg+1];
    memset(C, 0, sizeof(ui)*(max_deg+1));
    for(ui i = 0; i < n; i++) ++ C[degree[i]];
    for(ui i = 1; i <= max_deg; i++) C[i] += C[i-1];
    for(ui i = 0; i < n; i++) vs[--C[degree[i]]] = i;
        
    //count triangle
    ui * mark = new ui[n];
    memset(mark, 0, sizeof(ui)*n);
    ui * del = new ui[n];
    memset(del, 0, sizeof(ui)*n);
    
    init_hash();

    for(int i = n-1; i >= 0; i--){
        ui u = vs[i];
        for(ui j = p_pstart[u]; j < p_pstart[u+1]; j++) mark[p_edges[j]] = 1;
        for(ui j = n_pstart[u]; j < n_pstart[u+1]; j++) mark[n_edges[j]] = 2;
        for(ui j = p_pstart[u]; j < p_pstart[u+1]; j++){
            ui v = p_edges[j];
            if(!del[v]){
                for(ui k = p_pstart[v]; k < p_pstart[v+1]; k++){
                    ui w = p_edges[k];
                    if(!del[w] && mark[w] == 1 && w > v){
                        ++ sup_pp[u][v]; ++ sup_pp[v][u];
                        ++ sup_pp[u][w]; ++ sup_pp[w][u];
                        ++ sup_pp[v][w]; ++ sup_pp[w][v];
                    }
                }
                for(ui k = n_pstart[v]; k < n_pstart[v+1]; k++){
                    ui w = n_edges[k];
                    if(!del[w] && mark[w] == 2 && w > v){
                        ++ sup_nn[u][v]; ++ sup_nn[v][u];
                        ++ sup_np[u][w]; ++ sup_np[w][u];
                        ++ sup_np[v][w]; ++ sup_np[w][v];
                    }
                }
            }
        }
        for(ui j = n_pstart[u]; j < n_pstart[u+1]; j++){
            ui v = n_edges[j];
            if(!del[v]){
                for(ui k = p_pstart[v]; k < p_pstart[v+1]; k++){
                    ui w = p_edges[k];
                    if(!del[w] && mark[w] == 2 && w > v){
                        ++ sup_np[u][v]; ++ sup_np[v][u];
                        ++ sup_np[u][w]; ++ sup_np[w][u];
                        ++ sup_nn[v][w]; ++ sup_nn[w][v];
                    }
                }
                for(ui k = n_pstart[v]; k < n_pstart[v+1]; k++){
                    ui w = n_edges[k];
                    if(!del[w] && mark[w] == 1 && w > v){
                        ++ sup_np[u][v]; ++ sup_np[v][u];
                        ++ sup_np[v][w]; ++ sup_np[w][v];
                        ++ sup_nn[u][w]; ++ sup_nn[w][u];
                    }
                }
            }
        }
        for(ui j = p_pstart[u]; j < p_pstart[u+1]; j++) mark[p_edges[j]] = 0;
        for(ui j = n_pstart[u]; j < n_pstart[u+1]; j++) mark[n_edges[j]] = 0;
        del[u] = 1;
    }//u

    queue<pair<ui, ui>> Q;
    for(ui i = 0; i < n; i++){
        for(ui j = p_pstart[i]; j <p_pstart[i+1]; j++){
            ui v = p_edges[j];
            if(i < v){
                if(sup_pp[i][v] < pp_thre || sup_nn[i][v] < nn_thre){
                    Q.push(make_pair(i, v));
                }
            }
        }
        for(ui j = n_pstart[i]; j <n_pstart[i+1]; j++){
            ui v = n_edges[j];
            if(i < v){
                if(sup_np[i][v] < np_thre){
                    Q.push(make_pair(i, v));
                }
            }
        }
    }

    while (!Q.empty()) {
        pair<ui, ui> te = Q.front();
        ++ dele_count;
        Q.pop();
        ui u = te.first;
        ui v = te.second;
        
#ifdef _DEBUG_
        cout<<"delete edge ( "<<u<<" , "<<v<<" ) from G."<<endl;
#endif
        
        if(e_sign[u][v] == 1){ // e(u,v) is positive
            reord_deg(u, v);
            for(ui i = p_pstart[u]; i < p_pstart[u+1]; i++){
                ui w = p_edges[i];
                if(e_sign[v].find(w) != e_sign[v].end() && e_sign[v][w] == 1 && !e_del[u][w] && !e_del[v][w]){
                    if((sup_pp[u][w] --) == pp_thre && sup_nn[u][w] >= nn_thre)
                    {
                        Q.push(make_pair(u, w));
                    }
                    -- sup_pp[w][u];
                    if((sup_pp[v][w] --) == pp_thre && sup_nn[v][w] >= nn_thre)
                    {
                        Q.push(make_pair(v, w));
                    }
                    -- sup_pp[w][v];
                }
            }
            for(ui i = n_pstart[u]; i < n_pstart[u+1]; i++){
                ui w = n_edges[i];
                if(e_sign[v].find(w) != e_sign[v].end() && e_sign[v][w] == -1 && !e_del[u][w] && !e_del[v][w]){
                    if((sup_np[u][w] --) == np_thre)
                    {
                        Q.push(make_pair(u, w));
                    }
                    -- sup_np[w][u];
                    if((sup_np[v][w] --) == np_thre)
                    {
                        Q.push(make_pair(v, w));
                    }
                    -- sup_np[w][v];
                }
            }
        }
        else{ // e(u,v) is negative
            reord_deg(u, v);
            for(ui i = p_pstart[u]; i < p_pstart[u+1]; i++){
                ui w = p_edges[i];
                if(e_sign[v].find(w) != e_sign[v].end() && e_sign[v][w] == -1 && !e_del[u][w] && !e_del[v][w]){
                    if((sup_nn[u][w] --) == nn_thre && sup_pp[u][w] >= pp_thre)
                    {
                        Q.push(make_pair(u, w));
                    }
                    -- sup_nn[w][u];
                    if((sup_np[v][w] --) == np_thre)
                    {
                        Q.push(make_pair(v, w));
                    }
                    -- sup_np[w][v];
                }
            }
            for(ui i = n_pstart[u]; i < n_pstart[u+1]; i++){
                ui w = n_edges[i];
                if(e_sign[v].find(w) != e_sign[v].end() && e_sign[v][w] == 1 && !e_del[u][w] && !e_del[v][w]){
                    if((sup_np[u][w] --) == np_thre)
                    {
                        Q.push(make_pair(u, w));
                    }
                    -- sup_np[w][u];
                    if((sup_nn[v][w] --) == nn_thre && sup_pp[v][w] >= pp_thre)
                    {
                        Q.push(make_pair(v, w));
                    }
                    -- sup_nn[w][v];
                }
            }
        }
        e_del[u][v] = 1;
        e_del[v][u] = 1;
    }
    delete [] vs;
    delete [] C;
    delete [] mark;
    delete [] del;
//    t.restart();
}

void shrink_on_reduced_e()
{
    for(ui i = 0; i < n; i++){
        p_pend[i] = p_pstart[i];
        n_pend[i] = n_pstart[i];
    }
    for(ui i = 0; i < n; i++){
        for(ui j = p_pstart[i]; j < p_pstart[i+1]; j++){
            ui u = p_edges[j]; //e(i,u)
            if(!e_del[i][u]){
                p_edges[p_pend[i]++] = u;
            }
        }
    }
    for(ui i = 0; i < n; i++){
        for(ui j = n_pstart[i]; j < n_pstart[i+1]; j++){
            ui u = n_edges[j]; //e(i,u)
            if(!e_del[i][u]){
                n_edges[n_pend[i]++] = u;
            }
        }
    }
    
    //update degree
    for(ui i = 0; i < n; i++){
        p_degree[i] = p_pend[i] - p_pstart[i];
        n_degree[i] = n_pend[i] - n_pstart[i];
        degree[i] = p_degree[i] + n_degree[i];
    }
    
#ifdef _DEBUG_
    cout<<"after shrink_on_reduced_e(), : "<<endl;
    for(ui i = 0; i < n; i++){
        cout<<"for vertex "<<i<<" : "<<endl;
        cout<<"  + neis : ";
        for(ui j = p_pstart[i]; j < p_pend[i]; j++){
            cout<<p_edges[j]<<", ";
        }cout<<endl;
        cout<<"  - neis : ";
        for(ui j = n_pstart[i]; j < n_pend[i]; j++){
            cout<<n_edges[j]<<", ";
        }cout<<endl;
    }
#endif
}

void comp_process_order()
{
    max_core = 0;
    core = new ui[n];
    memset(core, 0, sizeof(ui)*n);
    ListLinearHeap *linear_heap = new ListLinearHeap(n, n-1);
    linear_heap->init(n, n-1, process_order, degree);
    for(ui i = 0; i < n; i ++) {
        ui u, key;
        linear_heap->pop_min(u, key);
        if(key > max_core) max_core = key;
        process_order[i] = u;
        core[u] = max_core;
        for(ui j = p_pstart[u];j < p_pend[u];j ++) if(core[p_edges[j]] == 0)
            linear_heap->decrement(p_edges[j], 1);
        for(ui j = n_pstart[u];j < n_pend[u];j ++) if(core[n_edges[j]] == 0)
            linear_heap->decrement(n_edges[j], 1);
    }
    delete linear_heap;
#ifdef _DEBUG_
    cout<<"peel sequence : ";
    for(ui i = 0; i< n; i++) cout<<process_order[i]<<", "; cout<<endl;
    cout<<"core  number  : ";
    for(ui i = 0; i< n; i++) cout<<core[process_order[i]]<<", ";cout<<endl;
#endif
}

void obtain_CL_CR_XL_XR(ui u, vector<ui> &CL, vector<ui> &CR, vector<ui> &XL, vector<ui> &XR)
{
    for(ui i = p_pstart[u]; i < p_pend[u]; i++) {
        ui v = p_edges[i]; //u's + nei
        if(ver_rank[v] > ver_rank[u]) {
            if(inCL[v] == 0) {CL.push_back(v); inCL[v] = 1;}
            if(inCR[v] == 0) {CR.push_back(v); inCR[v] = 1;}
            for(ui j = p_pstart[v]; j < p_pend[v]; j++) {
                ui w = p_edges[j]; //u's ++ nei
                if(ver_rank[w] > ver_rank[u]) {
                    if(inCL[w] == 0) {CL.push_back(w); inCL[w] = 1;}
                }
                else if (ver_rank[w] < ver_rank[u]) {
                    if(inXL[w] == 0) {XL.push_back(w); inXL[w] = 1;}
                }
            }
            for(ui j = n_pstart[v]; j < n_pend[v]; j++) {
                ui w = n_edges[j]; //u's +- nei
                if(ver_rank[w] > ver_rank[u]) {
                    if(inCR[w] == 0) {CR.push_back(w); inCR[w] = 1;}
                }
                else if (ver_rank[w] < ver_rank[u]) {
                    if(inXR[w] == 0) {XR.push_back(w); inXR[w] = 1;}
                }
            }
        }
        else {
            assert(ver_rank[v] < ver_rank[u]);
            if(inXL[v] == 0) {XL.push_back(v); inXL[v] = 1;}
            if(inXR[v] == 0) {XR.push_back(v); inXR[v] = 1;}
            for(ui j = p_pstart[v]; j < p_pend[v]; j++) {//unnecessary
                ui w = p_edges[j];
                if(ver_rank[w] < ver_rank[u] && inXL[w] == 0) {
                    XL.push_back(w);
                    inXL[w] = 1;
                }
            }
            for(ui j = n_pstart[v]; j < n_pend[v]; j++) {//unnecessary
                ui w = n_edges[j];
                if(ver_rank[w] < ver_rank[u] && inXR[w] == 0) {
                    XR.push_back(w);
                    inXR[w] = 1;
                }
            }
        }
    }
    for(ui i = n_pstart[u]; i < n_pend[u]; i++) {
        ui v = n_edges[i]; //u's - nei
        if(ver_rank[v] > ver_rank[u]) {
            if(inCL[v] == 0) {CL.push_back(v); inCL[v] = 1;}
            if(inCR[v] == 0) {CR.push_back(v); inCR[v] = 1;}
            for(ui j = p_pstart[v]; j < p_pend[v]; j++) {
                ui w = p_edges[j]; //u's -+ nei
                if(ver_rank[w] > ver_rank[u]) {
                    if(inCR[w] == 0) {CR.push_back(w); inCR[w] = 1;}
                }
                else if(ver_rank[w] < ver_rank[u]) {
                    if(inXR[w] == 0) {XR.push_back(w); inXR[w] = 1;}
                }
            }
            for(ui j = n_pstart[v]; j < n_pend[v]; j++) {
                ui w = n_edges[j]; //u's -- nei
                if(ver_rank[w] > ver_rank[u]) {
                    if(inCL[w] == 0) {CL.push_back(w); inCL[w] = 1;}
                }
                else if(ver_rank[w] < ver_rank[u]) {
                    if(inXL[w] == 0) {XL.push_back(w); inXL[w] = 1;}
                }
            }
        }
        else {
            assert(ver_rank[v] < ver_rank[u]);
            if(inXL[v] == 0) {XL.push_back(v); inXL[v] = 1;}
            if(inXR[v] == 0) {XR.push_back(v); inXR[v] = 1;}
            for(ui j = p_pstart[v]; j < p_pend[v]; j++) {//unnecessary
                ui w = p_edges[j];
                if(ver_rank[w] < ver_rank[u] && inXR[w] == 0) {
                    XR.push_back(w);
                    inXR[w] = 1;
                }
            }
            for(ui j = n_pstart[v]; j < n_pend[v]; j++) {//unnecessary
                ui w = n_edges[j];
                if(ver_rank[w] < ver_rank[u] && inXL[w] == 0) {
                    XL.push_back(w);
                    inXL[w] = 1;
                }
            }
        }
    }
#ifdef _DEBUG_
    cout<<"CL ("<<CL.size()<<"): "; for(auto e : CL) cout<<e<<","; cout<<endl;
    cout<<"CR ("<<CR.size()<<"): "; for(auto e : CR) cout<<e<<","; cout<<endl;
    cout<<"XL ("<<XL.size()<<"): "; for(auto e : XL) cout<<e<<","; cout<<endl;
    cout<<"XR ("<<XR.size()<<"): "; for(auto e : XR) cout<<e<<","; cout<<endl;
#endif
}

void construct_initial_subgraph(ui u, vector<ui>&CL, vector<ui>&CR,vector<ui>&P)
{
    vector<ui> S;
    S.push_back(u);
    assert(part[u] == 0);
    part[u] = 1;
    S.insert(S.end(), CL.begin(), CL.end());
    S.insert(S.end(), CR.begin(), CR.end());
    S.insert(S.end(), P.begin(), P.end());
    for(auto &v : S) {
        t_p_pend[v] = t_p_pstart[v];
        assert(p_pstart[v]==t_p_pstart[v]);
        for(ui i = p_pstart[v]; i < p_pend[v]; i++) if(part[p_edges[i]] != 0) {
            t_p_edges[t_p_pend[v]] = p_edges[i];
            t_p_pend[v]++;
        }
        t_n_pend[v] = t_n_pstart[v];
        assert(n_pstart[v]==t_n_pstart[v]);
        for(ui i = n_pstart[v]; i < n_pend[v]; i++) if(part[n_edges[i]] != 0) {
            t_n_edges[t_n_pend[v]] = n_edges[i];
            t_n_pend[v]++;
        }
        //be aware we use degree, p_degree, n_degree for t_...
        p_degree[v] = t_p_pend[v] - t_p_pstart[v];
        n_degree[v] = t_n_pend[v] - t_n_pstart[v];
        degree[v] = p_degree[v] + n_degree[v];
    }
    part[u] = 0;
#ifdef _DEBUG_
    cout<<"S:"; for(auto e : S) cout<<e<<","; cout<<endl;
    cout<<"the constructed initial subgraph:"<<endl;
    for(auto e : S) {
        cout<<e<<":"<<endl;
        cout<<"\tp nei :"; for(ui i = t_p_pstart[e]; i < t_p_pend[e]; i++) cout<<t_p_edges[i]<<","; cout<<endl;
        cout<<"\tn nei :"; for(ui i = t_n_pstart[e]; i < t_n_pend[e]; i++) cout<<t_n_edges[i]<<","; cout<<endl;
    }
#endif
#ifdef _CostlyAssert_
    for(ui i = 0; i < n; i++) assert(ver_del[i] == 0);
#endif
}

void construct_initial_subgraph(ui u, vector<ui>&CL, vector<ui>&CR,vector<ui>&P,vector<ui>&XL,vector<ui>&XR)
{
    vector<ui> S;
    S.push_back(u);
    ver_del[u] = 1;
    for(auto &e : CL) {
        S.push_back(e);
        ver_del[e] = 1;
    }
    for(auto &e : CR) {
        S.push_back(e);
        ver_del[e] = 1;
    }
    for(auto &e : P) {
        S.push_back(e);
        ver_del[e] = 1;
    }
    for(auto &e : XL) {
        S.push_back(e);
        ver_del[e] = 1;
    }
    for(auto &e : XR) if(ver_del[e] == 0){
        S.push_back(e);
        ver_del[e] = 1;
    }
    for(auto v : S) {
        t_p_pend[v] = t_p_pstart[v];
        for(ui i = p_pstart[v]; i < p_pend[v]; i++) if(ver_del[p_edges[i]]) {
            t_p_edges[t_p_pend[v]] = p_edges[i];
            t_p_pend[v]++;
        }
        t_n_pend[v] = t_n_pstart[v];
        for(ui i = n_pstart[v]; i < n_pend[v]; i++) if(ver_del[n_edges[i]]) {
            t_n_edges[t_n_pend[v]] = n_edges[i];
            t_n_pend[v]++;
        }
        //be aware we use degree, p_degree, n_degree for t_...
        p_degree[v] = t_p_pend[v] - t_p_pstart[v];
        n_degree[v] = t_n_pend[v] - t_n_pstart[v];
        degree[v] = p_degree[v] + n_degree[v];
    }
    for(auto e : S) ver_del[e] = 0;
#ifdef _DEBUG_
    cout<<"S:"; for(auto e : S) cout<<e<<","; cout<<endl;
    cout<<"the constructed initial subgraph:"<<endl;
    for(auto e : S) {
        cout<<e<<":"<<endl;
        cout<<"\tp nei :"; for(ui i = t_p_pstart[e]; i < t_p_pend[e]; i++) cout<<t_p_edges[i]<<","; cout<<endl;
        cout<<"\tn nei :"; for(ui i = t_n_pstart[e]; i < t_n_pend[e]; i++) cout<<t_n_edges[i]<<","; cout<<endl;
    }
#endif
#ifdef _CostlyAssert_
    for(ui i = 0; i < n; i++) assert(ver_del[i] == 0);
#endif
}

void build_M(ui u, vector<ui> &CL, vector<ui> &CR, vector<ui> &XL, vector<ui> &XR)
{
    vector<ui> S;
//    ui idx = 0;
//    tran[u] = idx ++;
    S.push_back(u);
    for(auto e : CL) {
//        tran[e] = idx ++;
        S.push_back(e);
    }
    for(auto e : XL) {
//        tran[e] = idx ++;
        S.push_back(e);
    }
    for(auto e : CR) if (inCL[e] == 0) {
//        tran[e] = idx ++;
        S.push_back(e);
    }
    for(auto e : XR) if (inXL[e] == 0) {
//        tran[e] = idx ++;
        S.push_back(e);
    }
    for(ui i = 0; i < n; i++) for(ui j = 0; j < n; j++) M[i][j] = 0;
    for(auto v : S) {
        for(ui j = p_pstart[v]; j < p_pend[v]; j++) {
            ui w = p_edges[j];
            if(inCL[w] || inCR[w] || inXL[w] || inXR[w]) {
                M[v][w] = 1; M[w][v] = 1;
            }
        }
        for(ui j = n_pstart[v]; j < n_pend[v]; j++) {
            ui w = n_edges[j];
            if(inCL[w] || inCR[w] || inXL[w] || inXR[w]) {
                M[v][w] = -1; M[w][v] = -1;
            }
        }
    }
#ifdef _DEBUG_
    cout<<"M:"<<endl;
    sort(S.begin(),S.end());
    cout<<"\t";
    for(auto e : S) cout<<e<<"\t"; cout<<endl;
    for(auto e : S) {
        cout<<e<<":\t"; for(auto x : S) cout<<M[e][x]<<"\t"; cout<<endl;
    } cout<<endl;
#endif
}

void build_M(ui u, vector<ui> &CL, vector<ui> &CR, vector<ui> &P, vector<ui> &XL, vector<ui> &XR)
{
    vector<ui> S;
//    ui idx = 0;
//    tran[u] = idx ++;
    S.push_back(u);
    for(auto &e : CL) {
//        tran[e] = idx ++;
        S.push_back(e);
    }
    for(auto &e : CR) {
//        tran[e] = idx ++;
        S.push_back(e);
    }
    for(auto &e : P) {
//        tran[e] = idx ++;
        S.push_back(e);
    }
    for(auto &e : XL) {
//        tran[e] = idx ++;
        S.push_back(e);
    }
    for(auto &e : XR) if (inXL[e] == 0) {
//        tran[e] = idx ++;
        S.push_back(e);
    }
    
#ifdef _CostlyAssert_
    if(S.size() > 1) {
        sort(S.begin(), S.end());
        for(ui i = 0; i < S.size() - 1; i ++) assert(S[i+1] > S[i]);
    }
#endif
    
//    for(ui i = 0; i < idx; i++) for(ui j = 0; j < idx; j++) M[i][j] = 0;
    for(ui i = 0; i < n; i++) for(ui j = 0; j < n; j++) M[i][j] = 0;
    for(auto v : S) {
        for(ui j = p_pstart[v]; j < p_pend[v]; j++) {
            ui w = p_edges[j];
            if(part[w] != 0 || inXL[w] || inXR[w]) {
                M[v][w] = 1; M[w][v] = 1;
            }
        }
        for(ui j = n_pstart[v]; j < n_pend[v]; j++) {
            ui w = n_edges[j];
            if(part[w] != 0 || inXL[w] || inXR[w]) {
                M[v][w] = -1; M[w][v] = -1;
            }
        }
    }
#ifdef _DEBUG_
    cout<<"M:"<<endl;
    sort(S.begin(),S.end());
    cout<<"\t";
    for(auto e : S) cout<<e<<"\t"; cout<<endl;
    for(auto e : S) {
        cout<<e<<":\t"; for(auto x : S) cout<<M[e][x]<<"\t"; cout<<endl;
    } cout<<endl;
#endif
}

bool obtain_partial_kplex_capacity(vector<ui> &KL, vector<ui> &KR)
{
    for(auto &e : KL) { p_degree[e] = 0; n_degree[e] = 0;} //for VR
    for(auto &e : KR) { p_degree[e] = 0; n_degree[e] = 0;} //for VR
    for(ui i = 0; i < KL.size(); i++) {
        ui u = KL[i];
        for(ui j = i + 1; j < KL.size(); j++) {
            ui v = KL[j];
            if(M[u][v] == 1) {
                ++ p_degree[u]; ++ p_degree[v];
            }
        }
    }
    for(ui i = 0; i < KR.size(); i++) {
        ui u = KR[i];
        for(ui j = i + 1; j < KR.size(); j++) {
            ui v = KR[j];
            if(M[u][v] == 1) {
                ++ p_degree[u]; ++ p_degree[v];
            }
        }
    }
    for(auto &u : KL) for(auto &v : KR) {
        if(M[u][v] == -1) {
            ++ n_degree[u]; ++ n_degree[v];
        }
    }
    ui Ksize = (ui)KL.size() + (ui)KR.size();
    for(auto &e : KL) {
        assert(Ksize > (p_degree[e] + n_degree[e]));
        int nneinum = Ksize - (p_degree[e] + n_degree[e]);
        partial_kplex_cap[e] = k - nneinum;
        if(partial_kplex_cap[e] < 0) return true;
    }
    for(auto &e : KR) {
        assert(Ksize > (p_degree[e] + n_degree[e]));
        int nneinum = Ksize - (p_degree[e] + n_degree[e]);
        partial_kplex_cap[e] = k - nneinum;
        if(partial_kplex_cap[e] < 0) return true;
    }
    return false;
}

ui get_mindeg_vertex(vector<ui>&KL, vector<ui>&KR, vector<ui>&CL, vector<ui>&CR, int & ustar_side, ui & ustar_pos)
{
#ifndef _VRinEnum_
    for(ui i = 0; i < CL.size(); i ++) {
        ui u = CL[i];
        for(ui j = i + 1; j < CL.size(); j ++) {
            ui v = CL[j];
            if(M[u][v]==1) {
                ++ p_degree[u]; ++ p_degree[v];
            }
        }
        for(auto &v : KL) if(M[u][v]==1) {
            ++ p_degree[u]; ++ p_degree[v];
        }
        for(auto &v : KR) if(M[u][v]==-1) {
            ++ n_degree[u]; ++ n_degree[v];
        }
    }
    for(ui i = 0; i < CR.size(); i ++) {
        ui u = CR[i];
        for(ui j = i + 1; j < CR.size(); j ++) {
            ui v = CR[j];
            if(M[u][v]==1) {
                ++ p_degree[u]; ++ p_degree[v];
            }
        }
        for(auto &v : KL) if(M[u][v]==-1) {
            ++ n_degree[u]; ++ n_degree[v];
        }
        for(auto &v : KR) if(M[u][v]==1) {
            ++ p_degree[u]; ++ p_degree[v];
        }
    }
    for(auto &u : CL) for(auto &v : CR) {
        if(M[u][v]==-1) {
            ++ n_degree[u]; ++ n_degree[v];
        }
    }
#endif
#ifdef _DEBUG_
    cout<<"degree vertex in KL KR CL CR : "<<endl;
    cout<<"KL : "<<endl;
    for(auto e : KL) cout<<"\t"<<e<<":"<<p_degree[e]<<","<<n_degree[e]<<endl;
    cout<<"KR : "<<endl;
    for(auto e : KR) cout<<"\t"<<e<<":"<<p_degree[e]<<","<<n_degree[e]<<endl;
    cout<<"CL : "<<endl;
    for(auto e : CL) cout<<"\t"<<e<<":"<<p_degree[e]<<","<<n_degree[e]<<endl;
    cout<<"CR : "<<endl;
    for(auto e : CR) cout<<"\t"<<e<<":"<<p_degree[e]<<","<<n_degree[e]<<endl;
#endif
    ui ustar = 0;
    int tmpdeg = INF;
    for(auto &u : KL) if(p_degree[u]+n_degree[u] < tmpdeg) {
        tmpdeg = p_degree[u]+n_degree[u];
        ustar = u;
        ustar_side = 1;
    }
    for(auto &u : KR) if(p_degree[u]+n_degree[u] < tmpdeg) {
        tmpdeg = p_degree[u]+n_degree[u];
        ustar = u;
        ustar_side = 2;
    }
    for(ui i = 0; i < CL.size(); i++) {
        ui u = CL[i];
        if(p_degree[u]+n_degree[u] < tmpdeg) {
            tmpdeg = p_degree[u]+n_degree[u];
            ustar = u;
            ustar_side = 3;
            ustar_pos = i;
        }
    }
    for(ui i = 0; i < CR.size(); i++) {
        ui u = CR[i];
        if(p_degree[u]+n_degree[u] < tmpdeg) {
            tmpdeg = p_degree[u]+n_degree[u];
            ustar = u;
            ustar_side = 4;
            ustar_pos = i;
        }
    }
    return ustar;
}

bool VR_in_Enum_baseline(vector<ui>&KL, vector<ui>&KR, vector<ui>&CL, vector<ui>&CR)
{
    vector<ui> Svec;
    for(auto &e : KL) {
        p_degree[e] = 0;
        n_degree[e] = 0;
        inQv[e] = 0;
        Svec.push_back(e);
    }
    for(auto &e : KR) {
        p_degree[e] = 0;
        n_degree[e] = 0;
        inQv[e] = 0;
        Svec.push_back(e);
    }
    for(auto &e : CL) {
        p_degree[e] = 0;
        n_degree[e] = 0;
        inQv[e] = 0;
        Svec.push_back(e);
        ver_del[e] = 1;
    }
    for(auto &e : CR) {
        p_degree[e] = 0;
        n_degree[e] = 0;
        inQv[e] = 0;
        if(ver_del[e]==0)
            Svec.push_back(e);
    }
    for(auto &e : CL) ver_del[e] = 0;

    for(ui i = 0; i < Svec.size(); i ++) {
        ui u = Svec[i];
        for(ui j = i + 1; j < Svec.size(); j ++) {
            ui v = Svec[j];
            if(M[u][v]==1) {
                ++ p_degree[u]; ++ p_degree[v];
            }
            if(M[u][v]==-1) {
                ++ n_degree[u]; ++ n_degree[v];
            }
        }
    }
    for(auto &e : KL) if(p_degree[e]+k<tau || n_degree[e]+k-1<tau || p_degree[e]+n_degree[e]+k<2*tau) return true;
    for(auto &e : KR) if(p_degree[e]+k<tau || n_degree[e]+k-1<tau || p_degree[e]+n_degree[e]+k<2*tau) return true;
    queue<ui> Q;
    for(auto &e : Svec) if(p_degree[e]+k<tau || n_degree[e]+k-1<tau || p_degree[e]+n_degree[e]+k<2*tau) {
        Q.push(e); inQv[e] = 1;
    }
    while (!Q.empty()) {
        ui u = Q.front();
        Q.pop();
        assert(inQv[u]==1);

        for(auto &e : Svec) {
            if(!inQv[e] && M[u][e]==1) {
                assert(p_degree[e]>0);
                --p_degree[e];
                if(p_degree[e]+k<tau || p_degree[e]+n_degree[e]+k<2*tau) {Q.push(e); inQv[e] = 1;}
            }
            if(!inQv[e] && M[u][e]==-1) {
                assert(n_degree[e]>0);
                --n_degree[e];
                if(n_degree[e]+k-1<tau || p_degree[e]+n_degree[e]+k<2*tau) {Q.push(e); inQv[e] = 1;}
            }
        }
    }
    for(auto &e : KL) if(inQv[e]) return true;
    for(auto &e : KR) if(inQv[e]) return true;
    vector<ui> tmpCL, tmpCR;
    for(auto &e : CL) if(inQv[e]==0) tmpCL.push_back(e);
    for(auto &e : CR) if(inQv[e]==0) tmpCR.push_back(e);
    CL=tmpCL;
    CR=tmpCR;
    return false;
}

bool VR_in_Enum(vector<ui>&KL, vector<ui>&KR, vector<ui>&CL, vector<ui>&CR)
{
    for(ui i = 0; i < CL.size(); i ++) {
        ui u = CL[i];
        for(ui j = i + 1; j < CL.size(); j ++) {
            ui v = CL[j];
            if(M[u][v]==1) {
                ++ p_degree[u]; ++ p_degree[v];
            }
        }
        for(auto &v : KL) if(M[u][v]==1) {
            ++ p_degree[u]; ++ p_degree[v];
        }
        for(auto &v : KR) if(M[u][v]==-1) {
            ++ n_degree[u]; ++ n_degree[v];
        }
    }
    for(ui i = 0; i < CR.size(); i ++) {
        ui u = CR[i];
        for(ui j = i + 1; j < CR.size(); j ++) {
            ui v = CR[j];
            if(M[u][v]==1) {
                ++ p_degree[u]; ++ p_degree[v];
            }
        }
        for(auto &v : KL) if(M[u][v]==-1) {
            ++ n_degree[u]; ++ n_degree[v];
        }
        for(auto &v : KR) if(M[u][v]==1) {
            ++ p_degree[u]; ++ p_degree[v];
        }
    }
    for(auto &u : CL) for(auto &v : CR) {
        if(M[u][v]==-1) {
            ++ n_degree[u]; ++ n_degree[v];
        }
    }
    queue<ui> Q;
    for(auto &u : KL) inQv[u]=0;
    for(auto &u : KR) inQv[u]=0;
    for(auto &u : CL) inQv[u]=0;
    for(auto &u : CR) inQv[u]=0;
    
    for(auto &e : KL) if(p_degree[e]+k<tau || n_degree[e]+k-1<tau || p_degree[e]+n_degree[e]+k<2*tau) return true;
    for(auto &e : KR) if(p_degree[e]+k<tau || n_degree[e]+k-1<tau || p_degree[e]+n_degree[e]+k<2*tau) return true;
    
    assert(Q.empty());
    for(auto &e : CL) if(p_degree[e]+k<tau || n_degree[e]+k-1<tau || p_degree[e]+n_degree[e]+k<2*tau) {
        Q.push(e); inQv[e] = 1;
    }
    for(auto &e : CR) if(p_degree[e]+k<tau || n_degree[e]+k-1<tau || p_degree[e]+n_degree[e]+k<2*tau) {
        Q.push(e); inQv[e] = 2;
    }
    while (!Q.empty()) {
        ui u = Q.front();
        Q.pop();
        ++ vr_prune_I;
        assert(inQv[u]==1 || inQv[u]==2);
        if(inQv[u]==1) {
            for(auto &e : KL) if(!inQv[e] && M[u][e]==1) {
                assert(p_degree[e]>0);
                --p_degree[e];
                if(p_degree[e]+k<tau || p_degree[e]+n_degree[e]+k<2*tau) {Q.push(e); inQv[e] = 1;}
            }
            for(auto &e : KR) if(!inQv[e] && M[u][e]==-1) {
                assert(n_degree[e]>0);
                --n_degree[e];
                if(n_degree[e]+k-1<tau || p_degree[e]+n_degree[e]+k<2*tau) {Q.push(e); inQv[e] = 2;}
            }
            for(auto &e : CL) if(!inQv[e] && M[u][e]==1) {
                assert(p_degree[e]>0);
                --p_degree[e];
                if(p_degree[e]+k<tau || p_degree[e]+n_degree[e]+k<2*tau) {Q.push(e); inQv[e] = 1;}
            }
            for(auto &e : CR) if(!inQv[e] && M[u][e]==-1) {
                assert(n_degree[e]>0);
                --n_degree[e];
                if(n_degree[e]+k-1<tau || p_degree[e]+n_degree[e]+k<2*tau) {Q.push(e); inQv[e] = 2;}
            }
        }
        else {
            for(auto &e : KL) if(!inQv[e] && M[u][e]==-1) {
                assert(n_degree[e]>0);
                --n_degree[e];
                if(n_degree[e]+k-1<tau || p_degree[e]+n_degree[e]+k<2*tau) {Q.push(e); inQv[e] = 1;}
            }
            for(auto &e : KR) if(!inQv[e] && M[u][e]==1) {
                assert(p_degree[e]>0);
                --p_degree[e];
                if(p_degree[e]+k<tau || p_degree[e]+n_degree[e]+k<2*tau) {Q.push(e); inQv[e] = 2;}
            }
            for(auto &e : CL) if(!inQv[e] && M[u][e]==-1) {
                assert(n_degree[e]>0);
                --n_degree[e];
                if(n_degree[e]+k-1<tau || p_degree[e]+n_degree[e]+k<2*tau) {Q.push(e); inQv[e] = 1;}
            }
            for(auto &e : CR) if(!inQv[e] && M[u][e]==1) {
                assert(p_degree[e]>0);
                --p_degree[e];
                if(p_degree[e]+k<tau || p_degree[e]+n_degree[e]+k<2*tau) {Q.push(e); inQv[e] = 2;}
            }
        }
    }
    for(auto &e : KL) if(inQv[e]) return true;
    for(auto &e : KR) if(inQv[e]) return true;
    vector<ui> tmpCL, tmpCR;
    for(auto &e : CL) if(inQv[e]==0) tmpCL.push_back(e);
    for(auto &e : CR) if(inQv[e]==0) tmpCR.push_back(e);
    CL=tmpCL;
    CR=tmpCR;
    return false;
}

void obtain_PNdeg_PNnei_in_Enum(vector<ui>&KL, vector<ui>&KR, vector<ui>&CL, vector<ui>&CR)
{
    vector<ui> S;
    S.insert(S.end(), KL.begin(), KL.end());
    S.insert(S.end(), KR.begin(), KR.end());
    S.insert(S.end(), CL.begin(), CL.end());
    S.insert(S.end(), CR.begin(), CR.end());
    for(auto &u : S) {
        Pnei[u].clear(); Pnnei[u].clear();
        Nnei[u].clear(); Nnnei[u].clear();
    }
    for(ui i = 0; i < S.size(); i ++) {
        ui u = S[i];
        for(ui j = i + 1; j < S.size(); j++) {
            ui v = S[j];
            if(part[u]==part[v]) { //u, v are from same side
                if(M[u][v]==1) { Pnei[u].push_back(v); Pnei[v].push_back(u); }
                else { Pnnei[u].push_back(v); Pnnei[v].push_back(u); }
            }
            else {
                if(M[u][v]==-1) { Nnei[u].push_back(v); Nnei[v].push_back(u); }
                else { Nnnei[u].push_back(v); Nnnei[v].push_back(u); }
            }
        }
    }
    for(auto &u : S) {
        p_degree[u] = (ui)Pnei[u].size();
        p_r_degree[u] = (ui)Pnnei[u].size();
        n_degree[u] = (ui)Nnei[u].size();
        n_r_degree[u] = (ui)Nnnei[u].size();
    }
#ifdef _DEBUG_
    cout<<"KCdeg and LRnei LRnnei: "<<endl;
    for(auto &u : S) {
        cout<<"vertex "<<u<<":"<<endl;
        cout<<"\tPnei ("<<Pnei[u].size()<<")   : "; for(auto &v : Pnei[u]) cout<<v<<","; cout<<endl;
        cout<<"\tPnnei ("<<Pnnei[u].size()<<") : "; for(auto &v : Pnnei[u]) cout<<v<<","; cout<<endl;
        cout<<"\tNnei ("<<Nnei[u].size()<<")   : "; for(auto &v : Nnei[u]) cout<<v<<","; cout<<endl;
        cout<<"\tNnnei ("<<Nnnei[u].size()<<") : "; for(auto &v : Nnnei[u]) cout<<v<<","; cout<<endl;
    }
#endif
}

void ER_in_Enum(vector<ui>&KL, vector<ui>&KR, vector<ui>&CL, vector<ui>&CR, vector<pair<ui, ui>>&rmpe, vector<pair<ui, ui>>&rmne)
{
    vector<ui> S;
    S.insert(S.end(), KL.begin(), KL.end());
    S.insert(S.end(), KR.begin(), KR.end());
    S.insert(S.end(), CL.begin(), CL.end());
    S.insert(S.end(), CR.begin(), CR.end());
    for(auto &u : S) {
        assert(inQv[u]==0);
        core[u] = p_degree[u]+n_degree[u];
    }
    sort(S.begin(), S.end(), mycomp_r);//decreasing order w.r.t. degree in KC

    //initialize edge support
    for(auto &u : S) {
        assert(inQv[u]==0);
        v_sta[u] = 0;
        for(auto &v : Pnei[u]) if(inQv[v]==0) {
            assert((part[u]==1 && part[v]==1) || (part[u]==2 && part[v]==2));
            assert(M[u][v] == 1);
            tri_pp[u][v] = 0; tri_nn[u][v] = 0;
        }
        for(auto &v : Nnei[u]) if(inQv[v]==0) {
            assert((part[u]==1 && part[v]==2) || (part[u]==2 && part[v]==1));
            assert(M[u][v] == -1);
            tri_np[u][v] = 0; tri_pn[u][v] = 0;
        }
    }
    //triangle count
    for(auto &u : S) {
        assert(inQv[u]==0 && v_sta[u]==0);
        for(auto &v : Pnei[u]) if(inQv[v]==0 && v_sta[v]==0) {
            for(auto &w : Pnei[v]) if(inQv[w]==0 && v_sta[w]==0) {
                if(v<w && M[u][w] == 1) {
                    ++tri_pp[u][v]; ++tri_pp[v][u];
                    ++tri_pp[u][w]; ++tri_pp[w][u];
                    ++tri_pp[v][w]; ++tri_pp[w][v];
                }
            }
            for(auto &w : Nnei[v]) if(inQv[w]==0 && v_sta[w]==0) {
                if(v<w && M[u][w] == -1) {
                    ++tri_nn[u][v]; ++tri_nn[v][u];
                    if(part[u]==1) {
                        ++tri_pn[v][w]; ++tri_pn[w][v];
                        ++tri_pn[u][w]; ++tri_pn[w][u];
                    }
                    else {
                        assert(part[u]==2);
                        ++tri_np[v][w]; ++tri_np[w][v];
                        ++tri_np[u][w]; ++tri_np[w][u];
                    }
                }
            }
        }
        for(auto &v : Nnei[u]) if(inQv[v]==0 && v_sta[v]==0) {
            for(auto &w : Pnei[v]) if(inQv[w]==0 && v_sta[w]==0) {
                if(v<w && M[u][w] == -1) {
                    ++tri_nn[w][v]; ++tri_nn[v][w];
                    if(part[u]==1) {
                        ++tri_np[v][u]; ++tri_np[u][v];
                        ++tri_np[u][w]; ++tri_np[w][u];
                    }
                    else {
                        assert(part[u]==2);
                        ++tri_pn[v][u]; ++tri_pn[u][v];
                        ++tri_pn[u][w]; ++tri_pn[w][u];
                    }
                }
            }
            for(auto &w : Nnei[v]) if(inQv[w]==0 && v_sta[w]==0) {
                if(v<w && M[u][w] == 1) {
                    ++tri_nn[u][w]; ++tri_nn[w][u];
                    if(part[u]==1) {
                        ++tri_pn[v][w]; ++tri_pn[w][v];
                        ++tri_pn[u][v]; ++tri_pn[v][u];
                    }
                    else {
                        assert(part[u]==2);
                        ++tri_np[v][w]; ++tri_np[w][v];
                        ++tri_np[u][v]; ++tri_np[v][u];
                    }
                }
            }
        }
        v_sta[u] = 1;
    }
#ifdef _DEBUG_
    cout<<"** edge support information : "<<endl;
    for(auto &u : S) {
        cout<<"for vertex "<<u<<endl;
        for(auto &v : Pnei[u]) if(inQv[v]==0) {
            cout<<"tri_pp["<<u<<"]["<<v<<"] = "<<tri_pp[u][v]<<endl;
            cout<<"tri_nn["<<u<<"]["<<v<<"] = "<<tri_nn[u][v]<<endl;
        }
        for(auto &v : Nnei[u]) if(inQv[v]==0) {
            cout<<"tri_pn["<<u<<"]["<<v<<"] = "<<tri_pn[u][v]<<endl;
            cout<<"tri_np["<<u<<"]["<<v<<"] = "<<tri_np[u][v]<<endl;
        }
    }
#endif
    //initialize Qe
//    vector<pair<ui, ui>> rmpe,rmne;
    queue<pair<ui, ui>> Qe;
    for(auto &u : S) {
        for(auto &v : Pnei[u]) if(inQv[v]==0 && v>u) {
            if(tri_pp[u][v]<tau-2*k || tri_nn[u][v]<tau-2*k+2 || tri_pp[u][v]+tri_nn[u][v]<2*tau-2*k) {
                Qe.push(make_pair(u, v));
                inQe[u][v] = 1;
                inQe[v][u] = 1;
            }
        }
        for(auto &v : Nnei[u]) if(inQv[v]==0 && v>u) {
            if(tri_pn[u][v]<tau-2*k+1 || tri_np[u][v]<tau-2*k+1 || tri_pn[u][v]+tri_np[u][v]<2*tau-2*k){
                Qe.push(make_pair(u, v));
                inQe[u][v] = 1;
                inQe[v][u] = 1;
            }
        }
    }
    while (!Qe.empty()) {
        ui u = Qe.front().first;
        ui v = Qe.front().second;
        Qe.pop();
        assert(inQv[u]==0 && inQv[v]==0);
        if(M[u][v] == 1) {
            assert(p_degree[u]>0 && p_degree[v]>0);
            --p_degree[u]; --p_degree[v];
            rmpe.push_back(make_pair(u, v));
            for(auto &w : Pnei[u]) if(inQv[w]==0 && M[v][w]==1){
                assert(tri_pp[u][w]>0 && tri_pp[w][u]>0); --tri_pp[u][w]; --tri_pp[w][u];
                if(!inQe[u][w] && (tri_pp[u][w]<tau-2*k || tri_nn[u][w]<tau-2*k+2 || tri_pp[u][w]+tri_nn[u][w]<2*tau-2*k) ) {
                    Qe.push(make_pair(u, w));
                    inQe[u][w] = 1; inQe[w][u] = 1;
                }
                assert(tri_pp[v][w]>0 &&tri_pp[w][v]>0); --tri_pp[v][w]; --tri_pp[w][v];
                if(!inQe[v][w] && (tri_pp[v][w]<tau-2*k || tri_nn[v][w]<tau-2*k+2 || tri_pp[v][w]+tri_nn[v][w]<2*tau-2*k) ) {
                    Qe.push(make_pair(v, w));
                    inQe[v][w] = 1; inQe[w][v] = 1;
                }
            }
            for(auto &w : Nnei[u]) if(inQv[w]==0 && M[v][w]==-1){
                if(part[u] == 1) {
                    assert(tri_pn[u][w]>0 && tri_pn[w][u]>0); --tri_pn[u][w]; --tri_pn[w][u];
                    assert(tri_pn[v][w]>0 && tri_pn[w][v]>0); --tri_pn[v][w]; --tri_pn[w][v];
                }
                else {
                    assert(part[u] == 2);
                    assert(tri_np[u][w]>0 && tri_np[w][u]>0); --tri_np[u][w]; --tri_np[w][u];
                    assert(tri_np[v][w]>0 && tri_np[w][v]>0); --tri_np[v][w]; --tri_np[w][v];
                }
                if(!inQe[u][w] && (tri_pn[u][w]<tau-2*k+1 || tri_np[u][w]<tau-2*k+1 || tri_pn[u][w]+tri_np[u][w]<2*tau-2*k) ) {
                    Qe.push(make_pair(u, w));
                    inQe[u][w] = 1; inQe[w][u] = 1;
                }
                if(!inQe[v][w] && (tri_pn[v][w]<tau-2*k+1 || tri_np[v][w]<tau-2*k+1 || tri_pn[v][w]+tri_np[v][w]<2*tau-2*k) ) {
                    Qe.push(make_pair(v, w));
                    inQe[v][w] = 1; inQe[w][v] = 1;
                }
            }
        }
        else {
            assert(M[u][v] == -1);
            assert(n_degree[u]>0 && n_degree[v]>0);
            --n_degree[u]; --n_degree[v];
            rmne.push_back(make_pair(u, v));
            for(auto &w : Pnei[u]) if(inQv[w]==0 && M[v][w]==-1){
                assert(tri_nn[u][w]>0 && tri_nn[w][u]>0); --tri_nn[u][w]; --tri_nn[w][u];
                if(!inQe[u][w] && (tri_pp[u][w]<tau-2*k || tri_nn[u][w]<tau-2*k+2 || tri_pp[u][w]+tri_nn[u][w]<2*tau-2*k) ) {
                    Qe.push(make_pair(u, w));
                    inQe[u][w] = 1; inQe[w][u] = 1;
                }
                if(part[u]==1) {
                    assert(tri_pn[v][w]>0 && tri_pn[w][v]>0); --tri_pn[v][w]; --tri_pn[w][v];
                }
                else {
                    assert(part[u]==2);
                    assert(tri_np[v][w]>0 && tri_np[w][v]>0); --tri_np[v][w]; --tri_np[w][v];
                }
                if(!inQe[v][w] && (tri_pn[w][v]<tau-2*k+1 || tri_np[w][v]<tau-2*k+1 || tri_pn[w][v]+tri_np[w][v]<2*tau-2*k) ) {
                    Qe.push(make_pair(v, w));
                    inQe[v][w] = 1; inQe[w][v] = 1;
                }
            }
            for(auto &w : Nnei[u]) if(inQv[w]==0 && M[v][w]==1){
                assert(tri_nn[v][w]>0 && tri_nn[w][v]>0); --tri_nn[v][w]; --tri_nn[w][v];
                if(!inQe[v][w] && (tri_pp[v][w]<tau-2*k || tri_nn[v][w]<tau-2*k+2 || tri_pp[v][w]+tri_nn[v][w]<2*tau-2*k) ) {
                    Qe.push(make_pair(v, w));
                    inQe[v][w] = 1; inQe[w][v] = 1;
                }
                if(part[u]==1) {
                    assert(tri_np[u][w]>0 && tri_np[w][u]>0); --tri_np[u][w]; --tri_np[w][u];
                }
                else {
                    assert(part[u]==2);
                    assert(tri_pn[u][w]>0 && tri_pn[w][u]>0); --tri_pn[u][w]; --tri_pn[w][u];
                }
                if(!inQe[u][w] && (tri_pn[u][w]<tau-2*k+1 || tri_np[u][w]<tau-2*k+1 || tri_pn[u][w]+tri_np[u][w]<2*tau-2*k) ) {
                    Qe.push(make_pair(u, w));
                    inQe[u][w] = 1; inQe[w][u] = 1;
                }
            }
        }
        M[u][v] = 0;
        M[v][u] = 0;
    }
    for(auto &e : rmpe) {
        assert(inQe[e.first][e.second] == 1 && inQe[e.second][e.first] == 1);
        inQe[e.first][e.second] = 0;
        inQe[e.second][e.first] = 0;
    }
    for(auto &e : rmne) {
        assert(inQe[e.first][e.second] == 1 && inQe[e.second][e.first] == 1);
        inQe[e.first][e.second] = 0;
        inQe[e.second][e.first] = 0;
    }
}

void refine_S_by_capacity(vector<ui> &t_S, vector<ui> &s_S, vector<ui> &KL, vector<ui> &KR, bool from_L)
{
    int sign;
    if(from_L == true) sign = 1;
    else sign = -1;
    for(auto &e : s_S) {
        bool flag = true;
        short tmp_capacity = k - 1;
        for(auto &x : KL) {
            if(M[e][x] != sign) {
                -- tmp_capacity;
                if(tmp_capacity < 0 || partial_kplex_cap[x] <= 0) {
                    flag = false;
                    break;
                }
            }
        }
        if(flag == false) continue;
        for(auto &x : KR) {
            if(M[e][x] != (-sign)) {
                -- tmp_capacity;
                if(tmp_capacity < 0 || partial_kplex_cap[x] <= 0) {
                    flag = false;
                    break;
                }
            }
        }
        if(flag == true) t_S.push_back(e);
    }
}

void refine_S_by_capacity_mdpivot(vector<ui> &t_S, vector<ui> &s_S, vector<ui> &KL, vector<ui> &KR, bool from_L)
{
    int sign;
    if(from_L == true) sign = 1;
    else sign = -1;
    for(auto &e : s_S) {
        p_degree[e] = 0 ; n_degree[e] = 0;
        bool flag = true;
        short tmp_capacity = k - 1;
        for(auto &x : KL) {
            if(M[e][x] != sign) {
                -- tmp_capacity;
                if(tmp_capacity < 0 || partial_kplex_cap[x] <= 0) {
                    flag = false;
                    break;
                }
            }
        }
        if(flag == false) continue;
        for(auto &x : KR) {
            if(M[e][x] != (-sign)) {
                -- tmp_capacity;
                if(tmp_capacity < 0 || partial_kplex_cap[x] <= 0) {
                    flag = false;
                    break;
                }
            }
        }
        if(flag == true) t_S.push_back(e);
    }
}

void refine_S_by_capacity(vector<ui> &t_S, vector<ui> &s_S, vector<ui> &KL, vector<ui> &KR, bool from_L, int collect_nnei_inK)
{
    int sign;
    if(from_L == true) sign = 1;
    else sign = -1;
    for(auto &e : s_S) {
        p_degree[e] = 0; n_degree[e] = 0;
        nneiKL[e].clear();
        nneiKR[e].clear();
        bool flag = true;
        short tmp_capacity = k - 1;
        for(auto &x : KL) {
            if(M[e][x] != sign) {
                nneiKL[e].push_back(x);
                if(flag){
                    -- tmp_capacity;
                    if(tmp_capacity < 0 || partial_kplex_cap[x] <= 0) {
                        flag = false;
                    }
                }
            }
        }
        for(auto &x : KR) {
            if(M[e][x] != (-sign)) {
                nneiKR[e].push_back(x);
                if(flag){
                    -- tmp_capacity;
                    if(tmp_capacity < 0 || partial_kplex_cap[x] <= 0) {
                        flag = false;
                    }
                }
            }
        }
        if(flag == true) t_S.push_back(e);
    }
    
#ifdef _DEBUG_
    for(auto &e : s_S) {
        cout<<e<<" : nnei KL : "; for(auto &x : nneiKL[e]) cout<<x<<","; cout<<endl;
        cout<<e<<" : nnei KR : "; for(auto &x : nneiKR[e]) cout<<x<<","; cout<<endl;
    }
#endif
}

bool early_termination(vector<ui> &KL, vector<ui> &KR, vector<ui> &CL, vector<ui> &CR, vector<ui> &XL, vector<ui> &XR)
{
    bool canbe_added = true;
    for(auto u : XL) {
        for(auto e : KL) {
            if(M[u][e] != 1) {
                canbe_added = false; break;
            }
        }
        if(canbe_added == false) continue;
        for(auto e : KR) {
            if(M[u][e] != -1) {
                canbe_added = false; break;
            }
        }
        if(canbe_added == false) continue;
        for(auto e : CL) {
            if(M[u][e] != 1) {
                canbe_added = false; break;
            }
        }
        if(canbe_added == false) continue;
        for(auto e : CR) {
            if(e == u) {
//                cout<<"overlap between XL & CR"<<endl;
                continue;
            }
            if(M[u][e] != -1) {
                canbe_added = false; break;
            }
        }
        if(canbe_added) return true;
    }
    canbe_added = true;
    for(auto u : XR) {
        for(auto e : KL) {
            if(M[u][e] != -1) {
                canbe_added = false; break;
            }
        }
        if(canbe_added == false) continue;
        for(auto e : KR) {
            if(M[u][e] != 1) {
                canbe_added = false; break;
            }
        }
        if(canbe_added == false) continue;
        for(auto e : CL) {
            if(e == u) {
//                cout<<"overlap between CL & XR"<<endl;
                continue;
            }
            if(M[u][e] != -1) {
                canbe_added = false; break;
            }
        }
        if(canbe_added == false) continue;
        for(auto e : CR) {
            if(M[u][e] != 1) {
                canbe_added = false; break;
            }
        }
        if(canbe_added) return true;
    }
    return false;
}

bool early_termination_adv_current_wrong(vector<ui> &KL, vector<ui> &KR, vector<ui> &CL, vector<ui> &CR, vector<ui> &XL, vector<ui> &XR)
{
    bool canbe_added = true;
    for(auto u : XL) {
        ui u_missed_num = 1;
        for(auto v : KL) {
            if(M[u][v] != 1) {
                ++ u_missed_num;
                if(u_missed_num > k) {canbe_added = false; break;}
                //u misses v \in KL, check if v is ample
                assert(partial_kplex_cap[v] >= 0 && partial_kplex_cap[v] <= k - 1);
                ui v_missed_num_in_KC = k - partial_kplex_cap[v];
                for(auto w : CL) if (M[v][w] != 1) {
                    ++ v_missed_num_in_KC;
                    if(v_missed_num_in_KC >= k ) {
                        canbe_added = false; break;
                    }
                }
                if(canbe_added == false) break;
                for(auto w : CR) if (M[v][w] != -1) {
                    ++ v_missed_num_in_KC;
                    if(v_missed_num_in_KC >= k ) {
                        canbe_added = false; break;
                    }
                }
                if(canbe_added == false) break;
            }
        }
        if(canbe_added == false) continue;
        for(auto v : KR) {
            if(M[u][v] != -1) {
                ++ u_missed_num;
                if(u_missed_num > k) {canbe_added = false; break;}
                //u misses v \in KR, check if v is ample
                assert(partial_kplex_cap[v] >= 0 && partial_kplex_cap[v] <= k - 1);
                ui v_missed_num_in_KC = k - partial_kplex_cap[v];
                for(auto w : CL) if (M[v][w] != -1) {
                    ++ v_missed_num_in_KC;
                    if(v_missed_num_in_KC >= k ) {
                        canbe_added = false; break;
                    }
                }
                if(canbe_added == false) break;
                for(auto w : CR) if (M[v][w] != 1) {
                    ++ v_missed_num_in_KC;
                    if(v_missed_num_in_KC >= k ) {
                        canbe_added = false; break;
                    }
                }
                if(canbe_added == false) break;
            }
        }
        if(canbe_added == false) continue;
        for(auto v : CL) {
            if(M[u][v] != 1) {
                ++ u_missed_num;
                if(u_missed_num > k) {canbe_added = false; break;}
                ui v_missed_num_in_KC = 1;
                for(auto w : KL) if(M[v][w] != 1) {
                    ++ v_missed_num_in_KC;
                    if(v_missed_num_in_KC >= k) {
                        canbe_added = false; break;
                    }
                }
                if(canbe_added == false) break;
                for(auto w : KR) if(M[v][w] != -1) {
                    ++ v_missed_num_in_KC;
                    if(v_missed_num_in_KC >= k) {
                        canbe_added = false; break;
                    }
                }
                if(canbe_added == false) break;
                for(auto w : CL) if( w != v && M[v][w] != 1) {
                    ++ v_missed_num_in_KC;
                    if(v_missed_num_in_KC >= k) {
                        canbe_added = false; break;
                    }
                }
                if(canbe_added == false) break;
                for(auto w : CR) if( w != v && M[v][w] != -1) {
                    ++ v_missed_num_in_KC;
                    if(v_missed_num_in_KC >= k) {
                        canbe_added = false; break;
                    }
                }
                if(canbe_added == false) break;
            }
        }
        if(canbe_added == false) continue;
        for(auto v : CR) {
            if(u == v) continue;
            if(M[u][v] != -1) {
                ++ u_missed_num;
                if(u_missed_num > k) {canbe_added = false; break;}
                ui v_missed_num_in_KC = 1;
                for(auto w : KL) if(M[v][w] != -1) {
                    ++ v_missed_num_in_KC;
                    if(v_missed_num_in_KC >= k) {
                        canbe_added = false; break;
                    }
                }
                if(canbe_added == false) break;
                for(auto w : KR) if(M[v][w] != 1) {
                    ++ v_missed_num_in_KC;
                    if(v_missed_num_in_KC >= k) {
                        canbe_added = false; break;
                    }
                }
                if(canbe_added == false) break;
                for(auto w : CL) if( w != v && M[v][w] != -1) {
                    ++ v_missed_num_in_KC;
                    if(v_missed_num_in_KC >= k) {
                        canbe_added = false; break;
                    }
                }
                if(canbe_added == false) break;
                for(auto w : CR) if( w != v && M[v][w] != 1) {
                    ++ v_missed_num_in_KC;
                    if(v_missed_num_in_KC >= k) {
                        canbe_added = false; break;
                    }
                }
                if(canbe_added == false) break;
            }
        }
        if(canbe_added) return true;
    }
    canbe_added = true;
    for(auto u : XR) {
        ui u_missed_num = 1;
        for(auto v : KL) {
            if(M[u][v] != -1) {
                ++ u_missed_num;
                if(u_missed_num > k) {canbe_added = false; break;}
                //u misses v \in KL, check if v is ample
                assert(partial_kplex_cap[v] >= 0 && partial_kplex_cap[v] <= k - 1);
                ui v_missed_num_in_KC = k - partial_kplex_cap[v];
                for(auto w : CL) if (M[v][w] != 1) {
                    ++ v_missed_num_in_KC;
                    if(v_missed_num_in_KC >= k ) {
                        canbe_added = false; break;
                    }
                }
                if(canbe_added == false) break;
                for(auto w : CR) if (M[v][w] != -1) {
                    ++ v_missed_num_in_KC;
                    if(v_missed_num_in_KC >= k ) {
                        canbe_added = false; break;
                    }
                }
                if(canbe_added == false) break;
            }
        }
        if(canbe_added == false) continue;
        for(auto v : KR) {
            if(M[u][v] != 1) {
                ++ u_missed_num;
                if(u_missed_num > k) {canbe_added = false; break;}
                //u misses v \in KR, check if v is ample
                assert(partial_kplex_cap[v] >= 0 && partial_kplex_cap[v] <= k - 1);
                ui v_missed_num_in_KC = k - partial_kplex_cap[v];
                for(auto w : CL) if (M[v][w] != -1) {
                    ++ v_missed_num_in_KC;
                    if(v_missed_num_in_KC >= k ) {
                        canbe_added = false; break;
                    }
                }
                if(canbe_added == false) break;
                for(auto w : CR) if (M[v][w] != 1) {
                    ++ v_missed_num_in_KC;
                    if(v_missed_num_in_KC >= k ) {
                        canbe_added = false; break;
                    }
                }
                if(canbe_added == false) break;
            }
        }
        if(canbe_added == false) continue;
        for(auto v : CL) {
            if(u == v) continue;
            if(M[u][v] != -1) {
                ++ u_missed_num;
                if(u_missed_num > k) {canbe_added = false; break;}
                ui v_missed_num_in_KC = 1;
                for(auto w : KL) if(M[v][w] != 1) {
                    ++ v_missed_num_in_KC;
                    if(v_missed_num_in_KC >= k) {
                        canbe_added = false; break;
                    }
                }
                if(canbe_added == false) break;
                for(auto w : KR) if(M[v][w] != -1) {
                    ++ v_missed_num_in_KC;
                    if(v_missed_num_in_KC >= k) {
                        canbe_added = false; break;
                    }
                }
                if(canbe_added == false) break;
                for(auto w : CL) if( w != v && M[v][w] != 1) {
                    ++ v_missed_num_in_KC;
                    if(v_missed_num_in_KC >= k) {
                        canbe_added = false; break;
                    }
                }
                if(canbe_added == false) break;
                for(auto w : CR) if( w != v && M[v][w] != -1) {
                    ++ v_missed_num_in_KC;
                    if(v_missed_num_in_KC >= k) {
                        canbe_added = false; break;
                    }
                }
                if(canbe_added == false) break;
            }
        }
        if(canbe_added == false) continue;
        for(auto v : CR) {
            if(M[u][v] != 1) {
                ++ u_missed_num;
                if(u_missed_num > k) {canbe_added = false; break;}
                ui v_missed_num_in_KC = 1;
                for(auto w : KL) if(M[v][w] != -1) {
                    ++ v_missed_num_in_KC;
                    if(v_missed_num_in_KC >= k) {
                        canbe_added = false; break;
                    }
                }
                if(canbe_added == false) break;
                for(auto w : KR) if(M[v][w] != 1) {
                    ++ v_missed_num_in_KC;
                    if(v_missed_num_in_KC >= k) {
                        canbe_added = false; break;
                    }
                }
                if(canbe_added == false) break;
                for(auto w : CL) if( w != v && M[v][w] != -1) {
                    ++ v_missed_num_in_KC;
                    if(v_missed_num_in_KC >= k) {
                        canbe_added = false; break;
                    }
                }
                if(canbe_added == false) break;
                for(auto w : CR) if( w != v && M[v][w] != 1) {
                    ++ v_missed_num_in_KC;
                    if(v_missed_num_in_KC >= k) {
                        canbe_added = false; break;
                    }
                }
                if(canbe_added == false) break;
            }
        }
        if(canbe_added) return true;
    }
    return false;
}

bool early_termination_adv(vector<ui> &KL, vector<ui> &KR, vector<ui> &CL, vector<ui> &CR, vector<ui> &XL, vector<ui> &XR)
{
#ifdef _DEBUG_
    cout<<endl<<"in ET: "<<endl;
    cout<<"KL ("<<KL.size()<<"): "; for(auto e : KL) cout<<e<<","; cout<<endl;
    cout<<"KR ("<<KR.size()<<"): "; for(auto e : KR) cout<<e<<","; cout<<endl;
    cout<<"CL ("<<CL.size()<<"): "; for(auto e : CL) cout<<e<<","; cout<<endl;
    cout<<"CR ("<<CR.size()<<"): "; for(auto e : CR) cout<<e<<","; cout<<endl;
    cout<<"XL ("<<XL.size()<<"): "; for(auto e : XL) cout<<e<<","; cout<<endl;
    cout<<"XR ("<<XR.size()<<"): "; for(auto e : XR) cout<<e<<","; cout<<endl;
#endif
    int KCsize = (int)KL.size() + (int)KR.size() + (int)CL.size() + (int)CR.size();
    bool canbe_added = true;
    for(auto &u : XL) {
        ui u_missed_num = 1;
        for(auto &v : KL) {
            if(M[u][v] != 1) {
                ++ u_missed_num;
                if(u_missed_num > k) {canbe_added = false; break;}
                assert(KCsize > (p_degree[v]+n_degree[v]));
                int v_missed_num = KCsize - (p_degree[v]+n_degree[v]);
                if(v_missed_num >= k) {canbe_added = false; break;}
            }
        }
        if(canbe_added == false) continue;
        for(auto &v : KR) {
            if(M[u][v] != -1) {
                ++ u_missed_num;
                if(u_missed_num > k) {canbe_added = false; break;}
                assert(KCsize > (p_degree[v]+n_degree[v]));
                int v_missed_num = KCsize - (p_degree[v]+n_degree[v]);
                if(v_missed_num >= k) {canbe_added = false; break;}
            }
        }
        if(canbe_added == false) continue;
        for(auto &v : CL) {
            if(M[u][v] != 1) {
                ++ u_missed_num;
                if(u_missed_num > k) {canbe_added = false; break;}
                assert(KCsize > (p_degree[v]+n_degree[v]));
                int v_missed_num_in_KC = KCsize - (p_degree[v]+n_degree[v]);
                if(v_missed_num_in_KC >= k) {canbe_added = false; break;}
            }
        }
        if(canbe_added == false) continue;
        for(auto &v : CR) {
            if(u == v) continue;
            if(M[u][v] != -1) {
                ++ u_missed_num;
                if(u_missed_num > k) {canbe_added = false; break;}
                assert(KCsize > (p_degree[v]+n_degree[v]));
                int v_missed_num_in_KC = KCsize - (p_degree[v]+n_degree[v]);
                if(v_missed_num_in_KC >= k) {canbe_added = false; break;}
            }
        }
        if(canbe_added) return true;
    }
    canbe_added = true;
    for(auto &u : XR) {
        ui u_missed_num = 1;
        for(auto v : KL) {
            if(M[u][v] != -1) {
                ++ u_missed_num;
                if(u_missed_num > k) {canbe_added = false; break;}
                assert(KCsize > (p_degree[v]+n_degree[v]));
                int v_missed_num = KCsize - (p_degree[v]+n_degree[v]);
                if(v_missed_num >= k) {canbe_added = false; break;}
            }
        }
        if(canbe_added == false) continue;
        for(auto &v : KR) {
            if(M[u][v] != 1) {
                ++ u_missed_num;
                if(u_missed_num > k) {canbe_added = false; break;}
                assert(KCsize > (p_degree[v]+n_degree[v]));
                int v_missed_num = KCsize - (p_degree[v]+n_degree[v]);
                if(v_missed_num >= k) {canbe_added = false; break;}
            }
        }
        if(canbe_added == false) continue;
        for(auto &v : CL) {
            if(u == v) continue;
            if(M[u][v] != -1) {
                ++ u_missed_num;
                if(u_missed_num > k) {canbe_added = false; break;}
                assert(KCsize > (p_degree[v]+n_degree[v]));
                int v_missed_num_in_KC = KCsize - (p_degree[v]+n_degree[v]);
                if(v_missed_num_in_KC >= k) {canbe_added = false; break;}
            }
        }
        if(canbe_added == false) continue;
        for(auto &v : CR) {
            if(M[u][v] != 1) {
                ++ u_missed_num;
                if(u_missed_num > k) {canbe_added = false; break;}
                assert(KCsize > (p_degree[v]+n_degree[v]));
                int v_missed_num_in_KC = KCsize - (p_degree[v]+n_degree[v]);
                if(v_missed_num_in_KC >= k) {canbe_added = false; break;}
            }
        }
        if(canbe_added) return true;
    }
    return false;
}

void get_skipvlist_by_pivoting(vector<ui> &KL, vector<ui> &KR, vector<ui> &CL, vector<ui> &CR, pair<vector<ui>, vector<ui>> & svs)
{
    assert(!CL.empty() || !CR.empty());
    vector<vector<ui>> CL_skipvlist_fromCL, CL_skipvlist_fromCR;
    vector<vector<ui>> CR_skipvlist_fromCL, CR_skipvlist_fromCR;
    int power = -1;
    for(auto u : CL) {
        vector<ui> tmp_skipvlist_fromCL, tmp_skipvlist_fromCR;
        for(ui i = 0; i < CL.size(); i++) if(M[u][CL[i]] == 1) {
            ui v = CL[i];
            bool cnn = false;
            for(auto &w : nneiKL[u]) {
                if(M[v][w]!=1) {
                    cnn = true; break;
                }
            }
            if(cnn) continue;
            for(auto &w : nneiKR[u]) {
                if(M[v][w]!=-1) {
                    cnn = true; break;
                }
            }
            if(cnn==false) tmp_skipvlist_fromCL.push_back(i);
        }
        for(ui i = 0; i < CR.size(); i++) if(M[u][CR[i]] == -1) {
            ui v = CR[i];
            bool cnn = false;
            for(auto &w : nneiKL[u]) {
                if(M[v][w]!=-1) {
                    cnn = true; break;
                }
            }
            if(cnn) continue;
            for(auto &w : nneiKR[u]) {
                if(M[v][w]!=1) {
                    cnn = true; break;
                }
            }
            if(cnn==false) tmp_skipvlist_fromCR.push_back(i);
        }
        int tmp_power = (int)tmp_skipvlist_fromCL.size()+(int)tmp_skipvlist_fromCR.size();
        if( tmp_power > power) {
            power = tmp_power;
            svs = make_pair(tmp_skipvlist_fromCL, tmp_skipvlist_fromCR);
        }
    }
    for(auto u : CR) {
        vector<ui> tmp_skipvlist_fromCL, tmp_skipvlist_fromCR;
        for(ui i = 0; i < CL.size(); i++) if(M[u][CL[i]] == -1) {
            ui v = CL[i];
            bool cnn = false;
            for(auto &w : nneiKL[u]) {
                if(M[v][w]!=1) {
                    cnn = true; break;
                }
            }
            if(cnn) continue;
            for(auto &w : nneiKR[u]) {
                if(M[v][w]!=-1) {
                    cnn = true; break;
                }
            }
            if(cnn==false) tmp_skipvlist_fromCL.push_back(i);
        }
        for(ui i = 0; i < CR.size(); i++) if(M[u][CR[i]] == 1) {
            ui v = CR[i];
            bool cnn = false;
            for(auto &w : nneiKL[u]) {
                if(M[v][w]!=-1) {
                    cnn = true; break;
                }
            }
            if(cnn) continue;
            for(auto &w : nneiKR[u]) {
                if(M[v][w]!=1) {
                    cnn = true; break;
                }
            }
            if(cnn==false) tmp_skipvlist_fromCR.push_back(i);
        }
        int tmp_power = (int)tmp_skipvlist_fromCL.size()+(int)tmp_skipvlist_fromCR.size();
        if( tmp_power > power) {
            power = tmp_power;
            svs = make_pair(tmp_skipvlist_fromCL, tmp_skipvlist_fromCR);
        }
    }
    assert(power != -1);
}

void get_skipvlist_by_pivoting_light(vector<ui> &KL, vector<ui> &KR, vector<ui> &CL, vector<ui> &CR, pair<vector<ui>, vector<ui>> & svs)
{
    assert(!CL.empty() || !CR.empty());
//    vector<ui> CL_power, CR_power;
    vector<vector<ui>> CL_skipvlist_fromCL, CL_skipvlist_fromCR;//?
    vector<vector<ui>> CR_skipvlist_fromCL, CR_skipvlist_fromCR;//?
    for(auto u : CL) {
        vector<ui> tmp_skipvlist_fromCL, tmp_skipvlist_fromCR;
        for(auto v : CL) if(M[u][v] == 1) {
            bool cnn = false;
            for(auto w : KL) if(M[u][w] != 1 && M[v][w] != 1) {
                cnn = true; break;
            }
            if(cnn) continue;
            for(auto w : KR) if(M[u][w] != -1 && M[v][w] != -1) {
                cnn = true; break;
            }
            if(cnn == false) tmp_skipvlist_fromCL.push_back(v);
        }
        for(auto v : CR) if(M[u][v] == -1) {
            bool cnn = false;
            for(auto w : KL) if(M[u][w] != 1 && M[v][w] != -1) {
                cnn = true; break;
            }
            if(cnn) continue;
            for(auto w : KR) if(M[u][w] != -1 && M[v][w] != 1) {
                cnn = true; break;
            }
            if(cnn == false) tmp_skipvlist_fromCR.push_back(v);
        }
        CL_skipvlist_fromCL.push_back(tmp_skipvlist_fromCL);
        CL_skipvlist_fromCR.push_back(tmp_skipvlist_fromCR);
    }
    for(auto u : CR) {
        vector<ui> tmp_skipvlist_fromCL, tmp_skipvlist_fromCR;
        for(auto v : CL) if(M[u][v] == -1) {
            bool cnn = false;
            for(auto w : KL) if(M[u][w] != -1 && M[v][w] != 1) {
                cnn = true; break;
            }
            if(cnn) continue;
            for(auto w : KR) if(M[u][w] != 1 && M[v][w] != -1) {
                cnn = true; break;
            }
            if(cnn == false) tmp_skipvlist_fromCL.push_back(v);
        }
        for(auto v : CR) if(M[u][v] == 1) {
            bool cnn = false;
            for(auto w : KL) if(M[u][w] != -1 && M[v][w] != -1) {
                cnn = true; break;
            }
            if(cnn) continue;
            for(auto w : KR) if(M[u][w] != 1 && M[v][w] != 1) {
                cnn = true; break;
            }
            if(cnn == false) tmp_skipvlist_fromCR.push_back(v);
        }
        CR_skipvlist_fromCL.push_back(tmp_skipvlist_fromCL);
        CR_skipvlist_fromCR.push_back(tmp_skipvlist_fromCR);
    }
    bool from_CL = true;
    ui ustar_idx = INF;
    int power = -1;
    for(ui i = 0; i < CL.size(); i++) {
        int tmp_power = (int)CL_skipvlist_fromCL[i].size() + (int)CL_skipvlist_fromCR[i].size();
        if(tmp_power > power) {
            ustar_idx = i;
            power = tmp_power;
        }
    }
    for(ui i = 0; i < CR.size(); i++) {
        int tmp_power = (int)CR_skipvlist_fromCL[i].size() + (int)CR_skipvlist_fromCR[i].size();
        if(tmp_power > power) {
            from_CL = false;
            ustar_idx = i;
            power = tmp_power;
        }
    }
    assert(power != -1 && ustar_idx != INF);
    if(from_CL) svs = make_pair(CL_skipvlist_fromCL[ustar_idx], CL_skipvlist_fromCR[ustar_idx]);
    else svs = make_pair(CR_skipvlist_fromCL[ustar_idx], CR_skipvlist_fromCR[ustar_idx]);
#ifdef _DEBUG_
    cout<<"pivot information:"<<endl;
    cout<<"[CL]:"<<endl;
    for(ui i = 0; i < CL.size(); i++) {
        cout<<CL[i]<<endl;
        cout<<"\tCL:"; for(auto e : CL_skipvlist_fromCL[i]) cout<<e<<","; cout<<endl;
        cout<<"\tCR:"; for(auto e : CL_skipvlist_fromCR[i]) cout<<e<<","; cout<<endl;
    }
    cout<<"[CR]:"<<endl;
    for(ui i = 0; i < CR.size(); i++) {
        cout<<CR[i]<<endl;
        cout<<"\tCL:"; for(auto e : CR_skipvlist_fromCL[i]) cout<<e<<","; cout<<endl;
        cout<<"\tCR:"; for(auto e : CR_skipvlist_fromCR[i]) cout<<e<<","; cout<<endl;
    }
    if(from_CL) cout<<"pivot = "<<CL[ustar_idx]<<", from CL."<<endl;
    else cout<<"pivot = "<<CR[ustar_idx]<<", from CR."<<endl;
    cout<<"svs : "; for(auto e : svs.first) cout<<e<<","; cout<<" | "; for(auto e : svs.second) cout<<e<<","; cout<<endl;
#endif
    
}

void get_skipvlist_by_pivoting_original(vector<ui> &KL, vector<ui> &KR, vector<ui> &CL, vector<ui> &CR, pair<vector<ui>, vector<ui>> & svs)
{
    assert(!CL.empty() || !CR.empty());
//    vector<ui> CL_power, CR_power;
    vector<vector<ui>> CL_skipvlist_fromCL, CL_skipvlist_fromCR;//?
    vector<vector<ui>> CR_skipvlist_fromCL, CR_skipvlist_fromCR;//?
    for(auto u : CL) {
        vector<ui> tmp_skipvlist_fromCL, tmp_skipvlist_fromCR;
        for(auto v : CL) if(M[u][v] == 1) {
            bool cnn = false;
            for(auto w : KL) if(M[u][w] != 1 && M[v][w] != 1) {
                cnn = true; break;
            }
            if(cnn) continue;
            for(auto w : KR) if(M[u][w] != -1 && M[v][w] != -1) {
                cnn = true; break;
            }
            if(cnn == false) tmp_skipvlist_fromCL.push_back(v);
        }
        for(auto v : CR) if(M[u][v] == -1) {
            bool cnn = false;
            for(auto w : KL) if(M[u][w] != 1 && M[v][w] != -1) {
                cnn = true; break;
            }
            if(cnn) continue;
            for(auto w : KR) if(M[u][w] != -1 && M[v][w] != 1) {
                cnn = true; break;
            }
            if(cnn == false) tmp_skipvlist_fromCR.push_back(v);
        }
        CL_skipvlist_fromCL.push_back(tmp_skipvlist_fromCL);
        CL_skipvlist_fromCR.push_back(tmp_skipvlist_fromCR);
    }
    for(auto u : CR) {
        vector<ui> tmp_skipvlist_fromCL, tmp_skipvlist_fromCR;
        for(auto v : CL) if(M[u][v] == -1) {
            bool cnn = false;
            for(auto w : KL) if(M[u][w] != -1 && M[v][w] != 1) {
                cnn = true; break;
            }
            if(cnn) continue;
            for(auto w : KR) if(M[u][w] != 1 && M[v][w] != -1) {
                cnn = true; break;
            }
            if(cnn == false) tmp_skipvlist_fromCL.push_back(v);
        }
        for(auto v : CR) if(M[u][v] == 1) {
            bool cnn = false;
            for(auto w : KL) if(M[u][w] != -1 && M[v][w] != -1) {
                cnn = true; break;
            }
            if(cnn) continue;
            for(auto w : KR) if(M[u][w] != 1 && M[v][w] != 1) {
                cnn = true; break;
            }
            if(cnn == false) tmp_skipvlist_fromCR.push_back(v);
        }
        CR_skipvlist_fromCL.push_back(tmp_skipvlist_fromCL);
        CR_skipvlist_fromCR.push_back(tmp_skipvlist_fromCR);
    }
    bool from_CL = true;
    ui ustar_idx = INF;
    int power = -1;
    for(ui i = 0; i < CL.size(); i++) {
        int tmp_power = (int)CL_skipvlist_fromCL[i].size() + (int)CL_skipvlist_fromCR[i].size();
        if(tmp_power > power) {
            ustar_idx = i;
            power = tmp_power;
        }
    }
    for(ui i = 0; i < CR.size(); i++) {
        int tmp_power = (int)CR_skipvlist_fromCL[i].size() + (int)CR_skipvlist_fromCR[i].size();
        if(tmp_power > power) {
            from_CL = false;
            ustar_idx = i;
            power = tmp_power;
        }
    }
    assert(power != -1 && ustar_idx != INF);
    if(from_CL) svs = make_pair(CL_skipvlist_fromCL[ustar_idx], CL_skipvlist_fromCR[ustar_idx]);
    else svs = make_pair(CR_skipvlist_fromCL[ustar_idx], CR_skipvlist_fromCR[ustar_idx]);
#ifdef _DEBUG_
    cout<<"pivot information:"<<endl;
    cout<<"[CL]:"<<endl;
    for(ui i = 0; i < CL.size(); i++) {
        cout<<CL[i]<<endl;
        cout<<"\tCL:"; for(auto e : CL_skipvlist_fromCL[i]) cout<<e<<","; cout<<endl;
        cout<<"\tCR:"; for(auto e : CL_skipvlist_fromCR[i]) cout<<e<<","; cout<<endl;
    }
    cout<<"[CR]:"<<endl;
    for(ui i = 0; i < CR.size(); i++) {
        cout<<CR[i]<<endl;
        cout<<"\tCL:"; for(auto e : CR_skipvlist_fromCL[i]) cout<<e<<","; cout<<endl;
        cout<<"\tCR:"; for(auto e : CR_skipvlist_fromCR[i]) cout<<e<<","; cout<<endl;
    }
    if(from_CL) cout<<"pivot = "<<CL[ustar_idx]<<", from CL."<<endl;
    else cout<<"pivot = "<<CR[ustar_idx]<<", from CR."<<endl;
    cout<<"svs : "; for(auto e : svs.first) cout<<e<<","; cout<<" | "; for(auto e : svs.second) cout<<e<<","; cout<<endl;
#endif
    
}

void enum_procedure(vector<ui> &KL, vector<ui> &KR, vector<ui> &CL, vector<ui> &CR, vector<ui> &XL, vector<ui> &XR, int depth)
{
    ++ branch_num;
    obtain_partial_kplex_capacity(KL, KR); //compute vertex degree in partial kplex <KL,KR>
    vector<ui> new_CL; refine_S_by_capacity(new_CL, CL, KL, KR, true, 1); //refine CL
    vector<ui> new_CR; refine_S_by_capacity(new_CR, CR, KL, KR, false, 1); //refine CR
    vector<ui> new_XL; refine_S_by_capacity(new_XL, XL, KL, KR, true); //refine XL
    vector<ui> new_XR; refine_S_by_capacity(new_XR, XR, KL, KR, false); //refine XR
    if(KL.size() + new_CL.size() < tau || KR.size() + new_CR.size() < tau) return;
#ifdef _ET_
    if(early_termination(KL,KR,new_CL,new_CR,new_XL,new_XR)) return;
    if(early_termination_adv(KL,KR,new_CL,new_CR,new_XL,new_XR)) return;
#endif
    if(new_CL.empty() && new_CR.empty() && new_XL.empty() && new_XR.empty()) {
#ifdef _StoreRes_
        vector<ui> k1, k2;
        for(auto e : KL) k1.push_back(original_id[e]);
        for(auto e : KR) k2.push_back(original_id[e]);
        Res.push_back(make_pair(k1, k2));
#endif
        int kplesize = (int)KL.size() + (int)KR.size();
        if( kplesize > max_kplex) max_kplex = kplesize ;
        if( kplesize < min_kplex) min_kplex = kplesize;
        ++ Res_num;
    }
    if(new_CL.empty() && new_CR.empty()) return;
#ifdef _VRinEnum_
    if(VR_in_Enum_baseline(KL,KR,new_CL,new_CR) || (new_CL.empty()&&new_CR.empty())) return;
#endif
    vector<ui> CL_unexpanded_vs, CR_unexpanded_vs;
    vector<bool> CL_flag(new_CL.size(),0), CR_flag(new_CR.size(),0);
#ifdef _PIVOT_
    pair<vector<ui>, vector<ui>> skip_vs;
    get_skipvlist_by_pivoting_original(KL,KR,new_CL,new_CR,skip_vs);
    for(auto e : skip_vs.first) skipv_CL[e] = 1;
    for(auto e : skip_vs.second) skipv_CR[e] = 1;
    for(ui i = 0; i < new_CL.size(); i++) if(skipv_CL[new_CL[i]] == 1) CL_flag[i] = 1;
    for(ui i = 0; i < new_CR.size(); i++) if(skipv_CR[new_CR[i]] == 1) CR_flag[i] = 1;
    for(auto e : skip_vs.first) skipv_CL[e] = 0;
    for(auto e : skip_vs.second) skipv_CR[e] = 0;
#endif
    for(ui i = 0; i < new_CL.size(); i++) {
        ui u = new_CL[i];
        if(CL_flag[i] == 1) { CL_unexpanded_vs.push_back(u); continue; }
        vector<ui> new_CL_prime;
        for(ui j = i + 1; j < new_CL.size(); j++) new_CL_prime.push_back(new_CL[j]);
        new_CL_prime.insert(new_CL_prime.end(), CL_unexpanded_vs.begin(), CL_unexpanded_vs.end());
        vector<ui> new_CR_prime;
        for(auto e : new_CR) if (e != u) new_CR_prime.push_back(e);
        vector<ui> next_KL = KL;
        next_KL.push_back(u);
        vector<ui> new_XR_prime;
        for(auto e : new_XR) if (e != u) new_XR_prime.push_back(e);
        enum_procedure(KR, next_KL, new_CR_prime, new_CL_prime, new_XR_prime, new_XL, ++depth);
        new_XL.push_back(u);
    }
    for(ui i = 0; i < new_CR.size(); i++) {
        ui u = new_CR[i];
        if(CR_flag[i] == 1) { CR_unexpanded_vs.push_back(u); continue; }
        vector<ui> new_CR_prime;
        for(ui j = i + 1; j < new_CR.size(); j++) new_CR_prime.push_back(new_CR[j]); //remove u from new_CR
        new_CR_prime.insert(new_CR_prime.end(), CR_unexpanded_vs.begin(), CR_unexpanded_vs.end());
        vector<ui> new_CL_prime;
        for(auto e : CL_unexpanded_vs) if (e != u) new_CL_prime.push_back(e);
        vector<ui> next_KR = KR;
        next_KR.push_back(u);
        vector<ui> new_XL_prime;
        for(auto e : new_XL) if (e != u) new_XL_prime.push_back(e);
        enum_procedure(next_KR, KL, new_CR_prime, new_CL_prime, new_XR, new_XL_prime, ++depth);
        new_XR.push_back(u);
    }
}

bool edge_prune_inEnum(vector<ui>&KL, vector<ui>&KR, vector<ui>&CL, vector<ui>&CR, vector<pair<ui,ui>>&rmpe, vector<pair<ui,ui>>&rmne)
{
#ifdef _DEBUG_
    cout<<"in edge_prune_inEnum()"<<endl;
    cout<<"KL: "; for(auto &u : KL) cout<<u<<","; cout<<endl;
    cout<<"KR: "; for(auto &u : KR) cout<<u<<","; cout<<endl;
    cout<<"CL: "; for(auto &u : CL) cout<<u<<","; cout<<endl;
    cout<<"CR: "; for(auto &u : CR) cout<<u<<","; cout<<endl;
#endif
#ifdef _CostlyAssert_
    for(auto &u : KL) assert(part[u]==1 || part[u]==2);
    for(auto &u : KR) assert(part[u]==1 || part[u]==2);
    for(auto &u : CL) assert(part[u]==1 || part[u]==2);
    for(auto &u : CR) assert(part[u]==1 || part[u]==2);
#endif
    vector<ui> S;
    S.insert(S.end(), KL.begin(), KL.end());
    S.insert(S.end(), KR.begin(), KR.end());
    S.insert(S.end(), CL.begin(), CL.end());
    S.insert(S.end(), CR.begin(), CR.end());
    for(auto &u : S) {
        Pnei[u].clear();
        Nnei[u].clear();
    }
    for(ui i = 0; i < S.size(); i ++) {
        ui u = S[i];
        for(ui j = i + 1; j < S.size(); j++) {
            ui v = S[j];
            if(part[u]==part[v]) { //u,v same side
                if(M[u][v]==1) { Pnei[u].push_back(v); Pnei[v].push_back(u); }
            }
            else { //u,v diff side
                if(M[u][v]==-1) { Nnei[u].push_back(v); Nnei[v].push_back(u); }
            }
        }
    }
#ifdef _DEBUG_
    cout<<"Pnei and Nnei : "<<endl;
    for(auto &u : S) {
        cout<<"vertex "<<u<<" : "<<endl;
        cout<<"\tPnei : "; for(auto &v : Pnei[u]) cout<<v<<","; cout<<endl;
        cout<<"\tNnei : "; for(auto &v : Nnei[u]) cout<<v<<","; cout<<endl;
    }
#endif
    for(auto &u : S) {
        for(auto &v : Pnei[u]) {
            tri_pp[u][v] = 0;
            tri_nn[u][v] = 0;
        }
        for(auto &v : Nnei[u]) {
            tri_pn[u][v] = 0;
            tri_np[u][v] = 0;
        }
        core[u] = (ui)Pnei[u].size() + (ui)Nnei[u].size(); // for triangle counting
        v_sta[u] = 0; // for triangle counting
    }
    sort(S.begin(), S.end(), mycomp_r);
#ifdef _CostlyAssert_
    for(auto &u : S) {
        for(auto &v : S) {
            if(u == v) continue;
            if(part[u]==part[v] && M[u][v]==1) {
                assert(tri_pp[u][v] == 0 && tri_pp[v][u] == 0);
                assert(tri_nn[u][v] == 0 && tri_nn[v][u] == 0);
            }
            else if (part[u] != part[v] && M[u][v]==-1) {
                assert(tri_pn[u][v] == 0 && tri_pn[v][u] == 0);
                assert(tri_np[u][v] == 0 && tri_np[v][u] == 0);
            }
        }
    }
#endif
#ifdef _DEBUG_
    cout<<"after sort, S :"<<endl; for(auto &u : S) cout<<u<<":"<<core[u]<<", "; cout<<endl;
#endif
//    cout<<"start count triangle"<<endl;
    //count triangle
    for(auto &u : S) {
        for(auto &v : Pnei[u]) {
            if(v_sta[v]) continue;
            for(auto &w : Pnei[v]) {
                if(v_sta[w]) continue;
                if(w==u)continue;
                assert(part[u]==part[v] && part[v]==part[w]);
                if(w>v && M[u][w]==1) {
                    ++ tri_pp[u][v]; ++ tri_pp[v][u];
                    ++ tri_pp[v][w]; ++ tri_pp[w][v];
                    ++ tri_pp[u][w]; ++ tri_pp[w][u];
                }
            }
            for(auto &w : Nnei[v]) {
                if(v_sta[w]) continue;
                assert(w != u);
                assert(part[u]==part[v] && part[v]!=part[w]);
                if(w>v && M[u][w]==-1) {
                    ++ tri_nn[u][v]; ++ tri_nn[v][u];
                    if(part[u]==1) {
                        ++ tri_pn[v][w]; ++ tri_pn[w][v];
                        ++ tri_pn[u][w]; ++ tri_pn[w][u];
                    }
                    else {
                        assert(part[u]==2);
                        ++ tri_np[v][w]; ++ tri_np[w][v];
                        ++ tri_np[u][w]; ++ tri_np[w][u];
                    }
                }
            }
        }
        for(auto &v : Nnei[u]) {
            if(v_sta[v]) continue;
            for(auto &w : Pnei[v]) {
                if(v_sta[w]) continue;
                assert(w != u);
                assert(part[u] != part[w] && part[v] == part[w]);
                if(w>v && M[u][w]==-1) {
                    ++ tri_nn[v][w]; ++ tri_nn[w][v];
                    if(part[u]==1) {
                        ++ tri_np[u][v]; ++ tri_np[v][u];
                        ++ tri_np[u][w]; ++ tri_np[w][u];
                    }
                    else {
                        assert(part[u]==2);
                        ++ tri_pn[u][v]; ++ tri_pn[v][u];
                        ++ tri_pn[u][w]; ++ tri_pn[w][u];
                    }
                }
            }
            for(auto &w : Nnei[v]) {
                if(v_sta[w]) continue;
                if(w==u) continue;
                assert(part[u]==part[w] && part[v] != part[w]);
                if(w>v && M[u][w]==1) {
                    ++ tri_nn[u][w]; ++ tri_nn[w][u];
                    if(part[u]==1) {
                        ++ tri_pn[u][v]; ++ tri_pn[v][u];
                        ++ tri_pn[v][w]; ++ tri_pn[w][v];
                    }
                    else {
                        assert(part[u]==2);
                        ++ tri_np[u][v]; ++ tri_np[v][u];
                        ++ tri_np[v][w]; ++ tri_np[w][v];
                    }
                }
            }
        }
        v_sta[u] = 1;
    }
#ifdef _DEBUG_
    cout<<"tri_pp >= "<<tau-2*k<<endl;
    cout<<"tri_nn >= "<<tau-2*k+2<<endl;
    cout<<"tri_pp + tri_nn >= "<<2*tau-2*k<<endl;
    cout<<"tri_pn >= "<<tau-2*k+1<<endl;
    cout<<"tri_np >= "<<tau-2*k+1<<endl;
    cout<<"tri_pn + tri_np >= "<<2*tau-2*k<<endl;
    cout<<"** edge support information : "<<endl;
    for(auto &u : S) {
        cout<<"for vertex "<<u<<endl;
        for(auto &v : Pnei[u]) {
            cout<<"\ttri_pp["<<u<<"]["<<v<<"] = "<<tri_pp[u][v]<<endl;
            cout<<"\ttri_nn["<<u<<"]["<<v<<"] = "<<tri_nn[u][v]<<endl;
        }
        for(auto &v : Nnei[u]) {
            cout<<"\ttri_pn["<<u<<"]["<<v<<"] = "<<tri_pn[u][v]<<endl;
            cout<<"\ttri_np["<<u<<"]["<<v<<"] = "<<tri_np[u][v]<<endl;
        }
    }
#endif
    //initialize Qe
//    vector<pair<ui, ui>> rmpe, rmne;
    queue<pair<ui, ui>> Qe;
    for(auto &u : S) {
        for(auto &v : Pnei[u]) if(v > u) {
            if(tri_pp[u][v]<tau-2*k || tri_nn[u][v]<tau-2*k+2 || tri_pp[u][v]+tri_nn[u][v]<2*tau-2*k) {
                Qe.push(make_pair(u, v));
                assert(inQe[u][v] == 0);
                assert(inQe[v][u] == 0);
                inQe[u][v] = 1;
                inQe[v][u] = 1;
            }
        }
        for(auto &v : Nnei[u]) if(v > u) {
            if(tri_pn[u][v]<tau-2*k+1 || tri_np[u][v]<tau-2*k+1 || tri_pn[u][v]+tri_np[u][v]<2*tau-2*k){
                Qe.push(make_pair(u, v));
                assert(inQe[u][v] == 0);
                assert(inQe[v][u] == 0);
                inQe[u][v] = 1;
                inQe[v][u] = 1;
            }
        }
    }
#ifndef _VRinEnum_
    cout<<"VR should be before ER!!!"<<endl;
    exit(1);
#endif
    queue<ui> Qv;
#ifdef _CostlyAssert_
    for(auto &u : S) assert(inQv[u]==0);
#endif
    while (!Qe.empty()) {
        pair<ui, ui> edge = Qe.front();
        Qe.pop();
        ui u = edge.first;
        ui v = edge.second;
        assert(inQe[u][v]==1 && inQe[v][u]==1);
        if(M[u][v]==1) { //(u,v): +
            M[u][v] = 0; //remove this edge
            M[v][u] = 0; //remove this edge
            rmpe.push_back(make_pair(u, v));
            assert(p_degree[u]>0 && p_degree[v]>0);
            --p_degree[u];
            --p_degree[v];
            if(inQv[u]==0 && (p_degree[u]+k<tau || p_degree[u]+n_degree[u]+k<2*tau)){
                Qv.push(u); inQv[u] = 1;
            }
            if(inQv[v]==0 && (p_degree[v]+k<tau || p_degree[v]+n_degree[v]+k<2*tau)){
                Qv.push(v); inQv[v] = 1;
            }
            for(auto &w : Pnei[u]) {
                if(M[u][w]==1 && M[v][w]==1) {
                    assert(tri_pp[u][w]>0 && tri_pp[w][u]>0);
                    -- tri_pp[u][w]; -- tri_pp[w][u];
                    if(inQe[u][w]==0 && (tri_pp[u][w]<tau-2*k || tri_pp[u][w]+tri_nn[u][w]<2*tau-2*k) ) {
                        Qe.push(make_pair(u, w));
                        inQe[u][w] = 1;
                        inQe[w][u] = 1;
                    }
                    assert(tri_pp[v][w]>0 && tri_pp[w][v]>0);
                    -- tri_pp[v][w]; -- tri_pp[w][v];
                    if(inQe[v][w]==0 && (tri_pp[v][w]<tau-2*k || tri_pp[v][w]+tri_nn[v][w]<2*tau-2*k) ) {
                        Qe.push(make_pair(v, w));
                        inQe[v][w] = 1;
                        inQe[w][v] = 1;
                    }
                }
            }
            for(auto &w : Nnei[u]) {
                if(M[u][w]==-1 && M[v][w]==-1) {
                    if(part[u]==1) {
                        assert(tri_pn[u][w]>0 && tri_pn[w][u]>0 && tri_pn[v][w]>0 && tri_pn[w][v]>0);
                        --tri_pn[u][w]; --tri_pn[w][u];
                        --tri_pn[v][w]; --tri_pn[w][v];
                    }
                    else {
                        assert(part[u]==2);
                        assert(tri_np[u][w]>0 && tri_np[w][u]>0 && tri_np[v][w]>0 && tri_np[w][v]>0);
                        --tri_np[u][w]; --tri_np[w][u];
                        --tri_np[v][w]; --tri_np[w][v];
                    }
                    if(inQe[u][w]==0 && (tri_pn[u][w]<tau-2*k+1 || tri_np[u][w]<tau-2*k+1 || tri_pn[u][w]+tri_np[u][w]<2*tau-2*k) ) {
                        Qe.push(make_pair(u, w));
                        inQe[u][w] = 1;
                        inQe[w][u] = 1;
                    }
                    if(inQe[v][w]==0 && (tri_pn[v][w]<tau-2*k+1 || tri_np[v][w]<tau-2*k+1 || tri_pn[v][w]+tri_np[v][w]<2*tau-2*k) ) {
                        Qe.push(make_pair(v, w));
                        inQe[v][w] = 1;
                        inQe[w][v] = 1;
                    }
                }
            }
        }
        else { //(u,v): -
            assert(M[u][v]==-1);
            M[u][v] = 0; //remove this edge
            M[v][u] = 0; //remove this edge
            rmne.push_back(make_pair(u, v));
            assert(n_degree[u]>0 && n_degree[v]>0);
            --n_degree[u];
            --n_degree[v];
            if(inQv[u]==0 && (n_degree[u]+k-1<tau || p_degree[u]+n_degree[u]+k<2*tau)){
                Qv.push(u); inQv[u] = 1;
            }
            if(inQv[v]==0 && (n_degree[v]+k-1<tau || p_degree[v]+n_degree[v]+k<2*tau)){
                Qv.push(v); inQv[v] = 1;
            }
            for(auto &w : Pnei[u]) {
                if(M[u][w]==1 && M[v][w]==-1) {
                    assert(tri_nn[u][w]>0 && tri_nn[w][u]>0);
                    -- tri_nn[u][w]; -- tri_nn[w][u];
                    if(inQe[u][w]==0 && (tri_nn[u][w]<tau-2*k+2 || tri_pp[u][w]+tri_nn[u][w]<2*tau-2*k) ) {
                        Qe.push(make_pair(u, w));
                        inQe[u][w] = 1;
                        inQe[w][u] = 1;
                    }
                    if(part[u]==1) {
                        assert(tri_pn[v][w]>0 && tri_pn[w][v]>0);
                        -- tri_pn[v][w]; -- tri_pn[w][v];
                    }
                    else {
                        assert(part[u]==2);
                        assert(tri_np[v][w]>0 && tri_np[w][v]>0);
                        -- tri_np[v][w]; -- tri_np[w][v];
                    }
                    if(inQe[v][w]==0 && (tri_pn[v][w]<tau-2*k+1 || tri_np[v][w]<tau-2*k+1 || tri_pn[v][w]+tri_np[v][w]<2*tau-2*k) ) {
                        Qe.push(make_pair(v, w));
                        inQe[v][w] = 1;
                        inQe[w][v] = 1;
                    }
                }
            }
            for(auto &w : Nnei[u]) {
                if(M[u][w]==-1 && M[v][w]==1) {
                    assert(tri_nn[v][w]>0 && tri_nn[w][v]>0);
                    -- tri_nn[v][w]; -- tri_nn[w][v];
                    if(inQe[v][w]==0 && (tri_nn[v][w]<tau-2*k+2 || tri_pp[v][w]+tri_nn[v][w]<2*tau-2*k) ) {
                        Qe.push(make_pair(v, w));
                        inQe[v][w] = 1;
                        inQe[w][v] = 1;
                    }
                    if(part[u]==1) {
                        assert(tri_np[u][w]>0 && tri_np[w][u]>0);
                        -- tri_np[u][w]; -- tri_np[w][u];
                    }
                    else {
                        assert(part[u]==2);
                        assert(tri_pn[u][w]>0 && tri_pn[w][u]>0);
                        -- tri_pn[u][w]; -- tri_pn[w][u];
                    }
                    if(inQe[u][w]==0 && (tri_pn[u][w]<tau-2*k+1 || tri_np[u][w]<tau-2*k+1 || tri_pn[u][w]+tri_np[u][w]<2*tau-2*k) ) {
                        Qe.push(make_pair(u, w));
                        inQe[u][w] = 1;
                        inQe[w][u] = 1;
                    }
                }
            }
        }
    }
    for(auto &e : rmpe) {
        ui u = e.first;
        ui v = e.second;
        assert(inQe[u][v]==1 && inQe[v][u]==1);
        inQe[u][v] = 0;
        inQe[v][u] = 0;
    }
    for(auto &e : rmne) {
        ui u = e.first;
        ui v = e.second;
        assert(inQe[u][v]==1 && inQe[v][u]==1);
        inQe[u][v] = 0;
        inQe[v][u] = 0;
    }
    er_round++;
    er_prune_cnt += (rmpe.size() + rmne.size());
    while (!Qv.empty()) {
        ui u = Qv.front();
#ifdef _DEBUG_
        cout<<"delete "<<u<<" in ER"<<endl;
#endif
        Qv.pop();
        ++ vr_prune_II;
        for(auto &v : Pnei[u]) {
            if(inQv[v]==0 && M[u][v]==1) {
                assert(p_degree[v]>0);
                --p_degree[v];
                if(p_degree[v]+k<tau || p_degree[v]+n_degree[v]+k<2*tau) {Qv.push(v); inQv[v] = 1;}
            }
        }
        for(auto &v : Nnei[u]) {
            if(inQv[v]==0 && M[u][v]==-1) {
                assert(n_degree[v]>0);
                --n_degree[v];
                if(n_degree[v]+k-1<tau || p_degree[v]+n_degree[v]+k<2*tau) {Qv.push(v); inQv[v] = 1;}
            }
        }
    }
    for(auto &e : KL) if(inQv[e]) return true;
    for(auto &e : KR) if(inQv[e]) return true;
    vector<ui> tmpCL, tmpCR;
    for(auto &e : CL) if(inQv[e]==0) tmpCL.push_back(e);
    for(auto &e : CR) if(inQv[e]==0) tmpCR.push_back(e);
    CL=tmpCL;
    CR=tmpCR;
    return false;
}

bool is_maximality_basic(vector<ui> &KL, vector<ui> &KR, vector<ui> &XL, vector<ui> &XR)
{
    for(auto &e : XL) {
        bool flag = true;
        short tmp_capacity = k - 1;
        for(auto &x : KL) {
            if(M[e][x] != 1) {
                -- tmp_capacity;
                if(tmp_capacity < 0 || partial_kplex_cap[x] <= 0) {
                    flag = false;
                    break;
                }
            }
        }
        if(flag == false) continue;
        for(auto &x : KR) {
            if(M[e][x] != (-1)) {
                -- tmp_capacity;
                if(tmp_capacity < 0 || partial_kplex_cap[x] <= 0) {
                    flag = false;
                    break;
                }
            }
        }
        if(flag == true) return false;
    }
    for(auto &e : XR) {
        bool flag = true;
        short tmp_capacity = k - 1;
        for(auto &x : KL) {
            if(M[e][x] != -1) {
                -- tmp_capacity;
                if(tmp_capacity < 0 || partial_kplex_cap[x] <= 0) {
                    flag = false;
                    break;
                }
            }
        }
        if(flag == false) continue;
        for(auto &x : KR) {
            if(M[e][x] != (1)) {
                -- tmp_capacity;
                if(tmp_capacity < 0 || partial_kplex_cap[x] <= 0) {
                    flag = false;
                    break;
                }
            }
        }
        if(flag == true) return false;
    }
    return true;
}

bool check_maximality(vector<ui> &KL, vector<ui> &KR, vector<ui> &CL, vector<ui> &CR, vector<ui> &XL, vector<ui> &XR, bool trueL)
{
    int Ssize = (int)KL.size() + (int)KR.size() + (int)CL.size() + (int)CR.size();
    for(auto &u : KL) {
        assert(Ssize > p_degree[u] + n_degree[u]);
        int nneinum = Ssize - (p_degree[u] + n_degree[u]);
        partial_kplex_cap[u] = k - nneinum;
        assert(partial_kplex_cap[u] >= 0);
    }
    for(auto &u : KR) {
        assert(Ssize > p_degree[u] + n_degree[u]);
        int nneinum = Ssize - (p_degree[u] + n_degree[u]);
        partial_kplex_cap[u] = k - nneinum;
        assert(partial_kplex_cap[u] >= 0);
    }
    for(auto &u : CL) {
        assert(Ssize > p_degree[u] + n_degree[u]);
        int nneinum = Ssize - (p_degree[u] + n_degree[u]);
        partial_kplex_cap[u] = k - nneinum;
        assert(partial_kplex_cap[u] >= 0);
    }
    for(auto &u : CR) {
        assert(Ssize > p_degree[u] + n_degree[u]);
        int nneinum = Ssize - (p_degree[u] + n_degree[u]);
        partial_kplex_cap[u] = k - nneinum;
        assert(partial_kplex_cap[u] >= 0);
    }
    
    /*/*/
    vector<ui> tKL, tKR;
    tKL.insert(tKL.begin(), KL.begin(), KL.end());
    tKL.insert(tKL.begin(), CL.begin(), CL.end());
    tKR.insert(tKR.begin(), KR.begin(), KR.end());
    tKR.insert(tKR.begin(), CR.begin(), CR.end());
    bool exist_full_v = false;
    for(auto &u: tKL) if(partial_kplex_cap[u]==0) {
        exist_full_v = true;
        for(ui i = p_pstart[u]; i < p_pend[u]; i ++) if(inXL[p_edges[i]]==1) {
            ui v = p_edges[i]; //test whether v can be added into K
            bool vcan = true;
            int vmissednum = 0;
            for(auto &w : tKL) if(M[v][w]!=1) {
                ++vmissednum;
                if(vmissednum>=k || partial_kplex_cap[w]==0) {
                    vcan = false; break;
                }
            }
            if(vcan==false) continue;
            for(auto &w : tKR) if(M[v][w]!=-1) {
                ++vmissednum;
                if(vmissednum>=k || partial_kplex_cap[w]==0) {
                    vcan = false; break;
                }
            }
            if(vcan) return false;
        }
        for(ui i = n_pstart[u]; i < n_pend[u]; i ++) if(inXR[n_edges[i]]==1) {
            ui v = n_edges[i]; //test whether v can be added into K
            bool vcan = true;
            int vmissednum = 0;
            for(auto &w : tKL) if(M[v][w]!=-1) {
                ++vmissednum;
                if(vmissednum>=k || partial_kplex_cap[w]==0) {
                    vcan = false; break;
                }
            }
            if(vcan==false) continue;
            for(auto &w : tKR) if(M[v][w]!=1) {
                ++vmissednum;
                if(vmissednum>=k || partial_kplex_cap[w]==0) {
                    vcan = false; break;
                }
            }
            if(vcan) return false;
        }
        return true;
    }
    for(auto &u: tKR) if(partial_kplex_cap[u]==0) {
        exist_full_v = true;
        for(ui i = p_pstart[u]; i < p_pend[u]; i ++) if(inXR[p_edges[i]]==1) {
            ui v = p_edges[i]; //test whether v can be added into K
            bool vcan = true;
            int vmissednum = 0;
            for(auto &w : tKL) if(M[v][w]!=-1) {
                ++vmissednum;
                if(vmissednum>=k || partial_kplex_cap[w]==0) {
                    vcan = false; break;
                }
            }
            if(vcan==false) continue;
            for(auto &w : tKR) if(M[v][w]!=1) {
                ++vmissednum;
                if(vmissednum>=k || partial_kplex_cap[w]==0) {
                    vcan = false; break;
                }
            }
            if(vcan) return false;
        }
        for(ui i = n_pstart[u]; i < n_pend[u]; i ++) if(inXL[n_edges[i]]==1) {
            ui v = n_edges[i]; //test whether v can be added into K
            bool vcan = true;
            int vmissednum = 0;
            for(auto &w : tKL) if(M[v][w]!=1) {
                ++vmissednum;
                if(vmissednum>=k || partial_kplex_cap[w]==0) {
                    vcan = false; break;
                }
            }
            if(vcan==false) continue;
            for(auto &w : tKR) if(M[v][w]!=-1) {
                ++vmissednum;
                if(vmissednum>=k || partial_kplex_cap[w]==0) {
                    vcan = false; break;
                }
            }
            if(vcan) return false;
        }
        return true;
    }
//    if(exist_full_v) return true;
    /*/*/
    for(auto &u : XL) {
        //test whether u can be added into k-plex
        int umissednum = 0;
        bool ucan = true;
        for(auto &v : KL) if(M[u][v]!=1) {
            umissednum++;
            if(umissednum >=k || partial_kplex_cap[v]==0) {
                ucan = false; break;
            }
        }
        if(ucan==false) continue;
        for(auto &v : KR) if(M[u][v]!=-1) {
            umissednum++;
            if(umissednum >=k || partial_kplex_cap[v]==0) {
                ucan = false; break;
            }
        }
        if(ucan==false) continue;
        for(auto &v : CL) if(M[u][v]!=1) {
            umissednum++;
            if(umissednum >=k || partial_kplex_cap[v]==0) {
                ucan = false; break;
            }
        }
        if(ucan==false) continue;
        for(auto &v : CR) if(M[u][v]!=-1) {
            umissednum++;
            if(umissednum >=k || partial_kplex_cap[v]==0) {
                ucan = false; break;
            }
        }
        if(ucan) return false;
    }
    for(auto &u : XR) {
        //test whether u can be added into k-plex
        int umissednum = 0;
        bool ucan = true;
        for(auto &v : KL) if(M[u][v]!=-1) {
            umissednum++;
            if(umissednum >=k || partial_kplex_cap[v]==0) {
                ucan = false; break;
            }
        }
        if(ucan==false) continue;
        for(auto &v : KR) if(M[u][v]!=1) {
            umissednum++;
            if(umissednum >=k || partial_kplex_cap[v]==0) {
                ucan = false; break;
            }
        }
        if(ucan==false) continue;
        for(auto &v : CL) if(M[u][v]!=-1) {
            umissednum++;
            if(umissednum >=k || partial_kplex_cap[v]==0) {
                ucan = false; break;
            }
        }
        if(ucan==false) continue;
        for(auto &v : CR) if(M[u][v]!=1) {
            umissednum++;
            if(umissednum >=k || partial_kplex_cap[v]==0) {
                ucan = false; break;
            }
        }
        if(ucan) return false;
    }
    return true;
}

bool check_maximality(vector<ui> &KL, vector<ui> &KR, vector<ui> &XL, vector<ui> &XR, bool trueL)
{
    vector<ui> tKL, tKR;
    if(trueL) {
        tKL = KL; tKR = KR;
    }
    else {
        tKL = KR; tKR = KL;
    }
    bool exist_full_v = false;
    for(auto &u: tKL) if(partial_kplex_cap[u]==0) {
        exist_full_v = true;
        for(ui i = p_pstart[u]; i < p_pend[u]; i ++) if(inXL[p_edges[i]]==1) {
            ui v = p_edges[i]; //test whether v can be added into K
            bool vcan = true;
            int vmissednum = 0;
            for(auto &w : tKL) if(M[v][w]!=1) {
                ++vmissednum;
                if(vmissednum>=k || partial_kplex_cap[w]==0) {
                    vcan = false; break;
                }
            }
            if(vcan==false) continue;
            for(auto &w : tKR) if(M[v][w]!=-1) {
                ++vmissednum;
                if(vmissednum>=k || partial_kplex_cap[w]==0) {
                    vcan = false; break;
                }
            }
            if(vcan) return false;
        }
        for(ui i = n_pstart[u]; i < n_pend[u]; i ++) if(inXR[n_edges[i]]==1) {
            ui v = n_edges[i]; //test whether v can be added into K
            bool vcan = true;
            int vmissednum = 0;
            for(auto &w : tKL) if(M[v][w]!=-1) {
                ++vmissednum;
                if(vmissednum>=k || partial_kplex_cap[w]==0) {
                    vcan = false; break;
                }
            }
            if(vcan==false) continue;
            for(auto &w : tKR) if(M[v][w]!=1) {
                ++vmissednum;
                if(vmissednum>=k || partial_kplex_cap[w]==0) {
                    vcan = false; break;
                }
            }
            if(vcan) return false;
        }
        return true;
    }
    for(auto &u: tKR) if(partial_kplex_cap[u]==0) {
        exist_full_v = true;
        for(ui i = p_pstart[u]; i < p_pend[u]; i ++) if(inXR[p_edges[i]]==1) {
            ui v = p_edges[i]; //test whether v can be added into K
            bool vcan = true;
            int vmissednum = 0;
            for(auto &w : tKL) if(M[v][w]!=-1) {
                ++vmissednum;
                if(vmissednum>=k || partial_kplex_cap[w]==0) {
                    vcan = false; break;
                }
            }
            if(vcan==false) continue;
            for(auto &w : tKR) if(M[v][w]!=1) {
                ++vmissednum;
                if(vmissednum>=k || partial_kplex_cap[w]==0) {
                    vcan = false; break;
                }
            }
            if(vcan) return false;
        }
        for(ui i = n_pstart[u]; i < n_pend[u]; i ++) if(inXL[n_edges[i]]==1) {
            ui v = n_edges[i]; //test whether v can be added into K
            bool vcan = true;
            int vmissednum = 0;
            for(auto &w : tKL) if(M[v][w]!=1) {
                ++vmissednum;
                if(vmissednum>=k || partial_kplex_cap[w]==0) {
                    vcan = false; break;
                }
            }
            if(vcan==false) continue;
            for(auto &w : tKR) if(M[v][w]!=-1) {
                ++vmissednum;
                if(vmissednum>=k || partial_kplex_cap[w]==0) {
                    vcan = false; break;
                }
            }
            if(vcan) return false;
        }
        return true;
    }
//    if(exist_full_v) return true;
    
    for(auto &u : XL) {
        //test whether u can be added into K
        int umissednum = 0;
        bool ucan = true;
        for(auto &v : KL) if(M[u][v]!=1) {
            umissednum++;
            assert(partial_kplex_cap[v]>0);
            if(umissednum >=k) {
                ucan = false; break;
            }
        }
        if(ucan==false) continue;
        for(auto &v : KR) if(M[u][v]!=-1) {
            umissednum++;
            assert(partial_kplex_cap[v]>0);
            if(umissednum >=k) {
                ucan = false; break;
            }
        }
        if(ucan) return false;
    }
    for(auto &u : XR) {
        //test whether u can be added into K
        int umissednum = 0;
        bool ucan = true;
        for(auto &v : KL) if(M[u][v]!=-1) {
            umissednum++;
            assert(partial_kplex_cap[v]>0);
            if(umissednum >=k) {
                ucan = false; break;
            }
        }
        if(ucan==false) continue;
        for(auto &v : KR) if(M[u][v]!=1) {
            umissednum++;
            assert(partial_kplex_cap[v]>0);
            if(umissednum >=k) {
                ucan = false; break;
            }
        }
        if(ucan) return false;
    }
    return true;
}

void enum_procedure_mindeg_pivot(vector<ui> &myKL, vector<ui> &myKR, vector<ui> &CL, vector<ui> &CR, vector<ui> &XL, vector<ui> &XR, int depth)
{
    ++ branch_num;
    vector<ui> KL=myKL, KR=myKR;
    bool tmpflag = obtain_partial_kplex_capacity(KL, KR);
    if(tmpflag) return; //compute vertex degree in partial kplex <KL,KR>
    vector<ui> new_CL; refine_S_by_capacity_mdpivot(new_CL, CL, KL, KR, true); //refine CL
    vector<ui> new_CR; refine_S_by_capacity_mdpivot(new_CR, CR, KL, KR, false); //refine CR
#ifdef _maintainX_
    vector<ui> new_XL; refine_S_by_capacity(new_XL, XL, KL, KR, true); //refine XL
    vector<ui> new_XR; refine_S_by_capacity(new_XR, XR, KL, KR, false); //refine XR
#endif
    if(KL.size() + new_CL.size() < tau || KR.size() + new_CR.size() < tau) return;
    if(new_CL.empty() && new_CR.empty() ) {
#ifdef _maintainX_
        if(new_XL.empty() && new_XR.empty()) {
            #ifdef _StoreRes_
            vector<ui> k1, k2;
            for(auto e : KL) k1.push_back(original_id[e]); for(auto e : KR) k2.push_back(original_id[e]);
            Res.push_back(make_pair(k1, k2));
            #endif
            int kplesize = (int)KL.size() + (int)KR.size();
            if( kplesize > max_kplex) max_kplex = kplesize ;
            if( kplesize < min_kplex) min_kplex = kplesize;
            ++ Res_num;
        }
#else
        bool ismax = check_maximality(KL, KR, XL, XR, true);
        if(ismax) {
            #ifdef _StoreRes_
            vector<ui> k1, k2;
            for(auto e : KL) k1.push_back(original_id[e]); for(auto e : KR) k2.push_back(original_id[e]);
            Res.push_back(make_pair(k1, k2));
            #endif
            int kplesize = (int)KL.size() + (int)KR.size();
            if( kplesize > max_kplex) max_kplex = kplesize ;
            if( kplesize < min_kplex) min_kplex = kplesize;
            ++ Res_num;
        }
#endif
        return;
    }
#ifdef _VRinEnum_
    if(VR_in_Enum(KL,KR,new_CL,new_CR) || (new_CL.empty()&&new_CR.empty())) return;
#endif
#ifdef _ERinEnum_
    vector<pair<ui, ui>> rmpe, rmne;
    if(depth < er_threshold){
        bool tmpbool = edge_prune_inEnum(KL, KR, new_CL, new_CR, rmpe, rmne);
        if(tmpbool || (new_CL.empty() && new_CR.empty()) ) {
            for(auto &e : rmpe) {
                ui u = e.first, v = e.second;
                assert(M[u][v]==0 && M[v][u]==0);
                M[u][v]=1; M[v][u]=1;
            }
            for(auto &e : rmne) {
                ui u = e.first, v = e.second;
                assert(M[u][v]==0 && M[v][u]==0);
                M[u][v]=-1; M[v][u]=-1;
            }
            return;
        }
    }
#endif
    int ustar_side = -1;
    ui ustar_pos;
    ui ustar = get_mindeg_vertex(KL, KR, new_CL, new_CR, ustar_side, ustar_pos);
    assert(ustar_side != -1);
    assert(KL.size()+KR.size()+new_CL.size()+new_CR.size() > p_degree[ustar]+n_degree[ustar]);
    if((KL.size()+KR.size()+new_CL.size()+new_CR.size()) - (p_degree[ustar]+n_degree[ustar]) <= k) {
#ifdef _maintainX_
        timer_in_enum.restart();
        bool ismaximal = check_maximality(KL, KR, new_CL, new_CR, new_XL, new_XR, true);
        T_check_max += timer_in_enum.elapsed();
#else
        bool ismaximal = check_maximality(KL, KR, new_CL, new_CR, XL, XR, true);
#endif
        if(ismaximal) {
            if(KL.size()+new_CL.size()>=tau && KR.size()+new_CR.size()>=tau) {
                #ifdef _StoreRes_
                vector<ui> k1, k2;
                for(auto e : KL) k1.push_back(original_id[e]); for(auto e : new_CL) k1.push_back(original_id[e]);
                for(auto e : KR) k2.push_back(original_id[e]); for(auto e : new_CR) k2.push_back(original_id[e]);
                Res.push_back(make_pair(k1, k2));
                #endif
                int kplesize = (int)KL.size() + (int)KR.size() + (int)new_CL.size() + (int)new_CR.size();
                if( kplesize > max_kplex) max_kplex = kplesize ;
                if( kplesize < min_kplex) min_kplex = kplesize;
                ++ Res_num;
            }
        }
    }
    else if (ustar_side == 1 || ustar_side == 2) {
        int kprime = partial_kplex_cap[ustar];
        assert(kprime >= 1);
        int tmpsign = 1;
        if(ustar_side == 2) tmpsign = -1;
        vector<ui> nonnei_inC;
        for(ui i = 0; i< new_CL.size(); i++) {
            ui v = new_CL[i];
            if(M[ustar][v] != tmpsign) {
                nonnei_inC.push_back(v);
                assert(part[v]==1);
            }
        }
        for(ui i = 0; i< new_CR.size(); i++) {
            ui v = new_CR[i];
            if(M[ustar][v] != (-tmpsign)) {
                nonnei_inC.push_back(v);
                assert(part[v]==2);
            }
        }
        assert(nonnei_inC.size() > kprime);
        for(int i = 0; i < kprime; i++) {
            ui v = nonnei_inC[i];
            if(part[v] == 1) {
                auto it = find(new_CL.begin(), new_CL.end(), v);
                assert(it != new_CL.end());
                new_CL.erase(it);
#ifdef _maintainX_
                new_XL.push_back(v);
                enum_procedure_mindeg_pivot(KL, KR, new_CL, new_CR, new_XL, new_XR, depth++);
                assert(new_XL.back()==v);
                new_XL.pop_back();
#else
                assert(inXL[v]==0);
                XL.push_back(v);
                inXL[v]=1;
                enum_procedure_mindeg_pivot(KL, KR, new_CL, new_CR, XL, XR, depth++);
                assert(XL.back()==v && inXL[v]==1);
                XL.pop_back();
                inXL[v]=0;
#endif
                KL.push_back(v);
            }
            else {
                auto it = find(new_CR.begin(), new_CR.end(), v);
                assert(it != new_CR.end());
                new_CR.erase(it);
#ifdef _maintainX_
                new_XR.push_back(v);
                enum_procedure_mindeg_pivot(KL, KR, new_CL, new_CR, new_XL, new_XR, depth++);
                assert(new_XR.back()==v);
                new_XR.pop_back();
#else
                assert(inXR[v]==0);
                XR.push_back(v);
                inXR[v]=1;
                enum_procedure_mindeg_pivot(KL, KR, new_CL, new_CR, XL, XR, depth++);
                assert(XR.back()==v && inXR[v]==1);
                XR.pop_back();
                inXR[v]=0;
#endif
                KR.push_back(v);
            }
        }
        for(ui i = kprime; i < nonnei_inC.size(); i++) {
            ui v = nonnei_inC[i];
            assert(part[v]==1 || part[v]==2);
            if(part[v]==1) {
                auto it = find(new_CL.begin(), new_CL.end(), v);
                assert(it != new_CL.end());
                new_CL.erase(it);
            }
            else {
                auto it = find(new_CR.begin(), new_CR.end(), v);
                assert(it != new_CR.end());
                new_CR.erase(it);
            }
        }
#ifdef _maintainX_
        enum_procedure_mindeg_pivot(KL, KR, new_CL, new_CR, new_XL, new_XR, depth++);
#else
        enum_procedure_mindeg_pivot(KL, KR, new_CL, new_CR, XL, XR, depth++);
#endif
    }
    else if (ustar_side == 3) {
        assert(new_CL[ustar_pos]==ustar && !new_CL.empty());
        if(ustar_pos != new_CL.size()-1) new_CL[ustar_pos] = new_CL[new_CL.size()-1];
        new_CL.pop_back();
#ifdef _maintainX_
        new_XL.push_back(ustar);
        enum_procedure_mindeg_pivot(KL, KR, new_CL, new_CR, new_XL, new_XR, depth++);
        assert(new_XL.back()==ustar);
        new_XL.pop_back();
#else
        assert(inXL[ustar]==0);
        XL.push_back(ustar);
        inXL[ustar] = 1;
        enum_procedure_mindeg_pivot(KL, KR, new_CL, new_CR, XL, XR, depth++);
        assert(XL.back()==ustar && inXL[ustar]==1);
        XL.pop_back();
        inXL[ustar] = 0;
#endif
        KL.push_back(ustar);
#ifdef _maintainX_
        enum_procedure_mindeg_pivot(KL, KR, new_CL, new_CR, new_XL, new_XR, depth++);
#else
        enum_procedure_mindeg_pivot(KL, KR, new_CL, new_CR, XL, XR, depth++);
#endif
    }
    else {
        assert(ustar_side == 4);
        assert(new_CR[ustar_pos]==ustar && !new_CR.empty());
        if(ustar_pos != new_CR.size()-1) new_CR[ustar_pos] = new_CR[new_CR.size()-1];
        new_CR.pop_back();
#ifdef _maintainX_
        new_XR.push_back(ustar);
        enum_procedure_mindeg_pivot(KL, KR, new_CL, new_CR, new_XL, new_XR, depth++);
        assert(new_XR.back()==ustar);
        new_XR.pop_back();
#else
        assert(inXR[ustar]==0);
        XR.push_back(ustar);
        inXR[ustar]=1;
        enum_procedure_mindeg_pivot(KL, KR, new_CL, new_CR, XL, XR, depth++);
        assert(XR.back()==ustar && inXR[ustar] ==1);
        XR.pop_back();
        inXR[ustar] = 0;
#endif
        KR.push_back(ustar);
#ifdef _maintainX_
        enum_procedure_mindeg_pivot(KL, KR, new_CL, new_CR, new_XL, new_XR, depth++);
#else
        enum_procedure_mindeg_pivot(KL, KR, new_CL, new_CR, XL, XR, depth++);
#endif
    }
#ifdef _ERinEnum_
    for(auto &e : rmpe) { ui u = e.first, v = e.second; assert(M[u][v]==0 && M[v][u]==0); M[u][v]=1; M[v][u]=1; }
    for(auto &e : rmne) { ui u = e.first, v = e.second; assert(M[u][v]==0 && M[v][u]==0); M[u][v]=-1; M[v][u]=-1; }
#endif
}

//for this enum procedure, CL and CR are not overlapped
void enum_procedure_adv_opt(vector<ui> &KL, vector<ui> &KR, vector<ui> &CL, vector<ui> &CR, vector<ui> &XL, vector<ui> &XR, int depth, bool trueL)
{
    ++ branch_num;
    bool tmpflag = obtain_partial_kplex_capacity(KL, KR);
    if(tmpflag) return; //compute vertex degree in partial kplex <KL,KR>
    vector<ui> new_CL; refine_S_by_capacity(new_CL, CL, KL, KR, true, 1); //refine CL
    vector<ui> new_CR; refine_S_by_capacity(new_CR, CR, KL, KR, false, 1); //refine CR
    
#ifdef _maintainX_
    vector<ui> new_XL; refine_S_by_capacity(new_XL, XL, KL, KR, true); //refine XL
    vector<ui> new_XR; refine_S_by_capacity(new_XR, XR, KL, KR, false); //refine XR
#endif
    if(KL.size() + new_CL.size() < tau || KR.size() + new_CR.size() < tau) return;
    if(new_CL.empty() && new_CR.empty() ) {
#ifdef _maintainX_
        if(new_XL.empty() && new_XR.empty()) {
            #ifdef _StoreRes_
            vector<ui> k1, k2;
            for(auto e : KL) k1.push_back(original_id[e]); for(auto e : KR) k2.push_back(original_id[e]);
            Res.push_back(make_pair(k1, k2));
            #endif
            int kplesize = (int)KL.size() + (int)KR.size();
            if( kplesize > max_kplex) max_kplex = kplesize ;
            if( kplesize < min_kplex) min_kplex = kplesize;
            ++ Res_num;
        }
#else
        //check maximality
        bool ismax = check_maximality(KL, KR, XL, XR, trueL);
        if(ismax) {
            #ifdef _StoreRes_
            vector<ui> k1, k2;
            for(auto e : KL) k1.push_back(original_id[e]); for(auto e : KR) k2.push_back(original_id[e]);
            Res.push_back(make_pair(k1, k2));
            #endif
            int kplesize = (int)KL.size() + (int)KR.size();
            if( kplesize > max_kplex) max_kplex = kplesize ;
            if( kplesize < min_kplex) min_kplex = kplesize;
            ++ Res_num;
        }
#endif
        return;
    }
#ifdef _VRinEnum_
    if(VR_in_Enum(KL,KR,new_CL,new_CR) || (new_CL.empty()&&new_CR.empty())) return;
#endif
    
#ifdef _ERinEnum_
    vector<pair<ui, ui>> rmpe, rmne;
    if(depth < er_threshold){
        bool tmpbool = edge_prune_inEnum(KL, KR, new_CL, new_CR, rmpe, rmne);
        if(tmpbool || (new_CL.empty() && new_CR.empty()) ) {
            for(auto &e : rmpe) {
                ui u = e.first, v = e.second;
                assert(M[u][v]==0 && M[v][u]==0);
                M[u][v]=1; M[v][u]=1;
            }
            for(auto &e : rmne) {
                ui u = e.first, v = e.second;
                assert(M[u][v]==0 && M[v][u]==0);
                M[u][v]=-1; M[v][u]=-1;
            }
            return;
        }
    }
#endif

#ifdef _maintainX_
#ifdef _ET_
    bool canET = false;
    timer_in_enum.restart();
    canET = early_termination_adv(KL,KR,new_CL,new_CR,new_XL,new_XR);
    T_ET += timer_in_enum.elapsed();
    if(canET) {
        for(auto &e : rmpe) {
            ui u = e.first, v = e.second;
            assert(M[u][v]==0 && M[v][u]==0);
            M[u][v]=1; M[v][u]=1;
        }
        for(auto &e : rmne) {
            ui u = e.first, v = e.second;
            assert(M[u][v]==0 && M[v][u]==0);
            M[u][v]=-1; M[v][u]=-1;
        }
        return;
    }
#endif
#endif
    
    for(auto &u : new_CL) core[u] = p_degree[u] + n_degree[u];
    for(auto &u : new_CR) core[u] = p_degree[u] + n_degree[u];
    sort(new_CL.begin(), new_CL.end(), mycomp);
    sort(new_CR.begin(), new_CR.end(), mycomp);
    vector<ui> CL_unexpanded_vs, CR_unexpanded_vs;
    vector<bool> CL_flag(new_CL.size(),0), CR_flag(new_CR.size(),0);
#ifdef _PIVOT_
    pair<vector<ui>, vector<ui>> skip_vs;
    get_skipvlist_by_pivoting(KL,KR,new_CL,new_CR,skip_vs);
    for(auto e : skip_vs.first) CL_flag[e] = 1;
    for(auto e : skip_vs.second) CR_flag[e] = 1;
#endif
    
#ifndef _maintainX_
    int mvtoXL=0, mvtoXR=0;
#endif
    
    for(ui i = 0; i < new_CL.size(); i++) {
        ui u = new_CL[i];
        if(CL_flag[i] == 1) { CL_unexpanded_vs.push_back(u); continue; }
        vector<ui> new_CL_prime;
        for(ui j = i + 1; j < new_CL.size(); j++) new_CL_prime.push_back(new_CL[j]);
        new_CL_prime.insert(new_CL_prime.end(), CL_unexpanded_vs.begin(), CL_unexpanded_vs.end());
        vector<ui> next_KL = KL;
        next_KL.push_back(u);
        
#ifdef _maintainX_
        enum_procedure_adv_opt(KR, next_KL, new_CR, new_CL_prime, new_XR, new_XL, ++depth, !trueL);
#else
    #ifdef _UB_
        if(UB(next_KL,KR,new_CL_prime,new_CR,u,true) >= 2*tau)
            enum_procedure_adv_opt(KR, next_KL, new_CR, new_CL_prime, XR, XL, ++depth, !trueL);
    #else
        enum_procedure_adv_opt(KR, next_KL, new_CR, new_CL_prime, XR, XL, ++depth, !trueL);
    #endif
            
#endif

#ifdef _maintainX_
        new_XL.push_back(u);
#else
        mvtoXL++;
        if(trueL) { assert(inXL[u]==0); inXL[u]=1; }
        else { assert(inXR[u]==0); inXR[u]=1; }
        XL.push_back(u);
#endif
    }
    
    for(ui i = 0; i < new_CR.size(); i++) {
        ui u = new_CR[i];
        if(CR_flag[i] == 1) { CR_unexpanded_vs.push_back(u); continue; }
        vector<ui> new_CR_prime;
        for(ui j = i + 1; j < new_CR.size(); j++) new_CR_prime.push_back(new_CR[j]); //remove u from new_CR
        new_CR_prime.insert(new_CR_prime.end(), CR_unexpanded_vs.begin(), CR_unexpanded_vs.end());
        vector<ui> next_KR = KR;
        next_KR.push_back(u);
                
#ifdef _maintainX_
        enum_procedure_adv_opt(next_KR, KL, new_CR_prime, CL_unexpanded_vs, new_XR, new_XL, ++depth, !trueL);
#else
    #ifdef _UB_
        if(UB(KL,next_KR,CL_unexpanded_vs,new_CR_prime,u,false) >= 2*tau)
            enum_procedure_adv_opt(next_KR, KL, new_CR_prime, CL_unexpanded_vs, XR, XL, ++depth, !trueL);
    #else
        enum_procedure_adv_opt(next_KR, KL, new_CR_prime, CL_unexpanded_vs, XR, XL, ++depth, !trueL);
    #endif
#endif

#ifdef _maintainX_
        new_XR.push_back(u);
#else
        mvtoXR++;
        if(trueL) { assert(inXR[u]==0); inXR[u]=1; }
        else { assert(inXL[u]==0); inXL[u]=1; }
        XR.push_back(u);
#endif
    }
    
#ifndef _maintainX_
    while (mvtoXL>0) {
        if(trueL) { assert(inXL[XL.back()]==1); inXL[XL.back()]=0; }
        else { assert(inXR[XL.back()]==1); inXR[XL.back()]=0; }
        XL.pop_back();
        --mvtoXL;
    }
    while (mvtoXR>0) {
        if(trueL) { assert(inXR[XR.back()]==1); inXR[XR.back()]=0; }
        else { assert(inXL[XR.back()]==1); inXL[XR.back()]=0; }
        XR.pop_back();
        --mvtoXR;
    }
#endif
    
#ifdef _ERinEnum_
    for(auto &e : rmpe) { ui u = e.first, v = e.second; assert(M[u][v]==0 && M[v][u]==0); M[u][v]=1; M[v][u]=1; }
    for(auto &e : rmne) { ui u = e.first, v = e.second; assert(M[u][v]==0 && M[v][u]==0); M[u][v]=-1; M[v][u]=-1; }
#endif
}

void shrink_CX(ui u, ui sidx, vector<ui>&S, vector<ui>&new_S, bool SisL, vector<ui>&unn_KL, vector<ui>&unn_KR, ui*LdK, ui*RdK, ui ks)
{
    int sign;
    ui * SdK;
    if(SisL) {
        sign = 1;
        SdK = LdK;
    }
    else {
        sign = -1;
        SdK = RdK;
    }
    for(ui i = sidx; i < S.size(); i++) {
        ui v = S[i];
        assert(SdK[v] + k >= ks);
        if(SdK[v] + k == ks) continue;
        bool skip = false;
        for(auto &e : unn_KL) {
            if(LdK[e] + k == ks && M[e][v] != sign) { // e is critical
                skip = true;
                break;
            }
        }
        if(skip) continue;
        for(auto &e : unn_KR) {
            if(RdK[e] + k == ks && M[e][v] != -sign) { // e is critical
                skip = true;
                break;
            }
        }
        if(skip) continue;
        if(v != u) new_S.push_back(v);
    }
}

void vertex_reduction(ui & del_cnt)
{
    ver_del = new bool[n];
    memset(ver_del, 0, sizeof(bool)*n);
    if(tau < k) {cout<<"tau < k"<<endl; exit(1);}
    ui pt = tau - k; //pt >= 0
    ui nt = tau - k + 1; //nt >= 1
    ui dt = 2 * tau - k;
    queue<ui> Q;
    for(ui i = 0; i < n; i++) if(p_degree[i] < pt || n_degree[i] < nt || degree[i] < dt) {
        Q.push(i);
        ver_del[i] = 1;
        ++ del_cnt;
    }
    while (!Q.empty()) {
        ui u = Q.front();
        Q.pop();
        for(ui i = p_pstart[u]; i < p_pstart[u+1]; i++) {
            ui v = p_edges[i];
            if(ver_del[v]) continue;
            assert(p_degree[v] > 0 && degree[v] > 0);
            -- p_degree[v];
            -- degree[v];
            if(p_degree[v] < pt || degree[v] < dt) {
                Q.push(v);
                ver_del[v] = 1;
                ++ del_cnt;
            }
        }
        for(ui i = n_pstart[u]; i < n_pstart[u+1]; i++) {
            ui v = n_edges[i];
            if(ver_del[v]) continue;
            assert(n_degree[v] > 0 && degree[v] > 0);
            -- n_degree[v];
            -- degree[v];
            if(n_degree[v] < nt || degree[v] < dt) {
                Q.push(v);
                ver_del[v] = 1;
                ++ del_cnt;
            }
        }
    }
}

void reconstruct_reduced_graph()
{
    ui * mapping = new ui[n];
    ui idx = 0;
    pm = 0;
    nm = 0;
    original_id = new ui[n];
    for(ui u = 0; u < n; u++) if(!ver_del[u]) {
        mapping[u] = idx;
        original_id[idx++] = u;
        for(ui i = p_pstart[u]; i < p_pstart[u+1]; i ++) if(!ver_del[p_edges[i]]) ++ pm;
        for(ui i = n_pstart[u]; i < n_pstart[u+1]; i ++) if(!ver_del[n_edges[i]]) ++ nm;
    }
    assert(pm%2==0 && nm%2==0);
#ifdef _DEBUG_
    cout<<"mapping:"<<endl; for(ui i = 0; i < n; i++) if(!ver_del[i]) cout<<i<<" -> "<<mapping[i]<<endl;
#endif
    ui original_n = n;
    n = idx;
    pm /= 2;
    nm /= 2;
    t_p_pstart = new ui[n+1];
    t_p_edges = new ui[2*pm];
    t_n_pstart = new ui[n+1];
    t_n_edges = new ui[2*nm];
    t_p_pend = new ui[n+1]; //for later usage
    t_n_pend = new ui[n+1]; //for later usage

    ui new_i = 0;
    t_p_pstart[0] = 0;
    for(ui i = 0; i < original_n; i++) if(!ver_del[i]){
        ui pos = t_p_pstart[new_i];
        for(ui j = p_pstart[i]; j < p_pstart[i+1]; j++) if(!ver_del[p_edges[j]]){
            t_p_edges[pos++] = mapping[p_edges[j]];
        }
        t_p_pstart[++new_i] = pos;
    }
    assert(new_i == n);
    new_i = 0;
    t_n_pstart[0] = 0;
    for(ui i = 0; i < original_n; i++) if(!ver_del[i]){
        ui pos = t_n_pstart[new_i];
        for(ui j = n_pstart[i]; j < n_pstart[i+1]; j++) if(!ver_del[n_edges[j]]){
            t_n_edges[pos++] = mapping[n_edges[j]];
        }
        t_n_pstart[++new_i] = pos;
    }
    assert(new_i == n);
    delete [] p_pstart;
    delete [] p_edges;
    delete [] p_pend;
    delete [] n_pstart;
    delete [] n_edges;
    delete [] n_pend;
    p_pstart = new ui[n+1];
    p_edges = new ui[2*pm];
    p_pend = new ui[n+1];
    n_pstart = new ui[n+1];
    n_edges = new ui[2*nm];
    n_pend = new ui[n+1];
    memcpy(p_pstart, t_p_pstart, sizeof(ui)*(n+1));
    memcpy(p_edges, t_p_edges, sizeof(ui)*(2*pm));
    memcpy(n_pstart, t_n_pstart, sizeof(ui)*(n+1));
    memcpy(n_edges, t_n_edges, sizeof(ui)*(2*nm));
    max_deg = 0;
    for(ui u = 0; u < n; u++) {
        p_pend[u] = p_pstart[u+1];
        n_pend[u] = n_pstart[u+1];
        degree[u] = (p_pstart[u+1] - p_pstart[u]) + (n_pstart[u+1] - n_pstart[u]);
        if(degree[u] > max_deg) max_deg = degree[u];
    }
    delete [] mapping;
#ifdef _DEBUG_
    cout<<"reconstructed graph:"<<endl;
    for(ui i = 0; i < n; i++) {
        cout<<i<<" ("<<degree[i]<<":"<<p_degree[i]<<"+"<<n_degree[i]<<"): "<<endl;
        cout<<"\tpositive neighbors: "; for(ui j = p_pstart[i]; j < p_pstart[i+1]; j++) cout<<p_edges[j]<<","; cout<<endl;
        cout<<"\tnegative neighbors: "; for(ui j = n_pstart[i]; j < n_pstart[i+1]; j++) cout<<n_edges[j]<<","; cout<<endl;
    }
#endif
}

void maximal_balanced_kplex_enum()
{
    CTprune();
    if(n == 0) return;
    process_order = new ui[n];
    for(ui i = 0; i < n; i++) process_order[i] = i;
    comp_process_order();
    ver_rank = new ui[n];
    for(ui i = 0; i < n; i++) ver_rank[process_order[i]] = i;
    inCL = new bool[n]; memset(inCL, 0, sizeof(bool)*n); //CL&CR may overlap
    inCR = new bool[n]; memset(inCR, 0, sizeof(bool)*n);
    inXL = new bool[n]; memset(inXL, 0, sizeof(bool)*n); //XL&XR may overlap
    inXR = new bool[n]; memset(inXR, 0, sizeof(bool)*n);
    inQv = new short[n]; //memset(inQv, 0, sizeof(short)*n);
    memset(ver_del, 0, sizeof(bool)*n); //for other use
    M_len = n;
    M = new short*[M_len];
    for(int i = 0; i < M_len; i++) M[i] = new short[M_len];
    partial_kplex_cap = new short[n];
    nneiKL.resize(n);
    nneiKR.resize(n);
    skipv_CL = new bool[n]; memset(skipv_CL, 0, sizeof(bool)*n);
    skipv_CR = new bool[n]; memset(skipv_CR, 0, sizeof(bool)*n);
    for(ui i = 0; i < n; i++) {
        ui u = process_order[i];
        vector<ui> KL, KR, CL, CR, XL, XR;
        KL.push_back(u);
        obtain_CL_CR_XL_XR(u, CL, CR, XL, XR);
        if(CL.size() + 1 < tau || CR.size() < tau) {
            for(auto e : CL) inCL[e] = 0;
            for(auto e : CR) inCR[e] = 0;
            for(auto e : XL) inXL[e] = 0;
            for(auto e : XR) inXR[e] = 0;
            continue;
        }
        build_M(u, CL, CR, XL, XR);
        sort(CL.begin(), CL.end(), mycomp);
        sort(CR.begin(), CR.end(), mycomp);
        for(auto e : CL) inCL[e] = 0;
        for(auto e : CR) inCR[e] = 0;
        for(auto e : XL) inXL[e] = 0;
        for(auto e : XR) inXR[e] = 0;
        enum_procedure(KL, KR, CL, CR, XL, XR, 1);
    }
}

//no overlap among CL CR P
//no overlapp between (CL+CR+P) and (XL+XR)
//XL and XR overlaps
void obtain_CL_CR_P_XL_XR(ui u, vector<ui>&CL, vector<ui>&CR, vector<ui>&P, vector<ui>&XL, vector<ui>&XR)
{
    //CL,CR
    for(ui j = p_pstart[u]; j < p_pend[u]; j++) {
        ui v = p_edges[j];
        if(ver_rank[v] > ver_rank[u]) {// u's +> neighbor
            CL.push_back(v); part[v] = 1;
        }
    }
    for(ui j = n_pstart[u]; j < n_pend[u]; j++) {
        ui v = n_edges[j];
        if(ver_rank[v] > ver_rank[u]) {// u's -> neighbor
            CR.push_back(v); part[v] = 2;
        }
    }
    //P
    for(ui j = p_pstart[u]; j < p_pend[u]; j++) {
        ui v = p_edges[j];
        if(ver_rank[v] > ver_rank[u]) {
            for(ui h = p_pstart[v]; h < p_pend[v]; h++) {
                ui w = p_edges[h]; //u's ++ nei
                if(ver_rank[w] > ver_rank[u] ) {
                    if(part[w] !=1 && part[w] != 2) {addtoKL[w]=1;}
                    if(part[w] == 0) {P.push_back(w); part[w] = 3;}
                }
                else if (ver_rank[w] < ver_rank[u]) {
                    if(inXL[w] == 0) {
                        XL.push_back(w); inXL[w] = 1;
                    }
                }
            }
            for(ui h = n_pstart[v]; h < n_pend[v]; h++) {
                ui w = n_edges[h]; //u's +- nei
                if(ver_rank[w] > ver_rank[u] ) {
                    if(part[w] !=1 && part[w] != 2) {addtoKR[w]=1;}
                    if(part[w] == 0) {P.push_back(w); part[w] = 3;}
                }
                else if (ver_rank[w] < ver_rank[u]) {
                    if(inXR[w] == 0) {
                        XR.push_back(w); inXR[w] = 1;
                    }
                }
            }
        }
        else {
            assert(ver_rank[v] < ver_rank[u]);
            if(inXL[v] == 0) {
                XL.push_back(v); inXL[v] = 1;
            }
            if(inXR[v] == 0) {
                XR.push_back(v); inXR[v] = 1;
            }
//            for(ui h = p_pstart[v]; h < p_pend[v]; h++) {
//                ui w = p_edges[h];
//                if(ver_rank[w] < ver_rank[u] && inXL[w] == 0) {
//                    XL.push_back(w); inXL[w] = 1;
//                }
//            }
//            for(ui h = n_pstart[v]; h < n_pend[v]; h++) {
//                ui w = n_edges[h];
//                if(ver_rank[w] < ver_rank[u] && inXR[w] == 0) {
//                    XR.push_back(w); inXR[w] = 1;
//                }
//            }
        }
    }
    for(ui j = n_pstart[u]; j < n_pend[u]; j++) {
        ui v = n_edges[j];
        if(ver_rank[v] > ver_rank[u]) {
            for(ui h = p_pstart[v]; h < p_pend[v]; h++) {
                ui w = p_edges[h]; //u's -+ nei
                if(ver_rank[w] > ver_rank[u]) {
                    if(part[w] !=1 && part[w] != 2) {addtoKR[w]=1;}
                    if(part[w] == 0) {P.push_back(w); part[w] = 3;}
                }
                else if (ver_rank[w] < ver_rank[u]) {
                    if(inXR[w] == 0) {
                        XR.push_back(w); inXR[w] = 1;
                    }
                }
            }
            for(ui h = n_pstart[v]; h < n_pend[v]; h++) {
                ui w = n_edges[h]; //u's -- nei
                if(ver_rank[w] > ver_rank[u]) {
                    if(part[w] !=1 && part[w] != 2) {addtoKL[w]=1;}
                    if(part[w] == 0) {P.push_back(w); part[w] = 3;}
                }
                else if (ver_rank[w] < ver_rank[u]) {
                    if(inXL[w] == 0) {
                        XL.push_back(w); inXL[w] = 1;
                    }
                }
            }
        }
        else {
            assert(ver_rank[v] < ver_rank[u]);
            if(inXL[v] == 0) {
                XL.push_back(v); inXL[v] = 1;
            }
            if(inXR[v] == 0) {
                XR.push_back(v); inXR[v] = 1;
            }
//            for(ui h = p_pstart[v]; h < p_pend[v]; h++) {
//                ui w = p_edges[h];
//                if(ver_rank[w] < ver_rank[u] && inXR[w] == 0) {
//                    XR.push_back(w); inXR[w] = 1;
//                }
//            }
//            for(ui h = n_pstart[v]; h < n_pend[v]; h++) {
//                ui w = n_edges[h];
//                if(ver_rank[w] < ver_rank[u] && inXL[w] == 0) {
//                    XL.push_back(w); inXL[w] = 1;
//                }
//            }
        }
    }

#ifdef _CostlyAssert_
    for(auto &e : CL) assert(part[e] == 1);
    for(auto &e : CR) assert(part[e] == 2);
    for(auto &e : P)  assert(part[e] == 3);
#endif
#ifdef _DEBUG_
    cout<<"CL ("<<CL.size()<<"): "; for(auto e : CL) cout<<e<<","; cout<<endl;
    cout<<"CR ("<<CR.size()<<"): "; for(auto e : CR) cout<<e<<","; cout<<endl;
    cout<<"P ("<<P.size()<<"): "; for(auto e : P) cout<<e<<","; cout<<endl;
    cout<<"XL ("<<XL.size()<<"): "; for(auto e : XL) cout<<e<<","; cout<<endl;
    cout<<"XR ("<<XR.size()<<"): "; for(auto e : XR) cout<<e<<","; cout<<endl;
#endif
}

void Cnm(int nn, int outLen, int starIdx, int mm, int * A, int AIdx, vector<vector<ui>> & combs)
{
    if(mm == 0){
        vector<ui> t_vec;
        for(int i = 0; i < outLen; i++){
            t_vec.push_back(A[i]);
        }
        combs.push_back(t_vec);
        return;
    }
    int endIdx = nn-mm+1;
    for(int i = starIdx; i < endIdx; i++){
        A[AIdx] = i;
        Cnm(nn, outLen, i+1, mm-1, A, AIdx+1, combs);
    }
}

void get_combs(ui dom_n, ui dom_m, vector<vector<ui>> & combs)
{
    int * A = new int[dom_n];
    memset(A, 0, sizeof(int)*dom_n);
    Cnm(dom_n, dom_m, 0, dom_m, A, 0, combs);
    delete [] A;
}

void enum_adv(vector<ui>&KL, vector<ui>&KR, vector<ui>&CL, vector<ui>&CR, vector<ui>&XL, vector<ui>&XR)
{
#ifdef _DEBUG_
    cout<<"\t\tin enum()."<<endl;
    cout<<"\t\tKL: "; for(auto e : KL) cout<<e<<","; cout<<endl;
    cout<<"\t\tKR: "; for(auto e : KR) cout<<e<<","; cout<<endl;
    cout<<"\t\tCL: "; for(auto e : CL) cout<<e<<","; cout<<endl;
    cout<<"\t\tCR: "; for(auto e : CR) cout<<e<<","; cout<<endl;
    cout<<"\t\tXL: "; for(auto e : XL) cout<<e<<","; cout<<endl;
    cout<<"\t\tXR: "; for(auto e : XR) cout<<e<<","; cout<<endl;
#endif
    
}

void vertex_reduction_on_initial_subgraph(ui u, vector<ui>&CL, vector<ui>&CR,vector<ui>&P,vector<ui>&XL,vector<ui>&XR)
{
    ui del_cnt = 0;
    ui pt = tau - k; //pt >= 0
    ui nt = tau - k + 1; //nt >= 1
    ui dt = 2 * tau - k;
    queue<ui> Q;
    if(p_degree[u] < pt || n_degree[u] < nt || degree[u] < dt) {
        Q.push(u);
        ver_del[u] = 1;
        ++ del_cnt;
    }
    for(auto &i : CL) if(p_degree[i] < pt || n_degree[i] < nt || degree[i] < dt) {
        Q.push(i);
        ver_del[i] = 1;
        ++ del_cnt;
    }
    for(auto &i : CR) if(p_degree[i] < pt || n_degree[i] < nt || degree[i] < dt) {
        Q.push(i);
        ver_del[i] = 1;
        ++ del_cnt;
    }
    for(auto &i : P) if(p_degree[i] < pt || n_degree[i] < nt || degree[i] < dt) {
        Q.push(i);
        ver_del[i] = 1;
        ++ del_cnt;
    }
    for(auto &i : XL) if(p_degree[i] < pt || n_degree[i] < nt || degree[i] < dt) {
        Q.push(i);
        ver_del[i] = 1;
        ++ del_cnt;
    }
    for(auto &i : XR) if(p_degree[i] < pt || n_degree[i] < nt || degree[i] < dt) {
        if(ver_del[i] == 0){
            Q.push(i);
            ver_del[i] = 1;
            ++ del_cnt;
        }
    }
    while (!Q.empty()) {
        ui u = Q.front();
        Q.pop();
        for(ui i = t_p_pstart[u]; i < t_p_pend[u]; i++) {
            ui v = t_p_edges[i];
            if(ver_del[v]) continue;
            assert(p_degree[v] > 0 && degree[v] > 0);
            -- p_degree[v];
            -- degree[v];
            if(p_degree[v] < pt || degree[v] < dt) {
                Q.push(v);
                ver_del[v] = 1;
                ++ del_cnt;
            }
        }
        for(ui i = t_n_pstart[u]; i < t_n_pend[u]; i++) {
            ui v = t_n_edges[i];
            if(ver_del[v]) continue;
            assert(n_degree[v] > 0 && degree[v] > 0);
            -- n_degree[v];
            -- degree[v];
            if(n_degree[v] < nt || degree[v] < dt) {
                Q.push(v);
                ver_del[v] = 1;
                ++ del_cnt;
            }
        }
    }
    vector<ui> newCL, newCR, newP, newXL, newXR;
    for(auto &e : CL) {
        if(ver_del[e]) {
            ver_del[e] = 0;
            inCL[e] = 0;
        }
        else {
            newCL.push_back(e);
        }
    }
    for(auto &e : CR) {
        if(ver_del[e]) {
            ver_del[e] = 0;
            inCR[e] = 0;
        }
        else {
            newCR.push_back(e);
        }
    }
    for(auto &e : P) {
        if(ver_del[e]) {
            ver_del[e] = 0;
            inP[e] = 0;
        }
        else {
            newP.push_back(e);
        }
    }
    for(auto &e : XL) {
        if(ver_del[e]) {
            inXL[e] = 0;
        }
        else {
            newXL.push_back(e);
        }
    }
    for(auto &e : XR) {
        if(ver_del[e]) {
            inXR[e] = 0;
        }
        else {
            newXR.push_back(e);
        }
    }
    for(auto &e : XL) ver_del[e] = 0;
    for(auto &e : XR) ver_del[e] = 0;
    CL = newCL;
    CR = newCR;
    P = newP;
    XL = newXL;
    XR = newXR;
}

void vertex_reduction_on_initial_subgraph(ui u, vector<ui>&CL, vector<ui>&CR,vector<ui>&P)
{
    ui del_cnt = 0;
    ui pt = tau - k; //pt >= 0
    ui nt = tau - k + 1; //nt >= 1
    ui dt = 2 * tau - k;
    queue<ui> Q;
    if(p_degree[u] < pt || n_degree[u] < nt || degree[u] < dt) {
        Q.push(u);
        ver_del[u] = 1;
        ++ del_cnt;
    }
    for(auto &i : CL) if(p_degree[i] < pt || n_degree[i] < nt || degree[i] < dt) {
        Q.push(i);
        ver_del[i] = 1;
        ++ del_cnt;
    }
    for(auto &i : CR) if(p_degree[i] < pt || n_degree[i] < nt || degree[i] < dt) {
        Q.push(i);
        ver_del[i] = 1;
        ++ del_cnt;
    }
    for(auto &i : P) if(p_degree[i] < pt || n_degree[i] < nt || degree[i] < dt) {
        Q.push(i);
        ver_del[i] = 1;
        ++ del_cnt;
    }
    while (!Q.empty()) {
        ui u = Q.front();
        Q.pop();
        for(ui i = t_p_pstart[u]; i < t_p_pend[u]; i++) {
            ui v = t_p_edges[i];
            if(ver_del[v]) continue;
            assert(p_degree[v] > 0 && degree[v] > 0);
            -- p_degree[v];
            -- degree[v];
            if(p_degree[v] < pt || degree[v] < dt) {
                Q.push(v);
                ver_del[v] = 1;
                ++ del_cnt;
            }
        }
        for(ui i = t_n_pstart[u]; i < t_n_pend[u]; i++) {
            ui v = t_n_edges[i];
            if(ver_del[v]) continue;
            assert(n_degree[v] > 0 && degree[v] > 0);
            -- n_degree[v];
            -- degree[v];
            if(n_degree[v] < nt || degree[v] < dt) {
                Q.push(v);
                ver_del[v] = 1;
                ++ del_cnt;
            }
        }
    }
    vector<ui> newCL, newCR, newP;
    for(auto &e : CL) {
        if(ver_del[e]) {
            ver_del[e] = 0;
            part[e] = 0;
        }
        else {
            newCL.push_back(e);
        }
    }
    for(auto &e : CR) {
        if(ver_del[e]) {
            ver_del[e] = 0;
            part[e] = 0;
        }
        else {
            newCR.push_back(e);
        }
    }
    for(auto &e : P) {
        if(ver_del[e]) {
            ver_del[e] = 0;
            part[e] = 0;
            addtoKL[e] = 0;
            addtoKR[e] = 0;
        }
        else {
            newP.push_back(e);
        }
    }
    CL = newCL;
    CR = newCR;
    P = newP;
}

//actually, KL={u} and KR={}, obtain degree in KLCL and in KRCR
void obtain_pnd_LRside(ui u, vector<ui>&CL, vector<ui>&CR, vector<ui>&P)
{
    assert(part[u] == 0);
    part[u] = 1;
    CL.push_back(u); //temporarily put u into CL, ease for programming
    for(auto &e : CL) {
        for(ui i = t_p_pstart[e]; i < t_p_pend[e]; i++) {
            if(part[t_p_edges[i]]==1) ++ pdL[e];
            else if (part[t_p_edges[i]]==2) ++pdR[e];
        }
        for(ui i = t_n_pstart[e]; i < t_n_pend[e]; i++) {
            if(part[t_n_edges[i]]==1) ++ ndL[e];
            else if (part[t_n_edges[i]]==2) ++ndR[e];
        }
    }
    for(auto &e : CR) {
        for(ui i = t_p_pstart[e]; i < t_p_pend[e]; i++) {
            if(part[t_p_edges[i]]==1) ++ pdL[e];
            else if (part[t_p_edges[i]]==2) ++pdR[e];
        }
        for(ui i = t_n_pstart[e]; i < t_n_pend[e]; i++) {
            if(part[t_n_edges[i]]==1) ++ ndL[e];
            else if (part[t_n_edges[i]]==2) ++ndR[e];
        }
    }
    for(auto &e : P) {
        for(ui i = t_p_pstart[e]; i < t_p_pend[e]; i++) {
            if(part[t_p_edges[i]]==1) ++ pdL[e];
            else if (part[t_p_edges[i]]==2) ++pdR[e];
        }
        for(ui i = t_n_pstart[e]; i < t_n_pend[e]; i++) {
            if(part[t_n_edges[i]]==1) ++ ndL[e];
            else if (part[t_n_edges[i]]==2) ++ndR[e];
        }
    }
    assert(CL.back() == u);
    CL.pop_back();
    part[u] = 0;
#ifdef _DEBUG_
    cout<<"pn degree LRside :"<<endl;
    cout<<"u : pdL= "<<pdL[u]<<", pdR= "<<pdR[u]<<", ndL= "<<ndL[u]<<", ndR= "<<ndR[u]<<endl;
    cout<<"CL : "<<endl; for(auto &e : CL) cout<<"\t"<<e<<" : pdL= "<<pdL[e]<<", pdR= "<<pdR[e]<<", ndL= "<<ndL[e]<<", ndR= "<<ndR[e]<<endl;
    cout<<"CR : "<<endl; for(auto &e : CR) cout<<"\t"<<e<<" : pdL= "<<pdL[e]<<", pdR= "<<pdR[e]<<", ndL= "<<ndL[e]<<", ndR= "<<ndR[e]<<endl;
    cout<<"P : "<<endl; for(auto &e : P) cout<<"\t"<<e<<" : pdL= "<<pdL[e]<<", pdR= "<<pdR[e]<<", ndL= "<<ndL[e]<<", ndR= "<<ndR[e]<<endl;
#endif
}

void obtain_pnd_LRside(ui u, vector<ui>&CL, vector<ui>&CR, vector<ui>&P, vector<ui>&XL, vector<ui>&XR)
{
    assert(!inCL[u]);
    inCL[u] = 1;
    CL.push_back(u); //temporarily put u into CL, ease for programming
    for(auto &e : CL) {
        for(ui i = t_p_pstart[e]; i < t_p_pend[e]; i++) {
            if(inCL[t_p_edges[i]]) ++ pdL[e];
            else if (inCR[t_p_edges[i]]) ++pdR[e];
        }
        for(ui i = t_n_pstart[e]; i < t_n_pend[e]; i++) {
            if(inCL[t_n_edges[i]]) ++ ndL[e];
            else if (inCR[t_n_edges[i]]) ++ndR[e];
        }
    }
    for(auto &e : CR) {
        for(ui i = t_p_pstart[e]; i < t_p_pend[e]; i++) {
            if(inCL[t_p_edges[i]]) ++ pdL[e];
            else if (inCR[t_p_edges[i]]) ++pdR[e];
        }
        for(ui i = t_n_pstart[e]; i < t_n_pend[e]; i++) {
            if(inCL[t_n_edges[i]]) ++ ndL[e];
            else if (inCR[t_n_edges[i]]) ++ndR[e];
        }
    }
    for(auto &e : P) {
        for(ui i = t_p_pstart[e]; i < t_p_pend[e]; i++) {
            if(inCL[t_p_edges[i]]) ++ pdL[e];
            else if (inCR[t_p_edges[i]]) ++pdR[e];
        }
        for(ui i = t_n_pstart[e]; i < t_n_pend[e]; i++) {
            if(inCL[t_n_edges[i]]) ++ ndL[e];
            else if (inCR[t_n_edges[i]]) ++ndR[e];
        }
    }
    for(auto &e : XL) {
        for(ui i = t_p_pstart[e]; i < t_p_pend[e]; i++) {
            if(inCL[t_p_edges[i]]) ++ XL_pdL[e];
            else if (inCR[t_p_edges[i]]) ++XL_pdR[e];
        }
        for(ui i = t_n_pstart[e]; i < t_n_pend[e]; i++) {
            if(inCL[t_n_edges[i]]) ++ XL_ndL[e];
            else if (inCR[t_n_edges[i]]) ++XL_ndR[e];
        }
    }
    for(auto &e : XR) {
        for(ui i = t_p_pstart[e]; i < t_p_pend[e]; i++) {
            if(inCL[t_p_edges[i]]) ++ XR_pdL[e];
            else if (inCR[t_p_edges[i]]) ++XR_pdR[e];
        }
        for(ui i = t_n_pstart[e]; i < t_n_pend[e]; i++) {
            if(inCL[t_n_edges[i]]) ++ XR_ndL[e];
            else if (inCR[t_n_edges[i]]) ++XR_ndR[e];
        }
    }
    assert(CL.back() == u);
    CL.pop_back();
    inCL[u] = 0;
#ifdef _DEBUG_
    cout<<"pn degree LRside :"<<endl;
    cout<<"u : pdL= "<<pdL[u]<<", pdR= "<<pdR[u]<<", ndL= "<<ndL[u]<<", ndR= "<<ndR[u]<<endl;
    cout<<"CL : "<<endl; for(auto &e : CL) cout<<"\t"<<e<<" : pdL= "<<pdL[e]<<", pdR= "<<pdR[e]<<", ndL= "<<ndL[e]<<", ndR= "<<ndR[e]<<endl;
    cout<<"CR : "<<endl; for(auto &e : CR) cout<<"\t"<<e<<" : pdL= "<<pdL[e]<<", pdR= "<<pdR[e]<<", ndL= "<<ndL[e]<<", ndR= "<<ndR[e]<<endl;
    cout<<"P : "<<endl; for(auto &e : P) cout<<"\t"<<e<<" : pdL= "<<pdL[e]<<", pdR= "<<pdR[e]<<", ndL= "<<ndL[e]<<", ndR= "<<ndR[e]<<endl;
#endif
}

void update_pnd_LRside(vector<ui>&SL, vector<ui>&SR, vector<ui>&ST, vector<ui>&SD, int f)
{
    if(f==0) { //make change
        for(auto &e : SL) {
            for(ui i = t_p_pstart[e]; i < t_p_pend[e]; i++) {
                ui u = t_p_edges[i];
                if(part[u]==1 || part[u]==2){
                    ++pdL[u];
                    assert(pdR[u]>0);
                    --pdR[u];
                }
            }
            for(ui i = t_n_pstart[e]; i < t_n_pend[e]; i++) {
                ui u = t_n_edges[i];
                if(part[u]==1 || part[u]==2){
                    ++ndL[u];
                    assert(ndR[u]>0);
                    --ndR[u];
                }
            }
        }
        for(auto &e : SR) {
            for(ui i = t_p_pstart[e]; i < t_p_pend[e]; i++) {
                ui u = t_p_edges[i];
                if(part[u]==1 || part[u]==2){
                    ++pdR[u];
                    assert(pdL[u]>0);
                    --pdL[u];
                }
            }
            for(ui i = t_n_pstart[e]; i < t_n_pend[e]; i++) {
                ui u = t_n_edges[i];
                if(part[u]==1 || part[u]==2){
                    ++ndR[u];
                    assert(ndL[u]>0);
                    --ndL[u];
                }
            }
        }
        for(auto &e : ST) {
            for(ui i = t_p_pstart[e]; i < t_p_pend[e]; i++) {
                ui u = t_p_edges[i];
                if(part[u]==1 || part[u]==2){
                    ++pdL[u];
                }
            }
            for(ui i = t_n_pstart[e]; i < t_n_pend[e]; i++) {
                ui u = t_n_edges[i];
                if(part[u]==1 || part[u]==2){
                    ++ndL[u];
                }
            }
        }
        for(auto &e : SD) {
            for(ui i = t_p_pstart[e]; i < t_p_pend[e]; i++) {
                ui u = t_p_edges[i];
                if(part[u]==1 || part[u]==2){
                    ++pdR[u];
                }
            }
            for(ui i = t_n_pstart[e]; i < t_n_pend[e]; i++) {
                ui u = t_n_edges[i];
                if(part[u]==1 || part[u]==2){
                    ++ndR[u];
                }
            }
        }
    }
    else { //restore
        for(auto &e : SL) {
            for(ui i = t_p_pstart[e]; i < t_p_pend[e]; i++) {
                ui u = t_p_edges[i];
                if(part[u]==1 || part[u]==2){
                    ++pdR[u];
                    assert(pdL[u]>0);
                    --pdL[u];
                }
            }
            for(ui i = t_n_pstart[e]; i < t_n_pend[e]; i++) {
                ui u = t_n_edges[i];
                if(part[u]==1 || part[u]==2){
                    ++ndR[u];
                    assert(ndL[u]>0);
                    --ndL[u];
                }
            }
        }
        for(auto &e : SR) {
            for(ui i = t_p_pstart[e]; i < t_p_pend[e]; i++) {
                ui u = t_p_edges[i];
                if(part[u]==1 || part[u]==2){
                    ++pdL[u];
                    assert(pdR[u]>0);
                    --pdR[u];
                }
            }
            for(ui i = t_n_pstart[e]; i < t_n_pend[e]; i++) {
                ui u = t_n_edges[i];
                if(part[u]==1 || part[u]==2){
                    ++ndL[u];
                    assert(ndR[u]>0);
                    --ndR[u];
                }
            }
        }
        for(auto &e : ST) {
            for(ui i = t_p_pstart[e]; i < t_p_pend[e]; i++) {
                ui u = t_p_edges[i];
                if(part[u]==1 || part[u]==2){
                    --pdL[u];
                }
            }
            for(ui i = t_n_pstart[e]; i < t_n_pend[e]; i++) {
                ui u = t_n_edges[i];
                if(part[u]==1 || part[u]==2){
                    --ndL[u];
                }
            }
        }
        for(auto &e : SD) {
            for(ui i = t_p_pstart[e]; i < t_p_pend[e]; i++) {
                ui u = t_p_edges[i];
                if(part[u]==1 || part[u]==2){
                    --pdR[u];
                }
            }
            for(ui i = t_n_pstart[e]; i < t_n_pend[e]; i++) {
                ui u = t_n_edges[i];
                if(part[u]==1 || part[u]==2){
                    --ndR[u];
                }
            }
        }
    }
}

bool refine_LRP_by_pndLR(ui u, vector<ui>&CL, vector<ui>&CR, vector<ui>&P)
{
    assert(part[u]==0);
    part[u] = 1;
    CL.push_back(u); //temporarily put u into CL, ease for programming
    
    queue<ui> Q;
    for(auto &e : CL) {
        if( ((pdL[e]<tau-2*k+1 || ndR[e]<tau-2*k+2) && (pdR[e]<tau-2*k+3 || ndL[e]<tau-2*k+3)) || (pdL[e]+ndR[e]<2*tau-2*k+1 && pdR[e]+ndL[e]<2*tau-2*k+2) ) {
            Q.push(e);
            ver_del[e] = 1;
        }
    }
    //
    for(auto &e : CR) {
        if( (( pdR[e]<tau-2*k+1 || ndL[e]<tau-2*k+2 ) && ( pdL[e]<tau-2*k+2 || ndR[e]<tau-2*k+4 )) || ( pdR[e]+ndL[e]<2*tau-2*k+1 && pdL[e]+ndR[e]<2*tau-2*k+2 ) ) {
            Q.push(e);
            ver_del[e] = 1;
        }
    }
    for(auto &e : P) {
        if( (( pdL[e]<tau-2*k+2 || ndR[e]<tau-2*k+4 ) && ( pdR[e]<tau-2*k+3 || ndL[e]<tau-2*k+3 )) || ( pdL[e]+ndR[e]<2*tau-2*k+2 && pdR[e]+ndL[e]<2*tau-2*k+2 ) ) {
            Q.push(e);
            ver_del[e] = 1;
        }
    }
    while (!Q.empty()) {
        ui v = Q.front();
        Q.pop();
        if (part[v] == 1) {
            for(ui i = t_p_pstart[v]; i < t_p_pend[v]; i ++) {
                ui w = t_p_edges[i];
                if(part[w]==1){
                    assert(pdL[w]>0);
                    --pdL[w];
                    if(ver_del[w]==1) continue;
                    if( ((pdL[w]<tau-2*k+1 || ndR[w]<tau-2*k+1) && (pdR[w]<tau-2*k+3 || ndL[w]<tau-2*k+3)) || (pdL[w]+ndR[w]<2*tau-2*k+1 && pdR[w]+ndL[w]<2*tau-2*k+2) ) {
                        Q.push(w);
                        ver_del[w] = 1;
                    }
                }
                else if (part[w]==2) {
                    assert(pdL[w]>0);
                    --pdL[w];
                    if(ver_del[w]==1) continue;
                    if( (( pdR[w]<tau-2*k+1 || ndL[w]<tau-2*k+2 ) && ( pdL[w]<tau-2*k+2 || ndR[w]<tau-2*k+4 )) || ( pdR[w]+ndL[w]<2*tau-2*k+1 && pdL[w]+ndR[w]<2*tau-2*k+2 ) ) {
                        Q.push(w);
                        ver_del[w] = 1;
                    }
                }
                else if (part[w]==3) {
                    assert(pdL[w]>0);
                    --pdL[w];
                    if(ver_del[w]==1) continue;
                    if( (( pdL[w]<tau-2*k+2 || ndR[w]<tau-2*k+4 ) && ( pdR[w]<tau-2*k+3 || ndL[w]<tau-2*k+3 )) || ( pdL[w]+ndR[w]<2*tau-2*k+2 && pdR[w]+ndL[w]<2*tau-2*k+2 ) ) {
                        Q.push(w);
                        ver_del[w] = 1;
                    }
                }
            }
            for(ui i = t_n_pstart[v]; i < t_n_pend[v]; i ++) {
                ui w = t_n_edges[i];
                if(part[w]==1){
                    assert(ndL[w]>0);
                    --ndL[w];
                    if(ver_del[w]==1) continue;
                    if( ((pdL[w]<tau-2*k+1 || ndR[w]<tau-2*k+1) && (pdR[w]<tau-2*k+3 || ndL[w]<tau-2*k+3)) || (pdL[w]+ndR[w]<2*tau-2*k+1 && pdR[w]+ndL[w]<2*tau-2*k+2) ) {
                        Q.push(w);
                        ver_del[w] = 1;
                    }
                }
                else if (part[w]==2) {
                    assert(ndL[w]>0);
                    --ndL[w];
                    if(ver_del[w]==1) continue;
                    if( (( pdR[w]<tau-2*k+1 || ndL[w]<tau-2*k+2 ) && ( pdL[w]<tau-2*k+2 || ndR[w]<tau-2*k+4 )) || ( pdR[w]+ndL[w]<2*tau-2*k+1 && pdL[w]+ndR[w]<2*tau-2*k+2 ) ) {
                        Q.push(w);
                        ver_del[w] = 1;
                    }
                }
                else if (part[w]==3) {
                    assert(ndL[w]>0);
                    --ndL[w];
                    if(ver_del[w]==1) continue;
                    if( (( pdL[w]<tau-2*k+2 || ndR[w]<tau-2*k+4 ) && ( pdR[w]<tau-2*k+3 || ndL[w]<tau-2*k+3 )) || ( pdL[w]+ndR[w]<2*tau-2*k+2 && pdR[w]+ndL[w]<2*tau-2*k+2 ) ) {
                        Q.push(w);
                        ver_del[w] = 1;
                    }
                }
            }
        }
        else if (part[v]==2) {
            for(ui i = t_p_pstart[v]; i < t_p_pend[v]; i ++) {
                ui w = t_p_edges[i];
                if(part[w]==1){
                    assert(pdR[w]>0);
                    --pdR[w];
                    if(ver_del[w]==1) continue;
                    if( ((pdL[w]<tau-2*k+1 || ndR[w]<tau-2*k+1) && (pdR[w]<tau-2*k+3 || ndL[w]<tau-2*k+3)) || (pdL[w]+ndR[w]<2*tau-2*k+1 && pdR[w]+ndL[w]<2*tau-2*k+2) ) {
                        Q.push(w);
                        ver_del[w] = 1;
                    }
                }
                else if (part[w]==2) {
                    assert(pdR[w]>0);
                    --pdR[w];
                    if(ver_del[w]==1) continue;
                    if( (( pdR[w]<tau-2*k+1 || ndL[w]<tau-2*k+2 ) && ( pdL[w]<tau-2*k+2 || ndR[w]<tau-2*k+4 )) || ( pdR[w]+ndL[w]<2*tau-2*k+1 && pdL[w]+ndR[w]<2*tau-2*k+2 ) ) {
                        Q.push(w);
                        ver_del[w] = 1;
                    }
                }
                else if (part[w]==3) {
                    assert(pdR[w]>0);
                    --pdR[w];
                    if(ver_del[w]==1) continue;
                    if( (( pdL[w]<tau-2*k+2 || ndR[w]<tau-2*k+4 ) && ( pdR[w]<tau-2*k+3 || ndL[w]<tau-2*k+3 )) || ( pdL[w]+ndR[w]<2*tau-2*k+2 && pdR[w]+ndL[w]<2*tau-2*k+2 ) ) {
                        Q.push(w);
                        ver_del[w] = 1;
                    }
                }
            }
            for(ui i = t_n_pstart[v]; i < t_n_pend[v]; i ++) {
                ui w = t_n_edges[i];
                if(part[w]==1){
                    assert(ndR[w]>0);
                    --ndR[w];
                    if(ver_del[w]==1) continue;
                    if( ((pdL[w]<tau-2*k+1 || ndR[w]<tau-2*k+1) && (pdR[w]<tau-2*k+3 || ndL[w]<tau-2*k+3)) || (pdL[w]+ndR[w]<2*tau-2*k+1 && pdR[w]+ndL[w]<2*tau-2*k+2) ) {
                        Q.push(w);
                        ver_del[w] = 1;
                    }
                }
                else if (part[w]==2) {
                    assert(ndR[w]>0);
                    --ndR[w];
                    if(ver_del[w]==1) continue;
                    if( (( pdR[w]<tau-2*k+1 || ndL[w]<tau-2*k+2 ) && ( pdL[w]<tau-2*k+2 || ndR[w]<tau-2*k+4 )) || ( pdR[w]+ndL[w]<2*tau-2*k+1 && pdL[w]+ndR[w]<2*tau-2*k+2 ) ) {
                        Q.push(w);
                        ver_del[w] = 1;
                    }
                }
                else if (part[w]==3) {
                    assert(ndR[w]>0);
                    --ndR[w];
                    if(ver_del[w]==1) continue;
                    if( (( pdL[w]<tau-2*k+2 || ndR[w]<tau-2*k+4 ) && ( pdR[w]<tau-2*k+3 || ndL[w]<tau-2*k+3 )) || ( pdL[w]+ndR[w]<2*tau-2*k+2 && pdR[w]+ndL[w]<2*tau-2*k+2 ) ) {
                        Q.push(w);
                        ver_del[w] = 1;
                    }
                }
            }
        }
        else {
            assert(part[v]==3);
        }
    }
    assert(CL.back()==u);
    CL.pop_back();
    part[u] = 0;
    
    vector<ui> newCL, newCR, newP;
    for(auto &e : CL) {
        if(ver_del[e]) {
            ver_del[e] = 0;
            part[e] = 0;
            pdL[e]=0; pdR[e]=0; ndL[e]=0; ndR[e]=0;
        }
        else {
            newCL.push_back(e);
        }
    }
    for(auto &e : CR) {
        if(ver_del[e]) {
            ver_del[e] = 0;
            part[e] = 0;
            pdL[e]=0; pdR[e]=0; ndL[e]=0; ndR[e]=0;
        }
        else {
            newCR.push_back(e);
        }
    }
    for(auto &e : P) {
        if(ver_del[e]) {
            ver_del[e] = 0;
            part[e] = 0;
            pdL[e]=0; pdR[e]=0; ndL[e]=0; ndR[e]=0;
            addtoKL[e] = 0;
            addtoKR[e] = 0;
        }
        else {
            newP.push_back(e);
        }
    }
    CL = newCL;
    CR = newCR;
    P = newP;
    if(ver_del[u]) {
        ver_del[u] = 0;
        pdL[u] = 0; pdR[u] = 0; ndL[u] = 0; ndR[u] = 0;
        return true;
    }
    else return false;
}

void refine_LRP_by_pndLR(ui u, vector<ui>&CL, vector<ui>&CR, vector<ui>&P, vector<ui>&XL, vector<ui>&XR)
{
    assert(!inCL[u]);
    inCL[u] = 1;
    CL.push_back(u); //temporarily put u into CL, ease for programming
    
    queue<ui> Q;
    for(auto &e : CL) {
        if( ((pdL[e]<tau-2*k+1 || ndR[e]<tau-2*k+2) && (pdR[e]<tau-2*k+3 || ndL[e]<tau-2*k+3)) || (pdL[e]+ndR[e]<2*tau-2*k+1 && pdR[e]+ndL[e]<2*tau-2*k+2) ) {
            Q.push(e);
            ver_del[e] = 1;
        }
    }
    //
    for(auto &e : CR) {
        if( (( pdR[e]<tau-2*k+1 || ndL[e]<tau-2*k+2 ) && ( pdL[e]<tau-2*k+2 || ndR[e]<tau-2*k+4 )) || ( pdR[e]+ndL[e]<2*tau-2*k+1 && pdL[e]+ndR[e]<2*tau-2*k+2 ) ) {
            Q.push(e);
            ver_del[e] = 1;
        }
    }
//    //
    for(auto &e : P) {
        if( (( pdL[e]<tau-2*k+2 || ndR[e]<tau-2*k+4 ) && ( pdR[e]<tau-2*k+3 || ndL[e]<tau-2*k+3 )) || ( pdL[e]+ndR[e]<2*tau-2*k+2 && pdR[e]+ndL[e]<2*tau-2*k+2 ) ) {
            Q.push(e);
            ver_del[e] = 1;
        }
    }
    //
    for(auto &e : XL) {
        if(XL_pdL[e]<tau-2*k+1 || XL_ndR[e]<tau-2*k+2 || XL_pdL[e]+XL_ndR[e]<2*tau-2*k+1) {
            Q.push(e);
            ver_del[e] = 1;
        }
    }
//
    for(auto &e : XR) {
        if(XR_pdR[e]<tau-2*k+1 || XR_ndL[e]<tau-2*k+2 || XR_pdR[e]+XR_ndL[e]<2*tau-2*k+1) {
            if(ver_del[e] == 0) {
                Q.push(e);
                ver_del[e] = 1;
            }
        }
    }
    while (!Q.empty()) {
        ui v = Q.front();
        Q.pop();
        if (inCL[v]) {
            for(ui i = t_p_pstart[v]; i < t_p_pend[v]; i ++) {
                ui w = t_p_edges[i];
                if(inCL[w]){
                    assert(pdL[w]>0);
                    --pdL[w];
                    if(ver_del[w]==1) continue;
                    if( ((pdL[w]<tau-2*k+1 || ndR[w]<tau-2*k+1) && (pdR[w]<tau-2*k+3 || ndL[w]<tau-2*k+3)) || (pdL[w]+ndR[w]<2*tau-2*k+1 && pdR[w]+ndL[w]<2*tau-2*k+2) ) {
                        Q.push(w);
                        ver_del[w] = 1;
                    }
                }
                else if (inCR[w]) {
                    assert(pdL[w]>0);
                    --pdL[w];
                    if(ver_del[w]==1) continue;
                    if( (( pdR[w]<tau-2*k+1 || ndL[w]<tau-2*k+2 ) && ( pdL[w]<tau-2*k+2 || ndR[w]<tau-2*k+4 )) || ( pdR[w]+ndL[w]<2*tau-2*k+1 && pdL[w]+ndR[w]<2*tau-2*k+2 ) ) {
                        Q.push(w);
                        ver_del[w] = 1;
                    }
                }
                else if (inP[w]) {
                    assert(pdL[w]>0);
                    --pdL[w];
                    if(ver_del[w]==1) continue;
                    if( (( pdL[w]<tau-2*k+2 || ndR[w]<tau-2*k+4 ) && ( pdR[w]<tau-2*k+3 || ndL[w]<tau-2*k+3 )) || ( pdL[w]+ndR[w]<2*tau-2*k+2 && pdR[w]+ndL[w]<2*tau-2*k+2 ) ) {
                        Q.push(w);
                        ver_del[w] = 1;
                    }
                }
            }
            for(ui i = t_n_pstart[v]; i < t_n_pend[v]; i ++) {
                ui w = t_n_edges[i];
                if(inCL[w]){
                    assert(ndL[w]>0);
                    --ndL[w];
                    if(ver_del[w]==1) continue;
                    if( ((pdL[w]<tau-2*k+1 || ndR[w]<tau-2*k+1) && (pdR[w]<tau-2*k+3 || ndL[w]<tau-2*k+3)) || (pdL[w]+ndR[w]<2*tau-2*k+1 && pdR[w]+ndL[w]<2*tau-2*k+2) ) {
                        Q.push(w);
                        ver_del[w] = 1;
                    }
                }
                else if (inCR[w]) {
                    assert(ndL[w]>0);
                    --ndL[w];
                    if(ver_del[w]==1) continue;
                    if( (( pdR[w]<tau-2*k+1 || ndL[w]<tau-2*k+2 ) && ( pdL[w]<tau-2*k+2 || ndR[w]<tau-2*k+4 )) || ( pdR[w]+ndL[w]<2*tau-2*k+1 && pdL[w]+ndR[w]<2*tau-2*k+2 ) ) {
                        Q.push(w);
                        ver_del[w] = 1;
                    }
                }
                else if (inP[w]) {
                    assert(ndL[w]>0);
                    --ndL[w];
                    if(ver_del[w]==1) continue;
                    if( (( pdL[w]<tau-2*k+2 || ndR[w]<tau-2*k+4 ) && ( pdR[w]<tau-2*k+3 || ndL[w]<tau-2*k+3 )) || ( pdL[w]+ndR[w]<2*tau-2*k+2 && pdR[w]+ndL[w]<2*tau-2*k+2 ) ) {
                        Q.push(w);
                        ver_del[w] = 1;
                    }
                }
            }
        }
        else if (inCR[v]) {
            for(ui i = t_p_pstart[v]; i < t_p_pend[v]; i ++) {
                ui w = t_p_edges[i];
                if(inCL[w]){
                    assert(pdR[w]>0);
                    --pdR[w];
                    if(ver_del[w]==1) continue;
                    if( ((pdL[w]<tau-2*k+1 || ndR[w]<tau-2*k+1) && (pdR[w]<tau-2*k+3 || ndL[w]<tau-2*k+3)) || (pdL[w]+ndR[w]<2*tau-2*k+1 && pdR[w]+ndL[w]<2*tau-2*k+2) ) {
                        Q.push(w);
                        ver_del[w] = 1;
                    }
                }
                else if (inCR[w]) {
                    assert(pdR[w]>0);
                    --pdR[w];
                    if(ver_del[w]==1) continue;
                    if( (( pdR[w]<tau-2*k+1 || ndL[w]<tau-2*k+2 ) && ( pdL[w]<tau-2*k+2 || ndR[w]<tau-2*k+4 )) || ( pdR[w]+ndL[w]<2*tau-2*k+1 && pdL[w]+ndR[w]<2*tau-2*k+2 ) ) {
                        Q.push(w);
                        ver_del[w] = 1;
                    }
                }
                else if (inP[w]) {
                    assert(pdR[w]>0);
                    --pdR[w];
                    if(ver_del[w]==1) continue;
                    if( (( pdL[w]<tau-2*k+2 || ndR[w]<tau-2*k+4 ) && ( pdR[w]<tau-2*k+3 || ndL[w]<tau-2*k+3 )) || ( pdL[w]+ndR[w]<2*tau-2*k+2 && pdR[w]+ndL[w]<2*tau-2*k+2 ) ) {
                        Q.push(w);
                        ver_del[w] = 1;
                    }
                }
            }
            for(ui i = t_n_pstart[v]; i < t_n_pend[v]; i ++) {
                ui w = t_n_edges[i];
                if(inCL[w]){
                    assert(ndR[w]>0);
                    --ndR[w];
                    if(ver_del[w]==1) continue;
                    if( ((pdL[w]<tau-2*k+1 || ndR[w]<tau-2*k+1) && (pdR[w]<tau-2*k+3 || ndL[w]<tau-2*k+3)) || (pdL[w]+ndR[w]<2*tau-2*k+1 && pdR[w]+ndL[w]<2*tau-2*k+2) ) {
                        Q.push(w);
                        ver_del[w] = 1;
                    }
                }
                else if (inCR[w]) {
                    assert(ndR[w]>0);
                    --ndR[w];
                    if(ver_del[w]==1) continue;
                    if( (( pdR[w]<tau-2*k+1 || ndL[w]<tau-2*k+2 ) && ( pdL[w]<tau-2*k+2 || ndR[w]<tau-2*k+4 )) || ( pdR[w]+ndL[w]<2*tau-2*k+1 && pdL[w]+ndR[w]<2*tau-2*k+2 ) ) {
                        Q.push(w);
                        ver_del[w] = 1;
                    }
                }
                else if (inP[w]) {
                    assert(ndR[w]>0);
                    --ndR[w];
                    if(ver_del[w]==1) continue;
                    if( (( pdL[w]<tau-2*k+2 || ndR[w]<tau-2*k+4 ) && ( pdR[w]<tau-2*k+3 || ndL[w]<tau-2*k+3 )) || ( pdL[w]+ndR[w]<2*tau-2*k+2 && pdR[w]+ndL[w]<2*tau-2*k+2 ) ) {
                        Q.push(w);
                        ver_del[w] = 1;
                    }
                }
            }
        }
        else {
            assert(inXL[v] || inXR[v] || inP[v]);
        }
    }
    assert(CL.back()==u);
    CL.pop_back();
    inCL[u] = 0;
    
    vector<ui> newCL, newCR, newP, newXL, newXR;
    for(auto &e : CL) {
        if(ver_del[e]) {
            ver_del[e] = 0;
            inCL[e] = 0;
            pdL[e]=0; pdR[e]=0; ndL[e]=0; ndR[e]=0;
        }
        else {
            newCL.push_back(e);
        }
    }
    for(auto &e : CR) {
        if(ver_del[e]) {
            ver_del[e] = 0;
            inCR[e] = 0;
            pdL[e]=0; pdR[e]=0; ndL[e]=0; ndR[e]=0;
        }
        else {
            newCR.push_back(e);
        }
    }
    for(auto &e : P) {
        if(ver_del[e]) {
            ver_del[e] = 0;
            inP[e] = 0;
            pdL[e]=0; pdR[e]=0; ndL[e]=0; ndR[e]=0;
        }
        else {
            newP.push_back(e);
        }
    }
    for(auto &e : XL) {
        if(ver_del[e]) {
            inXL[e] = 0;
            XL_pdL[e]=0; XL_pdR[e]=0; XL_ndL[e]=0; XL_ndR[e]=0;
        }
        else {
            newXL.push_back(e);
        }
    }
    for(auto &e : XR) {
        if(ver_del[e]) {
            inXR[e] = 0;
            XR_pdL[e]=0; XR_pdR[e]=0; XR_ndL[e]=0; XR_ndR[e]=0;
        }
        else {
            newXR.push_back(e);
        }
    }
    for(auto &e : XL) ver_del[e] = 0;
    for(auto &e : XR) ver_del[e] = 0;
    CL = newCL;
    CR = newCR;
    P = newP;
    XL = newXL;
    XR = newXR;
}

void refine_LR_by_pndLR(vector<ui>&newCL, vector<ui>&newCR, vector<ui>&rfL, vector<ui>&rfR, vector<ui>&rfCL, vector<ui>&rfCR)
{
    queue<ui> Q;
    for(auto &e : newCL) if(pdL[e] < tau - k || ndR[e] < tau - k + 1 || pdL[e] + ndR[e] < 2*tau - k) {
        Q.push(e);
        ver_del[e] = 1;
    }
    for(auto &e : newCR) if(pdR[e] < tau - k || ndL[e] < tau - k + 1 || pdR[e] + ndL[e] < 2*tau - k) {
        Q.push(e);
        ver_del[e] = 1;
    }

    while (!Q.empty()) {
        ui v = Q.front();
        Q.pop();
        if(part[v]==1){
            rfL.push_back(v);
            for(ui i = t_p_pstart[v]; i < t_p_pend[v]; i++) {
                ui w = t_p_edges[i];
                if(part[w]==1) {
                    assert(pdL[w] > 0);
                    --pdL[w];
                    if(ver_del[w] == 1) continue;
                    if(pdL[w] < tau - k || ndR[w] < tau - k + 1 || pdL[w] + ndR[w] < 2*tau - k) {
                        Q.push(w);
                        ver_del[w] = 1;
                    }
                }
                else if(part[w]==2) {
                    assert(pdL[w] > 0);
                    --pdL[w];
                }
            }
            for(ui i = t_n_pstart[v]; i < t_n_pend[v]; i++) {
                ui w = t_n_edges[i];
                if(part[w]==1) {
                    assert(ndL[w] > 0);
                    --ndL[w];
                }
                else if(part[w]==2) {
                    assert(ndL[w] > 0);
                    --ndL[w];
                    if(ver_del[w] == 1) continue;
                    if(pdR[w] < tau - k || ndL[w] < tau - k + 1 || pdR[w] + ndL[w] < 2*tau - k) {
                        Q.push(w);
                        ver_del[w] = 1;
                    }
                }
            }
        }
        else {
            assert(part[v]==2);
            rfR.push_back(v);
            for(ui i = t_p_pstart[v]; i < t_p_pend[v]; i++) {
                ui w = t_p_edges[i];
                if(part[w]==1) {
                    assert(pdR[w] > 0);
                    --pdR[w];
                }
                else if(part[w]==2) {
                    assert(pdR[w] > 0);
                    --pdR[w];
                    if(ver_del[w] == 1) continue;
                    if(pdR[w] < tau - k || ndL[w] < tau - k + 1 || pdR[w] + ndL[w] < 2*tau - k) {
                        Q.push(w);
                        ver_del[w] = 1;
                    }
                }
            }
            for(ui i = t_n_pstart[v]; i < t_n_pend[v]; i++) {
                ui w = t_n_edges[i];
                if(part[w]==1) {
                    assert(ndR[w] > 0);
                    --ndR[w];
                    if(ver_del[w] == 1) continue;
                    if(pdL[w] < tau - k || ndR[w] < tau - k + 1 || pdL[w] + ndR[w] < 2*tau - k) {
                        Q.push(w);
                        ver_del[w] = 1;
                    }
                }
                else if(part[w]==2) {
                    assert(ndR[w] > 0);
                    --ndR[w];
                }
            }
        }
    }
    for(ui i = 0; i < newCL.size(); ) {
        if(ver_del[newCL[i]]) {
            rfCL.push_back(newCL[i]);
            newCL[i] = newCL.back();
            newCL.pop_back();
        }
        else{
            ++ i;
        }
    }
    for(ui i = 0; i < newCR.size(); ) {
        if(ver_del[newCR[i]]) {
            rfCR.push_back(newCR[i]);
            newCR[i] = newCR.back();
            newCR.pop_back();
        }
        else{
            ++ i;
        }
    }
    
}

void recover_pndLR_by_rfLrfR(vector<ui>&rfL, vector<ui>&rfR)
{
    for(auto &e : rfL) {
        assert(ver_del[e] == 1);
        ver_del[e] = 0;
        for(ui i = t_p_pstart[e]; i < t_p_pend[e]; i++) {
            ui u = t_p_edges[i];
            if(part[u]==1 || part[u]==2) {
                ++ pdL[u];
            }
        }
        for(ui i = t_n_pstart[e]; i < t_n_pend[e]; i++) {
            ui u = t_n_edges[i];
            if(part[u]==1 || part[u]==2) {
                ++ ndL[u];
            }
        }
    }
    for(auto &e : rfR) {
        assert(ver_del[e] == 1);
        ver_del[e] = 0;
        for(ui i = t_p_pstart[e]; i < t_p_pend[e]; i++) {
            ui u = t_p_edges[i];
            if(part[u]==1 || part[u]==2) {
                ++ pdR[u];
            }
        }
        for(ui i = t_n_pstart[e]; i < t_n_pend[e]; i++) {
            ui u = t_n_edges[i];
            if(part[u]==1 || part[u]==2) {
                ++ ndR[u];
            }
        }
    }
}

int vr_rebuildg()
{
    memset(ver_del, 0, sizeof(bool)*n);
    ui pt = tau - k;
    ui nt = tau - k + 1;
    ui dt = 2 * tau - k;
    int del_cnt = 0;
    queue<ui> Q;
    for(ui i = 0; i < n; i++) if(p_degree[i] < pt || n_degree[i] < nt || degree[i] < dt) {
        Q.push(i);
        ver_del[i] = 1;
        ++ del_cnt;
    }
    while (!Q.empty()) {
        ui u = Q.front();
        Q.pop();
        for(ui i = p_pstart[u]; i < p_pend[u]; i++) {
            ui v = p_edges[i];
            if(ver_del[v]) continue;
            assert(p_degree[v] > 0 && degree[v] > 0);
            -- p_degree[v];
            -- degree[v];
            if(p_degree[v] < pt || degree[v] < dt) {
                Q.push(v);
                ver_del[v] = 1;
                ++ del_cnt;
            }
        }
        for(ui i = n_pstart[u]; i < n_pend[u]; i++) {
            ui v = n_edges[i];
            if(ver_del[v]) continue;
            assert(n_degree[v] > 0 && degree[v] > 0);
            -- n_degree[v];
            -- degree[v];
            if(n_degree[v] < nt || degree[v] < dt) {
                Q.push(v);
                ver_del[v] = 1;
                ++ del_cnt;
            }
        }
    }
    
    ui * mapping = new ui[n];
    ui idx = 0;
    pm = 0;
    nm = 0;
    for(ui u = 0; u < n; u++) if(!ver_del[u]) true_id[u] = original_id[u];
    for(ui u = 0; u < n; u++) if(!ver_del[u]) {
        mapping[u] = idx;
        original_id[idx++] = true_id[u];
        for(ui i = p_pstart[u]; i < p_pend[u]; i ++) if(!ver_del[p_edges[i]]) ++ pm;
        for(ui i = n_pstart[u]; i < n_pend[u]; i ++) if(!ver_del[n_edges[i]]) ++ nm;
    }
    assert(pm%2==0 && nm%2==0);
    ui original_n = n;
    n = idx;
    pm /= 2;
    nm /= 2;
    t_p_pstart = new ui[n+1];
    t_p_edges = new ui[2*pm];
    t_n_pstart = new ui[n+1];
    t_n_edges = new ui[2*nm];
    t_p_pend = new ui[n+1];
    t_n_pend = new ui[n+1];
    ui new_i = 0;
    t_p_pstart[0] = 0;
    for(ui i = 0; i < original_n; i++) if(!ver_del[i]){
        ui pos = t_p_pstart[new_i];
        for(ui j = p_pstart[i]; j < p_pend[i]; j++) if(!ver_del[p_edges[j]]){
            t_p_edges[pos++] = mapping[p_edges[j]];
        }
        t_p_pstart[++new_i] = pos;
    }
    assert(new_i == n);
    new_i = 0;
    t_n_pstart[0] = 0;
    for(ui i = 0; i < original_n; i++) if(!ver_del[i]){
        ui pos = t_n_pstart[new_i];
        for(ui j = n_pstart[i]; j < n_pend[i]; j++) if(!ver_del[n_edges[j]]){
            t_n_edges[pos++] = mapping[n_edges[j]];
        }
        t_n_pstart[++new_i] = pos;
    }
    assert(new_i == n);
    delete [] p_pstart;
    delete [] p_edges;
    delete [] p_pend;
    delete [] n_pstart;
    delete [] n_edges;
    delete [] n_pend;
    p_pstart = t_p_pstart;
    p_pend = t_p_pend;
    p_edges = t_p_edges;
    n_pstart = t_n_pstart;
    n_pend = t_n_pend;
    n_edges = t_n_edges;
    for(ui u = 0; u < n; u++) {
        p_pend[u] = p_pstart[u+1];
        n_pend[u] = n_pstart[u+1];
        degree[u] = (p_pstart[u+1] - p_pstart[u]) + (n_pstart[u+1] - n_pstart[u]);
        if(degree[u] > max_deg) max_deg = degree[u];
    }
    delete [] mapping;
    return del_cnt;
}

int er_rebuildg_deep()
{
    int dele_count = 0;
    ui * vs = new ui[n];
    ui * C = new ui[max_deg+1];
    memset(C, 0, sizeof(ui)*(max_deg+1));
    for(ui i = 0; i < n; i++) ++ C[degree[i]];
    for(ui i = 1; i <= max_deg; i++) C[i] += C[i-1];
    for(ui i = 0; i < n; i++) vs[--C[degree[i]]] = i;
    ui * mark = new ui[n];
    memset(mark, 0, sizeof(ui)*n);
    ui * del = new ui[n];
    memset(del, 0, sizeof(ui)*n);
    sup_pp = new dense_hash_map<ui, int>[n];
    sup_nn = new dense_hash_map<ui, int>[n];
    sup_np = new dense_hash_map<ui, int>[n];
    e_sign = new dense_hash_map<ui, int>[n];
    e_del = new dense_hash_map<ui, bool>[n];
    e_inQ = new dense_hash_map<ui, bool>[n];
    for(ui i = 0; i < n; i++){
        sup_pp[i].resize(degree[i]);
        sup_nn[i].resize(degree[i]);
        sup_np[i].resize(degree[i]);
        e_sign[i].resize(degree[i]);
        e_del[i].resize(degree[i]);
        e_inQ[i].resize(degree[i]);
    }
    for(ui i = 0; i < n; i++){
        sup_pp[i].set_empty_key(INF);
        sup_nn[i].set_empty_key(INF);
        sup_np[i].set_empty_key(INF);
        e_sign[i].set_empty_key(INF);
        e_del[i].set_empty_key(INF);
        e_inQ[i].set_empty_key(INF);
    }
    for(ui i = 0; i < n; i++){
        for(ui j = p_pstart[i]; j < p_pstart[i+1]; j++){
            sup_pp[i][p_edges[j]] = 0; sup_pp[p_edges[j]][i] = 0;
            sup_nn[i][p_edges[j]] = 0; sup_nn[p_edges[j]][i] = 0;
            sup_np[i][p_edges[j]] = 0; sup_np[p_edges[j]][i] = 0;
            e_sign[i][p_edges[j]] = 1; e_sign[p_edges[j]][i] = 1;
            e_del[i][p_edges[j]] = 0; e_del[p_edges[j]][i] = 0;
            e_inQ[i][p_edges[j]] = 0; e_inQ[p_edges[j]][i] = 0;
        }
        for(ui j = n_pstart[i]; j < n_pstart[i+1]; j++){
            sup_pp[i][n_edges[j]] = 0; sup_pp[n_edges[j]][i] = 0;
            sup_nn[i][n_edges[j]] = 0; sup_nn[n_edges[j]][i] = 0;
            sup_np[i][n_edges[j]] = 0; sup_np[n_edges[j]][i] = 0;
            e_sign[i][n_edges[j]] = -1; e_sign[n_edges[j]][i] = -1;
            e_del[i][n_edges[j]] = 0; e_del[n_edges[j]][i] = 0;
            e_inQ[i][n_edges[j]] = 0; e_inQ[n_edges[j]][i] = 0;
        }
    }
    for(int i = n-1; i >= 0; i--){
        ui u = vs[i];
        for(ui j = p_pstart[u]; j < p_pstart[u+1]; j++) mark[p_edges[j]] = 1;
        for(ui j = n_pstart[u]; j < n_pstart[u+1]; j++) mark[n_edges[j]] = 2;
        for(ui j = p_pstart[u]; j < p_pstart[u+1]; j++){
            ui v = p_edges[j];
            if(!del[v]){
                for(ui k = p_pstart[v]; k < p_pstart[v+1]; k++){
                    ui w = p_edges[k];
                    if(!del[w] && mark[w] == 1 && w > v){
                        ++ sup_pp[u][v]; ++ sup_pp[v][u];
                        ++ sup_pp[u][w]; ++ sup_pp[w][u];
                        ++ sup_pp[v][w]; ++ sup_pp[w][v];
                    }
                }
                for(ui k = n_pstart[v]; k < n_pstart[v+1]; k++){
                    ui w = n_edges[k];
                    if(!del[w] && mark[w] == 2 && w > v){
                        ++ sup_nn[u][v]; ++ sup_nn[v][u];
                        ++ sup_np[u][w]; ++ sup_np[w][u];
                        ++ sup_np[v][w]; ++ sup_np[w][v];
                    }
                }
            }
        }
        for(ui j = n_pstart[u]; j < n_pstart[u+1]; j++){
            ui v = n_edges[j];
            if(!del[v]){
                for(ui k = p_pstart[v]; k < p_pstart[v+1]; k++){
                    ui w = p_edges[k];
                    if(!del[w] && mark[w] == 2 && w > v){
                        ++ sup_np[u][v]; ++ sup_np[v][u];
                        ++ sup_np[u][w]; ++ sup_np[w][u];
                        ++ sup_nn[v][w]; ++ sup_nn[w][v];
                    }
                }
                for(ui k = n_pstart[v]; k < n_pstart[v+1]; k++){
                    ui w = n_edges[k];
                    if(!del[w] && mark[w] == 1 && w > v){
                        ++ sup_np[u][v]; ++ sup_np[v][u];
                        ++ sup_np[v][w]; ++ sup_np[w][v];
                        ++ sup_nn[u][w]; ++ sup_nn[w][u];
                    }
                }
            }
        }
        for(ui j = p_pstart[u]; j < p_pstart[u+1]; j++) mark[p_edges[j]] = 0;
        for(ui j = n_pstart[u]; j < n_pstart[u+1]; j++) mark[n_edges[j]] = 0;
        del[u] = 1;
    }//u
    queue<pair<ui, ui>> Q;
    for(ui i = 0; i < n; i++){
        for(ui j = p_pstart[i]; j <p_pstart[i+1]; j++){
            ui v = p_edges[j];
            if(i < v){
                if(sup_pp[i][v] < tau-2*k || sup_nn[i][v] < tau-2*k+2 || sup_pp[i][v]+sup_nn[i][v] < 2*tau-2*k){
                    Q.push(make_pair(i, v));
                    e_inQ[i][v] = 1; e_inQ[v][i] = 1;
                }
            }
        }
        for(ui j = n_pstart[i]; j <n_pstart[i+1]; j++){
            ui v = n_edges[j];
            if(i < v){
                if(sup_np[i][v] < 2*tau-2*k){
                    Q.push(make_pair(i, v));
                    e_inQ[i][v] = 1; e_inQ[v][i] = 1;
                }
            }
        }
    }
    while (!Q.empty()) {
        pair<ui, ui> te = Q.front();
        ++ dele_count;
        Q.pop();
        ui u = te.first;
        ui v = te.second;
        if(e_sign[u][v] == 1){ // e(u,v) is positive
            reord_deg(u, v);
            for(ui i = p_pstart[u]; i < p_pstart[u+1]; i++){
                ui w = p_edges[i];
                if(e_sign[v].find(w) != e_sign[v].end() && e_sign[v][w] == 1 && !e_del[u][w] && !e_del[v][w]){
                    assert(sup_pp[u][w]>0 && sup_pp[w][u]>0);
                    -- sup_pp[u][w]; -- sup_pp[w][u];
                    if(!e_inQ[u][w] && (sup_pp[u][w] < tau-2*k || sup_nn[u][w] < tau-2*k+2 || sup_pp[u][w]+sup_nn[u][w] < 2*tau-2*k) ) {
                        Q.push(make_pair(u, w));
                        e_inQ[u][w] = 1; e_inQ[w][u] = 1;
                    }
                    assert(sup_pp[v][w]>0 && sup_pp[w][v]>0);
                    -- sup_pp[v][w]; -- sup_pp[w][v];
                    if(!e_inQ[v][w] && (sup_pp[v][w] < tau-2*k || sup_nn[v][w] < tau-2*k+2 || sup_pp[v][w]+sup_nn[v][w] < 2*tau-2*k) ) {
                        Q.push(make_pair(v, w));
                        e_inQ[v][w] = 1; e_inQ[w][v] = 1;
                    }
                }
            }
            for(ui i = n_pstart[u]; i < n_pstart[u+1]; i++){
                ui w = n_edges[i];
                if(e_sign[v].find(w) != e_sign[v].end() && e_sign[v][w] == -1 && !e_del[u][w] && !e_del[v][w]){
                    assert(sup_np[u][w]>0 && sup_np[w][u]>0);
                    -- sup_np[u][w]; -- sup_np[w][u];
                    if(!e_inQ[u][w] && (sup_np[u][w] < 2*tau-2*k) ) {
                        Q.push(make_pair(u, w));
                        e_inQ[u][w] = 1; e_inQ[w][u] = 1;
                    }
                    assert(sup_np[v][w]>0 && sup_np[w][v]>0);
                    -- sup_np[v][w]; -- sup_np[w][v];
                    if(!e_inQ[v][w] && (sup_np[v][w] < 2*tau-2*k) ) {
                        Q.push(make_pair(v, w));
                        e_inQ[v][w] = 1; e_inQ[w][v] = 1;
                    }
                }
            }
        }
        else{ // e(u,v) is negative
            reord_deg(u, v);
            for(ui i = p_pstart[u]; i < p_pstart[u+1]; i++){
                ui w = p_edges[i];
                if(e_sign[v].find(w) != e_sign[v].end() && e_sign[v][w] == -1 && !e_del[u][w] && !e_del[v][w]){
                    assert(sup_nn[u][w]>0 && sup_nn[w][u]>0);
                    --sup_nn[u][w]; --sup_nn[w][u];
                    if(!e_inQ[u][w] && (sup_pp[u][w] < tau-2*k || sup_nn[u][w] < tau-2*k+2 || sup_pp[u][w]+sup_nn[u][w] < 2*tau-2*k) ) {
                        Q.push(make_pair(u, w));
                        e_inQ[u][w] = 1; e_inQ[w][u] = 1;
                    }
                    assert(sup_np[v][w]>0 && sup_np[w][v]>0);
                    -- sup_np[v][w]; -- sup_np[w][v];
                    if(!e_inQ[v][w] && (sup_np[v][w] < 2*tau-2*k) ) {
                        Q.push(make_pair(v, w));
                        e_inQ[v][w] = 1; e_inQ[w][v] = 1;
                    }
                }
            }
            for(ui i = n_pstart[u]; i < n_pstart[u+1]; i++){
                ui w = n_edges[i];
                if(e_sign[v].find(w) != e_sign[v].end() && e_sign[v][w] == 1 && !e_del[u][w] && !e_del[v][w]){
                    assert(sup_nn[v][w]>0 && sup_nn[w][v]>0);
                    --sup_nn[v][w]; --sup_nn[w][v];
                    if(!e_inQ[v][w] && (sup_pp[v][w] < tau-2*k || sup_nn[v][w] < tau-2*k+2 || sup_pp[v][w]+sup_nn[v][w] < 2*tau-2*k) ) {
                        Q.push(make_pair(v, w));
                        e_inQ[v][w] = 1; e_inQ[w][v] = 1;
                    }
                    assert(sup_np[u][w]>0 && sup_np[w][u]>0);
                    -- sup_np[u][w]; -- sup_np[w][u];
                    if(!e_inQ[u][w] && (sup_np[u][w] < 2*tau-2*k) ) {
                        Q.push(make_pair(u, w));
                        e_inQ[u][w] = 1; e_inQ[w][u] = 1;
                    }
                }
            }
        }
        e_del[u][v] = 1;
        e_del[v][u] = 1;
    }
    
    for(ui i = 0; i < n; i++){
        p_pend[i] = p_pstart[i];
        n_pend[i] = n_pstart[i];
    }
    for(ui i = 0; i < n; i++){
        for(ui j = p_pstart[i]; j < p_pstart[i+1]; j++){
            ui u = p_edges[j]; //e(i,u)
            if(!e_del[i][u]){
                p_edges[p_pend[i]++] = u;
            }
        }
    }
    for(ui i = 0; i < n; i++){
        for(ui j = n_pstart[i]; j < n_pstart[i+1]; j++){
            ui u = n_edges[j]; //e(i,u)
            if(!e_del[i][u]){
                n_edges[n_pend[i]++] = u;
            }
        }
    }
    
    //update degree
    for(ui i = 0; i < n; i++){
        p_degree[i] = p_pend[i] - p_pstart[i];
        n_degree[i] = n_pend[i] - n_pstart[i];
        degree[i] = p_degree[i] + n_degree[i];
    }
    
    delete [] vs;
    delete [] C;
    delete [] mark;
    delete [] del;
    delete [] sup_pp;
    delete [] sup_nn;
    delete [] sup_np;
    delete [] e_sign;
    delete [] e_del;
    delete [] e_inQ;
    return dele_count;
}

int er_rebuildg()
{
    int dele_count = 0;
    int pp_thre = max(tau - 2*k, 0);
    int nn_thre = max(tau - 2*k + 2, 0);
    int np_thre = max(tau - 2*k + 1, 0);
    ui * vs = new ui[n];
    ui * C = new ui[max_deg+1];
    memset(C, 0, sizeof(ui)*(max_deg+1));
    for(ui i = 0; i < n; i++) ++ C[degree[i]];
    for(ui i = 1; i <= max_deg; i++) C[i] += C[i-1];
    for(ui i = 0; i < n; i++) vs[--C[degree[i]]] = i;
    ui * mark = new ui[n];
    memset(mark, 0, sizeof(ui)*n);
    ui * del = new ui[n];
    memset(del, 0, sizeof(ui)*n);
    sup_pp = new dense_hash_map<ui, int>[n];
    sup_nn = new dense_hash_map<ui, int>[n];
    sup_np = new dense_hash_map<ui, int>[n];
    e_sign = new dense_hash_map<ui, int>[n];
    e_del = new dense_hash_map<ui, bool>[n];
    for(ui i = 0; i < n; i++){
        sup_pp[i].resize(degree[i]);
        sup_nn[i].resize(degree[i]);
        sup_np[i].resize(degree[i]);
        e_sign[i].resize(degree[i]);
        e_del[i].resize(degree[i]);
    }
    for(ui i = 0; i < n; i++){
        sup_pp[i].set_empty_key(INF);
        sup_nn[i].set_empty_key(INF);
        sup_np[i].set_empty_key(INF);
        e_sign[i].set_empty_key(INF);
        e_del[i].set_empty_key(INF);
    }
    for(ui i = 0; i < n; i++){
        for(ui j = p_pstart[i]; j < p_pstart[i+1]; j++){
            sup_pp[i][p_edges[j]] = 0; sup_pp[p_edges[j]][i] = 0;
            sup_nn[i][p_edges[j]] = 0; sup_nn[p_edges[j]][i] = 0;
            sup_np[i][p_edges[j]] = 0; sup_np[p_edges[j]][i] = 0;
            e_sign[i][p_edges[j]] = 1; e_sign[p_edges[j]][i] = 1;
            e_del[i][p_edges[j]] = 0; e_del[p_edges[j]][i] = 0;
        }
        for(ui j = n_pstart[i]; j < n_pstart[i+1]; j++){
            sup_pp[i][n_edges[j]] = 0; sup_pp[n_edges[j]][i] = 0;
            sup_nn[i][n_edges[j]] = 0; sup_nn[n_edges[j]][i] = 0;
            sup_np[i][n_edges[j]] = 0; sup_np[n_edges[j]][i] = 0;
            e_sign[i][n_edges[j]] = -1; e_sign[n_edges[j]][i] = -1;
            e_del[i][n_edges[j]] = 0; e_del[n_edges[j]][i] = 0;
        }
    }
    for(int i = n-1; i >= 0; i--){
        ui u = vs[i];
        for(ui j = p_pstart[u]; j < p_pstart[u+1]; j++) mark[p_edges[j]] = 1;
        for(ui j = n_pstart[u]; j < n_pstart[u+1]; j++) mark[n_edges[j]] = 2;
        for(ui j = p_pstart[u]; j < p_pstart[u+1]; j++){
            ui v = p_edges[j];
            if(!del[v]){
                for(ui k = p_pstart[v]; k < p_pstart[v+1]; k++){
                    ui w = p_edges[k];
                    if(!del[w] && mark[w] == 1 && w > v){
                        ++ sup_pp[u][v]; ++ sup_pp[v][u];
                        ++ sup_pp[u][w]; ++ sup_pp[w][u];
                        ++ sup_pp[v][w]; ++ sup_pp[w][v];
                    }
                }
                for(ui k = n_pstart[v]; k < n_pstart[v+1]; k++){
                    ui w = n_edges[k];
                    if(!del[w] && mark[w] == 2 && w > v){
                        ++ sup_nn[u][v]; ++ sup_nn[v][u];
                        ++ sup_np[u][w]; ++ sup_np[w][u];
                        ++ sup_np[v][w]; ++ sup_np[w][v];
                    }
                }
            }
        }
        for(ui j = n_pstart[u]; j < n_pstart[u+1]; j++){
            ui v = n_edges[j];
            if(!del[v]){
                for(ui k = p_pstart[v]; k < p_pstart[v+1]; k++){
                    ui w = p_edges[k];
                    if(!del[w] && mark[w] == 2 && w > v){
                        ++ sup_np[u][v]; ++ sup_np[v][u];
                        ++ sup_np[u][w]; ++ sup_np[w][u];
                        ++ sup_nn[v][w]; ++ sup_nn[w][v];
                    }
                }
                for(ui k = n_pstart[v]; k < n_pstart[v+1]; k++){
                    ui w = n_edges[k];
                    if(!del[w] && mark[w] == 1 && w > v){
                        ++ sup_np[u][v]; ++ sup_np[v][u];
                        ++ sup_np[v][w]; ++ sup_np[w][v];
                        ++ sup_nn[u][w]; ++ sup_nn[w][u];
                    }
                }
            }
        }
        for(ui j = p_pstart[u]; j < p_pstart[u+1]; j++) mark[p_edges[j]] = 0;
        for(ui j = n_pstart[u]; j < n_pstart[u+1]; j++) mark[n_edges[j]] = 0;
        del[u] = 1;
    }//u
    queue<pair<ui, ui>> Q;
    for(ui i = 0; i < n; i++){
        for(ui j = p_pstart[i]; j <p_pstart[i+1]; j++){
            ui v = p_edges[j];
            if(i < v){
                if(sup_pp[i][v] < pp_thre || sup_nn[i][v] < nn_thre){
                    Q.push(make_pair(i, v));
                }
            }
        }
        for(ui j = n_pstart[i]; j <n_pstart[i+1]; j++){
            ui v = n_edges[j];
            if(i < v){
                if(sup_np[i][v] < np_thre){
                    Q.push(make_pair(i, v));
                }
            }
        }
    }
    while (!Q.empty()) {
        pair<ui, ui> te = Q.front();
        ++ dele_count;
        Q.pop();
        ui u = te.first;
        ui v = te.second;
        if(e_sign[u][v] == 1){ // e(u,v) is positive
            reord_deg(u, v);
            for(ui i = p_pstart[u]; i < p_pstart[u+1]; i++){
                ui w = p_edges[i];
                if(e_sign[v].find(w) != e_sign[v].end() && e_sign[v][w] == 1 && !e_del[u][w] && !e_del[v][w]){
                    if((sup_pp[u][w] --) == pp_thre && sup_nn[u][w] >= nn_thre)
                    {
                        Q.push(make_pair(u, w));
                    }
                    -- sup_pp[w][u];
                    if((sup_pp[v][w] --) == pp_thre && sup_nn[v][w] >= nn_thre)
                    {
                        Q.push(make_pair(v, w));
                    }
                    -- sup_pp[w][v];
                }
            }
            for(ui i = n_pstart[u]; i < n_pstart[u+1]; i++){
                ui w = n_edges[i];
                if(e_sign[v].find(w) != e_sign[v].end() && e_sign[v][w] == -1 && !e_del[u][w] && !e_del[v][w]){
                    if((sup_np[u][w] --) == np_thre)
                    {
                        Q.push(make_pair(u, w));
                    }
                    -- sup_np[w][u];
                    if((sup_np[v][w] --) == np_thre)
                    {
                        Q.push(make_pair(v, w));
                    }
                    -- sup_np[w][v];
                }
            }
        }
        else{ // e(u,v) is negative
            reord_deg(u, v);
            for(ui i = p_pstart[u]; i < p_pstart[u+1]; i++){
                ui w = p_edges[i];
                if(e_sign[v].find(w) != e_sign[v].end() && e_sign[v][w] == -1 && !e_del[u][w] && !e_del[v][w]){
                    if((sup_nn[u][w] --) == nn_thre && sup_pp[u][w] >= pp_thre)
                    {
                        Q.push(make_pair(u, w));
                    }
                    -- sup_nn[w][u];
                    if((sup_np[v][w] --) == np_thre)
                    {
                        Q.push(make_pair(v, w));
                    }
                    -- sup_np[w][v];
                }
            }
            for(ui i = n_pstart[u]; i < n_pstart[u+1]; i++){
                ui w = n_edges[i];
                if(e_sign[v].find(w) != e_sign[v].end() && e_sign[v][w] == 1 && !e_del[u][w] && !e_del[v][w]){
                    if((sup_np[u][w] --) == np_thre)
                    {
                        Q.push(make_pair(u, w));
                    }
                    -- sup_np[w][u];
                    if((sup_nn[v][w] --) == nn_thre && sup_pp[v][w] >= pp_thre)
                    {
                        Q.push(make_pair(v, w));
                    }
                    -- sup_nn[w][v];
                }
            }
        }
        e_del[u][v] = 1;
        e_del[v][u] = 1;
    }
    
    for(ui i = 0; i < n; i++){
        p_pend[i] = p_pstart[i];
        n_pend[i] = n_pstart[i];
    }
    for(ui i = 0; i < n; i++){
        for(ui j = p_pstart[i]; j < p_pstart[i+1]; j++){
            ui u = p_edges[j]; //e(i,u)
            if(!e_del[i][u]){
                p_edges[p_pend[i]++] = u;
            }
        }
    }
    for(ui i = 0; i < n; i++){
        for(ui j = n_pstart[i]; j < n_pstart[i+1]; j++){
            ui u = n_edges[j]; //e(i,u)
            if(!e_del[i][u]){
                n_edges[n_pend[i]++] = u;
            }
        }
    }
    
    //update degree
    for(ui i = 0; i < n; i++){
        p_degree[i] = p_pend[i] - p_pstart[i];
        n_degree[i] = n_pend[i] - n_pstart[i];
        degree[i] = p_degree[i] + n_degree[i];
    }
    
    delete [] vs;
    delete [] C;
    delete [] mark;
    delete [] del;
    delete [] sup_pp;
    delete [] sup_nn;
    delete [] sup_np;
    delete [] e_sign;
    delete [] e_del;
    return dele_count;
}

//after running this function, must use p_pend and n_pend...
void CTprune()
{
    int round = 1;
    ver_del = new bool[n];
    memset(ver_del, 0, sizeof(bool)*n);
#ifdef _CTprune_
    while (round <= 10) {
        int dvnum = vr_rebuildg();
        int denum = er_rebuildg_deep();
        if(dvnum==0 && denum==0) break;
        ++ round;
    }
#endif
    pm = 0;
    nm = 0;
    for(ui u = 0; u < n; u++) {
        for(ui i = p_pstart[u]; i < p_pend[u]; i ++) ++ pm;
        for(ui i = n_pstart[u]; i < n_pend[u]; i ++) ++ nm;
    }
    assert(pm%2==0 && nm%2==0);
//    cout<<"n:"<<n<<", m:"<<pm+nm<<endl;
    pm /= 2;
    nm /= 2;
    t_p_pstart = new ui[n+1];
    t_p_edges = new ui[2*pm];
    t_n_pstart = new ui[n+1];
    t_n_edges = new ui[2*nm];
    t_p_pend = new ui[n+1]; //for later usage
    t_n_pend = new ui[n+1]; //for later usage
    memcpy(t_p_pstart, p_pstart, sizeof(ui)*(n+1));
    memcpy(t_p_pend, p_pend, sizeof(ui)*(n+1));
    memcpy(t_p_edges, p_edges, sizeof(ui)*(2*pm));
    memcpy(t_n_pstart, n_pstart, sizeof(ui)*(n+1));
    memcpy(t_n_pend, n_pend, sizeof(ui)*(n+1));
    memcpy(t_n_edges, n_edges, sizeof(ui)*(2*nm));
}

void reset_inCPX(vector<ui>&CL, vector<ui>&CR, vector<ui>&P, vector<ui>&XL, vector<ui>&XR)
{
    for(auto &e : CL) part[e] = 0;
    for(auto &e : CR) part[e] = 0;
    for(auto &e : P) {
        part[e] = 0;
        addtoKL[e] = 0;
        addtoKR[e] = 0;
    }
    for(auto &e : XL) inXL[e] = 0;
    for(auto &e : XR) inXR[e] = 0;
#ifdef _CostlyAssert_
    for(ui i = 0; i < n; i++) assert( part[i]==0 && !inXL[i] && !inXR[i] && !addtoKL[i] && !addtoKR[i]);
#endif
}

void reset_pnd_CP(vector<ui>&CL, vector<ui>&CR, vector<ui>&P)
{
    for(auto &e : CL) { pdL[e] = 0; pdR[e] = 0; ndL[e] = 0; ndR[e] = 0; }
    for(auto &e : CR) { pdL[e] = 0; pdR[e] = 0; ndL[e] = 0; ndR[e] = 0; }
    for(auto &e : P) { pdL[e] = 0; pdR[e] = 0; ndL[e] = 0; ndR[e] = 0; }
#ifdef _CostlyAssert_
    for(ui i = 0; i < n; i++) assert(pdL[i] == 0 && pdR[i] == 0 && ndL[i] == 0 && pdR[i] == 0);
#endif
}

//flag=1:   CL -> KR      -1: CL <- KR
//flag=2:   CR -> KL      -2: CR <- KL
//flag=3:   P -> KL       -3: P <- KL
//flag=4:   P -> KR       -4: P <- KR
void update_pndeg(ui v, int flag)
{
    if(flag == 1) {
        for(ui i = t_p_pstart[v]; i < t_p_pend[v]; i++) {
            ui u = t_p_edges[i];
            if(part[u] != 0){
                ++pdR[u];
                assert(pdL[u]>0);
                --pdL[u];
            }
        }
        for(ui i = t_n_pstart[v]; i < t_n_pend[v]; i++) {
            ui u = t_n_edges[i];
            if(part[u] != 0){
                ++ndR[u];
                assert(ndL[u]>0);
                --ndL[u];
            }
        }
    }
    else if(flag == -1) {
        for(ui i = t_p_pstart[v]; i < t_p_pend[v]; i++) {
            ui u = t_p_edges[i];
            if(part[u] != 0){
                assert(pdR[u]>0);
                --pdR[u];
                ++pdL[u];
            }
        }
        for(ui i = t_n_pstart[v]; i < t_n_pend[v]; i++) {
            ui u = t_n_edges[i];
            if(part[u] != 0){
                assert(ndR[u]>0);
                --ndR[u];
                ++ndL[u];
            }
        }
    }
    else if (flag == 2) {
        for(ui i = t_p_pstart[v]; i < t_p_pend[v]; i++) {
            ui u = t_p_edges[i];
            if(part[u] != 0){
                ++pdL[u];
                assert(pdR[u]>0);
                --pdR[u];
            }
        }
        for(ui i = t_n_pstart[v]; i < t_n_pend[v]; i++) {
            ui u = t_n_edges[i];
            if(part[u] != 0){
                ++ndL[u];
                assert(ndR[u]>0);
                --ndR[u];
            }
        }
    }
    else if (flag == -2) {
        for(ui i = t_p_pstart[v]; i < t_p_pend[v]; i++) {
            ui u = t_p_edges[i];
            if(part[u] != 0){
                assert(pdL[u]>0);
                --pdL[u];
                ++pdR[u];
            }
        }
        for(ui i = t_n_pstart[v]; i < t_n_pend[v]; i++) {
            ui u = t_n_edges[i];
            if(part[u] != 0){
                assert(ndL[u]>0);
                --ndL[u];
                ++ndR[u];
            }
        }
    }
    else if (flag == 3) {
        for(ui i = t_p_pstart[v]; i < t_p_pend[v]; i++) {
            ui u = t_p_edges[i];
            if(part[u] != 0){
                ++pdL[u];
            }
        }
        for(ui i = t_n_pstart[v]; i < t_n_pend[v]; i++) {
            ui u = t_n_edges[i];
            if(part[u] != 0){
                ++ndL[u];
            }
        }
    }
    else if (flag == -3) {
        for(ui i = t_p_pstart[v]; i < t_p_pend[v]; i++) {
            ui u = t_p_edges[i];
            if(part[u] != 0){
                assert(pdL[u] > 0);
                --pdL[u];
            }
        }
        for(ui i = t_n_pstart[v]; i < t_n_pend[v]; i++) {
            ui u = t_n_edges[i];
            if(part[u] != 0){
                assert(ndL[u] > 0);
                --ndL[u];
            }
        }
    }
    else if (flag == 4) {
        for(ui i = t_p_pstart[v]; i < t_p_pend[v]; i++) {
            ui u = t_p_edges[i];
            if(part[u] != 0){
                ++pdR[u];
            }
        }
        for(ui i = t_n_pstart[v]; i < t_n_pend[v]; i++) {
            ui u = t_n_edges[i];
            if(part[u] != 0){
                ++ndR[u];
            }
        }
    }
    else if (flag == -4) {
        for(ui i = t_p_pstart[v]; i < t_p_pend[v]; i++) {
            ui u = t_p_edges[i];
            if(part[u] != 0){
                assert(pdR[u]>0);
                --pdR[u];
            }
        }
        for(ui i = t_n_pstart[v]; i < t_n_pend[v]; i++) {
            ui u = t_n_edges[i];
            if(part[u] != 0){
                assert(ndR[u]>0);
                --ndR[u];
            }
        }
    }
}

bool prune_by_pndeg(ui u, int ex, vector<ui>&rfCL, vector<ui>&rfCR, vector<ui>&rfP)
{
    assert(ex >= 0);
    if(ex==0) return false;
    queue<ui> Q;
    assert(!ver_del[u]);
    if(pdL[u]+k+ex < tau || ndR[u]+k-1+ex < tau || pdL[u]+ndR[u]+k+ex < 2*tau) {
        Q.push(u);
        ver_del[u] = 1;
    }
    for(auto &e : A) if(!ver_del[e]) {
        if(part[e]==1) {
            if( (pdL[e]+k+ex < tau || ndR[e]+k-1+ex < tau || pdL[e]+ndR[e]+k+ex < 2*tau) && (pdR[e]+k-1+ex < tau || ndL[e]+k-1+ex < tau || pdR[e]+ndL[e]+k+ex < 2*tau) ) {
                Q.push(e);
                ver_del[e] = 1;
            }
        }
        else if (part[e]==2) {
            if( (pdR[e]+k+ex < tau || ndL[e]+k-1+ex < tau || pdR[e]+ndL[e]+k+ex < 2*tau) && (pdL[e]+k+ex < tau || ndR[e]+k-2+ex < tau || pdL[e]+ndR[e]+k+ex < tau) ) {
                Q.push(e);
                ver_del[e] = 1;
            }
        }
        else {
            assert(part[e]==3);
            if( (pdL[e]+k+ex < tau || ndR[e]+k-2+ex < tau || pdL[e]+ndR[e]+k+ex < 2*tau) && (pdR[e]+k-1+ex < tau || ndL[e]+k-1+ex < tau || pdR[e]+ndL[e]+k+ex < 2*tau) ) {
                Q.push(e);
                ver_del[e] = 1;
            }
        }
    }
    while (!Q.empty()) {
        ui v = Q.front();
        Q.pop();
        if(part[v]==1) {
            rfCL.push_back(v);
            for(ui i = t_p_pstart[v]; i < t_p_pend[v]; i++) if(part[t_p_edges[i]]!=0) {
                ui w = t_p_edges[i];
                assert(pdL[w] > 0);
                -- pdL[w];
                if(ver_del[w]) continue;
                if(part[w]==1 && ( (pdL[w]+k+ex < tau || pdL[w]+ndR[w]+k+ex < 2*tau) && (pdR[w]+k-1+ex < tau || ndL[w]+k-1+ex < tau || pdR[w]+ndL[w]+k+ex < 2*tau) ) ) {
                    Q.push(w);
                    ver_del[w] = 1;
                }
                //
                else if(part[w]==2 && ( (pdR[w]+k+ex < tau || ndL[w]+k-1+ex < tau || pdR[w]+ndL[w]+k+ex < 2*tau) && (pdL[w]+k+ex < tau || pdL[w]+ndR[w]+k+ex < tau) ) ) {
                    Q.push(w);
                    ver_del[w] = 1;
                }
                //
                else if (part[w]==3 && ( (pdL[w]+k+ex < tau || pdL[w]+ndR[w]+k+ex < 2*tau) && (pdR[w]+k-1+ex < tau || ndL[w]+k-1+ex < tau || pdR[w]+ndL[w]+k+ex < 2*tau) ) ) {
                    Q.push(w);
                    ver_del[w] = 1;
                }
            }
            for(ui i = t_n_pstart[v]; i < t_n_pend[v]; i++) if(part[t_n_edges[i]]!=0) {
                ui w = t_n_edges[i];
                assert(ndL[w] > 0);
                -- ndL[w];
                if(ver_del[w]) continue;
                if(part[w]==2 && ( ( ndL[w]+k-1+ex < tau || pdR[w]+ndL[w]+k+ex < 2*tau) && (pdL[w]+k+ex < tau || ndR[w]+k-2+ex < tau || pdL[w]+ndR[w]+k+ex < tau) ) ) {
                    Q.push(w);
                    ver_del[w] = 1;
                }
                //
                else if(part[w]==1 && ( (pdL[w]+k+ex < tau || ndR[w]+k-1+ex < tau || pdL[w]+ndR[w]+k+ex < 2*tau) && ( ndL[w]+k-1+ex < tau || pdR[w]+ndL[w]+k+ex < 2*tau) ) ) {
                    Q.push(w);
                    ver_del[w] = 1;
                }
                //
                else if (part[w]==3 && ( (pdL[w]+k+ex < tau || ndR[w]+k-2+ex < tau || pdL[w]+ndR[w]+k+ex < 2*tau) && ( ndL[w]+k-1+ex < tau || pdR[w]+ndL[w]+k+ex < 2*tau) ) ) {
                    Q.push(w);
                    ver_del[w] = 1;
                }
            }
        }
        else if (part[v]==2) {
            rfCR.push_back(v);
            for(ui i = t_p_pstart[v]; i < t_p_pend[v]; i++) if(part[t_p_edges[i]]!=0) {
                ui w = t_p_edges[i];
                assert(pdR[w] > 0);
                -- pdR[w];
                if(ver_del[w]) continue;
                if(part[w]==2 && ( (pdR[w]+k+ex < tau || pdR[w]+ndL[w]+k+ex < 2*tau) && (pdL[w]+k+ex < tau || ndR[w]+k-2+ex < tau || pdL[w]+ndR[w]+k+ex < tau) ) ) {
                    Q.push(w);
                    ver_del[w] = 1;
                }
                //
                else if(part[w]==1 && ( (pdL[w]+k+ex < tau || ndR[w]+k-1+ex < tau || pdL[w]+ndR[w]+k+ex < 2*tau) && (pdR[w]+k-1+ex < tau || pdR[w]+ndL[w]+k+ex < 2*tau) ) ) {
                    Q.push(w);
                    ver_del[w] = 1;
                }
                //
                else if (part[w]==3 && ( (pdL[w]+k+ex < tau || ndR[w]+k-2+ex < tau || pdL[w]+ndR[w]+k+ex < 2*tau) && (pdR[w]+k-1+ex < tau || pdR[w]+ndL[w]+k+ex < 2*tau) ) ) {
                    Q.push(w);
                    ver_del[w] = 1;
                }
            }
            for(ui i = t_n_pstart[v]; i < t_n_pend[v]; i++) if(part[t_n_edges[i]]!=0) {
                ui w = t_n_edges[i];
                assert(ndR[w] > 0);
                -- ndR[w];
                if(ver_del[w]) continue;
                if(part[w]==1 && ( (ndR[w]+k-1+ex < tau || pdL[w]+ndR[w]+k+ex < 2*tau) && (pdR[w]+k-1+ex < tau || ndL[w]+k-1+ex < tau || pdR[w]+ndL[w]+k+ex < 2*tau) ) ) {
                    Q.push(w);
                    ver_del[w] = 1;
                }
                //
                else if(part[w]==2 && ( (pdR[w]+k+ex < tau || ndL[w]+k-1+ex < tau || pdR[w]+ndL[w]+k+ex < 2*tau) && ( ndR[w]+k-2+ex < tau || pdL[w]+ndR[w]+k+ex < tau) ) ) {
                    Q.push(w);
                    ver_del[w] = 1;
                }
                //
                else if (part[w]==3 && ( (ndR[w]+k-2+ex < tau || pdL[w]+ndR[w]+k+ex < 2*tau) && (pdR[w]+k-1+ex < tau || ndL[w]+k-1+ex < tau || pdR[w]+ndL[w]+k+ex < 2*tau) ) ) {
                    Q.push(w);
                    ver_del[w] = 1;
                }
            }
        }
        else {
            assert(part[v]==3);
            rfP.push_back(v);
        }
    }
    for(auto &e : rfCL) if(inK[e]) return true;
    for(auto &e : rfCR) if(inK[e]) return true;
    for(auto &e : rfP) if(inK[e]) return true;
    return false;
}

void prune_by_pndeg_recover(vector<ui>&rfCL, vector<ui>&rfCR, vector<ui>&rfP)
{
    for(auto &e : rfCL) {
        assert(ver_del[e]);
        ver_del[e] = 0;
        for(ui i = t_p_pstart[e]; i < t_p_pend[e]; i++) if(part[t_p_edges[i]]!=0) {
            ++pdL[t_p_edges[i]];
        }
        for(ui i = t_n_pstart[e]; i < t_n_pend[e]; i++) if(part[t_n_edges[i]]!=0) {
            ++ndL[t_n_edges[i]];
        }
    }
    for(auto &e : rfCR) {
        assert(ver_del[e]);
        ver_del[e] = 0;
        for(ui i = t_p_pstart[e]; i < t_p_pend[e]; i++) if(part[t_p_edges[i]]!=0) {
            ++pdR[t_p_edges[i]];
        }
        for(ui i = t_n_pstart[e]; i < t_n_pend[e]; i++) if(part[t_n_edges[i]]!=0) {
            ++ndR[t_n_edges[i]];
        }
    }
    for(auto &e : rfP) {
        assert(ver_del[e]);
        ver_del[e] = 0;
    }
}

void seedg_rec(ui u, ui sidx, ui eidx, int kp)
{
    if(kp == 0) {
        if(sKL.size()+sKR.size()==k){
            vector<ui> newCL, newCR;
            for(auto &e : sCL) if(!inK[e] && !ver_del[e]) {
                newCL.push_back(e);
            }
            for(auto &e : sCR) if(!inK[e] && !ver_del[e]) {
                newCR.push_back(e);
            }
            invoke_enum++;
    #ifdef _mindegpivot_
            enum_procedure_mindeg_pivot(sKL, sKR, newCL, newCR, sXL, sXR, 1);
    #else
            enum_procedure_adv_opt(sKL, sKR, newCL, newCR, sXL, sXR, 1, true);
    #endif
            return;
        }
        else{
            int addtoXLnum = 0, addtoXRnum = 0;
            vector<ui> newCL, newCR;
            for(auto &e : sCL) if(!inK[e] && !ver_del[e]) {
                newCL.push_back(e);
                sXR.push_back(e);
                assert(inXR[e]==0);
                inXR[e]=1;
                ++ addtoXRnum;
            }
            for(auto &e : sCR) if(!inK[e] && !ver_del[e]) {
                newCR.push_back(e);
                sXL.push_back(e);
                assert(inXL[e]==0);
                inXL[e]=1;
                ++ addtoXLnum;
            }
            for(auto &e : sP) if(!inK[e] && !ver_del[e]) {
                if(addtoKL[e]) {
                    sXL.push_back(e);
                    assert(inXL[e]==0);
                    inXL[e]=1;
                    ++ addtoXLnum;
                }
                if(addtoKR[e]) {
                    sXR.push_back(e);
                    assert(inXR[e]==0);
                    inXR[e]=1;
                    ++ addtoXRnum;
                }
            }
            invoke_enum++;
    #ifdef _mindegpivot_
            enum_procedure_mindeg_pivot(sKL, sKR, newCL, newCR, sXL, sXR, 1);
    #else
            enum_procedure_adv_opt(sKL, sKR, newCL, newCR, sXL, sXR, 1, true);
    #endif
            while (addtoXLnum>0) {
                assert(inXL[sXL.back()]==1);
                inXL[sXL.back()]=0;
                sXL.pop_back();
                --addtoXLnum;
            }
            while (addtoXRnum>0) {
                assert(inXR[sXR.back()]==1);
                inXR[sXR.back()]=0;
                sXR.pop_back();
                --addtoXRnum;
            }
            return;
        }
        
    }
    for(ui i = sidx; i < eidx; i++) if(!ver_del[A[i]]) {
        ui v = A[i];
        if(part[v] == 1) {//v: CL -> KR
            sKR.push_back(v);
            assert(inK[v]==0);
            inK[v]=1;
            part[v]=2;
            bool core_pruned = false;
#ifdef _PP_
            update_pndeg(v,1);
            vector<ui> rfCL, rfCR, rfP;
            core_pruned = prune_by_pndeg(u,kp-1,rfCL,rfCR,rfP);
#endif
            if(!core_pruned) seedg_rec(u, i+1, eidx, kp-1);
#ifdef _PP_
            prune_by_pndeg_recover(rfCL,rfCR,rfP);
            update_pndeg(v,-1);
#endif
            assert(sKR.back()==v);
            sKR.pop_back();
            inK[v]=0;
            part[v] = 1;
        }
        else if (part[v] == 2) {//v: CR -> KL
            sKL.push_back(v);
            assert(inK[v]==0);
            inK[v]=1;
            part[v] = 1;
            bool core_pruned = false;
#ifdef _PP_
            update_pndeg(v,2);
            vector<ui> rfCL, rfCR, rfP;
            core_pruned = prune_by_pndeg(u,kp-1,rfCL,rfCR,rfP);
#endif
            if(!core_pruned) seedg_rec(u, i+1, eidx, kp-1);
#ifdef _PP_
            prune_by_pndeg_recover(rfCL,rfCR,rfP);
            update_pndeg(v,-2);
#endif
            assert(sKL.back()==v);
            sKL.pop_back();
            inK[v]=0;
            part[v] = 2;
        }
        else {
            assert(part[v] == 3);
            if(addtoKL[v]){
                sKL.push_back(v);//v: P -> KL
                assert(inK[v]==0);
                inK[v]=1;
                part[v]=1;
                bool core_pruned = false;
#ifdef _PP_
                update_pndeg(v,3);
                vector<ui> rfCL, rfCR, rfP;
                core_pruned = prune_by_pndeg(u,kp-1,rfCL,rfCR,rfP);
#endif
                if (!core_pruned) seedg_rec(u, i+1, eidx, kp-1);
#ifdef _PP_
                prune_by_pndeg_recover(rfCL,rfCR,rfP);
                update_pndeg(v,-3);
#endif
                assert(sKL.back()==v);
                sKL.pop_back();
                assert(inK[v]==1);
                inK[v]=0;
                part[v]=3;
            }
            if(addtoKR[v]){
                sKR.push_back(v);//v: P -> KR
                assert(inK[v]==0);
                inK[v] = 1;
                part[v] = 2;
                bool core_pruned = false;
#ifdef _PP_
                update_pndeg(v,4);
                vector<ui> rfCL, rfCR, rfP;
                core_pruned = prune_by_pndeg(u,kp-1,rfCL,rfCR,rfP);
#endif
                if (!core_pruned) seedg_rec(u, i+1, eidx, kp-1);
#ifdef _PP_
                prune_by_pndeg_recover(rfCL,rfCR,rfP);
                update_pndeg(v,-4);
#endif
                assert(sKR.back()==v);
                sKR.pop_back();
                assert(inK[v]==1);
                inK[v]=0;
                part[v] = 3;
            }
        }
    }
}

void maximal_balanced_kplex_enum_adv_opt()
{
    CTprune();
    process_order = new ui[n];
    for(ui i = 0; i < n; i++) process_order[i] = i;
    if(n==0) return;
    comp_process_order();
    ver_rank = new ui[n];
    for(ui i = 0; i < n; i++) ver_rank[process_order[i]] = i;
    part = new short[n]; memset(part, 0, sizeof(short)*n);
    inXL = new bool[n]; memset(inXL, 0, sizeof(bool)*n);
    inXR = new bool[n]; memset(inXR, 0, sizeof(bool)*n);
    inK = new bool[n]; memset(inK, 0, sizeof(bool)*n);
    pdL = new int[n]; memset(pdL, 0, sizeof(int)*n);
    pdR = new int[n]; memset(pdR, 0, sizeof(int)*n);
    ndL = new int[n]; memset(ndL, 0, sizeof(int)*n);
    ndR = new int[n]; memset(ndR, 0, sizeof(int)*n);
    partial_kplex_cap = new short[n];
    addtoKL = new bool[n]; memset(addtoKL, 0, sizeof(bool)*n);
    addtoKR = new bool[n]; memset(addtoKR, 0, sizeof(bool)*n);
    memset(ver_del, 0, sizeof(bool)*n); //for other use
    inQv = new short[n];
    M_len = n;
    M = new short*[M_len]; for(int i = 0; i < M_len; i++) M[i] = new short[M_len];
    nneiKL.resize(n);
    nneiKR.resize(n);

#ifdef _ERinEnum_
    Pnei.resize(n);
    Nnei.resize(n);
    tri_pp = new int*[n]; for(int i = 0; i < n; i++) tri_pp[i] = new int[n];
    tri_nn = new int*[n]; for(int i = 0; i < n; i++) tri_nn[i] = new int[n];
    tri_pn = new int*[n]; for(int i = 0; i < n; i++) tri_pn[i] = new int[n];
    tri_np = new int*[n]; for(int i = 0; i < n; i++) tri_np[i] = new int[n];
    inQe = new short*[n]; //remark edge is inQe or deleted
    for(int i = 0; i < n; i++) {
        inQe[i] = new short[n];
        memset(inQe[i], 0, sizeof(short)*n);
    }
    v_sta = new bool[n]; memset(v_sta, 0, sizeof(bool)*n); //for triangle counting
#endif
    for(ui i = 0; i < n; i++) {
        ui u = process_order[i];
        sKL.clear(); sKR.clear();
        sCL.clear(); sCR.clear(); sP.clear();
        sXL.clear(); sXR.clear();
        obtain_CL_CR_P_XL_XR(u, sCL, sCR, sP, sXL, sXR);
        assert(sCL.size()+sCR.size()<=max_core);
        if(sCL.size()+k < tau || sCR.size()+k-1 < tau || sCL.size()+sCR.size()+sP.size()+1 < 2*tau) {
            reset_inCPX(sCL, sCR, sP, sXL, sXR);
            continue;
        }
        construct_initial_subgraph(u, sCL, sCR, sP);
        vertex_reduction_on_initial_subgraph(u, sCL, sCR, sP);
        if(ver_del[u] || sCL.size()+k < tau || sCR.size()+k-1 < tau || sCL.size()+sCR.size()+sP.size()+1 < 2*tau) {
            ver_del[u] = 0;
            reset_inCPX(sCL, sCR, sP, sXL, sXR);
            continue;
        }
        obtain_pnd_LRside(u, sCL, sCR, sP);
#ifdef _ParVR_
        if(refine_LRP_by_pndLR(u, sCL, sCR, sP)) {
            reset_inCPX(sCL, sCR, sP, sXL, sXR);
            reset_pnd_CP(sCL, sCR, sP);
            continue;
        }
#endif
        for(auto &e : sCL) core[e] = pdL[e]+ndR[e]+pdR[e]+ndL[e]; //abuse core[]
        for(auto &e : sCR) core[e] = pdR[e]+ndL[e]+pdL[e]+ndR[e]; //abuse core[]
        sort(sCL.begin(), sCL.end(), mycomp);
        sort(sCR.begin(), sCR.end(), mycomp);
        build_M(u, sCL, sCR, sP, sXL, sXR);
        A.clear();
        A.insert(A.end(), sCL.begin(), sCL.end());
        A.insert(A.end(), sCR.begin(), sCR.end());
        A.insert(A.end(), sP.begin(), sP.end());
        assert(!inK[u]);
        inK[u] = 1;
        sKL.push_back(u);
        assert(part[u]==0);
        part[u] = 1; //abuse part[], assuming u is in Lside, ease for programming
        for(int kp = 0; kp < k ; kp++) seedg_rec(u, 0, (ui)A.size(), kp);
        assert(sKL.back()==u);
        part[u] = 0;
        sKL.pop_back();
        inK[u]=0;
        pdL[u] = 0; pdR[u] = 0; ndL[u] = 0; ndR[u] = 0;
        reset_inCPX(sCL, sCR, sP, sXL, sXR);
        reset_pnd_CP(sCL, sCR, sP);
#ifdef _CostlyAssert_
        for(ui i = 0; i < n; i++) assert( part[i]==0 && !inXL[i] && !inXR[i] && !addtoKL[i] && !addtoKR[i]);
        for(ui i = 0; i < n; i++) assert( inK[i]==0 && ver_del[i]==0);
        for(ui i = 0; i < n; i++) assert( pdL[i] == 0 && pdR[i] == 0 && ndL[i] == 0 && pdR[i] == 0);
#endif
    }
}

void display_results()
{
    for(auto &e : Res) {
        sort(e.first.begin(), e.first.end());
        sort(e.second.begin(), e.second.end());
        vector<ui> tmp_vec;
        if(e.first[0] > e.second[0]) {
            tmp_vec = e.first;
            e.first = e.second;
            e.second = tmp_vec;
        }
    }
    sort(Res.begin(),Res.end());
    for(auto e : Res) {
        cout<<"("; for(auto x : e.first) cout<<x<<",";
        cout<<"|"; for(auto x : e.second) cout<<x<<","; cout<<")"<<endl;
    }
}

void ShowCaseStudyResult(string input_graph)
{
    //read the original graph
    string buffer;
    map<ui, int> * s_G;
    ifstream input_file(input_graph, ios::in);
    if (!input_file.is_open()){cout << "cannot open file : "<<input_graph<<endl;exit(1);}
    else{
        input_file >> n >> m;
        s_G = new map<ui, int>[n];
        ui tu, tv;
        int flag;
        while (input_file >> tu >> tv >> flag){
            if(tu == tv) continue;
            assert(tu >= 0 && tu < n); assert(tv >= 0 && tv < n);
            assert(flag == 1 || flag == -1);
            s_G[tu].insert(make_pair(tv, flag));
            s_G[tv].insert(make_pair(tu, flag));
        }
        m = 0; pm = 0; nm = 0;
        for(ui i = 0; i < n; i++){
            const map<ui, int> & nei = s_G[i];
            for(auto e : nei){
                if(e.second == 1) ++ pm;
                else ++ nm;
            }
            m += nei.size();
        }
        assert(m%2 == 0); assert(pm%2 == 0); assert(nm%2 == 0);
        m /= 2; pm /= 2; nm /= 2;
        input_file.close();
    }
    
    vector<pair<int, pair<vector<ui>, vector<ui>>>> Rlt1;
    for(auto e : Res) {
        vector<ui> c1;
        for(auto x : e.first) c1.push_back(x);
        sort(c1.begin(),c1.end());
        
        vector<ui> c2;
        for(auto x : e.second) c2.push_back(x);
        sort(c2.begin(),c2.end());
        
        if(c1[0] < c2[0])
            Rlt1.push_back(make_pair(c1.size()+c2.size(), make_pair(c1, c2)));
        else
            Rlt1.push_back(make_pair(c1.size()+c2.size(), make_pair(c2, c1)));
    }
    sort(Rlt1.begin(), Rlt1.end());
    
    cout<<"Mapped Results: "<<endl;
    long long cnt = 1;
    for(auto &e : Rlt1) {
        cout<<"C"<<cnt++<<" : ("<<e.second.first.size()<<","<<e.second.second.size()<<"):";
        cout<<"\t( ";
        for(auto &x : e.second.first) cout<<id2str[x]<<", ";
        cout<<" | ";
        for(auto &x : e.second.second) cout<<id2str[x]<<", ";
        cout<<")."<<endl;
        int allmissed = 0;
        for(auto &x : e.second.first) {
            vector<ui> xnei, nnei;
            for(auto &w : e.second.first) {
                if(s_G[x].find(w) == s_G[x].end()) {
                    nnei.push_back(w);
                    ++ allmissed;
                }
                else if (s_G[x][w] == -1) {
                    xnei.push_back(w);
                    ++ allmissed;
                }
            }
            for(auto &w : e.second.second) {
                if(s_G[x].find(w) == s_G[x].end()) {
                    nnei.push_back(w);
                    ++ allmissed;
                }
                else if (s_G[x][w] == 1) {
                    xnei.push_back(w);
                    ++ allmissed;
                }
            }
            cout<<"\t\t "<<id2str[x]<<"["<<nnei.size()<<","<<xnei.size()<<"]: N(";
            for(auto s : nnei) cout<<id2str[s]<<","; cout<<"). X(";
            for(auto s : xnei) cout<<id2str[s]<<","; cout<<")."<<endl;
        }
        for(auto &x : e.second.second) {
            vector<ui> xnei, nnei;
            for(auto &w : e.second.first) {
                if(s_G[x].find(w) == s_G[x].end()) {
                    nnei.push_back(w);
                    ++ allmissed;
                }
                else if (s_G[x][w] == 1) {
                    xnei.push_back(w);
                    ++ allmissed;
                }
            }
            for(auto &w : e.second.second) {
                if(s_G[x].find(w) == s_G[x].end()) {
                    nnei.push_back(w);
                    ++ allmissed;
                }
                else if (s_G[x][w] == -1) {
                    xnei.push_back(w);
                    ++ allmissed;
                }
            }
            cout<<"\t\t "<<id2str[x]<<"["<<nnei.size()<<","<<xnei.size()<<"]: N(";
            for(auto s : nnei) cout<<id2str[s]<<","; cout<<"). X(";
            for(auto s : xnei) cout<<id2str[s]<<","; cout<<")."<<endl;
        }
        cout<<"\t\t *** all missed : "<<allmissed<<endl;
    }
    
    delete [] s_G;
}

int main(int argc, const char * argv[]) {
    if(argc < 5){ cout<<"\n Usage: [0]exe [1]input_graph [2]k [3]tau [4]algo \n"; exit(1); }
    k = atoi(argv[2]);
    if(k<2) {cout<<"k should be >= 2."<<endl; exit(1);}
    tau = atoi(argv[3]);
    if(tau<k) {cout<<"tau should be >= k."<<endl; exit(1);}
    algo = argv[4];
    cout<<"[G:"<<argv[1]<<", k:"<<k<<", tau:"<<tau<<", algo:"<<algo<<"]"<<endl;
    load_graph(argv[1]);
    Timer t;
    if(algo.compare("B") == 0)
        maximal_balanced_kplex_enum();
    else if (algo.compare("A") == 0)
        maximal_balanced_kplex_enum_adv_opt();
    else {
        cout<<"no matched algorithm!"<<endl; exit(1);
    }
    cout<<"Result Number: "<<Res_num<<" <"<<min_kplex<<"-"<<max_kplex<<">"<<endl;
    cout<<"Time Cost: "<<integer_to_string(t.elapsed())<<endl;
#ifdef _CheckRes_
    check_results(argv[1]);
#endif
#ifdef _PrintRes_
    display_results();
#endif
#ifdef _CaseStudy_
    #ifndef _StoreRes_
        cout<<"result must be stored before case study!"<<endl; exit(1);
    #endif
    ShowCaseStudyResult(argv[1]);
#endif
    delete_memo();
    return 0;
}
