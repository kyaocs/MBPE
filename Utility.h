//
//  Utility.h
//  maximal_balanced_kplex_enum
//
//  Created by kai on 2023/3/21.
//

#ifndef Utility_h
#define Utility_h

#include <iostream>
#include <fstream>
#include <map>
#include <unordered_map>
#include <set>
#include <unordered_set>
#include <vector>
#include <queue>
#define NDEBUG
#include <assert.h>
#include <algorithm>
#include <string.h>
#include <iomanip>
#include "Timer.h"
#include <sparsehash/dense_hash_map>

using namespace std;
using namespace google;
int * depth_distr;

//#define _StoreRes_
//#define _CheckRes_
//#define _PrintRes_
//#define _CostlyAssert_
//#define _DEBUG_
//#define _ShowProgress_
//#define _Statistic_
//#define _onlyShowG_
//#define _CaseStudy_
//#define _ET_
//#define _UB_
//#define _maintainX_
//#define _ERinEnum_

#define _CTprune_
#define _ParVR_
#define _PP_
#define _VRinEnum_
#define _PIVOT_
#define _mindegpivot_

typedef unsigned int ui;
const int INF = 1000000000;
#define miv(a,b) ((a)>(b)?(b):(a))
#define mav(a,b) ((a)<(b)?(b):(a))
int er_threshold = 2;
long long total_subg_inenum = 0;
double total_den = 0.0;

ui n; //number of vertices
ui m; //number of edges
ui pm; //number of positive edges
ui nm; //number of negative edges
ui max_deg;
ui max_core;
int tau; // size threshold
int k; //k-plex
string algo; //algorithm
unordered_map<ui, string> id2str;

ui * p_pstart;//positive
ui * p_edges;//positive
ui * p_pend;
ui * n_pstart;//negative
ui * n_edges;//negative
ui * n_pend;

ui * degree;
ui * p_degree;
ui * n_degree;

ui * process_order;
ui * core;
ui * ver_rank;
bool * ver_del;
short * inQv;
ui * original_id;
ui * true_id;
short * CNT;
bool * inKL;
bool * inKR;
bool * inCL;
bool * inCR;
bool * inXL;
bool * inXR;
bool * inK;
//ui * rdeg_inK;
short ** M;
//short * state;
short * part;
vector<ui> A;
ui M_len;
//ui * tran;
short * partial_kplex_cap;
bool * skipv_CL;
bool * skipv_CR;

bool * addtoKL;
bool * addtoKR;

ui * pos_inCL;
ui * pos_inCR;
ui * pos_inXL;
ui * pos_inXR;

ui * t_p_pstart;//positive
ui * t_p_edges;//positive
ui * t_p_pend;
ui * t_n_pstart;//negative
ui * t_n_edges;//negative
ui * t_n_pend;

ui * t_p_mark;
ui * t_n_mark;

ui * L_deg_inK;
ui * R_deg_inK;

int * pdL;
int * pdR;
int * ndL;
int * ndR;

int * XL_pdL;
int * XL_pdR;
int * XL_ndL;
int * XL_ndR;
int * XR_pdL;
int * XR_pdR;
int * XR_ndL;
int * XR_ndR;

int * XL_pdXL;
int * XL_ndXR;
int * XR_pdXR;
int * XR_ndXL;

bool * inSL;
bool * inSR;
bool * inSC;
bool * inP;

/* seed graph method*/
vector<ui> sKL, sKR, sCL, sCR, sP, sXL, sXR;
ui * pos;

dense_hash_map<ui, int> * sup_pp;
dense_hash_map<ui, int> * sup_nn;
dense_hash_map<ui, int> * sup_np;
dense_hash_map<ui, int> * e_sign;
dense_hash_map<ui, bool> * e_del;
dense_hash_map<ui, bool> * e_inQ;

long long total_rounds = 0;
long long total_C_size = 0;
long long total_X_size = 0;

long long total_high_rounds = 0;
long long total_high_C_size = 0;
long long total_high_X_size = 0;

long long total_mid_rounds = 0;
long long total_mid_C_size = 0;
long long total_mid_X_size = 0;

long long total_low_rounds = 0;
long long total_low_C_size = 0;
long long total_low_X_size = 0;

/*ER in Enum*/
vector<vector<ui>> Pnei, Nnei, Pnnei, Nnnei;
ui * p_r_degree;
ui * n_r_degree;
int ** tri_pp;
int ** tri_nn;
int ** tri_pn;
int ** tri_np;
short ** inQe;
bool * v_sta;
long long er_round = 0;
long long er_prune_cnt = 0;
long long vr_prune_I = 0;
long long vr_prune_II = 0;
//vector<vector<ui>> combs;

//optimize pivot
vector<vector<ui>> nneiKL, nneiKR;


vector<pair<vector<ui>, vector<ui>>> Res;
long long Res_num = 0;
int max_kplex = 0;
int min_kplex = INF;
long long branch_num = 0;
long long invoke_enum = 0;

Timer timer_in_enum;
long long T_u_toK = 0, T_update_deg = 0, T_shrink_CX = 0, T_recons_sg = 0, T_recover_sg = 0, T_recover_CX = 0, T_recover_deg = 0, T_u_toX = 0, T_final_recover = 0, T_prepare_CXK = 0, T_VR_in_Enum = 0, T_ER_in_Enum = 0, T_pivot = 0, T_ET = 0, T_check_max = 0;

long long T_init = 0, T_CX_M = 0, T_enum = 0;

bool mycomp(ui i, ui j)
{
    return core[i]<core[j];
}

bool mycomp_r(ui i, ui j)
{
    return core[i]>core[j];
}

string integer_to_string(long long number) {
    std::vector<ui> sequence;
    if(number == 0) sequence.push_back(0);
    while(number > 0) {
        sequence.push_back(number%1000);
        number /= 1000;
    }
    
    char buf[5];
    std::string res;
    for(unsigned int i = (ui)sequence.size();i > 0;i --) {
        if(i == sequence.size()) sprintf(buf, "%u", sequence[i-1]);
        else sprintf(buf, ",%03u", sequence[i-1]);
        res += std::string(buf);
    }
    return res;
}

void check_results(string input_graph)
{
    for(auto &e : Res) {
        sort(e.first.begin(), e.first.end());
        sort(e.second.begin(), e.second.end());
        vector<ui> tmp_vec;
//        cout<<e.first[0]<<","<<e.second[0]<<endl;
        if(e.first[0] > e.second[0]) {
            tmp_vec = e.first;
            e.first = e.second;
            e.second = tmp_vec;
        }
    }
    sort(Res.begin(),Res.end());
    //read the original graph
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
        input_file.close();
        
        if(p_pstart != nullptr) delete [] p_pstart;
        if(p_edges != nullptr) delete [] p_edges;
        if(n_pstart != nullptr) delete [] n_pstart;
        if(n_edges != nullptr) delete [] n_edges;
        if(inCL != nullptr) delete [] inCL;
        if(inCR != nullptr) delete [] inCR;
        if(ver_rank != nullptr) delete [] ver_rank;
        if(partial_kplex_cap != nullptr) delete [] partial_kplex_cap;
        
        p_pstart = new ui[n+1];
        p_edges = new ui[2*pm];
        n_pstart = new ui[n+1];
        n_edges = new ui[2*nm];
        inCL = new bool[n];
        inCR = new bool[n];
        ver_rank = new ui[n];
        partial_kplex_cap = new short[n];
        
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
        }
        assert(n_pstart[n] == 2*nm);
        delete [] s_G;
    }
    bool no_duplicate = true;
    int duplicate_num = 0;
    //check if there are duplicates
    for(ui i = 0; i < Res.size(); i++) {
        ui j = i + 1;
        if(j == Res.size()) break;
        auto K1 = Res[i]; auto K2 = Res[j];
//        cout<<"compare ";
//        cout<<"("; for(auto x : K1.first) cout<<x<<",";
//        cout<<"|"; for(auto x : K1.second) cout<<x<<","; cout<<")   vs   ";
//        cout<<"("; for(auto x : K2.first) cout<<x<<",";
//        cout<<"|"; for(auto x : K2.second) cout<<x<<","; cout<<")"<<endl;
        if(K1 == K2) {
            no_duplicate = false;
            ++ duplicate_num;
        }
    }
    cout<<"no duplicates:";
    if(no_duplicate) cout<<"YES"<<endl;
    else cout<<"NO("<<duplicate_num<<")"<<endl;
    
    //check correctness of each k-plex (ensure pstart&edges are original)
    bool each_kplex_is_correct = true;
    memset(inCL, 0, sizeof(bool)*n);
    memset(inCR, 0, sizeof(bool)*n);
    for(auto C : Res) {
        for(auto e : C.first) inCL[e] = 1;
        for(auto e : C.second) {
            if(inCL[e] == 1) {cout<<"inCL[e] == 1"<<endl; exit(1);} //ensure no same vertex in CL and CR
            inCR[e] = 1;
        }
        if(C.first.size() < tau) each_kplex_is_correct = false;
        if(C.second.size() < tau) each_kplex_is_correct = false;
        ui C_size = (ui)C.first.size() + (ui)C.second.size();
        for(auto u : C.first) {
//            cout<<" @ "<<u<<endl;
            ui tmp_deg = 0;
            for(ui i = p_pstart[u]; i < p_pstart[u+1]; i++) {
                ui v = p_edges[i];
//                cout<<"     cehck its pnei: "<<v<<endl;
                if(inCL[v] == 1) ++ tmp_deg;
            }
            for(ui i = n_pstart[u]; i < n_pstart[u+1]; i++) {
                ui v = n_edges[i];
//                cout<<"     cehck its nnei: "<<v<<endl;
                if(inCR[v] == 1) ++ tmp_deg;
            }
            if(tmp_deg >= C_size) {cout<<"tmp_deg >= C_size"<<endl; exit(1);}
            if(C_size - tmp_deg > k) {cout<<"gap > k"<<endl; each_kplex_is_correct = false;}
//            cout<<"Csize - tmp_deg = "<<C_size - tmp_deg<<endl;
        }
        for(auto u : C.second) {
//            cout<<" @ "<<u<<endl;
            ui tmp_deg = 0;
            for(ui i = p_pstart[u]; i < p_pstart[u+1]; i++) {
                ui v = p_edges[i];
//                cout<<"     cehck its pnei: "<<v<<endl;
                if(inCR[v] == 1) ++ tmp_deg;
            }
            for(ui i = n_pstart[u]; i < n_pstart[u+1]; i++) {
                ui v = n_edges[i];
//                cout<<"     cehck its nnei: "<<v<<endl;
                if(inCL[v] == 1) ++ tmp_deg;
            }
            if(tmp_deg >= C_size) {cout<<"tmp_deg >= C_size"<<endl; exit(1);}
            if(C_size - tmp_deg > k) {cout<<"gap > k"<<endl; each_kplex_is_correct = false;}
//            cout<<"Csize - tmp_deg = "<<C_size - tmp_deg<<endl;
        }
        for(auto e : C.first) inCL[e] = 0;
        for(auto e : C.second) inCR[e] = 0;
    }
    cout<<"correctness:";
    if(each_kplex_is_correct) cout<<"YES"<<endl;
    else cout<<"NO"<<endl;
    
    //check the maximality of each kplex
    bool each_kplex_is_maximal = true;
    memset(inCL, 0, sizeof(bool)*n);
    memset(inCR, 0, sizeof(bool)*n);
    memset(ver_rank, 0, sizeof(ui)*n); //reuse ver_rank
    for(auto C : Res) {
        
//        cout<<"checking ("; for(auto x : C.first) cout<<x<<",";
//        cout<<"|"; for(auto x : C.second) cout<<x<<","; cout<<")"<<endl;
        
        for(auto e : C.first) inCL[e] = 1;
        for(auto e : C.second) inCR[e] = 1;
        //check if a neighbor of k-plex C can be added in.
        unordered_set<ui> check_list;
        for(auto u : C.first) {
            for(ui i = p_pstart[u]; i < p_pstart[u+1]; i++)
                if(inCL[p_edges[i]] == 0 && inCR[p_edges[i]] == 0)
                    check_list.insert(p_edges[i]);
            for(ui i = n_pstart[u]; i < n_pstart[u+1]; i++)
                if(inCL[n_edges[i]] == 0 && inCR[n_edges[i]] == 0)
                    check_list.insert(n_edges[i]);
        }
        for(auto u : C.second) {
            for(ui i = p_pstart[u]; i < p_pstart[u+1]; i++)
                if(inCL[p_edges[i]] == 0 && inCR[p_edges[i]] == 0)
                    check_list.insert(p_edges[i]);
            for(ui i = n_pstart[u]; i < n_pstart[u+1]; i++)
                if(inCL[n_edges[i]] == 0 && inCR[n_edges[i]] == 0)
                    check_list.insert(n_edges[i]);
        }
//        cout<<"check_list : "; for(auto e : check_list) cout<<e<<",";cout<<endl;
        ui C_size = (ui)C.first.size() + (ui)C.second.size();
        for(auto u : C.first) {
            ui tmp_deg = 0;
            for(ui i = p_pstart[u]; i < p_pstart[u+1]; i++) {
                ui v = p_edges[i];
                if(inCL[v] == 1) ++ tmp_deg;
            }
            for(ui i = n_pstart[u]; i < n_pstart[u+1]; i++) {
                ui v = n_edges[i];
                if(inCR[v] == 1) ++ tmp_deg;
            }
            if(tmp_deg >= C_size) {cout<<"tmp_deg >= C_size"<<endl; exit(1);}
            if(C_size - tmp_deg > k) {cout<<"C_size - tmp_deg > k"<<endl; exit(1);}
            partial_kplex_cap[u] = k - (C_size - tmp_deg);
//            cout<<u<<" cap = "<< partial_kplex_cap[u]<<endl;
        }
        for(auto u : C.second) {
            ui tmp_deg = 0;
            for(ui i = p_pstart[u]; i < p_pstart[u+1]; i++) {
                ui v = p_edges[i];
                if(inCR[v] == 1) ++ tmp_deg;
            }
            for(ui i = n_pstart[u]; i < n_pstart[u+1]; i++) {
                ui v = n_edges[i];
                if(inCL[v] == 1) ++ tmp_deg;
            }
            if(tmp_deg >= C_size) {cout<<"tmp_deg >= C_size"<<endl; exit(1);}
            if(C_size - tmp_deg > k) {cout<<"C_size - tmp_deg > k"<<endl; exit(1);}
            partial_kplex_cap[u] = k - (C_size - tmp_deg);
//            cout<<u<<" cap = "<< partial_kplex_cap[u]<<endl;
        }
        for(auto u : check_list) {
            //record u's neighbors
            for(ui i = p_pstart[u]; i < p_pstart[u+1]; i++) ver_rank[p_edges[i]] = 1;
            for(ui i = n_pstart[u]; i < n_pstart[u+1]; i++) ver_rank[n_edges[i]] = 2;
//            cout<<"@ "<<u<<endl;
            //add u to CL
//            cout<<"add to CL"<<endl;
            ui tmp_missed = 1;
            bool can_be_added = true;
            for(auto e : C.first) {
                if(ver_rank[e] != 1) {
//                    cout<<e<<" in CL xxx"<<endl;
                    ++ tmp_missed;
                    if(partial_kplex_cap[e] <= 0) can_be_added = false;
                }
            }
            for(auto e : C.second) {
                if(ver_rank[e] != 2) {
//                    cout<<e<<" in CR xxx"<<endl;
                    ++ tmp_missed;
                    if(partial_kplex_cap[e] <= 0) can_be_added = false;
                }
            }
//            cout<<"tmp_missed = "<<tmp_missed<<endl;
            if(can_be_added == true && tmp_missed <= k) each_kplex_is_maximal = false;
            //add u to CR
//            cout<<"add to CR"<<endl;
            tmp_missed = 1;
            can_be_added = true;
            for(auto e : C.first) {
                if(ver_rank[e] != 2) {
//                    cout<<e<<" in CL xxx"<<endl;
                    ++ tmp_missed;
                    if(partial_kplex_cap[e] <= 0) can_be_added = false;
                }
            }
            for(auto e : C.second) {
                if(ver_rank[e] != 1) {
//                    cout<<e<<" in CR xxx"<<endl;
                    ++ tmp_missed;
                    if(partial_kplex_cap[e] <= 0) can_be_added = false;
                }
            }
//            cout<<"tmp_missed = "<<tmp_missed<<endl;
            if(can_be_added == true && tmp_missed <= k) each_kplex_is_maximal = false;
            for(ui i = p_pstart[u]; i < p_pstart[u+1]; i++) ver_rank[p_edges[i]] = 0;
            for(ui i = n_pstart[u]; i < n_pstart[u+1]; i++) ver_rank[n_edges[i]] = 0;
        }
        for(auto e : C.first) inCL[e] = 0;
        for(auto e : C.second) inCR[e] = 0;
    }//for each C in Res
    cout<<"maximality:";
    if(each_kplex_is_maximal) cout<<"YES"<<endl;
    else cout<<"NO"<<endl;
}

void delete_memo()
{
    if(p_pstart != nullptr) {
        delete [] p_pstart; p_pstart = nullptr;
    }
    if(p_edges != nullptr) {
        delete [] p_edges; p_edges = nullptr;
    }
    if(n_pstart != nullptr) {
        delete [] n_pstart; n_pstart = nullptr;
    }
    if(n_edges != nullptr) {
        delete [] n_edges; n_edges = nullptr;
    }
    if(degree != nullptr) {
        delete [] degree; degree = nullptr;
    }
    if(p_degree != nullptr) {
        delete [] p_degree; p_degree = nullptr;
    }
    if(n_degree != nullptr) {
        delete [] n_degree; n_degree = nullptr;
    }
    if(process_order != nullptr) {
        delete [] process_order; process_order = nullptr;
    }
    if(core != nullptr) {
        delete [] core; core = nullptr;
    }
    if(ver_rank != nullptr) {
        delete [] ver_rank; ver_rank = nullptr;
    }
    if(ver_del != nullptr) {
        delete [] ver_del; ver_del = nullptr;
    }
    if(original_id != nullptr) {
        delete [] original_id; original_id = nullptr;
    }
    if(inKL != nullptr) {
        delete [] inKL; inKL = nullptr;
    }
    if(inKR != nullptr) {
        delete [] inKR; inKR = nullptr;
    }
    if(inCL != nullptr) {
        delete [] inCL; inCL = nullptr;
    }
    if(inCR != nullptr) {
        delete [] inCR; inCR = nullptr;
    }
    if(inXL != nullptr) {
        delete [] inXL; inXL = nullptr;
    }
    if(inXR != nullptr) {
        delete [] inXR; inXR = nullptr;
    }
    for(int i = 0; i < M_len; i++) if(M[i] != nullptr){
        delete [] M[i];
        M[i] = nullptr;
    }
    if(partial_kplex_cap != nullptr) {
        delete [] partial_kplex_cap; partial_kplex_cap = nullptr;
    }
    if(skipv_CL != nullptr) {
        delete [] skipv_CL; skipv_CL = nullptr;
    }
    if(skipv_CR != nullptr) {
        delete [] skipv_CR; skipv_CR = nullptr;
    }
    if(pos_inCL != nullptr) {
        delete [] pos_inCL; pos_inCL = nullptr;
    }
    if(pos_inCR != nullptr) {
        delete [] pos_inCR; pos_inCR = nullptr;
    }
    if(pos_inXL != nullptr) {
        delete [] pos_inXL; pos_inXL = nullptr;
    }
    if(pos_inXR != nullptr) {
        delete [] pos_inXR; pos_inXR = nullptr;
    }
    if(t_p_pstart != nullptr) {
        delete [] t_p_pstart; t_p_pstart = nullptr;
    }
    if(t_p_edges != nullptr) {
        delete [] t_p_edges; t_p_edges = nullptr;
    }
    if(t_p_pend != nullptr) {
        delete [] t_p_pend; t_p_pend = nullptr;
    }
    if(t_n_pstart != nullptr) {
        delete [] t_n_pstart; t_n_pstart = nullptr;
    }
    if(t_n_edges != nullptr) {
        delete [] t_n_edges; t_n_edges = nullptr;
    }
    if(t_n_pend != nullptr) {
        delete [] t_n_pend; t_n_pend = nullptr;
    }
    if(t_p_mark != nullptr) {
        delete [] t_p_mark; t_p_mark = nullptr;
    }
    if(t_n_mark != nullptr) {
        delete [] t_n_mark; t_n_mark = nullptr;
    }
//    if(L_deg_inK != nullptr) {
//        delete [] L_deg_inK; L_deg_inK = nullptr;
//    }
//    if(R_deg_inK != nullptr) {
//        delete [] R_deg_inK; R_deg_inK = nullptr;
//    }
    if(CNT != nullptr) {
        delete [] CNT; CNT = nullptr;
    }
}


//    cout<<"depthdistr :"<<endl;
//    long long depth_total = 0;
//    long long idx = 0;
//    vector<pair<int , int>> tmpvec;
//    for(ui i = 1; i < x; i++) {
//        if(depth_distr[i]!=0) {
//            depth_total += (depth_distr[i]*i);
//            tmpvec.push_back(make_pair(i, depth_distr[i]));
//        }
//    }
//    cout<<"ave depth  = "<<(double)depth_total/branch_num<<endl;
//    for(auto e : tmpvec) {
//        cout<<e.first<<":"<<e.second<<", ";
//    }
//    cout<<"+++++++++++++++++++++++++++++++++++++"<<endl;
//    assert(total_rounds == branch_num);
//    cout<<"high portion = "<<(double)total_high_rounds/total_rounds<<endl;
//    cout<<"mid portion = "<<(double)total_mid_rounds/total_rounds<<endl;
//    cout<<"low portion = "<<(double)total_low_rounds/total_rounds<<endl;
//
//    if(total_high_rounds != 0) {
//        cout<<"ave C size in high: "<<total_high_C_size/total_high_rounds<<endl;
//        cout<<"ave X size in high: "<<total_high_X_size/total_high_rounds<<endl;
//
//    }
//    if(total_mid_rounds != 0) {
//        cout<<"ave C size in mid: "<<total_mid_C_size/total_mid_rounds<<endl;
//        cout<<"ave X size in mid: "<<total_mid_X_size/total_mid_rounds<<endl;
//
//    }
//    if(total_low_rounds != 0) {
//        cout<<"ave C size in low: "<<total_low_C_size/total_low_rounds<<endl;
//        cout<<"ave X size in low: "<<total_low_X_size/total_low_rounds<<endl;
//
//    }

#endif /* Utility_h */
