/**********************************************************************
 *  Name: delta_stepping_openmp.cpp                                   *
 *  Description: Delta-Stepping implementation in OpenMP              *
 *  Usage: ./delta_step_openmp <graph file> <aux file> <output file>  *
 *  Author: fengyanlin@pku.edu.cn                                     *
 **********************************************************************
 */
#include <omp.h>
#include <unordered_set>
#include <vector>
#include <queue>
#include <map>
#include <limits>
#include <iostream>
#include <algorithm>
#include <fstream>
#include <sstream>
#include <string>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
using namespace std;

#define HELLO { if (pid == 0) cout << "hello" << endl; }
#define DEBUG
// #define SERIAL
#define CHECK_DISCONNECTED

#define NTASK 80
#define SCHEDULE schedule(static, 2)

/*
 *  delta_stepping_openmp - solve the single source shotest path problem using the delta stepping algorithm
 *
 *  graph: the graph represented by a list of vertices
 *  source: ID of the source point
 *  
 *  returns: the sum of all the distances
 */

long long delta_stepping_openmp(vector<vector<pair<int, long long> > > & graph, vector<int> & nlight, int source,
                           vector<long long> & dist, long long delta, long long max_dist,
                           int nprocs, vector<unordered_set<int> > & B, vector<map<int, long long> > & R,
                           vector<unordered_set<int> > & S){
    int nbuckets = max_dist / delta + 1;
    int nnodes = graph.size();
    int next_bid[NTASK];
    int curr_bucket_empty_flag = 0;
    long long checksum = 0;

    fill(dist.begin(), dist.end(), numeric_limits<long long>::max());

#ifdef DEBUG
    cout << "running delta-stepping with [delta = " << delta << "] [max_dist = "
         << max_dist << "] [num_threads = " << nprocs << "] [src = " << source << "]" << endl;
#endif

#ifdef DEBUG
    int exit_flag = 0;
#pragma omp parallel for
    for (int tid = 0; tid < NTASK; tid++)
        for (int bi = 0; bi < nbuckets; bi++)
            if (!B[bi * NTASK + tid].empty()){
                exit_flag = 1;
            }
    if (exit_flag){
        cout << "bucket not empty" << endl;
        exit(0);
    }
#endif

    dist[source] = 0;
    B[source % NTASK].insert(source);
    int bid = -1;
    int nbid = -1;

#pragma omp parallel
{
    while (1){
#pragma omp single
{
        curr_bucket_empty_flag = 0;
        nbid = nbuckets;
}
#pragma omp for reduction(min:nbid) SCHEDULE
        for (int tid = 0; tid < NTASK; tid++){
            int next_bid = bid + 1;
            while (next_bid < nbuckets && B[next_bid * NTASK + tid].empty())
                next_bid++;
            if (next_bid < nbid)
                nbid = next_bid;
        }
        bid = nbid;
        if (bid >= nbuckets)
            break;

        while (!curr_bucket_empty_flag){
#pragma omp for SCHEDULE
            for (int tid = 0; tid < NTASK; tid++){
                for (auto v: B[bid * NTASK + tid]){
                    long long dv = dist[v];
                    int max_j = bid == nbuckets - 1 ? graph[v].size() : nlight[v];
                    for (int j = 0; j < max_j; j++){
                        int tv = graph[v][j].first;
                        long long w = graph[v][j].second;
                        long long dtv = dv + w;
#ifdef DEBUG
                        if (w < 0)
                            cout << "negative weight encountered" << endl;
#endif

                        int r_dest = (tid * NTASK) + (tv % NTASK);
                        auto it = R[r_dest].find(tv);
                        if (it == R[r_dest].end())
                            R[r_dest][tv] = dtv;
                        else if (dtv < it->second)
                            it->second = dtv;
                    }
                }
                if (bid != nbuckets - 1)
                    S[tid].insert(B[bid * NTASK + tid].begin(), B[bid * NTASK + tid].end());
                B[bid * NTASK + tid].clear();
            }
            // merge light requests
            // for (int i = 1; i < nprocs; i++){
            //     int to_merge = i * nprocs + pid;
            //     R[pid].insert(R[to_merge].begin(), R[to_merge].end());
            //     R[to_merge].clear();
            // }
#pragma omp for SCHEDULE
            for (int tid = 0; tid < NTASK; tid++){ 
                for (int i = 0; i < NTASK; i++){
                    for (pair<int, long long> edge : R[i * NTASK + tid]){
                        int v = edge.first;
                        long long dv = edge.second;
                        if (dv < dist[v]){
                            int from = min(dist[v] / delta, (long long)nbuckets - 1);
                            int to = min(dv / delta, (long long)nbuckets - 1);
#ifdef DEBUG
                            if (to < bid) ;
                                // cout << "[A] inserting into lower buckets" << endl;
#endif
                            from = from * NTASK + tid;  // tid == v % NTASK
                            to = to * NTASK + tid;
                            if (B[from].find(v) == B[from].end())
                                B[to].insert(v);
                            else if (from != to){
                                B[from].erase(v);
                                B[to].insert(v);
                            }
                            dist[v] = dv;
                        }
                    }
                    R[i * NTASK + tid].clear();
                }
            }

#pragma omp single
            curr_bucket_empty_flag = 1;
#pragma omp for reduction(min:curr_bucket_empty_flag) SCHEDULE
            for (int tid = 0; tid < NTASK; tid++)
                if (B[bid * NTASK + tid].empty() < curr_bucket_empty_flag)
                    curr_bucket_empty_flag = B[bid * NTASK + tid].empty();
        }
        // R[pid].clear();
        if (bid == nbuckets - 1)
            break;
#pragma omp for SCHEDULE
        for (int tid = 0; tid < NTASK; tid++){
            for (int v: S[tid]){
                long long dv = dist[v];
                for (int j = nlight[v]; j < graph[v].size(); j++){
                    int tv = graph[v][j].first;
                    long long w = graph[v][j].second;
                    long long dtv = w + dv;

                    int r_dest = (tid * NTASK) + (tv % NTASK);
                    auto it = R[r_dest].find(tv);
                    if (it == R[r_dest].end())
                        R[r_dest][tv] = dtv;
                    else if (dtv < it->second)
                        it->second = dtv;
                }
            }
            S[tid].clear();
        }
#pragma omp for SCHEDULE
        for (int tid = 0; tid < NTASK; tid++){
            for (int i = 0; i < NTASK; i++){
                for (pair<int, long long> edge : R[i * NTASK + tid]){
                    int v = edge.first;
                    long long dv = edge.second;
                    if (dv < dist[v]){
                        int from = min(dist[v] / delta, (long long)nbuckets - 1);
                        int to = min(dv / delta, (long long)nbuckets - 1);
#ifdef DEBUG
                        if (to <= bid)
                            cout << "[B] inserting into lower buckets " << bid << endl;
#endif
                        to = to * NTASK + tid;
                        if (B[from].find(v) == B[from].end())
                            B[to].insert(v);
                        else if (from != to){
                            B[from].erase(v);
                            B[to].insert(v);
                        }
                        dist[v] = dv;
                    }
                }
                R[i * NTASK + tid].clear();
            }
        }
    }
#pragma omp for reduction(+:checksum) SCHEDULE
    for (int v = 0; v < nnodes; v++)
#ifdef CHECK_DISCONNECTED
        if (dist[v] < numeric_limits<long long>::max())
            checksum += dist[v];
#else
        checksum += dist[v];
#endif
}
    return checksum;
}


int main(int argc, char * argv[]){
    string buf;
    size_t len = 0;
    int nnodes = 0, nedges = 0;
    ifstream grfile, ssfile;
    ofstream outfile;
    int src, trg;
    long long weight;

    long long delta = 20000;
    long long max_dist = 400000000;

    // disable sync with stdio
    ios_base::sync_with_stdio(false);

    if (argc != 5){
        printf("Usage: ./%s <graph file> <aux file> <output file> <delta>\n", argv[0]);
        exit(0);
    }

    delta = atoi(argv[4]);

    // set number of threads
    int nprocs = omp_get_num_procs();
#ifdef SERIAL
    nprocs = 1;
#endif
    omp_set_num_threads(nprocs);


    // read graph file
    grfile.open(argv[1]);
    while (getline(grfile, buf)){
        string op, tmp;
        if (buf[0] == 'c')
            continue;
        else if (buf[0] == 'p'){
            istringstream stream(buf);
            stream >> op >> tmp >> nnodes >> nedges;
        }
        else if (buf[0] == 'a'){
            istringstream stream(buf);
            stream >> op >> src >> trg >> weight;
            break;
        }
    }

    vector<vector<pair<int, long long> > > graph(nnodes);
    vector<long long> dist(nnodes);
    vector<int> nlight(nnodes, 0);

    graph[src - 1].push_back({trg - 1, weight});
    if (weight < delta)
        nlight[src - 1] += 1;
    for (int i = 0; i < nedges - 1; i++){
        string op;
        grfile >> op >> src >> trg >> weight;
        if (weight < delta){
            graph[src - 1].insert(graph[src - 1].begin() + nlight[src - 1], {trg - 1, weight});
            nlight[src - 1]++;
        }
        else
            graph[src - 1].push_back({trg - 1, weight});    
    }
    grfile.close();

    // initialize data structures
    if (nnodes < 1e6){
        max_dist = 1e7;
    }
    else if (nnodes < 1e7){
        max_dist = 4e7;
    }
    else {
        max_dist = 4e8;
    }
    int nbuckets = (max_dist / delta) + 1;
    vector<unordered_set<int> > B(nbuckets * NTASK);
    vector<map<int, long long> > R(NTASK * NTASK);
    vector<unordered_set<int> > S(NTASK);

    // read the .ss file line by line and solve the corresponding sssp problem
    ssfile.open(argv[2]);
    outfile.open(argv[3], fstream::app);
    while (getline(ssfile, buf)){
        string op;
        if (buf[0] == 's'){
            long long checksum = 0;
            istringstream stream(buf);
            stream >> op >> src;
            clock_t begin = clock();
            checksum = delta_stepping_openmp(graph, nlight, src - 1, dist, delta, max_dist, nprocs, B, R, S);
            clock_t end = clock();
            cout << "delta-stepping takes " << (double)(end - begin) / (CLOCKS_PER_SEC * nprocs) << " seconds" << endl;
            outfile << "d " << checksum << endl;
        }
    }
    ssfile.close();
    outfile.close();

    return 0;
}
