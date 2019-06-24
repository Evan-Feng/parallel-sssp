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
using namespace std;


#define DEBUG
// #define SERIAL
#define CHECK_DISCONNECTED


/*
 *  delta_stepping_openmp - solve the single source shotest path problem using the delta stepping algorithm
 *
 *  graph: the graph represented by a list of vertices
 *  source: ID of the source point
 *  
 *  returns: the distances between each vectex and the source point 
 */

void delta_stepping_openmp(vector<vector<pair<int, long long> > > & graph, vector<int> & nlight, int source, 
                           vector<long long> & dist, long long delta, long long max_dist,
                           int nprocs, vector<unordered_set<int> > & B, vector<map<int, long long> > & R,
                           vector<unordered_set<int> > & S, omp_lock_t * B_lock, omp_lock_t * R_lock){
    int nbuckets = max_dist / delta + 1;
    int nnodes = graph.size();
    int curr_bucket_empty[nprocs];
    int next_bid[nprocs];
    int curr_bucket_empty_flag = 0;

    fill(dist.begin(), dist.end(), numeric_limits<long long>::max());

#ifdef DEBUG
    cout << "running " << nprocs << " threads in parallel" << endl;
#endif

    dist[source] = 0;
    B[source % nprocs].insert(source);
    int bid = -1;

#pragma omp parallel
{
    int pid = omp_get_thread_num();

    while (1){
        next_bid[pid] = bid + 1;
        while (next_bid[pid] < nbuckets && B[next_bid[pid] * nprocs + pid].empty())
            next_bid[pid]++;
#pragma omp barrier
        if (pid == 0){
            curr_bucket_empty_flag = 0;
            bid = nbuckets;
        }
#pragma omp barrier
#pragma omp for reduction(min:bid)
        for (int i = 0; i < nprocs; i++)
            if (next_bid[i] < bid)
                bid = next_bid[i];
        if (bid >= nbuckets)
            break;

        while (!curr_bucket_empty_flag){
            for (auto v: B[bid * nprocs + pid]){
                long long dv = dist[v];
                int max_j = bid == nbuckets - 1 ? graph[v].size() : nlight[v];
                for (int j = 0; j < max_j; j++){
                    int tv = graph[v][j].first;
                    long long w = graph[v][j].second;
                    long long dtv = dv + w;

                    omp_set_lock(&R_lock[tv % nprocs]);
                    auto it = R[tv % nprocs].find(tv);
                    if (it == R[tv % nprocs].end())
                        R[tv % nprocs][tv] = dtv;
                    else if (dtv < it->second)
                        it->second = dtv;
                    omp_unset_lock(&R_lock[tv % nprocs]);
                }
            }
            S[pid].insert(B[bid * nprocs + pid].begin(), B[bid * nprocs + pid].end());
            B[bid * nprocs + pid].clear();

#pragma omp barrier

            for (pair<int, long long> edge : R[pid]){
                int v = edge.first;
                long long dv = edge.second;
                if (dv < dist[v]){
                    int from = min(dist[v] / delta, (long long)nbuckets - 1);
                    int to = min(dv / delta, (long long)nbuckets - 1);
                    from = from * nprocs + pid;  // pid == v % nprocs
                    to = to * nprocs + pid;
                    if (B[from].find(v) == B[from].end())
                        B[to].insert(v);
                    else if (from != to){
                        B[from].erase(v);
                        B[to].insert(v);
                    }
                    dist[v] = dv;
                }
            }
            R[pid].clear();

            if (pid == 0)
                curr_bucket_empty_flag = 1;
#pragma omp barrier
            curr_bucket_empty[pid] = B[bid * nprocs + pid].empty();
#pragma omp for reduction(min:curr_bucket_empty_flag)
            for (int i = 0; i < nprocs; i++)
                if (curr_bucket_empty[i] < curr_bucket_empty_flag)
                    curr_bucket_empty_flag = curr_bucket_empty[i];
        }
        R[pid].clear();
#pragma omp barrier
        for (int v: S[pid]){
            long long dv = dist[v];
            for (int j = nlight[v]; j < graph[v].size(); j++){
                int tv = graph[v][j].first;
                long long w = graph[v][j].second;
                long long dtv = w + dv;

                omp_set_lock(&R_lock[tv % nprocs]);
                auto it = R[tv % nprocs].find(tv);
                if (it == R[tv % nprocs].end())
                    R[tv % nprocs][tv] = dtv;
                else if (dtv < it->second)
                    it->second = dtv;
                omp_unset_lock(&R_lock[tv % nprocs]);
            }
        }
        S[pid].clear();
#pragma omp barrier
        for (pair<int, long long> edge : R[pid]){
            int v = edge.first;
            long long dv = edge.second;
            if (dv < dist[v]){
                int from = min(dist[v] / delta, (long long)nbuckets - 1);
                int to = min(dv / delta, (long long)nbuckets - 1);
                from = from * nprocs + pid;  // pid == v % nprocs
                to = to * nprocs + pid;
                if (B[from].find(v) == B[from].end())
                    B[to].insert(v);
                else if (from != to){
                    B[from].erase(v);
                    B[to].insert(v);
                }
                dist[v] = dv;
            }
        }
        R[pid].clear();
#pragma omp barrier
    }
}
}


int main(int argc, char * argv[]){
    string buf;
    size_t len = 0;
    int nnodes = 0, nedges = 0;
    ifstream grfile, ssfile;
    ofstream outfile;
    int src, trg;
    long long weight;

    long long delta = 8000;
    long long max_dist = 4000000;

    // disable sync with stdio
    ios_base::sync_with_stdio(false);

    if (argc != 4){
        printf("Usage: ./%s <graph file> <aux file> <output file>\n", argv[0]);
        exit(0);
    }

    // initialize data structures
    int nprocs = omp_get_num_procs();
#ifdef SERIAL
    nprocs = 1;
#endif
    omp_set_num_threads(nprocs);
    int nbuckets = (max_dist / delta) + 1;
    omp_lock_t B_lock[nbuckets * nprocs];
    omp_lock_t R_lock[nprocs];
    vector<unordered_set<int> > B(nbuckets * nprocs);
    vector<map<int, long long> > R(nprocs);
    vector<unordered_set<int> > S(nprocs);

#pragma omp parallel
{
    int pid = omp_get_thread_num();
    for (int bid = 0; bid < nbuckets; bid++)
        omp_init_lock(&B_lock[bid * nprocs + pid]);
    omp_init_lock(&R_lock[pid]);
}

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

    // read the .ss file line by line and solve the corresponding sssp problem
    ssfile.open(argv[2]);
    outfile.open(argv[3], fstream::app);
    while (getline(ssfile, buf)){
        string op;
        if (buf[0] == 's'){
            long long checksum = 0;
            istringstream stream(buf);
            stream >> op >> src;
            delta_stepping_openmp(graph, nlight, src - 1, dist, delta, max_dist, nprocs, B, R, S, B_lock, R_lock);
            for (long long n: dist){
               
#ifdef CHECK_DISCONNECTED
                if (n < numeric_limits<long long>::max())
                    checksum += n;
                else
                    cout << "inf dist encountered" << endl;
#else
                checksum += n;
#endif
            }
            outfile << "d " << checksum << endl;
        }
    }
    ssfile.close();
    outfile.close();

    return 0;
}
