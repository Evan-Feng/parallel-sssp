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

#define HELLO { if (omp_get_thread_num() == 0) cout << "hello" << endl; }
#define BID { if (omp_get_thread_num() == 0) cout << "bid = " << bid << endl; }
#define DEBUG
// #define CHECK_DISCONNECTED

#define NPROC 39
#define SCHEDULE schedule(static, 1000)

/*
 *  delta_stepping_openmp - solve the single source shotest path problem using the delta stepping algorithm
 *
 *  graph: the graph represented by a list of vertices
 *  source: ID of the source point
 *  
 *  returns: the sum of all the distances
 */
int nnodes, nedges;
long long delta, max_dist;
int curr_bucket_empty_flag = 0;
int local_bucket_empty[NPROC];
long long checksum = 0;
int bid = -1;
int exit_flag = 0;
vector<vector<pair<int, long long> > > graph;
int * B;
vector<vector<int> > R(NPROC);
vector<vector<int> > S(NPROC);
long long * dist;
vector<long long> nlight;
long long ll_inf = numeric_limits<long long>::max();
int int_inf = numeric_limits<int>::max();

// inline void delta_stepping_openmp(vector<vector<pair<int, long long> > > & graph, vector<int> & nlight, int source,
//                            vector<long long> & dist, long long delta, long long max_dist,
//                            int nprocs, vector<unordered_set<int> > & BB, vector<map<int, long long> > & RR,
//                            vector<unordered_set<int> > & SS){
void delta_stepping_openmp(int source){

    return;
}


int main(int argc, char * argv[]){
    string buf;
    size_t len = 0;
    ifstream grfile, ssfile;
    ofstream outfile;
    int src, trg;
    long long weight;
    int eof_flag;
    // long long delta = 20000;
    // long long max_dist = 400000000;

    // disable sync with stdio
    ios_base::sync_with_stdio(false);

    if (argc != 5){
        printf("Usage: ./%s <graph file> <aux file> <output file> <delta>\n", argv[0]);
        exit(0);
    }

    delta = atoi(argv[4]);

    // set number of threads
    omp_set_num_threads(NPROC);


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

    graph.resize(nnodes, vector<pair<int, long long> >(0));
    nlight.resize(nnodes, 0);
    dist = (long long *)malloc(nnodes * sizeof(long long));

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
    B = (int *)malloc(nnodes * sizeof(int));

#pragma omp parallel
{
    int pid = omp_get_thread_num();

    // read the .ss file line by line and solve the corresponding sssp problem
#pragma omp master
{
    ssfile.open(argv[2]);
    outfile.open(argv[3], fstream::app);
}
    while (1){
        string op;
#pragma omp barrier
#pragma omp master
{
        eof_flag = !getline(ssfile, buf);
}
#pragma omp barrier
        if (eof_flag)
            break;
        if (buf[0] == 's'){
            istringstream stream(buf);
            int src;
            stream >> op >> src;
            clock_t begin = clock();

/*------------------------------------------------------------------------*/
    int source = src - 1;
#pragma omp for SCHEDULE
    for (int v = 0; v < nnodes; v++){
        dist[v] = ll_inf;
        B[v] = int_inf;
    }

#ifdef DEBUG
#pragma omp single
{
    cout << "running parallel with [delta = " << delta << "] [max_dist = "
         << max_dist << "] [num_threads = " << NPROC << "] [src = " << source << "]" << endl;
}
#endif

#ifdef DEBUG
#pragma omp for SCHEDULE
    for (int pid = 0; pid < NPROC; pid++){
        if (!S[pid].empty()){
            exit_flag = 1;
        }
        if (!R[pid].empty()){
            exit_flag = 1;
        }
    }
    if (exit_flag){
        cout << "bucket not empty" << endl;
        exit(0);
    }
#endif

#pragma omp single
{
    dist[source] = 0;
    B[source] = 0;
    checksum = 0;
}

    // delta-stepping begins
    while (1){
#pragma omp single
{
        bid = int_inf;
        curr_bucket_empty_flag = 0;
}
#pragma omp for reduction(min:bid) SCHEDULE
        for (int v = 0; v < nnodes; v++)
            if (B[v] < bid)
                bid = B[v];

        if (bid == int_inf)
            break;

        // BID

        while (!curr_bucket_empty_flag){
#pragma omp for SCHEDULE
            for (int v = 0; v < nnodes; v++)
                if (B[v] == bid){
                    long long dv = dist[v];
                    for (int j = 0; j < nlight[v]; j++){
                        int tv = graph[v][j].first;
                        long long w = graph[v][j].second;
                        long long dtv = dv + w;
                        bool updated = false;
                        while (!updated){
                            long long old_dtv = dist[tv];
                            if (dtv < old_dtv)
                                updated = __sync_bool_compare_and_swap(&dist[tv], old_dtv, dtv);
                            else
                                break;
                        }
                        if (updated)
                            R[pid].push_back(tv);
                    }
                    B[v] = int_inf;
                    S[pid].push_back(v);
                }
            local_bucket_empty[pid] = 1;
            for (int v: R[pid]){
                long long dv = dist[v];
                int to_bid = dv / delta;
                if (to_bid >= bid)
                    B[v] = to_bid;
                if (to_bid == bid)
                    local_bucket_empty[pid] = 0;
            }
            R[pid].clear();
#pragma omp single
            curr_bucket_empty_flag = 1;
#pragma omp for reduction(min:curr_bucket_empty_flag)
            for (int pi = 0; pi < NPROC; pi++)
                if (local_bucket_empty[pi] < curr_bucket_empty_flag)
                    curr_bucket_empty_flag = local_bucket_empty[pi];
        }

        // handle heavy edges
        for (int v: S[pid]){
            long long dv = dist[v];
            for (int j = nlight[v]; j < graph[v].size(); j++){
                int tv = graph[v][j].first;
                long long w = graph[v][j].second;
                long long dtv = w + dv;
                bool updated = false;
                while (!updated){
                    long long old_dtv = dist[tv];
                    if (dtv < old_dtv)
                        updated = __sync_bool_compare_and_swap(&dist[tv], old_dtv, dtv);
                    else
                        break;
                }
                if (updated)
                    R[pid].push_back(tv);
            }
        }
        S[pid].clear();
#pragma omp barrier
        for (int v: R[pid]){
            long long dv = dist[v];
            int to_bid = dv / delta;
            if (to_bid > bid)
                B[v] = to_bid;
        }
        R[pid].clear();
    }
#pragma omp for reduction(+:checksum) SCHEDULE
    for (int v = 0; v < nnodes; v++)
#ifdef CHECK_DISCONNECTED
        if (dist[v] < ll_inf)
            checksum += dist[v];
#else
        checksum += dist[v];
#endif


/*------------------------------------------------------------------------*/

            clock_t end = clock();
#pragma omp master
{
            cout << "delta-stepping takes " << (double)(end - begin) / (CLOCKS_PER_SEC * NPROC) << " seconds" << endl;
            outfile << "d " << checksum << endl;
}
        }
    }
#pragma omp master
{
    ssfile.close();
    outfile.close();
}
}
    free(dist);
    free(B);

    return 0;
}
