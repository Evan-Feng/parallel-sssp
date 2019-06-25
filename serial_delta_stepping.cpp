/*********************************************************************
 *  Name: serial_delta_stepping.cpp                                  *
 *  Description: solve the single-source shortest path problem       *
 *               using the delta stepping algorithm                  *
 *  Usage: ./serial <graph file> <aux file> <output file>            *
 *  Author: fengyanlin@pku.edu.cn                                    *
 *********************************************************************
 */
#include <unordered_set>
#include <vector>
#include <queue>
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


#define DEBUG
// #define CHECK_DISCONNECTED


/*
 *  delta_stepping - solve the single source shotest path problem using the delta stepping algorithm
 *
 *  graph: the graph represented by a list of vertices
 *  source: ID of the source point
 *  
 *  returns: the distances between each vectex and the source point 
 */

void delta_stepping(vector<vector<pair<int, long long> > > & graph, vector<int> & nlight, int source, 
                    vector<long long> & dist, long long delta, long long max_dist, vector<vector<int> > & B, vector<int> & S, vector<int> & R){
    fill(dist.begin(), dist.end(), numeric_limits<long long>::max());
    int nbuckets = max_dist / delta;

#ifdef DEBUG
    cout << "running serial delta-stepping" << endl;
#endif

    dist[source] = 0;
    B[0].push_back(source);
    int bid = 0;
    while (1){
        while(bid < nbuckets && B[bid].empty())
            bid++;
        if (bid >= nbuckets)
            break;
        while (!B[bid].empty()){
            for (int v: B[bid]){
                long long dv = dist[v];
                int max_j = bid == nbuckets - 1 ? graph[v].size() : nlight[v];
                for (int j = 0; j < max_j; j++){
                    int tv = graph[v][j].first;
                    long long w = graph[v][j].second;
                    long long dtv = dv + w;
                    if (dtv < dist[tv]){
                        dist[tv] = dtv;
                        R.push_back(tv);
                    }
                }
                S.push_back(v);
            }
            B[bid].clear();

            for (int v : R){
                int dest = dist[v] / delta;
                if (dest >= nbuckets)
                    dest = nbuckets - 1;
                if (dest >= bid)
                    B[dest].push_back(v);
            }
            R.clear();
        }
        if (bid == nbuckets - 1){
            S.clear();
            break;
        }
        for (int v: S){
            long long dv = dist[v];
            for (int j = nlight[v]; j < graph[v].size(); j++){
                int tv = graph[v][j].first;
                long long w = graph[v][j].second;
                long long dtv = dv + w;
                if (dtv < dist[tv]){
                    dist[tv] = dtv;
                    R.push_back(tv);
                }
            }
        }
        S.clear();
        for (int v: R){
            int dest = dist[v] / delta;
            if (dest >= nbuckets)
                dest = nbuckets - 1;
            if (dest > bid)
                B[dest].push_back(v);
        }
        R.clear();
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

    long long delta = 20000;
    long long max_dist = 4e8;


    // disable sync with stdio
    ios_base::sync_with_stdio(false);

    if (argc != 4){
        printf("Usage: ./%s <graph file> <aux file> <output file>\n", argv[0]);
        exit(0);
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
    int nbuckets = max_dist / delta;
    vector<vector<int> > B(nbuckets);
    vector<int> S;
    vector<int> R;

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
            clock_t begin = clock();
            delta_stepping(graph, nlight, src - 1, dist, delta, max_dist, B, S, R);
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
            clock_t end = clock();
            cout << "serial delta-stepping takes " << (double)(end - begin) / CLOCKS_PER_SEC << " seconds" << endl;
            outfile << "d " << checksum << endl;
        }
    }
    ssfile.close();
    outfile.close();

    return 0;
}
