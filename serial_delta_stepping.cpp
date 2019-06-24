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
using namespace std;


// #define DEBUG
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
                    vector<long long> & dist, long long delta, long long max_dist){
    fill(dist.begin(), dist.end(), numeric_limits<long long>::max());

    int nbuckets = max_dist / delta + 1;
    vector<unordered_set<int> > B(nbuckets);

    dist[source] = 0;
    B[0].insert(source);
    for (int bid = 0; bid < nbuckets; bid++){
        unordered_set<int> S;
        vector<pair<int, long long> > R;
        while (!B[bid].empty()){
            for (int v: B[bid]){
                long long dv = dist[v];
                int max_j = bid == nbuckets - 1 ? graph[v].size() : nlight[v];
                for (int j = 0; j < max_j; j++){
                    int tv = graph[v][j].first;
                    long long w = graph[v][j].second;
                    R.push_back({tv, w + dv});
                }
            }
            S.insert(B[bid].begin(), B[bid].end());
            B[bid].clear();

            for (pair<int, long long> edge : R){
                int v = edge.first;
                long long dv = edge.second;

                if (dv < dist[v]){
                    dist[v] = dv;
                    int dest = dv / delta;
                    if (dest >= nbuckets)
                        dest = nbuckets - 1;
                    B[dest].insert(v);
                }
            }
            R.clear();
        }
        R.clear();
        for (int v: S){
            long long dv = dist[v];
            for (int j = nlight[v]; j < graph[v].size(); j++){
                int tv = graph[v][j].first;
                long long w = graph[v][j].second;

                R.push_back({tv, w + dv});
            }
        }
        for (pair<int, long long> edge : R){
            int v = edge.first;
            long long dv = edge.second;
            if (dv < dist[v]){
                int from = min(dist[v] / delta, (long long)nbuckets - 1);
                int to = min(dv / delta, (long long)nbuckets - 1);
                if (B[from].find(v) == B[from].end())
                    B[to].insert(v);
                else if (from != to){
                    B[from].erase(v);
                    B[to].insert(v);
                }
                dist[v] = dv;
            }
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
            delta_stepping(graph, nlight, src - 1, dist, delta, max_dist);
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
