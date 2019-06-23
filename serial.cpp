/*********************************************************************
 *  Name: serial.cpp                                                 *
 *  Description: solve the single-source shortest path problem       *
 *               using the dijkstras algorithm                       *
 *  Usage: ./serial <graph file> <aux file> <output file>            *
 *  Author: fengyanlin@pku.edu.cn                                    *
 *********************************************************************
 */
#include <vector>
#include <queue>
#include <limits>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <stdlib.h>
#include <stdio.h>
using namespace std;


/*
 *  dijkstras - solve the single source shotest path problem using the dijkstras algorithm
 *
 *  graph: the graph represented by a list of vertices
 *  source: ID of the source point
 *  
 *  returns: the distances between each vectex and the source point 
 */

void dijkstras(vector<vector<pair<int, long long> > > & graph, int source, vector<long long> & dist){
    vector<int> nodes;

    fill(dist.begin(), dist.end(), numeric_limits<long long>::max());
    dist[source] = 0;

    priority_queue<pair<long long, int>, vector<pair<long long, int> >, greater<pair<long long, int> > > q;
    q.push({0, source});
    while (!q.empty()){
        int v = q.top().second;
        long long dv = q.top().first;
        q.pop();
        if (dv != dist[v])
            continue;

        for (auto edge: graph[v]){
            int trg = edge.first;
            long long weight = edge.second;
            if (dist[v] + weight < dist[trg]){
                dist[trg] = dist[v] + weight;
                q.push({dist[trg], trg});
            }
        }
    }
}


int main(int argc, char * argv[]){
    string buf;
    int nnodes = 0, nedges = 0;
    ifstream grfile, ssfile;
    ofstream outfile;
    int src, trg;
    long long weight;

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
    graph[src - 1].push_back({trg - 1, weight});
    for (int i = 0; i < nedges - 1; i++){
        string op;
        grfile >> op >> src >> trg >> weight;
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
            dijkstras(graph, src - 1, dist);
            for (long long n: dist)
#ifdef CHECK_DISCONNECTED
                if (n < numeric_limits<long long>::max())
                    checksum += n;
                else
                    cout << "occurred INF dist" << endl;
#else
                checksum += n;
#endif
            outfile << "d " << checksum << endl;
        }
    }
    ssfile.close();
    outfile.close();

    return 0;
}
