/*********************************************************************
 *  Name: parallel.cpp                                               *
 *  Description: solve the single-source shortest path problem       *
 *               using the parallel delta stepping algorithm         *
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
#include <mpi.h>
using namespace std;


#define DEBUG
#define CHECK_DISCONNECTED


struct Edge{
    int trg;
    long long weight;
    Edge(int trg_, long long weight_): trg(trg_), weight(weight_){}
};


/*
 *  parallel_delta_stepping - solve the single source shotest path problem using the 
 *                            delta stepping algorithm
 *
 *  graph: the graph represented by a list of vertices
 *  source: ID of the source point
 *  
 *  returns: the distances between each vectex and the source point 
 */

void parallel_delta_stepping(int rank, int np, Edge * local_edges, int local_nnodes, int base_id, long long * dist,
                             int * local_offsets_light, int * local_offsets_heavy, int * local_nedges_each_vertex,
                             int source, long long * local_dist, long long delta, long long max_dist){
#ifdef DEBUG
    cout << "parallel_delta_stepping called by process " << rank << " with source " << source << endl;
    return;
#endif
#ifndef DEBUG
    fill(dist, dist + local_nnodes, numeric_limits<long long>::max());

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
                if (B[from].find(v) != B[from].end())
                    B[from].erase(v);
                B[to].insert(v);
                dist[v] = dv;
            }
        }
    }
#endif
}


int main(int argc, char * argv[]){
    string buf;
    int nnodes = 0, nedges = 0, local_nnodes = 0;
    ifstream grfile, ssfile;
    ofstream outfile;
    int src, trg;
    long long weight;
    long long delta = 1000;
    long long max_dist = 1000000;
    int np, rank, base_id;

    // disable sync with stdio
    ios_base::sync_with_stdio(false);

    if (argc != 4){
        printf("Usage: ./%s <graph file> <aux file> <output file>\n", argv[0]);
        exit(0);
    }

    MPI_Init(&argc, &argv); 
    MPI_Comm_size(MPI_COMM_WORLD, &np);  // get number of processes
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);  // get rank


    if (rank == 0){
        // proc 0 read graph file
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

        // broadcast nnodes and nedges
        MPI_Bcast(&nnodes, 1, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(&nedges, 1, MPI_INT, 0, MPI_COMM_WORLD);

        // compute local number of vertices
        local_nnodes = nnodes / np;
        if (rank == np - 1){
            local_nnodes += nnodes % np;
        }
        base_id = (nnodes / np) * rank;

        // convert graph to linear representation
        Edge * edges = (Edge *)malloc(nedges * sizeof(Edge));
        int * offsets_light = (int *)malloc(nnodes * sizeof(int));
        int * offsets_heavy = (int *)malloc(nnodes * sizeof(int));
        int * nedges_each_vertex = (int *)malloc(nnodes * sizeof(int));

        int local_nedges = offsets_light[1];
        Edge * local_edges = (Edge *)malloc(local_nedges * sizeof(Edge));
        int * local_offsets_light = (int *)malloc(local_nnodes * sizeof(int));
        int * local_offsets_heavy = (int *)malloc(local_nnodes * sizeof(int));
        int * local_nedges_each_vertex = (int *)malloc(local_nnodes * sizeof(int));
        long long * dist = (long long *)malloc(local_nnodes * sizeof(long long));

        Edge * ptr = edges;
        for (int v = 0; v < nnodes; v++){
            offsets_light[v] = ptr - edges;
            offsets_heavy[v] = ptr - edges + nlight[v];
            nedges_each_vertex[v] = graph[v].size();
            for (int j = 0; j < graph[v].size(); j++){
                ptr->trg = graph[v][j].first;
                ptr->weight = graph[v][j].second;
                ptr++;
            }
        }

        int * sendcounts_v = (int *)malloc(np * sizeof(int));
        int * displs_v = (int *)malloc(np * sizeof(int));
        int * sendcounts_e = (int *)malloc(np * sizeof(int));
        int * displs_e = (int *)malloc(np * sizeof(int));
        for (int pid = 0; pid < np; pid++){
            sendcounts_v[pid] = (pid == np - 1) ? ((nnodes / np) + (nnodes % np)) : (nnodes / np);
            displs_v[pid] = (nnodes / np) * pid;
        }
        for (int pid = 0; pid < np; pid++)
            displs_e[pid] = offsets_light[(nnodes / np) * pid] * sizeof(Edge);
        for (int pid = 0; pid < np; pid++)
            if (pid == np - 1)
                sendcounts_e[pid] = nedges * sizeof(Edge) - displs_e[pid];
            else
                sendcounts_e[pid] = (displs_e[pid + 1] - displs_e[pid]);

        MPI_Scatterv(nedges_each_vertex, sendcounts_v, displs_v, MPI_INT, local_nedges_each_vertex, local_nnodes, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Scatterv(offsets_light, sendcounts_v, displs_v, MPI_INT, local_offsets_light, local_nnodes, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Scatterv(offsets_heavy, sendcounts_v, displs_v, MPI_INT, local_offsets_heavy, local_nnodes, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Scatterv((char *)edges, sendcounts_e, displs_e, MPI_CHAR, (char *)local_edges, local_nedges * sizeof(Edge), MPI_CHAR, 0, MPI_COMM_WORLD);

        // proc 0 read the .ss file line by line and solve the corresponding sssp problem
        ssfile.open(argv[2]);
        outfile.open(argv[3], fstream::app);
        while (getline(ssfile, buf)){
            string op;
            if (buf[0] == 's'){
                long long checksum = 0;
                istringstream stream(buf);
                stream >> op >> src;
                MPI_Bcast(&src, 1, MPI_INT, 0, MPI_COMM_WORLD);
                parallel_delta_stepping(rank, np, local_edges, local_nnodes, base_id, dist, local_offsets_light,
                                        local_offsets_heavy, local_nedges_each_vertex, src - 1, dist, delta, max_dist);
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
                break;
            }
        }
        src = -1;
        MPI_Bcast(&src, 1, MPI_INT, 0, MPI_COMM_WORLD);
        ssfile.close();
        outfile.close();

        free(offsets_light);
        free(offsets_heavy);
        free(nedges_each_vertex);
        free(edges);
        free(sendcounts_v);
        free(sendcounts_e);
        free(displs_v);
        free(displs_e);
        free(local_offsets_light);
        free(local_offsets_heavy);
        free(local_nedges_each_vertex);
        free(local_edges);
        free(dist);
    }
    else {
        // broadcast nnodes and nedges
        MPI_Bcast(&nnodes, 1, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(&nedges, 1, MPI_INT, 0, MPI_COMM_WORLD);

        // compute local number of vertices
        local_nnodes = nnodes / np;
        if (rank == np - 1){
            local_nnodes += nnodes % np;
        }
        base_id = (nnodes / np) * rank;

        int * local_offsets_light = (int *)malloc(local_nnodes * sizeof(int));
        int * local_offsets_heavy = (int *)malloc(local_nnodes * sizeof(int));
        int * local_nedges_each_vertex = (int *)malloc(local_nnodes * sizeof(int));
        MPI_Scatterv(NULL, NULL, NULL, MPI_INT, local_nedges_each_vertex, local_nnodes, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Scatterv(NULL, NULL, NULL, MPI_INT, local_offsets_light, local_nnodes, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Scatterv(NULL, NULL, NULL, MPI_INT, local_offsets_heavy, local_nnodes, MPI_INT, 0, MPI_COMM_WORLD);
        int local_nedges = 0 ;
        for (int v = 0; v < local_nnodes; v++)
            local_nedges += local_nedges_each_vertex[v];
        Edge * local_edges = (Edge *)malloc(local_nedges * sizeof(Edge));
        long long * dist = (long long *)malloc(local_nnodes * sizeof(long long));

        MPI_Scatterv(NULL, NULL, NULL, MPI_CHAR, (char *)local_edges, local_nedges * sizeof(Edge), MPI_CHAR, 0, MPI_COMM_WORLD);

        cout << rank << ' ' << local_nnodes << endl;
        for (int v = 0; v < local_nedges_each_vertex[0]; v++)
            cout << local_edges[v].trg << ' ' << local_edges[v].weight << ' ' << endl;

        while (1){
            MPI_Bcast(&src, 1, MPI_INT, 0, MPI_COMM_WORLD);
            if (src == -1)
                break;
            parallel_delta_stepping(rank, np, local_edges, local_nnodes, base_id, dist, local_offsets_light,
                                    local_offsets_heavy, local_nedges_each_vertex, src - 1, dist, delta, max_dist);
        }
        free(local_offsets_light);
        free(local_offsets_heavy);
        free(local_nedges_each_vertex);
        free(local_edges);
        free(dist);
    }

    MPI_Finalize();

    return 0;
}
