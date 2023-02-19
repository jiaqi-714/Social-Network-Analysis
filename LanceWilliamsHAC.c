// Lance-Williams Algorithm for Hierarchical Agglomerative Clustering
// COMP2521 Assignment 2

#include <stdio.h>
#include <stdlib.h>

#include "Graph.h"
#include "LanceWilliamsHAC.h"

/**
 * Generates  a Dendrogram using the Lance-Williams algorithm (discussed
 * in the spec) for the given graph  g  and  the  specified  method  for
 * agglomerative  clustering. The method can be either SINGLE_LINKAGE or
 * COMPLETE_LINKAGE (you only need to implement these two methods).
 * 
 * The function returns a 'Dendrogram' structure.
 */

Dendrogram merge(Graph g, int nV, double **dist, Dendrogram curr[], double* **newdist, int method);
void findfirstmap(Graph g, double **dist);
void Singlelinkreassigndist(Graph g, double **dist, int nV, Dendrogram current, double* **newdist, int currentnode[], int numberofnode);
void combine(Dendrogram curr[], Dendrogram current, Dendrogram new, int *i, int currentnode[]);
void printall(double* **newdist, int nV);
void Completelinkreassigndist(Graph g, double **dist, int nV, Dendrogram current, double* **newdist, int currentnode[], int numberofnode);


//time complexity O(nV^3)
Dendrogram LanceWilliamsHAC(Graph g, int method) {
    int nV = GraphNumVertices(g);
    double **dist;
    dist = calloc(nV, sizeof(double));
    for (int i = 0; i < nV; i++) {
        dist[i] = calloc(nV, sizeof(double));
    }
    for (int source = 0; source < nV; source++) {
        for(int end = 0; end < nV; end++) {
            dist[source][end] = 10;
        }
    }

    findfirstmap(g, dist);
    double* **newdist;
    newdist = calloc(nV, sizeof(double*));
    for (int source = 0; source < nV; source++) {
        newdist[source] = calloc(nV, sizeof(double*));
    }
    for (int source = 0; source < nV; source++) {
        for(int end = 0; end < nV; end++) {
            newdist[source][end] = &dist[source][end];
        }
    }
    //printall(newdist, nV);
    Dendrogram curr[nV];
    for(int i = 0; i < nV; i++) {
        curr[i] = malloc(sizeof(Dendrogram));
        curr[i]->vertex = i;
        curr[i]->left = NULL;
        curr[i]->right = NULL;
    }
    Dendrogram current = NULL;
    for(int i = 0; i < nV - 1; i++)  {
        current = merge(g, nV, dist, curr, newdist, method);
    }
    
    for (int i = 0; i < nV; i++) {
        free(dist[i]);
        free(newdist[i]);
    }
    free(dist);
    free(newdist);
    return current;

}

//merge the graph, --time complexity O(nV^3)
Dendrogram merge(Graph g, int nV, double **dist, Dendrogram curr[], double* **newdist, int method) {
    double lowest = 10;
    int src;
    int ed;
    for (int source = 0; source < nV; source++) {
        for(int end = 0; end < nV; end++) {
            if (source != end) {
                if (lowest > *newdist[source][end]) {
                    lowest = *newdist[source][end];
                    src = source;
                    ed = end;
                }
            }
        }
    }
    //printf("-----merge %d %d - %.2f-----\n", src, ed, *newdist[src][ed]);
    int currentnode[nV];
    int numberofnode = 0;
    Dendrogram current = malloc(sizeof(Dendrogram));
    current->left = curr[src];
    current->right = curr[ed];
    combine(curr, current, current, &numberofnode, currentnode);
    if (method == 1) {
        Singlelinkreassigndist(g, dist, nV, current, newdist, currentnode, numberofnode);
    }
    else {
        Completelinkreassigndist(g, dist, nV, current, newdist, currentnode, numberofnode);
    }
    return current;
}

//replace the pointer of Dendrogram list by the cluster they belong and calcalate the number of nodes in "current" cluster
//for iteration in Singlelinkreassigndist and Completelinkreassigndist
void combine(Dendrogram curr[], Dendrogram current, Dendrogram new, int *i, int currentnode[]) {
    if (current == NULL) {
        return;
    }
    if (current->left == NULL && current->right == NULL) {
        curr[current->vertex] = new;
        currentnode[*i] = current->vertex;
        *i = *i + 1;
    }
    combine(curr, current->left, new, i, currentnode);
    combine(curr, current->right, new, i, currentnode);
}

//update the "newdist" graph according to "dist" by Singlelink --time complexity O(nV^3)
void Singlelinkreassigndist(Graph g, double **dist, int nV, Dendrogram current, double* **newdist, int currentnode[], int numberofnode) {
    for (int i = 0; i < numberofnode; i++) {
        for (int c = 0; c < numberofnode; c++) {
            dist[currentnode[i]][currentnode[c]] = 10;
            dist[currentnode[c]][currentnode[i]] = 10;
        }
    }
    for (int i = 0; i < numberofnode; i++) {
        for (int source = 0; source < nV; source++) {
            if(GraphIsAdjacent(g, source, currentnode[i])) {
                double min = dist[source][currentnode[i]];
                int low = currentnode[i];
                for (int c = 0; c < numberofnode; c++) {
                    if (GraphIsAdjacent(g, source, currentnode[c])) {
                        if (min > dist[source][currentnode[c]]) {
                            min = dist[source][currentnode[c]];
                            low = currentnode[c];
                        }
                    }
                }
                
                for (int c = 0; c < numberofnode; c++) {
                    newdist[source][currentnode[c]] = &dist[source][low];
                }
                
            }
            if(GraphIsAdjacent(g, currentnode[i], source)) {
                double min = dist[currentnode[i]][source];
                int low = currentnode[i];
                for (int c = 0; c < numberofnode; c++) {
                    if (GraphIsAdjacent(g, currentnode[c], source)) {
                        if (min > dist[currentnode[c]][source]) {
                            min = dist[currentnode[c]][source];
                            low = c;
                        }
                    }
                }
                //printf("low in %d %d, %dnodes\n",low, source, numberofnode);
                for (int c = 0; c < numberofnode; c++) {
                    //printf("change %d %d to %d %d\n", currentnode[c], source, low, source);
                    newdist[currentnode[c]][source] = &dist[low][source];
                }   
            }
            
        }
    }

}

//update the "newdist" graph according to "dist" by Completelink --time complexity O(nV^3)
void Completelinkreassigndist(Graph g, double **dist, int nV, Dendrogram current, double* **newdist, int currentnode[], int numberofnode) {
    for (int i = 0; i < numberofnode; i++) {
        for (int c = 0; c < numberofnode; c++) {
            dist[currentnode[i]][currentnode[c]] = 10;
            dist[currentnode[c]][currentnode[i]] = 10;
        }
    }
    for (int i = 0; i < numberofnode; i++) {
        for (int source = 0; source < nV; source++) {
            if(GraphIsAdjacent(g, source, currentnode[i])) {
                double lar = dist[source][currentnode[i]];
                int high = currentnode[i];
                for (int c = 0; c < numberofnode; c++) {
                    if (GraphIsAdjacent(g, source, currentnode[c])) {
                        if (lar < dist[source][currentnode[c]]) {
                            lar = dist[source][currentnode[c]];
                            high = currentnode[c];
                        }
                    }
                }
                for (int c = 0; c < numberofnode; c++) {
                    newdist[source][currentnode[c]] = &dist[source][high];
                }
            }
            if(GraphIsAdjacent(g, currentnode[i], source)) {
                double lar = dist[currentnode[i]][source];
                int high = currentnode[i];
                for (int c = 0; c < numberofnode; c++) {
                    if (GraphIsAdjacent(g, currentnode[c], source)) {
                        if (lar < dist[currentnode[c]][source]) {
                            lar = dist[currentnode[c]][source];
                            high = c;
                        }
                    }
                }
                //printf("high in %d %d, %dnodes\n",high, source, numberofnode);
                for (int c = 0; c < numberofnode; c++) {
                    newdist[currentnode[c]][source] = &dist[high][source];
                    //printf("change %d %d to %d %d\n", currentnode[c], source, high, source);
                }
            }
            
        }
    }
}

//print all things in dist --time complexity O(nV^2)
void printall(double** *dist, int nV) {
    for (int source = 0; source < nV; source++) {
        for(int end = 0; end < nV; end++) {
            printf("%d %d - %.2f ", source, end, *dist[source][end]);
        }
        printf("\n");
    }
    printf("\n");
}

//get the original dist map, --time complexity O(nV^2)
void findfirstmap(Graph g, double **dist) {
    int nV = GraphNumVertices(g);
    for (int src = 0; src < nV; src++) {
        AdjList curr = GraphOutIncident(g, src);
        if (curr != NULL) {
            while (curr->next != NULL){
                dist[src][curr->v] = (double) 1/curr->weight;
                curr= curr->next;
            }
            dist[src][curr->v] = (double) 1/curr->weight;
        }
    }
}

/**
 * Frees all memory associated with the given Dendrogram structure.
 */
void freeDendrogram(Dendrogram d) {
    if (d == NULL) {
		return;
	}
	freeDendrogram(d->left);
	freeDendrogram(d->right);
	free(d);
}
