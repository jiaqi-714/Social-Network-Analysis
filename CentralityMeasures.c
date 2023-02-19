#include "Dijkstra.h"
#include "Graph.h"
#include "CentralityMeasures.h"


#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
void findpath(PredNode **pred, PredNode *curr, double *i, int src);
void findpathno(PredNode **pred, PredNode *curr, double *i, double *t, int v, int src);


/**
 * Finds the closeness centrality for each vertex in the given graph and
 * returns the results in a NodeValues structure. --timecomplexity O(nV^2)
 */ 
NodeValues closenessCentrality(Graph g) {
    int nV = GraphNumVertices(g);
    NodeValues new;
    new.numNodes = nV;
    new.values = calloc(nV, sizeof(double));
    
    for(int i = 0; i < nV; i++){
        double n = 0;
        double total = 0;
        
        ShortestPaths curr_path = dijkstra(g,i);

        for(int c = 0; c < nV; c++) {
            if (c != i) {
                if (curr_path.pred[c] != NULL) {
                    total = total + curr_path.dist[c];
                    n++;
                }
            }
        }   
        if (total != 0) {
            new.values[i] = ((n)/(nV-1))*((n)/total);
        }
        else {
            new.values[i] = 0;
        }
        freeShortestPaths(curr_path);
    }
    return new;
}


/**
 * Finds  the  betweenness centrality for each vertex in the given graph
 * and returns the results in a NodeValues structure.--timecomplexity O(nV^3)
 */

NodeValues betweennessCentrality(Graph g) {
    double nV = GraphNumVertices(g);
    NodeValues new;
    new.numNodes = nV;
    new.values = calloc(nV, sizeof(double));
    
    for(int i = 0; i < nV; i++) { 
        double ratio = 0;
        //printf("------------%d----------\n", i);
        for(int src = 0; src < nV; src++) {
            ShortestPaths curr_path = dijkstra(g,src);
            for(int end = 0; end < nV; end++) {
                if (i != src && i != end && src != end) {
                    double total = 0;
                    double *tp = &total;
                    double no = 0;
                    double *np = &no;   
                    PredNode *current = curr_path.pred[end]; 
                    findpathno(curr_path.pred, current, np, tp, i, src);
                    if (total != 0) {
                        ratio = ratio + ((no) / total);
                    }
                }   
            }
            freeShortestPaths(curr_path);
        }  
        //printf("%lf\n" , ratio);
        new.values[i] = ratio;
        
    }
    
    return new;
}



// find how many shortest paths from one node to a src node, t = number of shortest paths from one node to a src node.
// i = the number of shortest paths will lose without certain node v --timecomplexity O(nV^2)
void findpathno(PredNode **pred, PredNode *curr, double *i, double *t, int v, int src) {
    if (curr != NULL) {
        findpathno(pred, curr->next, i, t, v, src);
        findpathno(pred, pred[curr->v], i, t, v, src);
        if (curr->v == v) {
            double counter = 0;
            double *cp = &counter;
            findpath(pred, pred[curr->v], cp, src);
            *i = *i + counter;
        }
        if (curr->v == src) {
            *t = *t + 1;
        }
    }
    else if (curr == NULL) {
        return;
    }
    
}

//find how many shortest paths from one node to "src" node will lose without certain node.
//i = the number of shortest paths will lose without certain node v --timecomplexity O(nV^2)
void findpath(PredNode **pred, PredNode *curr, double *i, int src) {
    
    if (curr != NULL) {
        findpath(pred, curr->next, i, src);
        findpath(pred, pred[curr->v], i, src);
        if (curr->v == src) {
            *i = *i + 1;
        }
    }

}


// normalise the betweennessCentrality, --timecomplexity O(nV)
NodeValues betweennessCentralityNormalised(Graph g) {
    NodeValues new = betweennessCentrality(g);
    int nV = new.numNodes;
    if (nV > 2) {
        for (int i = 0; i < nV; i++) {
            double total = new.values[i];
            new.values[i] = total/(((nV-2)*(nV-1)));
        }
    }


    return new;

}

// print all NodeValues --timecomplexity O(nV)
void showNodeValues(NodeValues nvs) {
    for (int i = 0;i < nvs.numNodes;i++) {
        printf("%d: <%lf>\n", i, nvs.values[i]);
    }
}
/**
 * Frees all memory associated with the given NodeValues structure. --timecomplexity O(nV)
 */
void freeNodeValues(NodeValues nvs){
    free(nvs.values);
}