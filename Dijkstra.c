#include "Dijkstra.h"
#include "Graph.h"
//#include "Graph.c"


#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

PredNode *insertPre(PredNode *current, int c);
int getlowestdistance(int *dist, int nV, int *visited);
void freepred(PredNode *i);

//Finds  all  shortest  paths  from  a given source vertex to all other
//vertices as discussed in the lectures. --timecomplexity O(nV^2)
ShortestPaths dijkstra(Graph g, Vertex src) {
    assert(g != NULL);
	ShortestPaths new;
	int nV = GraphNumVertices(g);
    new.numNodes = nV;
    new.src = src;
    new.dist = malloc(sizeof(int)*nV);
	new.pred = calloc(nV, sizeof(PredNode));
	int visited[nV];
	for (int i = 0; i < nV; i++) {
		new.dist[i] = 10000000;
		visited[i] = -1;
	}	
	new.dist[src] = 0;
    int c = getlowestdistance(new.dist, nV, visited);
	while (visited[c] == -1) {
		visited[c] = 1;
		int n = 0;
		while (n < nV) {
            if (GraphIsAdjacent(g, c, n) == 1) {
				AdjList i = GraphOutIncident(g, c);
				while (i->v != n) {
					i = i->next;
				}
				int weight = i->weight;
                if (new.dist[c] + weight == new.dist[n]) {
                    new.pred[n] = insertPre(new.pred[n], c);
                }
				if (new.dist[c] + weight < new.dist[n]) {
					if (new.pred[n] != NULL) {
						freepred(new.pred[n]);
						new.pred[n] = NULL;
					}
					new.dist[n] = new.dist[c] + weight;
					new.pred[n] = insertPre(new.pred[n], c);
				}
            }
			n++;
		}
		c = getlowestdistance(new.dist, nV, visited);
	}
	for(int i = 0; i < nV; i++) {
		if (new.dist[i] == 10000000 || new.dist[i] == 1000000) {
			new.dist[i] = 0;
		}
	}
	return new;
}

//insert the prenode, return the head of current prenode list --timecomplexity O(nV)
PredNode *insertPre(PredNode *current, int c) {
	if (current == NULL) {
		PredNode* new = malloc(sizeof(PredNode));
		new->v = c;
		new->next = NULL;
		return new;
	}
	else {
		PredNode* i = current;
		while (i->next != NULL){
			i = i->next;
		}
		PredNode* new = malloc(sizeof(PredNode));
		new->v = c;
		new->next = NULL;
		i->next = new;
		return current;
	}
}

//find a lowest unvisited one in dist[] and return it. --timecomplexity O(nV)
int getlowestdistance(int *dist, int nV, int *visited) {
	int i = 0;
	int current = dist[i];
	while (i < nV) {
		if (visited[i] == -1) {
			current = dist[i];
			break;
		}
		i++;
	}	
	i = 0;
	while (i < nV) {
		if (current > dist[i] && visited[i] == -1) {
			current = dist[i];
		}
		i++;
	}
	i = 0;
	while (i < nV) {
		if (current == dist[i] && visited[i] == -1) {
			
			return i;
		}
		i++;
	}
	return i;
}

// void showShortestPaths(ShortestPaths sps) {
// }

// Frees all memory associated with the given ShortestPaths structure. --timecomplexity O(nV^2)
void freeShortestPaths(ShortestPaths sps) {
	int c = sps.numNodes;
	free(sps.dist);
	int i = 0;
	while (i < c) {
		PredNode* c = sps.pred[i];
		freepred(c);
		i++;
	}
	free(sps.pred);
}

//Frees all memory associated with the given PredNode list. --timecomplexity O(nV)
void freepred(PredNode *i) {
	if (i != NULL) {
		freepred(i->next);
		free(i);
	}
}
