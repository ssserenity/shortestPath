#ifndef __ROUTE_H__
#define __ROUTE_H__

#include <vector>
#include<iostream>
#include<set> 
#include<queue>
#include<string>
#include<cstring>
#include <stdlib.h>     /* atoi */
#include<stack>

using namespace std;


/*--------------------------------------------------------------*/
#define MAX_VERTEX_NUM 600
#define INFINITY_DISTANCE 9999
/*Define the data type used in this file*/ 
typedef unsigned short WeightType;      // Weight type  
typedef unsigned short DataType;       //  Data type
typedef unsigned short Vertex;          // Vertex type 
typedef unsigned short IDType;          // the  data type of ID of edge 

/*Define the structure of adjacement element */

typedef struct AdjVNode *PtrToAdjVNode;
struct AdjVNode {
	Vertex AdjV;
	IDType edgeID;
	WeightType weight;
	PtrToAdjVNode next;
};

/*Define the adjacement list in the form of array*/ 

typedef struct VNode {
	PtrToAdjVNode firstEdge;
}AdjList[MAX_VERTEX_NUM];


/*Define the structure of graph in the form of link list*/

typedef struct GNode *PtrToGNode;
struct GNode {
	int numVertex;
	int numEdge;
	AdjList Graph;
};
typedef PtrToGNode LGraph;

/*Define the structure of edge*/

typedef struct ENode *ptrToENode;
struct ENode {
	IDType edgeID;
	Vertex v1, v2;
	WeightType weight ;
};
typedef ptrToENode Edge;

/*Define the structure of record that stores the weight and the path*/
typedef struct RecordNode *ptrToRecordNode;
struct RecordNode{
	WeightType weightSum;
	vector<IDType> path;
};
typedef ptrToRecordNode Record;


/*---------------------------------------------------------------------------------------------------*/

/*This function the main Entrance to find the shortest path */ 
void search_route(char *graph[5000], int edge_num, char *condition);

/*Build all the edges using information stored in the topo*/
vector<Edge> BuildAllEdges(char* topo[5000], int edge_num);

/*Convert the string which contains the information of edge to integer array*/
int*  CharToEdge(char *edge);

/*Create an Edge*/
Edge CreateEdge(int* array);



/*Build the Graph in link list form,store the topo informaiton in this graph*/
LGraph BuildGraph(int numVertex, vector<Edge> all_edges);

/*Initialize a graph with numV vertex but without any edge*/

LGraph CreateGraph(int numV);

/*Insert an edge to the graph*/
void InsertEdge(LGraph gra, Edge e);

/*Convert the string which contains the information of demand to integer array*/
vector<Vertex> ConvertDemandSet(char *demand);

/*Count on vertexes, and return the all vertexes set*/
set<Vertex> FindVertexSet(vector<Edge>all_edges);

/*Build visit vertex set that demanded,return the visit vertex set*/
set<Vertex> BuildVisitVertexSet(vector<Vertex> demandArray);

/*Build vertex set P_,P_ is all vertex set subtrace p_source and p_dest */
set<Vertex> BuildP_Set(set<Vertex> allVertexSet,Vertex p_source,Vertex p_dest);

/*Build not necessary visit vertex set*/
set<Vertex> BuildUnVisitVertexSet(set<Vertex> P_, set<Vertex> visit_set);

/*Use DFS algorithm to check if the source vertex and destination vertex are connected*/
bool DFS(LGraph gra,Vertex p_source,Vertex p_dest);
/*The recursion function part of DFS algorithm*/
void DFSRecur(LGraph gra,Vertex p_dest,Vertex v,bool check[],bool isIn[]);

/*Use Dijkstra Algorithm to find the shortest path from source vertex to destination vertex*/
Record Dijkstra(LGraph gra, Vertex p_source,Vertex p_dest);
/*Find the minum distance  vertex in the curent record in the Dijkstra function*/
Vertex FindMinDisVertex(set<Vertex>& uncheckedSet, WeightType disToSource[]);

/*Dynamic programmng recursion function,if find a feasible solution,add it to all_records*/
void DynamicRecursion(LGraph gra,Vertex v_prev,Vertex p_source,Vertex p_dest,set<Vertex> A,set<Vertex>B,Record pre_record,\
						vector<Record> & all_records,set<Vertex> &visit_set,set<Vertex> &un_visit_set);
						
/*Compare and find the shortes path in the records*/						
vector<IDType> GetShortestPath(vector<Record> &records);

/*Check if the final result is right,if right,return true*/
bool IsResultCorrect(vector<Edge> edges, vector<IDType> result,Vertex p_source,Vertex p_dest);
/*----------------------------------------------------------------------------------------------------------*/
#endif
