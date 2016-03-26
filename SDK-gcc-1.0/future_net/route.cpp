#include <vector>
#include<iostream>
#include<set> 
#include<queue>
#include<string>
#include<cstring>
#include <stdlib.h>     /* atoi */
#include<stack>

#include "route.h"
#include "lib_record.h"
#include <stdio.h>


/*---------------------------------------------------------------------------------------------------*/
//你要完成的功能总入口
void search_route(char *topo[5000], int edge_num, char *demand)
{
	/*Create all the edges*/
	vector<Edge> all_edges = BuildAllEdges(topo,edge_num);
	set<Vertex> allVertexSet = FindVertexSet(all_edges);
	int numVertex = allVertexSet.size();
	/*Build the graph*/
	LGraph gra = BuildGraph(numVertex,all_edges);
	/*Assign Id to the source vertex and destination vertex, and build the must visited vertex set S 
	and not necessary vistied vertex set  S_*/		
	vector<Vertex> demandArray = ConvertDemandSet(demand);
	Vertex p_source = demandArray[0];
	Vertex p_dest = demandArray[1];
	/* P_ is all vertex set subtrace p_source and p_dest */
	set<Vertex> P_ = BuildP_Set(allVertexSet,p_source,p_dest);
	/* A is to record already visited vertex in visit_set
	   B is to record already visited vertex in un_visit_set  */
	set<Vertex> visit_set,un_visit_set,A,B;
	visit_set = BuildVisitVertexSet(demandArray);
	un_visit_set = BuildUnVisitVertexSet(P_,visit_set);
	/*DFS to see if p_source and p_dest are connected or not, if not then  output no solution */
	bool connect = DFS(gra,p_source,p_dest);
	if (!connect){
		cout<<"No solution at all!!!"<<endl;
		return;
	}   
	else
	    cout<<"p_source and p_dest are connected "<<endl;    
	vector<Record> all_records;
	/*if the visit_set is empty,then use Dijkstra algorithm to find the shortest path*/
	if(visit_set.empty()){
		Record a_record = Dijkstra(gra,p_source,p_dest);
		all_records.push_back(a_record);
	}
	/* Use Dynamic programming method to find all paths that visit the specified vertex*/
	else{
		PtrToAdjVNode thisNode = gra->Graph[p_source].firstEdge;
		while (thisNode!= NULL) {
			// deep copy of A and B;
			set<Vertex> A1 = A;
			set<Vertex> B1 = B;
			Record a_record = new RecordNode;
			Vertex v = thisNode->AdjV;
			if(v == p_dest){
			}
			else{
				if(visit_set.count(v)){
					A1.insert(v);
				}
				else if(un_visit_set.count(v)){
					B1.insert(v);
				}
				a_record->weightSum = thisNode->weight;
				a_record->path.push_back(thisNode->edgeID);
				DynamicRecursion(gra,v,p_source,p_dest,A1,B1,a_record,all_records,visit_set,un_visit_set);
			}
			thisNode = thisNode->next;
		}
	}
	// because of some character encode problems, I write this comment in english
	// firstly, you need to initialize a vector to store the result like this:
	vector<IDType> result;
	
	/*If all_records is not empty*/
	if (all_records.size()){
		// then you use your algorithm to find a shortest route and put it into the "result". 	
		result = GetShortestPath(all_records);
		/*Check if the result is right*/
		if(!IsResultCorrect(all_edges, result,p_source,p_dest)){
			cout<<"The result is wrong!"<<endl;
		}else{
			cout<<"The result is right"<<endl;
		}
		// then print the results in the result.csv	
		for(std::vector<IDType>::iterator it = result.begin(); it != result.end(); ++it){
			record_result(*it);
			cout<<*it<<"|";
		}
	}else{
		cout<<"No solution at all!!!"<<endl;
	}
}
/*--------------------------------------------------------------------------------------------------------*/

/*Build all the edges using information stored in the topo*/
vector<Edge> BuildAllEdges(char *topo[5000], int edge_num)
{
	vector<Edge> all_edges ;
	for(int i=0;i<edge_num;i++){
		//cout<<topo[i]<<endl;
		int *array;
		array = CharToEdge(topo[i]);
		Edge a_edge = CreateEdge(array);
		all_edges.push_back(a_edge);
	}	
	return all_edges;
}


/*Convert the string which contains the information of edge to integer array*/
int* CharToEdge(char *edge)
{
	int *array = new int[4];
	string strArray[4] = {"","","",""};
	unsigned short count = 0;
	for(int i =0;i <strlen(edge) ;i++){
		if(edge[i] == ',')   {// the data was splited by','       
			array[count] = atoi(strArray[count].c_str());	//convert string to int
			count++;
		}else{
			strArray[count].push_back(edge[i]);
		}
		
	}
    array[count] = atoi(strArray[count].c_str());
    return array;
}

/*Convert the string which contains the information of demand to integer array*/
vector<Vertex> ConvertDemandSet(char *demand)
{
	vector<Vertex> array;
	string str = "";
	for(int i =0;i <strlen(demand) ;i++){
		if(demand[i] == ',' || demand[i] == '|'){   // the data was splited by','  or'|'     
		    int value = atoi(str.c_str());	//convert string to int
			array.push_back(value);
			str = "";
		}else{
			str.push_back(demand[i]);
		}
		
	}
    array.push_back(atoi(str.c_str())) ;
    return array;
}


/*Create an Edge*/
Edge CreateEdge(int* array)
{
		Edge a_edge = new ENode;
		a_edge->edgeID = *array;
		a_edge->v1 =  *(array+1);
		a_edge->v2 =  *(array+2);
		a_edge->weight =  *(array+3);
		delete[] array;
		return a_edge;
}
/*Count on vertexes, and return the all vertexes set*/
set<Vertex> FindVertexSet(vector<Edge>all_edges)
{
	set<Vertex> vertexSet;
	for(int i=0;i<all_edges.size();i++){
		vertexSet.insert(all_edges[i]->v1);
		vertexSet.insert(all_edges[i]->v2);
	}
	return vertexSet;
}

/*Build the Graph in link list form,store the topo informaiton in this graph*/
LGraph BuildGraph(int numV, vector<Edge> all_edges)
{
	LGraph gra;
	Edge e;
	Vertex v;
	gra = CreateGraph(numV);
	gra->numEdge = all_edges.size();
	if (gra->numEdge != 0) {
		for (int i = 0; i < gra->numEdge; i++) {
			e = all_edges[i];
			InsertEdge(gra, e);
		}
	}
	return gra;
}
/*Initialize a graph with numV vertex but without any edge*/
/*初始化一个有numV个顶点的但是没有边的图*/
LGraph CreateGraph(int numV)
{
	Vertex v, w;
	LGraph gra;
	gra = new GNode;
	gra->numVertex = numV;
	gra->numEdge = 0;
	for (int i = 0; i < numV; i++) {
		gra->Graph[i].firstEdge = NULL;
	}
	return gra;
}
/*Insert an edge to the graph*/
void InsertEdge(LGraph gra, Edge e)
{
	PtrToAdjVNode newNode;
	/*插入边<v1,v2>*/
	/*为v2建立邻接点*/
	newNode = new AdjVNode;
	newNode->AdjV = e->v2;
	newNode->weight = e->weight;
	newNode->edgeID = e->edgeID;
	/*将v2插入到v1的后面，并将v1后面的元素接到v2后面*/
	newNode->next = gra->Graph[e->v1].firstEdge;
	gra->Graph[e->v1].firstEdge = newNode;
}
/*------------------------------------------------------------------------------*/
/*Build visit vertex set that demanded,return the visit vertex set*/
set<Vertex> BuildVisitVertexSet(vector<Vertex> demandArray)
{
    set<Vertex> visit_set;
	for(int j=2;j<demandArray.size();j++){
		visit_set.insert(demandArray[j]);
	}
	return visit_set;
 } 

/*Build vertex set P_,P_ is all vertex set subtrace p_source and p_dest */
set<Vertex> BuildP_Set(set<Vertex> allVertexSet,Vertex p_source,Vertex p_dest)
{
	set<Vertex> P_ = allVertexSet;   // P_ is all vertex set subtrace p_source and p_dest;
	P_.erase(p_source);
	P_ .erase(p_dest); 
	return P_;
}

/*Build not necessary visit vertex set*/
set<Vertex> BuildUnVisitVertexSet(set<Vertex> P_, set<Vertex> visit_set)
{
	set<Vertex> un_visit_set = P_;
	for(set<Vertex>::iterator it = visit_set.begin(); it != visit_set.end(); ++it){
		Vertex v = *it;
		un_visit_set.erase(v);
	}
	return un_visit_set;
}

/*-------------------------------------------------------------------------------------*/ 

/*Use DFS algorithm to check if the source vertex and destination vertex are connected*/
bool DFS(LGraph graph,Vertex p_source,Vertex p_dest)   //从顶点p_source出发深度优先	搜索图
{
	
	bool *check = new bool[graph->numVertex]; /*标识符，表示该点是否已经检查*/
	bool *isIn = new bool[graph->numVertex];  /*标识符，表示该顶点是否与p_source连通 */
	for (int i = 0; i < graph->numVertex; i++) {  //初始化标识符
		check[i] = false;
		isIn[i] = false;
	}
	
    DFSRecur(graph,p_dest,p_source,check,isIn);
    if (isIn[p_dest]){  // the source vertex and destination vertex are connected
    	return true;
	} else{
	    return false;	
	}
}
/*The recursion function part of DFS algorithm*/
/*递归DFS*/
void DFSRecur(LGraph graph,Vertex p_dest,Vertex v,bool check[],bool isIn[])
{
	if (check[v]||graph->Graph[v].firstEdge == NULL) {
		return;
	}
	check[v] = true;
	PtrToAdjVNode currentNode = graph->Graph[v].firstEdge;
	while (currentNode!= NULL) {
		int i = currentNode->AdjV;
		isIn[i] = true;
		if (i == p_dest){
			return;
		}
		DFSRecur(graph, p_dest,i,check,isIn);
		currentNode = currentNode->next;
	}
}

/*-----------------------------------------------------------------------------*/

/*Use Dijkstra Algorithm to find the shortest path from source vertex to destination vertex*/
Record Dijkstra(LGraph gra, Vertex p_source,Vertex p_dest)
{
	int numV = gra->numVertex;
	set<Vertex> uncheckedSet;
	bool *checked = new bool[numV];
	WeightType *disToSource = new 	WeightType[numV];
	Vertex *prevV = new Vertex[numV];
	IDType edgePath[MAX_VERTEX_NUM][MAX_VERTEX_NUM];
	/*Initialize the distance*/
	for(int i = 0;i < numV; i++ ){
		uncheckedSet.insert(i);
		checked[i] = false;
		disToSource[i] = INFINITY_DISTANCE;// set the distance to infinity
	}
	disToSource[p_source] = 0; //distance from source to source
	while(!uncheckedSet.empty()){
		Vertex v = FindMinDisVertex(uncheckedSet, disToSource);
		uncheckedSet.erase(uncheckedSet.find(v));
		/* check the adjacement vertexes of v*/
		PtrToAdjVNode adVertex = gra->Graph[v].firstEdge;
		while(adVertex!= NULL){
			if(checked[adVertex->AdjV]){
			}
			else{
				WeightType altDis = disToSource[v]+adVertex->weight;
				if(altDis<disToSource[adVertex->AdjV]){
					disToSource[adVertex->AdjV] = altDis;
					 prevV[adVertex->AdjV] = v;
					 edgePath[v][adVertex->AdjV] = adVertex->edgeID; 
				}
			}
			adVertex = adVertex->next;
	    }
		checked[v] = true;
		if (checked[p_dest]) break;  // The destination have been checked.	
	}
	Record a_record = new RecordNode;
	a_record->weightSum = disToSource[p_dest];
	/*Retrieve the edge path information*/
	Vertex v_current = p_dest;
	Vertex w;
	stack<IDType> edges;
	while(v_current != p_source){
	    w = prevV[v_current];
	    edges.push(edgePath[w][v_current]);
	   v_current = w;
	}
	while(!edges.empty()){
		a_record->path.push_back(edges.top());
		edges.pop();
	}
	delete[] checked;
	delete[] disToSource;
	delete[] prevV;
    return a_record;	
}
/*Find the minum distance  vertex in the curent record in the Dijkstra function*/
Vertex FindMinDisVertex(set<Vertex>& uncheckedSet, WeightType disToSource[])
{
	Vertex index = MAX_VERTEX_NUM;
	WeightType dis = INFINITY_DISTANCE;
	set<Vertex>::iterator it = uncheckedSet.begin();
	while (it != uncheckedSet.end()) {
		if (disToSource[*it]<dis) {
			dis = disToSource[*it];
			index = *it;
		}
		it++;
	}
	return index;
}
/*Dynamic programmng recursion function*/
void DynamicRecursion(LGraph gra,Vertex v_prev,Vertex p_source,Vertex p_dest,set<Vertex> A,set<Vertex>B,Record pre_record,\
						vector<Record> & all_records,set<Vertex> &visit_set,set<Vertex> &un_visit_set)
{
	PtrToAdjVNode thisNode = gra->Graph[v_prev].firstEdge;
	if(!thisNode){
		delete pre_record;
	}
	while (thisNode!= NULL) {
		// deep copy of A,B and pre_record
		set<Vertex> A1 = A;
		set<Vertex> B1 = B;
		Record a_record = new RecordNode;
		*a_record = *pre_record;
		Vertex v = thisNode->AdjV;
		/*Base condition, terminate the recursion*/
		/*If the Vertex is already in A1 or B1,or it is the source vertex*/
		if(A1.count(v) ||B1.count(v) || (v == p_source)){
		}
		/*If it reaches the destination vertex and visited all the demand vertexes*/ 
		else if(v == p_dest){
			if(A1.size() == visit_set.size()){     // reaches the destination vertex
				a_record->weightSum += thisNode->weight;
				a_record->path.push_back(thisNode->edgeID);
				all_records.push_back(a_record);
			}
		}
		else{
			if(visit_set.count(v)){
				A1.insert(v);
				
			}
			else if(un_visit_set.count(v)){
				B1.insert(v);
			}
			a_record->weightSum += thisNode->weight;
			a_record->path.push_back(thisNode->edgeID);
			DynamicRecursion(gra,v,p_source, p_dest,A1,B1,a_record,all_records,visit_set,un_visit_set);
		}		
		thisNode = thisNode->next;		
	}	
}

/*Compare and find the shortes path in the records*/
vector<IDType> GetShortestPath(vector<Record> &records)
{
	WeightType minWeightSum = INFINITY_DISTANCE; 
	WeightType currentWeight;
	vector<IDType> shortestPath;
	int minIndex = -1;
	for(int i = 0;i<records.size();i++){
		currentWeight = records[i]->weightSum;
		if(currentWeight < minWeightSum){
			minWeightSum = currentWeight; 
			minIndex = i;
		}	
	}
	shortestPath = records[minIndex]->path;
	return shortestPath;
}
/*Check if the final result is right*/
bool IsResultCorrect(vector<Edge> edges, vector<IDType> result,Vertex p_source,Vertex p_dest)
{
	int num = result.size();
	if(edges[result[0]]->v1 != p_source || num > edges.size() ||\
	    edges[result[num-1]]->v2 != p_dest){
	    	return false;
		}
	for(int i = 0;i<num-1;i++){
		if(edges[result[i]]->v2 != edges[result[i+1]]->v1 ){
			return false;
		}
	}
	return true;
}












