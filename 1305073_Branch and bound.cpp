#include<bits/stdc++.h>
#include<stdio.h>
#include<fstream>
#include<iostream>
using namespace std;

#define NULL_VALUE -999999
#define INF 999999
#define MAX 1000000
#define WHITE 1
#define GREY 2
#define BLACK 3

int x_i[MAX],y_i[MAX];
clock_t start,stop;

class Queue
{
    int queueInitSize ;
    int queueMaxSize;
    int * data;
    int length;
    int front;
    int rear;
public:
    Queue();
    ~Queue();
    void enqueue(int item); //insert item in the queue
    int dequeue(); //returns the item according to FIFO
    bool empty(); //return true if Queue is empty
};

Queue::Queue()
{
    queueInitSize = 2 ;
    queueMaxSize = queueInitSize;
    data = new int[queueMaxSize] ; //allocate initial memory
    length = 0 ;
    front = 0;
    rear = 0;
}


void Queue::enqueue(int item)
{
	if (length == queueMaxSize)
	{
		int * tempData ;
		//allocate new memory space for tempList
		queueMaxSize = 2 * queueMaxSize ;
		tempData = new int[queueMaxSize] ;
		int i, j;
		j = 0;
		for( i = rear; i < length ; i++ )
		{
			tempData[j++] = data[i] ; //copy items from rear
		}
		for( i = 0; i < rear ; i++ )
		{
			tempData[j++] = data[i] ; //copy items before rear
		}
		rear = 0 ;
		front = length ;
		delete[] data ; //free the memory allocated before
		data = tempData ; //make list to point to new memory
	}

	data[front] = item ; //store new item
	front = (front + 1) % queueMaxSize ;
	length++ ;
}


bool Queue::empty()
{
	if(length == 0) return true ;
	else return false ;
}


int Queue::dequeue()
{
	if(length == 0) return NULL_VALUE ;
	int item = data[rear] ;
	rear = (rear + 1) % queueMaxSize ;
	length-- ;
	return item ;
}


Queue::~Queue()
{
    if(data) delete[] data; //deallocate memory
    data = 0; //set to NULL
}

//****************Queue class ends here************************

//****************Dynamic ArrayList class based************************
class ArrayList
{
	double * list;
	int length ;
	int listMaxSize ;
	int listInitSize ;
public:
	ArrayList() ;
	~ArrayList() ;
	int searchItem(double item) ;
    void insertItem(double item) ;
	void removeItem(double item) ;
	void removeItemAt(int item);
	double getItem(int position) ;
	int getLength();
	bool empty();
	void printList();
} ;


ArrayList::ArrayList()
{
	listInitSize = 2 ;
	listMaxSize = listInitSize ;
	list = new double[listMaxSize] ;
	length = 0 ;
}

void ArrayList::insertItem(double newitem)
{
	double * tempList ;
	if (length == listMaxSize)
	{
		//allocate new memory space for tempList
		listMaxSize = 2 * listMaxSize ;
		tempList = new double[listMaxSize] ;
		int i;
        for( i = 0; i < length ; i++ )
        {
            tempList[i] = list[i] ; //copy all items from list to tempList
        }
        delete[] list ; //free the memory allocated before
        list = tempList ; //make list to point to new memory
	};

	list[length] = newitem ; //store new item
	length++ ;
}

int ArrayList::searchItem(double item)
{
	int i = 0;
	for (i = 0; i < length; i++)
	{
		if( list[i] == item ) return i;
	}
	return NULL_VALUE;
}

void ArrayList::removeItemAt(int position) //do not preserve order of items
{
	if ( position < 0 || position >= length ) return ; //nothing to remove
	list[position] = list[length-1] ;
	length-- ;
}


void ArrayList::removeItem(double item)
{
	int position;
	position = searchItem(item) ;
	if ( position == NULL_VALUE ) return ; //nothing to remove
	removeItemAt(position) ;
}


double ArrayList::getItem(int position)
{
	if(position < 0 || position >= length) return NULL_VALUE ;
	//cout << "Return : " << list[position] << endl;
	return list[position] ;
}

int ArrayList::getLength()
{
	return length ;
}

bool ArrayList::empty()
{
    if(length==0)return true;
    else return false;
}

void ArrayList::printList()
{
    int i;
    for(i=0;i<length;i++)
        printf("%d ", list[i]);
    printf("Current size: %d, current length: %d\n", listMaxSize, length);
}

ArrayList::~ArrayList()
{
    if(list) delete [] list;
    list = 0 ;
}

//******************ArrayList class ends here*************************

//******************Graph class starts here**************************
class Graph
{
	int nVertices, nEdges ;
	bool directed ;
	ArrayList  * adjList ;
	ArrayList *weight;
	int *color;
	int *parent;
	int *distance;
	bool *visited;
	//define other variables required for bfs such as color, parent, and dist
	//you must use pointers and dynamic allocation

public:
    double tsp_cost;
    vector<int> tsp_path;
	Graph(bool dir = false);
	~Graph();
	void setnVertices(int n);
	void addEdge(int u, int v);
	void removeEdge(int u, int v, int w);
	bool isEdge(int u, int v, int w);
    int getDegree(int u);
    void printAdjVertices(int u);
    bool hasCommonAdjacent(int u, int v);
    int getDist(int u, int v);
    void printGraph();
	void bfs(int source); //will run bfs in the graph
	void dfs(int source); //will run dfs in the graph
	double firstMin(int i);
	double secondMin(int i);
	void TSP();
	void TSPbnb(double bound, double curr_weight,int level,vector<int> &current_path);
};


Graph::Graph(bool dir)
{
	nVertices = 0 ;
	nEdges = 0 ;
	adjList = 0 ;
	weight = 0;
	directed = dir ; //set direction of the graph
	color = new int [nVertices];
	visited = new bool[nVertices];
	parent = new int [nVertices];
	distance = new int [nVertices];
	//define other variables to be initialized
}

void Graph::setnVertices(int n)
{
	this->nVertices = n ;
	this->tsp_cost = DBL_MAX;
	if(adjList!=0) delete[] adjList ; //delete previous list
	adjList = new ArrayList[nVertices] ;
	if(weight!=0) delete[] weight;
	weight = new ArrayList [nVertices];

	for(int i=0;i<nVertices;i++) visited[i] = false;
}

void Graph::addEdge(int u, int v)
{
    double c;
    if(u<0 || v<0 || u>=nVertices || v>=nVertices) return; //vertex out of range
    this->nEdges++ ;
	adjList[u].insertItem(v) ;
	//cout << "(" << x_i[u] << "," << y_i[u] <<") (" << x_i[v] << "," << y_i[v] << ")"<< endl;
	c = sqrt(pow(x_i[u]-x_i[v],2) + pow(y_i[u]-y_i[v],2));
	//cout << "c: " << c << endl;
	weight[u].insertItem(c);
	if(!directed){
        adjList[v].insertItem(u) ;
        weight[v].insertItem(c) ;
        //cout << weight[v].getItem(0) << " ";
	}
}

void Graph::removeEdge(int u, int v,int w)
{
    //write this function
    int pos;
    if(isEdge(u,v,w)==1)nEdges--;

    adjList[u].removeItem(v);
    pos=adjList[u].searchItem(v);
    weight[u].removeItemAt(pos);

    if(!directed){

        adjList[v].removeItem(u);
        weight[v].removeItemAt(pos);

    }
}

bool Graph::isEdge(int u, int v, int w)
{
    //returns true if (u,v) is an edge, otherwise should return false
    for(int i=0;i<adjList[u].getLength() ;i++){
        if(adjList[u].getItem(i)==v && weight[u].getItem(i)==w) return true;
    }
    return false;
}

int Graph::getDegree(int u)
{
    //returns the degree of vertex u
    return adjList[u].getLength();
}

void Graph::printAdjVertices(int u)
{
    //prints all adjacent vertices of a vertex u
    for(int i=0;i<adjList[u].getLength();i++){
        printf("%d ",adjList[u].getItem(i));
    }
    printf("\n");
}

bool Graph::hasCommonAdjacent(int u, int v)
{
    //returns true if vertices u and v have common adjacent vertices
    for(int i=0;i<adjList[i].getLength();i++){
        for(int j=0;j<adjList[j].getLength();j++){
            if(adjList[u].getItem(i)==adjList[v].getItem(j)) return true;
        }
    }
    return false;

}

void Graph::bfs(int source)
{
    //complete this function
    //initialize BFS variables
    for(int i=0; i<nVertices; i++)
    {
        color[i] = WHITE ;
        parent[i] = -1 ;
        distance[i] = INF ;
    }
    Queue q ;
    color[source] = GREY;
    distance[source] = 0 ;
    parent[source] = -1;
    q.enqueue(source) ;
    while( !q.empty() )
    {
        //complete this part
        int u = q.dequeue();
        color[u] = GREY;
        for(int i=0;i<adjList[u].getLength();i++){
            int v = adjList[u].getItem(i);
            if(color[v]==WHITE){
                color[v] = GREY;
                q.enqueue(v);
                distance[v] = distance[u]+1;
                parent[v] = u;
            }
        }
        color[u] = BLACK;
    }
}

int Graph::getDist(int u, int v)
{
    //returns the shortest path distance from u to v
    //must call bfs using u as the source vertex, then use distance array to find the distance
    if(u>=0 && u<nVertices && v>=0 && v<nVertices){
        bfs(u);
        return distance[v];
    }
    return INF ;
}

void Graph::printGraph()
{
    printf("\nNumber of vertices: %d, Number of edges: %d\n", nVertices, nEdges);
    for(int i=0;i<nVertices;i++)
    {
        printf("%d th vertice: ", i);
        for(int j=0; j<adjList[i].getLength();j++)
        {
            //printf(" %d: %0.1f ", adjList[i].getItem(j),weight[i].getItem(j));
            cout << adjList[i].getItem(j) << ": " << weight[i].getItem(j) << " ";
        }
        printf("\n");
    }
}

Graph::~Graph()
{
    //write your destructor here
    if(adjList!=0) delete[] adjList;
    adjList = 0;
    if(weight!=0) delete[] weight;
    weight = 0;
    if(color) delete[] color;
    color = 0;
    if(distance) delete[] distance;
    distance = 0;
    if(parent) delete[] parent;
    parent = 0;
}


double Graph::firstMin(int i)
{
    double min = DBL_MAX;

    for(int j=0;j<weight[i].getLength();j++){
        double c = weight[i].getItem(j);
        if(c<min) min = c;
    }
    return min;
}

double Graph::secondMin(int i)
{
    double currentMax = weight[i].getItem(0);
    double secondMax=0;
    for( int j = 0; j< weight[i].getLength(); i++) {
        if (weight[i].getItem(j) > secondMax) secondMax = weight[i].getItem(j);

        if (weight[i].getItem(j) > currentMax){
            secondMax=currentMax;
            currentMax=weight[i].getItem(j);
        }

        return secondMax;
    }
}
void Graph::TSPbnb(double bound, double curr_weight,int level,vector<int> &current_path)
{

    if(level==nVertices){
        //cout << "Here\n";
        int ix = adjList[current_path[level-1]].searchItem(0);
        double curr_res = curr_weight + weight[current_path[level-1]].getItem(ix);

        if (curr_res < tsp_cost){
                tsp_path = current_path;
                tsp_cost = curr_res;
                current_path.clear();
                //cout << "Cost: " << tsp_cost <<endl;
        }
        return;
    }


    for(int i=0;i<adjList[current_path[level-1]].getLength();i++){

        int idx = adjList[current_path[level-1]].getItem(i);
        //cout << i << "," << idx << endl;

        if (visited[idx] == false){

            double temp = bound;
            curr_weight += weight[current_path[level-1]].getItem(i);

            if (level==1) bound -= ((firstMin(current_path[level-1]) + firstMin(idx))/2);
            else bound -= ((secondMin(current_path[level-1]) + firstMin(idx))/2);

            if (bound + curr_weight < tsp_cost){
                current_path.push_back(idx);
                visited[idx] = true;
                //cout << "curr_weight : " << curr_weight << endl;
                TSPbnb(bound, curr_weight, level+1,current_path);
            }

            //cout << "Pruning\n";
            curr_weight -= weight[current_path[level-1]].getItem(i);
            bound = temp;

            memset(visited, false, sizeof(visited));
            for (int j=0; j<=level-1; j++)
                visited[current_path[j]] = true;
        }
    }
}

void Graph::TSP()
{
    vector<int> current_path;

    double bound = 0;
    memset(visited, 0, sizeof(visited));

    for(int i=0;i<nVertices;i++) bound += firstMin(i) + secondMin(i);
    bound /= 2;
    //cout << bound << endl;

    visited[0] = true;
    current_path.push_back(0);

    TSPbnb(bound,0,1,current_path);
}


//**********************Graph class ends here******************************


//******main function to test your code*************************

int main(void)
{
    int n,x,y;
    int it=0;
    Graph g;
    ifstream input;
    input.open("1305073_input_8.txt");
    //input.open("input.txt");

    start = clock();
    input >> n;
    g.setnVertices(n);

    while(input >> x >> y){
        x_i[it] = x;
        y_i[it] = y;
        it++;
    }

    for(int i=0;i<n;i++){
        for(int j=i+1;j<n;j++){
            g.addEdge(i,j);
        }
    }

    while(1)
    {
        printf("1. TSP      2.Print Graph       3. Exit.\n");

        int ch;
        scanf("%d",&ch);
        if(ch==1)
        {
            start = clock();
            //cout << start << endl;
            g.TSP();

            //cout << stop << endl;

            cout << "TSP Cost: " << g.tsp_cost << endl;
            //cout << g.tsp_path.size() << endl;
            cout << "TSP Path: ";

            for(int i=0;i<g.tsp_path.size();i++) cout << g.tsp_path[i] << "-> " ;
            cout << 0 << endl;
            stop = clock();

            double et = (stop - start);
            //cout << et << endl;
            //ofstream output;
            //output.open("output2.csv",std::ios_base::app | std::ios_base::out);
            //output << n  << " , " << et << endl;

        }
        else if(ch==2)
        {
            g.printGraph();
        }
        else if(ch==3)
        {
            break;
        }
    }

}


