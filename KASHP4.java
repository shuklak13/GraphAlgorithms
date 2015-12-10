// @author Karan

import java.util.*;
import java.io.*; 

class GraphAlgorithms{
	public static void main(String [] args) throws IOException {
		Scanner sc = new Scanner(new File("inputFile.txt")); 
            //first line of input file should be the name of the data
            //the rest of the lines represent the edges
            //each line should have three integers, separated by a comma:
                //first and second: start and end node
                //third: weight of edge
		
        String title = sc.nextLine();
        int n = sc.nextInt();       //num vertices
        Graph g = new Graph(n);
        
        while(sc.hasNextInt())
            g.addEdge(sc.nextInt(), sc.nextInt(), sc.nextInt());

        //all the algo's are stored here
        GraphAlgo a = new GraphAlgo(g, 0);
            
        System.out.println(title);

		if(g.containsNegativeEdge()){
                if(a.bellmanFord())
                    a.printShortestPathsBF();
                else System.out.println("The graph contains a negative edge cycle");
        }
        else{
            a.dijkstra();
            a.printShortestPathsBF();
            System.out.println();
            a.floydWarshall();
            a.modifiedFloydWarshall();
        }
	}
}

class GraphAlgo{
    Graph g;
    int[] distances;     //distance to each vertex
    int[] predecessors;   //previous vertex to each vertex
    PriorityQueue unsettled;    //vertices unvisited, but adjacent to visited nodes (Dijkstra)
    Node[] settled; //vertices already visited (Dijkstra)
        int numSettled=0;
    int start;  //starting vertex, equals 0

    GraphAlgo(Graph grph, int s){
        start = s;
        g = grph;
        distances = new int[g.n];
        predecessors = new int[g.n];

		//set distances of everything to infinity (except start node)
        for(int i=0; i<g.n; i++){
            if(i==start){
            	distances[i] = 0;   //the distance between the start node and itself is 0
            	predecessors[i] = start+1;	//the start node is its own predecessor
            }
            else            distances[i] = Integer.MAX_VALUE;  //for other nodes, we don't know the distance yet
		}
    }
    
    //we only need to mess with the priority queue in Dijkstra's algorithm
    void relax(int newlySettled, boolean dijkstra){
        //adjacent nodes to the newlySettled node
        LinkedList neighbors = g.adjList[newlySettled-1];
        
        //iterate through every neighbor
        Node neighbor = neighbors.start;
        while(neighbor!=null){
            
            //if Dijkstra, only update a neighbor's distance if it is not already settled
            if(!dijkstra || !alreadySettled(neighbor.index)){
                int newDist = distances[newlySettled-1] + neighbor.weight;
                int oldDist = distances[neighbor.index-1];
                if(newDist < oldDist){
                    distances[neighbor.index-1] = newDist;  //update distance and predecessor
                    predecessors[neighbor.index-1] = newlySettled;  //if newDist shorter than oldDist
                    if(dijkstra)    unsettled.insert(neighbor); //add this node to the priorty queue
                }
            }
            
            neighbor = neighbor.next;
        }
        
        //remove anything in unsettled that is already settled
        //without this code, duplicates can stack in the unsettled queue
        if(dijkstra)
	        while(unsettled.size>0 && alreadySettled(unsettled.min().index))
	            unsettled.delMin();
    }
    
    boolean alreadySettled(int index){
        if(settled[index-1]!=null)  return true;
        return false;
    }
    
    //finds shortest path from start to all other vertices
    //does NOT work if edge weights are negative
    void dijkstra(){
        unsettled = new PriorityQueue(g.n * 2); //double size of priority queue so it will not get overfilled
        settled = new Node[g.n];
        
        //starting point is the first node
        unsettled.insert(g.adjList[0].start);
        
        //until the priority queue is empty...
        while(!unsettled.isEmpty()){
            Node next = unsettled.delMin(); //retrieve next closest neighbor
            settled[next.index-1] = next;   //add to settled list
            numSettled++;
            relax(next.index, true);
        }
    }
    
    //finds shortest path from start to all other vertices
    //DOES work if edge weights are negative
    boolean bellmanFord(){
       for(int i=1; i<=g.n; i++)
            for(int vertex=1; vertex<=g.n; vertex++)
                relax(vertex, false);
        
        return !negativeEdgeCycle();
    }
    
    //returns diameter of graph (uses an adjacency matrix)
    void floydWarshall(){
        int[][] sp = g.adjMatrix;   //sp = shortest path
        
        for(int inter=0; inter<g.n; inter++)
            for(int start=0; start<g.n; start++)
                for(int end=start; end<g.n; end++)
                    //for every start-inter-end combo, if start-inter and inter-end both exist...
                    if(sp[start][inter] != Integer.MAX_VALUE && sp[inter][end] != Integer.MAX_VALUE)
                        //and if (start-inter)+(start-inter) is shorter than (start-end), the we have a new shortest path
                        if(sp[start][end] > sp[start][inter] + sp[inter][end]){
                           sp[start][end] = sp[start][inter] + sp[inter][end];
                           sp[end][start] = sp[end][inter] + sp[inter][start];
                        }
        
        int diameter=0;
        int startD=0;
        int endD=0;
        for(int start=0; start<g.n; start++)
            for(int end=start; end<g.n; end++)
                if(sp[start][end]!=Integer.MAX_VALUE && sp[start][end] > diameter){   //update diameter if larger (but not infinite)
                   diameter = sp[start][end];
                   startD = start;  endD = end;
                }
        
        //the output value of both startD and endD is 1 greater than the stored value
        System.out.println("The diameter of the graph is " + diameter + " between "+ (++startD) +" and "+ (++endD));
    }
    int min(int x, int y){
        if(x<y) return x;
        else    return y;
    }
    int max(int x, int y){
        if(x>y) return x;
        else    return y;
    }
    //returns minimax of graph (uses an adjacency matrix)
    void modifiedFloydWarshall(){
        int[][] sm = g.adjMatrix;   //sm = smallest max
        
        for(int inter=0; inter<g.n; inter++)
            for(int start=0; start<g.n; start++)
                for(int end=start; end<g.n; end++)
                	//for every start-inter-end combo, if the edges start-inter and inter-end both exist...
                    if(sm[start][inter] != Integer.MAX_VALUE && sm[inter][end] != Integer.MAX_VALUE){
                    	//and if the max step of (start-inter)+(inter-end) is smaller than the max step of (start-end), the we have a new minimax
                        int oldmax = sm[start][end];
                        int newmax = max(sm[start][inter], sm[inter][end]);
                        int minimax = min(oldmax, newmax);
                        sm[start][end] = minimax;
                        sm[end][start] = minimax;
                    }

        //if g.n not reachable, there is no minimax
        if(sm[0][g.n-1]==Integer.MAX_VALUE)
			System.out.println("The minimax distance from 1 and "+g.n+" does not exist because "+g.n+" is not connected");
        
        //the output value of both startD and endD is 1 greater than the stored value
       else System.out.println("The minimax distance from 1 and "+ g.n +" is "+ sm[0][g.n-1]);
    }
    
    boolean negativeEdgeCycle(){
        //increment through every node
        for(int node=1; node<=g.n; node++){
            LinkedList neighbors = g.adjList[node-1];
            Node neighbor = neighbors.start;
            //increment through every neighbor
            while(neighbor!=null){
                //check if going through the edge again would further decrease the distance
                int newDist = distances[node-1] + neighbor.weight;
                if(newDist < distances[neighbor.index-1])
                    return true;    //if so, we have a NEGATIVE EDGE CYCLE!!!!!!! OH NOES!!!!!!
                neighbor = neighbor.next;
            }
        }
        return false;   //no negative edge cycles detected
    }
    
    void printShortestPathsBF(){
    	for(int i=1; i<=g.n; i++)
        	if(predecessors[i-1]==0)	//if no predecessor is defined, no path can be found
                System.out.println("No path from "+(start+1)+" to "+i);
            else 	//otherwise print out path
            	System.out.println(shortestPathBF(i));
    }
    int shortestPathBF(int to){
    	int distance;
        int predecessor = predecessors[to-1];
        if((predecessor-1)==start){ //stop once we hit the first node
            System.out.print(predecessor + " ");
            distance = 0;
        }
        else    //if not first node, add distance of previous node to distance of current node
        	distance = shortestPathBF(predecessor);
        System.out.print(to + " ");
        return distances[to-1];
    }
}

//retrieves lowest-key
class PriorityQueue{
    Node[] pq;
    int size=0;
    
    PriorityQueue(int capacity){
        pq = new Node[capacity * 2]; //double size of priority queue so it will not get overfilled
        //I might be a little bit paranoid... or not
    }
    
    void insert(Node n){
        if(size>=pq.length){    //in case queue run out of space, expand it
            Node[] temp = pq;
            pq = new Node[pq.length * 2];
            System.arraycopy(temp, 0, pq, 0, temp.length);
        }
        pq[size] = n;
        swim(size);
        size++;
    }
    
    Node min(){
        return pq[0];
    }
    
    Node delMin(){
        Node min = min();
        if(size!=0){  //priority queue has at least one element
            pq[0] = pq[size - 1];
            pq[size - 1] = null;
            sink(0);
            size--;
        }
        else
            System.out.println("Attempting to remove from empty priority queue");
        return min;
    }
    
    void swim(int n){
        if(pq[n].weight < pq[n/2].weight){
            swap(n, n/2);
            swim(n/2);
        }
    }
    
    void sink(int n){
        int left = n*2+1;
        int right = n*2+2;
        int smallerChild;
        
        if(pq[left]==null)  return;                                         //no children = no compare
        if(pq[right]==null)  smallerChild = left;                           //no right, so compare with left child
        else{                                                               //both children, compare smaller
            if(pq[left].weight < pq[right].weight)  smallerChild = left;
            else                                    smallerChild = right;
        }
        if(pq[n].weight > pq[smallerChild].weight){
            swap(n, smallerChild);
            sink(smallerChild);
        }
    }
    
    void swap(int n, int m){
        Node temp = pq[n];
        pq[n] = pq[m];
        pq[m] = temp;
    }
    
    boolean isEmpty(){
        if(pq[0]==null) return true;
        else return false;
    }
    
    void print(){
        System.out.println("Current unsettled queue:");
        print(0);
        System.out.println("Over");
    }
    
    //for testing purposes
    void print(int n){
        if(n>=size) return;
        System.out.print(pq[n].index+" ");
        print(n*2+1);
        print(n*2+2);
    }
}


//Contains an adjacency matrix and an adjacency list
class Graph {
        int n;
        int[][] adjMatrix;  LinkedList[] adjList;
    
	Graph(int m){
            n = m;
            adjMatrix = new int[n][n];      //we make the matrix bidirection (FW)
            adjList = new LinkedList[n];    //we make the adjacency list unidirectional (Dijkstra, BF)
            for(int i=0; i<n; i++){
                for(int j=i; j<n; j++)
                    if(i==j)adjMatrix[i][j]=0;
                    else{    adjMatrix[i][j]=Integer.MAX_VALUE;
                             adjMatrix[j][i]=Integer.MAX_VALUE; }
                adjList[i] = new LinkedList();
                adjList[i].addNode(new Node(i+1,i+1,0));
            }
        }
        
        void addEdge(int start, int end, int weight){
            adjMatrix[start-1][end-1] = weight;
            adjMatrix[end-1][start-1] = weight;
            adjList[start-1].addNode(new Node(start, end, weight));
        }
        
        boolean containsNegativeEdge(){
            //search for negative in adjList
            for(int i=0; i<n; i++){
                if(adjList[i].containsNegative())
                    return true;
            }
            //no negatives found
            return false;
        }
}

class LinkedList{
        Node start;
        Node last;
        
        LinkedList(){
        }
        
        void addNode(Node n){
            if(start==null){
                start = n;
                last = n;
            }else
                last = last.addNext(n);
        }
        
        boolean containsNegative(){
            if(start==null) return false;
            return start.isNegative();
        }
}

//These are the elements in the linked list
//Ex: the elements of the second linked lists are all of the vertices connected to the second vertex
class Node{
        int predecessor;    //node that "points" to this node
        int index;  //index of node this Node represents
        int weight; //distance from predecessor to this Node
        Node next; //points to the next object
        
        Node(int p, int n, int w){
            predecessor=p;
            index=n;
            weight=w;
        }
        
        Node addNext(Node n){
            next = n;
            return next;
        }
        
        //recursive function
        boolean isNegative(){
            if(weight<0)    return true;
            if(next==null)  return false;
            return next.isNegative();
        }

        void print(){
            System.out.println(predecessor + " " + index + " " + weight);
        }
}