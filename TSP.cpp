#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <unordered_set>
#include <queue>
#include <unordered_map>
#include <stack>
#include <limits.h>
#include "Matching.h"
#include "Graph.h"
#include <chrono>

using namespace std;

struct Point {
    double x, y;
    Point(){}
    Point(double X, double Y) : x(X), y(Y) {}
};

struct myEdge {
    int u, v;
    int distance;
    myEdge(int U, int V, int Distance) : u(U), v(V), distance(Distance) {}
};

int distance(const Point& a, const Point& b) {
    return round(hypot(a.x - b.x, a.y - b.y));
}

inline size_t key(int i,int j) {return (size_t) i << 32 | (unsigned int) j;}

struct myGraph {
    vector<Point> points;
    vector<myEdge> edges;
    unordered_map<size_t,vector<int>> pointIndexToEdgesIndex;
    vector<vector<double>> cost;
    int n;
    int k = 0;

    myGraph(int N): n(N), points(N){}
    void calDistance() {
        cost.resize(n,vector<double>(n,0));
        for(int i = 0; i < n; ++i) {
            for(int j = 0; j < n; ++j) {
                if(i != j) {
                    int dis = distance(points[i], points[j]);
                    edges.push_back(myEdge(i, j, dis));
                    pointIndexToEdgesIndex[key(j,i)].push_back(k);
                    pointIndexToEdgesIndex[key(i,j)].push_back(k++);
                    cost[i][j] = dis;
                }
            }
        }
    }
    vector<myEdge> findEdge(int u, int v) {
        vector<int> index = pointIndexToEdgesIndex[key(u,v)];
        vector<myEdge> e;
        for(int i : index) {
            e.push_back(edges[i]);
        }
        return e;
    }
    void addEdge(int u, int v) {
        int dis = cost[u][v];
        edges.push_back(myEdge(u, v, dis));
        pointIndexToEdgesIndex[key(u,v)].push_back(k);
        pointIndexToEdgesIndex[key(v,u)].push_back(k++);
    }
};


void reverseSection(vector<int>& tour, size_t i, size_t k) {
    while (i < k) {
        std::swap(tour[i], tour[k]);
        i++;
        k--;
    }
}

void perform2Opt(vector<int>& tour, myGraph& graph) {
    bool improvement = true;
    while (improvement) {
        improvement = false;
        int bestDelta = 0;
        size_t bestI = 0, bestK = 0;

        for (size_t i = 0; i < tour.size() - 1; ++i) {
            for (size_t k = i + 1; k < tour.size(); ++k) {
                // Calculate the change in distance if a 2-opt swap is performed
                int delta = graph.findEdge(i, k)[0].distance;

//                int delta = calculateDelta(tour, graph, i, k);

                // Check if this is the best improvement seen so far
                if (delta < bestDelta) {
                    bestDelta = delta;
                    bestI = i;
                    bestK = k;
                }
            }
        }

        // If a better route is found, perform the swap
        if (bestDelta < 0) {
            // Reversing the segment of the tour between bestI+1 and bestK
            std::reverse(tour.begin() + bestI + 1, tour.begin() + bestK + 1);
            improvement = true;
        }
    }
}

struct compareEdge {
    bool operator()(const myEdge& e1, const myEdge& e2) {
        return e1.distance > e2.distance;
    }
};

myGraph createMST(myGraph& graph) {
    int N = graph.n;
    unordered_set<int> visited;
    myGraph MST(N);
    MST.points = graph.points;

    //choose the points[0] as the first vertice to find MST
    int current = 0;
    visited.insert(0);

    priority_queue<myEdge, vector<myEdge>, compareEdge> pq;

    while(visited.size() < N) {
        int minDis = INT_MAX;
        for(int i = 0; i < N; ++i) {

            //skip the visited vertice
            if(visited.find(i) != visited.end()) {
                continue;
            }

            //add all adjacent edge of current vertice to priority queue
            myEdge curEdge = graph.findEdge(current, i)[0];
            pq.push(curEdge);
        }
        while(!pq.empty()) {
            myEdge minEdge = pq.top();
            pq.pop();
            int u = minEdge.u;
            int v = minEdge.v;
            if(visited.find(u) == visited.end()) {
                visited.insert(u);
                current = u;
                MST.edges.push_back(minEdge);
                break;
            }
            if(visited.find(v) == visited.end()) {
                visited.insert(v);
                current = v;
                MST.edges.push_back(minEdge);
                break;
            }
        }
    }
    return MST;
}

vector<int> findOdd(myGraph& graph) {
    int N = graph.n;
    vector<int> degrees(N, 0);

    for (myEdge edge : graph.edges) {
        degrees[edge.u]++;
        degrees[edge.v]++;
    }

    vector<int> odd;
    for (int i = 0; i < N; i++) {
        if (degrees[i] % 2 != 0) {
            odd.push_back(i);
        }
    }

    return odd;
}

vector<myEdge> perfectMatching_greedy(myGraph& graph, const vector<int>& oddVertices) {
    vector<myEdge> matching;
    unordered_set<int> matched; // To keep track of vertices that have already been matched

    for (int u : oddVertices) {
        if (matched.find(u) != matched.end()) {
            continue; // Skip if already matched
        }

        int minDistance = INT_MAX;
        int partner = -1;

        for (int v : oddVertices) {
            if (u != v && matched.find(v) == matched.end()) {
                int dist = graph.findEdge(u, v)[0].distance;
                if (dist < minDistance) {
                    minDistance = dist;
                    partner = v;
                }
            }
        }

        if (partner != -1) {
            matched.insert(u);
            matched.insert(partner);
            matching.push_back(myEdge(u, partner, minDistance));
        }
    }
    return matching;
}

vector<myEdge> perfectMatching(myGraph& graph, const vector<int>& oddVertices) {
    Graph O(oddVertices.size());
	vector<double> costO;
	for(int i = 0; i < (int)oddVertices.size(); i++)
	{
		for(int j = i+1; j < (int)oddVertices.size(); j++)
		{
            O.AddEdge(i, j);
            costO.push_back( graph.cost[oddVertices[i]][oddVertices[j]] );
		}
	}
    Matching M(O);
	auto p = M.SolveMinimumCostPerfectMatching(costO);
	list<int> matching = p.first;
    
    vector<myEdge> minWeightMatching;
    for(list<int>::iterator it = matching.begin(); it != matching.end(); it++) {
        pair<int, int> e = O.GetEdge( *it );
        int a = e.first;
        int b = e.second;
        myEdge e1 = graph.findEdge(oddVertices[a],oddVertices[b])[0];
        minWeightMatching.push_back(e1);
    }

    return minWeightMatching;
}
vector<int> findEulerianCircuit(myGraph& graph, const vector<myEdge>& matching) {
    // Add matching edges to the graph to make all degrees even
    for (const myEdge& e : matching) {
        graph.edges.push_back(e);
    }

    stack<int> stack;
    stack.push(0);
    unordered_set<int> usedEdges; // Use the edge index to mark as used

    vector<int> circuit;

    while (!stack.empty()) {
        int v = stack.top();

        bool found = false;
        for (size_t i = 0; i < graph.edges.size(); ++i) {
            myEdge& e = graph.edges[i];

            if ((e.u == v || e.v == v) && usedEdges.find(i) == usedEdges.end()) {
                stack.push(e.u == v ? e.v : e.u);
                usedEdges.insert(i);  // Mark edge as used by its index
                found = true;
                break;
            }
        }

        if (!found) {
            stack.pop();
            if (!stack.empty()) { // Avoid adding start vertex twice
                circuit.push_back(v);
            }
        }
    }
    if (circuit.size() == 0) {
        circuit.push_back(0);
    }
    return circuit; // The circuit might not need to be reversed depending on your requirements
}

vector<int> shortcutEulerianCircuit(const vector<int>& eulerCircuit) {
    vector<int> tour;
    unordered_set<int> visited;
    for (int vertex : eulerCircuit) {
        // If we have not visited this vertex before, add it to the Hamiltonian circuit
        if (visited.find(vertex) == visited.end()) {
            tour.push_back(vertex);
            visited.insert(vertex); // Mark the vertex as visited
        }
    }
    return tour;
}

vector<int> christofides_Test(myGraph& graph) {
    myGraph MST = createMST(graph);

    // cout<< "MST" << endl;
    // for(myEdge e : MST.edges) {
    //     cout<<e.u<<" "<<e.v<<" "<<e.distance<<endl;
    // }

    vector<int> oddDegreeVertices = findOdd(MST);
    //Create a graph with the odd degree vertices
	Graph O(oddDegreeVertices.size());
	vector<double> costO;
	for(int i = 0; i < (int)oddDegreeVertices.size(); i++)
	{
		for(int j = i+1; j < (int)oddDegreeVertices.size(); j++)
		{
            O.AddEdge(i, j);
            costO.push_back( graph.cost[oddDegreeVertices[i]][oddDegreeVertices[j]] );
		}
	}
    Matching M(O);
	auto p = M.SolveMinimumCostPerfectMatching(costO);
	list<int> matching = p.first;
    
    //pair<list<int>, double> solution = M.SolveMinimumCostPerfectMatching(graph.cost);
    // vector<myEdge> minWeightMatching = perfectMatching(graph, oddDegreeVertices);

    // cout<< "Matching" << endl;
    // for(int m : matching) {
    //     cout<< m<< endl;
    // }
    vector<myEdge> minWeightMatching;
    for(list<int>::iterator it = matching.begin(); it != matching.end(); it++) {
        pair<int, int> e = O.GetEdge( *it );
        int a = e.first;
        int b = e.second;
        myEdge e1 = graph.findEdge(a,b)[0];
        minWeightMatching.push_back(e1);
        // Debug output
    }
    // for(myEdge e : minWeightMatching) {
    //     cout<<e.u<<" "<<e.v<<" "<<e.distance<<endl;
    // }

    auto eulerCircuit = findEulerianCircuit(MST, minWeightMatching);
//    auto tour = shortcutEulerianCircuit(eulerCircuit);
//    return tour;
    return eulerCircuit;
}

vector<int> christofides(myGraph& graph) {
    myGraph MST = createMST(graph);
    vector<int> oddDegreeVertices = findOdd(MST);
    vector<myEdge> minWeightMatching = perfectMatching(graph, oddDegreeVertices);
    auto eulerCircuit = findEulerianCircuit(MST, minWeightMatching);
    auto tour = shortcutEulerianCircuit(eulerCircuit);
    return tour;
}

vector<int> christofides_bigInput(myGraph& graph) {
    myGraph MST = createMST(graph);
    vector<int> oddDegreeVertices = findOdd(MST);
    vector<myEdge> minWeightMatching = perfectMatching_greedy(graph, oddDegreeVertices);
    auto eulerCircuit = findEulerianCircuit(MST, minWeightMatching);
    auto tour = shortcutEulerianCircuit(eulerCircuit);
    return tour;
}

int main() {
    int N;
    cin >> N;
    myGraph graph(N);


    for (int i = 0; i < N; ++i) {
        cin >> graph.points[i].x >> graph.points[i].y;
    }
    graph.calDistance();

    vector<int> tour;
    if(N >= 960) {
        tour = christofides_bigInput(graph);
    } else {
        tour = christofides(graph);
    }
    //vector<int> tour = christofides_Test(graph);

    perform2Opt(tour, graph);
    // int tourlen = 0;
    // for(int i = 1; i < tour.size(); ++i) {
    //     tourlen += graph.cost[tour[i-1]][tour[i]];
    // }
    // tourlen += graph.cost[tour[tour.size()-1]][tour[0]];
    // cout<<"tour length is "<<tourlen<<endl;
    for (int index : tour) {
        cout << index << endl;
    }

    return 0;
}