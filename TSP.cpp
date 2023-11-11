#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <unordered_set>
#include <queue>
#include <unordered_map>
#include <stack>
#include <limits.h>

using namespace std;

struct Point {
    double x, y;
    Point(){}
    Point(double X, double Y) : x(X), y(Y) {}
};

struct Edge {
    int u, v;
    int distance;
    Edge(int U, int V, int Distance) : u(U), v(V), distance(Distance) {}
};

int distance(const Point& a, const Point& b) {
    return round(hypot(a.x - b.x, a.y - b.y));
}

inline size_t key(int i,int j) {return (size_t) i << 32 | (unsigned int) j;}

struct Graph {
    vector<Point> points;
    vector<Edge> edges;
    unordered_map<size_t,int> pointIndexToEdgesIndex;
    int n;

    Graph(int N): n(N), points(N){}
    void calDistance() {
        int k = 0;
        for(int i = 0; i < n; ++i) {
            for(int j = 0; j < n; ++j) {
                if(i != j) {
                    edges.push_back(Edge(i, j, distance(points[i], points[j])));
                    pointIndexToEdgesIndex[key(j,i)] = k;
                    pointIndexToEdgesIndex[key(i,j)] = k++;
                }
            }
        }
    }
    Edge findEdge(int u, int v) {
        int index = pointIndexToEdgesIndex[key(u,v)];
        return edges[index];
    }
};


void reverseSection(vector<int>& tour, size_t i, size_t k) {
    while (i < k) {
        std::swap(tour[i], tour[k]);
        i++;
        k--;
    }
}

void perform2Opt(vector<int>& tour, Graph& graph) {
    bool improvement = true;
    while (improvement) {
        improvement = false;
        int bestDelta = 0;
        size_t bestI = 0, bestK = 0;

        for (size_t i = 0; i < tour.size() - 1; ++i) {
            for (size_t k = i + 1; k < tour.size(); ++k) {
                // Calculate the change in distance if a 2-opt swap is performed
                int delta = graph.findEdge(i, k).distance;

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
    bool operator()(const Edge& e1, const Edge& e2) {
        return e1.distance > e2.distance;
    }
};

Graph createMST(Graph& graph) {
    int N = graph.n;
    unordered_set<int> visited;
    Graph MST(N);
    MST.points = graph.points;

    //choose the points[0] as the first vertice to find MST
    int current = 0;
    visited.insert(0);

    priority_queue<Edge, vector<Edge>, compareEdge> pq;

    while(visited.size() < N) {
        int minDis = INT_MAX;
        for(int i = 0; i < N; ++i) {

            //skip the visited vertice
            if(visited.find(i) != visited.end()) {
                continue;
            }

            //add all adjacent edge of current vertice to priority queue
            Edge curEdge = graph.findEdge(current, i);
            pq.push(curEdge);
        }
        while(!pq.empty()) {
            Edge minEdge = pq.top();
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

vector<int> findOdd(Graph& graph) {
    int N = graph.n;
    vector<int> degrees(N, 0);

    for (Edge edge : graph.edges) {
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

vector<Edge> perfectMatching(Graph& graph, const vector<int>& oddVertices) {
    vector<Edge> matching;
    unordered_set<int> matched; // To keep track of vertices that have already been matched

    for (int u : oddVertices) {
        if (matched.find(u) != matched.end()) {
            continue; // Skip if already matched
        }

        int minDistance = INT_MAX;
        int partner = -1;

        for (int v : oddVertices) {
            if (u != v && matched.find(v) == matched.end()) {
                int dist = graph.findEdge(u, v).distance;
                if (dist < minDistance) {
                    minDistance = dist;
                    partner = v;
                }
            }
        }

        if (partner != -1) {
            matched.insert(u);
            matched.insert(partner);
            matching.push_back(Edge(u, partner, minDistance));
        }
    }

    return matching;
}
vector<int> findEulerianCircuit(Graph& graph, const vector<Edge>& matching) {
    // Add matching edges to the graph to make all degrees even
    for (const Edge& e : matching) {
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
            Edge& e = graph.edges[i];

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

vector<int> christofides_Test(Graph& graph) {
    Graph MST = createMST(graph);

    cout<< "MST" << endl;
    for(Edge e : MST.edges) {
        cout<<e.u<<" "<<e.v<<" "<<e.distance<<endl;
    }

    vector<int> oddDegreeVertices = findOdd(MST);
    vector<Edge> minWeightMatching = perfectMatching(graph, oddDegreeVertices);

    cout<< "Matching" << endl;
    for(Edge e : minWeightMatching) {
        cout<<e.u<<" "<<e.v<<" "<<e.distance<<endl;
    }

    auto eulerCircuit = findEulerianCircuit(MST, minWeightMatching);
//    auto tour = shortcutEulerianCircuit(eulerCircuit);
//    return tour;
    return eulerCircuit;
}

vector<int> christofides(Graph& graph) {
    Graph MST = createMST(graph);
    vector<int> oddDegreeVertices = findOdd(MST);
    vector<Edge> minWeightMatching = perfectMatching(graph, oddDegreeVertices);
    auto eulerCircuit = findEulerianCircuit(MST, minWeightMatching);
    auto tour = shortcutEulerianCircuit(eulerCircuit);
    return tour;
}

int main() {
    int N;
    cin >> N;
    Graph graph(N);


    for (int i = 0; i < N; ++i) {
        cin >> graph.points[i].x >> graph.points[i].y;
    }
    graph.calDistance();

    vector<int> tour = christofides(graph);
//     vector<int> tour = christofides_Test(graph);

//     perform2Opt(tour, graph);
//     int tourlen = 0;
//     for(int i = 1; i < tour.size(); ++i) {
//         tourlen += graph.findEdge(tour[i-1],tour[i]).distance;
//     }
//     cout<<"tour length is "<<tourlen<<endl;
    for (int index : tour) {
        cout << index << endl;
    }

    return 0;
}