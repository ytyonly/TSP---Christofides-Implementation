#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <unordered_set>
#include <queue>
#include <unordered_map>
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
    return std::round(std::hypot(a.x - b.x, a.y - b.y));
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
                edges.push_back(Edge(i, j, distance(points[i], points[j])));
                pointIndexToEdgesIndex[key(i,j)] = k++;
            }
        }
    }
    Edge findEdge(int u, int v) {
        int index = pointIndexToEdgesIndex[key(u,v)];
        return edges[index];
    }
};

// vector<int> nearestNeighbor(const vector<Point>& points, vector<vector<int>> dis) {
//     int N = points.size();
//     vector<int> tour;
//     tour.reserve(N);
//     vector<bool> visited(N, false);
    
//     // Start from the first point
//     int current = 0;
//     tour.push_back(current);
//     visited[current] = true;


//     // Build the tour
//     for (int i = 1; i < N; ++i) {
//         int best = -1;
//         for (int j = 0; j < N; ++j) {
//             if (!visited[j] && (best == -1 || dis[current][j] < dis[current][best])) {
//                 best = j;
//             }
//         }
//         current = best;
//         tour.push_back(current);
//         visited[current] = true;
//     }

//     return tour;
// }

void reverseSection(vector<int>& tour, size_t i, size_t k) {
    while (i < k) {
        std::swap(tour[i], tour[k]);
        i++;
        k--;
    }
}

// Function to perform 2-opt optimization
void perform2Opt(vector<int>& tour, Graph graph) {
    bool improvement = true;
    while (improvement) {
        improvement = false;
        int bestDelta = 0;
        size_t bestI = 0, bestK = 0;
        for (size_t i = 0; i < tour.size() - 1; ++i) {
            for (size_t k = i + 1; k < tour.size(); ++k) {
                int delta = - graph.findEdge(tour[i], tour[i + 1]).distance - graph.findEdge(tour[k], tour[(k + 1) % tour.size()]).distance
                            + graph.findEdge(tour[i], tour[k]).distance - graph.findEdge(tour[i + 1], tour[(k + 1) % tour.size()]).distance;
                
                if (delta < bestDelta) {
                    bestDelta = delta;
                    bestI = i;
                    bestK = k;
                }
            }
        }
        if (bestDelta < 0) {
            reverseSection(tour, bestI + 1, bestK);
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

vector<int> christofides(Graph& graph) {
    auto mst = createMST(graph);
    // auto oddDegreeVertices = graph.findOddDegreeVertices();
    // auto minWeightMatching = graph.findMinWeightPerfectMatching(oddDegreeVertices);
    // auto multigraph = graph.combineMSTAndMatching(minWeightMatching);
    // auto eulerCircuit = graph.findEulerianCircuit();
    // auto tour = graph.shortcutEulerianCircuit(eulerCircuit);
    vector<int> tour = {};
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
    perform2Opt(tour, graph);
    // int tourlen = 0;
    // for(int i = 1; i < tour.size(); ++i) {
    //     tourlen += dis[tour[i-1]][tour[i]];
    // }
    // cout<<"tour length is "<<tourlen<<endl;
    for (int index : tour) {
        cout << index << endl;
    }

    return 0;
}