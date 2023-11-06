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

struct Blossom {
    int v;  // Vertex number
    int p;  // Parent of the blossom
    int base;  // Base of the blossom
    bool matched;  // True if the blossom is matched
};

// Function to find a perfect matching using the blossom algorithm
vector<Edge> perfectMatching(Graph& graph) {
    int n = graph.n;
    vector<Blossom> blossoms(n);
    vector<int> match(n, -1);  // Matching information
    vector<bool> used(n, false);  // Used in augmenting path
    vector<int> parent(n, -1);  // Parent in the tree
    queue<int> q;

    for (int u = 0; u < n; ++u) {
        if (match[u] == -1) {
            parent[u] = -1;
            used.assign(n, false);
            q.push(u);

            while (!q.empty()) {
                int v = q.front();
                q.pop();

                for (const Edge& edge : graph.edges) {
                    if (edge.u == v || edge.v == v) {
                        int to = (edge.u == v) ? edge.v : edge.u;

                        if (!used[to]) {
                            used[to] = true;
                            q.push(match[to]);

                            if (match[to] == -1) {
                                while (v != -1) {
                                    int p = parent[v];
                                    int pp = match[p];
                                    match[v] = p;
                                    match[p] = v;
                                    v = pp;
                                }
                            } else {
                                parent[to] = v;
                            }
                        } else if (blossoms[to].base != blossoms[v].base) {
                            int cur = v;
                            while (cur != -1) {
                                Blossom& b = blossoms[blossoms[cur].base];
                                b.matched = !b.matched;
                                cur = parent[b.base];
                            }
                            cur = to;
                            while (cur != -1) {
                                Blossom& b = blossoms[blossoms[cur].base];
                                b.matched = !b.matched;
                                cur = parent[b.base];
                            }
                        }
                    }
                }
            }
        }
    }

    vector<Edge> matching;
    for (int u = 0; u < n; ++u) {
        if (match[u] != -1 && u < match[u]) {
            matching.push_back(graph.findEdge(u, match[u]));
        }
    }

    return matching;
}

vector<int> christofides(Graph& graph) {
    Graph MST = createMST(graph);
    vector<int> oddDegreeVertices = findOdd(MST);
    // vector<Edge> minWeightMatching = findMinWeightPerfectMatching(oddDegreeVertices);
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