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

int main() {
    int N;
    cin >> N;
    Graph graph(N);


    for (int i = 0; i < N; ++i) {
        cin >> graph.points[i].x >> graph.points[i].y;
    }
    graph.calDistance();

    Graph MST = createMST(graph);
    for(Edge e : MST.edges) {
        cout<<e.u<<" "<<e.v<<" "<<e.distance<<endl;
    }
}