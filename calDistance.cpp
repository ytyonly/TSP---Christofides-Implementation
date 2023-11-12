#include<iostream>
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
    unordered_map<size_t,vector<int>> pointIndexToEdgesIndex;
    vector<vector<double>> cost; //for perfect matching function 
    int n;
    int k = 0;

    Graph(int N): n(N), points(N){}
    void calDistance() {
        cost.resize(n,vector<double>(n,0));
        for(int i = 0; i < n; ++i) {
            for(int j = 0; j < n; ++j) {
                if(i != j) {
                    int dis = distance(points[i], points[j]);
                    edges.push_back(Edge(i, j, dis));
                    pointIndexToEdgesIndex[key(j,i)].push_back(k);
                    pointIndexToEdgesIndex[key(i,j)].push_back(k++);
                    cost[i][j] = cost[j][i] = dis;
                }
            }
        }
    }
    vector<Edge> findEdge(int u, int v) {
        vector<int> index = pointIndexToEdgesIndex[key(u,v)];
        vector<Edge> e;
        for(int i : index) {
            e.push_back(edges[i]);
        }
        return e;
    }
    Edge addEdge(int u, int v) {
        int dis = cost[u][v];
        edges.push_back(Edge(u, v, dis));
        pointIndexToEdgesIndex[key(u,v)].push_back(k);
        pointIndexToEdgesIndex[key(v,u)].push_back(k++);
    }
};
int main() {
    int n;
    cin >> n;
    Graph graph(n);


    for (int i = 0; i < n; ++i) {
        cin >> graph.points[i].x >> graph.points[i].y;
    }
    graph.calDistance();
    vector<int> tour;
    for(int i = 0; i < n; ++i) {
        int t;
        cin >> t;
        tour.push_back(t);
    }
    int tourlen = 0;
    for(int i = 1; i < tour.size(); ++i) {
        tourlen += graph.cost[tour[i-1]][tour[i]];
    }
    tourlen += graph.cost[tour[tour.size() - 1]][tour[0]];
    cout<<"tour length is "<<tourlen<<endl;
}