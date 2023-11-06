#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
using namespace std;

struct Point {
    double x, y;
    Point(double x, double y) : x(x), y(y) {}
};

// Function to calculate the Euclidean distance and round it to the nearest integer
int distance(const Point& a, const Point& b) {
    return std::round(std::hypot(a.x - b.x, a.y - b.y));
}

vector<int> nearestNeighbor(const vector<Point>& points, vector<vector<int>> dis) {
    int N = points.size();
    vector<int> tour;
    tour.reserve(N);
    vector<bool> visited(N, false);
    
    // Start from the first point
    int current = 0;
    tour.push_back(current);
    visited[current] = true;


    // Build the tour
    for (int i = 1; i < N; ++i) {
        int best = -1;
        for (int j = 0; j < N; ++j) {
            if (!visited[j] && (best == -1 || dis[current][j] < dis[current][best])) {
                best = j;
            }
        }
        current = best;
        tour.push_back(current);
        visited[current] = true;
    }

    return tour;
}

void reverseSection(vector<int>& tour, size_t i, size_t k) {
    while (i < k) {
        std::swap(tour[i], tour[k]);
        i++;
        k--;
    }
}

// Function to perform 2-opt optimization
void perform2Opt(vector<int>& tour, const vector<Point>& points, vector<vector<int>> dis) {
    bool improvement = true;
    while (improvement) {
        improvement = false;
        int bestDelta = 0;
        size_t bestI = 0, bestK = 0;
        for (size_t i = 0; i < tour.size() - 1; ++i) {
            for (size_t k = i + 1; k < tour.size(); ++k) {
                int delta = - dis[tour[i]][tour[i + 1]] - dis[tour[k]][tour[(k + 1) % tour.size()]]
                            + dis[tour[i]][tour[k]] + dis[tour[i + 1]][tour[(k + 1) % tour.size()]];
                
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

int main() {
    int N;
    cin >> N;
    vector<Point> points;
    points.reserve(N);

    for (int i = 0; i < N; ++i) {
        double x, y;
        cin >> x >> y;
        points.emplace_back(x, y);
    }

    vector<vector<int>> dis(N, vector<int>(N, 0));
    for(int i = 0; i < N; ++i) {
        for(int j = 0; j < N; ++j) {
            dis[i][j] = distance(points[i], points[j]);
        }
    }
    

    vector<int> tour = nearestNeighbor(points, dis);
    perform2Opt(tour, points, dis);
    int tourlen = 0;
    for(int i = 1; i < tour.size(); ++i) {
        tourlen += dis[tour[i-1]][tour[i]];
    }
    cout<<"tour length is "<<tourlen<<endl;
    for (int index : tour) {
        cout << index << endl;
    }

    return 0;
}