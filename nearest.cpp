#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <cstdlib>
#include <ctime>
#include <climits>

using namespace std;

struct Point {
    double x, y;
    Point(double x, double y) : x(x), y(y) {}
};

struct Path {
    vector<int> tour;
    int distance = 0;
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

// Function to compute the delta for a 3-opt move
int compute3OptDelta(const vector<int>& tour, const vector<vector<int>>& dis, size_t i, size_t j, size_t k) {
    int oldCost = dis[tour[i]][tour[i + 1]] + dis[tour[j]][tour[j + 1]] + dis[tour[k]][tour[(k + 1) % tour.size()]];

    int newCost = dis[tour[i]][tour[k]] + dis[tour[j]][tour[i + 1]] + dis[tour[k]][tour[(j + 1) % tour.size()]];

    return newCost - oldCost;
}

// Function to apply a 3-opt move
void apply3OptMove(vector<int>& tour, size_t i, size_t j, size_t k) {
    reverseSection(tour, i + 1, j);
    reverseSection(tour, j + 1, k);
}

void perform3Opt(vector<int>& tour, const vector<Point>& points, vector<vector<int>> dis) {
    int trynum = 10000;
    while (trynum--) {
        // improvement = false;
        int bestDelta = 0;
        size_t bestI = 0, bestJ = 0, bestK = 0;
        for (size_t i = 0; i < tour.size() - 2; ++i) {
            for (size_t j = i + 1; j < tour.size() - 1; ++j) {
                for (size_t k = j + 1; k < tour.size(); ++k) {
                    int delta = compute3OptDelta(tour, dis, i, j, k);

                    if (delta < bestDelta) {
                        bestDelta = delta;
                        bestI = i;
                        bestJ = j;
                        bestK = k;
                    }
                }
            }
        }
        if (bestDelta < 0) {
            apply3OptMove(tour, bestI, bestJ, bestK);
            // improvement = true;
        }
    }
}

int calPathDis(vector<int> tour, vector<vector<int>> dis) {
    int tourlen = 0;
    for(int i = 1; i < tour.size(); ++i) {
        tourlen += dis[tour[i-1]][tour[i]];
    }
    tourlen += dis[tour[tour.size() - 1]][tour[0]];
    return tourlen;
}
vector<int> searchBetter(Path start, vector<vector<int>> dis) {
    //todo

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

    if(N == 1) {
        cout << 0 << endl;
        return 0;
    }

    vector<vector<int>> dis(N, vector<int>(N, 0));
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            dis[i][j] = distance(points[i], points[j]);
        }
    }

    // Initialize with nearest neighbor and perform 2-opt optimization
    vector<int> tour = nearestNeighbor(points, dis);
    //perform2Opt(tour, points, dis);

    Path start = {tour, calPathDis(tour, dis)};
    // if(N > 960) {
    //     for (int index : tour) {
    //         cout << index << endl;
    //     }
    //     return 0;
    // }
        
    // perform3Opt(tour, points, dis);

    // Output the final tour
    for (int index : tour) {
        cout << index << endl;
    }

    return 0;
}
