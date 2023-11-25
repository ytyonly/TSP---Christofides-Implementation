#include <iostream>
#include <utility>
#include <vector>
#include <cmath>
#include <algorithm>
#include <cstdlib>
#include <ctime>
#include <climits>
#include <random>

using namespace std;

time_t startTime;

struct Point {
    double x, y;
    Point(double x, double y) : x(x), y(y) {}
};

struct Path {
    vector<int> tour;
    int distance = 0;
    Path(vector<int> tour, int distance) : tour(move(tour)), distance(distance) {}
    Path() {}
};

// Function to calculate the Euclidean distance and round it to the nearest integer
int distance(const Point& a, const Point& b) {
    return round(hypot(a.x - b.x, a.y - b.y));
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
void perform2Opt(vector<int>& tour, vector<vector<int>> dis) {
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

void perform3Opt(vector<int>& tour, vector<vector<int>> dis) {
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
vector<int> searchBetter(Path start, vector<vector<int>> dis, int neighborNum) {
    //todo
    vector<Path> neighbors;
    for (int i = 0; i < neighborNum; ++i) {
        //generate neighbor of the start tour
        neighbors.push_back(shuffle(start.tour, dis, 2));
    }
    startTime = time(nullptr);
    while(time(nullptr) - startTime < 1.99) {
        //todo: didn't understand the logic
    }

}

Path shuffle(vector<int> tour, const vector<vector<int>> dis, int num) {
    int tourSize = tour.size();

    vector<int> newTour(tourSize);
    vector<int> keyPoints = tour;

    random_device rd;
    mt19937 g(rd());
    shuffle(keyPoints.begin(), keyPoints.end(), g);
    keyPoints.resize(num);

    for (int i = 0; i < tourSize; i++) {
        auto it = std::find(keyPoints.begin(), keyPoints.end(), i);
        if (it != keyPoints.end()) {
            int keyPointIndex = std::distance(keyPoints.begin(), it);
            int nextKeyPointIndex = (keyPointIndex + 1) % num;
            newTour[i] = tour[keyPoints[nextKeyPointIndex]];
        }
        else {
            newTour[i] = tour[i];
        }
    }
    perform2Opt(newTour, dis);
    Path path(newTour, calPathDis(newTour, dis));
//    for (int index : newTour) {
//        cout << index << endl;
//    }
    return path;
}

struct PathComparator {
    bool operator()(const Path& a, const Path& b) const {
        return a.distance < b.distance;
    }
};

vector<int> immuneAlgorithm(int population, int numofKeypoint, vector<int> &tour, vector<vector<int>> dis) {
    vector<Path> tours;
    Path bestPath(tour, calPathDis(tour, dis));
    PathComparator comp;

    for (int i = 0; i < population; i++) {
        tours.push_back(shuffle(bestPath.tour, dis, numofKeypoint));
    }
    sort(tours.begin(), tours.end(), comp);
    while (time(nullptr) - startTime < 1.99) {
        bestPath = tours[0];
        vector<Path> rest(tours.begin(), tours.begin() + int(population/2));
        tours.clear();
        tours.push_back(bestPath);
        int iteration = 2;

        for (int i = 0; i < int(population/2) * iteration; i++) {
            tours.push_back(shuffle(rest[i%int(population/2)].tour, dis, numofKeypoint));
        }
        sort(tours.begin(), tours.end(), comp);

        if (tours[0].distance < bestPath.distance) {
            bestPath = tours[0];
        }
        tours.resize(population);
    }
    return bestPath.tour;
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

    startTime = time(nullptr);
    // Initialize with nearest neighbor and perform 2-opt optimization
    vector<int> tour = nearestNeighbor(points, dis);
    perform2Opt(tour, dis);
//    vector<int> bestTour = immuneAlgorithm((int(log(N)) == 0) ? 1 : int(log(N)), 2, tour, points, dis);
//please use searchBetter instead
//    vector<int> bestTour = immuneAlgorithm(int(log(N)), 2, tour, dis);

    Path start = {tour, calPathDis(tour, dis)};
    vector<int> bestTour = searchBetter(start, dis, int(log(N)));
        
    // perform3Opt(tour, points, dis);

    // Output the final tour
    for (int index : bestTour) {
        cout << index << endl;
    }

    return 0;
}
