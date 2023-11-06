import sys
import matplotlib.pyplot as plt

def visualize(points, mst_edges):
    plt.figure(figsize=(8, 8))

    # Visualize the input points
    for x, y in points:
        plt.scatter(x, y, c='blue', marker='o')

    # Visualize the edges in the Minimum Spanning Tree
    for i, j in mst_edges:
        x1, y1 = points[i]
        x2, y2 = points[j]
        plt.plot([x1, x2], [y1, y2], 'r-')

    plt.title('Minimum Spanning Tree Visualization')
    plt.xlabel('X-axis')
    plt.ylabel('Y-axis')
    plt.grid()
    plt.show()

def main():
    # Read the input points as (x, y) coordinates
    num_points = int(input("Enter the number of points: "))
    points = []

    for _ in range(num_points):
        x, y = map(float, input().split())
        points.append((x, y))

    # Read the Minimum Spanning Tree as pairs of indexes (i, j)
    mst_edges = []
    num_edges = int(input("Enter the number of edges in the MST: "))
    for _ in range(num_edges):
        i, j = map(int, input().split())
        mst_edges.append((i, j))

    # Visualize the input points and the MST
    visualize(points, mst_edges)

if __name__ == "__main__":
    main()