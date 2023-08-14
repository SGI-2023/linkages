#include "mesh.h"

/*
 * Generate a mesh with a square linkage pattern, with dimensions N x M (where N = number of rows of squares, M = number
 * of columns.)
 *
 * Return the vertex positions and face vertex indices in the matrices V and F, respectively.
 */
void generateSquarePatternMesh(int N, int M, Eigen::MatrixXd& V, Eigen::MatrixXi& F) {
    const int numVertices{(N + 1) * (M + 1)};
    V.resize(numVertices, 3);
    const int numFaces{2 * N * M};
    F.resize(numFaces, 3);

    /**
     * (0, 0) (0, 1) (0, 2) ... (0, M - 1)
     * ...
     * (N - 1, 0) (N - 1, 1) (N - 1, 2) ... (N - 1, M - 1)
    */
    const auto linearizeIndex = [numCols = M + 1](int row, int col) {
        return (row * numCols) + col;
    };

    // Create mesh geometry
    const double rowOffset{1.0 / N};
    const double colOffset{1.0 / M};
    for (int row = 0; row < N + 1; ++row) {
        for (int col = 0; col < M + 1; ++col) {
            const int vIdx{linearizeIndex(row, col)};
            V.row(vIdx) = Eigen::Vector3d{col * colOffset - (1.0 / (2.0 * N)), row * rowOffset - (1.0 / (2.0 * M)), 0.0};
        }
    }
    
    // Create mesh topology
    // Note: could merge both loops into one for performance, just think the code is "cleaner" separately.
    int faceIdx{0};
    for (int row = 0; row < N; ++row) {
        for (int col = 0; col < M; ++col) {
            F.row(faceIdx) = Eigen::Vector3i{linearizeIndex(row, col),
                                             linearizeIndex(row, col + 1) , linearizeIndex(row + 1, col + 1)};
            F.row(faceIdx + 1) = Eigen::Vector3i{linearizeIndex(row, col), linearizeIndex(row + 1, col + 1),
                                                 linearizeIndex(row + 1, col)};
            faceIdx += 2;
        }
    }
    assert(faceIdx == numFaces);
}


/* TODO: Feel free to generate other linkage patterns, such as rectangle, parallelogram, etc.! */