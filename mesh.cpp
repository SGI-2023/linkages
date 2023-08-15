#include "mesh.h"
#include <igl/remove_unreferenced.h>

/*
 * Generate a mesh with a square linkage pattern, with dimensions N x M (where N = number of rows of squares, M = number
 * of columns.)
 *
 * Return the vertex positions and face vertex indices in the matrices V and F, respectively.
 */
void generateSquarePatternMesh(int N, int M, Eigen::MatrixXd& V, Eigen::MatrixXi& F) {
    N = 1; // TODO: FIX!
    F.resize(2 * N * M, 3); // Each square is composed of two triangles and there are N * M squares
    std::vector<Eigen::RowVector3d> positions;
    positions.reserve(4 * N * M); // Upper bound of the number of unique vertices, final number should be lower

    const auto linearizeIndex = [numCols = M](int row, int col) { return (row * numCols) + col; };

    const double rowOffset{1.0 / M};
    const double colOffset{1.0 / M};
    int vIdx{0};
    //V.row(vIdx) = Eigen::Vector3d{0.0, rowOffset, 0.0};
    positions.emplace_back(Eigen::RowVector3d{0.0, -rowOffset, 0.0});
    int fIdx{0};
    for (int row = 0; row < N; ++row) {
        for (int col = 0; col < M; ++col) {
            const int squareIdx{linearizeIndex(row, col)};

            if (squareIdx % 2 == 0) { // Even square pattern
                // vIdx is bottom-left
                // V.row(vIdx + 1) = Eigen::Vector3d{(col + 1) * colOffset, -(row + 1) * rowOffset, 0.0}; // Bottom-right
                // V.row(vIdx + 2) = Eigen::Vector3d{col * colOffset, row * rowOffset, 0.0};             // Top-left
                // V.row(vIdx + 3) = Eigen::Vector3d{(col + 1) * colOffset, row * rowOffset, 0.0};       // Top-right

                positions.emplace_back(Eigen::RowVector3d{(col + 1) * colOffset, -(row + 1) * rowOffset, 0.0}); // Bottom-right
                positions.emplace_back(Eigen::RowVector3d{col * colOffset, -row * rowOffset, 0.0});             // Top-left
                positions.emplace_back(Eigen::RowVector3d{(col + 1) * colOffset, -row * rowOffset, 0.0});       // Top-right

                F.row(fIdx) = Eigen::Vector3i{vIdx, vIdx + 1, vIdx + 2};
                F.row(fIdx + 1) = Eigen::Vector3i{vIdx + 2, vIdx + 1, vIdx + 3};
            } else { // Odd square pattern
                // vIdx is top-left
                // V.row(vIdx + 1) = Eigen::Vector3d{(col + 1) * colOffset, -row * rowOffset, 0.0};       // Top-right
                // V.row(vIdx + 2) = Eigen::Vector3d{col * colOffset, -(row + 1) * rowOffset, 0.0};       // Bottom-left
                // V.row(vIdx + 3) = Eigen::Vector3d{(col + 1) * colOffset, (row + 1) * rowOffset, 0.0}; // Bottom-right

                positions.emplace_back(Eigen::RowVector3d{(col + 1) * colOffset, -row * rowOffset, 0.0});
                positions.emplace_back(Eigen::RowVector3d{col * colOffset, -(row + 1) * rowOffset, 0.0});
                positions.emplace_back(Eigen::RowVector3d{(col + 1) * colOffset, -(row + 1) * rowOffset, 0.0});
                
                F.row(fIdx) = Eigen::Vector3i{vIdx, vIdx + 3, vIdx + 1};
                F.row(fIdx + 1) = Eigen::Vector3i{vIdx + 3, vIdx, vIdx + 2};
            }

            fIdx += 2;
            vIdx += 3;
        }
    }

    assert(vIdx + 1 == positions.size());
    V = Eigen::Map<Eigen::MatrixXd>(positions.front().data(), 3, positions.size());
    V.transposeInPlace();

    /*
    std::cout << "Pos=\n";
    for (int i = 0; i < vIdx + 1; ++i) {
        std::cout << "at " << i << ": " << positions[i] << "\n";
    }
    */
    std::cout << "Vertices =\n" << V << "\n";
    std::cout << "Faces =\n" << F << "\n";

}


/* TODO: Feel free to generate other linkage patterns, such as rectangle, parallelogram, etc.! */