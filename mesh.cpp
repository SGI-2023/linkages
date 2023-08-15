#include "mesh.h"

/*
 * Create |V| x |F| vertex-face adjacency matrix.
 */
Eigen::SparseMatrix<int> buildVertexFaceAdjacency(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F) {

    Eigen::SparseMatrix<int> A(V.rows(), F.rows());
    std::vector<Eigen::Triplet<int>> triplets;
    for (int fi = 0; fi < F.rows(); fi++) {
        for (int j = 0; j < 3; j++) {
            triplets.emplace_back(F(fi, j), fi, 1);
        }
    }
    A.setFromTriplets(triplets.begin(), triplets.end());
    return A;
}

/*
 * Create |F| x |V| face-vertex adjacency matrix (the transpose of the above, but still column-order.)
 */
Eigen::SparseMatrix<int> buildFaceVertexAdjacency(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F) {

    Eigen::SparseMatrix<int> A(F.rows(), V.rows());
    std::vector<Eigen::Triplet<int>> triplets;
    for (int fi = 0; fi < F.rows(); fi++) {
        for (int j = 0; j < 3; j++) {
            triplets.emplace_back(fi, F(fi, j), 1);
        }
    }
    A.setFromTriplets(triplets.begin(), triplets.end());
    return A;
}

/*
 * Return the set of oriented edges between vertices a and b. The sign indicates the orientation of the edge relative to
 * a->b.
 */
int edges(int a, int b, const Eigen::MatrixXi& E) {

    for (int i = 0; i < E.rows(); i++) {
        if (E(i, 0) == a && E(i, 1) == b) {
            return i;
        } else if (E(i, 0) == b && E(i, 1) == a) {
            return -i;
        }
    }
    throw std::logic_error("edge not found between vertices " + std::to_string(a) + " and " + std::to_string(b) + ".");
}

/*
 * Create |F| x |E| face-edge adjacency matrix.
 */
Eigen::SparseMatrix<int> buildFaceEdgeAdjacency(const Eigen::MatrixXi& E, const Eigen::MatrixXi& F) {

    Eigen::SparseMatrix<int> A(F.rows(), E.rows());
    std::vector<Eigen::Triplet<int>> triplets;
    for (int fi = 0; fi < F.rows(); fi++) {
        // Without more mesh structure, I'm not sure how to make this search operation faster.
        for (int j = 0; j < 3; j++) {
            int idx = edges(F(fi, j), F(fi, (j + 1) % 3), E);
            triplets.emplace_back(fi, abs(idx), idx > 0 ? 1 : -1);
        }
    }
    A.setFromTriplets(triplets.begin(), triplets.end());
    return A;
}

/*
 * Cut a mesh along an edge from vertex <a> to vertex <b>, duplicating vertex <a> but leaving <b> intact. Essentially
 * turning the edge between a and b into a hinge.
 */
void cutHinge(int a, int b, Eigen::MatrixXd& V, Eigen::MatrixXi& E, Eigen::MatrixXi& F, Eigen::SparseMatrix<int>& FE) {

    std::cerr << "cutting hinge" << std::endl;
    // Duplicate vertex <a>, adding the duplicate to the end of the vertex list.
    int nV = V.rows();
    V.conservativeResize(nV + 1, 3);
    V.row(nV) = V.row(a);

    // Find all faces with an edge between <a> and <b>, and edit topology.
    Eigen::VectorXi edgeVec = Eigen::VectorXi::Zero(E.rows());
    int eIdx = edges(a, b, E);
    edgeVec(abs(eIdx)) = (eIdx > 0) ? 1 : -1;
    Eigen::VectorXi faceVec = FE * edgeVec;
    for (int i = 0; i < F.cols(); i++) {
        if (faceVec(i) == 0) continue;
        // pick an arbitrary face whose edges connect to the new vertex
        for (int k = 0; k < 3; k++) {
            if (F(i, k) == a) {
                F(i, k) = nV;
                break;
            }
        }
        break;
    }

    // Recompute edges.
    igl::edges(F, E);
    FE = buildFaceEdgeAdjacency(E, F);
}


/*
 * Cut a mesh along a sequence of edges following the vertex indices in <cut>. Leaves both vertex endpoints of the cut
 * intact, creating a "slit"-like hole in the mesh.
 */
void cutSlit(const std::vector<int>& cut, Eigen::MatrixXd& V, Eigen::MatrixXi& E, Eigen::MatrixXi& F,
             Eigen::SparseMatrix<int>& FE) {

    std::cerr << "cutting slit" << std::endl;
    // Duplicate middle vertices, adding the duplicates to the end of the vertex list.
    int nCut = cut.size(); // number of cut vertices
    int nV = V.rows();
    V.conservativeResize(nV + nCut - 2, 3);
    for (int i = 1; i < nCut - 1; i++) V.row(nV + (i - 1)) = V.row(cut[i]);

    // Edit topology. This one is a bit trickier because we can't re-hook edges of an arbitrary face, we have to make
    // sure faces are consistently connected on either side of the slit. Argh, I don't like matrix-based meshes.
    Eigen::VectorXi edgeVec;
    for (int i = 1; i < nCut - 1; i++) {
        int eIdx = edges(cut[i], cut[i + 1], E);
        edgeVec = Eigen::VectorXi::Zero(E.rows());
        edgeVec(abs(eIdx)) = (eIdx > 0) ? 1 : -1;
        Eigen::VectorXi faceVec = FE * edgeVec;
        for (int fi = 0; fi < F.cols(); fi++) {
            // Only re-hook faces on the right side of the cut.
            if (faceVec(fi) <= 0) continue;
            for (int k = 0; k < 3; k++) {
                if (F(fi, k) == cut[i]) {
                    F(fi, k) = nV + (i - 1);
                    break;
                }
            }
            break;
        }
    }

    // Recompute edges.
    igl::edges(F, E);
    FE = buildFaceEdgeAdjacency(E, F);
}

/*
 * Generate a mesh with a square linkage pattern, with dimensions N x M (where N = number of rows of squares, M = number
 * of columns.)
 *
 * Return the vertex positions and face vertex indices in the matrices V and F, respectively.
 */
void generateSquarePatternMesh(int N, int M, Eigen::MatrixXd& V, Eigen::MatrixXi& F) {

    // Following the square pattern in the YouTube video "Auxetics & Metamaterials"
    // (https://www.youtube.com/watch?v=2tKj7wQtXxY).

    V.resize((N + 1) * (M + 1), 3);
    F.resize(2 * N * M, 3); // N x M squares, which each get cut up into two triangle

    // Write a regular grid first, with unit squares.
    double maxX = M * 0.5;
    double maxY = N * 0.5;
    for (int i = 0; i < N + 1; i++) {
        for (int j = 0; j < M + 1; j++) {
            Eigen::Vector3d p = {-maxX + j, maxY - i, 0};
            V.row(i * (M + 1) + j) = p;
        }
    }
    // Alternate diagonals. I don't really have a good reason for this.
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < M; j++) {
            // lower left triangle
            Eigen::Vector3i LL = {i * (M + 1) + j, (i + 1) * (M + 1) + j, (i + 1) * (M + 1) + (j + 1)};
            // upper right triangle
            Eigen::Vector3i UR = {i * (M + 1) + j, (i + 1) * (M + 1) + (j + 1), i * (M + 1) + (j + 1)};
            // upper left triangle
            Eigen::Vector3i UL = {i * (M + 1) + j, (i + 1) * (M + 1) + j, i * (M + 1) + (j + 1)};
            // lower right triangle
            Eigen::Vector3i LR = {i * (M + 1) + (j + 1), (i + 1) * (M + 1) + j, (i + 1) * (M + 1) + (j + 1)};

            bool sgn = (i % 2 == 0) == (j % 2 == 0);
            Eigen::Vector3i T1 = sgn ? LL : LR;
            Eigen::Vector3i T2 = sgn ? UR : UL;

            F(i * 2 * M + 2 * j, 0) = T1(0);
            F(i * 2 * M + 2 * j, 1) = T1(1);
            F(i * 2 * M + 2 * j, 2) = T1(2);

            F(i * 2 * M + (2 * j + 1), 0) = T2(0);
            F(i * 2 * M + (2 * j + 1), 1) = T2(1);
            F(i * 2 * M + (2 * j + 1), 2) = T2(2);
        }
    }

    // To speedup computation, create adjacency matrices.
    Eigen::MatrixXi E;
    igl::edges(F, E);
    Eigen::SparseMatrix<int> FE = buildFaceEdgeAdjacency(E, F);

    // Cut up the mesh! Create all horizontal cuts first.
    for (int i = 1; i < N; i++) {
        int j = 0;
        while (j < M) {
            if (i % 2 == 0) {
                if (j == 0) {
                    cutHinge(i * (M + 1) + j, i * (M + 1) + (j + 1), V, E, F, FE);
                    j += 1;
                } else if (j == M - 1) {
                    cutHinge(i * (M + 1) + (j + 1), i * (M + 1) + j, V, E, F, FE);
                    j += 1;
                } else if (j < M) {
                    std::vector<int> indices = {i * (M + 1) + j, i * (M + 1) + (j + 1), i * (M + 1) + (j + 2)};
                    cutSlit(indices, V, E, F, FE);
                    j += 2;
                }
            } else {
                if (j == M - 1) {
                    cutHinge(i * (M + 1) + (j + 1), i * (M + 1) + j, V, E, F, FE);
                    j += 1;
                } else if (j < M) {
                    std::vector<int> indices = {i * (M + 1) + j, i * (M + 1) + (j + 1), i * (M + 1) + (j + 2)};
                    cutSlit(indices, V, E, F, FE);
                    j += 2;
                }
            }
        }
    }

    // TODO: Now create vertical cuts.

    // Add random noise to vertices (debugging)
    for (int i = 0; i < V.rows(); i++) {
        V(i, 2) += ((double)rand() / (RAND_MAX));
    }
    std::cerr << "mesh created" << std::endl;
}


/* TODO: Feel free to generate other linkage patterns, such as rectangle, parallelogram, etc.! */