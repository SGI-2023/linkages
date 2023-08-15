#include "mesh.h"

/*
 * Create |V| x |F| vertex-face adjacency matrix.
 */
Eigen::SparseMatrix<int> buildVertexFaceAdjacency(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F) {

    Eigen::SparseMatrix<int> A(V.rows(), F.rows());
    std::vector<Eigen::Triplet<int>> triplets;
    for (int fi = 0; fi < F.rows(); fi++) {
        int d = F.cols();
        for (int j = 0; j < d; j++) {
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
        int d = F.cols();
        for (int j = 0; j < d; j++) {
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
        int d = F.cols();
        for (int j = 0; j < d; j++) {
            int idx = edges(F(fi, j), F(fi, (j + 1) % d), E);
            triplets.emplace_back(fi, abs(idx), idx > 0 ? 1 : -1);
        }
    }
    A.setFromTriplets(triplets.begin(), triplets.end());
    return A;
}

/*
 * Cut a quad mesh along an edge from vertex <a> to vertex <b>, duplicating vertex <a> but leaving <b> intact.
 * Essentially turning the edge between a and b into a hinge.
 */
void cutHinge(int a, int b, Eigen::MatrixXd& V, Eigen::MatrixXi& E, Eigen::MatrixXi& F, Eigen::SparseMatrix<int>& FE) {

    // Duplicate vertex <a>, adding the duplicate to the end of the vertex list.
    int nV = V.rows();
    V.conservativeResize(nV + 1, 3);
    V.row(nV) = V.row(a);

    // Find all faces with an edge between <a> and <b>, and edit topology.
    Eigen::VectorXi edgeVec = Eigen::VectorXi::Zero(E.rows());
    int eIdx = edges(a, b, E);
    edgeVec(abs(eIdx)) = (eIdx > 0) ? 1 : -1;
    Eigen::VectorXi faceVec = FE * edgeVec;
    if (faceVec.sum() == 0) {
        // don't cut along boundary edges
        for (int i = 0; i < F.rows(); i++) {
            if (faceVec(i) != 0) {
                // pick an arbitrary face whose edges connect to the new vertex
                for (int k = 0; k < F.cols(); k++) {
                    if (F(i, k) == a) {
                        F(i, k) = nV;
                        break;
                    }
                }
                break;
            }
        }
    }

    // Recompute edges.
    igl::edges(F, E);
    FE = buildFaceEdgeAdjacency(E, F);
}


/*
 * Cut a quad mesh along a sequence of edges following the vertex indices in <cut>. Leaves both vertex endpoints of the
 * cut intact, creating a "slit"-like hole in the mesh.
 */
void cutSlit(const std::vector<int>& cut, Eigen::MatrixXd& V, Eigen::MatrixXi& E, Eigen::MatrixXi& F,
             Eigen::SparseMatrix<int>& FE) {

    // Duplicate middle vertices, adding the duplicates to the end of the vertex list.
    int nCut = cut.size(); // number of cut vertices
    int nV = V.rows();
    V.conservativeResize(nV + nCut - 2, 3);
    for (int i = 1; i < nCut - 1; i++) V.row((nV - 1) + i) = V.row(cut[i]);

    // Edit topology. This one is a bit trickier because we can't re-hook edges of an arbitrary face, we have to make
    // sure faces are consistently connected on either side of the slit. Argh, I don't like matrix-based meshes.
    // First get all faces on one side of the cut, and record which cut vertices are incident.
    Eigen::VectorXi edgeVec;
    std::vector<int> cutFaces;
    std::vector<std::array<int, 2>> incVerts;
    for (int i = 1; i < nCut; i++) {
        int eIdx = edges(cut[i - 1], cut[i], E);
        edgeVec = Eigen::VectorXi::Zero(E.rows());
        edgeVec(abs(eIdx)) = (eIdx > 0) ? 1 : -1;
        Eigen::VectorXi faceVec = FE * edgeVec;
        if (faceVec.sum() != 0) continue; // don't cut along boundary edges
        for (int fi = 0; fi < F.rows(); fi++) {
            // Only re-hook faces on the right side of the cut.
            if (faceVec(fi) <= 0) continue;
            cutFaces.push_back(fi);
            incVerts.push_back({i - 1, i});
            break;
        }
    }
    // Now re-hook up the edges on the faces.
    for (int k = 0; k < cutFaces.size(); k++) {
        for (int ei = 0; ei < 2; ei++) {
            int cutIdx = incVerts[k][ei];
            if (cutIdx == 0) continue;
            if (cutIdx == nCut - 1) continue;
            for (int d = 0; d < F.cols(); d++) {
                if (F(cutFaces[k], d) == cut[cutIdx]) {
                    F(cutFaces[k], d) = (nV - 1) + cutIdx;
                    break;
                }
            }
        }
    }

    // Recompute edges.
    igl::edges(F, E);
    FE = buildFaceEdgeAdjacency(E, F);
}

/*
 * Generate a regular grid of N x M rectangles (where N = number of rows of squares, M = number of columns.)
 */
void buildRegularGrid(int N, int M, Eigen::MatrixXd& V, Eigen::MatrixXi& F, double width = 1.0, double height = 1.0) {

    V.resize((N + 1) * (M + 1), 3);
    F.resize(N * M, 4); // N x M quads

    // Write vertices (unit squares)
    double maxX = M * width * 0.5;
    double maxY = N * height * 0.5;
    for (int i = 0; i < N + 1; i++) {
        for (int j = 0; j < M + 1; j++) {
            Eigen::Vector3d p = {-maxX + j * width, maxY - i * height, 0};
            V.row(i * (M + 1) + j) = p;
        }
    }
    // Write quad faces, whose vertices are always in the same order:
    //
    // 0---3
    // |   |
    // 1---2
    //
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < M; j++) {
            F(i * M + j, 0) = i * (M + 1) + j;
            F(i * M + j, 1) = (i + 1) * (M + 1) + j;
            F(i * M + j, 2) = (i + 1) * (M + 1) + (j + 1);
            F(i * M + j, 3) = i * (M + 1) + (j + 1);
        }
    }
}

/*
 * Generate a mesh with a square linkage pattern, with dimensions N x M (where N = number of rows of squares, M = number
 * of columns.)
 *
 * Return the vertex positions and face vertex indices in the matrices V and F, respectively.
 */
void generateSquarePatternMesh(int N, int M, Eigen::MatrixXd& V, Eigen::MatrixXi& F) {

    buildRegularGrid(N, M, V, F);
    int nV = V.rows();
    int offset = 0;
    std::vector<int> dupVerts(nV); // store new index of duplicated vertex
    std::iota(begin(dupVerts), end(dupVerts), 0);

    // To speedup computation, create adjacency matrices.
    Eigen::MatrixXi E;
    igl::edges(F, E);
    Eigen::SparseMatrix<int> FE = buildFaceEdgeAdjacency(E, F);

    // Cut mesh! Make "horizontal" cuts first.
    // Iterate over each row of vertices, ignoring the top and bottom of boundary vertices.
    for (int i = 1; i < N; i++) {
        if ((M + i) % 2 == 0) {
            for (int j = 0; j < M; j += 2) {
                std::vector<int> indices = {i * (M + 1) + j, i * (M + 1) + (j + 1), i * (M + 1) + (j + 2)};
                cutSlit(indices, V, E, F, FE);

                dupVerts[i * (M + 1) + (j + 1)] = nV + offset;
                offset += 1;
            }
        } else {
            cutHinge(i * (M + 1) + 0, i * (M + 1) + 1, V, E, F, FE);
            dupVerts[i * (M + 1) + 0] = nV + offset;
            offset += 1;
            cutHinge(i * (M + 1) + M, i * (M + 1) + (M - 1), V, E, F, FE);
            dupVerts[i * (M + 1) + M] = nV + offset;
            offset += 1;
            for (int j = 1; j < M - 1; j += 2) {
                std::vector<int> indices = {i * (M + 1) + j, i * (M + 1) + (j + 1), i * (M + 1) + (j + 2)};
                cutSlit(indices, V, E, F, FE);

                dupVerts[i * (M + 1) + (j + 1)] = nV + offset;
                offset += 1;
            }
        }
    }
    // Now make vertical cuts. The logic is the same as for the horizontal cuts, just with the iteration order swapped
    // (iterate over columns instead of rows.)
    for (int j = 1; j < M; j++) {
        if ((N + j) % 2 != 0) {
            for (int i = 0; i < N; i += 2) {
                std::vector<int> indices = {i * (M + 1) + j, (i + 1) * (M + 1) + j, dupVerts[(i + 2) * (M + 1) + j]};
                cutSlit(indices, V, E, F, FE);
            }
        } else {
            cutHinge(0 * (M + 1) + j, dupVerts[1 * (M + 1) + j], V, E, F, FE);
            cutHinge(N * (M + 1) + j, (N - 1) * (M + 1) + j, V, E, F, FE);
            for (int i = 1; i < N - 1; i += 2) {
                std::vector<int> indices = {i * (M + 1) + j, (i + 1) * (M + 1) + j, dupVerts[(i + 2) * (M + 1) + j]};
                cutSlit(indices, V, E, F, FE);
            }
        }
    }

    // Triangulate quads with alternating diagonals. I don't really have a good reason for this.
    Eigen::MatrixXi tF(2 * F.rows(), 3); // tri faces
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < M; j++) {
            Eigen::Vector3i LL = {F(i * M + j, 0), F(i * M + j, 1), F(i * M + j, 2)}; // lower left triangle
            Eigen::Vector3i UR = {F(i * M + j, 0), F(i * M + j, 2), F(i * M + j, 3)}; // upper right triangle
            Eigen::Vector3i UL = {F(i * M + j, 0), F(i * M + j, 1), F(i * M + j, 3)}; // upper left triangle
            Eigen::Vector3i LR = {F(i * M + j, 3), F(i * M + j, 1), F(i * M + j, 2)}; // lower right triangle
            bool sgn = (i % 2 == 0) == (j % 2 == 0);
            Eigen::Vector3i T1 = sgn ? LL : LR;
            Eigen::Vector3i T2 = sgn ? UR : UL;
            tF(i * 2 * M + 2 * j, 0) = T1(0);
            tF(i * 2 * M + 2 * j, 1) = T1(1);
            tF(i * 2 * M + 2 * j, 2) = T1(2);

            tF(i * 2 * M + (2 * j + 1), 0) = T2(0);
            tF(i * 2 * M + (2 * j + 1), 1) = T2(1);
            tF(i * 2 * M + (2 * j + 1), 2) = T2(2);
        }
    }

    F = tF;
}


/* TODO: Feel free to generate other linkage patterns, such as rectangle, parallelogram, etc.! */