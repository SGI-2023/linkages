#include <igl/arap.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/opengl/glfw/imgui/ImGuiPlugin.h>
#include <igl/readOFF.h>
#include <igl/read_triangle_mesh.h>
#include <igl/unproject_onto_mesh.h>

#include "mesh.h"

// ==  geometric data
Eigen::MatrixXd VERTICES(0, 3); // |V| x 3 matrix; the i-th row contains the 3D position of the i-th vertex
Eigen::MatrixXi FACES(0, 3); // |F| x 3 matrix; the i-th row contains the vertex indices of the i-th face, in CCW order

// == viz parameters
Eigen::RowVector3d PURPLE(80.0 / 255.0, 64.0 / 255.0, 255.0 / 255.0);
Eigen::RowVector3d RED(1.0, 0.0, 0.0);
Eigen::RowVector3d GOLD(255.0 / 255.0, 228.0 / 255.0, 58.0 / 255.0);

// ==  solve data
std::vector<int> ANCHORS;  // indices of vertices which have been selected as anchors
std::vector<int> VIDX_MAP; // maps vertex index to index in ANCHORS
Eigen::MatrixXd HANDLES;
igl::ARAPData ARAP_DATA;
bool IS_DRAGGING = false;

// ==  program parameters
bool INTERNAL_MESH = false;
std::string MESH_FILEPATH;
Eigen::RowVector3f LAST_MOUSE;
bool SELECT_ANCHORS = true;
int LAST_VERTEX;
int N_ROWS = 4;
int N_COLS = 4;


void update_vis(igl::opengl::glfw::Viewer& viewer) {
    Eigen::MatrixXd points(ANCHORS.size(), 3);
    for (size_t i = 0; i < ANCHORS.size(); i++) points.row(i) = VERTICES.row(ANCHORS[i]);
    viewer.data().set_points(points, PURPLE);
}

void update_ARAP(igl::opengl::glfw::Viewer& viewer) {
    igl::arap_solve(HANDLES, ARAP_DATA, VERTICES);
    viewer.data().set_vertices(VERTICES);
    viewer.data().compute_normals();
}

void update_anchor_positions() {
    // Update positions of anchors
    HANDLES.resize(ANCHORS.size(), 3);
    VIDX_MAP.resize(VERTICES.rows());
    for (size_t i = 0; i < HANDLES.rows(); i++) {
        HANDLES.row(i) = VERTICES.row(ANCHORS[i]);
        VIDX_MAP[ANCHORS[i]] = i;
    }
}

bool key_pressed(igl::opengl::glfw::Viewer& viewer, unsigned char key, int mods) {
    switch (key) {
        case ' ': {
            SELECT_ANCHORS = !SELECT_ANCHORS;

            update_anchor_positions();

            break;
        }
        default: {
            return false;
        }
    }
    return true;
}

bool mouse_down(igl::opengl::glfw::Viewer& viewer, int button, int mods) {

    // Cast a ray in the view direction starting from the mouse position
    double x = viewer.current_mouse_x;
    double y = viewer.core().viewport(3) - viewer.current_mouse_y;
    LAST_MOUSE = Eigen::RowVector3f(x, y, 0);
    int fid;            // index of selected face
    Eigen::Vector3f bc; // barycentric coords in selected face
    if (igl::unproject_onto_mesh(Eigen::Vector2f(x, y), viewer.core().view, viewer.core().proj, viewer.core().viewport,
                                 VERTICES, FACES, fid, bc)) {

        // Get closest vertex. Could do this calculation purely intrinsically from barycentric coords and edge lengths,
        // but just use the extrinsic calculation.
        const Eigen::RowVector3d m3 = VERTICES.row(FACES(fid, 0)) * bc(0) + VERTICES.row(FACES(fid, 1)) * bc(1) +
                                      VERTICES.row(FACES(fid, 2)) * bc(2);
        int cid = 0;
        Eigen::Vector3d((VERTICES.row(FACES(fid, 0)) - m3).squaredNorm(),
                        (VERTICES.row(FACES(fid, 1)) - m3).squaredNorm(),
                        (VERTICES.row(FACES(fid, 2)) - m3).squaredNorm())
            .minCoeff(&cid);
        LAST_VERTEX = FACES(fid, cid);

        if (SELECT_ANCHORS) {
            auto it = std::find(ANCHORS.begin(), ANCHORS.end(), LAST_VERTEX);
            if (it == ANCHORS.end()) {
                ANCHORS.push_back(LAST_VERTEX);
                Eigen::VectorXi b(ANCHORS.size());
                for (size_t i = 0; i < ANCHORS.size(); i++) b(i) = ANCHORS[i];
                igl::arap_precomputation(VERTICES, FACES, VERTICES.cols(), b, ARAP_DATA);
            } else {
                ANCHORS.erase(it);
            }
        }

        IS_DRAGGING = true;
        update_vis(viewer);
        return true;
    }
    return false;
}

bool mouse_up(igl::opengl::glfw::Viewer& viewer, int button, int mods) {
    IS_DRAGGING = false;
    return false;
}

bool mouse_move(igl::opengl::glfw::Viewer& viewer, int button, int mods) {
    if (IS_DRAGGING && !SELECT_ANCHORS) {

        Eigen::RowVector3f drag_mouse(viewer.current_mouse_x, viewer.core().viewport(3) - viewer.current_mouse_y,
                                      LAST_MOUSE(2));
        Eigen::RowVector3f drag_scene, last_scene;
        igl::unproject(drag_mouse, viewer.core().view, viewer.core().proj, viewer.core().viewport, drag_scene);
        igl::unproject(LAST_MOUSE, viewer.core().view, viewer.core().proj, viewer.core().viewport, last_scene);
        LAST_MOUSE = drag_mouse;

        HANDLES.row(VIDX_MAP[LAST_VERTEX]) += (drag_scene - last_scene).cast<double>();

        update_ARAP(viewer);
        update_vis(viewer);
        return true;
    }
    return false;
}

void setup(igl::opengl::glfw::Viewer& viewer) {
    // ARAP precomputation
    ARAP_DATA.max_iter = 100;
    ARAP_DATA.with_dynamics = true;
    // ARAP_DATA.energy = igl::ARAP_ENERGY_TYPE_SPOKES;
    ARAP_DATA.energy = igl::ARAP_ENERGY_TYPE_SPOKES_AND_RIMS;
    // ARAP_DATA.energy = igl::ARAP_ENERGY_TYPE_ELEMENTS; // triangles or tets
    Eigen::VectorXi b(0);
    // igl::arap_precomputation(VERTICES, FACES, VERTICES.cols(), b, ARAP_DATA);

    viewer.data().clear();
    viewer.data().set_mesh(VERTICES, FACES);
}

int main(int argc, char* argv[]) {

    if (argc > 1) {
        MESH_FILEPATH = argv[1];
        INTERNAL_MESH = false;
        // Load in mesh
        igl::read_triangle_mesh(MESH_FILEPATH, VERTICES, FACES);
    } else {
        INTERNAL_MESH = true;
    }

    // Attach a menu plugin
    igl::opengl::glfw::Viewer viewer;
    igl::opengl::glfw::imgui::ImGuiPlugin plugin;
    viewer.plugins.push_back(&plugin);
    igl::opengl::glfw::imgui::ImGuiMenu menu;
    plugin.widgets.push_back(&menu);

    // Add content to the default menu
    menu.callback_draw_viewer_menu = [&]() {
        // Draw parent menu content
        menu.draw_viewer_menu();

        if (ImGui::CollapsingHeader("Mesh generation", ImGuiTreeNodeFlags_DefaultOpen)) {
            ImGui::InputInt("Rows", &N_ROWS);
            ImGui::InputInt("Cols", &N_COLS);
            if (ImGui::Button("Square linkage")) {
                generateSquarePatternMesh(N_ROWS, N_COLS, VERTICES, FACES);
                setup(viewer);
            }
        }

        ImGui::Text("Mode: %s", SELECT_ANCHORS ? "Anchor selection" : "ARAP deformation");
    };

    if (!INTERNAL_MESH) {
        setup(viewer);
    }

    // Set up callbacks
    viewer.data().set_face_based(true);
    viewer.core().is_animating = true;
    viewer.callback_key_pressed = &key_pressed;
    viewer.callback_mouse_down = &mouse_down;
    viewer.callback_mouse_up = &mouse_up;
    viewer.callback_mouse_move = &mouse_move;
    viewer.launch();
}
