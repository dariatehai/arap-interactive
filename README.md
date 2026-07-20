# ARAP Interactive Model

An interactive 3D mesh manipulation tool based on the As-Rigid-As-Possible (ARAP) algorithm. This project allows users to deform 3D meshes while maintaining their rigidity as much as possible.

The project contains:

- a custom local-global ARAP solver implemented with Eigen;
- a reference implementation based on libigl;
- an interactive GUI for loading meshes and manipulating one or more handles;
- runtime benchmarks over meshes with different resolutions;
- numerical comparison against libigl using relative RMS vertex error.

The project was developed for a 3D Scanning and Motion Capture course.

---

## Features

### Interactive deformation

- Load a triangular mesh from the GUI.
- Add one or more handle vertices by vertex index.
- Drag handle markers directly in the viewer.
- Remove the most recently added handle or clear all handles.
- Keep an anchored mesh region fixed during deformation.
- Recenter the camera.
- Reset the mesh to its undeformed state.

### ARAP implementations

The project includes two solver paths:

1. **Custom ARAP solver**
   - Cotangent-weighted deformation energy
   - Per-vertex local rotation estimation using SVD
   - Sparse global position solve
   - Reusable sparse factorization while handle positions change

2. **libigl ARAP**
   - Used as a reference implementation
   - Supports visual and numerical comparison with the custom solver

### Evaluation tools

- Runtime measurement for a fixed number of local-global iterations
- Comparison across mesh resolutions
- Relative RMS vertex error against libigl
- Optional per-iteration energy output for convergence inspection

---

## Building

From the project root:

```bash
cmake -S . -B build
cmake --build build --config Release
```

On multi-configuration generators such as Visual Studio, use:

```bash
cmake --build build --config Release
```

On single-configuration generators, configure the build type explicitly:

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build
```

---

## Troubleshooting

### CMake cannot find libigl

Check that libigl is installed or included as a submodule and that the project’s top-level `CMakeLists.txt` adds it correctly.

### Eigen headers are missing

Ensure Eigen is provided through libigl or installed separately and linked to the relevant targets.

### The viewer opens with no mesh

Use the GUI’s mesh-loading control and select a supported triangular mesh.

### A handle cannot be added

Check that:

- a mesh is loaded;
- the vertex index is valid;
- the vertex is not already a handle;
- the vertex does not belong to the fixed region.

### Dragging is slow

- Build in Release mode.
- Use a lower-resolution mesh.
- Reduce the number of local-global iterations per mouse update.
- Disable per-iteration console logging.
- Avoid changing the handle set during the drag.

### Sparse factorization fails

Possible causes include:

- insufficient positional constraints;
- disconnected components;
- degenerate geometry;
- problematic cotangent weights;
- an underconstrained free component.

Add suitable anchors and validate the input mesh.

---

## Authors

Developed as a group project for a 3D Scanning and Motion Capture course.

- Emil Alizada
- Shining Liu
- Daria Tehai
- Peter Zick​

