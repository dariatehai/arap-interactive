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
## Authors

Developed as a group project for a 3D Scanning and Motion Capture course.

- Emil Alizada
- Shining Liu
- Daria Tehai
- Peter Zick​

