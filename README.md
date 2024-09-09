# README
Simulates tensile tests of nanowires using Matplotlib.

### Available simulation parameters

| Parameter Name    | Description                                                                                                        |
|-------------------|--------------------------------------------------------------------------------------------------------------------|
| max_force         | Maximum magnitude of the force applied to the lattice                                                              |
| time_step         | How far the simulation advances each frame. Shorter steps are more accurate but cause the animation to run slower. |
| lattice_width     | Number of atoms along the x-axis of the lattice                                                                    |
| lattice_depth     | Number of atoms along the y-axis of the lattice                                                                    |
| lattice_height    | Number of atoms along the z-axis of the lattice                                                                    |
| simulation_length | Length of the final animation in seconds                                                                           |
| force_x_component | Component of the applied force along the x-axis. Value is relative to the other force components.                  |
| force_y_component | Component of the applied force along the y-axis. Value is relative to the other force components.                  |
| force_z_component | Component of the applied force along the z-axis. Value is relative to the other force components.                  |
| crystal_structure | Crystal structure of the lattice simulated. Accepts the following: "C", "FCC", "BCC", "HCP".                       |

### Supported crystal structures
* Simple cubic (C)
* Body-centered cubic (BCC)
* Face-centered cubic (FCC)
* Hexagonal close packed (HCP)
