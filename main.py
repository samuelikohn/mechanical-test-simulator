import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, FFMpegWriter
import numpy as np
from imageio_ffmpeg import get_ffmpeg_exe
plt.rcParams['animation.ffmpeg_path'] = get_ffmpeg_exe()

# Global list of atoms
atoms = []

# Simulation parameters
max_force = 2
time_step = 0.03
lattice_width = 4
lattice_depth = 4
lattice_height = 8
simulation_length = 3 # seconds
force_x_component = 0
force_y_component = 1
force_z_component = 0
crystal_structure = "FCC" # Supports: BCC, FCC, HCP, C

def main():
    """
    Main function. Creates and displays the animation with Matplotlib.

    Args:
        None

    Returns:
        None
    """

    # Verify that the chosen crystal structure is supported
    if crystal_structure not in ["BCC", "FCC", "C", "HCP"]:
        print("Invalid crystal structure")
        return

    fig = plt.figure(figsize=(8, 8))
    anim = FuncAnimation(fig, animate, interval=40, frames=generate_frames, repeat=True, cache_frame_data=False)
    anim.save('turn_4_animation.mp4', writer=FFMpegWriter(fps=25))


def animate(i):
    """
    Draws each frame of the animation.

    Args:
        i (int): Index for each frame. Received from the generator function 'generate_frames'.

    Returns:
        None
    """
    
    # Clear previous frame
    plt.clf()

    # Render animation in 3D
    ax = plt.axes(projection='3d')
    ax.set_axis_off()

    # Axis limits are based on atom positions
    x_min, y_min, z_min = place_atom(0, 0, 0)
    x_max, y_max, z_max = place_atom(lattice_width, lattice_depth, lattice_height)
    ax.set(xlim=[x_min - 1, x_max + 1], ylim=[y_min - 1, y_max + 1], zlim=[z_min - 1, z_max + 1])

    # Parameters for plotting spheres
    u, v = np.mgrid[0:2 * np.pi:15j, 0:np.pi:10j]

    # Time step. Smaller numbers are more accurate but animate slower
    tick(time_step)

    for atom in atoms:

        # Spheres are centered on each atom and have a radius of 0.5
        ax.plot_surface(atom.x + 0.5 * np.cos(u) * np.sin(v), atom.y + 0.5 * np.sin(u) * np.sin(v), atom.z + 0.5 * np.cos(v), color=atom.color)



def generate_frames():
    """
    Generator function for frame data. Generation stops when the animation reaches the desired number of frames.

    Args:
        None

    Returns:
        Iterator for frame data.
    """
    
    # Initialize atoms in the simulation with maximum applied force
    init_atoms(max_force, force_x_component, force_y_component, force_z_component)

    i = 0
    while i < (simulation_length * 25):
        yield i
        i += 1
        
        
def init_atoms(force, fx, fy, fz):
    """
    Initializes atoms used in the simulation.
    Atoms are arranged in a triangular lattice.
    Force is applied in equal and opposite magnitudes in the direction determined by fx and fy.
    Maximum force is applied on the furthest atom from the center of the lattice and force decreases linearly towards the center.
    Colors atoms based on the magnitude of the applied force.

    Args:
        force (float): Maximum force to apply.
        fx (float): Relative directional component of the applied force
        fy (float): Relative directional component of the applied force
        fz (float): Relative directional component of the applied force

    Returns:
        None
    """

    # Remove any existing atoms
    atoms.clear()

    # Find center
    x_min, y_min, z_min = place_atom(0, 0, 0)
    x_max, y_max, z_max = place_atom(lattice_width, lattice_depth, lattice_height)
    x_center = (x_min + x_max) / 2
    y_center = (y_min + y_max) / 2
    z_center = (z_min + z_max) / 2

    # Normalize force vector
    fx /= np.sqrt(fx ** 2 + fy ** 2 + fz ** 2)
    fy /= np.sqrt(fx ** 2 + fy ** 2 + fz ** 2)
    fz /= np.sqrt(fx ** 2 + fy ** 2 + fz ** 2)

    # Store values between loops
    x_positions = []
    y_positions = []
    z_positions = []
    distances = []

    for i in range(lattice_width):
        for j in range(lattice_depth):
            for k in range(lattice_height):

                # Get atom coordinates and add to list
                x_pos, y_pos, z_pos = place_atom(i, j, k)
                x_positions.append(x_pos)
                y_positions.append(y_pos)
                z_positions.append(z_pos)

                # Calculate distance vector from lattice center to atom
                dx = x_pos - x_center
                dy = y_pos - y_center
                dz = z_pos - z_center

                # Get magnitude of distance vector projected onto force vector and add to list
                dot_product = dx * fx + dy * fy + dz * fz
                distance = np.sqrt((dot_product * fx) ** 2 + (dot_product * fy) ** 2 + (dot_product * fz) ** 2)
                if dot_product < 0:
                    distance = -distance
                distances.append(distance)

    max_distance = max(abs(d) for d in distances)
    num_atoms = lattice_width * lattice_depth * lattice_height
    for i in range(num_atoms):
        
        # Force is distributed linearly across the tensile axis
        f_applied_x = fx * force * distances[i] / max_distance
        f_applied_y = fy * force * distances[i] / max_distance
        f_applied_z = fz * force * distances[i] / max_distance

        # R and B color components are fixed, magnitude of distance/force is proportional to G component
        color = (0.2, abs(distances[i]) / max_distance, 0.4)

        # Create atom
        atom(x_positions[i], y_positions[i], z_positions[i], f_applied_x, f_applied_y, f_applied_z, color)


def place_atom(i, j, k):
    """
    Calculates the position of atoms in a given crystal structure.

    Args:
        i (int): Row index of the atom
        j (int): Column index of the atom
        k (int): Height index of the atom

    Returns:
        x (float): x coordinate of the atom
        y (float): y coordinate of the atom
        z (float): z coordinate of the atom
    """

    match crystal_structure:
        case "BCC":
            x = 2 * i / np.sqrt(3) + (k % 2) / np.sqrt(3)
            y = 2 * j / np.sqrt(3) + (k % 2) / np.sqrt(3)
            z = k / np.sqrt(3)

        case "FCC":
            x = i / np.sqrt(2)
            y = j * np.sqrt(2) + (i % 2) / np.sqrt(2) + (k % 2) / np.sqrt(2)
            z = k / np.sqrt(2)

        case "C":
            x = i
            y = j
            z = k

        case "HCP":
            x = np.sqrt(3) / 2 * i + (k % 2) / (2 * np.sqrt(3))
            y = 0.5 * (i % 2) + j + (k % 2) / 2
            z = k * np.sqrt(2/3)

    return x, y, z


def tick(t):
    """
    Advances the simulation forward by a small amount of time.
    The position and force components of all atoms are updated accordingly.

    Args:
        t (float): Time step to apply to the system.

    Returns:
        None
    """

    # Update positions of all atoms
    for atom in atoms:
        atom.update_position(t)

    # Update force components of all atoms
    for atom in atoms:
        atom.update_force()


class atom:
    def add_force_vector(self, atom):
        """
        Adds the net force components caused by an atom to the force components of the current atom.

        Args:
            atom (atom): The atom that acts on the current atom.

        Returns:
            None
        """

        # Calculate total distance between atoms
        dx = atom.x - self.x
        dy = atom.y - self.y
        dz = atom.z - self.z
        distance = np.sqrt(dx ** 2 + dy ** 2 + dz ** 2)

        # Calculate net force based on the Lennard-Jones potential
        attractive_force = 1 / (distance ** 7)
        repulsive_force = abs(1 / (distance ** 13))
        force = attractive_force - repulsive_force

        # Calculate force components
        fx = force * dx / distance
        fy = force * dy / distance
        fz = force * dz / distance

        # Add to existing force components
        self.force_x += fx
        self.force_y += fy
        self.force_z += fz

    def update_force(self):
        """
        Recalculates the force components of an atom by adding the force vectors caused by all other atoms.

        Args:
            None

        Returns:
            None
        """

        # Set initial force
        self.force_x = self.applied_force_x
        self.force_y = self.applied_force_y
        self.force_z = self.applied_force_z

        # Add force from all other atoms
        for atom in atoms:
            if self != atom:
                self.add_force_vector(atom)

    def update_position(self, t):
        """
        Updates the position of the atom after a small amount of time.

        Args:
            t (float): Time step to advance the position by.

        Returns:
            None
        """

        # Update position
        self.x += self.force_x * t
        self.y += self.force_y * t
        self.z += self.force_z * t

    def __eq__(self, atom):
        """
        Equality operator for atom objects.

        Args:
            atom (atom): Atom that the current atom is being compared to.

        Returns:
            True if both x and y coordinates are the same.
            False otherwise
        """
        
        return self.x == atom.x and self.y == atom.y and self.z == atom.z
    
    def __init__(self, x, y, z, applied_force_x, applied_force_y, applied_force_z, color):
        """
        Constructor for the atom object.
        Assigns values to its x, y, and z coordinates and its color.
        Calculates x, y, and z force components from all existing atoms.

        Args:
            x (float): The x coordinate of the atom.
            y (float): The y coordinate of the atom.
            z (float): The z coordinate of the atom.
            applied_force_x (float): External force applied to the atom in the x direction.
            applied_force_y (float): External force applied to the atom in the y direction.
            applied_force_z (float): External force applied to the atom in the z direction.

        Returns:
            atom (atom): New atom object.
        """

        # Set initial position
        self.x = x
        self.y = y
        self.z = z

        # Set initial force
        self.force_x = applied_force_x
        self.force_y = applied_force_y
        self.force_z = applied_force_z

        # Stored for use in 'update_force' method
        self.applied_force_x = applied_force_x
        self.applied_force_y = applied_force_y
        self.applied_force_z = applied_force_z

        # Set color
        self.color = color

        # Add force of all existing atoms to self and from self to all existing atoms
        for atom in atoms:
            self.add_force_vector(atom)
            atom.add_force_vector(self)

        atoms.append(self)


if __name__ == "__main__":
    main()
