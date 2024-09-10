import numpy as np


class PortalFrameElement:
    def __init__(self, E, I, L, theta=0):
        """
        Initialize an element with properties.
        E: Young's Modulus
        I: Moment of Inertia
        L: Length of the element
        theta: Angle of the element in degrees (0 for horizontal, 90 for vertical)
        """
        self.E = E
        self.I = I
        self.L = L
        self.theta = np.radians(theta)

    def stiffness_matrix(self):
        """
        Returns the local stiffness matrix for a 2D frame element.
        """
        k = (self.E * self.I) / (self.L ** 3) * np.array([
            [12, 6 * self.L, -12, 6 * self.L],
            [6 * self.L, 4 * self.L * 2, -6 * self.L, 2 * self.L * 2],
            [-12, -6 * self.L, 12, -6 * self.L],
            [6 * self.L, 2 * self.L * 2, -6 * self.L, 4 * self.L * 2]
        ])
        return k


def assemble_global_stiffness_matrix(elements, global_dof):
    """
    Assemble the global stiffness matrix.
    elements: List of PortalFrameElement instances
    global_dof: Number of global degrees of freedom (size of the global matrix)
    """
    K_global = np.zeros((global_dof, global_dof))

    # For each element, add its stiffness matrix to the global stiffness matrix
    for element in elements:
        k_local = element.stiffness_matrix()
        # Here you need to map the element's local degrees of freedom to the global matrix
        # This requires knowledge of the node connectivity of the structure

        # Example: Adding the local stiffness matrix to the global matrix
        # Assuming a simple mapping of local DOFs to global DOFs
        # K_global[local_to_global_dof_indices] += k_local

    return K_global


def apply_boundary_conditions(K_global, boundary_conditions):
    """
    Apply boundary conditions to the global stiffness matrix.
    boundary_conditions: List of tuples indicating which DOFs are fixed (e.g., [(0, 'fixed'), (3, 'pinned')])
    """
    for dof, condition in boundary_conditions:
        if condition == 'fixed':
            K_global[dof, :] = 0
            K_global[:, dof] = 0
            K_global[dof, dof] = 1  # Keep diagonal as 1 to maintain the system solvable
        elif condition == 'pinned':
            K_global[dof, :] = 0
            K_global[:, dof] = 0
            K_global[dof, dof] = 1  # Similar logic for pinned conditions

    return K_global


def get_member_details(member_name):
    width = float(input(f"Enter width of {member_name} (in meters): "))
    height = float(input(f"Enter height (or depth) of {member_name} (in meters): "))
    length = float(input(f"Enter length of {member_name} (in meters): "))

    moment_of_inertia = (width * (height ** 3)) / 12
    youngs_modulus = float(input(f"Enter Young's modulus of {member_name} (in GPa): ")) * 1e9  # Convert GPa to Pa

    return width, height, length, moment_of_inertia, youngs_modulus


def get_load_details(member_name):
    load_present = input(f"Is there any load acting on {member_name}? (yes/no): ").lower()

    if load_present == 'yes':
        load_type = input(f"Enter load type for {member_name} (point/udl): ").lower()
        if load_type == 'point':
            load_magnitude = float(input(f"Enter magnitude of point load on {member_name} (in kN): "))
            load_location = float(input(f"Enter location of point load on {member_name} (in meters from one end): "))
            return load_type, load_magnitude, load_location
        elif load_type == 'udl':
            load_magnitude = float(input(f"Enter magnitude of UDL on {member_name} (in kN/m): "))
            return load_type, load_magnitude, None
    return None, None, None


def generate_matlab_code(columns, loads):
    matlab_code = """
% Portal Frame Analysis with User Input
clear;
clc;
members = {'Column1', 'Beam1', 'Column2'};
widths = zeros(1, 3);
heights = zeros(1, 3);
lengths = zeros(1, 3);
moments_of_inertia = zeros(1, 3);
youngs_modulus = zeros(1, 3);
load_type = cell(1, 3);
load_magnitude = zeros(1, 3);
load_location = zeros(1, 3);
"""

    for i, (width, height, length, moment_of_inertia, youngs_modulus) in enumerate(columns):
        matlab_code += f"widths({i + 1}) = {width};\n"
        matlab_code += f"heights({i + 1}) = {height};\n"
        matlab_code += f"lengths({i + 1}) = {length};\n"
        matlab_code += f"moments_of_inertia({i + 1}) = {moment_of_inertia};\n"
        matlab_code += f"youngs_modulus({i + 1}) = {youngs_modulus};\n"

    for i, (load_type, load_magnitude, load_location) in enumerate(loads):
        matlab_code += f"load_type{{{i + 1}}} = '{load_type}';\n"
        if load_magnitude is not None:
            matlab_code += f"load_magnitude({i + 1}) = {load_magnitude};\n"
        if load_location is not None:
            matlab_code += f"load_location({i + 1}) = {load_location};\n"

    matlab_code += """
% Plotting the Portal Frame
figure;
hold on;
for i = 1:2
    rectangle('Position', [i-1, 0, widths(i), heights(i)], 'FaceColor', 'b', 'EdgeColor', 'k');
end
rectangle('Position', [0, heights(1), widths(2), 0.1], 'FaceColor', 'r', 'EdgeColor', 'k');
xlim([-0.5, 2]);
ylim([0, max(heights) + 1]);
xlabel('Width (m)');
ylabel('Height (m)');
title('Portal Frame Diagram with Loads');
grid on;
hold off;
"""
    return matlab_code


def main():
    columns = []
    loads = []

    # Get details for the columns and beams
    for member_name in ['Column1', 'Beam1', 'Column2']:
        print(f"Enter details for {member_name}:")
        columns.append(get_member_details(member_name))
        loads.append(get_load_details(member_name))

    # Generate MATLAB code
    matlab_code = generate_matlab_code(columns, loads)
    with open('portal_frame_analysis_with_loads.m', 'w') as f:
        f.write(matlab_code)
    print("MATLAB code has been generated and saved to 'portal_frame_analysis_with_loads.m'.")

    # Create frame elements for stiffness matrix
    elements = [
        PortalFrameElement(E=columns[0][4], I=columns[0][3], L=columns[0][2], theta=90),
        PortalFrameElement(E=columns[1][4], I=columns[1][3], L=columns[1][2], theta=0),
        PortalFrameElement(E=columns[2][4], I=columns[2][3], L=columns[2][2], theta=90)
    ]

    global_dof = 6  # Number of global degrees of freedom, based on nodes
    K_global = assemble_global_stiffness_matrix(elements, global_dof)

    boundary_conditions = [(0, 'fixed'), (1, 'fixed'), (2, 'fixed')]
    K_global = apply_boundary_conditions(K_global, boundary_conditions)

    print("Global Stiffness Matrix after applying boundary conditions:")
    print(K_global)


if __name__ == "__main__":
    main()
