import numpy as np
import os

def read_clarky(filename):
    with open(filename, 'r') as f:
        lines = f.readlines()
    sections = []
    current_section = []
    for line in lines:
        stripped = line.strip()
        if not stripped:
            if current_section:
                sections.append(current_section)
                current_section = []
            continue
        parts = stripped.split()
        if len(parts) < 2:
            continue
        try:
            x, y = float(parts[0]), float(parts[1])
            current_section.append([x, y])
        except ValueError:
            pass
    if current_section:
        sections.append(current_section)
    upper = np.array(sections[0])
    lower = np.array(sections[1])
    if np.allclose(upper[-1], lower[0]):
        airfoil = np.vstack((upper, lower[1:]))
    else:
        airfoil = np.vstack((upper, lower))
    return airfoil

def rotate_points(points, angle_deg):
    angle_rad = np.deg2rad(angle_deg)
    rot = np.array([[np.cos(angle_rad), -np.sin(angle_rad)],
                    [np.sin(angle_rad), np.cos(angle_rad)]])
    return points @ rot.T

def generate_nonuniform_grid(xmin, xmax, nx, ymin, ymax, ny, refine_region):
    x1, x2 = refine_region
    n1 = nx // 3
    n2 = nx // 3
    n3 = nx - n1 - n2
    x_part1 = np.linspace(xmin, x1, n1, endpoint=False)
    x_part2 = np.linspace(x1, x2, n2, endpoint=False)
    x_part3 = np.linspace(x2, xmax, n3, endpoint=True)
    x_vals = np.concatenate([x_part1, x_part2, x_part3])
    y_vals = np.linspace(ymin, ymax, ny)
    XX, YY = np.meshgrid(x_vals, y_vals)
    return XX, YY

def generate_airfoil_mesh(angle_deg=-20.0, chord=1.5,
                          domain_size_x=(-0.5, 3.0),
                          domain_size_y=(-1.0, 1.0),
                          nx=600, ny=400, 
                          refine_region=(0.0, 1.5)):
    base_path = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
    data_path = os.path.join(base_path, "data")
    data_file = os.path.join(data_path, "airfoil.dat")
    
    airfoil_points = read_clarky(data_file)
    airfoil_points *= chord
    if angle_deg != 0.0:
        airfoil_points = rotate_points(airfoil_points, angle_deg)
    min_x = airfoil_points[:,0].min()
    airfoil_points[:,0] -= min_x
    
    xmin, xmax = domain_size_x
    ymin, ymax = domain_size_y
    XX, YY = generate_nonuniform_grid(xmin, xmax, nx, ymin, ymax, ny, refine_region)
    
    np.savetxt(os.path.join(data_path, "XX.txt"), XX)
    np.savetxt(os.path.join(data_path, "YY.txt"), YY)
    np.savetxt(os.path.join(data_path, "airfoil_points.txt"), airfoil_points)
    
    return XX, YY, airfoil_points

if __name__ == "__main__":
    XX, YY, airfoil = generate_airfoil_mesh(angle_deg=-20.0)
    print("Grid:", XX.shape)
    print("Number of profile points:", len(airfoil))









