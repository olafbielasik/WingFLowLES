import numpy as np
import matplotlib.pyplot as plt
import os
from scipy.ndimage import gaussian_filter

def compute_vorticity(u, v, dx, dy):
    omega = (np.gradient(v, axis=1) / dx - np.gradient(u, axis=0) / dy)
    return omega

def main():
    base_path = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
    data_path = os.path.join(base_path, "data")
    XX = np.loadtxt(os.path.join(data_path, "XX.txt"))
    YY = np.loadtxt(os.path.join(data_path, "YY.txt"))
    u = np.loadtxt(os.path.join(data_path, "u.txt"))
    v = np.loadtxt(os.path.join(data_path, "v.txt"))
    inside = np.loadtxt(os.path.join(data_path, "inside_mask.txt")).astype(bool)
    
    ny, nx = XX.shape
    dx = (XX[0, -1] - XX[0, 0]) / (nx - 1)
    dy = (YY[-1, 0] - YY[0, 0]) / (ny - 1)
    
    omega = compute_vorticity(u, v, dx, dy)
    omega = gaussian_filter(omega, sigma=1.0)
    omega[inside] = np.nan
    
    finite_vals = omega[np.isfinite(omega)]
    if finite_vals.size > 0:
        print("omega min =", finite_vals.min(), "omega max =", finite_vals.max())
    else:
        print("No finite values in omega.")
    
    plt.figure(figsize=(10, 6))
    levels = np.linspace(-15, 15, 61)
    cf = plt.contourf(XX, YY, omega, levels=levels, cmap='RdBu_r', extend='both')
    plt.colorbar(cf, label='Vorticity')
    plt.title("Flow vorticity around the wing (LES IB, Clark Y, nose up left)")
    plt.xlabel("x")
    plt.ylabel("y")
    plt.axis("equal")
    plt.savefig(os.path.join(data_path, "vorticity_plot.png"), dpi=300)
    plt.show()

if __name__ == "__main__":
    main()










