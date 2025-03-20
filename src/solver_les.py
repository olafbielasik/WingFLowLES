import numpy as np
import os
from numba import njit
import time

@njit
def point_in_polygon(px, py, polygon):
    inside = False
    n = len(polygon)
    for i in range(n-1):
        x1, y1 = polygon[i]
        x2, y2 = polygon[i+1]
        if ((y1 > py) != (y2 > py)):
            x_int = x1 + (py - y1) * (x2 - x1) / (y2 - y1)
            if x_int > px:
                inside = not inside
    return inside

@njit
def build_airfoil_mask_numba(XX, YY, polygon):
    ny, nx = XX.shape
    mask = np.zeros((ny, nx), dtype=np.bool_)
    for j in range(ny):
        for i in range(nx):
            px = XX[j, i]
            py = YY[j, i]
            if point_in_polygon(px, py, polygon):
                mask[j, i] = True
    return mask

@njit
def smagorinsky_nu_t(u, v, dx, dy, Cs):
    ny, nx = u.shape
    nu_t = np.zeros((ny, nx), dtype=np.float32)
    for j in range(1, ny-1):
        for i in range(1, nx-1):
            dudx = (u[j, i+1] - u[j, i-1]) / (2 * dx)
            dudy = (u[j+1, i] - u[j-1, i]) / (2 * dy)
            dvdx = (v[j, i+1] - v[j, i-1]) / (2 * dx)
            dvdy = (v[j+1, i] - v[j-1, i]) / (2 * dy)
            Sxx = dudx
            Syy = dvdy
            Sxy = 0.5 * (dudy + dvdx)
            S_mag = np.sqrt(Sxx**2 + Syy**2 + 2*(Sxy**2))
            delta = np.sqrt(dx * dy)
            nu_t[j, i] = (Cs * delta)**2 * S_mag
            if nu_t[j, i] > 0.1:
                nu_t[j, i] = 0.1
    return nu_t

@njit
def solve_poisson_pressure(p, u_star, v_star, dx, dy, dt, max_iter=100, tol=1e-4):
    ny, nx = p.shape
    beta = dx / dy
    for _ in range(max_iter):
        p_old = p.copy()
        for j in range(1, ny-1):
            for i in range(1, nx-1):
                div = ((u_star[j, i+1] - u_star[j, i-1]) / (2 * dx) +
                       (v_star[j+1, i] - v_star[j-1, i]) / (2 * dy))
                p[j, i] = ((p[j, i+1] + p[j, i-1]) +
                           beta**2 * (p[j+1, i] + p[j-1, i]) -
                           dx * dy * div / dt) / (2 + 2 * beta**2)
        for j in range(ny):
            p[j, 0] = p[j, 1]
            p[j, nx-1] = p[j, nx-2]
        for i in range(nx):
            p[0, i] = p[1, i]
            p[ny-1, i] = p[ny-2, i]
        err = np.sqrt(np.sum((p - p_old)**2))
        if err < tol:
            break
    return p

@njit
def step_fractional(u, v, p, inside_mask, dx, dy, dt, nu, Cs):
    ny, nx = u.shape
    nu_t = smagorinsky_nu_t(u, v, dx, dy, Cs)
    nu_eff = nu + nu_t

    u_star = u.copy()
    v_star = v.copy()
    for j in range(1, ny-1):
        for i in range(1, nx-1):
            if inside_mask[j, i]:
                u_star[j, i] = 0.0
                v_star[j, i] = 0.0
            else:
                dudx = (u[j, i+1] - u[j, i-1]) / (2 * dx)
                dudy = (u[j+1, i] - u[j-1, i]) / (2 * dy)
                dvdx = (v[j, i+1] - v[j, i-1]) / (2 * dx)
                dvdy = (v[j+1, i] - v[j-1, i]) / (2 * dy)
                conv_u = u[j, i] * dudx + v[j, i] * dudy
                conv_v = u[j, i] * dvdx + v[j, i] * dvdy
                lap_u = ((u[j, i+1] - 2 * u[j, i] + u[j, i-1]) / (dx * dx) +
                         (u[j+1, i] - 2 * u[j, i] + u[j-1, i]) / (dy * dy))
                lap_v = ((v[j, i+1] - 2 * v[j, i] + v[j, i-1]) / (dx * dx) +
                         (v[j+1, i] - 2 * v[j, i] + v[j-1, i]) / (dy * dy))
                u_star[j, i] = u[j, i] + dt * (-conv_u + nu_eff[j, i] * lap_u)
                v_star[j, i] = v[j, i] + dt * (-conv_v + nu_eff[j, i] * lap_v)

    for j in range(ny):
        u_star[j, 0] = 1.0
        v_star[j, 0] = 0.0
        u_star[j, nx-1] = u_star[j, nx-2]
        v_star[j, nx-1] = v_star[j, nx-2]
    for i in range(nx):
        u_star[0, i] = 0.0
        v_star[0, i] = 0.0
        u_star[ny-1, i] = 0.0
        v_star[ny-1, i] = 0.0

    p = solve_poisson_pressure(p, u_star, v_star, dx, dy, dt)

    for j in range(1, ny-1):
        for i in range(1, nx-1):
            if inside_mask[j, i]:
                u[j, i] = 0.0
                v[j, i] = 0.0
            else:
                dpdx = (p[j, i+1] - p[j, i-1]) / (2 * dx)
                dpdy = (p[j+1, i] - p[j-1, i]) / (2 * dy)
                u[j, i] = u_star[j, i] - dt * dpdx
                v[j, i] = v_star[j, i] - dt * dpdy

    for j in range(ny):
        u[j, 0] = 1.0
        v[j, 0] = 0.0
        u[j, nx-1] = u[j, nx-2]
        v[j, nx-1] = v[j, nx-2]
    for i in range(nx):
        u[0, i] = 0.0
        v[0, i] = 0.0
        u[ny-1, i] = 0.0
        v[ny-1, i] = 0.0

    return u, v, p

def main_solver_les(steps=60000, Re=2000, Cs=0.1, dt=0.0002): 
    base_path = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
    data_path = os.path.join(base_path, "data")
    frames_path = os.path.join(data_path, "video_frames")
    os.makedirs(frames_path, exist_ok=True)

    XX = np.loadtxt(os.path.join(data_path, "XX.txt"))
    YY = np.loadtxt(os.path.join(data_path, "YY.txt"))
    polygon = np.loadtxt(os.path.join(data_path, "airfoil_points.txt"))
    if not np.allclose(polygon[0], polygon[-1]):
        polygon = np.vstack([polygon, polygon[0]])

    ny, nx = XX.shape
    dx = (XX[0, -1] - XX[0, 0]) / (nx - 1)
    dy = (YY[-1, 0] - YY[0, 0]) / (ny - 1)

    nu = 1.0 / Re
    inside_mask = build_airfoil_mask_numba(XX, YY, polygon)

    rng = np.random.default_rng(1234)
    u = np.ones((ny, nx), dtype=np.float32) + 0.01 * rng.random((ny, nx)).astype(np.float32)
    v = 0.01 * rng.random((ny, nx)).astype(np.float32)
    p = np.zeros((ny, nx), dtype=np.float32)

    start_time = time.time()
    for n in range(steps):
        u, v, p = step_fractional(u, v, p, inside_mask, dx, dy, dt, nu, Cs)
        if n % 500 == 0: 
            print(f"Krok {n}/{steps}")
            np.savetxt(os.path.join(frames_path, f"u_{n}.txt"), u)
            np.savetxt(os.path.join(frames_path, f"v_{n}.txt"), v)

    total_time = time.time() - start_time
    print(f"Solver calculation time: {total_time:.2f} seconds")

    np.savetxt(os.path.join(data_path, "u.txt"), u)
    np.savetxt(os.path.join(data_path, "v.txt"), v)
    np.savetxt(os.path.join(data_path, "p.txt"), p)
    np.savetxt(os.path.join(data_path, "inside_mask.txt"), inside_mask.astype(int))

    print("LES IB simulation completed. Results recorded.")

if __name__ == "__main__":
    main_solver_les(steps=60000, Re=2000, Cs=0.1, dt=0.0002)















