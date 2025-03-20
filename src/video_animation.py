import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import os
import time
from scipy.ndimage import gaussian_filter

def compute_vorticity(u, v, dx, dy):
    omega = (np.gradient(v, axis=1) / dx - np.gradient(u, axis=0) / dy)
    return omega

base_path = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
data_path = os.path.join(base_path, "data")
frames_path = os.path.join(data_path, "video_frames")
os.makedirs(frames_path, exist_ok=True)

XX = np.loadtxt(os.path.join(data_path, "XX.txt"))
YY = np.loadtxt(os.path.join(data_path, "YY.txt"))
dx = (XX[0, -1] - XX[0, 0]) / (XX.shape[1] - 1)
dy = (YY[-1, 0] - YY[0, 0]) / (YY.shape[0] - 1)

frame_interval = 500
num_frames = 120 
frame_numbers = [n * frame_interval for n in range(num_frames)]

def load_frame(n):
    u = np.loadtxt(os.path.join(frames_path, f"u_{n}.txt"))
    v = np.loadtxt(os.path.join(frames_path, f"v_{n}.txt"))
    omega = compute_vorticity(u, v, dx, dy)
    omega = gaussian_filter(omega, sigma=1.0)
    return omega

start_time = time.time()
frames = [load_frame(n) for n in frame_numbers]
load_time = time.time() - start_time
print(f"Frame loading time: {load_time:.2f} seconds")

fig, ax = plt.subplots(figsize=(10, 6))
levels = np.linspace(-15, 15, 61)

def update(i):
    ax.clear()
    ax.contourf(XX, YY, frames[i], levels=levels, cmap='RdBu_r', extend='both')
    ax.set_title(f"LES simulation - Step {frame_numbers[i]}")
    ax.axis("equal")

ani = animation.FuncAnimation(fig, update, frames=num_frames, blit=False)
ani.save(os.path.join(data_path, "simulation.mp4"), writer="ffmpeg", fps=24, dpi=300) 
print(f"Animation generation time: {time.time() - start_time - load_time:.2f} seconds")

plt.show()





