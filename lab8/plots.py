import numpy as np
import matplotlib
matplotlib.use('Agg') 
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from mpl_toolkits.mplot3d import Axes3D
import os

plt.style.use('seaborn-v0_8-darkgrid')

def parse_snapshot_file(filename):
    snapshots = []
    current_data = []
    times = []

    if not os.path.exists(filename):
        return None, None, None, None

    print(f"Reading {filename}...")
    with open(filename, 'r') as f:
        lines = f.readlines()

    for line in lines:
        line = line.strip()
        if not line:
            continue

        if line.startswith("TIME"):
            if current_data:
                snapshots.append(np.array(current_data))
                current_data = []
            parts = line.split()
            if len(parts) > 1:
                times.append(float(parts[1]))
        else:
            vals = list(map(float, line.split()))
            if len(vals) == 3:
                current_data.append(vals)

    if current_data:
        snapshots.append(np.array(current_data))

    X_list, Y_list, Z_list = [], [], []

    for snap in snapshots:
        x_raw = snap[:, 0]
        y_raw = snap[:, 1]
        z_raw = snap[:, 2]

        unique_x = np.unique(x_raw)
        unique_y = np.unique(y_raw)
        Nx = len(unique_x)
        Ny = len(unique_y)

        X = x_raw.reshape((Nx, Ny))
        Y = y_raw.reshape((Nx, Ny))
        Z = z_raw.reshape((Nx, Ny))

        X_list.append(X)
        Y_list.append(Y)
        Z_list.append(Z)

    return X_list, Y_list, Z_list, times

def parse_exact_file(filename):
    if not os.path.exists(filename): return None
    print(f"Reading {filename}...")
    return np.loadtxt(filename)

def parse_error_file(filename):
    if not os.path.exists(filename): return None
    print(f"Reading {filename}...")
    max_errors = []
    with open(filename, 'r') as f:
        line = f.readline()
        errors = list(map(float, line.split()))
    return errors

def save_animation(X_list, Y_list, Z_list, times, errors, title, filename_out):
    if not Z_list: return

    print(f"Generating animation: {filename_out} (please wait)...")
    fig = plt.figure(figsize=(20, 10))
    ax1 = fig.add_subplot(121, projection='3d')
    ax2 = fig.add_subplot(122)
    ax2.plot(errors, label='Max Error', color='red', linewidth=2)
    ax2.set_xlabel('Time Step')
    ax2.set_ylabel('Max Absolute Error')
    ax2.set_title(f'Error Dynamics: {title}')
    ax2.set_yscale('log')
    ax2.legend()
    ax2.grid(True, which="both", ls="--", alpha=0.7)

    z_min = min(np.min(z) for z in Z_list)
    z_max = max(np.max(z) for z in Z_list)

    def update(frame):
        ax1.clear()
        ax1.set_title(f"{title}\nTime: {times[frame]:.4f}")
        ax1.set_xlabel("X")
        ax1.set_ylabel("Y")
        ax1.set_zlabel("U")
        ax1.set_zlim(z_min, z_max)
        surf = ax1.plot_surface(X_list[frame], Y_list[frame], Z_list[frame],
                               cmap='viridis', edgecolor='none')
        return surf

    ani = FuncAnimation(fig, update, frames=len(Z_list), interval=200)
    ani.save(filename_out, writer='pillow', fps=5)
    plt.close(fig)
    print(f"Done! Saved to {filename_out}")

def save_error_plot(errors, method_name, filename_out):
    if not errors: return
    print(f"Saving error plot: {filename_out}...")
    plt.figure(figsize=(10, 6))
    plt.plot(errors, label=f'Max Error ({method_name})', color='red', linewidth=2)
    plt.xlabel('Time Step')
    plt.ylabel('Max Absolute Error')
    plt.title(f'Error Dynamics: {method_name}')
    plt.yscale('log')
    plt.legend()
    plt.grid(True, which="both", ls="--", alpha=0.7)
    plt.savefig(filename_out, dpi=300)
    plt.close()
    print("Done!")

def save_exact_heatmap(Z, filename_out):
    if Z is None: return
    print(f"Saving heatmap: {filename_out}...")
    plt.figure(figsize=(8, 6))
    plt.title("Exact Solution (Heatmap)")
    plt.imshow(Z.T, origin='lower', cmap='plasma', interpolation='bilinear')
    plt.colorbar(label='U')
    plt.xlabel('Index X')
    plt.ylabel('Index Y')
    plt.savefig(filename_out, dpi=300)
    plt.close()
    print("Done!")


def main():
    if os.path.exists("result_fract_steps.txt"):
        err = parse_error_file("fract_steps_error.txt")
        X, Y, Z, T = parse_snapshot_file("result_fract_steps.txt")
        save_animation(X, Y, Z, T, err, "Fract steps", "result_fract_steps.gif")

    if os.path.exists("result_adi.txt"):
        err = parse_error_file("adi_error.txt")
        X, Y, Z, T = parse_snapshot_file("result_adi.txt")
        save_animation(X, Y, Z, T, err, "ADI Scheme", "adi_anim.gif")

    exact_Z = parse_exact_file("result_exact.txt")
    if exact_Z is not None:
        save_exact_heatmap(exact_Z, "exact_solution.png")

if __name__ == "__main__":
    main()