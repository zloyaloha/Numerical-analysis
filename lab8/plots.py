import numpy as np
import matplotlib
matplotlib.use('Agg') 
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from mpl_toolkits.mplot3d import Axes3D
import os

plt.style.use('seaborn-v0_8-darkgrid')

def parse_snapshot_file(filename):
    if not os.path.exists(filename):
        print(f"File not found: {filename}")
        return None, None, None, None

    print(f"Reading {filename}...")
    snapshots_X = []
    snapshots_Y = []
    snapshots_Z = []
    times = []

    current_data = []
    with open(filename, 'r') as f:
        lines = f.readlines()

    def process_block(data_block):
        if not data_block: return
        arr = np.array(data_block)
        x_raw = arr[:, 0]
        y_raw = arr[:, 1]
        z_raw = arr[:, 2]

        unique_x = np.unique(x_raw)
        unique_y = np.unique(y_raw)
        Nx = len(unique_x)
        Ny = len(unique_y)

        if len(x_raw) != Nx * Ny:
            print(f"Warning: Data block size {len(x_raw)} does not match grid {Nx}x{Ny}. Skipping frame.")
            return

        order = np.lexsort((y_raw, x_raw))
        X = x_raw[order].reshape((Nx, Ny))
        Y = y_raw[order].reshape((Nx, Ny))
        Z = z_raw[order].reshape((Nx, Ny))

        snapshots_X.append(X)
        snapshots_Y.append(Y)
        snapshots_Z.append(Z)

    for line in lines:
        line = line.strip()
        if not line: continue

        if line.startswith("TIME"):
            if current_data:
                process_block(current_data)
                current_data = []
            parts = line.split()
            if len(parts) > 1:
                times.append(float(parts[1]))
        else:
            try:
                vals = list(map(float, line.split()))
                if len(vals) == 3:
                    current_data.append(vals)
            except ValueError:
                continue

    if current_data:
        process_block(current_data)

    return snapshots_X, snapshots_Y, snapshots_Z, times

def parse_exact_file(filename):
    if not os.path.exists(filename): return None
    print(f"Reading {filename}...")
    return np.loadtxt(filename)

def parse_error_file(filename):
    if not os.path.exists(filename): return []
    print(f"Reading errors from {filename}...")
    with open(filename, 'r') as f:
        # Читаем все числа из файла, игнорируя переводы строк
        content = f.read()
        errors = list(map(float, content.split()))
    return errors

def save_animation(X_list, Y_list, Z_list, times, errors, title, filename_out):
    if not Z_list:
        print("No data to animate.")
        return

    print(f"Generating animation: {filename_out}...")
    fig = plt.figure(figsize=(16, 8))
    ax1 = fig.add_subplot(121, projection='3d')
    ax2 = fig.add_subplot(122)
    if errors and len(errors) > 0:
        ax2.plot(errors, label='Max Error', color='red')
        ax2.set_yscale('log')
        ax2.set_xlabel('Time Step')
        ax2.set_ylabel('Error')
        ax2.set_title('Error Dynamics')
        ax2.legend()
        ax2.grid(True)

    all_z = np.concatenate([z.flatten() for z in Z_list])
    z_min_global = np.min(all_z)
    z_max_global = np.max(all_z)
    if z_max_global > 1e5 or z_min_global < -1e5:
        print("Warning: Detected huge values (instability?). Clamping plot limits.")
        z_max_global = 2.0
        z_min_global = -0.5

    def update(frame):
        ax1.clear()
        ax1.set_title(f"{title}\nt = {times[frame]:.3f}")
        ax1.set_xlabel("X")
        ax1.set_ylabel("Y")
        ax1.set_zlabel("U")
        ax1.set_zlim(z_min_global, z_max_global)
        surf = ax1.plot_surface(X_list[frame], Y_list[frame], Z_list[frame],
                               cmap='viridis', edgecolor='none', rstride=1, cstride=1)
        return surf,

    ani = FuncAnimation(fig, update, frames=len(Z_list), interval=100, blit=False)
    ani.save(filename_out, writer='pillow', fps=10)
    plt.close(fig)
    print(f"Saved {filename_out}")

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