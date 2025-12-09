import numpy as np
import matplotlib.pyplot as plt

def main():
    # --- domain parameters ---
    lx = 11.0       # domain width
    ly = 2.0        # domain height
    Nx = 400        # resolution in x
    Ny = 200        # resolution in y

    # --- create grid ---
    XX, YY = np.meshgrid(
        np.linspace(0.0, lx, Nx + 1),
        np.linspace(0.0, ly, Ny + 1)
    )

    pi = np.pi

    # --- phi_mat formula (modify here as needed) ---
    phi_mat = -np.cos((4.0 * pi * (XX - 1.0))) * np.cos(4.0 * pi * YY) - 0.2 \
              + np.maximum(100.0 * (YY - ly + 0.05), 0.0)

    # --- plot ---
    plt.figure(figsize=(10, 4))
    plt.imshow(phi_mat, extent=[0, lx, 0, ly], origin='lower')
    plt.colorbar(label='phi value')
    plt.title("phi_mat Visualization")
    plt.xlabel("X")
    plt.ylabel("Y")
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    main()

