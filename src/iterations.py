import numpy as np
import matplotlib.pyplot as plt

def generate_phi(Nx=150, Ny=75, lx=4.0, ly=1.0):
    XX, YY = np.meshgrid(np.linspace(0.0, lx, Nx+1), np.linspace(0.0, ly, Ny+1))
    YY -= 0.25  # Shift the pattern down by 0.25
    phi_mat = -np.cos(4.0 /lx* np.pi * XX-1.0) * np.cos(4.0 * np.pi * YY) - 0.4 \
              + np.maximum(100.0 * (XX - YY - lx + 0.1), 0.0)
    return XX, YY, phi_mat

def plot_phi():
    XX, YY, phi_mat = generate_phi()
    plt.figure(figsize=(8, 4))
    plt.contourf(XX, YY, phi_mat, levels=50, cmap='viridis')
    plt.colorbar(label='Phi Value')
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.title('Phi Function Plot')
    plt.show()

if __name__ == "__main__":
    plot_phi()

