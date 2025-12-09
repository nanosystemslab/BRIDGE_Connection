import numpy as np
import matplotlib.pyplot as plt
from dolfin import *
from compliance import init  # Ensure compliance.py is in the same directory

# Run the initialization for 'cantilever_asymmetric'
Lag, Nx, Ny, lx, ly, XX, YY, phi_mat, Load, bcd = init("cantilever_asymmetric")

# Plot the first iteration of the material distribution
plt.figure(figsize=(8, 4))
plt.contourf(XX, YY, phi_mat, levels=100, cmap='viridis')
plt.colorbar(label="Level Set Function (φ)")
plt.xlabel("X-coordinate")
plt.ylabel("Y-coordinate")
plt.title("First Iteration: Cantilever Asymmetric")
plt.show()

