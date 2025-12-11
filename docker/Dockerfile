# Clean FEniCSx environment based on official dolfinx image
# This is the simplest and most reliable approach

FROM dolfinx/dolfinx:v0.7.2

# Install additional packages not in base image
RUN python3 -m pip install --no-cache-dir \
    jupyterlab \
    notebook \
    ipykernel \
    ipywidgets \
    gmsh \
    meshio \
    pyvista \
    matplotlib \
    seaborn \
    pandas \
    tqdm \
    rich

# Set working directory
WORKDIR /home/fenics/shared

# Expose Jupyter port
EXPOSE 8888

# Default command - start JupyterLab
CMD ["jupyter", "lab", "--ip=0.0.0.0", "--port=8888", "--no-browser", "--allow-root", "--NotebookApp.token=''"]