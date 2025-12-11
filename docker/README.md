# Docker Setup for BRIDGE Connection

This directory contains all Docker-related files for the BRIDGE Connection project.

## Directory Structure

```
docker/
├── Dockerfile                # Main Dockerfile based on official dolfinx/dolfinx:v0.7.2
├── .dockerignore            # Files to exclude from Docker build context
└── scripts/                 # Helper scripts for Docker operations
    ├── save-image.sh       # Save Docker image to tar file
    ├── load-image.sh       # Load Docker image from tar file
    └── push-to-hub.sh      # Push image to Docker Hub
```

## Quick Start

From the project root directory:

```bash
# Start the container
docker compose up -d

# Access Jupyter Lab
# Open browser to: http://localhost:8888
# No password required (development setup)

# Stop the container
docker compose down
```

## Building the Image

If you need to rebuild the image:

```bash
# From project root
docker compose build

# Or directly with Docker
docker build -f docker/Dockerfile -t bridge-connection:dolfinx .
```

## Image Details

- **Base Image**: `dolfinx/dolfinx:v0.7.2`
- **Size**: ~5.7GB
- **Key Components**:
  - FEniCSx/dolfinx v0.7.2
  - JupyterLab & Jupyter Notebook
  - Scientific Python stack (numpy, scipy, matplotlib, pandas)
  - Mesh tools (gmsh, meshio, pyvista)
  - Development tools (ipykernel, ipywidgets, tqdm, rich)

## Helper Scripts

### Save Image to File
```bash
./docker/scripts/save-image.sh bridge-connection:dolfinx
```

### Load Image from File
```bash
./docker/scripts/load-image.sh docker/images/bridge-connection.tar
```

### Push to Docker Hub
```bash
# Set your Docker Hub username
export DOCKER_HUB_USER=nanosystemslab
./docker/scripts/push-to-hub.sh bridge-connection:dolfinx
```

## Development Notes

- Authentication is disabled for easier development access
- Mount point: `/home/fenics/shared` maps to project root
- Jupyter config persists in named volume `jupyter_config`
- Container runs as root (development mode)

## Production Considerations

For production deployment:
1. Enable authentication (add token/password)
2. Run as non-root user
3. Use specific version tags instead of :latest
4. Consider using docker secrets for sensitive data