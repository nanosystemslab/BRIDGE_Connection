# BRIDGE Connection

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![Python](https://img.shields.io/badge/python-3.10-blue.svg)](https://www.python.org/downloads/)
[![FEniCSx](https://img.shields.io/badge/FEniCSx-0.7.2-orange.svg)](https://fenicsproject.org/)
[![Docker](https://img.shields.io/badge/docker-ready-brightgreen.svg)](https://hub.docker.com/r/nanosystemslab/bridge-connection)

## Overview

BRIDGE Connection is a computational mechanics research project focused on simulating the structural behavior of arch-block bridge systems using finite element methods (FEM). The project implements hyperelastic material models and contact mechanics to analyze the interaction between bridge components under various loading conditions.

## Key Features

- **Hyperelastic Material Modeling**: Implementation of nonlinear material behavior for accurate structural response
- **Contact Mechanics**: Simulation of contact interactions between arch and block components
- **FEniCSx/dolfinx Integration**: Leverages state-of-the-art FEM solver for efficient computation
- **Parametric Studies**: Framework for exploring different geometric and material configurations
- **Visualization Tools**: Integrated post-processing and visualization of results
- **Docker Environment**: Containerized setup for reproducible research

## Quick Start

### Prerequisites
- Docker Desktop installed
- Git for version control
- 8GB+ RAM recommended for simulations

### Installation

```bash
# Clone the repository
git clone https://github.com/nanosystemslab/BRIDGE_Connection.git
cd BRIDGE_Connection

# Pull the Docker image from Docker Hub
docker pull nanosystemslab/bridge-connection:latest

# Start the Docker container
docker compose up -d

# Access Jupyter Lab at http://localhost:8888
# (No password required in development mode)
```

## Project Structure

```
BRIDGE_Connection/
├── notebooks/
│   ├── active/           # Current research notebooks
│   │   └── working5.ipynb    # Main hyperelastic arch-block simulation
│   ├── examples/         # Example simulations and tutorials
│   └── experiments/      # Experimental notebooks and tests
├── src/
│   └── BRIDGE_Connection/   # Python package with core functionality
│       ├── materials/       # Material models (hyperelastic, etc.)
│       ├── geometry/        # Geometry generation and mesh utilities
│       ├── solvers/         # FEM solver configurations
│       └── visualization/   # Post-processing and plotting
├── data/                    # Input data and parameters
│   ├── meshes/             # Pre-generated mesh files
│   └── parameters/         # Configuration files
├── results/                 # Simulation outputs
│   ├── vtk/               # VTK files for ParaView
│   ├── figures/           # Generated plots and images
│   └── logs/              # Simulation logs
├── docker/                  # Docker configuration
└── tests/                   # Unit and integration tests
```

## Usage

### Running the Main Simulation

1. Open Jupyter Lab: http://localhost:8888
2. Navigate to `notebooks/active/working5.ipynb`
3. Run cells sequentially to:
   - Set up geometry parameters
   - Generate the arch-block mesh
   - Apply boundary conditions
   - Solve the hyperelastic problem
   - Visualize results

### Key Parameters

The simulation can be customized by modifying:
- **Geometry**: Arch radius, block dimensions, gap sizes
- **Materials**: Young's modulus, Poisson's ratio, hyperelastic parameters
- **Loading**: Applied forces, displacement boundary conditions
- **Solver**: Newton iteration parameters, convergence criteria

### Example Code

```python
import dolfinx
from BRIDGE_Connection import ArchBlockSystem

# Create system
system = ArchBlockSystem(
    arch_radius=10.0,
    block_width=2.0,
    material_E=1000.0,
    material_nu=0.3
)

# Apply loads and solve
system.apply_gravity()
system.solve()

# Visualize results
system.plot_displacement()
system.save_results("results/simulation_001")
```

## Development

### Building from Source

```bash
# Install in development mode
pip install -e .

# Run tests
pytest tests/

# Run specific notebook
jupyter lab notebooks/active/working5.ipynb
```

### Contributing

1. Fork the repository
2. Create a feature branch (`git checkout -b feature/amazing-feature`)
3. Commit changes (`git commit -m 'Add amazing feature'`)
4. Push to branch (`git push origin feature/amazing-feature`)
5. Open a Pull Request

See [CONTRIBUTING.md](./CONTRIBUTING.md) for detailed guidelines.

## Docker Environment

The project includes a fully configured Docker environment with:
- FEniCSx/dolfinx v0.7.2
- Scientific Python stack (numpy, scipy, matplotlib)
- Jupyter Lab for interactive development
- Mesh generation tools (gmsh, meshio)
- Visualization libraries (pyvista, paraview)

Docker image available at: [nanosystemslab/bridge-connection](https://hub.docker.com/r/nanosystemslab/bridge-connection)

## Research Applications

This project supports research in:
- Structural analysis of masonry arch bridges
- Contact mechanics in multi-body systems
- Hyperelastic material behavior under large deformations
- Seismic response of arch structures
- Optimization of arch-block configurations

## Publications

If you use this software in your research, please cite:
```bibtex
@software{bridge_connection_2024,
  title={BRIDGE Connection: FEM Analysis of Arch-Block Bridge Systems},
  author={Nanosystems Lab},
  year={2024},
  url={https://github.com/nanosystemslab/BRIDGE_Connection}
}
```

## Team

Developed by the [Nanosystems Lab](https://github.com/nanosystemslab) research group.

## License

This project is licensed under the GNU General Public License v3.0 - see [LICENSE](./LICENSE) for details.

## Support

- **Issues**: [GitHub Issues](https://github.com/nanosystemslab/BRIDGE_Connection/issues)
- **Discussions**: [GitHub Discussions](https://github.com/nanosystemslab/BRIDGE_Connection/discussions)

## Acknowledgments

- FEniCSx Project for the finite element framework
- DOLFINx developers for the computational backend

## Status

🚧 **Active Development** - This project is under active research and development. Features and APIs may change.

---

For detailed documentation, visit our [Wiki](https://github.com/nanosystemslab/BRIDGE_Connection/wiki)