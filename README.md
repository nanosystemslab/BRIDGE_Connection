# BRIDGE Connection

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![Python](https://img.shields.io/badge/python-3.10-blue.svg)](https://www.python.org/downloads/)
[![FEniCSx](https://img.shields.io/badge/FEniCSx-0.7.2-orange.svg)](https://fenicsproject.org/)
[![Docker](https://img.shields.io/badge/docker-ready-brightgreen.svg)](./docker/)

## Overview

BRIDGE Connection is a computational mechanics project for simulating hyperelastic arch-block systems using FEniCSx/dolfinx. This project focuses on finite element analysis of bridge-like structures with contact mechanics.

## Features

- Hyperelastic material modeling
- Contact mechanics simulation
- Arch-block system analysis
- FEniCSx/dolfinx based FEM solver
- Jupyter notebook workflow
- Docker containerized environment

## Quick Start

### Using Docker (Recommended)

```bash
# Clone the repository
git clone https://github.com/nanosystemslab/BRIDGE_Connection.git
cd BRIDGE_Connection

# Start the Docker container
docker compose up -d

# Access Jupyter Lab at http://localhost:8888
# (No password required in development mode)
```

See [docker/README.md](./docker/README.md) for detailed Docker setup instructions.

### Local Installation

```bash
# Install from source
pip install -e .

# Or install dependencies only
pip install -r requirements.txt
```

## Project Structure

```
BRIDGE_Connection/
├── notebooks/active/     # Active development notebooks
├── src/BRIDGE_Connection/ # Python source code
├── docker/               # Docker configuration
├── data/                 # Input data files
├── results/              # Output and results
└── tests/                # Test suite
```

See [STRUCTURE.md](./STRUCTURE.md) for detailed project organization.

## Usage

### Running Simulations

1. Start Jupyter Lab using Docker:
   ```bash
   docker compose up -d
   ```

2. Open browser to http://localhost:8888

3. Navigate to `notebooks/active/` for current work

4. Open `working5.ipynb` for the main hyperelastic simulation

### Key Notebooks

- `working5.ipynb` - Main hyperelastic arch-block simulation
- Examples in `notebooks/examples/` (coming soon)

## Development

### Docker Environment

All dependencies are containerized for consistency:
- FEniCSx/dolfinx v0.7.2
- Scientific Python stack
- Jupyter Lab
- Mesh generation tools

### Contributing

See [CONTRIBUTING.md](./CONTRIBUTING.md) for guidelines.

### Testing

```bash
# Run tests in Docker
docker compose exec bridge-connection pytest

# Or locally
pytest tests/
```

## Documentation

- [Docker Setup](./docker/README.md) - Container configuration
- [Project Structure](./STRUCTURE.md) - Code organization
- [API Documentation](https://bridge-connection.readthedocs.io/) (coming soon)

## Team

Developed by the Nanosystems Lab team.

## License

This project is licensed under the GPL v3 License - see [LICENSE](./LICENSE) for details.

## Citation

If you use this software in your research, please cite:
```bibtex
@software{bridge_connection,
  title={BRIDGE Connection},
  author={Nanosystems Lab},
  year={2024},
  url={https://github.com/nanosystemslab/BRIDGE_Connection}
}
```

## Issues

Please report issues at: https://github.com/nanosystemslab/BRIDGE_Connection/issues

## Acknowledgments

- FEniCSx Project for the finite element framework
- ME-672 course for initial Docker environment inspiration