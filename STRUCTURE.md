# Project Structure Guide

## Directory Organization

### 📓 `notebooks/`
- **`active/`** - Current working notebooks (e.g., working5.ipynb)
- **`experiments/`** - Experimental notebooks for testing new ideas
- **`examples/`** - Clean, documented example notebooks for reference

### 📚 `src/BRIDGE_Connection/`
Function library for reusable code. Import in notebooks as:
```python
from src.BRIDGE_Connection import geometry, materials, visualization
```

Planned modules:
- `geometry.py` - Mesh and geometry creation functions
- `materials.py` - Material models (hyperelastic, etc.)
- `visualization.py` - Plotting and VTK output utilities
- `solvers.py` - FEM solver configurations
- `boundary.py` - Boundary condition utilities
- `io_utils.py` - File I/O utilities

### 📁 `data/`
- **`meshes/`** - Saved mesh files (.msh, .geo)
- **`parameters/`** - Configuration files, parameter sets

### 📊 `results/`
- **`vtk/`** - VTK output files (.vtu, .pvd, .xdmf)
- **`figures/`** - Generated plots and visualizations
- **`logs/`** - Computation logs and solver output

### 🧪 `tests/`
Unit tests for library functions

### 📖 `docs/`
Documentation files

## Workflow

1. **Development**: Work in `notebooks/active/`
2. **Experiments**: Test new ideas in `notebooks/experiments/`
3. **Refactoring**: Extract common functions to `src/BRIDGE_Connection/`
4. **Results**: Save outputs to appropriate `results/` subdirectory
5. **Version Control**: The `.gitignore` is configured to exclude generated files

## Git Branches

- `main` - Main development branch
- `working-basic-boundary-condition` - Basic working implementation
- `modified-boundary-condition` - Modified boundary conditions
- `current-development` - Active development branch

## Notes

- The `merge/` directory is intentionally excluded from version control
- Generated files (.vtu, .msh, etc.) are automatically ignored by git
- Use .gitkeep files to maintain empty directory structure in git