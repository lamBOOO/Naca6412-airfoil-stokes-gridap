# Simulation

![alt text](./result.png)

## Tutorial

1. Install all dependencies
```bash
julia --project -e 'import Pkg; Pkg.instantiate()'
```

2. Execute the solver
```bash
julia --project solve.jl
```

## Changes
- Renameded the boundaries in the `geo`-file
- add physical points
- Comments in geo: Physical names were wrong; airfoil switched with wall
