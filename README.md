# Mora

`Mora` is a small Fortran solver for 1D electrostatic plasma expansion into vacuum. This repository appears to reproduce or adapt the model introduced in:

P. Mora, *Plasma Expansion into a Vacuum*, **Physical Review Letters** 90, 185002 (published May 7, 2003). DOI: [10.1103/PhysRevLett.90.185002](https://doi.org/10.1103/PhysRevLett.90.185002).

Note: the DOI above is a **PRL** paper, not PRE.

## Model

The current implementation advances a 1D ion expansion problem coupled to an electrostatic field. In normalized form, the code structure corresponds to:

```text
∂t ni + ∂x(ni ui) = 0
∂t ui + ui ∂x ui = Zi E
∂x E = ni - ne
E = -∂x φ
```

with a Boltzmann-like hot-electron closure used in the present source:

```text
ne = n0hot exp(-φ / Th)
```

The vacuum-side sheath is then extended beyond the ion front using the same electrostatic / Boltzmann closure used in the solver.

## Numerical method

From the current source layout:

- ion positions and velocities are advanced in time on a 1D ordered mesh / sheet representation
- ion density is rebuilt from the updated cell spacing
- the electrostatic potential is obtained from a tridiagonal linear solve
- the tridiagonal system is inverted by `gauss()` in `src/solver.f90`
- the solver iterates on `φ`, and optionally on `n0hot` and `Th`, to enforce electron-number and energy constraints

## Repository layout

- `src/`: Fortran sources and build files
- `tests/input.in`: sample runtime input file
- `test_ding_2015.py`: small parameter-check script
- `scripts/viz/`: generated figures
- `tools/`: notebooks for analysis / plotting

## Build

Use CMake from the repository root:

```bash
cmake -S src -B build
cmake --build build
```

This is intended to produce the executable `mora` in `build/`.

## Run

The program reads `input.in` from the current working directory. One simple way to run the sample case is:

```bash
cp tests/input.in build/input.in
cd build
./mora
```

Expected output files include:

- `conservation.txt`
- `historique.txt`
- `profil.txt`
- `spectres.txt`

## Input file

The current `read_input()` implementation expects `input.in` in this order:

```text
ncell nvacuum prog
itmax tmax
iter0 iter1 iter2 iter3
dti
n0hot
Thmax
lfini lmax
nb_cons En_cons
T_MeV LSS nLSS
```

See `tests/input.in` for an example.

## Current status

The repository has been cleaned up so the documentation, ignore rules, and sample input are more consistent. The current code path is now treated as a simplified single-species Mora solver with fixed hot-electron temperature, machine-readable outputs, and an integrated plotting notebook.

## Suggested wording

If you want a short repository description, this is the cleanest version:

> This code reproduces or adapts the 1D plasma-expansion model introduced in P. Mora, *Plasma Expansion into a Vacuum*, Phys. Rev. Lett. 90, 185002 (2003), DOI: 10.1103/PhysRevLett.90.185002.
