# IonBoost

`IonBoost` is a 1D Lagrangian solver for collisionless plasma expansion into a
vacuum with Boltzmann electrons — the mechanism behind TNSA-type ion
acceleration from short-pulse laser–solid interaction. It reproduces the model
and results of:

P. Mora, *Plasma Expansion into a Vacuum*, **Physical Review Letters** 90,
185002 (2003). DOI:
[10.1103/PhysRevLett.90.185002](https://doi.org/10.1103/PhysRevLett.90.185002).

The name comes from `IonBoost`, the original code written by Patrick Mora on
which this solver is based.

## Model

The code advances cold Lagrangian ion sheets coupled to an electrostatic field,
with electrons in isothermal Boltzmann equilibrium. In normalized units
(lengths in the initial Debye length, times in the inverse ion plasma
frequency, velocities in the ion-acoustic speed):

```text
∂t ni + ∂x(ni ui) = 0
∂t ui + ui ∂x ui = Zi E
∂x E = ni - ne
ne = n0hot exp(-φ / Th)
```

Beyond the ion front the field is continued analytically with the exact
Boltzmann vacuum-sheath solution.

## Numerical method

- ion sheets are advanced with a variable-step leapfrog scheme on an ordered
  1D mesh; the ion density is rebuilt from the sheet spacing
- the Poisson–Boltzmann equation for φ is solved by Newton iteration with a
  tridiagonal solve; the boundary condition at the ion front is the linearized
  `E = sqrt(2 pe)` relation (Eq. 8 of the paper)
- the vacuum sheath beyond the front uses the exact solution
  `φ = φ_f + 2 Th ln(1 + k x)` with `k = sqrt(n0 e^{-φ_f/Th} / 2 Th)`

## Validation

At `ωpi t = 50` (default notebook run) the solver matches the paper at the ~1%
level: the front field follows Eq. (9), the front velocity Eq. (10), the front
position `x/cs t = 5.59`, the front asymptotics Eqs. (14)/(17)/(18)/(20), and
the energy spectrum Eq. (21) with the Eq. (22) cutoff.

## Repository layout

- `src/`: Fortran sources and build files
- `sheath/`: sheath1d, a kinetic-ion Boltzmann-electron sheath solver
- `tests/input.in`: sample runtime input file
- `scripts/ionboost_plots.ipynb`: end-to-end notebook (build, run, and compare
  against the paper); figures are written to `scripts/viz/`
- `tools/`: notebooks for analysis / plotting

## Build

Use CMake from the repository root:

```bash
cmake -S src -B build
cmake --build build
```

This produces the executable `ionboost` in `build/`.

## Run

The program reads `input.in` from the current working directory:

```bash
cp tests/input.in build/input.in
cd build
./ionboost
```

Output files (CSV with a `#` header line):

- `conservation.txt`: time histories of particle numbers and energies
- `historique.txt`: time histories of front quantities
- `profil.txt`: final spatial profiles
- `spectres.txt`: final ion velocity and energy spectra

## Input file

`read_input()` expects `input.in` in this order:

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

See `tests/input.in` for an example. Note that the foil behaves as
semi-infinite (the regime of the 2003 paper) only while `tmax < lmax` in code
units; beyond that the run enters the thin-foil regime.
