"""
sheath1d — magnetized-sheath profile generator.

Kinetic ions (1D3V macroparticles, Boris push, oblique B) + Boltzmann
electrons closed by a nonlinear Poisson solve. Same solver family as the
Mora expansion code (Boltzmann closure + tridiagonal phi solve), with wall
and injection boundaries instead of vacuum expansion.

Normalized units:
  length  : lambda_D  (electron Debye length at the plasma entrance)
  potential: T_e / e   (phi <= 0 toward the wall)
  velocity: c_s = sqrt(T_e / m_i)
  time    : lambda_D / c_s   (inverse ion plasma frequency)

Model parameters:
  alpha_deg : angle between B and the WALL NORMAL (0 = normal incidence)
  tau       : T_i / T_e of the injected ions
  omega_hat : ion gyrofrequency in code units = Omega_i * lambda_D / c_s
              (= lambda_D / rho_s; 0 = unmagnetized)
  mu        : m_e / m_i (enters only the floating-wall condition)

Geometry: wall at s = 0, plasma entrance at s = L. B_hat = (cos a, 0, sin a).
Electrons: n_e = exp(phi), phi(L) = 0. Wall: floating (zero net current)
unless phi_wall is prescribed.

Output: steady-state phi(s), E(s), n_i(s), plus the wall IEAD sample.
"""

import numpy as np
from scipy.linalg import solve_banded


class Sheath1D:
    def __init__(self, alpha_deg=0.0, tau=1.0, omega_hat=0.0, mu=1.0/3672.0,
                 L=40.0, ncell=400, npart=200_000, u_inject=None,
                 phi_wall=None, seed=1, injection='chodura'):
        self.alpha = np.deg2rad(alpha_deg)
        self.tau = tau
        self.omega = omega_hat
        self.mu = mu
        self.L = L
        self.nc = ncell
        self.dx = L / ncell
        self.xe = np.linspace(0.0, L, ncell + 1)     # nodes
        self.rng = np.random.default_rng(seed)
        self.bhat = np.array([np.cos(self.alpha), 0.0, np.sin(self.alpha)])
        # injection drift along -B (toward the wall). Chodura-marginal by
        # default: parallel speed = sqrt(1 + tau) (generalized Bohm).
        self.u_inject = np.sqrt(1.0 + tau) if u_inject is None else u_inject
        self.injection = injection
        if injection == 'coulette':
            # chi(3) mean = 2*sqrt(2/pi) sigma -> <|v_par|> = 1.60 vth
            self.u_inject = 1.5958 * np.sqrt(tau)
        self.phi_wall_fixed = phi_wall
        self.phi_w = -3.0 if phi_wall is None else phi_wall
        self.phi = np.linspace(self.phi_w, 0.0, ncell + 1)
        # target macroparticle count sets the injection weight so that the
        # steady density at the entrance is ~1 in code units
        self.npart_target = npart
        self.x = np.empty(0)
        self.v = np.empty((0, 3))
        # flux of a drifting Maxwellian through the entrance (per unit area,
        # code units) -> particles injected per unit time
        self.wall_flux_acc = 0.0
        self.wall_time_acc = 0.0
        self.iead = []                                 # (E_impact, angle_deg)

    # ---- particle handling -------------------------------------------------

    def _inject(self, dt):
        """Inject ions at s=L from a drifting Maxwellian (drift -u along B)."""
        vth = np.sqrt(self.tau)
        # expected inward flux ~ n * u_drift for supersonic drift; sample the
        # inward-moving part of the shifted Maxwellian by rejection
        n_expect = self.flux0 * dt / self.w
        n_new = self.rng.poisson(n_expect)
        if n_new == 0:
            return
        v = np.empty((n_new, 3))
        # velocities in the (b, e1, e2) field-aligned frame
        if self.injection == 'coulette':
            # f ~ |v_par|^2 exp(-v^2/2vth^2): |v_par|/vth is chi(3)-distributed
            vpar = -vth * np.sqrt((self.rng.standard_normal((n_new, 3))**2).sum(1))
        else:
            vpar = -self.u_inject + vth * self.rng.standard_normal(n_new)
            vpar = np.where(vpar < 0.0, vpar, -np.abs(vpar))   # inward only
        vp1 = vth * self.rng.standard_normal(n_new)
        vp2 = vth * self.rng.standard_normal(n_new)
        e1 = np.array([-self.bhat[2], 0.0, self.bhat[0]])  # in-plane perp
        e2 = np.array([0.0, 1.0, 0.0])
        v = (vpar[:, None] * self.bhat + vp1[:, None] * e1 + vp2[:, None] * e2)
        # uniform placement across the entry the particle crosses in dt
        x = self.L + v[:, 0] * dt * self.rng.random(n_new)
        keep = x < self.L
        self.x = np.concatenate([self.x, x[keep]])
        self.v = np.concatenate([self.v, v[keep]])

    def _deposit(self):
        """CIC ion density on nodes, normalized to entrance density."""
        n = np.zeros(self.nc + 1)
        if len(self.x):
            f = self.x / self.dx
            i = np.clip(f.astype(int), 0, self.nc - 1)
            w = f - i
            np.add.at(n, i, 1.0 - w)
            np.add.at(n, i + 1, w)
        n *= self.w / self.dx
        n[0] *= 2.0; n[-1] *= 2.0                      # half-cells at ends
        return n

    def _efield(self):
        E = np.empty(self.nc + 1)
        E[1:-1] = -(self.phi[2:] - self.phi[:-2]) / (2 * self.dx)
        E[0] = -(self.phi[1] - self.phi[0]) / self.dx
        E[-1] = -(self.phi[-1] - self.phi[-2]) / self.dx
        return E

    def _push(self, dt):
        """Boris push with E = (Ex,0,0), B = omega * bhat (code units)."""
        if not len(self.x):
            return
        f = self.x / self.dx
        i = np.clip(f.astype(int), 0, self.nc - 1)
        w = f - i
        E = self._efield()
        Ex = (1.0 - w) * E[i] + w * E[i + 1]
        # half electric kick
        self.v[:, 0] += 0.5 * dt * Ex
        if self.omega > 0.0:
            t = self.bhat * np.tan(0.5 * dt * self.omega)
            s = 2.0 * t / (1.0 + t @ t)
            vp = self.v + np.cross(self.v, t)
            self.v = self.v + np.cross(vp, s)
        self.v[:, 0] += 0.5 * dt * Ex
        self.x += self.v[:, 0] * dt
        # wall absorption + IEAD tally
        hit = self.x <= 0.0
        if hit.any():
            vv = self.v[hit]
            ke = 0.5 * (vv ** 2).sum(1)                 # in T_e units (cs^2)
            ang = np.rad2deg(np.arccos(np.clip(-vv[:, 0] /
                       np.maximum(np.sqrt((vv ** 2).sum(1)), 1e-30), -1, 1)))
            self.iead.append(np.column_stack([ke, ang]))
            self.wall_flux_acc += hit.sum() * self.w
        gone = hit | (self.x >= self.L)
        self.x = self.x[~gone]
        self.v = self.v[~gone]

    # ---- field solve -------------------------------------------------------

    def _poisson(self, ni, newton_iters=30, tol=1e-10):
        """Solve phi'' = exp(phi) - ni with phi(L)=0, phi(0)=phi_w."""
        phi = self.phi.copy()
        phi[0] = self.phi_w
        phi[-1] = 0.0
        n = self.nc + 1
        for _ in range(newton_iters):
            r = np.zeros(n)
            r[1:-1] = ((phi[2:] - 2 * phi[1:-1] + phi[:-2]) / self.dx ** 2
                       - np.exp(phi[1:-1]) + ni[1:-1])
            ab = np.zeros((3, n))
            ab[0, 2:] = 1.0 / self.dx ** 2                     # upper
            ab[1, 1:-1] = -2.0 / self.dx ** 2 - np.exp(phi[1:-1])
            ab[1, 0] = ab[1, -1] = 1.0                          # Dirichlet
            ab[2, :-2] = 1.0 / self.dx ** 2                     # lower
            r[0] = r[-1] = 0.0
            d = solve_banded((1, 1), ab, -(-r))
            # solve J d = -r  (note residual sign convention above)
            phi[1:-1] -= d[1:-1]
            if np.max(np.abs(d)) < tol:
                break
        self.phi = phi

    def _update_floating_wall(self, relax=0.05):
        """Zero-net-current condition: Gamma_e(phi_w) = Gamma_i(wall)."""
        if self.phi_wall_fixed is not None or self.wall_time_acc <= 0.0:
            return
        gi = self.wall_flux_acc / self.wall_time_acc
        if gi <= 0.0:
            return
        # Boltzmann electron flux to wall: n_e(w) * sqrt(1/(2 pi mu)) with
        # n_e(w) = exp(phi_w) -> phi_w = ln( gi * sqrt(2 pi mu) )
        bx = max(abs(self.bhat[0]), 1e-6)
        target = np.log(gi * np.sqrt(2.0 * np.pi * self.mu) / bx)
        self.phi_w += relax * (target - self.phi_w)
        self.wall_flux_acc = 0.0
        self.wall_time_acc = 0.0

    # ---- driver ------------------------------------------------------------

    def run(self, t_end=None, dt=None, average_last=0.3, verbose=True):
        """March to steady state; returns dict with averaged profiles."""
        # entrance flux (code units): n0=1 times mean inward speed
        self.flux0 = self.u_inject
        # weight so steady content ~ npart_target: content ~ flux0 * L / u
        self.w = self.L / self.npart_target
        if dt is None:
            dt = 0.2 * self.dx / max(self.u_inject + 3 * np.sqrt(self.tau), 1.0)
            if self.omega > 0:
                dt = min(dt, 0.1 / self.omega)
        if t_end is None:
            t_end = 8.0 * self.L / min(self.u_inject, 1.0)
        nsteps = int(t_end / dt)
        n_avg_start = int(nsteps * (1.0 - average_last))
        acc_phi = np.zeros(self.nc + 1); acc_n = np.zeros(self.nc + 1); na = 0
        for step in range(nsteps):
            self._inject(dt)
            ni = self._deposit()
            self._poisson(ni)
            self._push(dt)
            self.wall_time_acc += dt
            if step % 50 == 0:
                self._update_floating_wall()
            if step >= n_avg_start:
                acc_phi += self.phi; acc_n += ni; na += 1
            if verbose and step % max(1, nsteps // 10) == 0:
                print(f"  step {step}/{nsteps}  N={len(self.x)}  "
                      f"phi_w={self.phi_w:+.3f}")
        phi = acc_phi / max(na, 1); ni = acc_n / max(na, 1)
        E = np.empty_like(phi)
        E[1:-1] = -(phi[2:] - phi[:-2]) / (2 * self.dx)
        E[0] = -(phi[1] - phi[0]) / self.dx; E[-1] = -(phi[-1] - phi[-2]) / self.dx
        iead = np.concatenate(self.iead) if self.iead else np.empty((0, 2))
        return dict(s=self.xe, phi=phi, E=E, ni=ni, ne=np.exp(phi),
                    phi_wall=self.phi_w, iead=iead,
                    params=dict(alpha_deg=np.rad2deg(self.alpha), tau=self.tau,
                                omega_hat=self.omega, mu=self.mu, L=self.L))


if __name__ == "__main__":
    # smoke: unmagnetized floating hydrogen sheath. Expect phi_wall ~ -2.8
    sh = Sheath1D(alpha_deg=0.0, tau=0.1, omega_hat=0.0, mu=1.0/1836.0,
                  L=30.0, ncell=300, npart=50_000)
    out = sh.run(verbose=True)
    print(f"phi_wall = {out['phi_wall']:.3f}  (H, tau=0.1; textbook ~ -2.8)")
    print(f"sheath e-fold width ~ few lambda_D: phi(5)={np.interp(5, out['s'], out['phi']):.3f}")
