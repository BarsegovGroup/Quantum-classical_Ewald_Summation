# Quantum–Classical Ewald Summation

A hybrid quantum–classical implementation of **Ewald summation** for the electrostatic energy of a periodic 3D lattice of point charges. The Ewald decomposition splits the total electrostatic energy into several terms; the long-range **reciprocal-space (Fourier) term** is evaluated on a quantum device using the **Quantum Fourier Transform (QFT)**, while the remaining terms are computed classically.

The notebook runs a complete, self-contained numerical experiment on a **64 × 64 × 64** grid (**18 qubits**) with a alternating-sign lattice of **30 752 point charges**, and validates the result against a direct O(N²) Coulomb sum.

---

## Reference

Mansur Ziiatdinov, Igor Novikov, Farid Ablayev, and Valeri Barsegov,
*Quantum algorithm for Ewald summation based computation of long-range electrostatics*,
arXiv:2512.20886 [quant-ph] (2025). <https://doi.org/10.48550/arXiv.2512.20886>

---

## Method

The total electrostatic energy is decomposed as

$E = E^L + E^S + E^{self} + E^{dipole}$

with the following components:

| Term | Where computed |
|------|----------------|
| $E^L$ | Reciprocal-space kernel applied to either a classical FFT (`fftn`) or the squared probability density read out of a quantum QFT circuit. |
| $E^S$ | Real-space sum with `erfc` screening, cut off at `RCUT`, using a 3D cell list and `n` periodic image shells. |
| $E^{self}$ | Closed-form self-interaction term (vectorized). |
| $E^{dipole}$ | Surface (dipole) correction, with surrounding-medium permittivity `epsp`. |
| $E_{Coulomb}$ | Direct `O(N²)` Coulomb sum over `m` image shells, used as ground truth. |

Here $\sigma$ is the Ewald split width, $V$ the supercell volume and $\mathbf{n}$ runs over periodic image shells. The full expressions for all four terms are given, with derivations, in the notebook itself.

### The quantum step

The charge distribution on the $M_x \times M_y \times M_z$ grid is encoded as the amplitudes of an $N_{qubits} = \log_2 M_x + \log_2 M_y + \log_2 M_z$ qubit state, with the flat index given by the bitstring $x \Vert y \Vert z$ ($x$ in the most significant bits). A **3D QFT** — one QFT per coordinate register, acting on disjoint qubits — is then applied.

The probability of measuring $\mathbf{k} = (k_x, k_y, k_z)$ is proportional to the squared structure factor, so $|\hat{\rho}(\mathbf{k})|^2$ is recovered from the measurement counts (rescaled by $N_{grid}$) and fed into the reciprocal-space Ewald kernel to give $E^L$.

The remaining terms ($E^S$ with a cell-list neighbour search, $E^{\mathrm{self}}$, $E^{\mathrm{dipole}}$) are classical and shared by both variants, so the quantum and classical totals differ **only** through $E^L$.

---

## What the notebook contains

| Section | Contents |
|---|---|
| Setup & notation | Grid, Ewald and sampling parameters; full glossary of every symbol used |
| Charge configuration | `place_box(...)` builds a centred, stride-2, alternating-sign (NaCl-like) lattice |
| Visualization | `plot_charges` (3D scatter of the supercell) and `plot_slice` (single lattice plane) |
| Quantum circuit | `initialize` → per-axis QFT (`X_REG`, `Y_REG`, `Z_REG`) → `measure_all` |
| Simulation | `AerSimulator` + `SamplerV2` sampling, or exact statevector probabilities |
| $E^L$ | Reciprocal-space Ewald kernel applied to either the measured density or a classical `fftn` |
| $E^S$ | Two-stage evaluation: cell-list neighbour build, then per-pair distance recomputation + `erfc` screening (mirroring an MD step) |
| $E^{\mathrm{self}}$, $E^{\mathrm{dipole}}$ | Vectorized closed-form terms |
| Ground truth | Direct, chunked O(N²) Coulomb sum over the unit cell plus image shells |
| Comparison | $E_{\mathrm{PME}}$ vs $E_{\mathrm{quantum}}$ vs $E_{\mathrm{Coulomb}}$, with relative RMS errors |
| Circuit complexity | Transpilation for IBM **FakeSherbrooke** (127 qubits); original vs compiled circuit depth |

---

## Requirements

Python **3.10+** and:

```
numpy
scipy
matplotlib
qiskit >= 2.0, < 3.0
qiskit-aer
qiskit-ibm-runtime
```

> The circuit uses `qiskit.circuit.library.QFT`, which is deprecated in Qiskit 2.1 and slated for removal in Qiskit 3.0 — hence the upper bound.

Install with:

```bash
pip install -r requirements.txt
```

or directly:

```bash
pip install numpy scipy matplotlib qiskit qiskit-aer qiskit-ibm-runtime
```

No IBM Quantum account or API token is needed — the notebook uses the local `AerSimulator` and the offline `FakeSherbrooke` backend model.

---

## Running

```bash
git clone https://github.com/IgorNovikov02/Quantum-Classical-Ewald-Summation.git
cd Quantum-Classical-Ewald-Summation
jupyter notebook Quantum-classical_Ewald-Summation.ipynb
```

Then run the cells top to bottom. The notebook is linear: every cell depends only on cells above it.

---

## Configuration

All knobs live in one cell near the top.

| Parameter | Default | Meaning |
|---|---|---|
| `EXACT_SIM` | `False` | `True` → exact statevector probabilities (no shot noise); `False` → sampling with `N_SHOTS` |
| `N_SHOTS` | `10_000_000` | Number of measurement samples |
| `MX`, `MY`, `MZ` | `64, 64, 64` | Grid points per axis; **each must be a power of 2** (⇒ 18 qubits) |
| `step` | `0.75e-10` m | Lattice spacing (box = 48 Å per side, $V = 1.1059\times10^{-25}\ \mathrm{m}^3$) |
| `SIGMA` | `0.880e-10` m | Ewald split width, balancing the real- and reciprocal-space contributions |
| `RCUT` | `5.0e-10` m | Real-space cutoff for $E^S$ |
| `epsp` | `1` | Relative permittivity of the surrounding medium (dipole term) |
| `n`, `m` | `1`, `1` | Periodic-image shells for $E^S$ and for the direct Coulomb reference |
| `CHUNK` | `100` | Row-block size for the direct Coulomb sum — affects peak memory only, not the result |

The charge configuration is set by `inp = place_box(32, 31, 31)` — 32 × 31 × 31 = 30 752 charges on the stride-2 sublattice, centred in the box, with signs alternating as $(-1)^{i_x+i_y+i_z}$. Any smaller box can be substituted by changing the three per-axis counts (each may run from 1 to 32); pass `centered=False` to anchor the box at the corner instead.

---

## Output

Running the notebook produces:

* the four Ewald energy terms and the totals $E_{\mathrm{PME}}$ (classical $E^L$) and $E_{\mathrm{quantum}}$ (measured $E^L$), in joules;
* the direct Coulomb reference energy and the relative RMS error of both variants against it;
* wall-clock timings for each term;
* the original and compiled circuit depth on IBM FakeSherbrooke;
* 3D and single-plane plots of the charge configuration.
