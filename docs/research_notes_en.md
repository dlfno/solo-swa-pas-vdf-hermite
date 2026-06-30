# Notas de investigacion (fuentes por tema)


# ===== swa-pas-data-product =====

I have gathered comprehensive primary-source data from the ESA SWA-PAS L2 Data User Guide (Fedorov 2020), the Owen et al. 2020 A&A instrument paper, and the online readthedocs documentation. Below are the complete notes.

---

# Solar Orbiter SWA‑PAS Level‑2 Velocity Distribution Function (VDF) — Technical Notes

## 0. Product identity and provenance

- **Mission / instrument:** Solar Orbiter, *Solar Wind Analyser* (SWA) suite, *Proton–Alpha Sensor* (PAS).
- **CDF dataset:** `solo_L2_swa-pas-vdf` (the 3‑D VDF product). Sister L2 products are `swa-pas-eflux` (1‑D omni‑directional differential energy flux) and `swa-pas-grnd-mom` (ground‑computed moments).
- **Primary documentation:** *SWA‑PAS L2 Data User Guide*, A. Fedorov (IRAP), V02, 16 Nov 2020 (ESA COSMOS). This is the authoritative source for variable names, units, frame conventions, and the velocity‑reconstruction formulae quoted below.
- **Instrument paper:** Owen, C. J., et al. 2020, "The Solar Orbiter Solar Wind Analyser (SWA) suite," *Astronomy & Astrophysics*, **642**, A16. DOI: `10.1051/0004-6361/201937259`.
- **Important citation correction:** there is **no** standalone "Livi et al. 2023 PAS in‑flight calibration" paper. *Livi, S., et al. 2023, A&A, 676, A36* (DOI `10.1051/0004-6361/202346304`) is **"First results from the Solar Orbiter Heavy Ion Sensor" (SWA‑HIS)** — a different SWA sensor. PAS in‑flight characterisation/operations are documented in **Owen et al. 2020** (build & calibration, Table 8), the **Fedorov L2 User Guide (2020)**, and **Louarn, P., et al. 2024, A&A, 682, A44** (SWA‑PAS solar‑wind moments / cross‑calibration; DOI `10.1051/0004-6361/202347874`). Cite those rather than a non‑existent Livi‑PAS‑2023.

---

## 1. The SWA‑PAS instrument: top‑hat electrostatic analyser

**Concept.** All three SWA sensors (EAS, PAS, HIS) are **top‑hat electrostatic analysers** (Owen et al. 2020, §1, §3.3). PAS is a quarter/partial‑spherical top‑hat ESA optimised for the narrow, strongly radial solar‑wind ion beam. It measures the full 3‑D VDF of solar‑wind ions **without mass or charge selection** (in practice protons H⁺ and alphas He²⁺) over **200 eV/q – 20 keV/q**, mounted behind a dedicated scalloped cut‑out in the spacecraft heat shield so that it looks roughly sunward (Owen et al. 2020, §3.3.1–3.3.2).

**Three measurement axes (how it forms the 3‑D matrix):**

1. **Energy / charge (E/q) — the ESA itself.** A spherical‑section electrostatic analyser transmits only ions in a narrow E/q passband set by the voltage on the inner hemisphere. The **analyser constant k ≈ 13–14 eV/V** (Owen Table 8). The energy index of a bin is fixed by the instantaneous analyser voltage. PAS steps through **96 logarithmically spaced energy levels** (92 used in a given sample) covering 200 eV/q–20 keV/q. **Relative energy resolution ΔE/E ≈ 5.5 %** (design), measured **3.0–9.3 %** (Owen Table 8).

2. **Elevation (polar) — entrance deflection system.** Two curved **deflector plates** plus a **"top‑cap"** electrode electrostatically steer incoming ions toward the ESA entrance. The **deflector‑to‑analyser voltage ratio (U_def/U_an)** selects the **elevation index**. PAS sweeps elevation continuously for each energy; the instantaneous elevation response width is ≈ 3°, binned into **9 elevation bins** spanning **−22.5° to +22.5°** (build cal: −20° to +23°), ~5° resolution (Owen §3.3.3, Table 8).

3. **Azimuth — anode/CEM array.** The azimuthal focal position on the detector encodes the **azimuth index**. PAS uses an array of **11 channel electron multipliers (CEMs / "channeltrons")**, each one an azimuth ("CEM") bin, all read out **simultaneously**. Azimuth coverage is **−24° to +42°** (a 66°‑wide fan, offset +9° from the Sun direction to capture solar‑wind aberration), ~5° resolution.

**Sweep scheme (Owen §3.3.3, Fedorov §2.2):** at fixed energy, PAS sweeps elevation while accumulating all 11 azimuth bins at once; on completing the elevation sweep it steps energy. A full **(96 energy × 9 elevation × 11 azimuth)** matrix is acquired in **~1 s**; basic accumulation per element ≈ **1 ms** (0.95 ms). Geometric factor per bin ≈ **5 × 10⁻⁶ cm² sr eV/eV** (sensitivity ≥ 4 × 10⁻⁶, build 4–6 × 10⁻⁶; Owen Table 8). A **peak‑tracking** algorithm recentres a reduced energy–elevation window on the beam to enable fast modes.

**Cadence / modes (Fedorov §1):**
- **Normal mode:** one full VDF every **4 s**.
- **Snapshot:** bursts of ~9 s every 300 s at K samplings/s (typ. 4 Hz).
- **Burst mode:** 300 s of continuous high‑cadence sampling, **up to ~15–20 Hz** (phase space reduced, e.g. 24 energies × 5 deflections).

---

## 2. The L2 `swa-pas-vdf` CDF product — variables, dimensions, units

From the L2 VDF CDF file (Fedorov 2020, Table 1; confirmed against the online product). **CDF variable names in `code` font.**

### Core data variable
| Variable | Dim. | Type | Units | Notes |
|---|---|---|---|---|
| `vdf` | **[time, 11, 9, 96]** | CDF_REAL4 | **`s^3 m^-6`** | Phase‑space density (the "ion number density in phase space"). `DEPEND_1=Azimuth`, `DEPEND_2=Elevation`, `DEPEND_3=Energy`. |

**Dimension ordering (critical):** the per‑record VDF matrix is **[azimuth = 11, elevation = 9, energy = 96]**; with the time record it is **time × azimuth × elevation × energy**. (Owen §3.3.3 writes it as a (96, 9, 11) = energy×elevation×azimuth matrix conceptually; the **CDF storage order is azimuth, elevation, energy** as given by `DEPEND_1/2/3`. Watch this — Owen's prose order is the reverse of the CDF axis order.)

**Units of `vdf`:** **`s^3 m^-6`** = s³·m⁻⁶ (phase‑space density f, such that n = ∫ f d³v with d³v in (m/s)³ giving m⁻³). This is the SI phase‑space density; to convert to the common CGS plasma unit s³·cm⁻⁶ multiply by 10⁻⁶.

### Coordinate / support variables
| Variable | Dim. | Units | Meaning |
|---|---|---|---|
| `Epoch` | [time] | ns (TT2000) | Record time, CDF_TIME_TT2000. |
| `Half_interval` | [time] | s | Acquisition half‑interval. |
| `SCET` | [time] | ticks | Onboard clock. |
| `Info` | [time] | — | Sampling category: 0 Ground, 1 Normal, 2 Snapshot, 3 Burst, 4 Engineering, 5 Calibration. **Use only Info = 1/2/3.** |
| `Energy` | [96] | **eV** | Centre of energy bins. `DELTA_PLUS_VAR=delta_p_Energy`, `DELTA_MINUS_VAR=delta_m_Energy`. |
| `delta_p_Energy` | [96] | eV | Upper energy half‑width. |
| `delta_m_Energy` | [96] | eV | Lower energy half‑width. |
| `Azimuth` | [11] | **deg** | Centre of CEM (azimuth) bins. `DELTA_PLUS_VAR=DELTA_MINUS_VAR=delta_Azimuth`. |
| `delta_Azimuth` | [11] | deg | Azimuth half‑width per CEM. |
| `Elevation` | [9] | **deg** | Centre of elevation bins. `DELTA_PLUS_VAR=DELTA_MINUS_VAR=delta_Elevation`. |
| `delta_Elevation` | [9] | deg | Elevation half‑width. |
| `Full_azimuth` | [11, 9] | deg | Full azimuth table az(iaz, iel) — azimuth depends on elevation. |
| `Full_elevation` | [11, 9] | deg | Full elevation table el(iaz, iel). |
| `Elevation_correction` | [96] | deg | Per‑energy elevation offset (the "saw‑tooth" bin‑boundary correction; sign alternates with energy parity). |
| `PAS_to_RTN` | **[time, 3, 3]** | None | Instrument→RTN rotation matrix (see §4). |
| `start_energy`,`nb_energy` | [time] | — | First energy index and count actually swept (window may be a subset, e.g. energies 10:73). |
| `start_elevation`,`nb_elevation` | [time] | — | Elevation window start/count (e.g. 3:7). |
| `start_CEM`,`nb_CEM` | [time] | — | Azimuth (CEM) window start/count. |
| `nb_K`,`K` | [time] | — | Sub‑samplings/s and current sub‑index (snapshot/burst). |

> Because peak‑tracking moves the active window, a record may populate only a sub‑block of the 11×9×96 array; the `start_*`/`nb_*` variables tell you which indices are valid.

---

## 3. Reconstructing 3‑D velocity from (E, elevation, azimuth)

### 3.1 Speed from energy

The ESA selects energy‑per‑charge; the `Energy` variable is in **eV**. The speed magnitude is

$$|\mathbf v| \;=\; \sqrt{\frac{2qE}{m}}\;=\;\sqrt{\frac{2E_{\rm eV}\,e}{m}}.$$

For **protons** (m = m_p, q = e) this reduces to the user‑guide constant **`E2V = 13.85`**:

$$|\mathbf v|\,[\mathrm{km/s}] \;=\; 13.85\,\sqrt{E_{\rm eV}}\qquad\big(\sqrt{2e/m_p}=1.385\times10^{4}\ \mathrm{m\,s^{-1}\,eV^{-1/2}}\big).$$

Per‑bin lower/upper speeds use the energy half‑widths (Fedorov §4):

$$v_{\min}[ie]=13.85\sqrt{E[ie]-\texttt{delta\_m\_Energy}},\qquad v_{\max}[ie]=13.85\sqrt{E[ie]+\texttt{delta\_p\_Energy}}.$$

(For **alphas**, same E/q step ⇒ q = 2e, m = 4m_p ⇒ |v|_α = |v|_p/√2; PAS does not mass‑separate, so the default conversion assumes protons.)

### 3.2 Spherical → Cartesian: the two conventions

**(a) "Lecture‑note" / standard‑spherical convention (all‑positive).** Treating azimuth `az` and elevation `el` as ordinary spherical angles gives the **look/arrival direction** unit vector:

$$\hat L_x=\cos(el)\cos(az),\quad \hat L_y=\sin(az)\cos(el),\quad \hat L_z=\sin(el).$$

This is the direction **from which** the ion arrived (where the instrument points). It is **not** the particle velocity.

**(b) Community / SWA convention (with minus signs).** Because the instrument measures the **look (arrival) direction**, the **particle velocity is the negative of the look direction**, $\mathbf v = -|\mathbf v|\,\hat L$. The SWA‑PAS User Guide folds this minus sign directly into the unit vectors. The **ion‑velocity unit vector** for a bin is (Fedorov 2020, §2.1):

**PAS analyser frame:**
$$\hat V^{\rm PAS}_X=-\cos(El)\cos(Az),\qquad \hat V^{\rm PAS}_Y=-\cos(El)\sin(Az),\qquad \hat V^{\rm PAS}_Z=+\sin(El).$$

**Spacecraft (SRF) frame:**
$$\boxed{\;\hat V^{\rm SRF}_X=-\cos(El)\cos(Az),\qquad \hat V^{\rm SRF}_Y=+\cos(El)\sin(Az),\qquad \hat V^{\rm SRF}_Z=-\sin(El)\;}$$

(PAS→SRF is a 180° rotation about the sun‑pointing X axis: X kept, Y and Z flipped.) The full velocity vector is

$$\mathbf v = |\mathbf v|\,\hat V = 13.85\sqrt{E_{\rm eV}}\;\big(\hat V_X,\hat V_Y,\hat V_Z\big)\ \ [\mathrm{km/s}].$$

**Which gives a positive radial bulk speed?** The **SWA minus‑sign convention** does. For a nominal beam (Az≈0, El≈0) it yields $\hat V^{\rm SRF}\approx(-1,0,0)$ — i.e. the ion velocity points **anti‑sunward** (−X_SRF), as physically required (solar wind flows away from the Sun). Combined with the SRF↔RTN relation $X_{\rm RTN}\!\approx\!-X_{\rm SC},\,Y_{\rm RTN}\!\approx\!-Y_{\rm SC}$ (Fedorov §2.1) and the per‑record `PAS_to_RTN`, the radial component comes out **V_R ≈ +|v| > 0** (e.g. +400 km/s). The all‑positive lecture‑note formula instead returns the **arrival direction** (≈ +X_SRF, sunward) → a spurious **negative** radial speed; the minus signs are mandatory to recover the outward‑flowing solar wind.

### 3.3 Per‑bin corner construction (for accurate moment cells)

The User Guide builds each velocity‑space cell from the **corner** angles `Full_azimuth`/`Full_elevation` (= az/el min, centre, max per (iaz, iel)) and the per‑energy `Elevation_correction` (Fedorov §4):

```
x_i = -cos(elArr_i) * cos(azArr_i)
y_i =  cos(elArr_i) * sin(azArr_i)
z_i = -sin(elArr_i)
if ((ie - se) % 2 == 0):  z_i += Elevation_correction[ie]
else:                     z_i -= Elevation_correction[ie]
P_i = [x_i, y_i, z_i]      # then scale by v_min/v_max for the cell corners
```

The parity‑dependent `Elevation_correction` produces the characteristic "saw‑tooth" bin boundaries in the V_x–V_z plane.

---

## 4. The `PAS_to_RTN` 3×3 rotation matrix

- **Variable:** `PAS_to_RTN`, dim **[time, 3, 3]**, CDF_REAL8, `UNITS = None`.
- **Metadata:** `COORDINATE_SYSTEM = SOLO_SWA_PAS`, `TARGET_SYSTEM = SOLO_SUN_RTN`. It maps the **PAS instrument xyz frame → heliocentric RTN** (R radially outward from Sun, T ≈ along orbital motion, N completing the triad). It is time‑dependent because the spacecraft attitude relative to RTN varies along the orbit.
- **Application:** build the velocity vector in the PAS frame, then rotate. With M = `PAS_to_RTN[t]` (rows = R, T, N components):

$$\mathbf v_{\rm RTN} = M\;\mathbf v_{\rm PAS},\qquad \begin{pmatrix}v_R\\v_T\\v_N\end{pmatrix}=M\begin{pmatrix}v_X^{\rm PAS}\\v_Y^{\rm PAS}\\v_Z^{\rm PAS}\end{pmatrix}.$$

  Apply per time record. M is a proper rotation (orthonormal, det = +1), so the inverse is its transpose (RTN→PAS = Mᵀ). The ground‑moment product separately supplies `V_SRF` and `V_RTN` (both km/s) and `P_SRF`/`P_RTN`, `TxTyTz_SRF`/`TxTyTz_RTN` for cross‑checking your own rotation. Note the SRF↔RTN sign relation $X_{\rm RTN}\!\approx\!-X_{\rm SRF}$, $Y_{\rm RTN}\!\approx\!-Y_{\rm SRF}$ holds only approximately ("most of the time"); use `PAS_to_RTN` for the exact, attitude‑correct transform.

---

## 5. Energy table, angular bins, and DELTA bin‑widths (moment integration)

**Energy table.** **96 channels** (92 used per sample), **logarithmically (quasi‑geometric) spaced** from 200 eV/q to 20 keV/q (`Energy` in eV, `VALIDMAX = 40000`). Because the spacing is log, each channel has **asymmetric half‑widths** `delta_p_Energy` and `delta_m_Energy` (both [96], eV) — they are **not** equal, which is why both DELTA_PLUS_VAR and DELTA_MINUS_VAR are provided.

**Angular bins.** 11 azimuth (CEM) bins over −24°…+42°; 9 elevation bins over ≈ −22.5°…+22.5°; nominal ~5–6° spacing. Bin half‑widths exist:
- `delta_Azimuth` [11] (DELTA_PLUS = DELTA_MINUS, symmetric),
- `delta_Elevation` [9] (symmetric),
- `Full_azimuth`/`Full_elevation` [11, 9] give exact min/centre/max corners (azimuth depends on elevation row).

**Do DELTA_PLUS/DELTA_MINUS variables exist? Yes** — for energy (asymmetric: `delta_p_Energy`, `delta_m_Energy`) and for angles (symmetric: `delta_Azimuth`, `delta_Elevation`), wired via the CDF `DELTA_PLUS_VAR`/`DELTA_MINUS_VAR` attributes. These are exactly what you need for the velocity‑space volume element.

**Moment integration (d³v = v² dv dΩ).** With dΩ = cos(El) dEl dAz (El is a latitude‑like polar angle measured from the X–Y plane), and using the per‑bin widths:

$$\int v^2\,dv = \frac{v_{\max}^3 - v_{\min}^3}{3},\qquad \int_{\rm bin} d\Omega = \big(\Delta Az_{\rm rad}\big)\,\big[\sin(El+\delta_{El}) - \sin(El-\delta_{El})\big].$$

So the discrete moments (m = proton mass) are:

$$n=\sum_{ie,iel,iaz} f\, v^2\,\Delta v\,\cos(El)\,\Delta El\,\Delta Az,\qquad \mathbf V=\frac{1}{n}\sum f\,\mathbf v\,v^2\,\Delta v\,\cos(El)\,\Delta El\,\Delta Az,$$
$$\mathsf{P}_{jk}=m\sum f\,(v_j-V_j)(v_k-V_k)\,v^2\,\Delta v\,\cos(El)\,\Delta El\,\Delta Az,$$

with $v_{\min/\max}=13.85\sqrt{E\mp\delta E}$, angles converted to radians, and $\mathbf v=|\mathbf v|\hat V$ using the §3.2(b) signed unit vectors (then optionally rotate to RTN via `PAS_to_RTN`).

**Caveats for moments (Fedorov §3):** VDF can be noisy/unstable in 300–400 eV (geometric factor falling) and is **unreliable below 300 eV**; sporadic low‑level "ghost" values appear at extreme elevation/azimuth bins; ground densities/pressures are perturbed for bulk speed 260–380 km/s and invalid below 260 km/s. In the ground‑moment product, use `validity` ≥ 2 (avoid 1).

---

## Sources

- Owen, C. J., et al. 2020, "The Solar Orbiter Solar Wind Analyser (SWA) suite," *A&A*, **642**, A16 — DOI `10.1051/0004-6361/201937259`. Full text: https://www.aanda.org/articles/aa/full_html/2020/10/aa37259-19/aa37259-19.html ; open PDF: https://discovery.ucl.ac.uk/id/eprint/10109649/1/SWA_Instrument_Paper_for_production.pdf ; ADS: https://ui.adsabs.harvard.edu/abs/2020A%26A...642A..16O/abstract
- Fedorov, A. 2020, *SWA‑PAS L2 Data User Guide*, V02 (16 Nov 2020), ESA COSMOS — PDF: https://www.cosmos.esa.int/documents/3689933/11863901/PAS_L2_Data_User_Guide_20201116.pdf (variable tables, units, frame & velocity‑reconstruction formulae, E2V=13.85, caveats).
- SWA‑PAS L2 Data User Guide (online, readthedocs) — product/variable listing incl. `vdf [time,11,9,96]`, `Energy/Azimuth/Elevation`, `delta_p_Energy`/`delta_m_Energy`, `PAS_to_RTN [time,3,3]`: https://swa-pas-data-user-guide.readthedocs.io/en/latest/products.html and home https://swa-pas-data-user-guide.readthedocs.io/en/latest/ ; mirror: https://rungk-om.github.io/pas-user-guide/
- Louarn, P., et al. 2024, "Multiscale … solar‑wind proton moments from SWA‑PAS," *A&A*, **682**, A44 — DOI `10.1051/0004-6361/202347874` (PAS moments / cross‑calibration in flight). PDF: https://verscharen.com/papers/louarn_etal24.pdf
- Livi, S., et al. 2023, "First results from the Solar Orbiter Heavy Ion Sensor," *A&A*, **676**, A36 — DOI `10.1051/0004-6361/202346304` (HIS, **not** PAS; cite only to disambiguate the "Livi 2023" reference): https://www.aanda.org/articles/aa/full_html/2023/08/aa46304-23/aa46304-23.html
- IRAP SWA/PAS instrument page: https://www.irap.omp.eu/en/project/solar-orbiter-pas-2/ ; CNES SWA page: https://solar-orbiter.cnes.fr/en/SOLO/GP_swa.htm

**Key local files produced for this task:** extracted instrument paper text at `/private/tmp/claude-501/-Users-alonsocervantesflores-Programaci-n-Melvin--Velocity-distribution-function/8d7892f7-2d12-4c46-bbfd-28b7902d9351/scratchpad/owen2020.txt` and the User Guide text at `/Users/alonsocervantesflores/.claude/projects/-Users-alonsocervantesflores-Programaci-n-Melvin--Velocity-distribution-function/8d7892f7-2d12-4c46-bbfd-28b7902d9351/tool-results/esa_userguide.txt`.

# ===== vdf-reconstruction-and-moments =====

# Reconstructing and Analyzing Measured Velocity Distribution Functions (VDFs) in Kinetic Plasma Physics

Notation: $f(\mathbf v)$ is the phase-space density (units $\mathrm{s^3\,m^{-6}}$ or $\mathrm{s^3\,cm^{-6}}$), $m$ and $q$ the species mass and charge, $\mathbf B$ the magnetic field, $\hat{\mathbf b}=\mathbf B/|\mathbf B|$. Subscripts $i,j,k$ index energy, elevation, and azimuth bins respectively. All frames are right-handed.

---

## 1. From instrument coordinates $(E,\theta,\phi)$ to a 3-D velocity scatter and to RTN

### 1.1 What an electrostatic analyzer (ESA) actually measures

A top-hat ESA such as **SPAN-i** on Parker Solar Probe selects particles by **energy-per-charge** $E/q$ (hemisphere voltage sweep, $k$-factor $E=k\,V_{\rm hemi}$, $k\approx16.7\ \mathrm{eV/V}$), by **elevation** $\theta$ (electrostatic deflectors, $\pm60^\circ$), and by **azimuth** $\phi$ (discrete anodes spanning $247.5^\circ$); SPAN-i has a $8\times32\times8$ (azimuth $\times$ energy $\times$ elevation) native resolution ([Livi et al. 2022, ApJ 938, 138](https://iopscience.iop.org/article/10.3847/1538-4357/ac93f5)). The measured quantity is a **count array** $C_{ijk}$ accumulated over $\Delta t$.

A particle counted in bin $(i,j,k)$ **arrives from** the look direction $\hat{\mathbf n}_{jk}$, so its velocity points the opposite way. With elevation $\theta$ measured from the instrument azimuthal ($x$–$y$) plane and azimuth $\phi$ in that plane:

$$
\hat{\mathbf n}_{jk}=\big(\cos\theta_j\cos\phi_k,\ \cos\theta_j\sin\phi_k,\ \sin\theta_j\big),
\qquad
\mathbf v_{ijk}=-\,v_i\,\hat{\mathbf n}_{jk}.
$$

The speed follows from the selected kinetic energy $E_i=q(E/q)_i$:

$$
\boxed{\,v_i=\sqrt{\dfrac{2E_i}{m}}\,},\qquad E_i=\tfrac12 m v_i^2 .
$$

(If a polar angle $\vartheta$ from the $z$-axis is used instead of elevation, replace $\cos\theta\to\sin\vartheta$, $\sin\theta\to\cos\vartheta$.)

### 1.2 Counts $\to$ phase-space density

The chain from counts to $f$ uses the differential directional number flux $j_N$ (particles $\mathrm{area^{-1}\,s^{-1}\,sr^{-1}\,eV^{-1}}$). Because the flux through an area element is $d\Phi=v\,f\,d^3v$ and $d^3v=v^2\,dv\,d\Omega$ with $dv=dE/(mv)$,

$$
j_N=\frac{d\Phi}{dE\,d\Omega}=\frac{v^2}{m}f=\frac{2E}{m^2}f
\;\Longrightarrow\;
\boxed{\,f=\frac{m^2}{2E^2}\,j_E=\frac{m^2}{2E}\,j_N\,}
$$

where $j_E=E\,j_N$ is the differential energy flux that the instrument pipeline reports. The counts are related to flux through the geometric factor $G$ (effective area $\times$ solid angle) and efficiency $\varepsilon$:

$$
C_{ijk}=j_N\,\varepsilon_{ijk}\,G_{ijk}\,\Delta E_i\,\Delta t
\;\Longrightarrow\;
f_{ijk}=\frac{m^2}{2E_i}\,\frac{C_{ijk}}{\varepsilon_{ijk}\,G_{ijk}\,\Delta E_i\,\Delta t}.
$$

This conversion ($f\propto E^{-2}\times$ energy-flux) is the standard ESA inversion; for SPAN-i, L2 products are delivered as differential energy flux from which $f$ is built this way ([Livi et al. 2022](https://iopscience.iop.org/article/10.3847/1538-4357/ac93f5)).

### 1.3 The velocity-space Jacobian

The spherical $\to$ Cartesian map $(v,\theta,\phi)\to(v_x,v_y,v_z)$ carries the volume element

$$
\boxed{\,d^3v=v^2\,dv\,d\Omega\,},\qquad
d\Omega=\cos\theta\,d\theta\,d\phi\ \ (\text{elevation convention}),
$$

(or $d\Omega=\sin\vartheta\,d\vartheta\,d\phi$ for a polar angle). The $v^2$ factor is the radial Jacobian and $\cos\theta$ (or $\sin\vartheta$) is the angular Jacobian — both are essential and a frequent source of error if dropped. This is exactly the weighting used when ESA samples on energy/angle shells are integrated to moments ([Wilson et al., arXiv:2505.09869](https://arxiv.org/pdf/2505.09869); [Broderick/Klein "Slepian reconstruction" I, arXiv:2501.17294](https://arxiv.org/html/2501.17294)).

### 1.4 Rotation into the physical (RTN) frame

Each Cartesian velocity is rotated from instrument $\to$ spacecraft (fixed mounting matrix $\mathsf R_{\rm I\to S}$) $\to$ RTN (time-dependent attitude, from SPICE quaternions):

$$
\mathbf v^{\rm RTN}=\mathsf R_{\rm S\to RTN}(t)\,\mathsf R_{\rm I\to S}\,\mathbf v^{\rm I}.
$$

Each measured $(i,j,k)$ bin thus becomes one point $\big(\mathbf v^{\rm RTN}_{ijk},\,f_{ijk}\big)$ — the **3-D velocity-space scatter** $(v_R,v_T,v_N)$. Rotations preserve $d^3v$ (orthogonal, $\det=1$), so moments are frame-independent.

---

## 2. Field-aligned (FA) coordinates and the gyrotropy assumption

### 2.1 Building the triad

Given $\mathbf B$ (RTN), set $\hat{\mathbf b}=\mathbf B/|\mathbf B|$. The parallel velocity is the projection

$$
v_\parallel=\mathbf v\cdot\hat{\mathbf b}.
$$

The perpendicular plane needs two unit vectors, and **the perp pair is non-unique**: any rigid rotation by a gyrophase angle $\psi$ about $\hat{\mathbf b}$ is an equally valid choice. One fixes the freedom with a reference vector $\hat{\mathbf r}$ (commonly the RTN-$\hat R$/sunward direction, the bulk velocity $\mathbf U$, or $\mathbf B\times\mathbf U$):

$$
\hat{\mathbf e}_{\perp 2}=\frac{\hat{\mathbf b}\times\hat{\mathbf r}}{|\hat{\mathbf b}\times\hat{\mathbf r}|},\qquad
\hat{\mathbf e}_{\perp 1}=\hat{\mathbf e}_{\perp 2}\times\hat{\mathbf b},\qquad
(\hat{\mathbf e}_{\perp 1},\hat{\mathbf e}_{\perp 2},\hat{\mathbf b})\ \text{right-handed}.
$$

Then $v_{\perp 1}=\mathbf v\cdot\hat{\mathbf e}_{\perp 1}$, $v_{\perp 2}=\mathbf v\cdot\hat{\mathbf e}_{\perp 2}$, and the gyrophase is $\psi=\arctan2(v_{\perp 2},v_{\perp 1})$.

### 2.2 Gyrotropy: $3\mathrm D\to2\mathrm D$ reduction

In a magnetized, low-$\beta$ plasma the strong background field is a symmetry axis: the distribution is invariant under rotation about $\hat{\mathbf b}$ (independent of $\psi$). This **gyrotropy** assumption lets one define

$$
\boxed{\,v_\perp=\sqrt{v_{\perp 1}^2+v_{\perp 2}^2}\,},\qquad
f(v_\parallel,v_{\perp 1},v_{\perp 2})\;\longrightarrow\;f(v_\parallel,v_\perp),
$$

obtained by averaging over gyrophase,

$$
f(v_\parallel,v_\perp)=\frac{1}{2\pi}\int_0^{2\pi}f(v_\parallel,v_\perp,\psi)\,d\psi .
$$

This 2-D **gyrotropic distribution function (GDF)** is the object the SPAN-i/MMS Slepian-basis reconstructions target, precisely because $\mathbf B$ "serves as a gyrotropic axis of symmetry" in low-$\beta$ plasma ([Broderick/Klein, arXiv:2501.17294](https://arxiv.org/html/2501.17294); [companion paper, ApJ 10.3847/1538-4357/ae1d71](https://iopscience.iop.org/article/10.3847/1538-4357/ae1d71)). It is also the natural frame for resolving field-aligned **proton beams**, which appear as a second population displaced in $v_\parallel$ ([Verniero et al. 2020, ApJS 248, 5](https://iopscience.iop.org/article/10.3847/1538-4365/ab86af)).

### 2.3 When gyrotropy fails: agyrotropy

Near shocks, reconnection sites, and electron/ion diffusion regions, the gyration is interrupted on a gyroperiod and $f$ depends on $\psi$ (**agyrotropy / non-gyrotropy**). The gyrophase average then discards real structure, so one must keep the full $f(v_\parallel,v_{\perp1},v_{\perp2})$ and quantify the departure with a scalar built from the pressure tensor rotated to the FA frame. The most-used measure is **Swisdak's $Q$**:

$$
\sqrt{Q}=\sqrt{\frac{P_{12}^2+P_{13}^2+P_{23}^2}{P_\perp^2+2P_\perp P_\parallel}}\in[0,1],
$$

with the field-aligned axis "1" along $\hat{\mathbf b}$ ([Swisdak 2016, GRL, 10.1002/2015GL066980](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1002/2015GL066980)). Related diagnostics are Scudder & Daughton's $A\!\varnothing_e$ (from the eigenvalues of the perpendicular $2\times2$ sub-block, $A\!\varnothing=2|\lambda_1-\lambda_2|/(\lambda_1+\lambda_2)$) and Aunai et al.'s $D_{\rm ng}$ (RMS of off-diagonal $P_{ij}$ normalized to thermal energy); all measure deviation of $\mathsf P$ from the gyrotropic diagonal form $\mathrm{diag}(P_\parallel,P_\perp,P_\perp)$ ([review of measures, search summary](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1002/2015GL066980)).

---

## 3. Plasma moments from a discrete VDF

### 3.1 Continuous definitions

The velocity moments (per species $s$) are ([Verscharen, Klein & Maruca 2019, *Living Rev. Solar Phys.* 16:5, Eqs. 27–29](https://pmc.ncbi.nlm.nih.gov/articles/PMC6934245/)):

$$
n=\int f\,d^3v
\tag{density}
$$
$$
\mathbf U=\frac1n\int \mathbf v\,f\,d^3v
\tag{bulk velocity}
$$
$$
\mathsf P=m\int(\mathbf v-\mathbf U)(\mathbf v-\mathbf U)\,f\,d^3v
\tag{pressure tensor}
$$

with temperature tensor $\mathsf T=\mathsf P/(n k_B)$ and heat flux $\mathbf q=\tfrac{m}{2}\int(\mathbf v-\mathbf U)|\mathbf v-\mathbf U|^2 f\,d^3v$. Equivalently $\mathsf P=m\!\int\!\mathbf v\mathbf v f\,d^3v-mn\,\mathbf U\mathbf U$ (raw second moment minus the dyad), which is numerically convenient.

### 3.2 Parallel/perpendicular temperatures

$$
T_\parallel=\hat{\mathbf b}\cdot\mathsf T\cdot\hat{\mathbf b}=\frac{p_\parallel}{n k_B},\qquad
T_\perp=\frac{\operatorname{Tr}\mathsf T-T_\parallel}{2}=\frac{p_\perp}{n k_B},
$$

(Verscharen et al. 2019, Eqs. 38–39), with scalar temperature $T=\tfrac13\operatorname{Tr}\mathsf T=\tfrac13(T_\parallel+2T_\perp)$, anisotropy $R=T_\perp/T_\parallel$, and plasma beta $\beta_s=8\pi n_s k_B T_s/B^2$ (their Eq. 15).

### 3.3 Discrete-sum form (the practical recipe)

Replace $d^3v\to v_i^2\,\Delta v_i\,\Delta\Omega_{jk}$ in each measured bin. Then

$$
\boxed{\,n=\sum_{i,j,k} f_{ijk}\;v_i^2\,\Delta v_i\,\Delta\Omega_{jk}\,}
$$
$$
n\,\mathbf U=\sum_{i,j,k}\mathbf v_{ijk}\,f_{ijk}\;v_i^2\,\Delta v_i\,\Delta\Omega_{jk},
\qquad \mathbf v_{ijk}=-v_i\hat{\mathbf n}_{jk}
$$
$$
\mathsf P=m\sum_{i,j,k}(\mathbf v_{ijk}-\mathbf U)(\mathbf v_{ijk}-\mathbf U)\,f_{ijk}\;v_i^2\,\Delta v_i\,\Delta\Omega_{jk}.
$$

**Velocity-bin width from energy-bin width.** From $E=\tfrac12 mv^2\Rightarrow dE=mv\,dv$:

$$
\Delta v_i=\frac{\Delta E_i}{m\,v_i}.
$$

For the usual logarithmically (geometrically) spaced energy channels with constant $\Delta E/E$, it is cleaner to use $\dfrac{dv}{v}=\tfrac12\dfrac{dE}{E}$, giving

$$
\Delta v_i=\tfrac12\,v_i\,\frac{\Delta E_i}{E_i},\qquad\text{so}\quad
v_i^2\,\Delta v_i=\tfrac12\,v_i^3\,\frac{\Delta E_i}{E_i}.
$$

**Solid-angle element from angular-bin widths.**

$$
\Delta\Omega_{jk}=\cos\theta_j\,\Delta\theta_j\,\Delta\phi_k
\quad(\text{elevation}),\qquad
\Delta\Omega_{jk}=\sin\vartheta_j\,\Delta\vartheta_j\,\Delta\phi_k\quad(\text{polar}).
$$

Equivalently, working directly from differential energy flux $j_E$ via $f=\tfrac{m^2}{2E^2}j_E$ turns the density sum into the familiar instrument form $n=\sum (m^2/2E_i^2)\,j_{E,ijk}\,v_i^2\Delta v_i\Delta\Omega_{jk}$; this is the discretization ESA pipelines use, and its accuracy is limited by finite energy/angle resolution and field-of-view truncation ([Wilson et al., arXiv:2505.09869](https://arxiv.org/pdf/2505.09869)). Caveats that bias these sums: truncated FOV (partial-sky coverage), finite low-/high-energy cutoffs, spacecraft potential shifting $E$, and one-count noise — each propagates differently into $n$, $\mathbf U$, and $\mathsf P$.

---

## 4. Thermal speeds and peculiar velocities

### 4.1 Thermal-speed definitions

Using the **"most probable speed"** convention standard in solar-wind kinetics ([Verscharen et al. 2019, Eqs. 4–5](https://pmc.ncbi.nlm.nih.gov/articles/PMC6934245/); [Wikipedia, *Thermal velocity*](https://en.wikipedia.org/wiki/Thermal_velocity)):

$$
\boxed{\,w=\sqrt{\frac{2k_B T}{m}}\,},\qquad
w_\parallel=\sqrt{\frac{2k_B T_\parallel}{m}},\qquad
w_\perp=\sqrt{\frac{2k_B T_\perp}{m}}.
$$

(Beware competing conventions: $\sqrt{k_BT/m}$, $\sqrt{3k_BT/m}$, $\sqrt{8k_BT/\pi m}$ all appear in the literature; the $\sqrt{2k_BT/m}$ form is the one that makes the Maxwellian exponent $-c^2/w^2$.)

### 4.2 Normalized peculiar velocity

$$
\xi_\parallel=\frac{v_\parallel-U_\parallel}{w_\parallel},\qquad
\xi_\perp=\frac{v_\perp-U_\perp}{w_\perp}.
$$

In the bulk frame the perpendicular drift is removed ($U_\perp\approx0$, the $\mathbf E\times\mathbf B$ frame), so $\xi_\perp\simeq v_\perp/w_\perp$. A **bi-Maxwellian** is then $f=\dfrac{n}{\pi^{3/2}w_\parallel w_\perp^2}\exp(-\xi_\parallel^2-\xi_\perp^2)$.

### 4.3 The cleaner "bulk-rest-frame peculiar velocity"

Define the peculiar (random) velocity by subtracting the bulk flow **before** projecting, which avoids having to define $U_\perp$ or a gyrophase:

$$
\mathbf c=\mathbf v-\mathbf U,\qquad
c_\parallel=\mathbf c\cdot\hat{\mathbf b},\qquad
\boxed{\,c_\perp=\big|\mathbf c-c_\parallel\hat{\mathbf b}\big|=\sqrt{|\mathbf c|^2-c_\parallel^2}\,}.
$$

Here $c_\perp\ge0$ by construction and the magnitude is rotation-invariant about $\hat{\mathbf b}$, so it is the robust quantity to histogram for a gyrotropic $f(c_\parallel,c_\perp)$ and to feed the pressure-tensor sum, $\mathsf P=m\sum \mathbf c\,\mathbf c\, f\,v^2\Delta v\,\Delta\Omega$, with $p_\parallel=m\sum c_\parallel^2 f\,(\cdots)$ and $p_\perp=\tfrac{m}{2}\sum c_\perp^2 f\,(\cdots)$.

---

## Key references

- [Livi, R., et al. 2022, *The Solar Probe ANalyzer—Ions on PSP*, ApJ 938, 138](https://iopscience.iop.org/article/10.3847/1538-4357/ac93f5) — SPAN-i instrument frame ($E/q$, deflection $\theta$, anode $\phi$), $k$-factor, counts→flux.
- [Verscharen, D., Klein, K. G. & Maruca, B. A. 2019, *The multi-scale nature of the solar wind*, Living Rev. Solar Phys. 16:5 (PMC mirror)](https://pmc.ncbi.nlm.nih.gov/articles/PMC6934245/) — moment definitions (Eqs. 27–29), $T_\parallel,T_\perp$ (38–39), thermal speeds (4–5), $\beta$ (15).
- [Broderick/Klein et al. 2025, *Recovering Ion Distribution Functions I: Slepian Reconstruction (MMS & Solar Orbiter)*, arXiv:2501.17294](https://arxiv.org/html/2501.17294) and [II: Gyrotropic Slepian Reconstruction, ApJ 10.3847/1538-4357/ae1d71](https://iopscience.iop.org/article/10.3847/1538-4357/ae1d71) — instrument→velocity reconstruction and the gyrotropic 2-D reduction.
- [Wilson et al. 2025, *How limited resolution of plasma analyzers affects moment accuracy*, arXiv:2505.09869](https://arxiv.org/pdf/2505.09869) — discrete moment sums over energy/angle bins, $v^2\,dv\,d\Omega$ discretization, resolution biases.
- [Verniero, J. L., et al. 2020, *PSP Observations of Proton Beams Simultaneous with Ion-scale Waves*, ApJS 248, 5](https://iopscience.iop.org/article/10.3847/1538-4365/ab86af) — field-aligned VDFs and beams.
- [Swisdak, M. 2016, *Quantifying gyrotropy in magnetic reconnection*, GRL, 10.1002/2015GL066980](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1002/2015GL066980) — agyrotropy measure $Q$ (and references to Scudder & Daughton 2008 $A\!\varnothing$, Aunai et al. 2013 $D_{\rm ng}$).
- [Wikipedia, *Thermal velocity*](https://en.wikipedia.org/wiki/Thermal_velocity) — thermal-speed conventions $w=\sqrt{2k_BT/m}$ and alternatives.

Sources:
- https://iopscience.iop.org/article/10.3847/1538-4357/ac93f5
- https://pmc.ncbi.nlm.nih.gov/articles/PMC6934245/
- https://arxiv.org/html/2501.17294
- https://iopscience.iop.org/article/10.3847/1538-4357/ae1d71
- https://arxiv.org/pdf/2505.09869
- https://iopscience.iop.org/article/10.3847/1538-4365/ab86af
- https://agupubs.onlinelibrary.wiley.com/doi/full/10.1002/2015GL066980
- https://en.wikipedia.org/wiki/Thermal_velocity

# ===== parametric-models =====

# Solar-Wind Ion Velocity Distribution Functions: Models, Fitting, and Model-Free Decomposition

The velocity distribution function (VDF) $f(\mathbf{v})$ gives the phase-space density of particles; all fluid moments (density, bulk velocity, pressure/temperature tensor, heat flux) are its velocity-space integrals. Because the solar wind is weakly collisional, its VDFs retain non-equilibrium structure that fluid theory averages away, which is precisely why kinetic modeling of the measured $f(\mathbf{v})$ matters (Marsch 2006; Verscharen, Klein & Maruca 2019).

---

## 1. Maxwellian / bi-Maxwellian model and the three proton-scale populations

### 1.1 Maxwellian = local thermodynamic equilibrium

A single drifting **Maxwellian** is the maximum-entropy solution of the Boltzmann equation and is the unique VDF of a collisionally relaxed gas in **thermal equilibrium**:

$$
f_M(\mathbf{v}) = \frac{n}{\pi^{3/2}\,w^3}\,\exp\!\left[-\frac{(\mathbf{v}-\mathbf{u})^2}{w^2}\right],\qquad w=\sqrt{\tfrac{2k_BT}{m}},
$$

where $w$ is the thermal speed and $T$ the (isotropic) temperature. It is isotropic in the bulk-flow frame; its only "structure" is density, drift, and one temperature.

### 1.2 Bi-Maxwellian = the simplest non-equilibrium (anisotropic) model

In a magnetized, weakly collisional plasma the magnetic field $\mathbf{B}$ breaks the isotropy: gyromotion thermalizes the two perpendicular directions but parallel and perpendicular temperatures need not be equal. The **bi-Maxwellian** (Verscharen, Klein & Maruca 2019, their Eq. 61) is

$$
f_{bM}(v_\parallel,v_\perp)=\frac{n}{\pi^{3/2}\,w_\parallel w_\perp^{2}}\,
\exp\!\left[-\frac{(v_\parallel-u_\parallel)^2}{w_\parallel^{2}}\right]
\exp\!\left[-\frac{(v_\perp-u_\perp)^2}{w_\perp^{2}}\right],
$$

with $w_{\parallel,\perp}=\sqrt{2k_BT_{\parallel,\perp}/m}$ and anisotropy $A=T_\perp/T_\parallel$. For a **gyrotropic** plasma $u_\perp=0$ and the perpendicular factor reduces to $\exp(-v_\perp^2/w_\perp^2)$; a nonzero $u_\perp$ encodes a perpendicular drift. The bi-Maxwellian "introduces temperature anisotropies with respect to the background magnetic field yet follows the Maxwellian behavior on any one-dimensional cut" (Verscharen et al. 2019). Fast-wind protons persistently show $T_{\perp p}>T_{\parallel p}$, which adiabatic expansion alone cannot produce ($T_\perp \ll T_\parallel$ would be expected), demanding active perpendicular heating — a kinetic signature fluid models miss.

### 1.3 The three proton-scale populations

Solar-wind proton VDFs measured since Helios are well organized with respect to $\mathbf{B}$ and are typically modeled as a **sum of bi-Maxwellians** (Marsch et al. 1982; Marsch 2006):

- **Proton core** — the global peak, containing ~90–95% of protons; close to a (an)isotropic Maxwellian and defining the bulk velocity. It is the densest, "coolest" component.

- **Proton beam** — "a second proton component streaming faster than the proton core component along the direction of the magnetic field with a relative speed $\gtrsim v_{Ap}$" (Verscharen et al. 2019), i.e. a field-aligned shoulder displaced sunward-outward, with drift $\Delta u \sim (1\text{–}2)\,v_{Ap}$ — corresponding to roughly $1.5\text{–}2\times$ the core thermal width along $\mathbf B$. The beam is "almost always directed away from the Sun and along the magnetic-field axis." A separate beam bi-Maxwellian (own $n_b,\,u_b,\,w_{\parallel b},\,w_{\perp b}$) is added to the core.

- **Alpha particles (He²⁺)** — a distinct ion species ($m_\alpha=4m_p$, charge $2e$, so $m/q=2$). They form $\lesssim 20\%$ of the mass density and "drift with respect to the proton core along the magnetic-field direction and away from the Sun with a typical drift speed $\lesssim v_{Ap}$" (Verscharen et al. 2019). Because electrostatic analyzers sort ions by **energy-per-charge** $E/q=\tfrac12 m v^2/q$, an alpha at the *same bulk speed* as the proton core appears at $E_\alpha/q_\alpha = 2\,(E_p/q_p)$ — i.e. at **double the proton $E/q$**. Equivalently, if that $E/q$ peak is (mis)converted to a speed assuming proton mass, alphas show up at $\sqrt{2}\times$ the proton speed: $v_{\rm app}=\sqrt{2(E/q)e/m_p}=\sqrt{2}\,v_p$. This $m/q=2$ offset is exactly what lets alphas be spectrally separated from protons.

---

## 2. Kappa and Gaussian/regularized-kappa distributions

### 2.1 The standard kappa VDF

Collisionless space plasmas almost universally show **suprathermal power-law tails** that a Maxwellian cannot fit. Vasyliūnas (1968) and Olbert (1968) introduced the **kappa distribution**, Maxwellian-like at low energy but a power law at high energy. In the Verscharen et al. (2019) convention (their Eq. 62), the isotropic kappa VDF is

$$
f_\kappa(\mathbf{v})=\frac{n}{\big[\pi(2\kappa-3)\big]^{3/2}\,w^{3}}\;
\frac{\Gamma(\kappa+1)}{\Gamma\!\left(\kappa-\tfrac12\right)}\;
\left[\,1+\frac{2(\mathbf{v}-\mathbf{u})^2}{(2\kappa-3)\,w^{2}}\,\right]^{-(\kappa+1)},
$$

where the $(2\kappa-3)$ factor is chosen so that $w$ remains the thermal speed and $T=mw^2/2k_B$ **independently of $\kappa$** (valid for $\kappa>3/2$). A common equivalent "Olbertian" form uses a thermal parameter $\theta$:

$$
f_\kappa(v)=\frac{n}{\pi^{3/2}\theta^{3}}\,\frac{\Gamma(\kappa+1)}{\kappa^{3/2}\,\Gamma(\kappa-\tfrac12)}
\left[1+\frac{v^2}{\kappa\theta^2}\right]^{-(\kappa+1)} .
$$

### 2.2 The kappa index $\kappa$: control of the "fat tails"

The spectral index $\kappa$ governs how non-Maxwellian the VDF is: tails are "more pronounced for smaller $\kappa$," and $f_\kappa\to f_M$ as $\kappa\to\infty$ (Verscharen et al. 2019). A low $\kappa$ means strong suprathermalization; a high $\kappa$ means a nearly thermalized population. The high-energy slope scales as $v^{-2(\kappa+1)}$. Physically, smaller $\kappa$ correlates with **weaker collisions** (lower density, larger heliocentric distance), so the tail departs further from equilibrium (Štverák et al. 2022).

### 2.3 Electron halo and the Gaussian/dual-kappa model

For solar-wind **electrons**, the standard decomposition is a **Maxwellian core + kappa halo** (a "Gaussian-kappa" or dual Maxwellian–kappa model; Pierrard & Lazar 2010; Lazar et al. 2017):

$$
f_e(\mathbf v)=f_{M}^{\rm core}(\mathbf v)+f_{\kappa}^{\rm halo}(\mathbf v),
$$

with the cold near-isotropic Maxwellian core fitting the $\lesssim 10$ s of eV bulk and the kappa halo capturing the $\sim 10^2\text{–}10^3$ eV suprathermal population (a sunward-aligned **strahl** is a third, field-aligned component). Fitted halo indices typically run $\kappa\sim 2\text{–}6$, decreasing with distance (e.g. $\kappa\approx 7.5$ at 0.4 AU, $5.4$ at 1 AU, $3.2$ at 3 AU; Štverák et al. 2022).

### 2.4 Regularized (Gaussian-cutoff) kappa

A known defect of the standard kappa is that velocity moments diverge for low $\kappa$ (and it implies unphysical superluminal particles). The **regularized kappa distribution (RKD)** of Scherer, Fichtner & Lazar (2017) multiplies the kappa core by a **Gaussian cutoff** $\exp(-\alpha^2 v^2/\theta^2)$ with $0<\alpha\ll1$:

$$
f_{\rm RKD}(v)\;\propto\;\left[1+\frac{v^2}{\kappa\theta^2}\right]^{-(\kappa+1)}\exp\!\left(-\alpha^{2}\frac{v^{2}}{\theta^{2}}\right),
$$

which keeps all moments finite, extends the allowed range to $\kappa>0$, and is now used for anisotropic in-situ electron halos and exospheric solar-wind models (Scherer et al. 2017; "regularized kappa-halos," Phys. Plasmas 2025).

---

## 3. Fitting: $\chi^2$ / Levenberg–Marquardt, GMM–EM, and "human bias"

### 3.1 Deterministic least-squares: Levenberg–Marquardt

The classical approach fits a chosen parametric model $f(\mathbf v;\boldsymbol\theta)$ (one or more bi-Maxwellians, or a kappa) by minimizing the (often weighted) chi-square between model and measured phase-space density over the velocity-space bins $i$:

$$
\chi^2(\boldsymbol\theta)=\sum_i \frac{\big[f_{\rm obs}(\mathbf v_i)-f(\mathbf v_i;\boldsymbol\theta)\big]^2}{\sigma_i^2}.
$$

The **Levenberg–Marquardt (LM)** algorithm interpolates between Gauss–Newton and gradient descent via a damping parameter $\lambda$,

$$
(\mathbf J^\top\mathbf J+\lambda\,\mathrm{diag}(\mathbf J^\top\mathbf J))\,\delta\boldsymbol\theta=\mathbf J^\top\mathbf r,
$$

($\mathbf J$ = Jacobian of residuals $\mathbf r$), giving robust convergence for nonlinear VDF fits. In heliophysics this is most commonly the **MPFIT** implementation (Markwardt 2009), used e.g. to fit core/beam/alpha bi-Maxwellians.

### 3.2 Unsupervised machine learning: GMM + Expectation–Maximization

De Marco et al. (2023, A&A) separate proton core, proton beam, and alphas in Solar Orbiter/PAS 3D VDFs with a **Gaussian Mixture Model (GMM)**: each measured velocity point is modeled as drawn from a mixture of $K$ Gaussians,

$$
P(\mathbf x_i\mid\boldsymbol\theta)=\sum_{k=1}^{K}w_k\,\mathcal N(\mathbf x_i\mid\boldsymbol\mu_k,\boldsymbol\Sigma_k),\qquad \sum_k w_k=1,
$$

whose parameters $\boldsymbol\theta=\{w_k,\boldsymbol\mu_k,\boldsymbol\Sigma_k\}$ are found by **Expectation–Maximization (EM)**: the E-step computes each point's *responsibility* (soft membership) $\gamma_{i,k}$, the M-step updates $\{w_k,\boldsymbol\mu_k,\boldsymbol\Sigma_k\}$ to maximize the log-likelihood, iterating to convergence. The authors stress that, unlike LM fitting, "there is no need to specify the magnetic field orientation or the expected location for beam and alphas to find the families," and each point gets a probabilistic membership rather than a forced category. The same GMM/EM machinery has been applied to magnetospheric (MMS) reconnection VDFs (Dupuis et al. 2020). EM is still anisotropic-Gaussian-based, so it inherits a bi-Maxwellian assumption, but it removes the manual seeding and region-masking bias of LM.

### 3.3 The "human bias" of model choice

Every parametric fit embeds physical *assumptions* that bias the result before any data are touched:

- **Thermal equilibrium / number of components** — choosing 1 vs 2 vs 3 Maxwellians, or Maxwellian-vs-kappa, predetermines whether a "beam" or "tail" exists. Overlapping core/beam/alpha curves make this "a challenging task due to the overlapping of the curves" (De Marco et al. 2023).
- **Collisionless tails** — using a Maxwellian forces $\kappa\to\infty$ and discards real suprathermal populations; using a kappa imposes a specific power-law form.
- **Anisotropy & symmetry** — bi-Maxwellians and gyrotropic kappas assume gyrotropy (azimuthal symmetry about $\mathbf B$) and reflection symmetry, discarding skewness, agyrotropy, and heat-flux distortions that are physically present. The "normal inverse-Gaussian" work (A&A 2024) shows real proton VDFs carry skewness/kurtosis no symmetric model captures.

Thus the recovered $n,\,\mathbf u,\,T_{\parallel,\perp},\,\kappa$ depend on the modeler's prior — different but "reasonable" models give different moments and instability thresholds from the *same* data.

---

## 4. Why a model-free decomposition is desirable

To see the VDF "as it is," one wants a **basis expansion** that is complete and physically natural but imposes no population count, no equilibrium assumption, and no symmetry — letting the data populate the coefficients.

A natural choice is the **Hermite expansion in velocity space** (with Fourier modes in physical space), because the Hermite functions are the eigenfunctions of the linearized collision/streaming operators and the lowest modes *are* the fluid moments:

$$
f(\mathbf v)=\sum_{n}\,c_{n}\,H_{n}\!\left(\tfrac{v-u}{w}\right)\,\frac{e^{-(v-u)^2/w^2}}{\sqrt{\pi}\,w},
$$

where $\{H_n\}$ are orthogonal under weight $e^{-v^2}$. The $n=0,1,2$ coefficients reproduce density, momentum, and pressure (a Maxwellian has only $c_0\neq0$); higher $n$ measure non-Maxwellianity, and the $\mathbf{B}$-parallel hierarchy turns phase mixing into simple coupling between adjacent Hermite moments (Schekochihin et al.; Servidio et al. 2017, MMS velocity-space cascade; Vencels et al. 2018, Fourier–Hermite spectral methods). Crucially, the **Hermite power spectrum** $|c_n|^2$ exposes a *velocity-space cascade* — structure at small velocity scales — that no finite sum of Maxwellians can represent.

Closely related model-free reconstructions include **Slepian / spherical-harmonic reconstruction** of MMS and Solar Orbiter ion VDFs (Verniero et al. / Broiles-style methods; "Recovering Ion Distribution Functions: Slepian Reconstruction," ApJ 2025), which optimally interpolate the measured counts on the spacecraft's angular–energy grid without assuming any analytic form.

The payoff: a decomposition (i) is unbiased by population count or equilibrium assumptions, (ii) keeps *all* the information (skewness, agyrotropy, fine velocity-space structure), (iii) yields a quantitative "distance from Maxwellian," and (iv) connects directly to kinetic theory (phase mixing, Landau damping, turbulent velocity-space cascade) — revealing the VDF as the measurement actually constrains it, rather than as a preselected model permits.

---

## Sources

- [Verscharen, Klein & Maruca (2019), *The multi-scale nature of the solar wind*, Living Rev. Sol. Phys. — bi-Maxwellian (Eq. 61), kappa (Eq. 62), core/beam/alpha drifts](https://pmc.ncbi.nlm.nih.gov/articles/PMC6934245/) · [arXiv](https://arxiv.org/abs/1902.03448)
- [Marsch (2006), *Kinetic Physics of the Solar Corona and Solar Wind*, Living Rev. Sol. Phys.](https://link.springer.com/article/10.12942/lrsp-2006-1)
- [De Marco et al. (2023), *Innovative technique for separating proton core, proton beam, and alpha particles…* (GMM + EM), A&A](https://www.aanda.org/articles/aa/full_html/2023/01/aa43719-22/aa43719-22.html)
- [Dupuis et al. (2020), *Characterizing magnetic reconnection regions using Gaussian mixture models on particle VDFs*, arXiv](https://arxiv.org/pdf/1910.10012)
- [Štverák et al. (2022), *Implications of Kappa Suprathermal Halo of the Solar Wind Electrons*, Front. Astron. Space Sci.](https://www.frontiersin.org/journals/astronomy-and-space-sciences/articles/10.3389/fspas.2022.892236/full)
- [Lazar et al. (2017), *Dual Maxwellian–Kappa modeling of the solar wind electrons*, A&A](https://www.aanda.org/articles/aa/full_html/2017/06/aa30194-16/aa30194-16.html)
- [Pierrard & Lazar (2010), *Kappa Distributions: Theory and Applications in Space Plasmas*, Solar Phys.](https://www.researchgate.net/publication/45906919_Kappa_Distributions_Theory_and_Applications_in_Space_Plasmas)
- [Scherer, Fichtner & Lazar — regularized kappa distributions (RKD) overview / applications](https://www.mdpi.com/2571-6182/6/3/36) · [Regularized kappa-halos with ALPS, arXiv](https://arxiv.org/html/2504.15955)
- [Servidio et al. (2017), *MMS observation of plasma velocity-space cascade: Hermite representation and theory*, arXiv](https://arxiv.org/pdf/1707.08180)
- [Vencels et al. (2018), *Spectral Approach to Plasma Kinetic Simulations Based on Hermite Decomposition*, Front. Astron. Space Sci.](https://www.frontiersin.org/journals/astronomy-and-space-sciences/articles/10.3389/fspas.2018.00027/full)
- [*Recovering Ion Distribution Functions: Slepian Reconstruction of VDFs from MMS and Solar Orbiter* (2025), ApJ](https://iopscience.iop.org/article/10.3847/1538-4357/adb6a0)
- [*Skewness and kurtosis of solar wind proton distribution functions* (normal inverse-Gaussian) (2024), A&A](https://www.aanda.org/articles/aa/full_html/2024/02/aa47874-23/aa47874-23.html)
- [Tu et al. (2004), *Dependence of the proton beam drift velocity on proton core plasma beta*, JGR](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2004JA010391)

# ===== hermite-laguerre-literature =====

I now have the primary sources fully extracted (Pezzi 2018 and Larosa 2025 PDFs read directly; Coburn 2024 and Servidio 2017 from publisher/abstract). I have enough to write the explainer.

One important correction surfaced in the research that I'll flag: the **Larosa et al. 2025 Hermite paper is on Parker Solar Probe (PSP), not Solar Orbiter** (arXiv:2512.01492). The method details the user described (HH, Cartesian-like gyrotropy, mirrored perpendicular axis, 50×50 grid) match that paper exactly.

---

# Spectral (Hermite / Hermite–Laguerre) Representations of Plasma Velocity Distributions — A Citation-Rich Explainer

## 0. The four papers at a glance

| Paper | Data / system | Basis | Gyrotropy | Grid | Key spectral result |
|---|---|---|---|---|---|
| **Servidio et al. 2017**, *PRL* **119**, 205101 ([arXiv:1707.08180](https://arxiv.org/abs/1707.08180), [APS](https://link.aps.org/doi/10.1103/PhysRevLett.119.205101)) | MMS ions, Earth's magnetosheath | 3D Hermite | Cartesian (3D $v_x,v_y,v_z$) | — | First observed power-law $P(m)$; $m^{-3/2}$ (weakly magnetized) vs $m^{-2}$ (magnetized) |
| **Pezzi et al. 2018**, *Phys. Plasmas* **25**, 060704 ([arXiv:1803.01633](https://arxiv.org/abs/1803.01633), [AIP](https://pubs.aip.org/aip/pop/article/25/6/060704/320076)) | Hybrid Vlasov–Maxwell (HVM) sim, 2.5D–3V | 3D Hermite | Cartesian + $\mathbf B_0\parallel \hat z$ | $N_m=100$/dir, ensemble of $64^3$ VDFs | Anisotropic v-space cascade: parallel $P(m_z)\sim m_z^{-2.01}$, perpendicular $P(m_\perp)\sim m_\perp^{-3.5}$ |
| **Coburn et al. 2024**, *ApJ* **964**, 100 ([IOP](https://iopscience.iop.org/article/10.3847/1538-4357/ad1329), [ADS](https://ui.adsabs.harvard.edu/abs/2024ApJ...964..100C)) | Solar Orbiter SWA electrons | **Hermite–Laguerre** | **Cylindrical** ($\mu=v_\perp^2$) | $60\times60$ | Spectral denoising / low-pass truncation at $m,l\approx14$ noise floor |
| **Larosa et al. 2025** ([arXiv:2512.01492](https://arxiv.org/abs/2512.01492)) | **Parker Solar Probe** SPAN-i ions (not Solar Orbiter) | **Hermite–Hermite** | **Cartesian-like**, $f(-v_\perp)=f(v_\perp)$ | $50\times50$ | Isotropic shell spectrum $P(m)\sim m^{-2}$ over $4\le m<12$; dual real/v-space cascade |

A foundational methods companion is **Roytershteyn & Delzanno 2018**, *Front. Astron. Space Sci.* **5**, 27 ([Frontiers](https://www.frontiersin.org/journals/astronomy-and-space-sciences/articles/10.3389/fspas.2018.00027/full)), and the gyrokinetic Hermite–Laguerre formulation **Mandell et al. 2018** ([arXiv:1708.04029](https://arxiv.org/pdf/1708.04029)).

> **Attribution caveat up front.** Two specifics in your prompt I could *not* find verbatim in any of these papers: the slope pair "$m^{-1.3}$, $n^{-3.1}$" and the local-slope formula $\alpha_m=(c_{m+2,0}^2-c_{m,0}^2)/(2c_{00}^2)$. They are consistent in *spirit* with the anisotropic-cascade results (shallow parallel, steep perpendicular) and local-exponent estimators in the literature, but appear to be your own working conventions or from a Solar-Orbiter-specific analysis I could not locate. I treat them as such below and explain the underlying concepts rigorously.

---

## 1. The Fourier analogy: a VDF as a spectrum of Hermite/Laguerre modes

In a collisionless plasma the VDF is generically *non-Maxwellian* — beams, rings, temperature anisotropy, "hammerhead" suprathermal tails. Just as a Fourier transform writes a signal as a sum of plane waves, the **Hermite transform writes the VDF as a Maxwellian times a series of Hermite polynomials**. Larosa et al. state the analogy explicitly: *"a single Fourier mode describes a plane wave, and likewise a single Hermite mode describes a Maxwellian; the more Hermite modes are non-zero the more the VDF is far from a Maxwellian"* ([Larosa et al. 2025](https://arxiv.org/abs/2512.01492)). The expansion is exact and invertible (a spectral basis), and **Parseval–Plancherel holds**, so the velocity-space "energy" (enstrophy) is preserved between physical and spectral space.

### 1.1 The orthonormal Hermite functions

In field-aligned coordinates $(\xi_\parallel,\xi_\perp)$ the gyrotropic VDF is expanded as

$$f(\xi_\parallel,\xi_\perp)=\sum_{m,n} c_{mn}\,\psi_m(\xi_\parallel)\,\psi_n(\xi_\perp),$$

with the **orthonormal Hermite functions** (Larosa Eq. 2; Pezzi Eq. 4):

$$\boxed{\;\psi_m(\xi)=\frac{H_m(\xi)}{\sqrt{2^m\,m!\,\sqrt\pi}}\;e^{-\xi^2/2}\;}$$

where $H_m$ are the *physicist's* Hermite polynomials, $H_m(\xi)=(-1)^m e^{\xi^2}\frac{d^m}{d\xi^m}e^{-\xi^2}$ (Larosa Eq. 3). These obey the orthonormality relation that makes the transform "Fourier-like":

$$\int_{-\infty}^{\infty}\psi_m(\xi)\,\psi_n(\xi)\,d\xi=\delta_{mn}.$$

Key normalization choices (Pezzi §; Larosa §2):
- The argument is the **normalized peculiar velocity** $\xi=(v-u)/v_{\rm th}$, i.e. the basis is *shifted to the local bulk frame* and *scaled by the thermal speed* in each direction.
- This means $m=0$ is the *local* drifting Maxwellian; the projection "quantifies high-order corrections to the particle DF, since the basis is shifted in the local fluid velocity frame, normalized to the ambient density and temperature" (Pezzi). Consequently large-scale fluctuations of density, bulk flow and temperature (1st/2nd moments) **do not contaminate the high-$m$ spectrum** (Larosa §2).
- $e^{-\xi^2/2}$ (not $e^{-\xi^2}$) makes the *functions* (not just polynomials) orthonormal — these are the "symmetrically weighted" Hermite functions / quantum-harmonic-oscillator eigenstates.

### 1.2 Two ways to handle the perpendicular plane — Cartesian (HH) vs cylindrical (HL)

Gyrotropy means $f$ is independent of gyrophase: $f=f(v_\parallel,v_\perp)$ with $v_\perp=\sqrt{v_{\perp1}^2+v_{\perp2}^2}$. There are two natural spectral encodings of that perpendicular plane, and this is the central methodological difference between Larosa and Coburn:

**(a) Hermite–Hermite (HH), "Cartesian-like" gyrotropy — Larosa et al. 2025.**
Treat $v_\perp$ as if it were a Cartesian axis and use a Hermite function $\psi_n(\xi_\perp)$ in the perpendicular direction too. Gyrotropy is imposed by **mirroring**: *"we impose $f(-v_\perp)=f(v_\perp)$ and extend the grid to negative values in the perpendicular direction consistently with the gyrotropy assumption … such a procedure implies null odd Hermite coefficients in the perpendicular direction"* (Larosa §2). So $c_{m,n}=0$ for odd $n$ — only even perpendicular modes carry power. This is exactly why their 1D spectrum shows an **oscillatory (even/odd) pattern** (Larosa Fig. 3): *"the oscillatory behavior … is due to the lack of power in the odd modes of the perpendicular spectrum caused by the gyrotropy assumption."*

**(b) Hermite–Laguerre (HL), cylindrical gyrotropy — Coburn et al. 2024.**
The *natural* basis for an azimuthally symmetric function in cylindrical velocity coordinates uses **associated Laguerre polynomials in $\mu=v_\perp^2$**. Coburn writes the electron VDF as

$$F_e(v_\parallel,v_\perp)=\sum_{m}\sum_{l} c_{ml}\,\psi_m(\hat v_\parallel)\,\Gamma_l^{1}(\hat v_\perp),\qquad \Gamma_l^{n}(\hat v_\perp)=e^{-\hat v_\perp^2/2}\,L_l^{n}(\hat v_\perp^2),$$

with $\hat v_\parallel=v_\parallel/v_{{\rm th},e}^\parallel$, $\hat v_\perp=v_\perp/v_{{\rm th},e}^\perp$, and the perpendicular weight index $n=1$. Here $L_l^{n}$ is the **associated Laguerre polynomial**, and the argument $\hat v_\perp^2 = \mu$ is the (normalized) magnetic-moment / energy variable. The $n=1$ weight absorbs the cylindrical Jacobian $v_\perp\,dv_\perp$, so the $\{\Gamma_l^1\}$ are orthonormal under the gyrotropic measure. This is the velocity-space basis of gyrokinetic theory (the "Laguerre–Hermite pseudo-spectral" formulation, [Mandell et al. 2018](https://arxiv.org/pdf/1708.04029); see also Frei/Jorge). In the gyrokinetic literature the general object is $\Gamma_n^{k}$ / $L_n^{|k|}$, where $k$ is the gyro-harmonic; for a gyrotropic (gyro-averaged) VDF only $k=0$ survives and the $n=1$ Jacobian weight is used as above.

**Why it matters:** HL with Laguerre is the *physically faithful* representation of a gyrotropic distribution (no spurious odd modes, a ring/shell is just a single low-$l$ Laguerre mode), whereas HH is simpler/uniform (same basis both axes) but pays the price of vanishing odd-$n$ modes and an oscillatory spectrum. A perpendicular Maxwellian is $\Gamma_0^1$ exactly (one HL mode) but requires many even Hermite modes in HH.

---

## 2. The coefficients $c_{mn}$ and their physical meaning

The spectral coefficient is the projection integral (Larosa Eq. 4; Coburn; Pezzi):

$$\boxed{\;c_{mn}=\int f(\xi_\parallel,\xi_\perp)\,\psi_m(\xi_\parallel)\,\psi_n(\xi_\perp)\,d\xi_\parallel\,d\xi_\perp\;}$$

(for HL, replace $\psi_n(\xi_\perp)\to\Gamma_l^1(\hat v_\perp)$ and integrate $d^3v$).

In practice the integral is done by **Gauss–Hermite (and Gauss–Laguerre) quadrature**: the measured VDF is interpolated onto the roots of the $(N_v{+}1)$-th Hermite polynomial *before* projecting (Larosa §2, following Servidio 2017 and Pezzi 2018). This exactly satisfies Parseval–Plancherel up to machine precision, *"avoiding spurious aliasing and convergence problems"* (Pezzi). Coburn's HL quadrature gives *"exact integrals for Hermite–Laguerre polynomials of order $2n-1$."*

**Physical hierarchy of moments** (the heart of the method — Pezzi §; Roytershteyn & Delzanno 2018):

| Index | Mode | Physical content |
|---|---|---|
| $c_{00}$ | Maxwellian amplitude | $\sim$ **density** (the local drifting Maxwellian); normalization reference |
| $m=1$ | 1st-order parallel | **bulk-flow / drift** fluctuations |
| $m=2$ | 2nd-order | **temperature** deformations (pressure anisotropy) |
| $m=3$ | 3rd-order | **heat flux** perturbations |
| low $m,n$ | — | fluid moments → "fluid–kinetic coupling is an intrinsic feature" (R&D 2018) |
| high $m,n$ | — | **fine-scale, non-thermal velocity-space structure** (beams, hammerheads, phase-space filamentation) |

Pezzi states it crisply: *"the index $m$ roughly corresponds to an order of the velocity moments: $m=1$ corresponds to bulk-flow fluctuations; $m=2$ to temperature deformations; $m=3$ to heat-flux perturbations, and so on."* Because the basis is centered on the local Maxwellian, **all the physics lives in $m,n\ge 1$**: a highly deformed $f$ produces large high-order coefficients.

**Enstrophy / free energy.** The summed squared coefficients give the velocity-space **enstrophy** (Pezzi Eq. 5; Larosa Eq. 7):

$$\Omega=\int \delta f^2\,d^3v=\sum_{m>0}\big[c_m\big]^2,$$

which is the cascading quantity. It is tied to a non-Maxwellianity indicator $\epsilon$ via $\Omega=\epsilon^2 n^2$, and is *"essentially the free energy in gyrokinetics"* (Pezzi, citing Schekochihin). Larosa pairs this with a Kaufmann–Paterson entropy measure $M_{KP}=(s_M-s)/(\tfrac32 k_B)$ to connect v-space structure to departure from LTE.

---

## 3. The 2D "spectrum", directional cuts, and the anisotropic velocity-space cascade

### 3.1 The 2D spectrum $\log(c_{mn}^2/c_{00}^2)$

Plotting $\log_{10}(c_{mn}^2/c_{00}^2)$ over the integer lattice $(m,n)$ is the velocity-space analogue of a 2D power spectrum (Pezzi Fig. 2; Larosa Fig. 3). Normalizing by $c_{00}^2$ (the Maxwellian power) makes it dimensionless and measures *fractional* departure from Maxwellian. The "DC term" $c_{00}$ is excluded from the cascade sum (it is the equilibrium, not a fluctuation).

### 3.2 Parallel and perpendicular cuts, and 1D shell-averaging

Two complementary reductions are used:

- **Directional cuts:** the parallel cut $c_{m0}^2/c_{00}^2$ (perpendicular index fixed at 0) measures structure *along* $\mathbf B$; the perpendicular cut $c_{0n}^2/c_{00}^2$ measures structure *across* $\mathbf B$. Pezzi forms reduced 1D spectra $P(m_z)$ and $P(m_\perp)$ this way (integrating over the orthogonal indices), and finds **different slopes in the two directions** (below).
- **Isotropic shell-averaged spectrum** (Servidio 2017; Pezzi; Larosa Eq. for $P(m)$): sum over concentric shells of unit width in Hermite space,
$$P(m)=\!\!\sum_{m-\frac12<|\mathbf m'|\le m+\frac12}\!\! c_{\mathbf m'}^2,\qquad \mathbf m=(m_\perp,m_\parallel),\ \ m=\sqrt{m_\perp^2+m_\parallel^2}.$$
This is the direct analogue of the omnidirectional $k$-spectrum in fluid turbulence.

### 3.3 Anisotropic cascade: parallel more extended than perpendicular

The central physical result (Pezzi 2018, confirming the Servidio 2017 phenomenology): **a background magnetic field makes the velocity-space cascade anisotropic — it extends to higher orders (shallower slope) along $\mathbf B$ and is suppressed (steeper slope) across $\mathbf B$.** In Pezzi's words: *"an anisotropy is revealed when considering the direction of $\mathbf B_0$ … spectra are stretched in the parallel direction … velocity gradients are stronger along the mean field and the cascade is inhibited across the mean field,"* which they call *analogous to the Shebalin effect* (the spatial-anisotropy analogue under a strong mean field). Mechanistically the parallel extension reflects **parallel phase mixing / Landau resonances** generating fine $v_\parallel$ structure.

**Quantitatively (Pezzi 2018, Fig. 3, the verified numbers):**
- Parallel reduced spectrum: $P(m_z)\sim m_z^{-2.01}$ — consistent with the magnetized prediction $m^{-2}$.
- Perpendicular reduced spectrum: $P(m_\perp)\sim m_\perp^{-3.5}$ — *"much lower in energy and steeper,"* confirming suppressed perpendicular cascade.

This shallow-parallel / steep-perpendicular pattern is exactly the qualitative behavior your "$m^{-1.3}$ (parallel), $n^{-3.1}$ (perpendicular)" pair encodes — parallel index small (extended cascade), perpendicular index large (suppressed). I could not source those *exact* HH values, but the inequality $|{\rm slope}_\parallel|<|{\rm slope}_\perp|$ is the robust, citable result.

### 3.4 The Servidio/Larosa spectral exponents and what they diagnose

The shell-averaged exponent is *physically diagnostic of the dominant process* (Servidio 2017; Larosa §1):
- $P(m)\sim m^{-3/2}$ — **weakly magnetized / isotropic** phase-space cascade, or a **phase-mixing / electric-field-dominated** regime.
- $P(m)\sim m^{-2}$ — **magnetized / magnetic-field-dominated** regime; a steeper exponent indicates **nonlinear-dominated** dynamics (cf. Parker et al. on ITG drift-kinetic turbulence; Celebre et al. on Vlasov–Poisson).

Larosa et al. 2025 measure $P(m)\approx m^{-2}$ over $4\le m<12$ for both their *wave* ($\sim28\,R_\odot$, superalfvénic, hammerhead beams) and *turbulent* ($\sim11\,R_\odot$, subalfvénic) streams — *"consistent with expectations for a low-$\beta$ plasma (Servidio et al. 2017),"* while cautioning the slope is *"very sensitive to the chosen range of $m$."* Notably their magnetic-field PSD slopes mirror the v-space story: the wave stream has a steeper $\sim f^{-2}$ magnetic spectrum, the turbulent stream $\sim f^{-3/2}$ — a real-space/v-space "dual cascade" correspondence.

### 3.5 Low-order vs high-order: beams vs cascade (Larosa's key contrast)

- **Low $m$ ($m=1,2,3$):** the *wave* stream has *more* power here than the turbulent stream — the dense proton **beam + hammerhead** structure projects onto just the first few Hermite modes. Larosa note synthetic tests: *"a dense beam at a few thermal speeds, and a much less dense beam at $\gtrsim5$ thermal speeds provide power only to the first few Hermite modes."* (A similar $1\le m\le4$ signature was seen in Solar Orbiter reconnection-exhaust VDFs.)
- **High $m$ ($8\lesssim m\lesssim15$):** the *turbulent* stream has *more* power — *"interpreted as an effect of the stronger velocity-space cascade … which induces fine-scale VDF distortions."* Its more-Maxwellian core also gives a larger $c_{00}$.

### 3.6 Spectral flattening, the noise floor, truncation, and reconstruction

At high order the measured spectrum **flattens to a plateau** — this is *instrumental noise / interpolation error*, not physics:
- **Larosa 2025:** *"we purposefully did not discuss the spectra for $m\gtrsim20$ because the flattening is likely due to interpolation errors"*; they validate that the signal lies above one-count and bi-Maxwellian noise estimates (their Appendix B).
- **Coburn 2024 (the explicit denoising application):** the HL coefficients *"decrease in relative strength until $m,l\approx14$ where the spectrum flattens,"* identified as the noise floor. They then **truncate**: zero all coefficients beyond the floor and **reconstruct** a *"low-pass-filtered VDF"* — a smooth, analytic, differentiable VDF free of instrumental shot noise. This is the spectral method's superpower: because the basis is analytic, the **denoised VDF and its velocity derivatives** can be computed exactly, enabling wave–particle-interaction analyses (Coburn; cf. Bowen et al. 2022). The grid size ($60\times60$ for Coburn; $50\times50$ for Larosa) is chosen so the *low-order* spectrum is converged and unaffected by the truncation.

This is the analogue of low-pass filtering a Fourier spectrum: keep the resolved large-scale modes, discard the high-"wavenumber" noise.

### 3.7 The local spectral slope $\alpha_m$

A single global power-law fit is fragile (slope depends on the fit range — Larosa explicitly warn of this). A **local slope** estimated mode-by-mode is more robust. Your stated estimator,

$$\alpha_m=\frac{c_{m+2,0}^2-c_{m,0}^2}{2\,c_{00}^2},$$

is a **two-mode-spaced, $c_{00}$-normalized finite-difference of the parallel-cut power**. The "$+2$" spacing is the natural step when same-parity modes must be compared (recall HH gyrotropy kills odd perpendicular modes; and even/odd Hermite modes can carry systematically different power, producing the oscillatory spectrum Larosa describe — so differencing $m\to m{+}2$ skips the parity oscillation). It returns a *local decrement* of velocity-space power along $\mathbf B$, i.e. a discrete proxy for $dP/dm$ that flags where the cascade steepens, where it transitions into the noise plateau, and how anisotropic it is relative to the perpendicular cut. I did not find this exact formula in the four papers — it reads as a sensible local-exponent/flux estimator built on their framework (Servidio/Pezzi shell spectra, Larosa parallel/perpendicular cuts) rather than a quoted equation; use it as your own diagnostic, cross-checked against the shell-averaged $P(m)$.

---

## 4. Paper-by-paper summary

### Servidio et al. 2017 — *PRL* 119, 205101 ([arXiv:1707.08180](https://arxiv.org/abs/1707.08180))
First **observational** detection of a velocity-space cascade. Using MMS high-resolution ion VDFs in Earth's turbulent magnetosheath, they apply a **3D Hermite transform** and find a broad-band **power-law** distribution of Hermite moments $P(m)$ — the v-space analogue of an inertial range. They propose a **Kolmogorov-like phenomenology** with two regimes: $P(m)\sim m^{-3/2}$ when weakly magnetized/isotropic, and $P(m)\sim m^{-2}$ when more strongly magnetized. The cascade is **anisotropic** with respect to $\mathbf B$ and **intermittent** (correlated with coherent structures). This paper sets the index conventions, the shell-averaging, and the phenomenology all later works inherit.

### Pezzi et al. 2018 — *Phys. Plasmas* 25, 060704 ([arXiv:1803.01633](https://arxiv.org/abs/1803.01633))
The **simulation** confirmation, using a 2.5D–3V **Hybrid Vlasov–Maxwell** run ($\delta B/B_0=1/3$, $\beta=0.5$, $\mathbf B_0\parallel\hat z$). Defines the 1D Hermite basis (their Eq. 4), $f=\sum_m f_m\psi_m$, coefficients $f_m=\int f\psi_m\,d^3v$, and enstrophy $\Omega=\sum_{m>0}f_m^2$ (Eq. 5). With $N_m=100$ modes/direction and an ensemble of $64^3$ VDFs they show: (1) a **magnetized anisotropic** v-space cascade — parallel $P(m_z)\sim m_z^{-2.01}$, perpendicular $P(m_\perp)\sim m_\perp^{-3.5}$; (2) the energy distribution follows the predicted $P(m)\sim m^{-2}$; (3) the cascade is **intermittent in real space**, enhanced at reconnecting current sheets (their "Hermite spectrogram" $P(\mathbf x,m)$). A spectral break near $m\simeq15$ marks numerical dissipation. Establishes the parallel-shallower / perpendicular-steeper anisotropy and the Shebalin-effect analogy.

### Coburn et al. 2024 — *ApJ* 964, 100 ([IOP](https://iopscience.iop.org/article/10.3847/1538-4357/ad1329))
Brings the spectral method to **electrons** and to **Solar Orbiter** (SWA), and uses it as a **denoising / analytic-VDF tool** rather than (only) a cascade diagnostic. Introduces the **Hermite–Laguerre** basis with **cylindrical gyrotropy**: $\psi_m(\hat v_\parallel)$ in parallel, **associated Laguerre $\Gamma_l^1(\hat v_\perp)=e^{-\hat v_\perp^2/2}L_l^1(\hat v_\perp^2)$** in perpendicular ($\mu=v_\perp^2$). On a $60\times60$ grid, coefficients fall off until a noise floor at $m,l\approx14$; **truncating there** yields a smooth, differentiable, noise-free electron VDF. They then study the regulation of the solar-wind **electron heat flux** by wave–particle interactions, exploiting the analytic VDF to compute velocity-space gradients.

### Larosa et al. 2025 — [arXiv:2512.01492](https://arxiv.org/abs/2512.01492) (Parker Solar Probe)
*Velocity-space turbulent cascade in the near-Sun solar wind.* Applies **Hermite–Hermite** decomposition with **Cartesian-like gyrotropy** ($f(-v_\perp)=f(v_\perp)$, perpendicular axis mirrored → odd perpendicular modes vanish, oscillatory spectrum) to **PSP SPAN-i** proton VDFs, $N_v=50$ ($50^2=2500$ Hermite points, matching the $8\times32\times8=2048$ instrument grid). Two streams: a superalfvénic **wave stream** ($\sim28\,R_\odot$, hammerhead beams, ion-cyclotron wave activity, low-$m$ dominated) and a subalfvénic **turbulent stream** ($\sim11\,R_\odot$, more Maxwellian core, more high-$m$ power). Both give shell-spectrum slope $\approx-2$ over $4\le m<12$. The Hermite analysis is cross-compared with field/VDF energy-conversion proxies (Local Energy Transfer, PVI, kurtosis, Kaufmann–Paterson $M_{KP}$), supporting a **dual cascade in real and velocity space**. Flattening at $m\gtrsim20$ is discarded as interpolation noise.
*(This is the paper matching your "Larosa et al. 2025 / HH / Cartesian-like / mirrored / $50\times50$" description — note it is **PSP**, not Solar Orbiter; the $60\times60$ grid in your prompt belongs to Coburn's HL analysis.)*

---

## 5. One-paragraph synthesis

A plasma VDF is decomposed onto an **orthonormal Hermite (parallel) and Hermite or associated-Laguerre (perpendicular) basis** centered on the local drifting Maxwellian — the velocity-space analogue of a Fourier transform, with $\psi_m=H_m(\xi)e^{-\xi^2/2}/\sqrt{2^m m!\sqrt\pi}$ and (cylindrical) $\Gamma_l^1=e^{-\hat v_\perp^2/2}L_l^1(\hat v_\perp^2)$. The coefficients $c_{mn}=\int f\,\psi_m\psi_n\,d\xi$ form a hierarchy: $c_{00}\!\sim$ density, $m=1,2,3\!\sim$ drift/temperature/heat-flux, high $m,n\!\sim$ fine non-thermal structure. The squared, $c_{00}^2$-normalized coefficients define a velocity-space "power spectrum" whose **shell-averaged** form is a power law $P(m)\sim m^{-3/2}$ (weak field) to $m^{-2}$ (magnetized), and whose **directional cuts** reveal an **anisotropic cascade — shallower/more extended along $\mathbf B$ (parallel phase mixing), steeper/suppressed across $\mathbf B$** (Pezzi: $m_z^{-2.01}$ vs $m_\perp^{-3.5}$). At high order the spectrum **flattens into an instrumental noise floor**, motivating **spectral truncation and reconstruction** of denoised, analytic, differentiable VDFs (Coburn). Local slope estimators (your $\alpha_m$) sharpen the measurement of where the cascade steepens and transitions to noise.

---

## Sources
- [Servidio et al. 2017, *PRL* 119, 205101 — APS](https://link.aps.org/doi/10.1103/PhysRevLett.119.205101) · [arXiv:1707.08180](https://arxiv.org/abs/1707.08180)
- [Pezzi et al. 2018, *Phys. Plasmas* 25, 060704 — AIP](https://pubs.aip.org/aip/pop/article/25/6/060704/320076/Velocity-space-cascade-in-magnetized-plasmas) · [arXiv:1803.01633](https://arxiv.org/abs/1803.01633)
- [Coburn et al. 2024, *ApJ* 964, 100 — IOPscience](https://iopscience.iop.org/article/10.3847/1538-4357/ad1329) · [ADS](https://ui.adsabs.harvard.edu/abs/2024ApJ...964..100C)
- [Larosa et al. 2025, "Velocity-space turbulent cascade in the near-Sun solar wind" (PSP) — arXiv:2512.01492](https://arxiv.org/abs/2512.01492) ([HTML](https://arxiv.org/html/2512.01492))
- [Roytershteyn & Delzanno 2018, *Front. Astron. Space Sci.* 5, 27](https://www.frontiersin.org/journals/astronomy-and-space-sciences/articles/10.3389/fspas.2018.00027/full)
- [Mandell et al. 2018, "Laguerre–Hermite pseudo-spectral velocity formulation of gyrokinetics" — arXiv:1708.04029](https://arxiv.org/pdf/1708.04029)

**Two flagged unverifiable specifics from the prompt:** the exact slope pair "$m^{-1.3}$, $n^{-3.1}$" and the local-slope formula $\alpha_m=(c_{m+2,0}^2-c_{m,0}^2)/(2c_{00}^2)$ do not appear verbatim in the sources above; they are consistent with (but not directly quotable from) the anisotropic-cascade results of Pezzi 2018 and the directional cuts of Larosa 2025. Also note the **mission correction**: the Larosa et al. 2025 Hermite paper is **Parker Solar Probe**, not Solar Orbiter — the Solar Orbiter spectral work is Coburn et al. 2024 (electrons, Hermite–Laguerre, $60\times60$).

The Larosa PSP PDF I downloaded for direct reading is at `/private/tmp/claude-501/-Users-alonsocervantesflores-Programaci-n-Melvin--Velocity-distribution-function/8d7892f7-2d12-4c46-bbfd-28b7902d9351/scratchpad/larosa_psp.pdf`.

# ===== quadrature-numerics =====

# Numerical Machinery for Gauss–Hermite / Gauss–Laguerre Quadrature and Hermite Decomposition of a VDF

---

## 1. Gaussian quadrature and why the nodes must be the polynomial roots

### 1.1 Setup

We approximate weighted integrals by an $N$-point rule

$$
\int_a^b p(x)\,f(x)\,dx \;\approx\; \sum_{i=1}^{N} w_i\,f(x_i),
$$

where $p(x)\ge 0$ is a fixed **weight function**, $\{x_i\}$ are the **nodes**, and $\{w_i\}$ are the **weights**. The two cases relevant to a velocity distribution function (VDF) decomposition are

| Family | Weight $p(x)$ | Interval | Orthogonal polynomial $P_N$ |
|---|---|---|---|
| Gauss–**Hermite** | $e^{-x^2}$ | $(-\infty,\infty)$ | Hermite $H_N$ |
| Gauss–**Laguerre** | $e^{-x}$ | $(0,\infty)$ | Laguerre $L_N$ |

An $N$-point rule has $2N$ free parameters ($N$ nodes + $N$ weights), so one can hope to make it **exact for all polynomials of degree $\le 2N-1$**. Gaussian quadrature achieves exactly this maximal order, and it does so precisely when the nodes are the $N$ roots of $P_N$.

### 1.2 Proof that the nodes must be the roots of $P_N$ (orthogonality / division argument)

Let $\{P_k\}_{k\ge0}$ be the family of polynomials orthogonal with respect to $p$:

$$
\langle P_j,P_k\rangle \;=\; \int_a^b p(x)\,P_j(x)P_k(x)\,dx \;=\; h_k\,\delta_{jk}.
$$

Take the nodes $x_1,\dots,x_N$ to be the $N$ (real, simple, interior) roots of $P_N$. Let $f$ be **any** polynomial of degree $\le 2N-1$. Divide $f$ by $P_N$:

$$
f(x) \;=\; q(x)\,P_N(x) \;+\; r(x), \qquad \deg q \le N-1,\quad \deg r \le N-1.
$$

Integrate against $p$:

$$
\int_a^b p\,f\,dx \;=\; \underbrace{\int_a^b p\,q\,P_N\,dx}_{=\,0} \;+\; \int_a^b p\,r\,dx .
$$

The first term vanishes by **orthogonality**: $q$ has degree $\le N-1$, so it is a linear combination of $P_0,\dots,P_{N-1}$, each orthogonal to $P_N$. Hence

$$
\int_a^b p\,f\,dx \;=\; \int_a^b p\,r\,dx .
$$

Now evaluate the quadrature sum. Because each node is a root, $P_N(x_i)=0$, so $f(x_i)=q(x_i)P_N(x_i)+r(x_i)=r(x_i)$:

$$
\sum_{i=1}^N w_i f(x_i) \;=\; \sum_{i=1}^N w_i r(x_i).
$$

If the weights are chosen so the rule integrates **every polynomial of degree $\le N-1$ exactly** (always possible with $N$ nodes — see §2), then $\sum_i w_i r(x_i)=\int_a^b p\,r\,dx$. Combining,

$$
\sum_{i=1}^N w_i f(x_i)=\int_a^b p\,r\,dx=\int_a^b p\,f\,dx .
$$

So the rule is exact up to degree $2N-1$. **The choice $x_i=\text{roots of }P_N$ is what makes the high-degree remainder $q\,P_N$ integrate to zero** — any other node set leaves a nonzero contribution and loses the extra $N$ orders. The roots are guaranteed real, simple, and inside $(a,b)$ precisely because $P_N$ is orthogonal with respect to a positive weight (a standard theorem of orthogonal polynomials). See [Wikipedia: Gauss–Legendre quadrature](https://en.wikipedia.org/wiki/Gauss%E2%80%93Legendre_quadrature) and the ACME lab notes [Gaussian Quadrature](https://acme.byu.edu/0000017a-17ef-d8b9-adfe-77ef21000001/vol2a-gaussianquadrature-fall2016-pdf).

---

## 2. Node–weight formula (Lagrange interpolation) and the Christoffel–Darboux closed form

### 2.1 Weights as moments of the Lagrange basis

Given the nodes, demand exactness for all polynomials of degree $\le N-1$. Build the **Lagrange basis**

$$
\ell_i(x) \;=\; \prod_{j\ne i}\frac{x-x_j}{x_i-x_j}, \qquad \ell_i(x_j)=\delta_{ij}, \quad \deg \ell_i = N-1 .
$$

Any degree-$\le N-1$ polynomial equals $\sum_i f(x_i)\ell_i(x)$. Integrating exactly forces

$$
\boxed{\,w_i \;=\; \int_a^b p(x)\,\ell_i(x)\,dx\,}.
$$

This is the defining "interpolatory" weight; with Gaussian nodes it automatically also gives exactness up to $2N-1$ by §1.2.

### 2.2 Reduction using $P_N$

Because $x_1,\dots,x_N$ are the roots of $P_N$,

$$
\prod_{j\ne i}(x-x_j) \;=\; \frac{P_N(x)}{(x-x_i)\,k_N}\cdot k_N\;=\;\frac{P_N(x)}{(x-x_i)}\cdot\frac{1}{c}, 
$$

and using $P_N'(x_i)=k_N\prod_{j\ne i}(x_i-x_j)$ (derivative of a product, only the surviving term) one obtains the standard closed form

$$
\ell_i(x)=\frac{P_N(x)}{P_N'(x_i)\,(x-x_i)}
\quad\Longrightarrow\quad
\boxed{\,w_i \;=\; \frac{1}{P_N'(x_i)}\int_a^b \frac{p(x)\,P_N(x)}{x-x_i}\,dx\,}.
$$

### 2.3 Christoffel–Darboux closed form

The Christoffel–Darboux identity for orthonormal $\{\hat P_k\}$,

$$
\sum_{k=0}^{N-1}\hat P_k(x)\hat P_k(y)
=\frac{a_N\bigl[\hat P_N(x)\hat P_{N-1}(y)-\hat P_{N-1}(x)\hat P_N(y)\bigr]}{x-y},
$$

evaluated on the kernel above collapses the integral and yields the practical weight formula

$$
\boxed{\,w_i=\frac{1}{\displaystyle\sum_{k=0}^{N-1}\hat P_k(x_i)^2}\,}
\qquad\text{(Christoffel numbers).}
$$

Specializing with the explicit derivative relations gives the textbook closed forms (Abramowitz & Stegun, Ch. 25), here for the **physicists'** Hermite and the standard Laguerre families:

$$
\textbf{Hermite (A\&S 25.4.46):}\quad
w_i=\frac{2^{\,N-1}\,N!\,\sqrt{\pi}}{N^2\,\bigl[H_{N-1}(x_i)\bigr]^2},
$$

$$
\textbf{Laguerre (A\&S 25.4.45):}\quad
w_i=\frac{x_i}{(N+1)^2\,\bigl[L_{N+1}(x_i)\bigr]^2}.
$$

The Hermite formula is confirmed in [Wolfram MathWorld: Hermite–Gauss Quadrature](https://mathworld.wolfram.com/Hermite-GaussQuadrature.html). NumPy uses the equivalent relation it documents as $w_k = c\,/\,\bigl(H'_N(x_k)\,H_{N-1}(x_k)\bigr)$ with $c$ fixed by normalizing $\int p\,dx$ ([numpy hermgauss docs](https://numpy.org/doc/stable/reference/generated/numpy.polynomial.hermite.hermgauss.html)).

---

## 3. Computing the nodes: Newton–Raphson vs. Golub–Welsch

### 3.1 Newton–Raphson root-finding

Directly solve $P_N(x_i)=0$ by iterating

$$
x^{(k+1)}=x^{(k)}-\frac{P_N(x^{(k)})}{P_N'(x^{(k)})},
$$

evaluating $P_N,P_N'$ via the three-term recurrence. **Drawbacks:** the roots cluster and the polynomial oscillates with huge dynamic range; without very good initial guesses (e.g. asymptotic estimates) Newton can **converge to an already-found root or skip a root entirely**, and overflow of $P_N$ for large $N$ is common. It works but is fragile and requires deflation / careful seeding.

### 3.2 Golub–Welsch algorithm (eigenvalues of the Jacobi matrix)

The robust modern method. The orthonormal polynomials obey a symmetric three-term recurrence (your notation)

$$
x\,P_m(x)=a_{m+1}P_{m+1}(x)+b_m P_m(x)+a_m P_{m-1}(x).
$$

Stack $m=0,\dots,N-1$: the recurrence is the matrix relation $x\,\mathbf P(x)=J\,\mathbf P(x)+a_N P_N(x)\,\mathbf e_N$, where $J$ is the **symmetric tridiagonal Jacobi matrix**

$$
J=\begin{pmatrix}
b_0 & a_1 & & \\
a_1 & b_1 & a_2 & \\
& a_2 & b_2 & \ddots \\
& & \ddots & \ddots & a_{N-1}\\
& & & a_{N-1} & b_{N-1}
\end{pmatrix}.
$$

At a node $x_i$ (root of $P_N$) the last term drops, so $J\mathbf P(x_i)=x_i\mathbf P(x_i)$:

- **Nodes** $x_i$ = **eigenvalues** of $J$ (real, simple, computed stably by the QR algorithm).
- **Weights** $w_i=\mu_0\,(v_{0,i})^2$, where $v_{0,i}$ is the **first component** of the $i$-th normalized eigenvector and $\mu_0=\int_a^b p(x)\,dx$ is the zeroth moment.

This is exactly the construction of [Golub & Welsch, *Math. Comp.* **23**(106), 221–230, 1969](https://bibbase.org/network/publication/golub-welsch-calculationofgaussquadraturerules-1969): "the eigenvalues are the sample points; the weights are $\mu_0$ times the square of the first entries of the corresponding eigenvectors," reducing quadrature construction to a symmetric tridiagonal eigenproblem solved by QR ([Reichel, *Computation of Gauss-type quadrature rules*](https://www.math.kent.edu/~reichel/publications/DC.pdf)).

### 3.3 Recurrence coefficients

Using **monic** coefficients $\alpha_m$ (diagonal $b_m$) and $\beta_m$ (off-diagonal $a_m=\sqrt{\beta_m}$):

**Hermite** (physicists', $H_{n+1}=2xH_n-2nH_{n-1}$, $\int e^{-x^2}H_mH_n=2^n n!\sqrt\pi\,\delta_{mn}$):

$$
b_m=0,\qquad a_m=\sqrt{\tfrac{m}{2}},\qquad \mu_0=\int_{-\infty}^{\infty}e^{-x^2}dx=\sqrt{\pi}.
$$

**Laguerre** ($(n+1)L_{n+1}=(2n+1-x)L_n-nL_{n-1}$, $\int_0^\infty e^{-x}L_mL_n=\delta_{mn}$):

$$
b_m=2m+1,\qquad a_m=m,\qquad \mu_0=\int_{0}^{\infty}e^{-x}dx=1.
$$

(equivalently monic $\beta_m=m^2$ for Laguerre, $\beta_m=m/2$ for Hermite).

### 3.4 Practical tools (NumPy)

```python
import numpy as np
x_H, w_H = np.polynomial.hermite.hermgauss(N)    # weight e^{-x^2} on (-inf, inf)
x_L, w_L = np.polynomial.laguerre.laggauss(N)    # weight e^{-x}   on (0, inf)
```

- [`numpy.polynomial.hermite.hermgauss(deg)`](https://numpy.org/doc/stable/reference/generated/numpy.polynomial.hermite.hermgauss.html): returns nodes/weights that "correctly integrate polynomials of degree `2*deg - 1` or less" with weight $e^{-x^2}$. Internally it builds the **companion (Jacobi) matrix**, takes its eigenvalues for the nodes (refined by one Newton step), and computes weights from the analytic formula of §2.3 — a modified Golub–Welsch for $N\lesssim 150$, with a stable asymptotic algorithm beyond.
- `numpy.polynomial.laguerre.laggauss(deg)`: the analogous routine for weight $e^{-x}$ on $(0,\infty)$.

(Note: NumPy's `hermgauss` is the *physicists'* convention $e^{-x^2}$; `hermegauss` is the *probabilists'* $e^{-x^2/2}$ — match the convention to your Hermite-function basis.)

---

## 4. Interpolation step: scattered VDF data onto the quadrature grid

The VDF is measured at **scattered points** $(\xi_{\parallel}^{(k)},\xi_{\perp}^{(k)})$ with values $f_k$, but the quadrature requires $f$ **at the node grid** $(\xi_{\parallel,i},\xi_{\perp,j})$ produced by Golub–Welsch (the nodes are fixed roots, almost never coincident with measurement locations). Use **inverse-distance weighting (IDW / Shepard interpolation)**:

$$
f(\xi_{\parallel,i},\xi_{\perp,j})
=\frac{\displaystyle\sum_k W_k\,f_k}{\displaystyle\sum_k W_k},
\qquad
W_k=\frac{1}{r_k^2+\varepsilon},
\qquad
r_k^2=(\xi_{\parallel,i}-\xi_{\parallel}^{(k)})^2+(\xi_{\perp,j}-\xi_{\perp}^{(k)})^2 .
$$

- The **convex-combination** form ($\sum_k W_k f_k/\sum_k W_k$) guarantees the interpolant lies within the data range and reproduces constants exactly.
- The regularizer $\varepsilon>0$ (a softening length) prevents the $1/r^2$ singularity when a node nearly coincides with a measurement, and controls smoothness: small $\varepsilon$ → sharp, nearly-nearest-neighbor; larger $\varepsilon$ → smoother but more diffusive. Choose $\varepsilon$ comparable to the squared local sample spacing.
- This interpolation error is the dominant error floor of the whole pipeline (see §5 convergence remark), independent of how many quadrature nodes you add.

---

## 5. The 3-step recipe and the Hermite-coefficient formula

### Recipe

1. **Nodes** $\{\xi_i\}$ (and the same set for $\perp$): eigenvalues of the Hermite Jacobi matrix via Golub–Welsch (`hermgauss(N)`).
2. **Weights** $\{w_H(\xi_i)\}$: from $\mu_0(v_{0,i})^2$ or the closed form of §2.3 (also returned by `hermgauss`).
3. **Interpolate** the measured VDF onto the tensor node grid $(\xi_i,\xi_j)$ via IDW (§4).

### Coefficient formula

Decompose the VDF in the orthonormal **Hermite-function** basis

$$
\psi_n(\xi)=\frac{1}{\sqrt{2^n n!\sqrt\pi}}\;e^{-\xi^2/2}\,H_n(\xi),
\qquad
\int_{-\infty}^{\infty}\psi_m(\xi)\psi_n(\xi)\,d\xi=\delta_{mn}
$$

(orthonormal with weight **1**, not $e^{-\xi^2}$). The projection coefficient is the plain integral

$$
c_{mn}=\iint f(\xi_\parallel,\xi_\perp)\,\psi_m(\xi_\parallel)\,\psi_n(\xi_\perp)\,d\xi_\parallel\,d\xi_\perp .
$$

To evaluate this with **Gauss–Hermite** (whose weight is $e^{-\xi^2}$) we insert $1=e^{-\xi^2}e^{+\xi^2}$ so the bracketed quantity is what the rule integrates:

$$
\int g(\xi)\,d\xi=\int e^{-\xi^2}\big[g(\xi)e^{+\xi^2}\big]d\xi\approx\sum_i w_H(\xi_i)\,e^{+\xi_i^2}\,g(\xi_i).
$$

Applied to both dimensions:

$$
\boxed{\;
c_{mn}=\sum_{i}\sum_{j}
f(\xi_i,\xi_j)\;
e^{\xi_i^2}\,w_H(\xi_i)\;
e^{\xi_j^2}\,w_H(\xi_j)\;
\psi_m(\xi_i)\,\psi_n(\xi_j)
\;}
$$

The factors $e^{\xi_i^2},\,e^{\xi_j^2}$ **exactly undo the $e^{-\xi^2}$ weight baked into the Gauss–Hermite weights**, converting the weighted quadrature into the unweighted projection integral. (Numerically, $e^{\xi^2}$ and $\psi$'s $e^{-\xi^2/2}$ partly cancel; for stability one combines them as $e^{\xi^2}\psi_m(\xi)=H_m(\xi)e^{-\xi^2/2}\cdot e^{\xi^2}/\sqrt{\cdots}$ and works with the $H_m$ form to avoid overflow at large $|\xi_i|$ — the "stable Hermite transform" idea, [arXiv:2604.02041](https://arxiv.org/pdf/2604.02041).)

### $(M,N)$ convergence: accuracy vs. cost

Two integers control the calculation:

- $M$ = number of Hermite **modes** retained ($m,n\le M$);
- $N$ = number of **quadrature nodes** per dimension.

**Quadrature requirement.** The integrand $f\,e^{\xi^2}\psi_m\psi_n$, when $f$ is itself a Hermite series of order $\le M$, is a polynomial times $e^{-\xi^2}$ of polynomial degree $\le 2M$ (the $e^{\xi^2}$ cancels the two Gaussian envelopes, leaving $\propto H_a H_b$). Gauss–Hermite is exact to degree $2N-1$, so **$N\ge M+1$ integrates every retained mode exactly.** Below that, low-order coefficients are aliased by high-order content (under-resolved quadrature).

**Why one still pushes $N$ higher.** A *measured/interpolated* VDF is **not** a finite Hermite polynomial; the integrand has tails and structure that need more nodes to resolve. As $N$ grows the coefficients $c_{mn}$ first change rapidly, then **plateau**: in practice the values **stop varying beyond $N\approx 60$**, because at that point the residual error is dominated not by quadrature but by the **IDW interpolation error** (§4), which $N$ cannot reduce. Adding nodes past the plateau only increases cost (the $N^2$ tensor grid, plus the $O(N^3)$ — or QR-accelerated — eigensolve) without improving accuracy.

**Practical balance.** Choose $M$ from the physical mode content you need to capture (drift, temperature anisotropy, skewness, heat flux ↔ low-order Hermite moments); choose $N$ just past the convergence plateau — typically $N\simeq 60$ here, comfortably satisfying $N\ge M+1$ for the modest $M$ of interest while sitting at the interpolation-limited floor. This gives the best accuracy-per-cost: smaller $N$ under-integrates, larger $N$ wastes work on error the interpolation already bounds.

---

## References

- G. H. Golub and J. H. Welsch, *Calculation of Gauss Quadrature Rules*, **Mathematics of Computation 23(106), 221–230 (1969)** — Jacobi-matrix eigenvalue construction of nodes/weights. [bibbase entry](https://bibbase.org/network/publication/golub-welsch-calculationofgaussquadraturerules-1969); algorithmic summary in [L. Reichel, *Computation of Gauss-type quadrature rules*](https://www.math.kent.edu/~reichel/publications/DC.pdf).
- M. Abramowitz and I. A. Stegun, *Handbook of Mathematical Functions*, Ch. 25 — formulas **25.4.45** (Gauss–Laguerre weights) and **25.4.46** (Gauss–Hermite weights); cross-checked at [Wolfram MathWorld: Hermite–Gauss Quadrature](https://mathworld.wolfram.com/Hermite-GaussQuadrature.html).
- NumPy: [`numpy.polynomial.hermite.hermgauss`](https://numpy.org/doc/stable/reference/generated/numpy.polynomial.hermite.hermgauss.html) (weight $e^{-x^2}$, exact to degree $2\,\text{deg}-1$, companion-matrix/Golub–Welsch nodes); `numpy.polynomial.laguerre.laggauss` (weight $e^{-x}$).
- General theory (orthogonality / roots argument, Christoffel–Darboux): [Gauss–Legendre quadrature (Wikipedia)](https://en.wikipedia.org/wiki/Gauss%E2%80%93Legendre_quadrature); [ACME Gaussian Quadrature lab](https://acme.byu.edu/0000017a-17ef-d8b9-adfe-77ef21000001/vol2a-gaussianquadrature-fall2016-pdf).
- Numerically stable Hermite transform (handling the $e^{\xi^2}\psi_n$ overflow at large nodes): [*Stable Hermite transforms via the Golub–Welsch algorithm*, arXiv:2604.02041](https://arxiv.org/pdf/2604.02041).