# Spectral Analysis of Solar Orbiter SWA-PAS Proton Velocity Distribution Functions: A Theory Reference

*A graduate-level reference for reconstructing the proton VDF in field-aligned coordinates, interpolating it onto a 60×60 Gauss–Hermite node grid, and computing the 2-D Hermite spectrum $\log_{10}(c_{mn}^2/c_{00}^2)$.*

---

## Table of contents

1. Solar wind, the SWA-PAS instrument, and the L2 VDF data product
2. From counts to phase-space density and velocity-space reconstruction
3. Plasma moments and the normalized peculiar velocity $\xi=(v-u)/w$
4. Gyrotropy and parametric VDF models (core/beam/alpha, kappa)
5. Polynomial decomposition: orthonormal Hermite functions, the Fourier analogy, $c_{mn}$, HH vs HL
6. Gaussian quadrature numerics
7. Interpolation onto the node grid (inverse-distance weighting)
8. The 2-D velocity-space spectrum, anisotropic cascade, and noise-floor flattening
9. References

> **Three attribution corrections carried from the source notes, stated once up front.**
> (i) There is **no** standalone "Livi et al. 2023 PAS in-flight calibration" paper; Livi et al. 2023 (A&A 676, A36) is the **SWA-HIS** Heavy Ion Sensor paper. PAS characterization is documented in Owen et al. (2020), the Fedorov L2 User Guide (2020), and Louarn et al. (2024).
> (ii) The **Larosa et al. 2025 Hermite–Hermite paper is Parker Solar Probe (SPAN-i), not Solar Orbiter.** The Solar Orbiter spectral work is Coburn et al. (2024) (electrons, Hermite–Laguerre, $60\times60$ grid).
> (iii) The exact parallel/perpendicular slope pair "$m^{-1.3}$, $n^{-3.1}$" and the local-slope estimator $\alpha_m=(c_{m+2,0}^2-c_{m,0}^2)/(2c_{00}^2)$ are **working conventions of this project**, not quoted from the cited literature; they are consistent in spirit with the anisotropic-cascade results of Pezzi et al. (2018) and the directional cuts of Larosa et al. (2025), and are treated as such in §8.

---

## 1. Solar wind, the SWA-PAS instrument, and the L2 VDF data product

### 1.1 Physical motivation

The solar wind is a weakly collisional, magnetized plasma. Because the proton–proton collision time is typically comparable to or longer than the expansion time, proton velocity distribution functions (VDFs) retain non-equilibrium structure — temperature anisotropy, field-aligned beams, suprathermal tails — that a fluid description averages away (Marsch 2006; Verscharen, Klein & Maruca 2019). The VDF $f(\mathbf v)$ is the phase-space density of particles; all fluid moments (density, bulk velocity, pressure/temperature tensor, heat flux) are velocity-space integrals of $f$. Recovering and decomposing the *measured* $f(\mathbf v)$ is therefore the natural way to study kinetic physics directly.

### 1.2 The SWA suite and the PAS sensor

The **Solar Wind Analyser (SWA)** suite aboard **Solar Orbiter** comprises three top-hat electrostatic analysers (EAS for electrons, PAS for protons/alphas, HIS for heavy ions) (Owen et al. 2020, §1, §3.3). The **Proton–Alpha Sensor (PAS)** is a quarter/partial-spherical top-hat ESA optimized for the narrow, strongly radial solar-wind ion beam. It measures the full 3-D VDF of solar-wind ions **without mass or charge selection** (in practice H⁺ and He²⁺) over **200 eV/q – 20 keV/q**, mounted behind a scalloped cut-out in the heat shield so that it looks roughly sunward (Owen et al. 2020, §3.3.1–3.3.2).

PAS forms its 3-D count matrix from three measurement axes:

1. **Energy-per-charge $E/q$ — the ESA itself.** A spherical-section analyser transmits only ions in a narrow $E/q$ passband set by the inner-hemisphere voltage, $E=k\,V_{\rm hemi}$ with **analyser constant $k\approx 13$–$14$ eV/V** (Owen Table 8). PAS steps through **96 logarithmically (quasi-geometric) spaced energy levels** (92 used per sample), **relative resolution $\Delta E/E\approx 5.5\%$ design** (3.0–9.3% measured).
2. **Elevation (polar angle) — entrance deflectors.** Two curved deflector plates plus a top-cap electrode steer ions by the deflector-to-analyser voltage ratio $U_{\rm def}/U_{\rm an}$. Binned into **9 elevation bins** spanning $\approx -22.5^\circ$ to $+22.5^\circ$ ($\sim 5^\circ$ resolution).
3. **Azimuth — CEM/anode array.** An array of **11 channel electron multipliers (CEMs)** read out simultaneously, covering **$-24^\circ$ to $+42^\circ$** (a $66^\circ$ fan offset $+9^\circ$ from the Sun direction to capture solar-wind aberration), $\sim 5^\circ$ resolution.

**Sweep scheme** (Owen §3.3.3; Fedorov §2.2): at fixed energy, PAS sweeps elevation while accumulating all 11 azimuth bins at once, then steps energy. A full $(96\times 9\times 11)$ matrix is acquired in $\sim 1$ s; per-element accumulation $\approx 1$ ms. Geometric factor per bin $\approx 5\times 10^{-6}\ \mathrm{cm^2\,sr\,eV/eV}$. A peak-tracking algorithm recentres a reduced energy–elevation window on the beam for fast modes.

**Cadence/modes** (Fedorov §1): **Normal** — one VDF every 4 s; **Snapshot** — $\sim 9$ s bursts every 300 s at $\sim$4 Hz; **Burst** — 300 s continuous, up to $\sim 15$–20 Hz with reduced phase space.

### 1.3 The L2 `swa-pas-vdf` CDF product

The 3-D product is the CDF dataset **`solo_L2_swa-pas-vdf`** (sister products: `swa-pas-eflux`, 1-D omni-directional differential energy flux; `swa-pas-grnd-mom`, ground-computed moments). The authoritative reference for variable names, units, and frame conventions is the **SWA-PAS L2 Data User Guide** (Fedorov 2020, V02).

**Core data variable.**

| Variable | Dim. | Units | Notes |
|---|---|---|---|
| `vdf` | **[time, 11, 9, 96]** | **$\mathrm{s^3\,m^{-6}}$** | Phase-space density $f$, with `DEPEND_1=Azimuth`, `DEPEND_2=Elevation`, `DEPEND_3=Energy`. |

**Dimension ordering (critical).** The per-record matrix is stored **[azimuth = 11, elevation = 9, energy = 96]**. Owen et al. (2020) write it conceptually as $(96,9,11)=$ energy×elevation×azimuth — *the reverse of the CDF axis order*. The CDF `DEPEND_1/2/3` order is authoritative.

**Units.** $\mathrm{s^3\,m^{-6}}$ (SI phase-space density, so that $n=\int f\,d^3v$ with $d^3v$ in $(\mathrm{m/s})^3$ yields $\mathrm{m^{-3}}$). To convert to the CGS plasma unit $\mathrm{s^3\,cm^{-6}}$, multiply by $10^{-6}$.

**Coordinate / support variables** (Fedorov 2020, Table 1).

| Variable | Dim. | Units | Meaning |
|---|---|---|---|
| `Epoch` | [time] | ns (TT2000) | Record time. |
| `Info` | [time] | — | Sampling category: 0 Ground, 1 Normal, 2 Snapshot, 3 Burst, 4 Eng., 5 Cal. **Use only Info = 1/2/3.** |
| `Energy` | [96] | eV | Bin centres; asymmetric half-widths `delta_p_Energy`, `delta_m_Energy`. |
| `Azimuth` | [11] | deg | CEM bin centres; half-width `delta_Azimuth` (symmetric). |
| `Elevation` | [9] | deg | Elevation bin centres; half-width `delta_Elevation` (symmetric). |
| `Full_azimuth` | [11, 9] | deg | Full azimuth table $az(i_{az},i_{el})$ — azimuth depends on elevation. |
| `Full_elevation` | [11, 9] | deg | Full elevation table $el(i_{az},i_{el})$. |
| `Elevation_correction` | [96] | deg | Per-energy "saw-tooth" offset; sign alternates with energy parity. |
| `PAS_to_RTN` | **[time, 3, 3]** | None | Instrument → RTN rotation matrix (§2.5). |
| `start_*`,`nb_*` | [time] | — | Window start index and count actually swept (peak-tracking subset). |

Because peak-tracking moves the active window, a record may populate only a sub-block of the $11\times 9\times 96$ array; the `start_*`/`nb_*` variables tell you which indices are valid. Both `DELTA_PLUS_VAR`/`DELTA_MINUS_VAR` are provided for energy (asymmetric, because of log spacing) and for angles (symmetric) — exactly what is needed for the velocity-space volume element (§2, §3).

---

## 2. From counts to phase-space density and velocity-space reconstruction

### 2.1 What an ESA measures, and counts → phase-space density

A top-hat ESA selects particles by energy-per-charge $E/q$, elevation $\theta$ (deflectors), and azimuth $\phi$ (anodes), accumulating a **count array** $C_{ijk}$ over an integration time $\Delta t$ (indices $i,j,k$ = energy, elevation, azimuth). For a representative instrument (PSP/SPAN-i), $E=k V_{\rm hemi}$ with $k\approx 16.7$ eV/V and an $8\times 32\times 8$ native resolution (Livi et al. 2022).

The conversion chain from counts to $f$ runs through the differential directional number flux $j_N$ (particles $\mathrm{area^{-1}\,s^{-1}\,sr^{-1}\,eV^{-1}}$). Because the flux through an area element is $d\Phi = v\,f\,d^3v$ with $d^3v=v^2\,dv\,d\Omega$ and $dv=dE/(mv)$,

$$
j_N=\frac{d\Phi}{dE\,d\Omega}=\frac{v^2}{m}f=\frac{2E}{m^2}f
\quad\Longrightarrow\quad
\boxed{\,f=\frac{m^2}{2E^2}\,j_E=\frac{m^2}{2E}\,j_N\,},
$$

where $j_E=E\,j_N$ is the differential energy flux that the L2 pipeline reports. Counts relate to flux through the geometric factor $G$ (effective area × solid angle) and efficiency $\varepsilon$:

$$
C_{ijk}=j_N\,\varepsilon_{ijk}\,G_{ijk}\,\Delta E_i\,\Delta t
\quad\Longrightarrow\quad
f_{ijk}=\frac{m^2}{2E_i}\,\frac{C_{ijk}}{\varepsilon_{ijk}\,G_{ijk}\,\Delta E_i\,\Delta t}.
$$

This $f\propto E^{-2}\times(\text{energy flux})$ inversion is the standard ESA reduction (Livi et al. 2022). For Solar Orbiter the L2 `swa-pas-vdf` product already delivers $f$ in $\mathrm{s^3\,m^{-6}}$, so the user starts from $f_{ijk}$ directly.

### 2.2 Speed from energy

The ESA selects $E/q$; with `Energy` in eV the speed magnitude is

$$
|\mathbf v|=\sqrt{\frac{2qE}{m}}=\sqrt{\frac{2E_{\rm eV}\,e}{m}} .
$$

For **protons** ($m=m_p,\ q=e$) this reduces to the User-Guide constant **`E2V = 13.85`** (Fedorov 2020):

$$
|\mathbf v|\,[\mathrm{km/s}]=13.85\,\sqrt{E_{\rm eV}},\qquad
\sqrt{2e/m_p}=1.385\times 10^{4}\ \mathrm{m\,s^{-1}\,eV^{-1/2}} .
$$

Per-bin lower/upper speeds use the energy half-widths,

$$
v_{\min}[i_e]=13.85\sqrt{E[i_e]-\texttt{delta\_m\_Energy}},\qquad
v_{\max}[i_e]=13.85\sqrt{E[i_e]+\texttt{delta\_p\_Energy}} .
$$

(For **alphas**, the same $E/q$ step gives $q=2e,\ m=4m_p$, so $|\mathbf v|_\alpha=|\mathbf v|_p/\sqrt 2$; PAS does not mass-separate, so the default conversion assumes protons. The $m/q=2$ offset is what later allows alphas to be spectrally separated — see §4.3.)

### 2.3 Spherical → Cartesian: arrival direction vs particle velocity

A particle counted in bin $(i,j,k)$ **arrives from** the look direction $\hat{\mathbf n}_{jk}$, so its velocity points the **opposite** way. With elevation measured from the instrument $x$–$y$ plane and azimuth in that plane:

$$
\hat{\mathbf n}_{jk}=\big(\cos\theta_j\cos\phi_k,\ \cos\theta_j\sin\phi_k,\ \sin\theta_j\big),
\qquad
\mathbf v_{ijk}=-\,v_i\,\hat{\mathbf n}_{jk}.
$$

(If a polar angle $\vartheta$ from the $z$-axis is used, replace $\cos\theta\to\sin\vartheta,\ \sin\theta\to\cos\vartheta$.)

**Two conventions, one physical requirement.**

*(a) Standard-spherical (all-positive) convention* gives the **look/arrival** unit vector $\hat L_x=\cos(el)\cos(az),\ \hat L_y=\sin(az)\cos(el),\ \hat L_z=\sin(el)$ — the direction from which the ion arrived, **not** the particle velocity.

*(b) SWA convention (with minus signs)* folds $\mathbf v=-|\mathbf v|\hat L$ directly into the unit vectors (Fedorov 2020, §2.1). In the PAS analyser frame:

$$
\hat V^{\rm PAS}_X=-\cos(El)\cos(Az),\quad
\hat V^{\rm PAS}_Y=-\cos(El)\sin(Az),\quad
\hat V^{\rm PAS}_Z=+\sin(El).
$$

In the spacecraft (SRF) frame (PAS→SRF is a $180^\circ$ rotation about the sun-pointing $X$ axis: $X$ kept, $Y$ and $Z$ flipped):

$$
\boxed{\;\hat V^{\rm SRF}_X=-\cos(El)\cos(Az),\quad
\hat V^{\rm SRF}_Y=+\cos(El)\sin(Az),\quad
\hat V^{\rm SRF}_Z=-\sin(El)\;}
$$

The full velocity vector is $\mathbf v=|\mathbf v|\hat V=13.85\sqrt{E_{\rm eV}}\,(\hat V_X,\hat V_Y,\hat V_Z)\ [\mathrm{km/s}]$.

**Why the minus signs are mandatory.** For a nominal beam ($Az\approx 0,\ El\approx 0$) the SWA convention yields $\hat V^{\rm SRF}\approx(-1,0,0)$, i.e. ion velocity pointing **anti-sunward** ($-X_{\rm SRF}$), as physically required for outward-flowing solar wind. Combined with the approximate relation $X_{\rm RTN}\approx-X_{\rm SC},\ Y_{\rm RTN}\approx-Y_{\rm SC}$ and the exact per-record `PAS_to_RTN`, the radial component comes out $V_R\approx +|v|>0$ (e.g. $+400$ km/s). The all-positive formula instead returns the arrival direction ($\approx +X_{\rm SRF}$, sunward) and a spurious **negative** $V_R$.

**Per-bin corner construction** (for accurate moment cells; Fedorov §4). Each velocity-space cell is built from the corner angles `Full_azimuth`/`Full_elevation` (min, centre, max per $(i_{az},i_{el})$) and the per-energy `Elevation_correction`:

```
x_i = -cos(elArr_i) * cos(azArr_i)
y_i =  cos(elArr_i) * sin(azArr_i)
z_i = -sin(elArr_i)
if ((ie - se) % 2 == 0):  z_i += Elevation_correction[ie]
else:                     z_i -= Elevation_correction[ie]
P_i = [x_i, y_i, z_i]   # then scale by v_min / v_max for the cell corners
```

The parity-dependent `Elevation_correction` produces the characteristic "saw-tooth" bin boundaries in the $V_x$–$V_z$ plane.

### 2.4 The velocity-space Jacobian

The spherical → Cartesian map $(v,\theta,\phi)\to(v_x,v_y,v_z)$ carries the volume element

$$
\boxed{\,d^3v=v^2\,dv\,d\Omega\,},\qquad
d\Omega=\cos\theta\,d\theta\,d\phi\ \ (\text{elevation convention}),
$$

(or $d\Omega=\sin\vartheta\,d\vartheta\,d\phi$ for a polar angle). The $v^2$ factor is the radial Jacobian, $\cos\theta$ the angular Jacobian; both are essential and a frequent source of error if dropped (Wilson et al. 2025; Broderick/Klein 2025).

### 2.5 Rotation into RTN, then field-aligned coordinates

Each Cartesian velocity is rotated instrument → spacecraft (fixed mounting matrix $\mathsf R_{\rm I\to S}$) → RTN (time-dependent attitude, from SPICE):

$$
\mathbf v^{\rm RTN}=\mathsf R_{\rm S\to RTN}(t)\,\mathsf R_{\rm I\to S}\,\mathbf v^{\rm I}.
$$

For Solar Orbiter this is encapsulated in the per-record `PAS_to_RTN` matrix $M$ (dim [time,3,3]; `COORDINATE_SYSTEM = SOLO_SWA_PAS`, `TARGET_SYSTEM = SOLO_SUN_RTN`), mapping the PAS frame to heliocentric **RTN** (R radially outward, T $\approx$ along orbital motion, N completing the triad):

$$
\begin{pmatrix}v_R\\v_T\\v_N\end{pmatrix}
= M\begin{pmatrix}v_X^{\rm PAS}\\v_Y^{\rm PAS}\\v_Z^{\rm PAS}\end{pmatrix}.
$$

$M$ is a proper rotation (orthonormal, $\det=+1$), so RTN→PAS $=M^{\mathsf T}$. It is time-dependent because attitude relative to RTN varies along the orbit; the approximate sign relation $X_{\rm RTN}\approx-X_{\rm SRF},\ Y_{\rm RTN}\approx-Y_{\rm SRF}$ holds only "most of the time," so use `PAS_to_RTN` for the exact transform. Because rotations preserve $d^3v$ (orthogonal, $\det=1$), moments are frame-independent. Each measured $(i,j,k)$ bin thus becomes one point $\big(\mathbf v^{\rm RTN}_{ijk},\,f_{ijk}\big)$ — the **3-D velocity-space scatter** $(v_R,v_T,v_N)$ that one scatter-plots.

**Field-aligned (FA) triad.** Given $\mathbf B$ in RTN, set $\hat{\mathbf b}=\mathbf B/|\mathbf B|$. The parallel velocity is the projection $v_\parallel=\mathbf v\cdot\hat{\mathbf b}$. The perpendicular pair is **non-unique** (any gyrophase rotation about $\hat{\mathbf b}$ is valid); one fixes it with a reference vector $\hat{\mathbf r}$ (commonly the RTN-$\hat R$ direction, the bulk velocity $\mathbf U$, or $\mathbf B\times\mathbf U$):

$$
\hat{\mathbf e}_{\perp 2}=\frac{\hat{\mathbf b}\times\hat{\mathbf r}}{|\hat{\mathbf b}\times\hat{\mathbf r}|},\qquad
\hat{\mathbf e}_{\perp 1}=\hat{\mathbf e}_{\perp 2}\times\hat{\mathbf b},\qquad
(\hat{\mathbf e}_{\perp 1},\hat{\mathbf e}_{\perp 2},\hat{\mathbf b})\ \text{right-handed},
$$

with $v_{\perp 1}=\mathbf v\cdot\hat{\mathbf e}_{\perp 1}$, $v_{\perp 2}=\mathbf v\cdot\hat{\mathbf e}_{\perp 2}$, and gyrophase $\psi=\arctan2(v_{\perp 2},v_{\perp 1})$.

---

## 3. Plasma moments and the normalized peculiar velocity

### 3.1 Continuous moment definitions

Per species $s$ (Verscharen, Klein & Maruca 2019, Eqs. 27–29):

$$
n=\int f\,d^3v,\qquad
\mathbf U=\frac1n\int \mathbf v\,f\,d^3v,\qquad
\mathsf P=m\int(\mathbf v-\mathbf U)(\mathbf v-\mathbf U)\,f\,d^3v,
$$

with temperature tensor $\mathsf T=\mathsf P/(n k_B)$ and heat flux $\mathbf q=\tfrac{m}{2}\int(\mathbf v-\mathbf U)|\mathbf v-\mathbf U|^2 f\,d^3v$. Numerically convenient is the raw-second-moment form $\mathsf P=m\!\int\!\mathbf v\mathbf v\,f\,d^3v-mn\,\mathbf U\mathbf U$.

### 3.2 Parallel / perpendicular temperatures and derived quantities

$$
T_\parallel=\hat{\mathbf b}\cdot\mathsf T\cdot\hat{\mathbf b}=\frac{p_\parallel}{n k_B},\qquad
T_\perp=\frac{\operatorname{Tr}\mathsf T-T_\parallel}{2}=\frac{p_\perp}{n k_B},
$$

(Verscharen et al. 2019, Eqs. 38–39), with scalar temperature $T=\tfrac13\operatorname{Tr}\mathsf T=\tfrac13(T_\parallel+2T_\perp)$, anisotropy $R=A=T_\perp/T_\parallel$, and plasma beta $\beta_s=8\pi n_s k_B T_s/B^2$ (their Eq. 15).

### 3.3 Discrete-sum moments (the SWA-PAS recipe)

Replace $d^3v\to v_i^2\,\Delta v_i\,\Delta\Omega_{jk}$ in each measured bin:

$$
\boxed{\,n=\sum_{i,j,k} f_{ijk}\;v_i^2\,\Delta v_i\,\Delta\Omega_{jk}\,},
$$
$$
n\,\mathbf U=\sum_{i,j,k}\mathbf v_{ijk}\,f_{ijk}\;v_i^2\,\Delta v_i\,\Delta\Omega_{jk},\qquad
\mathbf v_{ijk}=-v_i\hat{\mathbf n}_{jk},
$$
$$
\mathsf P=m\sum_{i,j,k}(\mathbf v_{ijk}-\mathbf U)(\mathbf v_{ijk}-\mathbf U)\,f_{ijk}\;v_i^2\,\Delta v_i\,\Delta\Omega_{jk}.
$$

**Velocity-bin width from energy-bin width.** From $E=\tfrac12 mv^2\Rightarrow dE=mv\,dv$:

$$
\Delta v_i=\frac{\Delta E_i}{m\,v_i}.
$$

For logarithmically (geometrically) spaced channels with constant $\Delta E/E$ it is cleaner to use $\dfrac{dv}{v}=\tfrac12\dfrac{dE}{E}$, giving

$$
\Delta v_i=\tfrac12\,v_i\,\frac{\Delta E_i}{E_i},
\qquad\text{so}\qquad
v_i^2\,\Delta v_i=\tfrac12\,v_i^3\,\frac{\Delta E_i}{E_i}.
$$

Equivalently, in the User-Guide corner form $\int v^2\,dv=(v_{\max}^3-v_{\min}^3)/3$ with $v_{\min/\max}=13.85\sqrt{E\mp\delta E}$.

**Solid-angle element.**

$$
\Delta\Omega_{jk}=\cos\theta_j\,\Delta\theta_j\,\Delta\phi_k\ \ (\text{elevation}),
\qquad
\int_{\rm bin}d\Omega=\Delta Az_{\rm rad}\,\big[\sin(El+\delta_{El})-\sin(El-\delta_{El})\big],
$$

so explicitly (Fedorov §3), with $m=$ proton mass and angles in radians,

$$
n=\sum f\,v^2\,\Delta v\,\cos(El)\,\Delta El\,\Delta Az,\qquad
\mathbf V=\frac1n\sum f\,\mathbf v\,v^2\,\Delta v\,\cos(El)\,\Delta El\,\Delta Az,
$$
$$
\mathsf P_{jk}=m\sum f\,(v_j-V_j)(v_k-V_k)\,v^2\,\Delta v\,\cos(El)\,\Delta El\,\Delta Az.
$$

Working directly from differential energy flux gives the familiar instrument form $n=\sum (m^2/2E_i^2)\,j_{E,ijk}\,v_i^2\Delta v_i\Delta\Omega_{jk}$ (Wilson et al. 2025).

**Caveats biasing the sums** (Fedorov §3; Wilson et al. 2025): truncated FOV (partial-sky coverage), finite low-/high-energy cutoffs, spacecraft-potential shift of $E$, one-count noise. The VDF is noisy in 300–400 eV (geometric factor falling) and **unreliable below 300 eV**; sporadic "ghost" values appear at extreme angular bins; ground densities/pressures are perturbed for bulk speed 260–380 km/s and invalid below 260 km/s. In `swa-pas-grnd-mom`, use `validity` $\ge 2$.

### 3.4 Thermal speeds and the normalized peculiar velocity $\xi$

Using the **most-probable-speed** convention standard in solar-wind kinetics (Verscharen et al. 2019, Eqs. 4–5):

$$
\boxed{\,w=\sqrt{\frac{2k_B T}{m}}\,},\qquad
w_\parallel=\sqrt{\frac{2k_B T_\parallel}{m}},\qquad
w_\perp=\sqrt{\frac{2k_B T_\perp}{m}} .
$$

(Competing conventions $\sqrt{k_BT/m},\ \sqrt{3k_BT/m},\ \sqrt{8k_BT/\pi m}$ exist; the $\sqrt{2k_BT/m}$ form is the one that makes the Maxwellian exponent $-c^2/w^2$.)

The **normalized peculiar velocity** — the argument of the Hermite basis (§5) — is

$$
\boxed{\;\xi=\frac{v-u}{w}\;},\qquad
\xi_\parallel=\frac{v_\parallel-U_\parallel}{w_\parallel},\qquad
\xi_\perp=\frac{v_\perp-U_\perp}{w_\perp}.
$$

In the bulk frame the perpendicular drift is removed ($U_\perp\approx 0$, the $\mathbf E\times\mathbf B$ frame), so $\xi_\perp\simeq v_\perp/w_\perp$. A bi-Maxwellian is then $f=\dfrac{n}{\pi^{3/2}w_\parallel w_\perp^2}\exp(-\xi_\parallel^2-\xi_\perp^2)$.

A cleaner **bulk-rest-frame peculiar velocity** subtracts the flow *before* projecting (avoiding any need to define $U_\perp$ or a gyrophase):

$$
\mathbf c=\mathbf v-\mathbf U,\qquad
c_\parallel=\mathbf c\cdot\hat{\mathbf b},\qquad
\boxed{\,c_\perp=\big|\mathbf c-c_\parallel\hat{\mathbf b}\big|=\sqrt{|\mathbf c|^2-c_\parallel^2}\,}.
$$

Here $c_\perp\ge 0$ by construction and is rotation-invariant about $\hat{\mathbf b}$ — the robust quantity to histogram for a gyrotropic $f(c_\parallel,c_\perp)$ and to feed the pressure-tensor sum, with $p_\parallel=m\sum c_\parallel^2 f\,(\cdots)$ and $p_\perp=\tfrac{m}{2}\sum c_\perp^2 f\,(\cdots)$.

---

## 4. Gyrotropy and parametric VDF models

### 4.1 Gyrotropy: the 3D → 2D reduction

In a magnetized, low-$\beta$ plasma the strong background field is a symmetry axis: $f$ is invariant under rotation about $\hat{\mathbf b}$ (independent of gyrophase $\psi$). This **gyrotropy** lets one define

$$
\boxed{\,v_\perp=\sqrt{v_{\perp 1}^2+v_{\perp 2}^2}\,},\qquad
f(v_\parallel,v_\perp)=\frac{1}{2\pi}\int_0^{2\pi}f(v_\parallel,v_\perp,\psi)\,d\psi,
$$

the 2-D **gyrotropic distribution function (GDF)** targeted by Slepian-basis reconstructions precisely because $\mathbf B$ "serves as a gyrotropic axis of symmetry" in low-$\beta$ plasma (Broderick/Klein 2025). It is also the natural frame for field-aligned proton beams (Verniero et al. 2020).

**When gyrotropy fails (agyrotropy).** Near shocks, reconnection sites, and diffusion regions, gyration is interrupted within a gyroperiod and $f$ depends on $\psi$; the gyrophase average then discards real structure. Departures are quantified by scalars built from the FA pressure tensor, the most-used being **Swisdak's $Q$**:

$$
\sqrt{Q}=\sqrt{\frac{P_{12}^2+P_{13}^2+P_{23}^2}{P_\perp^2+2P_\perp P_\parallel}}\in[0,1],
$$

with axis "1" along $\hat{\mathbf b}$ (Swisdak 2016). Related diagnostics are Scudder & Daughton's $A\!\varnothing_e=2|\lambda_1-\lambda_2|/(\lambda_1+\lambda_2)$ (eigenvalues of the perpendicular $2\times2$ sub-block) and Aunai et al.'s $D_{\rm ng}$ (RMS off-diagonal $P_{ij}$ normalized to thermal energy); all measure deviation of $\mathsf P$ from the gyrotropic diagonal form $\mathrm{diag}(P_\parallel,P_\perp,P_\perp)$.

### 4.2 Maxwellian and bi-Maxwellian

A single drifting **Maxwellian** is the maximum-entropy / thermal-equilibrium solution:

$$
f_M(\mathbf v)=\frac{n}{\pi^{3/2}\,w^3}\exp\!\Big[-\frac{(\mathbf v-\mathbf u)^2}{w^2}\Big],\qquad w=\sqrt{\tfrac{2k_BT}{m}} .
$$

In a magnetized weakly collisional plasma, gyromotion thermalizes the two perpendicular directions but $T_\parallel\ne T_\perp$ in general. The **bi-Maxwellian** (Verscharen et al. 2019, Eq. 61) is

$$
f_{bM}(v_\parallel,v_\perp)=\frac{n}{\pi^{3/2}\,w_\parallel w_\perp^{2}}
\exp\!\Big[-\frac{(v_\parallel-u_\parallel)^2}{w_\parallel^{2}}\Big]
\exp\!\Big[-\frac{(v_\perp-u_\perp)^2}{w_\perp^{2}}\Big],
$$

with $w_{\parallel,\perp}=\sqrt{2k_BT_{\parallel,\perp}/m}$ and anisotropy $A=T_\perp/T_\parallel$. For a gyrotropic plasma $u_\perp=0$. The bi-Maxwellian "introduces temperature anisotropies … yet follows the Maxwellian behavior on any one-dimensional cut" (Verscharen et al. 2019). Fast-wind protons persistently show $T_{\perp p}>T_{\parallel p}$, which adiabatic expansion alone cannot produce — a kinetic signature demanding active perpendicular heating.

### 4.3 The three proton-scale populations

Solar-wind proton VDFs are well organized with respect to $\mathbf B$ and are typically modeled as a **sum of bi-Maxwellians** (Marsch et al. 1982; Marsch 2006):

- **Proton core** — global peak, $\sim 90$–95% of protons, near-(an)isotropic Maxwellian, defines the bulk velocity; densest, coolest.
- **Proton beam** — "a second proton component streaming faster than the proton core … along the direction of the magnetic field with a relative speed $\gtrsim v_{Ap}$" (Verscharen et al. 2019), drift $\Delta u\sim(1$–$2)\,v_{Ap}$ ($\approx 1.5$–$2\times$ the core thermal width along $\mathbf B$), "almost always directed away from the Sun." Modeled by adding a separate beam bi-Maxwellian ($n_b,u_b,w_{\parallel b},w_{\perp b}$).
- **Alpha particles (He²⁺)** — distinct species ($m_\alpha=4m_p$, $q=2e$, $m/q=2$), $\lesssim 20\%$ of mass density, drifting along $\mathbf B$ away from the Sun at $\lesssim v_{Ap}$. Because ESAs sort by $E/q=\tfrac12 m v^2/q$, an alpha at the *same bulk speed* as the core appears at **double the proton $E/q$**; if (mis)converted assuming proton mass it appears at $\sqrt2\times$ the proton speed, $v_{\rm app}=\sqrt{2(E/q)e/m_p}=\sqrt2\,v_p$. This $m/q=2$ offset is exactly what lets alphas be spectrally separated.

### 4.4 Kappa and regularized-kappa distributions

Collisionless plasmas show suprathermal power-law tails a Maxwellian cannot fit; Vasyliūnas (1968) and Olbert (1968) introduced the **kappa distribution**. In the Verscharen et al. (2019, Eq. 62) convention,

$$
f_\kappa(\mathbf v)=\frac{n}{[\pi(2\kappa-3)]^{3/2}\,w^{3}}\,
\frac{\Gamma(\kappa+1)}{\Gamma(\kappa-\tfrac12)}\,
\Big[1+\frac{2(\mathbf v-\mathbf u)^2}{(2\kappa-3)\,w^{2}}\Big]^{-(\kappa+1)},
$$

where the $(2\kappa-3)$ factor keeps $w$ the thermal speed with $T=mw^2/2k_B$ **independently of $\kappa$** (for $\kappa>3/2$). An equivalent "Olbertian" form with thermal parameter $\theta$:

$$
f_\kappa(v)=\frac{n}{\pi^{3/2}\theta^{3}}\,\frac{\Gamma(\kappa+1)}{\kappa^{3/2}\,\Gamma(\kappa-\tfrac12)}
\Big[1+\frac{v^2}{\kappa\theta^2}\Big]^{-(\kappa+1)} .
$$

The index $\kappa$ controls the tails: more pronounced for smaller $\kappa$, with high-energy slope $\propto v^{-2(\kappa+1)}$, and $f_\kappa\to f_M$ as $\kappa\to\infty$ (Verscharen et al. 2019). Smaller $\kappa$ correlates with weaker collisions / larger heliocentric distance (Štverák et al. 2022).

For **electrons**, the standard decomposition is a **Maxwellian core + kappa halo** (dual Maxwellian–kappa; Pierrard & Lazar 2010; Lazar et al. 2017):

$$
f_e(\mathbf v)=f_M^{\rm core}(\mathbf v)+f_\kappa^{\rm halo}(\mathbf v),
$$

with fitted halo indices $\kappa\sim 2$–6 decreasing with distance (e.g. $\kappa\approx 7.5$ at 0.4 AU, 5.4 at 1 AU, 3.2 at 3 AU; Štverák et al. 2022). A field-aligned **strahl** is a third component.

The standard kappa has divergent moments for low $\kappa$ (and unphysical superluminal particles). The **regularized kappa distribution (RKD)** (Scherer, Fichtner & Lazar 2017) multiplies by a Gaussian cutoff:

$$
f_{\rm RKD}(v)\;\propto\;\Big[1+\frac{v^2}{\kappa\theta^2}\Big]^{-(\kappa+1)}\exp\!\Big(-\alpha^{2}\frac{v^{2}}{\theta^{2}}\Big),\qquad 0<\alpha\ll 1,
$$

keeping all moments finite and extending the allowed range to $\kappa>0$.

### 4.5 Fitting and the "human bias" of model choice — motivation for a model-free method

Parametric fits minimize a weighted chi-square,

$$
\chi^2(\boldsymbol\theta)=\sum_i\frac{[f_{\rm obs}(\mathbf v_i)-f(\mathbf v_i;\boldsymbol\theta)]^2}{\sigma_i^2},
$$

most commonly by **Levenberg–Marquardt** (interpolating Gauss–Newton and gradient descent via damping $\lambda$),

$$
\big(\mathbf J^{\mathsf T}\mathbf J+\lambda\,\mathrm{diag}(\mathbf J^{\mathsf T}\mathbf J)\big)\,\delta\boldsymbol\theta=\mathbf J^{\mathsf T}\mathbf r,
$$

(MPFIT implementation; Markwardt 2009). An unsupervised alternative is a **Gaussian Mixture Model + Expectation–Maximization** (De Marco et al. 2023, on Solar Orbiter/PAS), modeling each velocity point as $P(\mathbf x_i\mid\boldsymbol\theta)=\sum_k w_k\,\mathcal N(\mathbf x_i\mid\boldsymbol\mu_k,\boldsymbol\Sigma_k)$ with $\sum_k w_k=1$, giving each point a probabilistic membership and removing the manual seeding/region-masking bias of LM (though still bi-Maxwellian-based).

Every parametric fit embeds priors that bias the result *before any data are touched*: the number of components (1 vs 2 vs 3 Maxwellians, Maxwellian vs kappa) predetermines whether a "beam" or "tail" exists; bi-Maxwellians/gyrotropic kappas assume gyrotropy and reflection symmetry, discarding skewness/agyrotropy/heat-flux distortions that are physically present (the normal-inverse-Gaussian work, A&A 2024, shows real proton VDFs carry skewness/kurtosis no symmetric model captures). The recovered $n,\mathbf u,T_{\parallel,\perp},\kappa$ thus depend on the modeler's prior. This motivates a **model-free basis expansion** — the Hermite/Laguerre decomposition of §5 — which imposes no population count, no equilibrium assumption, and no symmetry, lets the data populate the coefficients, keeps *all* the information, yields a quantitative "distance from Maxwellian," and connects directly to kinetic theory (phase mixing, Landau damping, velocity-space cascade) (Schekochihin; Servidio et al. 2017; Vencels et al. 2018).

---

## 5. Polynomial decomposition: orthonormal Hermite functions, the Fourier analogy, $c_{mn}$

### 5.1 The Fourier analogy

Just as a Fourier transform writes a signal as a sum of plane waves, the **Hermite transform writes the VDF as a Maxwellian times a series of Hermite polynomials.** Larosa et al. (2025) state it explicitly: *"a single Fourier mode describes a plane wave, and likewise a single Hermite mode describes a Maxwellian; the more Hermite modes are non-zero the more the VDF is far from a Maxwellian."* The expansion is exact and invertible (a spectral basis), and **Parseval–Plancherel holds**, so the velocity-space "energy" (enstrophy) is preserved between physical and spectral space.

### 5.2 The orthonormal Hermite functions

In field-aligned coordinates $(\xi_\parallel,\xi_\perp)$ the gyrotropic VDF is expanded as

$$
f(\xi_\parallel,\xi_\perp)=\sum_{m,n}c_{mn}\,\psi_m(\xi_\parallel)\,\psi_n(\xi_\perp),
$$

with the **orthonormal Hermite functions** (Larosa Eq. 2; Pezzi Eq. 4):

$$
\boxed{\;\psi_m(\xi)=\frac{H_m(\xi)}{\sqrt{2^m\,m!\,\sqrt\pi}}\;e^{-\xi^2/2}\;}
$$

where $H_m$ are the *physicists'* Hermite polynomials, $H_m(\xi)=(-1)^m e^{\xi^2}\frac{d^m}{d\xi^m}e^{-\xi^2}$ (Larosa Eq. 3). These satisfy the "Fourier-like" orthonormality

$$
\int_{-\infty}^{\infty}\psi_m(\xi)\,\psi_n(\xi)\,d\xi=\delta_{mn}
$$

— note the weight is **1**, not $e^{-\xi^2}$: the $e^{-\xi^2/2}$ envelope makes the *functions* (not just the polynomials) orthonormal; these are the symmetrically weighted Hermite functions / quantum-harmonic-oscillator eigenstates.

Key normalization choices (Pezzi; Larosa §2):
- The argument is the **normalized peculiar velocity** $\xi=(v-u)/v_{\rm th}$ — the basis is *shifted to the local bulk frame* and *scaled by the thermal speed* in each direction.
- Hence $m=0$ is the *local* drifting Maxwellian; the projection "quantifies high-order corrections to the particle DF, since the basis is shifted in the local fluid velocity frame, normalized to the ambient density and temperature" (Pezzi). Consequently large-scale fluctuations of density, bulk flow and temperature (1st/2nd moments) **do not contaminate the high-$m$ spectrum** (Larosa §2).

### 5.3 The coefficients $c_{mn}$ and their physical meaning

The spectral coefficient is the projection integral (Larosa Eq. 4; Coburn; Pezzi):

$$
\boxed{\;c_{mn}=\iint f(\xi_\parallel,\xi_\perp)\,\psi_m(\xi_\parallel)\,\psi_n(\xi_\perp)\,d\xi_\parallel\,d\xi_\perp\;}
$$

(for HL, replace $\psi_n(\xi_\perp)\to\Gamma_l^1(\hat v_\perp)$ and integrate $d^3v$). In practice the integral is done by **Gauss–Hermite (and Gauss–Laguerre) quadrature** after interpolating the measured VDF onto the polynomial roots (§6, §7), which satisfies Parseval–Plancherel to machine precision, "avoiding spurious aliasing and convergence problems" (Pezzi).

**Physical hierarchy of moments** (the heart of the method; Pezzi; Roytershteyn & Delzanno 2018):

| Index | Physical content |
|---|---|
| $c_{00}$ | Maxwellian amplitude $\sim$ **density** (local drifting Maxwellian); normalization reference |
| $m=1$ | **bulk-flow / drift** fluctuations |
| $m=2$ | **temperature** deformations (pressure anisotropy) |
| $m=3$ | **heat-flux** perturbations |
| low $m,n$ | fluid moments → "fluid–kinetic coupling is an intrinsic feature" (R&D 2018) |
| high $m,n$ | **fine-scale, non-thermal velocity-space structure** (beams, hammerheads, phase-space filamentation) |

Pezzi: "the index $m$ roughly corresponds to an order of the velocity moments: $m=1$ corresponds to bulk-flow fluctuations; $m=2$ to temperature deformations; $m=3$ to heat-flux perturbations, and so on." Because the basis is centered on the local Maxwellian, **all the physics lives in $m,n\ge 1$**.

**Enstrophy / free energy.** The summed squared coefficients give the velocity-space enstrophy (Pezzi Eq. 5; Larosa Eq. 7):

$$
\Omega=\int \delta f^2\,d^3v=\sum_{m>0}\big[c_m\big]^2,
$$

the cascading quantity, tied to a non-Maxwellianity indicator via $\Omega=\epsilon^2 n^2$, and "essentially the free energy in gyrokinetics" (Pezzi, citing Schekochihin). Larosa pairs this with a Kaufmann–Paterson entropy measure $M_{KP}=(s_M-s)/(\tfrac32 k_B)$.

### 5.4 Two encodings of the perpendicular plane: Hermite–Hermite (HH) vs Hermite–Laguerre (HL)

Gyrotropy means $f=f(v_\parallel,v_\perp)$ with $v_\perp=\sqrt{v_{\perp1}^2+v_{\perp2}^2}$. There are two natural spectral encodings — the central methodological choice.

**(a) Hermite–Hermite (HH), "Cartesian-like" gyrotropy — Larosa et al. 2025 (PSP).** Treat $v_\perp$ as a Cartesian axis and use $\psi_n(\xi_\perp)$ perpendicular too. Gyrotropy is imposed by **mirroring**: *"we impose $f(-v_\perp)=f(v_\perp)$ and extend the grid to negative values in the perpendicular direction … such a procedure implies null odd Hermite coefficients in the perpendicular direction"* (Larosa §2). So $c_{m,n}=0$ for odd $n$ — only even perpendicular modes carry power, which is **why the 1-D perpendicular spectrum shows an oscillatory (even/odd) pattern** (Larosa Fig. 3): "the oscillatory behavior … is due to the lack of power in the odd modes of the perpendicular spectrum caused by the gyrotropy assumption." This is the basis used in the present project (HH, $60\times60$ grid via `numpy.polynomial.hermite.hermgauss`).

**(b) Hermite–Laguerre (HL), cylindrical gyrotropy — Coburn et al. 2024 (Solar Orbiter electrons).** The *natural* basis for an azimuthally symmetric function uses **associated Laguerre polynomials in $\mu=v_\perp^2$**:

$$
F_e(v_\parallel,v_\perp)=\sum_{m}\sum_{l}c_{ml}\,\psi_m(\hat v_\parallel)\,\Gamma_l^{1}(\hat v_\perp),
\qquad
\Gamma_l^{n}(\hat v_\perp)=e^{-\hat v_\perp^2/2}\,L_l^{n}(\hat v_\perp^2),
$$

with $\hat v_\parallel=v_\parallel/v_{{\rm th},e}^\parallel$, $\hat v_\perp=v_\perp/v_{{\rm th},e}^\perp$, and weight index $n=1$. The $n=1$ weight absorbs the cylindrical Jacobian $v_\perp\,dv_\perp$, so the $\{\Gamma_l^1\}$ are orthonormal under the gyrotropic measure. This is the velocity-space basis of gyrokinetics (the Laguerre–Hermite pseudo-spectral formulation; Mandell et al. 2018).

**Why it matters.** HL is the *physically faithful* representation of a gyrotropic distribution: a ring/shell or a perpendicular Maxwellian is a single low-$l$ Laguerre mode ($\Gamma_0^1$ exactly), with no spurious odd modes. HH is simpler/uniform (same basis both axes) but pays the price of vanishing odd-$n$ modes and an oscillatory spectrum, and a perpendicular Maxwellian requires many even Hermite modes.

---

## 6. Gaussian quadrature numerics

### 6.1 Setup and the two relevant families

Approximate weighted integrals by an $N$-point rule

$$
\int_a^b p(x)\,f(x)\,dx\;\approx\;\sum_{i=1}^N w_i\,f(x_i),
$$

with fixed weight $p(x)\ge 0$, nodes $\{x_i\}$, weights $\{w_i\}$.

| Family | Weight $p(x)$ | Interval | Orthogonal polynomial |
|---|---|---|---|
| Gauss–**Hermite** | $e^{-x^2}$ | $(-\infty,\infty)$ | $H_N$ |
| Gauss–**Laguerre** | $e^{-x}$ | $(0,\infty)$ | $L_N$ |

An $N$-point rule has $2N$ free parameters, so one can hope to be exact for all polynomials of degree $\le 2N-1$. Gaussian quadrature achieves this maximal order **precisely when the nodes are the $N$ roots of $P_N$.**

### 6.2 Proof sketch: the nodes must be the roots of $P_N$

Let $\{P_k\}$ be orthogonal with respect to $p$: $\langle P_j,P_k\rangle=\int_a^b p\,P_jP_k\,dx=h_k\delta_{jk}$. Take the nodes $x_1,\dots,x_N$ as the $N$ (real, simple, interior) roots of $P_N$. For any polynomial $f$ of degree $\le 2N-1$, divide by $P_N$:

$$
f(x)=q(x)P_N(x)+r(x),\qquad \deg q\le N-1,\ \deg r\le N-1.
$$

Integrate against $p$. The first term vanishes by **orthogonality** ($q$ is a combination of $P_0,\dots,P_{N-1}$, each orthogonal to $P_N$):

$$
\int_a^b p\,f\,dx=\underbrace{\int_a^b p\,q\,P_N\,dx}_{=0}+\int_a^b p\,r\,dx=\int_a^b p\,r\,dx.
$$

At each node, $P_N(x_i)=0\Rightarrow f(x_i)=r(x_i)$, so $\sum_i w_i f(x_i)=\sum_i w_i r(x_i)$. If the weights integrate every degree-$\le N-1$ polynomial exactly (§6.3), then $\sum_i w_i r(x_i)=\int_a^b p\,r\,dx$, hence

$$
\sum_{i=1}^N w_i f(x_i)=\int_a^b p\,f\,dx.
$$

The choice $x_i=\text{roots of }P_N$ is exactly what makes the high-degree remainder $q\,P_N$ integrate to zero; any other node set loses the extra $N$ orders. The roots are guaranteed real, simple, and interior because $P_N$ is orthogonal with respect to a positive weight.

### 6.3 Node weights: Lagrange form and Christoffel–Darboux closed form

Demand exactness for degree $\le N-1$ via the **Lagrange basis** $\ell_i(x)=\prod_{j\ne i}\frac{x-x_j}{x_i-x_j}$ ($\ell_i(x_j)=\delta_{ij}$):

$$
\boxed{\,w_i=\int_a^b p(x)\,\ell_i(x)\,dx\,}.
$$

Using the roots of $P_N$, $\ell_i(x)=\dfrac{P_N(x)}{P_N'(x_i)\,(x-x_i)}$, so

$$
\boxed{\,w_i=\frac{1}{P_N'(x_i)}\int_a^b\frac{p(x)\,P_N(x)}{x-x_i}\,dx\,}.
$$

The **Christoffel–Darboux identity** for orthonormal $\{\hat P_k\}$,

$$
\sum_{k=0}^{N-1}\hat P_k(x)\hat P_k(y)=\frac{a_N[\hat P_N(x)\hat P_{N-1}(y)-\hat P_{N-1}(x)\hat P_N(y)]}{x-y},
$$

collapses the integral into the **Christoffel-number** form

$$
\boxed{\,w_i=\frac{1}{\sum_{k=0}^{N-1}\hat P_k(x_i)^2}\,}.
$$

Specializing (Abramowitz & Stegun, Ch. 25):

$$
\textbf{Hermite (A\&S 25.4.46):}\quad
w_i=\frac{2^{\,N-1}\,N!\,\sqrt\pi}{N^2\,[H_{N-1}(x_i)]^2},
$$
$$
\textbf{Laguerre (A\&S 25.4.45):}\quad
w_i=\frac{x_i}{(N+1)^2\,[L_{N+1}(x_i)]^2}.
$$

NumPy uses the equivalent $w_k=c/\big(H'_N(x_k)\,H_{N-1}(x_k)\big)$ with $c$ fixed by $\int p\,dx$.

### 6.4 Computing the nodes: Newton–Raphson vs Golub–Welsch

**Newton–Raphson** solves $P_N(x_i)=0$ directly, iterating $x^{(k+1)}=x^{(k)}-P_N(x^{(k)})/P_N'(x^{(k)})$ via the three-term recurrence. Fragile: roots cluster, $P_N$ has huge dynamic range, so without good asymptotic seeds Newton can converge to an already-found root, skip a root, or overflow; it needs careful deflation/seeding.

**Golub–Welsch (1969)** — the robust modern method. The orthonormal polynomials obey

$$
x\,P_m(x)=a_{m+1}P_{m+1}(x)+b_m P_m(x)+a_m P_{m-1}(x),
$$

which stacks into $x\,\mathbf P(x)=J\,\mathbf P(x)+a_N P_N(x)\,\mathbf e_N$ with the **symmetric tridiagonal Jacobi matrix**

$$
J=\begin{pmatrix}
b_0 & a_1 & & \\
a_1 & b_1 & a_2 & \\
& a_2 & b_2 & \ddots \\
& & \ddots & \ddots & a_{N-1}\\
& & & a_{N-1} & b_{N-1}
\end{pmatrix}.
$$

At a node $x_i$ the last term drops, so $J\mathbf P(x_i)=x_i\mathbf P(x_i)$:

- **Nodes** $x_i$ = **eigenvalues** of $J$ (real, simple, computed stably by QR);
- **Weights** $w_i=\mu_0\,(v_{0,i})^2$, where $v_{0,i}$ is the first component of the $i$-th normalized eigenvector and $\mu_0=\int_a^b p\,dx$ is the zeroth moment.

**Recurrence coefficients** (monic $\beta_m=a_m^2$):

$$
\textbf{Hermite:}\quad b_m=0,\quad a_m=\sqrt{\tfrac m2},\quad \mu_0=\int_{-\infty}^\infty e^{-x^2}dx=\sqrt\pi;
$$
$$
\textbf{Laguerre:}\quad b_m=2m+1,\quad a_m=m,\quad \mu_0=\int_0^\infty e^{-x}dx=1.
$$

(Equivalently monic $\beta_m=m/2$ Hermite, $\beta_m=m^2$ Laguerre.)

### 6.5 Practical tools (NumPy)

```python
import numpy as np
x_H, w_H = np.polynomial.hermite.hermgauss(N)    # weight e^{-x^2} on (-inf, inf)
x_L, w_L = np.polynomial.laguerre.laggauss(N)    # weight e^{-x}   on (0, inf)
```

`hermgauss(deg)` returns nodes/weights that "correctly integrate polynomials of degree `2*deg - 1` or less" with weight $e^{-x^2}$; internally it builds the companion (Jacobi) matrix, takes its eigenvalues for the nodes (refined by one Newton step), and computes weights from the §6.3 closed form — a modified Golub–Welsch for $N\lesssim 150$. **Convention note:** `hermgauss` is the *physicists'* $e^{-x^2}$; `hermegauss` is the *probabilists'* $e^{-x^2/2}$ — match the convention to your basis. For the $60\times60$ grid of this project, $N=60$ in each direction.

---

## 7. Interpolation onto the node grid (inverse-distance weighting)

The measured VDF lives at **scattered points** $(\xi_\parallel^{(k)},\xi_\perp^{(k)})$ with values $f_k$, but the quadrature requires $f$ **at the fixed node grid** $(\xi_{\parallel,i},\xi_{\perp,j})$ — the Hermite roots, which almost never coincide with measurement locations. Use **inverse-distance weighting (IDW / Shepard interpolation)**:

$$
f(\xi_{\parallel,i},\xi_{\perp,j})=\frac{\sum_k W_k\,f_k}{\sum_k W_k},
\qquad
W_k=\frac{1}{r_k^2+\varepsilon},
\qquad
r_k^2=(\xi_{\parallel,i}-\xi_\parallel^{(k)})^2+(\xi_{\perp,j}-\xi_\perp^{(k)})^2 .
$$

- The **convex-combination** form ($\sum_k W_k f_k/\sum_k W_k$) guarantees the interpolant lies within the data range and reproduces constants exactly.
- The regularizer $\varepsilon>0$ (a softening length) prevents the $1/r^2$ singularity when a node nearly coincides with a measurement and controls smoothness: small $\varepsilon$ → sharp, nearly nearest-neighbor; larger $\varepsilon$ → smoother but more diffusive. Choose $\varepsilon$ comparable to the squared local sample spacing.
- This **interpolation error is the dominant error floor of the whole pipeline** (see §6.5 and §8.4): it is independent of how many quadrature nodes you add, and it is what makes the spectrum plateau beyond $N\approx 60$.

---

## 8. The 2-D velocity-space spectrum, anisotropic cascade, and noise-floor flattening

### 8.1 The Hermite-coefficient quadrature formula

The projection $c_{mn}=\iint f\,\psi_m(\xi_\parallel)\psi_n(\xi_\perp)\,d\xi_\parallel d\xi_\perp$ is an **unweighted** integral (weight 1), but Gauss–Hermite carries the weight $e^{-\xi^2}$. Insert $1=e^{-\xi^2}e^{+\xi^2}$ so the bracketed quantity is what the rule integrates:

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

The factors $e^{\xi_i^2},e^{\xi_j^2}$ exactly undo the $e^{-\xi^2}$ weight baked into the Gauss–Hermite weights, converting the weighted quadrature into the unweighted projection. **Numerical stability:** $e^{\xi^2}$ and the $e^{-\xi^2/2}$ inside $\psi$ partly cancel; for large $|\xi_i|$ combine them as $e^{\xi^2}\psi_m(\xi)=H_m(\xi)\,e^{-\xi^2/2}\cdot e^{\xi^2}/\sqrt{2^m m!\sqrt\pi}$ and work with the $H_m$ form to avoid overflow (the "stable Hermite transform"; arXiv:2604.02041).

**$(M,N)$ convergence.** Two integers control the calculation: $M$ = number of Hermite modes retained ($m,n\le M$), $N$ = quadrature nodes per dimension. If $f$ is itself a Hermite series of order $\le M$, the integrand $f\,e^{\xi^2}\psi_m\psi_n$ is (polynomial of degree $\le 2M$)$\times e^{-\xi^2}$; Gauss–Hermite is exact to degree $2N-1$, so **$N\ge M+1$ integrates every retained mode exactly** (below that, low-order coefficients are aliased by high-order content). A *measured/interpolated* VDF is not a finite Hermite polynomial, so one pushes $N$ higher to resolve tails — the coefficients change rapidly, then **plateau beyond $N\approx 60$**, because the residual error is then dominated by the IDW interpolation error (§7), which $N$ cannot reduce. Hence $N\simeq 60$ sits just past the convergence plateau at the interpolation-limited floor while comfortably satisfying $N\ge M+1$ for the modest $M$ of interest — the best accuracy-per-cost (smaller $N$ under-integrates; larger $N$ wastes the $N^2$ tensor grid and $O(N^3)$ eigensolve).

### 8.2 The 2-D spectrum and its reductions

Plotting $\log_{10}(c_{mn}^2/c_{00}^2)$ over the integer lattice $(m,n)$ is the velocity-space analogue of a 2-D power spectrum (Pezzi Fig. 2; Larosa Fig. 3). Normalizing by $c_{00}^2$ (the Maxwellian power) makes it dimensionless and measures *fractional* departure from Maxwellian; the "DC term" $c_{00}$ is excluded from the cascade sum (it is the equilibrium, not a fluctuation).

Two complementary reductions:

- **Directional cuts.** The parallel cut $c_{m0}^2/c_{00}^2$ (perpendicular index fixed at 0) measures structure *along* $\mathbf B$; the perpendicular cut $c_{0n}^2/c_{00}^2$ measures structure *across* $\mathbf B$. Pezzi forms reduced 1-D spectra $P(m_z)$ and $P(m_\perp)$ this way.
- **Isotropic shell-averaged spectrum** (Servidio 2017; Pezzi; Larosa): sum over concentric shells of unit width,
$$
P(m)=\!\!\sum_{m-\frac12<|\mathbf m'|\le m+\frac12}\!\!c_{\mathbf m'}^2,\qquad
\mathbf m=(m_\perp,m_\parallel),\ \ m=\sqrt{m_\perp^2+m_\parallel^2},
$$
the direct analogue of the omnidirectional $k$-spectrum in fluid turbulence.

### 8.3 The anisotropic velocity-space cascade and spectral slopes

The central physical result (Pezzi 2018, confirming Servidio 2017): **a background magnetic field makes the velocity-space cascade anisotropic — it extends to higher orders (shallower slope) along $\mathbf B$ and is suppressed (steeper slope) across $\mathbf B$.** Pezzi: "spectra are stretched in the parallel direction … velocity gradients are stronger along the mean field and the cascade is inhibited across the mean field" — an analogue of the Shebalin effect. Mechanistically, parallel extension reflects **parallel phase mixing / Landau resonances** generating fine $v_\parallel$ structure.

**Quantitatively (Pezzi 2018, Fig. 3, verified numbers):**

$$
P(m_z)\sim m_z^{-2.01}\ (\text{parallel}),\qquad
P(m_\perp)\sim m_\perp^{-3.5}\ (\text{perpendicular, much lower and steeper}).
$$

The shell-averaged exponent is *diagnostic of the dominant process* (Servidio 2017; Larosa §1):

- $P(m)\sim m^{-3/2}$ — **weakly magnetized / isotropic** phase-space cascade (phase-mixing / electric-field-dominated regime);
- $P(m)\sim m^{-2}$ — **magnetized / magnetic-field-dominated** regime; steeper indicates nonlinear-dominated dynamics.

Larosa et al. (2025, PSP) measure $P(m)\approx m^{-2}$ over $4\le m<12$ for both a *wave* stream ($\sim 28\,R_\odot$, superalfvénic, hammerhead beams) and a *turbulent* stream ($\sim 11\,R_\odot$, subalfvénic), "consistent with expectations for a low-$\beta$ plasma (Servidio et al. 2017)," while cautioning the slope is "very sensitive to the chosen range of $m$." Their magnetic-field PSD slopes mirror the v-space story (wave stream $\sim f^{-2}$, turbulent $\sim f^{-3/2}$) — a real-space/v-space dual-cascade correspondence.

**Low-order vs high-order — beams vs cascade** (Larosa's key contrast): at low $m$ ($1,2,3$) the wave stream has *more* power — the dense proton **beam + hammerhead** projects onto just the first few Hermite modes ("a dense beam at a few thermal speeds … provide power only to the first few Hermite modes"); at high $m$ ($8\lesssim m\lesssim 15$) the turbulent stream has *more* power, "interpreted as an effect of the stronger velocity-space cascade … fine-scale VDF distortions."

> **Project-specific working conventions (not in the cited literature).** This project's reported parallel/perpendicular slope pair "$m^{-1.3}$ (parallel), $n^{-3.1}$ (perpendicular)" and the local-slope estimator
> $$\alpha_m=\frac{c_{m+2,0}^2-c_{m,0}^2}{2\,c_{00}^2}$$
> are *working conventions*, not quoted from Pezzi/Servidio/Larosa/Coburn. They encode the same robust, citable inequality $|{\rm slope}_\parallel|<|{\rm slope}_\perp|$ (extended parallel cascade, suppressed perpendicular). The estimator is a two-mode-spaced, $c_{00}$-normalized finite difference of the parallel-cut power — a discrete proxy for $dP/dm$ along $\mathbf B$, flagging where the cascade steepens and transitions to the noise plateau. The **"+2" spacing is deliberate**: in HH gyrotropy odd perpendicular modes vanish and even/odd Hermite modes carry systematically different power (the oscillatory spectrum), so differencing $m\to m{+}2$ skips the parity oscillation. Cross-check it against the shell-averaged $P(m)$.

### 8.4 Spectral flattening (noise floor), truncation, and reconstruction

At high order the measured spectrum **flattens to a plateau** — this is **instrumental noise / interpolation error, not physics**:

- **Larosa 2025:** "we purposefully did not discuss the spectra for $m\gtrsim 20$ because the flattening is likely due to interpolation errors"; they validate that the signal lies above one-count and bi-Maxwellian noise estimates (their Appendix B).
- **Coburn 2024 (the explicit denoising application):** the HL coefficients "decrease in relative strength until $m,l\approx 14$ where the spectrum flattens," identified as the noise floor. They then **truncate** — zero all coefficients beyond the floor — and **reconstruct** a "low-pass-filtered VDF": smooth, analytic, differentiable, free of instrumental shot noise. Because the basis is analytic, the **denoised VDF and its velocity derivatives** can be computed exactly, enabling wave–particle-interaction analyses. The grid size ($60\times60$ for Coburn; $50\times50$ for Larosa) is chosen so the *low-order* spectrum is converged and unaffected by the truncation.

This is the direct analogue of low-pass filtering a Fourier spectrum: keep the resolved large-scale modes, discard the high-"wavenumber" noise. For the present project's $60\times60$ HH analysis, the same logic applies — read the converged low-order modes as physical (drift, anisotropy, heat flux, beam/hammerhead at low $m$; cascade structure at intermediate $m$), and treat the high-order plateau as the interpolation-limited noise floor.

### 8.5 One-paragraph synthesis

A solar-wind proton VDF measured by SWA-PAS is reconstructed from $(E,\text{elev},\text{azim})$ to $\mathbf v_{xyz}$ (with the mandatory SWA minus-sign convention so that $V_R>0$), rotated to RTN via `PAS_to_RTN`, then to field-aligned coordinates, and normalized to $\xi=(v-u)/w$ using the moment-derived $n,\mathbf U,T_{\parallel,\perp}$ and thermal speeds $w=\sqrt{2k_BT/m}$. Under the gyrotropy assumption it is decomposed on the orthonormal Hermite basis $\psi_m=H_m(\xi)e^{-\xi^2/2}/\sqrt{2^m m!\sqrt\pi}$ (Hermite–Hermite, with $f(-v_\perp)=f(v_\perp)$ killing odd perpendicular modes), the velocity-space analogue of a Fourier transform centered on the local Maxwellian. The coefficients $c_{mn}=\iint f\psi_m\psi_n\,d\xi$ are evaluated by Gauss–Hermite quadrature whose nodes are the Hermite roots (Golub–Welsch / `hermgauss`, $N=60$) after IDW interpolation of the scattered data onto the node grid. The squared, $c_{00}^2$-normalized spectrum $\log_{10}(c_{mn}^2/c_{00}^2)$ exposes the hierarchy $c_{00}\sim n$, $m=1,2,3\sim$ drift/temperature/heat-flux, high $m,n\sim$ fine non-thermal structure; its shell-averaged form is a power law ($m^{-3/2}$ weak-field to $m^{-2}$ magnetized), its directional cuts reveal an anisotropic cascade (shallow/extended parallel, steep/suppressed perpendicular; Pezzi $m_z^{-2.01}$ vs $m_\perp^{-3.5}$), and its high-order flattening marks an instrumental/interpolation noise floor that motivates spectral truncation and reconstruction of a denoised, analytic VDF (Coburn).

---

## 9. References

**Instrument and data product**
- Owen, C. J., et al. 2020, "The Solar Orbiter Solar Wind Analyser (SWA) suite," *A&A* **642**, A16. DOI: 10.1051/0004-6361/201937259. https://www.aanda.org/articles/aa/full_html/2020/10/aa37259-19/aa37259-19.html
- Fedorov, A. 2020, *SWA-PAS L2 Data User Guide*, V02 (16 Nov 2020), ESA COSMOS. https://www.cosmos.esa.int/documents/3689933/11863901/PAS_L2_Data_User_Guide_20201116.pdf
- SWA-PAS L2 Data User Guide (online). https://swa-pas-data-user-guide.readthedocs.io/en/latest/products.html ; mirror: https://rungk-om.github.io/pas-user-guide/
- Louarn, P., et al. 2024, "…solar-wind proton moments from SWA-PAS," *A&A* **682**, A44. DOI: 10.1051/0004-6361/202347874. https://verscharen.com/papers/louarn_etal24.pdf
- Livi, S., et al. 2023, "First results from the Solar Orbiter Heavy Ion Sensor," *A&A* **676**, A36 (SWA-HIS, **not** PAS). DOI: 10.1051/0004-6361/202346304. https://www.aanda.org/articles/aa/full_html/2023/08/aa46304-23/aa46304-23.html
- IRAP SWA/PAS page: https://www.irap.omp.eu/en/project/solar-orbiter-pas-2/ ; CNES SWA page: https://solar-orbiter.cnes.fr/en/SOLO/GP_swa.htm

**Counts → phase-space density, reconstruction, moments, gyrotropy**
- Livi, R., et al. 2022, "The Solar Probe ANalyzer—Ions on PSP," *ApJ* **938**, 138. https://iopscience.iop.org/article/10.3847/1538-4357/ac93f5
- Verscharen, D., Klein, K. G. & Maruca, B. A. 2019, "The multi-scale nature of the solar wind," *Living Rev. Solar Phys.* **16**:5. https://pmc.ncbi.nlm.nih.gov/articles/PMC6934245/ ; arXiv:1902.03448
- Broderick/Klein et al. 2025, "Recovering Ion Distribution Functions: Slepian Reconstruction (MMS & Solar Orbiter), I," arXiv:2501.17294, https://arxiv.org/html/2501.17294 ; II (Gyrotropic), *ApJ*, DOI: 10.3847/1538-4357/ae1d71, https://iopscience.iop.org/article/10.3847/1538-4357/ae1d71 (and DOI: 10.3847/1538-4357/adb6a0, https://iopscience.iop.org/article/10.3847/1538-4357/adb6a0)
- Wilson, L. B., et al. 2025, "How limited resolution of plasma analyzers affects moment accuracy," arXiv:2505.09869. https://arxiv.org/pdf/2505.09869
- Verniero, J. L., et al. 2020, "PSP Observations of Proton Beams Simultaneous with Ion-scale Waves," *ApJS* **248**, 5. https://iopscience.iop.org/article/10.3847/1538-4365/ab86af
- Swisdak, M. 2016, "Quantifying gyrotropy in magnetic reconnection," *GRL*. DOI: 10.1002/2015GL066980. https://agupubs.onlinelibrary.wiley.com/doi/full/10.1002/2015GL066980
- Wikipedia, "Thermal velocity." https://en.wikipedia.org/wiki/Thermal_velocity

**Parametric models and fitting**
- Marsch, E. 2006, "Kinetic Physics of the Solar Corona and Solar Wind," *Living Rev. Solar Phys.* https://link.springer.com/article/10.12942/lrsp-2006-1
- De Marco, R., et al. 2023, "Innovative technique for separating proton core, proton beam, and alpha particles" (GMM + EM), *A&A*. https://www.aanda.org/articles/aa/full_html/2023/01/aa43719-22/aa43719-22.html
- Dupuis, R., et al. 2020, "Characterizing magnetic reconnection regions using Gaussian mixture models on particle VDFs," arXiv:1910.10012. https://arxiv.org/pdf/1910.10012
- Štverák, Š., et al. 2022, "Implications of Kappa Suprathermal Halo of the Solar Wind Electrons," *Front. Astron. Space Sci.* https://www.frontiersin.org/journals/astronomy-and-space-sciences/articles/10.3389/fspas.2022.892236/full
- Lazar, M., et al. 2017, "Dual Maxwellian–Kappa modeling of the solar wind electrons," *A&A*. https://www.aanda.org/articles/aa/full_html/2017/06/aa30194-16/aa30194-16.html
- Pierrard, V. & Lazar, M. 2010, "Kappa Distributions: Theory and Applications in Space Plasmas," *Solar Phys.* https://www.researchgate.net/publication/45906919_Kappa_Distributions_Theory_and_Applications_in_Space_Plasmas
- Scherer, K., Fichtner, H. & Lazar, M. — regularized kappa distributions (RKD). https://www.mdpi.com/2571-6182/6/3/36 ; "Regularized kappa-halos with ALPS," arXiv:2504.15955, https://arxiv.org/html/2504.15955
- "Skewness and kurtosis of solar wind proton distribution functions" (normal inverse-Gaussian) 2024, *A&A*. https://www.aanda.org/articles/aa/full_html/2024/02/aa47874-23/aa47874-23.html
- Tu, C.-Y., et al. 2004, "Dependence of the proton beam drift velocity on proton core plasma beta," *JGR*. https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2004JA010391

**Hermite / Hermite–Laguerre spectral decomposition and velocity-space cascade**
- Servidio, S., et al. 2017, "MMS observation of plasma velocity-space cascade," *PRL* **119**, 205101. https://link.aps.org/doi/10.1103/PhysRevLett.119.205101 ; arXiv:1707.08180
- Pezzi, O., et al. 2018, "Velocity-space cascade in magnetized plasmas," *Phys. Plasmas* **25**, 060704. https://pubs.aip.org/aip/pop/article/25/6/060704/320076 ; arXiv:1803.01633
- Coburn, J., et al. 2024, *ApJ* **964**, 100 (Solar Orbiter electrons, Hermite–Laguerre, $60\times60$). https://iopscience.iop.org/article/10.3847/1538-4357/ad1329 ; ADS: https://ui.adsabs.harvard.edu/abs/2024ApJ...964..100C
- Larosa, A., et al. 2025, "Velocity-space turbulent cascade in the near-Sun solar wind" (**Parker Solar Probe**, SPAN-i, Hermite–Hermite, $50\times50$). arXiv:2512.01492. https://arxiv.org/abs/2512.01492 ; https://arxiv.org/html/2512.01492
- Roytershteyn, V. & Delzanno, G. L. 2018, *Front. Astron. Space Sci.* **5**, 27. https://www.frontiersin.org/journals/astronomy-and-space-sciences/articles/10.3389/fspas.2018.00027/full
- Mandell, N. R., et al. 2018, "Laguerre–Hermite pseudo-spectral velocity formulation of gyrokinetics," arXiv:1708.04029. https://arxiv.org/pdf/1708.04029
- Vencels, J., et al. 2018, "Spectral Approach to Plasma Kinetic Simulations Based on Hermite Decomposition," *Front. Astron. Space Sci.* https://www.frontiersin.org/journals/astronomy-and-space-sciences/articles/10.3389/fspas.2018.00027/full

**Quadrature numerics**
- Golub, G. H. & Welsch, J. H. 1969, "Calculation of Gauss Quadrature Rules," *Math. Comp.* **23**(106), 221–230. https://bibbase.org/network/publication/golub-welsch-calculationofgaussquadraturerules-1969
- Reichel, L., "Computation of Gauss-type quadrature rules." https://www.math.kent.edu/~reichel/publications/DC.pdf
- Abramowitz, M. & Stegun, I. A., *Handbook of Mathematical Functions*, Ch. 25 (formulas 25.4.45, 25.4.46).
- Wolfram MathWorld, "Hermite–Gauss Quadrature." https://mathworld.wolfram.com/Hermite-GaussQuadrature.html
- NumPy: `numpy.polynomial.hermite.hermgauss`, https://numpy.org/doc/stable/reference/generated/numpy.polynomial.hermite.hermgauss.html ; `numpy.polynomial.laguerre.laggauss`
- Wikipedia, "Gauss–Legendre quadrature." https://en.wikipedia.org/wiki/Gauss%E2%80%93Legendre_quadrature
- ACME lab notes, "Gaussian Quadrature." https://acme.byu.edu/0000017a-17ef-d8b9-adfe-77ef21000001/vol2a-gaussianquadrature-fall2016-pdf
- "Stable Hermite transforms via the Golub–Welsch algorithm," arXiv:2604.02041. https://arxiv.org/pdf/2604.02041