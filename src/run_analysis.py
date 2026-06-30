"""
Pipeline principal: reproduce las tres figuras de la tarea.

  Tarea 1 -> figures/task1_vdf_scatter.png   (slide 3 sin panel de energía)
  Tarea 2 -> figures/task2_hermite_grid.png  (slide 23, paneles a y b, HH)
  Tarea 3 -> figures/task3_hh_spectrum.png   (slide 25, panel izquierdo HH)
            figures/task3_spectrum_cuts.png  (diagnóstico, slide 26)

Uso:  .venv/bin/python src/run_analysis.py
"""
import os, sys, json
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import config as C
from cdf_io import read_vdf_timestep, read_mag_at, vdf_info
import physics as P
import hermite as H
import plots


def match_axes(vdf_shape, n_E, n_el, n_az):
    """Devuelve los índices de eje (ax_E, ax_el, ax_az) en la VDF según los tamaños."""
    axes = {}
    used = set()
    for label, size in (("E", n_E), ("el", n_el), ("az", n_az)):
        for ax, s in enumerate(vdf_shape):
            if s == size and ax not in used:
                axes[label] = ax; used.add(ax); break
    if len(axes) != 3:
        raise RuntimeError(f"No pude mapear ejes: vdf{vdf_shape} con (E={n_E},el={n_el},az={n_az})")
    return axes["E"], axes["el"], axes["az"]


def main():
    print("=" * 70)
    print("Reproducción VDF — Solar Orbiter SWA-PAS —", C.TARGET_TIME)
    print("=" * 70)

    # ---- 1. estructura del CDF de la VDF -------------------------------------
    print("\n[1] Estructura del CDF de la VDF:")
    for v, dims, nrec, dtype, units in vdf_info(C.VDF_CDF):
        print(f"    {v:16s} dims={str(dims):14s} nrec={nrec:<7d} {dtype:12s} {units}")

    # ---- 2. lee el paso temporal ---------------------------------------------
    d = read_vdf_timestep(C.VDF_CDF, C.TARGET_TIME)
    vdf = d["vdf"]; energy = d["energy"]; azim = d["azimuth"]; elev = d["elevation"]; M = d["M"]
    print(f"\n[2] Paso temporal usado: {d['t_used']}  (idx {d['idx']})  Info={d.get('info')}")
    print(f"    vdf shape={vdf.shape}  units={d['vdf_units']}")
    print(f"    Energy: {energy.size} ch  [{energy.min():.1f}, {energy.max():.1f}] eV")
    print(f"    Azimuth: {azim.size}  [{azim.min():.1f}, {azim.max():.1f}] deg")
    print(f"    Elevation: {elev.size}  [{elev.min():.1f}, {elev.max():.1f}] deg")
    print(f"    PAS_to_RTN det={np.linalg.det(M):+.3f}")

    n_E, n_el, n_az = energy.size, elev.size, azim.size
    axE, axel, axaz = match_axes(vdf.shape, n_E, n_el, n_az)
    print(f"    Ejes -> E:{axE}  el:{axel}  az:{axaz}")

    # ---- 3. mallas 3D de (E, el, az) en el orden de la VDF -------------------
    # anchos de bin reales del CDF (mejor que np.gradient): full = 2*half
    dE   = (d["delta_p_E"] + d["delta_m_E"]) if d["delta_p_E"] is not None else np.abs(np.gradient(energy))
    dele = (2.0 * d["delta_El"])            if d["delta_El"]  is not None else np.abs(np.gradient(elev))
    daz  = (2.0 * d["delta_Az"])            if d["delta_Az"]  is not None else np.abs(np.gradient(azim))

    # mallas 3D en el orden de ejes nativo de la VDF (indexing='ij')
    coords = [None, None, None]; widths = [None, None, None]
    coords[axE], widths[axE]   = energy, dE
    coords[axel], widths[axel] = elev, dele
    coords[axaz], widths[axaz] = azim, daz
    g  = np.meshgrid(coords[0], coords[1], coords[2], indexing="ij")
    gw = np.meshgrid(widths[0], widths[1], widths[2], indexing="ij")
    Emesh  = g[axE];  elmesh  = g[axel];  azmesh  = g[axaz]
    dEmesh = gw[axE]; delmesh = gw[axel]; dazmesh = gw[axaz]

    # ---- 4. reconstrucción de velocidades ------------------------------------
    # Convención SWA (Fedorov 2020): Vx=-cos·cos, Vy=-cos·sin, Vz=+sin  -> V_R>0
    def reconstruct(signs):
        vx, vy, vz = P.instrument_velocity(Emesh, azmesh, elmesh, signs=signs)
        return P.rotate_pas_to_rtn(vx, vy, vz, M)

    dV = P.velocity_volume(Emesh, dEmesh, elmesh, delmesh, dazmesh)
    signs = (-1.0, -1.0, +1.0)
    vR, vT, vN = reconstruct(signs)
    w = vdf * dV
    n_tot = w.sum()
    uR_test = (vR * w).sum() / n_tot
    if uR_test < 0:                       # respaldo: negación global del vector
        signs = (+1.0, +1.0, -1.0)
        vR, vT, vN = reconstruct(signs)
        uR_test = (vR * w).sum() / n_tot
    sign = signs
    print(f"\n[3] Convención de reconstrucción: signs={signs}  ->  u_R={uR_test/1e3:+.1f} km/s")

    # ---- 5. marco alineado al campo (b a lo largo del flujo) ----------------
    Braw, _, _ = read_mag_at(C.MAG_CDF, C.TARGET_TIME, window_s=4.0)
    bhat0 = Braw / np.linalg.norm(Braw)
    # momentos preliminares para fijar el signo de b a lo largo del flujo
    mom0 = P.compute_moments(vR, vT, vN, vdf, dV, bhat0)
    if mom0["u_par"] < 0:
        bhat0 = -bhat0
    b, e1, e2 = P.field_aligned_basis(bhat0 * np.linalg.norm(Braw))
    # asegura b alineado a lo largo del flujo de nuevo tras la base
    if (np.array([mom0["u"]]) @ b)[0] < 0:
        b = -b; e1 = -e1   # mantiene base derecha (e2 = b x e1)

    mom = P.compute_moments(vR, vT, vN, vdf, dV, b)
    n_cc = mom["n"] / 1e6                      # m^-3 -> cm^-3
    print(f"\n[4] Momentos del protón:")
    print(f"    n      = {n_cc:.2f} cm^-3")
    print(f"    |u|    = {np.linalg.norm(mom['u'])/1e3:.1f} km/s   u_RTN={np.round(mom['u']/1e3,1)}")
    print(f"    u_par  = {mom['u_par']/1e3:.1f} km/s")
    print(f"    T_par  = {mom['T_par']:.3e} K   T_perp = {mom['T_perp']:.3e} K")
    print(f"    w_par  = {mom['w_par']/1e3:.1f} km/s  w_perp = {mom['w_perp']/1e3:.1f} km/s")
    print(f"    T_perp/T_par = {mom['T_perp']/mom['T_par']:.2f}")

    # ---- 6. proyección FA + coordenadas normalizadas -------------------------
    v_par, v_p1, v_p2 = P.project_field_aligned(vR, vT, vN, b, e1, e2)
    # perp con deriva perpendicular removida (centra el anillo girotrópico en 0,
    # como slide 3 / 23a); paralelo se deja absoluto para conservar núcleo/haz.
    up1 = float(mom["u"] @ e1); up2 = float(mom["u"] @ e2)
    v_p1_dd = v_p1 - up1
    v_p2_dd = v_p2 - up2
    xi_par, xi_perp = P.normalized_coords(vR, vT, vN, mom["u"], b, mom["w_par"], mom["w_perp"])

    # aplana puntos válidos (vdf > 0)
    valid = (vdf > 0) & np.isfinite(vdf)
    f_v = vdf[valid]
    vpar_kms = v_par[valid] / 1e3
    vp1_kms  = v_p1_dd[valid] / 1e3
    vp2_kms  = v_p2_dd[valid] / 1e3
    xipar_v  = xi_par[valid]
    xiperp_v = xi_perp[valid]
    print(f"\n[5] Puntos válidos (vdf>0): {valid.sum()} de {vdf.size}")
    print(f"    PSD rango: [{f_v.min():.2e}, {f_v.max():.2e}] {d['vdf_units']}")

    # ================= TAREA 1 =================
    out1 = os.path.join(C.FIG_DIR, "task1_vdf_scatter.png")
    plots.plot_scatter_vdf(vpar_kms, vp1_kms, vp2_kms, f_v,
                           str(d["t_used"])[:19].replace("T", " "), out1)
    print(f"\n[Tarea 1] {out1}")

    # ================= TAREA 2 (HH) =================
    nodes, weights = H.gauss_hermite(C.N_NODES)
    logf_obs = np.log10(f_v)
    # refleja la VDF perpendicular (estilo Larosa: f(xi_perp)=f(-xi_perp))
    src_xi_par  = np.concatenate([xipar_v, xipar_v])
    src_xi_perp = np.concatenate([xiperp_v, -xiperp_v])
    src_logf    = np.concatenate([logf_obs, logf_obs])
    logf_grid = H.idw_interpolate(src_xi_par, src_xi_perp, src_logf,
                                  nodes, nodes, eps=C.IDW_EPS)
    out2 = os.path.join(C.FIG_DIR, "task2_hermite_grid.png")
    plots.plot_hermite_grid(xipar_v, xiperp_v, logf_obs, nodes, logf_grid, out2)
    print(f"[Tarea 2] {out2}")

    # ================= TAREA 3 (espectro HH) =================
    f_grid = 10.0 ** logf_grid                       # vuelve a lineal para c_mn
    c = H.hh_coefficients(f_grid, nodes, weights, C.M_ORDER, C.N_ORDER)
    out3 = os.path.join(C.FIG_DIR, "task3_hh_spectrum.png")
    plots.plot_hh_spectrum(c, out3)
    out3b = os.path.join(C.FIG_DIR, "task3_spectrum_cuts.png")
    plots.plot_spectrum_cuts(c, out3b)
    print(f"[Tarea 3] {out3}")
    print(f"          {out3b}")
    print(f"          c_00={c[0,0]:.4e}  log(c10^2/c00^2)={np.log10((c[1,0]**2)/(c[0,0]**2)):.2f}")

    # ================= BONUS: Hermite-Laguerre (HL, estilo Coburn) =================
    l_mu, l_w = H.gauss_laguerre(C.N_NODES)
    xiperp_hl = np.sqrt(l_mu)                          # nodos perp del grid HL (mu=xi_perp^2)
    # (display) IDW simple, como en la nota (slide 23c, con piso verde de IDW):
    logf_hl = H.idw_interpolate(xipar_v, xiperp_v, logf_obs, nodes, xiperp_hl, eps=C.IDW_EPS)
    # (espectro) los nodos de Laguerre llegan a xi_perp~15 (lejos de los datos): se
    # atenúa el campo lejano no medido para no inyectar un suelo de ruido artificial.
    logf_hl_spec = H.idw_interpolate(xipar_v, xiperp_v, logf_obs, nodes, xiperp_hl,
                                     eps=C.IDW_EPS, taper_d0=1.5)
    c_hl = H.hl_coefficients(10.0 ** logf_hl_spec, nodes, weights, l_mu, l_w, C.M_ORDER, C.N_ORDER)
    outA = os.path.join(C.FIG_DIR, "bonus_slide23_abc.png")
    plots.plot_slide23_abc(xipar_v, xiperp_v, logf_obs, nodes, logf_grid, xiperp_hl, logf_hl, outA)
    outB = os.path.join(C.FIG_DIR, "bonus_slide25_hh_hl.png")
    plots.plot_slide25_spectra(c, c_hl, outB)
    print(f"[Bonus  ] {outA}")
    print(f"          {outB}  (HL c_00={c_hl[0,0]:.3e})")

    # ---- guarda productos procesados ----
    np.savez(os.path.join(C.DATA_PROC, "vdf_processed.npz"),
             v_par=vpar_kms, v_perp1=vp1_kms, v_perp2=vp2_kms, psd=f_v,
             xi_par=xipar_v, xi_perp=xiperp_v, nodes=nodes, weights=weights,
             logf_grid=logf_grid, c_mn=c, c_mn_hl=c_hl, logf_hl=logf_hl,
             xiperp_hl=xiperp_hl, B=Braw, b_hat=b,
             u=mom["u"], w_par=mom["w_par"], w_perp=mom["w_perp"],
             n=mom["n"], T_par=mom["T_par"], T_perp=mom["T_perp"])
    meta = dict(time=str(d["t_used"]), sign=sign, n_cm3=float(n_cc),
                u_kms=list(np.round(mom["u"]/1e3, 2)), u_par_kms=float(mom["u_par"]/1e3),
                w_par_kms=float(mom["w_par"]/1e3), w_perp_kms=float(mom["w_perp"]/1e3),
                Tpar_K=float(mom["T_par"]), Tperp_K=float(mom["T_perp"]),
                B_nT=list(np.round(Braw, 3)), vdf_units=d["vdf_units"])
    with open(os.path.join(C.DATA_PROC, "moments.json"), "w") as fh:
        json.dump(meta, fh, indent=2)
    print("\n[ok] productos guardados en data/processed/")
    return meta


if __name__ == "__main__":
    main()
