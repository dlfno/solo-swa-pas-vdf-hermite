"""
Figuras del proyecto, con el estilo de la nota de B. Park:
  - Tarea 1: dispersión de la VDF en marco FA (slide 3, sin panel de energía)
  - Tarea 2: VDF interpolada en la malla de nodos de Hermite (slide 23, panel a y b)
  - Tarea 3: espectro 2D de Hermite log(c_mn^2/c_00^2) (slide 25, panel izquierdo)
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from config import PSD_CLIM, SPEC_CLIM, XI_LIM

plt.rcParams.update({
    "font.size": 11, "axes.linewidth": 1.0, "figure.dpi": 130,
    "mathtext.fontset": "cm", "font.family": "serif",
})
CMAP = "jet"


def _cell_edges(centers, first=None):
    """Bordes de celda (punto medio) para pcolormesh a partir de centros."""
    c = np.asarray(centers, float)
    e = np.empty(c.size + 1)
    e[1:-1] = 0.5 * (c[:-1] + c[1:])
    e[0]  = c[0]  - (c[1] - c[0]) / 2 if first is None else first
    e[-1] = c[-1] + (c[-1] - c[-2]) / 2
    return e


# ===========================================================================
# TAREA 1 — dispersión de la VDF (slide 3, sin panel de energía)
# ===========================================================================
def plot_scatter_vdf(vpar, vp1, vp2, f, time_str, outpath):
    """vpar,vp1,vp2 en km/s ; f en s^3/m^6 (>0)."""
    good = np.isfinite(f) & (f > 0)
    vpar, vp1, vp2, f = vpar[good], vp1[good], vp2[good], f[good]
    logf = np.log10(f)
    order = np.argsort(logf)            # dibuja primero los débiles, encima los fuertes
    vpar, vp1, vp2, logf = vpar[order], vp1[order], vp2[order], logf[order]

    fig = plt.figure(figsize=(15.0, 5.2))
    # columnas: [dispersión | barra de color | ESPACIADOR | v_par | v_perp1 | v_perp2]
    gs = GridSpec(2, 6, figure=fig, width_ratios=[1.3, 0.05, 0.42, 1, 1, 1],
                  hspace=0.08, wspace=0.42, left=0.055, right=0.985, top=0.90, bottom=0.13)

    # ---- panel izquierdo: dispersión coloreada por PSD ----
    ax_top = fig.add_subplot(gs[0, 0])
    ax_bot = fig.add_subplot(gs[1, 0], sharex=ax_top)
    for ax, vp, lbl in ((ax_top, vp1, r"$v_{\perp,1}\ [\mathrm{km\,s^{-1}}]$"),
                        (ax_bot, vp2, r"$v_{\perp,2}\ [\mathrm{km\,s^{-1}}]$")):
        sc = ax.scatter(vpar, vp, c=logf, cmap=CMAP, vmin=PSD_CLIM[0], vmax=PSD_CLIM[1],
                        s=26, alpha=0.65, edgecolors="none")
        ax.axhline(0, ls="--", color="0.4", lw=0.8)
        ax.set_ylabel(lbl)
        ax.grid(alpha=0.15)
    ax_bot.set_xlabel(r"$v_\parallel\ [\mathrm{km\,s^{-1}}]$")
    plt.setp(ax_top.get_xticklabels(), visible=False)
    cax = fig.add_subplot(gs[:, 1])
    cb = fig.colorbar(sc, cax=cax)
    cb.set_label(r"$\log(\mathrm{s^3\,m^{-6}})$", labelpad=2)

    # ---- paneles derechos: PSD vs cada componente ----
    titles = (r"$v_\parallel$", r"$v_{\perp,1}$", r"$v_{\perp,2}$")
    xs     = (vpar, vp1, vp2)
    xlabs  = (r"$v_\parallel\ [\mathrm{km\,s^{-1}}]$",
              r"$v_{\perp,1}\ [\mathrm{km\,s^{-1}}]$",
              r"$v_{\perp,2}\ [\mathrm{km\,s^{-1}}]$")
    f_lin = 10.0**logf
    for k in range(3):
        ax = fig.add_subplot(gs[:, 3 + k])
        ax.scatter(xs[k], f_lin, s=7, c="0.25", alpha=0.5, edgecolors="none")
        ax.set_yscale("log"); ax.set_ylim(1e-11, 1e-5)
        ax.set_title(titles[k]); ax.set_xlabel(xlabs[k]); ax.grid(alpha=0.2)
        if k == 0:
            ax.set_ylabel(r"PSD $[\mathrm{s^3\,m^{-6}}]$")
        else:
            plt.setp(ax.get_yticklabels(), visible=False)

    fig.suptitle(f"SolO PAS-VDF {time_str}", fontsize=14, y=0.985)
    fig.savefig(outpath, bbox_inches="tight")
    plt.close(fig)
    return outpath


# ===========================================================================
# TAREA 2 — VDF observada (a) e interpolada en la malla HH (b)  [slide 23]
# ===========================================================================
def plot_hermite_grid(xi_par_obs, xi_perp_obs, logf_obs,
                      nodes, logf_grid, outpath, mirror=True):
    """Panel a: dispersión observada (xi_par, xi_perp>=0).
       Panel b: malla H60xH60 interpolada (xi_perp reflejado, estilo Larosa)."""
    fig, (axa, axb) = plt.subplots(1, 2, figsize=(11, 4.6))

    # --- panel a: observado ---
    o = np.argsort(logf_obs)
    sc = axa.scatter(xi_par_obs[o], xi_perp_obs[o], c=logf_obs[o], cmap=CMAP,
                     vmin=PSD_CLIM[0], vmax=PSD_CLIM[1], s=22, alpha=0.7, edgecolors="none")
    axa.set_title("Observed VDF (a)")
    axa.set_xlabel(r"$\xi_\parallel$"); axa.set_ylabel(r"$\xi_\perp$")
    axa.set_xlim(-XI_LIM, XI_LIM); axa.set_ylim(0, XI_LIM)
    cb = fig.colorbar(sc, ax=axa); cb.set_label(r"$\log(\mathrm{s^3\,m^{-6}})$")

    # --- panel b: malla HH interpolada ---
    # nodos del eje perp (reflejados a ambos lados si mirror)
    xi_par = nodes
    xi_per = nodes
    # bordes de celda para pcolormesh
    def edges(c):
        c = np.sort(c); e = np.empty(c.size+1)
        e[1:-1] = 0.5*(c[:-1]+c[1:]); e[0] = c[0]-(c[1]-c[0])/2; e[-1] = c[-1]+(c[-1]-c[-2])/2
        return e
    ep, eperp = edges(xi_par), edges(xi_per)
    pm = axb.pcolormesh(ep, eperp, logf_grid.T, cmap=CMAP,
                        vmin=PSD_CLIM[0], vmax=PSD_CLIM[1], shading="auto")
    axb.set_title(r"HH quadrature $H_{60}\times H_{60}$ (b)")
    axb.set_xlabel(r"$\xi_\parallel$"); axb.set_ylabel(r"$\xi_\perp$")
    axb.set_xlim(-XI_LIM, XI_LIM); axb.set_ylim(-XI_LIM, XI_LIM)
    cb = fig.colorbar(pm, ax=axb); cb.set_label(r"$\log(\mathrm{s^3\,m^{-6}})$")

    fig.tight_layout()
    fig.savefig(outpath, bbox_inches="tight")
    plt.close(fig)
    return outpath


# ===========================================================================
# TAREA 3 — espectro 2D de Hermite (slide 25, panel izquierdo HH)
# ===========================================================================
# ===========================================================================
# BONUS — reproducción completa del slide 23 (a,b,c) y slide 25 (HH | HL)
# ===========================================================================
def plot_slide23_abc(xi_par_obs, xi_perp_obs, logf_obs,
                     h_nodes, logf_hh, xiperp_hl, logf_hl, outpath):
    """3 paneles como el slide 23: (a) observado, (b) HH (H60xH60, reflejado),
       (c) HL (H60xL60, cilíndrico mu=xi_perp^2)."""
    fig, axs = plt.subplots(1, 3, figsize=(16, 4.7))
    o = np.argsort(logf_obs)
    sc = axs[0].scatter(xi_par_obs[o], xi_perp_obs[o], c=logf_obs[o], cmap=CMAP,
                        vmin=PSD_CLIM[0], vmax=PSD_CLIM[1], s=20, alpha=0.7, edgecolors="none")
    axs[0].set_title("(a) Observed VDF — SWA-PAS"); axs[0].set_ylim(0, XI_LIM)
    eh = _cell_edges(h_nodes)
    axs[1].pcolormesh(eh, eh, logf_hh.T, cmap=CMAP, vmin=PSD_CLIM[0], vmax=PSD_CLIM[1], shading="auto")
    axs[1].set_title(r"(b) HH $H_{60}\times H_{60}$"); axs[1].set_ylim(-XI_LIM, XI_LIM)
    ep = _cell_edges(xiperp_hl, first=0.0)
    pm = axs[2].pcolormesh(eh, ep, logf_hl.T, cmap=CMAP, vmin=PSD_CLIM[0], vmax=PSD_CLIM[1], shading="auto")
    axs[2].set_title(r"(c) HL $H_{60}\times L_{60}$ (cyl.)"); axs[2].set_ylim(0, XI_LIM)
    for ax in axs:
        ax.set_xlabel(r"$\xi_\parallel$"); ax.set_ylabel(r"$\xi_\perp$"); ax.set_xlim(-XI_LIM, XI_LIM)
    cb = fig.colorbar(pm, ax=axs, fraction=0.025, pad=0.01); cb.set_label(r"$\log(\mathrm{s^3\,m^{-6}})$")
    fig.savefig(outpath, bbox_inches="tight"); plt.close(fig)
    return outpath


def plot_slide25_spectra(c_hh, c_hl, outpath):
    """2 paneles como el slide 25: espectro HH (izq) y HL (der)."""
    fig, axs = plt.subplots(1, 2, figsize=(12, 5.2))
    for ax, c, ttl in ((axs[0], c_hh, "HH quadrature"), (axs[1], c_hl, "HL quadrature")):
        M, N = c.shape
        with np.errstate(divide="ignore"):
            logr = np.log10(np.where(c**2 > 0, (c**2)/(c[0, 0]**2), np.nan))
        pm = ax.pcolormesh(np.arange(-0.5, M, 1.0), np.arange(-0.5, N, 1.0), logr.T,
                           cmap=CMAP, vmin=SPEC_CLIM[0], vmax=SPEC_CLIM[1], shading="auto")
        ax.set_title(ttl); ax.set_xlabel(r"$m$"); ax.set_ylabel(r"$n$")
        ax.set_xlim(0, M-1); ax.set_ylim(0, N-1); ax.set_aspect("equal")
        cb = fig.colorbar(pm, ax=ax); cb.set_label(r"$\log\,(c_{mn}^2/c_{00}^2)$")
    fig.tight_layout(); fig.savefig(outpath, bbox_inches="tight"); plt.close(fig)
    return outpath


def plot_hh_spectrum(c, outpath):
    """c : (M+1, N+1) coeficientes.  Grafica log10(c_mn^2 / c_00^2)."""
    M, N = c.shape
    ratio = (c**2) / (c[0, 0]**2)
    with np.errstate(divide="ignore"):
        logr = np.log10(np.where(ratio > 0, ratio, np.nan))
    fig, ax = plt.subplots(figsize=(6.4, 5.4))
    m_edges = np.arange(-0.5, M, 1.0)
    n_edges = np.arange(-0.5, N, 1.0)
    pm = ax.pcolormesh(m_edges, n_edges, logr.T, cmap=CMAP,
                       vmin=SPEC_CLIM[0], vmax=SPEC_CLIM[1], shading="auto")
    ax.set_title("HH quadrature")
    ax.set_xlabel(r"$m$"); ax.set_ylabel(r"$n$")
    ax.set_xlim(0, M-1); ax.set_ylim(0, N-1)
    ax.set_aspect("equal")
    cb = fig.colorbar(pm, ax=ax); cb.set_label(r"$\log\,(c_{mn}^2/c_{00}^2)$")
    fig.tight_layout()
    fig.savefig(outpath, bbox_inches="tight")
    plt.close(fig)
    return outpath


def plot_spectrum_cuts(c, outpath):
    """Cortes paralelo c_m0^2/c_00^2 (azul) y perpendicular c_0n^2/c_00^2 (rojo)
       en escala log-log (estilo slide 26).  En HH el perpendicular sólo tiene
       potencia en órdenes pares (mirroring f(xi_perp)=f(-xi_perp))."""
    M, N = c.shape
    par = (c[:, 0]**2) / (c[0, 0]**2)
    per = (c[0, :]**2) / (c[0, 0]**2)
    m = np.arange(M); n = np.arange(N)
    n_even = n[n % 2 == 0]                     # sólo órdenes pares (slide 26)
    fig, ax = plt.subplots(figsize=(6.4, 4.8))
    ax.loglog(m + 1, par, "o-", color="navy", ms=4, lw=1.2,
              label=r"$c_{m0}^2/c_{00}^2$ (parallel)")
    ax.loglog(n_even + 1, per[n_even], "o-", color="crimson", ms=4, lw=1.2,
              label=r"$c_{0n}^2/c_{00}^2$ (perp, even $n$)")
    # guías de ley de potencia indicativas (ajuste log-log en el rango inercial)
    for arr, idx, col, sl in ((par, m, "navy", None), (per, n_even, "crimson", None)):
        sel = (idx >= 1) & (idx <= 9) & (arr[idx] > 0)
        if sel.sum() >= 2:
            p = np.polyfit(np.log10(idx[sel] + 1), np.log10(arr[idx][sel]), 1)
            xx = np.array([1.3, 12.0])
            ax.plot(xx, 10**(p[1]) * xx**p[0], "--", color=col, lw=1.0, alpha=0.7,
                    label=fr"$\propto k^{{{p[0]:.1f}}}$")
    ax.set_xlabel(r"$m+1,\ n+1$"); ax.set_ylabel(r"$c^2/c_{00}^2$")
    ax.set_ylim(1e-6, 2); ax.legend(fontsize=9); ax.grid(alpha=0.3, which="both")
    fig.tight_layout(); fig.savefig(outpath, bbox_inches="tight"); plt.close(fig)
    return outpath
