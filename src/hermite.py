"""
Núcleo numérico de la descomposición de Hermite-Hermite (HH).

Implementa:
  - nodos y pesos de la cuadratura de Gauss-Hermite (algoritmo de Golub-Welsch,
    via numpy.polynomial.hermite.hermgauss -> autovalores de la matriz de Jacobi).
  - funciones de Hermite ORTONORMALES  psi_m(x) = (2^m m! sqrt(pi))^{-1/2} e^{-x^2/2} H_m(x)
    evaluadas con la recurrencia estable (sin overflow hasta m ~ cientos).
  - interpolación por peso inverso a la distancia (IDW) de la VDF dispersa
    sobre la malla de nodos.
  - coeficientes c_mn = sum_ij f(xi_i,xi_j) e^{xi_i^2} w_i e^{xi_j^2} w_j psi_m(xi_i) psi_n(xi_j).

Referencias: Golub & Welsch (1969); Servidio et al. (2017); Larosa et al. (2025).
"""
import numpy as np


# ---------------------------------------------------------------------------
# 1) Nodos y pesos de Gauss-Hermite  (Golub-Welsch)
# ---------------------------------------------------------------------------
def gauss_hermite(n):
    """Nodos x_i (raíces de H_n) y pesos w_i tales que
    \\int e^{-x^2} g(x) dx ~ sum_i w_i g(x_i).  Exacto para polinomios de grado <= 2n-1.

    numpy.polynomial.hermite.hermgauss usa la matriz companion (equivalente a
    Golub-Welsch): no requiere semillas iniciales ni deja raíces sin encontrar.
    """
    x, w = np.polynomial.hermite.hermgauss(n)
    return x, w


# ---------------------------------------------------------------------------
# 2) Funciones de Hermite ortonormales  psi_m(x)   (recurrencia estable)
# ---------------------------------------------------------------------------
def hermite_functions(x, mmax):
    """Devuelve un array de forma (len(x), mmax+1) con psi_0..psi_mmax.

    Recurrencia de las funciones de Hermite (acotadas, |psi_m| ~ pi^-1/4):
        psi_0 = pi^{-1/4} e^{-x^2/2}
        psi_1 = sqrt(2) x psi_0
        psi_{m+1} = sqrt(2/(m+1)) x psi_m - sqrt(m/(m+1)) psi_{m-1}
    """
    x = np.asarray(x, dtype=float)
    psi = np.empty((x.size, mmax + 1), dtype=float)
    psi[:, 0] = np.pi ** -0.25 * np.exp(-0.5 * x * x)
    if mmax >= 1:
        psi[:, 1] = np.sqrt(2.0) * x * psi[:, 0]
    for m in range(1, mmax):
        psi[:, m + 1] = (np.sqrt(2.0 / (m + 1)) * x * psi[:, m]
                         - np.sqrt(m / (m + 1.0)) * psi[:, m - 1])
    return psi


# ---------------------------------------------------------------------------
# 3) Interpolación por peso inverso a la distancia (IDW)
# ---------------------------------------------------------------------------
def idw_interpolate(xs, ys, vals, grid_x, grid_y, eps=1e-3, power=2,
                    taper_d0=None, floor=None):
    """Interpola valores `vals` definidos en puntos dispersos (xs, ys)
    sobre la malla cartesiana (grid_x[i], grid_y[j]).

        f(i,j) = sum_k w_k vals_k / sum_k w_k ,  w_k = 1 / (r_k^power + eps)
        r_k^2 = (xs_k - grid_x_i)^2 + (ys_k - grid_y_j)^2

    Si `taper_d0` se especifica, los nodos lejanos de TODO dato (distancia al
    vecino más cercano d_nn) se atenúan hacia `floor`:
        f <- f * t + floor * (1-t),   t = exp(-(d_nn/taper_d0)^2)
    (la VDF decae donde no se midió).  floor por defecto = min(vals).

    Devuelve un array (len(grid_x), len(grid_y))  [eje 0 = xi_par, eje 1 = xi_perp].
    """
    xs = np.asarray(xs); ys = np.asarray(ys); vals = np.asarray(vals)
    if floor is None:
        floor = float(np.min(vals))
    GX, GY = np.meshgrid(grid_x, grid_y, indexing="ij")     # (Ni, Nj)
    flatx, flaty = GX.ravel(), GY.ravel()
    res = np.empty(flatx.size, dtype=float)
    dnn = np.empty(flatx.size, dtype=float)
    block = 4096                                            # por bloques (memoria)
    for s in range(0, flatx.size, block):
        e = min(s + block, flatx.size)
        dx = flatx[s:e, None] - xs[None, :]
        dy = flaty[s:e, None] - ys[None, :]
        r2 = dx * dx + dy * dy
        w = 1.0 / (r2 ** (power / 2.0) + eps)
        res[s:e] = (w @ vals) / w.sum(axis=1)
        dnn[s:e] = np.sqrt(r2.min(axis=1))
    if taper_d0 is not None:
        t = np.exp(-(dnn / taper_d0) ** 2)
        res = res * t + floor * (1.0 - t)
    return res.reshape(GX.shape)


# ---------------------------------------------------------------------------
# 4) Coeficientes de la descomposición Hermite-Hermite
# ---------------------------------------------------------------------------
def hh_coefficients(f_grid, nodes, weights, mmax, nmax):
    """c_mn = sum_{i,j} f(xi_i,xi_j) [e^{xi_i^2} w_i] [e^{xi_j^2} w_j] psi_m(xi_i) psi_n(xi_j)

    f_grid : (Ni, Nj)  VDF interpolada en los nodos  (eje0 = paralelo, eje1 = perp)
    nodes, weights : nodos/pesos de Gauss-Hermite (mismos en ambos ejes)
    Devuelve c con forma (mmax+1, nmax+1).
    """
    # peso modificado  W_i = w_i e^{xi_i^2}  (estable en float64 hasta n~cientos)
    W = weights * np.exp(nodes * nodes)
    psi_par = hermite_functions(nodes, mmax)     # (Ni, mmax+1)
    psi_per = hermite_functions(nodes, nmax)     # (Nj, nmax+1)
    # c = psi_par^T  (W_i f_ij W_j)  psi_per
    A = (W[:, None] * f_grid) * W[None, :]       # (Ni, Nj)
    c = psi_par.T @ A @ psi_per                  # (mmax+1, nmax+1)
    return c


# ---------------------------------------------------------------------------
# 5) Hermite-Laguerre (HL):  Gauss-Laguerre + funciones de Laguerre Gamma_n^0
# ---------------------------------------------------------------------------
def gauss_laguerre(n):
    """Nodos mu_i (raíces de L_n) y pesos w_i:  \\int_0^inf e^{-x} g(x) dx ~ sum w_i g(mu_i)."""
    mu, w = np.polynomial.laguerre.laggauss(n)
    return mu, w


def laguerre_functions_k0(mu, nmax):
    """Funciones de Laguerre ortonormales para k=0:
        Gamma_n^0(mu) = e^{-mu/2} L_n(mu)        (con L_n estándar, normalizado)
    Recurrencia estable (acotada):
        phi_0 = e^{-mu/2},  phi_1 = e^{-mu/2}(1-mu)
        phi_{n+1} = ((2n+1-mu) phi_n - n phi_{n-1}) / (n+1)
    Devuelve array (len(mu), nmax+1)."""
    mu = np.asarray(mu, float)
    g = np.empty((mu.size, nmax + 1), dtype=float)
    e = np.exp(-0.5 * mu)
    g[:, 0] = e
    if nmax >= 1:
        g[:, 1] = e * (1.0 - mu)
    for n in range(1, nmax):
        g[:, n + 1] = ((2*n + 1 - mu) * g[:, n] - n * g[:, n - 1]) / (n + 1)
    return g


def hl_coefficients(f_grid, h_nodes, h_weights, l_nodes, l_weights, mmax, nmax):
    """Coeficientes Hermite-Laguerre:
       c_mn = sum_ij f(xi_i, sqrt(mu_j)) [e^{xi_i^2} w_H(xi_i)] [e^{mu_j} w_L(mu_j)]
                     psi_m(xi_i) Gamma_n^0(mu_j)
       f_grid : (Ni_par, Nj_perp)  con xi_perp_j = sqrt(mu_j).
       Devuelve c (mmax+1, nmax+1)."""
    WH = h_weights * np.exp(h_nodes * h_nodes)
    WL = l_weights * np.exp(l_nodes)
    psi = hermite_functions(h_nodes, mmax)       # (Ni, mmax+1)
    Gam = laguerre_functions_k0(l_nodes, nmax)   # (Nj, nmax+1)
    A = (WH[:, None] * f_grid) * WL[None, :]
    return psi.T @ A @ Gam


# ---------------------------------------------------------------------------
# Auto-test: ortonormalidad discreta y exactitud de la cuadratura
# ---------------------------------------------------------------------------
if __name__ == "__main__":
    n = 60
    x, w = gauss_hermite(n)
    print(f"Gauss-Hermite n={n}: {len(x)} nodos en [{x.min():.3f}, {x.max():.3f}]")

    # (a) cuadratura exacta:  \int e^{-x^2} x^2 dx = sqrt(pi)/2
    approx = np.sum(w * x**2)
    print(f"  \\int e^-x^2 x^2 dx : {approx:.12f}  vs  {np.sqrt(np.pi)/2:.12f}")

    # (b) ortonormalidad discreta de las funciones de Hermite:
    #     sum_i (w_i e^{x_i^2}) psi_m(x_i) psi_n(x_i) = delta_mn  (exacto para m+n <= 2n-1)
    W = w * np.exp(x*x)
    M = 40
    psi = hermite_functions(x, M)
    G = psi.T @ (W[:, None] * psi)          # (M+1, M+1) ~ Identidad
    off = np.abs(G - np.eye(M+1))
    print(f"  Ortonormalidad psi_m (m<= {M}):  max|G-I| = {off.max():.2e}")
    print(f"  diag(G)[:5] = {np.diag(G)[:5]}")

    # (c) las funciones de Hermite están acotadas (no hay overflow)
    psi_big = hermite_functions(x, 60)
    print(f"  max|psi_m| (m<=60) = {np.abs(psi_big).max():.4f}  (debe ser ~0.75)")

    # (d) Gauss-Laguerre + ortonormalidad de las funciones de Laguerre Gamma_n^0
    mu, wl = gauss_laguerre(n)
    print(f"\nGauss-Laguerre n={n}: {len(mu)} nodos en [{mu.min():.3f}, {mu.max():.3f}]")
    print(f"  min peso w_L = {wl.min():.2e}  (debe ser > 0, sin underflow)")
    WL = wl * np.exp(mu)
    print(f"  max W_L = w_L e^mu = {WL.max():.3e}  (finito)")
    ML = 40
    Gam = laguerre_functions_k0(mu, ML)
    GL = Gam.T @ (WL[:, None] * Gam)        # ~ Identidad
    offL = np.abs(GL - np.eye(ML+1))
    print(f"  Ortonormalidad Gamma_n^0 (n<= {ML}): max|G-I| = {offL.max():.2e}")
    Gam_big = laguerre_functions_k0(mu, 60)
    print(f"  max|Gamma_n^0| (n<=60) = {np.abs(Gam_big).max():.4f}  (acotada, sin overflow)")
