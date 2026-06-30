"""
Física de la reconstrucción de la VDF.

  - velocidad a partir de (energía, azimut, elevación)  ->  marco instrumento (xyz)
  - rotación PAS -> RTN
  - marco alineado al campo (field-aligned, FA):  v_par, v_perp1, v_perp2
  - momentos del plasma (n, u, T_par, T_perp, w_par, w_perp) con d^3v correcto
  - coordenadas normalizadas xi_par, xi_perp (velocidad peculiar en el marco de reposo)

Convenciones (nota de B. Park, slides 1-2):
  |v| = sqrt(2 q E / m)
  vx = |v| cos(el) cos(az)
  vy = |v| cos(el) sin(az)
  vz = |v| sin(el)
El instrumento mide la dirección de LLEGADA de las partículas; la velocidad
física lleva signo `sign` (=-1 en la convención de SWA, =+1 en la nota).
"""
import numpy as np
from config import Q_E, M_P, K_B


# ---------------------------------------------------------------------------
# Velocidad desde energía + ángulos
# ---------------------------------------------------------------------------
def speed_from_energy(E_eV):
    """Rapidez |v| [m/s] de un protón con energía E [eV]:  |v| = sqrt(2 q E / m)."""
    return np.sqrt(2.0 * Q_E * np.asarray(E_eV, float) / M_P)


def instrument_velocity(E_eV, azim_deg, elev_deg, signs=(-1.0, -1.0, +1.0)):
    """Vector velocidad en el marco PAS (xyz).  Acepta broadcasting.

    El instrumento mide la dirección de LLEGADA; la velocidad física es opuesta.
    Convención SWA / Guía de usuario PAS (Fedorov 2020, asimétrica):
        Vx = -cos(El) cos(Az)
        Vy = -cos(El) sin(Az)
        Vz = +sin(El)
    -> signs = (-1, -1, +1).  La convención de la nota (todo positivo) es signs=(+1,+1,+1).
    Con PAS_to_RTN da V_R > 0 (viento solar anti-solar).  Devuelve vx, vy, vz [m/s]."""
    v  = speed_from_energy(E_eV)
    az = np.deg2rad(np.asarray(azim_deg, float))
    el = np.deg2rad(np.asarray(elev_deg, float))
    sx, sy, sz = signs
    vx = sx * v * np.cos(el) * np.cos(az)
    vy = sy * v * np.cos(el) * np.sin(az)
    vz = sz * v * np.sin(el)
    return vx, vy, vz


def rotate_pas_to_rtn(vx, vy, vz, M):
    """Aplica la matriz 3x3 PAS_to_RTN a cada vector (vx,vy,vz).
    v_rtn[a] = sum_b M[a,b] v_xyz[b].  Devuelve vR, vT, vN (misma forma de entrada)."""
    vR = M[0, 0]*vx + M[0, 1]*vy + M[0, 2]*vz
    vT = M[1, 0]*vx + M[1, 1]*vy + M[1, 2]*vz
    vN = M[2, 0]*vx + M[2, 1]*vy + M[2, 2]*vz
    return vR, vT, vN


# ---------------------------------------------------------------------------
# Marco alineado al campo magnético
# ---------------------------------------------------------------------------
def field_aligned_basis(B):
    """Base ortonormal (b_hat, e_perp1, e_perp2) a partir de B (RTN).
    e_perp1 = b x R_hat (o T_hat si b || R_hat);  e_perp2 = b x e_perp1."""
    B = np.asarray(B, float)
    b = B / np.linalg.norm(B)
    ref = np.array([1.0, 0.0, 0.0])               # R_hat
    if abs(np.dot(b, ref)) > 0.99:
        ref = np.array([0.0, 1.0, 0.0])           # T_hat
    e1 = np.cross(b, ref); e1 /= np.linalg.norm(e1)
    e2 = np.cross(b, e1);  e2 /= np.linalg.norm(e2)
    return b, e1, e2


def project_field_aligned(vR, vT, vN, b, e1, e2):
    """Proyecta cada vector RTN sobre (b_hat, e_perp1, e_perp2)."""
    v_par   = vR*b[0]  + vT*b[1]  + vN*b[2]
    v_perp1 = vR*e1[0] + vT*e1[1] + vN*e1[2]
    v_perp2 = vR*e2[0] + vT*e2[1] + vN*e2[2]
    return v_par, v_perp1, v_perp2


# ---------------------------------------------------------------------------
# Momentos del plasma (integración con volumen de velocidad correcto)
# ---------------------------------------------------------------------------
def velocity_volume(E_eV, dE_eV, elev_deg, d_elev_deg, d_azim_deg):
    """Elemento de volumen de espacio de velocidades de cada bin:
        d^3v = |v|^2 d|v| dOmega ,  dOmega = cos(el) d(el) d(az)
        d|v| = sqrt(q/(2 m E)) dE
    Todas las entradas se difunden a la forma del bin.  Devuelve dV [m^3/s^3]."""
    v   = speed_from_energy(E_eV)
    dv  = np.sqrt(Q_E / (2.0 * M_P * np.asarray(E_eV, float))) * np.asarray(dE_eV, float)
    el  = np.deg2rad(np.asarray(elev_deg, float))
    dOm = np.cos(el) * np.deg2rad(np.asarray(d_elev_deg, float)) * np.deg2rad(np.asarray(d_azim_deg, float))
    return v*v * dv * dOm


def compute_moments(vR, vT, vN, f, dV, b):
    """Momentos del plasma a partir de la VDF discreta f [s^3/m^6] con volúmenes dV.
       n = sum f dV
       u = (1/n) sum v f dV
       P_ab = m sum (v-u)_a (v-u)_b f dV ;  T = P/(n k_B)
       T_par = b.T.b ;  T_perp = (tr T - T_par)/2
       w = sqrt(2 k_B T / m)
    Devuelve un diccionario con todo (unidades SI; velocidades en m/s)."""
    f = np.asarray(f, float); dV = np.asarray(dV, float)
    w = f * dV                                   # peso = densidad por bin
    n = np.sum(w)
    uR = np.sum(vR*w)/n; uT = np.sum(vT*w)/n; uN = np.sum(vN*w)/n
    u = np.array([uR, uT, uN])
    cR, cT, cN = vR-uR, vT-uT, vN-uN             # velocidad peculiar
    comp = [cR, cT, cN]
    P = np.empty((3, 3))
    for a in range(3):
        for bb in range(3):
            P[a, bb] = M_P * np.sum(comp[a]*comp[bb]*w)
    T = P / (n * K_B)                            # tensor de temperatura [K]
    T_par  = b @ T @ b
    T_perp = 0.5 * (np.trace(T) - T_par)
    w_par  = np.sqrt(2.0 * K_B * T_par  / M_P)
    w_perp = np.sqrt(2.0 * K_B * T_perp / M_P)
    u_par = u @ b
    return dict(n=n, u=u, u_par=u_par, T_par=T_par, T_perp=T_perp,
                w_par=w_par, w_perp=w_perp, P=P, T=T)


def normalized_coords(vR, vT, vN, u, b, w_par, w_perp):
    """Coordenadas normalizadas (xi_par, xi_perp) en el marco de reposo del plasma.
        c = v - u  (velocidad peculiar)
        c_par = c . b_hat ;  c_perp = |c - c_par b_hat|
        xi_par = c_par / w_par ;  xi_perp = c_perp / w_perp   (>= 0)
    """
    cR, cT, cN = vR-u[0], vT-u[1], vN-u[2]
    c_par = cR*b[0] + cT*b[1] + cN*b[2]
    c2    = cR*cR + cT*cT + cN*cN
    c_perp = np.sqrt(np.maximum(c2 - c_par*c_par, 0.0))
    return c_par/w_par, c_perp/w_perp
