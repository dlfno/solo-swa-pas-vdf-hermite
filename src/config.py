"""
Configuración global del proyecto: constantes físicas, rutas y parámetros.

Reproducción de la Función de Distribución de Velocidad (VDF) de protones
observada por Solar Orbiter / SWA-PAS el 2022-03-08 14:45:22, siguiendo la
nota de cátedra de B. Park (Instituto de Geofísica, UNAM).
"""
import os

# ---------------------------------------------------------------------------
# Constantes físicas (SI, CODATA 2018)
# ---------------------------------------------------------------------------
Q_E      = 1.602176634e-19      # carga elemental [C]
M_P      = 1.67262192369e-27    # masa del protón [kg]
K_B      = 1.380649e-23         # constante de Boltzmann [J/K]
EV_TO_J  = Q_E                  # 1 eV en Julios

# ---------------------------------------------------------------------------
# Rutas
# ---------------------------------------------------------------------------
HERE      = os.path.dirname(os.path.abspath(__file__))
ROOT      = os.path.dirname(HERE)
DATA_RAW  = os.path.join(ROOT, "data", "raw")
DATA_PROC = os.path.join(ROOT, "data", "processed")
FIG_DIR   = os.path.join(ROOT, "figures")
DOC_FIG   = os.path.join(ROOT, "docs", "figures")

VDF_CDF = os.path.join(DATA_RAW, "solo_L2_swa-pas-vdf_20220308_V02.cdf")
MAG_CDF = os.path.join(DATA_RAW, "solo_L2_mag-rtn-normal_20220308_V02.cdf")

# ---------------------------------------------------------------------------
# Parámetros del análisis
# ---------------------------------------------------------------------------
TARGET_TIME = "2022-03-08T14:45:22"   # instante a reproducir (slide 23)

# Descomposición de Hermite-Hermite (HH), estilo Larosa et al. (2025)
N_NODES = 60     # orden de la cuadratura de Gauss-Hermite (H_60)  -> nodos
M_ORDER = 60     # orden máximo paralelo  (m = 0..60)
N_ORDER = 60     # orden máximo perpendicular (n = 0..60)
IDW_EPS = 1e-3   # regularización del peso inverso a la distancia

# Escalas de color
PSD_CLIM  = (-11, -7)   # log10 PSD [s^3 m^-6]
SPEC_CLIM = (-8, 0)     # log10 (c_mn^2 / c_00^2)
XI_LIM    = 5.5         # rango de ejes en coordenadas normalizadas xi

for _d in (DATA_PROC, FIG_DIR, DOC_FIG):
    os.makedirs(_d, exist_ok=True)
