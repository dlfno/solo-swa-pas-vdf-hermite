"""
Lectura de los CDF de Solar Orbiter (SOAR):
  - MAG  L2 'mag-rtn-normal'  ->  B(RTN) en el instante objetivo
  - SWA-PAS L2 'swa-pas-vdf'  ->  VDF + rejillas (Energy, Azimuth, Elevation, PAS_to_RTN)
"""
import numpy as np
import cdflib


def _epoch_to_dt64(epoch):
    """Convierte un arreglo de épocas CDF (TT2000) a numpy.datetime64[ns]."""
    dt = cdflib.cdfepoch.to_datetime(epoch)
    return np.asarray(dt, dtype="datetime64[ns]")


def read_mag_at(path, target_iso, window_s=4.0):
    """Devuelve (B_vec[R,T,N] promedio [nT], t_used, n_samples) en una ventana
    de +/- window_s/2 alrededor de target_iso.  Si no hay muestras, usa la más cercana."""
    cdf = cdflib.CDF(path)
    epoch = cdf.varget("EPOCH")
    B = np.asarray(cdf.varget("B_RTN"), float)          # (Nt, 3)
    t = _epoch_to_dt64(epoch)
    t0 = np.datetime64(target_iso)
    half = np.timedelta64(int(window_s * 5e8), "ns")    # window_s/2 en ns
    sel = (t >= t0 - half) & (t <= t0 + half)
    # descarta relleno
    fill = ~np.all(np.isfinite(B), axis=1) | (np.abs(B).max(axis=1) > 1e6)
    sel &= ~fill
    if sel.sum() == 0:
        idx = np.argmin(np.abs(t - t0))
        return B[idx], t[idx], 1
    return B[sel].mean(axis=0), t0, int(sel.sum())


def vdf_info(path):
    """Imprime la estructura del CDF de la VDF (para inspección)."""
    cdf = cdflib.CDF(path)
    info = cdf.cdf_info()
    out = []
    for v in info.zVariables:
        vi = cdf.varinq(v)
        att = cdf.varattsget(v)
        units = att.get("UNITS", att.get("units", ""))
        out.append((v, tuple(vi.Dim_Sizes), vi.Last_Rec + 1, vi.Data_Type_Description, units))
    return out


def read_vdf_timestep(path, target_iso):
    """Lee el paso temporal más cercano a target_iso del producto swa-pas-vdf.

    Devuelve un diccionario con:
      vdf       : VDF 3D del bin temporal           forma (n_az, n_el, n_E) o la nativa
      energy    : tabla de energías [eV]            (n_E,)
      azimuth   : ángulos de azimut [deg]           (n_az,)
      elevation : ángulos de elevación [deg]        (n_el,)
      M         : matriz PAS_to_RTN 3x3
      t_used    : instante real del bin
      idx       : índice temporal
    Las dimensiones exactas se resuelven en tiempo de ejecución según el CDF.
    """
    cdf = cdflib.CDF(path)
    names = {n.lower(): n for n in cdf.cdf_info().zVariables}

    def pick(*cands):
        for c in cands:
            if c.lower() in names:
                return names[c.lower()]
        raise KeyError(f"ninguna de {cands} en {list(names.values())}")

    epoch_name = pick("Epoch", "EPOCH", "EPOCH_TIME")
    t = _epoch_to_dt64(cdf.varget(epoch_name))
    t0 = np.datetime64(target_iso)
    idx = int(np.argmin(np.abs(t - t0)))

    vdf_name = pick("vdf", "VDF")
    vdf = np.asarray(cdf.varget(vdf_name, startrec=idx, endrec=idx), float)[0]

    energy = np.asarray(cdf.varget(pick("Energy", "ENERGY")), float)
    if energy.ndim > 1:                      # a veces depende del tiempo
        energy = energy[idx] if energy.shape[0] == t.size else energy[0]
    azimuth = np.asarray(cdf.varget(pick("Azimuth", "AZIMUTH")), float)
    elevation = np.asarray(cdf.varget(pick("Elevation", "ELEVATION")), float)

    M_name = pick("PAS_to_RTN", "PAS_TO_RTN")
    M = np.asarray(cdf.varget(M_name, startrec=idx, endrec=idx), float).reshape(3, 3)

    def get(name, default=None):
        try:
            return np.asarray(cdf.varget(pick(name)), float)
        except KeyError:
            return default

    # anchos de bin (para el volumen de espacio de velocidades en los momentos)
    dpE = get("delta_p_Energy")          # half-width + (asimétrico)
    dmE = get("delta_m_Energy")          # half-width -
    dAz = get("delta_Azimuth")           # half-width azimut
    dEl = get("delta_Elevation")         # half-width elevación
    try:
        info_val = int(cdf.varget(pick("Info"), startrec=idx, endrec=idx)[0])
    except Exception:
        info_val = None

    return dict(vdf=vdf, energy=energy, azimuth=azimuth, elevation=elevation,
                M=M, t_used=t[idx], idx=idx, info=info_val,
                delta_p_E=dpE, delta_m_E=dmE, delta_Az=dAz, delta_El=dEl,
                vdf_units=cdf.varattsget(vdf_name).get("UNITS", "?"))
